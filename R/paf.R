# paf.R — Species-sensitivity-distribution PAF lookup
#
# Two public functions:
#   ssd_paf()   — full result list with optional CI
#   ssd_hc50()  — HC50 (for Concentration Addition in msPAF)
#
# Model files are fitted lazily on first request and cached in the
# per-user R cache directory (tools::R_user_dir("leachatetools", "cache")).
# Two method variants are cached separately:
#   method = "multi"   (default) — 6-distribution model average; current
#                                   best practice regardless of what ANZG used
#   method = "anzecc"            — per-analyte distribution chosen to best
#                                   replicate the original ANZG derivation
#                                   (see fit_dist column in analyte metadata)
#
# Chemistry normalisation (index-condition / bioavailability corrections):
#   ssd_paf() applies NONE of these — it evaluates the SSD at the concentration
#   you pass, as-is, for every analyte. Six analytes carry an index-condition
#   normalisation_formula in the metadata (NH3-N: pH/temperature un-ionised
#   fraction; Cu/Ni/Zn/Cd/Pb: hardness/DOC bioavailability). Those formulas are
#   applied automatically — and only — by add_amspaf(), which reads the required
#   co-analyte rows (pH, temperature, hardness, DOC, ...) from the data frame.
#
#   So: for a pipeline assessment, pass RAW measured concentrations plus the
#   co-analyte rows to add_amspaf() and let it correct once, internally. Do not
#   pre-correct. For a standalone single-sample ssd_paf() spot-check you may
#   pre-correct ammonia with correct_ammonia_ph_temp(conc, pH, temperature_C);
#   this is a convenience for the manual path only, never a pipeline step.
#
#   Water temperature is mandatory for any sample assessed for NH3-N (see
#   add_amspaf(require_temperature=)). Source it by measurement, or derive it
#   with estimate_water_temp() — optionally fed by get_silo_air_temp(), which
#   pulls daily mean air temperature from SILO for an Australian lat/lon.
#   Cr    : no external correction needed (freshwater subset already applied)
#   NO3-N : supply hardness + hardness_units; ssd_paf() blends the
#            soft/moderate/hard SSDs by hardness-class probability (see below),
#            smoothing the abrupt class boundaries. ssd_hc50()/the CA mixture
#            step still use a single representative class (.no3_class()).
#
# Probabilistic NO3-N hardness weighting (IMPLEMENTED; see .no3_weights()):
#   True hardness is assumed log-normally distributed around the measured
#   value with CV = hardness_cv (default 0.05, i.e. 5%).
#   Weights: p_soft = plnorm(30,  log(h), log(1+cv^2)^0.5)
#            p_hard = 1 - plnorm(150, log(h), log(1+cv^2)^0.5)
#            p_mod  = 1 - p_soft - p_hard
#   Expected PAF = p_soft*PAF_soft + p_mod*PAF_mod + p_hard*PAF_hard
#   This smooths the abrupt class transitions without requiring any
#   biological interpolation between the three distinct SSDs.

# ── Analyte registry ──────────────────────────────────────────────────────────

# Canonical dashboard name → safe file stem used in the .qs cache.
# NO3-N is intentionally absent: callers supply the hardness-specific variant
# (NO3-N_soft / NO3-N_mod / NO3-N_hard) or use the ssd_paf() hardness +
# hardness_units arguments which resolve the class automatically.
.SSD_NAME_MAP <- c(
  "NH3-N"                  = "NH3_N",
  "NO3-N_soft"             = "NO3_N_soft",
  "NO3-N_mod"              = "NO3_N_mod",
  "NO3-N_hard"             = "NO3_N_hard",
  "B"                      = "B",
  "Cr"                     = "Cr",
  "Cu"                     = "Cu",
  "Ni"                     = "Ni",
  "Zn"                     = "Zn",
  "Al"                     = "Al",
  # Unspeciated As → As(V) SSD. This is the conservative max-PAF ENVELOPE of
  # the two speciation SSDs, not a redox/oxygenation assumption: the As(V) SSD
  # yields a HIGHER PAF than the As(III) SSD at every environmentally realistic
  # concentration (verified 0.1–5000 µg/L; they only cross above ~10 mg/L
  # dissolved As). Although As(III) is more toxic to the single most-sensitive
  # diatom, the whole-assemblage As(V) SSD sits lower (its lower tail is pulled
  # down by corrected algal acute points), so As(V) is conservative regardless
  # of which species is actually present. Pinned by test-paf-as-envelope.R.
  "As"                     = "As_V",
  "As(V)"                  = "As_V",
  "As(III)"                = "As_III",
  "Cd"                     = "Cd",
  "Pb"                     = "Pb",
  "Mn"                     = "Mn",
  "Hg"                     = "Hg",
  "g-BHC"                  = "g_BHC",
  "Endrin"                 = "Endrin",
  "4,4-DDT"                = "4_4_DDT",
  "4.4-DDT"                = "4_4_DDT",
  "Malathion"              = "Malathion",
  "Parathion"              = "Parathion",
  "Azinphos-methyl"        = "Azinphos_Methyl",
  "Azinphos Methyl"        = "Azinphos_Methyl",   # metadata uses space; both accepted
  "Dimethoate"             = "Dimethoate",
  "Naphthalene"            = "Naphthalene",
  "Benzo(a)pyrene"         = "Benzo_a_pyrene",
  "Benzene"                = "Benzene",
  "o-Xylene"               = "o_Xylene",
  "Phenol"                 = "Phenol",
  "2-Chlorophenol"         = "2_Chlorophenol",
  "2,4,6-Trichlorophenol"  = "2_4_6_Trichlorophenol",
  "2.4.6-Trichlorophenol"  = "2_4_6_Trichlorophenol",  # metadata uses dots; both accepted
  "Pentachlorophenol"      = "Pentachlorophenol"
)

# Analytes with no SSD (insufficient data, QSAR-derived, or parked).
.SSD_NO_MODEL <- c(
  # Insufficient data / QSAR / low-reliability:
  "HCB", "Aldrin", "4,4-DDE", "Methoxychlor", "Demeton-S-methyl",
  "Phenanthrene", "Anthracene", "Fluoranthene", "Toluene", "Ethylbenzene",
  "Se",
  # Parked — large validation discrepancy; see data-raw/analyte_metadata_parked.csv:
  "Chlorpyrifos",   # ratio_95 = 3.13 — Warne exact algorithm not replicable
  "Diazinon"        # ratio_95 = 0.47 — fitting method difference
)

# ── Cache helpers ─────────────────────────────────────────────────────────────

.cache_dir <- function() {
  d <- tools::R_user_dir("leachatetools", which = "cache")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# Per-session in-memory cache (avoids repeated disk reads).
.paf_mem_cache <- new.env(parent = emptyenv())

.cache_key <- function(stem, method) paste0(stem, "__", method)

.load_or_fit <- function(analyte, stem, method, guideline_dir) {
  key <- .cache_key(stem, method)

  # 1. In-memory cache
  if (exists(key, envir = .paf_mem_cache))
    return(get(key, envir = .paf_mem_cache))

  # 2. Disk cache
  cache_file <- file.path(.cache_dir(), paste0(key, ".qs"))
  if (file.exists(cache_file)) {
    fit <- qs2::qs_read(cache_file)
    assign(key, fit, envir = .paf_mem_cache)
    return(fit)
  }

  # 3. Fit from raw data
  fit <- .fit_model(analyte, stem, method, guideline_dir)
  if (!is.null(fit)) {
    qs2::qs_save(fit, cache_file)
    assign(key, fit, envir = .paf_mem_cache)
  }
  fit
}

.fit_model <- function(analyte, stem, method, guideline_dir) {
  # Load the fitting infrastructure from ssd_fit.R
  meta_path <- system.file("extdata", "anzecc_analyte_metadata.csv",
                            package = "leachatetools")
  meta_all <- readr::read_csv(
    meta_path,
    col_types = readr::cols(
      acr                  = readr::col_double(),
      anzecc_hc_percentile = readr::col_integer(),
      trigger_divisor      = readr::col_double(),
      n_min_override       = readr::col_integer(),
      fit_dist             = readr::col_character(),
      .default             = readr::col_character()
    ),
    show_col_types = FALSE
  ) |> dplyr::mutate(ssd_available = as.logical(ssd_available))

  meta <- meta_all |> dplyr::filter(analyte == !!analyte)
  if (nrow(meta) == 0) {
    warning("No metadata for analyte '", analyte, "' \u2014 cannot fit model.")
    return(NULL)
  }

  # Resolve which distributions to fit
  dists <- if (method == "multi") {
    c("lnorm", "llogis", "burrIII3", "lgumbel", "gamma", "weibull")
  } else {
    # method == "anzecc": use per-analyte fit_dist from metadata
    .resolve_dists(meta$fit_dist)
  }

  # Load data and fit — delegates to the fitting infrastructure
  fit <- .fit_for_analyte(analyte, stem, meta, dists, guideline_dir)
  fit
}

.resolve_dists <- function(fd_raw) {
  dist_sets <- list(
    lnorm    = "lnorm",
    lgumbel  = "lgumbel",
    burrIII3 = "burrIII3",
    multi    = c("lnorm", "llogis", "burrIII3", "lgumbel", "gamma", "weibull")
  )
  if (is.na(fd_raw) || !nzchar(trimws(fd_raw))) return("lnorm")
  fd <- trimws(fd_raw)
  if (!is.null(dist_sets[[fd]])) dist_sets[[fd]]
  else trimws(strsplit(fd, ",")[[1]])
}

# NOTE: .fit_for_analyte() is defined in ssd_fit.R, which contains
# the full data-loading logic (XLSX readers, Warne 2000 CSV, etc.).
# It requires guideline_dir to point to the "guideline data/" folder
# containing the ANZG XLSX technical brief data files.
# Set via: options(leachatetools.guideline_dir = "/path/to/guideline data")

# ── NO3-N hardness class resolution ───────────────────────────────────────────

#' Resolve the appropriate NO3-N hardness-class analyte name.
#'
#' @param hardness_mg_L Numeric. Measured hardness in mg/L as CaCO3.
#' @param hardness_cv   Numeric. Accepted for signature parity with
#'   [ssd_paf()]; not used here (this helper applies hard cutoffs). The
#'   probabilistic blend lives in `.no3_weights()`.
#' @return Character: "NO3-N_soft", "NO3-N_mod", or "NO3-N_hard".
#' @keywords internal
.no3_class <- function(hardness_mg_L, hardness_cv = 0.05) {
  # Hard cutoffs per ANZG 2025 class boundaries. Used where a single
  # representative class is needed (ssd_hc50() / the CA mixture step).
  # ssd_paf() instead uses .no3_weights() for a smooth probabilistic blend.
  if (is.null(hardness_mg_L) || is.na(hardness_mg_L)) {
    warning("NO3-N: hardness_mg_L not supplied \u2014 defaulting to soft-water ",
            "class (most conservative). Supply hardness_mg_L for correct ",
            "class selection.")
    return("NO3-N_soft")
  }
  if      (hardness_mg_L <  30)  "NO3-N_soft"
  else if (hardness_mg_L <= 150) "NO3-N_mod"
  else                           "NO3-N_hard"
}

#' Probabilistic NO3-N hardness-class weights.
#'
#' The three nitrate SSDs (soft / moderate / hard) have abrupt class boundaries
#' at 30 and 150 mg/L CaCO3. Rather than snap a sample to one class, treat the
#' true hardness as log-normal around the measured value (CV = `hardness_cv`)
#' and weight the classes by the probability mass falling in each band. This
#' smooths the PAF across boundaries without biologically interpolating between
#' the distinct SSDs.
#'
#' @param hardness_mg_L Numeric. Measured hardness in mg/L as CaCO3.
#' @param hardness_cv   Numeric. CV of the hardness measurement (default 0.05).
#' @return Named numeric vector `c(soft, mod, hard)` summing to 1, or all-`NA`
#'   if `hardness_mg_L` is missing/non-positive.
#' @keywords internal
.no3_weights <- function(hardness_mg_L, hardness_cv = 0.05) {
  if (is.null(hardness_mg_L) || is.na(hardness_mg_L) || hardness_mg_L <= 0)
    return(c(soft = NA_real_, mod = NA_real_, hard = NA_real_))
  s  <- sqrt(log(1 + hardness_cv^2))          # log-scale sd; 0 -> hard cutoff
  ml <- log(hardness_mg_L)
  p_soft <- stats::plnorm(30,  meanlog = ml, sdlog = s)
  p_hard <- stats::plnorm(150, meanlog = ml, sdlog = s, lower.tail = FALSE)
  c(soft = p_soft, mod = 1 - p_soft - p_hard, hard = p_hard)
}

# ── Public API ────────────────────────────────────────────────────────────────

#' Estimate the fraction of species potentially affected at a concentration.
#'
#' @param analyte      Character. Analyte name (key in .SSD_NAME_MAP).
#'   Supply "NO3-N" together with `hardness` for automatic class
#'   selection, or supply the explicit class name ("NO3-N_soft" etc.).
#' @param conc         Numeric or `units` object. Concentration to evaluate.
#'   Bare numeric requires `conc_units`.
#' @param conc_units   Character. Unit of `conc` when it is bare numeric,
#'   e.g. `"ug/L"` or `"mg/L"`. Ignored when `conc` is a `units` object.
#' @param method       Character. "multi" (default) fits all 6 BCANZ
#'   distributions and model-averages; "anzecc" uses the per-analyte
#'   distribution that best matches the original ANZG derivation.
#' @param hardness     Numeric or `units` object, or `NULL`. Required for
#'   the NO3-N analyte.  Hardness at the time of measurement.  Bare numeric
#'   requires `hardness_units`.
#' @param hardness_units Character. Unit of `hardness` when it is bare
#'   numeric, e.g. `"mg/L"`.  Ignored when `hardness` is a `units` object
#'   or `NULL`.
#' @param hardness_cv  Numeric. CV of the hardness measurement, used to weight
#'   the three NO3-N hardness-class SSDs probabilistically (a log-normal blend
#'   that smooths the soft/moderate/hard boundaries). Default 0.05 (5%). Set
#'   to 0 to recover hard class cutoffs.
#' @param guideline_dir Character. Path to the "guideline data" folder
#'   containing ANZG XLSX files. Falls back to
#'   getOption("leachatetools.guideline_dir").
#' @param nboot        Integer. Bootstrap replicates for CI. 0 = no CI.
#' @param level        Numeric. Confidence level (default 0.95).
#'
#' @return Named list:
#'   $analyte        character
#'   $conc_ug_L      numeric — resolved concentration in µg/L
#'   $method         character
#'   $pct            numeric — % species affected, NA if no model
#'   $lower          numeric — lower CI %
#'   $upper          numeric — upper CI %
#'   $note           character vector — caveats
#' @examples
#' # Fraction of species affected by 9.3 mg/L total ammonia-N at the pH 7.0 /
#' # 20 degC index condition. Uses the package's bundled SSD data; set
#' # options(leachatetools.guideline_dir=) to fit from the ANZG XLSX files.
#' ssd_paf("NH3-N", conc = 9321, conc_units = "ug/L")
#'
#' # NO3-N needs hardness for automatic soft/moderate/hard SSD selection:
#' ssd_paf("NO3-N", conc = 50000, conc_units = "ug/L",
#'         hardness = 90, hardness_units = "mg/L")
#' @export
ssd_paf <- function(analyte,
                    conc,
                    conc_units     = NULL,
                    method         = c("multi", "anzecc"),
                    hardness       = NULL,
                    hardness_units = NULL,
                    hardness_cv    = 0.05,
                    guideline_dir  = getOption("leachatetools.guideline_dir"),
                    nboot          = 0L,
                    level          = 0.95) {
  method <- match.arg(method)
  stopifnot(is.character(analyte), length(analyte) == 1L)
  conc_ug_L    <- .resolve_to(conc, "ug/L", conc_units, "conc")
  hardness_mg_L <- if (!is.null(hardness))
    .resolve_to(hardness, "mg/L", hardness_units, "hardness") else NULL
  stopifnot(is.numeric(conc_ug_L), length(conc_ug_L) == 1L, conc_ug_L > 0)

  result <- list(analyte = analyte, conc_ug_L = conc_ug_L, method = method,
                 pct = NA_real_, lower = NA_real_, upper = NA_real_,
                 note = character(0))

  # NO3-N: blend the soft/moderate/hard SSDs by hardness-class probability
  # (smooth) rather than snapping to one class (abrupt). Falls back to the
  # single-class path when hardness is missing.
  if (analyte == "NO3-N") {
    w_full <- .no3_weights(hardness_mg_L, hardness_cv)
    if (!anyNA(w_full)) {
      classes <- c(soft = "NO3-N_soft", mod = "NO3-N_mod", hard = "NO3-N_hard")
      # Ignore negligible-weight classes so a far class never contributes NA.
      keep <- w_full > 1e-6
      w    <- w_full[keep] / sum(w_full[keep])
      sub  <- lapply(classes[keep], function(cl)
        ssd_paf(cl, conc_ug_L, conc_units = "ug/L", method = method,
                hardness = NULL, guideline_dir = guideline_dir,
                nboot = nboot, level = level))
      result$pct <- sum(w * vapply(sub, function(x) x$pct, numeric(1)))
      if (nboot > 0L) {
        result$lower <- sum(w * vapply(sub, function(x) x$lower, numeric(1)))
        result$upper <- sum(w * vapply(sub, function(x) x$upper, numeric(1)))
      }
      result$note <- sprintf(
        "NO3-N hardness blend (soft/mod/hard = %.2f/%.2f/%.2f)",
        w_full[["soft"]], w_full[["mod"]], w_full[["hard"]])
      return(result)
    }
    # hardness missing -> single representative class (with the soft default warning)
    analyte <- .no3_class(hardness_mg_L, hardness_cv)
    result$analyte <- analyte
  }

  if (analyte %in% .SSD_NO_MODEL) {
    result$note <- "No SSD available (QSAR/LR derivation, insufficient data, or parked)"
    return(result)
  }

  stem <- if (analyte %in% names(.SSD_NAME_MAP)) .SSD_NAME_MAP[[analyte]] else NULL
  if (is.null(stem)) {
    result$note <- paste0("Unknown analyte '", analyte, "'.")
    return(result)
  }

  fit <- tryCatch(
    .load_or_fit(analyte, stem, method, guideline_dir),
    error = function(e) { result$note <<- conditionMessage(e); NULL }
  )
  if (is.null(fit)) return(result)

  notes <- character(0)
  if (isTRUE(attr(fit, "acute_data")))
    notes <- c(notes, "SSD fitted to acute data (MR method); PAF is approximate")
  if (isTRUE(attr(fit, "acr_applied")))
    notes <- c(notes, paste0("ACR pre-applied (\u00f7", attr(fit, "trigger_divisor"),
                             "); conc_ug_L should be on chronic-equivalent scale"))

  ci_flag <- nboot > 0L
  hp_args <- list(x = fit, conc = conc_ug_L, ci = ci_flag,
                  level = level, proportion = TRUE)
  if (ci_flag) hp_args$nboot <- as.integer(nboot)

  hp_tbl <- tryCatch(
    do.call(ssdtools::ssd_hp, hp_args),
    error = function(e) {
      notes <<- c(notes, paste("ssd_hp error:", conditionMessage(e))); NULL
    }
  )

  if (!is.null(hp_tbl)) {
    result$pct   <- hp_tbl$est * 100
    if (ci_flag) {
      result$lower <- hp_tbl$lcl * 100
      result$upper <- hp_tbl$ucl * 100
    }
  }
  result$note <- notes
  result
}

#' Return HC50 for use in Concentration Addition (msPAF).
#'
#' HC50 is the concentration at which 50% of species are predicted to be
#' affected. Used as the denominator when computing Toxic Units for the
#' Concentration Addition combination step of msPAF.
#'
#' @inheritParams ssd_paf
#' @return Numeric scalar — HC50 in µg/L, or NA.
#' @examples
#' # HC50 (50% of species affected) for copper, used as the
#' # Toxic-Unit denominator in the msPAF concentration-addition step:
#' ssd_hc50("Cu")
#' @export
ssd_hc50 <- function(analyte,
                     method         = c("multi", "anzecc"),
                     hardness       = NULL,
                     hardness_units = NULL,
                     guideline_dir  = getOption("leachatetools.guideline_dir")) {
  method <- match.arg(method)
  hardness_mg_L <- if (!is.null(hardness))
    .resolve_to(hardness, "mg/L", hardness_units, "hardness") else NULL
  if (analyte == "NO3-N") analyte <- .no3_class(hardness_mg_L)
  if (analyte %in% .SSD_NO_MODEL) return(NA_real_)

  stem <- if (analyte %in% names(.SSD_NAME_MAP)) .SSD_NAME_MAP[[analyte]] else NULL
  if (is.null(stem)) return(NA_real_)

  fit <- tryCatch(
    .load_or_fit(analyte, stem, method, guideline_dir),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)

  tryCatch(
    ssdtools::ssd_hc(fit, proportion = 0.5, ci = FALSE)$est,
    error = function(e) NA_real_
  )
}

#' Quick point-estimate only (no CI, no bootstrapping).
#'
#' @inheritParams ssd_paf
#' @return Numeric — % species affected, or NA.
#' @examples
#' ssd_pct("Zn", conc = 30, conc_units = "ug/L")
#' @export
ssd_pct <- function(analyte, conc, conc_units = NULL, method = "multi",
                    hardness = NULL, hardness_units = NULL,
                    guideline_dir = getOption("leachatetools.guideline_dir")) {
  ssd_paf(analyte, conc, conc_units = conc_units, method = method,
          hardness = hardness, hardness_units = hardness_units,
          guideline_dir = guideline_dir, nboot = 0L)$pct
}

# ── Internal helpers for msPAF ────────────────────────────────────────────────

# Effective log-normal sigma derived from HC5 and HC50 of the fitted SSD.
# Used in compute_ca_group_mspaf() to build the concentration-addition
# mixture SSD without needing the raw distribution parameters.
#
# sigma = (log10(HC50) - log10(HC05)) / (-qnorm(0.05))
# This is exact for log-normal SSDs and a good approximation for others.
.ssd_sigma <- function(analyte, method, guideline_dir) {
  if (analyte == "NO3-N") analyte <- .no3_class(NULL)
  if (analyte %in% .SSD_NO_MODEL) return(NA_real_)
  stem <- if (analyte %in% names(.SSD_NAME_MAP)) .SSD_NAME_MAP[[analyte]] else NULL
  if (is.null(stem)) return(NA_real_)

  fit <- tryCatch(
    .load_or_fit(analyte, stem, method, guideline_dir),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)

  hc50 <- tryCatch(
    ssdtools::ssd_hc(fit, proportion = 0.5, ci = FALSE)$est,
    error = function(e) NA_real_
  )
  hc05 <- tryCatch(
    ssdtools::ssd_hc(fit, proportion = 0.05, ci = FALSE)$est,
    error = function(e) NA_real_
  )

  if (anyNA(c(hc50, hc05)) || hc05 <= 0 || hc50 <= 0) return(NA_real_)
  (log10(hc50) - log10(hc05)) / (-stats::qnorm(0.05))
}
