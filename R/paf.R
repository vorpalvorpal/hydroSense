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
# External corrections (apply before calling):
#   NH3-N : apply pH / temperature correction for un-ionised fraction
#   Cr    : no external correction needed (freshwater subset already applied)
#   NO3-N : supply hardness_mg_L to select the appropriate hardness class;
#            probabilistic class weighting (CV = 5%) is a TODO (see below)
#
# TODO: probabilistic NO3-N hardness weighting
#   True hardness is assumed log-normally distributed around the measured
#   value with CV = hardness_cv (default 0.05, i.e. 5%).
#   Weights: p_soft = plnorm(30,  log(h), log(1+cv^2)^0.5)
#            p_hard = 1 - plnorm(150, log(h), log(1+cv^2)^0.5)
#            p_mod  = 1 - p_soft - p_hard
#   Expected PAF = p_soft*PAF_soft + p_mod*PAF_mod + p_hard*PAF_hard
#   This smooths the abrupt class transitions without requiring any
#   biological interpolation between the three distinct SSDs.
#   Pending implementation — currently uses a hard class cutoff.

# ── Analyte registry ──────────────────────────────────────────────────────────

# Canonical dashboard name → safe file stem used in the .qs cache.
# NO3-N is intentionally absent: callers supply the hardness-specific variant
# (NO3-N_soft / NO3-N_mod / NO3-N_hard) or use the ssd_paf() hardness_mg_L
# argument which resolves the class automatically.
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
  "As"                     = "As",
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
  "Dimethoate"             = "Dimethoate",
  "Naphthalene"            = "Naphthalene",
  "Benzo(a)pyrene"         = "Benzo_a_pyrene",
  "Benzene"                = "Benzene",
  "o-Xylene"               = "o_Xylene",
  "Phenol"                 = "Phenol",
  "2-Chlorophenol"         = "2_Chlorophenol",
  "2,4,6-Trichlorophenol"  = "2_4_6_Trichlorophenol",
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
    fit <- qs::qread(cache_file)
    assign(key, fit, envir = .paf_mem_cache)
    return(fit)
  }

  # 3. Fit from raw data
  fit <- .fit_model(analyte, stem, method, guideline_dir)
  if (!is.null(fit)) {
    qs::qsave(fit, cache_file)
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
    warning("No metadata for analyte '", analyte, "' — cannot fit model.")
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
#' @param hardness_cv   Numeric. Coefficient of variation of hardness
#'   measurement (default 0.05 = 5%). Reserved for future probabilistic
#'   weighting — currently unused (hard cutoffs applied).
#' @return Character: "NO3-N_soft", "NO3-N_mod", or "NO3-N_hard".
.no3_class <- function(hardness_mg_L, hardness_cv = 0.05) {
  # TODO: implement probabilistic weighting using log-normal hardness
  # uncertainty (CV = hardness_cv). See module header for the formula.
  # For now: hard cutoffs per ANZG 2025 class boundaries.
  if (is.null(hardness_mg_L) || is.na(hardness_mg_L)) {
    warning("NO3-N: hardness_mg_L not supplied — defaulting to soft-water ",
            "class (most conservative). Supply hardness_mg_L for correct ",
            "class selection.")
    return("NO3-N_soft")
  }
  if      (hardness_mg_L <  30)  "NO3-N_soft"
  else if (hardness_mg_L <= 150) "NO3-N_mod"
  else                           "NO3-N_hard"
}

# ── Public API ────────────────────────────────────────────────────────────────

#' Estimate the fraction of species potentially affected at a concentration.
#'
#' @param analyte      Character. Analyte name (key in .SSD_NAME_MAP).
#'   Supply "NO3-N" together with `hardness_mg_L` for automatic class
#'   selection, or supply the explicit class name ("NO3-N_soft" etc.).
#' @param conc_ug_L    Numeric. Concentration in µg/L (after any external
#'   physicochemical corrections).
#' @param method       Character. "multi" (default) fits all 6 BCANZ
#'   distributions and model-averages; "anzecc" uses the per-analyte
#'   distribution that best matches the original ANZG derivation.
#' @param hardness_mg_L Numeric or NULL. Required for NO3-N analyte.
#'   Hardness in mg/L CaCO3 at the time of measurement.
#' @param hardness_cv  Numeric. CV of hardness measurement for probabilistic
#'   class weighting. Currently reserved — hard cutoffs used. Default 0.05.
#' @param guideline_dir Character. Path to the "guideline data" folder
#'   containing ANZG XLSX files. Falls back to
#'   getOption("leachatetools.guideline_dir").
#' @param nboot        Integer. Bootstrap replicates for CI. 0 = no CI.
#' @param level        Numeric. Confidence level (default 0.95).
#'
#' @return Named list:
#'   $analyte        character
#'   $conc_ug_L      numeric
#'   $method         character
#'   $pct            numeric — % species affected, NA if no model
#'   $lower          numeric — lower CI %
#'   $upper          numeric — upper CI %
#'   $note           character vector — caveats
#' @export
ssd_paf <- function(analyte,
                    conc_ug_L,
                    method         = c("multi", "anzecc"),
                    hardness_mg_L  = NULL,
                    hardness_cv    = 0.05,
                    guideline_dir  = getOption("leachatetools.guideline_dir"),
                    nboot          = 0L,
                    level          = 0.95) {
  method <- match.arg(method)
  stopifnot(is.character(analyte), length(analyte) == 1L)
  stopifnot(is.numeric(conc_ug_L), length(conc_ug_L) == 1L, conc_ug_L > 0)

  result <- list(analyte = analyte, conc_ug_L = conc_ug_L, method = method,
                 pct = NA_real_, lower = NA_real_, upper = NA_real_,
                 note = character(0))

  # Resolve NO3-N to a hardness-specific class
  if (analyte == "NO3-N") {
    analyte <- .no3_class(hardness_mg_L, hardness_cv)
    result$analyte <- analyte
  }

  if (analyte %in% .SSD_NO_MODEL) {
    result$note <- "No SSD available (QSAR/LR derivation, insufficient data, or parked)"
    return(result)
  }

  stem <- .SSD_NAME_MAP[[analyte]]
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
    notes <- c(notes, paste0("ACR pre-applied (÷", attr(fit, "trigger_divisor"),
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
#' @export
ssd_hc50 <- function(analyte,
                     method        = c("multi", "anzecc"),
                     hardness_mg_L = NULL,
                     guideline_dir = getOption("leachatetools.guideline_dir")) {
  method <- match.arg(method)
  if (analyte == "NO3-N") analyte <- .no3_class(hardness_mg_L)
  if (analyte %in% .SSD_NO_MODEL) return(NA_real_)

  stem <- .SSD_NAME_MAP[[analyte]]
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
#' @export
ssd_pct <- function(analyte, conc_ug_L, method = "multi",
                    hardness_mg_L = NULL,
                    guideline_dir = getOption("leachatetools.guideline_dir")) {
  ssd_paf(analyte, conc_ug_L, method = method,
          hardness_mg_L = hardness_mg_L,
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
  stem <- .SSD_NAME_MAP[[analyte]]
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
  (log10(hc50) - log10(hc05)) / (-qnorm(0.05))
}
