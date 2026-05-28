## ============================================================================
## Adjusted multi-substance Potentially Affected Fraction (AmsPAF)
## ============================================================================
##
## Estimates the fraction of aquatic species potentially affected by the
## mixture of toxicants at a monitoring location, adjusted for local
## geogenic background using the Added Risk Approach (ARA).
##
## ## Conceptual basis
##
## Standard msPAF (De Zwart & Posthuma 2005, Environmental Toxicology and
## Chemistry 24(10): 2665-2676) computes the fraction of species affected
## by a chemical mixture using Species Sensitivity Distributions (SSDs).
## SSDs are cumulative distributions of single-species toxicity thresholds
## (NOEC, EC10), fitted as log-normal curves. The fraction of species with
## thresholds below a measured concentration is the Potentially Affected
## Fraction (PAF) for that substance.
##
## For mixtures, substances acting by the same mechanism are combined via
## Concentration Addition (CA); dissimilar mechanisms are combined via
## Independent Action (IA):
##   - CA: msPAF = F_mix(sum(TU_i)) where TU_i = C_i / HC50_i
##   - IA: msPAF = 1 - prod(1 - PAF_i)
##   - Hybrid: CA within mode-of-action groups, IA across groups
##
## ## Local adjustment (ARA)
##
## In natural waters with elevated background metals (e.g. Zn in low-pH
## Sydney Basin sandstone groundwater), the resident community is
## pollution-induced-community-tolerant (PICT; Blanck et al. 1988) to
## geogenic concentrations. Assessing toxicity against raw measured
## concentrations would overstate the anthropogenic risk.
##
## The Added Risk Approach (Crommentuijn et al. 1997, RIVM Report 601501001;
## Struijs et al. 1997, Ecotoxicology and Environmental Safety 37(2):112-118)
## addresses this by assessing toxicity only on the anthropogenic increment
## above local background:
##
##   C_adj_i = max(C_norm_i - ref_norm_i, 0)
##
## where ref_norm_i is the normalised reference concentration (default:
## geometric mean of matched reference-site data, normalised to ANZG index
## conditions via the same chemistry normalisation applied to sample
## concentrations).  See `prepare_reference()` for alternatives.  Evaluating
## the SSD at C_adj_i is equivalent to evaluating a rightward-shifted SSD at
## C_norm_i, consistent with PICT theory.
##
## ## Chemistry normalisation
##
## Some ANZG SSDs are derived at specific water chemistry conditions (e.g.
## Cu at DOC = 0.5 mg/L, Zn at hardness 30 mg/L). Normalisation adjusts
## measured field concentrations to these index conditions before SSD lookup.
## Formulas are stored as R-expression strings in the analyte metadata CSV
## (`normalisation_formula` column) and parsed once at startup via
## `.parse_normalisation_formula()`. All formulas are stubs (identity) until
## populated from the ANZG technical briefs.
##
## ## SSD derivation
##
## SSD parameters are obtained from the fitted models in paf.R via
## ssd_hc50() (for the CA HC50 denominator) and ssd_paf() (for individual
## PAF evaluation). This uses the same BCANZ-validated SSD infrastructure
## as the single-substance PAF calculations, and is correct for any SSD
## distribution shape (not only log-normal).
##
## ## Mode of action groupings
##
## MOA group assignments are read from the `moa_group` column in the bundled
## analyte metadata CSV. Analytes with the same `moa_group` are combined by
## Concentration Addition (CA); groups are then combined via Independent Action
## (IA).
##
## Day-1 state: metals and inorganic nitrogen/sulphide species have
## `moa_group = "ionoregulatory"` (one CA group). Organics/pesticides have
## `moa_group = NA` → each gets a unique synthetic solo group `"_solo_<name>"`.
##
## ## Output interpretation
##
## AmsPAF is returned as a continuous risk metric (% of species potentially
## affected by the mixture, 0-100+).  Tier breaks are deliberately not
## provided by this package — msPAF semantics differ from single-substance
## guideline values (the "% affected" denominator depends on which
## substances are in the mixture).  See
## `vignette("chronic-amspaf-interpretation")` for a discussion of how to
## interpret msPAF values in an assessment context.
##
## ## References
##
## Blanck H, Wallin G, Waengberg SA (1988) Species-dependent variation in
##   algal sensitivity to chemical compounds. Ecotoxicology and Environmental
##   Safety 8(4):339-351.
## Crommentuijn T, Polder MD, van de Plassche EJ (1997) Maximum permissible
##   concentrations and negligible concentrations for metals, taking
##   background concentrations into account. RIVM Report 601501001.
## De Zwart D, Posthuma L (2005) Complex mixture toxicity for single and
##   multiple species: proposed methodologies. Environmental Toxicology and
##   Chemistry 24(10):2665-2676.
## Struijs J et al. (1997) Added risk approach to derive maximum permissible
##   concentrations for heavy metals. Ecotoxicology and Environmental Safety
##   37(2):112-118.
## Warne MStJ et al. (2018) Revised Method for Deriving Australian and New
##   Zealand Water Quality Guideline Values for Toxicants. ANZG.
##
## ============================================================================

## Analytes excluded from msPAF regardless of SSD availability.
##
##   NO3-N: GVs from NZ document of uncertain provenance, not ANZG
##   CH4:   Methane guideline is aesthetic/nuisance, not an ecotoxicity SSD
##   LHF:   Leachate indicator index — derived values, not a toxicity endpoint
##
.AMSPAF_EXCLUDED_ANALYTES <- c("NO3-N", "CH4", "LHF")


## ============================================================================
## add_amspaf
## ============================================================================

#' Compute the Adjusted multi-substance PAF (AmsPAF) for water quality samples
#'
#' Appends AmsPAF rows to a long-format water quality dataframe. AmsPAF
#' estimates the fraction of aquatic species potentially affected by the
#' combined toxicant mixture, adjusted for local geogenic background via the
#' Added Risk Approach. See the file-level header for full methodological
#' detail.
#'
#' The function accepts either per-sample or chronic-integrated chemistry (from
#' [time_weighted_aggregate()]). It does not need to know which — the
#' distinction is entirely in the input data. Similarly, `reference` may be a
#' raw long-format chemistry data frame or a pre-built [prepare_reference()]
#' object.
#'
#' @param df Long-format monitoring dataframe. Required columns:
#'   `sample_id`, `site_id`, `analyte`, `value` (concentrations in
#'   µg/L). Optional but recommended: `datetime` (propagated to AmsPAF
#'   rows if present), `detected` (assumed `TRUE` if absent), `imputed`
#'   (logical; if present, `n_analytes_imputed` is populated in output).
#'   Driver analytes needed for chemistry normalisation (e.g. pH, DOC) should
#'   be present as rows in `df`.
#' @param reference Background reference chemistry for the ARA adjustment.
#'   Accepts three forms:
#'   \itemize{
#'     \item A `prepared_reference` object from [prepare_reference()] —
#'       normalisation has already been applied; used directly.
#'     \item A raw long-format data frame (same schema as `df`) — will be
#'       passed to [prepare_reference()] internally.
#'     \item `NULL` (default) — no ARA adjustment; raw concentrations assessed
#'       directly against SSDs.
#'   }
#' @param analyte_metadata Data frame of analyte metadata, or `NULL` to load
#'   the bundled `inst/extdata/anzecc_analyte_metadata.csv`. Passed to
#'   [prepare_reference()] and [derive_ssd_params()].
#' @param method SSD method passed to [ssd_hc50()] and [ssd_paf()].
#'   `"multi"` (default) fits all 6 BCANZ distributions and model-averages;
#'   `"anzecc"` uses the per-analyte distribution matching the original ANZG
#'   derivation.
#' @param guideline_dir Path to the "guideline data" folder containing ANZG
#'   XLSX files. Falls back to `getOption("leachatetools.guideline_dir")`.
#' @param min_analytes Minimum number of analytes with fitted SSDs required
#'   to compute AmsPAF for a sample. Default `3`.
#' @param ref_summary Summary statistic for the reference distribution when
#'   `reference` is a raw data frame.  Passed through to
#'   [prepare_reference()].  Default `"geom_mean"` — the maximum-likelihood
#'   central tendency for log-normal concentrations, and a PICT-consistent
#'   estimate of the "typical" exposure the resident community has adapted
#'   to.  Other options: `"median"`, `"arith_mean"`, `"p80"`, `"p90"`,
#'   `"p95"`.
#'
#' @return The input `df` with AmsPAF rows appended. Each AmsPAF row carries:
#'   `value` (AmsPAF as a percentage, 0–100+), `detected = TRUE`,
#'   `analyte = "AmsPAF"`, `n_analytes_used` (integer),
#'   `n_analytes_imputed` (integer, 0 if `imputed` column absent),
#'   `dominant_analyte` (character), `max_paf` (numeric),
#'   `analyte_pafs` (list column of per-analyte diagnostic tibbles).
#'
#'   Tier breaks are not provided by this package — AmsPAF is a continuous
#'   risk metric and the threshold at which a community is "impacted"
#'   depends on the assessment context.  See
#'   `vignette("chronic-amspaf-interpretation")` for guidance.
#'
#' @seealso [ssd_paf()], [ssd_hc50()], [prepare_reference()],
#'   [time_weighted_aggregate()], [prescreen_analytes()], [impute_chemistry()]
#'
#' @references
#' De Zwart D, Posthuma L (2005) Environmental Toxicology and Chemistry
#' 24(10):2665-2676.
#'
#' @export
add_amspaf <- function(
    df,
    reference        = NULL,
    analyte_metadata = NULL,
    method           = c("multi", "anzecc"),
    guideline_dir    = getOption("leachatetools.guideline_dir"),
    min_analytes     = 3,
    ref_summary      = c("geom_mean", "median", "arith_mean",
                          "p80", "p90", "p95")
) {
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c("sample_id", "site_id", "analyte", "value")
  )
  method      <- match.arg(method)
  ref_summary <- match.arg(ref_summary)
  checkmate::assert_int(min_analytes, lower = 1L)

  ## ================================================================
  ## Step 1: Load analyte metadata and derive SSD parameters.
  ##
  ## Done once per function call — not per sample.
  ## ================================================================

  meta       <- .load_analyte_metadata(analyte_metadata)
  ssd_params <- derive_ssd_params(meta, method, guideline_dir)

  if (nrow(ssd_params) == 0L) {
    cli::cli_warn(
      "No analytes with fitted SSDs found. AmsPAF cannot be computed."
    )
    return(df)
  }

  ## ================================================================
  ## Step 2: Resolve the reference into a prepared_reference object.
  ## ================================================================

  if (inherits(reference, "prepared_reference")) {
    prep_ref <- reference
  } else if (is.null(reference)) {
    ## No ARA: empty ref_table → ref_norm joined as NA → coerced to 0
    prep_ref <- structure(
      list(
        ref_table = tibble::tibble(
          analyte  = character(0),
          ref_norm = numeric(0),
          n_obs    = integer(0)
        ),
        dropped = character(0),
        summary = ref_summary
      ),
      class = "prepared_reference"
    )
  } else {
    checkmate::assert_data_frame(reference)
    prep_ref <- prepare_reference(
      reference,
      analyte_metadata = meta,
      summary          = ref_summary
    )
  }

  ## Use only analyte+ref_norm columns for the join (avoid polluting tox_rows
  ## with n_obs / ref_lower / ref_upper if bootstrap_ci was used)
  ref_table <- if (nrow(prep_ref$ref_table) > 0L) {
    dplyr::select(prep_ref$ref_table, "analyte", "ref_norm")
  } else {
    prep_ref$ref_table
  }

  ## ================================================================
  ## Step 3: Compute AmsPAF per feature, per sample.
  ## ================================================================

  amspaf_df <-
    df |>
    dplyr::group_by(.data$site_id) |>
    dplyr::group_modify(\(.x, .y) {
      compute_amspaf_per_sample(
        sample_data   = .x,
        ref_table     = ref_table,
        ssd_params    = ssd_params,
        min_analytes  = min_analytes,
        method        = method,
        guideline_dir = guideline_dir
      )
    }) |>
    dplyr::ungroup()

  ## End-of-call summary of dropped analytes (single message rather than
  ## per-sample warnings — important for large datasets)
  .summarise_amspaf_diagnostics(amspaf_df, min_analytes)

  ## dropped_analytes is a diagnostic list-column; remove from the final
  ## output rows (still emitted in summary above)
  if ("dropped_analytes" %in% names(amspaf_df))
    amspaf_df <- dplyr::select(amspaf_df, -"dropped_analytes")

  amspaf_df <- dplyr::mutate(amspaf_df,
    analyte  = "AmsPAF",
    detected = TRUE
  )

  ## Tier breaks / regulatory interpretation are intentionally NOT provided
  ## by this package.  AmsPAF is a continuous risk metric (% of species
  ## potentially affected); the threshold at which a community is
  ## "impacted" depends on the assessment context (mixture composition,
  ## site-specific calibration, target protection level) and is a
  ## consumer-side decision.  See vignette("chronic-amspaf-interpretation").

  ## Propagate datetime from input df to AmsPAF rows (per-sample pipeline).
  if ("datetime" %in% names(df)) {
    sample_times <- dplyr::distinct(df, .data$sample_id, .data$datetime)
    amspaf_df    <- dplyr::left_join(amspaf_df, sample_times, by = "sample_id")
  }

  ## Propagate focal_date from input df to AmsPAF rows (chronic pipeline).
  if ("focal_date" %in% names(df)) {
    focal_times <- dplyr::distinct(df, .data$sample_id, .data$focal_date)
    amspaf_df   <- dplyr::left_join(amspaf_df, focal_times, by = "sample_id")
  }

  result <- dplyr::bind_rows(df, amspaf_df)
  if ("focal_date" %in% names(result)) {
    result <- dplyr::arrange(result, .data$focal_date)
  } else if ("datetime" %in% names(result)) {
    result <- dplyr::arrange(result, .data$datetime)
  }
  result
}


## ============================================================================
## derive_ssd_params
## ============================================================================

#' Derive SSD parameters for AmsPAF computation
#'
#' Reads analyte eligibility and mode-of-action groups from the bundled
#' metadata CSV, then calls [ssd_hc50()] and [.ssd_sigma()] to populate the
#' HC50 and effective sigma needed for Concentration Addition. Chemistry
#' normalisation formulas (currently all stubs) are parsed once and stored
#' as a list column for use in [compute_amspaf_per_sample()].
#'
#' Eligibility criteria:
#' \itemize{
#'   \item `ssd_available == TRUE` in the metadata
#'   \item `analyte` not in `.AMSPAF_EXCLUDED_ANALYTES`
#'   \item `ssd_hc50()` returns a non-NA value (model fits successfully)
#' }
#'
#' @param meta Analyte metadata tibble from [.load_analyte_metadata()].
#' @param method SSD method: `"multi"` or `"anzecc"`.
#' @param guideline_dir Path to ANZG guideline data folder.
#'
#' @return Tibble with columns `analyte`, `hc50`, `sigma`, `moa_group`,
#'   `parsed_formula` (list of language objects or NULLs),
#'   `coanalytes_req` (character).
#'
#' @keywords internal
derive_ssd_params <- function(meta, method, guideline_dir) {
  eligible <- meta |>
    dplyr::filter(.data$ssd_available == TRUE) |>
    dplyr::filter(!.data$analyte %in% .AMSPAF_EXCLUDED_ANALYTES)

  if (nrow(eligible) == 0L) {
    return(.empty_ssd_params())
  }

  params <- purrr::map_dfr(seq_len(nrow(eligible)), function(i) {
    nm  <- eligible$analyte[i]
    hc50 <- ssd_hc50(nm, method = method, guideline_dir = guideline_dir)
    if (is.na(hc50)) return(NULL)

    sigma <- .ssd_sigma(nm, method = method, guideline_dir = guideline_dir)

    ## MOA group from metadata column; NA/empty → unique solo group
    mg_raw <- eligible$moa_group[i]
    moa_group <- if (is.na(mg_raw) || !nzchar(mg_raw)) {
      paste0("_solo_", nm)
    } else {
      mg_raw
    }

    ## Parse normalisation formula once (cached in .normalise_cache)
    formula_str  <- eligible$normalisation_formula[i] %||% ""
    parsed_f     <- .parse_normalisation_formula(formula_str)
    coanalytes_r <- eligible$coanalytes_required[i] %||% ""

    tibble::tibble(
      analyte        = nm,
      hc50           = hc50,
      sigma          = sigma,
      moa_group      = moa_group,
      parsed_formula = list(parsed_f),
      coanalytes_req = coanalytes_r
    )
  })

  if (is.null(params) || nrow(params) == 0L) .empty_ssd_params() else params
}

.empty_ssd_params <- function() {
  tibble::tibble(
    analyte        = character(),
    hc50           = numeric(),
    sigma          = numeric(),
    moa_group      = character(),
    parsed_formula = list(),
    coanalytes_req = character()
  )
}


## ============================================================================
## compute_ca_group_mspaf
## ============================================================================

#' Compute msPAF for a single Concentration Addition group
#'
#' @param group_data Tibble with columns `C_adj`, `hc50`, `sigma`, `moa_group`;
#'   one row per analyte in the CA group.
#' @return msPAF as a proportion (not percentage).
#' @keywords internal
compute_ca_group_mspaf <- function(group_data) {
  group_data <- dplyr::filter(
    group_data,
    .data$C_adj > 0,
    !is.na(.data$hc50), .data$hc50 > 0,
    !is.na(.data$sigma), .data$sigma > 0
  )

  if (nrow(group_data) == 0L) return(0)

  group_data <- dplyr::mutate(group_data, TU = .data$C_adj / .data$hc50)
  TU_mix <- sum(group_data$TU)
  if (TU_mix <= 0) return(0)

  group_data <- dplyr::mutate(group_data, w = .data$TU / TU_mix)
  sigma_mix  <- sqrt(sum(group_data$w^2 * group_data$sigma^2))

  pnorm(log10(TU_mix) / sigma_mix)
}


## ============================================================================
## compute_amspaf_per_sample
## ============================================================================

#' Compute AmsPAF for each sample in a per-feature data block
#'
#' Internal workhorse called by [add_amspaf()] for each `uuid.feature` group.
#' Applies chemistry normalisation, ARA adjustment, and CA/IA mixture
#' combination per sample.
#'
#' @param sample_data Per-feature long-format df (may include co-analyte rows
#'   such as pH, DOC alongside toxicant rows).
#' @param ref_table Tibble `(analyte, ref_norm)` from `prep_ref$ref_table`.
#' @param ssd_params Tibble from [derive_ssd_params()].
#' @param min_analytes Minimum analytes required.
#' @param method SSD method.
#' @param guideline_dir Path to ANZG XLSX folder.
#'
#' @return Tibble with one row per sample that passes `min_analytes`, columns:
#'   `sample_id`, `value`, `n_analytes_used`, `n_analytes_imputed`,
#'   `dominant_analyte`, `max_paf`, `analyte_pafs`, `dropped_analytes`.
#' @keywords internal
compute_amspaf_per_sample <- function(
    sample_data,
    ref_table,
    ssd_params,
    min_analytes,
    method,
    guideline_dir
) {
  has_detected   <- "detected" %in% names(sample_data)
  has_imputed    <- "imputed"  %in% names(sample_data)

  if (!has_detected) {
    sample_data <- dplyr::mutate(sample_data, detected = TRUE)
  }

  sample_data |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::group_modify(\(.x, .y) {
      ## Build co-analyte lookup (all detected values in this sample)
      coanalyte_vals <- .x |>
        dplyr::filter(.data$detected) |>
        dplyr::select("analyte", "value") |>
        tibble::deframe()  # named numeric vector

      ## Filter to SSD-eligible analytes
      tox_rows <- .x |>
        dplyr::filter(.data$analyte %in% ssd_params$analyte) |>
        dplyr::left_join(
          dplyr::select(ssd_params, "analyte", "hc50", "sigma",
                        "moa_group", "parsed_formula", "coanalytes_req"),
          by = "analyte"
        ) |>
        dplyr::left_join(ref_table, by = "analyte") |>
        dplyr::mutate(
          ref_norm = tidyr::replace_na(.data$ref_norm, 0)
        )

      empty_row <- tibble::tibble(
        value              = numeric(0),
        n_analytes_used    = integer(0),
        n_analytes_imputed = integer(0),
        dominant_analyte   = character(0),
        max_paf            = numeric(0),
        analyte_pafs       = list(),
        dropped_analytes   = list()
      )

      if (nrow(tox_rows) < min_analytes) return(empty_row)

      ## ── Chemistry normalisation ────────────────────────────────────────
      ##
      ## BDL (quantified == FALSE): treated as zero exposure; C_norm = 0.
      ## Detected: apply normalisation formula using co-analyte values.
      ## Rows where normalisation returns NA (missing required co-analyte)
      ## are dropped and recorded (per-sample list column; summarised
      ## at the end of add_amspaf()).
      tox_rows <- dplyr::mutate(
        tox_rows,
        C_norm = purrr::pmap_dbl(
          list(
            q   = .data$detected,
            C   = .data$value,
            pf  = .data$parsed_formula,
            cr  = .data$coanalytes_req
          ),
          function(q, C, pf, cr) {
            if (!q) return(0)             # BDL → zero exposure
            co_names <- if (nzchar(cr %||% "")) {
              trimws(strsplit(cr, ",")[[1L]])
            } else character(0)
            co_names <- co_names[nzchar(co_names)]
            co_vals  <- coanalyte_vals[co_names[co_names %in% names(coanalyte_vals)]]
            .apply_normalisation(pf, C, co_vals)
          }
        )
      )

      ## Capture dropped analyte names + reason before filtering
      dropped <- dplyr::filter(tox_rows, is.na(.data$C_norm)) |>
        dplyr::transmute(
          .data$analyte,
          reason = "missing_co_analyte"
        )
      tox_rows <- dplyr::filter(tox_rows, !is.na(.data$C_norm))

      if (nrow(tox_rows) < min_analytes) return(empty_row)

      ## ARA shift
      tox_rows <- dplyr::mutate(
        tox_rows,
        C_adj = pmax(.data$C_norm - .data$ref_norm, 0)
      )

      ## Count imputed analytes (rows from impute_chemistry() that were BDL/missing)
      n_analytes_imputed <- if (has_imputed && "imputed" %in% names(tox_rows)) {
        sum(tox_rows$imputed, na.rm = TRUE)
      } else 0L

      ## ── Per-analyte PAF (diagnostic + dominant-analyte identification) ─
      tox_rows <- dplyr::mutate(tox_rows,
        PAF = purrr::map2_dbl(
          .data$C_adj, .data$analyte,
          function(c, a) {
            if (is.na(c) || c <= 0) return(0)
            paf_result <- tryCatch(
              ssd_paf(a, c, method = method, guideline_dir = guideline_dir, nboot = 0L),
              error = function(e) list(pct = NA_real_)
            )
            if (is.na(paf_result$pct)) NA_real_ else paf_result$pct / 100
          }
        )
      )

      ## ── CA msPAF per mode-of-action group ─────────────────────────────
      groups         <- unique(tox_rows$moa_group)
      mspaf_by_group <- vapply(
        groups,
        function(g) {
          compute_ca_group_mspaf(dplyr::filter(tox_rows, .data$moa_group == g))
        },
        numeric(1)
      )

      ## ── Combine groups via Independent Action ─────────────────────────
      amspaf <- 1 - prod(1 - mspaf_by_group)

      dominant <- if (any(!is.na(tox_rows$PAF))) {
        tox_rows$analyte[which.max(tox_rows$PAF)]
      } else NA_character_

      tibble::tibble(
        value              = amspaf * 100,
        n_analytes_used    = nrow(tox_rows),
        n_analytes_imputed = as.integer(n_analytes_imputed),
        dominant_analyte   = dominant,
        max_paf            = if (nrow(tox_rows) > 0L) max(tox_rows$PAF, na.rm = TRUE) else NA_real_,
        analyte_pafs       = list(
          dplyr::select(tox_rows, "analyte", "C_adj", "PAF", "moa_group")
        ),
        dropped_analytes   = list(dropped)
      )
    }) |>
    dplyr::ungroup()
}

#' Summarise per-sample dropped-analyte tally and emit a single cli message
#' @keywords internal
.summarise_amspaf_diagnostics <- function(amspaf_df, min_analytes) {
  if (!"dropped_analytes" %in% names(amspaf_df) || nrow(amspaf_df) == 0L)
    return(invisible(NULL))

  drop_long <- amspaf_df |>
    dplyr::select(dplyr::any_of(c("sample_id", "site_id")), "dropped_analytes") |>
    tidyr::unnest(cols = "dropped_analytes")

  if (nrow(drop_long) == 0L) return(invisible(NULL))

  per_analyte <- drop_long |>
    dplyr::count(.data$analyte, .data$reason, name = "n_samples")

  n_total_drops <- nrow(drop_long)
  n_samples_aff <- dplyr::n_distinct(drop_long$sample_id)

  cli::cli_inform(c(
    "i" = "add_amspaf: {n_total_drops} analyte row{?s} dropped across \\
           {n_samples_aff} sample{?s} (normalisation co-analyte missing).",
    paste0("    ", per_analyte$analyte, ": ", per_analyte$n_samples,
           " sample", ifelse(per_analyte$n_samples == 1L, "", "s"))
  ))
  invisible(per_analyte)
}


