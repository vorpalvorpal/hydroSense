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
## where ref_norm_i is the normalised reference concentration (80th percentile
## of matched reference site data, normalised to ANZG index conditions via the
## same chemistry normalisation applied to sample concentrations). Evaluating
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
## ## Tier breaks
##
## All four tier breaks are biologically anchored to ANZG species protection
## levels, expressed as msPAF percentages:
##   Tier 1 (Background):  AmsPAF < 1%   (99% species protection)
##   Tier 2 (Elevated):    1% - 5%       (95% species protection)
##   Tier 3 (Impacted):    5% - 10%      (90% species protection)
##   Tier 4 (Severe):      > 10% - 20%+  (80% species protection)
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

## Tier break constants (expressed as proportions internally, % in output)
.AMSPAF_BREAK_T1_T2 <- 0.01 ## 1%  — 99% species protection
.AMSPAF_BREAK_T2_T3 <- 0.05 ## 5%  — 95% species protection
.AMSPAF_BREAK_T3_T4 <- 0.10 ## 10% — 90% species protection
.AMSPAF_CEILING <- 0.20 ## 20% — 80% species protection (notional ceiling)

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
#' [compute_chronic_chemistry()]). It does not need to know which — the
#' distinction is entirely in the input data. Similarly, `reference` may be a
#' raw long-format chemistry data frame or a pre-built [prepare_reference()]
#' object.
#'
#' @param df Long-format monitoring dataframe. Required columns:
#'   `uuid.sample`, `uuid.feature`, `name.analyte`, `value` (concentrations in
#'   µg/L). Optional but recommended: `datetime.sample` (propagated to AmsPAF
#'   rows if present), `quantified` (assumed `TRUE` if absent), `imputed`
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
#' @param ref_percentile_for_anchor Percentile used when `reference` is a raw
#'   data frame. Default `0.80` (80th percentile, per ANZG convention).
#'
#' @return The input `df` with AmsPAF rows appended. Each AmsPAF row carries:
#'   `value` (AmsPAF as a percentage, 0–100+), `quantified = TRUE`,
#'   `name.analyte = "AmsPAF"`, `n_analytes_used` (integer),
#'   `n_analytes_imputed` (integer, 0 if `imputed` column absent),
#'   `dominant_analyte` (character), `max_paf` (numeric),
#'   `analyte_pafs` (list column of per-analyte diagnostic tibbles), and four
#'   guideline columns (`value/level_name/guideline/comments.guideline_1`
#'   through `_4`).
#'
#' @seealso [ssd_paf()], [ssd_hc50()], [prepare_reference()],
#'   [compute_chronic_chemistry()], [prescreen_analytes()], [impute_chemistry()]
#'
#' @references
#' De Zwart D, Posthuma L (2005) Environmental Toxicology and Chemistry
#' 24(10):2665-2676.
#'
#' @export
add_amspaf <- function(
    df,
    reference                 = NULL,
    analyte_metadata          = NULL,
    method                    = c("multi", "anzecc"),
    guideline_dir             = getOption("leachatetools.guideline_dir"),
    min_analytes              = 3,
    ref_percentile_for_anchor = 0.80
) {
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c("uuid.sample", "uuid.feature", "name.analyte", "value")
  )
  method <- match.arg(method)
  checkmate::assert_int(min_analytes, lower = 1L)
  checkmate::assert_number(ref_percentile_for_anchor, lower = 0.5, upper = 0.99)

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
    ## No ARA: empty quantiles → ref_norm treated as 0 inside compute_amspaf_per_sample()
    prep_ref <- structure(
      list(
        normalised_quantiles = tibble::tibble(
          name.analyte = character(0),
          ref_norm     = numeric(0)
        ),
        dropped    = character(0),
        percentile = ref_percentile_for_anchor
      ),
      class = "prepared_reference"
    )
  } else {
    checkmate::assert_data_frame(reference)
    prep_ref <- prepare_reference(
      reference,
      analyte_metadata = meta,
      percentile       = ref_percentile_for_anchor
    )
  }

  ref_quantiles <- prep_ref$normalised_quantiles

  ## ================================================================
  ## Step 3: Compute AmsPAF per feature, per sample.
  ## ================================================================

  amspaf_df <-
    df |>
    dplyr::group_by(uuid.feature) |>
    dplyr::group_modify(\(.x, .y) {
      compute_amspaf_per_sample(
        sample_data   = .x,
        ref_quantiles = ref_quantiles,
        ssd_params    = ssd_params,
        min_analytes  = min_analytes,
        method        = method,
        guideline_dir = guideline_dir
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      name.analyte = "AmsPAF",
      uuid.analyte = NA_character_,
      uuid = uuid::UUIDgenerate(use.time = FALSE, n = dplyr::n()),
      quantified = TRUE,
      RL_low = 0,

      ## ================================================================
      ## Guideline columns — tier breaks anchored to ANZG species
      ## protection levels.
      ## ================================================================

      value.guideline_1 = .AMSPAF_BREAK_T1_T2 * 100,
      level_name.guideline_1 = "Background",
      guideline.guideline_1 = "AmsPAF threshold",
      comments.guideline_1 = paste0(
        "1% of species potentially affected by the mixture. Corresponds to ",
        "the ANZG 99% species protection level. Concentrations are adjusted ",
        "for local geogenic background via the Added Risk Approach before ",
        "SSD evaluation, so background metals do not contribute to AmsPAF."
      ),

      value.guideline_2 = .AMSPAF_BREAK_T2_T3 * 100,
      level_name.guideline_2 = "Elevated",
      guideline.guideline_2 = "AmsPAF threshold",
      comments.guideline_2 = paste0(
        "5% of species potentially affected by the mixture. Corresponds to ",
        "the ANZG 95% species protection level (the standard default for ",
        "slightly-to-moderately disturbed ecosystems)."
      ),

      value.guideline_3 = .AMSPAF_BREAK_T3_T4 * 100,
      level_name.guideline_3 = "Impacted",
      guideline.guideline_3 = "AmsPAF threshold",
      comments.guideline_3 = paste0(
        "10% of species potentially affected by the mixture. Corresponds to ",
        "the ANZG 90% species protection level."
      ),

      value.guideline_4 = .AMSPAF_CEILING * 100,
      level_name.guideline_4 = "Severely impacted",
      guideline.guideline_4 = "AmsPAF threshold",
      comments.guideline_4 = paste0(
        "20% of species potentially affected by the mixture. Corresponds to ",
        "the ANZG 80% species protection level. Values may exceed this; the ",
        "boundary exists for consistent reporting and colour-scaling."
      )
    )

  ## Propagate datetime.sample from input df to AmsPAF rows.
  if ("datetime.sample" %in% names(df)) {
    sample_times <- dplyr::distinct(df, uuid.sample, datetime.sample)
    amspaf_df    <- dplyr::left_join(amspaf_df, sample_times, by = "uuid.sample")
  }

  result <- dplyr::bind_rows(df, amspaf_df)
  if ("datetime.sample" %in% names(result)) {
    result <- dplyr::arrange(result, datetime.sample)
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
#' @return Tibble with columns `name.analyte`, `hc50`, `sigma`, `moa_group`,
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
      name.analyte     = nm,
      hc50             = hc50,
      sigma            = sigma,
      moa_group        = moa_group,
      parsed_formula   = list(parsed_f),
      coanalytes_req   = coanalytes_r
    )
  })

  if (is.null(params) || nrow(params) == 0L) .empty_ssd_params() else params
}

.empty_ssd_params <- function() {
  tibble::tibble(
    name.analyte   = character(),
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
#' @param ref_quantiles Tibble `(name.analyte, ref_norm)` from
#'   `prep_ref$normalised_quantiles`.
#' @param ssd_params Tibble from [derive_ssd_params()].
#' @param min_analytes Minimum analytes required.
#' @param method SSD method.
#' @param guideline_dir Path to ANZG XLSX folder.
#'
#' @return Tibble with one row per sample that passes `min_analytes`, columns:
#'   `uuid.sample`, `value`, `n_analytes_used`, `n_analytes_imputed`,
#'   `dominant_analyte`, `max_paf`, `analyte_pafs`.
#' @keywords internal
compute_amspaf_per_sample <- function(
    sample_data,
    ref_quantiles,
    ssd_params,
    min_analytes,
    method,
    guideline_dir
) {
  has_quantified <- "quantified" %in% names(sample_data)
  has_imputed    <- "imputed"    %in% names(sample_data)

  if (!has_quantified) {
    sample_data <- dplyr::mutate(sample_data, quantified = TRUE)
  }

  sample_data |>
    dplyr::group_by(uuid.sample) |>
    dplyr::group_modify(\(.x, .y) {
      ## Build co-analyte lookup (all quantified values in this sample)
      coanalyte_vals <- .x |>
        dplyr::filter(.data$quantified) |>
        dplyr::select("name.analyte", "value") |>
        tibble::deframe()  # named numeric vector

      ## Filter to SSD-eligible analytes
      tox_rows <- .x |>
        dplyr::filter(.data$name.analyte %in% ssd_params$name.analyte) |>
        dplyr::left_join(
          dplyr::select(ssd_params, "name.analyte", "hc50", "sigma",
                        "moa_group", "parsed_formula", "coanalytes_req"),
          by = "name.analyte"
        ) |>
        dplyr::left_join(ref_quantiles, by = "name.analyte") |>
        dplyr::mutate(
          ref_norm = tidyr::replace_na(.data$ref_norm, 0)
        )

      if (nrow(tox_rows) < min_analytes) {
        return(tibble::tibble(
          value              = numeric(0),
          n_analytes_used    = integer(0),
          n_analytes_imputed = integer(0),
          dominant_analyte   = character(0),
          max_paf            = numeric(0),
          analyte_pafs       = list()
        ))
      }

      ## ── Chemistry normalisation ────────────────────────────────────────
      ##
      ## BDL (quantified == FALSE): treated as zero exposure; C_norm = 0.
      ## Detected: apply normalisation formula using co-analyte values.
      ## Rows where normalisation returns NA (missing required co-analyte)
      ## are dropped and reported.
      tox_rows <- dplyr::mutate(
        tox_rows,
        C_norm = purrr::pmap_dbl(
          list(
            q   = .data$quantified,
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

      ## Drop rows where normalisation failed (NA C_norm)
      n_dropped_norm <- sum(is.na(tox_rows$C_norm))
      if (n_dropped_norm > 0L) {
        cli::cli_inform(c(
          "!" = "Sample {.val {unique(.y$uuid.sample)}}: {n_dropped_norm} analyte \\
                 row{?s} dropped (normalisation returned NA — missing co-analyte)."
        ), .frequency = "always", .frequency_id = paste0(unique(.y$uuid.sample), "_norm"))
        tox_rows <- dplyr::filter(tox_rows, !is.na(.data$C_norm))
      }

      if (nrow(tox_rows) < min_analytes) {
        return(tibble::tibble(
          value              = numeric(0),
          n_analytes_used    = integer(0),
          n_analytes_imputed = integer(0),
          dominant_analyte   = character(0),
          max_paf            = numeric(0),
          analyte_pafs       = list()
        ))
      }

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
          .data$C_adj, .data$name.analyte,
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
        tox_rows$name.analyte[which.max(tox_rows$PAF)]
      } else NA_character_

      tibble::tibble(
        uuid.sample        = dplyr::first(.x$uuid.sample),
        value              = amspaf * 100,
        n_analytes_used    = nrow(tox_rows),
        n_analytes_imputed = as.integer(n_analytes_imputed),
        dominant_analyte   = dominant,
        max_paf            = if (nrow(tox_rows) > 0L) max(tox_rows$PAF, na.rm = TRUE) else NA_real_,
        analyte_pafs       = list(
          dplyr::select(tox_rows, "name.analyte", "C_adj", "PAF", "moa_group")
        )
      )
    }) |>
    dplyr::ungroup()
}


## ============================================================================
## Utility: classify AmsPAF value into tier label
## ============================================================================

#' Classify an AmsPAF value into a reporting tier
#'
#' @param amspaf_pct Numeric vector of AmsPAF values in percent.
#' @return Character vector of tier labels.
#' @export
classify_amspaf_tier <- function(amspaf_pct) {
  checkmate::assert_numeric(amspaf_pct, lower = 0)
  dplyr::case_when(
    amspaf_pct < .AMSPAF_BREAK_T1_T2 * 100 ~ "1_background",
    amspaf_pct < .AMSPAF_BREAK_T2_T3 * 100 ~ "2_elevated",
    amspaf_pct < .AMSPAF_BREAK_T3_T4 * 100 ~ "3_impacted",
    TRUE ~ "4_severely_impacted"
  )
}
