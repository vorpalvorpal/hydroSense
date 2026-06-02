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
## `.parse_normalisation_formula()`. Bioavailability/index-condition formulas
## are populated for Cu, Ni, Zn, Cd, Pb and NH3-N; analytes with an empty cell
## are treated as identity (no chemistry dependence).
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
## (IA). Users can override the bundled classification by supplying their own
## metadata CSV/data frame (see `.load_analyte_metadata()`).
##
## The bundled `moa_group` values follow the toxic-mode-of-action (TMoA)
## vocabulary of De Zwart & Posthuma (2005, Table 3/4), assigned from the
## mechanism of action stated in each substance's ANZG/Warne guidance document:
##   * Metals & inorganic ions — grouped by the mechanism stated in each ANZG
##     guidance doc, since De Zwart & Posthuma (2005) give only a principle
##     (CA within a shared *primary target receptor*; RA across) plus a single
##     two-metal example, not a metals lookup. Gill-ionoregulatory metals
##     (Al, Cd, Cu, Ni, Pb, Zn) share `"ionoregulatory (gill)"` (one CA group);
##     Hg and Se share `"sulfhydryl binding"`; As is `"arsenic"` (speciation-
##     specific); Cr, B, Mn each stay solo (no mechanism stated in their docs);
##     `"ammonia"` and `"nitrate"` are each their own group. Groups combine
##     across by IA.
##   * Organics are grouped by shared mechanism and CA-combined within a group:
##     "nonpolar narcosis" (PAHs, BTEX, HCB), "polar narcosis" (phenol,
##     mono-chlorophenol), "uncoupler oxidative phosphorylation" (poly-
##     chlorophenols), "AChE inhibition: organophosphate" (OP insecticides),
##     "neurotoxicant: cyclodiene" (aldrin/dieldrin/endrin/lindane), and
##     "neurotoxicant: DDT" (DDT/DDE/methoxychlor).
##
## Any analyte with `moa_group = NA` gets a unique synthetic solo group
## `"_solo_<name>"` (combined by IA, i.e. treated as its own TMoA).
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

#' Assert every ammonia-bearing sample carries a water temperature
#'
#' The NH3-N un-ionised-fraction normalisation requires water temperature.
#' Rather than silently dropping ammonia when temperature is absent, we fail
#' loudly: any `sample_id` that has an `NH3-N` row but no non-missing
#' `temperature` row triggers an error listing the offending samples.
#'
#' @param df Long-format chemistry data frame (`sample_id`, `analyte`, `value`).
#' @return Invisibly `TRUE`; called for its side effect (error on violation).
#' @keywords internal
.assert_temperature_present <- function(df) {
  if (!"NH3-N" %in% df$analyte) return(invisible(TRUE))

  amm_samples <- unique(df$sample_id[df$analyte == "NH3-N"])

  temp_ok <- df$analyte == "temperature" & !is.na(df$value)
  temp_samples <- unique(df$sample_id[temp_ok])

  missing <- setdiff(amm_samples, temp_samples)
  if (length(missing) == 0L) return(invisible(TRUE))

  n_missing <- length(missing)
  shown <- missing[seq_len(min(10L, n_missing))]
  more  <- if (n_missing > 10L) paste0(" (+", n_missing - 10L, " more)") else ""
  cli::cli_abort(c(
    "{n_missing} sample{?s} report NH3-N but have no water temperature.",
    "x" = "Affected: {.val {shown}}{more}",
    "i" = "Water temperature is mandatory for the ammonia pH/temperature correction. Supply a {.field temperature} row per sample, or derive one with {.fn estimate_water_temp} (optionally fed by {.fn get_silo_air_temp}).",
    "i" = "To assess a dataset that does not include ammonia, call {.code add_amspaf(..., require_temperature = FALSE)}."
  ))
}

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
#' @param require_temperature Logical (default `TRUE`). When `TRUE`, any sample
#'   that reports an `NH3-N` measurement **must** also carry a water
#'   `temperature` row (the ammonia un-ionised-fraction normalisation is
#'   undefined without it); a missing temperature is a hard error rather than a
#'   silent drop of ammonia. Supply temperature via direct measurement, or
#'   derive it with [estimate_water_temp()] (optionally fed by
#'   [get_silo_air_temp()]). Set `FALSE` only for datasets that do not assess
#'   ammonia.
#'
#' @return The input `df` with AmsPAF rows appended. Each AmsPAF row carries:
#'   `value` (AmsPAF as a percentage, 0–100+), `detected = TRUE`,
#'   `analyte = "AmsPAF"`, `n_analytes_used` (integer),
#'   `n_analytes_imputed` (integer, 0 if `imputed` column absent),
#'   `dominant_analyte` (character), `max_paf` (numeric),
#'   `analyte_pafs` (list column of per-analyte diagnostic tibbles, each with
#'   `analyte`, `C_adj`, `PAF`, `moa_group`, and `ref_source` — one of
#'   `"disabled"`, `"matched"`, `"unmatched"` recording how the ARA reference
#'   was resolved for that analyte).
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
#' @examples
#' \dontrun{
#' # Long-format monitoring data: one row per sample x analyte.
#' obs <- tibble::tibble(
#'   sample_id = c("S1", "S1", "S1"),
#'   site_id   = "downstream",
#'   analyte   = c("Cu", "Zn", "temperature"),
#'   value     = c(3.2, 18, 19)
#' )
#' ref <- tibble::tibble(
#'   sample_id = "R1", site_id = "upstream",
#'   analyte   = c("Cu", "Zn"), value = c(1.1, 6)
#' )
#' options(leachatetools.guideline_dir = "path/to/guideline data")
#' add_amspaf(obs, reference = ref)
#' }
#' @export
add_amspaf <- function(
    df,
    reference        = NULL,
    analyte_metadata = NULL,
    method           = c("multi", "anzecc"),
    guideline_dir    = getOption("leachatetools.guideline_dir"),
    min_analytes     = 3,
    ref_summary      = c("geom_mean", "median", "arith_mean",
                          "p80", "p90", "p95"),
    require_temperature = TRUE
) {
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c("sample_id", "site_id", "analyte", "value")
  )
  method      <- match.arg(method)
  ref_summary <- match.arg(ref_summary)
  checkmate::assert_int(min_analytes, lower = 1L)
  checkmate::assert_flag(require_temperature)

  ## Temperature is mandatory for ammonia: the NH3-N un-ionised-fraction
  ## normalisation is undefined without it.  Fail loudly (rather than silently
  ## dropping ammonia) for any sample that reports NH3-N but lacks a water
  ## temperature.  See estimate_water_temp() / get_silo_air_temp() for sourcing
  ## temperature.
  if (require_temperature) .assert_temperature_present(df)

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

  ## ARA is "enabled" whenever the caller supplied any reference (a
  ## prepared_reference or a raw data frame).  Only `reference = NULL`
  ## disables it.  This flag lets compute_amspaf_one_sample() distinguish
  ## "ARA deliberately off" from "ARA on but this analyte had no reference
  ## match" — both otherwise collapse to ref_norm = 0 (see ref_source).
  ara_enabled <- !is.null(reference)

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
        guideline_dir = guideline_dir,
        ara_enabled   = ara_enabled
      )
    }) |>
    dplyr::ungroup()

  ## End-of-call summary of dropped analytes (single message rather than
  ## per-sample warnings — important for large datasets)
  .summarise_amspaf_diagnostics(amspaf_df, min_analytes)

  ## End-of-call summary of analytes assessed without a reference match
  ## (ARA enabled but no background available — assessed against raw conc).
  .summarise_ara_coverage(amspaf_df, ara_enabled)

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
#'   `coanalytes_req` (character), and `fit` (list column of fitted SSD
#'   objects, one per analyte).  The fitted SSD is loaded **once per analyte
#'   here** (via the cached [.load_or_fit()] path) so that
#'   [compute_amspaf_per_sample()] can evaluate every sample's PAF in a single
#'   vectorised [ssdtools::ssd_hp()] call per analyte, rather than refitting /
#'   re-resolving the SSD inside a per-(sample × analyte) loop.
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

    ## Load the fitted SSD object once (cached in-memory + on disk by
    ## .load_or_fit()).  Stored as a list column for batched PAF evaluation.
    ## NO3-N never reaches here (it is in .AMSPAF_EXCLUDED_ANALYTES), so no
    ## hardness-class resolution is needed.
    stem <- .SSD_NAME_MAP[[nm]]
    fit  <- if (!is.null(stem)) {
      tryCatch(.load_or_fit(nm, stem, method, guideline_dir),
               error = function(e) NULL)
    } else NULL

    ## MOA group from metadata column; NA/empty → unique solo group.
    ## A user-supplied metadata CSV may omit `moa_group` entirely; treat an
    ## absent column as all-NA so every analyte falls into its own solo group
    ## (i.e. pure Response Addition across analytes).
    mg_raw <- if ("moa_group" %in% names(eligible)) eligible$moa_group[i] else NA_character_
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
      coanalytes_req = coanalytes_r,
      fit            = list(fit)
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
    coanalytes_req = character(),
    fit            = list()
  )
}


## ============================================================================
## compute_ca_group_mspaf
## ============================================================================

#' Compute msPAF for a single Concentration Addition group
#'
#' Implements the multispecies Concentration Addition model of De Zwart &
#' Posthuma (2005, Environ Toxicol Chem 24(10):2665-2676), eq. 6: hazard units
#' (here `TU = C_adj / HC50`) are summed across the components sharing a toxic
#' mode of action, and the combined proportion affected is the log-normal CDF
#' evaluated at `log10(ΣTU)` with a mixture slope equal to the **arithmetic
#' mean of the component slopes** (their `β̄_TMoA`; in the normal-CDF form the
#' slope is the standard deviation `σ̄ = mean(σ)`):
#'
#' \deqn{msPAF_{CA} = \Phi\!\left( \frac{\log_{10}(\sum TU)}{\bar{\sigma}} \right)}
#'
#' Note this is a plain (unweighted) average of the per-analyte sigmas, per the
#' primary source — NOT a TU-weighted or variance-style combination.  A
#' single-component group therefore reduces exactly to that component's own SSD.
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

  # De Zwart & Posthuma (2005) eq. 6: mixture slope = mean of component slopes.
  sigma_mix <- mean(group_data$sigma)

  pnorm(log10(TU_mix) / sigma_mix)
}


## ============================================================================
## compute_amspaf_per_sample
## ============================================================================

#' Compute AmsPAF for each sample in a per-feature data block
#'
#' Internal workhorse called by [add_amspaf()] for each `site_id` group.  Runs
#' in three phases so the (relatively expensive) SSD evaluation is **batched
#' across samples** rather than called once per (sample × analyte):
#' \enumerate{
#'   \item per-sample chemistry normalisation + ARA shift ([.amspaf_adjust()]);
#'   \item one vectorised [ssdtools::ssd_hp()] call per analyte across every
#'     sample ([.amspaf_add_paf()]);
#'   \item per-sample CA/IA mixture combination ([.amspaf_combine()]).
#' }
#'
#' @param sample_data Per-feature long-format df (may include co-analyte rows
#'   such as pH, DOC alongside toxicant rows).
#' @param ref_table Tibble `(analyte, ref_norm)` from `prep_ref$ref_table`.
#' @param ssd_params Tibble from [derive_ssd_params()] (incl. the `fit` column).
#' @param min_analytes Minimum analytes required.
#' @param method SSD method (used only for the rare NULL-fit fallback).
#' @param guideline_dir Path to ANZG XLSX folder (NULL-fit fallback only).
#' @param ara_enabled Logical; whether the caller supplied a reference.
#'   Controls the `ref_source` diagnostic (see [.amspaf_adjust()]).
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
    guideline_dir,
    ara_enabled = TRUE
) {
  if (!"detected" %in% names(sample_data)) {
    sample_data <- dplyr::mutate(sample_data, detected = TRUE)
  }
  has_imputed <- "imputed" %in% names(sample_data)

  sample_ids <- unique(sample_data$sample_id)

  ## ── Phase 1: per-sample normalisation + ARA ──────────────────────────────
  adj_list <- lapply(sample_ids, function(sid) {
    rows <- dplyr::filter(sample_data, .data$sample_id == .env$sid)
    c(list(sample_id = sid),
      .amspaf_adjust(rows, ref_table, ssd_params, ara_enabled))
  })

  ## Samples that pass min_analytes (after dropping missing-co-analyte rows)
  keep <- purrr::keep(adj_list, function(z) nrow(z$tox) >= min_analytes)
  if (length(keep) == 0L) return(.amspaf_empty_row(with_sample_id = TRUE))

  tox_keep <- purrr::map_dfr(keep, function(z) {
    dplyr::mutate(z$tox, sample_id = z$sample_id)
  })

  ## ── Phase 2: batched PAF (one ssd_hp() call per analyte) ─────────────────
  tox_keep <- .amspaf_add_paf(tox_keep, ssd_params, method, guideline_dir)

  ## Per-sample dropped-analyte tibbles (diagnostic list-column)
  dropped_lookup <- stats::setNames(
    lapply(keep, `[[`, "dropped"),
    vapply(keep, `[[`, character(1), "sample_id")
  )

  ## ── Phase 3: per-sample CA/IA combination ────────────────────────────────
  tox_keep |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::group_modify(\(.x, .y) {
      res <- .amspaf_combine(.x, has_imputed)
      res$dropped_analytes <- list(dropped_lookup[[.y$sample_id]])
      res
    }) |>
    dplyr::ungroup()
}


## ============================================================================
## compute_amspaf_one_sample
## ============================================================================

#' Compute AmsPAF for a single sample
#'
#' Thin wrapper over the shared AmsPAF helpers ([.amspaf_adjust()],
#' [.amspaf_add_paf()], [.amspaf_combine()]) that processes one sample in
#' isolation.  Retained so the normalisation / ARA / CA / IA steps can be
#' driven (and unit-tested) for a single sample; the batched
#' [compute_amspaf_per_sample()] path uses the same helpers so behaviour is
#' guaranteed identical.
#'
#' @param sample_rows Long-format chemistry rows for one sample (one row per
#'   analyte; may include co-analyte rows used for normalisation). Must carry
#'   a `detected` column.
#' @param ref_table Tibble `(analyte, ref_norm)` from `prep_ref$ref_table`.
#' @param ssd_params Tibble from [derive_ssd_params()] (incl. the `fit` column).
#' @param min_analytes Minimum analytes required.
#' @param method SSD method (used only for the rare NULL-fit fallback).
#' @param guideline_dir Path to ANZG XLSX folder (NULL-fit fallback only).
#' @param has_imputed Logical; whether the input carried an `imputed` column
#'   (controls `n_analytes_imputed` accounting).
#' @param ara_enabled Logical; whether the caller supplied a reference.
#'   Controls the `ref_source` diagnostic (see [.amspaf_adjust()]).
#'
#' @return A one-row tibble (or zero-row tibble if the sample fails
#'   `min_analytes`) with columns `value`, `n_analytes_used`,
#'   `n_analytes_imputed`, `dominant_analyte`, `max_paf`, `analyte_pafs`,
#'   `dropped_analytes`.
#' @keywords internal
compute_amspaf_one_sample <- function(
    sample_rows,
    ref_table,
    ssd_params,
    min_analytes,
    method,
    guideline_dir,
    has_imputed = FALSE,
    ara_enabled = TRUE
) {
  if (!"detected" %in% names(sample_rows)) {
    sample_rows <- dplyr::mutate(sample_rows, detected = TRUE)
  }

  adj <- .amspaf_adjust(sample_rows, ref_table, ssd_params, ara_enabled)

  if (nrow(adj$tox) < min_analytes) return(.amspaf_empty_row())

  tox <- .amspaf_add_paf(adj$tox, ssd_params, method, guideline_dir)
  res <- .amspaf_combine(tox, has_imputed)
  res$dropped_analytes <- list(adj$dropped)
  res
}


## ============================================================================
## Shared AmsPAF helpers
## ============================================================================
##
## These power both compute_amspaf_one_sample() (single sample) and
## compute_amspaf_per_sample() (batched across a feature block).  Keeping the
## normalisation / ARA / PAF / combine logic in one place guarantees the two
## entry points cannot diverge.

#' Normalise chemistry and apply the ARA shift for a block of sample rows
#'
#' Filters to SSD-eligible analytes, joins SSD params + the reference table,
#' applies the per-analyte chemistry normalisation (BDL → 0), records the
#' `ref_source` diagnostic, and computes the added-risk-adjusted concentration
#' `C_adj = max(C_norm - ref_norm, 0)`.  Operates on one *or many* samples; the
#' caller is responsible for any per-sample grouping downstream.
#'
#' `ref_source` distinguishes:
#' \itemize{
#'   \item `"disabled"` — ARA off (no reference supplied to [add_amspaf()]);
#'   \item `"matched"` — ARA on and a reference value was found for the analyte;
#'   \item `"unmatched"` — ARA on but no reference match (assessed against raw
#'     normalised concentration, i.e. `ref_norm = 0`).
#' }
#'
#' @param sample_rows Long-format chemistry rows (one or more samples). Must
#'   carry `detected`.
#' @param ref_table Tibble `(analyte, ref_norm)`.
#' @param ssd_params Tibble from [derive_ssd_params()].
#' @param ara_enabled Logical; whether the caller supplied a reference.
#' @return A list `list(tox, dropped)` where `tox` has the eligible,
#'   normalisation-resolved rows (with `C_norm`, `ref_norm`, `ref_source`,
#'   `C_adj`) and `dropped` is a `(analyte, reason)` tibble of rows dropped for
#'   a missing required co-analyte.
#' @keywords internal
.amspaf_adjust <- function(sample_rows, ref_table, ssd_params,
                           ara_enabled = TRUE) {
  ## Co-analyte lookup (all detected values present in this block).  When the
  ## block holds multiple samples this still works because the per-analyte
  ## normalisation only ever reads a sample's own co-analyte rows — but to be
  ## safe the batched caller passes one sample at a time into this helper.
  coanalyte_vals <- sample_rows |>
    dplyr::filter(.data$detected) |>
    dplyr::select("analyte", "value") |>
    tibble::deframe()

  empty_dropped <- tibble::tibble(
    analyte = character(0), reason = character(0)
  )

  tox <- sample_rows |>
    dplyr::filter(.data$analyte %in% ssd_params$analyte) |>
    dplyr::left_join(
      dplyr::select(ssd_params, "analyte", "hc50", "sigma",
                    "moa_group", "parsed_formula", "coanalytes_req"),
      by = "analyte"
    ) |>
    dplyr::left_join(ref_table, by = "analyte")

  if (nrow(tox) == 0L) {
    return(list(
      tox = tibble::tibble(
        analyte = character(0), value = numeric(0), detected = logical(0),
        hc50 = numeric(0), sigma = numeric(0), moa_group = character(0),
        C_norm = numeric(0), ref_norm = numeric(0),
        ref_source = character(0), C_adj = numeric(0)
      ),
      dropped = empty_dropped
    ))
  }

  ## ref_source records WHY ref_norm is what it is (see roxygen above).  This
  ## must be computed before ref_norm's NA is coerced to 0, otherwise
  ## "disabled" and "unmatched" become indistinguishable.
  if (!"ref_norm" %in% names(tox)) tox$ref_norm <- NA_real_
  tox <- dplyr::mutate(
    tox,
    ref_source = dplyr::case_when(
      !ara_enabled        ~ "disabled",
      !is.na(.data$ref_norm) ~ "matched",
      TRUE                ~ "unmatched"
    ),
    ref_norm = tidyr::replace_na(.data$ref_norm, 0)
  )

  ## ── Chemistry normalisation (BDL → 0; missing co-analyte → NA → dropped) ──
  tox <- dplyr::mutate(
    tox,
    C_norm = purrr::pmap_dbl(
      list(
        q  = .data$detected,
        C  = .data$value,
        pf = .data$parsed_formula,
        cr = .data$coanalytes_req
      ),
      function(q, C, pf, cr) {
        if (!q) return(0)
        co_names <- if (nzchar(cr %||% "")) {
          trimws(strsplit(cr, ",")[[1L]])
        } else character(0)
        co_names <- co_names[nzchar(co_names)]
        co_vals  <- coanalyte_vals[co_names[co_names %in% names(coanalyte_vals)]]
        .apply_normalisation(pf, C, co_vals)
      }
    )
  )

  dropped <- dplyr::filter(tox, is.na(.data$C_norm)) |>
    dplyr::transmute(.data$analyte, reason = "missing_co_analyte")

  tox <- tox |>
    dplyr::filter(!is.na(.data$C_norm)) |>
    dplyr::mutate(C_adj = pmax(.data$C_norm - .data$ref_norm, 0))

  list(tox = tox, dropped = dropped)
}

#' Add per-analyte PAF to adjusted tox rows, batched per analyte
#'
#' Evaluates each analyte's SSD once across **all** its `C_adj` values in a
#' single vectorised [ssdtools::ssd_hp()] call (via [.ssd_paf_vec()]), rather
#' than one [ssd_paf()] call per row.  The fitted SSD object is taken from the
#' `fit` list-column of `ssd_params` (loaded once in [derive_ssd_params()]).
#'
#' @param tox Adjusted tox rows from [.amspaf_adjust()] (needs `analyte`,
#'   `C_adj`); may span multiple samples.
#' @param ssd_params Tibble from [derive_ssd_params()] (uses `analyte`, `fit`).
#' @param method SSD method (NULL-fit fallback only).
#' @param guideline_dir Path to ANZG XLSX folder (NULL-fit fallback only).
#' @return `tox` with a numeric `PAF` column added (proportion, 0–1).
#' @keywords internal
.amspaf_add_paf <- function(tox, ssd_params, method, guideline_dir) {
  tox$PAF <- NA_real_
  if (nrow(tox) == 0L) return(tox)

  fit_lookup <- stats::setNames(ssd_params$fit, ssd_params$analyte)

  for (a in unique(tox$analyte)) {
    idx <- which(tox$analyte == a)
    tox$PAF[idx] <- .ssd_paf_vec(
      fit           = fit_lookup[[a]],
      conc          = tox$C_adj[idx],
      analyte       = a,
      method        = method,
      guideline_dir = guideline_dir
    )
  }
  tox
}

#' Vectorised SSD PAF lookup for one analyte
#'
#' Returns the proportion of species affected at each concentration in `conc`.
#' Concentrations that are `NA` or `<= 0` map to `PAF = 0`.  When a fitted SSD
#' object is available it is evaluated in a single [ssdtools::ssd_hp()] call
#' over the positive concentrations; otherwise it falls back to a per-value
#' [ssd_paf()] lookup (which re-resolves the model internally).
#'
#' @param fit Fitted SSD object (from the `fit` list-column) or `NULL`.
#' @param conc Numeric vector of adjusted concentrations.
#' @param analyte Analyte name (fallback `ssd_paf()` lookup only).
#' @param method SSD method (fallback only).
#' @param guideline_dir Path to ANZG XLSX folder (fallback only).
#' @return Numeric vector the same length as `conc`, PAF as a proportion (0–1).
#' @keywords internal
.ssd_paf_vec <- function(fit, conc, analyte, method, guideline_dir) {
  out <- numeric(length(conc))
  pos <- which(!is.na(conc) & conc > 0)
  if (length(pos) == 0L) return(out)

  if (is.null(fit)) {
    ## Fallback: scalar ssd_paf() per positive concentration.
    out[pos] <- vapply(conc[pos], function(c) {
      paf_result <- tryCatch(
        ssd_paf(analyte, c, method = method,
                guideline_dir = guideline_dir, nboot = 0L),
        error = function(e) list(pct = NA_real_)
      )
      if (is.na(paf_result$pct)) NA_real_ else paf_result$pct / 100
    }, numeric(1))
    return(out)
  }

  ## ssdtools::ssd_hp() expands duplicated `conc` values into a cross-join
  ## (e.g. 2 identical concentrations return 4 rows), so we cannot rely on the
  ## result being one row per input.  Evaluate on the UNIQUE concentrations and
  ## map back by value.  est is constant for a given conc (model-averaged), so
  ## the first row per conc is authoritative.
  uc <- unique(conc[pos])
  res <- tryCatch(
    ssdtools::ssd_hp(fit, conc = uc, ci = FALSE, proportion = TRUE),
    error = function(e) NULL
  )
  if (is.null(res)) {
    out[pos] <- NA_real_
    return(out)
  }
  res <- res[!duplicated(res$conc), c("conc", "est"), drop = FALSE]
  out[pos] <- res$est[match(conc[pos], res$conc)]
  out
}

#' Combine adjusted+PAF'd tox rows for one sample into a single AmsPAF row
#'
#' Performs Concentration Addition within each mode-of-action group and
#' Independent Action across groups, then assembles the one-row diagnostic
#' tibble.  Operates on a *single* sample's rows.
#'
#' @param tox Rows for one sample from [.amspaf_add_paf()] (needs `C_adj`,
#'   `hc50`, `sigma`, `moa_group`, `PAF`; `analyte`, `ref_source` carried into
#'   the `analyte_pafs` diagnostic if present).
#' @param has_imputed Logical; whether the input carried an `imputed` column.
#' @return A one-row tibble: `value`, `n_analytes_used`, `n_analytes_imputed`,
#'   `dominant_analyte`, `max_paf`, `analyte_pafs` (does NOT add `sample_id` or
#'   `dropped_analytes` — the caller attaches those).
#' @keywords internal
.amspaf_combine <- function(tox, has_imputed = FALSE) {
  n_analytes_imputed <- if (has_imputed && "imputed" %in% names(tox)) {
    sum(tox$imputed, na.rm = TRUE)
  } else 0L

  groups         <- unique(tox$moa_group)
  mspaf_by_group <- vapply(
    groups,
    function(g) compute_ca_group_mspaf(dplyr::filter(tox, .data$moa_group == g)),
    numeric(1)
  )
  amspaf <- 1 - prod(1 - mspaf_by_group)

  dominant <- if (any(!is.na(tox$PAF))) {
    tox$analyte[which.max(tox$PAF)]
  } else NA_character_

  pafs_cols <- intersect(
    c("analyte", "C_adj", "PAF", "moa_group", "ref_source"), names(tox)
  )

  tibble::tibble(
    value              = amspaf * 100,
    n_analytes_used    = nrow(tox),
    n_analytes_imputed = as.integer(n_analytes_imputed),
    dominant_analyte   = dominant,
    max_paf            = if (nrow(tox) > 0L) max(tox$PAF, na.rm = TRUE) else NA_real_,
    analyte_pafs       = list(dplyr::select(tox, dplyr::all_of(pafs_cols)))
  )
}

#' Construct an empty (zero-row) AmsPAF result tibble
#'
#' @param with_sample_id Logical; if `TRUE`, prepend a `sample_id` column (used
#'   by the batched [compute_amspaf_per_sample()] path).
#' @return A zero-row tibble with the AmsPAF output schema.
#' @keywords internal
.amspaf_empty_row <- function(with_sample_id = FALSE) {
  out <- tibble::tibble(
    value              = numeric(0),
    n_analytes_used    = integer(0),
    n_analytes_imputed = integer(0),
    dominant_analyte   = character(0),
    max_paf            = numeric(0),
    analyte_pafs       = list(),
    dropped_analytes   = list()
  )
  if (with_sample_id) {
    out <- tibble::add_column(out, sample_id = character(0), .before = 1L)
  }
  out
}

#' Summarise analytes assessed without a reference match (ARA coverage)
#'
#' When ARA is enabled but some analytes had no matching reference value, those
#' analytes were assessed against their raw normalised concentration
#' (`ref_norm = 0`).  This emits a single end-of-call cli message tallying how
#' many samples each such analyte affected, so the caller knows where the ARA
#' adjustment did *not* apply.  No-op when ARA is disabled.
#'
#' @param amspaf_df The assembled AmsPAF tibble (must carry `analyte_pafs`).
#' @param ara_enabled Logical; whether the caller supplied a reference.
#' @keywords internal
.summarise_ara_coverage <- function(amspaf_df, ara_enabled) {
  if (!ara_enabled) return(invisible(NULL))
  if (!"analyte_pafs" %in% names(amspaf_df) || nrow(amspaf_df) == 0L)
    return(invisible(NULL))

  pafs_long <- tryCatch(
    tidyr::unnest(
      dplyr::select(amspaf_df, dplyr::any_of(c("sample_id", "site_id")),
                    "analyte_pafs"),
      cols = "analyte_pafs"
    ),
    error = function(e) NULL
  )
  if (is.null(pafs_long) || !"ref_source" %in% names(pafs_long))
    return(invisible(NULL))

  unmatched <- dplyr::filter(pafs_long, .data$ref_source == "unmatched")
  if (nrow(unmatched) == 0L) return(invisible(NULL))

  per_analyte <- unmatched |>
    dplyr::count(.data$analyte, name = "n_samples")

  n_analytes <- nrow(per_analyte)

  cli::cli_inform(c(
    "i" = "add_amspaf: {n_analytes} analyte{?s} had no reference match \\
           (ARA enabled) and were assessed against raw normalised \\
           concentration:",
    paste0("    ", per_analyte$analyte, ": ", per_analyte$n_samples,
           " sample", ifelse(per_analyte$n_samples == 1L, "", "s"))
  ))
  invisible(per_analyte)
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


