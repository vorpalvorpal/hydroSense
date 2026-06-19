## ============================================================================
## Adjusted multi-substance Potentially Affected Fraction (msPAF)
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
## as the single-substance PAF calculations.
##
## ## The log-normal assumption in the CA step
##
## The package uses an SSD in two ways that make different distributional
## commitments — worth being explicit about:
##
##   * **Per-analyte PAF and the IA combination are shape-faithful.** ssd_paf()
##     reads each substance's PAF off the model-averaged "multi" curve (all 6
##     BCANZ distributions, weighted), so single-substance PAFs and the
##     across-group IA step (msPAF = 1 - prod(1 - PAF_i)) hold for any SSD
##     shape — IA needs only each component PAF, not a parametric form.
##   * **The within-group CA step is intrinsically parametric.** Concentration
##     Addition sums toxic units and reads the combined PAF off a SINGLE
##     log-normal CDF, msPAF_CA = Phi(log10(sum TU_i) / sigma_bar), with
##     sigma_bar the mean of the component SSD slopes (De Zwart & Posthuma
##     2005, eq. 6). CA treats co-acting substances as dilutions of one
##     "super-chemical", so it must express each as HC50 + a slope. This is the
##     standard msPAF-CA method and is EXACT when the component SSDs are
##     log-normal (or log-logistic, which has the same closed form); for a
##     model-averaged curve that is not itself log-normal it is an
##     APPROXIMATION, matched at HC50/HC05. compute_ca_group_mspaf() derives
##     each slope analytically from the fitted SSD's HC50 and HC05.
##
## In practice the approximation is small for the metals here (their SSDs are
## close to log-normal), but it is a real modelling choice, not an identity. A
## fully shape-faithful numerical CA (Monte-Carlo over species sensitivities
## under full rank-concordance, valid for any SSD shape) is planned for v2.0 as
## a selectable mspaf_method; see TODO sec. 7.
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
## msPAF is returned as a continuous risk metric (% of species potentially
## affected by the mixture, 0-100+).  Tier breaks are deliberately not
## provided by this package — msPAF semantics differ from single-substance
## guideline values (the "% affected" denominator depends on which
## substances are in the mixture).  See
## `vignette("chronic-mspaf-interpretation")` for a discussion of how to
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

## Session cache for SSD PAF spline closures. Keys are "method/analyte"
## (for shipped tables) or "method/analyte/guideline_dir" (for runtime builds).
.ssd_paf_lookup_env <- new.env(parent = emptyenv())

## Analytes excluded from msPAF regardless of SSD availability.
##
##   NO3-N: GVs from NZ document of uncertain provenance, not ANZG
##   CH4:   Methane guideline is aesthetic/nuisance, not an ecotoxicity SSD
##   LHF:   Leachate indicator index — derived values, not a toxicity endpoint
##
.MSPAF_EXCLUDED_ANALYTES <- c("NO3-N", "CH4", "LHF")


## ============================================================================
## add_mspaf
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
    "i" = "To assess a dataset that does not include ammonia, call {.code add_mspaf(..., require_temperature = FALSE)}."
  ))
}

#' Compute the Adjusted multi-substance PAF (msPAF) for water quality samples
#'
#' Appends msPAF rows to a long-format water quality dataframe. msPAF
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
#'   `sample_id`, `site_id`, `analyte`, `value`. Toxicant concentrations must
#'   ultimately be in µg/L for SSD lookup; supply them either via a
#'   `units.analyte` column (one unit string per row, e.g. `"mg/L"`) or via the
#'   `conc_units` argument (applied uniformly to all SSD-eligible rows). Optional
#'   but recommended: `datetime` (propagated to msPAF rows if present),
#'   `detected` (assumed `TRUE` if absent), `imputed` (logical; if present,
#'   `n_analytes_imputed` is populated in output). Driver analytes needed for
#'   chemistry normalisation (e.g. pH, DOC) should be present as rows in `df`.
#' @param reference Background reference chemistry for the ARA adjustment.
#'   Accepts four forms:
#'   \itemize{
#'     \item A `reference_model` from [fit_reference_model()] — contemporaneous
#'       (temporal) ARA; the model predicts what the reference site would show
#'       at each target sample's exact moment using hydrology and seasonality.
#'     \item A `prepared_reference` object from [prepare_reference()] —
#'       normalisation has already been applied; used directly (static ARA).
#'     \item A raw long-format data frame (same schema as `df`) — will be
#'       passed to [prepare_reference()] internally (static ARA).
#'     \item `NULL` (default) — no ARA adjustment; raw concentrations assessed
#'       directly against SSDs.
#'   }
#' @param tau,tau_units Exponential-decay half-life for chronic window
#'   integration when `reference` is a `reference_model`.  Default `90` days.
#' @param window,window_units Look-back window length for chronic integration.
#'   Default `365` days.
#' @param analyte_metadata Data frame of analyte metadata, or `NULL` to load
#'   the bundled `inst/extdata/anzecc_analyte_metadata.csv`. Passed to
#'   [prepare_reference()] and [derive_ssd_params()].
#' @param method SSD method passed to [ssd_hc50()] and [ssd_paf()].
#'   `"multi"` (default) fits all 6 BCANZ distributions and model-averages;
#'   `"anzecc"` uses the per-analyte distribution matching the original ANZG
#'   derivation.
#' @param guideline_dir Path to the "guideline data" folder containing ANZG
#'   XLSX files. Falls back to `getOption("hydroSense.guideline_dir")`.
#' @param min_analytes Minimum number of analytes with fitted SSDs required
#'   to compute msPAF for a sample. Default `3`.
#' @param ref_summary Summary statistic for the reference distribution when
#'   `reference` is a raw data frame.  Passed through to
#'   [prepare_reference()].  Default `"geom_mean"` — the maximum-likelihood
#'   central tendency for log-normal concentrations, and a PICT-consistent
#'   estimate of the "typical" exposure the resident community has adapted
#'   to.  Other options: `"median"`, `"arith_mean"`, `"p80"`, `"p90"`,
#'   `"p95"`.
#' @param conc_units Character. Unit string (e.g. `"mg/L"`, `"ug/L"`) applied
#'   uniformly to all SSD-eligible rows in `df` when `df` has no
#'   `units.analyte` column. Ignored when `df` already carries `units.analyte`.
#'   Required when `df` lacks `units.analyte` and toxicant concentrations are
#'   not already in µg/L.
#' @param require_temperature Logical (default `TRUE`). When `TRUE`, any sample
#'   that reports an `NH3-N` measurement **must** also carry a water
#'   `temperature` row (the ammonia un-ionised-fraction normalisation is
#'   undefined without it); a missing temperature is a hard error rather than a
#'   silent drop of ammonia. Supply temperature via direct measurement, or
#'   derive it with [estimate_water_temp()] (optionally fed by
#'   [get_silo_air_temp()]). Set `FALSE` only for datasets that do not assess
#'   ammonia.
#' @param return Output mode for draw-carrier input (see [summarise_draws()]).
#'   `"summary"` (default) collapses posterior draws to a central estimate plus
#'   a credible interval (`value`, `value_lower`, `value_upper`, `n_draws`);
#'   `"draws"` returns the raw per-draw msPAF rows (`draw_id 1..N`) for
#'   external risk models or further composition (e.g. into
#'   [time_weighted_aggregate()]).  For point (non-draw) input both modes
#'   return byte-identical output with no interval columns.
#' @param interval Width of the credible interval when `return = "summary"`.
#'   Default `0.90` (5th/95th percentile bounds).
#' @param central Central-tendency statistic when `return = "summary"`:
#'   `"median"` (default) or `"mean"`.
#'
#' @return The input `df` with msPAF rows appended. Each msPAF row carries:
#'   `value` (msPAF as a percentage, 0–100+), `detected = TRUE`,
#'   `analyte = "msPAF"`, `n_analytes_used` (integer),
#'   `n_analytes_imputed` (integer, 0 if `imputed` column absent),
#'   `dominant_analyte` (character), and `max_paf` (numeric).
#'
#'   The result carries two attributes (read them before further dplyr wrangling,
#'   which drops attributes): `"analyte_pafs"` — the per-analyte PAF breakdown as
#'   a flat tibble (`site_id`, `sample_id`, `draw_id` in draws mode, `analyte`,
#'   `C_adj`, `PAF`, `moa_group`, `ref_source`), retrieved with [analyte_pafs()];
#'   and `"ara_summary"` — per-(sample × analyte) ARA diagnostics, retrieved with
#'   [ara_summary()]. (`analyte_pafs` was formerly a per-row list-column; it is
#'   now a flat attribute — issue #30.)
#'
#'   Tier breaks are not provided by this package — msPAF is a continuous
#'   risk metric and the threshold at which a community is "impacted"
#'   depends on the assessment context.  See
#'   `vignette("chronic-mspaf-interpretation")` for guidance.
#'
#' @seealso [ssd_paf()], [ssd_hc50()], [prepare_reference()],
#'   [fit_reference_model()], [ara_summary()],
#'   [time_weighted_aggregate()], [prescreen_analytes()], [impute_chemistry()]
#'
#' @references
#' De Zwart D, Posthuma L (2005) Environmental Toxicology and Chemistry
#' 24(10):2665-2676.
#'
#' @examples
#' \donttest{
#' # Per-sample multi-substance PAF for the impacted site, with local
#' # background subtracted via the reference site. Uses the bundled SSD data.
#' demo <- leachate_demo()
#' ds  <- subset(demo, site_id == "downstream")
#' ref <- subset(demo, site_id == "reference")
#' out <- add_mspaf(ds, reference = ref)
#' subset(out, analyte == "msPAF", c("sample_id", "value"))
#' }
#' @export
add_mspaf <- function(
    df,
    reference        = NULL,
    analyte_metadata = NULL,
    method           = c("multi", "anzecc"),
    guideline_dir    = getOption("hydroSense.guideline_dir"),
    min_analytes     = 3,
    ref_summary      = c("geom_mean", "median", "arith_mean",
                          "p80", "p90", "p95"),
    conc_units       = NULL,
    require_temperature = TRUE,
    tau              = 90,
    tau_units        = "day",
    window           = 365,
    window_units     = "day",
    return           = c("summary", "draws"),
    interval         = 0.90,
    central          = c("median", "mean")
) {
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c("sample_id", "site_id", "analyte", "value")
  )
  method      <- match.arg(method)
  ref_summary <- match.arg(ref_summary)
  return      <- match.arg(return)
  central     <- match.arg(central)
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
      "No analytes with fitted SSDs found. msPAF cannot be computed."
    )
    return(df)
  }

  ## Convert toxicant concentrations to µg/L (SSD expectation).  Co-analyte
  ## rows (pH, DOC, hardness, temperature) are left untouched because their
  ## normalisation formulas reference them in their natural units.
  df <- .convert_df_tox_to_ugL(df, ssd_params$analyte, conc_units, "df")

  ## ================================================================
  ## Step 2: Resolve the reference into a prepared_reference object.
  ## ================================================================

  ## ARA is "enabled" whenever the caller supplied any reference (a
  ## prepared_reference or a raw data frame).  Only `reference = NULL`
  ## disables it.  This flag lets .mspaf_adjust() distinguish "ARA deliberately
  ## off" from "ARA on but this analyte had no reference match" — both otherwise
  ## collapse to ref_norm = 0 (see ref_source).
  ara_enabled <- !is.null(reference)

  if (inherits(reference, "reference_model")) {
    ## Temporal ARA path: resolve contemporaneous reference norms per sample
    tau_days    <- .resolve_to(tau,    "day", tau_units,    "tau")
    window_days <- .resolve_to(window, "day", window_units, "window")
    ref_table   <- .resolve_ref_norm(reference, df, tau_days, window_days)
    ## ref_table has columns: sample_id, analyte, ref_norm, ref_tier
  } else if (inherits(reference, "prepared_reference")) {
    prep_ref  <- reference
    ## Use only analyte+ref_norm columns for the join (avoid polluting tox_rows
    ## with n_obs / ref_lower / ref_upper if bootstrap_ci was used)
    ref_table <- if (nrow(prep_ref$ref_table) > 0L) {
      dplyr::select(prep_ref$ref_table, "analyte", "ref_norm")
    } else {
      prep_ref$ref_table
    }
  } else if (is.null(reference)) {
    ## No ARA: empty ref_table → ref_norm joined as NA → coerced to 0
    ref_table <- tibble::tibble(
      analyte  = character(0),
      ref_norm = numeric(0)
    )
  } else {
    checkmate::assert_data_frame(reference)
    prep_ref  <- prepare_reference(
      reference,
      analyte_metadata = meta,
      summary          = ref_summary
    )
    ref_table <- if (nrow(prep_ref$ref_table) > 0L) {
      dplyr::select(prep_ref$ref_table, "analyte", "ref_norm")
    } else {
      prep_ref$ref_table
    }
  }

  ## ================================================================
  ## Step 3: Compute msPAF per feature, per sample.
  ## ================================================================

  ## Per-site engine via split() + lapply (not group_modify): the per-site flat
  ## diagnostic attributes (analyte_pafs / dropped / ara_diag) would be stripped
  ## by group_modify's row-binding, so collect them explicitly. (issue #30)
  site_ids <- sort(unique(df$site_id))
  per_site <- lapply(site_ids, function(sid) {
    r <- compute_mspaf_per_sample(
      sample_data   = df[df$site_id == sid, , drop = FALSE],
      ref_table     = ref_table,
      ssd_params    = ssd_params,
      min_analytes  = min_analytes,
      method        = method,
      guideline_dir = guideline_dir,
      ara_enabled   = ara_enabled
    )
    r$site_id <- sid
    r
  })
  mspaf_df <- dplyr::bind_rows(per_site)

  collect_attr <- function(nm) {
    parts <- Map(function(r, sid) {
      a <- attr(r, nm)
      if (is.null(a) || nrow(a) == 0L) return(NULL)
      a$site_id <- sid
      a
    }, per_site, site_ids)
    dplyr::bind_rows(Filter(Negate(is.null), parts))
  }
  analyte_pafs_flat <- collect_attr("analyte_pafs")   # full per-(sample,draw)
  dropped_flat      <- collect_attr("dropped")
  ara_diag_flat     <- collect_attr("ara_diag")

  ## Structural diagnostics (dropped analytes, ref_source) are identical across
  ## draws; restrict to a representative draw for the end-of-call summaries and
  ## the ara_summary attribute (avoids ×N over-counting).
  rep_draw <- function(flat) {
    if (is.null(flat) || nrow(flat) == 0L || !"draw_id" %in% names(flat))
      return(flat)
    dplyr::filter(flat, .data$draw_id == min(flat$draw_id, na.rm = TRUE))
  }
  .summarise_mspaf_diagnostics(rep_draw(dropped_flat))
  .summarise_ara_coverage(rep_draw(analyte_pafs_flat), ara_enabled)

  ## ara_summary attribute: per-cell ARA diagnostics for the representative draw
  ## (sample_id + diagnostic columns, matching the historical shape).
  ara_diag_rep <- rep_draw(ara_diag_flat)
  ara_summary_out <- if (!is.null(ara_diag_rep) && nrow(ara_diag_rep) > 0L) {
    dplyr::select(ara_diag_rep, "sample_id", "analyte", "ref_norm", "C_norm",
                  "C_adj", "C_excess", "floor_fired", "ref_source", "ref_tier")
  } else NULL

  mspaf_df <- dplyr::mutate(mspaf_df,
    analyte  = "msPAF",
    detected = TRUE
  )

  ## Tier breaks / regulatory interpretation are intentionally NOT provided
  ## by this package.  msPAF is a continuous risk metric (% of species
  ## potentially affected); the threshold at which a community is
  ## "impacted" depends on the assessment context (mixture composition,
  ## site-specific calibration, target protection level) and is a
  ## consumer-side decision.  See vignette("chronic-mspaf-interpretation").

  ## Propagate datetime from input df to msPAF rows (per-sample pipeline).
  if ("datetime" %in% names(df)) {
    sample_times <- dplyr::distinct(df, .data$sample_id, .data$datetime)
    mspaf_df    <- dplyr::left_join(mspaf_df, sample_times, by = "sample_id")
  }

  ## Propagate focal_date from input df to msPAF rows (chronic pipeline).
  if ("focal_date" %in% names(df)) {
    focal_times <- dplyr::distinct(df, .data$sample_id, .data$focal_date)
    mspaf_df   <- dplyr::left_join(mspaf_df, focal_times, by = "sample_id")
  }

  result <- dplyr::bind_rows(df, mspaf_df)
  if ("focal_date" %in% names(result)) {
    result <- dplyr::arrange(result, .data$focal_date)
  } else if ("datetime" %in% names(result)) {
    result <- dplyr::arrange(result, .data$datetime)
  }

  ## Collapse draws to posterior median + CI (default), or return raw draws.
  ## Point input: summarise_draws is a no-op (identity).
  if (return == "summary") result <- summarise_draws(result, interval, central)

  ## Store ARA cell-level diagnostics as an attribute for ara_summary(), and the
  ## per-analyte PAF breakdown (full per-(sample, draw)) for analyte_pafs().
  attr(result, "ara_summary")  <- ara_summary_out
  attr(result, "analyte_pafs") <- analyte_pafs_flat

  result
}


## ============================================================================
## derive_ssd_params
## ============================================================================

#' Derive SSD parameters for msPAF computation
#'
#' Reads analyte eligibility and mode-of-action groups from the bundled
#' metadata CSV, then calls [ssd_hc50()] and `.ssd_sigma()` to populate the
#' HC50 and effective sigma needed for Concentration Addition. Chemistry
#' normalisation formulas (currently all stubs) are parsed once and stored
#' as a list column for use in [compute_mspaf_per_sample()].
#'
#' Eligibility criteria:
#' \itemize{
#'   \item `ssd_available == TRUE` in the metadata
#'   \item `analyte` not in `.MSPAF_EXCLUDED_ANALYTES`
#'   \item `ssd_hc50()` returns a non-NA value (model fits successfully)
#' }
#'
#' @param meta Analyte metadata tibble from `.load_analyte_metadata()`.
#' @param method SSD method: `"multi"` or `"anzecc"`.
#' @param guideline_dir Path to ANZG guideline data folder.
#'
#' @return Tibble with columns `analyte`, `hc50`, `sigma`, `moa_group`,
#'   `parsed_formula` (list of language objects or NULLs),
#'   `coanalytes_req` (character), and `fit` (list column of fitted SSD
#'   objects, one per analyte).  The fitted SSD is loaded **once per analyte
#'   here** (via the cached `.load_or_fit()` path) so that
#'   [compute_mspaf_per_sample()] can evaluate every sample's PAF in a single
#'   vectorised [ssdtools::ssd_hp()] call per analyte, rather than refitting /
#'   re-resolving the SSD inside a per-(sample × analyte) loop.
#'
#' @keywords internal
derive_ssd_params <- function(meta, method, guideline_dir) {
  eligible <- meta |>
    dplyr::filter(.data$ssd_available == TRUE) |>
    dplyr::filter(!.data$analyte %in% .MSPAF_EXCLUDED_ANALYTES)

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
    ## NO3-N never reaches here (it is in .MSPAF_EXCLUDED_ANALYTES), so no
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

  stats::pnorm(log10(TU_mix) / sigma_mix)
}


## ============================================================================
## compute_mspaf_per_sample
## ============================================================================

#' Compute msPAF for each sample in a per-feature data block
#'
#' Internal workhorse called by [add_mspaf()] for each `site_id` group.  Runs
#' in three phases so the (relatively expensive) SSD evaluation is **batched
#' across samples** rather than called once per (sample × analyte):
#' \enumerate{
#'   \item chemistry normalisation + ARA shift, vectorised across all rows
#'     ([.mspaf_adjust()]);
#'   \item one vectorised [ssdtools::ssd_hp()] call per analyte across every
#'     sample ([.mspaf_add_paf()]);
#'   \item CA/IA mixture combination via grouped reductions over
#'     (sample, draw, MOA group).
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
#'   Controls the `ref_source` diagnostic (see [.mspaf_adjust()]).
#'
#' @return Tibble with one row per sample that passes `min_analytes`, columns:
#'   `sample_id`, `value`, `n_analytes_used`, `n_analytes_imputed`,
#'   `dominant_analyte`, `max_paf`, `analyte_pafs`, `dropped_analytes`.
#' @keywords internal
compute_mspaf_per_sample <- function(
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

  ## Draw-carrier: broadcast exact cells so every row has a concrete draw_id.
  ## In the point case (no draw_id column, or all-NA) this assigns draw_id=1L
  ## everywhere; is_draws_mode=FALSE strips draw_id from the output at the end so
  ## the returned schema is byte-identical to pre-draws behaviour.
  is_draws_mode <- "draw_id" %in% names(sample_data) &&
    !all(is.na(sample_data[["draw_id"]]))
  draws       <- .draw_domain(sample_data)
  sample_data <- .broadcast_draws(sample_data, draws)

  ## ── Phase 1: normalisation + ARA, vectorised across the WHOLE frame ────────
  ## .mspaf_adjust() handles all (sample, draw) blocks in one pass (co-analytes
  ## joined wide; each analyte's formula evaluated once across its rows), so
  ## there is no per-block dplyr loop. (issue #30)
  adj <- .mspaf_adjust(sample_data, ref_table, ssd_params, ara_enabled)
  tox <- adj$tox
  tox$imp_n <- if (has_imputed && "imputed" %in% names(tox))
    as.integer(tox$imputed) else 0L

  ## Blocks (sample_id, draw_id) passing min_analytes after missing-co drops.
  n_by <- dplyr::count(tox, .data$sample_id, .data$draw_id, name = "nblk")
  keep_keys <- n_by |>
    dplyr::filter(.data$nblk >= min_analytes) |>
    dplyr::select("sample_id", "draw_id")
  if (nrow(keep_keys) == 0L) {
    empty <- .mspaf_empty_row(with_sample_id = TRUE)
    if (!is_draws_mode) empty <- dplyr::select(empty, -"draw_id")
    return(empty)
  }
  tox <- dplyr::semi_join(tox, keep_keys, by = c("sample_id", "draw_id"))

  ## ── Phase 2: batched PAF — one ssd_hp() per analyte across all rows. ──────
  tox <- .mspaf_add_paf(tox, ssd_params, method, guideline_dir)

  ## ── Phase 3: vectorised CA (per sample x draw x MOA group) then IA. ───────
  ## De Zwart & Posthuma (2005, Integr. Environ. Assess. Manag. 1:e1, eq.6):
  ## concentration addition within an MOA group (mixture slope = mean component
  ## slope; group msPAF = pnorm(log10(sum TU)/sigma_mix)); independent action
  ## across groups: msPAF = 1 - prod(1 - msPAF_group).
  ## NB: rows with NA moa_group contribute 0 here, reproducing the previous
  ## per-group behaviour (`filter(moa_group == NA)` selected nothing) — see the
  ## implementation note raised on issue #30.
  ca <- tox |>
    dplyr::filter(.data$C_adj > 0, !is.na(.data$moa_group),
                  is.finite(.data$hc50), .data$hc50 > 0,
                  is.finite(.data$sigma), .data$sigma > 0) |>
    dplyr::mutate(TU = .data$C_adj / .data$hc50) |>
    dplyr::group_by(.data$sample_id, .data$draw_id, .data$moa_group) |>
    dplyr::summarise(TU_mix = sum(.data$TU), sigma_mix = mean(.data$sigma),
                     .groups = "drop") |>
    dplyr::mutate(msPAF = ifelse(.data$TU_mix > 0,
                    stats::pnorm(log10(.data$TU_mix) / .data$sigma_mix), 0))
  ia <- ca |>
    dplyr::group_by(.data$sample_id, .data$draw_id) |>
    dplyr::summarise(mspaf = 1 - prod(1 - .data$msPAF), .groups = "drop")

  ## Per-block scalars over the kept (post-drop) tox rows.
  scal <- tox |>
    dplyr::group_by(.data$sample_id, .data$draw_id) |>
    dplyr::summarise(
      n_analytes_used    = dplyr::n(),
      n_analytes_imputed = as.integer(sum(.data$imp_n)),
      dominant_analyte   = if (any(!is.na(.data$PAF)))
                             .data$analyte[which.max(.data$PAF)] else NA_character_,
      max_paf            = max(.data$PAF, na.rm = TRUE),
      .groups = "drop"
    )

  result <- scal |>
    dplyr::left_join(ia, by = c("sample_id", "draw_id")) |>
    dplyr::mutate(value = dplyr::coalesce(.data$mspaf, 0) * 100) |>
    dplyr::select("sample_id", "draw_id", "value", "n_analytes_used",
                  "n_analytes_imputed", "dominant_analyte", "max_paf") |>
    dplyr::arrange(.data$sample_id, .data$draw_id)

  ## Flat diagnostics for the kept blocks, attached as attributes (no list-cols).
  brk <- tox |>
    dplyr::select(dplyr::any_of(c("sample_id", "draw_id", "analyte",
                                  "C_adj", "PAF", "moa_group", "ref_source"))) |>
    dplyr::arrange(.data$sample_id, .data$draw_id, .data$analyte)
  dropped_flat  <- dplyr::semi_join(adj$dropped, keep_keys,
                                    by = c("sample_id", "draw_id"))
  ara_diag_flat <- dplyr::semi_join(adj$ara_diag, keep_keys,
                                    by = c("sample_id", "draw_id"))

  if (!is_draws_mode) {
    result        <- dplyr::select(result, -"draw_id")
    brk           <- dplyr::select(brk, -dplyr::any_of("draw_id"))
    dropped_flat  <- dplyr::select(dropped_flat, -dplyr::any_of("draw_id"))
    ara_diag_flat <- dplyr::select(ara_diag_flat, -dplyr::any_of("draw_id"))
  }
  attr(result, "analyte_pafs") <- brk
  attr(result, "dropped")      <- dropped_flat
  attr(result, "ara_diag")     <- ara_diag_flat
  result
}


## ============================================================================
## Shared msPAF helpers
## ============================================================================
##
## Normalisation / ARA / PAF helpers used by compute_mspaf_per_sample(), each
## vectorised across all (sample x draw) rows in one call.

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
#'   \item `"disabled"` — ARA off (no reference supplied to [add_mspaf()]);
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
.mspaf_adjust <- function(sample_rows, ref_table, ssd_params,
                           ara_enabled = TRUE) {
  ## Per-(sample, draw) keys present in this frame. The batched caller passes the
  ## whole site frame (many samples x draws); normalisation reads each row's own
  ## sample/draw co-analyte values via a wide pivot (vectorised), so this helper
  ## no longer needs one-sample-at-a-time calls. (issue #30)
  keys <- intersect(c("sample_id", "draw_id"), names(sample_rows))

  empty_dropped <- tibble::tibble(
    analyte = character(0), reason = character(0)
  )

  ## A temporal reference is keyed by (sample_id, analyte); a static one by
  ## analyte only. Join on whichever keys ref_table carries (the reference is
  ## never drawn, so a per-sample ref broadcasts to all of that sample's draws).
  ref_keys <- intersect(c("sample_id", "analyte"), names(ref_table))
  tox <- sample_rows |>
    dplyr::filter(.data$analyte %in% ssd_params$analyte) |>
    dplyr::left_join(
      dplyr::select(ssd_params, "analyte", "hc50", "sigma",
                    "moa_group", "parsed_formula", "coanalytes_req"),
      by = "analyte"
    ) |>
    dplyr::left_join(ref_table, by = ref_keys)

  if (nrow(tox) == 0L) {
    return(list(
      tox = tibble::tibble(
        analyte = character(0), value = numeric(0), detected = logical(0),
        hc50 = numeric(0), sigma = numeric(0), moa_group = character(0),
        C_norm = numeric(0), ref_norm = numeric(0),
        ref_source = character(0), C_adj = numeric(0)
      ),
      dropped = empty_dropped,
      ara_diag = tibble::tibble(
        analyte = character(0), ref_norm = numeric(0), C_norm = numeric(0),
        C_adj = numeric(0), C_excess = numeric(0), floor_fired = logical(0),
        ref_source = character(0), ref_tier = character(0)
      )
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
  ## Co-analyte values, wide per (sample, draw): one column per detected analyte.
  ## Each analyte's formula is then evaluated ONCE across all its rows
  ## (.apply_normalisation vectorises over C + co-analyte vectors), instead of a
  ## per-row pmap. Missing/undetected co-analyte → NA propagates → row dropped.
  pre_cols <- names(tox)               # columns to keep (co-analyte join adds more)
  if (length(keys) > 0L) {
    co_wide <- sample_rows |>
      dplyr::filter(.data$detected) |>
      dplyr::select(dplyr::all_of(keys), "analyte", "value") |>
      dplyr::distinct(dplyr::across(dplyr::all_of(c(keys, "analyte"))),
                      .keep_all = TRUE) |>
      tidyr::pivot_wider(names_from = "analyte", values_from = "value")
    tox <- dplyr::left_join(tox, co_wide, by = keys, suffix = c("", ".co"))
    co_col <- function(nm, idx) {
      if (nm %in% names(tox)) as.numeric(tox[[nm]][idx])
      else if (paste0(nm, ".co") %in% names(tox)) as.numeric(tox[[paste0(nm, ".co")]][idx])
      else rep(NA_real_, length(idx))
    }
  } else {
    cav <- sample_rows |>
      dplyr::filter(.data$detected) |>
      dplyr::select("analyte", "value") |> tibble::deframe()
    co_col <- function(nm, idx) {
      if (nm %in% names(cav)) rep(as.numeric(cav[[nm]]), length(idx))
      else rep(NA_real_, length(idx))
    }
  }

  tox$C_norm <- NA_real_
  for (a in unique(tox$analyte)) {
    idx <- which(tox$analyte == a)
    pf  <- tox$parsed_formula[[idx[1L]]]
    cr  <- tox$coanalytes_req[idx[1L]]
    Cv  <- tox$value[idx]
    if (is.null(pf)) {
      cn <- Cv
    } else {
      co_names <- if (nzchar(cr %||% "")) trimws(strsplit(cr, ",")[[1L]]) else character(0)
      co_names <- co_names[nzchar(co_names)]
      co_list  <- stats::setNames(lapply(co_names, co_col, idx = idx), co_names)
      cn <- .apply_normalisation(pf, Cv, co_list)
      if (length(cn) != length(idx)) cn <- rep(cn[1L], length(idx))  # error -> NA
    }
    cn[!tox$detected[idx]] <- 0   # BDL → 0 (matches the per-row `if (!q) 0`)
    tox$C_norm[idx] <- cn
  }
  ## Drop the wide co-analyte helper columns added by the join.
  tox <- tox[, intersect(c(pre_cols, "C_norm"), names(tox)), drop = FALSE]

  ## dropped / ara_diag carry the (sample_id, draw_id) keys so the batched caller
  ## can attribute them per block without per-block lookups.
  dropped <- dplyr::filter(tox, is.na(.data$C_norm)) |>
    dplyr::transmute(dplyr::across(dplyr::any_of(c("sample_id", "draw_id"))),
                     .data$analyte, reason = "missing_co_analyte")

  tox <- tox |>
    dplyr::filter(!is.na(.data$C_norm)) |>
    dplyr::mutate(C_adj = pmax(.data$C_norm - .data$ref_norm, 0))

  ## ── ARA per-cell diagnostics ─────────────────────────────────────────────
  if (!"ref_tier" %in% names(tox)) tox$ref_tier <- NA_character_
  ara_diag <- dplyr::transmute(
    tox,
    dplyr::across(dplyr::any_of(c("sample_id", "draw_id"))),
    analyte     = .data$analyte,
    ref_norm    = .data$ref_norm,
    C_norm      = .data$C_norm,
    C_adj       = .data$C_adj,
    C_excess    = .data$C_norm - .data$ref_norm,
    floor_fired = .data$C_norm < .data$ref_norm,
    ref_source  = .data$ref_source,
    ref_tier    = .data$ref_tier
  )

  list(tox = tox, dropped = dropped, ara_diag = ara_diag)
}

#' Add per-analyte PAF to adjusted tox rows, batched per analyte
#'
#' Evaluates each analyte's SSD once across **all** its `C_adj` values in a
#' single vectorised [ssdtools::ssd_hp()] call (via [.ssd_paf_vec()]), rather
#' than one [ssd_paf()] call per row.  The fitted SSD object is taken from the
#' `fit` list-column of `ssd_params` (loaded once in [derive_ssd_params()]).
#'
#' @param tox Adjusted tox rows from [.mspaf_adjust()] (needs `analyte`,
#'   `C_adj`); may span multiple samples.
#' @param ssd_params Tibble from [derive_ssd_params()] (uses `analyte`, `fit`).
#' @param method SSD method (NULL-fit fallback only).
#' @param guideline_dir Path to ANZG XLSX folder (NULL-fit fallback only).
#' @return `tox` with a numeric `PAF` column added (proportion, 0–1).
#' @keywords internal
.mspaf_add_paf <- function(tox, ssd_params, method, guideline_dir) {
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

#' Resolve or build a spline-based PAF lookup closure for one analyte
#'
#' Returns a `stats::splinefun` closure that maps `log10(conc)` to PAF.
#' Shipped tables (`NULL guideline_dir` + known method/analyte) are loaded
#' from `inst/extdata/ssd_paf_lookup.qs2` and cached in the session.
#' Runtime tables are built adaptively and likewise cached.
#'
#' @param analyte Analyte name string.
#' @param method SSD method string (e.g. `"multi"`, `"anzecc"`).
#' @param fit `fitdists` object or `NULL`.  Required for runtime builds.
#' @param guideline_dir Path to ANZG XLSX folder, or `NULL` for shipped tables.
#' @return A `splinefun` closure, or `NULL` if the lookup cannot be built.
#' @keywords internal
.ssd_paf_lookup <- function(analyte, method, fit, guideline_dir) {
  key <- if (is.null(guideline_dir)) {
    paste(method, analyte, sep = "/")
  } else {
    paste(method, analyte, guideline_dir, sep = "/")
  }

  cached <- tryCatch(get(key, envir = .ssd_paf_lookup_env, inherits = FALSE),
                     error = function(e) NULL)
  if (!is.null(cached)) return(cached)

  if (is.null(guideline_dir)) {
    shipped_path <- system.file("extdata", "ssd_paf_lookup.qs2",
                                package = "hydroSense")
    if (nzchar(shipped_path) && file.exists(shipped_path)) {
      shipped_all <- tryCatch(qs2::qs_read(shipped_path), error = function(e) NULL)
      if (!is.null(shipped_all) && key %in% names(shipped_all)) {
        entry <- shipped_all[[key]]
        lg    <- seq(entry$log10_lo, entry$log10_hi, length.out = entry$n)
        spfun <- stats::splinefun(lg, entry$paf, method = "monoH.FC")
        assign(key, spfun, envir = .ssd_paf_lookup_env)
        return(spfun)
      }
    }
    ## Fall through to runtime build if key is missing from the shipped table.
  }

  ## Runtime build path.
  if (!inherits(fit, "fitdists")) return(NULL)

  .align_hp <- function(f, cc) {
    raw <- ssdtools::ssd_hp(f, conc = cc, ci = FALSE, proportion = TRUE)
    raw <- raw[!duplicated(raw$conc), ]
    raw$est[match(cc, raw$conc)]
  }

  ## Locate effective support via forward scan.
  scan_conc <- 10^seq(log10(1e-6), log10(1e9), length.out = 4000L)
  scan_paf  <- tryCatch(.align_hp(fit, scan_conc), error = function(e) NULL)
  if (is.null(scan_paf)) return(NULL)

  lo_idx    <- max(1L, which(scan_paf >= 1e-9)[1L] - 1L)
  hi_idx    <- min(length(scan_conc), utils::tail(which(scan_paf <= 1 - 1e-9), 1L) + 1L)
  log10_lo  <- log10(scan_conc[lo_idx])
  log10_hi  <- log10(scan_conc[hi_idx])

  ## Adaptive knot doubling until max|err| < 1e-8 on 8 000 check points.
  M <- 1025L
  spfun <- tryCatch({
    repeat {
      lg  <- seq(log10_lo, log10_hi, length.out = M)
      pg  <- .align_hp(fit, 10^lg)
      f   <- stats::splinefun(lg, pg, method = "monoH.FC")

      set.seed(99L)
      q_check <- 10^stats::runif(8000L, log10_lo, log10_hi)
      tr      <- .align_hp(fit, q_check)
      max_err <- max(abs(pmin(pmax(f(log10(q_check)), 0), 1) - tr),
                     na.rm = TRUE)
      if (max_err < 1e-8) break
      if (M >= 16384L) stop("accuracy budget not met")
      M <- M * 2L
    }
    f
  }, error = function(e) NULL)

  if (is.null(spfun)) return(NULL)

  assign(key, spfun, envir = .ssd_paf_lookup_env)
  spfun
}

#' Vectorised SSD PAF lookup for one analyte
#'
#' Returns the proportion of species affected at each concentration in `conc`.
#' Concentrations that are `NA`, `<= 0`, or non-finite map to `PAF = 0`.
#' When a fitted SSD object is available it is evaluated via a spline lookup
#' (fast path) or a single [ssdtools::ssd_hp()] call (exact fallback);
#' otherwise it falls back to a per-value [ssd_paf()] lookup.
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
  pos_idx <- which(!is.na(conc) & conc > 0 & is.finite(conc))
  if (length(pos_idx) == 0L) return(out)

  if (is.null(fit)) {
    out[pos_idx] <- vapply(conc[pos_idx], function(c) {
      paf_result <- tryCatch(
        ssd_paf(analyte, c, conc_units = "ug/L", method = method,
                guideline_dir = guideline_dir, nboot = 0L),
        error = function(e) list(pct = NA_real_)
      )
      if (is.na(paf_result$pct)) NA_real_ else paf_result$pct / 100
    }, numeric(1))
    return(out)
  }

  uc <- unique(conc[pos_idx])

  ## For a non-shipped (runtime-built) table, skip the spline build when the
  ## caller's unique concentration count is below the minimum knot threshold —
  ## building a 1025-knot spline to evaluate 3 points is never a win.
  use_lookup <- is.null(guideline_dir) || length(uc) >= 1025L

  if (use_lookup) {
    spfun <- .ssd_paf_lookup(analyte, method, fit, guideline_dir)
  } else {
    spfun <- NULL
  }

  if (!is.null(spfun)) {
    paf_uc      <- pmin(pmax(spfun(log10(uc)), 0), 1)
    out[pos_idx] <- paf_uc[match(conc[pos_idx], uc)]
  } else {
    ## ssdtools::ssd_hp() expands duplicated `conc` values into a cross-join
    ## (e.g. 2 identical concentrations return 4 rows), so we cannot rely on
    ## the result being one row per input.  Evaluate on the UNIQUE
    ## concentrations and map back by value.
    res <- tryCatch(
      ssdtools::ssd_hp(fit, conc = uc, ci = FALSE, proportion = TRUE),
      error = function(e) NULL
    )
    if (is.null(res)) {
      out[pos_idx] <- NA_real_
      return(out)
    }
    res <- res[!duplicated(res$conc), c("conc", "est"), drop = FALSE]
    out[pos_idx] <- res$est[match(conc[pos_idx], res$conc)]
  }
  out
}

#' Construct an empty (zero-row) msPAF result tibble
#'
#' @param with_sample_id Logical; if `TRUE`, prepend a `sample_id` column (used
#'   by the batched [compute_mspaf_per_sample()] path).
#' @return A zero-row tibble with the msPAF output schema.
#' @keywords internal
.mspaf_empty_row <- function(with_sample_id = FALSE) {
  out <- tibble::tibble(
    value              = numeric(0),
    n_analytes_used    = integer(0),
    n_analytes_imputed = integer(0),
    dominant_analyte   = character(0),
    max_paf            = numeric(0),
    draw_id            = integer(0)
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
#' @param pafs_long Flat per-analyte PAF breakdown (the `analyte_pafs` attribute).
#' @param ara_enabled Logical; whether the caller supplied a reference.
#' @keywords internal
.summarise_ara_coverage <- function(pafs_long, ara_enabled) {
  if (!ara_enabled) return(invisible(NULL))
  if (is.null(pafs_long) || nrow(pafs_long) == 0L ||
      !"ref_source" %in% names(pafs_long))
    return(invisible(NULL))

  unmatched <- dplyr::filter(pafs_long, .data$ref_source == "unmatched")
  if (nrow(unmatched) == 0L) return(invisible(NULL))

  per_analyte <- unmatched |>
    dplyr::count(.data$analyte, name = "n_samples")

  n_analytes <- nrow(per_analyte)

  cli::cli_inform(c(
    "i" = "add_mspaf: {n_analytes} analyte{?s} had no reference match \\
           (ARA enabled) and were assessed against raw normalised \\
           concentration:",
    paste0("    ", per_analyte$analyte, ": ", per_analyte$n_samples,
           " sample", ifelse(per_analyte$n_samples == 1L, "", "s"))
  ))
  invisible(per_analyte)
}

#' Summarise per-sample dropped-analyte tally and emit a single cli message
#' @keywords internal
.summarise_mspaf_diagnostics <- function(drop_long, min_analytes = NULL) {
  if (is.null(drop_long) || nrow(drop_long) == 0L) return(invisible(NULL))

  per_analyte <- drop_long |>
    dplyr::count(.data$analyte, .data$reason, name = "n_samples")

  n_total_drops <- nrow(drop_long)
  n_samples_aff <- dplyr::n_distinct(drop_long$sample_id)

  cli::cli_inform(c(
    "i" = "add_mspaf: {n_total_drops} analyte row{?s} dropped across \\
           {n_samples_aff} sample{?s} (normalisation co-analyte missing).",
    paste0("    ", per_analyte$analyte, ": ", per_analyte$n_samples,
           " sample", ifelse(per_analyte$n_samples == 1L, "", "s"))
  ))
  invisible(per_analyte)
}


