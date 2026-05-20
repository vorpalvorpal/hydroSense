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
##   C_adj_i = max(C_i - ref_i, 0)
##
## where ref_i is the local reference concentration (80th percentile of
## matched reference site data, per ANZG convention). Evaluating the SSD
## at C_adj_i is equivalent to evaluating a rightward-shifted SSD at C_i,
## consistent with PICT theory.
##
## ## SSD derivation
##
## SSD parameters are obtained from the fitted models in paf.R via
## ssd_hc50() (for the CA HC50 denominator) and ssd_paf() (for individual
## PAF evaluation). This uses the same BCANZ-validated SSD infrastructure
## as the single-substance PAF calculations, and is correct for any SSD
## distribution shape (not only log-normal).
##
## The effective log-normal sigma for the CA mixture SSD is derived from
## HC5 and HC50 of each fitted model via .ssd_sigma() (paf.R).
##
## ## Mode of action groupings
##
## Current implementation (hybrid CA/IA):
##   Group 1 (CA): Metals + NH3 — ionoregulatory / gill disruption
##   Combined across groups via IA (architecture ready for Group 2 addition)
##
## Reserved for future implementation:
##   Group 2 (CA): Phenols + BTEX — narcosis / baseline toxicity
##   To add Group 2: implement `compute_ca_group_mspaf()` calls for phenols/BTEX
##   analytes and add to the IA combination in `compute_amspaf_per_sample()`.
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

## Analytes excluded from msPAF regardless of type in analyte_types.
##
##   NO3-N: GVs from NZ document of uncertain provenance, not ANZG
##   CH4:   Methane guideline is aesthetic/nuisance, not an ecotoxicity SSD
##   LHF:   Leachate indicator index — derived values, not a toxicity endpoint
##
.AMSPAF_EXCLUDED_ANALYTES <- c("NO3-N", "CH4", "LHF")

## analyte_types$type values that participate in msPAF and their MOA group:
##
##   "metal"    → Group 1 CA (ionoregulatory / gill disruption).
##   "nitrogen" → Group 1 CA for NH3-N only.
##   "gas"      → Group 1 CA for H2S only.
##   "organic"  → Each analyte gets its own unique IA group.
##
.AMSPAF_INCLUDED_TYPES <- c("metal", "nitrogen", "gas", "organic")
.AMSPAF_GROUP1_NAMES <- c("NH3-N", "H2S") ## non-metal group 1 members


## ============================================================================
## add_amspaf
## ============================================================================
##
## Main entry point. Computes Adjusted msPAF for each sample in df that
## has sufficient analyte coverage, and binds the results back into df.
##
## Parameters:
##   df                        long-format monitoring dataframe.
##                             Required columns: uuid.sample, uuid.feature,
##                             name.analyte, value.
##                             Optional but recommended: datetime.sample
##                             (propagated to AmsPAF rows if present),
##                             quantified (assumed TRUE if absent).
##   analyte_types             data.frame with columns name (character) and
##                             type (character). Maps analyte names to their
##                             ecological type, which controls mode-of-action
##                             group assignment. Types that participate in
##                             AmsPAF: "metal", "nitrogen" (NH3-N only),
##                             "gas" (H2S only), "organic".
##   reference_data            Optional long-format dataframe in the same
##                             structure as df (columns: name.analyte, value,
##                             quantified). Used to derive per-analyte
##                             background reference concentrations for the
##                             ARA adjustment. When NULL (default), no ARA
##                             adjustment is applied and raw concentrations
##                             are used directly.
##   method                    SSD method passed to ssd_hc50() and ssd_paf().
##                             "multi" (default): 6-distribution model average.
##                             "anzecc": per-analyte distribution matching
##                             the original ANZG derivation.
##   guideline_dir             Path to the "guideline data" folder containing
##                             ANZG XLSX files. Falls back to
##                             getOption("leachatetools.guideline_dir").
##   min_analytes              Minimum number of analytes with fitted SSDs
##                             required to compute AmsPAF for a sample.
##   ref_percentile_for_anchor Percentile used to derive the ARA reference
##                             concentration from reference_data. Default 0.80
##                             (80th percentile per ANZG convention).
##
## ============================================================================

#' Compute the Adjusted multi-substance PAF (AmsPAF) for water quality samples
#'
#' Appends AmsPAF rows to a long-format water quality dataframe. AmsPAF
#' estimates the fraction of aquatic species potentially affected by the
#' combined toxicant mixture, adjusted for local geogenic background via the
#' Added Risk Approach. See the file-level header for full methodological
#' detail.
#'
#' @param df Long-format monitoring dataframe. Required columns:
#'   \code{uuid.sample}, \code{uuid.feature}, \code{name.analyte},
#'   \code{value} (concentrations in µg/L). Optional: \code{datetime.sample}
#'   (propagated to AmsPAF rows if present), \code{quantified} (assumed
#'   \code{TRUE} if absent).
#' @param analyte_types data.frame with columns \code{name} (character) and
#'   \code{type} (character) mapping analyte names to ecological types.
#'   Types included in AmsPAF: \code{"metal"}, \code{"nitrogen"} (NH3-N
#'   only), \code{"gas"} (H2S only), \code{"organic"} (each gets its own
#'   IA group). All other types are ignored.
#' @param reference_data Optional long-format dataframe in the same structure
#'   as \code{df} providing reference site chemistry for the ARA adjustment.
#'   Required columns: \code{name.analyte}, \code{value}, \code{quantified}.
#'   When \code{NULL} (default), no ARA adjustment is applied and raw
#'   concentrations are assessed against the SSDs directly.
#' @param method SSD method passed to \code{\link{ssd_hc50}} and
#'   \code{\link{ssd_paf}}. \code{"multi"} (default) fits all 6 BCANZ
#'   distributions and model-averages; \code{"anzecc"} uses the per-analyte
#'   distribution matching the original ANZG derivation.
#' @param guideline_dir Path to the "guideline data" folder containing ANZG
#'   XLSX files. Falls back to
#'   \code{getOption("leachatetools.guideline_dir")}.
#' @param min_analytes Minimum number of analytes with fitted SSDs required
#'   to compute AmsPAF for a sample. Default \code{3}.
#' @param ref_percentile_for_anchor Percentile used to summarise
#'   \code{reference_data} concentrations per analyte for the ARA adjustment.
#'   Default \code{0.80} (80th percentile, per ANZG convention).
#'
#' @return The input \code{df} with AmsPAF rows appended. Each AmsPAF row
#'   carries \code{value} (AmsPAF as a percentage, 0–100+), \code{quantified
#'   = TRUE}, \code{name.analyte = "AmsPAF"}, \code{n_analytes_used},
#'   \code{dominant_analyte}, \code{max_paf}, \code{analyte_pafs} (list
#'   column), and four guideline columns
#'   (\code{value/level_name/guideline/comments.guideline_1} through
#'   \code{_4}). Columns present in \code{df} but not in the AmsPAF rows
#'   are \code{NA}.
#'
#' @seealso \code{\link{ssd_paf}}, \code{\link{ssd_hc50}},
#'   \code{\link{derive_ssd_params}}, \code{\link{compute_amspaf_per_sample}}
#'
#' @references
#' De Zwart D, Posthuma L (2005) Environmental Toxicology and Chemistry
#' 24(10):2665-2676.
#'
#' @export

add_amspaf <- function(
  df,
  analyte_types,
  reference_data = NULL,
  method = c("multi", "anzecc"),
  guideline_dir = getOption("leachatetools.guideline_dir"),
  min_analytes = 3,
  ref_percentile_for_anchor = 0.80
) {
  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df), must.include = c("uuid.sample", "uuid.feature", "name.analyte", "value"))
  checkmate::assert_data_frame(analyte_types)
  checkmate::assert_names(names(analyte_types), must.include = c("name", "type"))
  checkmate::assert_data_frame(reference_data, null.ok = TRUE)
  method <- match.arg(method)
  checkmate::assert_int(min_analytes, lower = 1)
  checkmate::assert_number(ref_percentile_for_anchor, lower = 0.5, upper = 0.99)

  ## ================================================================
  ## Step 1: Derive SSD parameters using ssd_hc50() and .ssd_sigma().
  ##
  ## Done once per function call (O(n_analytes)), not per sample.
  ## ================================================================

  ssd_params <- derive_ssd_params(analyte_types, method, guideline_dir)

  if (nrow(ssd_params) == 0) {
    warning(
      "No analytes with fitted SSDs found in analyte_types. ",
      "AmsPAF cannot be computed."
    )
    return(df)
  }

  ## ================================================================
  ## Step 2: Pre-compute reference concentrations for ARA adjustment.
  ##
  ## If reference_data is NULL, ref_local = 0 for all analytes
  ## (no ARA adjustment; raw concentrations assessed directly).
  ## ================================================================

  ref_concs <- get_reference_concentrations(
    reference_data, ref_percentile_for_anchor, ssd_params$name.analyte
  )

  ## ================================================================
  ## Step 3: Compute AmsPAF per feature, per sample.
  ## ================================================================

  amspaf_df <-
    df |>
    dplyr::filter(name.analyte %in% ssd_params$name.analyte) |>
    dplyr::group_by(uuid.feature) |>
    dplyr::group_modify(\(.x, .y) {
      compute_amspaf_per_sample(
        sample_data = .x,
        reference_data = ref_concs,
        ssd_params = ssd_params,
        min_analytes = min_analytes,
        method = method,
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
      ## Guideline columns
      ##
      ## All four tier breaks are biologically anchored to ANZG species
      ## protection levels. Values are in percentage (e.g. 1.0 = 1%).
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
    amspaf_df <- dplyr::left_join(amspaf_df, sample_times, by = "uuid.sample")
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
##
## Builds the SSD parameter table for eligible analytes by calling
## ssd_hc50() and .ssd_sigma() from paf.R. No GV data or analyte metadata
## CSV is needed — all SSD information comes from the fitted models.
##
## Eligibility:
##   - analyte type is in .AMSPAF_INCLUDED_TYPES
##   - for "nitrogen" type: only NH3-N (other N species have no toxicity SSD)
##   - for "gas" type: only H2S
##   - name not in .AMSPAF_EXCLUDED_ANALYTES
##   - a fitted SSD exists (ssd_hc50() returns a non-NA value)
##
## Returns a tibble with columns:
##   name.analyte, hc50 (µg/L), sigma (effective log-normal sigma), moa_group
##
## moa_group assignment:
##   1 — ionoregulatory CA group (metals, NH3-N, H2S)
##   2, 3, ... — each organic gets a unique integer (never CA-combined)
##
## ============================================================================

#' Derive SSD parameters for AmsPAF computation
#'
#' Internal function called by \code{\link{add_amspaf}} to build the analyte
#' parameter table from fitted SSD models. For each eligible analyte,
#' \code{ssd_hc50()} provides the HC50 used as the CA denominator, and
#' \code{.ssd_sigma()} derives the effective log-normal sigma needed for the
#' concentration-addition mixture SSD.
#'
#' @param analyte_types data.frame with columns \code{name} and \code{type}.
#' @param method SSD method: \code{"multi"} or \code{"anzecc"}.
#' @param guideline_dir Path to ANZG guideline data folder.
#'
#' @return Tibble with columns \code{name.analyte}, \code{hc50},
#'   \code{sigma}, \code{moa_group}.
#'
#' @keywords internal

derive_ssd_params <- function(analyte_types, method, guideline_dir) {
  eligible <-
    analyte_types |>
    dplyr::filter(type %in% .AMSPAF_INCLUDED_TYPES) |>
    dplyr::filter(
      !(type == "nitrogen" & !name %in% .AMSPAF_GROUP1_NAMES),
      !(type == "gas" & !name %in% .AMSPAF_GROUP1_NAMES)
    ) |>
    dplyr::filter(!name %in% .AMSPAF_EXCLUDED_ANALYTES)

  if (nrow(eligible) == 0) {
    return(tibble::tibble(
      name.analyte = character(),
      hc50 = numeric(),
      sigma = numeric(),
      moa_group = integer()
    ))
  }

  params <- purrr::map_dfr(seq_len(nrow(eligible)), function(i) {
    analyte_name <- eligible$name[i]
    analyte_type <- eligible$type[i]

    hc50 <- ssd_hc50(analyte_name, method = method, guideline_dir = guideline_dir)
    if (is.na(hc50)) return(NULL)

    sigma <- .ssd_sigma(analyte_name, method = method, guideline_dir = guideline_dir)

    ## moa_group: 1 for ionoregulatory CA group; NA sentinel for organics
    ## (unique IDs assigned after the loop)
    moa_group <- if (analyte_type == "organic") NA_integer_ else 1L

    tibble::tibble(
      name.analyte = analyte_name,
      hc50 = hc50,
      sigma = sigma,
      moa_group = moa_group
    )
  })

  if (nrow(params) == 0) {
    return(tibble::tibble(
      name.analyte = character(),
      hc50 = numeric(),
      sigma = numeric(),
      moa_group = integer()
    ))
  }

  ## Assign unique group IDs to organics (sentinelled NA → 2, 3, ...).
  ## Group 1 is reserved for the ionoregulatory CA group.
  organic_idx <- which(is.na(params$moa_group))
  params$moa_group[organic_idx] <- seq_along(organic_idx) + 1L

  params
}


## ============================================================================
## get_reference_concentrations
## ============================================================================
##
## Computes per-analyte ARA reference concentrations from reference_data.
## Returns a tibble(name.analyte, ref_local) with ref_local = 0 when
## reference_data is NULL (no ARA adjustment).
##
## BDL handling: below-detection-limit measurements at reference sites are
## treated as zero (no natural background). The reasoning: if a metal is
## undetectable at reference, the natural background is genuinely near zero,
## and imputing a positive value would artificially inflate the ARA shift.
##
## ============================================================================

get_reference_concentrations <- function(reference_data, anchor_p, valid_analytes) {
  if (is.null(reference_data)) {
    return(tibble::tibble(
      name.analyte = valid_analytes,
      ref_local = 0
    ))
  }

  checkmate::assert_names(
    names(reference_data),
    must.include = c("name.analyte", "value", "quantified")
  )

  reference_data |>
    dplyr::filter(name.analyte %in% valid_analytes) |>
    dplyr::group_by(name.analyte) |>
    dplyr::summarise(
      ref_local = quantile(
        dplyr::if_else(quantified, value, 0),
        probs = anchor_p,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}


## ============================================================================
## compute_ca_group_mspaf
## ============================================================================
##
## Computes msPAF for a single CA group (set of analytes acting by the same
## mechanism) using Concentration Addition.
##
## Method (De Zwart & Posthuma 2005, eq. 5-9):
##   1. TU_i = C_adj_i / HC50_i  (toxic unit relative to SSD median)
##   2. TU_mix = sum(TU_i)
##   3. Evaluate mixture SSD at TU_mix using concentration-weighted sigma
##
## Mixture SSD under CA:
##   Evaluated as a log-normal with mu_mix = 0 (by construction) and
##   sigma_mix = sqrt(sum(w_i^2 * sigma_i^2)) where w_i = TU_i / TU_mix.
##   sigma_i is the effective log-normal sigma derived from HC5/HC50 of
##   each fitted SSD via .ssd_sigma().
##
## Parameters:
##   group_data  tibble with columns C_adj, hc50, sigma (effective log-normal
##               sigma from .ssd_sigma()), moa_group; one row per analyte
##
## Returns: msPAF_CA (proportion, not percentage)
##
## ============================================================================

compute_ca_group_mspaf <- function(group_data) {
  ## Drop analytes with zero ARA-adjusted concentration or missing parameters
  group_data <- dplyr::filter(
    group_data,
    C_adj > 0,
    !is.na(hc50), hc50 > 0,
    !is.na(sigma), sigma > 0
  )

  if (nrow(group_data) == 0) {
    return(0)
  }

  group_data <- dplyr::mutate(group_data, TU = C_adj / hc50)
  TU_mix <- sum(group_data$TU)

  if (TU_mix <= 0) {
    return(0)
  }

  group_data <- dplyr::mutate(group_data, w = TU / TU_mix)
  sigma_mix <- sqrt(sum(group_data$w^2 * group_data$sigma^2))

  mspaf_ca <- pnorm(log10(TU_mix) / sigma_mix)
  return(mspaf_ca)
}


## ============================================================================
## compute_amspaf_per_sample
## ============================================================================
##
## Computes AmsPAF for each sample in sample_data using pre-built reference
## concentrations and SSD parameters.
##
## Steps:
##   1. ARA shift: C_adj_i = max(C_i - ref_i, 0)
##   2. Individual PAF per analyte via ssd_paf() (for diagnostics)
##   3. Split analytes into mode-of-action groups
##   4. Compute msPAF_CA per group via compute_ca_group_mspaf()
##   5. Combine groups via Independent Action:
##        AmsPAF = 1 - prod(1 - msPAF_CA_k) for k in groups
##   6. Convert to percentage
##
## Returns one row per sample with: uuid.sample, value (AmsPAF %),
## and diagnostic columns.
##
## ============================================================================

compute_amspaf_per_sample <- function(
  sample_data,
  reference_data,
  ssd_params,
  min_analytes,
  method,
  guideline_dir
) {
  sample_data |>
    dplyr::group_by(uuid.sample) |>
    dplyr::filter(dplyr::n_distinct(name.analyte) >= min_analytes) |>
    dplyr::left_join(reference_data, by = "name.analyte") |>
    dplyr::left_join(ssd_params, by = "name.analyte") |>
    dplyr::filter(!is.na(hc50), !is.na(sigma)) |>
    dplyr::mutate(
      ## ARA shift. If ref_local is NA (analyte absent from reference_data),
      ## treat background as zero — conservative but avoids silent data loss.
      ref_local = tidyr::replace_na(ref_local, 0),
      C_adj = pmax(value - ref_local, 0)
    ) |>
    dplyr::group_modify(\(.x, .y) {
      if (nrow(.x) == 0) {
        return(tibble::tibble(
          value = numeric(0),
          n_analytes_used = integer(0),
          dominant_analyte = character(0),
          max_paf = numeric(0),
          analyte_pafs = list()
        ))
      }

      ## Per-analyte PAF via ssd_paf() for diagnostics and dominant-analyte
      ## identification. ssd_paf() uses the fitted model, so results are
      ## correct for any distribution shape, not only log-normal.
      .x <- dplyr::mutate(.x,
        PAF = purrr::map2_dbl(C_adj, name.analyte, function(c, a) {
          if (is.na(c)) return(NA_real_)
          if (c <= 0) return(0)
          paf_result <- ssd_paf(
            a, c,
            method = method,
            guideline_dir = guideline_dir,
            nboot = 0L
          )
          if (is.na(paf_result$pct)) NA_real_ else paf_result$pct / 100
        })
      )

      ## CA msPAF per mode-of-action group
      groups <- unique(.x$moa_group)
      mspaf_by_group <- vapply(
        groups,
        function(g) compute_ca_group_mspaf(dplyr::filter(.x, moa_group == g)),
        numeric(1)
      )

      ## Combine groups via Independent Action
      amspaf <- 1 - prod(1 - mspaf_by_group)

      dominant <- if (any(!is.na(.x$PAF))) {
        .x$name.analyte[which.max(.x$PAF)]
      } else {
        NA_character_
      }

      tibble::tibble(
        value = amspaf * 100,
        n_analytes_used = nrow(.x),
        dominant_analyte = dominant,
        max_paf = if (nrow(.x) > 0) max(.x$PAF, na.rm = TRUE) else NA_real_,
        analyte_pafs = list(dplyr::select(.x, name.analyte, C_adj, PAF, moa_group))
      )
    }) |>
    dplyr::ungroup()
}


## ============================================================================
## Utility: classify AmsPAF value into tier label
## ============================================================================

classify_amspaf_tier <- function(amspaf_pct) {
  checkmate::assert_numeric(amspaf_pct, lower = 0)
  dplyr::case_when(
    amspaf_pct < .AMSPAF_BREAK_T1_T2 * 100 ~ "1_background",
    amspaf_pct < .AMSPAF_BREAK_T2_T3 * 100 ~ "2_elevated",
    amspaf_pct < .AMSPAF_BREAK_T3_T4 * 100 ~ "3_impacted",
    TRUE ~ "4_severely_impacted"
  )
}
