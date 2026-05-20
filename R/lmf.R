## ============================================================================
## Leachate Mixing Fraction (LMF)
## ============================================================================
##
## PURPOSE
## -------
## The LMF answers the question: "Is leachate currently entering this sample
## point, and to what degree?"
##
## It is a SOURCE DETECTION index, not a toxicity index. It uses only
## conservative and quasi-conservative major-ion tracers. Non-conservative
## toxicants (metals, trace organics) are deliberately excluded because they
## have sediment storage and kinetic lag effects that would conflate "current
## leachate intrusion" with "historical plume residual".
##
## CONCEPTUAL BASIS
## ----------------
## The LMF is based on end-member mixing analysis (EMMA; Christophersen and
## Hooper 1992, Water Resources Research 28(1)). Two end-members are defined:
##
##   R (reference): clean background water at a matched reference site
##   L (leachate):  landfill leachate, characterised from leachate sampling
##                  infrastructure (ponds, sumps, treatment plant influent)
##
## For each measured ion, a per-ion mixing fraction is estimated:
##
##   f_i = (x_i - R_i) / (L_i - R_i)
##
## where x_i is the sample concentration, R_i is the reference mean, and
## L_i is the leachate end-member concentration.
##
## For a clean two-component mixture, every f_i should equal the true
## leachate mixing fraction. The LMF is the inverse-variance weighted mean
## of these per-ion estimates, giving more weight to ions that are more
## informative (large gradient L_i - R_i, low variability in reference).
##
## INVERSE-VARIANCE WEIGHTING
## ---------------------------
## Each ion i contributes a mixing-fraction estimate with variance:
##
##   sigma2_f_i ~= (sigma2_meas_i + sigma2_R_i) / (L_i - R_i)^2
##
## Precision is determined by:
##   - How well separated R and L are for this ion (large L - R = high
##     precision)
##   - How variable the reference is for this ion (large sigma_R = low
##     precision)
##   - Analytical measurement noise
##
## NH4 (near-zero in reference, high in leachate) therefore gets very high
## weight. Ca (variable in both) gets low weight. This is appropriate.
##
## LMF = sum(w_i * f_i) / sum(w_i)   where w_i = 1 / sigma2_f_i
##
## ION INFORMATIVENESS AND ADMISSION GATE
## ---------------------------------------
## Before computing LMF for any sample, each ion's intrinsic informativeness
## is assessed from the end-member calibration data alone:
##
##   informativeness_i = sigma_R_i / |L_i - R_i|
##
## Small = highly informative (large gradient, stable reference).
## Large = poorly informative (small gradient or noisy reference).
##
## Ions with informativeness_i <= informativeness_threshold are classified
## as "high-information". This classification is done once at calibration
## time and reflects stable chemical properties of the end-members, not
## per-sample missingness.
##
## A sample must have at least min_high_info_ions of these measured and
## available to receive an LMF value. Low-information ions are not excluded
## from the calculation — they contribute via their (small) inverse-variance
## weights.
##
## LEACHATE END-MEMBER STABILITY
## ------------------------------
## Raw leachate EC varies substantially (e.g., 5,000-25,000 uS/cm) between
## samples. However, the RATIOS of ions to Cl- are stable, reflecting
## leachate composition rather than dilution state.
##
## The leachate end-member L is defined via Cl-anchoring:
##   1. Compute mean ratio of each ion to Cl across leachate samples
##   2. Fix absolute scale at median Cl across leachate samples
##   3. L_i = mean_ratio_i * median_Cl
##
## Cl is used as the anchor because it is the most conservative groundwater
## tracer (no sorption, no speciation, no biological uptake).
##
## SPECIES COLLAPSING
## ------------------
## Ions that transform between species along flow paths while conserving
## total mass are collapsed before computing mixing fractions:
##
##   Total N = NH4-N + NO3-N + NO2-N
##     NH4 is nitrified to NO3 in oxic groundwater; collapsing to total N
##     makes the index robust across the redox gradient.
##     NOTE: Denitrification (N2 loss) breaks true N conservation,
##     introducing a small negative bias. No correction applied.
##
##   Total alkalinity = CO3 + HCO3 (as CaCO3 in meq/L)
##     CO3/HCO3 ratio is pH-dependent; total alkalinity is conserved
##     under most groundwater conditions.
##
## QUALITY CONTROL OUTPUTS
## -----------------------
## sigma_LMF: uncertainty on the mixing fraction estimate.
##   sigma_LMF = 1 / sqrt(sum(w_i))
##   High sigma_LMF = few informative ions measured. Returns NA if
##   sigma_LMF > max_sigma_lsi (expressed in percentage points).
##
## chi2_per_df: ion agreement diagnostic.
##   chi2 = sum(w_i * (f_i - LMF)^2),  df = n_ions - 1
##   chi2/df ~ 1 under the linear mixing model. High values indicate ion
##   disagreement: transformation activity (SO4 reduction, denitrification)
##   or alternative water source. Returns NA if chi2/df > max_chi2_per_df.
##
## NOTE: Full iterative per-species downweighting based on off-axis
## residuals is a planned future extension.
##
## TIER SYSTEM
## -----------
## LMF is reported as a percentage (0 = pure reference, 100 = pure leachate).
##
##   Tier 1 (Background):          LMF <= 1    (<=1% leachate equivalent)
##   Tier 2 (Trace impact):        LMF <= 5    (1-5%)
##   Tier 3 (Significant impact):  LMF <= 20   (5-20%)
##   Tier 4 (Severe impact):       LMF >  20   (>20%)
##
## Review these breaks against reference-site LMF distribution after
## deployment. The tier 1/2 break should sit near the 95th percentile
## of LMF on reference samples.
##
## CALIBRATION WINDOW
## ------------------
## Mirrors the MTUI calibration window logic (see add_metal_index.R):
##   - Derived from the input df's date range, not Sys.Date()
##   - Centred on the input data's date range where possible
##   - Shifted backwards if future data are unavailable
##   - Minimum sample counts applied separately to reference and leachate
##
## REFERENCES
## ----------
## Christophersen N, Hooper RP (1992) Water Resources Research 28(1):99-107.
## Christensen JB et al. (2001) Applied Geochemistry 16(7-8):659-718.
## Kjeldsen P et al. (2002) Crit Rev Environ Sci Technol 32(4).
## ============================================================================

## ============================================================================
## File-level constants
## ============================================================================

## Tier breaks in percentage units (0-100).
## Review and adjust after deployment once you have reference-site LMF data.
.LMF_BREAK_TRACE <- 1 ## 1% leachate equivalent
.LMF_BREAK_MODERATE <- 5 ## 5% leachate equivalent
.LMF_BREAK_SEVERE <- 20 ## 20% leachate equivalent
.LMF_BREAK_CEILING <- 100 ## pure leachate end-member

## Candidate analyte names in the meq-converted panel (trailing "_").
.LMF_ANALYTES_MEQ <- c(
  "Na_",
  "K_",
  "Ca_",
  "Mg_",
  "Cl_",
  "SO4\u00b2\u207b_",
  "F_",
  "NH3-N_",
  "NO3-N_",
  "NO2-N_",
  "CO3-CaCO3_",
  "HCO3-CaCO3_"
)

## Panel after species collapsing.
.LMF_PANEL_COLLAPSED <- c(
  "Cl_",
  "Na_",
  "K_",
  "Ca_",
  "Mg_",
  "total_N_",
  "total_alk_",
  "SO4\u00b2\u207b_",
  "F_"
)


## ============================================================================
## add_lmf() — main entry point
## ============================================================================

#' Compute the Leachate Mixing Fraction (LMF) for water quality samples
#'
#' Appends LMF rows to a long-format water quality dataframe. LMF is a
#' source-detection index estimating what fraction of a sample's chemistry
#' can be attributed to the leachate end-member under a two-component mixing
#' model. See the file-level header for full methodological detail.
#'
#' @param df Long-format dataframe (same structure as \code{data_df()} output).
#'   Required columns: \code{uuid.sample}, \code{uuid.feature},
#'   \code{name.analyte}, \code{value}, \code{quantified},
#'   \code{datetime.sample}.
#' @param calibration_window_years Number of years of calibration data to use.
#'   Window is centred on the input data's date range and shifted backwards if
#'   future data are unavailable. Default \code{5}.
#' @param min_leachate_total_n_mgl Minimum total N (mg/L, sum of NH4 + NO3 +
#'   NO2) for a leachate sample to qualify for end-member calibration. Excludes
#'   non-representative samples (dilution events, treatment upsets). Total N is
#'   used rather than NH4-N alone because aerated leachate features may have
#'   nitrified NH4 to NO3. Default \code{20}.
#' @param rsd_default Default relative analytical uncertainty, applied as
#'   \code{sigma_meas = rsd_default * |x|} with a floor at
#'   \code{rsd_default * |R|}. Default \code{0.05} (5\% RSD).
#' @param min_ref_samples Minimum reference samples required in the calibration
#'   window. Features matched to a reference site with fewer samples are
#'   skipped silently. Default \code{10}.
#' @param min_leachate_samples Minimum valid leachate samples required for
#'   end-member construction. The function stops with an informative error if
#'   not met. Default \code{10}.
#' @param max_sigma_lsi Maximum permitted \code{sigma_lmf} in percentage points
#'   (uncertainty on the mixing fraction estimate). Samples above this threshold
#'   return \code{NA} with reason code \code{"insufficient_precision"}.
#'   Default \code{10} (i.e., ±10 percentage points).
#' @param max_chi2_per_df Maximum permitted chi-squared per degree of freedom
#'   (computed on the \emph{original} inverse-variance weights, before robust
#'   reweighting). Retained as a diagnostic output in all cases. Default
#'   \code{Inf} (no hard gate — the robust reweighting handles outlier ions
#'   directly, making a hard chi2 gate redundant for most samples).
#' @param informativeness_threshold Threshold on \code{sigma_R / |L - R|} for
#'   classifying ions as high-information. Computed once from calibration data,
#'   not per-sample. Lower values are more stringent. Default \code{0.20}.
#' @param min_high_info_ions Minimum number of high-information ions that must
#'   be measured in a sample for LMF to be computed. Samples below this
#'   threshold return \code{NA} with reason code
#'   \code{"insufficient_high_info_ions"}. Default \code{3L}.
#' @param robust_iterations Number of Huber M-estimator reweighting passes.
#'   In each pass, ions whose residuals exceed \code{robust_threshold_k *
#'   MAD(residuals)} are downweighted proportionally. Three iterations is
#'   sufficient for convergence in practice. Set to \code{0L} to disable
#'   robust reweighting entirely. Default \code{3L}.
#' @param robust_threshold_k Threshold multiplier applied to the unweighted
#'   median absolute deviation (MAD) of per-ion residuals. Ions with
#'   \code{|r_i| > k * MAD} have their estimation weight reduced by a factor
#'   of \code{k * MAD / |r_i|}. Smaller values downweight more aggressively;
#'   \code{k = 1.5} is fairly tight, \code{k = 2.5} is close to the
#'   conventional Huber constant. Default \code{1.5}.
#' @param verbose If \code{TRUE}, prints an ion informativeness table via
#'   \code{cli::cli_inform()} showing each ion's score and high/low
#'   classification. Useful when tuning \code{informativeness_threshold}.
#'   Default \code{FALSE}.
#'
#' @param leachate_data Optional long-format dataframe in the same structure
#'   as \code{df} providing leachate end-member chemistry. Required columns:
#'   \code{uuid.sample}, \code{name.analyte}, \code{value}, \code{quantified}.
#'   When supplied, \emph{all} samples are used to build the leachate
#'   end-member regardless of date and regardless of total-N content (the
#'   total-N quality filter is bypassed — supply curated data). Samples
#'   missing the required LMF panel analytes are silently dropped.
#'   \code{NULL} (default) uses the standard leachate feature detection
#'   logic (requires dashboard infrastructure: \code{feature_df()},
#'   \code{data_df()}).
#' @param reference_data Optional long-format dataframe in the same structure
#'   as \code{df} providing reference site chemistry for end-member
#'   calibration. Required columns: \code{uuid.sample}, \code{name.analyte},
#'   \code{value}, \code{quantified}. When supplied, \emph{all} samples are
#'   used to build the reference end-member regardless of date; the
#'   calibration window is ignored. The same end-member is applied globally
#'   to all features. \code{NULL} (default) uses the standard per-feature
#'   matched reference site logic (requires dashboard infrastructure:
#'   \code{feature_sfc()}, \code{data_df()}).
#' @param lmf_analyte_uuid UUID assigned to LMF rows in the \code{uuid.analyte}
#'   column of the output. Set to a fixed UUID if integrating with a data
#'   system that tracks analytes by UUID. Default \code{NA_character_}.
#'
#' @return The input \code{df} with LMF rows appended. Each LMF row carries
#'   \code{value} (the robust LMF estimate as a percentage, 0 = pure reference,
#'   100 = pure leachate), \code{name.analyte = "LMF"}, \code{uuid.feature}
#'   (from the input), \code{lmf_naive} (the non-robust estimate for
#'   comparison), \code{lmf_reason} (\code{NA} on success, reason code on
#'   failure), \code{n_ions_used}, \code{n_ions_downweighted} (count of ions
#'   whose weight was reduced by robust reweighting), \code{sigma_lmf},
#'   \code{chi2_per_df} (diagnostic; computed on original weights), and four
#'   guideline columns (\code{value/level_name/guideline/comments.guideline_1}
#'   through \code{_4}). Columns present in \code{df} but not produced by
#'   the LMF computation are \code{NA} in the appended rows.
#'
#' @seealso \code{\link{build_leachate_endmember}},
#'   \code{\link{build_reference_endmember}},
#'   \code{\link{compute_lmf_for_sample}}
#'
#' @references
#' Christophersen N, Hooper RP (1992) Water Resources Research 28(1):99-107.
#' Christensen JB et al. (2001) Applied Geochemistry 16(7-8):659-718.
#' Kjeldsen P et al. (2002) Crit Rev Environ Sci Technol 32(4).
#'
#' @export

add_lmf <- function(
  df,
  leachate_data = NULL,
  reference_data = NULL,
  lmf_analyte_uuid = NA_character_,
  calibration_window_years = 5,
  min_leachate_total_n_mgl = 20,
  rsd_default = 0.05,
  min_ref_samples = 10,
  min_leachate_samples = 10,
  max_sigma_lsi = 10, ## percentage points
  max_chi2_per_df = Inf, ## diagnostic only; no hard gate
  informativeness_threshold = 0.20,
  min_high_info_ions = 3L,
  robust_iterations = 3L,
  robust_threshold_k = 1.5,
  verbose = FALSE
) {
  checkmate::assert_data_frame(df)
  checkmate::assert_data_frame(leachate_data, null.ok = TRUE)
  checkmate::assert_data_frame(reference_data, null.ok = TRUE)
  checkmate::assert_character(lmf_analyte_uuid, len = 1L)
  checkmate::assert_number(calibration_window_years, lower = 1)
  checkmate::assert_number(min_leachate_total_n_mgl, lower = 0)
  checkmate::assert_number(rsd_default, lower = 0, upper = 1)
  checkmate::assert_int(min_ref_samples, lower = 1)
  checkmate::assert_int(min_leachate_samples, lower = 1)
  checkmate::assert_number(max_sigma_lsi, lower = 0)
  checkmate::assert_number(max_chi2_per_df, lower = 1)
  checkmate::assert_number(informativeness_threshold, lower = 0)
  checkmate::assert_int(min_high_info_ions, lower = 1)
  checkmate::assert_int(robust_iterations, lower = 0)
  checkmate::assert_number(robust_threshold_k, lower = 0.1)
  checkmate::assert_flag(verbose)

  uuid_lmf <- lmf_analyte_uuid

  ## ================================================================
  ## Step 1: Derive calibration window from input df's date range.
  ##
  ## NOT derived from Sys.Date(). Anchored to the input data so
  ## results are reproducible and appropriate for the assessed period.
  ## The window is centred on the data's midpoint and shifted backwards
  ## per-site below if future calibration data are unavailable.
  ## ================================================================

  df_dates <- as.Date(df$datetime.sample)
  df_start <- min(df_dates, na.rm = TRUE)
  df_end <- max(df_dates, na.rm = TRUE)
  df_span <- as.numeric(difftime(df_end, df_start, units = "days")) / 365.25

  if (df_span >= calibration_window_years) {
    calibration_start <- df_start
    calibration_end <- df_end
  } else {
    df_centre <- df_start + (df_end - df_start) / 2
    half_window_days <- (calibration_window_years * 365.25) / 2
    calibration_start <- df_centre - lubridate::days(round(half_window_days))
    calibration_end <- df_centre + lubridate::days(round(half_window_days))
  }

  ## ================================================================
  ## Step 2: Prepare input data.
  ##
  ## Convert to meq/L (once only), extract LMF panel, verify no
  ## duplicates, apply BDL half-detection-limit replacement, collapse
  ## transforming species to conserved totals, pivot to wide format.
  ## ================================================================

  df_meq <-
    df |>
    to_meq() |>
    filter(name.analyte %in% .LMF_ANALYTES_MEQ)

  ## Verify no duplicate analyte values per sample.
  df_meq |>
    group_by(uuid.sample, name.analyte) |>
    mutate(n_values = n()) |>
    assertr::verify(
      n_values == 1,
      description = paste0(
        "ERROR: Multiple values of the same analyte in a single sample. ",
        "Resolve with the `multiple_values` argument of `data_df()` ",
        "before calling `add_lmf()`."
      )
    ) |>
    ungroup()

  ## BDL: replace detection-limit values with half the detection limit.
  ## Applied before collapsing so totals are correct.
  df_meq <-
    df_meq |>
    mutate(value = if_else(!quantified, value * 0.5, value))

  ## Pivot and collapse. The working wide panel has:
  ##   Cl_, Na_, K_, Ca_, Mg_   (conservative major ions)
  ##   total_N_                  (NH4 + NO3 + NO2 as N in meq)
  ##   total_alk_                (CO3 + HCO3 as CaCO3 in meq)
  ##   SO42-_                    (quasi-conservative, flagged by chi2)
  ##   F_                        (optional)
  df_wide <-
    df_meq |>
    select(uuid.sample, uuid.feature, datetime.sample, name.analyte, value) |>
    collapse_species(
      id_cols = c("uuid.sample", "uuid.feature", "datetime.sample")
    )

  ## ================================================================
  ## Step 3: Build the leachate end-member (once, shared all features).
  ##
  ## When leachate_data is provided, bypass the standard
  ## leachate feature detection and calibration window logic entirely.
  ## All samples in the override are used; the total-N quality filter
  ## is not applied (caller is responsible for data quality).
  ## ================================================================

  if (!is.null(leachate_data)) {
    if (verbose) {
      cli::cli_inform(
        "i" = paste0(
          "Using leachate_data (",
          n_distinct(leachate_data$uuid.sample),
          " samples, ",
          n_distinct(leachate_data$uuid.feature),
          " feature(s))."
        )
      )
    }
    leachate_em <- build_endmember_from_override(
      override_df = leachate_data,
      type = "leachate"
    )
  } else {
    leachate_em <- build_leachate_endmember(
      calibration_start = calibration_start,
      calibration_end = calibration_end,
      min_leachate_total_n_mgl = min_leachate_total_n_mgl,
      min_leachate_samples = min_leachate_samples
    )
  }

  ## ================================================================
  ## Step 4: Compute per-ion informativeness and classify ions.
  ##
  ## informativeness_i = sigma_R_i / |L_i - R_i|
  ##
  ## This is a system-wide property, computed from a pooled reference
  ## end-member, not per-feature. NH4 will reliably score as high-
  ## information at any site; Ca will reliably score as low.
  ##
  ## The high-information classification gates sample admission (must
  ## have >= min_high_info_ions measured). All ions still contribute
  ## to the weighted mean via their inverse-variance weights.
  ## ================================================================

  if (!is.null(reference_data)) {
    if (verbose) {
      cli::cli_inform(
        "i" = paste0(
          "Using reference_data for informativeness calibration (",
          n_distinct(reference_data$uuid.sample),
          " samples, ",
          n_distinct(reference_data$uuid.feature),
          " feature(s))."
        )
      )
    }
    pooled_ref <- build_endmember_from_override(
      override_df = reference_data,
      type = "reference"
    )$stats
  } else {
    pooled_ref <- build_pooled_reference_endmember(
      calibration_start = calibration_start,
      calibration_end = calibration_end
    )
  }

  informativeness_tbl <-
    pooled_ref |>
    inner_join(leachate_em$L_values |> select(ion, L), by = "ion") |>
    mutate(
      gradient = abs(L - R),
      informativeness = if_else(
        gradient > .Machine$double.eps * 1000,
        sigma_R / gradient,
        Inf
      ),
      high_info = informativeness <= informativeness_threshold
    ) |>
    arrange(informativeness)

  high_info_ions <- informativeness_tbl |>
    filter(high_info) |>
    pull(ion)

  ## ---------------------------------------------------------------
  ## Optionally report the informativeness table.
  ## ---------------------------------------------------------------
  if (verbose) {
    rows <- informativeness_tbl |>
      mutate(
        label = sprintf(
          "  %-18s  %.4f  [%s]",
          ion,
          informativeness,
          if_else(high_info, "HIGH", "low")
        )
      ) |>
      pull(label)

    cli::cli_inform(c(
      "i" = "Ion informativeness (threshold = {informativeness_threshold}):",
      rows,
      "i" = "High-information ions ({length(high_info_ions)}): {paste(high_info_ions, collapse=', ')}",
      "i" = "Minimum high-information ions required per sample: {min_high_info_ions}"
    ))
  }

  ## Stop early if the full panel can never meet the minimum — indicates
  ## a misconfigured threshold or broken end-members.
  if (length(high_info_ions) < min_high_info_ions) {
    stop(glue::glue(
      "Only {length(high_info_ions)} high-information ion(s) found ",
      "(informativeness <= {informativeness_threshold}), but ",
      "min_high_info_ions = {min_high_info_ions}. ",
      "Lower informativeness_threshold or min_high_info_ions, or check ",
      "end-member data. Run with verbose = TRUE to inspect."
    ))
  }

  ## ================================================================
  ## Step 5: Compute per-sample LMF for each feature in the input df.
  ##
  ## When reference_data is provided, the end-member is built once
  ## here (outside the per-feature loop) and applied to all features.
  ## When NULL, per-feature reference matching uses dashboard globals.
  ## ================================================================

  ## Pre-build the reference end-member if reference_data was provided.
  ## NULL when using standard per-feature matching.
  ref_endmember_override <- if (!is.null(reference_data)) {
    build_endmember_from_override(
      override_df = reference_data,
      type = "reference"
    )
  } else {
    NULL
  }

  lsi_results <-
    df_wide |>
    group_by(uuid.feature) |>
    group_modify(\(.x, .y) {
      feature_uuid <- .y[[1, 1]]
      if (!is.null(ref_endmember_override)) {
        ## Use the pre-built override end-member (same for all features).
        ref_endmember <- ref_endmember_override
      } else {
        ref_feature <- get_reference_site(feature_uuid)
        ref_endmember <- build_reference_endmember(
          reference_feature_uuid = ref_feature$uuid,
          cal_start = calibration_start,
          cal_end = calibration_end
        )
      }

      if (is.null(ref_endmember)) {
        message(glue::glue(
          "No reference data for feature {feature_uuid}. Skipping."
        ))
        return(tibble())
      }

      ## Apply minimum sample count only when NOT using an override.
      ## With an override the caller controls data quality directly.
      if (
        is.null(ref_endmember_override) &&
          ref_endmember$n_samples < min_ref_samples
      ) {
        message(glue::glue(
          "Only {ref_endmember$n_samples} reference samples for ",
          "feature {feature_uuid} (min {min_ref_samples}). Skipping."
        ))
        return(tibble())
      }

      ## Join R and L; exclude zero-gradient ions.
      endmembers <-
        ref_endmember$stats |>
        filter(ion %in% leachate_em$L_values$ion) |>
        left_join(leachate_em$L_values |> select(ion, L), by = "ion") |>
        filter(abs(L - R) > .Machine$double.eps * 1000)

      ## Compute LMF per sample.
      .x |>
        group_by(uuid.sample) |>
        group_modify(\(.s, .sid) {
          compute_lmf_for_sample(
            sample_wide = .s,
            endmembers = endmembers,
            high_info_ions = high_info_ions,
            min_high_info_ions = min_high_info_ions,
            rsd_default = rsd_default,
            max_sigma = max_sigma_lsi,
            max_chi2_df = max_chi2_per_df,
            uuid_lmf = uuid_lmf,
            ref_window_start = ref_endmember$window_start,
            ref_window_end = ref_endmember$window_end,
            leach_window_start = leachate_em$window_start,
            leach_window_end = leachate_em$window_end,
            cl_anchor = leachate_em$cl_anchor,
            robust_iterations = robust_iterations,
            robust_threshold_k = robust_threshold_k,
            verbose = verbose,
            datetime_sample = as.character(.s$datetime.sample[[1]])
          )
        }) |>
        ungroup()
    }) |>
    dplyr::ungroup()

  ## ================================================================
  ## Step 6: Bind LMF rows back into the input df.
  ## ================================================================

  lsi_results <- dplyr::mutate(lsi_results, name.analyte = "LMF")

  ## Propagate datetime.sample from input df to LMF rows.
  if ("datetime.sample" %in% names(df)) {
    sample_times <- dplyr::distinct(df, uuid.sample, datetime.sample)
    lsi_results <- dplyr::left_join(lsi_results, sample_times, by = "uuid.sample")
  }

  result <- dplyr::bind_rows(df, lsi_results)

  if ("datetime.sample" %in% names(result)) {
    result <- dplyr::arrange(result, datetime.sample)
  }
  result
}


## ============================================================================
## compute_lmf_for_sample() — core per-sample calculation
## ============================================================================
##
## All early-return paths call the internal fail() helper which delegates to
## make_lmf_row(), ensuring structurally identical output across all paths.

#' Compute LMF for a single sample
#'
#' Internal function called once per sample by \code{\link{add_lmf}}.
#' Computes the inverse-variance weighted mixing fraction and associated
#' quality metrics, applying admission and quality gates before returning.
#'
#' @param sample_wide One-row wide-format tibble with one column per measured
#'   (collapsed) panel ion.
#' @param endmembers Tibble with columns \code{ion}, \code{R},
#'   \code{sigma_R}, \code{n_ref}, \code{L}. One row per ion with usable
#'   end-member data for this feature.
#' @param high_info_ions Character vector of ion names classified as
#'   high-information at calibration time.
#' @param min_high_info_ions Minimum count of \code{high_info_ions} that must
#'   be present and non-NA in this sample.
#' @param rsd_default Default relative analytical uncertainty (fraction).
#' @param max_sigma Maximum permitted \code{sigma_lmf}.
#' @param max_chi2_df Maximum permitted chi-squared per degree of freedom.
#' @param uuid_lmf UUID of the LMF analyte entry in \code{analyteDF}.
#' @param ref_window_start,ref_window_end Calibration window dates for the
#'   reference end-member; stored in guideline comment strings.
#' @param leach_window_start,leach_window_end Calibration window dates for the
#'   leachate end-member; stored in guideline comment strings.
#' @param cl_anchor Cl concentration (meq/L) used to anchor the leachate
#'   end-member; stored in guideline comment strings.
#' @param robust_iterations Number of Huber reweighting passes. Inherited from
#'   \code{\link{add_lmf}}.
#' @param robust_threshold_k MAD multiplier for Huber downweighting. Inherited
#'   from \code{\link{add_lmf}}.
#' @param verbose If \code{TRUE}, emits a per-ion diagnostic table via
#'   \code{cli::cli_inform()} showing observed concentration, end-member
#'   values, per-ion mixing fraction, original and robust weights. Inherited
#'   from \code{\link{add_lmf}}. Default \code{FALSE}.
#' @param datetime_sample The sample datetime as a character string, used
#'   to label the per-ion diagnostic table when \code{verbose = TRUE}.
#'   Default \code{NA_character_}.
#'
#' @return A one-row tibble. On success: \code{value} is the robust LMF
#'   estimate, \code{quantified = TRUE}. On failure: \code{value = NA},
#'   \code{quantified = FALSE}, \code{lmf_reason} carries a descriptive
#'   reason code. All paths return identical column structure.
#'
#' @keywords internal

compute_lmf_for_sample <- function(
  sample_wide,
  endmembers,
  high_info_ions,
  min_high_info_ions,
  rsd_default,
  max_sigma,
  max_chi2_df,
  uuid_lmf,
  ref_window_start,
  ref_window_end,
  leach_window_start,
  leach_window_end,
  cl_anchor,
  robust_iterations = 3L,
  robust_threshold_k = 1.5,
  verbose = FALSE,
  datetime_sample = NA_character_
) {
  ## Helper: construct a failure row without repeating all columns.
  fail <- function(reason, n_ions = 0L, sigma = NA_real_, chi2 = NA_real_) {
    make_lmf_row(
      value = NA_real_,
      lmf_naive = NA_real_,
      reason = reason,
      n_ions = n_ions,
      n_downweighted = NA_integer_,
      sigma = sigma,
      chi2 = chi2,
      quantified = FALSE,
      uuid_lmf = uuid_lmf,
      ref_window_start = ref_window_start,
      ref_window_end = ref_window_end,
      leach_window_start = leach_window_start,
      leach_window_end = leach_window_end,
      cl_anchor = cl_anchor
    )
  }

  ## ---------------------------------------------------------------
  ## Admission gate: minimum high-information ions present in sample.
  ##
  ## "Available" means: in the endmembers table for this feature AND
  ## present as a non-NA column in this sample.
  ## ---------------------------------------------------------------

  ## Helper to safely extract a value from a wide row by column name.
  get_val <- function(nm) {
    if (nm %in% names(sample_wide)) sample_wide[[nm]] else NA_real_
  }

  hi_in_endmembers <- high_info_ions[high_info_ions %in% endmembers$ion]
  hi_vals <- purrr::map_dbl(hi_in_endmembers, get_val)
  available_hi <- hi_in_endmembers[!is.na(hi_vals)]
  n_hi_available <- length(available_hi)

  if (n_hi_available < min_high_info_ions) {
    return(fail(glue::glue(
      "insufficient_high_info_ions: {n_hi_available} of ",
      "{min_high_info_ions} required. ",
      "Available: [{paste(available_hi, collapse=', ')}]. ",
      "High-info ions in panel: [{paste(high_info_ions, collapse=', ')}]."
    )))
  }

  ## ---------------------------------------------------------------
  ## Collect all available panel ions (high and low informativeness)
  ## with non-NA values. All contribute via inverse-variance weights.
  ## ---------------------------------------------------------------

  all_vals <- purrr::map_dbl(endmembers$ion, get_val)
  available_mask <- !is.na(all_vals)

  ion_data <-
    endmembers[available_mask, ] |>
    mutate(
      x = all_vals[available_mask],

      ## Per-ion mixing fraction.
      ## f_i = (x_i - R_i) / (L_i - R_i)
      ## Values outside [0, 1] are permitted.
      f = (x - R) / (L - R),

      ## Analytical uncertainty: rsd_default * |x|, floor at
      ## rsd_default * |R| (approximates half the detection limit
      ## near reference concentrations).
      sigma_meas = pmax(rsd_default * abs(x), rsd_default * abs(R)),

      ## Per-ion variance on the mixing fraction:
      ## sigma2_f_i = (sigma2_meas + sigma2_R) / (L - R)^2
      var_f = (sigma_meas^2 + sigma_R^2) / (L - R)^2,
      weight = 1 / var_f
    ) |>
    filter(is.finite(weight), weight > 0)

  n_ions <- nrow(ion_data)

  if (n_ions == 0) {
    return(fail("no_valid_ions_after_variance_check"))
  }

  ## ---------------------------------------------------------------
  ## Pass 0: naive inverse-variance weighted mean.
  ##
  ## LMF = sum(w_i * f_i) / sum(w_i)  where w_i = 1 / sigma2_f_i
  ##
  ## This is the starting point for robust reweighting and is retained
  ## as lmf_naive in the output for comparison with the robust estimate.
  ## All arithmetic stays on the 0-1 scale until the final * 100 step.
  ## ---------------------------------------------------------------

  sum_w <- sum(ion_data$weight)
  lmf_naive <- sum(ion_data$weight * ion_data$f) / sum_w

  ## ---------------------------------------------------------------
  ## chi2/df on original weights: model-fit diagnostic.
  ##
  ## Computed on the naive estimate using original inverse-variance
  ## weights. Retained as a diagnostic output regardless of chi2 gate.
  ## Under the linear mixing model with all ions behaving conservatively,
  ## chi2/df ~ 1. High values flag ion disagreement (transformation
  ## activity, alternative sources) but do NOT suppress the result —
  ## the robust reweighting below handles outlier ions directly.
  ## ---------------------------------------------------------------

  df_resid <- n_ions - 1L
  chi2 <- sum(ion_data$weight * (ion_data$f - lmf_naive)^2)
  chi2_df <- if (df_resid > 0L) chi2 / df_resid else NA_real_

  ## ---------------------------------------------------------------
  ## Robust reweighting: Huber M-estimator.
  ##
  ## Iteratively downweights ions whose residuals are large relative
  ## to the spread of residuals across all ions, so that process-
  ## specific deviations (denitrification in total_N, sulphate
  ## reduction in SO4, young leachate SO4 spike, etc.) do not
  ## systematically bias the LMF estimate.
  ##
  ## At each iteration:
  ##   1. Compute residuals: r_i = f_i - LMF_current
  ##   2. Compute the unweighted MAD of |r_i| across all ions.
  ##      UNWEIGHTED because outlier detection is a question about
  ##      the distribution of residuals, not about estimation
  ##      reliability. An ion's informativeness (high or low weight)
  ##      is irrelevant to whether its residual is unusual.
  ##   3. Apply Huber downweighting:
  ##        w_robust_i = w_original_i * min(1, k * MAD / |r_i|)
  ##      Ions within k*MAD of LMF keep full weight. Ions beyond
  ##      are downweighted proportionally to their exceedance.
  ##      This is continuous — an ion twice the threshold gets half
  ##      its original weight, not zero.
  ##   4. Recompute LMF with the adjusted weights.
  ##
  ## Reference: Huber (1964) Annals of Mathematical Statistics 35(1).
  ## robust_threshold_k = 1.5 is fairly tight; 2.5 is closer to the
  ## conventional Huber constant used in regression.
  ## ---------------------------------------------------------------

  robust_weights <- ion_data$weight ## initialise at original weights
  lmf_robust <- lmf_naive ## initialise at naive estimate

  for (iter in seq_len(max(robust_iterations, 0L))) {
    residuals <- ion_data$f - lmf_robust

    ## Unweighted MAD of absolute residuals.
    ## Add a small floor to avoid division-by-zero when all ions agree
    ## perfectly (MAD = 0, which would downweight every non-zero residual
    ## to zero). The floor of 1e-6 is negligible on the 0-1 scale.
    mad_resid <- max(median(abs(residuals)), 1e-6)
    threshold <- robust_threshold_k * mad_resid

    ## Huber influence function: full weight within threshold,
    ## proportional downweight beyond.
    robust_weights <- ion_data$weight *
      pmin(1, threshold / pmax(abs(residuals), .Machine$double.eps))

    sum_w_robust <- sum(robust_weights)
    lmf_robust <- sum(robust_weights * ion_data$f) / sum_w_robust
  }

  ## Count how many ions had their weight meaningfully reduced.
  ## Threshold: weight reduced by more than 1% of original.
  n_downweighted <- sum(robust_weights < 0.99 * ion_data$weight)

  ## sigma_LMF from the robust weights.
  ## 1 / sqrt(sum(w_robust)) reflects the effective precision after
  ## downweighting outlier ions.
  sum_w_robust <- sum(robust_weights)
  sigma_lmf <- 1 / sqrt(sum_w_robust)

  ## ---------------------------------------------------------------
  ## Scale to percentage.
  ## All internal arithmetic was on the 0-1 scale; convert now.
  ## ---------------------------------------------------------------
  lmf_value <- lmf_robust * 100
  lmf_naive <- lmf_naive * 100
  sigma_lmf <- sigma_lmf * 100

  ## ---------------------------------------------------------------
  ## Verbose per-ion diagnostic table.
  ##
  ## Emitted once per sample when verbose = TRUE. Shows for each ion:
  ##   x     observed concentration (meq/L)
  ##   R     reference end-member mean (meq/L)
  ##   L     leachate end-member (meq/L)
  ##   f%    per-ion mixing fraction in percentage units
  ##   wt%   this ion's share of total inverse-variance weight
  ##   hi?   [H] = high-information, [-] = low-information
  ##
  ## f% values far from the aggregate LMF indicate an ion is not
  ## conforming to the two-component mixing model (e.g. SO4 under
  ## sulphate reduction, or N under denitrification). wt% shows
  ## which ions dominate the aggregate estimate.
  ## ---------------------------------------------------------------
  if (verbose) {
    ## Normalised robust weights for display (sum to 100%).
    sum_w_orig <- sum(ion_data$weight)

    diag_rows <-
      ion_data |>
      mutate(
        f_pct = f * 100,
        wt_orig_pct = weight / sum_w_orig * 100,
        wt_rob_pct = robust_weights / sum_w_robust * 100,
        hi_flag = if_else(ion %in% high_info_ions, "[H]", "[-]"),
        ## Flag ions that were meaningfully downweighted
        dw_flag = if_else(robust_weights < 0.99 * weight, "*", " ")
      ) |>
      arrange(desc(wt_rob_pct)) |>
      mutate(
        row_str = sprintf(
          "  %-14s  x=%7.3f  R=%7.3f  L=%8.3f  f%%=%6.1f  wt%%=%5.1f->%5.1f  %s%s",
          ion,
          x,
          R,
          L,
          f_pct,
          wt_orig_pct,
          wt_rob_pct,
          hi_flag,
          dw_flag
        )
      ) |>
      pull(row_str)

    dt_label <- if (!is.na(datetime_sample)) datetime_sample else "unknown date"

    cli::cli_inform(c(
      "i" = paste0(
        "Per-ion breakdown  [",
        dt_label,
        "]",
        "  LMF = ",
        round(lmf_value, 1),
        "%",
        "  (naive = ",
        round(lmf_naive, 1),
        "%)",
        "  sigma = ",
        round(sigma_lmf, 2),
        "  chi2/df = ",
        if (!is.na(chi2_df)) round(chi2_df, 2) else "NA",
        if (n_downweighted > 0) {
          paste0("  [", n_downweighted, " ion(s) downweighted *]")
        } else {
          ""
        }
      ),
      sprintf(
        "  %-14s  %9s  %9s  %9s  %7s  %13s  %s",
        "ion",
        "x(meq/L)",
        "R(meq/L)",
        "L(meq/L)",
        "f%",
        "wt%(orig->rob)",
        "hi?"
      ),
      diag_rows
    ))
  }

  ## ---------------------------------------------------------------
  ## Quality gates.
  ## ---------------------------------------------------------------

  if (!is.na(sigma_lmf) && sigma_lmf > max_sigma) {
    return(fail(
      glue::glue(
        "insufficient_precision: sigma_lmf = {round(sigma_lmf, 4)}, ",
        "threshold = {max_sigma}."
      ),
      n_ions = n_ions,
      sigma = sigma_lmf,
      chi2 = chi2_df
    ))
  }

  if (!is.na(chi2_df) && chi2_df > max_chi2_df) {
    return(fail(
      glue::glue(
        "poor_fit: chi2_per_df = {round(chi2_df, 2)}, ",
        "threshold = {max_chi2_df}. ",
        "Possible causes: SO4 reduction, partial denitrification, or ",
        "alternative water source. Investigate per-ion mixing fractions."
      ),
      n_ions = n_ions,
      sigma = sigma_lmf,
      chi2 = chi2_df
    ))
  }

  ## ---------------------------------------------------------------
  ## Successful result.
  ## ---------------------------------------------------------------

  make_lmf_row(
    value = lmf_value,
    lmf_naive = lmf_naive,
    reason = NA_character_,
    n_ions = n_ions,
    n_downweighted = n_downweighted,
    sigma = sigma_lmf,
    chi2 = chi2_df,
    quantified = TRUE,
    uuid_lmf = uuid_lmf,
    ref_window_start = ref_window_start,
    ref_window_end = ref_window_end,
    leach_window_start = leach_window_start,
    leach_window_end = leach_window_end,
    cl_anchor = cl_anchor
  )
}


## ============================================================================
## make_lmf_row() — standardised output constructor
## ============================================================================
##
## Factored out of compute_lmf_for_sample() so that every code path (success
## and all failure modes) produces structurally identical output, ensuring
## downstream bind_rows() calls are safe.

#' Construct a standardised LMF output row
#'
#' Internal helper that builds the one-row tibble returned by
#' \code{\link{compute_lmf_for_sample}} for both successful and failed
#' computations. All paths call this function so the output schema is
#' guaranteed consistent.
#'
#' @param value Numeric robust LMF estimate, or \code{NA_real_} on failure.
#' @param lmf_naive Numeric naive (non-robust) LMF estimate for comparison.
#' @param reason Character reason code, or \code{NA_character_} on success.
#' @param n_ions Integer count of ions used in the calculation.
#' @param n_downweighted Integer count of ions whose weight was meaningfully
#'   reduced (> 1\% reduction) by robust reweighting.
#' @param sigma Numeric \code{sigma_lmf} uncertainty estimate (from robust
#'   weights).
#' @param chi2 Numeric chi-squared per degree of freedom (on original weights;
#'   diagnostic only).
#' @param quantified Logical; \code{TRUE} on success, \code{FALSE} on failure.
#' @param uuid_lmf UUID of the LMF analyte entry in \code{analyteDF}.
#' @param ref_window_start,ref_window_end Reference calibration window dates.
#' @param leach_window_start,leach_window_end Leachate calibration window dates.
#' @param cl_anchor Cl anchor concentration (meq/L) for the leachate end-member.
#'
#' @return A one-row tibble with columns: \code{value} (robust LMF),
#'   \code{lmf_naive}, \code{lmf_reason}, \code{n_ions_used},
#'   \code{n_ions_downweighted}, \code{sigma_lmf}, \code{chi2_per_df},
#'   \code{uuid.analyte}, \code{uuid}, \code{quantified}, and four sets of
#'   guideline columns (\code{value/level_name/guideline/comments.guideline_1}
#'   through \code{_4}).
#'
#' @keywords internal

make_lmf_row <- function(
  value,
  lmf_naive,
  reason,
  n_ions,
  n_downweighted,
  sigma,
  chi2,
  quantified,
  uuid_lmf,
  ref_window_start,
  ref_window_end,
  leach_window_start,
  leach_window_end,
  cl_anchor
) {
  tibble(
    value = value,
    lmf_naive = lmf_naive,
    units.analyte = "%",
    lmf_reason = reason,
    n_ions_used = as.integer(n_ions),
    n_ions_downweighted = as.integer(n_downweighted),
    sigma_lmf = sigma,
    chi2_per_df = chi2,
    uuid.analyte = uuid_lmf,
    uuid = uuid::UUIDgenerate(),
    quantified = quantified,

    ## Tier breaks: mixing-fraction units, reviewed after deployment.
    value.guideline_1 = .LMF_BREAK_TRACE,
    level_name.guideline_1 = "Background",
    guideline.guideline_1 = "LMF threshold",
    comments.guideline_1 = glue::glue(
      "LMF <= {.LMF_BREAK_TRACE}%: chemistry indistinguishable from local ",
      "reference water. Reference end-member: {ref_window_start} to ",
      "{ref_window_end}."
    ),

    value.guideline_2 = .LMF_BREAK_MODERATE,
    level_name.guideline_2 = "Trace impact",
    guideline.guideline_2 = "LMF threshold",
    comments.guideline_2 = glue::glue(
      "LMF {.LMF_BREAK_TRACE}-{.LMF_BREAK_MODERATE}%: detectable leachate ",
      "signature. Monitor trend."
    ),

    value.guideline_3 = .LMF_BREAK_SEVERE,
    level_name.guideline_3 = "Significant impact",
    guideline.guideline_3 = "LMF threshold",
    comments.guideline_3 = glue::glue(
      "LMF {.LMF_BREAK_MODERATE}-{.LMF_BREAK_SEVERE}%: clear leachate ",
      "signature. Warrants investigation. Leachate end-member Cl anchor = ",
      "{round(cl_anchor, 1)} meq/L, {leach_window_start} to {leach_window_end}."
    ),

    value.guideline_4 = .LMF_BREAK_CEILING,
    level_name.guideline_4 = "Severe impact",
    guideline.guideline_4 = "LMF threshold",
    comments.guideline_4 = glue::glue(
      "LMF > {.LMF_BREAK_SEVERE}%: severe leachate impact. ",
      "Urgent investigation required. ",
      "Values > 100 indicate sample chemistry exceeds the leachate ",
      "end-member, possibly indicating proximity to a leachate source."
    )
  )
}


## ============================================================================
## collapse_species() — species aggregation helper
## ============================================================================
##
## Uses rowSums with na.rm = TRUE so partial N or alkalinity measurements are
## still used. If ALL constituent species for a group are NA in a row, the
## total is set back to NA (no data) rather than leaving a spurious zero.

#' Pivot and collapse transforming ion species to conserved totals
#'
#' Internal helper used by \code{\link{add_lmf}},
#' \code{\link{build_reference_endmember}}, and
#' \code{\link{build_leachate_endmember}} to avoid code duplication.
#' Pivots long-format meq data to wide format and replaces the individual N
#' species (NH4-N, NO3-N, NO2-N) with \code{total_N_} and the carbonate
#' species (CO3, HCO3) with \code{total_alk_}.
#'
#' @param df_meq Long-format dataframe containing at minimum columns
#'   \code{name.analyte} and \code{value}, plus any columns named in
#'   \code{id_cols}.
#' @param id_cols Character vector of column names to preserve as row
#'   identifiers during the wide pivot.
#'
#' @return Wide-format tibble with one row per unique combination of
#'   \code{id_cols} values. Individual N and carbonate species columns are
#'   replaced by \code{total_N_} and \code{total_alk_} respectively.
#'
#' @keywords internal

collapse_species <- function(df_meq, id_cols) {
  wide <-
    df_meq |>
    select(all_of(c(id_cols, "name.analyte", "value"))) |>
    pivot_wider(names_from = name.analyte, values_from = value)

  ## Collapse N species to total inorganic N.
  n_cols <- c("NH3-N_", "NO3-N_", "NO2-N_")
  n_present <- intersect(n_cols, names(wide))

  wide <- wide |>
    mutate(
      total_N_ = if (length(n_present) > 0) {
        rowSums(pick(all_of(n_present)), na.rm = TRUE)
      } else {
        NA_real_
      },
      total_N_ = if (length(n_present) > 0) {
        if_else(
          rowSums(!is.na(pick(all_of(n_present)))) == 0,
          NA_real_,
          total_N_
        )
      } else {
        NA_real_
      }
    )

  ## Collapse carbonate species to total alkalinity.
  alk_cols <- c("CO3-CaCO3_", "HCO3-CaCO3_")
  alk_present <- intersect(alk_cols, names(wide))

  wide <- wide |>
    mutate(
      total_alk_ = if (length(alk_present) > 0) {
        rowSums(pick(all_of(alk_present)), na.rm = TRUE)
      } else {
        NA_real_
      },
      total_alk_ = if (length(alk_present) > 0) {
        if_else(
          rowSums(!is.na(pick(all_of(alk_present)))) == 0,
          NA_real_,
          total_alk_
        )
      } else {
        NA_real_
      }
    )

  wide |>
    select(-any_of(c(n_cols, alk_cols)))
}


## ============================================================================
## build_endmember_from_override() — end-member builder for override data
## ============================================================================
##
## Processes a raw data_df()-structured dataframe into the same list format
## produced by build_reference_endmember() or build_leachate_endmember(),
## allowing callers to supply their own chemistry data for specific events
## instead of using the standard calibration window and feature-matching logic.
##
## All samples in the override dataframe are used regardless of date.
## Samples missing the LMF panel analytes are silently dropped (they would
## contribute nothing to the end-member anyway).
##
## The total-N quality filter is NOT applied for leachate overrides. The
## caller is responsible for supplying curated data.

#' Build an end-member from a caller-supplied override dataframe
#'
#' Internal helper called by \code{\link{add_lmf}} when
#' \code{reference_data} or \code{leachate_data} is
#' supplied. Processes the override data through the same meq conversion,
#' BDL substitution, and species-collapsing machinery as the standard
#' builders, but uses all available samples without date filtering.
#'
#' @param override_df Long-format dataframe in \code{data_df()} structure.
#'   All samples are used; samples missing the LMF panel are dropped.
#' @param type Character scalar, either \code{"reference"} or
#'   \code{"leachate"}, controlling which output list format is produced.
#'
#' @return For \code{type = "reference"}: a list with the same structure as
#'   \code{\link{build_reference_endmember}}: \code{stats} (tibble of ion,
#'   n_ref, R, sigma_R), \code{window_start}, \code{window_end},
#'   \code{n_samples}.
#'   For \code{type = "leachate"}: a list with the same structure as
#'   \code{\link{build_leachate_endmember}}: \code{L_values} (tibble of
#'   ion, mean_ratio, L), \code{cl_anchor}, \code{n_samples},
#'   \code{f_included}, \code{window_start}, \code{window_end}.
#'
#' @keywords internal

build_endmember_from_override <- function(override_df, type) {
  checkmate::assert_data_frame(override_df)
  checkmate::assert_choice(type, c("reference", "leachate"))

  ## Convert, apply BDL half-substitution, collapse species.
  ## Samples that don't survive the meq filter or lack any panel ions
  ## drop out naturally — no explicit check needed.
  override_meq <-
    override_df |>
    to_meq() |>
    filter(name.analyte %in% .LMF_ANALYTES_MEQ) |>
    mutate(value = if_else(!quantified, value * 0.5, value))

  if (nrow(override_meq) == 0) {
    stop(glue::glue(
      "build_endmember_from_override: no LMF panel analytes found in the ",
      "{type} override dataframe after meq conversion. Check that the ",
      "override contains the expected analytes: ",
      "{paste(.LMF_ANALYTES_MEQ, collapse = ', ')}."
    ))
  }

  override_wide <- collapse_species(
    df_meq = override_meq,
    id_cols = "uuid.sample"
  )

  ## Date range for use in window_start/window_end fields.
  ## These are informational only when using an override — they describe
  ## the span of the override data, not a calibration window.
  dates <- as.Date(override_df$datetime.sample)
  window_start <- min(dates, na.rm = TRUE)
  window_end <- max(dates, na.rm = TRUE)
  n_samples <- nrow(override_wide)

  if (type == "reference") {
    ## ---------------------------------------------------------------
    ## Reference end-member: mean and SD per ion.
    ## Require >= 3 observations per ion for a meaningful SD, consistent
    ## with build_reference_endmember().
    ## ---------------------------------------------------------------
    stats <-
      override_wide |>
      select(all_of(intersect(.LMF_PANEL_COLLAPSED, names(override_wide)))) |>
      pivot_longer(everything(), names_to = "ion", values_to = "value") |>
      filter(!is.na(value)) |>
      group_by(ion) |>
      summarise(
        n_ref = n(),
        R = mean(value, na.rm = TRUE),
        sigma_R = sd(value, na.rm = TRUE),
        .groups = "drop"
      ) |>
      filter(n_ref >= 3)

    if (nrow(stats) == 0) {
      stop(paste0(
        "build_endmember_from_override: no ions have >= 3 observations in ",
        "the reference override dataframe. Cannot build a reference ",
        "end-member. Check the override data."
      ))
    }

    return(list(
      stats = stats,
      window_start = window_start,
      window_end = window_end,
      n_samples = n_samples
    ))
  } else {
    ## ---------------------------------------------------------------
    ## Leachate end-member: Cl-anchored ratios.
    ## Mirrors build_leachate_endmember() but without the total-N
    ## quality filter (caller supplies curated data).
    ## ---------------------------------------------------------------
    cl_anchor <- median(override_wide$Cl_, na.rm = TRUE)

    if (is.na(cl_anchor) || cl_anchor <= 0) {
      stop(paste0(
        "build_endmember_from_override: cannot compute a Cl anchor from ",
        "the leachate override dataframe. Check that Cl is measured."
      ))
    }

    ratio_cols <- setdiff(
      intersect(.LMF_PANEL_COLLAPSED, names(override_wide)),
      "Cl_"
    )

    leachate_ratios <-
      override_wide |>
      filter(!is.na(Cl_), Cl_ > 0) |>
      mutate(across(all_of(ratio_cols), \(x) x / Cl_, .names = "ratio_{.col}"))

    ## F availability: same threshold as standard builder.
    f_col <- "ratio_F_"
    f_available <- f_col %in%
      names(leachate_ratios) &&
      sum(!is.na(leachate_ratios[[f_col]])) >= 3

    L_values <-
      leachate_ratios |>
      select(starts_with("ratio_")) |>
      summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) |>
      pivot_longer(everything(), names_to = "ion", values_to = "mean_ratio") |>
      mutate(
        ion = stringr::str_remove(ion, "^ratio_"),
        L = mean_ratio * cl_anchor
      ) |>
      filter(!is.na(L)) |>
      filter(!(ion == "F_" & !f_available))

    L_values <- bind_rows(
      L_values,
      tibble(ion = "Cl_", mean_ratio = 1, L = cl_anchor)
    )

    return(list(
      L_values = L_values,
      cl_anchor = cl_anchor,
      n_samples = n_samples,
      f_included = f_available,
      window_start = window_start,
      window_end = window_end
    ))
  }
}


## ============================================================================
## build_reference_endmember() — per-feature reference end-member
## ============================================================================
##
## Applies the same backwards-window-shift logic as MTUI: if reference data
## do not extend to the ideal forward edge of the calibration window (common
## when processing recent data), the window is shifted backwards so it ends at
## the latest available reference sample while retaining the full window width.
## Requires at least 3 observations per ion for a meaningful SD estimate.

#' Build a per-feature reference end-member for LMF calibration
#'
#' Computes mean (\code{R_i}) and standard deviation (\code{sigma_R_i}) for
#' each ion in the collapsed LMF panel from reference samples within the
#' calibration window. Called once per feature by \code{\link{add_lmf}}.
#'
#' @param reference_feature_uuid UUID of the matched reference feature.
#' @param cal_start Ideal calibration window start date.
#' @param cal_end Ideal calibration window end date. Shifted backwards if
#'   reference data do not extend this far.
#'
#' @return \code{NULL} if no reference data are found. Otherwise a named list
#'   with elements:
#'   \describe{
#'     \item{stats}{Tibble with columns \code{ion}, \code{n_ref}, \code{R},
#'       \code{sigma_R}. Only ions with >= 3 observations are included.}
#'     \item{window_start}{Actual calibration window start used.}
#'     \item{window_end}{Actual calibration window end used.}
#'     \item{n_samples}{Number of sample events contributing.}
#'   }
#'
#' @keywords internal

build_reference_endmember <- function(
  reference_feature_uuid,
  cal_start,
  cal_end
) {
  all_ref <-
    data_df() |>
    filter(uuid.feature == reference_feature_uuid) |>
    to_meq() |>
    filter(name.analyte %in% .LMF_ANALYTES_MEQ) |>
    mutate(value = if_else(!quantified, value * 0.5, value))

  if (nrow(all_ref) == 0) {
    return(NULL)
  }

  ## Backwards window shift if no future data.
  ref_max_date <- max(as.Date(all_ref$datetime.sample), na.rm = TRUE)
  window_start <- cal_start
  window_end <- cal_end

  if (ref_max_date < window_end) {
    shift_days <- as.numeric(difftime(window_end, ref_max_date, units = "days"))
    window_start <- window_start - lubridate::days(round(shift_days))
    window_end <- ref_max_date
  }

  ref_window <-
    all_ref |>
    filter(
      as.Date(datetime.sample) >= window_start,
      as.Date(datetime.sample) <= window_end
    )

  ref_wide <- collapse_species(ref_window, id_cols = "uuid.sample")

  R_stats <-
    ref_wide |>
    select(all_of(intersect(.LMF_PANEL_COLLAPSED, names(ref_wide)))) |>
    pivot_longer(everything(), names_to = "ion", values_to = "value") |>
    filter(!is.na(value)) |>
    group_by(ion) |>
    summarise(
      n_ref = n(),
      R = mean(value, na.rm = TRUE),
      sigma_R = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(n_ref >= 3)

  list(
    stats = R_stats,
    window_start = window_start,
    window_end = window_end,
    n_samples = nrow(ref_wide)
  )
}


## ============================================================================
## build_pooled_reference_endmember() — system-wide pooled reference
## ============================================================================
##
## Used only for computing per-ion informativeness in Step 4 of add_lmf().
## Draws on all reference features rather than a single matched site because
## informativeness reflects system-wide chemistry properties, not site-specific
## conditions.

#' Build a system-wide pooled reference end-member
#'
#' Pools data from all reference features to compute per-ion mean and standard
#' deviation. Used exclusively for computing ion informativeness scores in
#' \code{\link{add_lmf}}; per-feature reference end-members are built by
#' \code{\link{build_reference_endmember}}.
#'
#' @param calibration_start Calibration window start date.
#' @param calibration_end Calibration window end date. Shifted backwards if
#'   reference data do not extend this far.
#'
#' @return Tibble with columns \code{ion}, \code{n_ref}, \code{R},
#'   \code{sigma_R}. Only ions with >= 3 observations across all reference
#'   features are included.
#'
#' @keywords internal

build_pooled_reference_endmember <- function(
  calibration_start,
  calibration_end
) {
  ref_uuids <-
    feature_sfc() |>
    filter(reference == TRUE) |>
    pull(uuid)

  if (length(ref_uuids) == 0) {
    stop(
      "No reference features found. Check reference == TRUE in feature_sfc()."
    )
  }

  all_ref <-
    data_df() |>
    filter(uuid.feature %in% ref_uuids) |>
    to_meq() |>
    filter(name.analyte %in% .LMF_ANALYTES_MEQ) |>
    mutate(value = if_else(!quantified, value * 0.5, value))

  if (nrow(all_ref) == 0) {
    stop("No reference data found for any reference features.")
  }

  ref_max_date <- max(as.Date(all_ref$datetime.sample), na.rm = TRUE)
  window_start <- calibration_start
  window_end <- calibration_end

  if (ref_max_date < window_end) {
    shift_days <- as.numeric(difftime(window_end, ref_max_date, units = "days"))
    window_start <- window_start - lubridate::days(round(shift_days))
    window_end <- ref_max_date
  }

  ref_window <-
    all_ref |>
    filter(
      as.Date(datetime.sample) >= window_start,
      as.Date(datetime.sample) <= window_end
    )

  ref_wide <- collapse_species(ref_window, id_cols = "uuid.sample")

  ref_wide |>
    select(all_of(intersect(.LMF_PANEL_COLLAPSED, names(ref_wide)))) |>
    pivot_longer(everything(), names_to = "ion", values_to = "value") |>
    filter(!is.na(value)) |>
    group_by(ion) |>
    summarise(
      n_ref = n(),
      R = mean(value, na.rm = TRUE),
      sigma_R = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(n_ref >= 3)
}


## ============================================================================
## build_leachate_endmember() — Cl-anchored leachate end-member
## ============================================================================
##
## Leachate features are identified by ".L" + digit in their feature name
## (e.g., "B.L01"). F is included only if >= min_leachate_samples leachate
## samples have F measurements; otherwise it is excluded with a message.
## See the file-level header for full Cl-anchoring rationale.

#' Build the Cl-anchored leachate end-member for LMF calibration
#'
#' Constructs the leachate end-member L by computing the mean ratio of each
#' ion to Cl across valid leachate samples, then anchoring to the median Cl
#' concentration. This makes the end-member robust to EC variability between
#' leachate samples while preserving the stable compositional fingerprint.
#'
#' @param calibration_start Calibration window start date.
#' @param calibration_end Calibration window end date. Shifted backwards if
#'   leachate data do not extend this far.
#' @param min_leachate_total_n_mgl Minimum total N (mg/L) for a leachate
#'   sample to be considered valid. Applied in original mg/L units before meq
#'   conversion. Filters out anomalous samples not representing genuine
#'   leachate.
#' @param min_leachate_samples Minimum number of valid leachate samples
#'   required. Stops with an error if not met.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{L_values}{Tibble with columns \code{ion}, \code{mean_ratio},
#'       \code{L} (the end-member concentration in meq/L for each ion).}
#'     \item{cl_anchor}{Median Cl concentration (meq/L) used to anchor the
#'       end-member.}
#'     \item{n_samples}{Number of valid leachate samples used.}
#'     \item{f_included}{Logical; \code{TRUE} if F was included in the
#'       end-member.}
#'     \item{window_start}{Actual calibration window start used.}
#'     \item{window_end}{Actual calibration window end used.}
#'   }
#'
#' @keywords internal

build_leachate_endmember <- function(
  calibration_start,
  calibration_end,
  min_leachate_total_n_mgl,
  min_leachate_samples
) {
  leachate_features <-
    feature_df() |>
    filter(stringr::str_detect(name, "\\.L\\d"))

  if (nrow(leachate_features) == 0) {
    stop(paste0(
      "No leachate features found. Expected names matching '.L01', '.L02', ",
      "etc. (e.g., 'B.L01'). Check feature naming."
    ))
  }

  leachate_raw <-
    data_df() |>
    filter(uuid.feature %in% leachate_features$uuid)

  if (nrow(leachate_raw) == 0) {
    stop("No leachate data found for identified leachate features.")
  }

  leach_max_date <- max(as.Date(leachate_raw$datetime.sample), na.rm = TRUE)
  window_start <- calibration_start
  window_end <- calibration_end

  if (leach_max_date < window_end) {
    shift_days <- as.numeric(difftime(
      window_end,
      leach_max_date,
      units = "days"
    ))
    window_start <- window_start - lubridate::days(round(shift_days))
    window_end <- leach_max_date
  }

  leachate_window <-
    leachate_raw |>
    filter(
      as.Date(datetime.sample) >= window_start,
      as.Date(datetime.sample) <= window_end
    )

  ## ---------------------------------------------------------------
  ## Total-N quality filter in mg/L (original units, not meq).
  ## Using total N (NH4 + NO3 + NO2) not just NH4, because downstream
  ## aerated leachate features may have nitrified NH4 to NO3.
  ## ---------------------------------------------------------------

  valid_samples <-
    leachate_window |>
    filter(name.analyte %in% c("NH3-N", "NO3-N", "NO2-N")) |>
    group_by(uuid.sample) |>
    summarise(total_N_mgl = sum(value, na.rm = TRUE), .groups = "drop") |>
    filter(total_N_mgl >= min_leachate_total_n_mgl) |>
    pull(uuid.sample)

  leachate_valid <-
    leachate_window |>
    filter(uuid.sample %in% valid_samples)

  n_valid <- n_distinct(leachate_valid$uuid.sample)

  if (n_valid < min_leachate_samples) {
    stop(glue::glue(
      "Only {n_valid} valid leachate samples ({window_start} to ",
      "{window_end}, total N >= {min_leachate_total_n_mgl} mg/L). ",
      "Minimum: {min_leachate_samples}. ",
      "Extend the window, lower min_leachate_total_n_mgl, or add more ",
      "leachate sampling."
    ))
  }

  ## Convert, BDL, collapse.
  leachate_wide <-
    leachate_valid |>
    to_meq() |>
    filter(name.analyte %in% .LMF_ANALYTES_MEQ) |>
    mutate(value = if_else(!quantified, value * 0.5, value)) |>
    collapse_species(id_cols = "uuid.sample")

  ## Cl anchor: median Cl across valid leachate samples.
  cl_anchor <- median(leachate_wide$Cl_, na.rm = TRUE)

  if (is.na(cl_anchor) || cl_anchor <= 0) {
    stop(
      "Cannot compute Cl anchor. Check that Cl is measured in leachate samples."
    )
  }

  ## Compute per-sample ion/Cl ratios.
  ratio_cols <- setdiff(
    intersect(.LMF_PANEL_COLLAPSED, names(leachate_wide)),
    "Cl_"
  )

  leachate_ratios <-
    leachate_wide |>
    filter(!is.na(Cl_), Cl_ > 0) |>
    mutate(across(all_of(ratio_cols), \(x) x / Cl_, .names = "ratio_{.col}"))

  ## F availability check.
  f_available <- FALSE
  f_col <- "ratio_F_"

  if (f_col %in% names(leachate_ratios)) {
    n_f <- sum(!is.na(leachate_ratios[[f_col]]))
    if (n_f >= min_leachate_samples) {
      f_available <- TRUE
    } else {
      message(glue::glue(
        "F measured in only {n_f} leachate samples (min {min_leachate_samples}). ",
        "F excluded from LMF panel."
      ))
    }
  } else {
    message("F not found in leachate samples. F excluded from LMF panel.")
  }

  ## Average ratios and construct L_i values.
  L_values <-
    leachate_ratios |>
    select(starts_with("ratio_")) |>
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) |>
    pivot_longer(everything(), names_to = "ion", values_to = "mean_ratio") |>
    mutate(
      ion = stringr::str_remove(ion, "^ratio_"),
      L = mean_ratio * cl_anchor
    ) |>
    filter(!is.na(L)) |>
    filter(!(ion == "F_" & !f_available))

  ## Add Cl itself (ratio = 1, L_Cl = cl_anchor by definition).
  L_values <- bind_rows(
    L_values,
    tibble(ion = "Cl_", mean_ratio = 1, L = cl_anchor)
  )

  list(
    L_values = L_values,
    cl_anchor = cl_anchor,
    n_samples = n_valid,
    f_included = f_available,
    window_start = window_start,
    window_end = window_end
  )
}


## ============================================================================
## to_meq() — unit conversion
## ============================================================================
##
## Safe to call exactly once per pipeline. The "_" suffix on converted rows
## prevents collisions if bind_rows is used downstream. Rows missing
## valence.analyte or atomic_mass.analyte are silently omitted from the
## converted rows but remain in the originals returned by bind_rows.

#' Convert analyte concentrations to milliequivalents per litre
#'
#' Appends meq/L-converted rows to the input dataframe. Converted rows have
#' \code{name.analyte} suffixed with \code{"_"} to distinguish them from the
#' original rows, which are preserved unchanged.
#'
#' @param df Long-format dataframe. Required columns: \code{name.analyte},
#'   \code{value}, \code{units.analyte}, \code{valence.analyte},
#'   \code{atomic_mass.analyte}.
#'
#' @return The input \code{df} with converted rows appended via
#'   \code{bind_rows}. Converted rows have \code{name.analyte} suffixed with
#'   \code{"_"} and \code{value} in meq/L. Rows lacking \code{valence.analyte}
#'   or \code{atomic_mass.analyte} are excluded from conversion but retained
#'   in the original rows.
#'
#' @export

to_meq <- function(df) {
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c(
      "name.analyte",
      "value",
      "units.analyte",
      "valence.analyte",
      "atomic_mass.analyte"
    )
  )

  ## Split by units.analyte, convert each subset, then recombine.
  ##
  ## The original approach used group_by() + unique(units.analyte) inside
  ## mutate(), but set_units() requires a plain length-1 character scalar.
  ## unique() inside a mutate group can return a vector rather than a scalar
  ## in some contexts, causing:
  ##   "length(x) == 1 is not TRUE"
  ##
  ## The fix: split the dataframe by units, extract the unit string directly
  ## from the split key (a guaranteed scalar), convert each subset, then
  ## recombine. This is equivalent to the grouped approach but robust.

  to_convert <-
    df |>
    filter(!is.na(valence.analyte), !is.na(atomic_mass.analyte))

  if (nrow(to_convert) == 0) {
    return(bind_rows(df))
  }

  converted <-
    to_convert |>
    split(to_convert$units.analyte) |>
    purrr::imap(\(subset, unit_str) {
      subset |>
        mutate(
          name.analyte = glue::glue("{name.analyte}_"),
          value = units::set_units(value, unit_str, mode = "standard"),
          value = value / units::set_units(atomic_mass.analyte, g / mol),
          valence_abs = units::set_units(abs(valence.analyte), absValence),
          value = value * valence_abs,
          value = units::set_units(value, mEq / L),
          value = units::drop_units(value)
        ) |>
        select(-valence_abs)
    }) |>
    purrr::list_rbind()

  bind_rows(df, converted)
}
