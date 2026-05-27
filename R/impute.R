# ── Analyte group constants ───────────────────────────────────────────────────

#' WQ analytes used in the PCA pre-processing step
#'
#' All candidates; the actual variables used at any given site are the
#' intersection of this list with analytes that pass prescreen in the training
#' data.  ORP and DO are included because they are valuable at sites where they
#' are measured, but at most sites they will be filtered out by prescreen.
#'
#' Note: SO4 uses `\u00b2\u207b` (\u00b2\u207b) escape sequences to ensure correct
#' string matching regardless of how the source file is parsed by devtools.
#' @keywords internal
.WQ_BLOCK_CANDIDATES <- c(
  # Field parameters / redox
  "temperature", "ORP", "DO",
  # Major cations
  "Ca", "Mg", "Na", "K",
  # Major anions & alkalinity
  "Cl", "SO4\u00b2\u207b", "Alkalinity-total-CaCO3", "F",
  "Hardness-total-CaCO3",
  # Carbonate / alkalinity species
  "HCO3-CaCO3", "CO3-CaCO3", "OH-CaCO3",
  # Dissolved solids & suspended solids
  "TDS", "TSS",
  # Ionic balance totals
  "Anions-total", "Cations-total",
  # Organic carbon / oxygen demand (also used for organics hurdle)
  "DOC", "TOC", "BOD", "COD", "cBOD",
  # Nitrogen (excluding NH3-N which is a required driver)
  "NO2-N", "NO3-N", "NO2+NO3-N", "TKN-N", "N-total",
  # Sulfur
  "S",
  # Phosphorus
  "P-total", "P-reactive"
)

#' Analytes that must never enter the imputation model as response variables
#'
#' These are excluded from both the metals and organics groups.  They are
#' typically non-concentration measurements (counts, qualitative, physical)
#' for which a log-normal concentration model is inappropriate.
#' @keywords internal
.IMPUTE_EXCLUDED <- c(
  # Microbiological counts (colony-forming units — not concentrations)
  "Coliforms",
  "Escherichia coli",
  "Faecal Coliforms",
  "Heterotrophic Plate Count (22\u00b0C)",
  "Heterotrophic Plate Count (36\u00b0C)",
  "E. coli",
  # Qualitative / physical descriptors
  "Appearance",
  "Colour",
  "Turbidity",
  "Stage"
)

#' All metal-type analytes (used to define the metals imputation group)
#' @keywords internal
.METAL_ANALYTES <- c(
  "Al", "As", "B", "Ba", "Be", "Cd", "Co", "Cr", "Cr-6", "Cu",
  "Fe", "Hg", "Mn", "Mo", "Ni", "Pb", "Sb", "Se", "Sn", "Sr", "V", "Zn"
)

#' Analytes that satisfy the organic-carbon hurdle for the organics model
#' @keywords internal
.DOC_LIKE_ANALYTES <- c("DOC", "TOC", "BOD", "COD", "cBOD")


# ── fit_imputation_model() ────────────────────────────────────────────────────

#' Fit the Bayesian multivariate imputation model(s)
#'
#' Fits one or two `brms` multivariate GAMs — one for metals and one for
#' organics — using a PCA-compressed water-quality (WQ) block as additional
#' environmental predictors.  The returned model object is passed to
#' `impute_chemistry()` for prediction on new data.
#'
#' **Model structure**
#'
#' For each analyte group (metals / organics), the mean structure is:
#' ```
#' s(pH) + s(EC) + s(NH3-N) + s(WQ_PC1) + [s(WQ_PC2) + ...]
#' ```
#' where `WQ_PC*` are the leading principal components of the water-quality
#' block (see *WQ PCA* below).  All target analytes (metals or organics) are
#' modelled jointly with `rescor = TRUE`, so observed co-analytes at a given
#' sample condition the posterior of the missing ones through the residual
#' correlation matrix.
#'
#' **WQ PCA**
#'
#' WQ block variables (major ions, carbon/oxygen demand, nutrients, redox
#' indicators — see `.WQ_BLOCK_CANDIDATES`) that are present in `df` and pass
#' a detection-frequency check are pivoted to a wide matrix, median-imputed for
#' any missing cells, and submitted to `prcomp()`.  Principal components are
#' added until cumulative variance explained reaches `min_var_explained` or
#' `max_pcs` is reached (a warning is issued if the target cannot be met within
#' `max_pcs` axes).  A minimum of two PCs is always used.
#'
#' **Hurdles (applied at prediction time by `impute_chemistry()`)**
#'
#' - *Metals*: a sample is only imputed if at least one metal analyte is
#'   present (detected or BDL) in `df` for that sample.
#' - *Organics*: a sample is only imputed if at least one of
#'   {DOC, TOC, BOD, COD, cBOD} is present.
#'
#' **BDL drivers**
#'
#' When a required driver (pH, EC, NH3-N) is below the detection limit for a
#' sample, the stored detection-limit value is used as-is (conservative upper
#' bound).  A message is issued but the sample is retained — BDL-driver samples
#' represent genuine low-concentration conditions and are important for
#' calibrating the low end of the model.
#'
#' @param df Long-format chemistry data frame with columns `sample_id`,
#'   `site_id`, `datetime`, `analyte`, `value`, `detected`.
#' @param drivers Required driver analyte names. Default `c("pH", "EC",
#'   "NH3-N")`. Samples where any driver is entirely absent are dropped.
#' @param wq_candidates Candidate WQ analytes for PCA. Default
#'   `.WQ_BLOCK_CANDIDATES`.
#' @param metal_analytes Analyte names classified as metals.  Default
#'   `.METAL_ANALYTES`.
#' @param doc_like_analytes Analyte names used for the organics hurdle check
#'   (the "organic carbon present" requirement).  Default `.DOC_LIKE_ANALYTES`.
#' @param min_detect_freq Minimum detection frequency for a WQ analyte to be
#'   included in the PCA block.  Default `0.05`.
#' @param min_samples Minimum training samples after driver filtering.
#' @param min_var_explained Target cumulative variance for PCA axis selection.
#'   Default `0.75`.
#' @param max_pcs Maximum PCA axes to use.  Default `4L`.
#' @param family brms response family.  Must be `"gaussian"` (concentrations
#'   are log-transformed before fitting; residual correlations require
#'   Gaussian family).
#' @param iter,warmup,chains,cores brms MCMC settings.
#' @param save_dir If non-NULL, save the returned model object as a `.qs` file
#'   in this directory using `qs::qsave()`.
#' @param ... Additional arguments passed to `brms::brm()`.
#'
#' @return A named list of class `"imputation_model"`:
#'   - `$pca`: PCA fit + metadata (loadings, training medians, n_pcs, …)
#'   - `$metals`: list with `$fit` (brmsfit), `$analytes`, `$safe_names`
#'   - `$organics`: same structure, or `NULL` if no organics pass prescreen
#'   - `$drivers`, `$hurdle_metals`, `$hurdle_organics`: character vectors
#'   - `$fit_date`, `$n_samples`: metadata
#'   If `save_dir` is supplied, the path to the saved file is returned as
#'   `attr(result, "save_path")`.
#'
#' @seealso [impute_chemistry()]
#' @export
fit_imputation_model <- function(
    df,
    drivers           = c("pH", "EC", "NH3-N"),
    wq_candidates     = NULL,
    metal_analytes    = NULL,
    doc_like_analytes = NULL,
    min_detect_freq   = 0.05,
    min_samples       = 10L,
    min_var_explained = 0.75,
    max_pcs           = 4L,
    family            = "gaussian",
    iter              = 2000,
    warmup            = 1000,
    chains            = 4,
    cores             = parallel::detectCores(),
    save_dir          = NULL,
    ...
) {
  if (is.null(wq_candidates))     wq_candidates     <- .WQ_BLOCK_CANDIDATES
  if (is.null(metal_analytes))    metal_analytes    <- .METAL_ANALYTES
  if (is.null(doc_like_analytes)) doc_like_analytes <- .DOC_LIKE_ANALYTES

  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected"))
  checkmate::assert_character(drivers, min.len = 1L, any.missing = FALSE)
  checkmate::assert_count(min_samples)

  # ── 1. BDL driver handling ─────────────────────────────────────────────────
  # For samples where a driver is BDL, use the stored detection-limit value.
  # These are genuine low-concentration events and must not be silently dropped.
  n_bdl_drivers <- df |>
    dplyr::filter(.data$analyte %in% .env$drivers, !.data$detected) |>
    nrow()
  if (n_bdl_drivers > 0L) {
    cli::cli_inform(c(
      "i" = "{n_bdl_drivers} BDL row{?s} for driver analyte{?s} — using \\
             detection-limit value{?s} as conservative estimate."
    ))
  }
  # `value` already holds the DL for BDL rows; no transformation needed.
  # We just need to NOT filter out BDL driver rows.

  # ── 2. Drop samples missing any driver entirely ────────────────────────────
  samples_with_all_drivers <- df |>
    dplyr::filter(.data$analyte %in% .env$drivers) |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::summarise(n_drivers = dplyr::n_distinct(.data$analyte), .groups = "drop") |>
    dplyr::filter(.data$n_drivers == length(drivers)) |>
    dplyr::pull(.data$sample_id)

  n_dropped <- dplyr::n_distinct(df$sample_id) - length(samples_with_all_drivers)
  if (n_dropped > 0L) {
    cli::cli_inform(c(
      "!" = "{n_dropped} sample{?s} dropped: missing one or more drivers \\
             ({.val {drivers}}) entirely."
    ))
    df <- dplyr::filter(df, .data$sample_id %in% samples_with_all_drivers)
  }

  if (length(samples_with_all_drivers) < min_samples) {
    cli::cli_abort(c(
      "Only {length(samples_with_all_drivers)} sample{?s} remain after driver \\
       filtering — fewer than {.arg min_samples} = {min_samples}.",
      "i" = "Lower {.arg min_samples}, add more data, or choose different drivers."
    ))
  }

  # ── 3. Identify WQ block variables ────────────────────────────────────────
  # Keep WQ candidates that are present in df and meet detection-frequency threshold
  n_samples_total <- dplyr::n_distinct(df$sample_id)
  wq_present <- df |>
    dplyr::filter(.data$analyte %in% .env$wq_candidates) |>
    dplyr::group_by(.data$analyte) |>
    dplyr::summarise(
      detect_freq = dplyr::n_distinct(.data$sample_id) / n_samples_total,
      .groups = "drop"
    ) |>
    dplyr::filter(.data$detect_freq >= min_detect_freq) |>
    dplyr::pull(.data$analyte)

  cli::cli_inform(c(
    "i" = "WQ block: {length(wq_present)} variable{?s} available for PCA: \\
           {.val {sort(wq_present)}}."
  ))

  # ── 4. Identify analyte groups ─────────────────────────────────────────────
  all_analytes <- unique(df$analyte)
  # Exclude WQ block, drivers, and explicitly excluded analytes (microbiological
  # counts, qualitative/physical descriptors — not amenable to log-normal model)
  non_wq_non_driver <- setdiff(
    all_analytes,
    union(union(wq_candidates, drivers), .IMPUTE_EXCLUDED)
  )

  excl_present <- intersect(all_analytes, .IMPUTE_EXCLUDED)
  if (length(excl_present) > 0L) {
    cli::cli_inform(c(
      "i" = "{length(excl_present)} analyte{?s} excluded from imputation \\
             (non-concentration data): {.val {sort(excl_present)}}."
    ))
  }

  metals_in_df   <- intersect(non_wq_non_driver, metal_analytes)
  organics_in_df <- setdiff(non_wq_non_driver, metal_analytes)

  cli::cli_inform(c(
    "i" = "Metals group: {length(metals_in_df)} analyte{?s}: \\
           {.val {sort(metals_in_df)}}.",
    "i" = "Organics group: {length(organics_in_df)} analyte{?s}: \\
           {.val {sort(organics_in_df)}}."
  ))

  if (length(metals_in_df) == 0L && length(organics_in_df) == 0L) {
    cli::cli_warn(
      "No target analytes found outside the driver, WQ block, and excluded sets. \\
       Returning model with no fitted groups (imputation will be a no-op)."
    )
    return(structure(
      list(
        pca             = NULL,
        metals          = NULL,
        organics        = NULL,
        drivers         = drivers,
        wq_candidates   = wq_candidates,
        hurdle_metals   = metal_analytes,
        hurdle_organics = doc_like_analytes,
        fit_date        = Sys.Date(),
        n_samples       = length(samples_with_all_drivers)
      ),
      class = "imputation_model"
    ))
  }

  # ── 5. Fit PCA on WQ block ─────────────────────────────────────────────────
  pca_obj <- .prepare_wq_pca(
    df, wq_vars        = wq_present,
    min_var_explained  = min_var_explained,
    max_pcs            = max_pcs
  )

  cli::cli_inform(c(
    "i" = "WQ PCA: {pca_obj$n_pcs} axis/axes explain \\
           {round(100 * pca_obj$var_explained, 1)}% of WQ variance."
  ))

  # ── 6. Fit metals model ────────────────────────────────────────────────────
  metals_fit <- NULL
  if (length(metals_in_df) >= 1L) {
    cli::cli_inform(c("i" = "Fitting metals model …"))
    metals_fit <- .fit_group_model(
      df          = df,
      target_analytes = metals_in_df,
      drivers     = drivers,
      pca_obj     = pca_obj,
      family      = family,
      iter        = iter,
      warmup      = warmup,
      chains      = chains,
      cores       = cores,
      ...
    )
  }

  # ── 7. Fit organics model ──────────────────────────────────────────────────
  organics_fit <- NULL
  if (length(organics_in_df) >= 1L) {
    cli::cli_inform(c("i" = "Fitting organics model …"))
    organics_fit <- .fit_group_model(
      df          = df,
      target_analytes = organics_in_df,
      drivers     = drivers,
      pca_obj     = pca_obj,
      family      = family,
      iter        = iter,
      warmup      = warmup,
      chains      = chains,
      cores       = cores,
      ...
    )
  }

  # ── 8. Assemble result ────────────────────────────────────────────────────
  result <- structure(
    list(
      pca             = pca_obj,
      metals          = metals_fit,
      organics        = organics_fit,
      drivers         = drivers,
      wq_candidates   = wq_candidates,
      hurdle_metals   = metal_analytes,
      hurdle_organics = doc_like_analytes,
      fit_date        = Sys.Date(),
      n_samples       = length(samples_with_all_drivers)
    ),
    class = "imputation_model"
  )

  # ── 9. Save if requested ───────────────────────────────────────────────────
  if (!is.null(save_dir)) {
    if (!requireNamespace("qs", quietly = TRUE))
      cli::cli_abort("Package {.pkg qs} is required for saving models. \\
                      Install with: {.code install.packages('qs')}")
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    fname     <- sprintf("imputation_model_%s.qs", format(Sys.Date(), "%Y%m%d"))
    save_path <- file.path(save_dir, fname)
    qs::qsave(result, save_path)
    cli::cli_inform(c("v" = "Model saved to {.path {save_path}}"))
    attr(result, "save_path") <- save_path
  }

  result
}

#' @export
print.imputation_model <- function(x, ...) {
  cat(sprintf(
    "<imputation_model>  fitted %s | %d samples | %d WQ PCs (%.0f%% var)\n",
    x$fit_date, x$n_samples,
    x$pca$n_pcs, 100 * x$pca$var_explained
  ))
  if (!is.null(x$metals))
    cat(sprintf("  metals:   %d analytes\n", length(x$metals$analytes)))
  if (!is.null(x$organics))
    cat(sprintf("  organics: %d analytes\n", length(x$organics$analytes)))
  invisible(x)
}


# ── impute_chemistry() ────────────────────────────────────────────────────────

#' Impute missing and BDL chemistry using a fitted imputation model
#'
#' Applies the models fitted by [fit_imputation_model()] to `df`, returning
#' posterior mean estimates for missing and below-detection-limit (BDL)
#' observations in the metals and organics groups.
#'
#' **Hurdles**
#'
#' Imputed values are only returned for samples that meet the relevant hurdle:
#' - *Metals*: at least one metal analyte present (detected or BDL) at the
#'   sample.  Samples with no metals recorded are not imputed — a leachate
#'   metal pulse may simply not have arrived at this location yet.
#' - *Organics*: at least one of {DOC, TOC, BOD, COD, cBOD} present at the
#'   sample.
#'
#' Samples failing a hurdle pass through unchanged (non-imputed values are
#' preserved; BDL values remain flagged as BDL).
#'
#' @param df Long-format chemistry data frame (same schema as used for fitting).
#' @param model Fitted model from [fit_imputation_model()].
#' @param metal_hurdle Logical.  Apply metal-presence hurdle?  Default `TRUE`.
#' @param organic_hurdle Logical.  Apply DOC-like-presence hurdle?  Default
#'   `TRUE`.
#' @param bdl_cap Logical.  Cap imputed BDL values at the original detection
#'   limit?  Default `TRUE`.
#' @param return `"point"` (default) for posterior mean per cell; `"draws"` for
#'   one row per (sample × analyte × draw).
#'
#' @return `df` with BDL and missing cells in the metals/organics groups
#'   replaced by posterior mean estimates, plus columns:
#'   - `imputed` (logical) — `TRUE` for filled cells
#'   - `imputed_kind` — `"observed"`, `"censored_left"`, or `"missing"`
#'
#' @seealso [fit_imputation_model()]
#' @export
impute_chemistry <- function(
    df,
    model,
    metal_hurdle   = TRUE,
    organic_hurdle = TRUE,
    bdl_cap        = TRUE,
    return         = c("point", "draws")
) {
  return <- match.arg(return)
  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected"))
  if (!inherits(model, "imputation_model"))
    cli::cli_abort("{.arg model} must be an object returned by {.fn fit_imputation_model}.")

  # ── Compute WQ PC scores for new data ─────────────────────────────────────
  pca_scores <- .compute_pca_scores(df, model$pca)

  # ── Collect BDL detection limits ──────────────────────────────────────────
  all_targets <- c(
    if (!is.null(model$metals))   model$metals$analytes   else character(0),
    if (!is.null(model$organics)) model$organics$analytes else character(0)
  )
  dl_tbl <- df |>
    dplyr::filter(.data$analyte %in% .env$all_targets, !.data$detected) |>
    dplyr::select("sample_id", "analyte", detection_limit = "value")

  # ── Impute each group ──────────────────────────────────────────────────────
  result <- df  # start with original; overlay imputed values below

  if (!is.null(model$metals)) {
    eligible <- if (metal_hurdle) {
      df |>
        dplyr::filter(.data$analyte %in% model$hurdle_metals) |>
        dplyr::pull(.data$sample_id) |>
        unique()
    } else {
      unique(df$sample_id)
    }
    n_skip <- dplyr::n_distinct(df$sample_id) - length(eligible)
    if (n_skip > 0L)
      cli::cli_inform(c("i" = "Metals hurdle: skipping {n_skip} sample{?s} \\
                                (no metals present)."))

    result <- .predict_and_merge(
      df            = result,
      group         = model$metals,
      pca_scores    = pca_scores,
      drivers       = model$drivers,
      eligible_ids  = eligible,
      return        = return
    )
  }

  if (!is.null(model$organics)) {
    eligible_org <- if (organic_hurdle) {
      df |>
        dplyr::filter(.data$analyte %in% model$hurdle_organics) |>
        dplyr::pull(.data$sample_id) |>
        unique()
    } else {
      unique(df$sample_id)
    }
    n_skip_org <- dplyr::n_distinct(df$sample_id) - length(eligible_org)
    if (n_skip_org > 0L)
      cli::cli_inform(c("i" = "Organics hurdle: skipping {n_skip_org} sample{?s} \\
                                (no DOC-like variable present)."))

    result <- .predict_and_merge(
      df            = result,
      group         = model$organics,
      pca_scores    = pca_scores,
      drivers       = model$drivers,
      eligible_ids  = eligible_org,
      return        = return
    )
  }

  # Tag non-target rows (drivers, WQ vars, etc.) that were never imputed
  if (!"imputed" %in% names(result)) {
    result <- dplyr::mutate(result,
      imputed      = FALSE,
      imputed_kind = "observed"
    )
  } else {
    result <- dplyr::mutate(result,
      imputed      = dplyr::coalesce(.data$imputed, FALSE),
      imputed_kind = dplyr::coalesce(.data$imputed_kind, "observed")
    )
  }

  .check_bdl_imputed(result, dl_tbl, bdl_cap)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Fit PCA on the WQ block for training data
#' @keywords internal
.prepare_wq_pca <- function(df, wq_vars, min_var_explained = 0.75, max_pcs = 4L) {
  # Pivot WQ to wide (one row per sample); missing cells → NA
  wq_wide <- df |>
    dplyr::filter(.data$analyte %in% .env$wq_vars) |>
    dplyr::select("sample_id", "analyte", "value") |>
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      .by   = c("sample_id", "analyte")
    ) |>
    tidyr::pivot_wider(names_from = "analyte", values_from = "value")

  # Ensure every sample_id is present (even those with no WQ vars)
  all_samples <- tibble::tibble(sample_id = unique(df$sample_id))
  wq_wide     <- dplyr::left_join(all_samples, wq_wide, by = "sample_id")

  wq_matrix   <- as.matrix(dplyr::select(wq_wide, -"sample_id"))

  # Training medians (used to fill missing cells in new data)
  train_medians <- apply(wq_matrix, 2, stats::median, na.rm = TRUE)

  # Median-impute for PCA
  for (j in seq_len(ncol(wq_matrix))) {
    na_idx <- is.na(wq_matrix[, j])
    if (any(na_idx)) {
      wq_matrix[na_idx, j] <- if (!is.na(train_medians[j])) train_medians[j] else 0
    }
  }

  # Remove zero-variance columns (PCA would divide by 0)
  col_sds    <- apply(wq_matrix, 2, stats::sd)
  keep_cols  <- col_sds > 0
  if (!any(keep_cols))
    cli::cli_abort("All WQ block variables have zero variance — cannot fit PCA.")
  if (!all(keep_cols)) {
    dropped <- colnames(wq_matrix)[!keep_cols]
    cli::cli_inform(c("!" = "WQ PCA: dropping zero-variance variable{?s}: \\
                             {.val {dropped}}"))
    wq_matrix <- wq_matrix[, keep_cols, drop = FALSE]
  }

  pca_fit <- stats::prcomp(wq_matrix, center = TRUE, scale. = TRUE)
  cum_var <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)

  # Determine number of PCs
  n_needed <- which(cum_var >= min_var_explained)[1L]
  if (is.na(n_needed) || n_needed > max_pcs) {
    n_pcs <- min(max_pcs, length(cum_var))
    cli::cli_warn(c(
      "!" = "WQ PCA: first {n_pcs} axis/axes explain only \\
             {round(100 * cum_var[n_pcs], 1)}% of variance \\
             (target: {min_var_explained * 100}%).",
      "i" = "Consider adding more WQ variables or loosening \\
             {.arg min_var_explained}."
    ))
  } else {
    n_pcs <- max(2L, n_needed)
  }
  n_pcs <- min(n_pcs, ncol(wq_matrix))

  pc_scores <- tibble::as_tibble(pca_fit$x[, seq_len(n_pcs), drop = FALSE]) |>
    stats::setNames(paste0("WQ_PC", seq_len(n_pcs))) |>
    dplyr::mutate(sample_id = wq_wide$sample_id)

  list(
    fit           = pca_fit,
    medians       = train_medians,         # all WQ vars (before zero-var removal)
    active_vars   = colnames(wq_matrix),   # after zero-var removal
    n_pcs         = n_pcs,
    var_explained = cum_var[n_pcs],
    pc_scores     = pc_scores
  )
}

#' Project new data onto stored PCA
#' @keywords internal
.compute_pca_scores <- function(df, pca_obj) {
  wq_vars <- names(pca_obj$medians)

  wq_wide <- df |>
    dplyr::filter(.data$analyte %in% .env$wq_vars) |>
    dplyr::select("sample_id", "analyte", "value") |>
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      .by   = c("sample_id", "analyte")
    ) |>
    tidyr::pivot_wider(names_from = "analyte", values_from = "value")

  all_samples <- tibble::tibble(sample_id = unique(df$sample_id))
  wq_wide     <- dplyr::left_join(all_samples, wq_wide, by = "sample_id")

  # Ensure columns match active_vars (add NA for unseen WQ vars)
  for (v in pca_obj$active_vars) {
    if (!v %in% names(wq_wide)) wq_wide[[v]] <- NA_real_
  }

  wq_mat <- as.matrix(dplyr::select(wq_wide, dplyr::all_of(pca_obj$active_vars)))

  # Median-impute using training medians
  for (j in seq_len(ncol(wq_mat))) {
    col_nm  <- colnames(wq_mat)[j]
    med_val <- pca_obj$medians[col_nm]
    na_idx  <- is.na(wq_mat[, j])
    if (any(na_idx)) {
      wq_mat[na_idx, j] <- if (!is.na(med_val)) med_val else 0
    }
  }

  # Project onto training PCA space
  wq_scaled  <- scale(wq_mat,
                      center = pca_obj$fit$center,
                      scale  = pca_obj$fit$scale)
  scores_mat <- wq_scaled %*% pca_obj$fit$rotation[, seq_len(pca_obj$n_pcs), drop = FALSE]

  tibble::as_tibble(scores_mat) |>
    stats::setNames(paste0("WQ_PC", seq_len(pca_obj$n_pcs))) |>
    dplyr::mutate(sample_id = wq_wide$sample_id)
}

#' Fit a single brms group model (metals or organics)
#' @keywords internal
.fit_group_model <- function(df, target_analytes, drivers, pca_obj,
                              family, iter, warmup, chains, cores, ...) {
  eps_log <- 1e-9

  safe_analytes   <- stats::setNames(make.names(target_analytes), target_analytes)
  # safe_analytes: names = safe R names, values = original names
  # (reverse of the old convention to make lookup unambiguous)

  driver_col_names <- paste0(".drv_", make.names(drivers))

  # ── Wide drivers ───────────────────────────────────────────────────────────
  driver_wide <- df |>
    dplyr::filter(.data$analyte %in% .env$drivers) |>
    dplyr::select("sample_id", "analyte", "value") |>
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      .by   = c("sample_id", "analyte")
    ) |>
    tidyr::pivot_wider(
      names_from   = "analyte",
      values_from  = "value",
      names_prefix = ".drv_",
      names_repair = "universal"
    ) |>
    # Standardise column names: make.names applied after prefix
    dplyr::rename_with(~ paste0(".drv_", make.names(sub("^\\.drv_", "", .x))),
                       dplyr::starts_with(".drv_"))

  # ── Wide targets (log scale; NA for BDL and missing) ─────────────────────
  target_wide <- df |>
    dplyr::filter(.data$analyte %in% .env$target_analytes) |>
    dplyr::select("sample_id", "analyte", "value", "detected") |>
    dplyr::mutate(
      log_value = dplyr::if_else(
        .data$detected,
        log(pmax(.data$value, eps_log)),
        NA_real_
      )
    ) |>
    dplyr::summarise(
      log_value = mean(.data$log_value, na.rm = TRUE),
      .by        = c("sample_id", "analyte")
    ) |>
    tidyr::pivot_wider(names_from = "analyte", values_from = "log_value") |>
    # Rename to safe R column names using the safe_analytes map
    dplyr::rename(dplyr::any_of(
      stats::setNames(names(safe_analytes), unname(safe_analytes))
      # setNames(old_names, new_names): names = new safe names, values = old original names
    ))

  # ── Join everything ────────────────────────────────────────────────────────
  pc_cols    <- paste0("WQ_PC", seq_len(pca_obj$n_pcs))
  wide_df    <- driver_wide |>
    dplyr::left_join(
      dplyr::select(pca_obj$pc_scores, "sample_id", dplyr::all_of(pc_cols)),
      by = "sample_id"
    ) |>
    dplyr::left_join(target_wide, by = "sample_id")

  # ── brms formula ──────────────────────────────────────────────────────────
  rhs <- paste(
    c(paste0("s(", driver_col_names, ")"),
      paste0("s(", pc_cols, ")")),
    collapse = " + "
  )

  bf_list <- purrr::map(unname(safe_analytes), function(safe_nm) {
    brms::bf(stats::as.formula(paste0(safe_nm, " | mi() ~ ", rhs)))
  })
  brms_formula <- do.call(brms::mvbf, c(bf_list, list(rescor = TRUE)))

  cli::cli_inform(c(
    "i" = "brms: {length(target_analytes)} analyte{?s} × \\
           {nrow(wide_df)} sample{?s}. This may take several minutes."
  ))

  fit <- tryCatch(
    brms::brm(
      formula = brms_formula,
      data    = wide_df,
      family  = family,
      iter    = iter,
      warmup  = warmup,
      chains  = chains,
      cores   = cores,
      ...
    ),
    error = function(e) {
      cli::cli_abort(c(
        "brms model fitting failed.",
        "x" = "{conditionMessage(e)}",
        "i" = "If this is a Stan compilation error, check your Stan toolchain \\
               ({.url https://paul-buerkner.github.io/brms/})."
      ))
    }
  )

  list(
    fit          = fit,
    analytes     = target_analytes,        # original names
    safe_names   = safe_analytes,          # names=safe, values=original
    driver_cols  = driver_col_names,
    pc_cols      = pc_cols,
    wide_sample_ids = wide_df$sample_id
  )
}

#' Predict and merge imputed values for one analyte group
#' @keywords internal
.predict_and_merge <- function(df, group, pca_scores, drivers,
                                eligible_ids, return) {
  eps_log <- 1e-9

  target_analytes  <- group$analytes
  safe_analytes    <- group$safe_names    # names=safe, values=original
  driver_col_names <- group$driver_cols
  pc_cols          <- group$pc_cols

  # ── Build wide prediction df for eligible samples ─────────────────────────
  df_eligible <- dplyr::filter(df, .data$sample_id %in% .env$eligible_ids)
  if (nrow(df_eligible) == 0L) return(df)

  # Drivers wide
  driver_wide <- df_eligible |>
    dplyr::filter(.data$analyte %in% .env$drivers) |>
    dplyr::select("sample_id", "analyte", "value") |>
    dplyr::summarise(
      value = mean(.data$value, na.rm = TRUE),
      .by   = c("sample_id", "analyte")
    ) |>
    tidyr::pivot_wider(
      names_from   = "analyte",
      values_from  = "value",
      names_prefix = ".drv_",
      names_repair = "universal"
    ) |>
    dplyr::rename_with(~ paste0(".drv_", make.names(sub("^\\.drv_", "", .x))),
                       dplyr::starts_with(".drv_"))

  # Targets wide (log for detected, NA for BDL/missing → to be imputed)
  target_wide <- df_eligible |>
    dplyr::filter(.data$analyte %in% .env$target_analytes) |>
    dplyr::select("sample_id", "analyte", "value", "detected") |>
    dplyr::mutate(
      log_value = dplyr::if_else(
        .data$detected,
        log(pmax(.data$value, eps_log)),
        NA_real_
      )
    ) |>
    dplyr::summarise(
      log_value = mean(.data$log_value, na.rm = TRUE),
      .by        = c("sample_id", "analyte")
    ) |>
    tidyr::pivot_wider(names_from = "analyte", values_from = "log_value") |>
    dplyr::rename(dplyr::any_of(
      stats::setNames(names(safe_analytes), unname(safe_analytes))
    ))

  # PC scores
  pc_wide <- dplyr::filter(pca_scores, .data$sample_id %in% .env$eligible_ids) |>
    dplyr::select("sample_id", dplyr::all_of(pc_cols))

  # Ensure all training columns are present (add NA if analyte not in new data)
  all_target_safe <- names(safe_analytes)  # safe R names
  for (s in all_target_safe) {
    if (!s %in% names(target_wide)) target_wide[[s]] <- NA_real_
  }

  wide_new <- driver_wide |>
    dplyr::left_join(pc_wide,    by = "sample_id") |>
    dplyr::left_join(target_wide, by = "sample_id")

  # ── Posterior predictions ─────────────────────────────────────────────────
  if (return == "point") {
    epred <- tryCatch(
      brms::posterior_epred(group$fit, newdata = wide_new,
                            allow_new_levels = TRUE),
      error = function(e) cli::cli_abort(c(
        "brms::posterior_epred() failed during imputation.",
        "x" = "{conditionMessage(e)}"
      ))
    )
    pm_long <- .reshape_posterior_means(epred, wide_new$sample_id, safe_analytes)
  } else {
    post_draws <- brms::posterior_predict(group$fit, newdata = wide_new,
                                          allow_new_levels = TRUE)
    pm_long <- .reshape_posterior_draws(post_draws, wide_new$sample_id, safe_analytes)
  }

  # ── Tag imputation kind for each (sample, analyte) ─────────────────────────
  impute_kind <- df_eligible |>
    dplyr::filter(.data$analyte %in% .env$target_analytes) |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      sample_id,
      analyte,
      .imputed_kind = dplyr::if_else(.data$detected, "observed", "censored_left")
    )

  # Samples missing the analyte entirely → "missing"
  all_combos <- tidyr::expand_grid(
    sample_id = unique(df_eligible$sample_id),
    analyte   = target_analytes
  )
  impute_kind <- dplyr::left_join(all_combos, impute_kind,
                                  by = c("sample_id", "analyte")) |>
    dplyr::mutate(
      .imputed_kind = dplyr::if_else(is.na(.data$.imputed_kind),
                                     "missing", .data$.imputed_kind)
    )

  # ── Merge posterior values back into df ───────────────────────────────────
  # For eligible samples: replace BDL / missing target values with posterior
  # means; add imputed/imputed_kind columns.  Non-target rows are unchanged.

  val_col <- if (return == "point") ".post_mean" else ".post_value"

  if (!"imputed" %in% names(df)) {
    df <- dplyr::mutate(df, imputed = FALSE, imputed_kind = "observed")
  }

  target_rows_eligible <- df |>
    dplyr::filter(
      .data$sample_id %in% .env$eligible_ids,
      .data$analyte   %in% .env$target_analytes
    ) |>
    dplyr::left_join(impute_kind, by = c("sample_id", "analyte")) |>
    dplyr::mutate(
      imputed      = .data$.imputed_kind != "observed",
      imputed_kind = .data$.imputed_kind
    ) |>
    dplyr::select(-".imputed_kind")

  # Overlay posterior means onto imputed rows
  imputed_rows   <- dplyr::filter(target_rows_eligible, .data$imputed)
  observed_rows  <- dplyr::filter(target_rows_eligible, !.data$imputed)

  join_cols <- if (return == "point") {
    c("sample_id", "analyte", val_col)
  } else {
    c("sample_id", "analyte", "draw_id", val_col)
  }

  imputed_filled <- dplyr::left_join(
    imputed_rows,
    dplyr::select(pm_long, dplyr::all_of(join_cols)),
    by = c("sample_id", "analyte")
  ) |>
    dplyr::mutate(
      value    = dplyr::coalesce(.data[[val_col]], .data$value),
      detected = TRUE
    ) |>
    dplyr::select(-dplyr::all_of(val_col))

  # Rows for non-eligible samples (failed hurdle): keep unchanged
  non_eligible_target_rows <- df |>
    dplyr::filter(
      !(.data$sample_id %in% .env$eligible_ids),
      .data$analyte %in% .env$target_analytes
    )

  # Non-target rows
  non_target_rows <- df |>
    dplyr::filter(!(.data$analyte %in% .env$target_analytes))

  dplyr::bind_rows(
    non_target_rows,
    observed_rows,
    imputed_filled,
    non_eligible_target_rows
  ) |>
    dplyr::arrange(.data$sample_id, .data$analyte)
}


# ── Posterior reshape helpers (unchanged from original) ──────────────────────

#' @keywords internal
.reshape_posterior_means <- function(epred_draws, sample_ids, safe_analytes) {
  if (is.matrix(epred_draws)) {
    arr_list <- stats::setNames(list(epred_draws), unname(safe_analytes)[1L])
  } else {
    resp_nms <- dimnames(epred_draws)[[3L]]
    arr_list <- stats::setNames(
      lapply(seq_along(resp_nms), function(i) epred_draws[, , i]),
      resp_nms
    )
  }

  purrr::map2_dfr(arr_list, names(arr_list), function(mat, safe_nm) {
    orig_nm <- names(safe_analytes)[safe_analytes == safe_nm]
    if (length(orig_nm) == 0L) orig_nm <- safe_nm
    tibble::tibble(
      sample_id  = sample_ids,
      analyte    = orig_nm,
      .post_mean = exp(colMeans(mat))
    )
  })
}

#' @keywords internal
.reshape_posterior_draws <- function(post_draws, sample_ids, safe_analytes) {
  n_draws  <- dim(post_draws)[1L]
  resp_nms <- dimnames(post_draws)[[3L]]

  purrr::map2_dfr(seq_along(resp_nms), resp_nms, function(ri, safe_nm) {
    orig_nm <- names(safe_analytes)[safe_analytes == safe_nm]
    if (length(orig_nm) == 0L) orig_nm <- safe_nm
    mat <- post_draws[, , ri]
    tibble::tibble(
      sample_id   = rep(sample_ids, each = n_draws),
      analyte     = orig_nm,
      draw_id     = rep(seq_len(n_draws), times = length(sample_ids)),
      .post_value = exp(as.vector(t(mat)))
    )
  })
}


# ── BDL cap check (unchanged) ─────────────────────────────────────────────────

#' @keywords internal
.check_bdl_imputed <- function(result, dl_tbl, cap = TRUE) {
  if (nrow(dl_tbl) == 0L) return(result)

  bdl_rows <- dplyr::filter(result, .data$imputed_kind == "censored_left") |>
    dplyr::left_join(dl_tbl, by = c("sample_id", "analyte"))

  if (nrow(bdl_rows) == 0L || !"detection_limit" %in% names(bdl_rows))
    return(result)

  exceedances <- dplyr::filter(
    bdl_rows,
    !is.na(.data$detection_limit),
    .data$value > .data$detection_limit
  )

  if (nrow(exceedances) == 0L) return(result)

  n_ex        <- nrow(exceedances)
  analytes_ex <- unique(exceedances$analyte)

  cli::cli_warn(c(
    "!" = "{n_ex} imputed BDL value{?s} exceed the original detection limit.",
    "i" = "Affected analyte{?s}: {.val {analytes_ex}}.",
    "i" = "Using {.code mi()} instead of {.code cens('left')} (required for \\
           {.code rescor = TRUE}) means the left-censor constraint is not \\
           enforced during MCMC.",
    if (cap) "i" = "Values capped at DL ({.arg bdl_cap = TRUE})."
    else     "i" = "Values NOT capped ({.arg bdl_cap = FALSE})."
  ))

  if (!cap) return(result)

  exceedance_keys <- dplyr::select(exceedances, "sample_id", "analyte",
                                    cap_value = "detection_limit")
  dplyr::left_join(result, exceedance_keys, by = c("sample_id", "analyte")) |>
    dplyr::mutate(
      value = dplyr::if_else(!is.na(.data$cap_value), .data$cap_value, .data$value)
    ) |>
    dplyr::select(-"cap_value")
}
