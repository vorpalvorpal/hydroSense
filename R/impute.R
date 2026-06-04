# ── Analyte group constants ───────────────────────────────────────────────────

#' WQ analytes used in the PCA pre-processing step
#'
#' All candidates; the actual variables used at any given site are the
#' intersection of this list with analytes that pass prescreen in the training
#' data.  ORP and DO are included because they are valuable at sites where they
#' are measured, but at most sites they will be filtered out by prescreen.
#'
#' Note: the sulfate entry encodes its superscript charge (the SO4
#' two-minus symbol) with Unicode escape sequences in the source string, so
#' that matching works regardless of how the source file is parsed.
#' @keywords internal
.WQ_BLOCK_CANDIDATES <- c(
  # Field parameters / redox
  "temperature", "ORP", "DO",
  # Major cations
  "Ca", "Mg", "Na", "K",
  # Major anions & alkalinity
  "Cl", "SO4\u00b2\u207b", "Alkalinity-total-CaCO3", "F",
  # Hardness: total water hardness in mg/L as CaCO3 equivalents.
  # Callers must convert their data to this convention before passing in.
  "hardness",
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

#' Co-analytes required by ANZECC/ANZG metal normalisation formulas
#'
#' These are imputed separately (after metals/organics imputation) via
#' [impute_coanalytes()] so that [add_amspaf()] has values to normalise
#' against.  pH and EC are excluded — they are always present (required vars).
#' @keywords internal
.COANALYTE_TARGETS <- c("DOC", "Ca", "Mg", "hardness")

#' PCA variables that must NOT be log-transformed before the chemistry PCA
#'
#' Every other PCA variable is concentration-like — strictly positive and
#' strongly right-skewed, spanning orders of magnitude — so it is
#' `log10`-transformed before centring/scaling.  Without that, the PCA is
#' dominated by a handful of high-magnitude major ions (e.g. Cl, SO4, TDS) and
#' the leading axes mostly track absolute ionic strength rather than the
#' multiplicative covariance structure that drives metal/organic behaviour.
#'
#' The exclusions are the variables for which a log is meaningless or undefined:
#'   - `pH` — already a logarithmic scale (−log10 of H+ activity).
#'   - `temperature` — interval scale (°C); zero/negative values are valid.
#'   - `ORP` — redox potential (mV); routinely negative.
#'   - `DO` — dissolved oxygen (mg/L); legitimately ~0 in anoxic leachate
#'     plumes and only spans a narrow, near-linear range.
#' @keywords internal
.PCA_NO_LOG_VARS <- c("pH", "temperature", "ORP", "DO")


# ── brms availability guard ───────────────────────────────────────────────────

#' Stop with a friendly, actionable message if brms is not installed
#'
#' The Bayesian imputation step ([fit_imputation_model()] /
#' [impute_chemistry()]) is the only part of the package that needs
#' \pkg{brms}, so brms is an optional ("Suggests") dependency rather than a
#' hard requirement.  This keeps the package quick to install for users who
#' only need the LMF or AmsPAF tools.  When someone actually calls an
#' imputation function without brms installed, this guard explains — in plain
#' language — what to install and why.
#' @keywords internal
.require_brms <- function() {
  if (requireNamespace("brms", quietly = TRUE)) {
    return(invisible(TRUE))
  }
  cli::cli_abort(c(
    "The chemistry imputation step needs the {.pkg brms} package, which \\
     isn't installed yet.",
    "i" = "{.pkg brms} fits the Bayesian model that fills in missing and \\
           below-detection-limit results. It's optional, so it isn't \\
           installed automatically \u2014 only this imputation step uses it.",
    " " = "",
    "*" = "To install it, run this once at the R console:",
    " " = "{.code install.packages(\"brms\")}",
    " " = "",
    "i" = "{.pkg brms} also needs a working Stan engine (a C++ compiler). If \\
           the line above isn't enough, follow the short setup guide at \\
           {.url https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}.",
    "i" = "Once {.pkg brms} is installed, re-run this function \u2014 no other \\
           changes are needed."
  ))
}


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
#' s(PC1) + s(PC2) + ... + s(PCk)
#' ```
#' where `PC*` are the leading principal components of the unified chemistry
#' PCA (see *Chemistry PCA* below).  All target analytes (metals or organics)
#' are modelled jointly with `rescor = TRUE`, so observed co-analytes at a
#' given sample condition the posterior of the missing ones through the residual
#' correlation matrix.
#'
#' **Why `rescor = TRUE`** — the PCA captures the *instantaneous* chemical
#' covariance structure (what is measured together at a single moment), but
#' it cannot capture the *temporal-lag* covariance characteristic of
#' AMD/leachate-impacted aquifers, where conservative tracers move ahead of
#' redox-controlled metals.  At a post-pulse sample the PCA scores have
#' returned toward baseline but Cu/Pb/Zn/Mn remain elevated together — that
#' co-elevation is pure residual correlation with no predictor signal driving
#' it.  `rescor = TRUE` is the right machinery for this and is what makes
#' multivariate imputation borrow strength across analytes.
#'
#' **Costs of `rescor = TRUE`** — brms cannot combine `set_rescor(TRUE)` with
#' `cens("left")`, so this implementation uses `mi()` for BDL values and
#' applies a post-hoc cap (see [impute_chemistry()]).  The cap clips imputed
#' BDL cells to the original detection limit when the model predicts above
#' DL.  For sites where the chemistry context legitimately suggests high
#' concentrations the cap can fire frequently; results in that regime should
#' be inspected.  Three alternative configurations are worth benchmarking on
#' real hold-out data if predictive performance becomes a concern:
#'   - `rescor = TRUE` + `mi()` (current; expected to win on plume-affected
#'     groundwater because cross-analyte residual coupling captures plume
#'     dynamics that the predictor PCA misses).
#'   - `rescor = FALSE` + `cens("left")` (statistically clean for BDL; loses
#'     cross-analyte residual coupling).
#'   - `rescor = FALSE` + `cens("left")` + shared `(1 | sample_id)` (proper
#'     BDL handling with rank-1 latent-factor coupling across analytes).
#' Benchmark methodology: mask 10% of detected cells, fit each configuration,
#' compare hold-out RMSE / coverage.
#'
#' **Chemistry PCA**
#'
#' All `pca_vars` — major ions, pH, EC, NH3-N, DOC, nutrients, redox
#' indicators — that are present in `df` and pass a detection-frequency check
#' are submitted to `nipals::nipals()`, which handles within-sample missing
#' cells natively without prior imputation.  Using a single unified PCA (rather
#' than separate driver + WQ-block sets) eliminates predictor collinearity and
#' ensures normalisation co-analytes (DOC, Ca, Mg) influence the imputed metal
#' concentrations.  Principal components are added until cumulative variance
#' explained reaches `min_var_explained` or `max_pcs` is reached.  A minimum
#' of two PCs is always used.
#'
#' **Hurdles (applied at prediction time by `impute_chemistry()`)**
#'
#' - *Metals*: a sample is only imputed if at least one metal analyte is
#'   present (detected or BDL) in `df` for that sample.
#' - *Organics*: a sample is only imputed if at least one of
#'   DOC, TOC, BOD, COD or cBOD is present.
#'
#' **BDL required variables**
#'
#' When a `required_vars` analyte (pH or EC) is below the detection limit for a
#' sample, the stored detection-limit value is used as-is (conservative upper
#' bound).  A message is issued but the sample is retained.
#'
#' @param df Long-format chemistry data frame with columns `sample_id`,
#'   `site_id`, `datetime`, `analyte`, `value`, `detected`.
#' @param pca_vars Analyte names to include in the unified chemistry PCA (used
#'   as predictors for the brms model via PC scores).  Default: `c("pH", "EC",
#'   "NH3-N")` plus all `.WQ_BLOCK_CANDIDATES`.  Normalisation co-analytes
#'   (DOC, Ca, Mg) are included in the default set.
#' @param required_vars Analyte names that must be present in a sample for it
#'   to be retained in training and prediction.  Default `c("pH", "EC")`.
#'   Samples missing any of these are dropped entirely.
#' @param metal_analytes Analyte names classified as metals.  Default
#'   `.METAL_ANALYTES`.
#' @param doc_like_analytes Analyte names used for the organics hurdle check
#'   (the "organic carbon present" requirement).  Default `.DOC_LIKE_ANALYTES`.
#' @param min_target_detect_freq Minimum detection frequency (fraction of
#'   samples in which the analyte is *detected*) for a metal/organic to be
#'   included as an imputation target. Targets below this are dropped (they have
#'   too few detections to model and would otherwise inflate the model on
#'   near-all-BDL panels). Default `0.05`.
#' @param min_detect_freq Minimum detection frequency for a PCA variable to be
#'   retained.  Default `0.05`.  Required vars are always retained regardless.
#' @param min_samples Minimum training samples after required-var filtering.
#' @param min_var_explained Target cumulative variance for PCA axis selection.
#'   Default `0.75`.
#' @param max_pcs Maximum PCA axes to use.  Default `6L`.
#' @param family brms response family.  Must be `"gaussian"` (concentrations
#'   are log-transformed before fitting; residual correlations require
#'   Gaussian family).
#' @param impute_method How below-detection (BDL) values and cross-analyte
#'   coupling are handled. One of:
#'   \describe{
#'     \item{`"rescor_mi"`}{(default) Residual correlation across analytes
#'       (`rescor = TRUE`) with BDL/missing treated as imputable (`mi()`); the
#'       imputed BDL cells are capped at the detection limit post-hoc by
#'       [impute_chemistry()] (brms cannot combine `rescor` with `cens()`).}
#'     \item{`"cens"`}{Proper left-censoring of BDL at the detection limit
#'       (`cens("left")`), no residual correlation -- clean BDL handling but no
#'       cross-analyte coupling.}
#'     \item{`"cens_factor"`}{As `"cens"` plus a shared per-sample latent factor
#'       (`(1 | sample_id)` correlated across analytes), which re-introduces
#'       cross-analyte coupling while keeping proper censoring.}
#'   }
#'   See `vignette("imputation")` and the package benchmark for guidance on
#'   which to prefer.
#' @param iter,warmup,chains,cores brms MCMC settings.
#' @param save_dir If non-NULL, save the returned model object as a `.qs` file
#'   in this directory using `qs2::qs_save()`.
#' @param ... Additional arguments passed to `brms::brm()`.
#'
#' @return A named list of class `"imputation_model"`:
#'   - `$pca`: PCA fit + metadata (loadings, training medians, n_pcs, …)
#'   - `$metals`: list with `$fit` (brmsfit), `$analytes`, `$safe_names`
#'   - `$organics`: same structure, or `NULL` if no organics pass prescreen
#'   - `$required_vars`, `$pca_vars`, `$hurdle_metals`, `$hurdle_organics`
#'   - `$fit_date`, `$n_samples`: metadata
#'   If `save_dir` is supplied, the path to the saved file is returned as
#'   `attr(result, "save_path")`.
#'
#' @seealso [impute_chemistry()]
#' @examples
#' \dontrun{
#' # Requires a Stan toolchain (brms). Fit once, then reuse for imputation.
#' model <- fit_imputation_model(monitoring_long)
#' draws <- impute_chemistry(monitoring_long, model, return = "draws")
#' }
#' @export
fit_imputation_model <- function(
    df,
    pca_vars          = NULL,           # default built in body
    required_vars     = c("pH", "EC"),
    metal_analytes    = NULL,
    doc_like_analytes = NULL,
    min_detect_freq   = 0.05,
    min_target_detect_freq = 0.05,
    min_samples       = 10L,
    min_var_explained = 0.75,
    max_pcs           = 6L,
    family            = "gaussian",
    impute_method     = c("rescor_mi", "cens", "cens_factor"),
    iter              = 2000,
    warmup            = 1000,
    chains            = 4,
    cores             = parallel::detectCores(),
    save_dir          = NULL,
    ...
) {
  .require_brms()
  impute_method <- match.arg(impute_method)

  if (is.null(pca_vars))          pca_vars          <- c("pH", "EC", "NH3-N",
                                                          .WQ_BLOCK_CANDIDATES)
  if (is.null(metal_analytes))    metal_analytes    <- .METAL_ANALYTES
  if (is.null(doc_like_analytes)) doc_like_analytes <- .DOC_LIKE_ANALYTES

  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected"))
  checkmate::assert_character(required_vars, min.len = 1L, any.missing = FALSE)
  checkmate::assert_count(min_samples)

  # ── 1. BDL required-variable handling ────────────────────────────────────
  # For required vars (pH, EC) where a value is BDL, use the stored DL value.
  # These are genuine low-level conditions; the sample is retained.
  n_bdl_req <- df |>
    dplyr::filter(.data$analyte %in% .env$required_vars, !.data$detected) |>
    nrow()
  if (n_bdl_req > 0L) {
    cli::cli_inform(c(
      "i" = "{n_bdl_req} BDL row{?s} for required variable{?s} \u2014 using \\
             detection-limit value{?s} as conservative estimate."
    ))
  }
  # `value` already holds the DL for BDL rows; no transformation needed.

  # ── 2. Drop samples missing any required variable entirely ────────────────
  samples_with_required <- df |>
    dplyr::filter(.data$analyte %in% .env$required_vars) |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::summarise(n_req = dplyr::n_distinct(.data$analyte), .groups = "drop") |>
    dplyr::filter(.data$n_req == length(required_vars)) |>
    dplyr::pull(.data$sample_id)

  n_dropped <- dplyr::n_distinct(df$sample_id) - length(samples_with_required)
  if (n_dropped > 0L) {
    cli::cli_inform(c(
      "!" = "{n_dropped} sample{?s} dropped: missing one or more required \\
             variable{?s} ({.val {required_vars}}) entirely."
    ))
    df <- dplyr::filter(df, .data$sample_id %in% samples_with_required)
  }

  if (length(samples_with_required) < min_samples) {
    cli::cli_abort(c(
      "Only {length(samples_with_required)} sample{?s} remain after \\
       required-var filtering \u2014 fewer than {.arg min_samples} = {min_samples}.",
      "i" = "Lower {.arg min_samples}, add more data, or choose different \\
             {.arg required_vars}."
    ))
  }

  # ── 3. Filter pca_vars by presence frequency ──────────────────────────────
  # Required vars always pass regardless of presence frequency.
  # NB: this is a PRESENCE frequency (fraction of samples that have a row for
  # the analyte at all), not a detection frequency — a PCA variable is useful
  # as a predictor whether or not it was above the detection limit, so BDL
  # rows count towards retention here. (Contrast prescreen_analytes(), which
  # screens toxicants on true detection frequency.)
  n_samples_total <- dplyr::n_distinct(df$sample_id)
  pca_vars_present <- df |>
    dplyr::filter(.data$analyte %in% .env$pca_vars) |>
    dplyr::group_by(.data$analyte) |>
    dplyr::summarise(
      presence_freq = dplyr::n_distinct(.data$sample_id) / n_samples_total,
      .groups = "drop"
    ) |>
    dplyr::filter(.data$presence_freq >= min_detect_freq) |>
    dplyr::pull(.data$analyte)

  pca_vars_present <- union(
    intersect(required_vars, unique(df$analyte)),
    pca_vars_present
  )

  cli::cli_inform(c(
    "i" = "Chemistry PCA: {length(pca_vars_present)} variable{?s} available: \\
           {.val {sort(pca_vars_present)}}."
  ))

  # ── 4. Identify analyte groups ─────────────────────────────────────────────
  all_analytes <- unique(df$analyte)
  # Exclude all pca_vars and explicitly excluded analytes (microbiological
  # counts, qualitative/physical descriptors — not amenable to log-normal model)
  non_pca_non_excl <- setdiff(
    all_analytes,
    union(pca_vars, .IMPUTE_EXCLUDED)
  )

  excl_present <- intersect(all_analytes, .IMPUTE_EXCLUDED)
  if (length(excl_present) > 0L) {
    cli::cli_inform(c(
      "i" = "{length(excl_present)} analyte{?s} excluded from imputation \\
             (non-concentration data): {.val {sort(excl_present)}}."
    ))
  }

  metals_in_df   <- intersect(non_pca_non_excl, metal_analytes)
  organics_in_df <- setdiff(non_pca_non_excl, metal_analytes)

  # Drop target analytes detected too rarely to model. A brms regression needs
  # enough *detected* observations; near-/all-BDL analytes (e.g. ~100 organics
  # in a leachate panel) carry no signal and otherwise explode the model size.
  det_freq <- df |>
    dplyr::filter(.data$analyte %in% c(.env$metals_in_df, .env$organics_in_df)) |>
    dplyr::group_by(.data$analyte) |>
    dplyr::summarise(
      det_freq = dplyr::n_distinct(.data$sample_id[.data$detected]) / n_samples_total,
      .groups  = "drop"
    )
  keep_targets <- det_freq$analyte[det_freq$det_freq >= min_target_detect_freq]
  dropped <- setdiff(c(metals_in_df, organics_in_df), keep_targets)
  metals_in_df   <- intersect(metals_in_df,   keep_targets)
  organics_in_df <- intersect(organics_in_df, keep_targets)
  if (length(dropped) > 0L)
    cli::cli_inform(c(
      "i" = "Dropping {length(dropped)} target{?s} below \\
             min_target_detect_freq = {min_target_detect_freq}: \\
             {.val {sort(dropped)}}."
    ))

  cli::cli_inform(c(
    "i" = "Metals group: {length(metals_in_df)} analyte{?s}: \\
           {.val {sort(metals_in_df)}}.",
    "i" = "Organics group: {length(organics_in_df)} analyte{?s}: \\
           {.val {sort(organics_in_df)}}."
  ))

  if (length(metals_in_df) == 0L && length(organics_in_df) == 0L) {
    cli::cli_warn(
      "No target analytes found outside the PCA and excluded sets. \\
       Returning model with no fitted groups (imputation will be a no-op)."
    )
    return(structure(
      list(
        pca             = NULL,
        metals          = NULL,
        organics        = NULL,
        required_vars   = required_vars,
        pca_vars        = pca_vars,
        hurdle_metals   = metal_analytes,
        hurdle_organics = doc_like_analytes,
        fit_date        = Sys.Date(),
        n_samples       = length(samples_with_required)
      ),
      class = "imputation_model"
    ))
  }

  # ── 5. Fit unified chemistry PCA ──────────────────────────────────────────
  pca_obj <- .prepare_chem_pca(
    df, wq_vars        = pca_vars_present,
    min_var_explained  = min_var_explained,
    max_pcs            = max_pcs
  )

  cli::cli_inform(c(
    "i" = "Chemistry PCA: {pca_obj$n_pcs} axis/axes explain \\
           {round(100 * pca_obj$var_explained, 1)}% of variance."
  ))

  # ── 6. Fit metals model ────────────────────────────────────────────────────
  metals_fit <- NULL
  if (length(metals_in_df) >= 1L) {
    cli::cli_inform(c("i" = "Fitting metals model \u2026"))
    metals_fit <- .fit_group_model(
      df              = df,
      target_analytes = metals_in_df,
      pca_obj         = pca_obj,
      family          = family,
      iter            = iter,
      warmup          = warmup,
      chains          = chains,
      cores           = cores,
      impute_method   = impute_method,
      ...
    )
  }

  # ── 7. Fit organics model ──────────────────────────────────────────────────
  organics_fit <- NULL
  if (length(organics_in_df) >= 1L) {
    cli::cli_inform(c("i" = "Fitting organics model \u2026"))
    organics_fit <- .fit_group_model(
      df              = df,
      target_analytes = organics_in_df,
      pca_obj         = pca_obj,
      family          = family,
      iter            = iter,
      warmup          = warmup,
      chains          = chains,
      cores           = cores,
      impute_method   = impute_method,
      ...
    )
  }

  # ── 8. Assemble result ────────────────────────────────────────────────────
  result <- structure(
    list(
      pca             = pca_obj,
      metals          = metals_fit,
      organics        = organics_fit,
      required_vars   = required_vars,
      pca_vars        = pca_vars,
      hurdle_metals   = metal_analytes,
      hurdle_organics = doc_like_analytes,
      impute_method   = impute_method,
      fit_date        = Sys.Date(),
      n_samples       = length(samples_with_required)
    ),
    class = "imputation_model"
  )

  # ── 9. Save if requested ───────────────────────────────────────────────────
  if (!is.null(save_dir)) {
    if (!requireNamespace("qs2", quietly = TRUE))
      cli::cli_abort("Package {.pkg qs2} is required for saving models. \\
                      Install with: {.code install.packages('qs2')}")
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    fname     <- sprintf("imputation_model_%s.qs", format(Sys.Date(), "%Y%m%d"))
    save_path <- file.path(save_dir, fname)
    qs2::qs_save(result, save_path)
    cli::cli_inform(c("v" = "Model saved to {.path {save_path}}"))
    attr(result, "save_path") <- save_path
  }

  result
}

#' @export
print.imputation_model <- function(x, ...) {
  cat(sprintf(
    "<imputation_model>  fitted %s | %d samples | %d PCA vars | %d PCs (%.0f%% var)\n",
    x$fit_date, x$n_samples, length(x$pca_vars),
    x$pca$n_pcs, 100 * x$pca$var_explained
  ))
  if (!is.null(x$impute_method))
    cat(sprintf("  method:   %s\n", x$impute_method))
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
#' - *Organics*: at least one of DOC, TOC, BOD, COD or cBOD present at the
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
#' @examples
#' \dontrun{
#' model <- fit_imputation_model(monitoring_long)
#' # Point estimates (default), or full posterior draws with return = "draws":
#' imputed <- impute_chemistry(monitoring_long, model, return = "point")
#' }
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
  .require_brms()
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
      df           = result,
      group        = model$metals,
      pca_scores   = pca_scores,
      eligible_ids = eligible,
      return       = return
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
      df           = result,
      group        = model$organics,
      pca_scores   = pca_scores,
      eligible_ids = eligible_org,
      return       = return
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

  # The DL cap is the rescor_mi workaround (mi() does not enforce the censoring
  # bound during MCMC). The cens / cens_factor methods enforce it natively, so
  # their imputed values are not capped.
  if (identical(model$impute_method, "cens") ||
      identical(model$impute_method, "cens_factor")) {
    result
  } else {
    .check_bdl_imputed(result, dl_tbl, bdl_cap)
  }
}


# ── impute_coanalytes() ───────────────────────────────────────────────────────

#' Impute missing normalisation co-analytes from the fitted chemistry PCA
#'
#' Fits a univariate log-Gaussian GAM (`mgcv::gam`) for each target
#' co-analyte using the PC scores already computed by
#' [fit_imputation_model()].  Only samples where the co-analyte is entirely
#' absent are filled; BDL observations are left unchanged.
#'
#' This step belongs **after** [impute_chemistry()] and **before**
#' [time_weighted_aggregate()].  Imputed co-analyte values are never fed
#' back into the metals/organics model — the brms model ran on measured values
#' only and is already done.
#'
#' Using the chemistry PCA as the sole predictor set is appropriate because:
#' (a) the PCA already captures DOC/Ca/Mg variation in its axes; (b) a
#' univariate GAM on PC scores is unbiased and fast (no Stan required); (c)
#' the same PCA is used for the metals model so the co-analyte predictions
#' are conditioned on the same chemical environment summary.
#'
#' @param df Long-format chemistry data frame (same schema as
#'   [impute_chemistry()], with `imputed`/`imputed_kind` columns if
#'   [impute_chemistry()] has already been called).
#' @param model Fitted model from [fit_imputation_model()] (provides the PCA
#'   object and the list of `pca_vars`).
#' @param targets Co-analyte names to impute when missing.  Default
#'   `.COANALYTE_TARGETS` (`"DOC"`, `"Ca"`, `"Mg"`, `"hardness"`).  Only
#'   targets present in `model$pca_vars` are processed; others are
#'   skipped with a warning.
#' @param min_obs Minimum number of quantified observations required to fit a
#'   GAM for a target.  Targets with fewer observations are skipped.
#'   Default `10L`.
#'
#' @return `df` with missing co-analyte rows filled in, tagged with
#'   `imputed = TRUE` and `imputed_kind = "missing"`.  All other rows are
#'   unchanged.
#'
#' @seealso [fit_imputation_model()], [impute_chemistry()]
#' @examples
#' \dontrun{
#' # Deterministic GAM-based imputation of normalisation co-analytes
#' # (pH, DOC, hardness, ...) from the measured analyte suite.
#' impute_coanalytes(monitoring_long)
#' }
#' @export
impute_coanalytes <- function(
    df,
    model,
    targets = NULL,
    min_obs = 10L
) {
  if (is.null(targets)) targets <- .COANALYTE_TARGETS
  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected"))
  if (!inherits(model, "imputation_model"))
    cli::cli_abort(
      "{.arg model} must be an object returned by {.fn fit_imputation_model}."
    )
  if (is.null(model$pca))
    cli::cli_abort(
      "Model has no fitted PCA \u2014 did {.fn fit_imputation_model} find any \\
       target analytes?"
    )

  # ── Compute PC scores for all samples ─────────────────────────────────────
  pca_scores <- .compute_pca_scores(df, model$pca)
  pc_cols    <- paste0("PC", seq_len(model$pca$n_pcs))

  # Skip targets not represented in the PCA (they can't be predicted)
  targets_ok <- intersect(targets, model$pca_vars)
  skipped    <- setdiff(targets, model$pca_vars)
  if (length(skipped) > 0L)
    cli::cli_warn(c(
      "!" = "{length(skipped)} co-analyte target{?s} not in fitted \\
             {.arg pca_vars} \u2014 skipping: {.val {skipped}}."
    ))

  # Per-sample metadata for constructing new rows
  sample_meta <- df |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select("sample_id", "site_id", "datetime")

  result <- df

  for (tgt in targets_ok) {

    # "Present" means at least one row exists (detected or BDL)
    present_ids <- unique(result$sample_id[result$analyte == tgt])
    missing_ids <- setdiff(unique(result$sample_id), present_ids)

    if (length(missing_ids) == 0L) next   # already complete

    # Count quantified observations available for GAM fitting
    n_obs <- dplyr::n_distinct(
      result$sample_id[result$analyte == tgt & result$detected]
    )
    if (n_obs < min_obs) {
      cli::cli_warn(c(
        "!" = "Co-analyte {.val {tgt}}: only {n_obs} quantified sample{?s} \\
               (< {.arg min_obs} = {min_obs}) \u2014 skipping."
      ))
      next
    }

    cli::cli_inform(c(
      "i" = "Co-analyte {.val {tgt}}: imputing {length(missing_ids)} \\
             missing sample{?s} via GAM on {model$pca$n_pcs} PC score{?s}."
    ))

    # ── Fit GAM on quantified observations ──────────────────────────────────
    obs_vals <- result |>
      dplyr::filter(.data$analyte == tgt, .data$detected) |>
      dplyr::group_by(.data$sample_id) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::transmute(.data$sample_id,
                        log_tgt = log(pmax(.data$value, 1e-9)))

    gam_data <- pca_scores |>
      dplyr::filter(.data$sample_id %in% present_ids) |>
      dplyr::left_join(obs_vals, by = "sample_id") |>
      dplyr::filter(!is.na(.data$log_tgt))

    gam_formula <- stats::as.formula(
      paste("log_tgt ~",
            paste(paste0("s(", pc_cols, ")"), collapse = " + "))
    )

    gam_fit <- tryCatch(
      mgcv::gam(gam_formula, data = gam_data, family = stats::gaussian()),
      error = function(e) {
        cli::cli_warn(c(
          "!" = "GAM fit failed for {.val {tgt}}: {conditionMessage(e)}.",
          "i" = "Skipping imputation for this co-analyte."
        ))
        NULL
      }
    )
    if (is.null(gam_fit)) next

    # ── Predict for missing samples ──────────────────────────────────────────
    pred_data <- dplyr::filter(pca_scores, .data$sample_id %in% missing_ids)
    pred_vals <- exp(
      as.numeric(stats::predict(gam_fit, newdata = pred_data, type = "response"))
    )

    new_rows <- tibble::tibble(
      sample_id    = pred_data$sample_id,
      analyte      = tgt,
      value        = pred_vals,
      detected     = TRUE,
      imputed      = TRUE,
      imputed_kind = "missing"
    ) |>
      dplyr::left_join(sample_meta, by = "sample_id")

    result <- dplyr::bind_rows(result, new_rows)
  }

  # ── Ensure imputed/imputed_kind columns are populated on all rows ──────────
  if (!"imputed" %in% names(result)) {
    result <- dplyr::mutate(result, imputed = FALSE, imputed_kind = "observed")
  } else {
    result <- dplyr::mutate(result,
      imputed      = dplyr::coalesce(.data$imputed, FALSE),
      imputed_kind = dplyr::coalesce(.data$imputed_kind, "observed")
    )
  }

  dplyr::arrange(result, .data$sample_id, .data$analyte)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Log10-transform the concentration-like columns of a PCA matrix
#'
#' Applies `log10()` to every column whose name is not in [.PCA_NO_LOG_VARS],
#' leaving pH / temperature / ORP / DO on their native scale.  `NA` cells are
#' preserved (so NIPALS can still handle within-sample missingness), and a small
#' floor guards against `log10(0)` for genuine zeros.  Both the training PCA
#' (`.prepare_chem_pca()`) and the scoring projection (`.compute_pca_scores()`)
#' call this so the transform is identical on both paths.
#' @param mat Numeric matrix with named columns (samples × variables).
#' @param eps Positive floor applied before the log (default `1e-9`).
#' @keywords internal
.log_transform_pca <- function(mat, eps = 1e-9) {
  log_cols <- setdiff(colnames(mat), .PCA_NO_LOG_VARS)
  if (length(log_cols) > 0L) {
    mat[, log_cols] <- log10(pmax(mat[, log_cols, drop = FALSE], eps))
  }
  mat
}

#' Fit the unified chemistry PCA on training data
#'
#' This PCA spans the full unified chemistry predictor set (`pca_vars` in
#' `fit_imputation_model()`) — major ions, pH, EC, NH3-N, DOC, nutrients and
#' redox indicators.  Concentration-like variables are `log10`-transformed (see
#' [.log_transform_pca()]) before centring/scaling.  PC score columns are named
#' `PC1`, `PC2`, ….
#' @keywords internal
.prepare_chem_pca <- function(df, wq_vars, min_var_explained = 0.75, max_pcs = 4L) {
  # Pivot chemistry vars to wide (one row per sample); missing cells → NA
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

  # Training medians — kept as fallback for columns entirely absent in
  # scoring data.  Per-cell NAs within a column are handled by nipals.
  # Stored on the RAW scale: scoring fills absent columns with these medians
  # and then re-applies the same log transform below.
  train_medians <- apply(wq_matrix, 2, stats::median, na.rm = TRUE)

  # Log10-transform concentration-like variables (everything except
  # pH / temperature / ORP / DO).  Done before centring/scaling so the PCA
  # reflects multiplicative chemical variation rather than being dominated by
  # the highest-magnitude major ions.  NAs are preserved for NIPALS.
  wq_matrix <- .log_transform_pca(wq_matrix)

  # Remove zero-variance or all-NA columns
  col_sds   <- apply(wq_matrix, 2, stats::sd, na.rm = TRUE)
  keep_cols <- !is.na(col_sds) & col_sds > 0
  if (!any(keep_cols))
    cli::cli_abort("All PCA variables have zero variance \u2014 cannot fit PCA.")
  if (!all(keep_cols)) {
    dropped <- colnames(wq_matrix)[!keep_cols]
    cli::cli_inform(c("!" = "Chemistry PCA: dropping zero-variance variable{?s}: \\
                              {.val {dropped}}"))
    wq_matrix <- wq_matrix[, keep_cols, drop = FALSE]
  }

  # NIPALS PCA — handles missing cells without prior imputation
  ncomp   <- min(max_pcs, ncol(wq_matrix), nrow(wq_matrix) - 1L)
  pca_fit <- nipals::nipals(wq_matrix, ncomp = ncomp, center = TRUE, scale = TRUE)
  # nipals (>= 1.0) returns per-component proportions in `$R2`; older/other
  # builds expose a cumulative `$R2cum`. Accept either, deriving the cumulative
  # curve from `$R2` when needed.
  cum_var <- if (!is.null(pca_fit$R2cum)) {
    pca_fit$R2cum
  } else {
    cumsum(pca_fit$R2)
  }

  # Determine number of PCs
  n_needed <- which(cum_var >= min_var_explained)[1L]
  if (is.na(n_needed) || n_needed > max_pcs) {
    n_pcs <- min(max_pcs, length(cum_var))
    cli::cli_warn(c(
      "!" = "Chemistry PCA: first {n_pcs} axis/axes explain only \\
             {round(100 * cum_var[n_pcs], 1)}% of variance \\
             (target: {min_var_explained * 100}%).",
      "i" = "Consider adding more {.arg pca_vars} or loosening \\
             {.arg min_var_explained}."
    ))
  } else {
    n_pcs <- max(2L, n_needed)
  }
  n_pcs <- min(n_pcs, length(cum_var))

  pca_obj <- list(
    fit           = pca_fit,
    medians       = train_medians,         # all WQ vars (before zero-var removal)
    active_vars   = colnames(wq_matrix),   # after zero-var removal
    n_pcs         = n_pcs,
    var_explained = cum_var[n_pcs]
  )

  # Training scores MUST be produced by the same projection used at prediction
  # time (`.compute_pca_scores()`), otherwise the brms model is trained on one
  # score scale and predicted on another.  `nipals$scores` are NOT equal to the
  # regression projection used at scoring — they differ by a per-component
  # factor (the component eigenvalue) — so copying them here would silently
  # break imputation.  Deriving `pc_scores` via `.compute_pca_scores()` on the
  # training data guarantees train/predict consistency by construction.
  pca_obj$pc_scores <- .compute_pca_scores(df, pca_obj)

  pca_obj
}

#' Project new data onto stored NIPALS PCA axes
#'
#' Handles within-row missing values via NIPALS regression scoring: each
#' component score is estimated from observed variables only, then the residual
#' is deflated before the next component.  Columns entirely absent in `df` (not
#' measured at all, not just BDL) are filled with training medians.
#' @keywords internal
.compute_pca_scores <- function(df, pca_obj) {
  wq_vars <- pca_obj$active_vars  # columns used in training

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

  # Add entirely-absent training columns using training medians.
  # Per-cell NAs within a column are handled by the NIPALS scoring below.
  for (v in wq_vars) {
    if (!v %in% names(wq_wide)) wq_wide[[v]] <- pca_obj$medians[[v]]
  }

  wq_mat <- as.matrix(dplyr::select(wq_wide, dplyr::all_of(wq_vars)))

  # Fill columns that are entirely NA (not measured in this dataset at all)
  for (j in seq_len(ncol(wq_mat))) {
    if (all(is.na(wq_mat[, j]))) {
      col_nm       <- colnames(wq_mat)[j]
      wq_mat[, j]  <- pca_obj$medians[[col_nm]]
    }
  }

  # Apply the same log10 transform used at training time.  Done AFTER the
  # raw-scale median fills above so filled values are transformed identically
  # to how the training medians were, keeping centre/scale parameters valid.
  wq_mat <- .log_transform_pca(wq_mat)

  # Centre and scale using training parameters stored in the nipals object
  wq_scaled <- sweep(wq_mat,    2, pca_obj$fit$center, "-")
  wq_scaled <- sweep(wq_scaled, 2, pca_obj$fit$scale,  "/")

  # NIPALS regression scoring: per-row, per-component
  loadings   <- pca_obj$fit$loadings[, seq_len(pca_obj$n_pcs), drop = FALSE]
  scores_mat <- t(apply(wq_scaled, 1, .nipals_score_row,
                        loadings = loadings, n_pcs = pca_obj$n_pcs))
  colnames(scores_mat) <- paste0("PC", seq_len(pca_obj$n_pcs))

  tibble::as_tibble(scores_mat) |>
    dplyr::mutate(sample_id = wq_wide$sample_id)
}

#' NIPALS regression scoring for one centred/scaled observation
#'
#' Computes PC scores for a single row by projecting onto each loading vector
#' using only observed (non-NA) elements, then deflating the residual before
#' moving to the next component.  For a fully observed row this is identical to
#' the standard `x %*% loadings` projection.  For a row with missing values it
#' correctly down-weights the loading vectors to the observed subspace — without
#' the bias that zero/median imputation introduces.
#'
#' @param x Numeric vector (centred + scaled); `NA` marks missing variables.
#' @param loadings p × K loading matrix (unit-normalised columns from nipals).
#' @param n_pcs Number of components to score.
#' @keywords internal
.nipals_score_row <- function(x, loadings, n_pcs) {
  scores  <- numeric(n_pcs)
  x_resid <- x
  for (k in seq_len(n_pcs)) {
    lk  <- loadings[, k]
    obs <- !is.na(x_resid)
    if (any(obs)) {
      # Regression onto observed sub-vector of the loading; renormalise by
      # sum(lk[obs]^2) because the loading is unit-normalised over *all* p vars
      scores[k] <- sum(x_resid[obs] * lk[obs]) / sum(lk[obs]^2)
    }
    # Deflate: remove this component's contribution from observed elements
    x_resid[obs] <- x_resid[obs] - scores[k] * lk[obs]
  }
  scores
}

#' Fit a single brms group model (metals or organics)
#' @keywords internal
.fit_group_model <- function(df, target_analytes, pca_obj,
                              family, iter, warmup, chains, cores,
                              impute_method = "rescor_mi", ...) {
  eps_log <- 1e-9

  safe_analytes <- stats::setNames(make.names(target_analytes), target_analytes)
  # safe_analytes: names = safe R names, values = original names
  safe_vec <- unname(safe_analytes)
  pc_cols  <- paste0("PC", seq_len(pca_obj$n_pcs))
  rhs      <- paste(paste0("s(", pc_cols, ")"), collapse = " + ")
  pc_wide  <- dplyr::select(pca_obj$pc_scores, "sample_id", dplyr::all_of(pc_cols))

  # One row per (sample, analyte) with a safe column name and log value.
  base <- df |>
    dplyr::filter(.data$analyte %in% .env$target_analytes) |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      safe = unname(safe_analytes[.data$analyte]),
      lv   = log(pmax(.data$value, eps_log))
    )

  if (impute_method == "rescor_mi") {
    # Residual correlation + mi() for BDL/missing. BDL and missing are NA and
    # imputed; the post-hoc DL cap is applied in impute_chemistry().
    target_wide <- base |>
      dplyr::mutate(log_value = dplyr::if_else(.data$detected, .data$lv, NA_real_)) |>
      dplyr::select("sample_id", "safe", "log_value") |>
      tidyr::pivot_wider(names_from = "safe", values_from = "log_value")
    wide_df <- dplyr::left_join(pc_wide, target_wide, by = "sample_id")
    bf_list <- purrr::map(safe_vec, function(s)
      brms::bf(stats::as.formula(paste0(s, " | mi() ~ ", rhs))))
    brms_formula <- do.call(brms::mvbf, c(bf_list, list(rescor = TRUE)))

  } else {
    # cens / cens_factor: left-censor BDL at its detection limit (the BDL value),
    # rescor = FALSE. subset() lets each analyte use only the samples that
    # measured it; cens_factor adds a shared per-sample latent factor that
    # induces cross-analyte coupling.
    resp_wide <- base |>
      dplyr::select("sample_id", "safe", "lv") |>
      tidyr::pivot_wider(names_from = "safe", values_from = "lv")
    cens_wide <- base |>
      dplyr::mutate(cf = dplyr::if_else(.data$detected, "none", "left")) |>
      dplyr::select("sample_id", "safe", "cf") |>
      tidyr::pivot_wider(names_from = "safe", values_from = "cf",
                         names_glue = "cens_{safe}")
    wide_df <- pc_wide |>
      dplyr::left_join(resp_wide, by = "sample_id") |>
      dplyr::left_join(cens_wide, by = "sample_id")
    for (s in safe_vec) {
      wide_df[[paste0("sub_", s)]]  <- !is.na(wide_df[[s]])
      wide_df[[s]]                  <- dplyr::coalesce(wide_df[[s]], 0)
      wide_df[[paste0("cens_", s)]] <- dplyr::coalesce(wide_df[[paste0("cens_", s)]], "none")
    }
    grp <- if (impute_method == "cens_factor") " + (1 |q| sample_id)" else ""
    bf_list <- purrr::map(safe_vec, function(s)
      brms::bf(stats::as.formula(sprintf(
        "%s | cens(cens_%s) + subset(sub_%s) ~ %s%s", s, s, s, rhs, grp))))
    brms_formula <- do.call(brms::mvbf, c(bf_list, list(rescor = FALSE)))
  }

  cli::cli_inform(c(
    "i" = "brms ({impute_method}): {length(target_analytes)} analyte{?s} \u00d7 \\
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
    fit             = fit,
    analytes        = target_analytes,   # original names
    safe_names      = safe_analytes,     # names=safe, values=original
    pc_cols         = pc_cols,
    wide_sample_ids = wide_df$sample_id,
    impute_method   = impute_method
  )
}

#' Predict and merge imputed values for one analyte group
#' @keywords internal
.predict_and_merge <- function(df, group, pca_scores, eligible_ids, return) {
  eps_log <- 1e-9

  target_analytes <- group$analytes
  safe_analytes   <- group$safe_names   # names=safe, values=original
  pc_cols         <- group$pc_cols

  # ── Build wide prediction df for eligible samples ─────────────────────────
  df_eligible <- dplyr::filter(df, .data$sample_id %in% .env$eligible_ids)
  if (nrow(df_eligible) == 0L) return(df)

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

  # PC scores for eligible samples
  pc_wide <- dplyr::filter(pca_scores, .data$sample_id %in% .env$eligible_ids) |>
    dplyr::select("sample_id", dplyr::all_of(pc_cols))

  # Ensure all training analyte columns are present (add NA if absent)
  for (s in names(safe_analytes)) {
    if (!s %in% names(target_wide)) target_wide[[s]] <- NA_real_
  }

  wide_new <- pc_wide |>
    dplyr::left_join(target_wide, by = "sample_id")

  # cens / cens_factor models reference <safe>__cens and <safe>__sub columns in
  # the formula. For prediction we set them so every target cell is predicted
  # (cens "none", subset TRUE); the response placeholder is ignored.
  if (!is.null(group$impute_method) && group$impute_method != "rescor_mi") {
    for (s in unname(safe_analytes)) {
      if (!s %in% names(wide_new)) wide_new[[s]] <- 0
      wide_new[[s]]                 <- dplyr::coalesce(wide_new[[s]], 0)
      wide_new[[paste0("cens_", s)]] <- "none"
      wide_new[[paste0("sub_", s)]]  <- TRUE
    }
  }

  # ── Posterior predictions ─────────────────────────────────────────────────
  cens_method <- !is.null(group$impute_method) &&
    group$impute_method != "rescor_mi"
  if (cens_method) {
    # Models with subset() must be predicted one response at a time (brms
    # disallows joint prediction), so loop over analytes and assemble pm_long.
    pm_long <- purrr::map_dfr(unname(safe_analytes), function(s) {
      orig <- names(safe_analytes)[safe_analytes == s]
      if (return == "point") {
        ep <- brms::posterior_epred(group$fit, newdata = wide_new, resp = s,
                                    allow_new_levels = TRUE)
        tibble::tibble(sample_id = wide_new$sample_id, analyte = orig,
                       .post_mean = exp(colMeans(ep)))
      } else {
        pp <- brms::posterior_predict(group$fit, newdata = wide_new, resp = s,
                                      allow_new_levels = TRUE)
        nd <- nrow(pp)
        tibble::tibble(
          sample_id   = rep(wide_new$sample_id, each = nd),
          analyte     = orig,
          draw_id     = rep(seq_len(nd), times = ncol(pp)),
          .post_value = exp(as.vector(pp)))
      }
    })
  } else if (return == "point") {
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
