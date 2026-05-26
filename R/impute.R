#' Impute missing and below-detection-limit chemistry using a Bayesian
#' multivariate GAM
#'
#' Fits a multivariate log-normal GAM via `brms::brm()` with
#' `set_rescor(TRUE)` to capture cross-analyte covariance. Left-censored
#' observations (below-detection-limit, BDL) and missing analyte measurements
#' are treated as latent values and imputed from the multivariate joint
#' posterior using `mi()` (missing indicator). Missing analyte measurements at
#' a sample are imputed from the cross-analyte covariance structure and the
#' driver predictors.
#'
#' After fitting, the function returns a long-format data frame with all BDL
#' and missing cells replaced by posterior means (or draws), plus an `imputed`
#' column to flag which values were filled in.
#'
#' **Note on BDL handling:** `brms` does not support `cens()` (left-censored
#' likelihood) when `rescor = TRUE` is set — only `se`, `weights`, and `mi`
#' are supported addition arguments when residual correlations are estimated.
#' This function therefore uses `mi()` for both BDL and missing observations:
#' BDL values are set to `NA` before fitting and imputed from the multivariate
#' posterior. The strict left-censor constraint (concentration ≤ detection
#' limit) is not enforced during fitting. A post-hoc check (Option B) warns
#' and optionally caps imputed BDL values that exceed the original detection
#' limit.
#'
#' **Future improvement — Option D (custom Stan model):** A statistically
#' correct treatment would use a custom Stan likelihood that combines
#' left-censoring with cross-analyte residual correlations. This would enforce
#' the constraint `imputed_C ≤ DL` for BDL observations throughout the MCMC
#' chain rather than as a post-hoc correction. Option D is feasible but
#' requires bypassing brms and writing a custom Stan program. It is listed as
#' a planned future improvement.
#'
#' **Installation note:** `brms` requires a working Stan toolchain (rstan or
#' cmdstanr). If Stan is not installed, this function will abort with a message
#' pointing to the installation documentation.
#'
#' @param df Long-format chemistry data frame. Required columns:
#'   - `sample_id` (character)
#'   - `site_id` (character)
#'   - `datetime` (Date or POSIXct)
#'   - `analyte` (character)
#'   - `value` (numeric) — for BDL rows, set to the detection limit
#'   - `detected` (logical) — `FALSE` for BDL, `TRUE` for detected values
#'   Driver analytes (e.g. pH, EC, DOC) must be present as rows in `df` and
#'   must have `detected == TRUE` for every sample.
#' @param drivers Character vector of analyte names to use as GAM predictors.
#'   Default `c("pH", "EC", "DOC")`. These analytes are extracted from `df`,
#'   pivoted to wide columns, and used as spline predictors. They are
#'   **not** themselves imputed.
#' @param formula_template `brms::brmsformula` or `NULL`. If `NULL`, a default
#'   formula is constructed: `mvbind(analyte_1, ...) ~ s(d1) + s(d2) + ...`
#'   with one smooth per driver.
#' @param family brms family for the log-transformed response. Default
#'   `"gaussian"`. Concentrations are always log-transformed before fitting
#'   (i.e. the model is `log(conc) ~ ...`), so `"gaussian"` produces
#'   log-normal marginals. `"gaussian"` is required for `set_rescor(TRUE)` to
#'   be available in brms; `"lognormal"` is not supported with residual
#'   correlations. Posterior predictions are back-transformed with `exp()` so
#'   the returned `value` column is always on the original concentration scale.
#' @param iter Total MCMC iterations per chain. Default `2000`.
#' @param warmup Warmup iterations per chain. Default `1000`.
#' @param chains Number of MCMC chains. Default `4`.
#' @param cores Number of parallel cores. Default `parallel::detectCores()`.
#' @param bdl_cap Logical. If `TRUE` (default), imputed values for BDL
#'   observations that exceed the original detection limit are capped at the
#'   detection limit and a warning is issued. If `FALSE`, no capping is applied
#'   but a warning is still issued when imputed BDL values exceed the DL.
#'   Set to `FALSE` if you want to inspect the uncapped posterior means.
#' @param return Either `"point"` (default) to return posterior mean estimates
#'   for imputed cells, or `"draws"` to return one row per
#'   (sample × analyte × draw).
#' @param ... Additional arguments passed to `brms::brm()`.
#'
#' @return When `return = "point"`: the input `df` with imputed `value` for
#'   BDL and missing rows, plus two new columns:
#'   - `imputed` (logical) — `TRUE` for BDL or missing rows (now filled)
#'   - `imputed_kind` (character) — `"observed"`, `"censored_left"`, or
#'     `"missing"`
#'   The fitted `brmsfit` object is stored as `attr(result, "brmsfit")`.
#'
#'   When `return = "draws"`: as above but with an additional `draw_id`
#'   (integer) column and one row per (sample × analyte × draw). Heavier but
#'   propagates posterior uncertainty into `compute_chronic_chemistry()`.
#'
#' @examples
#' \dontrun{
#' imp <- impute_chemistry(chem_filtered, drivers = c("pH", "EC", "DOC"))
#' fit <- attr(imp, "brmsfit")
#' brms::pp_check(fit)
#' }
#'
#' @export
impute_chemistry <- function(
    df,
    drivers          = c("pH", "EC", "DOC"),
    formula_template = NULL,
    family           = "gaussian",
    iter             = 2000,
    warmup           = 1000,
    chains           = 4,
    cores            = parallel::detectCores(),
    bdl_cap          = TRUE,
    return           = c("point", "draws"),
    ...
) {
  return <- match.arg(return)

  # ── Input validation ───────────────────────────────────────────────────────
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected")
  )
  checkmate::assert_character(drivers, min.len = 1L, any.missing = FALSE)
  checkmate::assert_count(iter)
  checkmate::assert_count(warmup)
  checkmate::assert_count(chains)
  checkmate::assert_flag(bdl_cap)

  missing_drivers <- setdiff(drivers, unique(df$analyte))
  if (length(missing_drivers) > 0L) {
    cli::cli_abort(c(
      "Driver analyte{?s} not found in {.arg df}: {.val {missing_drivers}}.",
      "i" = "Drivers must be present as rows in {.arg df} with \\
             {.code analyte} matching the driver names."
    ))
  }

  # ── Separate targets from drivers ─────────────────────────────────────────
  driver_df  <- dplyr::filter(df, .data$analyte %in% .env$drivers)
  target_df  <- dplyr::filter(df, !(.data$analyte %in% .env$drivers))

  target_analytes <- unique(target_df$analyte)
  if (length(target_analytes) == 0L) {
    cli::cli_abort("No non-driver analytes found in {.arg df} to impute.")
  }

  # ── Check drivers are detected everywhere ─────────────────────────────────
  driver_missing <- driver_df |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::summarise(
      any_missing = any(!.data$detected),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$any_missing)

  if (nrow(driver_missing) > 0L) {
    bad <- unique(driver_missing$analyte)
    cli::cli_abort(c(
      "Driver analyte{?s} have missing or BDL values: {.val {bad}}.",
      "i" = "Drivers must be fully detected at every sample. \\
             Either remove BDL samples or choose different drivers."
    ))
  }

  # ── Save detection limits for BDL rows (Option B post-hoc check) ──────────
  # The `value` column stores the detection limit for BDL (detected = FALSE) rows.
  dl_tbl <- dplyr::filter(target_df, !.data$detected) |>
    dplyr::select("sample_id", "analyte", detection_limit = "value")

  # ── Pivot to wide format for brms ─────────────────────────────────────────
  # Drivers: one column each
  driver_wide <- driver_df |>
    dplyr::select("sample_id", "analyte", "value") |>
    tidyr::pivot_wider(
      names_from  = "analyte",
      values_from = "value",
      names_prefix = ".drv_"
    )

  driver_col_names <- paste0(".drv_", drivers)

  # Targets: log(value) for detected rows; NA for BDL and missing.
  #
  # We always fit on the log scale with family = "gaussian" and rescor = TRUE.
  # brms does not support cens() when rescor is estimated (only se/weights/mi
  # are allowed). Using mi() handles both BDL and truly-missing analytes as
  # latent missing values imputed from the multivariate joint distribution.
  #
  # Statistical consequence: we lose the strict left-censor constraint for BDL
  # values. In practice this is acceptable because:
  #   (a) the cross-analyte correlations constrain imputed values to plausible
  #       ranges conditioned on all co-analytes,
  #   (b) BDL values are typically far below the detected concentration range,
  #       so the posterior naturally stays low.
  # A post-hoc check (Option B) flags and caps exceedances below.
  eps_log <- 1e-9   # guard against log(0)
  target_wide_vals <- target_df |>
    dplyr::select("sample_id", "analyte", "value", "detected") |>
    dplyr::mutate(
      log_value = dplyr::if_else(
        .data$detected,
        log(pmax(.data$value, eps_log)),
        NA_real_   # BDL → NA → mi() imputes from multivariate posterior
      )
    ) |>
    dplyr::select("sample_id", "analyte", "log_value") |>
    tidyr::pivot_wider(
      names_from  = "analyte",
      values_from = "log_value"
    )

  # No separate censoring-indicator columns needed for mi().
  wide_df <- driver_wide |>
    dplyr::left_join(target_wide_vals, by = "sample_id")

  # ── Track imputation kind for each (sample, analyte) ─────────────────────
  impute_kind <- target_df |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      .imputed_kind = dplyr::case_when(
        !.data$detected ~ "censored_left",
        TRUE            ~ "observed"
      )
    ) |>
    dplyr::select("sample_id", "analyte", ".imputed_kind")

  # Samples with no measurement at all for an analyte: "missing"
  all_combos <- tidyr::expand_grid(
    sample_id = unique(target_df$sample_id),
    analyte   = target_analytes
  )
  impute_kind <- dplyr::left_join(all_combos, impute_kind,
                                  by = c("sample_id", "analyte")) |>
    dplyr::mutate(
      .imputed_kind = dplyr::if_else(
        is.na(.data$.imputed_kind), "missing", .data$.imputed_kind
      )
    )

  # ── Build brms formula ────────────────────────────────────────────────────
  if (is.null(formula_template)) {
    # Clean analyte names for use as R variable names in the formula
    # (brms requires syntactically valid column names)
    safe_analytes <- make.names(target_analytes)
    names(safe_analytes) <- target_analytes

    # Rename analyte columns in wide_df to safe R variable names
    rename_map <- stats::setNames(safe_analytes, target_analytes)
    wide_df    <- dplyr::rename(wide_df, dplyr::any_of(rename_map))

    driver_terms <- paste0("s(", driver_col_names, ")", collapse = " + ")

    bf_list <- purrr::map(unname(safe_analytes), function(a) {
      brms::bf(
        stats::as.formula(paste0(a, " | mi() ~ ", driver_terms))
      )
    })

    brms_formula <- do.call(brms::mvbf, c(bf_list, list(rescor = TRUE)))
  } else {
    brms_formula <- formula_template
    safe_analytes <- make.names(target_analytes)
    names(safe_analytes) <- target_analytes
  }

  # ── Fit the model ─────────────────────────────────────────────────────────
  cli::cli_inform(c(
    "i" = "Fitting brms multivariate GAM: \\
           {length(target_analytes)} analyte{?s} × \\
           {nrow(wide_df)} sample{?s}. \\
           This may take several minutes."
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
               ({.url https://paul-buerkner.github.io/brms/}).",
        "i" = "Try with {.code iter = 500, chains = 1} for a quick smoke test."
      ))
    }
  )

  # ── Extract posterior estimates for imputed cells ─────────────────────────
  if (return == "draws") {
    post_draws <- brms::posterior_predict(fit)  # [draws × obs × response]
    # Reshape to long format: (sample_id, analyte, draw_id, value)
    post_long <- .reshape_posterior_draws(
      post_draws, wide_df$sample_id, safe_analytes
    )

    result <- .build_imputed_df(df, post_long, impute_kind, drivers,
                                return = "draws")
    result <- .check_bdl_imputed(result, dl_tbl, bdl_cap)
    attr(result, "brmsfit") <- fit
    return(result)
  }

  # Point estimates: posterior expected values (mean of posterior draws).
  # brms::posterior_epred() is the modern API for E[Y | data]; it handles
  # mi() imputed cells correctly (returns posterior of the latent missing value).
  # Returns array [n_draws, n_obs, n_responses] on the log scale.
  epred_draws <- brms::posterior_epred(fit)
  pm_long <- .reshape_posterior_means(
    epred_draws, wide_df$sample_id, safe_analytes
  )

  result <- .build_imputed_df(df, pm_long, impute_kind, drivers,
                              return = "point")
  result <- .check_bdl_imputed(result, dl_tbl, bdl_cap)
  attr(result, "brmsfit") <- fit
  result
}

# ── Internal reshape helpers ──────────────────────────────────────────────────

#' Reshape brms::posterior_epred() output to long format of posterior means
#'
#' `posterior_epred()` returns `[n_draws, n_obs, n_responses]` on the log scale.
#' We average across draws and exp()-transform to get concentrations.
#' @keywords internal
.reshape_posterior_means <- function(epred_draws, sample_ids, safe_analytes) {
  # Handle both univariate (matrix [n_draws, n_obs]) and multivariate (3D array)
  if (is.matrix(epred_draws)) {
    # Single response — wrap in a named list
    arr_list <- list(epred_draws)
    names(arr_list) <- unname(safe_analytes)[1L]
  } else {
    # [n_draws, n_obs, n_responses]; response names in dim 3
    resp_nms <- dimnames(epred_draws)[[3L]]
    arr_list <- setNames(
      lapply(seq_along(resp_nms), function(i) epred_draws[, , i]),
      resp_nms
    )
  }

  purrr::map2_dfr(arr_list, names(arr_list), function(mat, safe_nm) {
    orig_nm <- names(safe_analytes)[safe_analytes == safe_nm]
    if (length(orig_nm) == 0L) orig_nm <- safe_nm  # fallback
    # mat is [n_draws, n_obs]; colMeans gives the posterior mean per observation
    tibble::tibble(
      sample_id    = sample_ids,
      analyte      = orig_nm,
      .post_mean   = exp(colMeans(mat))  # back-transform from log scale
    )
  })
}

#' Reshape brms::posterior_predict() output to long format
#' @keywords internal
.reshape_posterior_draws <- function(post_draws, sample_ids, safe_analytes) {
  # post_draws: [n_draws, n_obs, n_responses]
  n_draws <- dim(post_draws)[1L]
  resp_nms <- dimnames(post_draws)[[3L]]

  purrr::map2_dfr(
    seq_along(resp_nms), resp_nms,
    function(ri, safe_nm) {
      orig_nm <- names(safe_analytes)[safe_analytes == safe_nm]
      if (length(orig_nm) == 0L) orig_nm <- safe_nm
      mat <- post_draws[, , ri]  # [n_draws, n_obs]
      tibble::tibble(
        sample_id   = rep(sample_ids, each = n_draws),
        analyte     = orig_nm,
        draw_id     = rep(seq_len(n_draws), times = length(sample_ids)),
        .post_value = exp(as.vector(t(mat)))  # back-transform from log scale
      )
    }
  )
}

#' Merge posterior estimates back into the original long-format df
#' @keywords internal
.build_imputed_df <- function(df, post_long, impute_kind, drivers,
                              return = "point") {
  target_df  <- dplyr::filter(df, !(.data$analyte %in% .env$drivers))
  driver_df  <- dplyr::filter(df, .data$analyte %in% .env$drivers)

  # Tag existing observed rows
  target_df <- dplyr::left_join(
    target_df,
    impute_kind,
    by = c("sample_id", "analyte")
  ) |>
    dplyr::mutate(
      imputed      = .data$.imputed_kind != "observed",
      imputed_kind = .data$.imputed_kind
    ) |>
    dplyr::select(-".imputed_kind")

  if (return == "point") {
    val_col <- ".post_mean"
    # Replace BDL / missing values with posterior means
    imputed_rows  <- dplyr::filter(target_df, .data$imputed)
    observed_rows <- dplyr::filter(target_df, !.data$imputed)

    imputed_filled <- dplyr::left_join(
      imputed_rows,
      post_long |> dplyr::select("sample_id", "analyte", dplyr::all_of(val_col)),
      by = c("sample_id", "analyte")
    ) |>
      dplyr::mutate(
        value    = dplyr::coalesce(.data[[val_col]], .data$value),
        detected = TRUE
      ) |>
      dplyr::select(-dplyr::all_of(val_col))

    result_targets <- dplyr::bind_rows(observed_rows, imputed_filled)

  } else {
    # Draws: expand all target rows by n_draws, replace imputed with draws
    n_draws <- max(post_long$draw_id)
    imputed_rows  <- dplyr::filter(target_df, .data$imputed)
    observed_rows <- dplyr::filter(target_df, !.data$imputed)

    # Expand observed rows across all draws
    obs_expanded <- tidyr::crossing(
      observed_rows,
      draw_id = seq_len(n_draws)
    )

    # Imputed rows: join draws
    imp_expanded <- dplyr::left_join(
      tidyr::crossing(imputed_rows, draw_id = seq_len(n_draws)),
      post_long |>
        dplyr::select("sample_id", "analyte", "draw_id", ".post_value"),
      by = c("sample_id", "analyte", "draw_id")
    ) |>
      dplyr::mutate(
        value    = dplyr::coalesce(.data$.post_value, .data$value),
        detected = TRUE
      ) |>
      dplyr::select(-".post_value")

    result_targets <- dplyr::bind_rows(obs_expanded, imp_expanded)
  }

  # Recombine drivers (no imputation, no imputed columns)
  driver_df <- dplyr::mutate(
    driver_df,
    imputed      = FALSE,
    imputed_kind = "observed"
  )

  if (return == "draws") {
    driver_df <- tidyr::crossing(driver_df, draw_id = seq_len(max(post_long$draw_id)))
  }

  dplyr::bind_rows(result_targets, driver_df) |>
    dplyr::arrange(.data$sample_id, .data$analyte)
}

#' Post-hoc BDL check (Option B): warn/cap imputed values exceeding detection limit
#'
#' For BDL observations (`imputed_kind == "censored_left"`), if the imputed
#' posterior mean exceeds the original detection limit (stored in `dl_tbl`),
#' issue a warning. When `cap = TRUE`, the value is capped at the DL.
#'
#' This check is a minimum safeguard against the loss of the strict left-censor
#' constraint caused by using `mi()` instead of `cens()` (the latter is
#' incompatible with `rescor = TRUE` in brms). See the `impute_chemistry()`
#' documentation for details on Option D (custom Stan), which would handle this
#' correctly during fitting.
#'
#' @keywords internal
.check_bdl_imputed <- function(result, dl_tbl, cap = TRUE) {
  if (nrow(dl_tbl) == 0L) return(result)

  # Join detection limits onto BDL rows
  bdl_rows <- dplyr::filter(result, .data$imputed_kind == "censored_left") |>
    dplyr::left_join(dl_tbl, by = c("sample_id", "analyte"))

  if (nrow(bdl_rows) == 0L || !("detection_limit" %in% names(bdl_rows))) {
    return(result)
  }

  exceedances <- dplyr::filter(
    bdl_rows,
    !is.na(.data$detection_limit),
    .data$value > .data$detection_limit
  )

  if (nrow(exceedances) == 0L) return(result)

  n_ex <- nrow(exceedances)
  analytes_ex <- unique(exceedances$analyte)

  cli::cli_warn(c(
    "!" = "{n_ex} imputed BDL value{?s} exceed the original detection limit.",
    "i" = "Affected analyte{?s}: {.val {analytes_ex}}.",
    "i" = "This can occur because brms uses {.code mi()} (missing-indicator) \\
           instead of {.code cens('left')} for BDL observations when \\
           {.code rescor = TRUE} — the strict left-censor constraint is not \\
           enforced during MCMC.",
    if (cap) {
      "i" = "Values have been capped at the detection limit ({.arg bdl_cap = TRUE}). \\
             Set {.code bdl_cap = FALSE} to inspect uncapped values."
    } else {
      "i" = "Values have NOT been capped ({.arg bdl_cap = FALSE}). \\
             Use {.code bdl_cap = TRUE} to cap at the detection limit."
    }
  ))

  if (!cap) return(result)

  # Cap: replace imputed BDL values that exceed DL with the DL
  exceedance_keys <- dplyr::select(exceedances, "sample_id", "analyte",
                                    cap_value = "detection_limit")

  result <- dplyr::left_join(
    result, exceedance_keys, by = c("sample_id", "analyte")
  ) |>
    dplyr::mutate(
      value = dplyr::if_else(
        !is.na(.data$cap_value),
        .data$cap_value,
        .data$value
      )
    ) |>
    dplyr::select(-"cap_value")

  result
}
