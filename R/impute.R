#' Impute missing and below-detection-limit chemistry using a Bayesian
#' multivariate GAM
#'
#' Fits a multivariate log-normal GAM via `brms::brm()` with
#' `set_rescor(TRUE)` to capture cross-analyte covariance. Left-censored
#' observations (below-detection-limit, BDL) are modelled via `cens("left")`
#' in the brms likelihood — the latent concentration is inferred as somewhere
#' below the detection limit rather than assumed to be zero. Missing analyte
#' measurements at a sample are imputed from the cross-analyte covariance
#' structure and the driver predictors.
#'
#' After fitting, the function returns a long-format data frame with all BDL
#' and missing cells replaced by posterior means (or draws), plus an `imputed`
#' column to flag which values were filled in.
#'
#' **Installation note:** `brms` requires a working Stan toolchain (rstan or
#' cmdstanr). If Stan is not installed, this function will abort with a message
#' pointing to the installation documentation.
#'
#' @param df Long-format chemistry data frame. Required columns:
#'   - `uuid.sample` (character)
#'   - `uuid.feature` (character)
#'   - `datetime.sample` (Date or POSIXct)
#'   - `name.analyte` (character)
#'   - `value` (numeric) — for BDL rows, set to the detection limit
#'   - `quantified` (logical) — `FALSE` for BDL, `TRUE` for detected values
#'   Driver analytes (e.g. pH, EC, DOC) must be present as rows in `df` and
#'   must have `quantified == TRUE` for every sample.
#' @param drivers Character vector of analyte names to use as GAM predictors.
#'   Default `c("pH", "EC", "DOC")`. These analytes are extracted from `df`,
#'   pivoted to wide columns, and used as spline predictors. They are
#'   **not** themselves imputed.
#' @param formula_template `brms::brmsformula` or `NULL`. If `NULL`, a default
#'   formula is constructed: `mvbind(analyte_1, ...) ~ s(d1) + s(d2) + ...`
#'   with one smooth per driver.
#' @param family brms family for the response distribution. Default
#'   `"lognormal"` (appropriate for positive concentrations). Must be a family
#'   where the link is `log` so censoring on the log scale is handled
#'   correctly.
#' @param iter Total MCMC iterations per chain. Default `2000`.
#' @param warmup Warmup iterations per chain. Default `1000`.
#' @param chains Number of MCMC chains. Default `4`.
#' @param cores Number of parallel cores. Default `parallel::detectCores()`.
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
    family           = "lognormal",
    iter             = 2000,
    warmup           = 1000,
    chains           = 4,
    cores            = parallel::detectCores(),
    return           = c("point", "draws"),
    ...
) {
  return <- match.arg(return)

  # ── Dependency check ───────────────────────────────────────────────────────
  if (!requireNamespace("brms", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.pkg brms} is required for {.fn impute_chemistry} but is not installed.",
      "i" = "Install it with: {.code install.packages('brms')}",
      "i" = "brms also requires a Stan toolchain (rstan or cmdstanr).",
      "i" = "See {.url https://paul-buerkner.github.io/brms/} for instructions."
    ))
  }

  # ── Input validation ───────────────────────────────────────────────────────
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c("uuid.sample", "uuid.feature", "datetime.sample",
                     "name.analyte", "value", "quantified")
  )
  checkmate::assert_character(drivers, min.len = 1L, any.missing = FALSE)
  checkmate::assert_count(iter)
  checkmate::assert_count(warmup)
  checkmate::assert_count(chains)

  missing_drivers <- setdiff(drivers, unique(df$name.analyte))
  if (length(missing_drivers) > 0L) {
    cli::cli_abort(c(
      "Driver analyte{?s} not found in {.arg df}: {.val {missing_drivers}}.",
      "i" = "Drivers must be present as rows in {.arg df} with \\
             {.code name.analyte} matching the driver names."
    ))
  }

  # ── Separate targets from drivers ─────────────────────────────────────────
  driver_df  <- dplyr::filter(df, .data$name.analyte %in% .env$drivers)
  target_df  <- dplyr::filter(df, !(.data$name.analyte %in% .env$drivers))

  target_analytes <- unique(target_df$name.analyte)
  if (length(target_analytes) == 0L) {
    cli::cli_abort("No non-driver analytes found in {.arg df} to impute.")
  }

  # ── Check drivers are quantified everywhere ────────────────────────────────
  driver_missing <- driver_df |>
    dplyr::group_by(.data$uuid.sample, .data$name.analyte) |>
    dplyr::summarise(
      any_missing = any(!.data$quantified),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$any_missing)

  if (nrow(driver_missing) > 0L) {
    bad <- unique(driver_missing$name.analyte)
    cli::cli_abort(c(
      "Driver analyte{?s} have missing or BDL values: {.val {bad}}.",
      "i" = "Drivers must be fully quantified at every sample. \\
             Either remove BDL samples or choose different drivers."
    ))
  }

  # ── Pivot to wide format for brms ─────────────────────────────────────────
  # Drivers: one column each
  driver_wide <- driver_df |>
    dplyr::select("uuid.sample", "name.analyte", "value") |>
    tidyr::pivot_wider(
      names_from  = "name.analyte",
      values_from = "value",
      names_prefix = ".drv_"
    )

  driver_col_names <- paste0(".drv_", drivers)

  # Targets: value + censoring indicator per analyte
  target_wide_vals <- target_df |>
    dplyr::select("uuid.sample", "name.analyte", "value") |>
    tidyr::pivot_wider(
      names_from  = "name.analyte",
      values_from = "value"
    )

  # Censoring indicator: "left" if BDL (quantified == FALSE), "none" if detected
  # NA rows (analyte not measured at this sample) are left as NA in the wide df
  target_wide_cens <- target_df |>
    dplyr::mutate(
      .cens = dplyr::if_else(.data$quantified, "none", "left")
    ) |>
    dplyr::select("uuid.sample", "name.analyte", ".cens") |>
    tidyr::pivot_wider(
      names_from  = "name.analyte",
      values_from = ".cens",
      names_prefix = ".cens_"
    )

  wide_df <- driver_wide |>
    dplyr::left_join(target_wide_vals, by = "uuid.sample") |>
    dplyr::left_join(target_wide_cens, by = "uuid.sample")

  # ── Track imputation kind for each (sample, analyte) ─────────────────────
  impute_kind <- target_df |>
    dplyr::group_by(.data$uuid.sample, .data$name.analyte) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      .imputed_kind = dplyr::case_when(
        !.data$quantified ~ "censored_left",
        TRUE              ~ "observed"
      )
    ) |>
    dplyr::select("uuid.sample", "name.analyte", ".imputed_kind")

  # Samples with no measurement at all for an analyte: "missing"
  all_combos <- tidyr::expand_grid(
    uuid.sample  = unique(target_df$uuid.sample),
    name.analyte = target_analytes
  )
  impute_kind <- dplyr::left_join(all_combos, impute_kind,
                                  by = c("uuid.sample", "name.analyte")) |>
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

    # Rename analyte columns in wide_df to safe names
    old_names <- target_analytes
    new_names <- safe_analytes
    cens_old  <- paste0(".cens_", old_names)
    cens_new  <- paste0(".cens_", new_names)

    rename_map <- c(
      stats::setNames(new_names, old_names),
      stats::setNames(cens_new, cens_old)
    )
    wide_df <- dplyr::rename(wide_df, dplyr::any_of(rename_map))

    driver_terms <- paste0("s(", driver_col_names, ")", collapse = " + ")

    bf_list <- purrr::map(new_names, function(a) {
      cens_col <- paste0(".cens_", a)
      brms::bf(
        stats::as.formula(paste0(a, " | cens(", cens_col, ") ~ ", driver_terms))
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
      family  = brms::get_family(family),
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
    # Reshape to long format: (uuid.sample, name.analyte, draw_id, value)
    post_long <- .reshape_posterior_draws(
      post_draws, wide_df$uuid.sample, safe_analytes
    )

    result <- .build_imputed_df(df, post_long, impute_kind, drivers,
                                return = "draws")
    attr(result, "brmsfit") <- fit
    return(result)
  }

  # Point estimates: posterior mean
  post_means <- brms::fitted(fit, summary = TRUE, scale = "response")
  # post_means is [n_obs × (Estimate, Est.Error, ...) × response]
  # brms returns a 3D array with names in the third dim
  pm_long <- .reshape_posterior_means(
    post_means, wide_df$uuid.sample, safe_analytes
  )

  result <- .build_imputed_df(df, pm_long, impute_kind, drivers,
                              return = "point")
  attr(result, "brmsfit") <- fit
  result
}

# ── Internal reshape helpers ──────────────────────────────────────────────────

#' Reshape brms::fitted() output to long format
#' @keywords internal
.reshape_posterior_means <- function(post_means, sample_uids, safe_analytes) {
  # post_means is a 3-dim array [n_obs, n_stats, n_responses] or a matrix when
  # there's one response. Coerce to list of matrices, one per response.
  if (is.matrix(post_means)) {
    # single response (shouldn't happen in mvbind but handle gracefully)
    pm_list <- list(post_means)
    names(pm_list) <- names(safe_analytes)[1L]
  } else {
    pm_list <- lapply(dimnames(post_means)[[3L]], function(nm) post_means[, , nm])
    names(pm_list) <- dimnames(post_means)[[3L]]
  }

  purrr::map2_dfr(pm_list, names(pm_list), function(mat, safe_nm) {
    orig_nm <- names(safe_analytes)[safe_analytes == safe_nm]
    if (length(orig_nm) == 0L) orig_nm <- safe_nm  # fallback
    tibble::tibble(
      uuid.sample  = sample_uids,
      name.analyte = orig_nm,
      .post_mean   = mat[, "Estimate"]
    )
  })
}

#' Reshape brms::posterior_predict() output to long format
#' @keywords internal
.reshape_posterior_draws <- function(post_draws, sample_uids, safe_analytes) {
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
        uuid.sample  = rep(sample_uids, each = n_draws),
        name.analyte = orig_nm,
        draw_id      = rep(seq_len(n_draws), times = length(sample_uids)),
        .post_value  = as.vector(t(mat))
      )
    }
  )
}

#' Merge posterior estimates back into the original long-format df
#' @keywords internal
.build_imputed_df <- function(df, post_long, impute_kind, drivers,
                              return = "point") {
  target_df  <- dplyr::filter(df, !(.data$name.analyte %in% .env$drivers))
  driver_df  <- dplyr::filter(df, .data$name.analyte %in% .env$drivers)

  # Tag existing observed rows
  target_df <- dplyr::left_join(
    target_df,
    impute_kind,
    by = c("uuid.sample", "name.analyte")
  ) |>
    dplyr::mutate(
      imputed      = .data$.imputed_kind != "observed",
      imputed_kind = .data$.imputed_kind
    ) |>
    dplyr::select(-".imputed_kind")

  if (return == "point") {
    val_col <- ".post_mean"
    # Replace BDL / missing values with posterior means
    imputed_rows <- dplyr::filter(target_df, .data$imputed)
    observed_rows <- dplyr::filter(target_df, !.data$imputed)

    imputed_filled <- dplyr::left_join(
      imputed_rows,
      post_long |> dplyr::select("uuid.sample", "name.analyte", dplyr::all_of(val_col)),
      by = c("uuid.sample", "name.analyte")
    ) |>
      dplyr::mutate(
        value      = dplyr::coalesce(.data[[val_col]], .data$value),
        quantified = TRUE
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
        dplyr::select("uuid.sample", "name.analyte", "draw_id", ".post_value"),
      by = c("uuid.sample", "name.analyte", "draw_id")
    ) |>
      dplyr::mutate(
        value      = dplyr::coalesce(.data$.post_value, .data$value),
        quantified = TRUE
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
    dplyr::arrange(.data$uuid.sample, .data$name.analyte)
}
