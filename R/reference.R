#' Prepare reference chemistry for AmsPAF background subtraction
#'
#' Applies chemistry normalisation (if formulas are populated in the metadata)
#' and computes a per-analyte central-tendency summary from reference-site
#' data, optionally with a bootstrap confidence interval.  The resulting
#' object is passed as the `reference` argument to `add_amspaf()`.
#'
#' This is a pure function — it has no side effects and no internal cache.
#' In the chronic pipeline, call it once after computing chronic chemistry for
#' the reference feature(s), then pass the object into `add_amspaf()` for
#' every focal date.  In the per-sample pipeline, call it once on the raw
#' reference chemistry.
#'
#' **Summary statistic.**  The default is `"geom_mean"` — the geometric mean
#' of all detected observations.  This is preferred over a fixed quantile
#' because:
#' \itemize{
#'   \item it is the maximum-likelihood central tendency for log-normal
#'     concentrations (which is how aquatic concentration data typically
#'     distribute);
#'   \item it uses all observations rather than a single ranked point, so
#'     it is more robust to small reference datasets;
#'   \item it is PICT-consistent: the resident community has adapted to the
#'     integrated typical exposure over time, not to a particular upper
#'     quantile.
#' }
#' BDL observations contribute `0` to the geometric mean via an
#' \eqn{\epsilon}-shifted log: `exp(mean(log(value + eps)))`.  Other
#' summaries available: `"median"`, `"arith_mean"`, `"p80"`, `"p90"`,
#' `"p95"`.
#'
#' @param reference_data Long-format chemistry data frame for the reference
#'   (background) site(s). Same schema as the input to `add_amspaf()`:
#'   `analyte`, `value`, `detected`. Toxicant concentrations must be in µg/L
#'   before normalisation; supply them either via a `units.analyte` column or
#'   via the `conc_units` argument. BDL (`detected == FALSE`) observations
#'   contribute `0` to the summary statistic.
#' @param analyte_metadata Data frame of analyte metadata, or `NULL` to load
#'   the bundled `inst/extdata/anzecc_analyte_metadata.csv`. Accepts either a
#'   data frame or a file path string. Must contain columns `analyte`,
#'   `coanalytes_required`, and `normalisation_formula`.
#' @param summary Summary statistic for the reference distribution. One of
#'   `"geom_mean"` (default), `"median"`, `"arith_mean"`, `"p80"`, `"p90"`,
#'   `"p95"`.
#' @param bootstrap_ci Logical.  If `TRUE`, compute a 95% bootstrap CI on the
#'   reference summary statistic for each analyte (1,000 replicates by
#'   default).  Adds `ref_lower`, `ref_upper`, and `n_boot_valid` columns to
#'   `$ref_table`.  Default `FALSE`.
#' @param n_boot Number of bootstrap replicates if `bootstrap_ci = TRUE`.
#'   Default `1000L`.
#' @param conc_units Character. Unit string (e.g. `"mg/L"`, `"ug/L"`) applied
#'   uniformly to all SSD-eligible rows in `reference_data` when it has no
#'   `units.analyte` column. Ignored when `reference_data` carries
#'   `units.analyte`. Required when the data lacks `units.analyte` and toxicant
#'   concentrations are not already in µg/L.
#' @param eps Small positive guard added inside the log for geometric-mean
#'   computation, to handle BDL contributions of `0`.  Default `1e-9`.
#'
#' @return A list of class `"prepared_reference"` with elements:
#'   \itemize{
#'     \item `$ref_table`: tibble with columns `analyte`, `ref_norm`
#'       (normalised summary concentration), `n_obs` (count of observations
#'       contributing).  If `bootstrap_ci = TRUE`, also `ref_lower` /
#'       `ref_upper` (95% CI) and `n_boot_valid` (number of bootstrap draws
#'       that yielded a finite summary; CIs built from far fewer than `n_boot`
#'       draws are flagged with a warning).
#'     \item `$dropped`: character vector of analytes excluded due to no
#'       reference observations.
#'     \item `$summary`: the summary statistic used.
#'   }
#'
#' @examples
#' ref <- subset(leachate_demo(), site_id == "reference")
#'
#' # Default: geometric mean (recommended)
#' prep_ref <- prepare_reference(ref)
#' prep_ref$ref_table
#'
#' # Or a higher percentile of the local background distribution:
#' prepare_reference(ref, summary = "p80")$ref_table
#'
#' @export
prepare_reference <- function(
    reference_data,
    analyte_metadata = NULL,
    summary          = c("geom_mean", "median", "arith_mean",
                          "p80", "p90", "p95"),
    bootstrap_ci     = FALSE,
    n_boot           = 1000L,
    conc_units       = NULL,
    eps              = 1e-9
) {
  checkmate::assert_data_frame(reference_data)
  checkmate::assert_names(names(reference_data),
    must.include = c("analyte", "value", "detected"))
  checkmate::assert_flag(bootstrap_ci)
  checkmate::assert_int(n_boot, lower = 100L)
  checkmate::assert_number(eps, lower = 0)
  summary <- match.arg(summary)

  meta <- .load_analyte_metadata(analyte_metadata)

  ## Convert toxicant concentrations to µg/L before normalisation.
  ## Identify SSD-eligible analytes the same way derive_ssd_params() does so
  ## that co-analyte rows (pH, DOC, hardness, etc.) are left untouched.
  ssd_analytes <- meta$analyte[
    !is.na(meta$ssd_available) & meta$ssd_available == TRUE &
    !meta$analyte %in% .AMSPAF_EXCLUDED_ANALYTES
  ]
  reference_data <- .convert_df_tox_to_ugL(
    reference_data, ssd_analytes, conc_units, "reference_data"
  )

  # BDL reference observations are treated as 0 (not excluded).
  # Rationale: the summary represents what the local biota is adapted to.
  # If 99/100 reference samples are BDL, the background is genuinely near
  # zero — excluding BDLs would inflate the reference and understate risk.
  # Detected rows keep their measured value; BDL rows contribute 0.
  # Normalisation is applied only to detected rows (normalising 0 is a
  # no-op, but avoids formula edge cases).
  ref_q <- dplyr::mutate(
    reference_data,
    value = dplyr::if_else(.data$detected, .data$value, 0)
  )

  empty_result <- function() structure(
    list(
      ref_table = tibble::tibble(
        analyte  = character(0),
        ref_norm = numeric(0),
        n_obs    = integer(0)
      ),
      dropped = character(0),
      summary = summary
    ),
    class = "prepared_reference"
  )

  if (nrow(ref_q) == 0L) {
    cli::cli_warn(
      "No reference observations found; ARA subtraction will be zero for all analytes."
    )
    return(empty_result())
  }

  # Apply normalisation per analyte row (detected rows only; BDL rows stay 0)
  nq <- .normalise_ref_observations(ref_q, reference_data, meta)

  # Per-analyte summary statistic
  qnt <- nq |>
    dplyr::group_by(.data$analyte) |>
    dplyr::summarise(
      ref_norm = .ref_summary(.data$value_norm, summary, eps),
      n_obs    = sum(!is.na(.data$value_norm)),
      .groups  = "drop"
    )

  # Optional bootstrap CI per analyte
  if (bootstrap_ci) {
    ci_tbl <- nq |>
      dplyr::group_by(.data$analyte) |>
      dplyr::summarise(
        .ci = list(
          .ref_summary_bootstrap_ci(.data$value_norm, summary, eps, n_boot)
        ),
        .groups = "drop"
      ) |>
      tidyr::unnest_wider(".ci")

    qnt <- dplyr::left_join(qnt, ci_tbl, by = "analyte")

    ## Flag analytes where a large share of bootstrap draws were lost to
    ## non-finite summaries — the CI is built from too few effective draws to
    ## be trustworthy.
    if ("n_boot_valid" %in% names(qnt)) {
      shaky <- qnt$analyte[!is.na(qnt$n_boot_valid) &
                             qnt$n_boot_valid < 0.9 * n_boot]
      if (length(shaky) > 0L) {
        cli::cli_warn(c(
          "!" = "{length(shaky)} analyte{?s} lost > 10% of bootstrap draws to \\
                 non-finite values \u2014 CI may be unreliable: {.val {shaky}}.",
          "i" = "Inspect with the `n_boot_valid` column of `ref_table`."
        ))
      }
    }
  }

  # Analytes that ended up with no usable reference observations:
  # either never detected at the reference site, or normalisation returned
  # NA for every row (e.g. required co-analyte absent from reference data).
  all_analytes <- unique(reference_data$analyte)
  dropped      <- setdiff(all_analytes, qnt$analyte[qnt$n_obs > 0L])

  if (length(dropped) > 0L) {
    cli::cli_inform(c(
      "i" = "prepare_reference: {length(dropped)} analyte{?s} dropped \u2014 \\
             no usable normalised reference values: {.val {dropped}}."
    ))
  }

  # Warn about analytes with very few observations (CI may be unreliable)
  low_n <- qnt$analyte[qnt$n_obs > 0L & qnt$n_obs < 5L]
  if (length(low_n) > 0L) {
    cli::cli_warn(c(
      "!" = "{length(low_n)} analyte{?s} have < 5 reference observations \\
             \u2014 `ref_norm` estimate may be unreliable: {.val {low_n}}."
    ))
  }

  structure(
    list(
      ref_table = dplyr::filter(qnt, .data$n_obs > 0L),
      dropped   = dropped,
      summary   = summary
    ),
    class = "prepared_reference"
  )
}

# ── Shared normalisation helper ───────────────────────────────────────────────

#' Normalise reference observations to the SSD index condition
#'
#' Applies per-analyte bioavailability / physicochemical normalisation formulas
#' (hardness, pH, DOC) to a long-format chemistry data frame, returning the
#' frame with a `value_norm` column added.  BDL rows (`detected == FALSE`) are
#' normalised at their detection-limit value but callers are responsible for
#' setting `value = 0` for BDL rows *before* calling this if they want BDL to
#' contribute zero to downstream summaries (as `prepare_reference()` does).
#'
#' @param ref_df Long-format chemistry data frame with at minimum `analyte`,
#'   `value`, `detected`, and `sample_id`.
#' @param original_data The original reference chemistry frame used to look up
#'   co-analyte values (pH, DOC, hardness, Ca, Mg) per sample.
#' @param meta Analyte metadata tibble from `.load_analyte_metadata()`.
#' @return `ref_df` with a `value_norm` numeric column appended (NA where
#'   normalisation fails due to a missing required co-analyte).
#' @keywords internal
.normalise_ref_observations <- function(ref_df, original_data, meta) {
  has_sample_id <- "sample_id" %in% names(ref_df)

  ref_df |>
    dplyr::left_join(
      meta |>
        dplyr::select("analyte", "coanalytes_required", "normalisation_formula"),
      by = "analyte"
    ) |>
    dplyr::mutate(
      value_norm = purrr::pmap_dbl(
        list(
          det           = .data$detected,
          formula_str   = .data$normalisation_formula,
          C             = .data$value,
          sample_id_val = if (has_sample_id) .data$sample_id else NA_character_,
          co_req        = .data$coanalytes_required
        ),
        function(det, formula_str, C, sample_id_val, co_req) {
          if (!det) return(C)   # BDL: normalise at DL (caller sets value=0 for BDL→0 summary)
          parsed <- .parse_normalisation_formula(formula_str %||% "")
          if (is.null(parsed)) return(C)  # identity (no formula)
          coanalytes <- .extract_coanalytes_for_sample(
            original_data, sample_id_val, co_req %||% ""
          )
          .apply_normalisation(parsed, C, coanalytes)
        }
      )
    )
}


# ── Summary-statistic helpers ─────────────────────────────────────────────────

#' Compute a reference-distribution summary statistic
#' @keywords internal
.ref_summary <- function(x, summary, eps = 1e-9) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  switch(summary,
    geom_mean  = exp(mean(log(x + eps))) - eps,
    arith_mean = mean(x),
    median     = stats::median(x),
    p80        = stats::quantile(x, probs = 0.80, names = FALSE),
    p90        = stats::quantile(x, probs = 0.90, names = FALSE),
    p95        = stats::quantile(x, probs = 0.95, names = FALSE),
    cli::cli_abort("Unknown summary {.val {summary}}.")
  )
}

#' Bootstrap CI for a reference summary statistic
#' @keywords internal
.ref_summary_bootstrap_ci <- function(x, summary, eps = 1e-9, n_boot = 1000L) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3L) {
    return(list(ref_lower = NA_real_, ref_upper = NA_real_, n_boot_valid = 0L))
  }
  draws <- vapply(seq_len(n_boot), function(.) {
    .ref_summary(sample(x, size = n, replace = TRUE), summary, eps)
  }, numeric(1))
  ## Some draws can be non-finite (e.g. a resample that summarises to NA);
  ## quantile() drops them via na.rm but the caller needs to know how many
  ## draws actually contributed so it can flag unreliable CIs.
  list(
    ref_lower    = unname(stats::quantile(draws, 0.025, na.rm = TRUE)),
    ref_upper    = unname(stats::quantile(draws, 0.975, na.rm = TRUE)),
    n_boot_valid = sum(is.finite(draws))
  )
}

#' @export
print.prepared_reference <- function(x, ...) {
  has_ci <- "ref_lower" %in% names(x$ref_table)
  cat(sprintf(
    "<prepared_reference>  %d analytes | summary = %s%s | %d dropped\n",
    nrow(x$ref_table),
    x$summary,
    if (has_ci) " | bootstrap CI" else "",
    length(x$dropped)
  ))
  invisible(x)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Load and lightly validate the analyte metadata table
#'
#' Reads the bundled CSV when `meta` is NULL. Returns a tibble. Caches the
#' result in an internal environment to avoid re-reading on every
#' `prepare_reference()` / `add_amspaf()` call within a session.
#' @keywords internal
.meta_cache_env <- new.env(parent = emptyenv())

.load_analyte_metadata <- function(meta = NULL) {
  # Accept: NULL (load bundled CSV), a data frame, or a file path string.
  if (is.character(meta)) {
    if (length(meta) != 1L || !nzchar(meta)) {
      cli::cli_abort("{.arg analyte_metadata} must be a single file path, not {.val {meta}}.")
    }
    if (!file.exists(meta)) {
      cli::cli_abort("File not found: {.path {meta}}")
    }
    m <- readr::read_csv(meta, show_col_types = FALSE) |>
      dplyr::mutate(
        coanalytes_required  = dplyr::if_else(
          is.na(.data$coanalytes_required), "", .data$coanalytes_required),
        normalisation_formula = dplyr::if_else(
          is.na(.data$normalisation_formula), "", .data$normalisation_formula)
      )
    checkmate::assert_names(
      names(m),
      must.include = c("analyte", "coanalytes_required", "normalisation_formula")
    )
    return(m)
  }

  if (!is.null(meta)) {
    checkmate::assert_data_frame(meta)
    checkmate::assert_names(
      names(meta),
      must.include = c("analyte", "coanalytes_required", "normalisation_formula")
    )
    return(meta)
  }

  if (exists("meta", envir = .meta_cache_env, inherits = FALSE)) {
    return(get("meta", envir = .meta_cache_env, inherits = FALSE))
  }

  path <- system.file("extdata", "anzecc_analyte_metadata.csv",
                      package = "hydroSense")
  if (!nzchar(path)) {
    cli::cli_abort(
      "Cannot find {.file inst/extdata/anzecc_analyte_metadata.csv} inside the \\
       installed hydroSense package. Re-install or supply {.arg analyte_metadata} \\
       explicitly."
    )
  }

  m <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(
      coanalytes_required  = dplyr::if_else(
        is.na(.data$coanalytes_required), "", .data$coanalytes_required),
      normalisation_formula = dplyr::if_else(
        is.na(.data$normalisation_formula), "", .data$normalisation_formula)
    )

  assign("meta", m, envir = .meta_cache_env)
  m
}

#' Extract co-analyte values for a given sample from a long-format df
#'
#' @param df Long-format chemistry df
#' @param sample_id Sample identifier (may be NA for non-sample data)
#' @param coanalytes_str Comma-separated string of required co-analyte names
#' @return Named numeric vector; may be empty
#' @keywords internal
.extract_coanalytes_for_sample <- function(df, sample_id, coanalytes_str) {
  if (is.na(sample_id) || !nzchar(coanalytes_str %||% "")) {
    return(numeric(0))
  }
  required <- trimws(strsplit(coanalytes_str, ",")[[1L]])
  required <- required[nzchar(required)]
  if (length(required) == 0L) return(numeric(0))

  co <- dplyr::filter(df,
    .data$sample_id == .env$sample_id,
    .data$analyte %in% .env$required,
    .data$detected
  )
  if (nrow(co) == 0L) return(numeric(0))
  vals <- co$value
  names(vals) <- co$analyte
  vals[required[required %in% names(vals)]]
}

# Null-coalescing operator (base R doesn't have one before 4.4)
`%||%` <- function(x, y) if (is.null(x)) y else x
