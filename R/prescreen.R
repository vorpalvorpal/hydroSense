#' Pre-screen analytes by detection frequency
#'
#' Computes detection frequency (proportion of samples with `quantified == TRUE`)
#' per analyte and returns the names that meet a minimum threshold. Use this
#' before `impute_chemistry()` to drop analytes that were almost never detected —
#' imputing such analytes from near-zero priors adds noise without ecological
#' signal.
#'
#' @param df Long-format chemistry data frame with at least columns
#'   `name.analyte` (character) and `quantified` (logical). Typically also
#'   contains `uuid.sample` and `uuid.feature`.
#' @param k Minimum detection frequency (proportion, 0–1). Analytes with
#'   `n_quantified / n_samples < k` are excluded. Default `0.05` (5 %).
#' @param group_by_feature Logical. If `TRUE`, detection frequency is computed
#'   per `uuid.feature` and an analyte is included only if it passes in *all*
#'   features. If `FALSE` (default), frequency is pooled across all samples.
#' @param return Either `"vector"` (default) to return a character vector of
#'   included analyte names, or `"table"` to return a tibble with one row per
#'   analyte showing detection statistics and inclusion flag.
#'
#' @return When `return = "vector"`: a named character vector of passing analyte
#'   names. The vector carries an attribute `"excluded"` listing analytes that
#'   did not pass the threshold, so callers can record what was dropped.
#'
#'   When `return = "table"`: a tibble with columns `name.analyte`,
#'   `n_samples`, `n_quantified`, `detect_freq`, `included`.
#'
#' @examples
#' \dontrun{
#' included <- prescreen_analytes(chemistry, k = 0.05)
#' attr(included, "excluded")  # see what was dropped
#' chem_f <- dplyr::filter(chemistry, name.analyte %in% included)
#' }
#'
#' @export
prescreen_analytes <- function(
    df,
    k = 0.05,
    group_by_feature = FALSE,
    return = c("vector", "table")
) {
  return <- match.arg(return)

  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df), must.include = c("name.analyte", "quantified"))
  checkmate::assert_number(k, lower = 0, upper = 1)
  checkmate::assert_flag(group_by_feature)

  if (group_by_feature) {
    checkmate::assert_names(names(df), must.include = "uuid.feature",
      .var.name = "df (when group_by_feature = TRUE)")
  }

  if (!is.logical(df$quantified)) {
    cli::cli_abort(
      "{.arg df$quantified} must be logical, not {.cls {class(df$quantified)}}."
    )
  }

  if (group_by_feature) {
    tbl <- df |>
      dplyr::group_by(name.analyte, uuid.feature) |>
      dplyr::summarise(
        n_samples    = dplyr::n(),
        n_quantified = sum(.data$quantified, na.rm = TRUE),
        .groups      = "drop"
      ) |>
      dplyr::group_by(name.analyte) |>
      dplyr::summarise(
        n_samples    = sum(.data$n_samples),
        n_quantified = sum(.data$n_quantified),
        detect_freq  = min(.data$n_quantified / .data$n_samples),  # worst-case feature
        .groups      = "drop"
      ) |>
      dplyr::mutate(included = .data$detect_freq >= k)
  } else {
    tbl <- df |>
      dplyr::group_by(name.analyte) |>
      dplyr::summarise(
        n_samples    = dplyr::n(),
        n_quantified = sum(.data$quantified, na.rm = TRUE),
        detect_freq  = .data$n_quantified / .data$n_samples,
        .groups      = "drop"
      ) |>
      dplyr::mutate(included = .data$detect_freq >= k)
  }

  n_excluded <- sum(!tbl$included)
  if (n_excluded > 0) {
    cli::cli_inform(c(
      "i" = "prescreen_analytes: {n_excluded} analyte{?s} below k={k} detection \\
             frequency threshold: {.val {tbl$name.analyte[!tbl$included]}}."
    ))
  }

  if (return == "table") {
    return(tbl)
  }

  included_names <- tbl$name.analyte[tbl$included]
  excluded_names <- tbl$name.analyte[!tbl$included]
  structure(included_names, excluded = excluded_names)
}
