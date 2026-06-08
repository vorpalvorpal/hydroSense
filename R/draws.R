# ── Draw-carrier primitives ───────────────────────────────────────────────────
#
# Draw-carrier contract (v2.0 uncertainty propagation):
#
#   A long frame where each (sample_id, analyte) cell is EITHER:
#     * one exact row  (draw_id = NA)  — observation or any fixed input, OR
#     * N rows keyed by draw_id 1..N   — a posterior / predictive draw.
#
#   draw_id column absent entirely → legacy / point frame; treated as all-exact.
#
#   Cross-analyte alignment: within a fitted imputation group, draw_id d for
#   analyte A and draw_id d for analyte B came from the SAME posterior sample.
#   Never permute draw order for cells that share a (sample_id, group).
#
#   Temporal alignment: draw_id pairing across time (chronic window) is an
#   index-pairing approximation that assumes temporal independence.
#   OU/Kalman smoothing (#16) is the future upgrade; leave a hook, not code.


#' Active draw domain from a draw-carrier frame
#'
#' Returns the sorted integer vector of draw IDs present in `df`, or
#' `integer(0)` for a deterministic (all-exact or legacy) frame.
#'
#' Asserts the domain is contiguous `1..N` and that every draw-bearing
#' `(sample_id, analyte)` cell carries the full set (no ragged N).
#'
#' @param df Long-format data frame possibly carrying a `draw_id` column.
#' @return Sorted integer vector `1..N`, or `integer(0)`.
#' @keywords internal
.draw_domain <- function(df) {
  if (!"draw_id" %in% names(df)) return(integer(0))
  ids <- sort(unique(na.omit(as.integer(df[["draw_id"]]))))
  if (length(ids) == 0L) return(integer(0))

  expected <- seq_len(max(ids))
  if (!identical(ids, expected)) {
    cli::cli_abort(c(
      "Draw IDs must form a contiguous sequence 1..N.",
      "x" = "Found: {ids}"
    ))
  }

  # Ragged-N check: every draw-bearing (sample_id, analyte) cell must carry
  # the full 1..N set — partial joins or mid-pipeline resampling produce ragged
  # N and break cross-analyte alignment.
  if (all(c("sample_id", "analyte") %in% names(df))) {
    draw_rows <- df[!is.na(df[["draw_id"]]), , drop = FALSE]
    if (nrow(draw_rows) > 0L) {
      cell_key <- paste(draw_rows[["sample_id"]], draw_rows[["analyte"]],
                        sep = "\x01")
      ukey <- unique(cell_key)
      counts <- tabulate(match(cell_key, ukey))
      bad_n  <- sort(unique(counts[counts != length(ids)]))
      if (length(bad_n) > 0L) {
        cli::cli_abort(c(
          "Ragged draw counts: every draw-bearing (sample_id, analyte) cell \\
           must have exactly {length(ids)} draw{?s}.",
          "x" = "Found cell{?s} with {bad_n} draw{?s}."
        ))
      }
    }
  }

  ids
}


#' Summarise a draw-carrier frame to posterior median and credible interval
#'
#' Collapses per-cell posterior draws to a point estimate and credible
#' interval.  Returns `df` unchanged (identity) when the input carries no
#' draws (point frame), preserving the degradation guarantee — no CI columns
#' appear for callers that never produce draws.
#'
#' **Per-draw diagnostic columns** (`draw_id`, `dominant_analyte`, `max_paf`,
#' `analyte_pafs`) are dropped in draws mode; they are available by passing
#' `return = "draws"` to [add_amspaf()] or [time_weighted_aggregate()].
#'
#' Exact cells (`draw_id = NA`) each collapse to a degenerate interval:
#' `value_lower = value_upper = value`, `n_draws = 1`.
#'
#' @param df Long-format data frame following the draw-carrier contract.
#' @param interval Width of the credible interval.  Default `0.90` yields
#'   5th / 95th percentile bounds.
#' @param central Central-tendency statistic: `"median"` (default) or
#'   `"mean"`.
#'
#' @return `df` collapsed to one row per cell, with columns `value`
#'   (central estimate), `value_lower`, `value_upper` (interval bounds), and
#'   `n_draws`.  Returns `df` unchanged when the input carries no draws.
#'
#' @seealso [add_amspaf()], [time_weighted_aggregate()]
#' @examples
#' \dontrun{
#' # Collapse AmsPAF draws to median + 90 % CI
#' amspaf_draws |>
#'   add_amspaf(return = "draws") |>
#'   summarise_draws(interval = 0.90)
#' }
#' @export
summarise_draws <- function(df, interval = 0.90,
                             central = c("median", "mean")) {
  central <- match.arg(central)
  checkmate::assert_number(interval, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_data_frame(df)

  # Identity for point frames — no CI columns, schema unchanged
  if (length(.draw_domain(df)) == 0L) return(df)

  # Cell-identity keys: whichever canonical columns are present
  key_cols <- intersect(
    c("focal_date", "site_id", "sample_id", "analyte", "datetime"),
    names(df)
  )

  # Scalar columns invariant within a cell (structural, not draw-varying)
  scalar_cols <- c(
    "detected", "imputed", "imputed_kind",
    "n_samples_in_window", "n_imputed_in_window",
    "n_analytes_used", "n_analytes_imputed"
  )

  alpha       <- (1 - interval) / 2
  central_fn  <- if (central == "median") stats::median else base::mean

  # Compute value_lower / value_upper BEFORE value so that .data$value
  # always refers to the original column, not the newly-created summary.
  # (dplyr can expose freshly-created columns to subsequent expressions
  # within the same summarise() call.)
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) |>
    dplyr::summarise(
      value_lower = stats::quantile(.data$value, alpha,     na.rm = TRUE,
                                    names = FALSE),
      value_upper = stats::quantile(.data$value, 1 - alpha, na.rm = TRUE,
                                    names = FALSE),
      n_draws     = dplyr::n(),
      value       = central_fn(.data$value, na.rm = TRUE),
      dplyr::across(dplyr::any_of(scalar_cols), dplyr::first),
      .groups = "drop"
    )
}


#' Broadcast a draw-carrier frame to uniform draw coverage
#'
#' After this call every row has a concrete integer `draw_id`; downstream code
#' can `group_by(draw_id)` with no `NA`-special-casing.
#'
#' Exact cells (`draw_id = NA`) are replicated once per draw in `draws` with
#' `draw_id` filled.  Draw cells pass through unchanged, order preserved.
#'
#' **Call inside a per-site or per-sample scope** — not on the whole dataset —
#' because replication multiplies exact rows by N.
#'
#' @param df Long-format data frame following the draw-carrier contract.
#' @param draws Integer vector of active draw IDs (default: `.draw_domain(df)`).
#' @return `df` with `draw_id` filled for all rows; exact rows replicated N times.
#' @keywords internal
.broadcast_draws <- function(df, draws = .draw_domain(df)) {
  # Deterministic: no draws in frame → single pass, draw_id = 1L everywhere
  if (length(draws) == 0L) {
    df[["draw_id"]] <- 1L
    return(df)
  }

  if (!"draw_id" %in% names(df)) df[["draw_id"]] <- NA_integer_

  # Already fully broadcast (no NAs) → idempotent
  if (!anyNA(df[["draw_id"]])) return(df)

  exact_mask <- is.na(df[["draw_id"]])
  exact_rows <- df[ exact_mask, , drop = FALSE]
  draw_rows  <- df[!exact_mask, , drop = FALSE]

  if (nrow(exact_rows) == 0L) return(df)

  n_exact <- nrow(exact_rows)
  n_draws <- length(draws)

  # For each draw d: replicate every exact row with draw_id = d.
  # rep(1:n_exact, times=n_draws)  → 1,2,..,n_exact, 1,2,..,n_exact, ...
  # rep(draws,    each=n_exact)    → d1,d1,..,d2,d2,..
  idx       <- rep(seq_len(n_exact), times = n_draws)
  broadcast <- exact_rows[idx, , drop = FALSE]
  broadcast[["draw_id"]] <- rep(draws, each = n_exact)
  rownames(broadcast) <- NULL

  dplyr::bind_rows(draw_rows, broadcast)
}
