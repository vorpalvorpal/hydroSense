# ── interpolate_loglinear() ───────────────────────────────────────────────────

#' Log-linear daily interpolation of grab samples, with a span cap
#'
#' Fills a grid of days from a point's own grab samples ("anchors"). Between two
#' anchors the value is interpolated linearly in `log(value)` (so the midpoint
#' of a span is the geometric mean of its ends), or linearly in `value` when
#' `log = FALSE`. This is the toxicant interpolation used by
#' [mspaf_pipeline()] with `ara_form = "added_hazard"`.
#'
#' @section Span cap:
#' A day is only filled when the anchors that bracket it are close enough
#' together. For a grid day `g`, with `x_i` the last anchor on or before `g`:
#' * there is no anchor on or before `g`: the day is uncovered (leading days
#'   are never back-filled);
#' * `g` is an anchor date: the anchor's own value and detection flag;
#' * `g` is after the last anchor: covered only while `g - x_last <=
#'   span_cap`, holding the last anchor's value and flag;
#' * otherwise, with `x_j` the next anchor: covered only when `x_j - x_i <=
#'   span_cap`, and then interpolated between the two anchors.
#'
#' Uncovered days have `value = NA` and `detected = NA`. `span_cap = Inf`
#' disables the cap, so every day after the first anchor is covered.
#'
#' @section Below-detection anchors:
#' A below-detection anchor (`detected = FALSE`) enters the interpolation at
#' the value supplied for it (the caller decides: the detection limit, half of
#' it, ...). An interpolated day is flagged `detected = FALSE` only when
#' **both** bracketing anchors are below detection; one detected anchor is
#' enough to make the day detected. An anchor date, or a held day after the
#' last anchor, carries that anchor's own flag.
#'
#' @param dates `Date` vector of anchor (sample) dates. Must be unique; any
#'   order.
#' @param values Numeric anchor values, one per date.
#' @param detected Logical, one per date. `FALSE` marks a below-detection
#'   anchor; `NA` is treated as `TRUE`.
#' @param grid `Date` vector of the days to fill, in any order.
#' @param span_cap Maximum anchor-to-anchor span (and hold after the last
#'   anchor) in days. Default 120. `Inf` disables the cap.
#' @param log Logical. Interpolate in `log(value)` (default) or in `value`.
#'   With `log = TRUE` anchors with an `NA` or non-positive value are dropped;
#'   with `log = FALSE` only `NA` anchors are dropped.
#'
#' @return A data frame with one row per element of `grid`, in grid order, and
#'   columns `date` (= `grid`), `value` (numeric) and `detected` (logical).
#'   With no usable anchor every row is `NA`.
#'
#' @seealso [mspaf_pipeline()] (`ara_form = "added_hazard"`)
#' @examples
#' d0 <- as.Date("2024-01-01")
#' interpolate_loglinear(
#'   dates = d0 + c(0, 10, 200), values = c(2, 8, 4),
#'   detected = c(TRUE, FALSE, TRUE), grid = d0 + c(-1, 0, 5, 10, 100, 330)
#' )
#' @export
interpolate_loglinear <- function(dates, values, detected, grid,
                                  span_cap = 120, log = TRUE) {
  if (!inherits(dates, "Date")) dates <- as.Date(dates)
  if (!inherits(grid, "Date")) grid <- as.Date(grid)
  checkmate::assert_numeric(values)
  checkmate::assert_logical(detected)
  checkmate::assert_number(span_cap, lower = 0)
  checkmate::assert_flag(log)
  n <- length(dates)
  if (length(values) != n || length(detected) != n) {
    cli::cli_abort(c(
      "{.arg dates}, {.arg values} and {.arg detected} must have the same \\
       length.",
      "x" = "Lengths: {n}, {length(values)} and {length(detected)}."
    ), class = "hydroSense_input")
  }
  if (anyDuplicated(dates[!is.na(dates)])) {
    dup <- unique(dates[duplicated(dates) & !is.na(dates)])
    cli::cli_abort(c(
      "{.arg dates} must be unique.",
      "x" = "Duplicated: {.val {format(dup)}}.",
      "i" = "Reduce the samples to one value per day first."
    ), class = "hydroSense_input")
  }
  .interp_loglinear_num(as.numeric(dates), as.numeric(values), detected,
                        as.numeric(grid), span_cap, log, grid)
}

## Numeric core of interpolate_loglinear(): `x` anchor days (numeric), `g` grid
## days (numeric). `grid_dates` is returned as the date column.
.interp_loglinear_num <- function(x, v, det, g, span_cap, log, grid_dates) {
  det[is.na(det)] <- TRUE
  ok <- !is.na(x) & !is.na(v) & (!log | v > 0)
  ng <- length(g)
  val <- rep(NA_real_, ng)
  dt  <- rep(NA, ng)
  if (any(ok)) {
    o <- order(x[ok])
    x <- x[ok][o]
    y <- v[ok][o]
    q <- det[ok][o]
    yl <- if (log) base::log(y) else y
    nx <- length(x)
    i <- findInterval(g, x)                # number of anchors <= g
    has <- i > 0L & !is.na(g)
    ii <- i[has]
    gg <- g[has]
    exact <- x[ii] == gg
    last  <- !exact & ii == nx
    mid   <- !exact & ii < nx
    v_out <- rep(NA_real_, length(ii))
    d_out <- rep(NA, length(ii))
    ## anchor dates
    v_out[exact] <- y[ii[exact]]
    d_out[exact] <- q[ii[exact]]
    ## held after the last anchor
    held <- last & (gg - x[nx]) <= span_cap
    v_out[held] <- y[nx]
    d_out[held] <- q[nx]
    ## between two anchors
    if (any(mid)) {
      im <- ii[mid]
      jm <- im + 1L
      cov <- (x[jm] - x[im]) <= span_cap
      f <- (gg[mid] - x[im]) / (x[jm] - x[im])
      vm <- (1 - f) * yl[im] + f * yl[jm]
      if (log) vm <- exp(vm)
      dm <- q[im] | q[jm]
      vm[!cov] <- NA_real_
      dm[!cov] <- NA
      v_out[mid] <- vm
      d_out[mid] <- dm
    }
    val[has] <- v_out
    dt[has]  <- d_out
  }
  data.frame(date = grid_dates, value = val, detected = dt)
}
