#' Time-weighted chronic aggregation for any long-format value column
#'
#' For each (focal date × monitoring feature × analyte) combination, aggregates
#' values from the preceding `window` (in days) using exponential-decay temporal
#' weighting and forward-step duration weighting.  This is the chronic
#' exposure / response predictor used for downstream calibration against
#' biological community state.
#'
#' This function is value-agnostic: pass any long-format data frame with
#' `analyte` and `value` columns and it will compute one time-weighted value
#' per (focal_date × site_id × analyte).  Use it to aggregate raw chemistry
#' (`analyte` = chemical species, `summary = "geom_mean"`) or per-sample
#' msPAF values (`analyte = "msPAF"`, `summary = "arith_mean"`) into a
#' chronic predictor.
#'
#' **Forward-step duration weighting** treats each sample as representing the
#' period from its collection date to the next sample's date (or to
#' `focal_date` for the most recent sample).  This corrects for pulse-biased
#' sampling where storm events are sampled more frequently than base-flow
#' periods, which would otherwise over-weight episodic concentrations.
#'
#' **Exponential-decay temporal weighting** assigns higher weight to recent
#' samples, with a half-life of approximately `tau * log(2)` days (at the
#' default of 90 d).
#'
#' These two components are the only weighting scheme offered, by design:
#' forward-step duration weighting is the minimal correction for irregular /
#' pulse-biased sampling, and exponential decay is the standard memory kernel
#' for a community integrating a fluctuating exposure (one interpretable
#' parameter, `tau`, tied to the target biology's response time).  Richer
#' kernels would add parameters without a defensible way to fit them from
#' routine monitoring data, so they are left out — tune `tau` rather than
#' swapping the kernel.
#'
#' **Combined weight** for sample *i* is
#' `w_i = Δt_i × exp(-(focal_date - midpoint_i) / tau_d)`, where
#' `midpoint_i` is the midpoint of sample *i*'s representative interval.
#'
#' **Choosing `summary`:**
#'
#' - `"geom_mean"`: max-likelihood central tendency for log-normal data.
#'   Use for chemistry concentrations.  Caveat: with highly pulsed exposure
#'   and a non-linear SSD response, this underestimates the time-averaged
#'   risk metric (Jensen's inequality on the upper tail of the SSD).
#' - `"arith_mean"`: weighted arithmetic mean.  Use for bounded indices
#'   like msPAF percentages, or for chemistry when comparison against an
#'   arithmetic-mean compliance trigger is wanted.
#' - `"p90"`: 90th percentile (duration-weighted empirical CDF).  Diagnostic
#'   for "what's the upper-end exposure the community sees most of the
#'   time."
#'
#' For chronic msPAF specifically, the recommended pipeline is
#' `add_mspaf()` on per-sample chemistry, then `time_weighted_aggregate()`
#' on the resulting msPAF rows with `summary = "arith_mean"`.  This
#' computes the time-averaged msPAF, which integrates the toxic response
#' over time (consistent with how biological communities respond to
#' fluctuating exposure) — rather than the msPAF computed at a single
#' time-averaged chemistry, which can substantially under-state risk for
#' pulsed exposures.
#'
#' @param df Long-format data frame.  Required columns:
#'   - `sample_id` (character) — unique sample identifier
#'   - `site_id` (character) — monitoring feature identifier
#'   - `datetime` (Date or POSIXct) — sample collection date/time
#'   - `analyte` (character) — analyte / index name
#'   - `value` (numeric) — concentration or index value
#'   - `detected` (logical) — whether value is a detected observation
#'   Optional: `imputed` (logical) — propagated to `n_imputed_in_window`.
#' @param focal_dates Either a `Date` vector (applied to all features in `df`)
#'   or a data frame with columns `focal_date` (Date) and `site_id`
#'   (character) specifying per-feature focal dates.
#' @param tau Exponential-decay rate parameter.  Numeric or `units` object;
#'   bare numeric requires `tau_units`.  Default `NULL` → 90 days.  The
#'   effective half-life is `tau * log(2)` ≈ 62 days at the default.  Choose
#'   to match the response timescale of the downstream biology (algae:
#'   days–weeks; macroinvertebrates: weeks–months; fish: months–years).
#' @param tau_units Character unit string for `tau` when it is bare numeric,
#'   e.g. `"d"`.  Ignored when `tau` is a `units` object or `NULL`.
#' @param window Look-back window.  Numeric or `units` object; bare numeric
#'   requires `window_units`.  Default `NULL` → 365 days.
#' @param window_units Character unit string for `window` when it is bare
#'   numeric, e.g. `"d"`.  Ignored when `window` is a `units` object or
#'   `NULL`.
#' @param summary Aggregation method: `"geom_mean"` (default), `"arith_mean"`,
#'   `"p90"`.  See *Choosing `summary`* above.
#' @param anchor_outside_window Logical (default `TRUE`).  If `TRUE`, the most
#'   recent sample *before* the window start is included as an anchor to
#'   provide duration weighting at the leading edge of the window. Its
#'   duration is clipped to start at `window_start` so it only contributes
#'   the in-window portion of its interval.
#' @param eps Small positive guard added inside the log for geometric mean
#'   to avoid `log(0)`. Default `1e-9`.
#' @param return Output mode for draw-carrier input (see [summarise_draws()]).
#'   `"summary"` (default) collapses posterior draws to a central estimate plus
#'   a credible interval (`value`, `value_lower`, `value_upper`, `n_draws`);
#'   `"draws"` returns the raw per-draw chronic rows (`draw_id 1..N`).  For
#'   point (non-draw) input both modes return byte-identical output with no
#'   interval columns.  Draws are paired across time by `draw_id` (an
#'   index-pairing approximation that assumes temporal independence; see the
#'   OU/Kalman smoothing roadmap item).
#' @param interval Width of the credible interval when `return = "summary"`.
#'   Default `0.90` (5th/95th percentile bounds).
#' @param central Central-tendency statistic when `return = "summary"`:
#'   `"median"` (default) or `"mean"`.
#'
#' @return A long-format tibble with columns:
#'   - `focal_date` (Date)
#'   - `site_id` (character)
#'   - `sample_id` (character) — synthetic key `"chronic_<focal_date>_<site>"`
#'   - `analyte` (character)
#'   - `value` (numeric) — chronic time-weighted value
#'   - `detected` (logical) — always `TRUE`
#'   - `n_samples_in_window` (integer)
#'   - `n_imputed_in_window` (integer) — count of `imputed == TRUE` samples
#'     contributing; 0 if `imputed` column absent from `df`
#'
#' @examples
#' # Chronic (time-weighted geometric-mean) chemistry for one downstream
#' # analyte at two focal dates, from the bundled demo data.
#' cu <- subset(leachate_demo(), site_id == "downstream" & analyte == "Cu")
#' time_weighted_aggregate(
#'   cu,
#'   focal_dates = as.Date(c("2024-06-01", "2024-12-01")),
#'   tau = 90, tau_units = "d", window = 365, window_units = "d",
#'   summary = "geom_mean"
#' )
#'
#' @export
time_weighted_aggregate <- function(
    df,
    focal_dates,
    tau         = NULL,
    tau_units   = NULL,
    window      = NULL,
    window_units = NULL,
    summary     = c("geom_mean", "arith_mean", "p90"),
    anchor_outside_window = TRUE,
    eps         = 1e-9,
    return      = c("summary", "draws"),
    interval    = 0.90,
    central     = c("median", "mean")
) {
  summary <- match.arg(summary)
  return  <- match.arg(return)
  central <- match.arg(central)

  # ── Input validation ───────────────────────────────────────────────────────
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected")
  )
  tau_days    <- if (is.null(tau))    90  else .resolve_to(tau,    "d", tau_units,    "tau")
  window_days <- if (is.null(window)) 365 else .resolve_to(window, "d", window_units, "window")
  checkmate::assert_number(tau_days,    lower = 0.001)
  checkmate::assert_number(window_days, lower = 1)
  checkmate::assert_flag(anchor_outside_window)
  checkmate::assert_number(eps, lower = 0)

  has_imputed   <- "imputed" %in% names(df)
  is_draws_mode <- "draw_id" %in% names(df) && !all(is.na(df[["draw_id"]]))
  draws         <- .draw_domain(df)

  # Normalise datetime to Date for day-level arithmetic
  df <- dplyr::mutate(df,
    .date = as.Date(.data$datetime)
  )

  # ── Build the focal_dates × site_id grid ──────────────────────────────────
  if (is.data.frame(focal_dates)) {
    checkmate::assert_names(names(focal_dates),
      must.include = c("focal_date", "site_id"))
    focal_grid <- dplyr::mutate(focal_dates,
      focal_date = as.Date(.data$focal_date)
    )
  } else {
    focal_dates <- as.Date(focal_dates)
    features    <- unique(df$site_id)
    focal_grid  <- tidyr::expand_grid(
      focal_date = focal_dates,
      site_id    = features
    )
  }

  # ── Chronic values for every (focal_date, site_id) at once ────────────────
  # Vectorised over focal dates (see .twa_compute()): per point, the forward-
  # step intervals of the whole sample schedule are computed once and every
  # focal window is resolved from prefix sums, instead of re-filtering and
  # re-summarising the point's samples once per focal date.
  res <- .twa_compute(
    df, focal_grid,
    tau_days              = tau_days,
    window_days           = window_days,
    summary               = summary,
    anchor_outside_window = anchor_outside_window,
    eps                   = eps,
    has_imputed           = has_imputed,
    draws                 = draws
  )

  if (is.null(res)) {
    cli::cli_abort(c(
      "No values could be aggregated for any (focal_date, site_id).",
      "i" = "Check that {.arg focal_dates} and {.arg df$datetime} overlap \\
             within window = {window_days} days, and that \\
             {.arg df$site_id} matches the sites in {.arg focal_dates}."
    ))
  }

  # One row per (focal_date, site_id, analyte, draw): focal-grid row order,
  # then analyte (C locale) and draw_id within each row, as group_by() sorts.
  focal_out <- focal_grid$focal_date[res$r]
  site_out  <- focal_grid$site_id[res$r]
  out <- tibble::tibble(
    focal_date          = focal_out,
    site_id             = site_out,
    sample_id           = paste0("chronic_", focal_out, "_", site_out),
    analyte             = res$analyte,
    value               = res$value,
    detected            = TRUE,
    n_samples_in_window = res$n_in,
    n_imputed_in_window = res$n_imp
  )
  if (is_draws_mode) {
    # Same type as a broadcast draw_id column would carry.
    out$draw_id <- if (is.double(df$draw_id)) {
      as.double(res$draw)
    } else {
      as.integer(res$draw)
    }
  }

  if (return == "summary") out <- summarise_draws(out, interval, central)
  out
}

# ── expand_focal_dates ────────────────────────────────────────────────────────

#' Generate a sequence of focal dates for chronic msPAF computation
#'
#' A thin convenience wrapper around [base::seq.Date()] for generating the
#' `focal_dates` vector passed to [time_weighted_aggregate()].  The most
#' common use is a daily sequence for time-series analysis.
#'
#' @param start Start date (character `"YYYY-MM-DD"` or `Date`).
#' @param end   End date (character `"YYYY-MM-DD"` or `Date`).
#' @param by    Increment. Passed to [base::seq.Date()]. Common values:
#'   `"day"`, `"week"`, `"month"`. Default `"day"`.
#'
#' @return A `Date` vector from `start` to `end` at the specified increment.
#'
#' @examples
#' # Monthly focal dates across 2024, e.g. to feed time_weighted_aggregate().
#' expand_focal_dates("2024-01-01", "2024-12-31", by = "month")
#'
#' @export
expand_focal_dates <- function(start, end, by = "day") {
  start <- as.Date(start)
  end   <- as.Date(end)
  if (is.na(start)) cli::cli_abort("{.arg start} could not be parsed as a date.")
  if (is.na(end))   cli::cli_abort("{.arg end} could not be parsed as a date.")
  if (end < start)  cli::cli_abort("{.arg end} must be on or after {.arg start}.")
  seq.Date(start, end, by = by)
}


# ── Internal: vectorised chronic aggregation ──────────────────────────────────
#
# For one point, sort its sample schedule (one entry per sample_id x date) by
# date: d_1 <= ... <= d_n.  A focal date f with window start ws = f - window
# uses
#   * the in-window samples lo..hi (ws <= d_k <= f);
#   * the anchor lo - 1 (the latest sample strictly before ws), if requested.
# Sample k in lo..hi-1 represents [d_k, d_{k+1}] -- its forward-step interval
# over the WHOLE schedule, whatever the analyte -- so its weight
#   w_k = (d_{k+1} - d_k) * exp(-(f - m_k) / tau),  m_k the interval midpoint,
# factors as exp(-f / tau) * [(d_{k+1} - d_k) * exp(m_k / tau)], and the window
# sums of w and w * x over lo..hi-1 are differences of prefix sums.  Only the
# last in-window sample ([d_hi, f]) and the anchor ([ws, d_lo]) depend on f;
# they are added per focal date.  Prefix sums run per (analyte, draw) series
# and per block of focal dates re-centred on the block's first date, so that
# exp((m_k - centre) / tau) stays within double range.
#
# Matches the per-focal-date definition: zero-length intervals (same-day
# duplicates, a sample on the focal date) carry no weight but count in
# n_samples_in_window; the anchor never counts; a window whose values all have
# zero weight is NA; a (focal date, point) with no in-window sample is dropped.

## Exponent span covered by one prefix-sum block (exp(300) ~ 2e130).
.TWA_BLOCK_EXP <- 300
## exp(-x) underflows to exactly 0 in double precision beyond this x.
.TWA_EXP_UNDERFLOW <- 745.1332191019411
## Target size of the per-chunk (focal date x series) work vectors.
.TWA_CHUNK <- 2e6

#' Vectorised chronic aggregation over a focal grid
#'
#' @param df Long chemistry with a `.date` column (see
#'   [time_weighted_aggregate()]).
#' @param focal_grid Data frame `(focal_date, site_id)`; one output block per
#'   row, in row order.
#' @param tau_days,window_days,summary,anchor_outside_window,eps As in
#'   [time_weighted_aggregate()].
#' @param has_imputed Logical; `df` carries an `imputed` column.
#' @param draws Active draw domain ([.draw_domain()]); `integer(0)` for point
#'   input.
#' @return `NULL` when nothing could be aggregated, else a data frame
#'   `(r, analyte, draw, value, n_in, n_imp)` ordered by focal-grid row `r`,
#'   then analyte (C locale), then draw.
#' @keywords internal
.twa_compute <- function(df, focal_grid, tau_days, window_days, summary,
                         anchor_outside_window, eps, has_imputed, draws) {
  cols <- list(
    sample_id = df$sample_id,
    date      = as.numeric(df$.date),
    analyte   = df$analyte,
    value     = df$value,
    imputed   = if (has_imputed) {
      x <- df$imputed & TRUE
      !is.na(x) & x
    } else {
      rep(FALSE, nrow(df))
    },
    draw_id   = if (length(draws) > 0L) as.integer(df$draw_id) else NULL
  )

  grid_key <- as.character(focal_grid$site_id)
  grid_f   <- as.numeric(focal_grid$focal_date)
  rows_by_site <- split(seq_len(nrow(df)), as.character(df$site_id))
  sites <- intersect(unique(grid_key), names(rows_by_site))

  pieces <- lapply(sites, function(s) {
    rows <- rows_by_site[[s]]
    r_idx <- which(grid_key == s)
    .twa_site(
      lapply(cols, function(v) if (is.null(v)) NULL else v[rows]),
      f = grid_f[r_idx], r_idx = r_idx,
      tau_days = tau_days, window_days = window_days, summary = summary,
      anchor_outside_window = anchor_outside_window, eps = eps, draws = draws
    )
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1L))]
  if (length(pieces) == 0L) {
    return(NULL)
  }
  res <- do.call(rbind, pieces)
  if (nrow(res) == 0L) {
    return(NULL)
  }
  # Each grid row belongs to one point and arrives in (analyte, draw) order,
  # so a stable sort on the row index gives the final order.
  res[order(res$r, method = "radix"), , drop = FALSE]
}

#' Chronic aggregation for one point over all of its focal dates
#'
#' @param x List of the point's row vectors: `sample_id`, `date` (numeric
#'   days), `analyte`, `value`, `imputed` (logical), `draw_id` (or `NULL`).
#' @param f Numeric focal dates (days); `r_idx` their focal-grid rows.
#' @inheritParams .twa_compute
#' @return Data frame `(r, analyte, draw, value, n_in, n_imp)` or `NULL`.
#' @keywords internal
.twa_site <- function(x, f, r_idx, tau_days, window_days, summary,
                      anchor_outside_window, eps, draws) {
  ## ── Sample schedule: one entry per (sample_id, date), sorted by date ──────
  ## (stable, so same-day samples keep their order of appearance).
  sid_u <- unique(x$sample_id)
  first <- match(x$sample_id, sid_u)
  d_first <- x$date[match(sid_u, x$sample_id)]
  one_date <- identical(is.na(x$date), is.na(d_first[first])) &&
    all(x$date == d_first[first], na.rm = TRUE)
  if (one_date) {
    sch_d <- d_first
    k_row <- first
  } else {
    key <- paste(x$sample_id, x$date, sep = "\r")
    key_u <- unique(key)
    k_row <- match(key, key_u)
    sch_d <- x$date[match(key_u, key)]
  }
  o <- order(sch_d, method = "radix") # stable; NA dates last
  o <- o[!is.na(sch_d[o])]
  n <- length(o)
  if (n == 0L) {
    return(NULL)
  }
  k_row <- match(k_row, o) # schedule position; NA for undated samples
  sch_d <- sch_d[o]

  dlt <- c(diff(sch_d), NA_real_) # forward-step interval length
  mid <- sch_d + dlt / 2 # interval midpoint
  dpos <- c(dlt[-n] > 0, FALSE)

  ## ── Rows, broadcast across draws ──────────────────────────────────────────
  keep <- which(!is.na(k_row))
  pos <- keep
  dr <- if (length(draws) == 0L) rep(1L, length(keep)) else x$draw_id[keep]
  if (length(draws) > 0L && anyNA(dr)) {
    ## Exact cells (draw_id = NA) join every draw, after the draw rows
    ## (.broadcast_draws() order).
    ex <- keep[is.na(dr)]
    dk <- keep[!is.na(dr)]
    nd <- length(draws)
    pos <- c(dk, length(x$date) + rep(ex, times = nd))
    dr <- c(dr[!is.na(dr)], rep(as.integer(draws), each = length(ex)))
    keep <- c(dk, rep(ex, times = nd))
  }
  k_r <- k_row[keep]
  v_r <- x$value[keep]
  an_r <- x$analyte[keep]
  imp_r <- x$imputed[keep]

  ## Series = (analyte, draw), numbered in (analyte C-locale, draw) order.
  an_lev <- sort(unique(an_r), method = "radix", na.last = TRUE)
  nd1 <- max(dr) + 1
  g_code <- (match(an_r, an_lev) - 1) * nd1 + dr
  g_lev <- sort(unique(g_code))
  g_r <- match(g_code, g_lev)
  n_g <- length(g_lev)
  grp_an <- an_lev[(g_lev - 1) %/% nd1 + 1]
  grp_dr <- (g_lev - 1) %% nd1 + 1

  ## Rows sorted by (series, schedule position, row order); cells are the
  ## distinct (series, position) pairs.  Composite key = series * base + k.
  base <- n + 2
  ro <- order(g_r, k_r, pos, method = "radix")
  g_r <- g_r[ro]
  k_r <- k_r[ro]
  v_r <- v_r[ro]
  imp_r <- imp_r[ro]
  pos <- pos[ro]
  r_key <- g_r * base + k_r
  new_cell <- c(TRUE, r_key[-1L] != r_key[-length(r_key)])
  cell <- cumsum(new_cell)
  c_key <- r_key[new_cell]
  c_k <- k_r[new_cell]
  non_na <- !is.na(v_r)

  is_p90 <- summary == "p90"
  xv <- if (is_p90) {
    numeric(length(v_r))
  } else if (summary == "arith_mean") {
    as.numeric(v_r)
  } else {
    log(v_r + eps)
  }
  xv[!non_na] <- 0
  fin <- is.finite(xv)
  cm <- rowsum(
    cbind(
      rows  = 1,
      nonna = non_na,
      imp   = imp_r,
      x     = ifelse(fin, xv, 0),
      nan   = non_na & is.nan(xv),
      pinf  = non_na & !is.na(xv) & xv == Inf,
      ninf  = non_na & !is.na(xv) & xv == -Inf
    ),
    cell,
    reorder = FALSE
  )
  c_dpos <- dpos[c_k]
  cum0 <- function(z) c(0, cumsum(z))
  cc_rows <- cum0(cm[, "rows"])
  cc_imp <- cum0(cm[, "imp"])
  # Middle samples (lo..hi-1) carry weight only on a positive interval.
  cc_ok <- cum0(cm[, "nonna"] * c_dpos)
  cc_nan <- cum0(cm[, "nan"] * c_dpos)
  cc_pinf <- cum0(cm[, "pinf"] * c_dpos)
  cc_ninf <- cum0(cm[, "ninf"] * c_dpos)

  ## ── Focal windows ─────────────────────────────────────────────────────────
  ws <- f - window_days
  lo <- findInterval(ws, sch_d, left.open = TRUE) + 1L # first d_k >= ws
  hi <- findInterval(f, sch_d) # last d_k <= f
  valid <- !is.na(lo) & !is.na(hi) & lo <= hi
  if (!any(valid)) {
    return(NULL)
  }
  f <- f[valid]
  ws <- ws[valid]
  lo <- lo[valid]
  hi <- hi[valid]
  r_idx <- r_idx[valid]
  has_a <- anchor_outside_window & lo > 1L
  a_k <- lo - has_a # first schedule position in use

  # Anchor: [ws, d_lo], clipped to the window; last in-window sample: [d_hi, f].
  dlt_a <- sch_d[lo] - ws
  w_a <- ifelse(
    has_a, exp(-(f - (ws + dlt_a / 2)) / tau_days) * pmax(dlt_a, 0), 0
  )
  dlt_h <- f - sch_d[hi]
  w_h <- exp(-(f - (sch_d[hi] + dlt_h / 2)) / tau_days) * pmax(dlt_h, 0)

  # Middle samples whose weight exp(-(f - m_k) / tau) underflows to 0 carry
  # none (only possible when window / tau is very large).
  lo_v <- lo
  if (window_days / tau_days > 700 && n > 1L) {
    live <- findInterval(f - .TWA_EXP_UNDERFLOW * tau_days, mid[-n],
      left.open = TRUE
    ) + 1L
    lo_v <- pmax(lo, live)
  }

  ## ── Prefix sums of the middle-sample weights, per block x series ─────────
  if (!is_p90) {
    span <- .TWA_BLOCK_EXP * tau_days
    blk <- floor((f - min(f)) / span)
    ublk <- sort(unique(blk))
    b_i <- match(blk, ublk)
    centre <- min(f) + ublk * span
    n_b <- length(ublk)
    k_lo_b <- as.numeric(lapply(split(lo_v, b_i), min))
    k_hi_b <- as.numeric(lapply(split(hi, b_i), max)) - 1
    bb <- rep(seq_len(n_b), each = n_g)
    gg <- rep(seq_len(n_g), times = n_b)
    from <- findInterval(gg * base + k_lo_b[bb] - 1, c_key) + 1L
    to <- findInterval(gg * base + k_hi_b[bb], c_key)
    cnt <- pmax(to - from + 1L, 0L)
    e_cell <- sequence(cnt, from = from)
    e_seg <- rep((bb - 1L) * n_g + gg, cnt)
    e_b <- rep(bb, cnt)
    e_k <- c_k[e_cell]
    e_key <- e_seg * base + e_k
    e_w <- ifelse(
      c_dpos[e_cell], dlt[e_k] * exp((mid[e_k] - centre[e_b]) / tau_days), 0
    )
    cum_w <- if (length(e_w)) collapse::fcumsum(e_w * cm[e_cell, "nonna"], g = e_seg) else numeric(0)
    cum_wx <- if (length(e_w)) collapse::fcumsum(e_w * cm[e_cell, "x"], g = e_seg) else numeric(0)
    pref <- function(cum, seg, k) {
      p <- findInterval(seg * base + k, e_key)
      out <- numeric(length(p))
      hit <- p > 0L
      hit[hit] <- e_seg[p[hit]] == seg[hit]
      out[hit] <- cum[p[hit]]
      out
    }
  }

  ## ── Resolve (focal date x series) in chunks ───────────────────────────────
  n_f <- length(f)
  chunk_f <- if (is_p90) {
    rows_per_f <- length(k_r) / n_g * mean(hi - a_k + 1) / n
    max(1L, floor(.TWA_CHUNK / (n_g * max(1, rows_per_f))))
  } else {
    max(1L, floor(.TWA_CHUNK / n_g))
  }
  chunks <- split(seq_len(n_f), ceiling(seq_len(n_f) / chunk_f))

  out <- lapply(chunks, function(ii) {
    qi <- rep(ii, each = n_g)
    qg <- rep(seq_len(n_g), times = length(ii))
    gb <- qg * base
    rng <- function(cc, from_k, to_k) {
      cc[findInterval(gb + to_k, c_key) + 1L] -
        cc[findInterval(gb + from_k - 1, c_key) + 1L]
    }
    # A series is reported when it has any row in the samples in use.
    exist <- rng(cc_rows, a_k[qi], hi[qi]) > 0
    if (!any(exist)) {
      return(NULL)
    }
    qi <- qi[exist]
    qg <- qg[exist]
    gb <- gb[exist]
    n_in <- rng(cc_rows, lo[qi], hi[qi])
    n_imp <- rng(cc_imp, lo[qi], hi[qi])

    # Anchor and last-sample cells of each series (NA when absent).
    ia <- match(gb + a_k[qi], c_key)
    ia[!has_a[qi]] <- NA_integer_
    ih <- match(gb + hi[qi], c_key)

    value <- if (is_p90) {
      .twa_p90(
        qi, qg, gb, r_key, k_r, v_r, pos, mid, dlt, f, tau_days,
        a_k, hi, has_a, w_a, w_h
      )
    } else {
      wa <- w_a[qi]
      wh <- w_h[qi]
      ta <- !is.na(ia) & wa > 0
      th <- !is.na(ih) & wh > 0
      cell_sum <- function(col, idx, take) {
        z <- numeric(length(idx))
        z[take] <- cm[idx[take], col]
        z
      }
      has_mid <- hi[qi] - 1L >= lo_v[qi]
      mid_rng <- function(cc) {
        z <- rng(cc, lo_v[qi], hi[qi] - 1L)
        z[!has_mid] <- 0
        z
      }
      n_ok <- mid_rng(cc_ok) + cell_sum("nonna", ia, ta) +
        cell_sum("nonna", ih, th)
      n_nan <- mid_rng(cc_nan) + cell_sum("nan", ia, ta) + cell_sum("nan", ih, th)
      n_pinf <- mid_rng(cc_pinf) + cell_sum("pinf", ia, ta) +
        cell_sum("pinf", ih, th)
      n_ninf <- mid_rng(cc_ninf) + cell_sum("ninf", ia, ta) +
        cell_sum("ninf", ih, th)

      seg <- (b_i[qi] - 1L) * n_g + qg
      fac <- exp((centre[b_i[qi]] - f[qi]) / tau_days)
      mid_w <- fac * (pref(cum_w, seg, hi[qi] - 1L) - pref(cum_w, seg, lo_v[qi] - 1L))
      mid_wx <- fac * (pref(cum_wx, seg, hi[qi] - 1L) - pref(cum_wx, seg, lo_v[qi] - 1L))
      mid_w[!has_mid] <- 0
      mid_wx[!has_mid] <- 0
      sw <- mid_w + wa * cell_sum("nonna", ia, ta) + wh * cell_sum("nonna", ih, th)
      swx <- mid_wx + wa * cell_sum("x", ia, ta) + wh * cell_sum("x", ih, th)

      val <- swx / sw
      if (summary == "geom_mean") val <- exp(val)
      # Non-finite terms decide the sum, as they would in sum(w * x).
      val[n_ninf > 0] <- if (summary == "geom_mean") 0 else -Inf
      val[n_pinf > 0] <- Inf
      val[n_nan > 0 | (n_pinf > 0 & n_ninf > 0)] <- NaN
      val[n_ok == 0] <- NA_real_
      val
    }

    data.frame(
      r       = r_idx[qi],
      analyte = grp_an[qg],
      draw    = grp_dr[qg],
      value   = value,
      n_in    = as.integer(n_in),
      n_imp   = as.integer(n_imp),
      stringsAsFactors = FALSE
    )
  })
  out <- out[!vapply(out, is.null, logical(1L))]
  if (length(out) == 0L) {
    return(NULL)
  }
  do.call(rbind, out)
}

#' Duration-weighted 90th percentile for a set of (focal date x series) windows
#'
#' Expands each window to its rows, weights them exactly as the per-window
#' definition does, and takes the first value (in value order, ties in row
#' order) whose cumulative weight share reaches 0.9.
#' @keywords internal
.twa_p90 <- function(qi, qg, gb, r_key, k_r, v_r, pos, mid, dlt, f, tau_days,
                     a_k, hi, has_a, w_a, w_h) {
  n_q <- length(qi)
  from <- findInterval(gb + a_k[qi] - 1, r_key) + 1L
  to <- findInterval(gb + hi[qi], r_key)
  cnt <- pmax(to - from + 1L, 0L)
  pr <- sequence(cnt, from = from)
  pq <- rep(seq_len(n_q), cnt)
  k_p <- k_r[pr]
  q_p <- qi[pq]
  w_p <- exp(-(f[q_p] - mid[k_p]) / tau_days) * pmax(dlt[k_p], 0)
  is_a <- has_a[q_p] & k_p == a_k[q_p]
  is_h <- k_p == hi[q_p]
  w_p[is_a] <- w_a[q_p][is_a]
  w_p[is_h] <- w_h[q_p][is_h]
  v_p <- v_r[pr]
  ok <- !is.na(v_p) & !is.na(w_p) & w_p > 0

  value <- rep(NA_real_, n_q)
  if (!any(ok)) {
    return(value)
  }
  pq <- pq[ok]
  v_p <- v_p[ok]
  w_p <- w_p[ok]
  o <- order(pq, v_p, pos[pr][ok], method = "radix")
  pq <- pq[o]
  v_p <- v_p[o]
  w_p <- w_p[o]
  u_q <- unique(pq)
  tot <- collapse::fsum(w_p, g = pq, use.g.names = FALSE)
  share <- collapse::fcumsum(w_p, g = pq) / tot[match(pq, u_q)]
  hit <- which(share >= 0.9)
  hit <- hit[!duplicated(pq[hit])]
  value[pq[hit]] <- v_p[hit]
  value
}
