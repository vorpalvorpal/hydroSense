## Contract: `imputed` / `n_imputed_in_window` in the daily -> chronic path.
##
## * mspaf_daily() (interpolation = "model") returns a logical column `imputed`
##   (never NA).  A daily row is imputed when its point had NO grab sample on
##   that day (chemistry model-predicted / forward-filled); days with a grab
##   sample for that point (grab dates = as.Date(datetime), as in mspaf_daily's
##   daily grid) are not imputed.
## * mspaf_pipeline() (chronic = TRUE, point mode) reports
##   `n_imputed_in_window` = number of daily rows in [focal - window, focal]
##   with imputed == TRUE, hence
##   n_samples_in_window - n_imputed_in_window = number of measured grab days
##   in the window; 0 <= n_imputed <= n_samples; not identically 0 on sparse
##   grabs.  The chronic `value` is unchanged by this (see the pipeline golden).
## * The grab-sample path of time_weighted_aggregate() keeps propagating a
##   user-supplied `imputed` column.
##
## Synthetic fixtures only (helper-perf-fixtures.R).

.ni_target <- function() {
  tgt <- .pf_small_target()
  tgt[tgt$site_id %in% c("P1", "P2"), , drop = FALSE]
}

.ni_temperature <- function() {
  tmp <- .pf_scenario()$temperature
  tmp[tmp$datetime >= as.Date("2021-12-15") &
        tmp$datetime <= as.Date("2022-03-10"), ]
}

.ni_start <- as.Date("2022-01-01")
.ni_end   <- as.Date("2022-02-28")

.ni_daily <- function(target) {
  sc <- .pf_scenario()
  suppressWarnings(suppressMessages(mspaf_daily(
    target, temperature = .ni_temperature(), reference = sc$reference,
    reference_model = .pf_reference_model(), interpolation = "model",
    start = .ni_start, end = .ni_end
  )))
}

.ni_focal <- as.Date(c("2022-01-20", "2022-02-10", "2022-02-28"))

.ni_pipeline <- function(target, window = 30) {
  sc <- .pf_scenario()
  suppressWarnings(suppressMessages(mspaf_pipeline(
    target, reference = sc$reference, hydro = sc$rain, impute = FALSE,
    reference_args = list(
      min_obs_model = 20L,
      api_tau_bounds_short = c(3, 10),
      api_tau_bounds_long  = c(30, 90)
    ),
    daily_args = list(
      reference = sc$reference, temperature = .ni_temperature(),
      start = .ni_start, end = .ni_end
    ),
    chronic = TRUE, focal_dates = .ni_focal, tau = 90, window = window
  )))
}

## Expected `imputed` for a daily frame, derived independently from the target.
.ni_expected_imputed <- function(daily, target) {
  grab <- paste(target$site_id, as.Date(target$datetime))
  !(paste(daily$site_id, daily$date) %in% grab)
}

describe("mspaf_daily() imputed column", {
  tgt   <- .ni_target()
  daily <- .ni_daily(tgt)

  it("has a logical, non-NA `imputed` column", {
    expect_true("imputed" %in% names(daily))
    expect_type(daily$imputed, "logical")
    expect_false(anyNA(daily$imputed))
  })

  it("is FALSE exactly on days with a grab sample for that point", {
    skip_if_not("imputed" %in% names(daily), "imputed column missing")
    expect_equal(daily$imputed, .ni_expected_imputed(daily, tgt))
  })

  it("has both measured and imputed days on sparse grabs", {
    skip_if_not("imputed" %in% names(daily), "imputed column missing")
    expect_true(any(daily$imputed))
    expect_true(any(!daily$imputed))
  })
})

describe("mspaf_pipeline() n_imputed_in_window", {
  tgt <- .ni_target()
  out <- .ni_pipeline(tgt, window = 30)
  daily <- attr(out, "daily")

  it("is not identically zero on sparse grabs", {
    expect_true(any(out$n_imputed_in_window > 0))
  })

  it("lies in [0, n_samples_in_window]", {
    expect_true(all(out$n_imputed_in_window >= 0))
    expect_true(all(out$n_imputed_in_window <= out$n_samples_in_window))
  })

  it("equals n daily rows minus grab days in each focal window", {
    for (i in seq_len(nrow(out))) {
      f   <- out$focal_date[i]
      sid <- out$site_id[i]
      in_w <- daily$site_id == sid & daily$date >= f - 30 & daily$date <= f
      k <- sum(!.ni_expected_imputed(daily[in_w, ], tgt))
      n <- out$n_samples_in_window[i]
      expect_equal(n - out$n_imputed_in_window[i], k,
                   info = paste(sid, format(f)))
    }
  })
})

describe("grab-sample path", {
  it("still propagates a user-supplied imputed column", {
    tgt <- .ni_target()
    tgt <- tgt[tgt$analyte == "Cu", ]
    tgt$imputed <- FALSE
    tgt$imputed[1:3] <- TRUE
    res <- suppressWarnings(suppressMessages(time_weighted_aggregate(
      tgt, focal_dates = as.Date("2022-06-01"), tau = 90, tau_units = "d",
      window = 400, window_units = "d", return = "summary"
    )))
    expect_true(any(res$n_imputed_in_window > 0))
  })
})
