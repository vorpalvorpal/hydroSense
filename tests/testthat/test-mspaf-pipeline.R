## mspaf_pipeline() is an orchestrator: it sequences fit_imputation_model(),
## fit_reference_model() and mspaf_daily() and threads one imputation model
## through both downstream steps. These specs verify the *wiring* (the new
## behaviour) with mocked inner functions, so they are fast and Stan-free; the
## inner functions are covered by their own tests.

library(testthat)
library(hydroSense)

tgt <- tibble::tibble(
  sample_id = "t1", site_id = "target", datetime = as.Date("2024-01-01"),
  analyte = c("pH", "EC"), value = c(7, 200), detected = TRUE
)
ref <- tibble::tibble(
  sample_id = "r1", site_id = "reference", datetime = as.Date("2024-01-01"),
  analyte = c("pH", "EC"), value = c(7, 200), detected = TRUE
)
fake_daily <- function(...) {
  tibble::tibble(date = as.Date("2024-01-01"), mspaf = 1)
}

describe("mspaf_pipeline()", {
  it("threads the imputation model into both the reference and daily fits when impute = TRUE", {
    seen <- new.env()
    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list(metals = 1)), class = "imputation_model")
      },
      fit_reference_model = function(reference, ..., imputation_model = NULL) {
        seen$ref_im <- imputation_model
        "REFMODEL"
      },
      mspaf_daily = function(df, ..., reference_model = NULL,
                             imputation_model = NULL) {
        seen$daily_im <- imputation_model
        seen$daily_rm <- reference_model
        fake_daily()
      }
    )
    out <- mspaf_pipeline(tgt, ref, chronic = FALSE)
    expect_s3_class(seen$ref_im, "imputation_model")
    expect_s3_class(seen$daily_im, "imputation_model")
    expect_identical(seen$daily_rm, "REFMODEL")
    expect_true("mspaf" %in% names(out))
  })

  it("does not fit an imputation model and passes NULL downstream when impute = FALSE", {
    fitted <- FALSE
    local_mocked_bindings(
      fit_imputation_model = function(...) {
        fitted <<- TRUE
        stop("must not be called")
      },
      fit_reference_model = function(reference, ..., imputation_model = NULL) {
        expect_null(imputation_model)
        "RM"
      },
      mspaf_daily = function(df, ..., imputation_model = NULL) {
        expect_null(imputation_model)
        fake_daily()
      }
    )
    mspaf_pipeline(tgt, ref, impute = FALSE, chronic = FALSE)
    expect_false(fitted)
  })

  it("fits the imputation model on the reference by default and on the target when asked", {
    seen <- new.env()
    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        seen$site <- df$site_id[1]
        structure(list(groups = list()), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RM",
      mspaf_daily = function(df, ...) fake_daily()
    )
    mspaf_pipeline(tgt, ref, impute_on = "reference", chronic = FALSE)
    expect_identical(seen$site, "reference")
    mspaf_pipeline(tgt, ref, impute_on = "target", chronic = FALSE)
    expect_identical(seen$site, "target")
  })

  it("forwards daily_args and lets them override the interpolation default", {
    seen <- new.env()
    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list()), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RM",
      mspaf_daily = function(df, ..., interpolation = "model") {
        seen$interp <- interpolation
        seen$dots <- list(...)
        fake_daily()
      }
    )
    mspaf_pipeline(tgt, ref, chronic = FALSE, daily_args = list(
      interpolation = "forward_fill", require_temperature = FALSE
    ))
    expect_identical(seen$interp, "forward_fill")
    expect_false(seen$dots$require_temperature)
  })

  it("forwards reference_args to fit_reference_model", {
    seen <- new.env()
    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list()), class = "imputation_model")
      },
      fit_reference_model = function(reference, ..., min_obs_model = 20L) {
        seen$min_obs <- min_obs_model
        "RM"
      },
      mspaf_daily = function(df, ...) fake_daily()
    )
    mspaf_pipeline(tgt, ref, chronic = FALSE, reference_args = list(min_obs_model = 5L))
    expect_identical(seen$min_obs, 5L)
  })

  it("warns and proceeds without imputation when the imputation fit errors", {
    local_mocked_bindings(
      fit_imputation_model = function(...) stop("boom"),
      fit_reference_model = function(reference, ..., imputation_model = NULL) {
        expect_null(imputation_model)
        "RM"
      },
      mspaf_daily = function(df, ..., imputation_model = NULL) {
        expect_null(imputation_model)
        fake_daily()
      }
    )
    expect_warning(mspaf_pipeline(tgt, ref, chronic = FALSE), "Imputation skipped")
  })

  it("attaches the fitted reference and imputation models as attributes", {
    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list(m = 1)), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RMODEL",
      mspaf_daily = function(df, ...) fake_daily()
    )
    out <- mspaf_pipeline(tgt, ref, chronic = FALSE)
    expect_identical(attr(out, "reference_model"), "RMODEL")
    expect_s3_class(attr(out, "imputation_model"), "imputation_model")
  })

  it("rejects non-data-frame target or reference", {
    expect_error(mspaf_pipeline("nope", ref), "data")
    expect_error(mspaf_pipeline(tgt, 42), "data")
  })
})

## ── Chronic aggregation (final pipeline step) ──────────────────────────────
##
## These specs mock only the three "inner" steps (fit_imputation_model,
## fit_reference_model, mspaf_daily) and let the REAL time_weighted_aggregate()
## / .summarise_bracket() run, so the chronic reshape/wiring/summary added in
## .chronic_from_daily() is genuinely exercised end to end.

describe("mspaf_pipeline() chronic aggregation", {
  it("returns a chronic point frame (focal_date/site_id/value/n_samples_in_window, no mspaf) with the daily frame attached", {
    multi_day_daily <- tibble::tibble(
      date = as.Date("2024-01-01") + 0:4,
      site_id = "target",
      mspaf = c(0.1, 0.2, 0.3, 0.4, 0.5),
      n_analytes_used = 2L,
      dominant_analyte = "Cu",
      max_paf = c(0.1, 0.2, 0.3, 0.4, 0.5),
      n_measured_analytes = 2L,
      days_since_last_sample = 0
    )
    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list()), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RM",
      mspaf_daily = function(df, ...) multi_day_daily
    )

    out <- mspaf_pipeline(tgt, ref)

    expect_true(all(c("focal_date", "site_id", "value", "n_samples_in_window") %in% names(out)))
    expect_false("mspaf" %in% names(out))
    expect_identical(attr(out, "daily"), multi_day_daily)
    expect_identical(attr(out, "reference_model"), "RM")
    expect_s3_class(attr(out, "imputation_model"), "imputation_model")
  })

  it("computes the chronic value as the hand-computed exponential-decay-weighted arithmetic mean", {
    # Tiny controlled daily series: 3 consecutive days, single focal_date one
    # day after the last observation so every day gets a strictly positive
    # forward-step duration (delta_t = 1 for all three; the focal date itself
    # always gets delta_t = 0 under forward-step weighting, so picking a focal
    # date that coincides with the last data day would zero its weight).
    tiny_daily <- tibble::tibble(
      date = as.Date("2024-01-01") + 0:2,
      site_id = "target",
      mspaf = c(0.1, 0.4, 0.7)
    )
    tau_d <- 10
    window_d <- 30
    focal <- as.Date("2024-01-04")

    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list()), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RM",
      mspaf_daily = function(df, ...) tiny_daily
    )

    out <- mspaf_pipeline(
      tgt, ref,
      focal_dates = focal, tau = tau_d, window = window_d
    )

    # Hand computation, mirroring time_weighted_aggregate()'s forward-step
    # interval + midpoint exponential-decay weighting (R/chronic.R):
    #   interval_end   = c(dates[-1], focal_date)
    #   interval_start = dates
    #   delta_t        = interval_end - interval_start
    #   midpoint       = interval_start + delta_t / 2
    #   w              = delta_t * exp(-(focal_date - midpoint) / tau)
    #   value          = sum(w * mspaf) / sum(w)
    dates <- tiny_daily$date
    interval_end <- c(dates[-1], focal)
    interval_start <- dates
    delta_t <- as.numeric(interval_end - interval_start)
    midpoints <- interval_start + delta_t / 2
    w <- delta_t * exp(-as.numeric(focal - midpoints) / tau_d)
    expected_value <- sum(w * tiny_daily$mspaf) / sum(w)

    expect_equal(out$value[out$focal_date == focal & out$site_id == "target"],
      expected_value, tolerance = 1e-8)
  })

  it("spaces focal dates by focal_by (default daily, e.g. focal_by = 7 -> weekly)", {
    daily_range <- tibble::tibble(
      date = as.Date("2024-01-01") + 0:20,
      site_id = "target",
      mspaf = seq(0.1, 0.9, length.out = 21)
    )
    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list()), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RM",
      mspaf_daily = function(df, ...) daily_range
    )

    out_daily <- mspaf_pipeline(tgt, ref)
    daily_gaps <- diff(sort(unique(out_daily$focal_date)))
    expect_true(all(as.numeric(daily_gaps) == 1))

    out_weekly <- mspaf_pipeline(tgt, ref, focal_by = 7)
    weekly_gaps <- diff(sort(unique(out_weekly$focal_date)))
    expect_true(all(as.numeric(weekly_gaps) == 7))
  })

  it("returns the daily frame unchanged when chronic = FALSE, with model attrs but no daily attr", {
    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list(m = 1)), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RMODEL",
      mspaf_daily = function(df, ...) fake_daily()
    )

    out <- mspaf_pipeline(tgt, ref, chronic = FALSE)

    # Daily frame returned unchanged in content (the pipeline only adds the
    # reference_model / imputation_model attributes, so a full expect_identical
    # against a bare fake_daily() would spuriously fail on those attributes).
    expect_true("mspaf" %in% names(out))
    expect_identical(names(out), names(fake_daily()))
    expect_equal(out$mspaf, fake_daily()$mspaf)
    expect_equal(out$date, fake_daily()$date)
    expect_identical(attr(out, "reference_model"), "RMODEL")
    expect_s3_class(attr(out, "imputation_model"), "imputation_model")
    expect_null(attr(out, "daily"))
  })

  it("propagates draws-mode bracket uncertainty into chronic summary columns", {
    dates <- as.Date("2024-01-01") + 0:2
    focal <- as.Date("2024-01-04")
    draws_daily <- tidyr::expand_grid(date = dates, draw_id = 1:5) |>
      dplyr::mutate(
        site_id = "target",
        mspaf_ignorable   = stats::runif(dplyr::n(), 0.1, 0.9),
        mspaf_informative = stats::runif(dplyr::n(), 0.1, 0.9),
        n_analytes_used = 2L,
        dominant_analyte = "Cu",
        max_paf = 0.5,
        n_measured_analytes = 2L,
        days_since_last_sample = 0
      )

    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list()), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RM",
      mspaf_daily = function(df, ..., reference_model = NULL,
                             imputation_model = NULL, return = "summary",
                             ndraws = NULL, gap_uncertainty = "bracket",
                             interval = 0.9, central = "median") {
        draws_daily
      }
    )

    out <- mspaf_pipeline(
      tgt, ref,
      focal_dates = focal, tau = 10, window = 30,
      daily_args = list(ndraws = 5L, seed = 1L)
    )

    bracket_cols <- c(
      "median_ignorable", "lo_ignorable", "hi_ignorable",
      "median_informative", "lo_informative", "hi_informative",
      "precautionary_lo", "precautionary_hi"
    )
    expect_true(all(bracket_cols %in% names(out)))
    expect_identical(nrow(out), 1L) # one (focal_date x site) row
    expect_true(all(out$lo_ignorable <= out$median_ignorable))
    expect_true(all(out$median_ignorable <= out$hi_ignorable))
    expect_true(all(out$lo_informative <= out$median_informative))
    expect_true(all(out$median_informative <= out$hi_informative))
  })

  it("draws mode with gap_uncertainty = 'ignorable' returns only ignorable columns, no precautionary", {
    dates <- as.Date("2024-01-01") + 0:2
    focal <- as.Date("2024-01-04")
    draws_daily <- tidyr::expand_grid(date = dates, draw_id = 1:5) |>
      dplyr::mutate(
        site_id = "target",
        mspaf_ignorable = stats::runif(dplyr::n(), 0.1, 0.9)
      )

    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list()), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RM",
      mspaf_daily = function(df, ..., reference_model = NULL,
                             imputation_model = NULL, return = "summary",
                             ndraws = NULL, gap_uncertainty = "ignorable",
                             interval = 0.9, central = "median") {
        draws_daily
      }
    )

    out <- mspaf_pipeline(
      tgt, ref,
      focal_dates = focal, tau = 10, window = 30,
      daily_args = list(ndraws = 5L, seed = 1L, gap_uncertainty = "ignorable")
    )

    expect_true(all(c("median_ignorable", "lo_ignorable", "hi_ignorable") %in% names(out)))
    expect_false(any(c(
      "median_informative", "lo_informative", "hi_informative",
      "precautionary_lo", "precautionary_hi"
    ) %in% names(out)))
  })
})
