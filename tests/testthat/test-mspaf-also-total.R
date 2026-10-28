## also_total: total (non-ARA) msPAF alongside ARA msPAF from one run
## ============================================================================
##
## CONTRACT.
##
## In point mode (no `ndraws`) the static ARA reference passed as
## `reference` to mspaf_daily() (or via `daily_args = list(reference = )` in
## mspaf_pipeline()) only enters at the add_mspaf() call(s); imputation, the
## reference model, the target model, Kalman smoothing and daily prediction
## never read it.  A second add_mspaf(reference = NULL) pass over the same
## synthetic daily chemistry therefore yields total msPAF cheaply.
##
## A1. New argument `also_total = FALSE` on mspaf_daily() and mspaf_pipeline().
##     FALSE leaves output unchanged (no `mspaf_total` / `value_total`).
## A2. mspaf_daily(reference = <df>, also_total = TRUE) returns the same rows
##     as without it plus a column `mspaf_total`, equal (1e-10) to the `mspaf`
##     of a separate run with reference = NULL.  Its `mspaf` equals the
##     ARA-only run's.
## A3. mspaf_pipeline(daily_args = list(reference = <df>), also_total = TRUE)
##     returns the same rows plus `value_total`: the chronic aggregate of the
##     total daily series, equal to `value` of a separate pipeline run with no
##     daily_args$reference (same focal dates, tau, window).  `value` and
##     `n_samples_in_window` equal the ARA-only run's.  (mspaf_pipeline()'s own
##     `reference` argument is the reference-model fit data, not the static
##     ARA reference, so the daily_args route is the only one.)
## A4. also_total = TRUE with no reference: `mspaf_total` == `mspaf` and
##     `value_total` == `value`, no error.
## A5. also_total = TRUE in draws mode (ndraws given): classed error
##     `hydroSense_also_total_draws` (point mode only), from both functions.
## A6. also_total must be a single non-NA logical; anything else errors.

.at_sc <- function() .pf_scenario()

.at_tmp <- function() {
  tmp <- .at_sc()$temperature
  tmp[tmp$datetime >= as.Date("2021-12-15") &
        tmp$datetime <= as.Date("2022-03-10"), ]
}

.at_target <- function() {
  tgt <- .pf_small_target()
  tgt[tgt$site_id %in% c("P1", "P2"), , drop = FALSE]
}

.at_daily <- function(...) {
  sc <- .at_sc()
  args <- list(
    .at_target(), temperature = .at_tmp(),
    reference_model = .pf_reference_model(), interpolation = "model",
    start = as.Date("2022-01-01"), end = as.Date("2022-02-28")
  )
  args <- utils::modifyList(args, list(...), keep.null = TRUE)
  suppressWarnings(suppressMessages(do.call(mspaf_daily, args)))
}

.at_pipeline <- function(daily_args = list(), ...) {
  sc <- .at_sc()
  args <- list(
    .at_target(), reference = sc$reference, hydro = sc$rain, impute = FALSE,
    reference_args = list(
      min_obs_model = 20L,
      api_tau_bounds_short = c(3, 10),
      api_tau_bounds_long  = c(30, 90)
    ),
    daily_args = utils::modifyList(list(
      temperature = .at_tmp(),
      start = as.Date("2022-01-01"), end = as.Date("2022-02-28")
    ), daily_args),
    chronic = TRUE,
    focal_dates = as.Date(c("2022-01-20", "2022-02-10", "2022-02-28")),
    tau = 90, window = 365
  )
  args <- utils::modifyList(args, list(...), keep.null = TRUE)
  out <- suppressWarnings(suppressMessages(do.call(mspaf_pipeline, args)))
  attributes(out)[c("daily", "reference_model", "imputation_model")] <- NULL
  out
}

.at_strip <- function(x) {
  for (a in setdiff(names(attributes(x)),
                    c("names", "row.names", "class"))) attr(x, a) <- NULL
  x
}

## Memoised runs shared across tests.
.at_cache <- new.env()
.at_get <- function(key, f) {
  if (is.null(.at_cache[[key]])) .at_cache[[key]] <- f()
  .at_cache[[key]]
}
.at_d_ara   <- function() .at_get("d_ara",   function() .at_daily(reference = .at_sc()$reference))
.at_d_tot   <- function() .at_get("d_tot",   function() .at_daily(reference = .at_sc()$reference, also_total = TRUE))
.at_d_none  <- function() .at_get("d_none",  function() .at_daily())
.at_p_ara   <- function() .at_get("p_ara",   function() .at_pipeline(list(reference = .at_sc()$reference)))
.at_p_tot   <- function() .at_get("p_tot",   function() .at_pipeline(list(reference = .at_sc()$reference), also_total = TRUE))
.at_p_none  <- function() .at_get("p_none",  function() .at_pipeline())

describe("A2: mspaf_daily(also_total = TRUE)", {
  skip_on_cran()

  it("adds mspaf_total and leaves the other columns and rows unchanged", {
    ara <- .at_d_ara(); tot <- .at_d_tot()
    expect_true("mspaf_total" %in% names(tot))
    expect_false("mspaf_total" %in% names(ara))
    expect_identical(nrow(tot), nrow(ara))
    expect_equal(.at_strip(tot[, names(ara)]), .at_strip(ara), tolerance = 1e-10)
  })

  it("mspaf_total equals mspaf of a reference = NULL run", {
    tot <- .at_d_tot(); none <- .at_d_none()
    expect_identical(tot$site_id, none$site_id)
    expect_identical(tot$date, none$date)
    expect_equal(tot$mspaf_total, none$mspaf, tolerance = 1e-10)
  })
})

describe("A3: mspaf_pipeline(also_total = TRUE)", {
  skip_on_cran()

  it("adds value_total and leaves the ARA columns unchanged", {
    ara <- .at_p_ara(); tot <- .at_p_tot()
    expect_true("value_total" %in% names(tot))
    expect_false("value_total" %in% names(ara))
    expect_identical(nrow(tot), nrow(ara))
    expect_equal(tot[, names(ara)], ara, tolerance = 1e-10)
    expect_equal(tot$n_samples_in_window, ara$n_samples_in_window)
  })

  it("value_total equals value of a run with no static reference", {
    tot <- .at_p_tot(); none <- .at_p_none()
    expect_equal(tot$value_total, none$value, tolerance = 1e-10)
  })
})

describe("A4: also_total = TRUE with no reference", {
  skip_on_cran()

  it("mspaf_daily(): mspaf_total equals mspaf", {
    d <- .at_daily(also_total = TRUE)
    expect_true("mspaf_total" %in% names(d))
    expect_equal(d$mspaf_total, d$mspaf, tolerance = 1e-10)
  })

  it("mspaf_pipeline(): value_total equals value", {
    p <- .at_pipeline(also_total = TRUE)
    expect_true("value_total" %in% names(p))
    expect_equal(p$value_total, p$value, tolerance = 1e-10)
  })
})

describe("A5: draws mode is rejected", {
  it("mspaf_daily() aborts with hydroSense_also_total_draws", {
    expect_error(
      .at_daily(also_total = TRUE, ndraws = 3L, seed = 1L),
      class = "hydroSense_also_total_draws"
    )
  })

  it("mspaf_pipeline() aborts with hydroSense_also_total_draws", {
    expect_error(
      .at_pipeline(list(ndraws = 3L, seed = 1L), also_total = TRUE),
      class = "hydroSense_also_total_draws"
    )
  })
})

## A validation error must come from checking also_total, not from R rejecting
## an argument the function does not have yet.
.at_expect_validation_error <- function(expr, info) {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  expect_false(is.null(err), info = info)
  if (!is.null(err)) {
    msg <- paste(conditionMessage(err), collapse = " ")
    expect_false(grepl("unused argument", msg), info = info)
    expect_match(msg, "also_total", info = info)
  }
}

describe("A6: also_total validation", {
  bad <- list(NA, "yes", 1, c(TRUE, FALSE), NULL, logical(0))

  it("mspaf_daily() rejects non-flag also_total", {
    for (b in bad) {
      .at_expect_validation_error(.at_daily(also_total = b), deparse(b))
    }
  })

  it("mspaf_pipeline() rejects non-flag also_total", {
    for (b in bad) {
      .at_expect_validation_error(.at_pipeline(also_total = b), deparse(b))
    }
  })
})
