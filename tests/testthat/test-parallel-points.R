## ============================================================================
## Parallel per-point work in mspaf_daily() / mspaf_pipeline()
## ============================================================================
##
## CONTRACT.
##
## mspaf_daily(parallel = TRUE) dispatches the per-point ("site") loop through
## future.apply::future_lapply(); each point's target-model fit
## (.fit_daily_target) and daily prediction (.predict_daily_tox) run INSIDE
## that parallel map.  Draws run sequentially within a point.
##
## P1. Point mode (no ndraws): output with parallel = TRUE under
##     future::plan(multisession, workers = 2) is IDENTICAL to serial, via
##     mspaf_pipeline(daily_args = list(parallel = TRUE)), for: static ARA
##     reference, no reference, also_total = TRUE, impute = TRUE; the chronic
##     result and the "daily" attribute alike.
## P2. mspaf_daily(parallel = TRUE) is identical to serial, including the
##     analyte_pafs / ara_summary / impact_tiers attributes (row order too).
## P3. Structural: with a forking plan, the per-point fit runs in worker
##     processes (pid differs from the parent), once per point.
## P4. Globals guard: with future.globals.maxSize = 200 MB, the parallel map
##     does not error, i.e. no future ships more than 200 MB of globals.
## P5. Draws mode: parallel + draws is supported; per-point seeds come from
##     future.apply's L'Ecuyer streams, so results are NOT required to equal
##     serial (the stream differs by design), but two parallel runs with the
##     same seed are identical.
##
## NB: multisession workers load the INSTALLED hydroSense, not a
## devtools::load_all() tree; P1/P2/P4/P5 therefore test the installed build.
## ============================================================================

.pp_skip <- function() {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_not_installed("withr")
}

.pp_plan <- function(strategy = future::multisession, workers = 2L,
                     env = parent.frame()) {
  old <- future::plan()
  withr::defer(future::plan(old), envir = env)
  withr::defer(future::plan(future::sequential), envir = env)
  future::plan(strategy, workers = workers)
}

## multisession workers load the INSTALLED hydroSense.  If that differs from the
## code under test (e.g. a devtools::load_all() tree, or a stale install) an
## equality check would compare two different implementations, so skip.
.pp_skip_if_stale_install <- function() {
  here <- paste(deparse(hydroSense:::mspaf_daily), collapse = "\n")
  there <- future::value(future::future(
    paste(deparse(hydroSense:::mspaf_daily), collapse = "\n"),
    globals = FALSE
  ))
  skip_if(!identical(here, there),
          "multisession workers load a hydroSense that differs from the code under test")
}

.pp_sc  <- function() .pf_scenario()
.pp_tmp <- function() {
  tmp <- .pp_sc()$temperature
  tmp[tmp$datetime >= as.Date("2021-12-15") &
        tmp$datetime <= as.Date("2022-03-10"), ]
}

.pp_pipeline <- function(par, daily_args = list(), ...) {
  sc <- .pp_sc()
  args <- list(
    .pf_small_target(), reference = sc$reference, hydro = sc$rain,
    impute = FALSE, impute_seed = 1L,
    reference_args = list(
      min_obs_model = 20L,
      api_tau_bounds_short = c(3, 10),
      api_tau_bounds_long  = c(30, 90)
    ),
    daily_args = utils::modifyList(list(
      temperature = .pp_tmp(), parallel = par,
      start = as.Date("2022-01-01"), end = as.Date("2022-02-28")
    ), daily_args),
    chronic = TRUE,
    focal_dates = as.Date(c("2022-01-20", "2022-02-10", "2022-02-28")),
    tau = 90, window = 365
  )
  args <- utils::modifyList(args, list(...), keep.null = TRUE)
  suppressWarnings(suppressMessages(do.call(mspaf_pipeline, args)))
}

.pp_expect_same_pipeline <- function(...) {
  ser <- .pp_pipeline(FALSE, ...)
  par <- .pp_pipeline(TRUE, ...)
  expect_gt(length(unique(ser$site_id)), 2L)  # >= 3 points really split
  expect_identical(
    attr(par, "daily"), attr(ser, "daily"),
    label = "daily attribute (parallel)", expected.label = "serial"
  )
  keep <- c("daily", "reference_model", "imputation_model")
  strip <- function(x) { attributes(x)[keep] <- NULL; x }
  expect_equal(strip(par), strip(ser), tolerance = 0)
  expect_identical(strip(par), strip(ser))
  invisible(list(ser = ser, par = par))
}

test_that("P1a: pipeline parallel == serial, static ARA reference", {
  .pp_skip()
  .pp_plan()
  .pp_skip_if_stale_install()
  .pp_expect_same_pipeline(daily_args = list(reference = .pp_sc()$reference))
})

test_that("P1b: pipeline parallel == serial, no reference", {
  .pp_skip()
  .pp_plan()
  .pp_skip_if_stale_install()
  .pp_expect_same_pipeline()
})

test_that("P1c: pipeline parallel == serial, also_total = TRUE", {
  .pp_skip()
  .pp_plan()
  .pp_skip_if_stale_install()
  r <- .pp_expect_same_pipeline(
    daily_args = list(reference = .pp_sc()$reference), also_total = TRUE
  )
  expect_true("value_total" %in% names(r$par))
})

test_that("P1d: pipeline parallel == serial, impute = TRUE", {
  .pp_skip()
  .pp_plan()
  .pp_skip_if_stale_install()
  .pp_expect_same_pipeline(
    impute = TRUE, daily_args = list(reference = .pp_sc()$reference)
  )
})

.pp_daily <- function(par, ...) {
  args <- list(
    .pf_small_target(), temperature = .pp_tmp(),
    reference_model = .pf_reference_model(), interpolation = "model",
    reference = .pp_sc()$reference,
    start = as.Date("2022-01-01"), end = as.Date("2022-02-28"),
    parallel = par
  )
  args <- utils::modifyList(args, list(...), keep.null = TRUE)
  suppressWarnings(suppressMessages(do.call(mspaf_daily, args)))
}

test_that("P2: mspaf_daily parallel == serial incl. attributes and order", {
  .pp_skip()
  .pp_plan()
  .pp_skip_if_stale_install()
  ser <- .pp_daily(FALSE)
  par <- .pp_daily(TRUE)
  expect_gt(length(unique(ser$site_id)), 2L)
  expect_identical(par, ser)
  for (a in c("analyte_pafs", "ara_summary", "impact_tiers")) {
    expect_identical(attr(par, a), attr(ser, a), label = paste0("attr ", a))
  }
  expect_setequal(
    intersect(names(attributes(ser)), c("analyte_pafs", "ara_summary")),
    intersect(names(attributes(par)), c("analyte_pafs", "ara_summary"))
  )
})

test_that("P3: per-point fit runs inside the parallel map (worker pids)", {
  .pp_skip()
  skip_if_not(future::supportsMulticore(), "forking not supported here")
  .pp_plan(future::multicore, workers = 2L)
  log <- tempfile("pp-pids-")
  withr::defer(unlink(log))
  parent <- Sys.getpid()
  orig <- hydroSense:::.fit_daily_target
  testthat::local_mocked_bindings(
    .fit_daily_target = function(...) {
      cat(Sys.getpid(), "\n", file = log, append = TRUE)
      orig(...)
    },
    .package = "hydroSense"
  )
  invisible(.pp_daily(TRUE))
  pids <- as.integer(readLines(log))
  expect_length(pids, 3L)
  expect_false(
    any(pids == parent),
    info = paste("per-point .fit_daily_target() ran in the parent process;",
                 "fit must happen inside the future_lapply() map")
  )
})

test_that("P4: parallel map ships < 200 MB globals per future", {
  .pp_skip()
  .pp_plan()
  .pp_skip_if_stale_install()
  withr::local_options(future.globals.maxSize = 200 * 1024^2)
  expect_no_error(par <- .pp_daily(TRUE))
  expect_s3_class(par, "data.frame")
  expect_no_error(.pp_pipeline(
    TRUE, impute = TRUE, daily_args = list(reference = .pp_sc()$reference)
  ))
})

test_that("P5: parallel draws reproducible for a fixed seed", {
  .pp_skip()
  .pp_plan()
  .pp_skip_if_stale_install()
  run <- function() {
    .pp_daily(TRUE, ndraws = 4L, seed = 7L, return = "summary",
              reference = NULL)
  }
  a <- run()
  b <- run()
  expect_identical(a, b)
  expect_gt(length(unique(a$site_id)), 2L)
  ## Different seed changes the stochastic output.
  d <- .pp_daily(TRUE, ndraws = 4L, seed = 8L, return = "summary",
                 reference = NULL)
  expect_false(isTRUE(all.equal(a, d)))
})

## Forking variants run the in-tree code in the workers, so they hold even when
## the installed hydroSense is stale.
test_that("P2-fork: mspaf_daily parallel == serial under multicore", {
  .pp_skip()
  skip_if_not(future::supportsMulticore(), "forking not supported here")
  .pp_plan(future::multicore, workers = 2L)
  ser <- .pp_daily(FALSE)
  par <- .pp_daily(TRUE)
  expect_gt(length(unique(ser$site_id)), 2L)
  expect_identical(par, ser)
  for (a in c("analyte_pafs", "ara_summary", "impact_tiers")) {
    expect_identical(attr(par, a), attr(ser, a), label = paste0("attr ", a))
  }
})

test_that("P1-fork: pipeline parallel == serial under multicore (impute)", {
  .pp_skip()
  skip_if_not(future::supportsMulticore(), "forking not supported here")
  .pp_plan(future::multicore, workers = 2L)
  .pp_expect_same_pipeline(
    impute = TRUE, daily_args = list(reference = .pp_sc()$reference)
  )
})

test_that("P5-fork: parallel draws reproducible under multicore", {
  .pp_skip()
  skip_if_not(future::supportsMulticore(), "forking not supported here")
  .pp_plan(future::multicore, workers = 2L)
  run <- function(s) .pp_daily(TRUE, ndraws = 4L, seed = s,
                               return = "summary", reference = NULL)
  expect_identical(run(7L), run(7L))
  expect_false(isTRUE(all.equal(run(7L), run(8L))))
})
