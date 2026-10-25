## Structural (vectorisation) tests for the msPAF performance refactor
## ============================================================================
##
## CONTRACT.  The refactor is a pure performance change: results must stay
## identical to the pre-refactor code within 1e-10 relative (pinned by
## test-perf-goldens.R).  THIS file pins the SHAPE of the work instead, with
## call counts rather than wall-clock timings (no flaky timing assertions):
##
##   S1  .resolve_ref_norm_instant() evaluates the antecedent precipitation
##       index in one pass over all target dates: at most 1 + 2 * n_model
##       calls of .compute_api() per resolver call (one event-API pass for
##       the tier-1 gate, one short- and one long-store pass per tier-2
##       analyte model), independent of the number of dates.
##   S2  .resolve_ref_norm_instant() predicts each tier-2 reference GAM once
##       per resolver call (at most n_model predict() calls), not once per
##       date.
##   S3  mspaf_daily() resolves the reference norm on the daily grid ONCE per
##       call, shared by every monitoring point, instead of once per point.
##   S4  The number of .compute_api() calls and reference-GAM predictions in
##       one mspaf_daily() call does not grow with the length of the daily
##       grid (O(1) per point, not one per day).
##   S5  time_weighted_aggregate() does not loop over focal dates with
##       purrr::pmap() and handles 20,000 focal dates x 3 points.
##   S6  time_weighted_aggregate() makes the same number of dplyr verb calls
##       (and per-group .aggregate_weighted() calls, if that helper still
##       exists) for 60 focal dates as for 6: its work is vectorised over
##       focal dates rather than iterated per focal date.
##
## SEAMS — what the implementer may and may not rename:
##   MUST KEEP (the tests mock or count through these):
##     * the public functions mspaf_daily(), mspaf_pipeline(),
##       time_weighted_aggregate() and their existing arguments;
##     * .compute_api(hydro_values, hydro_dates, target_dates, tau) as the
##       single entry point for API evaluation (new optional arguments are
##       fine; inlining the recursion elsewhere is not);
##     * .resolve_ref_norm_instant(ref_model, target_df) as the single entry
##       point for instantaneous reference-norm resolution, returning
##       (sample_id, analyte, ref_norm, ref_tier) (new optional arguments are
##       fine);
##     * reference-GAM predictions dispatched through the predict() S3
##       generic (stats::predict(gam, ...)), not mgcv::predict.gam() called
##       directly — S2/S4 count them through a class-based S3 method.
##   FREE TO CHANGE: every other private helper, including adding a
##     precomputed-reference argument to .fit_daily_target() or
##     fit_target_model() (the existing `ref_q` pattern of
##     .resolve_target_impact() is the model), caching within one call,
##     renaming or removing .aggregate_weighted() (S6 counts it only if it
##     exists), and using base R / data.table / dplyr inside
##     time_weighted_aggregate() as long as no per-focal-date loop remains.
##
## On the pre-refactor code every test here FAILS (by design).  The counting
## shims "trip" as soon as a bound is exceeded and then fail fast, so the
## failures are quick even though the old code is slow.

library(testthat)


## ---------------------------------------------------------------------------
## Counting helpers
## ---------------------------------------------------------------------------

## When a counting shim exceeds its bound it "trips": it records the fact and
## signals a `pf_structural_trip` condition.  That class deliberately does NOT
## inherit from "error", so the package's own tryCatch(error = ...) fallbacks
## (around fit_target_model(), the daily-grid resolver, GAM prediction) cannot
## swallow it; .pf_catch_trip() unwinds the whole call and returns the message
## so the test fails fast with a clear reason.
.pf_trip <- function(msg) {
  cond <- structure(class = c("pf_structural_trip", "condition"),
                    list(message = msg, call = NULL))
  signalCondition(cond)
  stop(msg, call. = FALSE) # no handler installed: fail loudly anyway
}

.pf_catch_trip <- function(expr) {
  tryCatch({
    force(expr)
    NULL
  }, pf_structural_trip = function(c) conditionMessage(c))
}

## Reference-GAM prediction counter: GAMs of a reference model are given an
## extra S3 class whose predict() method counts, optionally trips at a limit,
## and then defers to predict.gam() via NextMethod().
.pf_gam_ctr <- new.env()
.pf_gam_ctr$n <- 0L
.pf_gam_ctr$limit <- Inf
.pf_gam_ctr$tripped <- FALSE
registerS3method("predict", "pf_counted_gam", function(object, ...) {
  .pf_gam_ctr$n <- .pf_gam_ctr$n + 1L
  if (.pf_gam_ctr$n > .pf_gam_ctr$limit) {
    .pf_gam_ctr$tripped <- TRUE
    .pf_trip(sprintf("STRUCTURAL: reference GAM predicted more than %s times",
                     format(.pf_gam_ctr$limit)))
  }
  NextMethod()
}, envir = asNamespace("stats"))

.pf_reset_gam_ctr <- function(limit = Inf) {
  .pf_gam_ctr$n <- 0L
  .pf_gam_ctr$limit <- limit
  .pf_gam_ctr$tripped <- FALSE
}

.pf_counted_reference_model <- function() {
  rm <- .pf_reference_model()
  for (nm in names(rm$models)) {
    g <- rm$models[[nm]]$gamm_fit
    if (!is.null(g)) class(rm$models[[nm]]$gamm_fit) <- c("pf_counted_gam", class(g))
  }
  rm
}

.pf_n_model_tier <- function(rm) {
  sum(vapply(rm$models, function(m) {
    identical(m$tier, "model") && !is.null(m$gamm_fit)
  }, logical(1)))
}

## Wrap a function with a counter that trips (see .pf_trip()) once `limit()` is exceeded.  `weight(...)` lets a call count
## for more than one unit (e.g. the number of dates it resolves).
.pf_counting_shim <- function(ctr, orig, limit = function() Inf,
                              weight = function(...) 1L, what = "call") {
  force(orig)
  function(...) {
    ctr$calls <- ctr$calls + 1L
    ctr$units <- ctr$units + weight(...)
    if (ctr$units > limit()) {
      ctr$tripped <- TRUE
      .pf_trip(sprintf("STRUCTURAL: %s count %d exceeded the bound %s",
                       what, ctr$units, format(limit())))
    }
    orig(...)
  }
}

.pf_new_ctr <- function() {
  e <- new.env()
  e$calls <- 0L
  e$units <- 0L
  e$tripped <- FALSE
  e
}

## mspaf_daily() on the small 3-point target over [start, end], with an
## external temperature series covering exactly the grid (so the per-point
## date set equals the grid).
.pf_small_daily <- function(rm, start, end) {
  sc <- .pf_scenario()
  tmp <- sc$temperature
  tmp <- tmp[tmp$datetime >= start & tmp$datetime <= end, ]
  suppressWarnings(suppressMessages(mspaf_daily(
    .pf_small_target(), temperature = tmp, reference = sc$reference,
    reference_model = rm, interpolation = "model", start = start, end = end
  )))
}


## ---------------------------------------------------------------------------
## S1 / S2: the resolver works in one pass over its target dates
## ---------------------------------------------------------------------------

describe(".resolve_ref_norm_instant() is vectorised over target dates", {
  it("S1: calls .compute_api() at most 1 + 2 * n_model times for 200 dates", {
    rm <- .pf_reference_model()
    n_model <- .pf_n_model_tier(rm)
    expect_gte(n_model, 1L) # fixture sanity: tier-2 models exist
    # One event-API pass (tier-1 hydro gate, tau = match_window_days) plus a
    # short- and a long-store pass per tier-2 model: 1 + 2 * n_model (= 5).
    # The pre-refactor code makes 1 call per date plus 2 per (unmatched date
    # x tier-2 analyte): ~ 200 + 2 * 2 * 190 here.
    bound <- 1L + 2L * n_model
    qd <- seq(as.Date("2022-01-01"), by = "day", length.out = 200L)
    ctr <- .pf_new_ctr()
    local_mocked_bindings(.compute_api = .pf_counting_shim(
      ctr, hydroSense:::.compute_api, limit = function() bound,
      what = ".compute_api() call"
    ))
    res <- NULL
    trip <- .pf_catch_trip(res <- hydroSense:::.resolve_ref_norm_instant(
      rm, tibble::tibble(sample_id = as.character(qd), datetime = qd)
    ))
    expect(is.null(trip), paste(trip, "while resolving 200 dates (API evaluated per date)"))
    if (is.null(trip)) {
      expect_lte(ctr$units, bound)
      expect_identical(nrow(res), length(qd) * length(rm$models))
    }
  })

  it("S2: predicts each tier-2 reference GAM at most once for 200 dates", {
    rm <- .pf_counted_reference_model()
    n_model <- .pf_n_model_tier(rm)
    qd <- seq(as.Date("2022-01-01"), by = "day", length.out = 200L)
    .pf_reset_gam_ctr(limit = n_model)
    withr::defer(.pf_reset_gam_ctr())
    res <- NULL
    trip <- .pf_catch_trip(res <- hydroSense:::.resolve_ref_norm_instant(
      rm, tibble::tibble(sample_id = as.character(qd), datetime = qd)
    ))
    expect(is.null(trip), paste(trip, "while resolving 200 dates (GAM predicted per date)"))
    if (is.null(trip)) {
      expect_lte(.pf_gam_ctr$n, n_model)
      expect_identical(nrow(res), length(qd) * length(rm$models))
    }
  })
})


## ---------------------------------------------------------------------------
## S3: one daily-grid reference resolution per mspaf_daily() call
## ---------------------------------------------------------------------------

describe("mspaf_daily() shares the reference-norm resolution across points", {
  it("S3: resolves each daily-grid date once, not once per point", {
    rm    <- .pf_reference_model()
    tgt   <- .pf_small_target()
    start <- as.Date("2022-03-01")
    end   <- as.Date("2022-04-29")
    grid  <- seq(start, end, by = "day")
    n_pts <- length(unique(tgt$site_id))
    expect_identical(n_pts, 3L)

    # The reference norm depends only on (date, reference_model) — and the
    # reference model carries the single hydrology series, which every point
    # shares — so nothing about it is point-specific.  Allowed work:
    #   * the daily grid, resolved ONCE: length(grid) date-evaluations.  (In
    #     general the per-point date set is the grid PLUS any extra days of
    #     the external temperature series — .fill_external_temperature() does
    #     not clip to start/end — so a shared resolution must cover the union
    #     of the points' date sets; here the temperature series covers
    #     exactly the grid, so the union is the grid.)
    #   * fit_target_model()'s own resolution of each point's grab samples
    #     (one per distinct (sample_id, date) of that point) — fit_target_model
    #     is public and may keep resolving its own anchors.
    # Pre-refactor: each point re-resolves the whole grid (3 x 60) on top of
    # its grabs, i.e. ~ 3 * 60 + 30 = 210 > 60 + 30 = 90.
    n_grab <- nrow(dplyr::distinct(tgt, .data$site_id, .data$sample_id,
                                   d = as.Date(.data$datetime)))
    bound_dates <- length(grid) + n_grab
    # At most one daily-grid call plus one grab call per point (or a single
    # call if the grab dates are folded into the shared resolution).
    bound_calls <- 1L + n_pts

    ctr <- .pf_new_ctr()
    local_mocked_bindings(.resolve_ref_norm_instant = .pf_counting_shim(
      ctr, hydroSense:::.resolve_ref_norm_instant,
      limit  = function() bound_dates,
      weight = function(ref_model, target_df, ...) {
        nrow(dplyr::distinct(tibble::tibble(
          s = target_df$sample_id, d = as.Date(target_df$datetime)
        )))
      },
      what = "reference-norm date-evaluation"
    ))
    out <- NULL
    trip <- .pf_catch_trip(out <- .pf_small_daily(rm, start, end))

    expect(is.null(trip), paste(trip, sprintf(
      "(%d-day grid x %d points: the grid is re-resolved for every point)",
      length(grid), n_pts)))
    if (is.null(trip)) {
      expect_lte(ctr$units, bound_dates)
      expect_lte(ctr$calls, bound_calls)
      expect_setequal(unique(out$site_id), c("P1", "P2", "P3"))
    }
  })
})


## ---------------------------------------------------------------------------
## S4: API / reference-GAM work does not grow with the daily grid
## ---------------------------------------------------------------------------

describe("mspaf_daily() evaluates hydrology in O(1) passes per point", {
  it("S4: .compute_api() calls and reference-GAM predictions do not grow with the grid", {
    rm <- .pf_counted_reference_model()
    start <- as.Date("2022-03-01")

    # Mocks live in a function frame so they are undone when it returns.
    run_counted <- function(ctr, end, limit = function() Inf) {
      local_mocked_bindings(.compute_api = .pf_counting_shim(
        ctr, hydroSense:::.compute_api, limit = limit,
        what = ".compute_api() call"
      ))
      .pf_catch_trip(.pf_small_daily(rm, start, end))
    }

    # Run A: 30-day grid.  Count everything.
    ctr_a <- .pf_new_ctr()
    .pf_reset_gam_ctr()
    withr::defer(.pf_reset_gam_ctr())
    expect_null(run_counted(ctr_a, start + 29L))
    api_a <- ctr_a$units
    gam_a <- .pf_gam_ctr$n

    # Run B: 150-day grid, same grabs.  Must not need more calls than A.
    # Everything outside the reference resolver (fit_target_model()'s tau
    # search, the per-analyte hydro features of .fit_daily_target(), the
    # smoother's hydro z-score) already makes a grid-independent number of
    # .compute_api() calls (54 on this fixture pre-refactor, identical for
    # 40- and 170-day grids); only the per-date resolver loop grows.  After
    # the refactor the whole count is grid-independent, so B == A (B < A is
    # allowed for caching).  Reference-GAM predictions come only from the
    # resolver, so they must not grow either.
    ctr_b <- .pf_new_ctr()
    .pf_reset_gam_ctr(limit = gam_a)
    trip <- run_counted(ctr_b, start + 149L, limit = function() api_a)

    expect(is.null(trip), paste(trip, sprintf(
      "on a 150-day grid; a 30-day grid needed %d .compute_api() calls and %d reference-GAM predictions",
      api_a, gam_a)))
    if (is.null(trip)) {
      expect_lte(ctr_b$units, api_a)
      expect_lte(.pf_gam_ctr$n, gam_a)
    }
  })
})


## ---------------------------------------------------------------------------
## S5 / S6: time_weighted_aggregate() is vectorised over focal dates
## ---------------------------------------------------------------------------

## Three points sampled monthly for 55 years, one analyte (the shape of the
## chronic msPAF step on grab data).
.pf_long_monthly <- function() {
  set.seed(81L)
  d0 <- seq(as.Date("1970-01-15"), as.Date("2024-12-15"), by = "month")
  dplyr::bind_rows(lapply(c("P1", "P2", "P3"), function(s) {
    d <- d0 + sample(-5:5, length(d0), replace = TRUE)
    tibble::tibble(
      sample_id = paste0(s, "_", seq_along(d)), site_id = s, datetime = d,
      analyte = "msPAF", value = stats::rlnorm(length(d), log(20), 0.5),
      detected = TRUE, imputed = stats::runif(length(d)) < 0.1
    )
  }))
}

describe("time_weighted_aggregate() is vectorised over focal dates", {
  it("S5: handles 20,000 focal dates x 3 points without purrr::pmap()", {
    df <- .pf_long_monthly()
    focal <- seq(as.Date("1970-01-01"), by = "day", length.out = 20000L)
    # The pre-refactor code iterates (focal_date x site_id) with
    # purrr::pmap(); 60,000 iterations of dplyr calls is what makes the
    # chronic step slow.  Any per-focal-date loop is forbidden; pmap is the
    # one the old code uses.  (Mock lives in a function frame.)
    without_pmap <- function(expr) {
      local_mocked_bindings(
        pmap = function(...) {
          stop("STRUCTURAL: purrr::pmap() must not be called by time_weighted_aggregate()",
               call. = FALSE)
        },
        .package = "purrr"
      )
      expr
    }
    big <- without_pmap(time_weighted_aggregate(
      df, focal, tau = 90, tau_units = "d", window = 365, window_units = "d",
      summary = "arith_mean"
    ))
    # Every focal date on/after a point's first sample has data in its window
    # (monthly sampling, 365-day window).
    first <- tapply(df$datetime, df$site_id, min)
    n_expected <- sum(vapply(first, function(f) sum(focal >= f), integer(1)))
    expect_identical(nrow(big), n_expected)
    expect_true(all(big$n_samples_in_window >= 1L))
    # The big call agrees with small calls on a sample of focal dates.
    pick <- focal[seq(400L, 20000L, by = 997L)]
    small <- time_weighted_aggregate(df, pick, tau = 90, tau_units = "d",
                                     window = 365, window_units = "d",
                                     summary = "arith_mean")
    sub <- big[big$focal_date %in% pick, , drop = FALSE]
    sub <- sub[order(sub$focal_date, sub$site_id, method = "radix"), ]
    small <- small[order(small$focal_date, small$site_id, method = "radix"), ]
    expect_equal(as.data.frame(sub), as.data.frame(small), tolerance = 1e-10)
  })

  it("S6: dplyr / per-group helper calls do not grow with the number of focal dates", {
    sc <- .pf_scenario()
    tgt <- sc$target
    verbs <- c("filter", "distinct", "arrange", "left_join", "inner_join",
               "group_by", "summarise", "mutate", "select", "bind_rows")
    count_calls <- function(focal) {
      ctrs <- stats::setNames(lapply(verbs, function(v) .pf_new_ctr()), verbs)
      shims <- stats::setNames(lapply(verbs, function(v) {
        .pf_counting_shim(ctrs[[v]], getExportedValue("dplyr", v), what = v)
      }), verbs)
      agg <- .pf_new_ctr()
      ns <- asNamespace("hydroSense")
      has_agg <- exists(".aggregate_weighted", envir = ns, inherits = FALSE)
      run <- function() {
        do.call(local_mocked_bindings, c(shims, list(.package = "dplyr")))
        if (has_agg) {
          local_mocked_bindings(.aggregate_weighted = .pf_counting_shim(
            agg, get(".aggregate_weighted", envir = ns), what = "aggregate"
          ))
        }
        time_weighted_aggregate(tgt, focal, tau = 90, tau_units = "d",
                                window = 365, window_units = "d",
                                summary = "arith_mean")
      }
      run()
      c(vapply(ctrs, function(e) e$calls, integer(1)),
        .aggregate_weighted = agg$calls)
    }
    few  <- count_calls(seq(as.Date("2022-01-01"), by = "week", length.out = 6L))
    many <- count_calls(seq(as.Date("2022-01-01"), by = "week", length.out = 60L))
    # Pre-refactor: ~7 dplyr verbs and 11 .aggregate_weighted() calls (one
    # per analyte) for every (focal date x point), so `many` ~ 10 x `few`.
    # A per-group R-level helper call (e.g. .aggregate_weighted() inside a
    # grouped summarise()) is still a loop over focal dates and is counted.
    # Checked for summary = "arith_mean" (the chronic msPAF summary); the
    # "p90" path is pinned only by the goldens.
    grew <- names(many)[many > few]
    expect(length(grew) == 0L, sprintf(
      "call counts grew with focal dates (6 -> 60): %s",
      paste(sprintf("%s %d -> %d", grew, few[grew], many[grew]), collapse = ", ")
    ))
  })
})
