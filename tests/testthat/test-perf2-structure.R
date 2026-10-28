## Structural (vectorisation / caching) tests for the second msPAF
## performance pass (perf/vectorise-daily-2)
## ============================================================================
##
## CONTRACT.  Results must stay identical to the pre-refactor code (pinned by
## test-perf2-goldens.R).  THIS file pins the SHAPE of the work, mostly with
## call counts; the one wall-clock test (T1b) compares two sizes with a
## generous bound and an absolute floor, and is skipped on CRAN.
##
##   T1  .compute_daily_diag() is linear in the number of days:
##       a) no per-day purrr::map_dfr()/tibble construction -- the number of
##          purrr / tibble / dplyr calls does not grow with the grid;
##       b) 4x the days costs < 6x the time (the per-day scan of every row
##          makes the pre-refactor code ~O(days^2): > 8x here).
##   T2  .interpolate_analyte() has no R-level loop over target dates (body
##       inspection: no for / while / repeat, apply family or purrr map).
##   T3  .build_daily_chem() subsets the site's rows once, not with a
##       dplyr::filter() + dplyr::arrange() pair per analyte: the number of
##       those calls does not grow with the number of analytes.
##   T4  .pivot_chem_wide() (and so .compute_pca_scores()) reshape without
##       tidyr::pivot_wider() or a grouped dplyr::summarise(): neither is
##       called, and neither appears in the function body.
##   T5  During one fit_target_model() call, the antecedent precipitation
##       index is evaluated at most ONCE per (rain series, tau): every
##       .compute_api() call outside the reference-norm resolver has a
##       distinct (hydro values, hydro dates, tau) key.  Checked for the
##       pooled tau search and for the per-analyte (pool = FALSE) path.
##   T6  The same holds for .fit_daily_target()'s own work after the target
##       model is fitted (per-analyte static hydro features, the smoother's
##       hydro z-score): pooled analytes share one tau pair, so they share
##       one API evaluation per tau.
##   T7  The pooled tau search does not rebuild its per-analyte frames on
##       every optimiser evaluation: the number of tibble / dplyr / purrr
##       calls made by .fit_pooled_impact_response() with auto_select = TRUE
##       (dozens of AIC evaluations) exceeds that with auto_select = FALSE
##       (no search) by at most a small constant.
##
## SEAMS -- what the implementer may and may not rename (in addition to the
## seams of test-perf-structure.R):
##   MUST KEEP: .compute_api(hydro_values, hydro_dates, target_dates, tau) as
##     the single entry point that EVALUATES the API (T5/T6 count it; a cache
##     must sit in front of it -- e.g. evaluate once over every date the fit
##     needs and index into that -- not inside it, and the recursion must not
##     be inlined elsewhere); .resolve_ref_norm_instant() (its own API calls
##     are excluded from T5/T6); fit_target_model(); .fit_daily_target();
##     .fit_pooled_impact_response(anchors_all, target_analytes, hydro,
##     hydro_type, tau_bounds_short, tau_bounds_long, auto_select,
##     min_obs_model, eps) (new optional arguments are fine);
##     .compute_daily_diag(), .interpolate_analyte(), .build_daily_chem(),
##     .pivot_chem_wide(), .compute_pca_scores() with their current
##     signatures; .optimise_tau_pair() and its stats::optimize() search.
##   FREE TO CHANGE: every other private helper; caching state that lives for
##     one fit (an environment created per fit_target_model() /
##     .fit_daily_target() call and passed down); base R or vctrs instead of
##     dplyr/tidyr/purrr inside the named functions.
##
## On the pre-refactor code every test here FAILS (by design) except the body
## inspections that the old code already satisfies (none, at the time of
## writing).  Counting shims "trip" as soon as a bound is exceeded and then
## fail fast.

library(testthat)


## ---------------------------------------------------------------------------
## Counting helpers (as in test-perf-structure.R, prefixed .pf2s_)
## ---------------------------------------------------------------------------

## A trip signals a non-error condition so the package's own
## tryCatch(error = ...) fallbacks cannot swallow it.
.pf2s_trip <- function(msg) {
  cond <- structure(class = c("pf_structural_trip", "condition"),
                    list(message = msg, call = NULL))
  signalCondition(cond)
  stop(msg, call. = FALSE)
}

.pf2s_catch_trip <- function(expr) {
  tryCatch({
    force(expr)
    NULL
  }, pf_structural_trip = function(c) conditionMessage(c))
}

.pf2s_new_ctr <- function() {
  e <- new.env()
  e$calls <- 0L
  e
}

## Counting wrappers for exported functions of `pkg`, mocked for the duration
## of `run()`.  Returns the named call counts.  `limit` (named, optional)
## trips a function's counter once exceeded.
.pf2s_count_pkg_calls <- function(pkg, fns, run, limit = NULL) {
  ns <- asNamespace(pkg)
  fns <- fns[vapply(fns, exists, logical(1), envir = ns, inherits = FALSE)]
  ctrs <- stats::setNames(lapply(fns, function(f) .pf2s_new_ctr()), fns)
  shims <- stats::setNames(lapply(fns, function(f) {
    orig <- get(f, envir = ns)
    ctr <- ctrs[[f]]
    lim <- if (!is.null(limit) && f %in% names(limit)) limit[[f]] else Inf
    force(orig)
    function(...) {
      ctr$calls <- ctr$calls + 1L
      if (ctr$calls > lim) {
        .pf2s_trip(sprintf("STRUCTURAL: %s::%s() called more than %s times",
                           pkg, f, format(lim)))
      }
      orig(...)
    }
  }), fns)
  (function() {
    do.call(testthat::local_mocked_bindings, c(shims, list(.package = pkg)))
    run()
  })()
  vapply(ctrs, function(e) e$calls, integer(1))
}

## Several packages at once: list(pkg = fns).  Mocks nest via a recursive
## function frame so all are active together.
.pf2s_count_calls <- function(spec, run, limits = list()) {
  out <- list()
  go <- function(i) {
    if (i > length(spec)) return(run())
    pkg <- names(spec)[i]
    out[[pkg]] <<- .pf2s_count_pkg_calls(pkg, spec[[i]], function() go(i + 1L),
                                         limit = limits[[pkg]])
  }
  go(1L)
  unlist(unname(lapply(names(out), function(p) {
    stats::setNames(out[[p]], paste0(p, "::", names(out[[p]])))
  })))
}

.pf2s_spec <- list(
  purrr = c("map", "map_dfr", "map_dfc", "map_dbl", "map_int", "map_lgl",
            "map_chr", "map2", "map2_dfr", "pmap", "pmap_dfr", "imap"),
  tibble = c("tibble", "as_tibble", "tribble"),
  dplyr = c("filter", "arrange", "mutate", "select", "summarise", "group_by",
            "bind_rows", "bind_cols", "left_join", "inner_join", "distinct",
            "count", "if_else", "first")
)

## Function names called anywhere in a function's body (including nested
## function definitions); `pkg::fn` is reported as "pkg::fn" and as "fn".
.pf2s_called <- function(fn) {
  found <- character()
  walk <- function(e) {
    if (is.call(e)) {
      h <- e[[1L]]
      if (is.name(h)) {
        found <<- c(found, as.character(h))
      } else if (is.call(h) && identical(h[[1L]], as.name("::")) ||
                 is.call(h) && identical(h[[1L]], as.name(":::"))) {
        found <<- c(found, paste0(as.character(h[[2L]]), "::",
                                  as.character(h[[3L]])),
                    as.character(h[[3L]]))
      }
      for (a in as.list(e)) if (!missing(a)) walk(a)
    } else if (is.pairlist(e)) {
      for (a in as.list(e)) if (!missing(a)) walk(a)
    }
  }
  walk(body(fn))
  unique(found)
}

.pf2s_fn <- function(name) get(name, envir = asNamespace("hydroSense"))


## ---------------------------------------------------------------------------
## T1: .compute_daily_diag() is linear in the number of days
## ---------------------------------------------------------------------------

describe(".compute_daily_diag() is linear in the number of days", {
  it("T1a: makes no per-day purrr / tibble / dplyr calls", {
    tox <- paste0("T", 1:4)
    small <- .pf2_diag_long(40L, k = 4L)
    big   <- .pf2_diag_long(400L, k = 4L)
    run <- function(x) function() {
      hydroSense:::.compute_daily_diag(x, tox, "S")
    }
    n_small <- .pf2s_count_calls(.pf2s_spec, run(small))
    n_big   <- .pf2s_count_calls(.pf2s_spec, run(big))
    grew <- names(n_big)[n_big > n_small]
    expect(length(grew) == 0L, sprintf(
      "call counts grew with the number of days (40 -> 400): %s",
      paste(sprintf("%s %d -> %d", grew, n_small[grew], n_big[grew]),
            collapse = ", ")))
    expect_false("map_dfr" %in% .pf2s_called(.pf2s_fn(".compute_daily_diag")),
                 label = ".compute_daily_diag() body calls purrr::map_dfr()")
  })

  it("T1b: 4x the days costs < 6x the time", {
    skip_on_cran()
    ## Many tox rows per day so the pre-refactor per-day scan of every row
    ## (O(days x rows) = O(days^2)) dominates its per-day overhead.
    k <- 150L
    tox <- paste0("T", seq_len(k))
    small <- .pf2_diag_long(600L, k = k)
    big   <- .pf2_diag_long(2400L, k = k)
    f <- hydroSense:::.compute_daily_diag
    time_it <- function(x, reps) {
      min(vapply(seq_len(reps), function(i) {
        system.time(f(x, tox, "S"))[["elapsed"]]
      }, numeric(1)))
    }
    invisible(f(small, tox, "S")) # warm up
    t_small <- max(time_it(small, 3L), 1e-3)
    t_big <- time_it(big, 1L)
    if (t_big < 1) t_big <- min(t_big, time_it(big, 2L))
    ratio <- t_big / t_small
    ## The absolute floor keeps a fast (linear) implementation from failing
    ## on timer noise when both runs take a few milliseconds.
    expect(ratio < 6 || t_big < 0.1, sprintf(
      "600 -> 2400 days took %.3f s -> %.3f s (%.1fx); linear work would be ~4x",
      t_small, t_big, ratio))
  })
})


## ---------------------------------------------------------------------------
## T2: .interpolate_analyte() has no per-date R loop
## ---------------------------------------------------------------------------

describe(".interpolate_analyte() is vectorised over target dates", {
  it("T2: its body has no for / while / repeat, apply-family or purrr loop", {
    called <- .pf2s_called(.pf2s_fn(".interpolate_analyte"))
    loops <- c("for", "while", "repeat", "vapply", "sapply", "lapply",
               "mapply", "Map", "apply", "Vectorize", "map", "map_dbl",
               "map_lgl", "map_dfr", "pmap", "map2")
    hit <- intersect(loops, called)
    expect(length(hit) == 0L, sprintf(
      ".interpolate_analyte() loops over target dates in R: %s",
      paste(hit, collapse = ", ")))
  })
})


## ---------------------------------------------------------------------------
## T3: .build_daily_chem() subsets the site rows once
## ---------------------------------------------------------------------------

describe(".build_daily_chem() does not filter + arrange once per analyte", {
  it("T3: dplyr::filter()/arrange() calls do not grow with the analytes", {
    mk <- function(n_an) {
      set.seed(94L)
      d <- seq(as.Date("2022-01-03"), by = "20 days", length.out = 12L)
      an <- paste0("X", seq_len(n_an))
      g <- expand.grid(datetime = d, analyte = an, KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)
      tibble::tibble(sample_id = paste0("s", as.integer(g$datetime)),
                     site_id = "S", datetime = g$datetime,
                     analyte = g$analyte, value = stats::runif(nrow(g), 1, 9),
                     detected = TRUE, units.analyte = "ug/L")
    }
    grid <- seq(as.Date("2022-01-01"), as.Date("2022-09-30"), by = "day")
    run <- function(x) function() {
      hydroSense:::.build_daily_chem(x, grid, "forward_fill", "drop",
                                     c("X1", "X2"))
    }
    spec <- list(dplyr = c("filter", "arrange"))
    few  <- .pf2s_count_calls(spec, run(mk(3L)))
    many <- .pf2s_count_calls(spec, run(mk(30L)))
    grew <- names(many)[many > few]
    expect(length(grew) == 0L, sprintf(
      "per-analyte subsetting (3 -> 30 analytes): %s",
      paste(sprintf("%s %d -> %d", grew, few[grew], many[grew]),
            collapse = ", ")))
  })
})


## ---------------------------------------------------------------------------
## T4: .pivot_chem_wide() without pivot_wider / grouped summarise
## ---------------------------------------------------------------------------

describe(".pivot_chem_wide() reshapes without tidyr::pivot_wider()", {
  it("T4a: neither pivot_wider() nor summarise() is called", {
    x <- .pf2_pivot_input()
    pca <- .pf2_imputation_model()$pca
    wq <- c("pH", "EC", "DOC", "Ca", "Mg", "hardness", "temperature", "NH3-N")
    spec <- list(tidyr = c("pivot_wider", "pivot_longer", "spread"),
                 dplyr = c("summarise", "summarize", "group_by"))
    limits <- list(tidyr = c(pivot_wider = 0, pivot_longer = 0, spread = 0),
                   dplyr = c(summarise = 0, summarize = 0, group_by = 0))
    trip <- .pf2s_catch_trip(.pf2s_count_calls(spec, function() {
      hydroSense:::.pivot_chem_wide(x, wq, hydroSense:::.PCA_NO_LOG_VARS)
      hydroSense:::.compute_pca_scores(x, pca)
    }, limits = limits))
    expect(is.null(trip), paste(trip, "inside .pivot_chem_wide() / .compute_pca_scores()"))
  })

  it("T4b: the body of .pivot_chem_wide() names neither", {
    called <- .pf2s_called(.pf2s_fn(".pivot_chem_wide"))
    hit <- intersect(c("tidyr::pivot_wider", "dplyr::summarise",
                       "dplyr::summarize", "pivot_wider"), called)
    expect(length(hit) == 0L, paste(".pivot_chem_wide() calls:",
                                    paste(hit, collapse = ", ")))
  })
})


## ---------------------------------------------------------------------------
## T5 / T6: one API evaluation per (rain series, tau) per fit
## ---------------------------------------------------------------------------

## Count .compute_api() calls by (hydro values, hydro dates, tau) key,
## skipping calls made while .resolve_ref_norm_instant() is on the stack (and,
## if `skip_fit`, while fit_target_model() is).  Trips on the first repeat.
.pf2s_api_once <- function(run, skip_fit = FALSE) {
  st <- new.env()
  st$seen <- character()
  st$depth <- 0L
  st$n <- 0L
  api0 <- hydroSense:::.compute_api
  res0 <- hydroSense:::.resolve_ref_norm_instant
  ftm0 <- hydroSense::fit_target_model
  guard <- function(orig) {
    force(orig)
    function(...) {
      st$depth <- st$depth + 1L
      on.exit(st$depth <- st$depth - 1L)
      orig(...)
    }
  }
  api_shim <- function(hydro_values, hydro_dates, target_dates, tau, ...) {
    if (st$depth == 0L) {
      key <- paste(rlang::hash(list(as.numeric(hydro_values),
                                    as.numeric(as.Date(hydro_dates)))),
                   sprintf("%.17g", tau))
      st$n <- st$n + 1L
      if (key %in% st$seen) {
        .pf2s_trip(sprintf(
          "STRUCTURAL: .compute_api() evaluated twice for the same rain series and tau = %.6g (after %d distinct evaluations)",
          tau, length(st$seen)))
      }
      st$seen <- c(st$seen, key)
    }
    api0(hydro_values, hydro_dates, target_dates, tau, ...)
  }
  mocks <- list(.compute_api = api_shim,
                .resolve_ref_norm_instant = guard(res0))
  if (skip_fit) mocks$fit_target_model <- guard(ftm0)
  trip <- (function() {
    do.call(testthat::local_mocked_bindings,
            c(mocks, list(.package = "hydroSense")))
    .pf2s_catch_trip(run())
  })()
  list(trip = trip, n = st$n, distinct = length(st$seen))
}

describe("the API is evaluated once per (rain series, tau) during a fit", {
  it("T5a: pooled fit_target_model() (tau search over a factor-smooth GAM)", {
    rm <- .pf_reference_model()
    p1 <- .pf2_p1()
    r <- .pf2s_api_once(function() {
      suppressWarnings(suppressMessages(fit_target_model(
        p1, rm, api_tau_bounds_short = c(3, 10),
        api_tau_bounds_long = c(30, 90)
      )))
    })
    expect(is.null(r$trip), paste(r$trip,
      "-- the pooled tau search recomputes the API per analyte and per evaluation"))
    if (is.null(r$trip)) expect_gt(r$n, 0L)
  })

  it("T5b: per-analyte fit_target_model() (pool = FALSE)", {
    rm <- .pf_reference_model()
    p1 <- .pf2_p1()
    r <- .pf2s_api_once(function() {
      suppressWarnings(suppressMessages(fit_target_model(
        p1, rm, pool = FALSE, api_tau_bounds_short = c(3, 10),
        api_tau_bounds_long = c(30, 90)
      )))
    })
    expect(is.null(r$trip), paste(r$trip,
      "-- the per-analyte tau searches recompute the API at the same tau"))
    if (is.null(r$trip)) expect_gt(r$n, 0L)
  })

  it("T6: .fit_daily_target() shares one API per tau across pooled analytes", {
    inp <- .pf2_daily_model_inputs()
    rm <- .pf_reference_model()
    r <- .pf2s_api_once(function() {
      fdm <- suppressWarnings(suppressMessages(hydroSense:::.fit_daily_target(
        site_rows = inp$site_rows, reference_model = rm,
        imputation_model = NULL, conc_units = NULL, meta = inp$meta,
        tox_analytes = inp$tox, daily_long = inp$daily_long
      )))
      expect_false(is.null(fdm))
      expect_true(any(vapply(fdm$tm$models, function(m) isTRUE(m$pooled),
                             logical(1))))
    }, skip_fit = TRUE)
    expect(is.null(r$trip), paste(r$trip,
      "-- per-analyte static features / smoother z-scores recompute the API"))
    if (is.null(r$trip)) expect_gt(r$n, 0L)
  })
})


## ---------------------------------------------------------------------------
## T7: the pooled tau search does not rebuild frames per evaluation
## ---------------------------------------------------------------------------

## Synthetic pooled anchors: three analytes on 24 shared dates whose g
## responds to the 6-day and 60-day API (so the pooled fit is "useful").
.pf2s_pool_anchors <- function() {
  sc <- .pf_scenario()
  set.seed(95L)
  d <- sort(sample(seq(as.Date("2021-02-01"), as.Date("2023-10-01"),
                       by = "day"), 24L))
  hs <- .pf2_ref_api(sc$rain$value, sc$rain$date, d, 6)
  hl <- .pf2_ref_api(sc$rain$value, sc$rain$date, d, 60)
  do.call(rbind, lapply(c("Cu", "Zn", "Ni"), function(a) {
    g <- 0.02 * hs + 0.004 * hl + stats::rnorm(length(d), 0, 0.2)
    tibble::tibble(date = d, analyte = a, I = g, scale_c = NA_real_, g = g)
  }))
}

describe("the pooled tau search reuses its per-analyte frames", {
  it("T7: tibble / dplyr / purrr calls do not grow with optimiser evaluations", {
    sc <- .pf_scenario()
    anch <- .pf2s_pool_anchors()
    hydro <- tibble::tibble(date = sc$rain$date, value = sc$rain$value)
    run <- function(auto) function() {
      suppressWarnings(hydroSense:::.fit_pooled_impact_response(
        anch, c("Cu", "Zn", "Ni"), hydro, "rainfall",
        c(3, 10), c(30, 90), auto, 12L, 1e-9
      ))
    }
    ## Baseline: no tau search (one fit at the default taus).
    base <- .pf2s_count_calls(.pf2s_spec, run(FALSE))
    ## With the search: dozens of AIC evaluations.  Allow 2 extra calls per
    ## counted function (one-off setup), trip beyond that.
    limits <- stats::setNames(lapply(names(.pf2s_spec), function(p) {
      fns <- .pf2s_spec[[p]]
      b <- unname(base[paste0(p, "::", fns)])
      stats::setNames(ifelse(is.na(b), 0, b) + 2, fns)
    }), names(.pf2s_spec))
    trip <- .pf2s_catch_trip(.pf2s_count_calls(.pf2s_spec, run(TRUE),
                                               limits = limits))
    expect(is.null(trip), paste(trip,
      sprintf("(baseline without the tau search: %s)",
              paste(sprintf("%s %d", names(base)[base > 0], base[base > 0]),
                    collapse = ", "))))
  })
})
