## Characterisation ("golden master") tests for the second msPAF performance
## pass (perf/vectorise-daily-2)
## ============================================================================
##
## CONTRACT.  The vectorisation of the point-mode daily hot spots
##   .compute_daily_diag()      (O(days^2) per-day scan -> linear),
##   .interpolate_analyte() / .build_daily_chem()   (per-date / per-analyte
##                                                    loops),
##   .pivot_chem_wide() / .compute_pca_scores()     (pivot_wider + grouped
##                                                    summarise),
##   .compute_hydro_features() / .compute_api()     (API recomputed per
##                                                    analyte and per optimiser
##                                                    evaluation),
##   fit_target_model()'s tau search (anchors_for / build_df rebuilt per
##                                     evaluation),
##   .fit_daily_target() / .predict_daily_tox() / .daily_tox_from_model(),
## is a PURE PERFORMANCE change.  Every result below must stay identical to
## the pre-refactor code (origin/main bdb0740):
##   * numeric (double) values within 1e-10 RELATIVE, element by element (an
##     absolute floor of 1e-12 applies to values at or near zero), with the
##     same NA pattern;
##   * keys, dates, counts, character / logical / integer columns, names,
##     column names AND order, classes and row order EXACT.
##   * The tau optimisation keeps the same optimiser: the sequence of
##     (tau_short, tau_long) pairs that .optimise_tau_pair() evaluates, and
##     the AIC at each, are pinned (tau_trace), as is the chosen tau.
##     Caching hydro features may only remove repeated work, never change
##     what is evaluated.
##
## SEAMS the goldens call directly (keep their signatures; new optional
## arguments are fine): .compute_daily_diag(daily_long, tox_analytes, site),
## .interpolate_analyte(obs_dates, obs_values, obs_detected, target_dates,
## interpolation, leading_edge, log_space), .build_daily_chem(site_rows,
## dates, interpolation, leading_edge, tox_analytes), .pivot_chem_wide(df,
## wq_vars, no_log_vars), .compute_pca_scores(df, pca_obj),
## .compute_hydro_features(hydro, target_dates, tau_short, tau_long,
## hydro_type), .optimise_tau_pair(aic_fn, tau_bounds_short,
## tau_bounds_long, ...) (the tau search must go through it with an aic_fn),
## .resolve_target_impact(target_model, query), .fit_daily_target(...),
## .predict_daily_tox(fdm), .daily_tox_from_model(...), and the fields of a
## target_model / fdm scaffold read by .pf2_tm_view() / .pf2_fdm_view().
##
## PUBLIC OUTPUT.  mspaf_pipeline() point mode (chronic = FALSE, impute =
## TRUE) is compared with and without a static ARA reference.  Concurrent
## feature branches may ADD columns to the public output, so there (only)
## every golden column must be present, in the same relative order, and
## match; extra columns are ignored.
##
## PLATFORM.  The pure-data goldens are portable.  Anything that fits a
## target model is not: the tau search optimises a REML/AIC score that can
## jump under ulp-level BLAS/libm differences (see test-perf-goldens.R), so
## those goldens run only on the platform that generated them
## (.pf_golden_platform: aarch64-apple-darwin20, R 4.4.3, reference BLAS)
## and are skipped elsewhere, including CI.
##
## The goldens in _goldens/perf2_*.rds come from data-raw/make_perf2_goldens.R
## run on the PRE-refactor code with the fixtures of helper-perf2-fixtures.R.
## Never regenerate them with refactored code.
##
## Structural (vectorisation / caching) expectations live in
## test-perf2-structure.R.

library(testthat)

.pf2_read_golden <- function(name) {
  path <- .pf2_golden_path(name)
  if (!file.exists(path)) {
    stop("Golden file missing: ", path,
         " (generate it with data-raw/make_perf2_goldens.R on the PRE-refactor code)")
  }
  g <- readRDS(path)
  attr(g, "golden_meta") <- NULL
  g
}

## Element-wise relative comparison of two double vectors.
.pf2_expect_num <- function(actual, golden, label, rel = 1e-10) {
  expect_type(actual, "double")
  expect_identical(length(actual), length(golden), label = paste(label, "length"))
  if (length(actual) != length(golden)) return(invisible())
  expect_identical(names(actual), names(golden), label = paste(label, "names"))
  expect_identical(is.na(actual), is.na(golden),
                   label = paste(label, "NA pattern"))
  ok <- !is.na(golden) & !is.na(actual)
  a <- as.numeric(actual)[ok]
  g <- as.numeric(golden)[ok]
  same  <- a == g # covers equal infinities
  diff  <- abs(a - g)
  scale <- pmax(abs(a), abs(g))
  bad <- !same & !(diff <= rel * scale | diff <= 1e-12)
  i <- which(bad)[1L]
  expect(
    !any(bad),
    sprintf("%s: %d of %d value(s) differ by more than %g relative; first at [%d]: actual %.17g, golden %.17g",
            label, sum(bad), length(g), rel, which(ok)[i], a[i], g[i])
  )
}

## Recursive comparison following the CONTRACT.  `extra_cols = TRUE` (public
## output only) lets `actual` carry columns the golden lacks.
.pf2_expect_same <- function(actual, golden, label, rel = 1e-10,
                             extra_cols = FALSE) {
  if (is.data.frame(golden)) {
    expect_s3_class(actual, "data.frame")
    if (!is.data.frame(actual)) return(invisible())
    expect_identical(class(actual), class(golden), label = paste(label, "class"))
    if (extra_cols) {
      missing <- setdiff(names(golden), names(actual))
      expect(length(missing) == 0L,
             sprintf("%s: golden column(s) missing: %s", label,
                     paste(missing, collapse = ", ")))
      expect_identical(intersect(names(actual), names(golden)), names(golden),
                       label = paste(label, "column order"))
      if (length(missing)) return(invisible())
    } else {
      expect_identical(names(actual), names(golden), label = paste(label, "names"))
      if (!identical(names(actual), names(golden))) return(invisible())
    }
    expect_identical(nrow(actual), nrow(golden), label = paste(label, "nrow"))
    if (!identical(nrow(actual), nrow(golden))) return(invisible())
    for (nm in names(golden)) {
      .pf2_expect_same(actual[[nm]], golden[[nm]], paste0(label, "$", nm),
                       rel = rel)
    }
    return(invisible())
  }
  if (is.list(golden)) {
    expect_type(actual, "list")
    expect_identical(names(actual), names(golden), label = paste(label, "names"))
    expect_identical(length(actual), length(golden), label = paste(label, "length"))
    if (!identical(names(actual), names(golden)) ||
        length(actual) != length(golden)) {
      return(invisible())
    }
    for (i in seq_along(golden)) {
      .pf2_expect_same(actual[[i]], golden[[i]],
                       paste0(label, "[[", names(golden)[i] %||% i, "]]"),
                       rel = rel, extra_cols = extra_cols)
    }
    return(invisible())
  }
  if (is.double(golden) && !inherits(golden, c("Date", "POSIXt", "difftime"))) {
    expect_identical(class(actual), class(golden), label = paste(label, "class"))
    expect_identical(attributes(unname(actual)), attributes(unname(golden)),
                     label = paste(label, "attributes"))
    .pf2_expect_num(actual, golden, label, rel = rel)
    return(invisible())
  }
  expect_identical(actual, golden, label = label)
}

## One it() per case of a golden set.
.pf2_check_set <- function(golden, actual, what, extra_cols = FALSE) {
  it(paste(what, "covers the same cases"), {
    expect_identical(names(actual), names(golden))
  })
  for (nm in names(golden)) {
    local({
      nm <- nm
      it(paste(what, "matches the pre-refactor output:", nm), {
        expect_true(nm %in% names(actual))
        .pf2_expect_same(actual[[nm]], golden[[nm]], paste(what, nm),
                         extra_cols = extra_cols)
      })
    })
  }
}


## ---------------------------------------------------------------------------
## Pure-data goldens (portable)
## ---------------------------------------------------------------------------

describe(".compute_daily_diag() golden", {
  .pf2_check_set(.pf2_read_golden("daily_diag"), .pf2_diag_cases(),
                 ".compute_daily_diag()")
})

describe(".interpolate_analyte() golden", {
  .pf2_check_set(.pf2_read_golden("interpolate_analyte"), .pf2_interp_cases(),
                 ".interpolate_analyte()")
})

describe(".build_daily_chem() golden", {
  .pf2_check_set(.pf2_read_golden("build_daily_chem"), .pf2_build_chem_cases(),
                 ".build_daily_chem()")
})

describe(".pivot_chem_wide() golden", {
  .pf2_check_set(.pf2_read_golden("pivot_chem_wide"), .pf2_pivot_cases(),
                 ".pivot_chem_wide()")
})

describe(".compute_pca_scores() golden", {
  .pf2_check_set(.pf2_read_golden("pca_scores"), .pf2_pca_cases(),
                 ".compute_pca_scores()")
})

describe(".compute_hydro_features() golden", {
  .pf2_check_set(.pf2_read_golden("hydro_features"), .pf2_hydro_cases(),
                 ".compute_hydro_features()")
})


## ---------------------------------------------------------------------------
## API caching must be keyed on the rain series AND tau AND the target dates
## ---------------------------------------------------------------------------

describe("hydro features never come from a stale cache", {
  ## A cache keyed only on tau (or on the rain series' length/dates) would
  ## hand back the first series' API for the second.  Each call must agree
  ## with an independent reference evaluation (.pf2_ref_api()).
  it("depends on the rain values, the tau and the target dates", {
    sc  <- .pf_scenario()
    chf <- hydroSense:::.compute_hydro_features
    api <- .pf2_ref_api
    rain_a <- sc$rain
    rain_b <- sc$rain
    rain_b$value <- rev(rain_b$value) # same dates, different values
    td1 <- seq(as.Date("2021-01-01"), as.Date("2021-12-31"), by = "7 days")
    td2 <- td1 + 3L
    calls <- list(
      list(rain_a, td1, 5, 60), list(rain_b, td1, 5, 60),
      list(rain_a, td2, 5, 60), list(rain_a, td1, 6, 60),
      list(rain_b, td1, 5, 61), list(rain_a, td1, 5, 60)
    )
    for (k in seq_along(calls)) {
      cl <- calls[[k]]
      got <- chf(cl[[1]], cl[[2]], cl[[3]], cl[[4]], "rainfall")
      expect_identical(got$date, cl[[2]])
      .pf2_expect_num(got$hydro_short,
                      api(cl[[1]]$value, cl[[1]]$date, cl[[2]], cl[[3]]),
                      paste("call", k, "hydro_short"))
      .pf2_expect_num(got$hydro_long,
                      api(cl[[1]]$value, cl[[1]]$date, cl[[2]], cl[[4]]),
                      paste("call", k, "hydro_long"))
    }
  })
})


## ---------------------------------------------------------------------------
## Model-fitting goldens (platform-specific)
## ---------------------------------------------------------------------------

describe("fit_target_model() golden (tau search, pooled / per-analyte / WQ)", {
  .pf_skip_unless_golden_platform()
  .pf2_check_set(.pf2_read_golden("fit_target_model"), .pf2_fit_cases(),
                 "fit_target_model()")
})

describe(".fit_daily_target() / .predict_daily_tox() / .daily_tox_from_model() golden", {
  .pf_skip_unless_golden_platform()
  .pf2_check_set(.pf2_read_golden("daily_model"), .pf2_daily_model_cases(),
                 "daily model")
})

describe("mspaf_pipeline() point-mode golden (impute = TRUE)", {
  .pf_skip_unless_golden_platform()
  .pf2_check_set(.pf2_read_golden("pipeline_point"),
                 .pf2_pipeline_point_cases(), "mspaf_pipeline()",
                 extra_cols = TRUE)
})


describe("fit_target_model() does not reuse hydro features across fits", {
  ## Runs after the fit goldens above, so any cache that outlived those fits
  ## would be visible here.  The same target is refitted on a rain series
  ## shifted by a year; every model-tier anchor's hydro features must be the
  ## API of THAT series at the chosen tau.
  it("anchors carry the API of the series they were fitted on", {
    .pf_skip_unless_golden_platform()
    sc <- .pf_scenario()
    rain2 <- sc$rain
    rain2$value <- c(rain2$value[-(1:365)], rain2$value[1:365])
    tm2 <- suppressWarnings(suppressMessages(fit_target_model(
      .pf2_p1(), .pf_reference_model(), hydro = rain2, pool = FALSE,
      api_tau_bounds_short = c(3, 10), api_tau_bounds_long = c(30, 90)
    )))
    api <- .pf2_ref_api
    golden <- .pf2_read_golden("fit_target_model")$per_analyte
    model_an <- names(Filter(function(m) m$tier == "model", tm2$models))
    expect_gt(length(model_an), 0L)
    for (nm in model_an) {
      m <- tm2$models[[nm]]
      .pf2_expect_num(m$anchors$hydro_short,
                      api(rain2$value, rain2$date, m$anchors$date, m$tau_short),
                      paste(nm, "hydro_short"))
      .pf2_expect_num(m$anchors$hydro_long,
                      api(rain2$value, rain2$date, m$anchors$date, m$tau_long),
                      paste(nm, "hydro_long"))
      g <- golden[[paste0("anchors:", nm)]]
      if (!is.null(g)) {
        expect_false(isTRUE(all.equal(m$anchors$hydro_short, g$hydro_short)),
                     label = paste(nm, "differs from the original-rain fit"))
      }
    }
  })
})
