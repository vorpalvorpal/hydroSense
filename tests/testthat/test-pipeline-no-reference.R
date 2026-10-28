## mspaf_pipeline() without a reference site (reference = NULL)
## ============================================================================
##
## CONTRACT.
##
## A downstream caller must run msPAF for monitoring points that have NO
## reference site of their own water type, and must not borrow one.
##
## N1. mspaf_pipeline(target, reference = NULL, ...) is allowed.
##     * fit_reference_model() is never called;
##     * `ara` defaults to "none"; ara = "static" or "contemporaneous" with
##       reference = NULL is a classed error `hydroSense_ara_needs_reference`;
##     * mspaf_daily() is called with interpolation = "linear" unless
##       daily_args$interpolation says otherwise;
##     * daily_args$interpolation = "model" is a classed error
##       `hydroSense_model_needs_reference`.
## N2. Imputation: with impute = TRUE the model is fitted on the target.  A
##     defaulted impute_on becomes "target"; an explicit
##     impute_on = "reference" is a classed error
##     `hydroSense_impute_needs_reference`.
## N3. Output has the structure of a normal point-mode run (chronic frame,
##     `daily` attribute, `reference_model` attribute NULL, `imputation_model`
##     attribute present when impute = TRUE) and equals a manual replication:
##     fit_imputation_model(target) (if impute); mspaf_daily(target,
##     interpolation = "linear", reference = NULL, reference_model = NULL,
##     imputation_model = ...); then the chronic aggregation of the daily
##     `mspaf` line as .chronic_from_daily() does (tol 1e-10).
## N4. also_total = TRUE with reference = NULL: value_total equals value.
## N5. Every error above is raised before any model is fitted (neither
##     fit_reference_model nor fit_imputation_model is called).
## N6. Behaviour with a reference is unchanged: a normal run equals a run on
##     the same inputs with the arguments spelled out (impute_on, ara
##     explicit), and still fits the reference model once.
##
## All data are synthetic (helper-perf-fixtures.R).

.nr_sc <- function() .pf_scenario()

.nr_target <- function() {
  tgt <- .pf_small_target()
  tgt[tgt$site_id %in% c("P1", "P2"), , drop = FALSE]
}

.nr_tmp <- function() {
  tmp <- .nr_sc()$temperature
  tmp[tmp$datetime >= as.Date("2021-12-15") &
        tmp$datetime <= as.Date("2022-03-10"), ]
}

.nr_daily_base <- function() {
  list(temperature = .nr_tmp(),
       start = as.Date("2022-01-01"), end = as.Date("2022-02-28"))
}

.nr_focal <- as.Date(c("2022-01-20", "2022-02-10", "2022-02-28"))

.nr_ref_args <- function() {
  list(min_obs_model = 20L, api_tau_bounds_short = c(3, 10),
       api_tau_bounds_long = c(30, 90))
}

## Run mspaf_pipeline(); counts model fits.  Returns list(out, err, n_ref,
## n_imp, daily_calls).
.nr_run <- function(..., daily_args = list(), keep_attrs = TRUE) {
  counts <- new.env()
  counts$ref <- 0L; counts$imp <- 0L; counts$daily <- list()
  real_ref <- hydroSense:::fit_reference_model
  real_imp <- hydroSense:::fit_imputation_model
  real_daily <- hydroSense:::mspaf_daily
  testthat::local_mocked_bindings(
    fit_reference_model = function(...) {
      counts$ref <- counts$ref + 1L; real_ref(...)
    },
    fit_imputation_model = function(...) {
      counts$imp <- counts$imp + 1L; real_imp(...)
    },
    mspaf_daily = function(...) {
      counts$daily[[length(counts$daily) + 1L]] <- list(...)[
        c("interpolation", "reference", "reference_model")]
      real_daily(...)
    },
    .package = "hydroSense"
  )
  args <- list(
    target = .nr_target(), hydro = .nr_sc()$rain, impute = FALSE,
    daily_args = utils::modifyList(.nr_daily_base(), daily_args),
    chronic = TRUE, focal_dates = .nr_focal, tau = 90, window = 365
  )
  args <- utils::modifyList(args, list(...), keep.null = TRUE)
  err <- NULL
  out <- tryCatch(
    suppressWarnings(suppressMessages(do.call(mspaf_pipeline, args))),
    error = function(e) { err <<- e; NULL }
  )
  if (!is.null(out) && !keep_attrs) {
    attributes(out)[c("daily", "reference_model", "imputation_model")] <- NULL
  }
  list(out = out, err = err, n_ref = counts$ref, n_imp = counts$imp,
       daily_calls = counts$daily)
}

.nr_cache <- new.env()
.nr_get <- function(key, f) {
  if (is.null(.nr_cache[[key]])) .nr_cache[[key]] <- f()
  .nr_cache[[key]]
}

.nr_noimp <- function() .nr_get("noimp", function()
  .nr_run(reference = NULL, also_total = TRUE))
.nr_imp <- function() .nr_get("imp", function()
  .nr_run(reference = NULL, impute = TRUE))

.nr_chronic <- function(daily, ycol = "mspaf") {
  long <- tibble::tibble(
    sample_id = paste0("daily_", daily$site_id, "_", daily$date),
    site_id = daily$site_id, datetime = daily$date,
    analyte = "msPAF", detected = TRUE, value = daily[[ycol]]
  )
  time_weighted_aggregate(
    long, focal_dates = .nr_focal, tau = 90, tau_units = "d",
    window = 365, window_units = "d", summary = "arith_mean",
    return = "summary"
  )
}

.nr_manual <- function(impute) {
  imodel <- NULL
  if (impute) {
    imodel <- suppressWarnings(suppressMessages(
      fit_imputation_model(.nr_target(), required_vars = c("pH", "EC"),
                           groups = NULL)))
  }
  args <- c(list(.nr_target(), interpolation = "linear", reference = NULL,
                 reference_model = NULL, imputation_model = imodel),
            .nr_daily_base())
  daily <- suppressWarnings(suppressMessages(do.call(mspaf_daily, args)))
  list(daily = daily, chronic = .nr_chronic(daily), imodel = imodel)
}

.nr_key <- function(x, cols = "value") {
  x <- as.data.frame(x)[, c("focal_date", "site_id", cols), drop = FALSE]
  x <- x[order(x$site_id, x$focal_date, method = "radix"), , drop = FALSE]
  rownames(x) <- NULL
  x
}

describe("N1: reference = NULL", {
  skip_on_cran()

  it("runs without fitting a reference model", {
    r <- .nr_noimp()
    expect_null(r$err)
    expect_identical(r$n_ref, 0L)
  })

  it("calls mspaf_daily() with linear interpolation and no reference", {
    r <- .nr_noimp()
    expect_null(r$err)
    expect_length(r$daily_calls, 1L)
    expect_identical(r$daily_calls[[1]]$interpolation, "linear")
    expect_null(r$daily_calls[[1]]$reference)
    expect_null(r$daily_calls[[1]]$reference_model)
  })

  it("honours daily_args$interpolation = 'forward_fill'", {
    r <- .nr_run(reference = NULL,
                 daily_args = list(interpolation = "forward_fill"))
    expect_null(r$err)
    expect_identical(r$daily_calls[[1]]$interpolation, "forward_fill")
  })

  it("rejects ara = 'static' / 'contemporaneous' with a classed error", {
    for (a in c("static", "contemporaneous")) {
      r <- .nr_run(reference = NULL, ara = a)
      expect_s3_class(r$err, "hydroSense_ara_needs_reference")
    }
  })

  it("accepts ara = 'none' explicitly", {
    r <- .nr_run(reference = NULL, ara = "none")
    expect_null(r$err)
  })

  it("rejects daily_args$interpolation = 'model' with a classed error", {
    r <- .nr_run(reference = NULL, daily_args = list(interpolation = "model"))
    expect_s3_class(r$err, "hydroSense_model_needs_reference")
  })
})

describe("N2: imputation without a reference", {
  skip_on_cran()

  it("defaulted impute_on becomes 'target': fitted on the target", {
    r <- .nr_imp()
    expect_null(r$err)
    expect_identical(r$n_imp, 1L)
    expect_identical(r$n_ref, 0L)
    expect_false(is.null(attr(r$out, "imputation_model")))
  })

  it("explicit impute_on = 'reference' is a classed error", {
    r <- .nr_run(reference = NULL, impute = TRUE, impute_on = "reference")
    expect_s3_class(r$err, "hydroSense_impute_needs_reference")
  })

  it("explicit impute_on = 'target' is accepted", {
    r <- .nr_run(reference = NULL, impute = TRUE, impute_on = "target")
    expect_null(r$err)
  })

  it("impute_on = 'reference' is fine when impute = FALSE", {
    r <- .nr_run(reference = NULL, impute = FALSE, impute_on = "reference")
    expect_null(r$err)
  })
})

describe("N3: output structure and equality", {
  skip_on_cran()

  it("has the point-mode structure, a NULL reference_model attribute", {
    out <- .nr_noimp()$out
    expect_s3_class(out, "data.frame")
    expect_true(all(c("focal_date", "site_id", "value") %in% names(out)))
    expect_s3_class(attr(out, "daily"), "data.frame")
    expect_true("mspaf" %in% names(attr(out, "daily")))
    expect_true("reference_model" %in% names(attributes(out)) ||
                  is.null(attr(out, "reference_model")))
    expect_null(attr(out, "reference_model"))
    expect_null(attr(out, "imputation_model"))
  })

  it("has the same columns as a normal run with a reference", {
    normal <- .nr_get("normal", function() .nr_run(
      reference = .nr_sc()$reference, reference_args = .nr_ref_args(),
      also_total = TRUE))
    expect_null(normal$err)
    expect_identical(names(.nr_noimp()$out), names(normal$out))
  })

  it("equals the manual replication (no imputation)", {
    m <- .nr_manual(FALSE)
    out <- .nr_noimp()$out
    expect_equal(.nr_key(out), .nr_key(m$chronic), tolerance = 1e-10)
    d <- function(x) {
      x <- as.data.frame(x)[, c("date", "site_id", "mspaf")]
      x[order(x$site_id, x$date, method = "radix"), ]
    }
    expect_equal(d(attr(out, "daily")), d(m$daily), tolerance = 1e-10,
                 check.attributes = FALSE)
  })

  it("equals the manual replication (imputation on the target)", {
    m <- .nr_manual(TRUE)
    out <- .nr_imp()$out
    expect_false(is.null(m$imodel))
    expect_equal(.nr_key(out), .nr_key(m$chronic), tolerance = 1e-10)
  })
})

describe("N4: also_total without a reference", {
  skip_on_cran()

  it("value_total equals value", {
    out <- .nr_noimp()$out
    expect_true("value_total" %in% names(out))
    expect_equal(out$value_total, out$value, tolerance = 1e-10)
  })
})

describe("N5: errors precede any model fitting", {
  skip_on_cran()

  it("fits nothing when the arguments are rejected", {
    cases <- list(
      list(args = list(ara = "static", impute = TRUE),
           cls = "hydroSense_ara_needs_reference"),
      list(args = list(ara = "contemporaneous", impute = TRUE),
           cls = "hydroSense_ara_needs_reference"),
      list(args = list(impute = TRUE, impute_on = "reference"),
           cls = "hydroSense_impute_needs_reference"),
      list(args = list(impute = TRUE,
                       daily_args = list(interpolation = "model")),
           cls = "hydroSense_model_needs_reference")
    )
    for (cs in cases) {
      r <- do.call(.nr_run, c(list(reference = NULL), cs$args))
      expect_s3_class(r$err, cs$cls)
      expect_identical(r$n_ref, 0L)
      expect_identical(r$n_imp, 0L)
      expect_length(r$daily_calls, 0L)
    }
  })
})

describe("N6: behaviour with a reference is unchanged", {
  skip_on_cran()

  it("a default run equals one with impute_on / ara spelled out", {
    base <- .nr_run(reference = .nr_sc()$reference,
                    reference_args = .nr_ref_args(), keep_attrs = FALSE)
    spelled <- .nr_run(reference = .nr_sc()$reference,
                       reference_args = .nr_ref_args(),
                       impute_on = "reference", ara = "none",
                       keep_attrs = FALSE)
    expect_null(base$err)
    expect_null(spelled$err)
    expect_identical(base$n_ref, 1L)
    expect_identical(base$daily_calls[[1]]$interpolation, "model")
    expect_equal(.nr_key(base$out), .nr_key(spelled$out), tolerance = 1e-10)
  })

  it("impute_on stays 'reference' by default when a reference is given", {
    r <- .nr_run(reference = .nr_sc()$reference, impute = TRUE,
                 reference_args = .nr_ref_args())
    expect_null(r$err)
    expect_identical(r$n_imp, 1L)
    expect_identical(r$n_ref, 1L)
  })
})
