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
    out <- mspaf_pipeline(tgt, ref)
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
    mspaf_pipeline(tgt, ref, impute = FALSE)
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
    mspaf_pipeline(tgt, ref, impute_on = "reference")
    expect_identical(seen$site, "reference")
    mspaf_pipeline(tgt, ref, impute_on = "target")
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
    mspaf_pipeline(tgt, ref, daily_args = list(
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
    mspaf_pipeline(tgt, ref, reference_args = list(min_obs_model = 5L))
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
    expect_warning(mspaf_pipeline(tgt, ref), "Imputation skipped")
  })

  it("attaches the fitted reference and imputation models as attributes", {
    local_mocked_bindings(
      fit_imputation_model = function(df, ...) {
        structure(list(groups = list(m = 1)), class = "imputation_model")
      },
      fit_reference_model = function(reference, ...) "RMODEL",
      mspaf_daily = function(df, ...) fake_daily()
    )
    out <- mspaf_pipeline(tgt, ref)
    expect_identical(attr(out, "reference_model"), "RMODEL")
    expect_s3_class(attr(out, "imputation_model"), "imputation_model")
  })

  it("rejects non-data-frame target or reference", {
    expect_error(mspaf_pipeline("nope", ref), "data")
    expect_error(mspaf_pipeline(tgt, 42), "data")
  })
})
