## Characterisation ("golden master") tests for the msPAF performance refactor
## ============================================================================
##
## CONTRACT.  The vectorisation of the three chronic-msPAF hot spots
##   1. .resolve_ref_norm_instant()  (reference-norm resolution),
##   2. .compute_api()               (antecedent precipitation index),
##   3. time_weighted_aggregate()    (chronic aggregation),
## and of their callers (mspaf_daily(), mspaf_pipeline()) is a PURE
## PERFORMANCE change.  Every result below must stay identical to the
## pre-refactor code (commit 86891ab):
##   * numeric columns within 1e-10 RELATIVE, element by element (an absolute
##     floor of 1e-12 applies to values at or near zero), with the same NA
##     pattern;
##   * keys, dates, counts (n_samples_in_window, n_imputed_in_window, ...),
##     character columns, column names/order, classes and row order EXACT.
##   * Exception: the row order of .resolve_ref_norm_instant() is not part of
##     its contract (every consumer joins or looks up by key), so that golden
##     is compared after sorting by (sample_id, analyte).
##
## The goldens in _goldens/*.rds were produced by data-raw/make_perf_goldens.R
## running the PRE-refactor code on the synthetic fixtures defined in
## helper-perf-fixtures.R.  Never regenerate them with refactored code.
##
## Structural (vectorisation) expectations live in test-perf-structure.R.
##
## PLATFORM.  The component goldens are portable. The end-to-end
## mspaf_pipeline() golden is not: the target model chooses tau by optimising
## a fit score that jumps (a REML/GCV optimum switch), so ulp-level BLAS or
## libm differences between platforms change the chosen tau and move msPAF
## by ~0.5%. That golden therefore runs only on the platform that generated
## it (.pf_golden_platform below) and is skipped elsewhere, including CI.

library(testthat)

.pf_read_golden <- function(name) {
  path <- .pf_golden_path(name)
  if (!file.exists(path)) {
    stop("Golden file missing: ", path,
         " (generate it with data-raw/make_perf_goldens.R on the PRE-refactor code)")
  }
  g <- readRDS(path)
  attr(g, "golden_meta") <- NULL
  g
}

## Element-wise relative comparison of two numeric vectors.
.pf_expect_num <- function(actual, golden, label) {
  expect_type(actual, "double")
  expect_identical(length(actual), length(golden), label = label)
  expect_identical(is.na(actual), is.na(golden),
                   label = paste(label, "NA pattern"))
  ok <- !is.na(golden) & !is.na(actual)
  diff  <- abs(actual[ok] - golden[ok])
  scale <- pmax(abs(golden[ok]), abs(actual[ok]))
  bad <- diff > 1e-10 * scale & diff > 1e-12
  expect(
    !any(bad),
    sprintf("%s: %d value(s) differ by more than 1e-10 relative; first at [%d]: actual %.17g, golden %.17g",
            label, sum(bad), which(ok)[which(bad)[1L]],
            actual[ok][which(bad)[1L]], golden[ok][which(bad)[1L]])
  )
}

## Full data-frame comparison (see CONTRACT above).
.pf_expect_frame <- function(actual, golden, label) {
  expect_identical(class(actual), class(golden), label = paste(label, "class"))
  expect_identical(names(actual), names(golden), label = paste(label, "names"))
  expect_identical(nrow(actual), nrow(golden), label = paste(label, "nrow"))
  if (!identical(names(actual), names(golden)) ||
      !identical(nrow(actual), nrow(golden))) {
    return(invisible())
  }
  for (nm in names(golden)) {
    a <- actual[[nm]]
    g <- golden[[nm]]
    lab <- paste0(label, "$", nm)
    if (is.double(g) && !inherits(g, c("Date", "POSIXt"))) {
      .pf_expect_num(as.numeric(a), as.numeric(g), lab)
      expect_identical(class(a), class(g), label = paste(lab, "class"))
    } else {
      expect_identical(a, g, label = lab)
    }
  }
}


describe(".compute_api() golden", {
  golden <- .pf_read_golden("compute_api")
  actual <- .pf_api_cases()

  it("covers the same cases", {
    expect_identical(names(actual), names(golden))
  })
  for (nm in names(golden)) {
    local({
      nm <- nm
      it(paste("matches the pre-refactor values:", nm), {
        .pf_expect_num(actual[[nm]], golden[[nm]], nm)
      })
    })
  }
  it("still rejects a non-positive or non-scalar tau", {
    sc <- .pf_scenario()
    api <- hydroSense:::.compute_api
    expect_error(api(sc$rain$value, sc$rain$date, as.Date("2022-01-01"), 0))
    expect_error(api(sc$rain$value, sc$rain$date, as.Date("2022-01-01"), c(1, 2)))
  })
})


describe(".resolve_ref_norm_instant() golden", {
  golden <- .pf_read_golden("resolve_ref_norm_instant")
  actual <- lapply(.pf_resolver_cases(), .pf_sort_resolver)

  it("covers matched, model and static tiers", {
    expect_setequal(unique(golden$grid$ref_tier), c("matched", "model", "static"))
  })
  for (nm in names(golden)) {
    local({
      nm <- nm
      it(paste("matches the pre-refactor reference norms:", nm), {
        .pf_expect_frame(actual[[nm]], golden[[nm]], nm)
      })
    })
  }
})


describe("time_weighted_aggregate() golden", {
  golden <- .pf_read_golden("time_weighted_aggregate")
  actual <- .pf_twa_cases()

  it("covers the same cases", {
    expect_identical(names(actual), names(golden))
  })
  for (nm in names(golden)) {
    local({
      nm <- nm
      it(paste("matches the pre-refactor chronic values:", nm), {
        .pf_expect_frame(actual[[nm]], golden[[nm]], nm)
      })
    })
  }
  it("still aborts when no (focal_date, site_id) has in-window data", {
    sc <- .pf_scenario()
    expect_error(
      time_weighted_aggregate(sc$target, as.Date("2000-01-01"),
                              tau = 90, tau_units = "d",
                              window = 365, window_units = "d"),
      "No values could be aggregated"
    )
  })
})


describe("mspaf_daily() point-mode golden", {
  golden <- .pf_read_golden("mspaf_daily_point")
  actual <- .pf_daily_cases()

  for (case in names(golden)) {
    for (part in names(golden[[case]])) {
      local({
        case <- case
        part <- part
        it(paste("matches the pre-refactor output:", case, part), {
          expect_true(part %in% names(actual[[case]]))
          .pf_expect_frame(actual[[case]][[part]], golden[[case]][[part]],
                           paste(case, part))
        })
      })
    }
  }
})


describe("mspaf_pipeline() chronic golden", {
  .pf_skip_unless_golden_platform()
  golden <- .pf_read_golden("mspaf_pipeline_chronic")
  actual <- .pf_pipeline_case()

  for (part in names(golden)) {
    local({
      part <- part
      it(paste("matches the pre-refactor pipeline output:", part), {
        .pf_expect_frame(actual[[part]], golden[[part]], part)
      })
    })
  }
})
