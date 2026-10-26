## Review findings against the msPAF performance refactor
## ============================================================================
##
## CONTRACT (as in test-perf-goldens.R).  The refactor is a pure performance
## change: results stay identical to the pre-refactor code (commit 86891ab)
## within 1e-10 RELATIVE, with the same NA pattern; keys, counts, tiers and
## other non-double columns exact.  A reviewer reproduced two breaks of that
## contract; this file pins the pre-refactor behaviour for both.
##
##   R1  Batch prediction failure must not demote every date.  When
##       predict.gam() errors (it does on an Inf hydro feature; it returns NA
##       for an NA one), the pre-refactor resolver lost only the failing
##       dates to the static tier, because it predicted one date at a time.
##       A batched predict() whose error demotes the whole batch breaks this.
##       R1c keeps the batching: with no error, at most one predict() per
##       tier-2 model per resolver call; with an error, extra calls are
##       allowed but never more than one per date (per model) on top of the
##       batch.
##   R2  Prefix-sum cancellation in time_weighted_aggregate(summary =
##       "arith_mean").  When earlier values dwarf those in the current
##       window, the difference of two prefix sums loses digits.  These
##       values are legitimately tiny (a 1e-7 msPAF background is a real
##       result), so they are compared PURELY relatively: the 1e-12 absolute
##       floor of test-perf-goldens.R would hide exactly this error.  The
##       geometric-mean path is pinned for the same inputs.
##
## The goldens in _goldens/review_*.rds were produced by
## data-raw/make_review_goldens.R running the PRE-refactor R/ sources on the
## synthetic fixtures of helper-perf-fixtures.R and helper-perf-review.R.
## Never regenerate them with refactored code.

library(testthat)

.pfr_read_golden <- function(name) {
  path <- .pf_golden_path(name)
  if (!file.exists(path)) {
    stop("Golden file missing: ", path,
         " (generate it with data-raw/make_review_goldens.R on the PRE-refactor code)")
  }
  g <- readRDS(path)
  attr(g, "golden_meta") <- NULL
  g
}

## Element-wise PURELY relative comparison (no absolute floor).
.pfr_expect_rel <- function(actual, golden, label, tol = 1e-10) {
  expect_type(actual, "double")
  expect_identical(length(actual), length(golden), label = label)
  if (length(actual) != length(golden)) return(invisible())
  expect_identical(is.na(actual), is.na(golden),
                   label = paste(label, "NA pattern"))
  ok <- !is.na(golden) & !is.na(actual)
  a <- actual[ok]
  g <- golden[ok]
  same <- a == g # also covers equal infinities
  diff  <- abs(a - g)
  scale <- pmax(abs(a), abs(g))
  bad <- !same & !(diff <= tol * scale)
  i <- which(bad)[1L]
  expect(
    !any(bad),
    sprintf("%s: %d of %d value(s) differ by more than %g relative; first at [%d]: actual %.17g, golden %.17g (rel %.3g)",
            label, sum(bad), length(g), tol, which(ok)[i], a[i], g[i],
            diff[i] / scale[i])
  )
}

.pfr_expect_frame <- function(actual, golden, label) {
  actual <- as.data.frame(actual, stringsAsFactors = FALSE)
  golden <- as.data.frame(golden, stringsAsFactors = FALSE)
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
      .pfr_expect_rel(as.numeric(a), as.numeric(g), lab)
    } else {
      expect_identical(class(a), class(g), label = paste(lab, "class"))
      neq <- which(!(a == g) | (is.na(a) != is.na(g)))
      neq <- union(neq, which(is.na(a == g) & !(is.na(a) & is.na(g))))
      expect(identical(a, g), sprintf(
        "%s: %d of %d value(s) differ; first at [%d]: actual %s, golden %s",
        lab, length(neq), length(g), neq[1L], format(a[neq[1L]]),
        format(g[neq[1L]])
      ))
    }
  }
}

.pfr_tier_table <- function(x) {
  t <- table(factor(x$ref_tier, c("matched", "model", "static")))
  paste(names(t), as.integer(t), sep = " = ", collapse = ", ")
}


## ---------------------------------------------------------------------------
## R1: a failing predict() demotes only the failing dates
## ---------------------------------------------------------------------------

describe(".resolve_ref_norm_instant() with an Inf hydro value", {
  golden <- .pfr_read_golden("review_resolver_inf")

  for (ht in c("rainfall", "stage")) {
    local({
      ht <- ht
      it(sprintf("R1 (%s): tiers and ref_norm match the per-date code", ht), {
        res <- hydroSense:::.resolve_ref_norm_instant(
          .pfr_inf_model(ht), .pfr_inf_target()
        )
        res <- .pf_sort_resolver(res)
        g <- golden[[ht]]
        # Fixture sanity: the golden has dates on every tier.
        expect_gt(sum(g$ref_tier == "model"), 0L)
        expect_gt(sum(g$ref_tier == "static" & g$analyte %in% c("Cu", "Zn")), 0L)
        expect_identical(
          .pfr_tier_table(res), .pfr_tier_table(g),
          label = sprintf("tier counts (%s hydrology, Inf at %s)", ht,
                          format(.pfr_inf_date))
        )
        .pfr_expect_frame(res, g, paste0("resolver_inf_", ht))
      })
    })
  }
})

## Per-model predict() counter: GAMs carry an extra S3 class whose method
## counts calls under the analyte's name, then defers to predict.gam().
.pfr_pred_ctr <- new.env()
registerS3method("predict", "pfr_counted_gam", function(object, ...) {
  nm <- attr(object, "pfr_analyte")
  .pfr_pred_ctr[[nm]] <- (.pfr_pred_ctr[[nm]] %||% 0L) + 1L
  NextMethod()
}, envir = asNamespace("stats"))

.pfr_counted <- function(rm) {
  for (nm in names(rm$models)) {
    g <- rm$models[[nm]]$gamm_fit
    if (!is.null(g) && identical(rm$models[[nm]]$tier, "model")) {
      class(g) <- c("pfr_counted_gam", class(g))
      attr(g, "pfr_analyte") <- nm
      rm$models[[nm]]$gamm_fit <- g
    }
  }
  rm
}

.pfr_count_predicts <- function(rm, tdf) {
  rm <- .pfr_counted(rm)
  rm(list = ls(.pfr_pred_ctr), envir = .pfr_pred_ctr)
  hydroSense:::.resolve_ref_norm_instant(rm, tdf)
  nms <- names(rm$models)[vapply(rm$models, function(m) {
    inherits(m$gamm_fit, "pfr_counted_gam")
  }, logical(1))]
  stats::setNames(vapply(nms, function(n) .pfr_pred_ctr[[n]] %||% 0L,
                         integer(1)), nms)
}

describe(".resolve_ref_norm_instant() keeps prediction batched", {
  it("R1c: one predict() per model without errors; at most one more per date with them", {
    tdf <- .pfr_inf_target()
    n_dates <- length(unique(as.Date(tdf$datetime)))

    clean <- .pfr_count_predicts(.pf_reference_model(), tdf)
    expect_gte(length(clean), 1L) # fixture sanity: tier-2 models exist
    expect(all(clean <= 1L), sprintf(
      "finite hydrology: predict() called %s times per model over %d dates (expected at most 1)",
      paste(names(clean), clean, sep = " ", collapse = ", "), n_dates))

    for (ht in c("rainfall", "stage")) {
      bad <- .pfr_count_predicts(.pfr_inf_model(ht), tdf)
      expect(all(bad <= 1L + n_dates), sprintf(
        "Inf in %s hydrology: predict() called %s times per model over %d dates (bound 1 + %d)",
        ht, paste(names(bad), bad, sep = " ", collapse = ", "), n_dates,
        n_dates))
    }
  })
})


## ---------------------------------------------------------------------------
## R2: arithmetic-mean TWA is not degraded by large earlier values
## ---------------------------------------------------------------------------

describe("time_weighted_aggregate() keeps full precision after large values", {
  golden <- .pfr_read_golden("review_prefix_sums")
  inp <- .pfr_prefix_inputs()
  labels <- c(
    a = "(a) 50 then 1e-9, tau = window = 365",
    b = "(b) 1e3 then 1e-12, tau = window = 30",
    c = "(c) U(0,100) then 1e-10*U, tau 90, window 365",
    d = "(d) msPAF-like pulses on a 1e-7 background, tau = window = 365"
  )

  for (nm in names(inp)) {
    local({
      nm <- nm
      x <- inp[[nm]]
      for (s in c("arith_mean", "geom_mean")) {
        local({
          s <- s
          it(sprintf("R2%s %s: %s matches the pre-refactor code", nm, s,
                     labels[[nm]]), {
            res <- suppressMessages(time_weighted_aggregate(
              x$df, x$focal, tau = x$tau, tau_units = "d",
              window = x$window, window_units = "d", summary = s
            ))
            .pfr_expect_frame(res, golden[[paste0(nm, "_", s)]],
                              paste0("prefix_", nm, "_", s))
          })
        })
      }
    })
  }
})
