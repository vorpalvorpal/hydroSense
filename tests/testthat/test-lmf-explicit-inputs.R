## ============================================================================
## add_lmf(): end-member data are explicit inputs — no dashboard fallback
## ============================================================================
##
## Contract (TDD): add_lmf() used to fall back, when leachate_data or
## reference_data was NULL, to host-environment dashboard globals
## (feature_df(), data_df(), feature_sfc(), get_reference_site()). That
## fallback is removed. With either argument NULL (or missing), add_lmf()
## errors with a message naming the missing argument, in both calibration
## modes, and never calls a dashboard global.
##
## The globals are stubbed to fail with a sentinel message so that reaching
## the old fallback cannot masquerade as the required error.

library(testthat)
library(hydroSense)

.ex_data <- function() {
  set.seed(101)
  list(
    ref = .pi_samples(rep("2022-06-01", 12), 0, "ref", "refsite"),
    lea = .pi_samples(rep("2022-06-01", 12), 1, "lea", "leasite"),
    q   = .pi_samples("2022-06-01", 0.5, "q", "query", seed = 7)
  )
}

for (mode in c("window", "trailing")) {
  test_that(paste0("[", mode, "] NULL leachate_data errors naming leachate_data"), {
    d <- .ex_data()
    .pi_stub_dashboard_globals()
    err <- expect_error(
      add_lmf(d$q, leachate_data = NULL, reference_data = d$ref, calibration = mode),
      "leachate_data"
    )
    expect_no_match(conditionMessage(err), "SENTINEL")
  })

  test_that(paste0("[", mode, "] NULL reference_data errors naming reference_data"), {
    d <- .ex_data()
    .pi_stub_dashboard_globals()
    err <- expect_error(
      add_lmf(d$q, leachate_data = d$lea, reference_data = NULL, calibration = mode),
      "reference_data"
    )
    expect_no_match(conditionMessage(err), "SENTINEL")
  })

  test_that(paste0("[", mode, "] omitted end-member arguments error naming them"), {
    d <- .ex_data()
    .pi_stub_dashboard_globals()
    expect_error(add_lmf(d$q, reference_data = d$ref, calibration = mode), "leachate_data")
    expect_error(add_lmf(d$q, leachate_data = d$lea, calibration = mode), "reference_data")
  })
}

test_that("add_lmf with explicit data never calls a dashboard global", {
  d <- .ex_data()
  .pi_stub_dashboard_globals()
  expect_no_error(add_lmf(d$q, leachate_data = d$lea, reference_data = d$ref))
  expect_no_error(add_lmf(d$q, leachate_data = d$lea, reference_data = d$ref,
                          calibration = "trailing"))
})

test_that("the package no longer declares the dashboard globals", {
  ## globalVariables() registrations for host-environment functions are what
  ## kept the fallback quiet under R CMD check; with the fallback gone they
  ## must go too.
  gv <- utils::globalVariables(package = asNamespace("hydroSense"))
  expect_false(any(c("data_df", "feature_df", "feature_sfc",
                     "get_reference_site") %in% gv))
})
