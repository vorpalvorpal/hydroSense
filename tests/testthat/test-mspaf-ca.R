## Unit tests for compute_ca_group_mspaf() — the within-TMoA Concentration
## Addition combination.
##
## Pins the De Zwart & Posthuma (2005) eq. 6 form:
##   msPAF_CA = Phi( log10(sum TU) / mean(sigma) ),  TU = C_adj / HC50
## and guards against the previous (incorrect) variance-style combination
## sigma_mix = sqrt(sum(w^2 * sigma^2)).

library(testthat)
library(leachatetools)

ca <- function(...) leachatetools:::compute_ca_group_mspaf(tibble::tibble(...))

test_that("a single-component group reduces to that component's own SSD", {
  # TU = 1 (C_adj = HC50) -> log10(1) = 0 -> Phi(0) = 0.5
  expect_equal(
    ca(C_adj = 10, hc50 = 10, sigma = 0.7, moa_group = "g"),
    0.5
  )
  # TU = 10 -> Phi(1 / sigma)
  expect_equal(
    ca(C_adj = 100, hc50 = 10, sigma = 0.7, moa_group = "g"),
    pnorm(log10(10) / 0.7)
  )
})

test_that("a multi-component group uses mean(sigma) and summed TU (eq. 6)", {
  C_adj <- c(50, 30)
  hc50  <- c(10, 20)
  sigma <- c(0.6, 1.0)
  TU    <- C_adj / hc50            # 5, 1.5  -> sum 6.5
  expected <- pnorm(log10(sum(TU)) / mean(sigma))

  got <- ca(C_adj = C_adj, hc50 = hc50, sigma = sigma,
            moa_group = c("g", "g"))
  expect_equal(got, expected)

  # Regression guard: must NOT match the old sqrt(sum(w^2 * sigma^2)) form.
  w        <- TU / sum(TU)
  old_form <- pnorm(log10(sum(TU)) / sqrt(sum(w^2 * sigma^2)))
  expect_false(isTRUE(all.equal(got, old_form)))
})

test_that("zero / invalid groups return 0", {
  expect_equal(ca(C_adj = 0,  hc50 = 10, sigma = 0.7, moa_group = "g"), 0)
  expect_equal(
    ca(C_adj = 10, hc50 = NA_real_, sigma = 0.7, moa_group = "g"),
    0
  )
})
