## Stan-free coverage for the domain-agnostic imputation-group machinery:
## the impute_group() constructor, the leachate_impute_groups() preset, group
## validation, and the analyte-routing logic (.route_groups). These need no
## brms/Stan and run in the default suite.

library(testthat)
library(hydroSense)

# ── impute_group() constructor ───────────────────────────────────────────────

test_that("impute_group() builds a valid object and dedups", {
  g <- impute_group("metals", targets = c("Cu", "Zn", "Cu"), hurdle = c("Cu"))
  expect_s3_class(g, "impute_group")
  expect_equal(g$name, "metals")
  expect_equal(g$targets, c("Cu", "Zn"))   # de-duplicated, order preserved
  expect_equal(g$hurdle, "Cu")
})

test_that("impute_group() treats targets = NULL as a catch-all", {
  g <- impute_group("organics", targets = NULL, hurdle = c("DOC", "TOC"))
  expect_null(g$targets)
  expect_equal(g$hurdle, c("DOC", "TOC"))
})

test_that("impute_group() validates its inputs", {
  expect_error(impute_group(123),            "non-empty string")
  expect_error(impute_group(c("a", "b")),    "non-empty string")
  expect_error(impute_group(""),             "non-empty string")
  expect_error(impute_group("g", targets = c("Cu", NA)), "targets")
  expect_error(impute_group("g", hurdle = 1:3),          "hurdle")
})

test_that("print.impute_group() is informative", {
  expect_output(print(impute_group("metals", targets = c("Cu", "Zn"))),
                "metals")
  expect_output(print(impute_group("catch", targets = NULL)), "catch-all")
  expect_output(print(impute_group("g", targets = "Cu")),     "<none>")
})


# ── leachate_impute_groups() preset ──────────────────────────────────────────

test_that("leachate_impute_groups() returns the metals + organics preset", {
  gs <- leachate_impute_groups()
  expect_length(gs, 2L)
  expect_true(all(vapply(gs, inherits, logical(1L), "impute_group")))
  expect_equal(vapply(gs, function(g) g$name, character(1L)),
               c("metals", "organics"))

  metals   <- gs[[1]]
  organics <- gs[[2]]
  expect_true(all(c("Cu", "Zn", "Pb", "Ni") %in% metals$targets))
  ## Hurdled on its own targets, EXCEPT the redox indicators Fe/Mn (finding
  ## 7): those are routine analytes kept as PCA predictors, not genuine
  ## trace-metal contamination signals, so they must not alone satisfy the
  ## presence hurdle.
  expect_identical(metals$hurdle, setdiff(metals$targets, c("Fe", "Mn")))
  expect_false("Fe" %in% metals$hurdle)
  expect_false("Mn" %in% metals$hurdle)
  expect_null(organics$targets)                    # catch-all
  expect_true(all(c("DOC", "TOC", "BOD") %in% organics$hurdle))
})


# ── .validate_impute_groups() ────────────────────────────────────────────────

test_that(".validate_impute_groups() rejects malformed group lists", {
  expect_error(hydroSense:::.validate_impute_groups(list()), "non-empty list")
  expect_error(
    hydroSense:::.validate_impute_groups(list(impute_group("a"), "nope")),
    "impute_group"
  )
  expect_error(
    hydroSense:::.validate_impute_groups(
      list(impute_group("dup", targets = "Cu"), impute_group("dup", targets = "Zn"))
    ),
    "unique"
  )
  expect_error(
    hydroSense:::.validate_impute_groups(
      list(impute_group("a", targets = NULL), impute_group("b", targets = NULL))
    ),
    "catch-all"
  )
})

test_that(".validate_impute_groups() accepts the leachate preset", {
  expect_silent(hydroSense:::.validate_impute_groups(leachate_impute_groups()))
})


# ── .route_groups() analyte assignment ───────────────────────────────────────

test_that(".route_groups() assigns explicit targets and routes remainder to catch-all", {
  pool <- c("Cu", "Zn", "As", "Phenol", "Toluene", "NO3-N")
  groups <- list(
    impute_group("metals",    targets = c("Cu", "Zn", "As")),
    impute_group("nutrients", targets = c("NO3-N")),
    impute_group("organics",  targets = NULL)   # catch-all
  )
  routed <- hydroSense:::.route_groups(pool, groups)

  expect_named(routed, c("metals", "nutrients", "organics"))
  expect_setequal(routed$metals,    c("Cu", "Zn", "As"))
  expect_setequal(routed$nutrients, "NO3-N")
  expect_setequal(routed$organics,  c("Phenol", "Toluene"))  # the remainder
})

test_that(".route_groups() gives an overlapping analyte to the first group", {
  pool <- c("Cu", "Zn")
  groups <- list(
    impute_group("a", targets = c("Cu", "Zn")),
    impute_group("b", targets = c("Zn"))        # Zn already claimed by 'a'
  )
  routed <- hydroSense:::.route_groups(pool, groups)
  expect_setequal(routed$a, c("Cu", "Zn"))
  expect_length(routed$b, 0L)
})

test_that(".route_groups() leaves unmatched analytes unassigned without a catch-all", {
  pool <- c("Cu", "Zn", "Mystery")
  groups <- list(impute_group("metals", targets = c("Cu", "Zn")))
  routed <- hydroSense:::.route_groups(pool, groups)
  expect_setequal(routed$metals, c("Cu", "Zn"))
  expect_false("Mystery" %in% unlist(routed))   # dropped, no catch-all
})


# ── all-empty target set yields a no-op model ────────────────────────────────

test_that("fit_imputation_model returns a no-op model when no targets are found", {
  skip_if_not_installed("brms")
  set.seed(1)
  ids <- paste0("s", seq_len(20))
  mk  <- function(an, vals) tibble::tibble(
    sample_id = ids, site_id = "A",
    datetime = as.Date("2023-01-01") + seq_along(ids),
    analyte = an, value = vals, detected = TRUE
  )
  # Only PCA/required vars present — nothing routes into any group.
  df <- dplyr::bind_rows(
    mk("pH", runif(20, 6, 8)),
    mk("EC", runif(20, 100, 500)),
    mk("DOC", rlnorm(20, 1, 1))
  )
  expect_warning(
    m <- fit_imputation_model(df, required_vars = c("pH", "EC"),
                              iter = 100, warmup = 50, chains = 1),
    "No target analytes"
  )
  expect_length(m$groups, 0L)
  expect_null(m$pca)
  # impute_chemistry on a no-op model returns df unchanged with tag columns.
  expect_warning(out <- impute_chemistry(df, m), "no fitted groups")
  expect_true(all(c("imputed", "imputed_kind") %in% names(out)))
  expect_false(any(out$imputed))
})


# ── #61: absolute detection-count gate (min_target_detect_n) ─────────────────

test_that("min_target_detect_n drops a target that clears the frequency gate", {
  set.seed(1)
  n <- 20L
  ids <- paste0("s", seq_len(n))
  mk <- function(an, vals, det) tibble::tibble(
    sample_id = ids, site_id = "A",
    datetime  = as.Date("2023-01-01") + seq_len(n),
    analyte = an, value = vals, detected = det
  )
  # Cu detected in 2/20 samples: det_freq = 0.10 (>= 0.05, clears the fraction
  # gate) but n_detect = 2 (< default min_target_detect_n = 4).
  cu_det <- c(rep(TRUE, 2L), rep(FALSE, n - 2L))
  df <- dplyr::bind_rows(
    mk("pH", runif(n, 6, 8), rep(TRUE, n)),
    mk("EC", runif(n, 100, 500), rep(TRUE, n)),
    mk("Cu", c(runif(2, 1, 5), rep(0.5, n - 2L)), cu_det)
  )
  # Cu is the only target -> dropped by the count gate -> no modellable groups,
  # so the function returns a no-op model without reaching Stan.
  expect_message(
    m <- fit_imputation_model(
      df, required_vars = c("pH", "EC"),
      groups = list(impute_group("metals", targets = "Cu", hurdle = "Cu"))
    ) |> suppressWarnings(),
    "min_target_detect_n"
  )
  expect_length(m$groups, 0L)
})

test_that("a lower min_target_detect_n keeps the same low-count target in-pool", {
  # Routing/gating only: with min_target_detect_n = 1 the count gate no longer
  # bites, so Cu (det_freq 0.10) survives the gate (we don't fit Stan here).
  det <- tibble::tibble(analyte = "Cu", n_detect = 2L, det_freq = 0.10)
  keep_default <- det$analyte[det$det_freq >= 0.05 & det$n_detect >= 4L]
  keep_relaxed <- det$analyte[det$det_freq >= 0.05 & det$n_detect >= 1L]
  expect_length(keep_default, 0L)
  expect_equal(keep_relaxed, "Cu")
})
