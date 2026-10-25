## ============================================================================
## hydrosense_provenance(): what computed a derived product
## ============================================================================
##
## Contract (TDD): an exported hydrosense_provenance() returns a named list
## with at least
##   package_version  character(1), == as.character(packageVersion("hydroSense"))
##   git_sha          character(1); NA_character_ allowed (e.g. CRAN install),
##                    otherwise a 7-40 character lowercase hex SHA
##   r_version        character(1), == as.character(getRversion())
##   computed_at      POSIXct(1) in UTC, the time of the call
## so that a derived store (LMF / chronic msPAF) can record what produced it.

library(testthat)
library(hydroSense)

test_that("hydrosense_provenance() is exported", {
  expect_true("hydrosense_provenance" %in% getNamespaceExports("hydroSense"))
})

test_that("hydrosense_provenance() returns the required named fields", {
  t0 <- Sys.time()
  p  <- hydrosense_provenance()
  t1 <- Sys.time()

  expect_type(p, "list")
  expect_true(all(c("package_version", "git_sha", "r_version", "computed_at") %in% names(p)))

  expect_type(p$package_version, "character")
  expect_length(p$package_version, 1L)
  expect_identical(p$package_version, as.character(utils::packageVersion("hydroSense")))

  expect_type(p$git_sha, "character")
  expect_length(p$git_sha, 1L)
  if (!is.na(p$git_sha)) expect_match(p$git_sha, "^[0-9a-f]{7,40}$")

  expect_type(p$r_version, "character")
  expect_length(p$r_version, 1L)
  expect_identical(p$r_version, as.character(getRversion()))

  expect_s3_class(p$computed_at, "POSIXct")
  expect_length(p$computed_at, 1L)
  expect_identical(attr(p$computed_at, "tzone"), "UTC")
  expect_gte(as.numeric(p$computed_at), floor(as.numeric(t0)))
  expect_lte(as.numeric(p$computed_at), ceiling(as.numeric(t1)))
})

test_that("hydrosense_provenance() fields are all named (no unnamed elements)", {
  p <- hydrosense_provenance()
  expect_false(is.null(names(p)))
  expect_true(all(nzchar(names(p))))
})
