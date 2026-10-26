## Generate the golden outputs pinned by tests/testthat/test-perf-goldens.R
##
## Purpose: characterisation ("golden master") fixtures for the msPAF
## performance refactor (branch perf/mspaf-vectorise).  The goldens record what
## the PRE-REFACTOR code returns on fully synthetic inputs, so the vectorised
## code can be held to "identical within 1e-10 relative".
##
## The goldens were generated with the R/ sources of commit 86891ab
## (origin/main before the refactor; R/ last changed in 0ddb0b2).  Each golden
## records this in attr(, "golden_meta").  DO NOT regenerate them with refactored code: that would make
## the tests compare the new code with itself.  Regenerate only if the
## fixtures (tests/testthat/helper-perf-fixtures.R) change on purpose, and then
## from a checkout of the pre-refactor R/ sources.
##
## All inputs are synthetic (fixed seeds, invented site codes); no monitoring
## data is read.
##
## Usage (from the package root; does not install the package):
##   Rscript data-raw/make_perf_goldens.R

pkg <- normalizePath(".")
suppressMessages(devtools::load_all(pkg, quiet = TRUE))
library(testthat)
source(file.path(pkg, "tests", "testthat", "helper-perf-fixtures.R"))

out_dir <- file.path(pkg, "tests", "testthat", "_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sha <- tryCatch(
  system2("git", c("-C", pkg, "rev-parse", "--short", "HEAD"), stdout = TRUE),
  error = function(e) NA_character_
)
r_sha <- tryCatch(
  system2("git", c("-C", pkg, "log", "-1", "--format=%h", "--", "R/"),
          stdout = TRUE),
  error = function(e) NA_character_
)
dirty <- tryCatch(
  length(system2("git", c("-C", pkg, "status", "--porcelain", "R/"),
                 stdout = TRUE)) > 0L,
  error = function(e) NA
)
stamp <- list(
  generated_from = sha, r_sources_last_commit = r_sha,
  r_sources_dirty = dirty,
  r_version = R.version.string, generated_on = as.character(Sys.Date())
)

save_golden <- function(x, name) {
  attr(x, "golden_meta") <- stamp
  saveRDS(x, file.path(out_dir, paste0(name, ".rds")), version = 3)
  message("wrote ", name)
}

## Generate under the same collation testthat imposes on every test
## (testthat::local_reproducible_output() sets LC_COLLATE = "C"), so any
## locale-sensitive ordering inside the code under test matches the test run.
withr::local_collate("C")

t0 <- Sys.time()
save_golden(.pf_api_cases(), "compute_api")
save_golden(lapply(.pf_resolver_cases(), .pf_sort_resolver), "resolve_ref_norm_instant")
save_golden(.pf_twa_cases(), "time_weighted_aggregate")
save_golden(.pf_daily_cases(), "mspaf_daily_point")
save_golden(.pf_pipeline_case(), "mspaf_pipeline_chronic")
message("done in ", format(Sys.time() - t0))
