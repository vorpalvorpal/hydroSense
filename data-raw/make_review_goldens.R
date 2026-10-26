## Generate the golden outputs pinned by tests/testthat/test-perf-review.R
##
## Purpose: pin two review findings against the msPAF performance refactor
## (branch perf/mspaf-vectorise) to what the PRE-REFACTOR code returns:
##   * review_resolver_inf.rds : .resolve_ref_norm_instant() when one hydro
##     value is Inf (predict.gam() errors on some dates only);
##   * review_prefix_sums.rds  : time_weighted_aggregate() on series whose
##     earlier values dwarf the current window's (prefix-sum cancellation).
##
## The R/ code is loaded from a checkout of the pre-refactor commit 86891ab;
## the fixtures (tests/testthat/helper-perf-fixtures.R and
## helper-perf-review.R) and the output directory are taken from THIS
## checkout.  Each golden records its provenance in attr(, "golden_meta").
## DO NOT run it against refactored R/ sources.
##
## All inputs are synthetic (fixed seeds, invented site codes).
##
## Usage (from the package root; does not install the package):
##   git worktree add --detach <old> 86891ab
##   Rscript data-raw/make_review_goldens.R <old>
##   git worktree remove <old>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("usage: Rscript data-raw/make_review_goldens.R <pre-refactor checkout>")
}
old_pkg <- normalizePath(args[1])
pkg <- normalizePath(".")

old_sha <- system2("git", c("-C", old_pkg, "rev-parse", "--short", "HEAD"),
                   stdout = TRUE)
if (!identical(old_sha, "86891ab")) {
  stop("expected a checkout of 86891ab, got ", old_sha)
}
dirty <- length(system2("git", c("-C", old_pkg, "status", "--porcelain", "R/"),
                        stdout = TRUE)) > 0L
if (dirty) stop("R/ in ", old_pkg, " has local changes")

suppressMessages(devtools::load_all(old_pkg, quiet = TRUE))
library(testthat)
source(file.path(pkg, "tests", "testthat", "helper-perf-fixtures.R"))
source(file.path(pkg, "tests", "testthat", "helper-perf-review.R"))

out_dir <- file.path(pkg, "tests", "testthat", "_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- list(
  r_sources_commit = old_sha, r_sources_dirty = dirty,
  fixtures_from = system2("git", c("-C", pkg, "rev-parse", "--short", "HEAD"),
                          stdout = TRUE),
  r_version = R.version.string, generated_on = as.character(Sys.Date())
)

save_golden <- function(x, name) {
  attr(x, "golden_meta") <- stamp
  saveRDS(x, file.path(out_dir, paste0(name, ".rds")), version = 3)
  message("wrote ", name)
}

## Same collation as the test run (testthat sets LC_COLLATE = "C").
withr::local_collate("C")

t0 <- Sys.time()
save_golden(lapply(.pfr_resolver_inf_cases(), .pf_sort_resolver),
            "review_resolver_inf")
save_golden(.pfr_prefix_cases(), "review_prefix_sums")
message("done in ", format(Sys.time() - t0))
