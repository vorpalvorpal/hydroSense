## Generate the golden outputs pinned by tests/testthat/test-perf2-goldens.R
##
## Purpose: characterisation ("golden master") fixtures for the second msPAF
## performance pass (branch perf/vectorise-daily-2): .compute_daily_diag(),
## .interpolate_analyte(), .build_daily_chem(), .pivot_chem_wide(),
## .compute_pca_scores(), .compute_hydro_features(), the fit_target_model()
## tau search (pooled and per-analyte, with the .optimise_tau_pair() AIC
## trace), .fit_daily_target() / .predict_daily_tox() /
## .daily_tox_from_model(), and mspaf_pipeline() in point mode.  The goldens
## record what the PRE-REFACTOR code returns on fully synthetic inputs, so the
## vectorised code can be held to "identical within 1e-10 relative".
##
## They were generated from origin/main at bdb0740 (R/ sources unchanged).
## Each golden records its source commit in attr(, "golden_meta").  DO NOT
## regenerate them with refactored code: that would make the tests compare
## the new code with itself.  Regenerate only if the fixtures
## (tests/testthat/helper-perf2-fixtures.R) change on purpose, and then from a
## checkout of the pre-refactor R/ sources.
##
## The fit-level and pipeline goldens are platform-specific (see PLATFORM in
## test-perf2-goldens.R); generate them on that platform.
##
## All inputs are synthetic (fixed seeds, invented site codes); no monitoring
## data is read.
##
## Usage (from the package root, with R/ committed; does not install):
##   Rscript data-raw/make_perf2_goldens.R

pkg <- normalizePath(".")
dirty <- length(system2("git", c("-C", pkg, "status", "--porcelain", "R/"),
                        stdout = TRUE)) > 0L
if (dirty) stop("R/ has uncommitted changes; commit them before regenerating")
suppressMessages(devtools::load_all(pkg, quiet = TRUE))
library(testthat)
source(file.path(pkg, "tests", "testthat", "helper-perf-fixtures.R"))
source(file.path(pkg, "tests", "testthat", "helper-perf2-fixtures.R"))

out_dir <- file.path(pkg, "tests", "testthat", "_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp <- list(
  generated_from = system2("git", c("-C", pkg, "rev-parse", "--short", "HEAD"),
                           stdout = TRUE),
  r_sources_last_commit = system2("git", c("-C", pkg, "log", "-1",
                                           "--format=%h", "--", "R/"),
                                  stdout = TRUE),
  r_sources_dirty = FALSE,
  r_version = R.version.string,
  platform = R.version$platform,
  generated_on = as.character(Sys.Date())
)

save_golden <- function(x, name) {
  attr(x, "golden_meta") <- stamp
  saveRDS(x, file.path(out_dir, paste0("perf2_", name, ".rds")), version = 3)
  message("wrote perf2_", name)
}

## Same collation testthat imposes on every test (LC_COLLATE = "C").
withr::local_collate("C")

t0 <- Sys.time()
save_golden(.pf2_diag_cases(), "daily_diag")
save_golden(.pf2_interp_cases(), "interpolate_analyte")
save_golden(.pf2_build_chem_cases(), "build_daily_chem")
save_golden(.pf2_pivot_cases(), "pivot_chem_wide")
save_golden(.pf2_pca_cases(), "pca_scores")
save_golden(.pf2_hydro_cases(), "hydro_features")
save_golden(.pf2_fit_cases(), "fit_target_model")
save_golden(.pf2_daily_model_cases(), "daily_model")
save_golden(.pf2_pipeline_point_cases(), "pipeline_point")
message("done in ", format(Sys.time() - t0))
