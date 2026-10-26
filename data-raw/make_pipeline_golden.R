## Regenerate the end-to-end mspaf_pipeline() golden after a deliberate model
## change (first used for the msPAF stability fix: fsc basis for the pooled
## tau GAM, REML and a PC gate for the marginal imputation GAMs).
##
## Only _goldens/mspaf_pipeline_chronic.rds is written. Every other golden
## still comes from the pre-refactor code (data-raw/make_perf_goldens.R) and
## must not be regenerated. The pipeline golden is platform-specific; see the
## PLATFORM note in tests/testthat/test-perf-goldens.R.
##
## Usage (from the package root, on the golden platform, with R/ committed):
##   Rscript data-raw/make_pipeline_golden.R

pkg <- normalizePath(".")
dirty <- length(system2("git", c("-C", pkg, "status", "--porcelain", "R/"),
                        stdout = TRUE)) > 0L
if (dirty) stop("R/ has uncommitted changes; commit them before regenerating")
suppressMessages(devtools::load_all(pkg, quiet = TRUE))
library(testthat)
source(file.path(pkg, "tests", "testthat", "helper-perf-fixtures.R"))

stamp <- list(
  generated_from = system2("git", c("-C", pkg, "rev-parse", "--short", "HEAD"),
                           stdout = TRUE),
  r_sources_last_commit = system2("git", c("-C", pkg, "log", "-1", "--format=%h",
                                           "--", "R/"), stdout = TRUE),
  r_sources_dirty = FALSE,
  r_version = R.version.string, generated_on = as.character(Sys.Date()),
  reason = "msPAF stability fix (fsc pooled basis; REML + PC gate imputation)"
)

withr::local_collate("C")
x <- .pf_pipeline_case()
attr(x, "golden_meta") <- stamp
saveRDS(x, file.path(pkg, "tests", "testthat", "_goldens",
                     "mspaf_pipeline_chronic.rds"), version = 3)
message("wrote mspaf_pipeline_chronic")
