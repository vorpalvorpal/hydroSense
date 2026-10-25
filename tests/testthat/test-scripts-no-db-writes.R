## ============================================================================
## scripts/run_mspaf_pipeline.R must not write back to a measurement database
## ============================================================================
##
## Contract (TDD): LMF and msPAF are derived products and are never written
## into a live measurement DB. The old helper script wrote interpolated
## temperature and msPAF rows back into the monitoring DuckDB (and recorded
## imputation models there). If the script still exists it must contain no DB
## writes; deleting the script also satisfies this spec.

library(testthat)

.script_path <- function() {
  file.path(testthat::test_path("..", ".."), "scripts", "run_mspaf_pipeline.R")
}

## Code lines with whole-line comments removed (roxygen / prose may mention
## the history; only executable code is checked).
.code_lines <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x[!grepl("^\\s*#", x)]
}

test_that("run_mspaf_pipeline.R has no database write-backs (or is gone)", {
  path <- .script_path()
  if (!file.exists(path)) {
    succeed("scripts/run_mspaf_pipeline.R has been removed")
    return(invisible())
  }
  code <- .code_lines(path)

  forbidden <- c(
    "dbWriteTable"                  = "dbWriteTable",
    "dbAppendTable"                 = "dbAppendTable",
    "INSERT INTO"                   = "(?i)\\bINSERT\\s+(OR\\s+\\w+\\s+)?INTO\\b",
    "DELETE FROM"                   = "(?i)\\bDELETE\\s+FROM\\b",
    "UPDATE ... SET"                = "(?i)\\bUPDATE\\s+\\w+\\s+SET\\b",
    "CREATE TABLE"                  = "(?i)\\bCREATE\\s+(OR\\s+REPLACE\\s+)?TABLE\\b",
    ".write_mspaf_to_db"            = "\\.write_mspaf_to_db",
    ".fill_temperature_gaps"        = "\\.fill_temperature_gaps",
    "read_only = FALSE"             = "read_only\\s*=\\s*FALSE"
  )
  for (nm in names(forbidden)) {
    hits <- grep(forbidden[[nm]], code, perl = TRUE, value = TRUE)
    expect_equal(length(hits), 0L,
                 info = paste0("found ", nm, ": ", trimws(hits[1])))
  }
})
