#!/usr/bin/env Rscript
# Benchmark .compute_api() — issue #49 (recursive linear reservoir).
#
# The pre-#49 implementation was a per-target truncated windowed sum
# (O(n_hydro * n_targets)); #49 replaces it with the exact recursive reservoir
# evaluated as a single C-level stats::filter() pass on a daily grid
# (O(n_hydro + n_targets)), with no truncation horizon.  This benchmark records
# the new cost across realistic hydro-series lengths and target-set sizes.
# Outside R CMD check (lives under bench/).
#
# Usage: Rscript bench/bench-compute-api.R [label]
suppressMessages({
  library(bench)
  devtools::load_all(".", quiet = TRUE)
})

label <- (commandArgs(trailingOnly = TRUE)[1]) %||% "after"
set.seed(49L)

api <- hydroSense:::.compute_api

grid <- expand.grid(
  n_hydro   = c(2000L, 5000L, 10000L),
  n_target  = c(30L, 200L, 1000L),
  tau       = c(7, 180)
)

results <- purrr::pmap_dfr(grid, function(n_hydro, n_target, tau) {
  hd <- seq(as.Date("2000-01-01"), by = "day", length.out = n_hydro)
  hv <- pmax(0, rnorm(n_hydro, 2, 4))
  td <- hd[round(seq(1, n_hydro, length.out = n_target))]
  b <- bench::mark(api(hv, hd, td, tau), iterations = 25, check = FALSE)
  tibble::tibble(
    n_hydro = n_hydro, n_target = n_target, tau = tau,
    median = b$median, mem_alloc = b$mem_alloc
  )
})

print(results)
saveRDS(results, sprintf("bench/results-compute-api-%s.rds", label))
