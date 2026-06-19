#!/usr/bin/env Rscript
# Benchmark impute_chemistry() over a sparse fixture with many entirely-absent
# target cells (issue #53).
#
# The #53 fix makes the merge fabricate new rows for eligible sample x target
# cells that are absent from the input, instead of discarding their (already
# computed) predictions. The predictions themselves are unchanged, so the only
# new cost is dplyr join/bind work. The plan flags one risk: the draws-mode
# join must expand per-draw solely from `pm_long` and must NOT introduce a
# many-to-many blow-up. This benchmark records the wall-clock + memory cost of
# impute_chemistry() in point and draws mode on a deliberately sparse frame so a
# regression in that claim is visible when comparing `before` vs `after`.
#
# Outside R CMD check (lives under bench/).
# Usage: Rscript bench/bench-impute-fabricate.R [label]   (label: before|after)
suppressMessages({
  library(bench)
  devtools::load_all(".", quiet = TRUE)
})

label <- (commandArgs(trailingOnly = TRUE)[1]) %||% "after"
set.seed(53L)

# Sparse fixture: every sample carries the drivers (pH/EC/NH3-N/DOC) so all are
# eligible, but each metal is present in only ~half of samples. The absent
# metal cells are exactly what the #53 fix fabricates, so they dominate the
# merge's new workload.
make_sparse_chem <- function(n) {
  samples <- paste0("s", seq_len(n))
  metals <- c("Cu", "Zn", "Ni")
  drivers <- c("pH", "EC", "NH3-N", "DOC")

  driver_rows <- tidyr::expand_grid(sample_id = samples, analyte = drivers) |>
    dplyr::mutate(
      value = dplyr::case_when(
        analyte == "pH" ~ runif(dplyr::n(), 6.5, 8.5),
        analyte == "EC" ~ runif(dplyr::n(), 100, 500),
        analyte == "NH3-N" ~ runif(dplyr::n(), 0.01, 0.5),
        TRUE ~ runif(dplyr::n(), 0.2, 5)
      )
    )
  # Each metal present for a random ~50% of samples (the rest are fabricated).
  metal_rows <- purrr::map_dfr(metals, function(m) {
    keep <- samples[as.logical(rbinom(n, 1L, 0.5))]
    tibble::tibble(
      sample_id = keep, analyte = m,
      value = exp(rnorm(length(keep), log(2), 0.5))
    )
  })

  dplyr::bind_rows(driver_rows, metal_rows) |>
    dplyr::mutate(
      site_id = "f1",
      datetime = as.Date("2023-01-01") + (match(sample_id, samples) - 1L),
      detected = TRUE
    )
}

# One realistic-but-fast size; fit the model once (fitting is brms-dominated and
# is NOT what changed), then time only the impute_chemistry() calls.
n <- 150L
df <- make_sparse_chem(n)
model <- suppressMessages(fit_imputation_model(
  df,
  required_vars = c("pH", "EC"),
  iter = 500, warmup = 250, chains = 1, cores = 1
))

results <- bench::mark(
  point = suppressMessages(
    impute_chemistry(df, model, return = "point")
  ),
  draws = suppressMessages(
    impute_chemistry(df, model, return = "draws", ndraws = 50L)
  ),
  iterations = 5, check = FALSE
)

print(results[, c("expression", "median", "mem_alloc")])
saveRDS(results, sprintf("bench/results-impute-fabricate-%s.rds", label))
