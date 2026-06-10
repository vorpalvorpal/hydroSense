#!/usr/bin/env Rscript
# Benchmark add_amspaf() draws-mode scaling — issue #30 (vectorise the engine).
# No prior benchmark coverage existed for this engine; this fills that gap and
# provides the before/after baseline. Outside R CMD check (lives under bench/).
#
# Usage: Rscript bench/bench-amspaf.R [label]   (label tags the saved results,
#        e.g. "baseline" before the change, "after" once vectorised)
suppressMessages({ library(bench); devtools::load_all(".", quiet = TRUE) })

label <- (commandArgs(trailingOnly = TRUE)[1]) %||% "baseline"
set.seed(30L)

CO <- c(pH = 7.5, DOC = 2.0, Ca = 6.0, Mg = 4.0, hardness = 30.0)
.METALS <- c("Cu", "Zn", "Ni", "Cd", "Pb", "As", "Cr", "Mn", "Co", "Se",
             "Ag", "Al", "Fe", "Hg", "B")

# Build a draw-bearing chemistry frame: n_samples x n_draws, n_analytes metals.
make_frame <- function(n_samples, n_draws, n_analytes) {
  metals <- .METALS[seq_len(n_analytes)]
  dates  <- as.Date("2020-01-01") + seq_len(n_samples) - 1L
  base   <- stats::setNames(stats::runif(n_analytes, 0.5, 50), metals)
  md <- purrr::map_dfr(seq_len(n_samples), function(i) {
    purrr::map_dfr(metals, function(a) tibble::tibble(
      sample_id = paste0("s", i), site_id = "f1", datetime = dates[i],
      analyte = a, value = base[[a]] * exp(stats::rnorm(n_draws, 0, 0.25)),
      detected = TRUE, draw_id = seq_len(n_draws)))
  })
  co <- purrr::map_dfr(seq_len(n_samples), function(i) tibble::tibble(
    sample_id = paste0("s", i), site_id = "f1", datetime = dates[i],
    analyte = names(CO), value = unname(CO), detected = TRUE,
    draw_id = NA_integer_))
  dplyr::bind_rows(md, co)
}

grid <- rbind(
  expand.grid(n_samples = c(30L, 100L), n_draws = c(1L, 8L), n_analytes = 9L),
  expand.grid(n_samples = 30L, n_draws = 8L, n_analytes = c(3L, 15L))
)

res <- purrr::pmap_dfr(grid, function(n_samples, n_draws, n_analytes) {
  df <- make_frame(n_samples, n_draws, n_analytes)
  b  <- bench::mark(
    add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws"),
    iterations = 2L, check = FALSE, filter_gc = FALSE
  )
  tibble::tibble(
    n_samples = n_samples, n_draws = n_draws, n_analytes = n_analytes,
    n_blocks = n_samples * n_draws,
    median_s = as.numeric(b$median), mem_mb = as.numeric(b$mem_alloc) / 1024^2
  )
})
print(res, n = Inf)

dir.create("bench", showWarnings = FALSE)
out <- file.path("bench", sprintf("results-%s.rds", label))
saveRDS(res, out)
cat("WROTE", out, "\n")
