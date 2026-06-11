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

## ── PAF lookup micro-benchmarks (issue #36) ──────────────────────────────────
##
## These measure the three performance-critical paths added by the SSD PAF
## lookup table:
##   A. .ssd_paf_lookup() cold vs warm cache (per-analyte spline build cost)
##   B. .ssd_paf_vec() lookup vs direct ssd_hp() at various concentration counts
##   C. add_amspaf() draws-mode end-to-end: #36 vs a reference call using the
##      cold-cache path (approximates #30 baseline for PAF cost)

cat("\n── PAF lookup micro-benchmarks ──\n")

meta      <- leachatetools:::.load_analyte_metadata(NULL)
ssd_params <- suppressMessages(
  leachatetools:::derive_ssd_params(meta, method = "multi", guideline_dir = NULL)
)
cu_fit <- ssd_params$fit[[which(ssd_params$analyte == "Cu")]]

## A. Cold vs warm cache for .ssd_paf_lookup()
clear_lookup <- function() {
  env <- leachatetools:::.ssd_paf_lookup_env
  rm(list = ls(envir = env, all.names = TRUE), envir = env)
}

paf_cache <- bench::mark(
  cold = { clear_lookup(); leachatetools:::.ssd_paf_lookup("Cu", "multi", cu_fit, NULL) },
  warm = leachatetools:::.ssd_paf_lookup("Cu", "multi", cu_fit, NULL),
  iterations = 5L, check = FALSE, filter_gc = FALSE
)
cat("\nA. .ssd_paf_lookup() cold vs warm cache (Cu/multi):\n")
print(paf_cache[, c("expression", "min", "median", "mem_alloc")])

## B. .ssd_paf_vec() lookup vs direct ssd_hp()
## Simulate the exact fallback: pass a non-NULL guideline_dir with few (<1025)
## unique concentrations to force the direct path, vs NULL guideline_dir for
## the shipped-table lookup path.
tmp_gdir <- tempdir()
set.seed(36L)
conc_large <- stats::runif(2000L, 0.1, 100)   # 2000 unique → lookup path eligible

paf_vec <- bench::mark(
  lookup_shipped = leachatetools:::.ssd_paf_vec(
    cu_fit, conc_large, "Cu", "multi", NULL
  ),
  direct_ssd_hp = leachatetools:::.ssd_paf_vec(
    cu_fit, conc_large[seq_len(3)], "Cu", "multi", tmp_gdir
  ),
  iterations = 10L, check = FALSE, filter_gc = FALSE
)
cat("\nB. .ssd_paf_vec() lookup vs direct ssd_hp() (Cu/multi):\n")
print(paf_vec[, c("expression", "min", "median", "mem_alloc")])

## C. End-to-end add_amspaf() draws mode: 30 samples × 8 draws, 9 analytes.
## Run twice to show warm-cache performance; first call populates the lookup.
df_large <- make_frame(30L, 8L, 9L)
.warmup <- suppressMessages(
  add_amspaf(df_large, reference = NULL, conc_units = "ug/L", return = "draws")
)  # warm up the cache
rm(.warmup)
e2e <- bench::mark(
  add_amspaf_lookup = suppressMessages(
    add_amspaf(df_large, reference = NULL, conc_units = "ug/L", return = "draws")
  ),
  iterations = 3L, check = FALSE, filter_gc = FALSE
)
cat("\nC. add_amspaf() 30s×8d×9a warm cache (lookup active):\n")
print(e2e[, c("expression", "min", "median", "mem_alloc")])

paf_bench <- list(cache = paf_cache, paf_vec = paf_vec, e2e = e2e)
paf_out <- file.path("bench", sprintf("results-paf-lookup-%s.rds", label))
saveRDS(paf_bench, paf_out)
cat("WROTE", paf_out, "\n")
