#!/usr/bin/env Rscript
# Issue #15 crossover-scale experiment (B.S01, ARA): compare the daily msPAF
# band under c = HC5 / HC1 / HC0.1 (SSD proportions 0.05 / 0.01 / 0.001) against
# the deterministic centre. c sets where asinh(I/c) switches additive ->
# multiplicative; a lower c puts more of the impact range in the proportional
# regime -> stronger baseline variance-stabilisation. We override .analyte_c per
# run (no permanent code change) and save the summaries for plotting.
#
#   Rscript dev/crossover_experiment.R   (slow: 3 x draws regen)

suppressMessages({ library(dplyr); devtools::load_all(".", quiet = TRUE) })
options(hydroSense.guideline_dir = "guideline data")
cc <- qs2::qs_read("test data/bs01_v3_cache.qs2")
da <- cc$daily_args
N <- 20L; SEED <- 42L; INTERVAL <- 0.5; TOX_RSD <- 0.15
GAP <- as.Date(c("2023-12-18", "2024-02-17"))   # Jan-2024 baseline gap
EVT <- as.Date(c("2024-08-15", "2024-10-15"))   # 2024 event

orig_analyte_c <- hydroSense:::.analyte_c
set_prop <- function(p) {
  f <- function(fit) {
    hc <- tryCatch(ssdtools::ssd_hc(fit, proportion = p, ci = FALSE)$est,
                   error = function(e) NA_real_)
    if (length(hc) != 1L || !is.finite(hc) || hc <= 0) {   # numerical guard
      hc <- tryCatch(ssdtools::ssd_hc(fit, proportion = 0.05, ci = FALSE)$est,
                     error = function(e) NA_real_)
    }
    if (length(hc) != 1L || !is.finite(hc) || hc <= 0) {
      cli::cli_abort("no positive HC for transform scale")
    }
    hc
  }
  assignInNamespace(".analyte_c", f, "hydroSense")
}

alpha <- (1 - INTERVAL) / 2     # NB: not named `lo` -> would shadow the lo column
band <- function(dr) dr |> dplyr::group_by(date) |>
  dplyr::summarise(median = median(mspaf, na.rm = TRUE),
                   lo = quantile(mspaf, alpha, names = FALSE, na.rm = TRUE),
                   hi = quantile(mspaf, 1 - alpha, names = FALSE, na.rm = TRUE),
                   .groups = "drop")
win <- function(d, w) dplyr::filter(d, date >= w[1], date <= w[2])

props <- c(HC5 = 0.05, HC1 = 0.01, "HC0.1" = 0.001)
res <- list()
for (nm in names(props)) {
  cat("running", nm, "(p =", props[[nm]], ") ...\n")
  set_prop(props[[nm]])
  dr  <- suppressMessages(do.call(mspaf_daily, c(da, list(
    reference = da$reference_model, ndraws = N, seed = SEED,
    return = "draws", grab_cv = TOX_RSD))))
  det <- suppressMessages(do.call(mspaf_daily, c(da, list(
    reference = da$reference_model))))[, c("date", "mspaf")]
  ap  <- attr(dr, "analyte_pafs")
  nh  <- ap[ap$analyte == "NH3-N" & ap$date > GAP[1] & ap$date < GAP[2], ]
  res[[nm]] <- list(band = band(dr), det = det,
                    nh3_gap_q99 = quantile(nh$PAF, .99, names = FALSE))
  qs2::qs_save(res, "dev/crossover_experiment.qs2")   # incremental: survive a crash
}
assignInNamespace(".analyte_c", orig_analyte_c, "hydroSense")

qs2::qs_save(res, "dev/crossover_experiment.qs2")

cat("\n== ARA summary by crossover (baseline gap / 2024 event) ==\n")
for (nm in names(res)) {
  b <- res[[nm]]$band
  bg <- win(b, GAP); ev <- win(b, EVT)
  cat(sprintf("%-6s baseline: med %.1f IQR %.1f | event: med %.1f IQR %.1f | NH3-N gap q99 %.3f\n",
              nm, mean(bg$median), mean(bg$hi - bg$lo),
              mean(ev$median), mean(ev$hi - ev$lo), res[[nm]]$nh3_gap_q99))
}
cat("WROTE dev/crossover_experiment.qs2\n")
