#!/usr/bin/env Rscript
# Before/after comparison for the issue #15 asinh transform, on B.S01.
#
# "Before" = the frozen additive-model results in dev/before_asinh/ (captured by
# the freeze step on the pre-#15 code). "After" = the current working caches
# dev/bs01_kalman_*.qs2 regenerated on the #15 branch (run dev/bs01_plot_kalman.R
# and dev/loo_coverage_bs01.R first).
#
# Reports two things the user asked for:
#   1. INVARIANTS that must NOT change (engine + anchor-exact centre).
#   2. The new uncertainty band vs the old band and the deterministic centre --
#      is the change actually an improvement (baseline tightens, events hold)?
#
#   Rscript dev/compare_asinh.R

suppressMessages({ library(dplyr); devtools::load_all(".", quiet = TRUE) })
options(hydroSense.guideline_dir = "guideline data")
B <- "dev/before_asinh"
GAP <- as.Date(c("2023-12-18", "2024-02-17"))   # the Jan-2024 baseline gap
EVT <- as.Date(c("2024-08-15", "2024-10-15"))   # the 2024 impact event

ok  <- function(x) cat(sprintf("  [OK]   %s\n", x))
bad <- function(x) cat(sprintf("  [DIFF] %s\n", x))

cat("== 1. INVARIANTS (must be unchanged) ==\n")

## 1a. Engine: add_mspaf on the grabs (R/mspaf.R is untouched by #15).
inv <- qs2::qs_read(file.path(B, "engine_invariant.qs2"))
cc  <- qs2::qs_read("test data/bs01_v3_cache.qs2")
eng_now <- suppressMessages(add_mspaf(cc$target_chem, reference = NULL)) |>
  dplyr::filter(.data$analyte == "msPAF") |>
  dplyr::select(sample_id, value) |> dplyr::arrange(sample_id)
d_eng <- max(abs(eng_now$value - inv$engine_mspaf_total$value), na.rm = TRUE)
if (isTRUE(d_eng < 1e-9)) ok(sprintf("per-sample add_mspaf identical (max|d|=%.2g)", d_eng)) else
  bad(sprintf("per-sample add_mspaf CHANGED (max|d|=%.3g)", d_eng))

## 1b. Anchor-exactness: deterministic centre at grab dates unchanged.
c_old <- qs2::qs_read(file.path(B, "bs01_kalman_centre.qs2"))
c_new <- qs2::qs_read("dev/bs01_kalman_centre.qs2")
for (which in c("ara", "total")) {
  j <- dplyr::inner_join(
    dplyr::rename(c_old[[which]], old = "mspaf"),
    dplyr::rename(c_new[[which]], new = "mspaf"), by = "date") |>
    dplyr::filter(date %in% inv$grab_dates)
  d <- max(abs(j$old - j$new), na.rm = TRUE)
  if (isTRUE(d < 1e-2)) ok(sprintf("deterministic centre at grabs (%s) unchanged (max|d|=%.3g)", which, d)) else
    bad(sprintf("deterministic centre at grabs (%s) CHANGED (max|d|=%.3g)", which, d))
}

cat("\n== 2. IMPROVEMENT (band vs old band vs deterministic) ==\n")

band_summary <- function(path) {
  qs2::qs_read(path) |> dplyr::group_by(date) |>
    dplyr::summarise(median = median(mspaf, na.rm = TRUE),
                     mean = mean(mspaf, na.rm = TRUE),
                     iqr = quantile(mspaf, .75, names = FALSE, na.rm = TRUE) -
                           quantile(mspaf, .25, names = FALSE, na.rm = TRUE),
                     .groups = "drop")
}
win <- function(df, w) dplyr::filter(df, date >= w[1], date <= w[2])

for (which in c("ara", "tot")) {
  cat(sprintf("\n-- %s --\n", toupper(which)))
  old <- band_summary(file.path(B, sprintf("bs01_kalman_draws_%s.qs2", which)))
  new <- band_summary(sprintf("dev/bs01_kalman_draws_%s.qs2", which))
  det <- (if (which == "tot") c_new$total else c_new$ara) |> dplyr::rename(det = "mspaf")
  for (lab in c("baseline gap", "2024 event")) {
    w <- if (lab == "baseline gap") GAP else EVT
    o <- win(old, w); n <- win(new, w); dd <- win(det, w)
    cat(sprintf("  %-13s: median %.1f->%.1f  mean %.1f->%.1f  IQR %.1f->%.1f  (det %.1f)\n",
                lab, mean(o$median), mean(n$median), mean(o$mean), mean(n$mean),
                mean(o$iqr), mean(n$iqr), mean(dd$det, na.rm = TRUE)))
  }
}

## Per-analyte headline: NH3-N PAF upper tail in the baseline gap.
cat("\n-- NH3-N baseline-gap PAF upper tail (the #39 smoking gun) --\n")
nh <- function(path) {
  ap <- attr(qs2::qs_read(path), "analyte_pafs")
  ap <- ap[ap$analyte == "NH3-N" & ap$date > GAP[1] & ap$date < GAP[2], ]
  c(median = median(ap$PAF), mean = mean(ap$PAF),
    q99 = quantile(ap$PAF, .99, names = FALSE))
}
o <- nh(file.path(B, "bs01_kalman_draws_ara.qs2")); n <- nh("dev/bs01_kalman_draws_ara.qs2")
cat(sprintf("  median %.3f->%.3f  mean %.3f->%.3f  q99 %.3f->%.3f  (want q99 down)\n",
            o["median"], n["median"], o["mean"], n["mean"], o["q99"], n["q99"]))

## LOO coverage (run dev/loo_coverage_bs01.R after A to refresh the working file).
cat("\n-- LOO pooled coverage (nominal 0.90) --\n")
lo <- qs2::qs_read(file.path(B, "loo_coverage_bs01.qs2"))$pooled
ln <- tryCatch(qs2::qs_read("dev/loo_coverage_bs01.qs2")$pooled, error = function(e) NULL)
cat("  before:\n"); print(as.data.frame(lo[lo$block == 1, ]), row.names = FALSE)
if (!is.null(ln)) { cat("  after:\n"); print(as.data.frame(ln[ln$block == 1, ]), row.names = FALSE) } else
  cat("  after: (re-run dev/loo_coverage_bs01.R)\n")
