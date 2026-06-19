#!/usr/bin/env Rscript
# Freeze the CURRENT (pre-rework) deterministic daily msPAF centre line for
# B.S01, so the Kalman-smoother rework can be sanity-checked against it.
#
# Why this exists: the rework REPLACES the `.interp_residual` centre line with a
# state-space posterior. The new centre should differ from this baseline (that's
# the point) but NOT wildly. This script saves the old centre line + a plot-free
# data object the rework can diff/plot against.
#
# It uses POINT MODE only (interpolation = "model", no ndraws), which does NOT
# touch the OU/draws uncertainty machinery (that path has bugs and is being
# replaced). Heavy setup (DB pull, SILO, reference model) is read from the v3
# cache, so no DB/network is needed.
#
# Run from the package root:  Rscript dev/baseline_bs01_centreline.R

suppressMessages({
  library(dplyr)
  devtools::load_all(".", quiet = TRUE)
})

CACHE_V3 <- "test data/bs01_v3_cache.qs2"
OUT_QS2  <- "dev/baseline_bs01_centreline.qs2"
GUIDE    <- "guideline data"   # SSD/guideline tables (mspaf_daily needs these)

stopifnot(file.exists(CACHE_V3))
options(hydroSense.guideline_dir = GUIDE)

cc <- qs2::qs_read(CACHE_V3)
da <- cc$daily_args            # df, reference_model, temperature, start, end,
                               # interpolation = "model", require_temperature, min_analytes

cat("Loaded v3 cache. daily_args fields:", paste(names(da), collapse = ", "), "\n")

## ── Point-mode centre lines (deterministic; no ndraws) ───────────────────────
# ARA mode: reference = the fitted reference_model (impact = C_norm - ref_norm,
#   which cancels ref in C_excess -> msPAF reflects the impact only).
# non-ARA (total concentration) mode: reference = NULL (no background subtracted)
#   but reference_model is still supplied so interpolation = "model" can fit the
#   season-blind impact model used to interpolate between grabs.

run_point <- function(reference) {
  do.call(mspaf_daily, c(da, list(reference = reference)))
}

cat("Running ARA point-mode centre line...\n")
ara_now <- run_point(reference = da$reference_model)

cat("Running non-ARA (total) point-mode centre line...\n")
tot_now <- run_point(reference = NULL)

pick <- function(x) x[, intersect(c("date", "site_id", "mspaf"), names(x)), drop = FALSE]
ara_now <- pick(ara_now)
tot_now <- pick(tot_now)

## ── Cross-check against the cached ARA series (made at v3 time) ───────────────
# If the point-mode centre line has drifted since the cache was built, report it
# so we know the "current" baseline is the freshly regenerated one.
if (!is.null(cc$daily_ara)) {
  cmp <- dplyr::inner_join(
    dplyr::rename(pick(cc$daily_ara), mspaf_cached = "mspaf"),
    dplyr::rename(ara_now,           mspaf_now    = "mspaf"),
    by = intersect(c("date", "site_id"), names(ara_now))
  )
  if (nrow(cmp) > 0) {
    md <- max(abs(cmp$mspaf_cached - cmp$mspaf_now), na.rm = TRUE)
    rr <- suppressWarnings(stats::cor(cmp$mspaf_cached, cmp$mspaf_now,
                                      use = "complete.obs"))
    cat(sprintf("ARA cached-vs-regenerated: n=%d  max|diff|=%.3g  cor=%.5f\n",
                nrow(cmp), md, rr))
  }
}

baseline <- list(
  ara         = ara_now,         # date, site_id, mspaf  (ARA / added-risk)
  total       = tot_now,         # date, site_id, mspaf  (non-ARA / total conc)
  source      = "dev/baseline_bs01_centreline.R",
  generated   = Sys.time(),
  point_mode  = TRUE,
  note        = paste("Pre-rework .interp_residual centre line for B.S01.",
                      "Compare new state-space centre against this:",
                      "expect a difference (mostly mid-gap, ~0 at anchors),",
                      "but not a wild one.")
)
qs2::qs_save(baseline, OUT_QS2)
cat("WROTE", OUT_QS2, "\n")
cat(sprintf("  ARA   rows: %d  (mspaf range %.4g .. %.4g)\n",
            nrow(ara_now), min(ara_now$mspaf, na.rm = TRUE),
            max(ara_now$mspaf, na.rm = TRUE)))
cat(sprintf("  total rows: %d  (mspaf range %.4g .. %.4g)\n",
            nrow(tot_now), min(tot_now$mspaf, na.rm = TRUE),
            max(tot_now$mspaf, na.rm = TRUE)))
