#!/usr/bin/env Rscript
# Issue #31: LOO anchor-coverage of the daily residual smoother on real B.S01
# data, swept over ou_scale (= gamma multiplier) and block size.
#
# The harness (.loo_anchor_coverage / .loo_coverage_series, R/kalman_bridge.R)
# holds out each anchor (block=1) or each run of `block` consecutive anchors and
# checks whether the held-out residual falls inside the state-posterior interval.
# Coverage is in LATENT residual space (what the bridge models), so this isolates
# the bridge-variance calibration from the convex-SSD/floor transform (#15):
#   * coverage ~ nominal at scale=1  -> bridge variance OK; AmsPAF inflation is
#     the transform (#15), not the variance.
#   * coverage > nominal (over-covered) -> bridge too wide -> lower ou_scale.
# block >= 2 probes mid-gap width (plain LOO mostly tests the pinch zone).
#
#   Rscript dev/loo_coverage_bs01.R

suppressMessages({ library(dplyr); devtools::load_all(".", quiet = TRUE) })
options(hydroSense.guideline_dir = "guideline data")

cc <- qs2::qs_read("test data/bs01_v3_cache.qs2")
da <- cc$daily_args

cat("Fitting target model for B.S01 ...\n")
tm <- suppressMessages(fit_target_model(da$df, da$reference_model))

NOMINAL <- 0.90
scales  <- c(0.25, 0.5, 1, 2)
blocks  <- c(1L, 3L)

cat(sprintf("\nPooled LOO coverage (nominal %.0f%%), by ou_scale x block:\n",
            100 * NOMINAL))
grid <- expand.grid(scale = scales, block = blocks)
pooled <- purrr::pmap_dfr(grid, function(scale, block) {
  cv <- .loo_anchor_coverage(tm, interval = NOMINAL, block = block, scale = scale)
  p  <- dplyr::filter(cv, .data$analyte == "(pooled)")
  tibble::tibble(ou_scale = scale, block = block,
                 pooled_coverage = round(p$coverage, 3),
                 mean_width = round(p$mean_width, 3), n = p$n)
})
print(as.data.frame(pooled))

cat("\nPer-analyte coverage at ou_scale = 1, block = 1 ",
    "(coverage >> nominal = band too wide for that analyte):\n", sep = "")
pa <- .loo_anchor_coverage(tm, interval = NOMINAL, block = 1L, scale = 1) |>
  dplyr::filter(!is.na(.data$coverage)) |>
  dplyr::arrange(dplyr::desc(.data$coverage))
print(as.data.frame(pa), row.names = FALSE)

qs2::qs_save(list(pooled = pooled, per_analyte = pa),
             "dev/loo_coverage_bs01.qs2")
cat("\nWROTE dev/loo_coverage_bs01.qs2\n")
