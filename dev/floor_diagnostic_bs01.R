#!/usr/bin/env Rscript
# Issue #15: is the ARA additive-model floor (max(C_norm - ref_norm, 0)) a
# material distortion that the started-ratio would fix, or cosmetic?
#
# Uses the per-analyte breakdown carried on the cached daily draws
# (attr(draws, "analyte_pafs"): draw_id, analyte, C_adj, PAF, date). C_adj == 0
# is exactly a floor-fire (impact clipped to zero). We ask:
#   1. How often does the floor fire (overall, per analyte)?
#   2. Is it concentrated at baseline (low deterministic centre)?
#   3. Does the floor drive the MEAN inflation, or is that the upper tail of the
#      convex SSD? (The started-ratio reshapes the near-crossover/lower tail; it
#      cannot change the upper-tail convex amplification.)
#
#   Rscript dev/floor_diagnostic_bs01.R

suppressMessages(library(dplyr))
dr  <- qs2::qs_read("dev/bs01_kalman_draws_ara.qs2")
ap  <- attr(dr, "analyte_pafs")
ctr <- qs2::qs_read("dev/bs01_kalman_centre.qs2")$ara |>
  dplyr::rename(det = "mspaf") |>
  dplyr::mutate(regime = dplyr::case_when(det < 10 ~ "baseline",
                                          det >= 20 ~ "pulse",
                                          TRUE ~ "mid"))
ap <- ap |>
  dplyr::mutate(floored = .data$C_adj <= 0) |>
  dplyr::left_join(ctr[, c("date", "det", "regime")], by = "date")

cat(sprintf("Total per-(day,analyte,draw) cells: %d\n", nrow(ap)))
cat(sprintf("Overall floor-firing rate (C_adj==0): %.1f%%\n",
            100 * mean(ap$floored)))

cat("\n[1] Floor-firing rate + mean PAF per analyte:\n")
ap |> group_by(analyte) |>
  summarise(n = n(), floor_rate = round(100 * mean(floored), 1),
            mean_PAF = round(mean(PAF, na.rm = TRUE), 2), .groups = "drop") |>
  arrange(desc(n)) |> as.data.frame() |> print(row.names = FALSE)

cat("\n[2] Floor-firing rate by day regime (deterministic ARA centre):\n")
ap |> filter(!is.na(regime)) |> group_by(regime) |>
  summarise(n = n(), floor_rate = round(100 * mean(floored), 1),
            mean_PAF = round(mean(PAF, na.rm = TRUE), 2), .groups = "drop") |>
  as.data.frame() |> print(row.names = FALSE)

cat("\n[3] Is the mean inflation floor-driven or convex-SSD (upper tail)?\n")
thr <- stats::quantile(ap$PAF, 0.9, na.rm = TRUE)
ap |> mutate(grp = ifelse(PAF >= thr, "top-decile PAF (mean drivers)",
                          "lower 90% PAF")) |>
  group_by(grp) |>
  summarise(n = n(), floor_rate = round(100 * mean(floored), 1),
            mean_C_adj = round(mean(C_adj, na.rm = TRUE), 2),
            mean_PAF = round(mean(PAF, na.rm = TRUE), 3), .groups = "drop") |>
  as.data.frame() |> print(row.names = FALSE)

cat("\n[4] Near-crossover occupancy (where started-ratio differs most):\n")
cat("    Of NON-floored cells, fraction with tiny positive C_adj (< 1 ug/L norm):\n")
nz <- ap |> filter(!floored)
cat(sprintf("      %.1f%% of non-floored cells have C_adj < 1\n",
            100 * mean(nz$C_adj < 1, na.rm = TRUE)))
cat(sprintf("      median non-floored C_adj = %.2f\n",
            stats::median(nz$C_adj, na.rm = TRUE)))
