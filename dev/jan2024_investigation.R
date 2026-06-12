#!/usr/bin/env Rscript
# Investigation: why does the daily ARA median "spike" around Jan 2024 while the
# deterministic centre stays flat? (B.S01)
#
# Findings (see issue threads #15, #31, #34, #39):
#  * Jan 2024 sits in a ~63-day sampling gap (last grab 2023-12-17, next
#    2024-02-18). The deterministic centre interpolates flat ~8.5%; the draw
#    median climbs to ~26%.
#  * The combined median ~= combined mean ~= IA(per-analyte means) ~27%: the
#    combination of ~7 skewed INDEPENDENT analytes symmetrises, so the median
#    does NOT escape the inflation (it is not a mean-vs-median issue at the
#    combined level).
#  * Root cause: HOMOSCEDASTIC ADDITIVE process variance. NH3-N draw spread is
#    ~constant in absolute terms across regimes (IQR ~1374 baseline vs ~1697
#    event) while the level changes, so baseline draws are grossly over-wide
#    (NH3-N impact draws reach PAF q99 ~0.95 during a 0.02-0.21 mg/L period).
#    Through the convex SSD this inflates E[PAF] and the combined median.
#  * Verdict: the flat deterministic is the defensible central estimate; the
#    spiked median is an over-dispersion artefact. The fix is a proportional
#    (multiplicative) variance: the started-ratio transform (#15 con b) or
#    state-dependent variance (#34) -- NOT lowering ou_scale (#31 pooled
#    coverage is ~nominal because it averages over regimes).
#
#   Rscript dev/jan2024_investigation.R

suppressMessages(library(dplyr))
cc  <- qs2::qs_read("test data/bs01_v3_cache.qs2")
dr  <- qs2::qs_read("dev/bs01_kalman_draws_ara.qs2")
ap  <- attr(dr, "analyte_pafs")
ctr <- qs2::qs_read("dev/bs01_kalman_centre.qs2")$ara |> dplyr::rename(det = "amspaf")
GAP <- as.Date(c("2023-12-18", "2024-02-17"))
apg <- ap |> filter(date > GAP[1], date < GAP[2])
drg <- dr |> filter(date > GAP[1], date < GAP[2])
ia  <- function(x) 100 * (1 - prod(1 - x))

cat("[1] combined vs per-analyte composition in the gap:\n")
pa <- apg |> group_by(analyte) |>
  summarise(mean_i = mean(PAF), median_i = median(PAF), .groups = "drop") |>
  filter(mean_i > 0)
cat(sprintf("    combined: mean=%.1f%% median=%.1f%% | IA(means)=%.1f%% IA(medians)=%.1f%%\n",
            mean(drg$amspaf), median(drg$amspaf), ia(pa$mean_i), ia(pa$median_i)))

cat("[2] NH3-N spread baseline gap vs 2024 event (homoscedasticity test):\n")
w <- function(a, b) ap |> filter(analyte == "NH3-N", date >= as.Date(a), date <= as.Date(b))
for (r in list(c("baseline", "2023-12-18", "2024-02-17"),
               c("event",    "2024-08-15", "2024-10-15"))) {
  d <- w(r[2], r[3])
  cat(sprintf("    %-9s C_adj median=%7.1f IQR=%8.1f | PAF median=%.3f mean=%.3f q90=%.3f\n",
              r[1], median(d$C_adj), IQR(d$C_adj), median(d$PAF), mean(d$PAF),
              quantile(d$PAF, .9, names = FALSE)))
}

cat("[3] median/deterministic ratio, sampled vs mid-gap (pervasive, not gap-only):\n")
dr |> group_by(date) |>
  summarise(median = median(amspaf, na.rm = TRUE),
            dsl = first(days_since_last_sample), .groups = "drop") |>
  left_join(ctr, by = "date") |>
  filter(date >= as.Date("2023-09-01"), date <= as.Date("2024-05-31")) |>
  mutate(grp = ifelse(dsl <= 10, "<=10d since grab", ">10d (mid-gap)")) |>
  group_by(grp) |>
  summarise(det = round(mean(det), 1), median = round(mean(median), 1),
            ratio = round(mean(median) / mean(det), 1), n = n(), .groups = "drop") |>
  as.data.frame() |> print(row.names = FALSE)
