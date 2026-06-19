#!/usr/bin/env Rscript
# Issue #32: visual before/after comparison of independent vs coupled 90% credible
# bands on B.S01 (ARA path).
#
# Two stacked panels (A: independent draws, B: coupled draws), grey-dashed
# deterministic centre, 90% credible ribbon, and median line.  A caption reports
# the mean CI width ratio (coupled/independent).
#
#   Rscript dev/bs01_coupling_compare.R  (first run: slow — two draw sets)

suppressMessages({ library(dplyr); library(ggplot2); devtools::load_all(".", quiet = TRUE) })
options(hydroSense.guideline_dir = "guideline data")

CACHE        <- "test data/bs01_v3_cache.qs2"
DRAWS_INDEP  <- "dev/bs01_kalman_draws_ara_indep.qs2"
DRAWS_COUP   <- "dev/bs01_kalman_draws_ara_coupled.qs2"
CENTRE       <- "dev/bs01_kalman_centre.qs2"
OUT_PNG      <- "dev/bs01_coupling_compare.png"

N_DRAWS <- 20L
SEED    <- 42L

cc <- qs2::qs_read(CACHE)
da <- cc$daily_args

## ── Generate draws (cached) ───────────────────────────────────────────────────
run_draws <- function(couple) {
  args                  <- da
  args$ndraws           <- N_DRAWS
  args$seed             <- SEED
  args$return           <- "draws"
  args$couple_residuals <- couple
  args$grab_cv          <- 0.15
  suppressMessages(do.call(amspaf_daily, args))
}

if (!file.exists(DRAWS_INDEP)) {
  cat("Independent draws ...\n")
  qs2::qs_save(run_draws(FALSE), DRAWS_INDEP)
  gc()
}
if (!file.exists(DRAWS_COUP)) {
  cat("Coupled draws ...\n")
  qs2::qs_save(run_draws(TRUE), DRAWS_COUP)
  gc()
}

## ── Deterministic centre ──────────────────────────────────────────────────────
if (file.exists(CENTRE)) {
  ctr_cache <- qs2::qs_read(CENTRE)
  centre    <- ctr_cache$ara |> dplyr::select("date", "amspaf")
} else {
  cat("Computing deterministic centre (point mode) ...\n")
  centre <- suppressMessages(do.call(amspaf_daily, da)) |>
    dplyr::select("date", "amspaf")
}

## ── Summarise draws at 90% interval ──────────────────────────────────────────
summarise_draws <- function(path) {
  qs2::qs_read(path) |>
    dplyr::group_by(date) |>
    dplyr::summarise(
      amspaf_median = median(amspaf, na.rm = TRUE),
      amspaf_lower  = quantile(amspaf, 0.05, na.rm = TRUE),
      amspaf_upper  = quantile(amspaf, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

summ_indep <- summarise_draws(DRAWS_INDEP)
summ_coup  <- summarise_draws(DRAWS_COUP)

## ── Width ratio (for caption) ─────────────────────────────────────────────────
width_indep <- mean(summ_indep$amspaf_upper - summ_indep$amspaf_lower, na.rm = TRUE)
width_coup  <- mean(summ_coup$amspaf_upper  - summ_coup$amspaf_lower,  na.rm = TRUE)
ratio       <- width_coup / width_indep
shift       <- mean(summ_coup$amspaf_median - summ_indep$amspaf_median, na.rm = TRUE)

## ── Join centre ───────────────────────────────────────────────────────────────
df_indep <- dplyr::left_join(summ_indep, dplyr::rename(centre, amspaf_centre = "amspaf"), by = "date")
df_coup  <- dplyr::left_join(summ_coup,  dplyr::rename(centre, amspaf_centre = "amspaf"), by = "date")

ymax <- max(c(df_indep$amspaf_upper, df_coup$amspaf_upper), na.rm = TRUE)
ylim <- ggplot2::ylim(0, ymax)

START <- min(centre$date, na.rm = TRUE)
END   <- max(centre$date, na.rm = TRUE)
xlim  <- ggplot2::scale_x_date(limits = c(START, END))
thm   <- ggplot2::theme_minimal(base_size = 10)
ylab  <- "% species affected (ARA)"
cap   <- sprintf("Mean CI width: independent = %.2f, coupled = %.2f — ratio %.3f (%.1f%%)",
                 width_indep, width_coup, ratio, 100 * ratio)

pA <- ggplot(df_indep, aes(date)) +
  geom_ribbon(aes(ymin = amspaf_lower, ymax = amspaf_upper), fill = "steelblue", alpha = 0.22) +
  geom_line(aes(y = amspaf_centre), colour = "grey35", linetype = "21", linewidth = 0.5) +
  geom_line(aes(y = amspaf_median), colour = "steelblue4", linewidth = 0.7) +
  xlim + ylim +
  labs(title = "A. Independent draws (90% CI)", x = NULL, y = ylab) +
  thm

pB <- ggplot(df_coup, aes(date)) +
  geom_ribbon(aes(ymin = amspaf_lower, ymax = amspaf_upper), fill = "steelblue3", alpha = 0.22) +
  geom_line(aes(y = amspaf_centre), colour = "grey35", linetype = "21", linewidth = 0.5) +
  geom_line(aes(y = amspaf_median), colour = "steelblue", linewidth = 0.7) +
  xlim + ylim +
  labs(title = "B. Coupled draws (90% CI)", x = NULL, y = ylab, caption = cap) +
  thm

g <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(pA, pB, ncol = 1)
} else {
  gridExtra::grid.arrange(pA, pB, ncol = 1)
}
ggsave(OUT_PNG, g, width = 12, height = 7, dpi = 120)
cat("WROTE", OUT_PNG, "\n")

## ── Console report ────────────────────────────────────────────────────────────
cat(sprintf("\nMean CI width (independent): %5.1f\n", width_indep))
cat(sprintf("Mean CI width (coupled):     %5.1f\n", width_coup))
cat(sprintf("Width ratio (coupled/indep): %5.3f\n", ratio))
cat(sprintf("Median shift (coupled centre - indep centre): %5.1f (should be ~0 by construction)\n", shift))
