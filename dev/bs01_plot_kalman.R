#!/usr/bin/env Rscript
# B.S01 comparison plot for the Kalman daily-uncertainty rework (issue #16).
#
# Six stacked panels (shared x-axis):
#   A. OLD centre line, no ARA (total)        — pre-rework .interp_residual
#   B. OLD centre line, ARA
#   C. NEW centre + credible interval, no ARA  — state-space smoother (old centre
#      overlaid dashed for comparison)
#   D. NEW centre + credible interval, ARA     (old centre overlaid dashed)
#   E. Leachate mixing fraction (existing, grab samples)
#   F. Sampling events — a dot per grab; fill alpha encodes how complete the
#      sampling was that day (full metals + WQ = solid black; one/two analytes =
#      light grey). A guide for reading the panels above, not a calibrated scale.
#
# Heavy setup (DB pull, SILO, reference model, deterministic daily series) is read
# from the v3 cache; the OLD centre lines from the frozen baseline. The NEW draws
# are cached separately (slow). Run from the package root:
#   Rscript dev/bs01_plot_kalman.R

suppressMessages({
  library(dplyr); library(ggplot2)
  devtools::load_all(".", quiet = TRUE)
})

CACHE_V3  <- "test data/bs01_v3_cache.rds"
BASELINE  <- "dev/baseline_bs01_centreline.rds"
CACHE_ARA <- "dev/bs01_kalman_cache_ara.rds"   # separate caches = independent
CACHE_TOT <- "dev/bs01_kalman_cache_tot.rds"   # checkpoints, lighter memory
GUIDE     <- "guideline data"
N_DRAWS   <- 20L   # add_amspaf over N draws x 1800 days dominates runtime
SEED      <- 42L
INTERVAL  <- 0.9
TOX_RSD   <- 0.15
options(leachatetools.guideline_dir = GUIDE)

stopifnot(file.exists(CACHE_V3), file.exists(BASELINE))
cc   <- readRDS(CACHE_V3)
base <- readRDS(BASELINE)
da   <- cc$daily_args
START <- da$start; END <- da$end
pick <- function(x) x[, intersect(c("date","amspaf","amspaf_lower","amspaf_upper"),
                                  names(x)), drop = FALSE]

## ── NEW: state-space daily AmsPAF with credible interval (ARA + non-ARA) ──────
## Each call is cached separately and memory released between them.
run_draws <- function(reference) suppressMessages(do.call(amspaf_daily, c(da, list(
  reference = reference, ndraws = N_DRAWS, seed = SEED, return = "summary",
  interval = INTERVAL, grab_cv = TOX_RSD))))
if (!file.exists(CACHE_ARA)) {
  cat("ARA draws ...\n"); saveRDS(run_draws(da$reference_model), CACHE_ARA); gc()
}
if (!file.exists(CACHE_TOT)) {
  cat("total (non-ARA) draws ...\n"); saveRDS(run_draws(NULL), CACHE_TOT); gc()
}
new_ara <- pick(readRDS(CACHE_ARA)); new_tot <- pick(readRDS(CACHE_TOT))

## ── Sampling-events completeness (per grab date) ─────────────────────────────
tc <- cc$target_chem
ev <- tc |>
  dplyr::filter(is.na(.data$detected) | .data$detected) |>
  dplyr::group_by(date = as.Date(.data$datetime)) |>
  dplyr::summarise(n_analytes = dplyr::n_distinct(.data$analyte), .groups = "drop")
# Alpha ~ completeness; a rough guide only (not a calibrated scale).
full <- stats::quantile(ev$n_analytes, 0.9, names = FALSE)
ev$alpha <- pmax(0.12, pmin(1, ev$n_analytes / full))
ev$y <- 1

## ── Panels ───────────────────────────────────────────────────────────────────
xlim <- ggplot2::scale_x_date(limits = c(START, END))
thm  <- ggplot2::theme_minimal(base_size = 10)
ylab <- "% species affected"

old_line <- function(d, ttl, sub) ggplot(d, aes(date, amspaf)) +
  geom_line(colour = "grey25", linewidth = 0.6) + xlim +
  labs(title = ttl, subtitle = sub, x = NULL, y = ylab) + thm

new_band <- function(dnew, dold, ttl, sub) ggplot() +
  geom_ribbon(data = dnew, aes(date, ymin = amspaf_lower, ymax = amspaf_upper),
              fill = "steelblue", alpha = 0.25) +
  geom_line(data = dold, aes(date, amspaf), colour = "grey55",
            linetype = "21", linewidth = 0.5) +
  geom_line(data = dnew, aes(date, amspaf), colour = "steelblue4",
            linewidth = 0.7) +
  xlim + labs(title = ttl, subtitle = sub, x = NULL, y = ylab) + thm

pA <- old_line(base$total, "A. Old centre line — no ARA (total mixture)",
               "Pre-rework deterministic interpolation")
pB <- old_line(base$ara,   "B. Old centre line — Added Risk (ARA)",
               "Pre-rework deterministic interpolation")
pC <- new_band(new_tot, base$total,
               "C. New centre + 90% CI — no ARA",
               "State-space smoother (blue) vs old centre (grey dashed)")
pD <- new_band(new_ara, base$ara,
               "D. New centre + 90% CI — Added Risk (ARA)",
               "State-space smoother (blue) vs old centre (grey dashed)")
pE <- ggplot(cc$lmf_ts, aes(datetime, lmf)) +
  geom_line(colour = "grey60") + geom_point(size = 1, colour = "grey20") +
  xlim + labs(title = "E. Leachate mixing fraction (grab samples)",
              x = NULL, y = "LMF (%)") + thm
pF <- ggplot(ev, aes(date, y)) +
  geom_hline(yintercept = 1, colour = "grey80", linewidth = 0.4) +
  geom_point(aes(alpha = alpha), shape = 21, fill = "black",
             colour = "grey30", size = 3, stroke = 0.4) +
  scale_alpha_identity() +
  xlim + scale_y_continuous(limits = c(0.8, 1.2), breaks = NULL) +
  labs(title = "F. Sampling events (fill alpha = sampling completeness)",
       subtitle = "Solid = full metals + WQ suite; faint = only one or two analytes",
       x = NULL, y = NULL) + thm

g <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(pA, pB, pC, pD, pE, pF, ncol = 1)
} else {
  gridExtra::grid.arrange(pA, pB, pC, pD, pE, pF, ncol = 1)
}
ggsave("dev/bs01_kalman_compare.png", g, width = 10, height = 18, dpi = 120)
cat("WROTE dev/bs01_kalman_compare.png\n")
