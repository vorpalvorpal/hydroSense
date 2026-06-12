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

CACHE_V3   <- "test data/bs01_v3_cache.qs2"
BASELINE   <- "dev/baseline_bs01_centreline.qs2"
DRAWS_ARA  <- "dev/bs01_kalman_draws_ara.qs2"  # raw per-(day,draw) AmsPAF — LONG to build
DRAWS_TOT  <- "dev/bs01_kalman_draws_tot.qs2"  #   (the add_amspaf cost; see issue #30)
CENTRE     <- "dev/bs01_kalman_centre.qs2"     # deterministic centre lines — fast (point mode)
CHRONIC_ARA <- "dev/bs01_kalman_chronic_ara.qs2" # chronic TWA of the daily draws (cached)
CHRONIC_TOT <- "dev/bs01_kalman_chronic_tot.qs2"
GUIDE      <- "guideline data"
TAU        <- 90    # chronic TWA decay (days); half-life ~ tau*log(2) ~ 62 d
WINDOW     <- 365   # chronic look-back window (days)
N_DRAWS    <- 20L   # raise once #30 makes draws cheap
SEED       <- 42L
INTERVAL   <- 0.5   # 25-75% band. Change this freely: re-summarising the cached
                    # draws is INSTANT. Only deleting the DRAWS_* caches forces
                    # the long re-run.
TOX_RSD    <- 0.15
options(leachatetools.guideline_dir = GUIDE)

stopifnot(file.exists(CACHE_V3), file.exists(BASELINE))
cc   <- qs2::qs_read(CACHE_V3)
base <- qs2::qs_read(BASELINE)
da   <- cc$daily_args
START <- da$start; END <- da$end

## ── NEW: state-space daily AmsPAF (ARA + non-ARA) ────────────────────────────
## Cache the RAW per-(day, draw) AmsPAF once (return = "draws"); this is the slow
## part (add_amspaf over N x days — issue #30). Any credible interval is then an
## instant re-summarise of the cache. The deterministic centre line (which is
## interval-independent) comes from a fast point-mode run, cached separately.
draws_run <- function(reference) suppressMessages(do.call(amspaf_daily, c(da, list(
  reference = reference, ndraws = N_DRAWS, seed = SEED, return = "draws",
  grab_cv = TOX_RSD))))
if (!file.exists(DRAWS_ARA)) {
  cat("ARA draws (long) ...\n");   qs2::qs_save(draws_run(da$reference_model), DRAWS_ARA); gc()
}
if (!file.exists(DRAWS_TOT)) {
  cat("total draws (long) ...\n"); qs2::qs_save(draws_run(NULL),               DRAWS_TOT); gc()
}
if (!file.exists(CENTRE)) {
  cat("deterministic centre lines (point mode, fast) ...\n")
  ctr <- function(reference)
    suppressMessages(do.call(amspaf_daily, c(da, list(reference = reference))))
  qs2::qs_save(list(ara   = ctr(da$reference_model)[, c("date", "amspaf")],
                    total = ctr(NULL)[,              c("date", "amspaf")]), CENTRE)
}

## Summarise the cached draws at INTERVAL (instant) and attach the deterministic
## centre — equivalent to amspaf_daily(return = "summary", interval = INTERVAL),
## but the interval is now a cheap post-hoc choice. (summarise_draws() would give
## the same lower/upper from this per-(day,draw) frame.)
ctr_cache <- qs2::qs_read(CENTRE)
band <- function(draws_path, centre_df) {
  lo <- (1 - INTERVAL) / 2
  q  <- qs2::qs_read(draws_path) |>
    dplyr::group_by(.data$date) |>
    dplyr::summarise(
      amspaf_mean  = mean(.data$amspaf, na.rm = TRUE),
      amspaf_lower = stats::quantile(.data$amspaf, lo,     names = FALSE, na.rm = TRUE),
      amspaf_upper = stats::quantile(.data$amspaf, 1 - lo, names = FALSE, na.rm = TRUE),
      .groups = "drop")
  ## centre_df carries the deterministic point-mode centre (renamed amspaf_det)
  dplyr::left_join(dplyr::rename(centre_df, amspaf_det = "amspaf"), q, by = "date")
}
new_ara <- band(DRAWS_ARA, ctr_cache$ara)
new_tot <- band(DRAWS_TOT, ctr_cache$total)

## ── Chronic (time-weighted) AmsPAF ───────────────────────────────────────────
## The chronic step is a weighted ARITHMETIC mean of the daily AmsPAF (the only
## linear summary for a bounded index — geom_mean would add a second Jensen gap
## inside the chronic step). It is run PER DRAW (each draw is a coherent daily
## path), then summarised. Because the aggregation is linear:
##   * the chronic MEAN = time-average of the daily means → the decision endpoint
##     E[chronic AmsPAF], and it is robust to the #23 temporal-pairing choice
##     (correlation changes a sum's variance, never its expectation);
##   * the chronic INTERVAL still inherits the #23 index-pairing approximation.
## We plot the chronic mean as the centre line and overlay the DETERMINISTIC
## chronic (TWA of the point-mode daily centre) as a dashed reference — the gap
## between them is the residual nonlinearity bias, much smaller than daily
## because time-averaging many days pulls the distribution back toward symmetry.
chronic_summary <- function(draws_path, centre_df, cache_path) {
  if (file.exists(cache_path)) return(qs2::qs_read(cache_path))
  draws_df <- qs2::qs_read(draws_path)
  site     <- draws_df$site_id[[1]]
  to_long  <- function(d, val, id_prefix, with_draw)
    dplyr::transmute(
      d,
      sample_id = paste0(id_prefix, .data$date),
      site_id   = if ("site_id" %in% names(d)) .data$site_id else site,
      datetime  = .data$date,
      analyte   = "AmsPAF",
      value     = {{ val }},
      detected  = TRUE,
      draw_id   = if (with_draw) .data$draw_id else NA_integer_
    )
  inp <- to_long(draws_df, .data$amspaf, "d", TRUE) |>
    dplyr::filter(is.finite(.data$value))
  focal <- sort(unique(inp$datetime))

  ## draws → chronic mean (central = "mean") + credible band at INTERVAL
  cw <- time_weighted_aggregate(
    inp, focal_dates = focal, tau = TAU, tau_units = "d",
    window = WINDOW, window_units = "d",
    summary = "arith_mean", return = "summary",
    interval = INTERVAL, central = "mean"
  ) |>
    dplyr::transmute(date = .data$focal_date, amspaf = .data$value,
                     amspaf_lower = .data$value_lower,
                     amspaf_upper = .data$value_upper)

  ## deterministic daily centre → chronic (reference, no draws)
  det <- to_long(centre_df, .data$amspaf, "c", FALSE) |>
    dplyr::filter(is.finite(.data$value)) |>
    dplyr::select(-"draw_id")
  cd <- time_weighted_aggregate(
    det, focal_dates = focal, tau = TAU, tau_units = "d",
    window = WINDOW, window_units = "d", summary = "arith_mean"
  ) |>
    dplyr::transmute(date = .data$focal_date, amspaf_det = .data$value)

  out <- dplyr::left_join(cw, cd, by = "date")
  qs2::qs_save(out, cache_path)
  out
}
chr_tot <- chronic_summary(DRAWS_TOT, ctr_cache$total, CHRONIC_TOT)
chr_ara <- chronic_summary(DRAWS_ARA, ctr_cache$ara,   CHRONIC_ARA)

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

## Daily panels: blue = draw MEAN (E[AmsPAF]); grey dashed = deterministic
## centre (point mode); band = credible interval. Showing both makes the
## uncertainty-driven inflation of the mean directly visible.
new_band <- function(dnew, ttl, sub) ggplot(dnew, aes(date)) +
  geom_ribbon(aes(ymin = amspaf_lower, ymax = amspaf_upper),
              fill = "steelblue", alpha = 0.22) +
  geom_line(aes(y = amspaf_det,  colour = "deterministic centre"),
            linetype = "21", linewidth = 0.55) +
  geom_line(aes(y = amspaf_mean, colour = "mean  E[AmsPAF]"),
            linewidth = 0.7) +
  scale_colour_manual(values = c("deterministic centre" = "grey35",
                                 "mean  E[AmsPAF]" = "steelblue4"),
                      name = NULL) +
  xlim + labs(title = ttl, subtitle = sub, x = NULL, y = ylab) + thm +
  ggplot2::theme(legend.position = "top")

ci_lab <- sprintf("%g%% CI (%g-%g pct)", 100 * INTERVAL,
                  100 * (1 - INTERVAL) / 2, 100 * (1 - (1 - INTERVAL) / 2))
pC <- new_band(new_tot,
               sprintf("C. Daily AmsPAF + %s — no ARA (total mixture)", ci_lab),
               "Blue = mean E[AmsPAF]; grey dashed = deterministic centre")
pD <- new_band(new_ara,
               sprintf("D. Daily AmsPAF + %s — Added Risk (ARA)", ci_lab),
               "Blue = mean E[AmsPAF]; grey dashed = deterministic centre")

## Chronic panels: green = chronic MEAN (E[chronic AmsPAF], the decision
## endpoint); grey dashed = deterministic chronic (TWA of the point-mode centre).
chronic_band <- function(d, ttl, sub) ggplot(d, aes(date)) +
  geom_ribbon(aes(ymin = amspaf_lower, ymax = amspaf_upper),
              fill = "seagreen", alpha = 0.22) +
  geom_line(aes(y = amspaf_det), colour = "grey55", linetype = "21",
            linewidth = 0.5) +
  geom_line(aes(y = amspaf), colour = "seagreen4", linewidth = 0.7) +
  xlim + labs(title = ttl, subtitle = sub, x = NULL, y = ylab) + thm

chr_sub <- sprintf("Chronic mean (green) + %s vs deterministic chronic (grey dashed); tau=%gd, window=%gd",
                   ci_lab, TAU, WINDOW)
pE <- chronic_band(chr_tot,
               "E. Chronic (time-weighted) AmsPAF — no ARA (total mixture)",
               chr_sub)
pF <- chronic_band(chr_ara,
               "F. Chronic (time-weighted) AmsPAF — Added Risk (ARA)",
               chr_sub)

pG <- ggplot(cc$lmf_ts, aes(datetime, lmf)) +
  geom_line(colour = "grey60") + geom_point(size = 1, colour = "grey20") +
  xlim + labs(title = "G. Leachate mixing fraction (grab samples)",
              x = NULL, y = "LMF (%)") + thm
pH <- ggplot(ev, aes(date, y)) +
  geom_hline(yintercept = 1, colour = "grey80", linewidth = 0.4) +
  geom_point(aes(alpha = alpha), shape = 21, fill = "black",
             colour = "grey30", size = 3, stroke = 0.4) +
  scale_alpha_identity() +
  xlim + scale_y_continuous(limits = c(0.8, 1.2), breaks = NULL) +
  labs(title = "H. Sampling events (fill alpha = sampling completeness)",
       subtitle = "Solid = full metals + WQ suite; faint = only one or two analytes",
       x = NULL, y = NULL) + thm

g <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(pC, pD, pE, pF, pG, pH, ncol = 1)
} else {
  gridExtra::grid.arrange(pC, pD, pE, pF, pG, pH, ncol = 1)
}
ggsave("dev/bs01_kalman_compare.png", g, width = 10, height = 19, dpi = 120)
cat("WROTE dev/bs01_kalman_compare.png\n")
