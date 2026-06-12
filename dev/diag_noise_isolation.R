#!/usr/bin/env Rscript
# Noise-source isolation for the daily AmsPAF credible band (diagnostic).
#
# Symptom (B.S01): the daily band does NOT pinch at grab anchors (IQR ~constant
# with distance to a grab) and E[AmsPAF] is inflated ~state-independently, which
# washes the event signal out of the chronic mean.
#
# Hypothesis: the per-draw impact/WQ GAM coefficient perturbation
# (.perturb_target_model -> beta ~ N(coef, Vp)) is added on top of a residual
# bridge whose anchors were built ONCE from the unperturbed fit, so dGAM is
# never re-absorbed at the grabs -> the reconstructed grab value floats by dGAM
# -> no pinch, uniform inflation (amplified one-sided by the convex SSD + floor).
#
# This script regenerates draws on a baseline+pulse window with one noise source
# disabled at a time and reports grab-day vs mid-gap band width and the mean.
#   Rscript dev/diag_noise_isolation.R

suppressMessages({ library(dplyr); devtools::load_all(".", quiet = TRUE) })

GUIDE  <- "guideline data"
options(leachatetools.guideline_dir = GUIDE)
cc <- qs2::qs_read("test data/bs01_v3_cache.qs2")
da <- cc$daily_args
## Reduce to a baseline (late 2023, deterministic ~8) + pulse (2024-09, ~44).
da$start <- as.Date("2023-06-01"); da$end <- as.Date("2024-12-31")
N <- 20L; SEED <- 42L

## Swap the GAM perturbation for an identity to disable the S1-S3 trend draws.
ident <- function(tm, perturb_reference = FALSE) tm
orig  <- leachatetools:::.perturb_target_model
gam_off <- function() assignInNamespace(".perturb_target_model", ident, "leachatetools")
gam_on  <- function() assignInNamespace(".perturb_target_model", orig,  "leachatetools")

run_draws <- function(ou_scale, grab_cv, gam) {
  if (gam) gam_on() else gam_off()
  on.exit(gam_on(), add = TRUE)
  suppressMessages(do.call(amspaf_daily, c(da, list(
    reference = da$reference_model, ndraws = N, seed = SEED, return = "draws",
    ou_scale = ou_scale, kappa = 0.5, grab_cv = grab_cv))))
}

## Deterministic centre (config-independent; point mode).
det <- suppressMessages(do.call(amspaf_daily, c(da, list(
  reference = da$reference_model))))[, c("date", "amspaf")]

configs <- tibble::tribble(
  ~name,          ~ou_scale, ~grab_cv, ~gam,
  "FULL",         1,         0.15,     TRUE,
  "no_GAM",       1,         0.15,     FALSE,
  "bridge_only",  1,         NA,       FALSE,
  "GAM_only",     0,         NA,       TRUE,
  "meas_only",    0,         0.15,     FALSE
)

summ <- function(dr, nm) {
  per <- dr |> group_by(date) |>
    summarise(dsl = first(days_since_last_sample),
              mn = mean(amspaf, na.rm = TRUE),
              med = median(amspaf, na.rm = TRUE),
              w = quantile(amspaf, .75, names = FALSE, na.rm = TRUE) -
                  quantile(amspaf, .25, names = FALSE, na.rm = TRUE),
              .groups = "drop")
  grab <- per |> filter(dsl == 0)
  gap  <- per |> filter(dsl >= 31, dsl <= 90)
  tibble::tibble(
    config       = nm,
    grabday_IQR  = round(median(grab$w), 1),
    midgap_IQR   = round(median(gap$w), 1),
    pinch_ratio  = round(median(grab$w) / median(gap$w), 2),
    mean_PAF     = round(mean(per$mn), 1),
    median_PAF   = round(median(per$med), 1)
  )
}

res <- purrr::pmap_dfr(configs, function(name, ou_scale, grab_cv, gam) {
  cat("running", name, "...\n")
  dr <- run_draws(ou_scale, if (is.na(grab_cv)) NULL else grab_cv, gam)
  summ(dr, name)
})

cat(sprintf("\nWindow %s .. %s | deterministic mean PAF = %.1f\n",
            da$start, da$end, mean(det$amspaf, na.rm = TRUE)))
cat("(pinch_ratio << 1 = band collapses at grabs = healthy; ~1 = no pinch)\n\n")
print(as.data.frame(res))
qs2::qs_save(res, "dev/diag_noise_isolation.qs2")
