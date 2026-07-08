# Route C validation harness (dev/plan-route-c.md, section "Validation").
#
# Comparisons against the real B.S01 leachate panel, reusing the 3-seed masking
# hold-out that produced the rescor_mi/cens/cens_factor numbers recorded in
# memory (`finding_impute_method_benchmark.md`):
#
#   A. impute_method = "factor" vs rescor_mi vs the MARGINAL BASELINE, on the
#      SAME masked hold-out cells: point recovery (RMSE/bias), calibration
#      (90%/50% coverage), convergence.
#      - Convergence for "factor" uses the ROTATION-INVARIANT implied covariance
#        Rhat (.route_c_convergence()), NOT max(rhat(Lambda)): a factor model is
#        only identified up to rotation of Lambda, so raw-Lambda Rhat is not a
#        meaningful gate; Sigma = Lambda Lambda' + Psi is what prediction uses.
#      - The MARGINAL BASELINE (each metal in its own single-analyte group, so
#        every analyte fits Stage-1 GAM + an INDEPENDENT censored residual, no
#        shared factor) is the always-runnable "no-borrowing" reference that
#        isolates the finding-3 benefit: factor beating marginal on the
#        conditioning check IS the borrowing payoff.
#   B. External baseline zCompositions::lrEM is attempted best-effort but is
#      INFEASIBLE on the full B.S01 panel (it needs a fully-complete column and
#      none exists); it skips gracefully. The marginal baseline in (A) is the
#      real comparison. lrEM code retained for a reduced/completable sub-panel.
#
# Lives under dev/ (gitignored via a forced add, matching plan-route-c.md's
# precedent) because it is coupled to the private "test data/monitoring.duckdb"
# and is a validation/decision script, not package code.
#
#   Usage: Rscript dev/bench-route-c.R [quick] [seeds=N]
#
# `quick` uses iter=800/chains=2 (fast sanity check, wider intervals);
# omit for the production iter=2000/chains=4 config used for the recorded
# rescor_mi/cens/cens_factor numbers. Default seeds: 1 (quick) / 3 (full).

suppressMessages({
  library(DBI); library(duckdb); library(dplyr); library(tidyr)
  devtools::load_all(".", quiet = TRUE)
})

args   <- commandArgs(TRUE)
QUICK  <- "quick" %in% args
SEEDS  <- {
  s <- sub("^seeds=", "", grep("^seeds=", args, value = TRUE))
  if (length(s)) as.integer(s) else if (QUICK) 1L else 3L
}
ITER   <- if (QUICK) 800 else 2000
WARMUP <- ITER / 2
CHAINS <- if (QUICK) 2L else 4L
CORES  <- CHAINS
BATCH  <- 15L
NDRAWS <- min(1000L, CHAINS * (ITER - WARMUP))
DB     <- "test data/monitoring.duckdb"

.NO_CONVERT <- c("pH", "Temperature", "Hardness-total-CaCO3", "ORP", "DO", "EC",
                 "SAR", "Stage", "Water Height-10min", "Precipitation-10min",
                 "Precipitation-1hr", "Precipitation-24hr", "Temperature-10min",
                 "Wind Speed-10min", "Wind Direction-10min", "temperature")
.DERIVED_ANALYTES <- c("AmsPAF", "msPAF", "LMF")

# ── Load B.S01 chemistry (long) — identical to test data/bench_impute_method.R
load_chem <- function(feature = "B.S01") {
  con <- dbConnect(duckdb(), DB, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  uuid <- dbGetQuery(con, sprintf(
    "SELECT uuid FROM feature WHERE name = '%s'", feature))$uuid
  raw <- dbGetQuery(con, sprintf("
    SELECT an.uuid_sample AS sample_id, s.datetime AS datetime, a.name AS analyte,
           an.value AS value, an.quantified AS detected, an.rl_low AS rl_low
    FROM analysis an
    JOIN sample s ON an.uuid_sample = s.uuid
    JOIN lab_method lm ON an.uuid_lab = lm.uuid
    JOIN analyte a ON lm.uuid_analyte = a.uuid
    WHERE s.uuid_feature = '%s'", uuid))
  raw |>
    group_by(sample_id, analyte) |> slice(1L) |> ungroup() |>
    mutate(
      value    = if_else(!analyte %in% .NO_CONVERT, value * 1000, value),
      rl_low   = if_else(!analyte %in% .NO_CONVERT, rl_low * 1000, rl_low),
      detected = as.logical(detected),
      datetime = as.Date(datetime),
      analyte  = case_when(analyte == "Temperature" ~ "temperature",
                           analyte == "Hardness-total-CaCO3" ~ "hardness",
                           TRUE ~ analyte),
      site_id  = feature
    ) |>
    filter(!analyte %in% .DERIVED_ANALYTES)
}

chem0 <- load_chem("B.S01")
ROUTINE <- intersect(c("As", "Cr", "Cu", "Ni", "Pb", "Zn"),
                     intersect(.METAL_ANALYTES, unique(chem0$analyte)))

keep_analytes <- union(c("pH", "EC", "NH3-N", .WQ_BLOCK_CANDIDATES), ROUTINE)
chem0 <- dplyr::filter(chem0, analyte %in% keep_analytes)
rsamp <- chem0 |> filter(analyte %in% ROUTINE) |> distinct(sample_id) |> pull(sample_id)
chem0 <- chem0 |> filter(sample_id %in% rsamp)

GROUPS <- list(impute_group("metals", targets = ROUTINE, hurdle = ROUTINE))

cat(sprintf("B.S01: %d samples, routine metals: %s\n",
            length(rsamp), paste(ROUTINE, collapse = ",")))

# ── Masking (identical scheme to the recorded rescor_mi/cens/cens_factor run) ─
mask_chem <- function(seed) {
  set.seed(seed)
  det  <- chem0 |> filter(analyte %in% ROUTINE, detected)
  mask <- det |> slice_sample(prop = 0.10) |>
    transmute(sample_id, analyte, truth = log10(pmax(value, 1e-6)))
  mkey <- paste(mask$sample_id, mask$analyte)
  masked <- chem0 |>
    mutate(.mk = paste(sample_id, analyte) %in% mkey,
           value    = if_else(.mk, value * 2, value),     # DL headroom
           detected = if_else(.mk, FALSE, detected)) |>
    select(-.mk)
  list(masked = masked, mask = mask)
}

# ── A. impute_method comparison (brms methods + factor) ──────────────────────
run_one <- function(meth, seed) {
  mc <- mask_chem(seed)
  # "marginal" is now a first-class impute_method: per-analyte censored GAM with
  # a Student-t posterior predictive (parameter + residual-variance uncertainty),
  # no cross-analyte borrowing. (It fits per-analyte GAMs internally, so it takes
  # the normal metals group, not one group per analyte.)
  is_marg    <- meth == "marginal"
  groups_use <- GROUPS
  method_use <- meth

  t_fit <- system.time(
    m <- fit_imputation_model(mc$masked, required_vars = c("pH", "EC"),
                              groups = groups_use, impute_method = method_use,
                              iter = ITER, warmup = WARMUP,
                              chains = CHAINS, cores = CORES)
  )[["elapsed"]]

  # Convergence: for "factor" use the rotation-INVARIANT Sigma Rhat (raw Lambda
  # Rhat is not a meaningful gate); brms methods use full Rhat; the all-
  # degenerate marginal baseline has no Stan fit to diagnose.
  if (is_marg) {
    max_rhat <- NA_real_; ndiv <- NA_real_
  } else if (method_use == "factor") {
    conv     <- hydroSense:::.route_c_convergence(m$groups$metals)
    max_rhat <- conv$sigma_rhat
    ndiv <- tryCatch(sum(m$groups$metals$fit$diagnostic_summary(quiet = TRUE)$num_divergent),
                     error = function(e) NA_real_)
  } else {
    fit      <- m$groups$metals$fit
    max_rhat <- suppressWarnings(max(brms::rhat(fit), na.rm = TRUE))
    ndiv <- tryCatch(sum(subset(brms::nuts_params(fit),
                        Parameter == "divergent__")$Value), error = function(e) NA_real_)
  }
  imp <- impute_chemistry(mc$masked, m, return = "draws", bdl_cap = FALSE,
                          ndraws = NDRAWS, batch_size = BATCH)
  s <- imp |> inner_join(mc$mask, by = c("sample_id", "analyte")) |>
    mutate(pl = log10(pmax(value, 1e-6))) |>
    group_by(sample_id, analyte, truth) |>
    summarise(pm   = mean(pl),
              lo90 = quantile(pl, 0.05), hi90 = quantile(pl, 0.95),
              lo50 = quantile(pl, 0.25), hi50 = quantile(pl, 0.75),
              .groups = "drop")
  tibble(config = meth, seed = seed, n = nrow(s),
         t_fit_s = round(t_fit, 1),
         max_rhat = round(max_rhat, 3), n_divergent = ndiv,
         err   = s$pm - s$truth,
         in90  = s$truth >= s$lo90 & s$truth <= s$hi90,
         in50  = s$truth >= s$lo50 & s$truth <= s$hi50,
         w90   = s$hi90 - s$lo90)
}

# ── B. Route D benchmark: zCompositions on Stage-1 residuals ─────────────────
# lrEM/lrDA carry no WQ-predictor mean, so run them on the SAME Stage-1 GAM
# residuals the factor model conditions on, for a like-for-like comparison
# (plan-route-c.md "Route D benchmark (zCompositions)").
run_route_d <- function(seed) {
  if (!requireNamespace("zCompositions", quietly = TRUE)) {
    message("zCompositions not installed — skipping Route D benchmark.")
    return(NULL)
  }
  mc <- mask_chem(seed)
  # A factor-method fit gives us the Stage-1 gams for free.
  m <- fit_imputation_model(mc$masked, required_vars = c("pH", "EC"),
                            groups = GROUPS, impute_method = "factor",
                            iter = ITER, warmup = WARMUP,
                            chains = CHAINS, cores = CORES)
  grp <- m$groups$metals
  pca_scores <- hydroSense:::.compute_pca_scores(mc$masked, m$pca)

  # Residual matrix: rows = samples that measured >=1 routine metal, columns =
  # ROUTINE metals, entries = lv - mu_j(X_i) (detected) with BDL cells at their
  # censoring bound (log(DL) - mu_j), matching Stage 1 exactly.
  target_vals <- mc$masked |>
    dplyr::filter(analyte %in% ROUTINE) |>
    dplyr::select(sample_id, analyte, value, detected) |>
    dplyr::group_by(sample_id, analyte) |>
    dplyr::slice(1L) |> dplyr::ungroup() |>
    dplyr::mutate(lv = log(pmax(value, grp$log_floors[analyte])))

  resid_mat <- matrix(NA_real_, nrow = length(rsamp), ncol = length(ROUTINE),
                      dimnames = list(rsamp, ROUTINE))
  mu_mat    <- matrix(NA_real_, nrow = length(rsamp), ncol = length(ROUTINE),
                      dimnames = list(rsamp, ROUTINE))
  cens_mat  <- matrix(FALSE, nrow = length(rsamp), ncol = length(ROUTINE),
                      dimnames = list(rsamp, ROUTINE))
  for (a in ROUTINE) {
    rows_a <- target_vals[target_vals$analyte == a, ]
    pc_a   <- dplyr::filter(pca_scores, sample_id %in% rows_a$sample_id)
    mu_a   <- as.numeric(stats::predict(grp$gams[[a]], newdata = pc_a, type = "response"))
    idx    <- match(pc_a$sample_id, rsamp)
    lv_a   <- rows_a$lv[match(pc_a$sample_id, rows_a$sample_id)]
    resid_mat[idx, a] <- lv_a - mu_a
    mu_mat[idx, a]    <- mu_a
    cens_mat[idx, a]  <- !rows_a$detected[match(pc_a$sample_id, rows_a$sample_id)]
  }
  keep <- rowSums(!is.na(resid_mat)) >= 2L   # zCompositions needs >=2 parts/row
  resid_mat <- resid_mat[keep, , drop = FALSE]
  mu_mat    <- mu_mat[keep, , drop = FALSE]
  cens_mat  <- cens_mat[keep, , drop = FALSE]

  # lrEM/lrDA are compositional: they expect strictly positive "parts" and do
  # their own log-ratio transform internally, so they cannot take the (signed)
  # log-scale residual directly. Exponentiate first (parts = exp(residual), a
  # positive multiplicative deviation from the GAM mean); the censoring bound
  # goes through the same transform. Convert the imputed part back via log()
  # + mu_j to land back on the natural-log concentration scale for scoring.
  dl_mat <- matrix(0, nrow = nrow(resid_mat), ncol = ncol(resid_mat),
                   dimnames = dimnames(resid_mat))
  X <- exp(resid_mat)
  X[cens_mat] <- NA_real_                    # lrEM wants the *unobserved* value as NA
  dl_mat[cens_mat] <- exp(resid_mat[cens_mat])   # ... and the bound in a parallel matrix
  # True-missing (never measured in this masked draw, not part of the mask):
  # median-fill so lrEM's log-ratio geometry has a complete matrix to work on.
  fill <- stats::median(X, na.rm = TRUE)
  X[is.na(X) & !cens_mat] <- fill

  # B.S01's ragged co-observation (median ~6 metals/sample — the very reason
  # the plan rejects a dense full-Sigma model) means several columns/rows are
  # >80% BDL; zCompositions' default z.delete would silently drop them. Keep
  # everything (z.delete = FALSE) so the comparison covers the same cells the
  # factor model was scored on, even where lrEM's own ALR geometry strains.
  point <- tryCatch(
    zCompositions::lrEM(X, label = NA, dl = dl_mat, rob = FALSE, ini.cov = "multRepl",
                        z.warning = 1, z.delete = FALSE),
    error = function(e) { message("lrEM failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(point)) return(NULL)

  mc$mask |>
    dplyr::filter(sample_id %in% rownames(point)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      pred_log10 = (mu_mat[sample_id, analyte] + log(point[sample_id, analyte])) / log(10)
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(config = "route_d_lrEM", seed = seed,
                     err = pred_log10 - truth)
}

# ── Run everything, pool, summarise ───────────────────────────────────────────
grid <- expand.grid(meth = c("rescor_mi", "marginal", "factor"),
                    seed = seq_len(SEEDS), stringsAsFactors = FALSE)
cells <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  meth <- grid$meth[i]; seed <- grid$seed[i]
  cat(sprintf("\n#### %s | seed %d ####\n", meth, seed))
  tryCatch(run_one(meth, seed),
           error = function(e) { message(meth, "/", seed, " FAILED: ",
                                          conditionMessage(e)); NULL })
}))

route_d <- bind_rows(lapply(seq_len(SEEDS), function(sd) {
  cat(sprintf("\n#### route_d_lrEM | seed %d ####\n", sd))
  tryCatch(run_route_d(sd), error = function(e) {
    message("route_d/", sd, " FAILED: ", conditionMessage(e)); NULL
  })
}))

summ <- cells |>
  group_by(config) |>
  summarise(n_cells = n(), seeds = n_distinct(seed),
            rmse = sqrt(mean(err^2)), mae = mean(abs(err)), bias = mean(err),
            cov90 = mean(in90), cov50 = mean(in50), width90 = mean(w90),
            max_rhat = max(max_rhat), n_div_max = max(n_divergent),
            t_fit_s = mean(t_fit_s), .groups = "drop") |>
  arrange(rmse) |>
  mutate(across(c(rmse, mae, bias, cov90, cov50, width90, t_fit_s), \(x) round(x, 3)))

summ_d <- if (nrow(route_d) > 0) {
  route_d |> group_by(config) |>
    summarise(n_cells = n(), seeds = n_distinct(seed),
              rmse = sqrt(mean(err^2)), mae = mean(abs(err)), bias = mean(err),
              .groups = "drop") |>
    mutate(across(c(rmse, mae, bias), \(x) round(x, 3)))
} else NULL

cat("\n================ Route C validation — B.S01 ================\n")
cat("\n-- A. impute_method (point recovery / calibration / convergence) --\n")
print(as.data.frame(summ), row.names = FALSE)
cat("\n-- B. Route D benchmark (zCompositions::lrEM on Stage-1 residuals) --\n")
if (!is.null(summ_d)) print(as.data.frame(summ_d), row.names = FALSE) else
  cat("(skipped — zCompositions not installed or lrEM failed)\n")
cat("\nPooled over", SEEDS, "mask seed(s). Lower RMSE/MAE = better recovery;",
    "cov90 target ~0.90. NOTE: for 'factor' max_rhat is the rotation-INVARIANT",
    "Sigma Rhat (.route_c_convergence), not raw-Lambda Rhat; 'marginal' has no",
    "Stan fit (NA). Win condition (plan-route-c.md): factor rmse ~= rescor_mi,",
    "cov90 ~= 0.90, Sigma Rhat clean (<1.1), no high divergence rate, and factor",
    "beats the marginal baseline on the finding-3 conditioning check.\n")

saveRDS(list(summ = summ, cells = cells, summ_d = summ_d, route_d = route_d,
             routine = ROUTINE, seeds = SEEDS, iter = ITER, chains = CHAINS),
        sprintf("dev/bench_route_c_results%s.rds", if (QUICK) "_quick" else ""))
