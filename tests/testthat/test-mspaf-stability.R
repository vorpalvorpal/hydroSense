## Numerical stability of msPAF under ulp-level input changes
## ============================================================================
##
## CONTRACT.
##
## Chronic msPAF on the synthetic pipeline fixture used to move by ~1% (daily)
## and ~0.3% (chronic) when the input concentrations were multiplied by
## (1 + 1e-12 * U(-1, 1)).  Two mechanisms amplified ulp noise:
##   (i)  the pooled tau GAM in .fit_pooled_impact_response() (fit_pool) used
##        mgcv's factor smooth bs = "fs".  Its constructor runs nat.param() and
##        gives each of the two null-space columns its own ridge penalty; those
##        columns span a 2-fold degenerate eigenspace, so their orientation
##        rotates arbitrarily under ulp changes in x.  The profiled AIC(tau)
##        surface is therefore discontinuous (a 1e-9 change in tau moved AIC by
##        up to ~2.4 on this fixture, ~6.6 in the design study) and the Brent
##        tau optimisation (.optimise_tau_pair) jumps.
##   (ii) the marginal imputation GAMs (.fit_group_model_marginal(),
##        lv ~ s(PC1) + ... + s(PC5)) were fitted by GCV and saturate when few
##        values are detected (fixture Ni: n_det = 11, edf = 11), so their
##        predictions move ~1e-3 for a 1e-12 input change.
##
## S1. Smooth class "fsc" (registered as an S3 method of
##     mgcv::smooth.construct, so s(x, f, bs = "fsc") works in plain mgcv
##     calls): the fs factor smooth with the null space spanned in a fixed,
##     canonical orientation r * [1, (x - mean(x)) / sd(x)], with P updated to
##     match so prediction via the fs Predict.matrix reproduces the fit.
##       * model matrix stable: x -> x * (1 + 1e-12) changes the fsc model
##         matrix by < 1e-8 max abs (plain fs: up to ~4.6 on the same data);
##       * predict(fit, newdata = training data) == fitted(fit) to 1e-10;
##       * same column space as fs: near-unpenalised (sp = 1e-12) fsc and fs
##         fits agree to 1e-8 (fitted and predicted), and every fsc column is
##         in the span of the fs columns;
##       * same penalty structure as fs (one wiggliness + two null-space
##         ridges per term), so the number of smoothing parameters matches;
##       * fallback: when the base smooth's null space is not 2-D (tp with
##         m = 1 -> 1-D, m = 3 -> 3-D) fsc is plain fs: same X and S.
## S2. The pooled profiled AIC(tau_short, tau_long) is continuous: on the
##     fixture's P1 pooled call |AIC(ts, tl) - AIC(ts + d1, tl + d2)| < 1e-6
##     for steps of 1e-9 .. 4e-9 in either tau, at several anchors (today up
##     to ~2.4).
## S3. Marginal imputation GAMs are fitted by REML with a PC-count gate
##     m = min(P, floor(n_det / 3) - 1) using the FIRST m PCs (PC1..PCm);
##     m = 0 gives lv ~ 1.
##       * fixture Ni (n_det = 11, so m = 2) imputed values move < 1e-6
##         relative under a 1e-12 input perturbation (today ~8e-4);
##       * every marginal GAM has edf <= n_det - 2;
##       * n_det = 4 or 5 (m = 0) fits lv ~ 1 without error (today it aborts:
##         "fewer unique covariate combinations than ... degrees of freedom");
##       * n_det = 6..9 (m = 1 or 2) also fits without error, which needs the
##         basis dimension k of each s(PCj) kept below n_det.
## S4. End-to-end (skip_on_cran, ~3 pipeline runs): .pf_pipeline_case()'s
##     configuration with seed 1, eps = 1e-12, target and reference (and the
##     daily static-ARA reference, the same perturbed frame) values multiplied
##     by (1 + eps * U(-1, 1)): daily and chronic msPAF max relative change
##     <= 1e-6 (values with |x| > 1e-3), every tau chosen by
##     .optimise_tau_pair() moves <= 1e-6 relative, identical NA pattern.
##     Today: daily 1.0e-2, chronic 3.2e-3, tau 1.7e-1.
## S5. Invariants (guards; pass today): every optimised tau lies within its
##     bounds and tau_long >= 1.5 * tau_short (the separation
##     .optimise_tau_pair() imposes); daily and chronic msPAF lie in [0, 100];
##     results are invariant to row shuffles of target and reference to 1e-10
##     relative with the same NA pattern.  (Not bit-identical: today shuffled
##     runs differ by ~7e-16, a summation-order effect, which is why the bar
##     is 1e-10, the D2 bar of test-mspaf-robust.R.)
##
## SEAMS relied on (internal names; update here if they move):
##   * hydroSense:::.fit_pooled_impact_response(...) -- mocked once inside
##     fit_target_model() only to capture its arguments for P1;
##   * hydroSense:::.optimise_tau_pair(aic_fn, tau_bounds_short,
##     tau_bounds_long, ...) -- mocked to capture the pooled AIC closure
##     (pool_aic) and, in S4/S5, wrapped to log every optimised tau pair;
##   * hydroSense:::.fit_group_model_marginal(target_analytes, safe_analytes,
##     base, pc_wide, pc_cols, log_floors, group_name) -- called directly with
##     synthetic frames; returns $gams, one mgcv gam per analyte;
##   * fit_imputation_model(...)$groups[[g]]$gams[[analyte]] (mgcv gam).
##
## GOLDENS.  This fix deliberately changes results that depend on the pooled
## tau GAM or the marginal imputation: the end-to-end mspaf_pipeline() chronic
## golden (test-perf-goldens.R) and the shuffled-reference comparison against
## it (test-mspaf-robust.R D2) are regenerated once with the new code and
## reviewed.  Component goldens that use neither -- .compute_api(),
## time_weighted_aggregate(), .resolve_ref_norm_instant(), and the static /
## temporal-ARA mspaf_daily() point-mode golden (its 10-grab target is below
## min_obs_model, so no pooled fit runs; no imputation) -- must stay
## IDENTICAL to their goldens.
##
## All data are synthetic (helper-perf-fixtures.R and the frames built here).

library(testthat)

.ms_q <- function(expr) suppressWarnings(suppressMessages(expr))

## Relative change over values that are not ~0 (|x| > floor).
.ms_rel <- function(a, b, floor = 1e-3) {
  ok <- is.finite(a) & is.finite(b) & pmax(abs(a), abs(b)) > floor
  if (!any(ok)) return(0)
  max(abs(a[ok] - b[ok]) / abs(b[ok]))
}


# ── S1: the "fsc" factor-smooth basis ─────────────────────────────────────────

## Three factor levels, 30 points each, x well away from 0 (like hydro_long).
## Seed 2 is one where plain fs rotates its null space under x * (1 + 1e-12):
## max |dX| = 4.58 for fs vs 3e-14 for a canonical null space (mgcv 1.9-3,
## reference BLAS; other seeds give 0.27 .. 4.5, or ~1e-14 when the rotation
## happens not to trigger).
.ms_fs_data <- function(seed = 2L, n = 30L, nl = 3L) {
  withr::with_seed(seed, {
    f  <- factor(rep(letters[seq_len(nl)], each = n))
    x  <- stats::rgamma(n * nl, 2, 0.05) + 40
    x2 <- stats::runif(n * nl, 0, 30)
    y  <- sin(x / 20) + 0.5 * cos(x2 / 8) + as.numeric(f) / 3 +
      stats::rnorm(n * nl, 0, 0.3)
    data.frame(x = x, x2 = x2, f = f, y = y)
  })
}

.ms_smooth <- function(d, bs, k = 4L, m = NA) {
  spec <- mgcv::s(x, f, bs = bs, k = k, m = m)
  mgcv::smoothCon(spec, data = d, knots = NULL)[[1L]]
}

.ms_fit2 <- function(d, bs, sp = NULL) {
  form <- stats::as.formula(sprintf(
    "y ~ s(x, f, bs = '%s', k = 4) + s(x2, f, bs = '%s', k = 4)", bs, bs))
  .ms_q(mgcv::gam(form, data = d, method = "REML", sp = sp))
}

describe("S1: the fsc factor-smooth basis", {

  it("is usable in a plain mgcv::gam() call", {
    d <- .ms_fs_data()
    fit <- .ms_q(mgcv::gam(y ~ s(x, f, bs = "fsc"), data = d))
    expect_s3_class(fit, "gam")
    # Same penalty structure as fs: wiggliness + one ridge per null-space column.
    ref <- .ms_q(mgcv::gam(y ~ s(x, f, bs = "fs"), data = d))
    expect_identical(length(fit$sp), length(ref$sp))
    expect_identical(length(stats::coef(fit)), length(stats::coef(ref)))
  })

  it("has a model matrix stable to a 1e-12 relative change in x", {
    # Plain fs on seed 2: max |dX| ~ 4.6 (the null-space pair rotates).
    for (seed in 1:6) {
      d  <- .ms_fs_data(seed)
      d2 <- d
      d2$x <- d$x * (1 + 1e-12)
      dX <- max(abs(.ms_smooth(d, "fsc")$X - .ms_smooth(d2, "fsc")$X))
      expect_lt(dX, 1e-8, label = sprintf("fsc max |dX| (seed %d)", seed))
    }
  })

  it("predicts its own training data exactly (P matches the basis)", {
    d <- .ms_fs_data()
    fit <- .ms_fit2(d, "fsc")
    expect_lt(max(abs(stats::predict(fit, newdata = d) - stats::fitted(fit))),
              1e-10)
  })

  it("spans the same column space as fs", {
    d <- .ms_fs_data()
    Xs <- .ms_smooth(d, "fs")$X
    Xc <- .ms_smooth(d, "fsc")$X
    expect_identical(dim(Xc), dim(Xs))
    # Every fsc column is a linear combination of the fs columns.
    res <- stats::qr.resid(qr(Xs), Xc)
    expect_lt(max(abs(res)) / max(abs(Xc)), 1e-8)

    # Near-unpenalised fits are the same projection (differences O(sp)).
    f0 <- .ms_fit2(d, "fs",  sp = rep(1e-12, 6L))
    f1 <- .ms_fit2(d, "fsc", sp = rep(1e-12, 6L))
    expect_lt(max(abs(stats::fitted(f0) - stats::fitted(f1))), 1e-8)
    nd <- d[withr::with_seed(5L, sample.int(nrow(d), 20L)), ]
    nd$x <- nd$x * 1.05
    expect_lt(max(abs(stats::predict(f0, nd) - stats::predict(f1, nd))), 1e-8)
  })

  it("falls back to plain fs when the null space is not 2-D", {
    d <- .ms_fs_data()
    # tp penalty order m = 1 -> null space {1} (1-D); m = 3 -> {1, x, x^2}.
    for (m in c(1L, 3L)) {
      s_fs  <- .ms_smooth(d, "fs",  k = 5L, m = m)
      s_fsc <- .ms_smooth(d, "fsc", k = 5L, m = m)
      expect_false(ncol(s_fs$Xb) - s_fs$base$rank == 2L)
      expect_equal(s_fsc$X, s_fs$X, tolerance = 1e-12,
                   label = sprintf("fsc X (m = %d)", m))
      expect_equal(s_fsc$S, s_fs$S, tolerance = 1e-12,
                   label = sprintf("fsc S (m = %d)", m))
    }
  })
})


# ── S2: continuity of the pooled AIC(tau) objective ───────────────────────────

## Reproduce the pipeline's P1 pooled call exactly (the design study's
## cap_args): canonicalised inputs, imputation model on the reference (seed 1),
## the reference model with the pipeline's reference_args, then the P1 target
## model -- whose .fit_pooled_impact_response() call is captured and replayed
## with .optimise_tau_pair() mocked to hand back the pool_aic closure.
.ms_pool_aic <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    sc <- .pf_scenario()
    cz <- hydroSense:::.canonicalise_mspaf_input
    ref <- cz(sc$reference)
    tgt <- cz(sc$target)
    im <- .ms_q(fit_imputation_model(ref, seed = 1L))
    rm <- .ms_q(fit_reference_model(
      ref, hydro = sc$rain, imputation_model = im, min_obs_model = 20L,
      api_tau_bounds_short = c(3, 10), api_tau_bounds_long = c(30, 90)
    ))
    cap <- new.env()
    local({
      local_mocked_bindings(.fit_pooled_impact_response = function(...) {
        cap$args <- list(...)
        stop("captured")
      })
      try(.ms_q(fit_target_model(tgt[tgt$site_id == "P1", ], rm,
                                 imputation_model = im)), silent = TRUE)
    })
    stopifnot(!is.null(cap$args))
    pooled <- hydroSense:::.fit_pooled_impact_response
    local({
      local_mocked_bindings(.optimise_tau_pair = function(aic_fn, ...) {
        cap$aic <- aic_fn
        list(tau_short = 5, tau_long = 60)
      })
      .ms_q(do.call(pooled, cap$args))
    })
    stopifnot(is.function(cap$aic))
    cache <<- cap$aic
    cache
  }
})

describe("S2: the pooled profiled AIC(tau) is continuous in tau", {

  it("moves < 1e-6 for 1e-9 .. 4e-9 steps in either tau", {
    skip_on_cran()
    aic <- .ms_pool_aic()
    # Today the AIC range over these steps is 0.14, 0.88, 1.18, 0.20 and
    # 2.39 at the five anchors (the first is the design study's chosen tau).
    anchors <- list(c(23.758482, 156.36037), c(5, 60), c(2, 30), c(10, 120),
                    c(28, 340))
    for (a in anchors) {
      a0 <- aic(a[1], a[2])
      expect_true(is.finite(a0))
      steps <- (1:4) * 1e-9
      d_long  <- vapply(steps, function(h) aic(a[1], a[2] + h), numeric(1)) - a0
      d_short <- vapply(steps, function(h) aic(a[1] + h, a[2]), numeric(1)) - a0
      expect_lt(max(abs(d_long)), 1e-6,
                label = sprintf("|dAIC| in tau_long at (%g, %g)", a[1], a[2]))
      expect_lt(max(abs(d_short)), 1e-6,
                label = sprintf("|dAIC| in tau_short at (%g, %g)", a[1], a[2]))
    }
  })
})


# ── S3: marginal imputation GAMs ──────────────────────────────────────────────

.ms_perturb <- function(df, seed, eps = 1e-12) {
  withr::with_seed(seed, {
    df$value <- df$value * (1 + eps * stats::runif(nrow(df), -1, 1))
  })
  df
}

.ms_impute <- function(ref, tgt) {
  m <- .ms_q(fit_imputation_model(ref, seed = 1L))
  it <- withr::with_seed(7L, .ms_q(impute_chemistry(tgt, m)))
  list(model = m,
       imputed = it[order(it$sample_id, it$analyte, method = "radix"), ])
}

## Every marginal GAM of a fitted imputation model, with its n_det.
.ms_marginal_gams <- function(model) {
  out <- list()
  for (g in model$groups) {
    if (!identical(g$impute_method, "marginal")) next
    for (a in names(g$gams)) {
      out[[a]] <- list(gam = g$gams[[a]], n_det = nrow(g$gams[[a]]$model),
                       P = length(g$pc_cols), pc_cols = g$pc_cols)
    }
  }
  out
}

.ms_gate <- function(n_det, P) as.integer(min(P, floor(n_det / 3) - 1L))

## PCs a fitted marginal GAM actually uses (none for lv ~ 1).
.ms_pcs_used <- function(fit) {
  sort(grep("^PC", all.vars(stats::formula(fit)[[3L]]), value = TRUE))
}

## Synthetic single-analyte inputs for .fit_group_model_marginal().
.ms_marginal_direct <- function(n_det, P = 5L, n_bdl = 3L, seed = 1L) {
  n <- n_det + n_bdl
  sid <- sprintf("s%03d", seq_len(n))
  pcw <- withr::with_seed(seed + 1L, {
    m <- matrix(stats::rnorm(n * P), n, P,
                dimnames = list(NULL, paste0("PC", seq_len(P))))
    tibble::as_tibble(m)
  })
  pcw <- dplyr::mutate(pcw, sample_id = sid, .before = 1L)
  base <- withr::with_seed(seed, tibble::tibble(
    sample_id = sid, analyte = "Xx",
    detected = c(rep(TRUE, n_det), rep(FALSE, n_bdl)),
    lv = 0.4 * pcw$PC1 + stats::rnorm(n, 0, 0.5)
  ))
  .ms_q(hydroSense:::.fit_group_model_marginal(
    target_analytes = "Xx", safe_analytes = c(Xx = "Xx"), base = base,
    pc_wide = pcw, pc_cols = paste0("PC", seq_len(P)),
    log_floors = c(Xx = -10), group_name = "test"
  ))$gams$Xx
}

describe("S3: marginal imputation GAMs", {

  it("imputes the thin analyte (fixture Ni) stably under 1e-12 perturbation", {
    sc <- .pf_scenario()
    b <- .ms_impute(sc$reference, sc$target)
    p <- .ms_impute(.ms_perturb(sc$reference, 1L),
                    .ms_perturb(sc$target, 101L))
    expect_identical(p$imputed[c("sample_id", "analyte", "imputed")],
                     b$imputed[c("sample_id", "analyte", "imputed")])
    k <- b$imputed$imputed %in% TRUE & b$imputed$analyte == "Ni"
    expect_gt(sum(k), 10L)
    # Today: 7.7e-4.
    expect_lt(max(abs(p$imputed$value[k] - b$imputed$value[k]) /
                    abs(b$imputed$value[k])), 1e-6)
  })

  it("fits every fixture marginal GAM by REML, gated, with edf <= n_det - 2", {
    gams <- .ms_marginal_gams(.ms_impute(.pf_scenario()$reference,
                                         .pf_scenario()$target)$model)
    expect_true(all(c("Cu", "Ni", "Zn") %in% names(gams)))
    for (a in names(gams)) {
      g <- gams[[a]]
      # Today: GCV, and Ni (n_det = 11) has edf 11.0.
      expect_identical(g$gam$method, "REML", label = paste(a, "method"))
      expect_lte(sum(g$gam$edf), g$n_det - 2, label = paste(a, "edf"))
      m <- .ms_gate(g$n_det, g$P)
      expect_identical(.ms_pcs_used(g$gam), sort(g$pc_cols[seq_len(m)]),
                       label = sprintf("%s PCs used (n_det %d, P %d)",
                                       a, g$n_det, g$P))
    }
    expect_identical(.ms_pcs_used(gams$Ni$gam), c("PC1", "PC2"))
  })

  it("fits lv ~ 1 when the gate gives m = 0 (n_det = 4, 5)", {
    for (n in 4:5) {
      fit <- .ms_marginal_direct(n)
      expect_s3_class(fit, "gam")
      expect_identical(.ms_pcs_used(fit), character(0))
      expect_length(stats::coef(fit), 1L)
      expect_identical(fit$method, "REML")
    }
  })

  it("uses the first m = min(P, floor(n_det/3) - 1) PCs (n_det >= 10)", {
    cases <- rbind(c(11, 5), c(12, 5), c(15, 5), c(17, 5), c(18, 5),
                   c(44, 5), c(44, 2), c(30, 3))
    for (i in seq_len(nrow(cases))) {
      n <- cases[i, 1]
      P <- cases[i, 2]
      fit <- .ms_marginal_direct(n, P = P)
      m <- .ms_gate(n, P)
      expect_identical(.ms_pcs_used(fit), sort(paste0("PC", seq_len(m))),
                       label = sprintf("PCs used (n_det %d, P %d)", n, P))
      expect_identical(fit$method, "REML")
      expect_lte(sum(fit$edf), n - 2, label = sprintf("edf (n_det %d)", n))
    }
  })

  it("fits the gated model without error for n_det = 6..9 (m = 1, 2)", {
    # s(PCj) at mgcv's default k = 10 cannot be fitted to < 10 points, so the
    # basis dimension must be kept below n_det here.
    for (n in 6:9) {
      fit <- .ms_marginal_direct(n)
      expect_identical(.ms_pcs_used(fit),
                       paste0("PC", seq_len(.ms_gate(n, 5L))),
                       label = sprintf("PCs used (n_det %d)", n))
      expect_lte(sum(fit$edf), n - 2, label = sprintf("edf (n_det %d)", n))
    }
  })
})


# ── S4 / S5: end-to-end perturbation stability and invariants ─────────────────

## One pipeline run in .pf_pipeline_case()'s configuration, with the target
## and reference transformed by target_fn / ref_fn (the transformed reference
## is also the daily static-ARA reference).  Every .optimise_tau_pair() result
## is logged together with its bounds.
.ms_pipeline_run <- function(target_fn = identity, ref_fn = identity) {
  args <- .pf_pipeline_args()
  args$target <- target_fn(args$target)
  args$reference <- ref_fn(args$reference)
  args$daily_args$reference <- args$reference
  args$daily_args$start <- as.Date("2023-06-01")
  args$daily_args$end   <- as.Date("2023-12-12")
  tmp <- args$daily_args$temperature
  args$daily_args$temperature <- tmp[tmp$datetime >= as.Date("2023-04-20"), ]

  log <- new.env()
  log$v <- list()
  orig <- hydroSense:::.optimise_tau_pair
  out <- local({
    local_mocked_bindings(.optimise_tau_pair = function(aic_fn, tau_bounds_short,
                                                        tau_bounds_long, ...) {
      r <- orig(aic_fn, tau_bounds_short, tau_bounds_long, ...)
      log$v[[length(log$v) + 1L]] <- c(
        ts = r$tau_short, tl = r$tau_long,
        bs_lo = tau_bounds_short[1L], bs_hi = tau_bounds_short[2L],
        bl_lo = tau_bounds_long[1L], bl_hi = tau_bounds_long[2L]
      )
      r
    })
    .ms_q(do.call(mspaf_pipeline, args))
  })
  daily <- as.data.frame(attr(out, "daily"))
  chronic <- as.data.frame(out)
  list(
    tau = do.call(rbind, log$v),
    daily = daily[order(daily$site_id, daily$date, method = "radix"), ],
    chronic = chronic[order(chronic$site_id, chronic$focal_date,
                            method = "radix"), ]
  )
}

.ms_runs <- local({
  cache <- list()
  function(which) {
    if (is.null(cache[[which]])) {
      pert <- function(salt) function(df) .ms_perturb(df, 1L * 1000L + salt)
      shuf <- function(seed) function(df) {
        df[withr::with_seed(seed, sample.int(nrow(df))), , drop = FALSE]
      }
      cache[[which]] <<- switch(which,
        base     = .ms_pipeline_run(),
        perturb  = .ms_pipeline_run(pert(1L), pert(2L)),
        shuffled = .ms_pipeline_run(shuf(11L), shuf(12L))
      )
    }
    cache[[which]]
  }
})

describe("S4: mspaf_pipeline() is stable under 1e-12 input perturbation", {

  it("moves daily and chronic msPAF and tau by <= 1e-6 relative (seed 1)", {
    skip_on_cran()
    b <- .ms_runs("base")
    p <- .ms_runs("perturb")
    expect_identical(p$daily[c("site_id", "date")], b$daily[c("site_id", "date")])
    expect_identical(p$chronic[c("site_id", "focal_date")],
                     b$chronic[c("site_id", "focal_date")])
    expect_identical(is.na(p$daily$mspaf), is.na(b$daily$mspaf))
    expect_identical(is.na(p$chronic$value), is.na(b$chronic$value))
    # Today: daily 1.0e-2, chronic 3.2e-3, tau 1.7e-1.
    expect_lte(.ms_rel(p$daily$mspaf, b$daily$mspaf), 1e-6)
    expect_lte(.ms_rel(p$chronic$value, b$chronic$value), 1e-6)
    expect_identical(dim(p$tau), dim(b$tau))
    tau_cols <- c("ts", "tl")
    expect_lte(max(abs(p$tau[, tau_cols] - b$tau[, tau_cols]) /
                     b$tau[, tau_cols]), 1e-6)
  })
})

describe("S5: invariants of the pipeline result (guards)", {

  it("keeps every optimised tau within bounds with tau_long >= 1.5 tau_short", {
    skip_on_cran()
    tau <- .ms_runs("base")$tau
    expect_gt(nrow(tau), 0L)
    expect_true(all(tau[, "ts"] >= tau[, "bs_lo"] & tau[, "ts"] <= tau[, "bs_hi"]))
    expect_true(all(tau[, "tl"] >= tau[, "bl_lo"] & tau[, "tl"] <= tau[, "bl_hi"]))
    expect_true(all(tau[, "tl"] >= 1.5 * tau[, "ts"]))
  })

  it("keeps daily and chronic msPAF within [0, 100]", {
    skip_on_cran()
    b <- .ms_runs("base")
    v <- c(b$daily$mspaf, b$chronic$value)
    v <- v[!is.na(v)]
    expect_gt(length(v), 0L)
    expect_true(all(v >= 0 & v <= 100))
  })

  it("is invariant to row shuffles of target and reference", {
    skip_on_cran()
    b <- .ms_runs("base")
    s <- .ms_runs("shuffled")
    expect_identical(s$daily[c("site_id", "date")], b$daily[c("site_id", "date")])
    expect_identical(s$chronic[c("site_id", "focal_date")],
                     b$chronic[c("site_id", "focal_date")])
    expect_identical(is.na(s$daily$mspaf), is.na(b$daily$mspaf))
    expect_identical(is.na(s$chronic$value), is.na(b$chronic$value))
    expect_lte(.ms_rel(s$daily$mspaf, b$daily$mspaf, floor = 0), 1e-10)
    expect_lte(.ms_rel(s$chronic$value, b$chronic$value, floor = 0), 1e-10)
    expect_identical(dim(s$tau), dim(b$tau))
    expect_lte(max(abs(s$tau[, c("ts", "tl")] - b$tau[, c("ts", "tl")]) /
                     b$tau[, c("ts", "tl")]), 1e-10)
  })
})
