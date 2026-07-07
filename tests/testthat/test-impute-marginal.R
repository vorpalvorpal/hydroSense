## ============================================================================
## BDD specs for impute_method = "marginal"
## ============================================================================
##
## The marginal method: a per-analyte left-censored GAM (Stage-1 spline mean on
## PC scores + posterior-predictive residual), with NO cross-analyte factor. It
## is the "no-borrowing" method — accurate and fast on panels where cross-metal
## correlation is weak/ragged (e.g. B.S01), and it needs NEITHER brms NOR
## cmdstanr (mgcv only). Its uncertainty is a proper posterior predictive
## (GAM parameter uncertainty via Vp + residual variance), not the under-
## dispersed plug-in the factor method's degenerate path uses.
##
## Written test-first; drives the implementation red -> green.

library(testthat)
library(hydroSense)


describe("impute_method = 'marginal' fitting", {

  it("fits per-analyte GAMs needing neither brms nor cmdstanr", {
    ## The whole point: no Stan backend. Mock BOTH engine guards to abort; a
    ## marginal fit must still succeed (it only uses mgcv).
    local_mocked_bindings(
      .require_brms     = function() cli::cli_abort("brms must not be required"),
      .require_cmdstanr = function() cli::cli_abort("cmdstanr must not be required")
    )
    df <- .imp_chem(n = 40, seed = 1)
    model <- fit_imputation_model(
      df, impute_method = "marginal",
      groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                                 hurdle = c("Cu", "Zn")))
    )
    grp <- model$groups[["metals"]]
    expect_identical(model$impute_method, "marginal")
    expect_type(grp$gams, "list")
    expect_setequal(names(grp$gams), c("Cu", "Zn"))
    expect_true(all(vapply(grp$gams, inherits, logical(1), "gam")))
    expect_null(grp$fit)   # no Stan fit object
  })
})


describe("impute_method = 'marginal' prediction", {

  it("truncates imputed BDL cells at the detection limit, with no cap", {
    df <- .imp_chem(n = 60, seed = 2)
    dl <- 0.5
    bdl_ids <- paste0("s", 1:10)
    df$detected[df$analyte == "Zn" & df$sample_id %in% bdl_ids] <- FALSE
    df$value[df$analyte == "Zn" & df$sample_id %in% bdl_ids]    <- dl

    model <- fit_imputation_model(
      df, impute_method = "marginal",
      groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                                 hurdle = c("Cu", "Zn")))
    )
    out <- impute_chemistry(df, model, return = "draws", ndraws = 500)
    zn_bdl <- out[out$analyte == "Zn" & out$sample_id %in% bdl_ids &
                    out$imputed_kind == "censored_left", ]
    expect_true(nrow(zn_bdl) > 0)
    expect_true(all(zn_bdl$value <= dl + 1e-9))
    expect_null(bdl_cap_summary(out))   # truncated by construction, cap never fires
  })

  it("caps point-mode BDL estimates at the DL too", {
    df <- .imp_chem(n = 60, seed = 2)
    dl <- 0.5
    bdl_ids <- paste0("s", 1:10)
    df$detected[df$analyte == "Zn" & df$sample_id %in% bdl_ids] <- FALSE
    df$value[df$analyte == "Zn" & df$sample_id %in% bdl_ids]    <- dl
    model <- fit_imputation_model(
      df, impute_method = "marginal",
      groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                                 hurdle = c("Cu", "Zn")))
    )
    out <- impute_chemistry(df, model, return = "point")
    zn_bdl <- out[out$analyte == "Zn" & out$sample_id %in% bdl_ids &
                    out$imputed_kind == "censored_left", ]
    expect_true(nrow(zn_bdl) > 0)
    expect_true(all(zn_bdl$value <= dl + 1e-9))
  })

  it("does NOT borrow across analytes (the anti-finding-3)", {
    ## Marginal is independent by construction: a co-measured metal must not
    ## change the imputed one (contrast the factor method, where it does).
    set.seed(3); n <- 80
    df <- .imp_chem(n = n, seed = 3)
    fac <- stats::rnorm(n)   # even with real Cu-Zn coupling in TRAINING...
    df$value[df$analyte == "Cu"] <- exp(1.0 * fac + stats::rnorm(n, 0, 0.4))
    df$value[df$analyte == "Zn"] <- exp(1.0 * fac + stats::rnorm(n, 0, 0.4))
    model <- fit_imputation_model(
      df, impute_method = "marginal",
      groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                                 hurdle = c("Cu", "Zn")))
    )
    ctx <- dplyr::bind_rows(
      .imp_rows("q1", "pH", 7), .imp_rows("q1", "EC", 400),
      .imp_rows("q1", "Cl", 500)
    )
    lo <- dplyr::bind_rows(ctx, .imp_rows("q1", "Cu", 1))
    hi <- dplyr::bind_rows(ctx, .imp_rows("q1", "Cu", 100))
    v_lo <- impute_chemistry(lo, model, return = "point") |>
      (\(o) o$value[o$sample_id == "q1" & o$analyte == "Zn"])()
    v_hi <- impute_chemistry(hi, model, return = "point") |>
      (\(o) o$value[o$sample_id == "q1" & o$analyte == "Zn"])()
    ## ...the imputed Zn is identical: it depends only on q1's PC scores, which
    ## are the same in lo and hi. (Cu is not a PCA predictor.)
    expect_equal(v_lo, v_hi, tolerance = 1e-6)
  })

  it("propagates GAM parameter uncertainty, not just residual noise", {
    ## The calibration upgrade: predictive variance for an imputed cell must
    ## equal Xp Vp Xp' (mean uncertainty) + sig2 (residual), NOT sig2 alone.
    df <- .imp_chem(n = 40, seed = 5)
    ## Drop Zn for a handful of samples so they need a (truncation-free) missing
    ## imputation — clean variance decomposition, no truncation to distort it.
    miss_ids <- paste0("s", 1:8)
    df <- df[!(df$analyte == "Zn" & df$sample_id %in% miss_ids), ]
    model <- fit_imputation_model(
      df, impute_method = "marginal",
      groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                                 hurdle = c("Cu", "Zn")))
    )
    out <- impute_chemistry(df, model, return = "draws", ndraws = 4000)

    sid  <- miss_ids[1]
    draws <- out$value[out$sample_id == sid & out$analyte == "Zn" & out$imputed]
    expect_gt(length(draws), 1000)          # a distribution
    v <- stats::var(log(draws))             # draws are exp(eta + eps); log => eta + eps

    gam_zn <- model$groups[["metals"]]$gams[["Zn"]]
    sig2   <- gam_zn$sig2
    pc  <- hydroSense:::.compute_pca_scores(df, model$pca)
    pc1 <- pc[pc$sample_id == sid, , drop = FALSE]
    Xp  <- stats::predict(gam_zn, newdata = pc1, type = "lpmatrix")
    param_var <- as.numeric(Xp %*% stats::vcov(gam_zn, unconditional = TRUE) %*% t(Xp))

    expect_gt(param_var, 0)                                 # mean uncertainty is real
    expect_equal(v, param_var + sig2, tolerance = 0.2)      # draws = param + residual
  })
})
