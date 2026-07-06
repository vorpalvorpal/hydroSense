## ============================================================================
## BDD specifications for Route C — low-rank censored factor imputation
## ============================================================================
##
## Executable contract for dev/plan-route-c.md. Written test-first: every spec
## describes target behaviour and is guarded with `.skip_route_c()` until the
## corresponding piece is built. Implement in the order below — the prediction
## kernels (§A) are Stan-free and encode findings 1-3, so they can be driven
## green before any Stan is written.
##
## Target seams (created by the implementation; see the plan):
##   .factor_condition(y, mu, Lambda, psi)
##       -> list(mean, cov) : Gaussian conditional of the missing analytes
##          (y == NA) given the observed ones, under N(mu, Lambda Lambda' + Psi).
##          May use the cheap k x k Woodbury form; the ANSWER must match the
##          brute-force full-Sigma reference (.rc_ref_conditional).
##   .factor_condition_draw(y, mu, Lambda, psi, ndraws, upper = NULL)
##       -> ndraws x m matrix of draws for the missing block; entries with a
##          finite `upper` (a BDL target's log detection limit) are drawn from
##          the left-truncated Normal so no draw exceeds the bound.
##   .fit_group_model(..., impute_method = "factor")
##       -> group with $Lambda (JxK), $Psi (length J), $gams (named), $k, and
##          the existing $safe_names / $analytes fields.
##   impute_chemistry(df, model)  with model$impute_method == "factor"
##       -> end-to-end; the DL cap is a no-op (draws respect the bound).

library(testthat)
library(hydroSense)


## ----------------------------------------------------------------------------
## §A  Prediction kernel — the finding-3 conditioning (Stan-free)
## ----------------------------------------------------------------------------
describe("factor conditional prediction kernel (.factor_condition)", {

  it("matches the full-Sigma Gaussian conditional (findings 1 & 3)", {
    .skip_route_c(".factor_condition()")

    fm  <- .rc_factor(J = 4L, k = 2L, seed = 1L)
    ## Observe Cu (idx 1) and Cd (idx 3); impute Zn (2) and Pb (4).
    y <- c(fm$mu[["Cu"]] + 0.9, NA, fm$mu[["Cd"]] - 0.4, NA)

    got <- hydroSense:::.factor_condition(y, fm$mu, fm$Lambda, fm$psi)
    ref <- .rc_ref_conditional(y, fm$mu, fm$Lambda, fm$psi)

    expect_equal(as.numeric(got$mean), as.numeric(ref$mean), tolerance = 1e-8)
    expect_equal(matrix(got$cov, 2, 2), matrix(ref$cov, 2, 2), tolerance = 1e-8)
  })

  it("returns the marginal when no analyte is observed (|O| = 0)", {
    .skip_route_c(".factor_condition()")

    fm <- .rc_factor(J = 4L, k = 2L, seed = 2L)
    y  <- rep(NA_real_, 4L)   # nothing observed

    got   <- hydroSense:::.factor_condition(y, fm$mu, fm$Lambda, fm$psi)
    Sigma <- fm$Lambda %*% t(fm$Lambda) + diag(fm$psi)

    expect_equal(as.numeric(got$mean), as.numeric(fm$mu), tolerance = 1e-8)
    expect_equal(unname(matrix(got$cov, 4, 4)), unname(Sigma), tolerance = 1e-8)
  })

  it("is invariant to rotation of the loadings", {
    .skip_route_c(".factor_condition()")

    fm <- .rc_factor(J = 4L, k = 2L, seed = 3L)
    y  <- c(fm$mu[["Cu"]] + 0.7, NA, NA, fm$mu[["Pb"]] + 0.3)

    theta <- 0.6
    R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
    Lambda_rot <- fm$Lambda %*% R   # ΛR: same ΛΛ', so predictions must be identical

    a <- hydroSense:::.factor_condition(y, fm$mu, fm$Lambda,     fm$psi)
    b <- hydroSense:::.factor_condition(y, fm$mu, Lambda_rot,    fm$psi)

    expect_equal(unname(as.numeric(a$mean)), unname(as.numeric(b$mean)), tolerance = 1e-8)
    expect_equal(unname(as.matrix(a$cov)),   unname(as.matrix(b$cov)),   tolerance = 1e-8)
  })

  it("is stable and correct when only one analyte is observed (|O| = 1)", {
    .skip_route_c(".factor_condition()")

    fm <- .rc_factor(J = 4L, k = 2L, seed = 4L)
    y  <- c(fm$mu[["Cu"]] + 1.2, NA, NA, NA)  # only Cu observed

    got <- hydroSense:::.factor_condition(y, fm$mu, fm$Lambda, fm$psi)
    ref <- .rc_ref_conditional(y, fm$mu, fm$Lambda, fm$psi)

    expect_true(all(is.finite(got$mean)))
    expect_equal(as.numeric(got$mean), as.numeric(ref$mean), tolerance = 1e-8)
  })

  it("moves an unobserved analyte in the direction of the observed residual", {
    .skip_route_c(".factor_condition()")

    ## Single positive-loading factor: every metal loads +1, so a high observed
    ## Cu residual must raise the conditional mean of the unobserved Zn.
    analytes <- c("Cu", "Zn")
    Lambda <- matrix(c(1, 1), nrow = 2, dimnames = list(analytes, "F1"))
    psi    <- stats::setNames(c(0.1, 0.1), analytes)
    mu     <- stats::setNames(c(0, 0), analytes)

    lo <- hydroSense:::.factor_condition(c(-1, NA), mu, Lambda, psi)$mean
    hi <- hydroSense:::.factor_condition(c( 2, NA), mu, Lambda, psi)$mean

    expect_gt(hi, lo)
  })
})


## ----------------------------------------------------------------------------
## §A' Draws + BDL truncation (Stan-free)
## ----------------------------------------------------------------------------
describe("factor conditional draws (.factor_condition_draw)", {

  it("never lets a BDL target draw exceed its detection limit (finding 2)", {
    .skip_route_c(".factor_condition_draw()")

    fm <- .rc_factor(J = 4L, k = 2L, seed = 5L)
    y  <- c(fm$mu[["Cu"]] + 0.8, NA, NA, NA)   # observe Cu; impute Zn, Cd, Pb
    ## `upper` is ordered by the missing block (Zn, Cd, Pb): Zn is a BDL target
    ## with a tight detection limit; Cd/Pb are truly missing (unbounded).
    upper <- c(log(0.5), NA, NA)

    draws <- hydroSense:::.factor_condition_draw(
      y, fm$mu, fm$Lambda, fm$psi, ndraws = 500L, upper = upper
    )
    expect_equal(dim(draws), c(500L, 3L))
    expect_true(all(draws[, 1] <= log(0.5) + 1e-9))   # Zn respects the DL
  })

  it("centres unbounded draws on the conditional mean", {
    .skip_route_c(".factor_condition_draw()")

    fm <- .rc_factor(J = 4L, k = 2L, seed = 6L)
    y  <- c(fm$mu[["Cu"]] + 0.5, NA, fm$mu[["Cd"]] - 0.2, NA)

    cmean <- hydroSense:::.factor_condition(y, fm$mu, fm$Lambda, fm$psi)$mean
    draws <- hydroSense:::.factor_condition_draw(
      y, fm$mu, fm$Lambda, fm$psi, ndraws = 4000L, upper = c(NA, NA)
    )
    expect_equal(colMeans(draws), as.numeric(cmean), tolerance = 0.1)
  })
})


## ----------------------------------------------------------------------------
## §B  Fit — the two-stage censored factor model (needs cmdstanr)
## ----------------------------------------------------------------------------
describe("factor model fit (.fit_group_model, impute_method = 'factor')", {

  it("returns Lambda, Psi, per-analyte gams and k", {
    skip_if_not_installed("cmdstanr")
    .skip_route_c("the 'factor' branch of .fit_group_model()")

    df    <- .imp_chem(n = 60, seed = 1)
    model <- fit_imputation_model(
      df, impute_method = "factor",
      groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                                 hurdle = c("Cu", "Zn"))),
      iter = 400, warmup = 200, chains = 2
    )
    grp <- model$groups[["metals"]]

    expect_equal(nrow(grp$Lambda), length(grp$analytes))
    expect_true(is.numeric(grp$Psi) && length(grp$Psi) == length(grp$analytes))
    expect_type(grp$gams, "list")
    expect_true(all(vapply(grp$gams, inherits, logical(1), "gam")))
    expect_true(is.integer(grp$k) || is.numeric(grp$k))
  })

  it("rejects analyte names that collide under make.names()", {
    skip_if_not_installed("cmdstanr")
    .skip_route_c("the 'factor' branch of .fit_group_model()")

    df <- .imp_chem(n = 40, seed = 2)
    expect_error(
      fit_imputation_model(
        df, impute_method = "factor",
        groups = list(impute_group("g", targets = c("Cr-6", "Cr.6"),
                                   hurdle = c("Cr-6", "Cr.6")))
      ),
      regexp = "collid|unique|safe", ignore.case = TRUE
    )
  })

  it("converges with benign geometry (no rescor_mi-style funnel)", {
    skip_if_not_installed("cmdstanr")
    .skip_route_c("the 'factor' branch of .fit_group_model()")

    df    <- .imp_chem(n = 80, seed = 3)
    model <- fit_imputation_model(
      df, impute_method = "factor",
      groups = list(impute_group("metals", targets = c("Cu", "Zn", "Cd"),
                                 hurdle = c("Cu", "Zn", "Cd"))),
      iter = 800, warmup = 400, chains = 2
    )
    ## The low-rank structure exists precisely to avoid rescor_mi's R-hat ~ 1.6.
    rhats <- brms::rhat(model$groups[["metals"]]$fit)
    expect_lt(max(rhats, na.rm = TRUE), 1.1)
  })
})


## ----------------------------------------------------------------------------
## §C  End-to-end impute_chemistry() with the factor method (needs cmdstanr)
## ----------------------------------------------------------------------------
describe("impute_chemistry(impute_method = 'factor')", {

  it("imputes BDL cells without exceeding the DL, and needs no cap (findings 1-2)", {
    skip_if_not_installed("cmdstanr")
    .skip_route_c("end-to-end factor prediction")

    df <- .imp_chem(n = 60, seed = 1)
    ## Make some Zn cells BDL at a known limit.
    dl <- 0.5
    bdl_ids <- paste0("s", 1:10)
    df$detected[df$analyte == "Zn" & df$sample_id %in% bdl_ids] <- FALSE
    df$value[df$analyte == "Zn" & df$sample_id %in% bdl_ids]    <- dl

    model <- fit_imputation_model(
      df, impute_method = "factor",
      groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                                 hurdle = c("Cu", "Zn"))),
      iter = 400, warmup = 200, chains = 2
    )
    out <- impute_chemistry(df, model, return = "draws")

    imp_zn_bdl <- out[out$analyte == "Zn" & out$sample_id %in% bdl_ids &
                        out$imputed_kind == "censored_left", ]
    expect_true(all(imp_zn_bdl$value <= dl + 1e-9))
    ## The cap should never have fired for the factor method.
    expect_null(bdl_cap_summary(out))
  })

  it("lets co-measured metals inform an imputed missing metal (finding 3)", {
    skip_if_not_installed("cmdstanr")
    .skip_route_c("end-to-end factor prediction")

    df <- .imp_chem(n = 80, seed = 7)
    model <- fit_imputation_model(
      df, impute_method = "factor",
      groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                                 hurdle = c("Cu", "Zn"))),
      iter = 500, warmup = 250, chains = 2
    )

    ## One new sample, Zn absent (to be imputed), observed Cu at two levels.
    ctx <- dplyr::bind_rows(
      .imp_rows("q1", "pH", 7), .imp_rows("q1", "EC", 400),
      .imp_rows("q1", "Cl", 500)
    )
    lo <- dplyr::bind_rows(ctx, .imp_rows("q1", "Cu", 1))
    hi <- dplyr::bind_rows(ctx, .imp_rows("q1", "Cu", 100))

    zn_lo <- impute_chemistry(lo, model, return = "point")
    zn_hi <- impute_chemistry(hi, model, return = "point")
    v_lo <- zn_lo$value[zn_lo$sample_id == "q1" & zn_lo$analyte == "Zn"]
    v_hi <- zn_hi$value[zn_hi$sample_id == "q1" & zn_hi$analyte == "Zn"]

    ## If Cu and Zn co-vary positively, the higher-Cu sample must impute a
    ## higher Zn — the conditioning rescor_mi never delivered at predict time.
    expect_false(isTRUE(all.equal(v_lo, v_hi)))
  })

  it("propagates uncertainty through return = 'draws'", {
    skip_if_not_installed("cmdstanr")
    .skip_route_c("end-to-end factor prediction")

    df <- .imp_chem(n = 60, seed = 9)
    model <- fit_imputation_model(
      df, impute_method = "factor",
      groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                                 hurdle = c("Cu", "Zn"))),
      iter = 400, warmup = 200, chains = 2
    )
    out <- impute_chemistry(df, model, return = "draws")
    imp <- out[out$imputed & out$analyte == "Zn", ]
    per_cell <- tapply(imp$value, imp$sample_id, function(v) length(unique(v)))
    expect_true(all(per_cell > 1))   # a distribution, not a point
  })
})
