## ============================================================================
## BDD specifications derived from the imputation-engine code review
## ============================================================================
##
## Scope: the review findings that SURVIVE the planned "Route C" rewrite of the
## fitting/prediction path (a low-rank censored factor model). Route C resolves
## findings 1-3 (the default method, the per-draw DL cap, and prediction-time
## cross-analyte conditioning) by construction, so they are NOT encoded here —
## their contracts belong with the Route C model. Everything below lives in the
## shared predictor-building / hurdle / co-analyte machinery that Route C keeps,
## and should be hardened first (findings 5 and 6 also improve the PC-score
## inputs Route C will consume).
##
## These are written test-first: each spec encodes the *target* (post-fix)
## behaviour and is expected to FAIL against the current implementation, so each
## is guarded by `.skip_tdd()` (from helper-lmf-review.R) to keep the suite
## green until the fix lands. To drive a fix red-green: delete the `.skip_tdd()`
## line, watch it fail, implement, watch it pass.
##
## Every spec here is brms/Stan-free (it exercises deterministic helpers), so it
## runs in the default suite rather than being gated behind a Stan toolchain.

library(testthat)
library(hydroSense)


## ----------------------------------------------------------------------------
## Finding 5 — PCA predictor cells that are below detection are held at the FULL
## detection limit, with no half-DL substitution (the LMF path uses DL/2). BDL
## nutrients/DOC are PCA variables, so the leading axes are biased upward in
## proportion to censoring. Fix: treat BDL predictor cells at DL/2 (or censored)
## for train/predict parity with the LMF path.
## ----------------------------------------------------------------------------
describe("PCA predictor treatment of below-detection cells (finding 5)", {

  it("scores a BDL-at-DL cell differently from a detected-at-DL cell", {
    train <- .imp_chem(n = 40, seed = 1)
    pca   <- hydroSense:::.prepare_chem_pca(
      train, wq_vars = .imp_pca_vars(), min_var_explained = 0.75, max_pcs = 4L
    )

    ## One query sample; DOC reported exactly at a detection limit of 100.
    q_detected <- .imp_rows("q1", "DOC", 100, detected = TRUE)
    q_bdl      <- .imp_rows("q1", "DOC", 100, detected = FALSE)
    ctx <- dplyr::bind_rows(
      .imp_rows("q1", "pH", 7), .imp_rows("q1", "EC", 400),
      .imp_rows("q1", "Cl", 500)
    )

    sc_detected <- hydroSense:::.compute_pca_scores(
      dplyr::bind_rows(ctx, q_detected), pca)
    sc_bdl <- hydroSense:::.compute_pca_scores(
      dplyr::bind_rows(ctx, q_bdl), pca)

    ## A non-detect at the limit must not score identically to a detection at
    ## the limit: DL/2 makes the BDL value strictly smaller before the log.
    expect_false(isTRUE(all.equal(
      as.numeric(dplyr::select(sc_detected, dplyr::starts_with("PC"))),
      as.numeric(dplyr::select(sc_bdl,      dplyr::starts_with("PC")))
    )))
  })
})


## ----------------------------------------------------------------------------
## Finding 6 — the eps = 1e-9 log floor is a single absolute value shared across
## analytes of wildly different scales. A genuine zero in a high-magnitude
## column (e.g. Cl ~ 1e5) maps to log10(1e-9) = -9, an extreme outlier ~14
## orders below the column's real values, dragging the PCA. Fix: a floor tied to
## the column's own scale (e.g. half its smallest positive value).
## ----------------------------------------------------------------------------
describe("scale-aware log floor for the PCA transform (finding 6)", {

  it("does not map a zero to a column-scale outlier", {
    ## A high-magnitude column with one genuine zero.
    mat <- matrix(c(1e4, 1e5, 5e4, 0), ncol = 1,
                  dimnames = list(NULL, "Cl"))
    out <- hydroSense:::.log_transform_pca(mat)

    transformed_zero <- out[4, 1]
    min_positive     <- min(out[-4, 1])

    ## The floored zero should sit near the column's smallest real value, not
    ## ~13 orders of magnitude below it (log10(1e-9) = -9 vs log10(1e4) = 4).
    expect_lt(min_positive - transformed_zero, 3)
  })
})


## ----------------------------------------------------------------------------
## Finding 7 — the metals presence hurdle is .METAL_ANALYTES, which includes Fe
## and Mn. Those are the routine redox indicators (removed from targets because
## they are PCA predictors), so a sample that measured ONLY Fe/Mn passes the
## trace-metals hurdle and receives fabricated Cd/Hg/Pb/etc. Fix: the hurdle
## should exclude the redox indicators (or be intersected with fitted targets).
## ----------------------------------------------------------------------------
describe("metals presence hurdle excludes redox indicators (finding 7)", {

  it("does not let Fe or Mn alone satisfy the trace-metals hurdle", {
    ## leachate_impute_groups() returns an UNNAMED list keyed by each group's
    ## $name field, so select by name rather than by list name.
    grps <- leachate_impute_groups()
    grp  <- Filter(function(g) g$name == "metals", grps)[[1L]]
    expect_false("Fe" %in% grp$hurdle)
    expect_false("Mn" %in% grp$hurdle)
    ## The hurdle should still gate on genuine trace metals.
    expect_true(all(c("Cd", "Pb", "Zn") %in% grp$hurdle))
  })
})


## ----------------------------------------------------------------------------
## Finding 8a — the PCA wide-pivot collapses duplicate (sample, analyte) rows
## with mean() on the RAW scale, then log-transforms; targets are logged first.
## For a log-normal quantity mean-of-raw != exp(mean-of-log). Fix: collapse
## duplicates on the log (geometric-mean) scale, consistently with targets.
## ----------------------------------------------------------------------------
describe("duplicate collapse is consistent for log-normal analytes (finding 8a)", {

  it("collapses duplicate cells on the log scale (geometric mean)", {
    train <- .imp_chem(n = 40, seed = 2)
    pca   <- hydroSense:::.prepare_chem_pca(
      train, wq_vars = .imp_pca_vars(), min_var_explained = 0.75, max_pcs = 4L
    )

    ctx <- dplyr::bind_rows(
      .imp_rows("q1", "pH", 7), .imp_rows("q1", "EC", 400)
    )
    ## Duplicate Cl rows two orders of magnitude apart: geometric mean = 1000,
    ## arithmetic mean = 5050.
    q_dupe <- dplyr::bind_rows(
      ctx, .imp_rows("q1", "Cl", 100), .imp_rows("q1", "Cl", 10000)
    )
    q_geom <- dplyr::bind_rows(ctx, .imp_rows("q1", "Cl", 1000))

    sc_dupe <- hydroSense:::.compute_pca_scores(q_dupe, pca)
    sc_geom <- hydroSense:::.compute_pca_scores(q_geom, pca)

    ## The duplicated sample must score like its geometric-mean single value.
    expect_equal(
      as.numeric(dplyr::select(sc_dupe, dplyr::starts_with("PC"))),
      as.numeric(dplyr::select(sc_geom, dplyr::starts_with("PC"))),
      tolerance = 1e-6
    )
  })
})


## ----------------------------------------------------------------------------
## Bug B2 — the safe-name inverse lookup names(safe)[safe == safe_nm] is a
## value-equality match. If two analytes collide under make.names() (e.g.
## "Cr-6" and "Cr.6"), it silently returns two originals and mis-maps. Fix: a
## guard that detects the collision and errors early.
## ----------------------------------------------------------------------------
describe("safe-name mapping rejects make.names collisions (bug B2)", {

  it("errors informatively when two analytes share a safe name", {
    ## Target API introduced by the fix: a bijectivity guard used wherever the
    ## engine builds make.names() identifiers for responses.
    expect_error(
      hydroSense:::.assert_safe_analyte_names(c("Cr-6", "Cr.6")),
      regexp = "collid|unique|safe",
      ignore.case = TRUE
    )
    ## A non-colliding set must pass silently.
    expect_no_error(
      hydroSense:::.assert_safe_analyte_names(c("Cr-6", "Zn", "Cu"))
    )
  })
})


## ----------------------------------------------------------------------------
## Bug B3 — mean(x, na.rm = TRUE) over an all-NA (sample, analyte) group returns
## NaN (not NA), which can reach the NIPALS scoring. The zero-variance/all-NA
## guards mostly cover it, so it is latent rather than active. Fix: coerce NaN
## cells so scoring never sees them; scores stay finite.
## ----------------------------------------------------------------------------
describe("all-NA cells never inject NaN into PC scores (bug B3)", {

  it("returns finite scores for a sample whose WQ cell is all-NA", {
    train <- .imp_chem(n = 40, seed = 3)
    pca   <- hydroSense:::.prepare_chem_pca(
      train, wq_vars = .imp_pca_vars(), min_var_explained = 0.75, max_pcs = 4L
    )

    ## Query sample with a present-but-NA DOC row plus real context.
    q <- dplyr::bind_rows(
      .imp_rows("q1", "pH", 7), .imp_rows("q1", "EC", 400),
      .imp_rows("q1", "Cl", 500),
      .imp_rows("q1", "DOC", NA_real_)
    )
    sc <- hydroSense:::.compute_pca_scores(q, pca)
    pc <- as.numeric(dplyr::select(sc, dplyr::starts_with("PC")))
    expect_true(all(is.finite(pc)))
  })
})


## ----------------------------------------------------------------------------
## Findings documented but not encoded as failing specs (no crisp deterministic
## contract; recorded so the fix isn't lost):
##
##   - Finding 4 (co-analyte imputation leakage + covariate shift): DOC/Ca/Mg
##     are predicted from PC axes partly built from themselves, and missing
##     samples use median-filled scores for those very cells. The fix is a
##     design change (hold each co-analyte out of its own scoring axes); its
##     effect is distributional, not a single-value invariant.
##
##   - Finding 8b (co-analyte GAM draws ignore smoothing-parameter uncertainty,
##     so intervals run slightly narrow): a calibration property best checked by
##     coverage on held-out data, not a pointwise assertion.
##
##   - Bug B1 (the documented 2-PC floor is silently violable when only one
##     component exists): behaviour is correct for the only reachable input; the
##     defect is that the docstring overpromises. A documentation fix, so there
##     is no distinguishing runtime contract to encode.
## ----------------------------------------------------------------------------
describe("documented-only findings (4, 8b, B1)", {
  it("finding 4: co-analyte imputation avoids self-prediction leakage")
  it("finding 8b: co-analyte draws include smoothing-parameter uncertainty")
  it("bug B1: n_pcs never contradicts the documented minimum-PC promise")
})
