## Tests for the multi-output daily uncertainty bracket (issue #50).
##
## mspaf_daily() gains a `gap_uncertainty = c("bracket","ignorable",
## "informative")` selector. Of all propagated uncertainty sources only the
## latent Kalman residual is missingness-dependent (it balloons in gaps). The
## INFORMATIVE (lower) envelope freezes that residual at its posterior MEAN on
## in-gap days while reusing every other per-draw source; the IGNORABLE (upper)
## envelope keeps the simulation-smoother draw (today's behaviour). The two are
## nested, coincide at observation days, and diverge only across gaps. A
## deterministic point line and a precautionary composite [lo_informative,
## hi_ignorable] (a decision bound, NOT a calibrated CI) round out the frame.
##
## Functions under test:
##   .summarise_bracket(draws_df, interval, central, gap_uncertainty)  [pure]
##   mspaf_daily(..., gap_uncertainty=)                                [public]
##
## Scientific basis: Rubin (1976) ignorable vs informative missingness; Durbin &
## Koopman (2002) simulation smoother (frozen mean = x_hat, the deviation term
## dropped); Jensen convex multi-stressor ordering (#39/#42).

library(testthat)
library(hydroSense)


## ── Shared fixtures ──────────────────────────────────────────────────────────

make_chem_b <- function(site, dates, mult = 1, seed = 1L) {
  set.seed(seed)
  analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
  purrr::map_dfr(dates, function(d) {
    tibble::tibble(
      sample_id = paste0(site, format(d, "%Y%m%d")),
      site_id = site,
      datetime = d,
      analyte = analytes,
      value = c(
        exp(stats::rnorm(1, log(0.5), 0.3)) * mult,
        exp(stats::rnorm(1, log(5), 0.4)) * mult,
        exp(stats::rnorm(1, log(0.3), 0.3)) * mult,
        stats::runif(1, 6.5, 8), stats::runif(1, 1, 5), stats::runif(1, 20, 60),
        stats::runif(1, 4, 12),  stats::runif(1, 2, 8)
      ),
      detected = TRUE,
      units.analyte = dplyr::case_when(
        analyte %in% c("Cu", "Zn", "Ni") ~ "ug/L",
        TRUE ~ NA_character_
      )
    )
  })
}

make_hydro_b <- function(n = 700L, seed = 99L) {
  set.seed(seed)
  tibble::tibble(
    date  = seq(as.Date("2020-07-01"), by = "day", length.out = n),
    value = pmax(0, stats::rnorm(n, 2, 4))
  )
}

## One reference model fitted once; biweekly grabs leave real 14-day gaps.
.bf <- local({
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 14L)
  hydro <- make_hydro_b()
  ref <- make_chem_b("reference", dates, seed = 1L)
  tgt <- make_chem_b("target", dates, mult = 5, seed = 2L)
  rm <- tryCatch(
    fit_reference_model(ref,
      hydro = hydro, conc_units = "ug/L",
      min_obs_model = 10L,
      api_tau_bounds_short = c(7, 7),
      api_tau_bounds_long = c(30, 30)
    ),
    error = function(e) NULL
  )
  list(rm = rm, tgt = tgt, dates = dates)
})

## Run mspaf_daily with bracket defaults, silencing the rainfall CLI warning.
run_daily <- function(..., gap_uncertainty = "bracket") {
  suppressWarnings(suppressMessages(
    mspaf_daily(
      .bf$tgt,
      reference_model = .bf$rm,
      interpolation = "model",
      require_temperature = FALSE,
      conc_units = "ug/L",
      gap_uncertainty = gap_uncertainty,
      ...
    )
  ))
}

capture_warns <- function(expr) {
  w <- character()
  withCallingHandlers(
    suppressMessages(force(expr)),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    }
  )
  w
}


## ═══════════════════════════════════════════════════════════════════════════
## A. .summarise_bracket() — pure summariser, known-answer oracles
## ═══════════════════════════════════════════════════════════════════════════

describe(".summarise_bracket()", {
  ## Two envelopes of per-draw msPAF for a single (site, date); informative is a
  ## deliberately narrower spread than ignorable.
  toy <- tibble::tibble(
    date = as.Date("2021-02-01"),
    site_id = "S",
    draw_id = 1:5,
    mspaf_ignorable = c(1, 2, 3, 4, 5),
    mspaf_informative = c(2, 2, 3, 4, 4)
  )
  lo_p <- 0.1
  hi_p <- 0.9 # interval = 0.8

  it("bracket mode returns both envelopes + precautionary composite", {
    out <- hydroSense:::.summarise_bracket(
      toy,
      interval = 0.8, central = "median", gap_uncertainty = "bracket"
    )
    expect_true(all(c(
      "date", "site_id",
      "median_informative", "lo_informative", "hi_informative",
      "median_ignorable", "lo_ignorable", "hi_ignorable",
      "precautionary_lo", "precautionary_hi"
    ) %in% names(out)))
    expect_equal(nrow(out), 1L)

    # ignorable envelope: median + (1-interval)/2 quantiles of c(1..5)
    expect_equal(out$median_ignorable, stats::median(1:5))
    expect_equal(
      out$lo_ignorable,
      stats::quantile(c(1, 2, 3, 4, 5), lo_p, names = FALSE)
    )
    expect_equal(
      out$hi_ignorable,
      stats::quantile(c(1, 2, 3, 4, 5), hi_p, names = FALSE)
    )
    # informative envelope on the narrower draws
    expect_equal(out$median_informative, stats::median(c(2, 2, 3, 4, 4)))
    expect_equal(
      out$lo_informative,
      stats::quantile(c(2, 2, 3, 4, 4), lo_p, names = FALSE)
    )
    expect_equal(
      out$hi_informative,
      stats::quantile(c(2, 2, 3, 4, 4), hi_p, names = FALSE)
    )
    # precautionary composite = [lo_informative, hi_ignorable]
    expect_equal(out$precautionary_lo, out$lo_informative)
    expect_equal(out$precautionary_hi, out$hi_ignorable)
  })

  it("ignorable mode returns only the ignorable columns (no informative, no precautionary)", {
    out <- hydroSense:::.summarise_bracket(
      toy,
      interval = 0.8, central = "median", gap_uncertainty = "ignorable"
    )
    expect_true(all(c("median_ignorable", "lo_ignorable", "hi_ignorable") %in%
      names(out)))
    expect_false(any(grepl("informative", names(out))))
    expect_false(any(grepl("precautionary", names(out))))
  })

  it("informative mode returns only the informative columns", {
    out <- hydroSense:::.summarise_bracket(
      toy,
      interval = 0.8, central = "median", gap_uncertainty = "informative"
    )
    expect_true(all(c("median_informative", "lo_informative", "hi_informative") %in%
      names(out)))
    expect_false(any(grepl("ignorable", names(out))))
    expect_false(any(grepl("precautionary", names(out))))
  })

  it("central='mean' uses the per-day draw mean for the centre", {
    out <- hydroSense:::.summarise_bracket(
      toy,
      interval = 0.8, central = "mean", gap_uncertainty = "bracket"
    )
    expect_equal(out$median_ignorable, mean(1:5)) # column keeps its name
    expect_equal(out$median_informative, mean(c(2, 2, 3, 4, 4)))
  })

  it("nested bounds hold within each row (informative ⊆ ignorable)", {
    out <- hydroSense:::.summarise_bracket(
      toy,
      interval = 0.8, central = "median", gap_uncertainty = "bracket"
    )
    expect_lte(out$lo_ignorable, out$lo_informative)
    expect_lte(out$hi_informative, out$hi_ignorable)
  })
})


## ═══════════════════════════════════════════════════════════════════════════
## B. mspaf_daily() bracket integration
## ═══════════════════════════════════════════════════════════════════════════

describe("mspaf_daily(gap_uncertainty=)", {
  ## ── Schema ────────────────────────────────────────────────────────────────

  it("bracket summary returns deterministic + both envelopes + precautionary", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    out <- run_daily(ndraws = 20L, seed = 1L, return = "summary")
    expect_true(all(c(
      "date", "site_id", "deterministic",
      "median_informative", "lo_informative", "hi_informative",
      "median_ignorable", "lo_ignorable", "hi_ignorable",
      "precautionary_lo", "precautionary_hi"
    ) %in% names(out)))
    # the old single-envelope columns are gone (pre-v1, no shim)
    expect_false(any(c("mspaf", "mspaf_lower", "mspaf_upper") %in% names(out)))
    expect_equal(
      nrow(out),
      nrow(dplyr::distinct(dplyr::select(out, "date", "site_id")))
    )
  })

  it("ignorable summary exposes only ignorable envelope columns", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    out <- run_daily(
      ndraws = 20L, seed = 1L, return = "summary",
      gap_uncertainty = "ignorable"
    )
    expect_true(all(c(
      "deterministic", "median_ignorable",
      "lo_ignorable", "hi_ignorable"
    ) %in% names(out)))
    expect_false(any(grepl("informative|precautionary", names(out))))
  })

  it("informative summary exposes only informative envelope columns", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    out <- run_daily(
      ndraws = 20L, seed = 1L, return = "summary",
      gap_uncertainty = "informative"
    )
    expect_true(all(c(
      "deterministic", "median_informative",
      "lo_informative", "hi_informative"
    ) %in% names(out)))
    expect_false(any(grepl("ignorable|precautionary", names(out))))
  })

  it("bracket draws mode carries mspaf_ignorable and mspaf_informative per draw", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    out <- run_daily(ndraws = 6L, seed = 1L, return = "draws")
    expect_true(all(c(
      "date", "site_id", "draw_id",
      "mspaf_ignorable", "mspaf_informative"
    ) %in% names(out)))
    counts <- out |>
      dplyr::group_by(.data$date, .data$site_id) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    expect_true(all(counts$n == 6L))
  })

  ## ── Deterministic line ────────────────────────────────────────────────────

  it("the deterministic column equals the point-mode (ndraws=NULL) estimate", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    pt <- run_daily(ndraws = NULL) |> dplyr::arrange(.data$site_id, .data$date)
    br <- run_daily(ndraws = 20L, seed = 1L, return = "summary") |>
      dplyr::arrange(.data$site_id, .data$date)
    j <- dplyr::inner_join(
      dplyr::select(pt, "date", "site_id", mspaf_pt = "mspaf"),
      dplyr::select(br, "date", "site_id", "deterministic"),
      by = c("date", "site_id")
    )
    expect_equal(j$deterministic, j$mspaf_pt,
      tolerance = 1e-6,
      label = "deterministic == grabs-exact point line"
    )
  })

  ## ── Regression: ignorable reproduces the pre-#50 draw summary numerically ──

  it("the ignorable envelope equals the quantiles/median of its own draws", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    interval <- 0.9
    lo_p <- 0.05
    hi_p <- 0.95
    summ <- run_daily(
      ndraws = 30L, seed = 5L, return = "summary",
      interval = interval, central = "median"
    ) |>
      dplyr::arrange(.data$site_id, .data$date)
    draws <- run_daily(ndraws = 30L, seed = 5L, return = "draws")
    ref <- draws |>
      dplyr::group_by(.data$date, .data$site_id) |>
      dplyr::summarise(
        med = stats::median(.data$mspaf_ignorable),
        lo = stats::quantile(.data$mspaf_ignorable, lo_p, names = FALSE),
        hi = stats::quantile(.data$mspaf_ignorable, hi_p, names = FALSE),
        .groups = "drop"
      ) |>
      dplyr::arrange(.data$site_id, .data$date)
    expect_equal(summ$median_ignorable, ref$med, tolerance = 1e-8)
    expect_equal(summ$lo_ignorable, ref$lo, tolerance = 1e-8)
    expect_equal(summ$hi_ignorable, ref$hi, tolerance = 1e-8)
  })

  ## ── RNG-neutrality: adding the informative path does not perturb ignorable ─

  it("ignorable draws are byte-identical whether or not the informative path runs", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    br <- run_daily(
      ndraws = 8L, seed = 13L, return = "draws",
      gap_uncertainty = "bracket"
    ) |>
      dplyr::arrange(.data$site_id, .data$date, .data$draw_id)
    ig <- run_daily(
      ndraws = 8L, seed = 13L, return = "draws",
      gap_uncertainty = "ignorable"
    ) |>
      dplyr::arrange(.data$site_id, .data$date, .data$draw_id)
    expect_equal(br$mspaf_ignorable, ig$mspaf_ignorable,
      tolerance = 1e-12,
      label = "informative path is RNG-neutral"
    )
  })

  ## ── Nesting and coincidence ───────────────────────────────────────────────

  it("the informative band nests inside the ignorable band", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    out <- run_daily(
      ndraws = 40L, seed = 7L, return = "summary",
      interval = 0.9
    )
    tol <- 1e-6
    # Structural nesting (informative removes one independent noise source);
    # allow a small Monte-Carlo slack per day, require it to hold on the mean.
    expect_gte(mean(out$lo_informative - out$lo_ignorable), -tol)
    expect_lte(mean(out$hi_informative - out$hi_ignorable), tol)
    expect_true(mean(out$lo_ignorable <= out$lo_informative + tol) > 0.9)
    expect_true(mean(out$hi_informative <= out$hi_ignorable + tol) > 0.9)
  })

  it("the two envelopes coincide at observation days and diverge across gaps", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    out <- run_daily(
      ndraws = 40L, seed = 7L, return = "summary",
      interval = 0.9
    )
    is_grab <- out$date %in% .bf$dates
    # between-envelope half-width gap (ignorable upper minus informative upper)
    gap_w <- (out$hi_ignorable - out$hi_informative)
    expect_lt(mean(gap_w[is_grab], na.rm = TRUE),
      mean(gap_w[!is_grab], na.rm = TRUE),
      label = "divergence is larger in gaps than at grabs"
    )
    # at grabs the two medians are close (within at-anchor residual variance)
    med_gap <- abs(out$median_ignorable - out$median_informative)
    expect_lt(
      mean(med_gap[is_grab], na.rm = TRUE),
      mean(med_gap[!is_grab], na.rm = TRUE) + 1e-9
    )
  })

  it("the bracket collapses (informative == ignorable) under dense daily sampling", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    daily_dates <- seq(as.Date("2021-03-01"), by = "day", length.out = 30L)
    tgt_dense <- make_chem_b("dense", daily_dates, mult = 5, seed = 4L)
    out <- suppressWarnings(suppressMessages(mspaf_daily(
      tgt_dense,
      reference_model = .bf$rm, interpolation = "model",
      require_temperature = FALSE, conc_units = "ug/L",
      gap_uncertainty = "bracket", ndraws = 25L, seed = 3L,
      return = "summary", interval = 0.9
    )))
    skip_if(nrow(out) == 0L, "dense run produced no rows")
    # no gaps -> the residual is pinned everywhere -> bands ~coincide
    expect_equal(out$lo_informative, out$lo_ignorable, tolerance = 0.05 *
      (max(out$hi_ignorable) - min(out$lo_ignorable) + 1e-9))
    expect_equal(out$hi_informative, out$hi_ignorable, tolerance = 0.05 *
      (max(out$hi_ignorable) - min(out$lo_ignorable) + 1e-9))
  })

  ## ── Central-tendency ordering (Jensen, #39/#42) ───────────────────────────

  it("ignorable median is the largest central tendency; deterministic ~ informative (Jensen)", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    out <- run_daily(
      ndraws = 60L, seed = 9L, return = "summary",
      central = "median"
    )
    tol <- 1e-6
    # Jensen on the convex multi-stressor combine: more input variance -> higher
    # msPAF central tendency. The IGNORABLE envelope carries the most variance
    # (residual drawn across every gap), so its median is the largest -- robustly
    # above both the informative median and the deterministic line. These two
    # gaps are ~0.1 msPAF and hold across platforms.
    expect_gte(mean(out$median_ignorable), mean(out$median_informative) - tol)
    expect_gte(mean(out$median_ignorable), mean(out$deterministic) - tol)
    # Deterministic and informative nearly COINCIDE: they differ only by the
    # in-gap residual freeze and the per-draw trend perturbation, a second-order
    # Jensen effect that acts on the handful of observed days only (~0.004 msPAF
    # here) and sits below Monte-Carlo noise. So no strict ordering is asserted
    # between them (it inverts across BLAS/RNG platforms); we assert only that
    # they agree to within a small tolerance. Per-day gap structure is covered by
    # the nesting and coincidence tests.
    expect_equal(mean(out$median_informative), mean(out$deterministic),
      tolerance = 0.05
    )
    # The ignorable envelope sits at or above the deterministic line on a clear
    # majority of days.
    expect_gt(mean(out$deterministic <= out$median_ignorable + tol), 0.6)
  })

  ## ── Precautionary composite ───────────────────────────────────────────────

  it("the precautionary composite is exactly [lo_informative, hi_ignorable]", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    out <- run_daily(ndraws = 20L, seed = 1L, return = "summary")
    expect_equal(out$precautionary_lo, out$lo_informative, tolerance = 1e-12)
    expect_equal(out$precautionary_hi, out$hi_ignorable, tolerance = 1e-12)
  })

  ## ── CLI warning ───────────────────────────────────────────────────────────

  it("warns about ignorable gap treatment for rainfall hydrology + uncertainty", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    w <- capture_warns(mspaf_daily(
      .bf$tgt,
      reference_model = .bf$rm, interpolation = "model",
      require_temperature = FALSE, conc_units = "ug/L",
      gap_uncertainty = "bracket", ndraws = 5L, seed = 1L
    ))
    expect_true(any(grepl("ignorable", w, ignore.case = TRUE)))
  })

  it("does not emit the ignorable-gap warning in point mode (no uncertainty)", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    w <- capture_warns(mspaf_daily(
      .bf$tgt,
      reference_model = .bf$rm, interpolation = "model",
      require_temperature = FALSE, conc_units = "ug/L",
      gap_uncertainty = "bracket", ndraws = NULL
    ))
    expect_false(any(grepl("ignorable", w, ignore.case = TRUE)))
  })

  ## ── Edge cases ────────────────────────────────────────────────────────────

  it("point mode (ndraws=NULL) ignores gap_uncertainty and returns the point schema", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    out <- run_daily(ndraws = NULL)
    expect_true(all(c("date", "site_id", "mspaf") %in% names(out)))
    expect_false(any(grepl(
      "informative|ignorable|precautionary|deterministic",
      names(out)
    )))
  })

  it("a single tox analyte still yields a defined informative envelope", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    tgt1 <- dplyr::filter(.bf$tgt, .data$analyte %in%
      c("Cu", "pH", "DOC", "hardness", "Ca", "Mg"))
    out <- suppressWarnings(suppressMessages(mspaf_daily(
      tgt1,
      reference_model = .bf$rm, interpolation = "model",
      require_temperature = FALSE, conc_units = "ug/L", min_analytes = 1L,
      gap_uncertainty = "bracket", ndraws = 10L, seed = 1L,
      return = "summary"
    )))
    skip_if(nrow(out) == 0L, "single-analyte run produced no rows")
    expect_true(all(is.finite(out$median_informative)))
    expect_true(all(is.finite(out$lo_informative)))
    expect_true(all(is.finite(out$hi_informative)))
  })

  it("is reproducible: same seed gives identical bracket draws", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    r1 <- run_daily(ndraws = 6L, seed = 77L, return = "draws")
    r2 <- run_daily(ndraws = 6L, seed = 77L, return = "draws")
    expect_equal(r1$mspaf_ignorable, r2$mspaf_ignorable)
    expect_equal(r1$mspaf_informative, r2$mspaf_informative)
  })
})


## ═══════════════════════════════════════════════════════════════════════════
## D. Graceful degradation when the daily target model cannot be fitted
## ═══════════════════════════════════════════════════════════════════════════
## When .fit_daily_target() returns NULL in draws mode (e.g. every tox analyte
## has fewer than min_obs_model anchors), the site must contribute exact
## forward-fill rows with zero-width intervals — the bracket collapses. The
## fallback previously omitted the `draw_id` carrier, so the draws assembler
## crashed ("Column `draw_id` doesn't exist") instead of degrading. This guards
## that regression.

describe("mspaf_daily() draws-mode graceful degradation (NULL target fit)", {
  ## All metals below detection -> no impact anchors -> no daily impact model
  ## fits -> .fit_daily_target() returns NULL for the site.
  bdl_target <- function() {
    tgt <- make_chem_b("target", .bf$dates, mult = 5, seed = 7L)
    tgt$detected[tgt$analyte %in% c("Cu", "Zn", "Ni")] <- FALSE
    tgt
  }

  it("returns collapsed zero-width bracket instead of erroring", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")

    warns <- character()
    out <- withCallingHandlers(
      suppressMessages(
        mspaf_daily(bdl_target(),
          reference_model = .bf$rm, interpolation = "model",
          require_temperature = FALSE, conc_units = "ug/L",
          gap_uncertainty = "bracket", ndraws = 8L, seed = 1L,
          return = "summary"
        )
      ),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    expect_true(any(grepl("Target model fit failed", warns)))
    expect_s3_class(out, "tbl_df")
    expect_gt(nrow(out), 0L)
    ## Bracket collapsed: both envelopes zero-width and equal to the centre.
    expect_true(all(out$hi_ignorable - out$lo_ignorable < 1e-8, na.rm = TRUE))
    expect_true(all(out$hi_informative - out$lo_informative < 1e-8, na.rm = TRUE))
    expect_equal(out$median_ignorable, out$deterministic, tolerance = 1e-6)
    expect_equal(out$median_informative, out$median_ignorable, tolerance = 1e-6)
  })

  it("carries a draw_id and degrades without erroring in return = 'draws' mode", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")

    out <- suppressWarnings(suppressMessages(
      mspaf_daily(bdl_target(),
        reference_model = .bf$rm, interpolation = "model",
        require_temperature = FALSE, conc_units = "ug/L",
        gap_uncertainty = "bracket", ndraws = 8L, seed = 1L,
        return = "draws"
      )
    ))
    expect_s3_class(out, "tbl_df")
    expect_true("draw_id" %in% names(out))
    expect_gt(nrow(out), 0L)
  })
})


## ═══════════════════════════════════════════════════════════════════════════
## E. Tier-aware ou_scale (#62) + public loo_anchor_coverage()
## ═══════════════════════════════════════════════════════════════════════════

describe("tier-aware ou_scale + loo_anchor_coverage()", {
  it("loo_anchor_coverage() returns a per-analyte + pooled table from a target fit", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    tm <- suppressWarnings(suppressMessages(
      fit_target_model(.bf$tgt, reference_model = .bf$rm, conc_units = "ug/L")
    ))
    cov <- loo_anchor_coverage(tm, interval = 0.9)
    expect_s3_class(cov, "tbl_df")
    expect_true(all(c("analyte", "tier", "n", "coverage", "mean_width") %in%
                      names(cov)))
    expect_true("(pooled)" %in% cov$analyte)
    fin <- cov$coverage[!is.na(cov$coverage)]
    expect_true(all(fin >= 0 & fin <= 1))
  })

  it("a tier-named ou_scale covering all tiers equals the scalar form", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    # The fixture's tox analytes are all model tier, so naming every tier with
    # the same value must reproduce the scalar ou_scale exactly (same seed).
    s_scalar <- run_daily(ndraws = 6L, seed = 5L, ou_scale = 2)
    s_named  <- run_daily(ndraws = 6L, seed = 5L,
                          ou_scale = c(model = 2, bridge = 2))
    expect_equal(s_named$lo_ignorable, s_scalar$lo_ignorable)
    expect_equal(s_named$hi_ignorable, s_scalar$hi_ignorable)
  })

  it("rejects an unnamed multi-element ou_scale", {
    skip_if(is.null(.bf$rm), "Reference model not fitted")
    expect_error(
      suppressWarnings(run_daily(ndraws = 4L, seed = 1L, ou_scale = c(1, 2))),
      "name"
    )
  })
})
