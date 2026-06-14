## Behaviour specification for the variance-stabilising transform of the daily
## impact residual (issue #15 <U+2014> fixes the baseline over-dispersion of #39).
##
## Plan (issue #15 comments): smooth the ARA impact I = C_norm - ref_norm in a
## variance-stabilising space g = asinh(I / c), with per-analyte scale c = HC5
## (the SSD 5% hazard concentration). ref stays ONLY in the ARA difference; the
## transform involves c (from the SSD) alone. asinh (not log(I+c)) because the
## impact is signed (C < ref is real precipitation).
##
## New code (does not exist yet <U+2014> every it() starts with skip(), so the suite is
## PENDING until the implement skill lands it):
##   R/transform.R (or R/target_model.R):
##     .g_transform(I, c)   -> asinh(I / c)
##     .g_inverse(g, c)     -> c * sinh(g)
##     .analyte_c(fit)      -> ssdtools::ssd_hc(fit, proportion = 0.05)$est  (HC5)
##   Wiring (R/target_model.R / R/amspaf_daily.R): anchors and the impact-GAM
##   residual are built in g-space; reconstruction inverse-transforms back to I.
##
## Units: HC5 from the fitted SSD is in the SAME normalised concentration space
## as C_norm / ref_norm (normalisation maps measured C onto the SSD scale), so
## c and I share units and I/c is dimensionless.
##
## NOTE: the real-data baseline-tightening acceptance (NH3-N PAF q99 collapses
## on B.S01 without breaking event coverage) is a Stage-4 dev validation
## (dev/jan2024_investigation.R, dev/loo_coverage_bs01.R), not a unit test. The
## wiring specs below assert the mathematical behaviour that validation relies on.

library(testthat)
library(leachatetools)

PENDING <- "pending: #15 <U+2014> asinh variance-stabilising transform"

## <U+2500><U+2500> Pure transform helpers <U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500>

describe(".g_transform() / .g_inverse()", {
  it("round-trips: g_inverse(g_transform(I, c), c) == I for signed I", {
    c <- 12.5
    I <- c(-5000, -100, -1, -1e-6, 0, 1e-6, 1, 100, 5000)
    expect_equal(
      leachatetools:::.g_inverse(leachatetools:::.g_transform(I, c), c),
      I,
      tolerance = 1e-10
    )
  })

  it("maps zero impact to zero in both directions", {
    expect_identical(leachatetools:::.g_transform(0, 7), 0)
    expect_identical(leachatetools:::.g_inverse(0, 7), 0)
  })

  it("is sign-preserving and strictly increasing in I", {
    c <- 3
    I <- sort(c(-1000, -10, -0.1, 0, 0.1, 10, 1000))
    g <- leachatetools:::.g_transform(I, c)
    expect_true(all(diff(g) > 0)) # strictly increasing
    expect_identical(sign(g), sign(I)) # sign preserved
  })

  it("is additive for |I| << c: g ~= I/c (slope 1/c at the origin)", {
    # asinh(x) -> x as x -> 0, so g_transform(I, c) -> I/c for |I| << c.
    c <- 50
    I <- c * c(-1e-3, -1e-4, 1e-4, 1e-3)
    expect_equal(leachatetools:::.g_transform(I, c), I / c, tolerance = 1e-6)
  })

  it("is logarithmic for |I| >> c: g ~= sign(I)*log(2|I|/c)", {
    # asinh(x) ~ sign(x)*log(2|x|) for |x| >> 1, so g ~ sign(I)*log(2|I|/c).
    c <- 4
    I <- c(-1e7, 1e7)
    expect_equal(leachatetools:::.g_transform(I, c),
      sign(I) * log(2 * abs(I) / c),
      tolerance = 1e-6
    )
  })

  it("propagates NA and maps +/-Inf to +/-Inf", {
    expect_identical(leachatetools:::.g_transform(NA_real_, 5), NA_real_)
    expect_identical(leachatetools:::.g_inverse(NA_real_, 5), NA_real_)
    expect_identical(leachatetools:::.g_transform(c(-Inf, Inf), 5), c(-Inf, Inf))
  })

  it("errors on a non-positive scale c", {
    expect_snapshot(leachatetools:::.g_transform(1, 0), error = TRUE)
    expect_snapshot(leachatetools:::.g_transform(1, -3), error = TRUE)
  })

  it("is vectorised over I with a scalar c", {
    c <- 9
    I <- runif(100, -200, 200)
    g <- leachatetools:::.g_transform(I, c)
    expect_length(g, 100L)
    expect_equal(leachatetools:::.g_inverse(g, c), I, tolerance = 1e-10)
  })
})

## <U+2500><U+2500> Per-analyte scale c = HC5 <U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500>

describe(".analyte_c()", {
  ## A real fitted SSD for the oracle (multi method, no guideline_dir).
  get_fit <- function(analyte = "Cu") {
    meta <- leachatetools:::.load_analyte_metadata(NULL)
    sp <- suppressMessages(
      leachatetools:::derive_ssd_params(meta,
        method = "multi",
        guideline_dir = NULL
      )
    )
    sp$fit[[which(sp$analyte == analyte)]]
  }

  it("returns the SSD 5% hazard concentration (HC5)", {
    fit <- get_fit("Cu")
    # Oracle: same call the package uses for HC5 (R/paf.R).
    hc5 <- ssdtools::ssd_hc(fit, proportion = 0.05, ci = FALSE)$est
    expect_equal(leachatetools:::.analyte_c(fit), hc5, tolerance = 1e-8)
  })

  it("returns a single finite positive scale for a normal fit", {
    cc <- leachatetools:::.analyte_c(get_fit("Zn"))
    expect_length(cc, 1L)
    expect_true(is.finite(cc) && cc > 0)
  })

  it("errors or returns NA (never silently invalid) for a NULL / unusable fit", {
    # Exact degraded behaviour (error vs NA) is the implementer's per-plan
    # choice; the invariant is that it never returns a non-positive scale that
    # would break asinh(I / c).
    out <- tryCatch(leachatetools:::.analyte_c(NULL),
      error = function(e) NA_real_
    )
    expect_true(is.na(out) || (is.finite(out) && out > 0))
  })
})

## <U+2500><U+2500> Transform wiring into the smoother (the math the target-model uses) <U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500><U+2500>
## These exercise the existing residual smoother fed transformed anchors, which
## is exactly what fit_target_model() will do internally once wired. They
## specify: (1) anchor round-trip exactness (centre unchanged at grabs),
## (2) geometric mid-gap interpolation, (3) baseline draw tightening because
## gamma is re-estimated in the compressed g-space.

describe("asinh transform wiring (daily impact smoother)", {
  it("reproduces the measured impact at grab anchors (centre unchanged at grabs)", {
    dates <- as.Date("2021-01-01") + (seq_len(12) - 1L) * 30L
    I <- c(1, 2, 50, 3, 100, 2, 1, 40, 2, 3, 80, 1)
    cc <- 20
    tdates <- seq(min(dates), max(dates), by = "day")
    g <- leachatetools:::.g_transform(I, cc)
    sm <- leachatetools:::.residual_smoother(dates, g, tdates)
    gi <- match(dates, sm$grid_dates)
    I_hat <- leachatetools:::.g_inverse(sm$mean[gi], cc)
    # smoother pins ~exactly at anchors (tiny anchor obs-noise), so the
    # back-transformed deterministic centre recovers the measured impact.
    expect_equal(I_hat, I, tolerance = 1e-2)
  })

  it("interpolates geometrically between anchors, well below the linear value", {
    dates <- as.Date(c("2021-01-01", "2021-03-02")) # 60-day gap
    I <- c(1, 1000)
    cc <- 10
    tdates <- seq(min(dates), max(dates), by = "day")
    g <- leachatetools:::.g_transform(I, cc)
    sm <- leachatetools:::.residual_smoother(dates, g, tdates)
    mid <- as.Date("2021-01-31")
    gi <- match(mid, sm$grid_dates)
    I_mid <- leachatetools:::.g_inverse(sm$mean[gi], cc)
    linear_mid <- stats::approx(as.numeric(dates), I,
      xout = as.numeric(mid)
    )$y # ~500
    # geometric (g-space) interpolation pulls the mid-gap impact far below the
    # arithmetic midpoint.
    expect_lt(I_mid, linear_mid / 2)
  })

  it("bounds the baseline draw spread at the c-scale, unlike the additive smoother", {
    withr::local_seed(42)
    dates <- as.Date("2021-01-01") + (seq_len(24) - 1L) * 30L
    I <- rep(2, length(dates))
    I[c(10, 11)] <- 6000 # baseline + 1 event
    cc <- 100
    tdates <- seq(min(dates), max(dates), by = "day")
    base_day <- as.Date("2021-02-15") # baseline MID-GAP day (between anchors)

    # additive smoother: global gamma is inflated by the event spike
    sm_add <- leachatetools:::.residual_smoother(dates, I, tdates)
    dr_add <- leachatetools:::.kalman_draw(sm_add$model, 500L)
    add_q99 <- stats::quantile(
      dr_add[match(base_day, sm_add$grid_dates), ], 0.99,
      names = FALSE
    )

    # g-space smoother: gamma re-estimated in the compressed space
    g <- leachatetools:::.g_transform(I, cc)
    sm_g <- leachatetools:::.residual_smoother(dates, g, tdates)
    dr_g <- leachatetools:::.g_inverse(
      leachatetools:::.kalman_draw(sm_g$model, 500L), cc
    )
    g_q99 <- stats::quantile(
      dr_g[match(base_day, sm_g$grid_dates), ], 0.99,
      names = FALSE
    )

    # the transform dramatically tightens the baseline upper tail ...
    expect_lt(g_q99, add_q99 / 5)
    # ... and keeps it near the c-scale rather than the event scale (6000).
    expect_lt(g_q99, 10 * cc)
  })
})

## ── Stage 2 integration: fit_target_model() uses the transform ────────────────
## End-to-end wiring through the real fit_target_model() + .resolve_target_impact()
## (heavy synthetic fits, so each it() builds its own model AFTER skip()).

describe("fit_target_model() asinh wiring (issue #15)", {
  make_chem <- function(site, dates, mult = 1, seed = 1) {
    set.seed(seed)
    analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
    purrr::map_dfr(dates, function(d) {
      tibble::tibble(
        sample_id = paste0(site, format(d, "%Y%m%d")), site_id = site,
        datetime = d, analyte = analytes,
        value = c(
          exp(stats::rnorm(1, log(0.5), 0.3)) * mult,
          exp(stats::rnorm(1, log(5), 0.4)) * mult,
          exp(stats::rnorm(1, log(0.3), 0.3)) * mult,
          stats::runif(1, 6.5, 8), stats::runif(1, 1, 5), stats::runif(1, 20, 60),
          stats::runif(1, 4, 12), stats::runif(1, 2, 8)
        ),
        detected = TRUE
      )
    })
  }
  make_tm <- function() {
    dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
    hydro <- tibble::tibble(
      date = seq(as.Date("2020-07-01"), by = "day", length.out = 700),
      value = pmax(0, stats::rnorm(700, 2, 4))
    )
    ref <- make_chem("reference", dates, seed = 1)
    tgt <- make_chem("target", dates, mult = 5, seed = 2)
    rm <- suppressMessages(fit_reference_model(
      ref, hydro = hydro, conc_units = "ug/L", min_obs_model = 10L,
      api_windows_short = 7L, api_windows_long = 30L
    ))
    suppressMessages(fit_target_model(
      tgt, rm, conc_units = "ug/L", min_obs_model = 10L,
      api_windows_short = 7L, api_windows_long = 30L
    ))
  }

  it("stores a finite positive per-analyte HC5 transform scale on each model", {
    tm <- make_tm()
    sc <- vapply(tm$models, function(m) m$scale_c %||% NA_real_, numeric(1L))
    expect_true(all(is.finite(sc) & sc > 0))
    # Cu's stored scale matches the SSD HC5 oracle.
    cu_fit <- {
      meta <- leachatetools:::.load_analyte_metadata(NULL)
      sp <- suppressMessages(leachatetools:::derive_ssd_params(
        meta, method = "multi", guideline_dir = NULL))
      sp$fit[[which(sp$analyte == "Cu")]]
    }
    expect_equal(unname(sc["Cu"]),
                 ssdtools::ssd_hc(cu_fit, proportion = 0.05, ci = FALSE)$est,
                 tolerance = 1e-6)
  })

  it("reconstructs the measured impact at grab anchors (end-to-end anchor-exact)", {
    tm <- make_tm()
    nm <- "Zn"
    anch <- tm$models[[nm]]$anchors
    res <- leachatetools:::.resolve_target_impact(
      tm, tibble::tibble(date = anch$date), analytes = nm)
    j <- dplyr::inner_join(
      dplyr::select(anch, date, I),
      dplyr::select(res, date, impact), by = "date")
    # transform round-trip is anchor-exact: reconstructed impact == measured I
    # (to the smoother's tiny anchor obs-noise, ~1e-5 rel).
    expect_equal(j$impact, j$I, tolerance = 1e-3)
  })
})

## ── Stage 3: S6 grab measurement error mapped into g-space (delta method) ──────

describe(".s6_var_to_g() (issue #15)", {
  it("maps observation variance by the squared transform slope", {
    # delta method: Var_g = Var_I * g'(I)^2 = Var_I / (I^2 + c^2)
    var_i <- 4; impact <- 30; cc <- 10
    expect_equal(leachatetools:::.s6_var_to_g(var_i, impact, cc),
                 var_i / (impact^2 + cc^2), tolerance = 1e-12)
  })

  it("turns multiplicative grab error into ~constant g-space noise at high impact", {
    # multiplicative measurement error var_I = (cv*I)^2 -> Var_g -> cv^2 as I>>c,
    # and shrinks toward 0 as I->0 (so baseline anchors are not over-noised).
    cv <- 0.15; cc <- 10
    vg_hi <- leachatetools:::.s6_var_to_g((cv * 1e4)^2, 1e4, cc)
    vg_lo <- leachatetools:::.s6_var_to_g((cv * 0.1)^2, 0.1, cc)
    expect_equal(vg_hi, cv^2, tolerance = 1e-3)   # plateau at high impact
    expect_lt(vg_lo, cv^2 / 100)                  # vanishes at baseline
  })
})

## ── transform parameter: pseudo_log vs additive ──────────────────────────────
## Reuses the make_tm() helper defined in the "fit_target_model() asinh wiring"
## describe block above.  Here we rebuild it inline to keep each it() self-
## contained (the fit is ~5 s on synthetic data).

describe("fit_target_model() transform parameter", {
  make_chem2 <- function(site, dates, mult = 1, seed = 1) {
    set.seed(seed)
    analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
    purrr::map_dfr(dates, function(d) {
      tibble::tibble(
        sample_id = paste0(site, format(d, "%Y%m%d")), site_id = site,
        datetime = d, analyte = analytes,
        value = c(
          exp(stats::rnorm(1, log(0.5), 0.3)) * mult,
          exp(stats::rnorm(1, log(5), 0.4)) * mult,
          exp(stats::rnorm(1, log(0.3), 0.3)) * mult,
          stats::runif(1, 6.5, 8), stats::runif(1, 1, 5), stats::runif(1, 20, 60),
          stats::runif(1, 4, 12), stats::runif(1, 2, 8)
        ),
        detected = TRUE
      )
    })
  }
  make_tm2 <- function(transform = "pseudo_log") {
    dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
    hydro <- tibble::tibble(
      date  = seq(as.Date("2020-07-01"), by = "day", length.out = 700),
      value = pmax(0, stats::rnorm(700, 2, 4))
    )
    ref <- make_chem2("reference", dates, seed = 1)
    tgt <- make_chem2("target", dates, mult = 5, seed = 2)
    rm <- suppressMessages(fit_reference_model(
      ref, hydro = hydro, conc_units = "ug/L", min_obs_model = 10L,
      api_windows_short = 7L, api_windows_long = 30L
    ))
    suppressMessages(fit_target_model(
      tgt, rm, conc_units = "ug/L", min_obs_model = 10L,
      api_windows_short = 7L, api_windows_long = 30L,
      transform = transform
    ))
  }

  it("pseudo_log (default) stores finite positive scale_c on every model", {
    tm <- make_tm2("pseudo_log")
    sc <- vapply(tm$models, function(m) m$scale_c %||% NA_real_, numeric(1L))
    expect_true(all(is.finite(sc) & sc > 0))
  })

  it("additive sets scale_c = NA on every model (pre-#15 additive path)", {
    tm <- make_tm2("additive")
    sc <- vapply(tm$models, function(m) m$scale_c %||% NA_real_, numeric(1L))
    expect_true(all(is.na(sc)))
  })

  it("both transforms produce finite anchor reconstructions (end-to-end)", {
    tm_pl <- make_tm2("pseudo_log")
    tm_ad <- make_tm2("additive")
    nm <- intersect(names(tm_pl$models), names(tm_ad$models))[[1L]]
    anch_pl <- tm_pl$models[[nm]]$anchors
    anch_ad <- tm_ad$models[[nm]]$anchors
    res_pl <- leachatetools:::.resolve_target_impact(
      tm_pl, tibble::tibble(date = anch_pl$date), analytes = nm)
    res_ad <- leachatetools:::.resolve_target_impact(
      tm_ad, tibble::tibble(date = anch_ad$date), analytes = nm)
    expect_true(all(is.finite(res_pl$impact)))
    expect_true(all(is.finite(res_ad$impact)))
  })

  it("the two transforms produce different centre-line predictions", {
    # With high-impact synthetic data, pseudo_log compresses the g-space whereas
    # additive does not.  The predicted impact should differ (not identical).
    tm_pl <- make_tm2("pseudo_log")
    tm_ad <- make_tm2("additive")
    nm <- intersect(names(tm_pl$models), names(tm_ad$models))[[1L]]
    dates <- seq(as.Date("2021-01-15"), by = "30 days", length.out = 5)
    res_pl <- leachatetools:::.resolve_target_impact(
      tm_pl, tibble::tibble(date = dates), analytes = nm)
    res_ad <- leachatetools:::.resolve_target_impact(
      tm_ad, tibble::tibble(date = dates), analytes = nm)
    # They should differ for at least some interpolated dates
    expect_false(isTRUE(all.equal(res_pl$impact, res_ad$impact)))
  })
})
