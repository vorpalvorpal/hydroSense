## Tests for .fit_daily_target() + .predict_daily_tox() (issue #16, Chunk C).
##
## These helpers split the formerly-monolithic .daily_tox_from_model() into a
## fit-once scaffold and a per-draw predictor.  Properties tested:
##
##   D1. .fit_daily_target returns a list with required fields
##   D2. fdm$ou factors have n_target == length(fdm$qdates)
##   D3. degenerate OU when anchor count < 2
##   D4. point-mode .predict_daily_tox returns finite C_raw values
##   D5. eps_paths = NULL and eps_paths = zeros give identical output
##   D6. non-zero eps changes C_raw output (variance injected)
##   D7. .daily_tox_from_model wrapper produces expected output

library(testthat)
library(leachatetools)


## ── Shared setup (fitted once at file-parse time) ─────────────────────────────

make_chem <- function(site, dates, mult = 1, seed = 1) {
  set.seed(seed)
  analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
  purrr::map_dfr(dates, function(d) tibble::tibble(
    sample_id = paste0(site, format(d, "%Y%m%d")),
    site_id   = site,
    datetime  = d,
    analyte   = analytes,
    value     = c(
      exp(stats::rnorm(1, log(0.5), 0.3)) * mult,
      exp(stats::rnorm(1, log(5),   0.4)) * mult,
      exp(stats::rnorm(1, log(0.3), 0.3)) * mult,
      stats::runif(1, 6.5, 8), stats::runif(1, 1, 5), stats::runif(1, 20, 60),
      stats::runif(1, 4, 12),  stats::runif(1, 2, 8)
    ),
    detected = TRUE
  ))
}

make_hydro <- function(n = 700, seed = 99) {
  set.seed(seed)
  tibble::tibble(
    date  = seq(as.Date("2020-07-01"), by = "day", length.out = n),
    value = pmax(0, stats::rnorm(n, 2, 4))
  )
}

## One-time expensive setup: fits reference model on synthetic data
.td <- local({
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
  hydro <- make_hydro()
  ref   <- make_chem("reference", dates, seed = 1)
  tgt   <- make_chem("target",    dates, mult = 5, seed = 2)
  rm    <- fit_reference_model(ref, hydro = hydro, conc_units = "ug/L",
                                min_obs_model = 10L,
                                api_windows_short = 7L,
                                api_windows_long  = 30L)
  all_dates <- seq(dates[1], dates[length(dates)], by = "day")
  meta      <- leachatetools:::.load_analyte_metadata(NULL)
  tox_nms   <- meta$analyte[!is.na(meta$ssd_available) & meta$ssd_available]
  daily_long <- leachatetools:::.build_daily_chem(
    site_rows     = tgt,
    dates         = all_dates,
    interpolation = "forward_fill",
    leading_edge  = "drop",
    tox_analytes  = tox_nms
  )
  list(rm = rm, tgt = tgt, all_dates = all_dates, daily_long = daily_long,
       meta = meta, tox_nms = tox_nms)
})


## ── D1: .fit_daily_target structure ──────────────────────────────────────────

test_that("D1: .fit_daily_target returns a list with required fields", {
  fdm <- leachatetools:::.fit_daily_target(
    site_rows        = .td$tgt,
    reference_model  = .td$rm,
    imputation_model = NULL,
    conc_units       = "ug/L",
    meta             = .td$meta,
    tox_analytes     = .td$tox_nms,
    daily_long       = .td$daily_long
  )
  expect_type(fdm, "list")
  expect_true(all(c("tm", "modelled", "qdates", "co_split",
                    "wq_long", "fac_lookup", "measured_key", "ou") %in% names(fdm)))
  expect_s3_class(fdm$tm, "target_model")
  expect_true(length(fdm$modelled) > 0L)
  expect_true(all(fdm$modelled %in% names(fdm$tm$models)))
  expect_true(length(fdm$qdates) == length(.td$all_dates))
})


## ── D2: OU factors have correct length ───────────────────────────────────────

test_that("D2: fdm$ou[[nm]]$factors$n_target == length(fdm$qdates)", {
  fdm <- leachatetools:::.fit_daily_target(
    site_rows        = .td$tgt,
    reference_model  = .td$rm,
    imputation_model = NULL,
    conc_units       = "ug/L",
    meta             = .td$meta,
    tox_analytes     = .td$tox_nms,
    daily_long       = .td$daily_long
  )
  for (nm in fdm$modelled) {
    expect_equal(fdm$ou[[nm]]$factors$n_target, length(fdm$qdates),
                 info = paste("analyte", nm))
  }
})


## ── D3: degenerate OU when anchors < 2 ───────────────────────────────────────

test_that("D3: degenerate OU params when analyte has < 2 anchor observations", {
  ## Build minimal single-obs target chemistry (1 grab per analyte)
  single_date <- as.Date("2021-07-01")
  one_grab    <- make_chem("target", single_date, mult = 3, seed = 5)
  daily_1     <- leachatetools:::.build_daily_chem(
    site_rows     = one_grab,
    dates         = seq(single_date, single_date + 30, by = "day"),
    interpolation = "forward_fill",
    leading_edge  = "drop",
    tox_analytes  = .td$tox_nms
  )
  fdm <- leachatetools:::.fit_daily_target(
    site_rows        = one_grab,
    reference_model  = .td$rm,
    imputation_model = NULL,
    conc_units       = "ug/L",
    meta             = .td$meta,
    tox_analytes     = .td$tox_nms,
    daily_long       = daily_1
  )
  if (!is.null(fdm)) {
    for (nm in fdm$modelled) {
      expect_true(fdm$ou[[nm]]$params$degenerate,
                  info = paste("expected degenerate OU for", nm, "with 1 grab"))
    }
  }
})


## ── D4: point-mode prediction returns finite values ──────────────────────────

test_that("D4: .predict_daily_tox (point mode) returns finite C_raw for all rows", {
  fdm <- leachatetools:::.fit_daily_target(
    site_rows        = .td$tgt,
    reference_model  = .td$rm,
    imputation_model = NULL,
    conc_units       = "ug/L",
    meta             = .td$meta,
    tox_analytes     = .td$tox_nms,
    daily_long       = .td$daily_long
  )
  model_rows <- leachatetools:::.predict_daily_tox(fdm)
  expect_true(!is.null(model_rows))
  expect_true(nrow(model_rows) > 0L)
  expect_true(all(is.finite(model_rows$value)))
  expect_true(all(model_rows$value >= 0))
  expect_true(all(fdm$modelled %in% model_rows$analyte))
})


## ── D5: eps_paths = NULL ≡ eps_paths = all-zeros ─────────────────────────────

test_that("D5: NULL eps_paths and zero eps_paths give identical output", {
  fdm <- leachatetools:::.fit_daily_target(
    site_rows        = .td$tgt,
    reference_model  = .td$rm,
    imputation_model = NULL,
    conc_units       = "ug/L",
    meta             = .td$meta,
    tox_analytes     = .td$tox_nms,
    daily_long       = .td$daily_long
  )
  rows_null <- leachatetools:::.predict_daily_tox(fdm, eps_paths = NULL)
  zero_eps  <- stats::setNames(
    lapply(fdm$modelled, function(nm) rep(0, length(fdm$qdates))),
    fdm$modelled
  )
  rows_zero <- leachatetools:::.predict_daily_tox(fdm, eps_paths = zero_eps)
  expect_equal(rows_null$value, rows_zero$value)
})


## ── D6: non-zero eps changes output ──────────────────────────────────────────

test_that("D6: positive eps_paths shifts C_raw upwards for at least one analyte", {
  fdm <- leachatetools:::.fit_daily_target(
    site_rows        = .td$tgt,
    reference_model  = .td$rm,
    imputation_model = NULL,
    conc_units       = "ug/L",
    meta             = .td$meta,
    tox_analytes     = .td$tox_nms,
    daily_long       = .td$daily_long
  )
  rows_base <- leachatetools:::.predict_daily_tox(fdm, eps_paths = NULL)

  ## Inject a large positive constant ε on the first modelled analyte
  nm        <- fdm$modelled[[1L]]
  big_eps   <- stats::setNames(
    list(rep(1000, length(fdm$qdates))),
    nm
  )
  rows_eps  <- leachatetools:::.predict_daily_tox(fdm, eps_paths = big_eps)

  base_nm <- rows_base$value[rows_base$analyte == nm]
  eps_nm  <- rows_eps$value[rows_eps$analyte == nm]
  expect_true(all(eps_nm >= base_nm - .Machine$double.eps),
              label = "positive eps should shift value >= base")
  expect_true(any(eps_nm > base_nm),
              label = "at least some values strictly higher")
})


## ── D7: .daily_tox_from_model wrapper end-to-end ─────────────────────────────

test_that("D7: .daily_tox_from_model replaces modelled-analyte rows and attaches impact_tiers", {
  out <- leachatetools:::.daily_tox_from_model(
    daily_long       = .td$daily_long,
    site_rows        = .td$tgt,
    reference_model  = .td$rm,
    imputation_model = NULL,
    conc_units       = "ug/L",
    meta             = .td$meta,
    tox_analytes     = .td$tox_nms
  )
  ## Still a data frame with the expected columns
  expect_true(is.data.frame(out))
  expect_true(all(c(".date", "value", "analyte", "detected") %in% names(out)))
  ## impact_tiers attribute is attached
  tiers <- attr(out, "impact_tiers")
  expect_true(!is.null(tiers))
  expect_true(all(c("analyte", "impact_tier") %in% names(tiers)))
  expect_true(all(tiers$impact_tier %in% c("model", "bridge", "wq")))
  ## All values finite
  tox_rows <- out[out$analyte %in% .td$tox_nms, , drop = FALSE]
  expect_true(all(is.finite(tox_rows$value)))
})
