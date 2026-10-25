## ============================================================================
## Synthetic fixtures for the msPAF performance-refactor characterisation tests
## ============================================================================
##
## Used by test-perf-goldens.R, test-perf-structure.R and the golden generator
## data-raw/make_perf_goldens.R.  Everything here is SYNTHETIC: site codes
## ("P1".."P3", "R1"), dates, chemistry, rainfall and temperature are drawn
## from fixed seeds.  No real monitoring data is used.
##
## Design (chosen so every branch of the three hot functions is exercised):
##   * three target points with different sampling regimes:
##       P1  dense, irregular, with storm clusters, two samples sharing a
##           date (zero-length forward-step interval), BDL Cu/Ni;
##       P2  regular ~45-day sampling that starts later (leading edge);
##       P3  sparse, with a > 365-day gap (focal dates with no in-window
##           samples; windows that only see the anchor).
##   * a reference point R1 sampled ~monthly: Cu and Zn are rain-driven
##     (tier-2 "model" GAMs), Ni is sparse (tier-3 "static"); sampling dates
##     coincide with some target dates (tier-1 "matched" candidates).
##   * a daily rainfall series with NA values AND missing date ranges, which
##     starts after some of the API query dates.
##   * a daily water-temperature series with a gap.

.pf_dates_cluster <- function(seed) {
  set.seed(seed)
  base <- seq(as.Date("2021-01-12"), as.Date("2023-12-10"), by = "31 days")
  jit  <- base + sample(-12:12, length(base), replace = TRUE)
  storms <- as.Date(c("2021-06-03", "2021-06-05", "2021-06-09",
                      "2022-02-20", "2022-02-23",
                      "2023-03-01", "2023-03-04"))
  dup <- as.Date("2022-08-15") # sampled twice on the same day (see below)
  sort(c(jit, storms, dup, dup))
}

## Long-format chemistry for one point.  `mult` scales the toxicants (impact).
.pf_point_chem <- function(site, dates, mult, seed, ni_every = 2L,
                           bdl_frac = 0.15, rain = NULL) {
  set.seed(seed)
  n  <- length(dates)
  # Rain-driven signal: mean rainfall over the 7 days before each sample.
  rain_sig <- if (is.null(rain)) rep(0, n) else vapply(dates, function(d) {
    r <- rain$value[rain$date > d - 7 & rain$date <= d]
    if (length(r) == 0L) 0 else mean(r, na.rm = TRUE)
  }, numeric(1))
  rain_sig[!is.finite(rain_sig)] <- 0
  doy <- as.integer(format(dates, "%j"))
  sid <- sprintf("%s_%03d", site, seq_len(n))

  rows <- list(
    Cu = mult * exp(log(1.2) + 0.08 * rain_sig + stats::rnorm(n, 0, 0.25)),
    Zn = mult * exp(log(8) + 0.06 * rain_sig + stats::rnorm(n, 0, 0.30)),
    Ni = mult * exp(log(1.5) + stats::rnorm(n, 0, 0.30)),
    `NH3-N` = mult * exp(log(250) + stats::rnorm(n, 0, 0.5)),
    pH  = stats::runif(n, 6.6, 7.8),
    EC  = stats::runif(n, 150, 450),
    DOC = stats::runif(n, 1, 5),
    hardness = stats::runif(n, 25, 90),
    Ca  = stats::runif(n, 4, 14),
    Mg  = stats::runif(n, 2, 9),
    temperature = 15 + 6 * sin(2 * pi * (doy - 15) / 365.25) +
      stats::rnorm(n, 0, 0.7)
  )
  out <- dplyr::bind_rows(lapply(names(rows), function(a) {
    tibble::tibble(
      sample_id = sid, site_id = site, datetime = dates, analyte = a,
      value = rows[[a]], detected = TRUE
    )
  }))

  # Ni measured only on every `ni_every`-th sample (ragged analyte coverage).
  ni_keep <- seq_len(n) %% ni_every == 1L
  out <- out[!(out$analyte == "Ni" & !out$sample_id %in% sid[ni_keep]), ]
  # NH3-N missing on a few samples.
  out <- out[!(out$analyte == "NH3-N" & out$sample_id %in% sid[seq(3, n, 7)]), ]

  # Below-detection Cu / Ni: value is the detection limit, detected = FALSE.
  for (a in c("Cu", "Ni")) {
    idx <- which(out$analyte == a)
    bdl <- idx[stats::runif(length(idx)) < bdl_frac]
    out$value[bdl]    <- if (a == "Cu") 0.5 else 0.8
    out$detected[bdl] <- FALSE
  }

  out$units.analyte <- ifelse(out$analyte %in% c("Cu", "Zn", "Ni", "NH3-N"),
                              "ug/L", NA_character_)
  out[order(out$datetime, out$sample_id, out$analyte, method = "radix"), ]
}

## Daily rainfall with NA values and two missing date ranges.
.pf_rain <- function(seed = 11L) {
  set.seed(seed)
  d <- seq(as.Date("2020-09-01"), as.Date("2023-12-31"), by = "day")
  wet <- stats::runif(length(d)) < 0.35
  v <- ifelse(wet, stats::rgamma(length(d), shape = 0.8, scale = 9), 0)
  v[stats::runif(length(d)) < 0.04] <- NA_real_ # NA rain days (0 input)
  rain <- tibble::tibble(date = d, value = v)
  # Missing rows (record gaps): no row at all for these dates.
  gap1 <- d >= as.Date("2021-08-10") & d <= as.Date("2021-08-25")
  gap2 <- d >= as.Date("2022-11-01") & d <= as.Date("2022-11-04")
  rain[!(gap1 | gap2), ]
}

## Daily water temperature (datetime, value) with a 3-week gap.
.pf_temperature <- function(seed = 12L) {
  set.seed(seed)
  d <- seq(as.Date("2021-01-01"), as.Date("2023-12-31"), by = "day")
  doy <- as.integer(format(d, "%j"))
  v <- 15 + 6 * sin(2 * pi * (doy - 15) / 365.25) + stats::rnorm(length(d), 0, 0.4)
  keep <- !(d >= as.Date("2022-06-01") & d <= as.Date("2022-06-21"))
  tibble::tibble(datetime = d[keep], value = v[keep])
}

## The full synthetic scenario.  Memoised per R session (cheap anyway).
.pf_scenario <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    rain <- .pf_rain()
    p1_dates <- .pf_dates_cluster(21L)
    p2_dates <- seq(as.Date("2021-04-20"), as.Date("2023-11-30"), by = "45 days")
    p3_dates <- as.Date(c("2021-02-02", "2021-03-30", "2021-06-15",
                          "2021-09-01", "2021-10-20",
                          # > 365-day gap
                          "2022-12-05", "2023-02-14", "2023-05-02",
                          "2023-07-19", "2023-10-03", "2023-12-12"))
    r1_dates <- sort(unique(c(
      seq(as.Date("2020-10-05"), as.Date("2023-12-20"), by = "26 days"),
      p1_dates[c(5, 12, 20)] # same-day reference grabs -> tier-1 matches
    )))

    target <- dplyr::bind_rows(
      .pf_point_chem("P1", p1_dates, mult = 4, seed = 31L, rain = rain),
      .pf_point_chem("P2", p2_dates, mult = 2, seed = 32L, rain = rain),
      .pf_point_chem("P3", p3_dates, mult = 6, seed = 33L, ni_every = 1L,
                     rain = rain)
    )
    reference <- .pf_point_chem("R1", r1_dates, mult = 1, seed = 41L,
                                ni_every = 4L, rain = rain)

    cache <<- list(
      target      = target,
      reference   = reference,
      rain        = rain,
      temperature = .pf_temperature()
    )
    cache
  }
})

## Reference model fitted on the scenario (no imputation).  Tau bounds are
## narrowed so the fit is fast but still selects within an interval.
.pf_reference_model <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    sc <- .pf_scenario()
    cache <<- suppressWarnings(suppressMessages(fit_reference_model(
      sc$reference, hydro = sc$rain,
      min_obs_model = 20L,
      api_tau_bounds_short = c(3, 10),
      api_tau_bounds_long  = c(30, 90)
    )))
    cache
  }
})

## Arguments for the production-shaped chronic pipeline call.
.pf_pipeline_args <- function() {
  sc <- .pf_scenario()
  list(
    target    = sc$target,
    reference = sc$reference,
    hydro     = sc$rain,
    impute    = TRUE,
    impute_seed = 1L,
    reference_args = list(
      min_obs_model = 20L,
      api_tau_bounds_short = c(3, 10),
      api_tau_bounds_long  = c(30, 90)
    ),
    daily_args = list(
      reference   = sc$reference, # static ARA (prepare_reference inside add_mspaf)
      temperature = sc$temperature
    ),
    chronic  = TRUE,
    focal_by = "day",
    tau      = 90,
    window   = 365
  )
}

.pf_golden_path <- function(name) {
  testthat::test_path("_goldens", paste0(name, ".rds"))
}

## ---------------------------------------------------------------------------
## Golden case builders.  Each returns a NAMED LIST of results computed by
## whatever code is loaded; the generator saves them, the tests recompute and
## compare.  Keep these deterministic (fixed seeds, no Sys.time()).
## ---------------------------------------------------------------------------

## time_weighted_aggregate(): many settings over the grab-sample fixture.
.pf_twa_cases <- function() {
  sc  <- .pf_scenario()
  tgt <- sc$target
  # imputed flag with some TRUE and some NA (NA must not count as imputed).
  set.seed(51L)
  tgt$imputed <- stats::runif(nrow(tgt)) < 0.2
  tgt$imputed[stats::runif(nrow(tgt)) < 0.05] <- NA
  # A few NA values (dropped inside the weighted summary).
  tgt$value[which(tgt$analyte == "Zn")[c(4, 9)]] <- NA_real_

  # Focal dates: before the record, on the first sample, on sample dates
  # (zero-length final interval), around the duplicate same-day pair, inside
  # P3's > 365-day gap, after the record, plus a weekly sweep.
  edge <- as.Date(c(
    "2020-06-01", "2020-12-31", "2021-01-12", "2021-02-02", "2021-06-05",
    "2021-06-06", "2022-08-15", "2022-08-16", "2022-10-25", "2022-12-04",
    "2022-12-05", "2023-03-04", "2023-12-12", "2024-03-01", "2025-01-20"
  ))
  sweep <- seq(as.Date("2020-12-01"), as.Date("2024-02-01"), by = "3 weeks")
  # Every sample date of P1 is also a focal date (window edges hit exactly).
  p1_dates <- sort(unique(as.Date(tgt$datetime[tgt$site_id == "P1"])))
  focal <- sort(unique(c(edge, sweep)))

  twa <- function(df, ...) {
    suppressMessages(time_weighted_aggregate(df, ...))
  }

  out <- list()
  out$geom_90_365 <- twa(tgt, focal, tau = 90, tau_units = "d",
                         window = 365, window_units = "d")
  out$arith_30_90 <- twa(tgt, focal, tau = 30, tau_units = "d",
                         window = 90, window_units = "d",
                         summary = "arith_mean")
  out$p90_365_730 <- twa(tgt, focal, tau = 365, tau_units = "d",
                         window = 730, window_units = "d", summary = "p90")
  # Defaults (tau = 90 d, window = 365 d) and `units` objects.
  out$defaults <- twa(tgt, edge)
  out$units_weeks <- twa(tgt, edge,
                         tau = units::set_units(6, "week"),
                         window = units::set_units(26, "week"))
  # No anchor: the leading edge sees only in-window samples.
  out$no_anchor <- twa(tgt, focal, tau = 90, tau_units = "d",
                       window = 120, window_units = "d",
                       anchor_outside_window = FALSE)
  # Window shorter than the inter-sample spacing: many windows hold only the
  # anchor-clipped interval or a single sample; eps changes the geomean guard.
  out$short_window <- twa(tgt, focal, tau = 5, tau_units = "d",
                          window = 10, window_units = "d", eps = 1e-3)
  # Focal dates exactly on P1's sample dates (sample at window_start and at
  # focal_date), window = 60 d.
  out$on_sample_dates <- twa(tgt[tgt$site_id == "P1", ], p1_dates,
                             tau = 20, tau_units = "d",
                             window = 60, window_units = "d",
                             summary = "arith_mean")
  # Per-site focal-date frame, including a site absent from df and repeated
  # rows (output follows the frame's row order).
  fd <- tibble::tibble(
    focal_date = as.Date(c("2022-03-01", "2021-07-01", "2023-06-30",
                           "2022-03-01", "2023-01-15", "2023-01-15")),
    site_id    = c("P1", "P2", "P3", "P1", "PX", "P3")
  )
  out$focal_frame <- twa(tgt, fd, tau = 90, tau_units = "d",
                         window = 365, window_units = "d")
  # POSIXct (UTC) datetimes, no imputed column.
  tgt_ct <- sc$target
  tgt_ct$datetime <- as.POSIXct(paste(tgt_ct$datetime, "10:30:00"), tz = "UTC")
  out$posixct_no_imputed <- twa(tgt_ct, edge, tau = 90, tau_units = "d",
                                window = 365, window_units = "d")
  # Daily-series shape (as in mspaf_pipeline's chronic step): one analyte,
  # one sample per day with gaps, a leading edge and NA values.
  set.seed(52L)
  dd <- seq(as.Date("2022-01-01"), as.Date("2023-06-30"), by = "day")
  dd <- dd[!(dd >= as.Date("2022-07-01") & dd <= as.Date("2022-07-20"))]
  daily <- dplyr::bind_rows(lapply(c("P1", "P2"), function(s) {
    keep <- if (s == "P2") dd >= as.Date("2022-04-01") else rep(TRUE, length(dd))
    v <- 20 + 10 * sin(seq_along(dd[keep]) / 30) + stats::rnorm(sum(keep))
    v[sample.int(length(v), 5)] <- NA
    tibble::tibble(sample_id = paste0("daily_", s, "_", dd[keep]),
                   site_id = s, datetime = dd[keep], analyte = "msPAF",
                   detected = TRUE, value = v)
  }))
  out$daily_series <- twa(daily,
                          seq(as.Date("2022-01-01"), as.Date("2023-07-15"),
                              by = "3 days"),
                          tau = 90, tau_units = "d", window = 365,
                          window_units = "d", summary = "arith_mean")
  # Draws mode: a draw-bearing msPAF series (summary and raw draws).
  set.seed(53L)
  dr <- daily[daily$site_id == "P1" & daily$datetime <= as.Date("2022-03-31"), ]
  dr <- dplyr::bind_rows(lapply(1:4, function(d) {
    dplyr::mutate(dr, value = .data$value * exp(stats::rnorm(dplyr::n(), 0, 0.1)),
                  draw_id = d)
  }))
  fdr <- seq(as.Date("2022-01-10"), as.Date("2022-04-10"), by = "10 days")
  out$draws_summary <- twa(dr, fdr, tau = 30, tau_units = "d", window = 60,
                           window_units = "d", summary = "arith_mean")
  out$draws_raw <- twa(dr, fdr, tau = 30, tau_units = "d", window = 60,
                       window_units = "d", summary = "arith_mean",
                       return = "draws")
  out$draws_mean_80 <- twa(dr, fdr, tau = 30, tau_units = "d", window = 60,
                           window_units = "d", summary = "geom_mean",
                           interval = 0.8, central = "mean")
  # Mixed exact + draw cells: exact rows are broadcast across draws.
  mixed <- dplyr::bind_rows(
    dr,
    dplyr::mutate(dr[dr$draw_id == 1L, ], analyte = "pH", value = 7,
                  draw_id = NA_integer_)
  )
  out$draws_mixed <- twa(mixed, fdr, tau = 30, tau_units = "d", window = 60,
                         window_units = "d", return = "draws")
  out
}

## .compute_api(): dates before the rain record, inside missing-row gaps, on NA
## days, after the record, duplicated and unsorted targets; unsorted hydro with
## duplicate-date rows; several tau.
.pf_api_cases <- function() {
  sc   <- .pf_scenario()
  rain <- sc$rain
  api  <- hydroSense:::.compute_api
  targets <- as.Date(c(
    "2020-01-01", "2020-08-31", "2020-09-01", "2020-09-02", "2021-03-15",
    "2021-08-09", "2021-08-10", "2021-08-17", "2021-08-25", "2021-08-26",
    "2022-11-02", "2022-11-05", "2023-12-31", "2024-01-01", "2024-06-30",
    "2021-03-15", "2020-01-01"
  ))
  na_days <- rain$date[is.na(rain$value)][1:5]
  targets <- c(targets, na_days, na_days + 1L)
  daily_grid <- seq(as.Date("2020-07-01"), as.Date("2024-02-01"), by = "day")

  # Unsorted hydro with duplicate-date rows (duplicates sum).
  set.seed(61L)
  perm <- sample.int(nrow(rain))
  hv <- c(rain$value[perm], 5, 7.5, NA)
  hd <- c(rain$date[perm], as.Date(c("2021-03-10", "2021-03-10", "2021-03-11")))

  list(
    tau7      = api(rain$value, rain$date, targets, 7),
    tau60     = api(rain$value, rain$date, targets, 60),
    tau_frac  = api(rain$value, rain$date, targets, 2.5),
    grid_tau5 = api(rain$value, rain$date, daily_grid, 5),
    grid_tau90 = api(rain$value, rain$date, daily_grid, 90),
    unsorted_dups = api(hv, hd, targets, 7),
    empty_hydro = api(numeric(0), as.Date(character(0)), targets, 7),
    single_target = api(rain$value, rain$date, as.Date("2022-02-22"), 30),
    no_targets = api(rain$value, rain$date, as.Date(character(0)), 30)
  )
}

## .resolve_ref_norm_instant(): matched / model / static tiers, dates before the
## hydro record and in rain gaps, POSIXct datetimes, duplicate rows.
.pf_resolver_cases <- function() {
  sc <- .pf_scenario()
  rm <- .pf_reference_model()
  ref_dates <- sort(unique(as.Date(sc$reference$datetime)))
  d <- sort(unique(c(
    as.Date(c("2020-06-01", "2020-08-31", "2020-09-01", "2021-08-12",
              "2021-08-26", "2022-11-03", "2024-01-15")),
    ref_dates[c(3, 10, 25)] + c(0L, 2L, -5L),  # within +-5 d of a ref grab
    ref_dates[c(4, 30)] + 6L,                  # just outside the match window
    seq(as.Date("2022-01-01"), as.Date("2022-03-31"), by = "day")
  )))
  tdf <- tibble::tibble(sample_id = paste0("q", seq_along(d)), datetime = d)
  # Duplicate (sample_id, date) rows collapse; POSIXct datetimes are coerced.
  tdf2 <- tibble::tibble(
    sample_id = c("a", "a", "b", "c"),
    datetime  = as.POSIXct(c("2022-05-01 09:00", "2022-05-01 15:00",
                             "2020-07-01 00:00", "2023-02-14 12:00"),
                           tz = "UTC"),
    extra = 1:4
  )
  list(
    grid = hydroSense:::.resolve_ref_norm_instant(rm, tdf),
    dups_posixct = hydroSense:::.resolve_ref_norm_instant(rm, tdf2)
  )
}

## Row-order-free view of a resolver result (the resolver's row order is not
## part of its contract: every consumer joins or looks up by key).
.pf_sort_resolver <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[order(x$sample_id, x$analyte, method = "radix"), , drop = FALSE] |>
    `rownames<-`(NULL)
}

## The production-shaped chronic pipeline.  The daily grid is restricted to
## ~6 months (start/end) so the test stays fast; the models are still fitted
## on the full record and the chronic windows run off the start of the daily
## series.  Returns the chronic frame plus the daily (mspaf_daily point-mode)
## frame and its per-analyte PAF breakdown.
.pf_pipeline_case <- function() {
  args <- .pf_pipeline_args()
  args$daily_args$start <- as.Date("2023-06-01")
  args$daily_args$end   <- as.Date("2023-12-12")
  # The external temperature series deliberately extends beyond the daily grid
  # at both ends: mspaf_daily() adds those days to the per-site date set
  # (.fill_external_temperature() does not clip to start/end), and that set
  # drives the smoother grid and the hydro z-score standardisation.
  tmp <- args$daily_args$temperature
  args$daily_args$temperature <- tmp[tmp$datetime >= as.Date("2023-04-20"), ]
  out <- suppressWarnings(suppressMessages(do.call(mspaf_pipeline, args)))
  daily <- attr(out, "daily")
  apafs <- attr(daily, "analyte_pafs")
  apafs <- as.data.frame(apafs)
  apafs <- apafs[order(apafs$site_id, apafs$date, apafs$analyte,
                       method = "radix"), , drop = FALSE]
  rownames(apafs) <- NULL
  attributes(out)[c("daily", "reference_model", "imputation_model")] <- NULL
  attr(daily, "analyte_pafs") <- NULL
  attr(daily, "ara_summary")  <- NULL
  list(chronic = out, daily = daily, analyte_pafs = apafs)
}

## Small 3-point target used by the structural (call-count) tests: ten grabs
## per point in 2022, so fit_target_model() is quick (bridge tier) while the
## daily-grid reference resolution still runs over every day.
.pf_small_target <- function() {
  sc <- .pf_scenario()
  mk <- function(site, start, step, mult, seed) {
    d <- seq(as.Date(start), by = paste(step, "days"), length.out = 10L)
    .pf_point_chem(site, d, mult = mult, seed = seed, ni_every = 1L,
                   bdl_frac = 0, rain = sc$rain)
  }
  dplyr::bind_rows(
    mk("P1", "2022-01-05", 30L, 4, 71L),
    mk("P2", "2022-01-20", 33L, 2, 72L),
    mk("P3", "2022-02-02", 28L, 6, 73L)
  )
}

## mspaf_daily() point mode on the small target: the daily grid starts before
## every point's first grab (leading-edge drop) and the external temperature
## series overhangs the grid at both ends.  Static ARA (raw reference frame,
## prepared inside add_mspaf()) and temporal ARA (the reference_model itself,
## resolved per synthetic daily sample inside add_mspaf()).
.pf_daily_cases <- function() {
  sc  <- .pf_scenario()
  rm  <- .pf_reference_model()
  tgt <- .pf_small_target()
  tmp <- sc$temperature
  tmp <- tmp[tmp$datetime >= as.Date("2021-12-15") &
               tmp$datetime <= as.Date("2022-04-10"), ]
  run <- function(reference) {
    out <- suppressWarnings(suppressMessages(mspaf_daily(
      tgt, temperature = tmp, reference = reference, reference_model = rm,
      interpolation = "model", start = as.Date("2022-01-01"),
      end = as.Date("2022-03-31")
    )))
    apafs <- as.data.frame(attr(out, "analyte_pafs"))
    apafs <- apafs[order(apafs$site_id, apafs$date, apafs$analyte,
                         method = "radix"), , drop = FALSE]
    rownames(apafs) <- NULL
    ara <- attr(out, "ara_summary")
    attr(out, "analyte_pafs") <- NULL
    attr(out, "ara_summary")  <- NULL
    res <- list(daily = out, analyte_pafs = apafs)
    if (!is.null(ara) && "ref_tier" %in% names(ara) &&
        any(!is.na(ara$ref_tier))) {
      ara <- as.data.frame(ara)
      ara <- ara[order(ara$sample_id, ara$analyte, method = "radix"), ,
                 drop = FALSE]
      rownames(ara) <- NULL
      res$ara_summary <- ara
    }
    res
  }
  list(static_ara = run(sc$reference), temporal_ara = run(rm))
}
