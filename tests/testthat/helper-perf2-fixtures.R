## ============================================================================
## Synthetic fixtures for the second msPAF performance pass ("perf2")
## ============================================================================
##
## Used by test-perf2-goldens.R, test-perf2-structure.R and the golden
## generator data-raw/make_perf2_goldens.R.  Builds on helper-perf-fixtures.R
## (the synthetic scenario, rainfall, temperature and reference model).
## Everything here is SYNTHETIC: fixed seeds, invented site codes ("P1",
## "P3", "R1", "S"), no monitoring data.
##
## The hot spots of a real point-mode mspaf_pipeline() run that this pass
## vectorises, and the case builders that pin each one:
##   .compute_daily_diag()                 .pf2_diag_cases()
##   .interpolate_analyte()                .pf2_interp_cases()
##   .build_daily_chem()                   .pf2_build_chem_cases()
##   .pivot_chem_wide()                    .pf2_pivot_cases()
##   .compute_pca_scores()                 .pf2_pca_cases()
##   .compute_hydro_features()             .pf2_hydro_cases()
##   fit_target_model() (tau search,       .pf2_fit_cases()
##     .fit_pooled_impact_response(),
##     .fit_impact_response(),
##     .optimise_tau_pair() trace,
##     .resolve_target_impact())
##   .fit_daily_target() / .predict_daily_tox() /
##     .daily_tox_from_model()             .pf2_daily_model_cases()
##   mspaf_pipeline() point mode           .pf2_pipeline_point_cases()

## ---------------------------------------------------------------------------
## Small shared inputs
## ---------------------------------------------------------------------------

## SSD-eligible analytes (log-space interpolation, diag counting).
.pf2_tox <- function() {
  meta <- hydroSense:::.load_analyte_metadata(NULL)
  meta$analyte[!is.na(meta$ssd_available) & meta$ssd_available == TRUE]
}

## P1 grabs before 2022-10-01: ~27 samples (pooled tau search runs for Cu,
## Zn and NH3-N; Ni has 11 anchors and stays a flat bridge), including the
## duplicate same-day pair of 2022-08-15 and BDL Cu/Ni.
.pf2_p1 <- function() {
  sc <- .pf_scenario()
  sc$target[sc$target$site_id == "P1" &
              sc$target$datetime < as.Date("2022-10-01"), ]
}

## Imputation model fitted on the synthetic reference (PCA + groups).
.pf2_imputation_model <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    sc <- .pf_scenario()
    cache <<- suppressWarnings(suppressMessages(fit_imputation_model(
      sc$reference, required_vars = c("pH", "EC"), seed = 1L
    )))
    cache
  }
})

## Forward-filled daily chemistry for P1 over 2022 (the "model" path's input
## to .daily_tox_from_model()), with the external temperature series filled
## in.  The grid starts before nothing (P1 is sampled from 2021) but extends
## past the last grab (2022-09-xx) to exercise the trailing carry-forward.
.pf2_daily_long <- function() {
  sc <- .pf_scenario()
  site_rows <- .pf2_p1()
  site_rows$datetime <- as.Date(site_rows$datetime)
  dates <- seq(as.Date("2022-01-01"), as.Date("2022-11-30"), by = "day")
  dl <- hydroSense:::.build_daily_chem(site_rows, dates, "forward_fill",
                                       "drop", .pf2_tox())
  tmp <- sc$temperature
  tmp <- tmp[tmp$datetime >= as.Date("2021-12-20") &
               tmp$datetime <= as.Date("2022-12-10"), ]
  hydroSense:::.fill_external_temperature(dl, tmp)
}

## Independent reference API (does not call package code, so a stale cache
## inside .compute_api() cannot also corrupt the expectation): NA rain is 0,
## duplicate-date rows sum, first-order recursion on a complete daily grid.
.pf2_ref_api <- function(values, dates, targets, tau) {
  dates <- as.Date(dates)
  targets <- as.Date(targets)
  if (length(dates) == 0L) return(rep(0, length(targets)))
  values[is.na(values)] <- 0
  g0 <- min(dates)
  span <- as.integer(max(c(dates, targets)) - g0) + 1L
  rain <- numeric(span)
  pos <- as.integer(dates - g0) + 1L
  for (i in seq_along(pos)) rain[pos[i]] <- rain[pos[i]] + values[i]
  state <- as.numeric(stats::filter(rain, exp(-1 / tau), method = "recursive"))
  tp <- as.integer(targets - g0) + 1L
  out <- numeric(length(targets))
  ok <- tp >= 1L & tp <= span
  out[ok] <- state[tp[ok]]
  out
}

.pf2_golden_path <- function(name) {
  testthat::test_path("_goldens", paste0("perf2_", name, ".rds"))
}

## Row-order-free view of a frame keyed by the given columns.
.pf2_sort <- function(x, keys) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  o <- do.call(order, c(unname(as.list(x[keys])), list(method = "radix")))
  x <- x[o, , drop = FALSE]
  rownames(x) <- NULL
  x
}

## ---------------------------------------------------------------------------
## .compute_daily_diag()
## ---------------------------------------------------------------------------

## Synthetic daily_long: n days x k tox analytes plus one co-analyte, grabs
## every `every` days (all tox measured on a grab day).  Used by the golden
## cases (small n) and the scaling tests (large n).
.pf2_diag_long <- function(n, k = 4L, every = 9L, start = "2001-03-01") {
  d  <- seq(as.Date(start), by = "day", length.out = n)
  an <- c(paste0("T", seq_len(k)), "pH")
  g  <- expand.grid(.date = d, analyte = an, KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
  g$value     <- 1
  g$detected  <- TRUE
  g$.measured <- (as.integer(g$.date - d[1L]) %% every) == 3L
  tibble::as_tibble(g)
}

.pf2_diag_cases <- function() {
  diag <- hydroSense:::.compute_daily_diag
  tox  <- .pf2_tox()
  out  <- list()

  ## Real daily chemistry (forward-filled + temperature fill; rows appended
  ## out of date order by the temperature fill).
  dl <- .pf2_daily_long()
  out$p1_forward_fill <- diag(dl, tox, "P1")

  ## Leading edge backfilled (days before the first grab: NA days_since).
  sc <- .pf_scenario()
  p3 <- sc$target[sc$target$site_id == "P3", ]
  p3$datetime <- as.Date(p3$datetime)
  dl3 <- hydroSense:::.build_daily_chem(
    p3, seq(as.Date("2020-12-01"), as.Date("2022-12-31"), by = "day"),
    "linear", "backfill", tox
  )
  out$p3_backfill <- diag(dl3, tox, "P3")

  ## Synthetic: shuffled rows, duplicated measured (date, analyte) rows (n
  ## counts ROWS), NA .measured, a measured co-analyte that must be ignored,
  ## and tox analytes passed as a superset of those present.
  set.seed(91L)
  sx <- .pf2_diag_long(60L, k = 3L, every = 7L)
  sx$analyte[sx$analyte == "T1"] <- "Cu"
  sx$analyte[sx$analyte == "T2"] <- "Zn"
  sx$analyte[sx$analyte == "T3"] <- "Ni"
  sx$.measured[sx$analyte == "pH" & as.integer(format(sx$.date, "%d")) == 1L] <- TRUE
  dup <- sx[sx$analyte == "Cu" & sx$.measured, ][1:2, ]
  sx <- rbind(sx, dup)
  sx$.measured[which(sx$analyte == "Zn")[c(5, 12, 40)]] <- NA
  sx <- sx[sample.int(nrow(sx)), ]
  out$synthetic_messy <- diag(sx, c("Cu", "Zn", "Ni", "Cd"), "S")

  ## Tox rows but no grab anywhere: every days_since is NA.
  nm <- .pf2_diag_long(15L, k = 2L)
  nm$analyte[nm$analyte == "T1"] <- "Cu"
  nm$analyte[nm$analyte == "T2"] <- "Zn"
  nm$.measured <- FALSE
  out$no_grabs <- diag(nm, c("Cu", "Zn"), "S")

  ## No tox rows at all: the empty schema.
  out$empty <- diag(nm[nm$analyte == "pH", ], c("Cu", "Zn"), "S")
  out
}

## ---------------------------------------------------------------------------
## .interpolate_analyte()
## ---------------------------------------------------------------------------

.pf2_interp_cases <- function() {
  ia <- hydroSense:::.interpolate_analyte
  od <- as.Date(c("2022-03-10", "2022-01-05", "2022-02-01", "2022-02-01",
                  "2022-04-20", "2022-01-20"))
  ov <- c(4, 2, 0, 9, 1.5, 3)       # includes a 0 (log-space floor) and dup
  odt <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
  td <- seq(as.Date("2021-12-25"), as.Date("2022-05-10"), by = "day")
  td_messy <- c(td[c(80, 3, 3, 40, 120)], as.Date("2022-01-20"),
                as.Date("2021-01-01"))
  out <- list()
  for (interp in c("forward_fill", "linear")) {
    for (edge in c("drop", "backfill")) {
      for (lg in c(FALSE, TRUE)) {
        nm <- paste(interp, edge, if (lg) "log" else "lin", sep = "_")
        out[[nm]] <- ia(od, ov, odt, td, interp, edge, log_space = lg)
      }
    }
  }
  out$messy_targets <- ia(od, ov, odt, td_messy, "linear", "backfill", TRUE)
  out$single_obs <- ia(as.Date("2022-02-01"), 5, FALSE, td, "linear",
                       "backfill", TRUE)
  out$empty_obs <- ia(as.Date(character()), numeric(), logical(), td,
                      "linear", "drop", FALSE)
  out$no_targets <- ia(od, ov, odt, as.Date(character()), "linear",
                       "drop", FALSE)
  out
}

## ---------------------------------------------------------------------------
## .build_daily_chem()
## ---------------------------------------------------------------------------

.pf2_build_chem_cases <- function() {
  bdc <- hydroSense:::.build_daily_chem
  sc  <- .pf_scenario()
  tox <- .pf2_tox()
  rows <- sc$target[sc$target$site_id %in% c("P1", "P3"), ]
  rows <- rows[rows$site_id == "P1" | rows$datetime < as.Date("2022-01-01"), ]
  rows$site_id <- "P1"
  rows$datetime <- as.Date(rows$datetime)
  ## Every Ni grab below detection on P3's dates -> falls back to all rows
  ## for that analyte only when there is no detected row at all: make a
  ## separate all-BDL analyte.
  extra <- rows[rows$analyte == "Ni", ]
  extra$analyte  <- "Cd"
  extra$detected <- FALSE
  extra$value    <- 0.2
  rows <- rbind(rows, extra)
  rows$valence.analyte <- ifelse(rows$analyte %in% c("Ca", "Mg"), 2, NA_real_)
  set.seed(92L)
  rows <- rows[sample.int(nrow(rows)), ]
  dates <- seq(as.Date("2020-12-20"), as.Date("2024-01-10"), by = "day")
  out <- list()
  for (interp in c("forward_fill", "linear")) {
    for (edge in c("drop", "backfill")) {
      out[[paste(interp, edge, sep = "_")]] <-
        bdc(rows, dates, interp, edge, tox)
    }
  }
  ## A grid that only partly overlaps the record, and no pass-through cols.
  bare <- rows[, c("sample_id", "site_id", "datetime", "analyte", "value",
                   "detected")]
  out$bare_partial <- bdc(bare, seq(as.Date("2023-06-01"),
                                    as.Date("2023-12-31"), by = "day"),
                          "linear", "drop", tox)
  out
}

## ---------------------------------------------------------------------------
## .pivot_chem_wide() and .compute_pca_scores()
## ---------------------------------------------------------------------------

.pf2_pivot_input <- function() {
  sc <- .pf_scenario()
  x <- sc$target[sc$target$site_id == "P1", ]
  set.seed(93L)
  ## Duplicate (sample_id, analyte) rows: concentration analytes collapse by
  ## geometric mean, no-log analytes by arithmetic mean.
  dup <- x[x$analyte %in% c("EC", "pH", "Ca") & x$sample_id %in% x$sample_id[1:6], ]
  dup$value <- dup$value * 1.3
  ## A cell whose every row is BDL (concentration -> NA, not NaN), a zero
  ## value (floored before the log), and BDL pH (kept: no-log variable).
  bdl <- x[x$analyte == "Mg" & x$sample_id == x$sample_id[10], ]
  x$detected[x$analyte == "Mg" & x$sample_id == x$sample_id[10]] <- FALSE
  bdl$detected <- FALSE
  x$value[x$analyte == "DOC" & x$sample_id == x$sample_id[3]] <- 0
  x$detected[x$analyte == "pH" & x$sample_id == x$sample_id[4]] <- FALSE
  x <- rbind(x, dup, bdl)
  x[sample.int(nrow(x)), ]
}

.pf2_pivot_cases <- function() {
  pcw <- hydroSense:::.pivot_chem_wide
  x <- .pf2_pivot_input()
  nolog <- hydroSense:::.PCA_NO_LOG_VARS
  wq <- c("pH", "EC", "DOC", "Ca", "Mg", "hardness", "temperature", "NH3-N",
          "Alk")                              # Alk absent from the data
  nodet <- x[, setdiff(names(x), "detected")]
  ## Daily WQ shape (sample_id = date string, no detected column), as built
  ## by .fit_daily_target() for the WQ layer.
  dl <- .pf2_daily_long()
  co <- dl[!dl$analyte %in% .pf2_tox(), ]
  daily_wq <- tibble::tibble(sample_id = as.character(co$.date),
                             analyte = co$analyte, value = co$value)
  list(
    mixed      = pcw(x, wq, nolog),
    no_detected = pcw(nodet, wq, nolog),
    subset_vars = pcw(x, c("Mg", "pH"), nolog),
    daily_wq   = pcw(daily_wq, wq, nolog),
    no_match   = pcw(x, c("Alk", "Fe"), nolog)
  )
}

.pf2_pca_cases <- function() {
  cps <- hydroSense:::.compute_pca_scores
  pca <- .pf2_imputation_model()$pca
  sc  <- .pf_scenario()
  x   <- .pf2_pivot_input()
  dl  <- .pf2_daily_long()
  co  <- dl[!dl$analyte %in% .pf2_tox(), ]
  daily_wq <- tibble::tibble(sample_id = as.character(co$.date),
                             analyte = co$analyte, value = co$value)
  ## A frame where one active variable is entirely absent (training median
  ## fill) and one is present but all NA.
  partial <- sc$target[sc$target$site_id == "P2" & sc$target$analyte != "DOC", ]
  partial$value[partial$analyte == "Mg"] <- NA_real_
  list(
    target     = cps(sc$target, pca),
    messy      = cps(x, pca),
    daily_wq   = cps(daily_wq, pca),
    partial    = cps(partial, pca)
  )
}

## ---------------------------------------------------------------------------
## .compute_hydro_features()
## ---------------------------------------------------------------------------

.pf2_hydro_cases <- function() {
  chf <- hydroSense:::.compute_hydro_features
  sc  <- .pf_scenario()
  td  <- c(seq(as.Date("2020-08-15"), as.Date("2024-01-20"), by = "5 days"),
           as.Date(c("2021-08-12", "2021-08-12", "2020-01-01")))
  stage <- sc$rain
  stage$value <- 1 + cumsum(dplyr::coalesce(stage$value, 0)) / 1000
  list(
    rain_short_long = chf(sc$rain, td, 4.5, 80, "rainfall"),
    rain_equal_tau  = chf(sc$rain, td, 30, 30, "rainfall"),
    stage           = chf(stage, td, 5, 40, "stage"),
    rain_empty_dates = chf(sc$rain, as.Date(character()), 3, 30, "rainfall")
  )
}

## ---------------------------------------------------------------------------
## fit_target_model(): tau search, pooled / per-analyte / WQ-layer fits
## ---------------------------------------------------------------------------

## Golden view of a target_model: scalars per analyte in one frame, anchors /
## d_anchors frames, impact / WQ GAM coefficients and AICs.
.pf2_tm_view <- function(tm) {
  nms <- names(tm$models)
  num_or_na <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
  scal <- data.frame(
    analyte     = nms,
    tau_short   = vapply(nms, function(n) tm$models[[n]]$tau_short, numeric(1)),
    tau_long    = vapply(nms, function(n) tm$models[[n]]$tau_long, numeric(1)),
    tier        = vapply(nms, function(n) tm$models[[n]]$tier, character(1)),
    n_obs       = vapply(nms, function(n) as.integer(tm$models[[n]]$n_obs), integer(1)),
    scale_c     = vapply(nms, function(n) num_or_na(tm$models[[n]]$scale_c), numeric(1)),
    pooled      = vapply(nms, function(n) isTRUE(tm$models[[n]]$pooled), logical(1)),
    pool_levels = vapply(nms, function(n) {
      paste(tm$models[[n]]$pool_levels %||% character(), collapse = "|")
    }, character(1)),
    pool_center = vapply(nms, function(n) num_or_na(tm$models[[n]]$pool_center), numeric(1)),
    pool_scale  = vapply(nms, function(n) num_or_na(tm$models[[n]]$pool_scale), numeric(1)),
    has_impact_fit = vapply(nms, function(n) !is.null(tm$models[[n]]$impact_fit), logical(1)),
    has_wq_fit  = vapply(nms, function(n) !is.null(tm$models[[n]]$wq_fit), logical(1)),
    stringsAsFactors = FALSE
  )
  rownames(scal) <- NULL
  out <- list(models = scal, pc_cols = tm$pc_cols, hydro = tm$hydro)
  for (n in nms) {
    m <- tm$models[[n]]
    out[[paste0("anchors:", n)]] <- m$anchors
    if (!is.null(m$d_anchors)) out[[paste0("d_anchors:", n)]] <- m$d_anchors
    if (!is.null(m$impact_fit)) {
      out[[paste0("impact_coef:", n)]] <- stats::coef(m$impact_fit)
      out[[paste0("impact_aic:", n)]]  <- stats::AIC(m$impact_fit)
    }
    if (!is.null(m$wq_fit)) {
      out[[paste0("wq_coef:", n)]] <- stats::coef(m$wq_fit)
      out[[paste0("wq_aic:", n)]]  <- stats::AIC(m$wq_fit)
    }
  }
  out
}

## Run `expr` while recording every AIC evaluation the tau optimiser makes:
## .optimise_tau_pair() is wrapped so its aic_fn logs (call, tau_short,
## tau_long, aic).  Returns list(value, trace).
.pf2_with_tau_trace <- function(expr_fn) {
  log <- new.env()
  log$rows <- list()
  log$call <- 0L
  orig <- hydroSense:::.optimise_tau_pair
  testthat::local_mocked_bindings(
    .optimise_tau_pair = function(aic_fn, tau_bounds_short, tau_bounds_long, ...) {
      log$call <- log$call + 1L
      k <- log$call
      wrapped <- function(ts, tl) {
        a <- aic_fn(ts, tl)
        log$rows[[length(log$rows) + 1L]] <- data.frame(
          call = k, tau_short = ts, tau_long = tl, aic = a
        )
        a
      }
      orig(wrapped, tau_bounds_short, tau_bounds_long, ...)
    },
    .package = "hydroSense"
  )
  value <- expr_fn()
  trace <- if (length(log$rows)) do.call(rbind, log$rows) else
    data.frame(call = integer(), tau_short = numeric(), tau_long = numeric(),
               aic = numeric())
  list(value = value, trace = trace)
}

## Query grid for .resolve_target_impact() on the fitted models.
.pf2_query <- function() {
  tibble::tibble(date = seq(as.Date("2021-01-01"), as.Date("2022-10-15"),
                            by = "3 days"))
}

.pf2_fit_cases <- function() {
  rm  <- .pf_reference_model()
  p1  <- .pf2_p1()
  fit <- function(...) {
    suppressWarnings(suppressMessages(fit_target_model(p1, rm, ...)))
  }
  rti <- function(tm) {
    suppressWarnings(hydroSense:::.resolve_target_impact(tm, .pf2_query()))
  }
  out <- list()

  ## A: pooled factor-smooth fit, default tau bounds (the production call of
  ##    .fit_daily_target()).
  a <- .pf2_with_tau_trace(function() fit())
  out$pooled_default <- c(.pf2_tm_view(a$value),
                          list(tau_trace = a$trace, resolve = rti(a$value)))

  ## B: per-analyte fits (pool = FALSE), narrowed bounds.
  b <- .pf2_with_tau_trace(function() {
    fit(pool = FALSE, api_tau_bounds_short = c(3, 10),
        api_tau_bounds_long = c(30, 90))
  })
  out$per_analyte <- c(.pf2_tm_view(b$value),
                       list(tau_trace = b$trace, resolve = rti(b$value)))

  ## C: impute-first + WQ layer (PCA), pooled, narrowed bounds.
  c3 <- .pf2_with_tau_trace(function() {
    fit(imputation_model = .pf2_imputation_model(),
        api_tau_bounds_short = c(3, 10), api_tau_bounds_long = c(30, 90))
  })
  out$pooled_wq <- c(.pf2_tm_view(c3$value),
                     list(tau_trace = c3$trace, resolve = rti(c3$value)))
  out
}

## ---------------------------------------------------------------------------
## .fit_daily_target() / .predict_daily_tox() / .daily_tox_from_model()
## ---------------------------------------------------------------------------

## Golden view of an fdm scaffold (everything except model objects).
.pf2_fdm_view <- function(fdm) {
  out <- list(
    modelled     = fdm$modelled,
    qdates       = fdm$qdates,
    measured_key = fdm$measured_key,
    ref_q        = .pf2_sort(fdm$ref_q, c("date", "analyte")),
    fac_lut      = fdm$fac_lut,
    tm_models    = .pf2_tm_view(fdm$tm)$models
  )
  for (nm in fdm$modelled) {
    s <- fdm$static[[nm]]
    out[[paste0("static:", nm)]] <- data.frame(
      hydro_short = s$hydro_short, hydro_long = s$hydro_long,
      ref_vec = s$ref_vec, fac = as.numeric(s$fac)
    )
    out[[paste0("static_lp:", nm)]] <- if (is.null(s$lp)) NULL else unname(c(s$lp))
    sm <- fdm$smoothers[[nm]]
    out[[paste0("smoother:", nm)]] <- data.frame(
      grid_dates = sm$grid_dates, mean = sm$mean, var = sm$var
    )
    out[[paste0("smoother_anchor_dates:", nm)]] <- sm$anchor_dates
  }
  out[["co_grab_map"]] <- fdm$co_grab_map
  out[["co_split_n"]] <- vapply(fdm$co_split, nrow, integer(1))
  out
}

.pf2_daily_model_inputs <- function() {
  site_rows <- .pf2_p1()
  site_rows$datetime <- as.Date(site_rows$datetime)
  list(
    site_rows = site_rows,
    daily_long = .pf2_daily_long(),
    meta = hydroSense:::.load_analyte_metadata(NULL),
    tox = .pf2_tox()
  )
}

.pf2_daily_model_cases <- function() {
  inp <- .pf2_daily_model_inputs()
  rm  <- .pf_reference_model()
  ## (a) fit-once scaffold with the imputation model (WQ tier, PCA scoring of
  ##     the daily WQ frame) and its point-mode prediction.
  fdm <- suppressWarnings(suppressMessages(hydroSense:::.fit_daily_target(
    site_rows = inp$site_rows, reference_model = rm,
    imputation_model = .pf2_imputation_model(), conc_units = NULL,
    meta = inp$meta, tox_analytes = inp$tox, daily_long = inp$daily_long
  )))
  pred <- suppressWarnings(hydroSense:::.predict_daily_tox(fdm))
  tiers <- attr(pred, "impact_tiers")
  attr(pred, "impact_tiers") <- NULL
  ## (b) the point-mode wrapper without an imputation model.
  dtm <- suppressWarnings(suppressMessages(hydroSense:::.daily_tox_from_model(
    inp$daily_long, inp$site_rows, rm, NULL, NULL, inp$meta, inp$tox
  )))
  dtm_tiers <- attr(dtm, "impact_tiers")
  attr(dtm, "impact_tiers") <- NULL
  list(
    fdm_wq   = .pf2_fdm_view(fdm),
    predict_wq = pred,
    predict_wq_tiers = as.data.frame(tiers),
    daily_tox = dtm,
    daily_tox_tiers = as.data.frame(dtm_tiers)
  )
}

## ---------------------------------------------------------------------------
## mspaf_pipeline() point mode (chronic = FALSE), impute = TRUE
## ---------------------------------------------------------------------------

## Targets P1 (dense: pooled model tier) and P3 (sparse: bridge tier).
.pf2_pipeline_point_args <- function(ara) {
  args <- .pf_pipeline_args()
  sc <- .pf_scenario()
  args$target <- sc$target[sc$target$site_id %in% c("P1", "P3"), ]
  args$chronic <- FALSE
  args$focal_by <- NULL
  args$tau <- NULL
  args$window <- NULL
  args$daily_args$reference <- if (ara) sc$reference else NULL
  args$daily_args$start <- as.Date("2023-05-01")
  args$daily_args$end   <- as.Date("2023-09-30")
  tmp <- args$daily_args$temperature
  args$daily_args$temperature <- tmp[tmp$datetime >= as.Date("2023-04-20"), ]
  args[!vapply(args, is.null, logical(1))]
}

.pf2_pipeline_point_case <- function(ara) {
  args <- .pf2_pipeline_point_args(ara)
  out <- suppressWarnings(suppressMessages(do.call(mspaf_pipeline, args)))
  apafs <- as.data.frame(attr(out, "analyte_pafs"))
  apafs <- apafs[order(apafs$site_id, apafs$date, apafs$analyte,
                       method = "radix"), , drop = FALSE]
  rownames(apafs) <- NULL
  tiers <- attr(out, "impact_tiers")
  ara_s <- attr(out, "ara_summary")
  daily <- out
  attributes(daily)[c("analyte_pafs", "ara_summary", "impact_tiers",
                      "reference_model", "imputation_model")] <- NULL
  res <- list(daily = daily, analyte_pafs = apafs)
  if (!is.null(tiers)) res$impact_tiers <- as.data.frame(tiers)
  if (!is.null(ara_s) && nrow(ara_s) > 0L) {
    keys <- intersect(c("sample_id", "site_id", "analyte"), names(ara_s))
    res$ara_summary <- .pf2_sort(ara_s, keys)
  }
  res
}

.pf2_pipeline_point_cases <- function() {
  list(ara = .pf2_pipeline_point_case(TRUE),
       total = .pf2_pipeline_point_case(FALSE))
}
