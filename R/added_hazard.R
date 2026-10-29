# ── Chronic added-hazard msPAF (mspaf_pipeline(ara_form = "added_hazard")) ────
#
# Every point (each target, and the reference) is built identically from its
# own samples on a daily grid: toxicants log-linear with a span cap
# (interpolate_loglinear()), co-analytes interpolated without a cap, and each
# toxicant normalised with that point's own co-analytes. Per (target, day) the
# added hazard D = max((msPAF_T - msPAF_R) / (1 - msPAF_R), 0) is formed and
# time-weighted with exp(-age / tau) over a look-back window, renormalised over
# the covered days. Everything is vectorised as (days x analytes) matrices and
# the time-weighting is a causal convolution, so a 25-year daily grid runs in
# seconds. See the roxygen of mspaf_pipeline() for the user-facing definition.

#' Validate the added-hazard arguments of mspaf_pipeline()
#'
#' Raises every argument error of the `ara_form = "added_hazard"` path before
#' anything is computed, and resolves `ara` and `impute`.
#'
#' @return A list `(ara, reference)`; `reference` is `NULL` when `ara =
#'   "none"`.
#' @keywords internal
.ah_validate <- function(reference, impute, impute_missing, daily_args,
                         chronic_summary, ara, ara_missing, min_coverage,
                         tau, window) {
  conflict <- function(msg, hint) {
    cli::cli_abort(c(msg, "i" = hint), class = "hydroSense_ara_conflict",
                   call = rlang::caller_env(2))
  }
  interp <- daily_args$interpolation %||% "loglinear"
  if (!identical(interp, "loglinear")) {
    conflict(
      "{.code ara_form = \"added_hazard\"} interpolates log-linearly.",
      "Drop {.code daily_args$interpolation} or set it to \\
       {.val loglinear}."
    )
  }
  if (!ara_missing) {
    ara <- match.arg(ara, c("static", "contemporaneous", "none"))
    if (ara == "static") {
      conflict(
        "{.code ara = \"static\"} is not available with \\
         {.code ara_form = \"added_hazard\"}.",
        "Use {.code ara = \"contemporaneous\"} (with a {.arg reference}) or \\
         {.code ara = \"none\"}."
      )
    }
    if (ara == "contemporaneous" && is.null(reference)) {
      conflict(
        "{.code ara = \"contemporaneous\"} needs a {.arg reference}.",
        "Supply a reference site, or use {.code ara = \"none\"}."
      )
    }
  } else {
    ara <- if (is.null(reference)) "none" else "contemporaneous"
  }
  if (!impute_missing && isTRUE(impute)) {
    conflict(
      "{.code ara_form = \"added_hazard\"} never imputes.",
      "Use {.code impute = FALSE}."
    )
  }
  if (!identical(chronic_summary, "arith_mean")) {
    conflict(
      "{.code ara_form = \"added_hazard\"} needs \\
       {.code chronic_summary = \"arith_mean\"}.",
      "The added hazard is time-weighted as an arithmetic mean."
    )
  }
  if (!is.null(reference)) {
    checkmate::assert_names(names(reference), must.include = "site_id",
                            .var.name = "reference")
    n_ref <- length(unique(reference$site_id))
    if (n_ref != 1L) {
      cli::cli_abort(c(
        "{.arg reference} must contain exactly one {.field site_id}.",
        "x" = "It contains {n_ref}."
      ), class = "hydroSense_reference_sites", call = rlang::caller_env())
    }
  }
  checkmate::assert_number(min_coverage, lower = 0, upper = 1)
  checkmate::assert_number(tau, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_number(window, lower = 0, finite = TRUE, null.ok = TRUE)
  span_cap <- daily_args$span_cap %||% 120
  checkmate::assert_number(span_cap, lower = 0, .var.name = "daily_args$span_cap")
  if (!is.null(daily_args$temperature)) {
    checkmate::assert_data_frame(daily_args$temperature,
                                 .var.name = "daily_args$temperature")
    checkmate::assert_names(names(daily_args$temperature),
                            must.include = c("datetime", "value"),
                            .var.name = "daily_args$temperature")
  }
  list(ara = ara, reference = if (ara == "none") NULL else reference)
}

#' Canonicalise, unit-convert and reduce a long chemistry frame
#'
#' One value per (site, analyte, date): when any row that day is detected,
#' the mean of the detected values (`detected = TRUE`); otherwise the smallest
#' reported limit (`detected = FALSE`). Toxicants are converted to ug/L.
#'
#' @return `list(red, dates)`: the reduced data frame (`site_id`, `analyte`,
#'   `date` (numeric day), `value`, `detected`) and the numeric dates of every
#'   input row (for the grid).
#' @keywords internal
.ah_reduce_input <- function(df, ssd_analytes, conc_units, arg) {
  checkmate::assert_names(names(df),
    must.include = c("site_id", "datetime", "analyte", "value"),
    .var.name = arg
  )
  df <- as.data.frame(.canonicalise_mspaf_input(df))
  if (!"detected" %in% names(df)) df$detected <- TRUE
  df$detected <- as.logical(df$detected)
  df$detected[is.na(df$detected)] <- TRUE
  df$value <- as.numeric(df$value)
  df <- .convert_df_tox_to_ugL(df, ssd_analytes, conc_units, arg)
  day <- as.numeric(as.Date(df$datetime))
  all_days <- day[!is.na(day)]
  keep <- !is.na(day) & !is.na(df$value)
  d <- data.frame(site_id = as.character(df$site_id[keep]),
                  analyte = as.character(df$analyte[keep]),
                  date = day[keep], value = df$value[keep],
                  detected = df$detected[keep], stringsAsFactors = FALSE)
  if (!nrow(d)) {
    return(list(red = d, dates = all_days))
  }
  g <- collapse::GRP(d, c("site_id", "analyte", "date"), sort = TRUE,
                     return.groups = TRUE)
  n_det <- collapse::fsum(as.integer(d$detected), g, use.g.names = FALSE)
  m_det <- collapse::fmean(ifelse(d$detected, d$value, NA_real_), g,
                           use.g.names = FALSE, na.rm = TRUE)
  m_min <- collapse::fmin(d$value, g, use.g.names = FALSE)
  red <- data.frame(
    site_id  = as.character(g$groups$site_id),
    analyte  = as.character(g$groups$analyte),
    date     = as.numeric(g$groups$date),
    value    = ifelse(n_det > 0, m_det, m_min),
    detected = n_det > 0,
    stringsAsFactors = FALSE
  )
  list(red = red, dates = all_days)
}

## Interpolate a co-analyte without a cap: detected positive anchors only,
## linear in value (or log(value)), constant beyond the first/last anchor; one
## anchor is a constant, none is NA.
.ah_co_series <- function(x, y, g, log) {
  if (!length(x)) return(rep(NA_real_, length(g)))
  o <- order(x)
  x <- x[o]
  y <- y[o]
  if (length(x) == 1L) return(rep(y, length(g)))
  yy <- if (log) base::log(y) else y
  r <- stats::approx(x, yy, xout = g, rule = 2, ties = mean)$y
  if (log) exp(r) else r
}

## Daily co-analyte series of one point, as a named list of grid vectors.
.ah_coanalytes <- function(red_site, g, extra = character(0)) {
  pos <- red_site[red_site$detected & red_site$value > 0, , drop = FALSE]
  one <- function(a, log) {
    x <- pos[pos$analyte == a, , drop = FALSE]
    .ah_co_series(x$date, x$value, g, log)
  }
  ## Hardness: measured where it exists, else 2.497 Ca + 4.118 Mg on dates
  ## where both are detected.
  hd <- pos[pos$analyte == "hardness", c("date", "value"), drop = FALSE]
  ca <- pos[pos$analyte == "Ca", c("date", "value"), drop = FALSE]
  mg <- pos[pos$analyte == "Mg", c("date", "value"), drop = FALSE]
  both <- setdiff(intersect(ca$date, mg$date), hd$date)
  if (length(both)) {
    hd <- rbind(hd, data.frame(
      date = both,
      value = 2.497 * ca$value[match(both, ca$date)] +
        4.118 * mg$value[match(both, mg$date)]
    ))
  }
  co <- list(
    pH          = one("pH", FALSE),
    temperature = one("temperature", FALSE),
    DOC         = one("DOC", TRUE),
    hardness    = .ah_co_series(hd$date, hd$value, g, TRUE),
    Ca          = one("Ca", TRUE),
    Mg          = one("Mg", TRUE)
  )
  for (a in setdiff(extra, names(co))) co[[a]] <- one(a, FALSE)
  co
}

## Daily series of the global temperature for one point (NULL when the series
## has nothing for it): daily mean, linear between days, constant beyond.
.ah_temperature <- function(temperature, site, g) {
  if (is.null(temperature)) return(NULL)
  tt <- as.data.frame(temperature)
  if ("site_id" %in% names(tt)) tt <- tt[tt$site_id == site, , drop = FALSE]
  d <- as.numeric(as.Date(tt$datetime))
  v <- as.numeric(tt$value)
  ok <- !is.na(d) & is.finite(v)
  if (!any(ok)) return(NULL)
  m <- tapply(v[ok], d[ok], mean)
  x <- as.numeric(names(m))
  y <- as.numeric(m)
  if (length(x) == 1L) return(rep(y, length(g)))
  stats::approx(x, y, xout = g, rule = 2)$y
}

## Normalise a vector of concentrations of analyte `k` (row of `ssd`).
.ah_normalise <- function(ssd, k, C, co) {
  pf <- ssd$parsed_formula[[k]]
  if (is.null(pf)) return(C)
  out <- .apply_normalisation(pf, C, co)
  if (length(out) != length(C)) out <- rep(out[1L], length(C))
  as.numeric(out)
}

## Build one point: raw / detected / normalised (days x toxicants) matrices
## and its co-analyte series. Below-detection anchors enter at
## `bdl_mult * LOR`.
.ah_build_point <- function(red_site, ssd, g, span_cap, bdl_mult, temp,
                            co_extra) {
  tox <- ssd$analyte
  n <- length(g)
  raw <- matrix(NA_real_, n, length(tox), dimnames = list(NULL, tox))
  det <- matrix(NA, n, length(tox), dimnames = list(NULL, tox))
  co <- .ah_coanalytes(red_site, g, co_extra)
  if (!is.null(temp)) co$temperature <- temp
  for (k in seq_along(tox)) {
    o <- red_site[red_site$analyte == tox[k], , drop = FALSE]
    if (!nrow(o)) next
    v <- ifelse(o$detected, o$value, bdl_mult * o$value)
    s <- .interp_loglinear_num(o$date, v, o$detected, g, span_cap, TRUE, g)
    raw[, k] <- s$value
    det[, k] <- s$detected
  }
  norm <- raw
  for (k in seq_along(tox)) {
    if (all(is.na(raw[, k]))) next
    norm[, k] <- .ah_normalise(ssd, k, raw[, k], co)
  }
  list(raw = raw, det = det, norm = norm, co = co)
}

## Static reference summary (raw ug/L) per toxicant for the gap-filler:
## max(geometric mean with below-detection results at LOR/2,
##     median(LOR of the below-detection results) / 2).
.ah_static_summary <- function(red_ref, tox) {
  vapply(tox, function(a) {
    o <- red_ref[red_ref$analyte == a, , drop = FALSE]
    if (!nrow(o)) return(NA_real_)
    v <- ifelse(o$detected, o$value, o$value / 2)
    v <- v[v > 0]
    gm <- if (length(v)) exp(mean(log(v))) else NA_real_
    if (all(o$detected)) return(gm)
    max(gm, stats::median(o$value[!o$detected]) / 2, na.rm = TRUE)
  }, numeric(1))
}

## msPAF (proportion) of each row of a (days x analytes) matrix of normalised
## concentrations: CA within mode-of-action group over components with C > 0,
## IA across groups. NA where no analyte is present.
.ah_mspaf_rows <- function(M, ssd) {
  n <- nrow(M)
  present <- rowSums(!is.na(M)) > 0L
  surv <- rep(1, n)
  for (grp in unique(ssd$moa_group)) {
    cols <- which(ssd$moa_group == grp)
    Mg <- M[, cols, drop = FALSE]
    pos <- !is.na(Mg) & Mg > 0
    Mz <- Mg
    Mz[!pos] <- 0
    tu <- rowSums(sweep(Mz, 2L, ssd$hc50[cols], "/"))
    ns <- rowSums(pos)
    sg <- rowSums(sweep(pos * 1, 2L, ssd$sigma[cols], "*")) / ns
    p <- numeric(n)
    hit <- tu > 0
    p[hit] <- stats::pnorm(log10(tu[hit]) / sg[hit])
    surv <- surv * (1 - p)
  }
  out <- 1 - surv
  out[!present] <- NA_real_
  out
}

## Single-substance PAF (proportion) of each cell; 0 for C <= 0, NA for NA.
.ah_paf_cells <- function(M, ssd) {
  P <- M
  for (k in seq_len(ncol(M))) {
    C <- M[, k]
    p <- rep(0, length(C))
    hit <- !is.na(C) & C > 0
    p[hit] <- stats::pnorm(log10(C[hit] / ssd$hc50[k]) / ssd$sigma[k])
    p[is.na(C)] <- NA_real_
    P[, k] <- p
  }
  P
}

## Causal exponentially weighted window sums of the columns of X at grid
## positions `k` (1-based; may fall outside the grid, where X is taken as 0):
## sum over a = 0..window of w[a + 1] * X[k - a, ].
.ah_window_sums <- function(X, k, w) {
  X <- as.matrix(X)
  n <- nrow(X)
  L <- length(w)
  out <- matrix(0, length(k), ncol(X))
  inside <- !is.na(k) & k >= 1L & k - (L - 1L) <= n
  if (!any(inside) || n == 0L) return(out)
  kmax <- max(k[inside])
  pad_end <- max(0L, kmax - n)
  Xp <- rbind(matrix(0, L - 1L, ncol(X)), X, matrix(0, pad_end, ncol(X)))
  Y <- stats::filter(Xp, w, method = "convolution", sides = 1L)
  Y <- matrix(as.numeric(Y), nrow = nrow(Xp))
  out[inside, ] <- Y[(L - 1L) + k[inside], , drop = FALSE]
  out
}

#' Run the chronic added-hazard msPAF
#'
#' Engine of `mspaf_pipeline(ara_form = "added_hazard")`; arguments are
#' validated by [.ah_validate()] first.
#'
#' @return See [mspaf_pipeline()].
#' @keywords internal
.added_hazard_pipeline <- function(target, reference, daily_args, chronic,
                                   focal_dates, focal_by, tau, window,
                                   also_total, contributions, min_coverage) {
  span_cap <- daily_args$span_cap %||% 120
  tau      <- tau %||% 90
  window   <- floor(window %||% 365)
  method   <- match.arg(daily_args$method %||% "multi", c("multi", "anzecc"))
  guideline_dir <- if ("guideline_dir" %in% names(daily_args)) {
    daily_args$guideline_dir
  } else {
    getOption("hydroSense.guideline_dir")
  }
  conc_units  <- daily_args$conc_units
  temperature <- daily_args$temperature

  ## SSD parameters: toxicants are the analytes that have them.
  meta <- .load_analyte_metadata(daily_args$analyte_metadata)
  ssd_all <- derive_ssd_params(meta, method, guideline_dir)
  ssd_all <- ssd_all[is.finite(ssd_all$hc50) & ssd_all$hc50 > 0 &
                       is.finite(ssd_all$sigma) & ssd_all$sigma > 0, ,
                     drop = FALSE]

  tg <- .ah_reduce_input(target, ssd_all$analyte, conc_units, "target")
  rf <- if (is.null(reference)) NULL else
    .ah_reduce_input(reference, ssd_all$analyte, conc_units, "reference")
  red_t <- tg$red
  red_r <- if (is.null(rf)) NULL else rf$red
  sites <- unique(as.character(target$site_id))

  ## Toxicants present in the data, in SSD-table order.
  present <- unique(c(red_t$analyte, red_r$analyte))
  ssd <- ssd_all[ssd_all$analyte %in% present, , drop = FALSE]
  tox <- ssd$analyte
  co_extra <- unique(unlist(lapply(ssd$coanalytes_req, function(cr) {
    if (is.na(cr) || !nzchar(cr)) character(0) else trimws(strsplit(cr, ",")[[1L]])
  })))
  co_extra <- co_extra[nzchar(co_extra)]

  ## ── Grid ──────────────────────────────────────────────────────────────────
  all_days <- c(tg$dates, if (!is.null(rf)) rf$dates)
  if (!length(tg$dates)) {
    cli::cli_abort("{.arg target} has no dated rows.")
  }
  start <- if (!is.null(daily_args$start)) as.numeric(as.Date(daily_args$start))
           else min(all_days)
  end <- if (!is.null(daily_args$end)) {
    as.numeric(as.Date(daily_args$end))
  } else {
    hi <- max(tg$dates)
    if (is.finite(span_cap)) hi + span_cap else hi
  }
  if (end < start) cli::cli_abort("The daily grid is empty: {.arg end} is before {.arg start}.")
  g <- seq(start, end)
  N <- length(g)
  grid <- as.Date(g, origin = "1970-01-01")

  ## Global temperature (no site_id column): one series for every point.
  temp_global <- if (!is.null(temperature) &&
                     !"site_id" %in% names(temperature)) {
    .ah_temperature(temperature, NULL, g)
  } else NULL
  temp_for <- function(site) {
    if (is.null(temperature)) return(NULL)
    if (!"site_id" %in% names(temperature)) return(temp_global)
    .ah_temperature(temperature, site, g)
  }

  ## ── Reference point (built once) ─────────────────────────────────────────
  Rp <- NULL
  Rfill <- NULL
  ref_sample_days <- numeric(0)
  if (!is.null(red_r)) {
    ref_site <- unique(as.character(reference$site_id))
    rsite <- red_r[red_r$site_id == ref_site, , drop = FALSE]
    Rp <- .ah_build_point(rsite, ssd, g, span_cap, 0.5, temp_for(ref_site),
                          co_extra)
    rstat <- .ah_static_summary(rsite, tox)
    Rfill <- Rp$norm
    for (k in seq_along(tox)) {
      Rfill[, k] <- if (is.na(rstat[k])) NA_real_ else
        .ah_normalise(ssd, k, rep(rstat[k], N), Rp$co)
    }
    ref_sample_days <- unique(rsite$date[rsite$analyte %in% ssd_all$analyte])
  }

  ## ── Focal dates per site ─────────────────────────────────────────────────
  focal_for <- function(site) {
    if (is.null(focal_dates)) {
      return(as.numeric(expand_focal_dates(grid[1L], grid[N], by = focal_by)))
    }
    if (is.data.frame(focal_dates)) {
      fd <- focal_dates
      if ("site_id" %in% names(fd)) fd <- fd[fd$site_id == site, , drop = FALSE]
      return(sort(unique(as.numeric(as.Date(fd$focal_date)))))
    }
    sort(unique(as.numeric(as.Date(focal_dates))))
  }

  w <- exp(-(0:window) / tau)
  w_all <- sum(w)

  chronic_out <- list()
  daily_out <- list()
  ctb_out <- list()
  dctb_out <- list()
  sdays_out <- list()

  for (s in sites) {
    tsite <- red_t[red_t$site_id == s, , drop = FALSE]
    Tp <- .ah_build_point(tsite, ssd, g, span_cap, 1, temp_for(s), co_extra)
    Tn <- Tp$norm
    bd <- !is.na(Tn) & !Tp$det
    bd[is.na(bd)] <- FALSE
    if (!is.null(Rp)) {
      Rn <- Rp$norm
      need <- !is.na(Tn) & is.na(Rn)
      Rn[need] <- Rfill[need]
      Tn[bd] <- Rn[bd]
      Tn[is.na(Rn)] <- NA_real_
      Rn[is.na(Tn)] <- NA_real_
    } else {
      Tn[bd] <- 0
      Rn <- Tn
      Rn[!is.na(Rn)] <- 0
    }
    msT <- .ah_mspaf_rows(Tn, ssd)
    if (!is.null(Rp)) {
      msR <- .ah_mspaf_rows(Rn, ssd)
      D <- ifelse(is.na(msT) | is.na(msR), NA_real_,
                  ifelse(msR >= 1 - 1e-12, 0,
                         pmax((msT - msR) / (1 - msR), 0)))
    } else {
      msR <- ifelse(is.na(msT), NA_real_, 0)
      D <- msT
    }
    n_an <- as.integer(rowSums(!is.na(Tn)))
    cov <- !is.na(D)

    daily_out[[s]] <- tibble::tibble(
      date = grid, site_id = s, mspaf = 100 * D, mspaf_total = 100 * msT,
      mspaf_reference = 100 * msR, n_analytes = n_an
    )

    ## Per-analyte apportioning of each day's D (and msT).
    cA <- cT <- NULL
    if (contributions) {
      pT <- .ah_paf_cells(Tn, ssd)
      pR <- .ah_paf_cells(Rn, ssd)
      added <- pmax(pT - pR, 0)
      sa <- rowSums(added, na.rm = TRUE)
      cA <- added * ifelse(sa > 0, D / sa, 0)
      cA[is.na(cA)] <- 0
      if (also_total) {
        sT <- rowSums(pT, na.rm = TRUE)
        cT <- pT * ifelse(sT > 0, msT / sT, 0)
        cT[is.na(cT)] <- 0
      }
      if (!chronic) {
        keep <- which(!is.na(Tn), arr.ind = TRUE)
        mk <- function(C, variant) {
          x <- tibble::tibble(
            site_id = s, date = grid[keep[, 1L]], analyte = tox[keep[, 2L]],
            group = ssd$moa_group[keep[, 2L]], contribution = 100 * C[keep]
          )
          if (also_total) x <- tibble::add_column(x, variant = variant, .before = 1L)
          x
        }
        dctb_out[[s]] <- dplyr::bind_rows(mk(cA, "ara"),
                                          if (also_total) mk(cT, "total"))
      }
    }

    ## Target sample days (any toxicant result) and the pairing flag.
    sdays <- sort(unique(tsite$date[tsite$analyte %in% ssd_all$analyte]))
    sdays_out[[s]] <- tibble::tibble(
      site_id = s, date = as.Date(sdays, origin = "1970-01-01"),
      paired = sdays %in% ref_sample_days
    )

    if (!chronic) next

    ## ── Time-weighting at the focal dates ───────────────────────────────────
    fd <- focal_for(s)
    k <- as.integer(fd - start) + 1L
    covn <- as.numeric(cov)
    D0 <- ifelse(cov, D, 0)
    T0 <- ifelse(cov, msT, 0)
    is_sd <- as.numeric(cov & (g %in% sdays))
    sums <- .ah_window_sums(cbind(covn, D0 * covn, T0 * covn), k, w)
    cnt <- .ah_window_sums(cbind(covn, is_sd), k, rep(1, window + 1L))
    den <- sums[, 1L]
    coverage <- pmin(pmax(den / w_all, 0), 1)
    ok <- den > 0 & coverage >= min_coverage
    value <- ifelse(ok, 100 * sums[, 2L] / den, NA_real_)
    n_samp <- findInterval(fd, sdays) -
      findInterval(fd - window, sdays, left.open = TRUE)
    fdate <- as.Date(fd, origin = "1970-01-01")
    ch <- tibble::tibble(
      focal_date = fdate, site_id = s,
      sample_id = paste0("chronic_", fdate, "_", s),
      analyte = "msPAF", value = value
    )
    if (also_total) ch$value_total <- ifelse(ok, 100 * sums[, 3L] / den, NA_real_)
    ch$detected <- TRUE
    ch$n_samples_in_window <- as.integer(n_samp)
    ch$n_imputed_in_window <- as.integer(round(cnt[, 1L] - cnt[, 2L]))
    ch$coverage <- coverage
    chronic_out[[s]] <- ch

    if (contributions) {
      one_variant <- function(C, variant) {
        cs <- .ah_window_sums(C * covn, k, w)
        cs <- 100 * cs / den
        cs[!ok, ] <- NA_real_
        x <- tibble::tibble(
          site_id = s,
          focal_date = rep(fdate, times = length(tox)),
          analyte = rep(tox, each = length(fd)),
          group = rep(ssd$moa_group, each = length(fd)),
          contribution = as.numeric(cs)
        )
        if (also_total) x <- tibble::add_column(x, variant = variant, .before = 1L)
        x
      }
      ctb_out[[s]] <- dplyr::bind_rows(one_variant(cA, "ara"),
                                       if (also_total) one_variant(cT, "total"))
    }
  }

  daily <- dplyr::bind_rows(daily_out)
  sample_days <- dplyr::bind_rows(sdays_out)
  if (!length(sdays_out)) {
    sample_days <- tibble::tibble(site_id = character(), date = as.Date(character()),
                                  paired = logical())
  }
  if (!chronic) {
    if (contributions) {
      attr(daily, "contributions") <- .ah_arrange_ctb(dplyr::bind_rows(dctb_out),
                                                      "date")
    }
    attr(daily, "sample_days") <- sample_days
    return(daily)
  }
  out <- dplyr::bind_rows(chronic_out)
  if (nrow(out)) out <- out[order(out$site_id, out$focal_date), , drop = FALSE]
  if (contributions) {
    attr(out, "contributions") <- .ah_arrange_ctb(dplyr::bind_rows(ctb_out),
                                                  "focal_date")
  }
  attr(out, "daily") <- daily
  attr(out, "sample_days") <- sample_days
  out
}

## Order a contributions frame by (variant,) site, date, analyte.
.ah_arrange_ctb <- function(x, date_col) {
  if (!nrow(x)) return(x)
  keys <- list(x$site_id, x[[date_col]], x$analyte)
  if ("variant" %in% names(x)) keys <- c(list(x$variant), keys)
  x[do.call(order, keys), , drop = FALSE]
}
