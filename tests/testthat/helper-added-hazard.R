## Test helpers for the chronic added-hazard msPAF path
## (mspaf_pipeline(ara_form = "added_hazard"), interpolate_loglinear()).
##
## Everything here is SYNTHETIC.  The ORACLE below is a compact, independent
## re-statement of the specification maths (it does not call any new package
## code); the only package internals it uses are derive_ssd_params() (hc50,
## sigma, moa_group, parsed normalisation formula) and analyte_aliases().
##
## Conventions: days are integer offsets from `.ah_base`; concentrations in
## ug/L; below-detection rows carry detected = FALSE and value = LOR.

.ah_base <- as.Date("2024-01-01")
.ah_day <- function(d) .ah_base + d

## ── SSD parameters ───────────────────────────────────────────────────────────

.ah_ssd <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      meta  <- hydroSense:::.load_analyte_metadata(NULL)
      cache <<- hydroSense:::derive_ssd_params(meta, "multi", NULL)
    }
    cache
  }
})
.ah_par <- function(a, what) {
  p <- .ah_ssd()
  p[[what]][match(a, p$analyte)]
}
.ah_hc50  <- function(a) .ah_par(a, "hc50")
.ah_sigma <- function(a) .ah_par(a, "sigma")
.ah_moa   <- function(a) .ah_par(a, "moa_group")

## Single-substance PAF (proportion): 0 for C <= 0.
.ah_paf1 <- function(C, a) {
  ifelse(is.finite(C) & C > 0,
         stats::pnorm(log10(pmax(C, 1e-300) / .ah_hc50(a)) / .ah_sigma(a)), 0)
}

## msPAF (proportion) of a named vector of normalised concentrations (one
## day): CA within mode-of-action group over components with C > 0, IA across.
.ah_ms1 <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  surv <- 1
  for (g in unique(.ah_moa(names(x)))) {
    xx <- x[.ah_moa(names(x)) == g]
    xx <- xx[xx > 0]
    if (!length(xx)) next
    tu <- sum(xx / .ah_hc50(names(xx)))
    surv <- surv * (1 - stats::pnorm(log10(tu) / mean(.ah_sigma(names(xx)))))
  }
  1 - surv
}

## Added hazard from msT, msR (proportions).
.ah_D <- function(msT, msR) {
  ifelse(is.na(msT) | is.na(msR), NA_real_,
         ifelse(msR >= 1 - 1e-12, 0, pmax((msT - msR) / (1 - msR), 0)))
}

## ── data builders ────────────────────────────────────────────────────────────

## spec: named list; each element is a data.frame(day, value, detected) or a
## list(day =, value =, detected =).  Returns the long input frame.
.ah_rows <- function(site, spec, base = .ah_base) {
  tox <- c("Cu", "Zn", "NH3-N", "B", "Cr")
  out <- lapply(names(spec), function(a) {
    s <- spec[[a]]
    n <- length(s$day)
    det <- if (is.null(s$detected)) rep(TRUE, n) else rep_len(s$detected, n)
    tibble::tibble(
      sample_id = sprintf("%s_%s_%s", site, a, s$day),
      site_id = site, datetime = base + s$day, analyte = a,
      value = rep_len(s$value, n), detected = det,
      units.analyte = if (a %in% tox) "ug/L" else NA_character_
    )
  })
  dplyr::bind_rows(out)
}

## Constant-concentration rows: every analyte at `vals[[a]]` on every day.
.ah_const <- function(site, days, vals, base = .ah_base, detected = list()) {
  spec <- lapply(names(vals), function(a) {
    list(day = days, value = vals[[a]],
         detected = if (is.null(detected[[a]])) TRUE else detected[[a]])
  })
  names(spec) <- names(vals)
  .ah_rows(site, spec, base)
}

## ── the main scenario (two targets + one reference) ──────────────────────────
##
## T1: monthly then a 150-day target gap (uncovered), BDL Cu/Zn, measured
##     hardness on some dates and Ca+Mg on others, a co-analyte-only date
##     (day 200), an analyte (Cr) the reference never measures.
## T2: leading BDL anchor, duplicated rows on one day, alias analyte names,
##     no temperature (so NH3-N is dropped), starts before the reference.
## REF: a 210-day gap (days 90..300) so the static gap-filler is used.

.ah_scn_T1 <- function() {
  d <- c(0, 30, 60, 90, 150, 300, 330, 360)
  tox <- .ah_rows("T1", list(
    Cu = list(day = d, value = c(3, 0.5, 0.5, 5, 8, 12, 6, 0.5),
              detected = c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)),
    Zn = list(day = d, value = c(30, 55, 80, 60, 20, 90, 70, 100),
              detected = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)),
    B = list(day = d, value = c(3000, 4000, 9000, 7000, 5000, 12000, 8000, 6000)),
    `NH3-N` = list(day = d, value = c(400, 900, 300, 1500, 600, 2500, 800, 700)),
    Cr = list(day = d, value = c(50, 80, 120, 60, 200, 90, 70, 40)),
    pH = list(day = c(d, 200), value = c(7.1, 7.4, 6.9, 7.6, 7.2, 7.0, 7.5, 7.3, 7.2)),
    DOC = list(day = c(d, 200), value = c(1.2, 2.5, 4, 0.8, 3.3, 1.9, 2.2, 5, 2.8)),
    temperature = list(day = d, value = c(12, 15, 18, 22, 16, 14, 11, 13)),
    hardness = list(day = c(0, 30, 90, 300), value = c(40, 55, 70, 35)),
    Ca = list(day = c(0, 60, 150, 330, 360), value = c(30, 12, 15, 10, 11)),
    Mg = list(day = c(0, 60, 150, 330, 360), value = c(20, 5, 6, 4, 4.5))
  ))
  tox
}

.ah_scn_T2 <- function() {
  d <- c(-20, 30, 45, 100)
  main <- .ah_rows("T2", list(
    Cu = list(day = c(-20, 30, 45, 100), value = c(0.5, 10, 0.5, 4),
              detected = c(FALSE, TRUE, FALSE, TRUE)),
    Zn = list(day = d, value = c(25, 40, 15, 65)),
    B = list(day = d, value = c(2500, 5000, 3000, 8000)),
    `NH3-N` = list(day = d, value = c(300, 700, 500, 900)),
    pH = list(day = d, value = c(7.0, 7.3, 7.1, 7.4)),
    DOC = list(day = d, value = c(2, 3, 1.5, 2.5)),
    `Hardness-total-CaCO3` = list(day = d, value = c(35, 60, 45, 65))
  ))
  extra <- .ah_rows("T2", list(
    Cu = list(day = c(30, 45, 100), value = c(11, 0.3, 0.5),
              detected = c(TRUE, FALSE, FALSE))
  ))
  ## day 30: two detected rows (mean 10.5); day 45: two BDL rows (min LOR 0.3);
  ## day 100: one detected + one BDL row (detected mean 4).
  extra$sample_id <- paste0(extra$sample_id, "_dup")
  dplyr::bind_rows(main, extra)
}

.ah_scn_REF <- function() {
  d <- c(0, 30, 60, 90, 300, 330, 360)
  .ah_rows("REF", list(
    Cu = list(day = d, value = c(0.8, 0.6, 0.8, 1.0, 0.7, 0.9, 0.6),
              detected = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)),
    Zn = list(day = d, value = c(5, 8, 6, 10, 4, 7, 9)),
    B = list(day = d, value = c(500, 600, 400, 700, 550, 450, 650)),
    `NH3-N` = list(day = d, value = c(100, 150, 90, 200, 120, 110, 130)),
    pH = list(day = d, value = c(7.0, 7.2, 7.1, 7.3, 6.9, 7.2, 7.0)),
    DOC = list(day = d, value = c(1, 1.5, 2, 1.2, 0.9, 1.8, 1.1)),
    temperature = list(day = d, value = c(11, 14, 17, 21, 15, 13, 10)),
    hardness = list(day = d, value = c(30, 35, 40, 45, 32, 38, 36))
  ))
}

.ah_scn_target <- function() dplyr::bind_rows(.ah_scn_T1(), .ah_scn_T2())

## ── ORACLE ───────────────────────────────────────────────────────────────────

.ah_canon <- function(d) {
  al <- analyte_aliases()
  hit <- match(d$analyte, al$alias)
  d$analyte <- ifelse(is.na(hit), d$analyte, al$canonical[hit])
  d
}

## One value per (site, analyte, date).
.ah_reduce <- function(d) {
  d <- .ah_canon(as.data.frame(d))
  d$date <- as.Date(d$datetime)
  key <- paste(d$site_id, d$analyte, as.numeric(d$date), sep = "|")
  parts <- lapply(split(seq_len(nrow(d)), key), function(i) {
    x <- d[i, , drop = FALSE]
    if (any(x$detected)) {
      v <- mean(x$value[x$detected]); q <- TRUE
    } else {
      v <- min(x$value); q <- FALSE
    }
    data.frame(site_id = x$site_id[1], analyte = x$analyte[1], date = x$date[1],
               value = v, detected = q, stringsAsFactors = FALSE)
  })
  do.call(rbind, unname(parts))
}

## Oracle interpolation (loop form, independent of the package).
.ah_interp <- function(dates, v, det, grid, cap = 120, log = TRUE) {
  n <- length(grid)
  val <- rep(NA_real_, n); dt <- rep(NA, n)
  ok <- !is.na(v) & (!log | v > 0)
  if (any(ok)) {
    o <- order(dates[ok])
    x <- as.numeric(dates[ok])[o]; y <- v[ok][o]; q <- det[ok][o]
    for (k in seq_len(n)) {
      g <- as.numeric(grid[k]); i <- sum(x <= g)
      if (i == 0L) next
      if (x[i] == g) { val[k] <- y[i]; dt[k] <- q[i]; next }
      if (i == length(x)) {
        if (g - x[i] <= cap) { val[k] <- y[i]; dt[k] <- q[i] }
        next
      }
      j <- i + 1L
      if (x[j] - x[i] <= cap) {
        f <- (g - x[i]) / (x[j] - x[i])
        val[k] <- if (log) exp((1 - f) * log(y[i]) + f * log(y[j]))
                  else (1 - f) * y[i] + f * y[j]
        dt[k] <- q[i] || q[j]
      }
    }
  }
  data.frame(date = grid, value = val, detected = dt)
}

## Co-analyte daily series for one point (no cap, constant beyond the ends).
.ah_co <- function(red, site, grid) {
  x <- red[red$site_id == site, , drop = FALSE]
  g <- as.numeric(grid)
  anchors <- function(df) df[order(df$date), , drop = FALSE]
  series <- function(df, log) {
    if (!nrow(df)) return(rep(NA_real_, length(g)))
    df <- anchors(df)
    if (nrow(df) == 1L) return(rep(df$value, length(g)))
    y <- if (log) log(df$value) else df$value
    r <- stats::approx(as.numeric(df$date), y, xout = g, rule = 2)$y
    if (log) exp(r) else r
  }
  one <- function(a, log) {
    series(x[x$analyte == a & x$detected & x$value > 0, , drop = FALSE], log)
  }
  hd <- x[x$analyte == "hardness" & x$detected & x$value > 0,
          c("date", "value"), drop = FALSE]
  ca <- x[x$analyte == "Ca" & x$detected & x$value > 0, c("date", "value")]
  mg <- x[x$analyte == "Mg" & x$detected & x$value > 0, c("date", "value")]
  both <- intersect(as.numeric(ca$date), as.numeric(mg$date))
  both <- setdiff(both, as.numeric(hd$date))
  hc <- if (length(both)) {
    data.frame(date = as.Date(both, origin = "1970-01-01"),
               value = 2.497 * ca$value[match(both, as.numeric(ca$date))] +
                 4.118 * mg$value[match(both, as.numeric(mg$date))])
  } else hd[0, ]
  hh <- rbind(hd, hc)
  data.frame(pH = one("pH", FALSE), temperature = one("temperature", FALSE),
             DOC = one("DOC", TRUE), hardness = series(hh, TRUE),
             Ca = one("Ca", TRUE), Mg = one("Mg", TRUE))
}

.ah_norm <- function(a, C, co) {
  f <- .ah_par(a, "parsed_formula")[[1]]
  if (is.null(f)) return(C)
  env <- list2env(as.list(co), parent = baseenv())
  env$C <- C
  out <- eval(f, env)
  if (length(out) == 1L) out <- rep(out, length(C))
  out
}

## Point build: raw / detected / normalised matrices (days x toxicants).
.ah_point <- function(red, site, tox, grid, cap, bdl_mult, temp = NULL) {
  co <- .ah_co(red, site, grid)
  if (!is.null(temp)) co$temperature <- temp
  raw <- det <- matrix(NA_real_, length(grid), length(tox),
                       dimnames = list(NULL, tox))
  for (a in tox) {
    o <- red[red$site_id == site & red$analyte == a, , drop = FALSE]
    if (!nrow(o)) next
    s <- .ah_interp(o$date, ifelse(o$detected, o$value, bdl_mult * o$value),
                    o$detected, grid, cap, log = TRUE)
    raw[, a] <- s$value; det[, a] <- as.numeric(s$detected)
  }
  nrm <- raw
  for (a in tox) nrm[, a] <- .ah_norm(a, raw[, a], co)
  list(raw = raw, det = det, norm = nrm, co = co)
}

## Static reference summary (raw ug/L) for the gap-filler.
.ah_static <- function(red, ref, tox) {
  vapply(tox, function(a) {
    o <- red[red$site_id == ref & red$analyte == a, , drop = FALSE]
    if (!nrow(o)) return(NA_real_)
    v <- ifelse(o$detected, o$value, o$value / 2)
    gm <- exp(mean(log(v)))
    if (all(o$detected)) gm else max(gm, stats::median(o$value[!o$detected]) / 2)
  }, numeric(1))
}

## msPAF per row of a (days x analytes) matrix; NA row when no analyte present.
.ah_ms_mat <- function(M) {
  apply(M, 1, function(r) .ah_ms1(r))
}

## Time-weighting of one daily series x (NA = uncovered) at grid index i.
.ah_tw <- function(i, tau = 90, window = 365) {
  a <- 0:window
  list(idx = i - a, w = exp(-a / tau))
}

#' Oracle for mspaf_pipeline(ara_form = "added_hazard").
#'
#' target / reference: long frames (site_id, datetime, analyte, value,
#' detected).  Returns list(chronic, daily, contributions, sample_days).
.ah_oracle <- function(target, reference = NULL, span_cap = 120, tau = 90,
                       window = 365, min_coverage = 0.2, focal_dates = NULL,
                       grid = NULL, temp = NULL) {
  tg <- .ah_reduce(target)
  rf <- if (is.null(reference)) NULL else .ah_reduce(reference)
  red <- rbind(tg, rf)
  tox_all <- .ah_ssd()$analyte
  tox <- intersect(tox_all, unique(red$analyte))
  sites <- unique(tg$site_id)
  ref_site <- if (is.null(rf)) NULL else unique(rf$site_id)
  allraw <- .ah_canon(as.data.frame(rbind(target, reference)))
  if (is.null(grid)) {
    lo <- min(as.Date(allraw$datetime))
    hi <- max(as.Date(target$datetime))
    hi <- if (is.finite(span_cap)) hi + span_cap else hi
    grid <- seq(lo, hi, by = "day")
  }
  ng <- length(grid)
  if (is.null(focal_dates)) focal_dates <- grid

  R <- NULL; Rstat <- NULL
  if (!is.null(rf)) {
    R <- .ah_point(red, ref_site, tox, grid, span_cap, 0.5, temp)
    Rstat <- .ah_static(red, ref_site, tox)
  }
  wfull <- exp(-(0:window) / tau)
  out <- list(chronic = list(), daily = list(), contributions = list(),
              sample_days = list())

  for (s in sites) {
    Tp <- .ah_point(red, s, tox, grid, span_cap, 1, temp)
    Tn <- Tp$norm
    if (!is.null(R)) {
      Rn <- R$norm
      for (a in tox) {
        need <- !is.na(Tn[, a]) & is.na(Rn[, a])
        if (any(need)) {
          Rn[need, a] <- if (is.na(Rstat[a])) NA_real_ else
            .ah_norm(a, rep(Rstat[a], sum(need)), R$co[need, , drop = FALSE])
        }
      }
      bd <- Tp$det[, tox, drop = FALSE] == 0 & !is.na(Tn)
      bd[is.na(bd)] <- FALSE
      T0 <- Tn; T0[bd] <- 0          # msT0: target alone, BDL days at 0
      Tn[bd] <- Rn[bd]
      Tn[is.na(Rn)] <- NA; Rn[is.na(Tn)] <- NA
      T0[is.na(Tn)] <- NA
    } else {
      Rn <- Tn; Rn[!is.na(Rn)] <- 0
      bd <- Tp$det[, tox, drop = FALSE] == 0 & !is.na(Tn)
      bd[is.na(bd)] <- FALSE
      Tn[bd] <- 0
      T0 <- Tn
    }
    msT <- .ah_ms_mat(Tn)
    msT0 <- .ah_ms_mat(T0)
    msR <- if (is.null(R)) ifelse(is.na(msT), NA_real_, 0) else .ah_ms_mat(Rn)
    D <- if (is.null(R)) msT else .ah_D(msT, msR)
    nan <- rowSums(!is.na(Tn))

    ## per-analyte apportioning
    pT <- pR <- Tn; pT[] <- NA_real_; pR[] <- NA_real_
    for (a in tox) { pT[, a] <- .ah_paf1(Tn[, a], a); pR[, a] <- .ah_paf1(Rn[, a], a) }
    pT[is.na(Tn)] <- NA; pR[is.na(Tn)] <- NA
    added <- pmax(pT - pR, 0)
    sa <- rowSums(added, na.rm = TRUE)
    pT0 <- T0; pT0[] <- NA_real_
    for (a in tox) pT0[, a] <- .ah_paf1(T0[, a], a)
    pT0[is.na(T0)] <- NA
    sT <- rowSums(pT0, na.rm = TRUE)
    cA <- added * ifelse(sa > 0, D / sa, 0)
    cA[is.na(cA)] <- 0
    cT <- pT0 * ifelse(sT > 0, msT0 / sT, 0)
    cT[is.na(cT)] <- 0

    ## time-weighting
    tgd <- tg[tg$site_id == s & tg$analyte %in% tox_all, ]
    sdates <- sort(unique(tgd$date))
    nf <- length(focal_dates)
    ch <- data.frame(focal_date = focal_dates, site_id = s, value = NA_real_,
                     value_total = NA_real_, coverage = NA_real_,
                     n_samples_in_window = NA_integer_,
                     n_imputed_in_window = NA_integer_)
    ca_rows <- list()
    for (k in seq_len(nf)) {
      f <- focal_dates[k]
      i <- match(as.numeric(f), as.numeric(grid))
      a <- 0:window
      idx <- i - a
      okidx <- idx >= 1L
      cov <- rep(FALSE, length(a)); cov[okidx] <- !is.na(D[idx[okidx]])
      den <- sum(wfull[cov])
      cvg <- den / sum(wfull)
      ch$coverage[k] <- cvg
      ch$n_samples_in_window[k] <- sum(sdates >= f - window & sdates <= f)
      cd <- grid[idx[okidx]][cov[okidx]]
      ch$n_imputed_in_window[k] <- sum(!(cd %in% sdates))
      if (den > 0 && cvg >= min_coverage) {
        wi <- wfull[cov]; ii <- idx[cov]
        ch$value[k] <- 100 * sum(wi * D[ii]) / den
        ch$value_total[k] <- 100 * sum(wi * msT0[ii]) / den
        cAi <- 100 * colSums(cA[ii, , drop = FALSE] * wi) / den
        cTi <- 100 * colSums(cT[ii, , drop = FALSE] * wi) / den
      } else {
        cAi <- cTi <- rep(NA_real_, length(tox)); names(cAi) <- names(cTi) <- tox
      }
      ca_rows[[k]] <- rbind(
        data.frame(variant = "ara", site_id = s, focal_date = f, analyte = tox,
                   group = .ah_moa(tox), contribution = unname(cAi)),
        data.frame(variant = "total", site_id = s, focal_date = f, analyte = tox,
                   group = .ah_moa(tox), contribution = unname(cTi)))
    }
    out$chronic[[s]] <- ch
    out$contributions[[s]] <- do.call(rbind, ca_rows)
    out$daily[[s]] <- data.frame(
      date = grid, site_id = s, mspaf = 100 * D, mspaf_total = 100 * msT0,
      mspaf_reference = 100 * msR, n_analytes = nan)
    ## pairing flag
    tdates <- sort(unique(tgd$date))
    rdates <- if (is.null(rf)) as.Date(character()) else
      unique(rf$date[rf$analyte %in% tox_all])
    out$sample_days[[s]] <- data.frame(site_id = s, date = tdates,
                                       paired = tdates %in% rdates)
  }
  lapply(out, function(z) { r <- do.call(rbind, unname(z)); rownames(r) <- NULL; r })
}

## ── comparison helpers ───────────────────────────────────────────────────────

## Merge two frames on `by`, returning the joined frame with x/y suffixes.
.ah_join <- function(got, exp, by) {
  got <- as.data.frame(got); exp <- as.data.frame(exp)
  merge(got, exp, by = by, suffixes = c(".got", ".exp"), all = TRUE)
}

## Expect two numeric vectors equal (NA positions identical) within `tol`.
.ah_expect_close <- function(got, exp, tol = 1e-6, label = NULL) {
  testthat::expect_identical(is.na(got), is.na(exp), label = label)
  ok <- !is.na(got) & !is.na(exp)
  testthat::expect_true(all(abs(got[ok] - exp[ok]) <= tol),
                        label = paste0(if (is.null(label)) "values" else label,
                                       " within ", format(tol)))
}
