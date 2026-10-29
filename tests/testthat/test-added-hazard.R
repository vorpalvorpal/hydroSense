## Chronic added-hazard msPAF: mspaf_pipeline(ara_form = "added_hazard")
## ============================================================================
##
## CONTRACT (condensed from the specification; all data synthetic).
##
## mspaf_pipeline(..., ara_form = c("excess", "added_hazard"), min_coverage = 0.2)
##
## * "excess" (default) is the existing behaviour, untouched.
## * "added_hazard": every point (each target, and the reference) is
##   interpolated from its own samples onto a daily grid (log-linear for
##   toxicants, span cap `daily_args$span_cap`, default 120 d), detection-limit
##   rules are applied, and per day
##       D = max((msPAF_T - msPAF_R) / (1 - msPAF_R), 0)
##   is time-weighted with w_a = exp(-a / tau), a = 0..window (tau 90,
##   window 365), renormalised over covered days.
##   Argument errors (all before any model is fitted):
##     daily_args$interpolation other than "loglinear" (default), or
##       "loglinear" with ara_form = "excess"  -> hydroSense_ara_conflict
##     ara = "static"                          -> hydroSense_ara_conflict
##     impute = TRUE                           -> hydroSense_ara_conflict
##     reference with != 1 site_id             -> hydroSense_reference_sites
##   `ara` unspecified: "contemporaneous" with a reference, else "none".
## * Below detection: reference anchors at LOR/2, target anchors at LOR; a
##   target below detection on a day gives T := R (zero excess); the static
##   gap-filler s = max(geomean with BDL at LOR/2, median(LOR)/2); an analyte
##   the reference never measured is dropped.  With no reference, target BDL
##   days enter at concentration 0 and D = msPAF_T.
## * Output (chronic = TRUE): one row per (site, focal date) with focal_date,
##   site_id, sample_id, analyte, value, [value_total], detected,
##   n_samples_in_window, n_imputed_in_window, coverage; value is NA when
##   coverage < min_coverage.  attr "daily" (date, site_id, mspaf [=100 D],
##   mspaf_total, mspaf_reference, n_analytes), attr "sample_days" (site_id,
##   date, paired), attr "contributions" (see roxygen: apportioned by added
##   single-substance PAF; sums to value / value_total).
##
## Toxicants used: Cu (DOC-normalised), Zn (pH/hardness/DOC), NH3-N
## (pH/temperature), B, Cr: all have SSD parameters from the shipped metadata,
## so none of these tests needs the (gitignored) guideline data.

library(testthat)

.ah_cache <- new.env()
.ah_get <- function(key, f) {
  if (is.null(.ah_cache[[key]])) .ah_cache[[key]] <- f()
  .ah_cache[[key]]
}

## Run the pipeline on the added-hazard path.
.ah_run <- function(target, reference = NULL, ..., daily_args = list()) {
  args <- list(target = target, reference = reference, impute = FALSE,
               ara_form = "added_hazard", daily_args = daily_args)
  args <- utils::modifyList(args, list(...), keep.null = TRUE)
  suppressWarnings(suppressMessages(do.call(mspaf_pipeline, args)))
}

## Make any model fit an error, to prove validation happens first.
.ah_block_fits <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    fit_reference_model = function(...) stop("fit_reference_model was called"),
    fit_imputation_model = function(...) stop("fit_imputation_model was called"),
    .package = "hydroSense", .env = env
  )
}

## ── scenario runs (memoised) ─────────────────────────────────────────────────

.ah_main <- function() .ah_get("main", function() {
  .ah_run(.ah_scn_target(), .ah_scn_REF(), also_total = TRUE,
          contributions = TRUE)
})
.ah_main_oracle <- function() .ah_get("main_oracle", function() {
  .ah_oracle(.ah_scn_target(), .ah_scn_REF())
})

## Constant two-toxicant scenario (hand-computable): Cu in the ionoregulatory
## group, B in the boron group; DOC 0.5 makes the Cu normalisation the identity.
.ah_sb_days <- seq(0, 390, by = 30)
.ah_sb <- function(tvals = c(Cu = 10, B = 8000, DOC = 0.5),
                   rvals = c(Cu = 2, B = 2000, DOC = 0.5),
                   tdet = list(), rdet = list(), days = .ah_sb_days,
                   rdays = days) {
  list(target = .ah_const("T1", days, as.list(tvals), detected = tdet),
       reference = .ah_const("REF", rdays, as.list(rvals), detected = rdet))
}
.ah_sb_focal <- .ah_day(390)

## ═════════════════════════════════════════════════════════════════════════════
## Oracle self-checks (these do not touch the new package code)
## ═════════════════════════════════════════════════════════════════════════════

describe("test oracle: hand-computed numbers", {

  it("interpolates log-linearly, caps spans and holds after the last sample", {
    g <- .ah_day(c(-1, 0, 5, 10, 130, 131))
    o <- .ah_interp(.ah_day(c(0, 10)), c(2, 8), c(TRUE, TRUE), g)
    expect_equal(o$value, c(NA, 2, 4, 8, 8, NA), tolerance = 1e-12)
    o <- .ah_interp(.ah_day(c(0, 121)), c(2, 8), c(TRUE, TRUE), .ah_day(60))
    expect_true(is.na(o$value))
    o <- .ah_interp(.ah_day(c(0, 10)), c(2, 8), c(FALSE, TRUE), .ah_day(5))
    expect_identical(o$detected, TRUE)
  })

  it("matches the single-substance PAF and the CA/IA combination by hand", {
    p <- .ah_ssd()
    hc <- p$hc50[p$analyte == "Cu"]; sg <- p$sigma[p$analyte == "Cu"]
    expect_equal(.ah_paf1(10, "Cu"), stats::pnorm(log10(10 / hc) / sg))
    pb <- .ah_paf1(8000, "B")
    expect_equal(.ah_ms1(c(Cu = 10, B = 8000)), 1 - (1 - .ah_paf1(10, "Cu")) * (1 - pb))
    expect_equal(.ah_ms1(c(Cu = 10, B = NA)), .ah_paf1(10, "Cu"))
    expect_true(is.na(.ah_ms1(c(Cu = NA_real_))))
    expect_equal(.ah_D(0.5, 0.2), 0.375)
    expect_equal(.ah_D(0.2, 0.5), 0)
  })

  it("reduces duplicate rows and applies the static summary rule", {
    d <- data.frame(site_id = "S", datetime = .ah_day(c(0, 0, 1, 1, 2, 2)),
                    analyte = "Cu", value = c(2, 4, 0.5, 0.3, 4, 0.5),
                    detected = c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE))
    r <- .ah_reduce(d)
    r <- r[order(r$date), ]
    expect_equal(r$value, c(3, 0.3, 4))
    expect_identical(r$detected, c(TRUE, FALSE, TRUE))
    ref <- data.frame(site_id = "REF", analyte = "Cu", date = .ah_day(0:2),
                      value = c(2, 2, 2), detected = FALSE)
    expect_equal(.ah_static(ref, "REF", "Cu"), c(Cu = 1))
    ref$detected <- c(TRUE, FALSE, FALSE); ref$value <- c(8, 2, 2)
    expect_equal(.ah_static(ref, "REF", "Cu"),
                 c(Cu = max(exp(mean(log(c(8, 1, 1)))), 1)))
  })

  it("oracle chronic value equals the hand kernel on a constant series", {
    s <- .ah_sb()
    o <- .ah_oracle(s$target, s$reference, focal_dates = .ah_sb_focal)
    msT <- .ah_ms1(c(Cu = 10, B = 8000)); msR <- .ah_ms1(c(Cu = 2, B = 2000))
    D <- (msT - msR) / (1 - msR)
    expect_equal(o$chronic$value, 100 * D, tolerance = 1e-10)
    expect_equal(o$chronic$value_total, 100 * msT, tolerance = 1e-10)
    expect_equal(o$chronic$coverage, 1)
    ct <- o$contributions[o$contributions$variant == "ara", ]
    expect_equal(sum(ct$contribution), 100 * D, tolerance = 1e-10)
  })

  it("oracle exercises every branch in the main scenario", {
    o <- .ah_main_oracle()
    expect_setequal(unique(o$chronic$site_id), c("T1", "T2"))
    d1 <- o$daily[o$daily$site_id == "T1", ]
    expect_true(any(!is.na(d1$mspaf)))
    expect_true(any(is.na(d1$mspaf)))              # target gap 150..300 uncovered
    expect_true(any(o$chronic$coverage > 0.2 & !is.na(o$chronic$value)))
    expect_true(any(is.na(o$chronic$value)))
    expect_false(any(o$sample_days$paired[o$sample_days$site_id == "T2"][-2]))
  })
})

## ═════════════════════════════════════════════════════════════════════════════
## Argument validation (before any fitting)
## ═════════════════════════════════════════════════════════════════════════════

describe("mspaf_pipeline(ara_form =): argument validation", {
  s <- .ah_sb(days = c(0, 30, 60), rdays = c(0, 30, 60))

  it("rejects daily_args$interpolation other than 'loglinear'", {
    .ah_block_fits()
    for (m in c("linear", "model")) {
      expect_error(
        mspaf_pipeline(s$target, s$reference, impute = FALSE,
                       ara_form = "added_hazard",
                       daily_args = list(interpolation = m)),
        class = "hydroSense_ara_conflict"
      )
    }
  })

  it("rejects interpolation = 'loglinear' with ara_form = 'excess'", {
    .ah_block_fits()
    expect_error(
      mspaf_pipeline(s$target, s$reference, impute = FALSE, ara_form = "excess",
                     daily_args = list(interpolation = "loglinear")),
      class = "hydroSense_ara_conflict"
    )
    expect_error(
      mspaf_pipeline(s$target, s$reference, impute = FALSE,
                     daily_args = list(interpolation = "loglinear")),
      class = "hydroSense_ara_conflict"
    )
  })

  it("rejects ara = 'static'", {
    .ah_block_fits()
    expect_error(
      mspaf_pipeline(s$target, s$reference, impute = FALSE,
                     ara_form = "added_hazard", ara = "static",
                     daily_args = list(reference = s$reference)),
      class = "hydroSense_ara_conflict"
    )
    expect_error(
      mspaf_pipeline(s$target, s$reference, impute = FALSE,
                     ara_form = "added_hazard", ara = "static"),
      class = "hydroSense_ara_conflict"
    )
  })

  it("rejects impute = TRUE", {
    .ah_block_fits()
    expect_error(
      mspaf_pipeline(s$target, s$reference, impute = TRUE,
                     ara_form = "added_hazard"),
      class = "hydroSense_ara_conflict"
    )
  })

  it("requires a single-site reference", {
    .ah_block_fits()
    two <- dplyr::bind_rows(s$reference,
                            .ah_const("REF2", c(0, 30, 60), list(Cu = 2, DOC = 0.5)))
    expect_error(
      mspaf_pipeline(s$target, two, impute = FALSE, ara_form = "added_hazard"),
      class = "hydroSense_reference_sites"
    )
  })

  it("rejects an unknown ara_form, a non-arith_mean summary (hydroSense_ara_conflict), and ara = 'contemporaneous' without a reference", {
    .ah_block_fits()
    expect_error(
      mspaf_pipeline(s$target, s$reference, impute = FALSE, ara_form = "bogus")
    )
    expect_error(
      mspaf_pipeline(s$target, s$reference, impute = FALSE,
                     ara_form = "added_hazard", chronic_summary = "p90"),
      class = "hydroSense_ara_conflict"
    )
    expect_error(
      mspaf_pipeline(s$target, NULL, impute = FALSE, ara_form = "added_hazard",
                     ara = "contemporaneous"),
      class = "hydroSense_ara_conflict"
    )
  })
})

describe("mspaf_pipeline(ara_form = 'excess'): default behaviour is unchanged", {
  it("a call without ara_form equals one with ara_form = 'excess'", {
    skip_on_cran()
    sc <- .pf_scenario()
    tgt <- .pf_small_target()
    tgt <- tgt[tgt$site_id == "P1", , drop = FALSE]
    tmp <- sc$temperature
    tmp <- tmp[tmp$datetime >= as.Date("2021-12-15") &
                 tmp$datetime <= as.Date("2022-03-10"), ]
    base <- list(
      target = tgt, reference = NULL, impute = FALSE,
      daily_args = list(temperature = tmp, start = as.Date("2022-01-01"),
                        end = as.Date("2022-02-28")),
      focal_dates = as.Date(c("2022-01-20", "2022-02-28")),
      tau = 90, window = 365
    )
    run <- function(...) suppressWarnings(suppressMessages(
      do.call(mspaf_pipeline, c(base, list(...)))))
    a <- run()
    b <- run(ara_form = "excess")
    attributes(a)[c("daily", "reference_model", "imputation_model")] <- NULL
    attributes(b)[c("daily", "reference_model", "imputation_model")] <- NULL
    expect_equal(a, b, tolerance = 1e-12)
    expect_false("coverage" %in% names(a))
  })
})

## ═════════════════════════════════════════════════════════════════════════════
## Equivalence with the oracle
## ═════════════════════════════════════════════════════════════════════════════

describe("added_hazard: output structure", {
  skip_on_cran()

  it("has the documented columns, one row per (site, focal date), all grid days", {
    out <- as.data.frame(.ah_main())
    expect_true(all(c("focal_date", "site_id", "sample_id", "analyte", "value",
                      "value_total", "detected", "n_samples_in_window",
                      "n_imputed_in_window", "coverage") %in% names(out)))
    expect_equal(anyDuplicated(out[c("site_id", "focal_date")]), 0L)
    o <- .ah_main_oracle()$chronic
    ## focal_by = "day" (default): every grid day, every site
    expect_equal(nrow(out), nrow(o))
    expect_equal(range(out$focal_date), range(o$focal_date))
    expect_setequal(out$site_id, c("T1", "T2"))
    expect_identical(out$sample_id,
                     paste0("chronic_", out$focal_date, "_", out$site_id))
    expect_true(all(out$analyte == "msPAF"))
    expect_true(all(out$detected))
  })

  it("keeps coverage in [0, 1] and gates value at min_coverage = 0.2", {
    out <- as.data.frame(.ah_main())
    expect_true(all(out$coverage >= 0 & out$coverage <= 1))
    expect_identical(is.na(out$value), out$coverage < 0.2 | out$coverage == 0)
    expect_identical(is.na(out$value_total), is.na(out$value))
  })
})

describe("added_hazard: equivalence with the oracle", {
  skip_on_cran()

  it("matches value, value_total, coverage and window counts (1e-6)", {
    out <- as.data.frame(.ah_main())
    o <- .ah_main_oracle()$chronic
    m <- .ah_join(out, o, c("site_id", "focal_date"))
    expect_equal(nrow(m), nrow(o))
    expect_false(anyNA(m$coverage.got))
    .ah_expect_close(m$value.got, m$value.exp, label = "value")
    .ah_expect_close(m$value_total.got, m$value_total.exp, label = "value_total")
    .ah_expect_close(m$coverage.got, m$coverage.exp, label = "coverage")
    expect_equal(as.numeric(m$n_samples_in_window.got),
                 as.numeric(m$n_samples_in_window.exp))
    expect_equal(as.numeric(m$n_imputed_in_window.got),
                 as.numeric(m$n_imputed_in_window.exp))
  })

  it("matches the daily added hazard, total and reference series", {
    d <- as.data.frame(attr(.ah_main(), "daily"))
    expect_true(all(c("date", "site_id", "mspaf", "mspaf_total",
                      "mspaf_reference", "n_analytes") %in% names(d)))
    o <- .ah_main_oracle()$daily
    m <- .ah_join(d, o, c("site_id", "date"))
    ## every day the oracle covers is present in the package's frame
    cov_exp <- !is.na(m$mspaf.exp)
    expect_false(anyNA(m$mspaf.got[cov_exp]))
    ## days the package returns beyond the oracle's are uncovered
    expect_true(all(is.na(m$mspaf.got[is.na(m$mspaf.exp)])))
    .ah_expect_close(m$mspaf.got[cov_exp], m$mspaf.exp[cov_exp], label = "mspaf")
    .ah_expect_close(m$mspaf_total.got[cov_exp], m$mspaf_total.exp[cov_exp],
                     label = "mspaf_total")
    .ah_expect_close(m$mspaf_reference.got[cov_exp], m$mspaf_reference.exp[cov_exp],
                     label = "mspaf_reference")
    expect_equal(as.numeric(m$n_analytes.got[cov_exp]),
                 as.numeric(m$n_analytes.exp[cov_exp]))
  })

  it("exposes the pairing flag per target sample day", {
    sd <- as.data.frame(attr(.ah_main(), "sample_days"))
    expect_true(all(c("site_id", "date", "paired") %in% names(sd)))
    o <- .ah_main_oracle()$sample_days
    m <- .ah_join(sd, o, c("site_id", "date"))
    expect_equal(nrow(m), nrow(o))          # co-analyte-only day 200 excluded
    expect_identical(as.logical(m$paired.got), as.logical(m$paired.exp))
    ## spot checks: T1 day 150 is unpaired, day 30 paired; T2 day 30 paired
    f <- function(s, d) sd$paired[sd$site_id == s & sd$date == .ah_day(d)]
    expect_false(f("T1", 150)); expect_true(f("T1", 30)); expect_true(f("T2", 30))
    expect_false(f("T2", 45))
    expect_equal(sum(sd$date[sd$site_id == "T1"] == .ah_day(200)), 0L)
  })

  it("is independent per target site (joint call = separate calls)", {
    tg <- .ah_scn_target(); ref <- .ah_scn_REF()
    fd <- .ah_day(seq(10, 220, by = 30))
    run <- function(sites) as.data.frame(.ah_run(
      tg[tg$site_id %in% sites, ], ref, also_total = TRUE, focal_dates = fd))
    joint <- run(c("T1", "T2"))
    for (s in c("T1", "T2")) {
      alone <- run(s)
      a <- alone[order(alone$focal_date), ]
      j <- joint[joint$site_id == s, ]; j <- j[order(j$focal_date), ]
      expect_equal(nrow(j), length(fd))
      .ah_expect_close(j$value, a$value, tol = 1e-9, label = paste(s, "value"))
      .ah_expect_close(j$value_total, a$value_total, tol = 1e-9)
      .ah_expect_close(j$coverage, a$coverage, tol = 1e-9)
    }
  })
})

describe("added_hazard: no reference (ara = 'none')", {
  skip_on_cran()

  it("matches the oracle and defaults ara to 'none' without a reference", {
    t1 <- .ah_scn_T1()
    out <- as.data.frame(.ah_run(t1, NULL, contributions = TRUE))
    o <- .ah_oracle(t1, NULL)
    m <- .ah_join(out, o$chronic, c("site_id", "focal_date"))
    expect_equal(nrow(m), nrow(o$chronic))
    .ah_expect_close(m$value.got, m$value.exp, label = "value")
    .ah_expect_close(m$coverage.got, m$coverage.exp, label = "coverage")
    ctb <- as.data.frame(attr(.ah_run(t1, NULL, contributions = TRUE), "contributions"))
    s <- tapply(ctb$contribution, ctb$focal_date, sum)
    v <- out$value[match(names(s), as.character(out$focal_date))]
    ok <- !is.na(v)
    expect_true(all(abs(s[ok] - v[ok]) < 1e-6))
  })

  it("equals an explicit ara = 'none'", {
    s <- .ah_sb(days = c(0, 30, 60, 90))
    a <- .ah_run(s$target, NULL, focal_dates = .ah_day(90))
    b <- .ah_run(s$target, NULL, ara = "none", focal_dates = .ah_day(90))
    expect_equal(as.data.frame(a)$value, as.data.frame(b)$value, tolerance = 1e-12)
  })

  it("equals an explicit ara = 'contemporaneous' when a reference is given", {
    s <- .ah_sb(days = c(0, 30, 60, 90))
    a <- .ah_run(s$target, s$reference, focal_dates = .ah_day(90))
    b <- .ah_run(s$target, s$reference, ara = "contemporaneous",
                 focal_dates = .ah_day(90))
    expect_equal(as.data.frame(a)$value, as.data.frame(b)$value, tolerance = 1e-12)
  })

  it("is the target's own msPAF on a hand case; BDL target days enter at 0", {
    s <- .ah_sb()
    out <- as.data.frame(.ah_run(s$target, NULL, focal_dates = .ah_sb_focal))
    expect_equal(out$value, 100 * .ah_ms1(c(Cu = 10, B = 8000)), tolerance = 1e-8)
    ## every toxicant below detection: msPAF 0 (covered, not NA)
    z <- .ah_sb(tvals = c(Cu = 1, B = 100, DOC = 0.5),
                tdet = list(Cu = FALSE, B = FALSE))
    out0 <- as.data.frame(.ah_run(z$target, NULL, focal_dates = .ah_sb_focal))
    expect_false(is.na(out0$value))
    expect_equal(out0$value, 0, tolerance = 1e-12)
  })
})

describe("added_hazard: ara = 'none' with a reference supplied", {
  skip_on_cran()
  it("is honoured: no ARA, so value equals value_total and the no-reference run", {
    s <- .ah_sb(days = c(0, 30, 60, 90))
    a <- as.data.frame(.ah_run(s$target, s$reference, ara = "none",
                               also_total = TRUE, focal_dates = .ah_day(90)))
    b <- as.data.frame(.ah_run(s$target, NULL, focal_dates = .ah_day(90)))
    expect_equal(a$value, a$value_total, tolerance = 1e-10)
    expect_equal(a$value, b$value, tolerance = 1e-10)
  })
})

describe("added_hazard: daily_args$temperature", {
  skip_on_cran()
  it("is used for every point, overriding own temperature results (NH3-N)", {
    tg <- .ah_scn_T1(); ref <- .ah_scn_REF()
    fd <- .ah_day(c(60, 100, 330))
    tmp <- data.frame(datetime = seq(.ah_day(-30), .ah_day(520), by = "day"),
                      value = 25)
    with_t <- as.data.frame(.ah_run(tg, ref, focal_dates = fd,
                                    daily_args = list(temperature = tmp)))
    without <- as.data.frame(.ah_run(tg, ref, focal_dates = fd))
    o <- .ah_oracle(tg, ref, focal_dates = fd, temp = 25)$chronic
    .ah_expect_close(with_t$value, o$value, label = "value with temperature")
    expect_gt(max(abs(with_t$value - without$value)), 1e-6)
    o0 <- .ah_oracle(tg, ref, focal_dates = fd)$chronic
    .ah_expect_close(without$value, o0$value, label = "value without")
  })
})

describe("added_hazard: value_total is the target alone (msT0), BDL days at 0", {
  skip_on_cran()

  it("a target analyte below detection every day adds nothing to the total", {
    s <- .ah_sb(tvals = c(Cu = 1, B = 8000, DOC = 0.5), tdet = list(Cu = FALSE))
    args <- list(also_total = TRUE, contributions = TRUE, focal_dates = .ah_sb_focal)
    r <- do.call(.ah_run, c(list(s$target, s$reference), args))
    no_cu <- s$target[s$target$analyte != "Cu", ]
    r2 <- do.call(.ah_run, c(list(no_cu, s$reference), args))
    expect_equal(as.data.frame(r)$value_total, as.data.frame(r2)$value_total,
                 tolerance = 1e-10)
    expect_equal(as.data.frame(r)$value_total, 100 * .ah_paf1(8000, "B"),
                 tolerance = 1e-8)
    ctb <- as.data.frame(attr(r, "contributions"))
    tt <- ctb[ctb$variant == "total", ]
    expect_equal(sum(tt$contribution[tt$analyte == "Cu"]), 0)
    expect_equal(sum(tt$contribution), as.data.frame(r)$value_total, tolerance = 1e-6)
  })

  it("equals the ara = 'none' value_total on the analytes the reference has", {
    t1 <- .ah_scn_T1(); ref <- .ah_scn_REF()
    fd <- .ah_day(seq(60, 480, by = 30))
    a <- as.data.frame(.ah_run(t1, ref, also_total = TRUE, focal_dates = fd))
    t1r <- t1[t1$analyte != "Cr", ]
    b <- as.data.frame(.ah_run(t1r, NULL, ara = "none", also_total = TRUE,
                               focal_dates = fd))
    .ah_expect_close(a$value_total, b$value_total, tol = 1e-8, label = "value_total")
    expect_true(sum(!is.na(a$value_total)) > 5)
  })
})

## ═════════════════════════════════════════════════════════════════════════════
## Hand-computed checks
## ═════════════════════════════════════════════════════════════════════════════

describe("added_hazard: the D formula on a constant series", {
  skip_on_cran()
  s <- .ah_sb()
  out <- .ah_get("sb", function() as.data.frame(.ah_run(
    s$target, s$reference, also_total = TRUE, contributions = TRUE,
    focal_dates = .ah_sb_focal)))
  msT <- .ah_ms1(c(Cu = 10, B = 8000)); msR <- .ah_ms1(c(Cu = 2, B = 2000))
  D <- max((msT - msR) / (1 - msR), 0)

  it("gives 100 * max((msT - msR) / (1 - msR), 0) with full coverage", {
    expect_equal(nrow(out), 1L)
    expect_equal(out$value, 100 * D, tolerance = 1e-8)
    expect_equal(out$coverage, 1, tolerance = 1e-12)
    expect_gt(D, 0.05)    # a non-trivial case
  })

  it("gives value_total = 100 * msT", {
    expect_equal(out$value_total, 100 * msT, tolerance = 1e-8)
  })

  it("floors D at zero when the target is below the reference", {
    z <- .ah_sb(tvals = c(Cu = 1, B = 1000, DOC = 0.5))
    o <- as.data.frame(.ah_run(z$target, z$reference, focal_dates = .ah_sb_focal))
    expect_equal(o$value, 0, tolerance = 1e-12)
  })

  it("apportions D by added single-substance PAF (variant 'ara')", {
    ctb <- as.data.frame(attr(.ah_run(s$target, s$reference, also_total = TRUE,
                                      contributions = TRUE,
                                      focal_dates = .ah_sb_focal),
                              "contributions"))
    expect_true(all(c("variant", "site_id", "focal_date", "analyte", "group",
                      "contribution") %in% names(ctb)))
    expect_setequal(unique(ctb$variant), c("ara", "total"))
    add <- c(Cu = .ah_paf1(10, "Cu") - .ah_paf1(2, "Cu"),
             B = .ah_paf1(8000, "B") - .ah_paf1(2000, "B"))
    exp_a <- 100 * D * add / sum(add)
    a <- ctb[ctb$variant == "ara", ]
    for (an in c("Cu", "B")) {
      expect_equal(sum(a$contribution[a$analyte == an]), unname(exp_a[an]),
                   tolerance = 1e-8)
    }
    expect_equal(sum(a$contribution), out$value, tolerance = 1e-6)
    expect_identical(unique(a$group[a$analyte == "Cu"]), "ionoregulatory (gill)")
    expect_identical(unique(a$group[a$analyte == "B"]), "boron")
  })

  it("apportions msT by PAF_i(T) (variant 'total') and sums to value_total", {
    ctb <- as.data.frame(attr(.ah_run(s$target, s$reference, also_total = TRUE,
                                      contributions = TRUE,
                                      focal_dates = .ah_sb_focal),
                              "contributions"))
    tt <- ctb[ctb$variant == "total", ]
    pT <- c(Cu = .ah_paf1(10, "Cu"), B = .ah_paf1(8000, "B"))
    exp_t <- 100 * msT * pT / sum(pT)
    for (an in c("Cu", "B")) {
      expect_equal(sum(tt$contribution[tt$analyte == an]), unname(exp_t[an]),
                   tolerance = 1e-8)
    }
    expect_equal(sum(tt$contribution), out$value_total, tolerance = 1e-6)
  })
})

describe("added_hazard: detection-limit rules (hand-computed)", {
  skip_on_cran()

  one <- function(s, ...) as.data.frame(.ah_run(
    s$target, s$reference, focal_dates = .ah_sb_focal, ...))

  it("enters a below-detection reference result at LOR/2", {
    ## reference Cu is below detection (LOR 2) on every date -> R = 1
    s <- .ah_sb(tvals = c(Cu = 10, DOC = 0.5), rvals = c(Cu = 2, DOC = 0.5),
                rdet = list(Cu = FALSE))
    msT <- .ah_paf1(10, "Cu"); msR <- .ah_paf1(1, "Cu")
    expect_equal(one(s)$value, 100 * (msT - msR) / (1 - msR), tolerance = 1e-8)
  })

  it("never lets the static gap-filler fall below LOR/2 (geometric mean with 0 would be ~0)", {
    ## reference sampled on days 0 and 30 only, both below detection (LOR 2);
    ## the target is sampled long after the reference's 120-day hold has
    ## ended, so every target day uses the gap-filler s = max(1, 2/2) = 1.
    tg <- .ah_const("T1", c(300, 330, 360, 390, 420), list(Cu = 10, DOC = 0.5))
    rf <- .ah_const("REF", c(0, 30), list(Cu = 2, DOC = 0.5),
                    detected = list(Cu = FALSE))
    out <- as.data.frame(.ah_run(tg, rf, focal_dates = .ah_day(420)))
    msT <- .ah_paf1(10, "Cu"); msR <- .ah_paf1(1, "Cu")
    D <- (msT - msR) / (1 - msR)
    expect_false(is.na(out$value))
    expect_equal(out$value, 100 * D, tolerance = 1e-8)
    ## the buggy summary (geometric mean with BDL as 0) would give D ~ msT
    expect_gt(abs(out$value - 100 * msT), 1)
    d <- as.data.frame(attr(.ah_run(tg, rf, focal_dates = .ah_day(420)), "daily"))
    row <- d[d$date == .ah_day(400) & d$site_id == "T1", ]
    expect_equal(row$mspaf_reference, 100 * msR, tolerance = 1e-8)
  })

  it("uses the static summary max(geomean, median(LOR)/2) with mixed results", {
    ## reference Cu: 8 (det), <2, <2 -> geomean(8,1,1) = 2 > LOR/2 = 1
    tg <- .ah_const("T1", c(300, 330, 360, 390, 420), list(Cu = 10, DOC = 0.5))
    rf <- .ah_rows("REF", list(
      Cu = list(day = c(0, 30, 60), value = c(8, 2, 2),
                detected = c(TRUE, FALSE, FALSE)),
      DOC = list(day = 0, value = 0.5)))
    out <- as.data.frame(.ah_run(tg, rf, focal_dates = .ah_day(420)))
    s <- max(exp(mean(log(c(8, 1, 1)))), 1)
    msT <- .ah_paf1(10, "Cu"); msR <- .ah_paf1(s, "Cu")
    expect_equal(out$value, 100 * (msT - msR) / (1 - msR), tolerance = 1e-8)
  })

  it("gives a target below detection zero excess and zero ARA contribution", {
    ## Cu below detection on every target date; B elevated
    s <- .ah_sb(tvals = c(Cu = 1, B = 8000, DOC = 0.5),
                tdet = list(Cu = FALSE))
    out <- as.data.frame(.ah_run(s$target, s$reference, also_total = TRUE,
                                 contributions = TRUE, focal_dates = .ah_sb_focal))
    ctb <- as.data.frame(attr(.ah_run(s$target, s$reference, also_total = TRUE,
                                      contributions = TRUE,
                                      focal_dates = .ah_sb_focal),
                              "contributions"))
    a <- ctb[ctb$variant == "ara", ]
    expect_equal(sum(a$contribution[a$analyte == "Cu"]), 0)
    ## Cu behaves as T := R, so the hazard comes from B alone
    cuR <- c(Cu = 2)
    msT <- .ah_ms1(c(Cu = 2, B = 8000)); msR <- .ah_ms1(c(Cu = 2, B = 2000))
    expect_equal(out$value, 100 * (msT - msR) / (1 - msR), tolerance = 1e-8)
    expect_equal(sum(a$contribution[a$analyte == "B"]), out$value, tolerance = 1e-6)
  })

  it("drops an analyte the reference never measured (both sides)", {
    s <- .ah_sb()
    cr <- .ah_const("T1", .ah_sb_days, list(Cr = 500))
    base <- one(s)
    withcr <- as.data.frame(.ah_run(dplyr::bind_rows(s$target, cr), s$reference,
                                    contributions = TRUE,
                                    focal_dates = .ah_sb_focal))
    expect_equal(withcr$value, base$value, tolerance = 1e-10)
    cc <- as.data.frame(attr(.ah_run(dplyr::bind_rows(s$target, cr), s$reference,
                                     contributions = TRUE,
                                     focal_dates = .ah_sb_focal), "contributions"))
    expect_equal(sum(cc$contribution[cc$analyte == "Cr"]), 0)
  })
})

describe("added_hazard: time-weighting and coverage gating (hand-computed)", {
  skip_on_cran()

  ## Target sampled every day 0..199, Cu 10 for days 0..99 and 40 from day 100;
  ## reference Cu = 2 throughout; DOC 0.5 on both sides.
  step <- function() .ah_get("step", function() {
    days <- 0:199
    cu <- ifelse(days < 100, 10, 40)
    tg <- dplyr::bind_rows(
      .ah_rows("T1", list(Cu = list(day = days, value = cu),
                          DOC = list(day = c(0, 199), value = 0.5))))
    rf <- .ah_rows("REF", list(Cu = list(day = c(0, 100, 199), value = 2),
                               DOC = list(day = 0, value = 0.5)))
    list(target = tg, reference = rf)
  })
  pT1 <- .ah_paf1(10, "Cu"); pT2 <- .ah_paf1(40, "Cu"); pR <- .ah_paf1(2, "Cu")
  D1 <- (pT1 - pR) / (1 - pR); D2 <- (pT2 - pR) / (1 - pR)
  w <- exp(-(0:365) / 90)

  run <- function(focal, ...) {
    s <- step()
    as.data.frame(.ah_run(s$target, s$reference, also_total = TRUE,
                          focal_dates = .ah_day(focal), ...))
  }

  it("weights by exp(-age / 90), renormalised over covered days", {
    out <- run(199)
    ## ages 0..99 are the D2 days (days 199..100), ages 100..199 the D1 days
    num <- sum(w[1:100]) * D2 + sum(w[101:200]) * D1
    expect_equal(out$value, 100 * num / sum(w[1:200]), tolerance = 1e-8)
    numT <- sum(w[1:100]) * pT2 + sum(w[101:200]) * pT1
    expect_equal(out$value_total, 100 * numT / sum(w[1:200]), tolerance = 1e-8)
    expect_equal(out$coverage, sum(w[1:200]) / sum(w), tolerance = 1e-10)
  })

  it("reports window counts (distinct sample days; covered days that are not samples)", {
    out <- run(c(30, 199))
    expect_equal(as.numeric(out$n_samples_in_window), c(31, 200))
    expect_equal(as.numeric(out$n_imputed_in_window), c(0, 0))
  })

  it("exposes the daily D, total and reference series (in %)", {
    s <- step()
    d <- as.data.frame(attr(.ah_run(s$target, s$reference, also_total = TRUE,
                                    focal_dates = .ah_day(199)), "daily"))
    r <- d[d$date == .ah_day(150) & d$site_id == "T1", ]
    expect_equal(r$mspaf, 100 * D2, tolerance = 1e-8)
    expect_equal(r$mspaf_total, 100 * pT2, tolerance = 1e-8)
    expect_equal(r$mspaf_reference, 100 * pR, tolerance = 1e-8)
    expect_equal(as.numeric(r$n_analytes), 1)
  })

  it("gates at the default min_coverage of 0.2", {
    out <- run(c(5, 30, 199))
    cv <- vapply(c(5, 30, 199), function(k) sum(w[1:(k + 1)]) / sum(w), 0)
    expect_equal(out$coverage, cv, tolerance = 1e-10)
    expect_true(is.na(out$value[1]))
    expect_false(is.na(out$value[2]))
    expect_false(is.na(out$value[3]))
    expect_true(is.na(out$value_total[1]))
    ## coverage is reported even where value is gated
    expect_false(anyNA(out$coverage))
  })

  it("honours a custom min_coverage (value NA exactly below it)", {
    out <- run(c(5, 30, 199), min_coverage = 0.5)
    expect_identical(is.na(out$value), c(TRUE, TRUE, FALSE))
    out <- run(c(5, 30, 199), min_coverage = 0.95)
    expect_identical(is.na(out$value), c(TRUE, TRUE, TRUE))
    out <- run(c(5, 30, 199), min_coverage = 0)
    expect_false(anyNA(out$value))
    expect_equal(out$value[1], 100 * D1, tolerance = 1e-8)  # only D1 days covered
  })
})

## ═════════════════════════════════════════════════════════════════════════════
## Contributions against the oracle (main scenario)
## ═════════════════════════════════════════════════════════════════════════════

describe("added_hazard: contributions", {
  skip_on_cran()

  ctb <- function() as.data.frame(attr(.ah_main(), "contributions"))

  it("has the documented shape", {
    x <- ctb()
    expect_true(all(c("variant", "site_id", "focal_date", "analyte", "group",
                      "contribution") %in% names(x)))
    expect_setequal(unique(x$variant), c("ara", "total"))
    expect_true(all(x$group == .ah_moa(x$analyte)))
    expect_false(any(x$analyte %in% c("pH", "DOC", "hardness", "temperature")))
  })

  it("matches the oracle per analyte (1e-6)", {
    x <- ctb()
    o <- .ah_main_oracle()$contributions
    by <- c("variant", "site_id", "focal_date", "analyte")
    m <- .ah_join(x[c(by, "contribution")], o[c(by, "contribution")], by)
    has_exp <- !is.na(m$contribution.exp)
    got <- m$contribution.got
    got[is.na(got) & has_exp] <- 0     # an absent row means zero contribution
    expect_true(all(abs(got[has_exp] - m$contribution.exp[has_exp]) < 1e-6))
    ## where value is NA the contributions are NA (or absent)
    expect_true(all(is.na(m$contribution.got[!has_exp])))
  })

  it("sums to value (variant 'ara') and value_total (variant 'total')", {
    x <- ctb(); out <- as.data.frame(.ah_main())
    for (v in c("ara", "total")) {
      xv <- x[x$variant == v, ]
      s <- stats::aggregate(contribution ~ site_id + focal_date, xv, sum)
      m <- merge(s, out, by = c("site_id", "focal_date"))
      ref <- if (v == "ara") m$value else m$value_total
      expect_gt(nrow(m), 100)
      expect_true(all(abs(m$contribution - ref) < 1e-6), label = v)
    }
  })

  it("has no contribution rows with a value where value is NA", {
    x <- ctb(); out <- as.data.frame(.ah_main())
    m <- merge(x, out[c("site_id", "focal_date", "value")],
               by = c("site_id", "focal_date"))
    expect_true(all(is.na(m$contribution[is.na(m$value)])))
  })

  it("group contribution = sum over the group's analytes", {
    x <- ctb(); x <- x[x$variant == "ara" & !is.na(x$contribution), ]
    g <- stats::aggregate(contribution ~ site_id + focal_date + group, x, sum)
    out <- as.data.frame(.ah_main())
    s <- stats::aggregate(contribution ~ site_id + focal_date, g, sum)
    m <- merge(s, out, by = c("site_id", "focal_date"))
    expect_true(all(abs(m$contribution - m$value) < 1e-6))
  })
})

## ═════════════════════════════════════════════════════════════════════════════
## Regression facts
## ═════════════════════════════════════════════════════════════════════════════

describe("regression: the 120-day cap stops a 16-year-old sample being carried", {
  skip_on_cran()

  ## one elevated metals sample in 2005, nothing for 16 years, then monthly
  b05 <- as.Date("2005-01-01")
  gap <- as.numeric(as.Date("2021-01-01") - b05)
  later <- gap + 30 * (0:12)
  mk <- function() {
    tg <- .ah_rows("T1", list(
      Cu = list(day = c(0, later), value = c(60, rep(2.5, 13))),
      Zn = list(day = c(0, later), value = c(600, rep(30, 13))),
      DOC = list(day = c(0, later), value = 1.5),
      pH = list(day = c(0, later), value = 7),
      hardness = list(day = c(0, later), value = 50)), base = b05)
    rf <- .ah_rows("REF", list(
      Cu = list(day = c(0, later), value = 1),
      Zn = list(day = c(0, later), value = 5),
      DOC = list(day = c(0, later), value = 1.5),
      pH = list(day = c(0, later), value = 7),
      hardness = list(day = c(0, later), value = 50)), base = b05)
    list(target = tg, reference = rf)
  }
  fd <- as.Date(c("2005-02-01", "2010-06-01", "2021-07-01"))

  it("is NA with zero coverage years after the old sample under the default cap", {
    s <- mk()
    out <- as.data.frame(.ah_run(s$target, s$reference, focal_dates = fd))
    ## The old sample sits in a 16-year gap, so only its own day is covered:
    ## a month later the covered kernel weight is far below the gate.
    expect_true(is.na(out$value[out$focal_date == fd[1]]))
    expect_true(is.na(out$value[out$focal_date == fd[2]]))
    expect_equal(out$coverage[out$focal_date == fd[2]], 0)
    expect_false(is.na(out$value[out$focal_date == fd[3]]))
  })

  it("carries it with span_cap = Inf", {
    s <- mk()
    out <- as.data.frame(.ah_run(s$target, s$reference, focal_dates = fd,
                                 daily_args = list(span_cap = Inf)))
    expect_false(is.na(out$value[out$focal_date == fd[2]]))
    expect_gt(out$coverage[out$focal_date == fd[2]], 0.2)
  })
})

describe("regression: a one-day x30 pulse does not move the chronic value", {
  skip_on_cran()

  ## Quiet period: target Cu 3, reference Cu 0.5 (DOC 0.5).  Monthly samples
  ## on both sides; pulse (x30) on the monthly sample day P.
  days <- seq(0, 420, by = 30); P <- 240
  run <- function(pulse, bracket) {
    cu <- rep(3, length(days))
    if (pulse) cu[days == P] <- 90
    tdays <- days; tcu <- cu
    if (bracket) { tdays <- c(days, P - 1, P + 1); tcu <- c(cu, 3, 3) }
    tg <- .ah_rows("T1", list(Cu = list(day = tdays, value = tcu),
                              DOC = list(day = c(0, 420), value = 0.5)))
    rf <- .ah_rows("REF", list(Cu = list(day = days, value = 0.5),
                               DOC = list(day = c(0, 420), value = 0.5)))
    o <- .ah_run(tg, rf, focal_dates = .ah_day(P + c(30, 90)))
    list(out = as.data.frame(o), daily = as.data.frame(attr(o, "daily")))
  }

  it("moves D by < 1 point at +30 d and +90 d when the pulse is bracketed daily", {
    base <- run(FALSE, TRUE); pul <- run(TRUE, TRUE)
    ## the pulse day's own D is large (the test is meaningful)
    expect_gt(pul$daily$mspaf[pul$daily$date == .ah_day(P)], 20)
    expect_lt(base$daily$mspaf[base$daily$date == .ah_day(P)], 60)
    d <- pul$out$value - base$out$value
    expect_true(all(abs(d) < 1))
    expect_true(all(d > 0))                   # it does register, just briefly
  })

  it("does move it by more than a point when the pulse is interpolated across the month", {
    base <- run(FALSE, FALSE); pul <- run(TRUE, FALSE)
    expect_gt(pul$out$value[1] - base$out$value[1], 1)
  })
})

## Regression (c), receptor below detection gives zero contribution, is covered
## by 'gives a target below detection zero excess and zero ARA contribution'.
