## Per-analyte and per-group contributions to msPAF
## ============================================================================
##
## CONTRACT.
##
## msPAF is assembled per sample (or per daily synthetic sample) as
##   C_adj_a = pmax(C_norm_a - ref_norm_a, 0)   (ARA; C_norm_a without ARA)
##   TU_a    = C_adj_a / HC50_a
##   msPAF_g = pnorm(log10(sum_{a in g} TU_a) / mean_{a in g}(sigma_a))
##             (concentration addition within a mode-of-action group g)
##   msPAF   = 100 * (1 - prod_g (1 - msPAF_g))   (independent action)
##
## Apportioning rule, per sample / day:
##   p_g = msPAF_g (proportion), clamped to <= 1 - 1e-12
##   L_g = -log(1 - p_g),  L = sum_g L_g
##   c_g = msPAF * L_g / L          (0 when L = 0)
##   share_a = TU_a / sum_{a in g} TU_a  (equal split if sum TU = 0, c_g > 0)
##   c_a = c_g * share_a
## so sum_a c_a = sum_g c_g = msPAF exactly (tolerance 1e-8).
##
## Chronic: each analyte's daily contribution series is aggregated with exactly
## the same weights, window and rows as the msPAF itself, so for every
## (site_id, focal_date) sum_a chronic contribution = chronic `value` (1e-8).
##
## API.
##  1. add_mspaf(..., contributions = FALSE).  TRUE attaches
##     attr(out, "contributions"): long frame with sample_id, site_id (and
##     datetime when the input carries it), analyte, group (mode-of-action
##     group), c_adj, tu, paf (single-analyte PAF, %), group_mspaf (%),
##     contribution (%).  Rows only for analytes that entered the msPAF.
##  2. mspaf_daily(..., contributions = FALSE).  TRUE attaches the daily long
##     breakdown (same columns, keyed by site_id + date) for the msPAF returned.
##  3. mspaf_pipeline(..., contributions = FALSE).  TRUE attaches the chronic
##     long breakdown: site_id, focal_date, analyte, group, contribution (%);
##     sum over analytes = value per (site_id, focal_date).
##  4. Default FALSE: no attribute; output otherwise identical, and TRUE does
##     not change the main output either.
##  5. Point mode only: TRUE with draws -> classed error
##     `hydroSense_contributions_draws`.  The flag is validated like other
##     logical arguments (error mentions `contributions`, and is never an
##     "unused argument" error).
##  6. With also_total = TRUE (skipped until that argument exists), the chronic
##     frame gains a `variant` column ("ara" / "total") summing to value and
##     value_total respectively.
##
## All data are synthetic (frames built here and helper-perf-fixtures.R).

library(testthat)

.ct_tol <- 1e-8

# ── add_mspaf() fixtures ──────────────────────────────────────────────────────

## Three analytes with no normalisation formula, in two mode-of-action groups:
## As(III) and As(V) ("arsenic"), B ("boron").  Concentrations in ug/L.
##   S1: all three well above HC50-ish levels
##   S2: modest
##   S3: tiny (msPAF near zero but positive)
##   S4: saturating boron (group msPAF -> 1; exercises the clamp)
.ct_point <- function() {
  vals <- list(
    S1 = c(`As(III)` = 2000, `As(V)` = 1500, B = 8000),
    S2 = c(`As(III)` = 300,  `As(V)` = 900,  B = 1200),
    S3 = c(`As(III)` = 1,    `As(V)` = 2,    B = 5),
    S4 = c(`As(III)` = 500,  `As(V)` = 400,  B = 1e9)
  )
  dts <- as.Date("2024-03-01") + 0:3
  dplyr::bind_rows(lapply(seq_along(vals), function(i) {
    tibble::tibble(
      sample_id = names(vals)[i], site_id = "H1", datetime = dts[i],
      analyte = names(vals[[i]]), value = unname(vals[[i]]),
      detected = TRUE, units.analyte = "ug/L"
    )
  }))
}

## Reference for ARA: constant concentrations, so ref_norm = that value under
## any summary.  S1/S2 sit above it, S5 (added below) sits at or below it.
.ct_reference <- function() {
  d <- as.Date("2023-01-01") + 30 * (0:5)
  dplyr::bind_rows(lapply(seq_along(d), function(i) {
    tibble::tibble(
      sample_id = paste0("R", i), site_id = "REF", datetime = d[i],
      analyte = c("As(III)", "As(V)", "B"), value = c(100, 200, 600),
      detected = TRUE, units.analyte = "ug/L"
    )
  }))
}

## S5: every analyte at or below the reference -> all C_adj = 0 -> msPAF 0.
.ct_point_with_zero <- function() {
  dplyr::bind_rows(
    .ct_point(),
    tibble::tibble(
      sample_id = "S5", site_id = "H1", datetime = as.Date("2024-03-05"),
      analyte = c("As(III)", "As(V)", "B"), value = c(50, 200, 100),
      detected = TRUE, units.analyte = "ug/L"
    )
  )
}

.ct_add <- function(df, reference = NULL, ...) {
  suppressWarnings(suppressMessages(
    add_mspaf(df, reference = reference, require_temperature = FALSE, ...)
  ))
}

.ct_mspaf_rows <- function(x) {
  as.data.frame(x[x$analyte == "msPAF", , drop = FALSE])
}

.ct_contrib_cols <- c("sample_id", "site_id", "analyte", "group", "c_adj",
                      "tu", "paf", "group_mspaf", "contribution")

## SSD parameters used to build the hand-computed expectations.
.ct_ssd <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      meta  <- hydroSense:::.load_analyte_metadata(NULL)
      cache <<- hydroSense:::derive_ssd_params(meta, "multi", NULL)
    }
    cache
  }
})

## Expected per-analyte rows for one sample from the formula in the header.
.ct_expected <- function(c_adj_named) {
  p <- .ct_ssd()
  p <- p[match(names(c_adj_named), p$analyte), ]
  d <- data.frame(analyte = p$analyte, group = p$moa_group,
                  c_adj = unname(c_adj_named), hc50 = p$hc50,
                  sigma = p$sigma, stringsAsFactors = FALSE)
  d$tu <- d$c_adj / d$hc50
  grp <- split(d, d$group)
  pg  <- vapply(grp, function(g) {
    in_ca <- g$c_adj > 0
    if (!any(in_ca)) return(0)
    stats::pnorm(log10(sum(g$tu[in_ca])) / mean(g$sigma[in_ca]))
  }, numeric(1))
  mspaf <- 100 * (1 - prod(1 - pg))
  Lg <- -log(1 - pmin(pg, 1 - 1e-12))
  cg <- if (sum(Lg) > 0) mspaf * Lg / sum(Lg) else 0 * Lg
  d$group_mspaf <- 100 * pg[d$group]
  d$c_group     <- cg[d$group]
  d$share <- ave(d$tu, d$group, FUN = function(t) t / sum(t))
  d$contribution <- d$c_group * d$share
  list(rows = d, mspaf = mspaf, c_g = cg)
}


describe("add_mspaf(contributions = TRUE): hand-computed two-group case", {

  it("matches c_g and c_a from the apportioning formula (no ARA)", {
    df  <- .ct_point()
    out <- .ct_add(df, contributions = TRUE)
    ctb <- attr(out, "contributions")
    expect_s3_class(ctb, "data.frame")
    expect_true(all(.ct_contrib_cols %in% names(ctb)))

    ms <- .ct_mspaf_rows(out)
    for (sid in c("S1", "S2", "S3")) {
      x   <- df[df$sample_id == sid, ]
      exp <- .ct_expected(stats::setNames(x$value, x$analyte))
      got <- as.data.frame(ctb[ctb$sample_id == sid, ])
      got <- got[match(exp$rows$analyte, got$analyte), ]

      expect_equal(ms$value[ms$sample_id == sid], exp$mspaf, tolerance = 1e-10)
      expect_identical(got$group, exp$rows$group)
      expect_equal(got$c_adj, exp$rows$c_adj, tolerance = 1e-10)
      expect_equal(got$tu, exp$rows$tu, tolerance = 1e-10)
      expect_equal(got$group_mspaf, exp$rows$group_mspaf, tolerance = 1e-10)
      expect_equal(got$contribution, exp$rows$contribution, tolerance = 1e-10)

      ## Group contributions: sum of the group's analyte rows = c_g.
      cg_got <- tapply(got$contribution, got$group, sum)
      expect_equal(as.vector(cg_got[names(exp$c_g)]), unname(exp$c_g),
                   tolerance = 1e-10)
    }
  })

  it("within a group, analyte shares are proportional to TU", {
    out <- .ct_add(.ct_point(), contributions = TRUE)
    ctb <- as.data.frame(attr(out, "contributions"))
    ars <- ctb[ctb$group == "arsenic" & ctb$sample_id %in% c("S1", "S2"), ]
    for (sid in c("S1", "S2")) {
      g <- ars[ars$sample_id == sid, ]
      expect_equal(nrow(g), 2L)
      expect_equal(g$contribution / sum(g$contribution), g$tu / sum(g$tu),
                   tolerance = 1e-10)
    }
  })

  it("carries the single-analyte PAF (%) and the analyte_pafs diagnostics", {
    out <- .ct_add(.ct_point(), contributions = TRUE)
    ctb <- as.data.frame(attr(out, "contributions"))
    ap  <- as.data.frame(attr(out, "analyte_pafs"))
    m   <- merge(ctb, ap, by = c("sample_id", "site_id", "analyte"))
    expect_equal(nrow(m), nrow(ctb))
    expect_equal(m$paf, 100 * m$PAF, tolerance = 1e-10)
    expect_equal(m$c_adj, m$C_adj, tolerance = 1e-12)
    expect_identical(m$group, m$moa_group)
    expect_true(all(m$paf >= 0 & m$paf <= 100))
  })

  it("stays finite and exact when a group saturates (clamp)", {
    out <- .ct_add(.ct_point(), contributions = TRUE)
    ctb <- as.data.frame(attr(out, "contributions"))
    s4  <- ctb[ctb$sample_id == "S4", ]
    ms  <- .ct_mspaf_rows(out)
    expect_true(all(is.finite(s4$contribution)))
    expect_true(all(s4$contribution >= 0))
    expect_equal(sum(s4$contribution), ms$value[ms$sample_id == "S4"],
                 tolerance = .ct_tol)
    ## Boron dominates the saturated sample.
    expect_gt(s4$contribution[s4$analyte == "B"], 0.9 * sum(s4$contribution))
  })
})


describe("add_mspaf(contributions = TRUE): exact-sum identity and keys", {

  for (ara in c(FALSE, TRUE)) {
    it(paste0("sum of contributions = msPAF per sample (ARA = ", ara, ")"), {
      df  <- .ct_point_with_zero()
      ref <- if (ara) .ct_reference() else NULL
      out <- .ct_add(df, reference = ref, contributions = TRUE)
      ms  <- .ct_mspaf_rows(out)
      ctb <- as.data.frame(attr(out, "contributions"))
      expect_true(all(ctb$contribution >= 0))
      tot <- tapply(ctb$contribution, ctb$sample_id, sum)
      got <- as.vector(tot[ms$sample_id])
      got[is.na(got)] <- 0          # a sample without rows contributes 0
      expect_equal(got, ms$value, tolerance = .ct_tol)

      ## Group level: sum over groups of (sum over the group's analytes) = msPAF.
      by_g <- stats::aggregate(contribution ~ sample_id + group, ctb, sum)
      tot_g <- tapply(by_g$contribution, by_g$sample_id, sum)
      expect_equal(as.vector(tot_g[names(tot)]), as.vector(tot),
                   tolerance = .ct_tol)
    })
  }

  it("ARA: C_adj is the reference-adjusted concentration", {
    out <- .ct_add(.ct_point(), reference = .ct_reference(),
                   contributions = TRUE)
    ctb <- as.data.frame(attr(out, "contributions"))
    exp <- .ct_expected(c(`As(III)` = 2000 - 100, `As(V)` = 1500 - 200,
                          B = 8000 - 600))
    s1  <- ctb[ctb$sample_id == "S1", ]
    s1  <- s1[match(exp$rows$analyte, s1$analyte), ]
    expect_equal(s1$c_adj, exp$rows$c_adj, tolerance = 1e-8)
    expect_equal(s1$contribution, exp$rows$contribution, tolerance = 1e-10)
  })

  it("a zero-msPAF sample has zero contributions", {
    out <- .ct_add(.ct_point_with_zero(), reference = .ct_reference(),
                   contributions = TRUE)
    ms  <- .ct_mspaf_rows(out)
    expect_equal(ms$value[ms$sample_id == "S5"], 0)
    ctb <- as.data.frame(attr(out, "contributions"))
    s5  <- ctb$contribution[ctb$sample_id == "S5"]
    expect_true(all(s5 == 0))
    expect_equal(sum(s5), 0)
  })

  it("keys match the msPAF rows (sample_id, site_id, datetime)", {
    out <- .ct_add(.ct_point_with_zero(), reference = .ct_reference(),
                   contributions = TRUE)
    ms  <- .ct_mspaf_rows(out)
    ctb <- as.data.frame(attr(out, "contributions"))
    expect_true("datetime" %in% names(ctb))
    k_ctb <- unique(ctb[, c("sample_id", "site_id", "datetime")])
    k_ms  <- ms[, c("sample_id", "site_id", "datetime")]
    ## every contribution key is an msPAF row, with the same datetime
    m <- merge(k_ctb, k_ms, by = c("sample_id", "site_id", "datetime"))
    expect_equal(nrow(m), nrow(k_ctb))
    ## every msPAF row with a positive value has contributions
    pos <- ms$sample_id[ms$value > 0]
    expect_true(all(pos %in% ctb$sample_id))
    ## one row per (sample, analyte)
    expect_false(anyDuplicated(ctb[, c("sample_id", "analyte")]) > 0)
  })

  it("rows only for analytes that entered that sample's msPAF", {
    df <- .ct_point()
    ## Only S1 reports Cr (its own group "chromium"); a non-SSD row
    ## (Standing water level) never appears.
    extra <- tibble::tibble(
      sample_id = c("S1", "S1"), site_id = "H1",
      datetime = as.Date("2024-03-01"),
      analyte = c("Cr", "Standing water level"), value = c(40, 3),
      detected = TRUE, units.analyte = c("ug/L", "m")
    )
    out <- .ct_add(dplyr::bind_rows(df, extra), contributions = TRUE)
    ctb <- as.data.frame(attr(out, "contributions"))
    ap  <- as.data.frame(attr(out, "analyte_pafs"))
    expect_false("Standing water level" %in% ctb$analyte)
    expect_true("Cr" %in% ctb$analyte[ctb$sample_id == "S1"])
    expect_false("Cr" %in% ctb$analyte[ctb$sample_id != "S1"])
    for (sid in unique(ctb$sample_id)) {
      used <- ap$analyte[ap$sample_id == sid]
      pos  <- ap$analyte[ap$sample_id == sid & ap$C_adj > 0]
      have <- ctb$analyte[ctb$sample_id == sid]
      expect_true(all(have %in% used))
      expect_true(all(pos %in% have))
    }
  })
})


describe("add_mspaf(): default, draws and validation", {

  it("default FALSE: no attribute; TRUE leaves the main output unchanged", {
    df   <- .ct_point()
    base <- .ct_add(df)
    expect_null(attr(base, "contributions"))
    expect_null(attr(.ct_add(df, contributions = FALSE), "contributions"))
    with <- .ct_add(df, contributions = TRUE)
    attr(with, "contributions") <- NULL
    expect_identical(with, base)
  })

  it("draw-carrier input with contributions = TRUE is a classed error", {
    df <- .ct_point()
    dr <- dplyr::bind_rows(lapply(1:2, function(i) {
      x <- df
      x$draw_id <- i
      x$value <- x$value * (1 + 0.1 * i)
      x
    }))
    expect_error(
      suppressMessages(add_mspaf(dr, require_temperature = FALSE,
                                 return = "draws", contributions = TRUE)),
      class = "hydroSense_contributions_draws"
    )
  })

  it("validates the flag like other logical arguments", {
    for (bad in list("yes", NA, c(TRUE, FALSE), 1)) {
      err <- expect_error(
        suppressMessages(add_mspaf(.ct_point(), require_temperature = FALSE,
                                   contributions = bad))
      )
      expect_match(conditionMessage(err), "contributions")
      expect_no_match(conditionMessage(err), "unused argument")
    }
  })
})


# ── mspaf_daily() ─────────────────────────────────────────────────────────────

## Two points (P1, P2) of the small structural target over two months; groups:
## ionoregulatory (gill) {Cu, Zn, Ni} and ammonia {NH3-N}.
.ct_daily_target <- function() {
  tgt <- .pf_small_target()
  tgt[tgt$site_id %in% c("P1", "P2"), , drop = FALSE]
}

.ct_temperature <- function() {
  tmp <- .pf_scenario()$temperature
  tmp[tmp$datetime >= as.Date("2021-12-15") &
        tmp$datetime <= as.Date("2022-03-10"), ]
}

.ct_daily <- function(reference, ...) {
  suppressWarnings(suppressMessages(mspaf_daily(
    .ct_daily_target(), temperature = .ct_temperature(),
    reference = reference, reference_model = .pf_reference_model(),
    interpolation = "model",
    start = as.Date("2022-01-01"), end = as.Date("2022-02-28"), ...
  )))
}

## Memoised daily runs (each ~1-3 s).
.ct_daily_cache <- local({
  cache <- list()
  function(key) {
    if (is.null(cache[[key]])) {
      ref <- if (grepl("^ara", key)) .pf_scenario()$reference else NULL
      cache[[key]] <<- if (grepl("_ctb$", key)) {
        .ct_daily(ref, contributions = TRUE)
      } else {
        .ct_daily(ref)
      }
    }
    cache[[key]]
  }
})

describe("mspaf_daily(contributions = TRUE)", {

  for (ara in c("none", "ara")) {
    it(paste0("sum over analytes = daily mspaf for every (site, date) [", ara, "]"), {
      out <- .ct_daily_cache(paste0(ara, "_ctb"))
      ctb <- as.data.frame(attr(out, "contributions"))
      expect_s3_class(attr(out, "contributions"), "data.frame")
      expect_true(all(c("site_id", "date", "analyte", "group", "c_adj", "tu",
                        "paf", "group_mspaf", "contribution") %in% names(ctb)))
      expect_s3_class(ctb$date, "Date")
      expect_gt(nrow(ctb), 0L)
      expect_true(all(is.finite(ctb$contribution) & ctb$contribution >= 0))

      tot <- stats::aggregate(contribution ~ site_id + date, ctb, sum)
      m <- merge(as.data.frame(out[, c("site_id", "date", "mspaf")]), tot,
                 by = c("site_id", "date"), all.x = TRUE)
      m$contribution[is.na(m$contribution)] <- 0
      expect_equal(nrow(m), nrow(out))
      expect_equal(m$contribution, m$mspaf, tolerance = .ct_tol)

      ## Both groups present and each one's analytes split by TU.
      expect_setequal(unique(ctb$group),
                      c("ionoregulatory (gill)", "ammonia"))
      g <- ctb[ctb$group == "ionoregulatory (gill)", ]
      g$share_c  <- ave(g$contribution, g$site_id, g$date,
                        FUN = function(v) v / sum(v))
      g$share_tu <- ave(g$tu, g$site_id, g$date, FUN = function(v) v / sum(v))
      expect_equal(g$share_c, g$share_tu, tolerance = 1e-8)
    })
  }

  it("keys match the returned rows and the analyte_pafs breakdown", {
    out <- .ct_daily_cache("ara_ctb")
    ctb <- as.data.frame(attr(out, "contributions"))
    ap  <- as.data.frame(attr(out, "analyte_pafs"))
    k_ctb <- unique(ctb[, c("site_id", "date")])
    k_out <- as.data.frame(out[, c("site_id", "date")])
    expect_equal(nrow(merge(k_ctb, k_out)), nrow(k_ctb))
    expect_false(anyDuplicated(ctb[, c("site_id", "date", "analyte")]) > 0)
    m <- merge(ctb, ap, by = c("site_id", "date", "analyte"))
    expect_equal(nrow(m), nrow(ctb))
    expect_equal(m$c_adj, m$C_adj, tolerance = 1e-12)
    expect_equal(m$paf, 100 * m$PAF, tolerance = 1e-10)
    expect_identical(m$group, m$moa_group)
  })

  it("default FALSE: no attribute; TRUE leaves the main output unchanged", {
    base <- .ct_daily_cache("ara")
    expect_null(attr(base, "contributions"))
    with <- .ct_daily_cache("ara_ctb")
    attr(with, "contributions") <- NULL
    expect_identical(with, base)
  })

  it("ndraws with contributions = TRUE is a classed error", {
    expect_error(
      suppressWarnings(suppressMessages(mspaf_daily(
        .ct_daily_target(), temperature = .ct_temperature(),
        reference_model = .pf_reference_model(), interpolation = "model",
        start = as.Date("2022-01-01"), end = as.Date("2022-01-10"),
        ndraws = 2L, seed = 1L, contributions = TRUE
      ))),
      class = "hydroSense_contributions_draws"
    )
  })

  it("validates the flag like other logical arguments", {
    for (bad in list("yes", NA, c(TRUE, FALSE))) {
      err <- expect_error(suppressMessages(mspaf_daily(
        .ct_daily_target(), require_temperature = FALSE, contributions = bad
      )))
      expect_match(conditionMessage(err), "contributions")
      expect_no_match(conditionMessage(err), "unused argument")
    }
  })
})


# ── mspaf_pipeline() (chronic) ────────────────────────────────────────────────

.ct_pipeline <- function(ara, ...) {
  sc <- .pf_scenario()
  daily_args <- list(temperature = .ct_temperature(),
                     start = as.Date("2022-01-01"), end = as.Date("2022-02-28"))
  if (ara) daily_args$reference <- sc$reference
  suppressWarnings(suppressMessages(mspaf_pipeline(
    .ct_daily_target(), reference = sc$reference, hydro = sc$rain,
    impute = FALSE,
    reference_args = list(
      min_obs_model = 20L,
      api_tau_bounds_short = c(3, 10),
      api_tau_bounds_long  = c(30, 90)
    ),
    daily_args = daily_args,
    chronic = TRUE,
    focal_dates = as.Date(c("2022-01-20", "2022-02-10", "2022-02-28")),
    tau = 90, window = 365, ...
  )))
}

.ct_strip_pipeline <- function(x) {
  attributes(x)[c("daily", "reference_model", "imputation_model")] <- NULL
  x
}

## Chronic sum check: for every (site, focal date) with a value, the
## contributions sum to it.
.ct_expect_chronic_sum <- function(out, ctb, value_col = "value") {
  tot <- stats::aggregate(contribution ~ site_id + focal_date, ctb, sum)
  main <- as.data.frame(out[, c("site_id", "focal_date", value_col)])
  main <- main[!is.na(main[[value_col]]), ]
  m <- merge(main, tot, by = c("site_id", "focal_date"), all.x = TRUE)
  expect_equal(nrow(m), nrow(main))
  m$contribution[is.na(m$contribution)] <- 0
  expect_equal(m$contribution, m[[value_col]], tolerance = .ct_tol)
}

describe("mspaf_pipeline(contributions = TRUE): chronic breakdown", {

  for (ara in c(FALSE, TRUE)) {
    it(paste0("sum over analytes = chronic value per (site, focal_date) [ARA = ",
              ara, "]"), {
      out <- .ct_pipeline(ara, contributions = TRUE)
      ctb <- attr(out, "contributions")
      expect_s3_class(ctb, "data.frame")
      ctb <- as.data.frame(ctb)
      expect_true(all(c("site_id", "focal_date", "analyte", "group",
                        "contribution") %in% names(ctb)))
      expect_s3_class(ctb$focal_date, "Date")
      expect_false(anyDuplicated(ctb[, c("site_id", "focal_date",
                                         "analyte")]) > 0)
      fin <- ctb$contribution[!is.na(ctb$contribution)]
      expect_gt(length(fin), 0L)
      expect_true(all(fin >= 0))
      .ct_expect_chronic_sum(out, ctb)

      ## Keys are a subset of the chronic frame's (site, focal_date).
      k_ctb <- unique(ctb[, c("site_id", "focal_date")])
      k_out <- as.data.frame(out[, c("site_id", "focal_date")])
      expect_equal(nrow(merge(k_ctb, k_out)), nrow(k_ctb))

      ## Group-level totals are derivable by summing; no duplicate group rows.
      expect_false(any(is.na(ctb$analyte)))
      expect_true(all(ctb$analyte %in% c("Cu", "Zn", "Ni", "NH3-N")))

      ## TRUE leaves the chronic frame and its daily attribute unchanged.
      base <- .ct_pipeline(ara)
      expect_null(attr(base, "contributions"))
      expect_identical(.ct_strip_pipeline(
        structure(out, contributions = NULL)), .ct_strip_pipeline(base))
      expect_identical(attr(out, "daily")[names(attr(base, "daily"))],
                       attr(base, "daily"))
    })
  }

  it("ndraws with contributions = TRUE is a classed error", {
    sc <- .pf_scenario()
    expect_error(
      suppressWarnings(suppressMessages(mspaf_pipeline(
        .ct_daily_target(), reference = sc$reference, hydro = sc$rain,
        impute = FALSE,
        reference_args = list(min_obs_model = 20L,
                              api_tau_bounds_short = c(3, 10),
                              api_tau_bounds_long  = c(30, 90)),
        daily_args = list(temperature = .ct_temperature(), ndraws = 2L,
                          seed = 1L, start = as.Date("2022-01-01"),
                          end = as.Date("2022-01-10")),
        focal_dates = as.Date("2022-01-10"),
        contributions = TRUE
      ))),
      class = "hydroSense_contributions_draws"
    )
  })

  it("validates the flag like other logical arguments", {
    sc <- .pf_scenario()
    for (bad in list("yes", NA, c(TRUE, FALSE))) {
      err <- expect_error(suppressMessages(mspaf_pipeline(
        .ct_daily_target(), reference = sc$reference, impute = FALSE,
        contributions = bad
      )))
      expect_match(conditionMessage(err), "contributions")
      expect_no_match(conditionMessage(err), "unused argument")
    }
  })

  it("also_total = TRUE: a `variant` column, each variant summing to its value", {
    skip_if_not("also_total" %in% names(formals(mspaf_pipeline)),
                "mspaf_pipeline() has no also_total argument yet")
    out <- .ct_pipeline(TRUE, contributions = TRUE, also_total = TRUE)
    ctb <- as.data.frame(attr(out, "contributions"))
    expect_true("variant" %in% names(ctb))
    expect_setequal(unique(ctb$variant), c("ara", "total"))
    expect_false(anyDuplicated(ctb[, c("variant", "site_id", "focal_date",
                                       "analyte")]) > 0)
    .ct_expect_chronic_sum(out, ctb[ctb$variant == "ara", ], "value")
    .ct_expect_chronic_sum(out, ctb[ctb$variant == "total", ], "value_total")
  })
})
