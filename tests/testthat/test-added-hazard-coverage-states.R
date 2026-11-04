## Added-hazard msPAF: per-toxicant coverage states
## ============================================================================
##
## CONTRACT (condensed; all data synthetic).
##
## mspaf_pipeline(ara_form = "added_hazard") takes daily_args$measured_within
## (a single non-negative number of days, may be Inf; default 30) and returns,
## on both the chronic (chronic = TRUE) and the daily (chronic = FALSE)
## output, attribute "analyte_coverage": columns site_id, side, analyte, state,
## date_from, date_to (character x4, Date x2).  Per (site_id, side, analyte)
## the rows are the maximal runs of constant state over the whole daily grid
## (contiguous, no overlap, first date_from = grid start, last date_to = grid
## end, adjacent states differ, NA is a state).
##
## State of a point / toxicant / day:
##   own normalised value not NA: "measured" if anchor_distance <=
##     measured_within, else "modelled" (interpolated or held; own series
##     only, never the gap-filler);
##   reference only, own value NA but static gap-filler not NA:
##     "background_fill";
##   otherwise NA.
## The target side never has "background_fill".
## Ordering: site_id, side, analyte, date_from (C locale).
## The setting changes no other output and is ignored by ara_form = "excess".
##
## Fixtures use B, Cr and Mn (no co-analytes in their normalisation) so the
## days of each state can be worked out by hand; Cu (needs DOC) is used for the
## co-analyte cases.  Grid of every fixture below: day 0 .. 50 (first sample
## day to last target sample day 40 + span_cap 10).

library(testthat)

.cs_cache <- new.env()
.cs_get <- function(key, f) {
  if (is.null(.cs_cache[[key]])) .cs_cache[[key]] <- f()
  .cs_cache[[key]]
}

.cs_run <- function(target, reference = NULL, ..., daily_args = list()) {
  da <- utils::modifyList(list(span_cap = 10, measured_within = 2), daily_args,
                          keep.null = TRUE)
  args <- list(target = target, reference = reference, impute = FALSE,
               ara_form = "added_hazard", daily_args = da)
  args <- utils::modifyList(args, list(...), keep.null = TRUE)
  suppressWarnings(suppressMessages(do.call(mspaf_pipeline, args)))
}

.cs_attr <- function(x) attr(x, "analyte_coverage", exact = TRUE)

## Expected frame from a compact spec: list of c(state, from_day, to_day).
.cs_exp <- function(site, side, analyte, runs) {
  data.frame(
    site_id = site, side = side, analyte = analyte,
    state = vapply(runs, function(r) as.character(r[[1L]]), character(1)),
    date_from = .ah_day(vapply(runs, function(r) as.numeric(r[[2L]]), 0)),
    date_to = .ah_day(vapply(runs, function(r) as.numeric(r[[3L]]), 0)),
    stringsAsFactors = FALSE
  )
}
.cs_runs <- function(...) list(...)
.cs_pick <- function(cov, site, side, analyte) {
  x <- as.data.frame(cov)
  x <- x[x$site_id == site & x$side == side & x$analyte == analyte, , drop = FALSE]
  rownames(x) <- NULL
  x
}
.cs_expect_runs <- function(cov, site, side, analyte, runs) {
  got <- .cs_pick(cov, site, side, analyte)
  exp <- .cs_exp(site, side, analyte, runs)
  expect_equal(got, exp, ignore_attr = TRUE, label = paste(site, side, analyte))
}

## Target T1: B at 0, 8 (below detection), 40; Cr at 0, 40.
.cs_t1 <- function(bdet = c(TRUE, TRUE, TRUE)) {
  .ah_rows("T1", list(
    B = list(day = c(0, 8, 40), value = c(3000, 4000, 5000), detected = bdet),
    Cr = list(day = c(0, 40), value = c(50, 80))
  ))
}
## Target T2: B at 10 only.
.cs_t2 <- function() .ah_rows("T2", list(B = list(day = 10, value = 3500)))
## Reference: B at 4, 12, 30; Mn at 5, 15; never Cr.
.cs_ref <- function() {
  r <- .ah_rows("REF", list(
    B = list(day = c(4, 12, 30), value = c(500, 600, 700)),
    Mn = list(day = c(5, 15), value = c(300, 400))
  ))
  r$units.analyte <- "ug/L"   # the helper only sets units for its own toxicants
  r
}

## ═════════════════════════════════════════════════════════════════════════════
## Validation
## ═════════════════════════════════════════════════════════════════════════════

describe("daily_args$measured_within: validation", {
  tg <- .cs_t1()

  it("accepts 0, a positive number and Inf", {
    for (v in list(0, 7, 30.5, Inf)) {
      expect_no_error(.cs_run(tg, daily_args = list(measured_within = v)))
    }
  })

  it("errors on negative, NA, non-numeric or non-scalar values", {
    for (v in list(-1, NA, NA_real_, "30", c(1, 2), numeric(0), -Inf)) {
      expect_error(.cs_run(tg, daily_args = list(measured_within = v)),
                   regexp = "measured_within",
                   label = paste(deparse(v), collapse = ""))
    }
  })

  it("validates before any model is fitted or anything is computed", {
    reg <- .cs_ref()
    testthat::local_mocked_bindings(
      fit_reference_model = function(...) stop("fit_reference_model was called"),
      .package = "hydroSense"
    )
    expect_error(.cs_run(tg, reg, daily_args = list(measured_within = -1)),
                 regexp = "measured_within")
  })
})

## ═════════════════════════════════════════════════════════════════════════════
## Structure
## ═════════════════════════════════════════════════════════════════════════════

describe("analyte_coverage: presence, columns and types", {
  tg <- .cs_t1()
  chk <- function(out) {
    cov <- .cs_attr(out)
    expect_false(is.null(cov))
    expect_s3_class(cov, "data.frame")
    expect_named(cov, c("site_id", "side", "analyte", "state", "date_from",
                        "date_to"))
    expect_type(cov$site_id, "character")
    expect_type(cov$side, "character")
    expect_type(cov$analyte, "character")
    expect_type(cov$state, "character")
    expect_s3_class(cov$date_from, "Date")
    expect_s3_class(cov$date_to, "Date")
    expect_true(all(cov$state %in% c("measured", "modelled", "background_fill", NA)))
    expect_true(all(cov$side %in% c("target", "reference")))
  }

  it("is on the chronic output (no reference)", {
    chk(.cs_run(tg, chronic = TRUE))
  })
  it("is on the daily output (no reference)", {
    out <- .cs_run(tg, chronic = FALSE)
    chk(out)
  })
  it("is on both outputs with a reference", {
    chk(.cs_run(tg, .cs_ref(), chronic = TRUE))
    chk(.cs_run(tg, .cs_ref(), chronic = FALSE))
  })
  it("is the same on the chronic and the daily output", {
    a <- .cs_attr(.cs_run(tg, .cs_ref(), chronic = TRUE))
    b <- .cs_attr(.cs_run(tg, .cs_ref(), chronic = FALSE))
    expect_equal(as.data.frame(a), as.data.frame(b), ignore_attr = TRUE)
  })
})

describe("analyte_coverage: partition properties", {
  out <- .cs_get("part", function() {
    .cs_run(dplyr::bind_rows(.cs_t1(), .cs_t2()), .cs_ref())
  })
  cov <- as.data.frame(.cs_attr(out))
  grid <- as.Date(attr(out, "daily")$date)
  grp <- split(cov, list(cov$site_id, cov$side, cov$analyte), drop = TRUE)

  it("has runs that are contiguous, ordered and cover the grid exactly", {
    expect_gt(length(grp), 0L)
    for (g in grp) {
      expect_true(all(g$date_from <= g$date_to))
      expect_identical(g$date_from[1L], min(grid))
      expect_identical(g$date_to[nrow(g)], max(grid))
      if (nrow(g) > 1L) {
        expect_identical(g$date_from[-1L], g$date_to[-nrow(g)] + 1)
      }
    }
  })

  it("has maximal runs: adjacent runs differ in state (NA is a state)", {
    expect_gt(length(grp), 0L)
    for (g in grp) {
      if (nrow(g) > 1L) {
        s <- g$state
        expect_true(all(is.na(s[-1L]) != is.na(s[-length(s)]) |
                          (!is.na(s[-1L]) & !is.na(s[-length(s)]) &
                             s[-1L] != s[-length(s)])))
      }
    }
  })

  it("is ordered by site_id, side, analyte, date_from (C locale)", {
    o <- order(cov$site_id, cov$side, cov$analyte, cov$date_from,
               method = "radix")
    expect_identical(o, seq_len(nrow(cov)))
  })

  it("has one set of rows per site, side and toxicant of the run", {
    expect_setequal(unique(cov$site_id), c("T1", "T2"))
    expect_setequal(unique(cov$side), c("target", "reference"))
    expect_setequal(unique(cov$analyte), c("B", "Cr", "Mn"))
    expect_length(grp, 2L * 2L * 3L)
  })
})

## ═════════════════════════════════════════════════════════════════════════════
## Hand-computed states: target only
## ═════════════════════════════════════════════════════════════════════════════
## span_cap 10, measured_within 2.
## B (T1) anchors 0, 8, 40.  0..8 covered (gap 8): distances 0 1 2 3 4 3 2 1 0;
##   9..39 uncovered (gap 32); 40..50 held (distance 0..10).
## Cr (T1) anchors 0, 40: gap 40 > cap, so 1..39 uncovered (only the last
##   anchor is held); 40..50 held.

describe("analyte_coverage: states of a target alone", {
  out <- .cs_get("tonly", function() .cs_run(.cs_t1()))
  cov <- .cs_attr(out)

  it("gives measured / modelled / NA runs for an interpolated toxicant", {
    .cs_expect_runs(cov, "T1", "target", "B", .cs_runs(
      c("measured", 0, 2), c("modelled", 3, 5), c("measured", 6, 8),
      c(NA, 9, 39), c("measured", 40, 42), c("modelled", 43, 50)))
  })

  it("holds after the last anchor and is NA before/beyond the span cap", {
    .cs_expect_runs(cov, "T1", "target", "Cr", .cs_runs(
      c("measured", 0, 0), c(NA, 1, 39),
      c("measured", 40, 42), c("modelled", 43, 50)))
  })

  it("has no reference rows without a reference", {
    expect_false("reference" %in% cov$side)
  })

  it("treats distance == measured_within as measured and +1 as modelled", {
    # [decision] B day 2 (distance 2) is measured, day 3 (distance 3) modelled
    # with measured_within = 2: asserted by the runs above; here with 3.
    o3 <- .cs_run(.cs_t1(), daily_args = list(measured_within = 3))
    .cs_expect_runs(.cs_attr(o3), "T1", "target", "B", .cs_runs(
      c("measured", 0, 3), c("modelled", 4, 4), c("measured", 5, 8),
      c(NA, 9, 39), c("measured", 40, 43), c("modelled", 44, 50)))
  })

  it("measured_within = 0 makes only anchor days measured", {
    o0 <- .cs_run(.cs_t1(), daily_args = list(measured_within = 0))
    .cs_expect_runs(.cs_attr(o0), "T1", "target", "B", .cs_runs(
      c("measured", 0, 0), c("modelled", 1, 7), c("measured", 8, 8),
      c(NA, 9, 39), c("measured", 40, 40), c("modelled", 41, 50)))
  })

  it("measured_within = Inf makes every covered day measured", {
    oi <- .cs_run(.cs_t1(), daily_args = list(measured_within = Inf))
    .cs_expect_runs(.cs_attr(oi), "T1", "target", "B", .cs_runs(
      c("measured", 0, 8), c(NA, 9, 39), c("measured", 40, 50)))
  })

  it("defaults measured_within to 30", {
    # anchors 0 and 40, span_cap Inf: day 20 is 20 from either (measured),
    # with a default of 30 all days of the 0..40 stretch are within 20.
    tg <- .ah_rows("T1", list(B = list(day = c(0, 40), value = c(1000, 2000))))
    d <- .cs_run(tg, daily_args = list(span_cap = 100, measured_within = NULL))
    .cs_expect_runs(.cs_attr(d), "T1", "target", "B", .cs_runs(
      c("measured", 0, 70), c("modelled", 71, 140)))
  })

  it("counts below-detection anchors as results", {
    ob <- .cs_run(.cs_t1(bdet = c(FALSE, FALSE, FALSE)))
    .cs_expect_runs(.cs_attr(ob), "T1", "target", "B", .cs_runs(
      c("measured", 0, 2), c("modelled", 3, 5), c("measured", 6, 8),
      c(NA, 9, 39), c("measured", 40, 42), c("modelled", 43, 50)))
  })

  it("is NA throughout when a required co-analyte is missing", {
    tg <- dplyr::bind_rows(
      .ah_rows("T1", list(Cu = list(day = c(0, 4), value = c(5, 6)))),
      .ah_rows("T1", list(B = list(day = c(0, 40), value = c(1000, 2000))))
    )
    o <- .cs_run(tg)  # no DOC: Cu cannot be normalised; grid 0..50
    .cs_expect_runs(.cs_attr(o), "T1", "target", "Cu",
                    .cs_runs(c(NA, 0, 50)))
  })

  it("is measured/modelled for a normalised toxicant once its co-analyte exists", {
    tg <- .ah_rows("T1", list(
      Cu = list(day = c(0, 4), value = c(5, 6)),
      DOC = list(day = 0, value = 1)
    ))
    o <- .cs_run(tg, daily_args = list(measured_within = 1))  # grid 0..14
    .cs_expect_runs(.cs_attr(o), "T1", "target", "Cu", .cs_runs(
      c("measured", 0, 1), c("modelled", 2, 2), c("measured", 3, 5),
      c("modelled", 6, 14)))
  })
})

## ═════════════════════════════════════════════════════════════════════════════
## Hand-computed states: with a reference
## ═════════════════════════════════════════════════════════════════════════════
## Reference B anchors 4, 12, 30 (cap 10, mw 2): own distances
##   4..12: 0 1 2 3 4 3 2 1 0; 13..29 uncovered; 30..40 held (0..10);
##   0..3 and 41..50 uncovered.  The static gap-filler (B was measured) turns
##   every uncovered day into "background_fill".
## Reference Mn anchors 5, 15 (gap 10): 5..15 distances 0 1 2 3 4 5 4 3 2 1 0;
##   16..25 held (1..10, own "modelled"); the rest uncovered -> gap-filled
##   "background_fill".
## Reference Cr: never measured -> NA over the whole grid.
## Target Mn: never measured -> NA over the whole grid.

describe("analyte_coverage: states with a reference", {
  out <- .cs_get("withref", function() {
    .cs_run(dplyr::bind_rows(.cs_t1(), .cs_t2()), .cs_ref())
  })
  cov <- .cs_attr(out)

  it("gives the reference's own measured/modelled days and gap-fills the rest", {
    .cs_expect_runs(cov, "T1", "reference", "B", .cs_runs(
      c("background_fill", 0, 3), c("measured", 4, 6), c("modelled", 7, 9),
      c("measured", 10, 12), c("background_fill", 13, 29),
      c("measured", 30, 32), c("modelled", 33, 40),
      c("background_fill", 41, 50)))
  })

  it("marks static gap-fill days 'background_fill' even where the target is uncovered", {
    # target T1 Mn is NA everywhere, reference Mn is gap-filled / own.
    # modelled (held, own series) 18..25 then background_fill from 26.
    .cs_expect_runs(cov, "T1", "reference", "Mn", .cs_runs(
      c("background_fill", 0, 4), c("measured", 5, 7), c("modelled", 8, 12),
      c("measured", 13, 17), c("modelled", 18, 25),
      c("background_fill", 26, 50)))
  })

  it("keeps interpolated and held reference days 'modelled', not background_fill", {
    r <- .cs_pick(cov, "T1", "reference", "B")
    state_on <- function(i) r$state[r$date_from <= .ah_day(i) &
                                      .ah_day(i) <= r$date_to]
    for (i in c(7:9, 33:40)) expect_identical(state_on(i), "modelled")
    for (i in c(0:3, 13:29, 41:50)) expect_identical(state_on(i), "background_fill")
  })

  it("never has background_fill on the target side", {
    expect_false("background_fill" %in% cov$state[cov$side == "target"])
    expect_true("background_fill" %in% cov$state[cov$side == "reference"])
  })

  it("splits measured -> modelled -> background_fill into separate maximal runs", {
    r <- .cs_pick(cov, "T1", "reference", "Mn")
    i <- which(r$state == "modelled" & r$date_from == .ah_day(18))
    expect_length(i, 1L)
    expect_identical(r$state[i + 1L], "background_fill")
    expect_identical(r$state[i - 1L], "measured")
    expect_identical(r$date_from[i + 1L], r$date_to[i] + 1)
  })

  it("is NA on the reference side for a toxicant it never measured", {
    .cs_expect_runs(cov, "T1", "reference", "Cr", .cs_runs(c(NA, 0, 50)))
  })

  it("is all NA for the target for a toxicant only the reference has", {
    .cs_expect_runs(cov, "T1", "target", "Mn", .cs_runs(c(NA, 0, 50)))
    .cs_expect_runs(cov, "T2", "target", "Mn", .cs_runs(c(NA, 0, 50)))
  })

  it("keeps target-side runs identical to the no-reference run", {
    .cs_expect_runs(cov, "T1", "target", "B", .cs_runs(
      c("measured", 0, 2), c("modelled", 3, 5), c("measured", 6, 8),
      c(NA, 9, 39), c("measured", 40, 42), c("modelled", 43, 50)))
  })

  it("repeats the reference rows under each target site_id", {
    r1 <- .cs_pick(cov, "T1", "reference", "B")
    r2 <- .cs_pick(cov, "T2", "reference", "B")
    expect_gt(nrow(r2), 0L)
    expect_equal(r1[, -1], r2[, -1], ignore_attr = TRUE)
    for (a in c("Cr", "Mn")) {
      expect_equal(.cs_pick(cov, "T1", "reference", a)[, -1],
                   .cs_pick(cov, "T2", "reference", a)[, -1],
                   ignore_attr = TRUE)
    }
  })

  it("gives the second target its own runs", {
    .cs_expect_runs(cov, "T2", "target", "B", .cs_runs(
      c(NA, 0, 9), c("measured", 10, 12), c("modelled", 13, 20),
      c(NA, 21, 50)))
    .cs_expect_runs(cov, "T2", "target", "Cr", .cs_runs(c(NA, 0, 50)))
  })

  it("has the same set of toxicants on both sides", {
    expect_setequal(cov$analyte[cov$side == "reference"],
                    cov$analyte[cov$side == "target"])
  })

  it("is NA on the reference side for a toxicant with no reference result at all", {
    o <- .cs_run(.cs_t1(), .ah_rows("REF", list(
      B = list(day = c(4, 12), value = c(500, 600)))))
    .cs_expect_runs(.cs_attr(o), "T1", "reference", "Cr", .cs_runs(c(NA, 0, 50)))
  })

  it("counts reference below-detection anchors as results", {
    ref <- .ah_rows("REF", list(
      B = list(day = c(4, 12, 30), value = c(500, 600, 700),
               detected = c(FALSE, FALSE, FALSE))))
    o <- .cs_run(.cs_t1(), ref)
    .cs_expect_runs(.cs_attr(o), "T1", "reference", "B", .cs_runs(
      c("background_fill", 0, 3), c("measured", 4, 6), c("modelled", 7, 9),
      c("measured", 10, 12), c("background_fill", 13, 29),
      c("measured", 30, 32), c("modelled", 33, 40),
      c("background_fill", 41, 50)))
  })

  it("never labels the target background_fill when the target is below detection", {
    o <- .cs_run(.cs_t1(bdet = c(FALSE, FALSE, FALSE)), .cs_ref())
    cv <- .cs_attr(o)
    expect_false("background_fill" %in% cv$state[cv$side == "target"])
    .cs_expect_runs(cv, "T1", "target", "B", .cs_runs(
      c("measured", 0, 2), c("modelled", 3, 5), c("measured", 6, 8),
      c(NA, 9, 39), c("measured", 40, 42), c("modelled", 43, 50)))
  })

  it("gives the same states on the daily output", {
    od <- .cs_run(dplyr::bind_rows(.cs_t1(), .cs_t2()), .cs_ref(),
                  chronic = FALSE)
    cd <- .cs_attr(od)
    .cs_expect_runs(cd, "T1", "reference", "B", .cs_runs(
      c("background_fill", 0, 3), c("measured", 4, 6), c("modelled", 7, 9),
      c("measured", 10, 12), c("background_fill", 13, 29),
      c("measured", 30, 32), c("modelled", 33, 40),
      c("background_fill", 41, 50)))
    .cs_expect_runs(cd, "T2", "reference", "Cr", .cs_runs(c(NA, 0, 50)))
    expect_false("background_fill" %in% cd$state[cd$side == "target"])
  })

  it("is background_fill only where the own series is absent (measured_within = Inf)", {
    oi <- .cs_run(.cs_t1(), .cs_ref(), daily_args = list(measured_within = Inf))
    .cs_expect_runs(.cs_attr(oi), "T1", "reference", "B", .cs_runs(
      c("background_fill", 0, 3), c("measured", 4, 12),
      c("background_fill", 13, 29), c("measured", 30, 40),
      c("background_fill", 41, 50)))
  })

  it("has no background_fill anywhere without a reference", {
    o <- .cs_run(dplyr::bind_rows(.cs_t1(), .cs_t2()))
    expect_false("background_fill" %in% .cs_attr(o)$state)
    o2 <- .cs_run(dplyr::bind_rows(.cs_t1(), .cs_t2()), chronic = FALSE)
    expect_false("background_fill" %in% .cs_attr(o2)$state)
  })

  it("leaves the numeric outputs unchanged by the background_fill label", {
    # Values worked out on the behaviour before the relabel (the gap-filler
    # feeds the reference exactly as before); only the state label changes.
    o <- .cs_run(.cs_t1(), .cs_ref(), also_total = TRUE)
    d <- attr(o, "daily")
    i <- c(1, 5, 10, 41, 51)   # grid days 0, 4, 9, 40, 50
    expect_equal(d$mspaf[i], c(24.837654573370592, 29.564511669289161, NA,
                               40.40577505056244, 40.935032671569495),
                 tolerance = 1e-8)
    expect_equal(d$mspaf_total[i], c(26.575739640402496, 30.70628129435644, NA,
                                     42.300875330091699, 42.300875330091699),
                 tolerance = 1e-8)
    # the reference msPAF on a gap-filled day equals the static summary's
    expect_equal(d$mspaf_reference[i], c(2.312441232596385, 1.6210147073963732,
                                         NA, 3.1800065881168127,
                                         2.312441232596385), tolerance = 1e-8)
    expect_identical(d$n_analytes[i], c(1L, 1L, 0L, 1L, 1L))
  })
})

## ═════════════════════════════════════════════════════════════════════════════
## Edge cases and invariances
## ═════════════════════════════════════════════════════════════════════════════

describe("analyte_coverage: no toxicant", {
  it("is a zero-row frame with the six columns and types", {
    tg <- .ah_rows("T1", list(pH = list(day = c(0, 10), value = c(7, 7.2))))
    out <- tryCatch(.cs_run(tg), error = function(e) e)
    expect_false(inherits(out, "error"))
    cov <- .cs_attr(out)
    expect_equal(nrow(cov), 0L)
    expect_named(cov, c("site_id", "side", "analyte", "state", "date_from",
                        "date_to"))
    expect_type(cov$site_id, "character")
    expect_type(cov$side, "character")
    expect_type(cov$analyte, "character")
    expect_type(cov$state, "character")
    expect_s3_class(cov$date_from, "Date")
    expect_s3_class(cov$date_to, "Date")
  })
})

describe("daily_args$measured_within: changes no other output", {
  tgt <- .ah_scn_target()
  ref <- .ah_scn_REF()
  a <- .cs_get("inv_a", function() {
    .cs_run(tgt, ref, also_total = TRUE, contributions = TRUE,
            daily_args = list(span_cap = 120, measured_within = 0))
  })
  b <- .cs_get("inv_b", function() {
    .cs_run(tgt, ref, also_total = TRUE, contributions = TRUE,
            daily_args = list(span_cap = 120, measured_within = Inf))
  })

  it("leaves value, value_total and coverage identical", {
    for (col in c("value", "value_total", "coverage")) {
      expect_identical(as.data.frame(a)[[col]], as.data.frame(b)[[col]])
    }
  })
  it("leaves contributions, daily and sample_days identical", {
    for (nm in c("contributions", "daily", "sample_days")) {
      expect_equal(attr(a, nm), attr(b, nm))
    }
  })
  it("changes the coverage states themselves", {
    expect_false(isTRUE(all.equal(.cs_attr(a), .cs_attr(b))))
  })
  it("works on the daily output too", {
    da <- .cs_run(tgt, ref, chronic = FALSE,
                  daily_args = list(span_cap = 120, measured_within = 0))
    db <- .cs_run(tgt, ref, chronic = FALSE,
                  daily_args = list(span_cap = 120, measured_within = Inf))
    expect_equal(as.data.frame(da), as.data.frame(db), ignore_attr = "analyte_coverage")
    expect_equal(attr(da, "contributions"), attr(db, "contributions"))
    expect_false(isTRUE(all.equal(.cs_attr(da), .cs_attr(db))))
  })
})

describe("ara_form = 'excess': measured_within is ignored", {
  it("neither errors nor changes the result, whatever the value", {
    skip_on_cran()
    sc <- .pf_scenario()
    tgt <- .pf_small_target()
    tgt <- tgt[tgt$site_id == "P1", , drop = FALSE]
    tmp <- sc$temperature
    tmp <- tmp[tmp$datetime >= as.Date("2021-12-15") &
                 tmp$datetime <= as.Date("2022-03-10"), ]
    run <- function(mw) {
      da <- list(temperature = tmp, start = as.Date("2022-01-01"),
                 end = as.Date("2022-02-28"))
      if (!is.null(mw)) da$measured_within <- mw
      suppressWarnings(suppressMessages(mspaf_pipeline(
        tgt, NULL, impute = FALSE, ara_form = "excess", daily_args = da,
        focal_dates = as.Date(c("2022-01-20", "2022-02-28")),
        tau = 90, window = 365)))
    }
    base <- run(NULL)
    for (v in list(0, 5, Inf, -3, "x")) {
      other <- run(v)
      attributes(base)[c("daily", "reference_model", "imputation_model")] <- NULL
      attributes(other)[c("daily", "reference_model", "imputation_model")] <- NULL
      expect_equal(other, base, tolerance = 1e-12)
    }
    expect_null(attr(run(NULL), "analyte_coverage"))
  })
})
