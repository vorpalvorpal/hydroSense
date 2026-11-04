## Added-hazard msPAF: daily_args$fill_reference
## ============================================================================
##
## CONTRACT (condensed; all data synthetic).
##
## mspaf_pipeline(ara_form = "added_hazard", daily_args = list(fill_reference =
## <data frame shaped like reference>)): the reference's static gap-filler
## summary (used on "background_fill" days, where the target is covered and the
## reference's own series is not) is computed per analyte from the
## fill_reference rows instead of the reference rows.  Nothing else changes:
## the reference's own series, the grid, the coverage states and the target
## all still come from `reference`.  An analyte with no fill_reference rows has
## no fill.  NULL / absent = today's behaviour.  fill_reference with another
## site_id than the reference's is an error, as is a non-data-frame.  Ignored
## when there is no reference.
##
## Fixture (B only unless stated; B has no co-analyte normalisation):
##   target T1: B = 5000 every 5 days, day 0..40 (covered on the whole grid
##     0..50 with span_cap 10);
##   reference REF: B at days 0, 2, 4 (4000, 2000, 1000), so its own series
##     covers days 0..14 and days 15..50 are "background_fill".
##   old static summary = geometric mean of 4000, 2000, 1000 = 2000.

library(testthat)

.fr_run <- function(target, reference, fill = NULL, chronic = FALSE, ...) {
  da <- list(span_cap = 10)
  if (!is.null(fill)) da$fill_reference <- fill
  suppressWarnings(suppressMessages(
    mspaf_pipeline(target, reference, impute = FALSE, chronic = chronic,
                   ara_form = "added_hazard", daily_args = da, ...)
  ))
}

.fr_target <- function() {
  .ah_rows("T1", list(B = list(day = seq(0, 40, by = 5), value = 5000)))
}
.fr_ref <- function() {
  .ah_rows("REF", list(B = list(day = c(0, 2, 4), value = c(4000, 2000, 1000))))
}
## A lower subset-like fill set (not rows of the reference).
.fr_fill <- function(vals = c(600, 400), detected = TRUE) {
  .ah_rows("REF", list(B = list(day = c(100, 110)[seq_along(vals)],
                                value = vals, detected = detected)))
}

## Expected daily mspaf (percent) of a fill day for a static B of `s`.
.fr_exp_D <- function(s, tgt = 5000) {
  100 * .ah_D(.ah_ms1(c(B = tgt)), .ah_ms1(c(B = unname(s))))
}
.fr_fill_days <- function(d) d$date >= .ah_day(15) & d$date <= .ah_day(50)

describe("daily_args$fill_reference: fill-day values", {
  it("uses the fill_reference summary on background_fill days", {
    fill <- .fr_fill()
    base <- .fr_run(.fr_target(), .fr_ref())
    got <- .fr_run(.fr_target(), .fr_ref(), fill)
    s_old <- .ah_static(.ah_reduce(.fr_ref()), "REF", "B")
    s_new <- .ah_static(.ah_reduce(fill), "REF", "B")
    expect_equal(unname(s_old), 2000)
    expect_equal(unname(s_new), sqrt(600 * 400))
    fd <- .fr_fill_days(got)
    expect_true(all(abs(base$mspaf[fd] - .fr_exp_D(s_old)) < 1e-6))
    expect_true(all(abs(got$mspaf[fd] - .fr_exp_D(s_new)) < 1e-6))
    expect_true(all(abs(got$mspaf_reference[fd] -
                          100 * .ah_ms1(c(B = unname(s_new)))) < 1e-6))
    ## lower background => more added hazard
    expect_true(all(got$mspaf[fd] > base$mspaf[fd]))
  })

  it("gives zero added hazard when the fill summary equals the target", {
    fill <- .ah_rows("REF", list(B = list(day = c(100, 110), value = 5000)))
    got <- .fr_run(.fr_target(), .fr_ref(), fill)
    fd <- .fr_fill_days(got)
    expect_true(all(abs(got$mspaf[fd]) < 1e-9))
  })

  it("applies the below-detection rules to the fill_reference rows", {
    fill <- .ah_rows("REF", list(B = list(day = c(100, 110, 120),
                                          value = c(600, 400, 800),
                                          detected = c(FALSE, TRUE, FALSE))))
    s_new <- .ah_static(.ah_reduce(fill), "REF", "B")
    got <- .fr_run(.fr_target(), .fr_ref(), fill)
    fd <- .fr_fill_days(got)
    expect_true(all(abs(got$mspaf[fd] - .fr_exp_D(s_new)) < 1e-6))
  })

  it("moves the chronic value the same way", {
    fill <- .fr_fill()
    base <- .fr_run(.fr_target(), .fr_ref(), chronic = TRUE,
                    focal_dates = .ah_day(c(40, 50)))
    got <- .fr_run(.fr_target(), .fr_ref(), fill, chronic = TRUE,
                   focal_dates = .ah_day(c(40, 50)))
    expect_true(all(got$value > base$value))
  })
})

describe("daily_args$fill_reference: what stays the same", {
  it("leaves days with a reference value of its own unchanged", {
    base <- .fr_run(.fr_target(), .fr_ref())
    got <- .fr_run(.fr_target(), .fr_ref(), .fr_fill())
    own <- !.fr_fill_days(got)
    expect_equal(as.data.frame(got)[own, ], as.data.frame(base)[own, ])
  })

  it("leaves the coverage states unchanged", {
    base <- .fr_run(.fr_target(), .fr_ref())
    got <- .fr_run(.fr_target(), .fr_ref(), .fr_fill())
    expect_equal(attr(got, "analyte_coverage"), attr(base, "analyte_coverage"))
    st <- as.data.frame(attr(got, "analyte_coverage"))
    st <- st[st$side == "reference" & st$analyte == "B", ]
    expect_true("background_fill" %in% st$state)
  })

  it("is identical to no fill_reference when fill_reference equals reference", {
    ref <- .fr_ref()
    base <- .fr_run(.fr_target(), ref)
    got <- .fr_run(.fr_target(), ref, ref)
    expect_equal(got, base)
    expect_equal(attr(got, "analyte_coverage"), attr(base, "analyte_coverage"))
  })
})

describe("daily_args$fill_reference: analytes without fill rows", {
  it("gives no fill for an analyte absent from fill_reference", {
    tg <- dplyr::bind_rows(
      .fr_target(),
      .ah_rows("T1", list(Cr = list(day = seq(0, 40, by = 5), value = 100)))
    )
    rf <- dplyr::bind_rows(
      .fr_ref(),
      .ah_rows("REF", list(Cr = list(day = c(0, 2, 4), value = 20)))
    )
    base <- .fr_run(tg, rf)
    got <- .fr_run(tg, rf, .fr_fill())   # fill has B only
    fd <- .fr_fill_days(got)
    ## base: both analytes on fill days; with B-only fill, Cr drops out
    expect_true(all(base$n_analytes[fd & base$date <= .ah_day(50)] == 2L))
    expect_true(all(got$n_analytes[fd] == 1L))
    s_new <- .ah_static(.ah_reduce(.fr_fill()), "REF", "B")
    expect_true(all(abs(got$mspaf[fd] - .fr_exp_D(s_new)) < 1e-6))
    ## reference Cr has no background_fill state any more
    st <- as.data.frame(attr(got, "analyte_coverage"))
    st <- st[st$side == "reference" & st$analyte == "Cr", ]
    expect_false("background_fill" %in% st$state)
  })
})

describe("daily_args$fill_reference: validation and no reference", {
  it("errors on a fill_reference from another site", {
    bad <- .fr_fill()
    bad$site_id <- "OTHER"
    expect_error(.fr_run(.fr_target(), .fr_ref(), bad))
    mixed <- dplyr::bind_rows(.fr_fill(), bad)
    expect_error(.fr_run(.fr_target(), .fr_ref(), mixed))
  })

  it("errors on a fill_reference that is not a data frame", {
    expect_error(.fr_run(.fr_target(), .fr_ref(), "REF"))
    expect_error(.fr_run(.fr_target(), .fr_ref(), c(B = 500)))
  })

  it("is ignored when there is no reference", {
    base <- .fr_run(.fr_target(), NULL)
    got <- .fr_run(.fr_target(), NULL, .fr_fill())
    expect_equal(got, base)
  })
})
