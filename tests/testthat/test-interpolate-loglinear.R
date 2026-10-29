## interpolate_loglinear(): log-linear daily interpolation with a span cap
## ============================================================================
##
## CONTRACT.
##
##   interpolate_loglinear(dates, values, detected, grid, span_cap = 120,
##                         log = TRUE)
##
## * dates (Date, unique, any order), values (numeric), detected (logical;
##   FALSE = below detection, `values` then holds the anchor value to use,
##   e.g. the detection limit), grid (Date vector).
## * Returns a data.frame with columns `date` (= grid), `value`, `detected`,
##   one row per grid day, in grid order.
## * log = TRUE drops anchors with NA or non-positive value; log = FALSE drops
##   only NA anchors.  Duplicated dates: error of class `hydroSense_input`.
## * For a grid day g with x_i the last anchor date <= g:
##     - no anchor on or before g: value NA, detected NA;
##     - g equals an anchor date: that anchor's value and detected flag;
##     - g after the last anchor: covered iff g - x_last <= span_cap; value and
##       flag held from the last anchor;
##     - otherwise, with the next anchor x_j: covered iff x_j - x_i <= span_cap;
##       value linear in log(value) between the two anchors (in value when
##       log = FALSE); detected is FALSE only when BOTH bracketing anchors are
##       below detection, else TRUE;
##     - uncovered days: value NA and detected NA.
## * span_cap = Inf disables the cap.  No anchors: all NA.
##
## All data are synthetic.

library(testthat)

.il_d0 <- as.Date("2024-01-01")
.il <- function(off, v, det = rep(TRUE, length(off)), grid_off, ...) {
  interpolate_loglinear(.il_d0 + off, v, det, .il_d0 + grid_off, ...)
}

describe("interpolate_loglinear(): return shape", {
  it("returns a data.frame with date, value, detected, one row per grid day", {
    out <- .il(c(0, 10), c(2, 8), grid_off = -1:12)
    expect_s3_class(out, "data.frame")
    expect_named(out, c("date", "value", "detected"))
    expect_equal(nrow(out), 14L)
    expect_identical(out$date, .il_d0 + (-1:12))
    expect_type(out$value, "double")
    expect_type(out$detected, "logical")
  })

  it("keeps grid order when the grid is not sorted, and accepts unsorted anchors", {
    g <- c(5, 0, 12, -1, 10)
    out <- .il(c(10, 0), c(8, 2), grid_off = g)
    expect_identical(out$date, .il_d0 + g)
    expect_equal(out$value, c(4, 2, 8, NA, 8), tolerance = 1e-12)
  })
})

describe("interpolate_loglinear(): values", {
  it("returns the anchor value on an anchor date", {
    out <- .il(c(0, 10, 25), c(2, 8, 3), grid_off = c(0, 10, 25))
    expect_equal(out$value, c(2, 8, 3), tolerance = 1e-12)
    expect_true(all(out$detected))
  })

  it("interpolates linearly in log(value): midpoint is the geometric mean", {
    out <- .il(c(0, 10), c(2, 8), grid_off = 5)
    expect_equal(out$value, 4, tolerance = 1e-12)
    out2 <- .il(c(0, 10), c(2, 8), grid_off = c(2, 7))
    expect_equal(out2$value, exp(log(2) + c(0.2, 0.7) * (log(8) - log(2))),
                 tolerance = 1e-12)
  })

  it("interpolates linearly in value when log = FALSE", {
    out <- .il(c(0, 10), c(2, 8), grid_off = c(5, 2), log = FALSE)
    expect_equal(out$value, c(5, 3.2), tolerance = 1e-12)
  })

  it("is NA with detected NA before the first anchor", {
    out <- .il(c(3, 10), c(2, 8), grid_off = c(-5, 0, 2, 3))
    expect_true(all(is.na(out$value[1:3])))
    expect_true(all(is.na(out$detected[1:3])))
    expect_equal(out$value[4], 2)
  })
})

describe("interpolate_loglinear(): span cap", {
  it("covers an interior gap of exactly span_cap days but not span_cap + 1", {
    ## anchors at 0, 120 (gap 120) and 120, 241 (gap 121)
    out <- .il(c(0, 120, 241), c(1, 4, 16), grid_off = c(0, 60, 120, 180, 241))
    expect_equal(out$value[2], 2, tolerance = 1e-12)       # geometric mean
    expect_true(is.na(out$value[4]))                        # gap of 121 days
    expect_true(is.na(out$detected[4]))
    expect_equal(out$value[c(1, 3, 5)], c(1, 4, 16))        # anchors survive
  })

  it("respects a custom span_cap", {
    out <- .il(c(0, 40), c(1, 4), grid_off = 20, span_cap = 40)
    expect_equal(out$value, 2, tolerance = 1e-12)
    out <- .il(c(0, 40), c(1, 4), grid_off = 20, span_cap = 39)
    expect_true(is.na(out$value))
  })

  it("holds the last sample for span_cap days after it, then NA", {
    out <- .il(c(0, 10), c(2, 8), grid_off = c(10, 11, 130, 131))
    expect_equal(out$value[1:3], c(8, 8, 8))
    expect_true(is.na(out$value[4]))
    expect_true(is.na(out$detected[4]))
    out <- .il(c(0, 10), c(2, 8), grid_off = c(40, 41), span_cap = 30)
    expect_equal(out$value[1], 8)
    expect_true(is.na(out$value[2]))
  })

  it("disables the cap with span_cap = Inf", {
    out <- .il(c(0, 5000), c(1, 4), grid_off = c(2500, 9000), span_cap = Inf)
    expect_equal(out$value[1], 2, tolerance = 1e-12)
    expect_equal(out$value[2], 4)
    expect_true(out$detected[2])
  })
})

describe("interpolate_loglinear(): detection flags", {
  it("flags a day below detection only when BOTH bracketing anchors are", {
    out <- .il(c(0, 10), c(2, 8), det = c(FALSE, FALSE), grid_off = 5)
    expect_identical(out$detected, FALSE)
    out <- .il(c(0, 10), c(2, 8), det = c(FALSE, TRUE), grid_off = 5)
    expect_identical(out$detected, TRUE)
    out <- .il(c(0, 10), c(2, 8), det = c(TRUE, FALSE), grid_off = 5)
    expect_identical(out$detected, TRUE)
  })

  it("uses the anchor value of a below-detection anchor in the interpolation", {
    out <- .il(c(0, 10), c(2, 8), det = c(FALSE, TRUE), grid_off = 5)
    expect_equal(out$value, 4, tolerance = 1e-12)
  })

  it("reports an exact below-detection date as detected = FALSE with its value", {
    out <- .il(c(0, 10, 20), c(2, 0.5, 8), det = c(TRUE, FALSE, TRUE),
               grid_off = 10)
    expect_identical(out$detected, FALSE)
    expect_equal(out$value, 0.5)
  })

  it("holds the below-detection flag after the last anchor", {
    out <- .il(c(0, 10), c(2, 0.5), det = c(TRUE, FALSE), grid_off = c(10, 50))
    expect_identical(out$detected, c(FALSE, FALSE))
    expect_equal(out$value, c(0.5, 0.5))
  })
})

describe("interpolate_loglinear(): anchor dropping", {
  it("drops NA, zero and negative anchors when log = TRUE", {
    out <- .il(c(0, 3, 5, 7, 10), c(2, NA, 0, -1, 8), grid_off = c(3, 5, 7))
    ## remaining anchors: (0, 2) and (10, 8)
    expect_equal(out$value, exp(log(2) + c(0.3, 0.5, 0.7) * log(4)),
                 tolerance = 1e-12)
  })

  it("keeps zero and negative anchors when log = FALSE, drops only NA", {
    out <- .il(c(0, 3, 5, 10), c(2, NA, 0, 8), grid_off = c(3, 5, 8), log = FALSE)
    ## anchors: (0, 2), (5, 0), (10, 8)
    expect_equal(out$value, c(2 + (0 - 2) * 3 / 5, 0, 4.8), tolerance = 1e-12)
  })

  it("returns all NA when every anchor is dropped, or there are none", {
    out <- .il(c(0, 10), c(0, NA), grid_off = 0:3)
    expect_true(all(is.na(out$value)))
    expect_true(all(is.na(out$detected)))
    out <- interpolate_loglinear(as.Date(character()), numeric(), logical(),
                                 .il_d0 + 0:3)
    expect_equal(nrow(out), 4L)
    expect_true(all(is.na(out$value)))
    expect_true(all(is.na(out$detected)))
  })
})

describe("interpolate_loglinear(): input validation", {
  it("errors with class hydroSense_input on duplicated dates", {
    expect_error(
      interpolate_loglinear(.il_d0 + c(0, 5, 5), c(1, 2, 3), rep(TRUE, 3),
                            .il_d0 + 0:6),
      class = "hydroSense_input"
    )
  })
})
