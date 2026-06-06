## Regression guard for the unspeciated-arsenic SSD default.
##
## `.SSD_NAME_MAP[["As"]]` resolves unspeciated As to the As(V) SSD on the
## grounds that this is the conservative *max-PAF envelope* of the two
## speciation SSDs — i.e. the As(V) SSD yields a higher (more protective) PAF
## than the As(III) SSD at every environmentally realistic concentration, even
## though As(III) is mechanistically more toxic to the single most-sensitive
## diatom (the whole-assemblage As(V) SSD sits lower because its lower tail is
## pulled down by the corrected algal acute points).
##
## This is an EMPIRICAL property of the current Warne 2000 datasets, not a
## theorem. If either SSD's underlying data is ever revised so that As(III)
## overtakes As(V) within the realistic range, the `As -> As_V` default would
## silently stop being conservative. These tests fail loudly if that happens,
## prompting a switch to a true runtime envelope (or an As(III) default).

library(testthat)
library(leachatetools)

gd <- function() {
  d <- getOption("leachatetools.guideline_dir")
  if (is.null(d)) "guideline data" else d
}

test_that("unspeciated 'As' resolves to the As(V) SSD", {
  expect_identical(leachatetools:::.SSD_NAME_MAP[["As"]], "As_V")
})

test_that("As(V) PAF >= As(III) PAF across the realistic concentration range", {
  skip_on_cran()
  # 0.1 µg/L (well below any DGV) up to 5000 µg/L (5 mg/L dissolved As — already
  # an extreme leachate value). The SSDs do not cross until ~10 mg/L.
  grid <- c(0.1, 0.5, 1, 3, 5, 7, 10, 13, 24, 50, 100, 200, 360, 1000, 2000, 5000)

  paf_v <- vapply(grid, function(c)
    ssd_paf("As(V)",   c, conc_units = "ug/L", guideline_dir = gd())$pct, numeric(1))
  paf_i <- vapply(grid, function(c)
    ssd_paf("As(III)", c, conc_units = "ug/L", guideline_dir = gd())$pct, numeric(1))

  ok <- !is.na(paf_v) & !is.na(paf_i)
  skip_if(sum(ok) < length(grid) / 2, "SSD fits unavailable in this environment")

  # The whole point: As(V) is the conservative envelope at every realistic conc.
  expect_true(
    all(paf_v[ok] >= paf_i[ok] - 1e-9),
    info = paste0(
      "As(III) overtook As(V) at conc(s): ",
      paste(grid[ok][paf_i[ok] > paf_v[ok] + 1e-9], collapse = ", "),
      " µg/L — the unspeciated 'As -> As_V' default is no longer the ",
      "conservative max-PAF envelope and must be revisited."
    )
  )
})

test_that("the unspeciated-As default equals the max-PAF envelope of both SSDs", {
  skip_on_cran()
  # ssd_paf('As') uses the As_V stem; confirm that equals max(As_V, As_III) at a
  # representative leachate concentration, i.e. the As_V default IS the envelope.
  conc <- 50
  paf_as   <- ssd_paf("As",      conc, conc_units = "ug/L", guideline_dir = gd())$pct
  paf_v    <- ssd_paf("As(V)",   conc, conc_units = "ug/L", guideline_dir = gd())$pct
  paf_i    <- ssd_paf("As(III)", conc, conc_units = "ug/L", guideline_dir = gd())$pct
  skip_if(anyNA(c(paf_as, paf_v, paf_i)), "SSD fits unavailable")

  expect_equal(paf_as, max(paf_v, paf_i), tolerance = 1e-8)
})
