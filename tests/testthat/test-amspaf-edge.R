## Edge cases for the vectorised AmsPAF engine (issue #30). The equivalence
## golden already covers BDL -> 0, missing-co-analyte drop, and the
## <min_analytes drop; these add unmatched-reference and all-BDL behaviour.

library(testthat)
library(hydroSense)

co <- c(pH = 7.5, DOC = 2.0, Ca = 6.0, Mg = 4.0, hardness = 30.0)
mk <- function(sid, cu, zn, ni, det = TRUE) {
  tibble::tibble(
    sample_id = sid, site_id = "f1", datetime = as.Date("2024-01-01"),
    analyte = c("Cu", "Zn", "Ni", names(co)),
    value = c(cu, zn, ni, unname(co)),
    detected = c(rep(det, 3L), rep(TRUE, length(co))))
}

test_that("unmatched reference: analyte assessed against raw conc (ref_norm = 0)", {
  df  <- dplyr::bind_rows(mk("s1", 5, 10, 0.3), mk("s2", 50, 20, 0.1))
  ## Reference only carries Zn + co-analytes -> Cu and Ni are 'unmatched'.
  ref <- dplyr::bind_rows(mk("r1", NA, 3, NA), mk("r2", NA, 4, NA)) |>
    dplyr::filter(!is.na(.data$value))
  pref <- suppressMessages(prepare_reference(ref, conc_units = "ug/L"))
  out  <- suppressMessages(add_amspaf(df, reference = pref, conc_units = "ug/L"))

  bd <- analyte_pafs(out)
  expect_identical(sort(unique(bd$ref_source[bd$analyte %in% c("Cu", "Ni")])),
                   "unmatched")
  expect_identical(unique(bd$ref_source[bd$analyte == "Zn"]), "matched")
  ## unmatched analytes: C_adj == C_norm (nothing subtracted) -> positive
  expect_true(all(bd$C_adj[bd$analyte == "Cu"] > 0))
})

test_that("all-BDL metals: C_adj = 0 -> AmsPAF = 0 (co-analytes still present)", {
  df  <- mk("s1", 5, 10, 0.3, det = FALSE)   # Cu/Zn/Ni below detection
  out <- suppressMessages(add_amspaf(df, reference = NULL, conc_units = "ug/L"))
  amspaf <- dplyr::filter(out, .data$analyte == "AmsPAF")
  expect_identical(nrow(amspaf), 1L)
  expect_equal(amspaf$value, 0)
  bd <- analyte_pafs(out)
  expect_true(all(bd$C_adj == 0))
  expect_true(all(bd$PAF == 0))
})

test_that("ARA off: ref_source is 'disabled' for every analyte", {
  df  <- mk("s1", 5, 10, 0.3)
  out <- suppressMessages(add_amspaf(df, reference = NULL, conc_units = "ug/L"))
  bd  <- analyte_pafs(out)
  expect_identical(unique(bd$ref_source), "disabled")
})
