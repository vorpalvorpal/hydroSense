## API: analyte_pafs() accessor + attribute (issue #30 — replaced the former
## per-row list-column with a flat per-(sample, draw, analyte) attribute).

library(testthat)
library(leachatetools)

co <- c(pH = 7.5, DOC = 2.0, Ca = 6.0, Mg = 4.0, hardness = 30.0)
mk <- function(sid, cu, zn, ni, draw_id = NA_integer_) {
  tibble::tibble(
    sample_id = sid, site_id = "f1", datetime = as.Date("2024-01-01"),
    analyte = c("Cu", "Zn", "Ni", names(co)),
    value = c(cu, zn, ni, unname(co)), detected = TRUE,
    draw_id = draw_id)
}

test_that("point mode: analyte_pafs() returns a flat breakdown, not a column", {
  df  <- dplyr::bind_rows(mk("s1", 5, 10, 0.3), mk("s2", 50, 20, 0.1))
  out <- suppressMessages(add_amspaf(df, reference = NULL, conc_units = "ug/L"))

  expect_false("analyte_pafs" %in% names(out))          # no list-column
  bd <- analyte_pafs(out)
  expect_s3_class(bd, "data.frame")
  expect_true(all(c("site_id", "sample_id", "analyte", "C_adj", "PAF",
                    "moa_group", "ref_source") %in% names(bd)))
  expect_false("draw_id" %in% names(bd))                # point mode: no draw_id
  # one row per assessed (sample, analyte): Cu/Zn/Ni for 2 samples
  expect_setequal(unique(bd$analyte), c("Cu", "Zn", "Ni"))
  expect_identical(nrow(bd), 6L)
  expect_true(all(bd$PAF >= 0 & bd$PAF <= 1))
})

test_that("draws mode: analyte_pafs() is keyed per (sample, draw, analyte)", {
  metals <- purrr::map_dfr(1:3, function(d)
    dplyr::bind_rows(
      tibble::tibble(sample_id = "s1", site_id = "f1",
                     datetime = as.Date("2024-01-01"),
                     analyte = c("Cu", "Zn", "Ni"),
                     value = c(5, 10, 0.3) * d, detected = TRUE, draw_id = d)))
  cofix <- mk("s1", 0, 0, 0)[4:8, ] |> dplyr::mutate(draw_id = NA_integer_)
  df <- dplyr::bind_rows(metals, cofix)
  out <- suppressMessages(add_amspaf(df, reference = NULL, conc_units = "ug/L",
                                     return = "draws"))
  bd <- analyte_pafs(out)
  expect_true("draw_id" %in% names(bd))
  expect_setequal(unique(bd$draw_id), 1:3)
  # 3 metals x 3 draws
  expect_identical(nrow(bd), 9L)
})

test_that("analyte_pafs() returns NULL (with message) when attribute absent", {
  expect_message(res <- analyte_pafs(tibble::tibble(x = 1)), "analyte_pafs")
  expect_null(res)
})
