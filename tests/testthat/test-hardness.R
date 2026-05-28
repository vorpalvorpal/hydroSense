## Tests for derive_hardness()

library(testthat)
library(leachatetools)

base_chem <- function() {
  tibble::tribble(
    ~sample_id, ~site_id, ~datetime,            ~analyte,   ~value, ~detected,
    "A",        "X",       as.Date("2024-01-01"), "Ca",       20,     TRUE,
    "A",        "X",       as.Date("2024-01-01"), "Mg",        5,     TRUE,
    "B",        "X",       as.Date("2024-01-02"), "Ca",       10,     TRUE,
    "B",        "X",       as.Date("2024-01-02"), "hardness", 50,     TRUE,
    "C",        "X",       as.Date("2024-01-03"), "Ca",       20,     TRUE,
    "C",        "X",       as.Date("2024-01-03"), "Mg",        5,     TRUE,
    "C",        "X",       as.Date("2024-01-03"), "hardness", 70.535, TRUE,
    "D",        "X",       as.Date("2024-01-04"), "Ca",       20,     TRUE,
    "D",        "X",       as.Date("2024-01-04"), "Mg",        5,     TRUE,
    "D",        "X",       as.Date("2024-01-04"), "hardness",200,     TRUE,
    "E",        "X",       as.Date("2024-01-05"), "Ca",       30,     TRUE
  )
}

test_that("Ca + Mg → derive hardness exactly", {
  out  <- derive_hardness(base_chem(), verbose = FALSE)
  hA   <- out$value[out$sample_id == "A" & out$analyte == "hardness"]
  expect_length(hA, 1L)
  expect_equal(hA, 2.497 * 20 + 4.118 * 5)
})

test_that("Ca + hardness → derive Mg exactly", {
  out  <- derive_hardness(base_chem(), verbose = FALSE)
  mgB  <- out$value[out$sample_id == "B" & out$analyte == "Mg"]
  expect_length(mgB, 1L)
  expect_equal(mgB, (50 - 2.497 * 10) / 4.118)
})

test_that("Mg + hardness → derive Ca exactly", {
  df <- tibble::tribble(
    ~sample_id, ~site_id, ~datetime,            ~analyte,   ~value, ~detected,
    "F",        "X",       as.Date("2024-01-06"), "Mg",        4,     TRUE,
    "F",        "X",       as.Date("2024-01-06"), "hardness", 30,     TRUE
  )
  out <- derive_hardness(df, verbose = FALSE)
  ca  <- out$value[out$sample_id == "F" & out$analyte == "Ca"]
  expect_equal(ca, (30 - 4.118 * 4) / 2.497)
})

test_that("all three present + consistent → no new row, no warning", {
  df <- dplyr::filter(base_chem(), sample_id == "C")
  expect_silent(out <- derive_hardness(df, verbose = FALSE))
  # Only the original 3 rows should remain (no derived row appended)
  expect_equal(nrow(out), 3L)
})

test_that("all three present + inconsistent → warn but don't modify", {
  df <- dplyr::filter(base_chem(), sample_id == "D")
  expect_warning(
    out <- derive_hardness(df, verbose = FALSE),
    regexp = "disagreement"
  )
  # Original 3 rows preserved
  expect_equal(nrow(out), 3L)
})

test_that("only one of three present → no change", {
  df <- dplyr::filter(base_chem(), sample_id == "E")
  out <- derive_hardness(df, verbose = FALSE)
  expect_equal(nrow(out), nrow(df))
})

test_that("negative derived Mg dropped (hardness < 2.497*Ca)", {
  df <- tibble::tribble(
    ~sample_id, ~site_id, ~datetime,            ~analyte,   ~value, ~detected,
    "G",        "X",       as.Date("2024-01-07"), "Ca",      100,     TRUE,
    "G",        "X",       as.Date("2024-01-07"), "hardness", 10,     TRUE   # impossible
  )
  expect_warning(
    out <- derive_hardness(df, verbose = FALSE),
    regexp = "negative derived Mg"
  )
  # No Mg row was appended
  expect_false("Mg" %in% out$analyte)
})

test_that("hardness ≈ 0 handled without divide-by-zero", {
  df <- tibble::tribble(
    ~sample_id, ~site_id, ~datetime,            ~analyte,   ~value, ~detected,
    "H",        "X",       as.Date("2024-01-08"), "Ca",       0.01,   TRUE,
    "H",        "X",       as.Date("2024-01-08"), "Mg",       0.01,   TRUE,
    "H",        "X",       as.Date("2024-01-08"), "hardness", 0,      TRUE
  )
  # hardness = 0 but calc says ~0.066 → relative err > tolerance
  expect_warning(
    out <- derive_hardness(df, tolerance = 0.05, verbose = FALSE),
    regexp = "disagreement"
  )
})

test_that("idempotent when no derivations needed", {
  df  <- dplyr::filter(base_chem(), sample_id == "C")
  out1 <- derive_hardness(df, verbose = FALSE)
  out2 <- derive_hardness(out1, verbose = FALSE)
  expect_equal(nrow(out1), nrow(out2))
})

test_that("derived rows tagged with imputed = TRUE, imputed_kind = derived", {
  df  <- dplyr::filter(base_chem(), sample_id == "A")
  out <- derive_hardness(df, verbose = FALSE)
  derived_row <- out[out$analyte == "hardness", ]
  expect_true(derived_row$imputed)
  expect_equal(derived_row$imputed_kind, "derived")
})

test_that("output preserves site_id and datetime metadata on derived rows", {
  df  <- dplyr::filter(base_chem(), sample_id == "A")
  out <- derive_hardness(df, verbose = FALSE)
  derived_row <- out[out$analyte == "hardness", ]
  expect_equal(derived_row$site_id, "X")
  expect_equal(derived_row$datetime, as.Date("2024-01-01"))
})
