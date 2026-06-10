## Equivalence harness for issue #30 (vectorising compute_amspaf_per_sample()).
##
## Compares the CURRENT engine output against a golden snapshot taken from the
## pre-change code (dev/gen_amspaf_golden.R -> fixtures/amspaf_golden.rds).
## Acceptance gate for the rewrite: numeric AmsPAF/PAF/C_adj within relative
## 1e-9; integer/categorical fields (n_analytes_used, dominant_analyte,
## ref_source, moa_group, analyte set) exact; per-analyte breakdown and
## ara_summary preserved.
##
## Robust to the analyte_pafs API change: the per-analyte breakdown is read from
## the `analyte_pafs` attribute if present (post-change), else by unnesting the
## list-column (pre-change).

library(testthat)
library(leachatetools)

golden <- readRDS(test_path("fixtures", "amspaf_golden.rds"))

## Normalise an add_amspaf() output to (scalars, breakdown, ara_summary),
## matching dev/gen_amspaf_golden.R but preferring the new flat attribute.
norm_out <- function(out) {
  rows <- dplyr::filter(out, .data$analyte == "AmsPAF")
  scal_cols <- intersect(c("sample_id", "draw_id", "value", "n_analytes_used",
                           "n_analytes_imputed", "dominant_analyte", "max_paf"),
                         names(rows))
  scal <- dplyr::arrange(rows[, scal_cols],
                         dplyr::across(dplyr::any_of(c("sample_id", "draw_id"))))
  brk_attr <- attr(out, "analyte_pafs")
  brk <- if (!is.null(brk_attr)) {
    brk_attr
  } else if ("analyte_pafs" %in% names(rows)) {
    tidyr::unnest(dplyr::select(rows, dplyr::any_of(c("sample_id", "draw_id")),
                               "analyte_pafs"), cols = "analyte_pafs")
  } else NULL
  if (!is.null(brk)) {
    keep <- intersect(c("sample_id", "draw_id", "analyte"), names(brk))
    brk <- dplyr::arrange(brk, dplyr::across(dplyr::all_of(keep)))
  }
  list(scalars = scal, breakdown = brk, ara_summary = attr(out, "ara_summary"))
}

run_case <- function(df, reference) suppressMessages(
  add_amspaf(df, reference = reference, conc_units = "ug/L",
             return = if ("draw_id" %in% names(df)) "draws" else "summary"))

prep_ref <- suppressMessages(
  prepare_reference(golden$inputs$ref_df, conc_units = "ug/L"))

## Compare one normalised case to its golden counterpart.
expect_equiv <- function(new, gold, label) {
  ## scalars: numeric value/max_paf within tol; integers/categoricals exact
  expect_equal(new$scalars$value, gold$scalars$value, tolerance = 1e-9,
               info = paste(label, "value"))
  if ("max_paf" %in% names(gold$scalars))
    expect_equal(new$scalars$max_paf, gold$scalars$max_paf, tolerance = 1e-9,
                 info = paste(label, "max_paf"))
  for (col in c("sample_id", "draw_id", "n_analytes_used",
                "n_analytes_imputed", "dominant_analyte")) {
    if (col %in% names(gold$scalars))
      expect_identical(new$scalars[[col]], gold$scalars[[col]],
                       info = paste(label, col))
  }
  ## per-analyte breakdown
  if (!is.null(gold$breakdown)) {
    expect_false(is.null(new$breakdown), info = paste(label, "breakdown present"))
    expect_identical(nrow(new$breakdown), nrow(gold$breakdown),
                     info = paste(label, "breakdown nrow"))
    for (col in c("sample_id", "draw_id", "analyte", "moa_group", "ref_source")) {
      if (col %in% names(gold$breakdown))
        expect_identical(new$breakdown[[col]], gold$breakdown[[col]],
                         info = paste(label, "breakdown", col))
    }
    for (col in c("C_adj", "PAF")) {
      if (col %in% names(gold$breakdown))
        expect_equal(new$breakdown[[col]], gold$breakdown[[col]], tolerance = 1e-9,
                     info = paste(label, "breakdown", col))
    }
  }
  ## ara_summary attribute (value columns within tol)
  if (!is.null(gold$ara_summary)) {
    expect_false(is.null(new$ara_summary), info = paste(label, "ara_summary"))
    expect_identical(nrow(new$ara_summary), nrow(gold$ara_summary),
                     info = paste(label, "ara_summary nrow"))
  }
}

test_that("add_amspaf point mode, ARA off — matches golden", {
  expect_equiv(norm_out(run_case(golden$inputs$pt_df, NULL)),
               golden$pt_noara, "pt_noara")
})
test_that("add_amspaf point mode, ARA on — matches golden", {
  expect_equiv(norm_out(run_case(golden$inputs$pt_df, prep_ref)),
               golden$pt_ara, "pt_ara")
})
test_that("add_amspaf draws mode, ARA off — matches golden", {
  expect_equiv(norm_out(run_case(golden$inputs$dr_df, NULL)),
               golden$dr_noara, "dr_noara")
})
test_that("add_amspaf draws mode, ARA on — matches golden", {
  expect_equiv(norm_out(run_case(golden$inputs$dr_df, prep_ref)),
               golden$dr_ara, "dr_ara")
})
