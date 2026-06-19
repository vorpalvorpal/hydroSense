## Stan-free coverage for impute.R: the deterministic helpers and the
## mgcv-based co-analyte imputation. These need no brms/Stan, so they run in
## the default suite and cover the parts of impute.R that the (gated) brms
## smoke test does not.

library(testthat)
library(hydroSense)

# Chemistry with the standard PCA panel; DOC missing for some samples so
# impute_coanalytes() has something to fill.
make_coan_chem <- function(n = 30, drop_doc_for = 6L, seed = 1) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  mk <- function(an, vals) {
    tibble::tibble(sample_id = ids, site_id = "A",
                   datetime = as.Date("2023-01-01") + seq_len(n),
                   analyte = an, value = vals, detected = TRUE)
  }
  df <- dplyr::bind_rows(
    mk("pH",          runif(n, 6, 8)),
    mk("Cl",          rlnorm(n, 5, 1)),
    mk("SO4²⁻", rlnorm(n, 4, 1.5)),
    mk("Ca",          rlnorm(n, 3, 1)),
    mk("temperature", runif(n, 5, 20)),
    mk("DOC",         rlnorm(n, 1, 1))
  )
  # Remove whole DOC rows for the last `drop_doc_for` samples.
  drop_ids <- tail(ids, drop_doc_for)
  df[!(df$analyte == "DOC" & df$sample_id %in% drop_ids), ]
}

pca_vars_used <- c("pH", "Cl", "SO4²⁻", "Ca", "temperature", "DOC")

# Build an imputation_model carrying only a fitted PCA (no brms models needed
# for impute_coanalytes()).
make_pca_model <- function(df) {
  pca <- hydroSense:::.prepare_chem_pca(
    df, wq_vars = pca_vars_used, min_var_explained = 0.75, max_pcs = 4L
  )
  structure(
    list(pca = pca, pca_vars = pca_vars_used, fit_date = "2023-01-01",
         n_samples = dplyr::n_distinct(df$sample_id), groups = list()),
    class = "imputation_model"
  )
}

test_that("impute_coanalytes fills missing co-analyte rows via the GAM path", {
  df    <- make_coan_chem(n = 30, drop_doc_for = 6L)
  model <- make_pca_model(df)

  out <- suppressMessages(impute_coanalytes(df, model, targets = "DOC"))

  # Every sample now has a DOC row.
  doc <- out[out$analyte == "DOC", ]
  expect_equal(dplyr::n_distinct(doc$sample_id), 30L)
  # The six previously-missing samples are flagged as imputed.
  imp <- doc[doc$imputed, ]
  expect_equal(nrow(imp), 6L)
  expect_true(all(imp$imputed_kind == "missing"))
  expect_true(all(is.finite(imp$value) & imp$value > 0))
  # Observed rows are not flagged.
  expect_false(any(doc$imputed[!(doc$sample_id %in% imp$sample_id)]))
})

test_that("impute_coanalytes skips a target with too few observations", {
  df    <- make_coan_chem(n = 30, drop_doc_for = 25L)  # only 5 quantified DOC
  model <- make_pca_model(df)
  expect_warning(
    out <- impute_coanalytes(df, model, targets = "DOC", min_obs = 10L),
    "quantified sample"
  )
  # No imputed DOC rows added (target skipped).
  doc <- out[out$analyte == "DOC", ]
  expect_false(any(doc$imputed %in% TRUE))
})

test_that("impute_coanalytes errors without a fitted PCA", {
  df    <- make_coan_chem()
  model <- structure(list(pca = NULL, pca_vars = pca_vars_used),
                     class = "imputation_model")
  expect_error(impute_coanalytes(df, model), "no fitted PCA")
})

test_that(".check_bdl_imputed caps imputed BDL values above the detection limit", {
  result <- tibble::tibble(
    sample_id    = c("s1", "s2", "s3"),
    analyte      = "Cu",
    value        = c(5, 0.5, 2),
    imputed_kind = c("censored_left", "censored_left", "observed")
  )
  dl <- tibble::tibble(sample_id = c("s1", "s2"), analyte = "Cu",
                       detection_limit = c(1, 1))

  # s1 is BDL with value 5 > DL 1 -> exceedance -> warn + cap to 1.
  expect_warning(
    capped <- hydroSense:::.check_bdl_imputed(result, dl, cap = TRUE), "exceed"
  )
  expect_equal(capped$value[capped$sample_id == "s1"], 1)  # capped
  expect_equal(capped$value[capped$sample_id == "s2"], 0.5) # untouched (below DL)
  expect_equal(capped$value[capped$sample_id == "s3"], 2)   # observed, untouched
})

test_that(".check_bdl_imputed attaches an auditable per-cell cap summary", {
  result <- tibble::tibble(
    sample_id    = c("s1", "s2", "s3"),
    analyte      = c("Cu", "Zn", "Cu"),
    value        = c(5, 4, 0.5),
    imputed_kind = "censored_left"
  )
  dl <- tibble::tibble(sample_id = c("s1", "s2", "s3"), analyte = c("Cu", "Zn", "Cu"),
                       detection_limit = c(1, 1, 1))
  suppressWarnings(capped <- hydroSense:::.check_bdl_imputed(result, dl, cap = TRUE))

  s <- bdl_cap_summary(capped)
  expect_s3_class(s, "tbl_df")
  expect_setequal(s$analyte, c("Cu", "Zn"))           # s1/Cu and s2/Zn exceeded; s3 did not
  expect_true(all(s$capped))
  expect_equal(s$max_ratio[s$analyte == "Cu"], 5)     # 5 / DL 1
  # accessor on a frame with no activations returns NULL invisibly + a message
  expect_message(out <- bdl_cap_summary(result), "No detection-limit cap")
  expect_null(out)
})

test_that(".check_bdl_imputed caps per-row without duplicating posterior draws", {
  # return = "draws" gives several rows per cell; the cap must clip each row
  # against the single DL without many-to-many row duplication.
  result <- tibble::tibble(
    sample_id    = "s1",
    analyte      = "Cu",
    .draw        = 1:3,
    value        = c(5, 0.5, 8),        # draws 1 and 3 exceed DL = 1
    imputed_kind = "censored_left"
  )
  dl <- tibble::tibble(sample_id = "s1", analyte = "Cu", detection_limit = 1)
  suppressWarnings(out <- hydroSense:::.check_bdl_imputed(result, dl, cap = TRUE))

  expect_equal(nrow(out), 3L)                # no duplication
  expect_equal(out$value, c(1, 0.5, 1))      # only exceeding draws capped
  expect_equal(bdl_cap_summary(out)$n_rows, 2L)
})

test_that(".check_bdl_imputed warns but does not cap when cap = FALSE", {
  result <- tibble::tibble(sample_id = "s1", analyte = "Cu", value = 5,
                           imputed_kind = "censored_left")
  dl <- tibble::tibble(sample_id = "s1", analyte = "Cu", detection_limit = 1)
  expect_warning(
    out <- hydroSense:::.check_bdl_imputed(result, dl, cap = FALSE), "NOT capped"
  )
  expect_equal(out$value, 5)
})

test_that(".check_bdl_imputed is a no-op with an empty DL table or no exceedance", {
  result <- tibble::tibble(sample_id = "s1", analyte = "Cu", value = 0.4,
                           imputed_kind = "censored_left")
  empty_dl <- tibble::tibble(sample_id = character(), analyte = character(),
                             detection_limit = numeric())
  expect_identical(hydroSense:::.check_bdl_imputed(result, empty_dl), result)

  # Present DL but value below it -> no exceedance, unchanged, no warning.
  dl <- tibble::tibble(sample_id = "s1", analyte = "Cu", detection_limit = 1)
  expect_identical(hydroSense:::.check_bdl_imputed(result, dl), result)
})

test_that("print.imputation_model summarises the fit", {
  m <- structure(
    list(fit_date = "2023-01-01", n_samples = 30L,
         pca_vars = c("pH", "Cl", "DOC"),
         pca = list(n_pcs = 2L, var_explained = 0.8),
         impute_method = "rescor_mi",
         groups = list(metals = list(analytes = c("Cu", "Zn")))),
    class = "imputation_model"
  )
  expect_output(print(m), "imputation_model")
  expect_output(print(m), "metals:")
})
