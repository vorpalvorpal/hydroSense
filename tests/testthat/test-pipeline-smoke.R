## End-to-end smoke tests for the chronic AmsPAF pipeline
##
## Structure:
##   § 1  Always-run: prescreen → chronic → prepare_reference (no external deps)
##   § 2  Gated on guideline_dir: full pipeline including add_amspaf
##   § 3  Gated on BRMS_SMOKE_TEST=1: full pipeline including impute_chemistry
##
## The pipeline under test:
##   prescreen_analytes() → [impute_chemistry()] →
##   compute_chronic_chemistry() → prepare_reference() → add_amspaf()

library(testthat)
library(leachatetools)

# ── Shared synthetic dataset ──────────────────────────────────────────────────
#
# Two sites: "ref" (reference) and "ds" (downstream).
# Three metal analytes (Cu, Zn, Ni) and three drivers (pH, EC, DOC).
# 40 samples per site spread over 3 years so focal dates in year 3 have
# a full 365-day window to draw from.

make_pipeline_chem <- function(seed = 42) {
  set.seed(seed)
  n_per_site <- 40
  sites    <- c("ref", "ds")
  # Analytes: three metals, three imputation drivers, two Ni/Zn co-analytes
  # Ca and Mg are required by the Ni normalisation formula.
  # hardness (as Ca hardness proxy) is required by Zn.
  analytes <- c("Cu", "Zn", "Ni",
                "pH", "EC", "DOC",
                "Ca", "Mg", "hardness")

  dates <- sort(sample(
    seq(as.Date("2023-01-01"), as.Date("2025-12-31"), by = "day"),
    n_per_site
  ))

  purrr::map_dfr(sites, function(site) {
    tidyr::expand_grid(
      sample_id = paste0(site, "_s", seq_len(n_per_site)),
      analyte   = analytes
    ) |>
      dplyr::mutate(
        site_id  = site,
        datetime = dates[match(
          sub(paste0(site, "_s"), "", sample_id),
          as.character(seq_len(n_per_site))
        )],
        value = dplyr::case_when(
          analyte == "pH"       ~ runif(dplyr::n(), 6.8, 7.8),
          analyte == "EC"       ~ runif(dplyr::n(), 150, 400),
          analyte == "DOC"      ~ runif(dplyr::n(), 0.5, 5),
          analyte == "Ca"       ~ runif(dplyr::n(), 4, 12),     # mg/L (Ni valid: 2–70)
          analyte == "Mg"       ~ runif(dplyr::n(), 2, 8),      # mg/L (Ni valid: 1.6–54)
          analyte == "hardness" ~ runif(dplyr::n(), 25, 80),    # mg/L CaCO3 (Zn valid: 20–440)
          analyte == "Cu"       ~ exp(rnorm(dplyr::n(), log(2), 0.6)),
          analyte == "Zn"       ~ exp(rnorm(dplyr::n(), log(5), 0.7)),
          analyte == "Ni"       ~ exp(rnorm(dplyr::n(), log(3), 0.5))
        ),
        # A few BDL rows in the downstream site
        detected = !(site == "ds" & analyte == "Cu" &
                       sample_id %in% paste0("ds_s", 1:3)),
        imputed  = FALSE
      )
  })
}


# ─────────────────────────────────────────────────────────────────────────────
# § 1  Core pipeline: prescreen → chronic → prepare_reference
# No external dependencies; always runs.
# ─────────────────────────────────────────────────────────────────────────────

test_that("§1a prescreen_analytes returns expected analytes and protects drivers", {
  chem     <- make_pipeline_chem()
  included <- prescreen_analytes(chem, k = 0.05)
  expect_type(included, "character")
  # All analytes should pass (high detection frequency in synthetic data)
  expect_true("Cu"  %in% included)
  expect_true("Zn"  %in% included)
  expect_true("Ni"  %in% included)
  # Drivers are always present and detected → in included
  expect_true("pH"  %in% included)
  expect_true("DOC" %in% included)
})

test_that("§1b compute_chronic_chemistry produces one row per (date × site × analyte)", {
  chem       <- make_pipeline_chem()
  focal      <- as.Date(c("2025-04-01", "2025-10-01"))
  chr_chem   <- compute_chronic_chemistry(chem, focal_dates = focal,
                                          tau_days = 90, window_days = 365)

  # Each focal date × each site × each analyte should have exactly one row
  expected_rows <- length(focal) * length(unique(chem$site_id)) *
    length(unique(chem$analyte))
  expect_equal(nrow(chr_chem), expected_rows)

  # Required columns present
  expected_cols <- c("focal_date", "site_id", "sample_id",
                     "analyte", "value", "detected",
                     "n_samples_in_window", "n_imputed_in_window")
  expect_true(all(expected_cols %in% names(chr_chem)))

  # All values finite and positive
  expect_true(all(is.finite(chr_chem$value)))
  expect_true(all(chr_chem$value > 0))
})

test_that("§1c prepare_reference produces a valid prepared_reference object", {
  chem    <- make_pipeline_chem()
  ref_raw <- dplyr::filter(chem, site_id == "ref")
  focal   <- as.Date("2025-04-01")
  chr_ref <- compute_chronic_chemistry(ref_raw, focal_dates = focal,
                                       tau_days = 90, window_days = 365)
  prep    <- prepare_reference(chr_ref)

  expect_s3_class(prep, "prepared_reference")
  expect_true(is.data.frame(prep$normalised_quantiles))
  expect_true(all(c("analyte", "ref_norm") %in% names(prep$normalised_quantiles)))
  expect_true(all(is.finite(prep$normalised_quantiles$ref_norm)))
  expect_equal(prep$percentile, 0.80)
})

test_that("§1d expand_focal_dates integrates with compute_chronic_chemistry", {
  chem       <- make_pipeline_chem()
  focal      <- expand_focal_dates("2025-01-01", "2025-03-31", by = "week")
  expect_gte(length(focal), 13L)  # roughly 13 weekly intervals

  chr_chem <- compute_chronic_chemistry(chem, focal_dates = focal,
                                        tau_days = 90, window_days = 365)
  # One row per (date × site × analyte)
  n_sites    <- length(unique(chem$site_id))
  n_analytes <- length(unique(chem$analyte))
  expect_equal(nrow(chr_chem), length(focal) * n_sites * n_analytes)
})

test_that("§1e chronic reference baseline is lower for reference than downstream", {
  # Reference site should have lower chronic concentrations for metals on average
  chem <- make_pipeline_chem()

  # Add an elevated downstream signal
  chem <- dplyr::mutate(chem,
    value = dplyr::if_else(
      site_id == "ds" & analyte == "Cu",
      value * 5,  # 5× elevated Cu at downstream
      value
    )
  )

  focal    <- as.Date("2025-06-01")
  chr_chem <- compute_chronic_chemistry(chem, focal_dates = focal,
                                        tau_days = 90, window_days = 365)

  ref_cu <- chr_chem$value[chr_chem$site_id == "ref" & chr_chem$analyte == "Cu"]
  ds_cu  <- chr_chem$value[chr_chem$site_id == "ds"  & chr_chem$analyte == "Cu"]
  expect_gt(ds_cu, ref_cu)
})


# ─────────────────────────────────────────────────────────────────────────────
# § 2  Full pipeline including add_amspaf
# Requires leachatetools.guideline_dir option to be set.
# ─────────────────────────────────────────────────────────────────────────────

test_that("§2 full chronic AmsPAF pipeline runs end-to-end", {
  skip_if(
    is.null(getOption("leachatetools.guideline_dir")),
    "Skipping: leachatetools.guideline_dir not set (no ANZG XLSX available)"
  )

  chem <- make_pipeline_chem()

  # Step 1: prescreen
  included <- prescreen_analytes(chem, k = 0.05)
  chem_f   <- dplyr::filter(chem, analyte %in% included)

  # Step 2: chronic integration for downstream and reference
  focal <- as.Date(c("2025-04-01", "2025-10-01"))

  ds_chem  <- dplyr::filter(chem_f, site_id == "ds")
  ref_chem <- dplyr::filter(chem_f, site_id == "ref")

  chr_ds  <- compute_chronic_chemistry(ds_chem, focal_dates = focal)
  chr_ref <- compute_chronic_chemistry(ref_chem, focal_dates = focal)

  # Step 3: prepare reference
  prep_ref <- prepare_reference(chr_ref)

  # Step 4: AmsPAF
  out <- add_amspaf(chr_ds, reference = prep_ref)

  # Should have AmsPAF rows appended
  amspaf_rows <- dplyr::filter(out, analyte == "AmsPAF")
  expect_gte(nrow(amspaf_rows), 1L)

  # AmsPAF values should be non-negative finite numbers
  expect_true(all(is.finite(amspaf_rows$value)))
  expect_true(all(amspaf_rows$value >= 0))

  # n_analytes_used should be populated
  expect_true("n_analytes_used" %in% names(amspaf_rows))
  expect_true(all(amspaf_rows$n_analytes_used >= 1L))

  # focal_date column should be preserved in output
  expect_true("focal_date" %in% names(amspaf_rows))
  expect_equal(sort(unique(amspaf_rows$focal_date)), sort(focal))
})


# ─────────────────────────────────────────────────────────────────────────────
# § 3  Full pipeline including impute_chemistry (brms smoke test)
# Gated on BRMS_SMOKE_TEST=1 and brms/Stan installation.
# ─────────────────────────────────────────────────────────────────────────────

test_that("§3 full pipeline with imputation runs end-to-end", {
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping: brms not installed or BRMS_SMOKE_TEST != '1'"
  )
  skip_if(
    is.null(getOption("leachatetools.guideline_dir")),
    "Skipping: leachatetools.guideline_dir not set"
  )

  chem <- make_pipeline_chem()

  # Step 1: prescreen (all analytes pass in synthetic data)
  included <- prescreen_analytes(chem, k = 0.05)
  chem_f   <- dplyr::filter(chem, analyte %in% included)

  # Step 2: impute (fast settings for CI)
  imp <- impute_chemistry(
    chem_f,
    drivers = c("pH", "EC", "DOC"),
    iter    = 500,
    warmup  = 250,
    chains  = 1,
    cores   = 1
  )

  expect_true(all(is.finite(imp$value)))
  expect_true("imputed" %in% names(imp))

  # Step 3: chronic
  focal <- as.Date(c("2025-04-01", "2025-10-01"))

  chr_ds  <- compute_chronic_chemistry(
    dplyr::filter(imp, site_id == "ds"), focal_dates = focal
  )
  chr_ref <- compute_chronic_chemistry(
    dplyr::filter(imp, site_id == "ref"), focal_dates = focal
  )

  # Step 4: prepare reference + AmsPAF
  prep_ref <- prepare_reference(chr_ref)
  out      <- add_amspaf(chr_ds, reference = prep_ref)

  amspaf_rows <- dplyr::filter(out, analyte == "AmsPAF")
  expect_gte(nrow(amspaf_rows), 1L)
  expect_true(all(is.finite(amspaf_rows$value)))

  # n_analytes_imputed should be propagated from imp through chronic to AmsPAF
  expect_true("n_analytes_imputed" %in% names(amspaf_rows))
})
