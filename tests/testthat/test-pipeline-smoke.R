## End-to-end smoke tests for the chronic AmsPAF pipeline
##
## Structure:
##   § 1  Always-run: prescreen → chronic → prepare_reference (no external deps)
##   § 2  Gated on guideline_dir: full pipeline including add_amspaf
##   § 3  Gated on BRMS_SMOKE_TEST=1: full pipeline including impute_chemistry
##
## The pipeline under test:
##   prescreen_analytes() → [impute_chemistry()] →
##   time_weighted_aggregate() → prepare_reference() → add_amspaf()

library(testthat)
library(hydroSense)

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
  included <- prescreen_analytes(chem, k = 0.05, conc_units = "ug/L")
  expect_type(included, "character")
  # All analytes should pass (high detection frequency in synthetic data)
  expect_true("Cu"  %in% included)
  expect_true("Zn"  %in% included)
  expect_true("Ni"  %in% included)
  # Drivers are always present and detected → in included
  expect_true("pH"  %in% included)
  expect_true("DOC" %in% included)
})

test_that("§1b time_weighted_aggregate produces one row per (date × site × analyte)", {
  chem       <- make_pipeline_chem()
  focal      <- as.Date(c("2025-04-01", "2025-10-01"))
  chr_chem   <- time_weighted_aggregate(chem, focal_dates = focal,
                                          tau = 90, tau_units = "d",
                                          window = 365, window_units = "d")

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
  chr_ref <- time_weighted_aggregate(ref_raw, focal_dates = focal,
                                       tau = 90, tau_units = "d",
                                       window = 365, window_units = "d")
  prep    <- prepare_reference(chr_ref, conc_units = "ug/L")

  expect_s3_class(prep, "prepared_reference")
  expect_true(is.data.frame(prep$ref_table))
  expect_true(all(c("analyte", "ref_norm", "n_obs") %in% names(prep$ref_table)))
  expect_true(all(is.finite(prep$ref_table$ref_norm)))
  expect_equal(prep$summary, "geom_mean")
})

test_that("§1d expand_focal_dates integrates with time_weighted_aggregate", {
  chem       <- make_pipeline_chem()
  focal      <- expand_focal_dates("2025-01-01", "2025-03-31", by = "week")
  expect_gte(length(focal), 13L)  # roughly 13 weekly intervals

  chr_chem <- time_weighted_aggregate(chem, focal_dates = focal,
                                        tau = 90, tau_units = "d",
                                        window = 365, window_units = "d")
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
  chr_chem <- time_weighted_aggregate(chem, focal_dates = focal,
                                        tau = 90, tau_units = "d",
                                        window = 365, window_units = "d")

  ref_cu <- chr_chem$value[chr_chem$site_id == "ref" & chr_chem$analyte == "Cu"]
  ds_cu  <- chr_chem$value[chr_chem$site_id == "ds"  & chr_chem$analyte == "Cu"]
  expect_gt(ds_cu, ref_cu)
})


# ─────────────────────────────────────────────────────────────────────────────
# § 2  Full pipeline including add_amspaf
# Requires hydroSense.guideline_dir option to be set.
# ─────────────────────────────────────────────────────────────────────────────

test_that("§2 full chronic AmsPAF pipeline runs end-to-end", {

  chem <- make_pipeline_chem()

  # Step 1: prescreen
  included <- prescreen_analytes(chem, k = 0.05, conc_units = "ug/L")
  chem_f   <- dplyr::filter(chem, analyte %in% included)

  # Step 2: chronic integration for downstream and reference
  focal <- as.Date(c("2025-04-01", "2025-10-01"))

  ds_chem  <- dplyr::filter(chem_f, site_id == "ds")
  ref_chem <- dplyr::filter(chem_f, site_id == "ref")

  chr_ds  <- time_weighted_aggregate(ds_chem, focal_dates = focal)
  chr_ref <- time_weighted_aggregate(ref_chem, focal_dates = focal)

  # Step 3: prepare reference
  prep_ref <- prepare_reference(chr_ref, conc_units = "ug/L")

  # Step 4: AmsPAF
  out <- add_amspaf(chr_ds, reference = prep_ref, conc_units = "ug/L")

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

test_that("§3 full pipeline with imputation runs end-to-end (Path B)", {
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping: brms not installed or BRMS_SMOKE_TEST != '1'"
  )

  chem <- make_pipeline_chem()

  # Step 1: prescreen
  included <- prescreen_analytes(chem, k = 0.05,
                                  protect = c("pH", "EC"), conc_units = "ug/L")
  chem_f   <- dplyr::filter(chem, analyte %in% included)

  # Step 2: fit + impute (fast settings)
  model <- fit_imputation_model(
    chem_f,
    required_vars = c("pH", "EC"),
    iter    = 500,
    warmup  = 250,
    chains  = 1,
    cores   = 1
  )
  imp <- impute_chemistry(chem_f, model)

  expect_true(all(is.finite(imp$value)))
  expect_true("imputed" %in% names(imp))

  # Step 3: per-sample AmsPAF
  prep_ref <- prepare_reference(dplyr::filter(imp, site_id == "ref"),
                                conc_units = "ug/L")
  ps       <- add_amspaf(dplyr::filter(imp, site_id == "ds"),
                         reference = prep_ref, conc_units = "ug/L")
  ps_amspaf <- dplyr::filter(ps, analyte == "AmsPAF")
  expect_gte(nrow(ps_amspaf), 1L)

  # Step 4: chronic Path B — time-weighted arithmetic mean of per-sample AmsPAF
  focal     <- as.Date(c("2025-04-01", "2025-10-01"))
  chr_amspaf <- time_weighted_aggregate(
    ps_amspaf, focal_dates = focal, summary = "arith_mean"
  )

  expect_gte(nrow(chr_amspaf), 1L)
  expect_true(all(is.finite(chr_amspaf$value)))
})
