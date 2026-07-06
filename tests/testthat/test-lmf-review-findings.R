## ============================================================================
## TDD specifications derived from the LMF code review
## ============================================================================
##
## Each test below encodes the *target* (post-fix) behaviour for a confirmed
## review finding. They are written test-first: with the current
## implementation they are expected to FAIL, so each is guarded by
## `.skip_tdd()` to keep the suite green until the corresponding fix lands.
##
## To drive a fix red-green: delete the `.skip_tdd(...)` line at the top of the
## relevant test, watch it fail, implement the fix in R/lmf.R, watch it pass.
##
## Findings without a crisp behavioural contract are documented but not
## encoded here:
##   - Issue 7 (LMF is "fraction of median-strength leachate", not volumetric):
##     roxygen/doc wording only.
##   - Smaller note (leachate total-N filter runs on raw pre-BDL values in
##     implicit mg/L): internal consistency, no observable contract.
##   - Issue 2 full fix (GLS with the reference covariance matrix): the
##     directional overdispersion contract is captured below; the full
##     covariance solve is an implementation choice.
##   - Bug B3 (rename max_sigma_lsi -> max_sigma_lmf): API rename with a
##     deprecation shim; better verified by a dedicated deprecation test once
##     the shim exists.

library(testthat)
library(hydroSense)


## ----------------------------------------------------------------------------
## Issue 1 — the variance model omits uncertainty in the leachate end-member L.
## var_f = (sigma2_meas + sigma2_R) / (L - R)^2 treats L as exact. The missing
## first-order term is f^2 * sigma2_L / (L - R)^2, largest at high LMF.
## ----------------------------------------------------------------------------

test_that("[issue 1] sigma_lmf responds to leachate end-member spread", {
  ref <- .rev_endmembers()$reference

  ## Two leachate end-members with (approximately) the SAME mean composition
  ## but very different sample-to-sample spread.
  set.seed(2); lea_tight  <- .rev_samples(rep(1, 12), "lea", "leasite", noise_sd = 0.01)
  set.seed(2); lea_spread <- .rev_samples(rep(1, 12), "lea", "leasite", noise_sd = 0.25)

  ## High mixing fraction so the f^2 * sigma2_L term is material.
  q <- .rev_samples(0.8, "q", "query", noise_sd = 0)

  r_tight  <- .rev_lmf_row(add_lmf(q, leachate_data = lea_tight,
                                   reference_data = ref, max_sigma_lsi = 1e6))
  r_spread <- .rev_lmf_row(add_lmf(q, leachate_data = lea_spread,
                                   reference_data = ref, max_sigma_lsi = 1e6))

  ## A noisier leachate end-member must widen the mixing-fraction error bar.
  expect_gt(r_spread$sigma_lmf, r_tight$sigma_lmf * 1.15)
})


## ----------------------------------------------------------------------------
## Issue 2 — cross-ion error correlation / overdispersion is ignored. When ions
## genuinely disagree (chi2/df > 1) the inverse-variance sigma_lmf is optimistic.
## Minimal fix: inflate sigma_lmf by sqrt(chi2/df) (Birge ratio) when > 1; full
## fix: GLS with the reference covariance. Either satisfies this directional
## contract.
## ----------------------------------------------------------------------------

test_that("[issue 2] sigma_lmf is inflated when ions disagree (chi2/df > 1)", {
  em <- .rev_endmembers()

  q_clean <- .rev_samples(0.5, "q", "query", noise_sd = 0)

  ## Push ions off the mixing line in opposing directions: the mean fraction
  ## stays ~0.5 (weights barely move) but chi2/df climbs well above 1.
  q_conflict <- q_clean
  q_conflict <- .rev_scale_ion(q_conflict, "Na", 1.25)
  q_conflict <- .rev_scale_ion(q_conflict, "K",  1.25)
  q_conflict <- .rev_scale_ion(q_conflict, "Ca", 0.75)
  q_conflict <- .rev_scale_ion(q_conflict, "Mg", 0.75)

  r_clean <- .rev_lmf_row(add_lmf(q_clean, leachate_data = em$leachate,
                                  reference_data = em$reference, max_sigma_lsi = 1e6))
  r_conf  <- .rev_lmf_row(add_lmf(q_conflict, leachate_data = em$leachate,
                                  reference_data = em$reference, max_sigma_lsi = 1e6))

  ## Precondition: the conflict sample really is overdispersed.
  expect_gt(r_conf$chi2_per_df, 3 * max(r_clean$chi2_per_df, 1e-6))
  ## Contract: overdispersion widens the reported uncertainty.
  expect_gt(r_conf$sigma_lmf, r_clean$sigma_lmf * 1.3)
})


## ----------------------------------------------------------------------------
## Issue 3 — robust reweighting uses raw-residual MAD with k = 1.5. Against a
## MAD over <=9 ions, ~half sit at >=1x MAD, so downweighting fires on nearly
## every sample even under a perfect mixing model (a "pseudo-median" pull).
## Fix: studentize residuals (r_i / sigma_f_i) before the MAD and/or use a less
## aggressive, scale-consistent threshold. Contract: a clean, well-fitting
## sample keeps full weight.
## ----------------------------------------------------------------------------

test_that("[issue 3] a well-fitting sample triggers no robust downweighting", {
  em <- .rev_endmembers()
  set.seed(4)
  ## Only ordinary measurement noise — no ion genuinely departs the mixing line.
  q <- .rev_samples(0.5, "q", "query", noise_sd = 0.02)

  r <- .rev_lmf_row(add_lmf(q, leachate_data = em$leachate,
                            reference_data = em$reference))

  expect_equal(r$n_ions_downweighted, 0L)
})


## ----------------------------------------------------------------------------
## Issue 4 — none of the tracers is leachate-specific, and max_chi2_per_df
## defaults to Inf, so the Huber loop can quietly downweight the dissenting
## high-information ions and converge confidently on a self-consistent Na-Cl
## subset (road salt / saline intrusion masquerading as leachate).
## ----------------------------------------------------------------------------

test_that("[issue 4] an alternative Na-Cl source is surfaced, not silently scored", {
  em <- .rev_endmembers()

  ## Only Na and Cl look leachate-like; every other tracer says reference.
  q <- .rev_samples(0.0, "q", "query", noise_sd = 0)
  q <- .rev_set_f(q, "Na", 0.6)
  q <- .rev_set_f(q, "Cl", 0.6)

  r <- .rev_lmf_row(add_lmf(q, leachate_data = em$leachate,
                            reference_data = em$reference, max_sigma_lsi = 1e6))

  ## The disagreement is visible in chi2 ...
  expect_gt(r$chi2_per_df, 5)
  ## ... and, per the proposed fix, surfaced via a diagnostic that flags the
  ## "robust estimate carried by the low-info majority against the sensitive
  ## tracers" pattern. Proposed contract: an n_high_info_downweighted column.
  expect_gt(r$n_high_info_downweighted, 0L)
})

test_that("[issue 4] a finite max_chi2_per_df rejects the mis-attributed sample", {
  em <- .rev_endmembers()
  q <- .rev_samples(0.0, "q", "query", noise_sd = 0)
  q <- .rev_set_f(q, "Na", 0.6)
  q <- .rev_set_f(q, "Cl", 0.6)

  r <- .rev_lmf_row(add_lmf(q, leachate_data = em$leachate,
                            reference_data = em$reference,
                            max_chi2_per_df = 5, max_sigma_lsi = 1e6))

  expect_true(is.na(r$value))
  expect_match(r$lmf_reason, "poor_fit")
})


## ----------------------------------------------------------------------------
## Issue 5 — BDL handling overstates precision where it matters most. A
## non-detect is substituted at DL/2 but then given sigma_meas ~ 5% RSD. For a
## clean sample the highest-weight ion (total N) is typically a non-detect, so
## the clean end is anchored on an imputed value carrying fictitious precision.
## Honest uncertainty of a censored value is ~ DL/2 (or DL/sqrt(12)).
## ----------------------------------------------------------------------------

test_that("[issue 5] censoring the highest-weight ion lowers precision", {
  em <- .rev_endmembers()

  ## Clean sample (f = 0). Same chemistry twice; in one, the N species are
  ## reported as non-detects (as they typically are downstream).
  base <- .rev_samples(0.0, "q", "query", noise_sd = 0)
  n_species <- c("NH3-N", "NO3-N", "NO2-N")

  q_detected  <- base
  q_nondetect <- base
  q_nondetect$detected[q_nondetect$analyte %in% n_species] <- FALSE

  sig_det <- .rev_lmf_row(add_lmf(q_detected, leachate_data = em$leachate,
                                  reference_data = em$reference,
                                  max_sigma_lsi = 1e6))$sigma_lmf
  sig_nd  <- .rev_lmf_row(add_lmf(q_nondetect, leachate_data = em$leachate,
                                  reference_data = em$reference,
                                  max_sigma_lsi = 1e6))$sigma_lmf

  ## Non-detect total N must reduce confidence, not leave it ~unchanged.
  expect_gt(sig_nd, sig_det * 1.3)
})


## ----------------------------------------------------------------------------
## Issue 6 — collapse_species() sums N species with rowSums(na.rm = TRUE), so a
## sample missing NH4 (often the dominant leachate N species) gets a biased-low
## total_N that then enters as the highest-weight ion at full confidence.
## Contract: a partial total that omits the dominant species is not reported as
## a confident total (NA when NH4 is absent).
## ----------------------------------------------------------------------------

test_that("[issue 6] a partial N sum missing NH4 is not reported at full confidence", {
  meq <- tibble::tribble(
    ~sample_id, ~analyte,   ~value,
    "A",        "NH3-N_",   3.0,
    "A",        "NO3-N_",   0.5,
    "A",        "NO2-N_",   0.1,
    "A",        "Cl_",      2.0,
    "B",        "NO3-N_",   0.5,   # NH4 missing
    "B",        "Cl_",      2.0
  )

  wide <- hydroSense:::collapse_species(meq, id_cols = "sample_id")

  a <- wide$total_N_[wide$sample_id == "A"]
  b <- wide$total_N_[wide$sample_id == "B"]

  expect_equal(a, 3.6)             # complete sum is unchanged
  expect_true(is.na(b))            # partial sum without NH4 -> NA
})


## ----------------------------------------------------------------------------
## Smaller note — the leachate end-member averages per-sample ion:Cl ratios
## (mean-of-ratios). When Cl and the ion carry independent noise this is biased
## HIGH by ~ (1 + CV_Cl^2) relative to ratio-of-means (Jensen). A median ratio
## or ratio-of-means is unbiased.
## ----------------------------------------------------------------------------

test_that("[smaller note] leachate ion:Cl ratio is unbiased under independent noise", {
  set.seed(21)
  n <- 60
  ## True meq ratio Na:Cl = 2. Independent 30% noise on each => mean-of-ratios
  ## inflates by ~1 + 0.3^2 = 1.09; ratio-of-means recovers ~2.
  cl_mg <- pmax(1000 * (1 + stats::rnorm(n, 0, 0.30)), 1)
  na_mg <- pmax(1297 * (1 + stats::rnorm(n, 0, 0.30)), 1)

  mk <- function(a, v, val, mass) tibble::tibble(
    sample_id = paste0("l", seq_len(n)), analyte = a, value = val,
    detected = TRUE, datetime = as.Date("2022-01-01"),
    units.analyte = "mg/L", valence.analyte = v, atomic_mass.analyte = mass
  )
  lea <- dplyr::bind_rows(
    mk("Cl", 1, cl_mg, 35.45),
    mk("Na", 1, na_mg, 22.99)
  )

  em <- hydroSense:::build_endmember_from_override(lea, type = "leachate")
  na_ratio <- em$L_values$mean_ratio[em$L_values$ion == "Na_"]

  expect_lt(na_ratio, 2.10)        # unbiased estimator stays near 2.0
})


## ----------------------------------------------------------------------------
## Bug B1 — if every feature is skipped (no/insufficient reference data),
## group_modify() returns a frame with only site_id, and the subsequent
## left_join(..., by = "sample_id") errors on the missing column instead of
## returning df unchanged.
## ----------------------------------------------------------------------------

test_that("[bug B1] add_lmf returns df unchanged when every feature is skipped", {
  em <- .rev_endmembers()
  q  <- .rev_samples(0.5, "q", "query", noise_sd = 0)

  ## get_reference_site is a host-environment global (not a package binding);
  ## stub it in the global env so the non-override path can run.
  assign("get_reference_site", function(feature_id) list(uuid = "none"),
         envir = globalenv())
  withr::defer(suppressWarnings(rm("get_reference_site", envir = globalenv())))

  ## Valid pooled reference (so informativeness/admission succeed) but every
  ## per-feature reference build returns NULL -> every feature is skipped.
  testthat::local_mocked_bindings(
    build_pooled_reference_endmember = function(...) {
      tibble::tibble(
        ion     = c("Na_", "Cl_", "total_N_"),
        n_ref   = 10L,
        R       = c(0.4, 0.4, 0.05),
        sigma_R = c(0.02, 0.02, 0.005)
      )
    },
    build_reference_endmember = function(...) NULL,
    .package = "hydroSense"
  )

  expect_no_error(
    out <- suppressMessages(
      add_lmf(q, leachate_data = em$leachate, reference_data = NULL)
    )
  )
  expect_false("LMF" %in% out$analyte)   # nothing computed
  expect_equal(nrow(out), nrow(q))       # input preserved
})


## ----------------------------------------------------------------------------
## Bug B2 — build_endmember_from_override() admits F with >= 3 measured samples
## while the standard build_leachate_endmember() requires >= min_leachate_samples
## (default 10). The two paths can disagree on panel composition for the same
## data. Contract: the override applies the same F-inclusion threshold.
## ----------------------------------------------------------------------------

test_that("[bug B2] override leachate applies the standard F-inclusion threshold", {
  set.seed(12)
  ## 5 leachate samples with F measured in all of them: below the standard
  ## threshold (10), so F must be EXCLUDED to match build_leachate_endmember().
  lea <- .rev_samples(rep(1, 5), "lea", "leasite", noise_sd = 0.05)

  em <- hydroSense:::build_endmember_from_override(lea, type = "leachate")

  expect_false("F_" %in% em$L_values$ion)
})
