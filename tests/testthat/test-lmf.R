## End-to-end test of the LMF (leachate mixing fraction) pillar via add_lmf(),
## using the caller-supplied end-member override path (leachate_data /
## reference_data) so no dashboard infrastructure is needed.
##
## Construction: a two-component linear mixing model. Query samples are known
## mixtures of a reference and a leachate end-member; add_lmf() should recover
## the mixing fraction (as a 0-100 percentage).

library(testthat)
library(leachatetools)

# Major-ion end-member concentrations (mg/L). Leachate >> reference so the
# ions are informative. valence/atomic mass drive the meq conversion.
.lmf_ions <- function() {
  tibble::tribble(
    ~analyte,      ~valence, ~mass,   ~R,    ~L,
    "Na",           1,       22.99,   10,    1000,
    "K",            1,       39.10,   2,     200,
    "Ca",           2,       40.08,   20,    400,
    "Mg",           2,       24.31,   5,     150,
    "Cl",           1,       35.45,   15,    2000,
    "F",            1,       19.00,   0.2,   5,
    "SO4²⁻", 2,    96.06,   10,    300,
    "NH3-N",        1,       14.01,   0.1,   200,
    "NO3-N",        1,       14.01,   0.5,   5,
    "NO2-N",        1,       14.01,   0.05,  1,
    "CO3-CaCO3",    2,       100.09,  5,     100,
    "HCO3-CaCO3",   2,       100.09,  50,    800
  )
}

# Build long-format samples at mixing fractions `fvec` (0 = reference,
# 1 = leachate), with small multiplicative noise.
.lmf_samples <- function(fvec, prefix, site, noise_sd = 0.02) {
  ions <- .lmf_ions()
  dplyr::bind_rows(lapply(seq_along(fvec), function(i) {
    f    <- fvec[i]
    conc <- ((1 - f) * ions$R + f * ions$L) *
            (1 + stats::rnorm(nrow(ions), 0, noise_sd))
    tibble::tibble(
      sample_id           = paste0(prefix, "_", i),
      site_id             = site,
      analyte             = ions$analyte,
      value               = pmax(conc, 1e-6),
      detected            = TRUE,
      datetime            = as.Date("2022-06-01"),
      units.analyte       = "mg/L",
      valence.analyte     = ions$valence,
      atomic_mass.analyte = ions$mass
    )
  }))
}

.lmf_endmembers <- function() {
  set.seed(101)
  list(
    reference = .lmf_samples(rep(0, 12), "ref", "refsite"),
    leachate  = .lmf_samples(rep(1, 12), "lea", "leasite")
  )
}

test_that("add_lmf recovers the two-component mixing fraction", {
  em    <- .lmf_endmembers()
  fvec  <- c(0, 0.25, 0.5, 0.75, 1.0)
  set.seed(7)
  query <- .lmf_samples(fvec, "q", "query")

  out <- add_lmf(query, leachate_data = em$leachate, reference_data = em$reference)
  lmf <- out[out$analyte == "LMF", ]
  lmf <- lmf[order(lmf$sample_id), ]

  expect_equal(nrow(lmf), length(fvec))
  expect_true(all(is.na(lmf$lmf_reason)))            # all quantified
  # Recovered LMF tracks 100 * mixing fraction.
  expect_equal(lmf$value, 100 * fvec, tolerance = 3)
})

test_that("LMF is monotonic increasing with leachate fraction", {
  em    <- .lmf_endmembers()
  set.seed(11)
  query <- .lmf_samples(seq(0, 1, by = 0.1), "q", "query")
  out   <- add_lmf(query, leachate_data = em$leachate, reference_data = em$reference)
  lmf   <- out[out$analyte == "LMF", ]
  lmf   <- lmf[order(as.integer(sub("q_", "", lmf$sample_id))), ]
  expect_true(all(diff(lmf$value) > 0))
})

test_that("add_lmf appends LMF rows and preserves the input unchanged", {
  em    <- .lmf_endmembers()
  set.seed(7)
  query <- .lmf_samples(c(0, 1), "q", "query")
  out   <- add_lmf(query, leachate_data = em$leachate, reference_data = em$reference)

  # Original rows are all retained.
  expect_equal(sum(out$analyte != "LMF"), nrow(query))
  lmf <- out[out$analyte == "LMF", ]
  expect_equal(nrow(lmf), 2L)                        # one per sample
  expect_true(all(lmf$units.analyte == "%"))
  expect_true(all(lmf$detected))
  expect_false(any(is.na(lmf$datetime)))             # datetime propagated
  expect_true(all(is.finite(lmf$sigma_lmf)))
  expect_true(all(lmf$n_ions_used > 0))
})

test_that("a sample with too few high-information ions fails gracefully", {
  em    <- .lmf_endmembers()
  set.seed(7)
  # Keep only Ca + Mg (low-information major ions): below min_high_info_ions.
  query <- .lmf_samples(0.5, "q", "query")
  query <- query[query$analyte %in% c("Ca", "Mg"), ]
  out   <- add_lmf(query, leachate_data = em$leachate, reference_data = em$reference)
  lmf   <- out[out$analyte == "LMF", ]
  expect_equal(nrow(lmf), 1L)
  expect_true(is.na(lmf$value))
  expect_match(lmf$lmf_reason, "insufficient_high_info_ions")
})

test_that("robust reweighting downweights a single grossly outlying ion", {
  em    <- .lmf_endmembers()
  set.seed(7)
  query <- .lmf_samples(0.5, "q", "query")
  # Spike Cl far above its mixing-consistent value: a transformation/outlier ion.
  query$value[query$analyte == "Cl"] <- query$value[query$analyte == "Cl"] * 6
  out   <- add_lmf(query, leachate_data = em$leachate, reference_data = em$reference)
  lmf   <- out[out$analyte == "LMF", ]
  expect_gt(lmf$n_ions_downweighted, 0)
  # Robust estimate should remain near the true 50% despite the outlier.
  expect_equal(lmf$value, 50, tolerance = 12)
})
