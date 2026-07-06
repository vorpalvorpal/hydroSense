## ============================================================================
## Shared builders for the LMF code-review TDD specs (test-lmf-review-findings.R)
## ============================================================================
##
## Mirrors the two-component mixing construction in test-lmf.R but exposes the
## extra knobs (noise level, non-detects, dropped/perturbed ions) needed to
## exercise each individual review finding in isolation.
##
## All concentrations are mg/L; valence + atomic mass drive the meq conversion
## inside add_lmf(). Leachate >> reference so the ions are informative.

## Major-ion end-member concentrations (mg/L).
.rev_ions <- function() {
  tibble::tribble(
    ~analyte,          ~valence, ~mass,   ~R,    ~L,
    "Na",               1,       22.99,   10,    1000,
    "K",                1,       39.10,   2,     200,
    "Ca",               2,       40.08,   20,    400,
    "Mg",               2,       24.31,   5,     150,
    "Cl",               1,       35.45,   15,    2000,
    "F",                1,       19.00,   0.2,   5,
    "SO4²⁻",  2,       96.06,   10,    300,
    "NH3-N",            1,       14.01,   0.1,   200,
    "NO3-N",            1,       14.01,   0.5,   5,
    "NO2-N",            1,       14.01,   0.05,  1,
    "CO3-CaCO3",        2,       100.09,  5,     100,
    "HCO3-CaCO3",       2,       100.09,  50,    800
  )
}

## Long-format samples at mixing fractions `fvec` (0 = reference, 1 = leachate)
## with optional multiplicative gaussian noise on each concentration.
.rev_samples <- function(fvec, prefix, site, noise_sd = 0.0,
                         datetime = as.Date("2022-06-01")) {
  ions <- .rev_ions()
  dplyr::bind_rows(lapply(seq_along(fvec), function(i) {
    f     <- fvec[i]
    noise <- if (noise_sd > 0) stats::rnorm(nrow(ions), 0, noise_sd) else 0
    conc  <- ((1 - f) * ions$R + f * ions$L) * (1 + noise)
    tibble::tibble(
      sample_id           = paste0(prefix, "_", i),
      site_id             = site,
      analyte             = ions$analyte,
      value               = pmax(conc, 1e-6),
      detected            = TRUE,
      datetime            = datetime,
      units.analyte       = "mg/L",
      valence.analyte     = ions$valence,
      atomic_mass.analyte = ions$mass
    )
  }))
}

## 12-sample reference + leachate end-members at a given noise level.
.rev_endmembers <- function(noise_sd = 0.02, seed = 101) {
  set.seed(seed)
  list(
    reference = .rev_samples(rep(0, 12), "ref", "refsite", noise_sd),
    leachate  = .rev_samples(rep(1, 12), "lea", "leasite", noise_sd)
  )
}

## Set a single ion in a long sample df to an exact mixing fraction `f`
## (linear in concentration, so f is preserved through the meq conversion).
.rev_set_f <- function(df, ion, f) {
  ions <- .rev_ions()
  r <- ions$R[ions$analyte == ion]
  l <- ions$L[ions$analyte == ion]
  df$value[df$analyte == ion] <- r + f * (l - r)
  df
}

## Multiply a single ion's value (used to push an ion off the mixing line).
.rev_scale_ion <- function(df, ion, factor) {
  df$value[df$analyte == ion] <- df$value[df$analyte == ion] * factor
  df
}

## Pull the single LMF output row from an add_lmf() result.
.rev_lmf_row <- function(out) {
  out[out$analyte == "LMF", , drop = FALSE]
}

## Mark a pending TDD spec: these encode the *target* (post-fix) behaviour and
## are expected to FAIL against the current implementation. Remove the
## .skip_tdd() line at the top of a test to activate it for red-green TDD.
.skip_tdd <- function(issue) {
  testthat::skip(paste0("TDD pending — ", issue,
                        " (remove .skip_tdd() to activate)"))
}
