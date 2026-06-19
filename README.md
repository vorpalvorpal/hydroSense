
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hydroSense

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/vorpalvorpal/hydroSense/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vorpalvorpal/hydroSense/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/vorpalvorpal/hydroSense/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/vorpalvorpal/hydroSense/actions/workflows/pkgdown.yaml)
[![brms-smoke](https://github.com/vorpalvorpal/hydroSense/actions/workflows/brms-smoke.yaml/badge.svg)](https://github.com/vorpalvorpal/hydroSense/actions/workflows/brms-smoke.yaml)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**hydroSense** analyses water quality in freshwater systems
potentially impacted by landfill leachate. It turns routine monitoring
chemistry into ecologically interpretable, mixture-aware risk estimates,
with explicit handling of below-detection data and local background
conditions.

The package has three pillars:

- **Chemistry imputation.** A domain-agnostic Bayesian multivariate
  model (`fit_imputation_model()` / `impute_chemistry()`) imputes
  below-detection and missing analyte concentrations while preserving
  cross-analyte correlation structure, and can return full posterior
  draws. Analyte groups are configured with `impute_group()`;
  `leachate_impute_groups()` is the bundled leachate preset.
- **Leachate detection (LMF).** End-member mixing analysis (`add_lmf()`)
  estimates the *leachate-mixing fraction* of each sample from its
  major-ion signature, distinguishing leachate influence from natural
  variation.
- **Multi-substance toxicity (msPAF).** Species sensitivity
  distributions (`ssd_paf()`) feed a multi-substance
  potentially-affected-fraction (`add_mspaf()`) that combines toxicants
  by concentration- and response-addition and subtracts the *local*
  background via the Added Risk Approach. The background can be a static
  reference or a *contemporaneous* one (`fit_reference_model()`) that
  predicts what the reference site would have shown at each sample’s
  moment from hydrology and season; a season-blind companion
  (`fit_target_model()`) predicts the impacted site between grabs, so
  `mspaf_daily()` can return a continuous daily risk series.

Toxicity is grounded in the ANZECC/ANZG freshwater guideline datasets,
with bioavailability and physicochemical normalisations (hardness, pH,
DOC, temperature) applied so that every sample is compared on the same
index-condition basis.

## Installation

You can install the development version from
[GitHub](https://github.com/vorpalvorpal/hydroSense):

``` r
# install.packages("pak")
pak::pak("vorpalvorpal/hydroSense")
```

**brms** is a hard dependency (`Imports`) and installs with the package,
but the Bayesian imputation step additionally needs a working
[Stan](https://mc-stan.org/) toolchain, which is not pulled in
automatically; the toxicity and LMF pillars do not need Stan.

## Usage

A typical msPAF assessment takes long-format monitoring data (one row
per sample × analyte) and returns a per-sample multi-substance PAF,
adjusted for local background:

``` r
library(hydroSense)

# Point the toxicity engine at the ANZG guideline-data folder once per session.
options(hydroSense.guideline_dir = "path/to/guideline data")

mspaf <- add_mspaf(
  df        = monitoring_long,   # sample_id, site_id, analyte, value, ...
  reference = reference_long     # local upstream / background chemistry
)
```

Single-substance PAF and SSD quantities are available directly:

``` r
# Fraction of species potentially affected by 9.3 mg/L total ammonia-N
# (at the pH 7.0 / 20 degC index condition):
ssd_paf("NH3-N", conc_ug_L = 9321)

# Concentration affecting 50% of species (used as the Toxic-Unit denominator):
ssd_hc50("Cu")
```

## Learn more

- [Leachate Mixing
  Fraction](https://vorpalvorpal.github.io/hydroSense/articles/leachate-mixing-fraction.html)
  — what fraction of a sample’s chemistry can be explained as a mixture
  of local reference water and landfill leachate.
- [Imputation](https://vorpalvorpal.github.io/hydroSense/articles/imputation.html)
  — how missing and below-detection values are imputed.
- [Analyte
  normalisation](https://vorpalvorpal.github.io/hydroSense/articles/normalisation.html)
  — how bioavailability and physicochemical corrections bring
  concentrations to the SSD index condition.
- [msPAF](https://vorpalvorpal.github.io/hydroSense/articles/chronic-mspaf-interpretation.html)
  — what the returned PAF values mean for an environmental assessment.
- [Function
  reference](https://vorpalvorpal.github.io/hydroSense/reference/index.html)
  — all exported functions, grouped by pillar.

## License

MIT © Robin Shannon. Toxicity data are derived from the ANZECC/ANZG
freshwater guideline materials; see those sources for their own
licensing terms.
