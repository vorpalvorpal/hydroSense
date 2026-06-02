# leachatetools

**leachatetools** analyses water quality in freshwater systems
potentially impacted by landfill leachate. It turns routine monitoring
chemistry into ecologically interpretable, mixture-aware risk estimates,
with explicit handling of below-detection data and local background
conditions.

The package has three pillars:

- **Leachate detection (LMF).** End-member mixing analysis
  ([`add_lmf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_lmf.md))
  estimates the *leachate-mixing fraction* of each sample from its
  major-ion signature, distinguishing leachate influence from natural
  variation.
- **Chemistry imputation.** A Bayesian multivariate model
  ([`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)
  /
  [`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md))
  imputes below-detection and missing analyte concentrations while
  preserving cross-analyte correlation structure, and can return full
  posterior draws.
- **Multi-substance toxicity (AmsPAF).** Species sensitivity
  distributions
  ([`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md))
  feed a multi-substance potentially-affected-fraction
  ([`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md))
  that combines toxicants by concentration- and response-addition and
  subtracts the *local* background via the Added Risk Approach.

Toxicity is grounded in the ANZECC/ANZG freshwater guideline datasets,
with bioavailability and physicochemical normalisations (hardness, pH,
DOC, temperature) applied so that every sample is compared on the same
index-condition basis.

## Installation

You can install the development version from
[GitHub](https://github.com/vorpalvorpal/leachatetools):

``` r

# install.packages("pak")
pak::pak("vorpalvorpal/leachatetools")
```

Bayesian imputation additionally requires a working
[Stan](https://mc-stan.org/) toolchain via **brms** (a `Suggests`
dependency); the toxicity and LMF pillars do not.

## Usage

A typical AmsPAF assessment takes long-format monitoring data (one row
per sample × analyte) and returns a per-sample multi-substance PAF,
adjusted for local background:

``` r

library(leachatetools)

# Point the toxicity engine at the ANZG guideline-data folder once per session.
options(leachatetools.guideline_dir = "path/to/guideline data")

amspaf <- add_amspaf(
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

- [Analyte
  normalisation](https://vorpalvorpal.github.io/leachatetools/articles/normalisation.html)
  — how bioavailability and physicochemical corrections bring
  concentrations to the SSD index condition.
- [Interpreting chronic AmsPAF
  outputs](https://vorpalvorpal.github.io/leachatetools/articles/chronic-amspaf-interpretation.html)
  — what the returned PAF values mean for an environmental assessment.
- [Function
  reference](https://vorpalvorpal.github.io/leachatetools/reference/index.html)
  — all exported functions, grouped by pillar.

## License

MIT © Robin Shannon. Toxicity data are derived from the ANZECC/ANZG
freshwater guideline materials; see those sources for their own
licensing terms.
