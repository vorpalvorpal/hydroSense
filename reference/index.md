# Package index

## Example data

- [`leachate_demo()`](https://vorpalvorpal.github.io/leachatetools/reference/leachate_demo.md)
  : Synthetic leachate-impacted water-quality data (demo)

## Leachate detection (LMF)

End-member mixing analysis to quantify the leachate-mixing fraction of a
sample.

- [`add_lmf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_lmf.md)
  : Compute the Leachate Mixing Fraction (LMF) for water quality samples
- [`to_meq()`](https://vorpalvorpal.github.io/leachatetools/reference/to_meq.md)
  : Convert analyte concentrations to milliequivalents per litre

## Chemistry imputation

Bayesian multivariate imputation of below-detection and missing analyte
concentrations, plus deterministic co-analyte imputation.

- [`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)
  : Fit the Bayesian multivariate imputation model(s)
- [`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
  : Impute missing and BDL chemistry using a fitted imputation model
- [`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md)
  : Impute missing normalisation co-analytes from the fitted chemistry
  PCA

## Toxicity & single-substance PAF

Species sensitivity distributions and the potentially affected fraction
for an individual analyte.

- [`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
  : Estimate the fraction of species potentially affected at a
  concentration.
- [`ssd_pct()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_pct.md)
  : Quick point-estimate only (no CI, no bootstrapping).
- [`ssd_hc50()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_hc50.md)
  : Return HC50 for use in Concentration Addition (msPAF).

## Multi-substance PAF (AmsPAF)

Concentration- and response-addition mixture PAF with Added-Risk local
background adjustment.

- [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
  : Compute the Adjusted multi-substance PAF (AmsPAF) for water quality
  samples
- [`prepare_reference()`](https://vorpalvorpal.github.io/leachatetools/reference/prepare_reference.md)
  : Prepare reference chemistry for AmsPAF background subtraction

## Water chemistry & normalisation

Bioavailability and physicochemical adjustments that bring measured
concentrations to the SSD index condition.

- [`correct_ammonia_ph_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/correct_ammonia_ph_temp.md)
  : Correct total ammonia-N to the ANZG reference pH and temperature
- [`derive_hardness()`](https://vorpalvorpal.github.io/leachatetools/reference/derive_hardness.md)
  : Three-way reconciliation between Ca, Mg, and hardness
- [`estimate_water_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/estimate_water_temp.md)
  : Estimate water temperature from air temperature measurements
- [`get_silo_air_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/get_silo_air_temp.md)
  : Fetch daily mean air temperature from SILO for an Australian
  location

## Screening & temporal aggregation

Detection-frequency prescreening and chronic time-weighted averaging.

- [`prescreen_analytes()`](https://vorpalvorpal.github.io/leachatetools/reference/prescreen_analytes.md)
  : Pre-screen analytes by detection frequency
- [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/leachatetools/reference/time_weighted_aggregate.md)
  : Time-weighted chronic aggregation for any long-format value column
- [`expand_focal_dates()`](https://vorpalvorpal.github.io/leachatetools/reference/expand_focal_dates.md)
  : Generate a sequence of focal dates for chronic AmsPAF computation
