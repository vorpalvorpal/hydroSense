# Pivot long chemistry rows to a wide per-sample predictor frame

Shared by
[`.prepare_chem_pca()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-prepare_chem_pca.md)
(training) and
[`.compute_pca_scores()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_pca_scores.md)
(scoring) so predictor construction is identical on both paths.
Below-detection cells (`detected == FALSE`) for concentration-like
analytes are set to `NA` before collapsing, so a BDL predictor cell is
treated as missing (scored by NIPALS from the sample's observed
predictors) rather than substituted at its detection limit or DL/2 —
this package avoids substituting magic numbers for non-detects. The
interval-scale `no_log_vars` (pH, temperature, ORP, DO) are kept as-is
when BDL: "below DL" is a concentration idea and they are essentially
never BDL. Duplicate `(sample, analyte)` rows collapse via geometric
mean for concentration-like analytes (matching how targets are logged
before collapsing) and arithmetic mean for `no_log_vars`. An all-`NA`
collapse yields `NaN` from
[`mean()`](https://rdrr.io/r/base/mean.html)/`exp(NaN)`; coerced back to
`NA` so it reaches NIPALS scoring as missing, not a bogus numeric value
(bug B3). Callers without a `detected` column (e.g. the hydro-layer WQ
frames in `R/target_model.R`, which carry no BDL concept) are treated as
fully detected — no halving.

## Usage

``` r
.pivot_chem_wide(df, wq_vars, no_log_vars)
```
