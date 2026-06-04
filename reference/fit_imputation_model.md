# Fit the Bayesian multivariate imputation model(s)

Fits one or two `brms` multivariate GAMs — one for metals and one for
organics — using a PCA-compressed water-quality (WQ) block as additional
environmental predictors. The returned model object is passed to
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
for prediction on new data.

## Usage

``` r
fit_imputation_model(
  df,
  pca_vars = NULL,
  required_vars = c("pH", "EC"),
  metal_analytes = NULL,
  doc_like_analytes = NULL,
  min_detect_freq = 0.05,
  min_samples = 10L,
  min_var_explained = 0.75,
  max_pcs = 6L,
  family = "gaussian",
  impute_method = c("rescor_mi", "cens", "cens_factor"),
  iter = 2000,
  warmup = 1000,
  chains = 4,
  cores = parallel::detectCores(),
  save_dir = NULL,
  ...
)
```

## Arguments

- df:

  Long-format chemistry data frame with columns `sample_id`, `site_id`,
  `datetime`, `analyte`, `value`, `detected`.

- pca_vars:

  Analyte names to include in the unified chemistry PCA (used as
  predictors for the brms model via PC scores). Default:
  `c("pH", "EC", "NH3-N")` plus all `.WQ_BLOCK_CANDIDATES`.
  Normalisation co-analytes (DOC, Ca, Mg) are included in the default
  set.

- required_vars:

  Analyte names that must be present in a sample for it to be retained
  in training and prediction. Default `c("pH", "EC")`. Samples missing
  any of these are dropped entirely.

- metal_analytes:

  Analyte names classified as metals. Default `.METAL_ANALYTES`.

- doc_like_analytes:

  Analyte names used for the organics hurdle check (the "organic carbon
  present" requirement). Default `.DOC_LIKE_ANALYTES`.

- min_detect_freq:

  Minimum detection frequency for a PCA variable to be retained. Default
  `0.05`. Required vars are always retained regardless.

- min_samples:

  Minimum training samples after required-var filtering.

- min_var_explained:

  Target cumulative variance for PCA axis selection. Default `0.75`.

- max_pcs:

  Maximum PCA axes to use. Default `6L`.

- family:

  brms response family. Must be `"gaussian"` (concentrations are
  log-transformed before fitting; residual correlations require Gaussian
  family).

- impute_method:

  How below-detection (BDL) values and cross-analyte coupling are
  handled. One of:

  `"rescor_mi"`

  :   (default) Residual correlation across analytes (`rescor = TRUE`)
      with BDL/missing treated as imputable (`mi()`); the imputed BDL
      cells are capped at the detection limit post-hoc by
      [`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
      (brms cannot combine `rescor` with `cens()`).

  `"cens"`

  :   Proper left-censoring of BDL at the detection limit
      (`cens("left")`), no residual correlation – clean BDL handling but
      no cross-analyte coupling.

  `"cens_factor"`

  :   As `"cens"` plus a shared per-sample latent factor
      (`(1 | sample_id)` correlated across analytes), which
      re-introduces cross-analyte coupling while keeping proper
      censoring.

  See
  [`vignette("imputation")`](https://vorpalvorpal.github.io/leachatetools/articles/imputation.md)
  and the package benchmark for guidance on which to prefer.

- iter, warmup, chains, cores:

  brms MCMC settings.

- save_dir:

  If non-NULL, save the returned model object as a `.qs` file in this
  directory using
  [`qs2::qs_save()`](https://rdrr.io/pkg/qs2/man/qs_save.html).

- ...:

  Additional arguments passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

## Value

A named list of class `"imputation_model"`:

- `$pca`: PCA fit + metadata (loadings, training medians, n_pcs, …)

- `$metals`: list with `$fit` (brmsfit), `$analytes`, `$safe_names`

- `$organics`: same structure, or `NULL` if no organics pass prescreen

- `$required_vars`, `$pca_vars`, `$hurdle_metals`, `$hurdle_organics`

- `$fit_date`, `$n_samples`: metadata If `save_dir` is supplied, the
  path to the saved file is returned as `attr(result, "save_path")`.

## Details

**Model structure**

For each analyte group (metals / organics), the mean structure is:

    s(PC1) + s(PC2) + ... + s(PCk)

where `PC*` are the leading principal components of the unified
chemistry PCA (see *Chemistry PCA* below). All target analytes (metals
or organics) are modelled jointly with `rescor = TRUE`, so observed
co-analytes at a given sample condition the posterior of the missing
ones through the residual correlation matrix.

**Why `rescor = TRUE`** — the PCA captures the *instantaneous* chemical
covariance structure (what is measured together at a single moment), but
it cannot capture the *temporal-lag* covariance characteristic of
AMD/leachate-impacted aquifers, where conservative tracers move ahead of
redox-controlled metals. At a post-pulse sample the PCA scores have
returned toward baseline but Cu/Pb/Zn/Mn remain elevated together — that
co-elevation is pure residual correlation with no predictor signal
driving it. `rescor = TRUE` is the right machinery for this and is what
makes multivariate imputation borrow strength across analytes.

**Costs of `rescor = TRUE`** — brms cannot combine `set_rescor(TRUE)`
with `cens("left")`, so this implementation uses `mi()` for BDL values
and applies a post-hoc cap (see
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)).
The cap clips imputed BDL cells to the original detection limit when the
model predicts above DL. For sites where the chemistry context
legitimately suggests high concentrations the cap can fire frequently;
results in that regime should be inspected. Three alternative
configurations are worth benchmarking on real hold-out data if
predictive performance becomes a concern:

- `rescor = TRUE` + `mi()` (current; expected to win on plume-affected
  groundwater because cross-analyte residual coupling captures plume
  dynamics that the predictor PCA misses).

- `rescor = FALSE` + `cens("left")` (statistically clean for BDL; loses
  cross-analyte residual coupling).

- `rescor = FALSE` + `cens("left")` + shared `(1 | sample_id)` (proper
  BDL handling with rank-1 latent-factor coupling across analytes).
  Benchmark methodology: mask 10% of detected cells, fit each
  configuration, compare hold-out RMSE / coverage.

**Chemistry PCA**

All `pca_vars` — major ions, pH, EC, NH3-N, DOC, nutrients, redox
indicators — that are present in `df` and pass a detection-frequency
check are submitted to
[`nipals::nipals()`](https://kwstat.github.io/nipals/reference/nipals.html),
which handles within-sample missing cells natively without prior
imputation. Using a single unified PCA (rather than separate driver +
WQ-block sets) eliminates predictor collinearity and ensures
normalisation co-analytes (DOC, Ca, Mg) influence the imputed metal
concentrations. Principal components are added until cumulative variance
explained reaches `min_var_explained` or `max_pcs` is reached. A minimum
of two PCs is always used.

**Hurdles (applied at prediction time by
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md))**

- *Metals*: a sample is only imputed if at least one metal analyte is
  present (detected or BDL) in `df` for that sample.

- *Organics*: a sample is only imputed if at least one of DOC, TOC, BOD,
  COD or cBOD is present.

**BDL required variables**

When a `required_vars` analyte (pH or EC) is below the detection limit
for a sample, the stored detection-limit value is used as-is
(conservative upper bound). A message is issued but the sample is
retained.

## See also

[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires a Stan toolchain (brms). Fit once, then reuse for imputation.
model <- fit_imputation_model(monitoring_long)
draws <- impute_chemistry(monitoring_long, model, return = "draws")
} # }
```
