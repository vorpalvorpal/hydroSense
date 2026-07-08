# Fit the Bayesian multivariate imputation model(s)

Fits a `brms` multivariate GAM for each **imputation group** (see
[`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md)),
using a PCA-compressed water-quality (WQ) block as additional
environmental predictors. The engine itself is domain-agnostic: the
leachate-specific groups (a `metals` group and a catch-all `organics`
group) are supplied by the default `groups = leachate_impute_groups()`
and can be swapped for any other chemistry by passing your own list of
[`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md)
objects. The returned model object is passed to
[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)
for prediction on new data.

## Usage

``` r
fit_imputation_model(
  df,
  pca_vars = NULL,
  required_vars = c("pH", "EC"),
  groups = NULL,
  exclude = NULL,
  no_log_vars = NULL,
  min_detect_freq = 0.05,
  min_target_detect_freq = 0.05,
  min_target_detect_n = 4L,
  min_samples = 10L,
  min_var_explained = 0.75,
  max_pcs = 6L,
  family = "gaussian",
  impute_method = c("marginal", "rescor_mi", "cens", "cens_factor", "factor"),
  iter = 2000,
  warmup = 1000,
  chains = 4,
  cores = parallel::detectCores(),
  k = NULL,
  seed = NULL,
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

- groups:

  A list of
  [`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md)
  objects describing which analytes to impute and how each group is
  hurdled. Default `NULL` uses
  [`leachate_impute_groups()`](https://vorpalvorpal.github.io/hydroSense/reference/leachate_impute_groups.md)
  (a `metals` group plus a catch-all `organics` group). Supply your own
  list to impute a different chemistry.

- exclude:

  Analyte names that must never be modelled as response variables in any
  group (e.g. counts, qualitative descriptors). Default `NULL` uses
  [.IMPUTE_EXCLUDED](https://vorpalvorpal.github.io/hydroSense/reference/dot-IMPUTE_EXCLUDED.md).

- no_log_vars:

  Analyte names that must **not** be log-transformed before the
  chemistry PCA (interval-scale or already-logarithmic variables such as
  pH and temperature). Default `NULL` uses
  [.PCA_NO_LOG_VARS](https://vorpalvorpal.github.io/hydroSense/reference/dot-PCA_NO_LOG_VARS.md).

- min_detect_freq:

  Minimum detection frequency for a PCA variable to be retained. Default
  `0.05`. Required vars are always retained regardless.

- min_target_detect_freq:

  Minimum detection frequency (fraction of samples in which the analyte
  is *detected*) for a metal/organic to be included as an imputation
  target. Targets below this are dropped (they have too few detections
  to model and would otherwise inflate the model on near-all-BDL
  panels). Default `0.05`. Combined with `min_target_detect_n` (both
  gates must pass).

- min_target_detect_n:

  Minimum **absolute** number of distinct samples in which a
  metal/organic is detected for it to be included as an imputation
  target. The fraction gate above already implies a count floor of
  `min_target_detect_freq * n_samples`, but that scales with dataset
  size and collapses on small datasets; this absolute floor guarantees
  enough anchors to constrain the fit regardless of dataset size.
  Default `4L` (non-binding on typical panels, where the fraction gate
  dominates).

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
  handled. Defaults to `"marginal"` — the fastest, best- calibrated
  method on the package benchmark; the others are retained for
  comparison and for panels with genuinely strong cross-analyte
  coupling. See
  [`vignette("imputation")`](https://vorpalvorpal.github.io/hydroSense/articles/imputation.md)
  for the assumptions, trade-offs and benchmark. One of:

  `"marginal"`

  :   **(default)** Per-analyte left-censored
      [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (spline mean
      on the PC scores) with an independent Student-t posterior
      predictive – **no** cross-analyte borrowing. Uses mgcv only (no
      brms, no cmdstanr) and is ~400x faster than `"rescor_mi"`. Best
      MAE and well-calibrated 90% coverage on the B.S01 benchmark; on
      that panel the cross-metal borrowing the other methods attempt
      does not help (it mis-conditions sparse analytes). Uncertainty is
      a proper posterior predictive: GAM parameter uncertainty
      (`beta ~ N(coef, Vp)`, unconditional) plus residual-variance
      uncertainty (`sigma^2` drawn from its scaled-inverse-chi^2
      posterior, giving the t-tails), with BDL cells truncated at the
      detection limit (no post-hoc cap).

  `"rescor_mi"`

  :   Residual correlation across analytes (`rescor = TRUE`) with
      BDL/missing treated as imputable (`mi()`); the imputed BDL cells
      are capped at the detection limit post-hoc by
      [`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)
      (brms cannot combine `rescor` with `cens()`). Slightly better
      hold-out RMSE than `"marginal"` but at ~400x the cost and with
      over-wide, uninformative intervals; the `mi()` + correlation
      geometry is funnel-prone, so `adapt_delta = 0.95` and an `lkj(2)`
      prior on the residual correlation are set by default (override via
      `control` / `prior` in `...`). The geometry stays hard (tree-depth
      saturation, low E-BFMI, worst-case R̂ ≈ 1.6 on a hard mask): trust
      the **point estimate**, but check
      [`brms::rhat()`](https://mc-stan.org/posterior/reference/rhat.html)
      before relying on the **draws**.

  `"cens"`

  :   Proper left-censoring of BDL at the detection limit
      (`cens("left")`), no residual correlation – clean BDL handling but
      no cross-analyte coupling.

  `"cens_factor"`

  :   Proper left-censoring **with** cross-analyte coupling. Fitted as a
      single long-format model with a shared per-sample latent factor
      `(1 | sample_id)` common to all analytes, so an observed metal
      informs the unobserved/BDL ones at that sample. The factor is
      well-identified (each sample contributes several analyte
      observations); `adapt_delta = 0.95` is set by default to clear the
      factor's mild funnel (override via `control` in `...`).

  `"factor"`

  :   **(Route C)** Low-rank left-censored latent factor model
      (`Sigma = Lambda Lambda' + Psi`, rank `k = 2` by default),
      resolving findings 1-3 at the source: BDL cells are censored in
      the likelihood (no post-hoc cap), and the latent factor is
      inferred from a sample's *observed* metals at prediction time, so
      measured metals genuinely condition the missing/BDL ones. Fitted
      in two stages: a per-analyte
      [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) mean on the
      PC scores, then a Stan factor model on the residuals (needs
      cmdstanr). See `dev/plan-route-c.md`.

  See
  [`vignette("imputation")`](https://vorpalvorpal.github.io/hydroSense/articles/imputation.md)
  and the package benchmark for guidance on which to prefer.

- iter, warmup, chains, cores:

  brms MCMC settings.

- k:

  Number of latent factors, only used when `impute_method = "factor"`
  (Route C). `NULL` (default) uses `min(2, J - 1)` per group, where `J`
  is that group's number of target analytes (capped below `J` for
  identifiability — see `dev/plan-route-c.md`). Choosing `k = 2` vs
  `k = 3` by held-out coverage is a deliberate validation step; raise it
  here once you have evidence a group's analytes need a second/third
  shared axis. Ignored (with a message) for a single-analyte group,
  which always falls back to a Stage-1-only marginal fit regardless of
  `k`. Ignored by the brms methods.

- seed:

  Optional integer seed, only used when `impute_method = "factor"`
  (Route C): passed to the Stage-2 Stan sampler (`cmdstanr`'s
  `$sample(seed = ...)`) for reproducible factor fits. `NULL` (default)
  samples with a random seed each call. The brms methods are already
  seedable via `seed` in `...` (passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html));
  this argument only affects the factor method.

- save_dir:

  If non-NULL, save the returned model object as a `.qs` file in this
  directory using
  [`qs2::qs_save()`](https://rdrr.io/pkg/qs2/man/qs_save.html).

- ...:

  Additional arguments passed to the sampler. For the three brms methods
  they go to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html); the
  Stan **`backend`** defaults to `"cmdstanr"` when the cmdstanr package
  and a CmdStan install are both available (cached compiled binaries +
  faster sampling, statistically equivalent to rstan), otherwise
  `"rstan"`; pass `backend = ...` here or set
  `options(hydroSense.brms_backend = ...)` to override. For
  `impute_method = "factor"` (Route C) there is no brms call — `...` is
  passed to cmdstanr's `$sample()` instead (e.g. `adapt_delta`,
  `max_treedepth`, `parallel_chains`), so brms-only arguments do not
  apply.

## Value

A named list of class `"imputation_model"`:

- `$pca`: PCA fit + metadata (loadings, training medians, n_pcs,
  `no_log_vars`, …)

- `$groups`: a named list (one entry per fitted group, named by the
  group's `name`); each entry has `$fit` (brmsfit), `$analytes`,
  `$safe_names`, `$name`, `$hurdle`. Empty if no group had any
  modellable analytes.

- `$group_specs`: the input list of
  [`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md)
  objects

- `$required_vars`, `$pca_vars`, `$exclude`, `$impute_method`

- `$fit_date`, `$n_samples`: metadata If `save_dir` is supplied, the
  path to the saved file is returned as `attr(result, "save_path")`.

## Details

**Model structure**

For each group, the mean structure is:

    s(PC1) + s(PC2) + ... + s(PCk)

where `PC*` are the leading principal components of the unified
chemistry PCA (see *Chemistry PCA* below). A group's target analytes are
modelled jointly with `rescor = TRUE`, so observed analytes at a given
sample condition the posterior of the missing ones through the residual
correlation matrix. Separate groups are fitted as separate models and do
not share residual correlation.

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
[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)).
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
of `min(2, available components)` PCs is used — i.e. the floor of two
only applies when at least two components exist.

**Hurdles (applied at prediction time by
[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md))**

Each group may carry a *presence hurdle* (see
[`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md)):
a sample is only imputed for that group if it carries at least one of
the hurdle analytes (detected or BDL). Under the leachate preset the
metals group is hurdled on metal presence and the organics group on
DOC-like presence.

**BDL required variables**

When a `required_vars` analyte (pH or EC) is below the detection limit
for a sample, the stored detection-limit value is used as-is
(conservative upper bound). A message is issued but the sample is
retained.

## See also

[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md),
[`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md),
[`leachate_impute_groups()`](https://vorpalvorpal.github.io/hydroSense/reference/leachate_impute_groups.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires a Stan toolchain (brms). Fit once, then reuse for imputation.
model <- fit_imputation_model(monitoring_long)
draws <- impute_chemistry(monitoring_long, model, return = "draws")

# A different domain: two custom groups instead of the leachate preset.
model2 <- fit_imputation_model(
  monitoring_long,
  groups = list(
    impute_group("trace_metals", targets = c("As", "Cd", "Pb"),
                 hurdle = c("As", "Cd", "Pb")),
    impute_group("nutrients", targets = NULL, hurdle = c("NO3-N", "P-total"))
  )
)
} # }
```
