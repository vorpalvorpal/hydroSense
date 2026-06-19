# hydroSense (development)

## Rename, dependency, and a one-call pipeline (breaking)

* **Package renamed `leachatetools` → `hydroSense`.** Update `library()` calls
  and the option `leachatetools.guideline_dir` → `hydroSense.guideline_dir`.
* **`AmsPAF` renamed to `msPAF` throughout** — function names
  (`add_amspaf()` → `add_mspaf()`, `amspaf_daily()` → `mspaf_daily()`,
  `classify_amspaf_tier()` → `classify_mspaf_tier()`), the emitted analyte
  label (`"AmsPAF"` → `"msPAF"`), and result columns (`amspaf*` → `mspaf*`).
  No backwards-compatibility shims.
* **`brms` moved from `Suggests` to `Imports`** — it now installs with the
  package (a working Stan toolchain is still required separately for fitting).
* **New `mspaf_pipeline()`** bundles `fit_imputation_model()` →
  `fit_reference_model()` → `mspaf_daily()` into one call, with imputation
  **on by default**. It threads one imputation model through both downstream
  steps. Imputation moves the msPAF result (a MAR/MNAR modelling choice), so
  pass `impute = FALSE` for the non-imputed calculation.

## Imputation completes the analyte panel (issue #53)

* `impute_chemistry()` now imputes analytes that are **entirely absent** for an
  eligible sample, not just BDL/censored cells on rows that already exist. Each
  fabricated cell becomes a new row tagged `imputed_kind = "missing"`,
  `detected = TRUE`, carrying the sample's `site_id`/`datetime`. This activates
  the documented "anchor a sparse analyte from a well-sampled one" behaviour
  (e.g. Zn-only samples gain a model-predicted Cu row) that the merge had
  silently dropped. Prefer `return = "draws"` when these values feed a further
  model, so the imputation uncertainty propagates.

## Daily gap-uncertainty bracket: informative vs ignorable envelopes (issue #50)

* `mspaf_daily()` gains a `gap_uncertainty` argument
  (`"bracket"` (default), `"ignorable"`, `"informative"`). In observation gaps
  the latent Kalman residual reverts to its marginal variance, widening the
  band — the honest posterior only under **ignorable (MAR)** missingness. Field
  sampling is often **informative (MNAR)** (gaps exist *because* the system was
  judged quiescent), so the band is now **bracketed** between two extremes: the
  ignorable envelope keeps the simulation-smoother draw (previous behaviour),
  while the informative envelope freezes the residual at its posterior mean on
  in-gap days. The two are nested and coincide at observation days; their
  separation is the cost of the ignorability assumption (Rubin 1976).
* **Breaking (pre-v1):** draws-mode output schema changes. `"summary"` now
  returns a `deterministic` centre line plus `median_*`/`lo_*`/`hi_*` columns
  per envelope (and `precautionary_lo`/`precautionary_hi` in `"bracket"` mode —
  a labelled *decision* bound `[lo_informative, hi_ignorable]`, not a calibrated
  interval). `"draws"` returns `mspaf_ignorable`/`mspaf_informative` per draw.
  The old `mspaf`/`mspaf_lower`/`mspaf_upper` columns are removed.
* `mspaf_daily()` emits a one-time CLI warning, when uncertainty is requested
  against rainfall hydrology, that gaps are treated as ignorable and the
  intervals over-state uncertainty in gaps known to be quiescent.

## Free recession constant in the rainfall API (issue #49)

* The Antecedent Precipitation Index is now an **exact recursive linear
  reservoir** (`S_t = k^dt * S_{t-1} + P_t`, `k = exp(-1/tau)`; Maillet 1905,
  Kohler & Linsley 1951) with no truncation horizon, evaluated in C via
  `stats::filter()`. The recession constant `tau` is selected per analyte by
  profiled GAM AIC (golden-section, with a `tau_long >= 1.5*tau_short`
  separation and a delta-AIC >= 2 adoption gate over a parsimonious default).
* **Breaking:** `fit_reference_model()` and `fit_target_model()` arguments
  `api_windows_short`/`api_windows_long` are replaced by
  `api_tau_bounds_short`/`api_tau_bounds_long` (length-2 `c(lo, hi)` day
  ranges; defaults `c(1, 30)` and `c(20, 365)`). Per-analyte model fields
  `window_short`/`window_long` are renamed `tau_short`/`tau_long`.

## msPAF engine vectorisation + `analyte_pafs` API change (issue #30)

* `compute_mspaf_per_sample()` is vectorised across all (sample × draw) blocks
  in one pass (normalisation/ARA with a wide co-analyte join and per-analyte
  formula evaluation; CA/IA via grouped reductions) instead of a per-block
  dplyr/purrr loop. Output is mathematically identical (verified by a before/after
  equivalence harness across `add_mspaf()`, `time_weighted_aggregate()`, and
  `mspaf_daily()`, point and draws modes). Removes the per-block overhead that
  dominated draws-mode runtime; the residual cost is now `ssdtools::ssd_hp()`.
* **API change (pre-v1):** the per-analyte PAF breakdown is no longer a per-row
  `analyte_pafs` **list-column** — it is a flat tibble attribute, retrieved with
  the new exported `analyte_pafs()` accessor (mirroring `ara_summary()`), keyed
  by `(site_id, sample_id, draw_id, analyte)`. Affects `add_mspaf()` and
  `mspaf_daily()` output schemas.

## Daily msPAF uncertainty via a state-space (Kalman) residual smoother (issue #16)

`mspaf_daily()` propagates uncertainty through the daily msPAF time series via
`ndraws`, `seed`, `return`, `interval`, `central`, `grab_cv`, `ou_scale`, and
`kappa`.  The latent residual impact state is modelled as a continuous-time
AR(1)/Ornstein–Uhlenbeck process fitted by a **Kalman filter + smoother** (KFAS):
the smoother posterior mean is the centre line and its posterior covariance gives
coherent draws, so the centre line and the credible band come from **one** model.
This replaces the earlier deterministic interpolation + bolted-on OU bridge.

Uncertainty sources:

* **S1–S3** GAM posteriors (Vp matrix) — reference, impact, WQ-layer fits
* **S4/S5** State-space residual smoother for the impact residual *S* / WQ
  residual *d*: posterior variance pinches at grabs, balloons mid-gap, and (via
  `kappa`) balloons faster across gaps spanning high flow
* **S6** Grab measurement error, entering as the smoother's observation noise
* **S7** Co-analyte grab measurement error (coherent per-grab multipliers)

Key properties: θ is fitted by 1-D MLE (γ by moments) with a Brownian/degenerate
fallback ladder for sparse analytes; the daily grid is clipped per analyte to its
grab span (no extrapolation beyond the grabs); residual draws are precomputed
once per analyte (`simulateSSM`) — the single point where future cross-analyte /
inter-site coupling would correlate the draws.  Requires the **KFAS** package.

New return modes:

* `return = "summary"` (default): one row per (date, site_id) with columns
  `mspaf`, `mspaf_lower`, `mspaf_upper`.  The centre line (`mspaf`) is the
  deterministic smoother posterior mean; bounds are empirical quantiles of the
  stochastic draws.
* `return = "draws"`: one row per (date, site_id, draw_id) with `draw_id`
  running 1..N.

The `seed` argument saves and restores `.Random.seed` so that calling
`mspaf_daily()` with a seed does not alter the caller's RNG state.

A new `parallel = FALSE` argument opts in to parallel draw execution via
`future.apply::future_lapply()`, honouring whatever `future::plan()` the
caller has established.  The **future.apply** package is listed in `Suggests`.

A new vignette `vignette("daily-uncertainty")` describes the model,
parameters, and how to interpret credible envelopes.

## End-to-end uncertainty propagation (issue #2)

Added a draw-carrier contract throughout the pipeline: a long-format chemistry
frame where each `(sample_id, analyte)` cell is either one exact row
(`draw_id = NA`, broadcasts as a constant) or N rows keyed `draw_id = 1..N`.

New exports:

* `summarise_draws()` — collapses a draw-bearing frame to median + credible
  interval; identity on point frames.
* `draw_measurement_error()` — applies per-measurement lognormal draws on
  observed detected values.
* `add_mspaf(..., return = "draws")` — returns per-draw msPAF rows.
* `time_weighted_aggregate(..., return = "draws")` — propagates draws through
  the 365-day time-weighted chronic aggregate.
* `impute_coanalytes(..., return = "draws")` — full GAM posterior-predictive
  draws for co-analyte imputation.
