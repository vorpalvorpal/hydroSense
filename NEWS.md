# leachatetools (development)

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

## AmsPAF engine vectorisation + `analyte_pafs` API change (issue #30)

* `compute_amspaf_per_sample()` is vectorised across all (sample × draw) blocks
  in one pass (normalisation/ARA with a wide co-analyte join and per-analyte
  formula evaluation; CA/IA via grouped reductions) instead of a per-block
  dplyr/purrr loop. Output is mathematically identical (verified by a before/after
  equivalence harness across `add_amspaf()`, `time_weighted_aggregate()`, and
  `amspaf_daily()`, point and draws modes). Removes the per-block overhead that
  dominated draws-mode runtime; the residual cost is now `ssdtools::ssd_hp()`.
* **API change (pre-v1):** the per-analyte PAF breakdown is no longer a per-row
  `analyte_pafs` **list-column** — it is a flat tibble attribute, retrieved with
  the new exported `analyte_pafs()` accessor (mirroring `ara_summary()`), keyed
  by `(site_id, sample_id, draw_id, analyte)`. Affects `add_amspaf()` and
  `amspaf_daily()` output schemas.

## Daily AmsPAF uncertainty via a state-space (Kalman) residual smoother (issue #16)

`amspaf_daily()` propagates uncertainty through the daily AmsPAF time series via
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
  `amspaf`, `amspaf_lower`, `amspaf_upper`.  The centre line (`amspaf`) is the
  deterministic smoother posterior mean; bounds are empirical quantiles of the
  stochastic draws.
* `return = "draws"`: one row per (date, site_id, draw_id) with `draw_id`
  running 1..N.

The `seed` argument saves and restores `.Random.seed` so that calling
`amspaf_daily()` with a seed does not alter the caller's RNG state.

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
* `add_amspaf(..., return = "draws")` — returns per-draw AmsPAF rows.
* `time_weighted_aggregate(..., return = "draws")` — propagates draws through
  the 365-day time-weighted chronic aggregate.
* `impute_coanalytes(..., return = "draws")` — full GAM posterior-predictive
  draws for co-analyte imputation.
