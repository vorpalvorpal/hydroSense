# leachatetools (development)

## Daily AmsPAF with uncertainty quantification (issue #16)

`amspaf_daily()` now accepts `ndraws`, `seed`, `return`, `interval`, `central`,
`grab_cv`, and `ou_scale` arguments to propagate seven sources of uncertainty
through the daily AmsPAF time series.

Seven uncertainty sources are propagated:

* **S1** Reference GAM posterior (Vp matrix, multivariate normal draws)
* **S2** Impact GAM posterior (Vp matrix)
* **S3** WQ-layer GAM posterior (Vp matrix)
* **S4** OU bridge on normalised residual state S (between-grab extrapolation)
* **S5** OU bridge on WQ score d (leading/trailing zones)
* **S6** Anchor grab measurement error (lognormal, CV = `grab_cv`)
* **S7** Co-analyte grab measurement error (coherent per-grab multipliers)

The OU bridge uses method-of-moments parameter estimation and an exact Cholesky
factorisation that converges to the Brownian bridge as θ → 0.  Cholesky factors
are precomputed once per fit (`fit_once` architecture); each of the N draws
requires only a matrix–vector product.

New return modes:

* `return = "summary"` (default): one row per (date, site_id) with columns
  `amspaf`, `amspaf_lower`, `amspaf_upper`.  The centre line (`amspaf`) is the
  deterministic heuristic-blend point prediction; bounds are empirical quantiles
  of the stochastic draws.
* `return = "draws"`: one row per (date, site_id, draw_id) with `draw_id`
  running 1..N.

The `seed` argument saves and restores `.Random.seed` so that calling
`amspaf_daily()` with a seed does not alter the caller's RNG state.

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
