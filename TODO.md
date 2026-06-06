# leachatetools — project TODO

Short list of the remaining open work. The larger v2.0 features live as GitHub
issues under the [v2.0 milestone](https://github.com/vorpalvorpal/leachatetools/milestone/1):
uncertainty propagation (#2), imputation generalisation (#3), numerical
Monte-Carlo CA (#1), and the macroinvertebrate calibration layer (#4).

The conceptual-review decisions (mixture-toxicity correctness, ammonia
normalisation, ARA reference statistic, the two SSD representations, chronic
aggregation, organics-under-BDL, the BDL/imputation-method benchmark, etc.) were
all resolved and are recorded in git history and in the code/roxygen/vignettes —
they have been removed from this file now that they are done.

---

## Open items

### 1. Imputation-method follow-ups

- **`cens_factor` rebuilt as a long-format shared-factor model — RESOLVED
  2026-06-06.** The old wide `(1 |q| sample_id)` form silently produced no
  coupling (no `cor` parameter under `subset()`) and would not converge (R̂ up to
  2.1, 554 divergences). Replaced with a single long-format univariate model
  `lv | cens(cf) ~ 0 + safe + s(PC, by = safe) + (1 | sample_id)` (per-analyte
  residual SD), where `(1 | sample_id)` is a genuine shared per-sample latent
  factor. Validated on B.S01 6-metal masked data: R̂ 1.013, 75 divergences (now
  with `adapt_delta = 0.95` default), **shared factor SD 0.54 [0.33, 0.76]** (the
  coupling is real), RMSE 0.445; draws path works. Reworked the `cens_factor`
  branches of `.fit_group_model()` + `.predict_and_merge()` (+ new
  `.predict_factor_long()`); the reshape is fully contained in the brms
  internals, `cens`/`rescor_mi` untouched. Pinned by a `VarCorr` assertion in
  `test-impute-smoke.R`.
- **Repeated cross-validation.** The §10b benchmark uses a single random mask /
  one seed (n = 19 held-out cells, noisy). Repeat across multiple seeds / folds
  for a stable accuracy + calibration comparison. (Benchmark scripts live in the
  gitignored `test data/`.)
