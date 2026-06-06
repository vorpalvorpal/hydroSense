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

- **`cens_factor` is structurally broken (diagnosed 2026-06-06).** The intended
  shared per-sample latent factor (`(1 |q| sample_id)` correlated across the
  metals) is **never created**: because each response uses `subset()` (a
  different subset of rows), brms cannot estimate cross-response random-effect
  correlations, so `prior_summary()` shows **no `cor` parameter**. The method
  therefore adds an independent per-(sample, analyte) random intercept with one
  observation per cell — weakly identified vs the residual (→ the 554 divergences
  / R̂ up to 2.1 in the sweep) — while delivering **none** of the advertised
  coupling. A prior/adapt_delta sweep confirmed tuning cannot fix it (divergences
  fall but R̂ worsens to 2.1). The correct fix is a **long-format restructure**:
  one univariate model `value | cens(cf) ~ 0 + analyte + s(PC, by = analyte) +
  (1 | sample_id)` where `(1 | sample_id)` is a single per-sample latent shared
  across analytes (the genuine coupling), well-identified because each sample
  contributes several analyte observations. Reworks the cens_factor branch of
  `.fit_group_model()` **and** `.predict_and_merge()`. Decision pending: rebuild
  (this restructure) vs remove cens_factor vs keep-but-hard-gate.
- **Repeated cross-validation.** The §10b benchmark used a single random mask /
  one seed (n = 19 held-out cells, noisy). Repeat across multiple seeds / folds
  for a stable accuracy + calibration comparison once the method question above
  is settled. (Benchmark scripts live in the gitignored `test data/`.)
