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

### 1. Run the brms path in CI (Stan-in-CI)

`BRMS_SMOKE_TEST=1` still does not run in CI: rstan cannot compile a model on the
Ubuntu runner (toolchain). The brms smoke tests pass locally, so the imputation
path is currently only exercised on a developer machine. Needs a reliable
Stan-in-CI setup — most likely `cmdstanr` (precompiled), or a dedicated job —
so regressions in `impute.R` are caught automatically.

### 2. Imputation-method follow-ups

From the §10b BDL/imputation-method benchmark (`impute_method` =
`rescor_mi`/`cens`/`cens_factor`; `rescor_mi` confirmed as default):

- **`cens_factor` reparameterisation.** At 2000 iterations the shared-latent-factor
  model did not converge (R̂ 1.147, 554 divergences). Try a **non-centred**
  parameterisation of the factor (and/or tighter priors) before trusting it; it
  is currently experimental.
- **Repeated cross-validation.** The benchmark used a single random mask / one
  seed (n = 19 held-out cells, noisy). Repeat across multiple seeds / folds for a
  stable accuracy + calibration comparison. (Benchmark scripts live in the
  gitignored `test data/`.)
