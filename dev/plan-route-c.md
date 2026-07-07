# Route C — low-rank censored factor model for chemistry imputation

Replaces the brms `rescor_mi` / `cens` / `cens_factor` fitting-and-prediction
path for the metals/organics groups with a **left-censored latent factor
model**, resolving review findings 1–3 at the source:

- **Finding 1** (weakest default + post-hoc cap): censoring handled in the
  likelihood; no cap.
- **Finding 2** (per-draw DL cap corrupts the posterior): BDL cells are bounded
  latent parameters; prediction draws are truncated at the DL by construction.
- **Finding 3** (no prediction-time cross-analyte conditioning): the latent
  factor is inferred from a new sample's *observed* metals, so measured metals
  genuinely inform the missing/BDL ones — for in-sample **and** new samples.

> **Executable contract.** This plan is paired with a BDD suite,
> `tests/testthat/test-route-c.R` (+ `helper-route-c.R`). Every numbered
> acceptance criterion below maps to an `it(...)` spec there. Implement
> red→green: the specs are written test-first and skipped via
> `.skip_route_c()` until the code exists. **Do the prediction specs first**
> (they are Stan-free — see §Prediction) — they encode findings 1–3 and can be
> driven green before any Stan is written, using a frozen factor model.

## Status of the prerequisite (DONE)

The `claude/impute-review-hardening` pass has **landed** (PR #69). Consequences
this plan depends on:

- **PC-score inputs are hardened.** `.compute_pca_scores()` /
  `.prepare_chem_pca()` now share `.pivot_chem_wide()`, which: treats a BDL
  *concentration predictor* cell as **missing (NA)** — not full-DL, not DL/2
  (finding 5); collapses duplicate `(sample, analyte)` cells by geometric mean
  (finding 8a); and coerces all-NA collapses back to `NA` (bug B3).
- **Scale-aware log floors exist.** `.scale_aware_log_floor(x, eps)` (half the
  smallest observed positive value). Per-column floors are stored on the PCA
  object as `pca_obj$log_floors`; per-analyte floors are stored on each fitted
  group as `group$log_floors`. **Reuse these — do not reintroduce `1e-9`.**
- **Safe-name guard exists.** `.assert_safe_analyte_names(targets)` errors on a
  `make.names()` collision. Call it in the `factor` fit branch too.
- **Hurdle fixed.** `leachate_impute_groups()` drops Fe/Mn from the metals
  hurdle (finding 7). No action needed here.
- **Predictor-side censoring is explicitly OUT of scope.** BDL predictors are
  NA (unbiased, substitution-free). Quantitatively recovering the "it's low"
  signal via a *censored predictor* pre-stage is a **separate, gated follow-up:
  issue #68** — do not fold it into this work.

## Why a low-rank factor covariance (vs dense `rescor`)

With ~20 metals a full residual covariance is ~190 correlations; on ragged
leachate panels (median co-observation ≈ 6 metals) most are unidentified — the
cause of `rescor_mi`'s R̂ ≈ 1.6 geometry. A rank-`k` factor structure
`Σ = ΛΛ' + Ψ` (k = 2–3) has ~75 parameters, every sample with ≥2 metals informs
the shared loadings, the geometry is benign, and the conditional prediction
inverts a `k×k` matrix. It is the natural home for the `#48` factor-model
follow-up. `cens_factor` is the rank-1 special case that (a) can't separate
metal families onto different axes and (b) draws the factor from the population
for new samples instead of inferring it.

## Model

For sample `i`, analyte `j`, on the **natural-log** concentration scale
(matches the existing target scale: `.fit_group_model()` uses `lv = log(...)`,
**not** `log10`):

```
y_ij = mu_j(X_i) + (Λ f_i)_j + eps_ij
f_i    ~ N(0, I_k)          # latent per-sample factor scores (k = 2..3)
eps_ij ~ N(0, psi_j)        # idiosyncratic per-analyte noise
```

- `mu_j(X_i)` = spline mean on the (hardened) PC scores — the existing predictor
  structure, unchanged.
- `Λ` (J×k) loadings; `Ψ = diag(psi_j)`; induced `Σ = ΛΛ' + Ψ`.
- **Detected** cell: Gaussian likelihood on `y_ij`.
- **BDL** cell: a latent parameter with `upper = log(DL_ij)` (data
  augmentation) — proper left-censoring, jointly with the factor and Σ.

### Two-stage implementation (do this; not the one-stage joint fit)
Decouple the well-trodden spline mean from the novel censored-covariance part:

1. **Stage 1 (mean).** Per target analyte, fit
   `mgcv::gam(lv ~ s(PC1)+…+s(PCk))` on **detected** observations, where
   `lv = log(pmax(value, group$log_floors[analyte]))`. Store each GAM. Compute:
   - detected residual `r_ij = lv_ij - mu_j(X_i)`;
   - for a BDL cell, the residual censoring bound
     `b_ij = log(DL_ij) - mu_j(X_i)` (the residual is left-censored at `b_ij`).
2. **Stage 2 (censored factor model).** Fit `r_i = Λ f_i + eps_i` in Stan with
   BDL residuals as upper-bounded latents. Small, fast, benign geometry.

Known limitation (accept it): Stage-1 mean uncertainty is not propagated into
Stage 2 (minor where detections exist). A one-stage joint fit is a later
refinement, not part of this task.

### Prediction (closed form, per posterior draw, in R — NO Stan, NO refit)
**Implement and test this first.** For a new sample with detected target set
`O` (residuals `r_O`) and to-impute set `M`, using observed-analyte loadings
`Λ_O` (rows of Λ for O) and `Ψ_O` (their idiosyncratic variances), per posterior
draw `d`:

```
V = (I_k + Λ_O' Ψ_O⁻¹ Λ_O)⁻¹            # k×k factor posterior covariance
m = V Λ_O' Ψ_O⁻¹ r_O                      # k factor posterior mean
y_M | y_O  ~  N( mu_M + Λ_M m ,  Λ_M V Λ_M' + Ψ_M )
```

Sample `y_M`; for any BDL **target** cell, draw from that Normal **truncated to
`(-∞, log(DL)]`** (`truncnorm::rtruncnorm`). This does findings 1–3 in one
operation: missing metals conditioned on observed ones, BDL metals conditioned
**and** bounded, full per-draw uncertainty. The `k×k` inverse is stable even
from a single observed metal (`|O| = 1`).

Edge cases the spec pins down:
- `|O| = 0` (no target analyte observed for the sample): `V = I_k`, `m = 0`, so
  `y_M ~ N(mu_M, Λ_M Λ_M' + Ψ_M)` — the marginal predictive. Must not error.
- Rotation invariance: replacing `Λ` with `Λ R` (`R` orthogonal) leaves all
  predictions unchanged (only `ΛΛ'` enters). The spec asserts this.

### Identifiability
Fix `Λ`'s top `k×k` block lower-triangular with positive diagonal (standard
factor-analysis constraint), or fit free and rotate post hoc — rotation is
irrelevant to prediction. Choose `k` by held-out coverage / PSIS-LOO: start at
`k = 2`, test whether `k = 3` earns its keep.

## Concrete implementation targets

New/edited symbols (all in `R/impute.R` unless noted):

- `impute_method` gains the value `"factor"` in `fit_imputation_model()`'s
  `match.arg(...)` and its roxygen `@param` enumeration. Keep the three brms
  methods during a deprecation window.
- `.fit_group_model(..., impute_method = "factor")` — new branch. Returns a list
  with, in addition to the existing fields:
  - `gams` — named list of Stage-1 `mgcv::gam` objects, keyed by original
    analyte name;
  - `Lambda` — J×k matrix (rows keyed by safe analyte name via `rownames`);
  - `Psi` — length-J named numeric (idiosyncratic variances);
  - `k` — integer;
  - `fit` — the Stage-2 `CmdStanMCMC` (or a light wrapper) for diagnostics;
  - `impute_method = "factor"`.
  Call `.assert_safe_analyte_names(target_analytes)` at entry (as the brms
  branches now do).
- `.predict_factor_conditional(group, pc_wide, df_eligible, return, ndraws,
  batch_size)` — NEW. (Takes `df_eligible` in addition to
  `.predict_factor_long()`'s args, because it must read each sample's own
  observed target values.) Returns the same `pm_long` shape (`sample_id,
  analyte, .post_mean` for `point`; `sample_id, analyte, draw_id, .post_value`
  for `draws`) so `.predict_and_merge()`'s merge/fabricate logic is untouched.
  Internally: compute `mu_j` from `group$gams` on `pc_wide`; read the observed
  target residuals for each eligible sample from `df_eligible`; apply the
  closed-form update above per posterior draw of `(Λ, Ψ)`; truncate BDL target
  draws at `log(DL)`. Emits rows only for cells that need prediction
  (BDL/missing); observed cells are left to the caller.
- **Draw-domain consistency:** every group's per-draw prediction must emit the
  same number of draws so the draw carrier stays coherent. Stan groups use the
  posterior draw count `chains*(iter-warmup)`; a degenerate single-analyte group
  stores the same figure as `post_ndraws` and replicates its plug-in `(Λ,Ψ)`
  that many times (fresh Gaussian noise per draw). Point mode collapses the
  degenerate group to a single evaluation.
- **Engine gating (not brms):** the factor method uses cmdstanr + mgcv and must
  NOT trigger `.require_brms()`. `fit_imputation_model()` and
  `impute_chemistry()` call `.require_cmdstanr()` when
  `impute_method == "factor"`, `.require_brms()` otherwise.
- `.predict_and_merge()` — add a `factor` branch that dispatches to
  `.predict_factor_conditional()` (parallel to the existing `cens_factor`
  dispatch to `.predict_factor_long()`).
- `impute_chemistry()` — the DL cap is a **no-op** for `factor` (draws already
  respect the bound). Leave `.check_bdl_imputed()` in place for the brms methods;
  short-circuit it when `model$impute_method == "factor"` and document in
  `bdl_cap_summary()` that the factor method needs no cap.
- The Stan program lives at `inst/stan/factor_censored.stan` (compiled via
  `cmdstanr::cmdstan_model()`; cache the compiled model like `.brms_backend()`
  implies).

### Stan program sketch (`inst/stan/factor_censored.stan`)
```stan
data {
  int<lower=1> N;                 // samples
  int<lower=1> J;                 // analytes
  int<lower=1> K;                 // factors
  int<lower=0> N_obs;             // detected residual cells
  int<lower=0> N_cens;            // BDL residual cells
  array[N_obs] int  obs_row;      // 1..N   for each detected cell
  array[N_obs] int  obs_col;      // 1..J
  vector[N_obs]     r_obs;        // detected residuals
  array[N_cens] int cen_row;
  array[N_cens] int cen_col;
  vector[N_cens]    b_cens;       // residual upper bounds = log(DL) - mu
}
parameters {
  matrix[N, K] f;                          // factor scores (non-centred)
  matrix[J, K] Lambda_free;                // loadings (constrain top block below)
  vector<lower=0>[J] psi;                  // idiosyncratic sd
  vector<upper=0>[N_cens] r_cens_raw;      // r_cens - b_cens  <= 0  (left-censored)
}
transformed parameters {
  matrix[J, K] Lambda = Lambda_free;
  for (a in 1:K) for (b in (a+1):K) Lambda[a, b] = 0;   // lower-triangular top block
  // (diagonal positivity: put lower=0 on the K diagonal entries, omitted for brevity)
}
model {
  to_vector(f) ~ std_normal();
  to_vector(Lambda_free) ~ std_normal();
  psi ~ student_t(3, 0, 1);
  for (n in 1:N_obs)
    r_obs[n] ~ normal(dot_product(Lambda[obs_col[n]], f[obs_row[n]]), psi[obs_col[n]]);
  for (n in 1:N_cens) {
    real rc = b_cens[n] + r_cens_raw[n];   // <= b_cens
    rc ~ normal(dot_product(Lambda[cen_col[n]], f[cen_row[n]]), psi[cen_col[n]]);
  }
}
```
(`Ψ = psi^2`. Diagonal-positivity on the top `K×K` block is the one fiddly
identifiability detail — a `positive_ordered` or `lower=0` diagonal is fine.)

## Packages
- `cmdstanr` + `posterior` — sample Stage 2, wrangle draws (matches the current
  cmdstanr default).
- `mgcv` — Stage-1 PC-score smooths (already a dependency).
- `truncnorm` — truncated BDL draws at prediction (add to `Suggests`).
- base / `Matrix` — the `k×k` factor update.
- `gllvm` was considered and **rejected**: no left-censored Gaussian family.
- `zCompositions::lrEM` was intended as the Route D benchmark but **cannot run
  on B.S01** (needs a fully-complete column; none exists). Replaced by an
  internal marginal/no-factor baseline (see Validation).

## Validation

> **Phase-D reality check.** Two things learned attempting this on B.S01:
> (1) `zCompositions::lrEM` **cannot run** on the B.S01 panel — it needs at
> least one fully-complete column and there is none, so the external Route D
> baseline as originally specified is infeasible. (2) The k=2 fit shows
> apparent non-convergence, but this is **rotation non-identifiability of
> `Lambda`, not a sampling failure**: the rotation-invariant implied covariance
> `Σ = ΛΛ' + Ψ` — the only thing prediction uses — converges cleanly (R̂ ≈ 1.0
> while a raw `Lambda` element sits ~1.14). So the convergence question must be
> asked of invariant functionals, and the external benchmark must be replaced.
> The offline validation below is `dev/bench-route-c.R`; it is NOT a unit test
> (the unit tests are smoke + invariant-sanity only).

### Convergence — measure invariant functionals, not raw `Lambda`
Use `.route_c_convergence(group)`: max R̂ over the implied `Σ` entries, over
`psi`, and `lp__`. These are rotation-invariant and are what imputation
depends on. Report R̂ / ESS / divergences / E-BFMI on **`Σ` and `psi`** (not
`Lambda`) at production panel size and budget. On real B.S01 metals (more
analytes than the 2–4 analyte unit fixtures → k identified better), decide
whether any residual `Σ` non-convergence remains after a reasonable budget; if
so that is a genuine follow-up (more data / longer runs), separate from this
task.

### Baseline — internal "no-borrowing" marginal (replaces the infeasible lrEM)
Since `lrEM` can't run, the informative comparison is the **marginal / no-factor
baseline**: the Stage-1 GAM mean + an *independent* left-censored residual per
analyte (i.e. the degenerate path applied to every analyte, k=0, no shared
factor). This always runs and directly isolates the finding-3 borrowing
benefit: factor vs marginal = "does conditioning on co-measured metals help?".
Also compare against the already-benchmarked `rescor_mi` / `cens_factor` / `cens`
numbers from the impute-method benchmark. (Optional: run `lrEM` on a reduced,
completable sub-panel purely as an external sanity point — not the main
comparison.)

### B.S01 hold-out (reuse the existing harness) — `dev/bench-route-c.R`
Reuse the 3-seed masking hold-out (memory: `calibration/impute-method-benchmark`)
that produced the `rescor_mi` RMSE 0.30 vs 0.75/0.92 numbers. Per masked cell,
for **factor vs marginal-baseline vs rescor_mi vs cens_factor**:
- **Point recovery:** RMSE / bias vs held-out truth. Target ≈ `rescor_mi`
  (~0.30) — no accuracy loss.
- **Calibration:** 90% interval coverage ≈ nominal (`rescor_mi` fails this;
  `cens_factor` half-achieves it).
- **Convergence:** `.route_c_convergence()` on `Σ`/`psi` + divergences/E-BFMI.
- **Conditioning (the finding-3 win):** mask Zn, vary the *observed* Cu/Cd/Pb in
  the same sample, confirm imputed Zn moves and **in the right direction vs the
  marginal baseline (which cannot move at all)**. Quantify effect size to settle
  `k = 2` vs `k = 3` by held-out predictive log-density / coverage.

**Win condition:** point accuracy ≈ `rescor_mi`, coverage ≈ nominal, `Σ`
convergence clean, and imputed metals demonstrably respond to co-measured metals
(beat the marginal baseline on conditioning) — collapsing the current "accurate
OR trustworthy" trade-off into one method.

## Rollout
1. ~~Land `claude/impute-review-hardening`~~ **DONE (PR #69).**
2. **Prediction first (Stan-free).** Implement `.predict_factor_conditional()`
   against a frozen factor model; drive the §Prediction BDD specs green. This
   demonstrates the finding-3 conditioning on B.S01 (≈ a day; also fixes `k`).
3. Write `inst/stan/factor_censored.stan` + the `factor` branch of
   `.fit_group_model()`; drive the fit BDD specs green (gated on `cmdstanr`).
4. Wire `impute_method = "factor"` end-to-end; drive the end-to-end BDD specs
   green (draws respect the DL; imputed metals condition on observed).
5. Benchmark on B.S01 (`dev/bench-route-c.R`): factor vs the marginal baseline
   vs `rescor_mi`/`cens_factor`, reporting RMSE, coverage, `Σ`-convergence and
   the finding-3 conditioning effect. If the win condition holds, flip the
   default and deprecate the brms methods. Fold in findings 4/8b if desired.
6. Only if the #68 sensitivity check says predictor BDL matters, revisit
   censored predictors — separately.
7. Open follow-up: real-data `Σ`-convergence at k=2 under a larger
   data/iteration budget than a unit test affords (the residual small-N
   rotation difficulty), if the B.S01 benchmark still shows any invariant-
   functional non-convergence.
