# leachatetools — project TODO

Living list of outstanding work, gathered from the conceptual review (2026-06)
and in-code `TODO` markers. Grouped by theme; each item notes **why** it
matters and **where** it lives. Roughly ordered most- to least-impactful within
each section.

---

## 0. Expand test coverage and documentation across the whole package (cross-cutting)

Two recent bugs (PCA `R2cum`-NULL crash and the train/predict score-scale
mismatch, §9a) lived undetected because the imputation path is only covered by
brms smoke tests that skip when Stan is unavailable. Treat that as a signal:
**every module needs brms-independent unit tests for its core numerics**, not
just end-to-end smoke tests. Audit each `R/*.R` file for:
- untested internal helpers (especially anything doing maths or reshaping);
- functions whose only coverage is a skip-prone smoke test;
- public functions lacking worked-example documentation / a vignette.
Target: a `covr::package_coverage()` baseline, then raise it deliberately.
Pair each new test with a docstring/example so coverage and docs grow together.
Run `BRMS_SMOKE_TEST=1` in CI so the brms path is actually exercised somewhere.

**Progress (2026-06-03) — substantially done.**
- Brms-independent tests added across the package: `test-paf-values`,
  `test-watertemp`, `test-to-meq`, `test-lmf` (LMF end-to-end), `test-impute-helpers`
  (`impute_coanalytes`, `.check_bdl_imputed`, `print`), plus existing
  `test-impute-pca`, `test-silo`, `test-temperature-mandatory`, `test-mspaf-ca`,
  `test-paf-as-envelope`. Suite now 242 pass / 0 fail / 2 skip.
- Coverage baseline established with covr: **~56% (≈66% with the brms path on)**;
  up from ~26%/41%. Big remaining gaps are `impute.R` (Stan-gated) and `ssd_fit.R`
  (needs the gitignored guideline XLSX + warm cache).
- Docs: all 17 exports have title/@return/@examples and every argument documented;
  13/17 examples are runnable (remaining \dontrun are Stan-only + network-only).
  Added a bundled example dataset (`leachate_demo()`, qs2 qdata format) and an LMF
  vignette; imputation + normalisation + chronic-AmsPAF vignettes also present.
- Writing the LMF test surfaced and fixed three latent namespace bugs in `lmf.R`
  (unqualified dplyr verbs, unqualified tidyr pivots, unregistered units in
  `to_meq`); R CMD check now passes with 0 notes/warnings.
- **Still open:** running `BRMS_SMOKE_TEST=1` in CI. Attempted on the Ubuntu
  runner but reverted — rstan cannot compile a model there (toolchain). The brms
  smoke tests pass locally; a separate, reliable Stan-in-CI setup (e.g. cmdstanr,
  or a dedicated job) is needed to exercise the brms path in CI.

**Also completed this iteration (infrastructure, not originally listed here):**
pkgdown website + GitHub Pages (now at vorpalvorpal.github.io/leachatetools, after
removing a stale account-level CNAME), R-CMD-check + pkgdown CI workflows, README,
LICENSE, `.Rbuildignore`; migrated archived `qs` → `qs2`; declared previously
undeclared dependencies; fixed an unknown-analyte crash in `ssd_paf()`; unified the
codebase on qualified `pkg::` calls.

---

## 1. Uncertainty propagation (cross-cutting — highest priority)

**Goal:** carry uncertainty end-to-end — through imputation, normalisation, SSD
lookup, CA/RA mixture combination, and temporal aggregation — and expose it at
the package boundary so downstream consumers (beyond this package) can keep
propagating it.

Current state: uncertainty exists in pieces but is collapsed too early.
- `impute_chemistry(return = "draws")` already emits posterior draws, but the
  chronic/AmsPAF pipeline consumes point estimates (`return = "point"`).
- `ssd_paf()` can return CIs (`ci = TRUE`) but the AmsPAF path uses point PAFs.
- The analytic CA combination (`compute_ca_group_mspaf`) uses a single sigma per
  group with no uncertainty on HC50/sigma themselves.

Design questions to settle before implementing:
- **Representation:** posterior draws (Monte Carlo, composable, heavy) vs.
  parametric moments (mean + variance, light, needs a propagation rule at each
  step). Draws are the honest default given imputation is already Bayesian.
- **Plumbing:** thread a `draw_id` dimension through `time_weighted_aggregate()`
  and `add_amspaf()` so a single posterior draw flows coherently from chemistry
  → PAF → msPAF → time-average, then summarise at the very end.
- **SSD parameter uncertainty:** draw HC50/HC05 (hence sigma) from the fitted
  `ssdtools` model rather than treating them as fixed.
- **Boundary output:** return per-endpoint posterior draws (or quantiles) so an
  external risk model can ingest them.

This is a substantial refactor — do it deliberately, likely after the
imputation generalisation (§6) lands so the draws interface is stable.

---

## 2. Mixture-toxicity method correctness (ecotox/stats)

- **CA mixture-slope (sigma) formula — RESOLVED & FIXED (2026-06).** Confirmed
  against the primary source (De Zwart & Posthuma 2005, eq. 6,
  `references/Zwart.Posthuma2005.pdf` p.2672). The canonical CA combination is
  `msPAF_CA = Φ(log10(ΣTU) / σ̄)` with `σ̄ = mean(σ)` — the **plain arithmetic
  average of the component SSD slopes** (their `β̄_TMoA`; `β = √3·σ/π`), NOT a
  TU-weighted mean and NOT the old `sqrt(sum(w²σ²))` variance form. HUs/TUs are
  summed (non-log) then the combined PAF is read off one log-normal CDF.
  `compute_ca_group_mspaf()` now uses `sigma_mix <- mean(sigma)`; pinned by
  `tests/testthat/test-mspaf-ca.R`.

- **TMoA / mode-of-action grouping — RESOLVED (hybrid, 2026-06-02).** The
  bundled `moa_group` column is now populated from each substance's stated
  mechanism of action (per its ANZG/Warne guidance doc), mapped onto the De
  Zwart & Posthuma (2005, Table 3/4) TMoA vocabulary:
  - **Metals & inorganic ions — grouped by stated mechanism** (the primary
    source gives only a principle — CA within a shared *primary target
    receptor*, RA across — plus a single two-metal example, NOT a per-metal
    rule; the ANZG docs don't use TMoA language but their mechanistic text
    points to shared modes). Gill-ionoregulatory metals Al/Cd/Cu/Ni/Pb/Zn →
    one CA group `ionoregulatory (gill)`; Hg/Se → `sulfhydryl binding`; As →
    `arsenic` (speciation-specific); Cr/B/Mn each solo.
    Evidence basis: Al/Cd/Pb/Mn/Hg/Se/As from the older Warne HTML factsheets;
    Cu/Ni/Zn/Cr/B from the NEWER ANZG PDF tech briefs (§2.1 Mechanism). Those
    briefs confirm Cu (Na+ regulation, gill primary target), Ni (ionoregulatory
    + respiratory + ROS, gill-dependent), Zn (Ca2+ uptake disruption at gill) →
    ionoregulatory group; and explicitly state Cr(III) has "no evidence of a
    specific mode of action" and B acts by passive membrane diffusion with
    species-specific effects → both correctly solo. Ni is the least-pure member
    (respiratory/ROS alongside ionoregulatory) but ionoregulatory is primary.
    `ammonia` and `nitrate` each their own group. Ammonia/nitrate checked
    against the NEW ANZG tech briefs (PDF, not Warne HTML): ammonia §2.1 =
    gill structural damage + uncoupling oxidative phosphorylation + CNS ATP
    depletion + osmoregulatory disruption (NOT metal-style gill ion-uptake
    competition); nitrate §2.1 = methaemoglobinaemia (O2-pigment conversion)
    + salinity-driven osmoregulatory effects, gills low-permeability to NO3.
    Both confirmed distinct from the ionoregulatory-metal mechanism → solo.
  - **Organics CA-grouped by shared mechanism:** `nonpolar narcosis` (PAHs,
    BTEX, HCB), `polar narcosis` (phenol, 2-chlorophenol), `uncoupler
    oxidative phosphorylation` (2,4,6-TCP, pentachlorophenol), `AChE
    inhibition: organophosphate` (parathion, azinphos-methyl, dimethoate,
    demeton-S-methyl), `neurotoxicant: cyclodiene` (aldrin, dieldrin, endrin,
    g-BHC/lindane), `neurotoxicant: DDT` (DDT, DDE, methoxychlor).
  Users can override the whole scheme by supplying their own metadata CSV.
  Doc comment in `R/mspaf.R` updated; classification lives in
  `inst/extdata/anzecc_analyte_metadata.csv`. Hg+Se are CA-combined as a
  shared sulfhydryl-binding TMoA (both docs cite sulfhydryl/enzyme disruption).
  Possible future refinement: Cr/B/Mn currently stay solo because their ANZG
  docs state no clear primary mode — revisit if better mechanism data emerge
  (e.g. whether Cr(III) belongs with the ionoregulatory group).

- **As speciation default — RESOLVED (2026-06-02).** `.SSD_NAME_MAP` keeps
  unspeciated `"As" → "As_V"`, but the rationale is now corrected and verified:
  this is the conservative **max-PAF envelope** of the two speciation SSDs, NOT
  a redox/oxygenation assumption. Empirically (checked 0.1–5000 µg/L) the As(V)
  SSD yields a higher PAF than the As(III) SSD at every realistic concentration
  — they cross only above ~10 mg/L dissolved As. So even though landfill
  leachate is reducing (As(III) dominant), As(III) oxidises only slowly in
  freshwater, and As(III) is more toxic to the single most-sensitive diatom,
  the whole-assemblage As(V) SSD sits lower (lower tail pulled down by corrected
  algal acute points) → As(V) is conservative regardless of speciation.
  Comment in `R/paf.R` + CSV note updated; pinned by
  `tests/testthat/test-paf-as-envelope.R` (fails loudly if a future SSD-data
  revision ever lets As(III) overtake As(V), at which point switch to a true
  runtime max-of-both envelope). Users can still pass As(III)/As(V) explicitly.

---

## 3. Ammonia (NH3-N) normalisation — trace & de-risk double correction — RESOLVED 2026-06-02

Traced both paths and **found + fixed an inverted normalisation formula**.

- **SSD basis confirmed:** the metadata declares NH3-N as total ammonia-N at the
  ANZG index condition pH 7.0 / 20 °C (`dgv_conditions`), toxicity driven by the
  un-ionised NH3 fraction `f = 1/(1+10^(pKa−pH))`, `pKa = 0.09018 + 2729.92/T(K)`.
- **Direction bug (fixed):** the metadata formula was `C * f_ref / f_sample`,
  which normalises a high-pH sample *downward* — backwards. A 900 µg/L sample at
  pH 8.5 holds 100 µg/L un-ionised NH3 (28× the DGV-equivalent ~3.5 µg/L) yet the
  old formula scaled it to 31.9 µg/L (appears harmless). Corrected to
  `C * f_sample / f_ref`; the same sample now normalises to ~25,400 µg/L. This
  understated ammonia risk at high pH and overstated it at low pH — substantial,
  so the corrected direction is pinned by tests in `test-normalise.R`.
- **Two-path asymmetry documented, not unified:** `ssd_paf()` is the manual path
  (caller pre-corrects); `add_amspaf()` auto-applies the metadata formula. Added
  the exported helper **`correct_ammonia_ph_temp(conc, pH, temperature_C)`** so
  `ssd_paf()` callers can do the external step with the *identical* maths
  (`test-normalise.R` asserts helper == metadata formula at several pH/T points).
  Roxygen on the helper + the `R/paf.R` header now warn loudly NOT to pre-correct
  before `add_amspaf()` (double-application).
- **Incidental fix:** repaired CSV encoding corruption (`<c2><b5>`→µ ×8,
  `<e2><80><94>`→— ×1) and a malformed As row (unquoted commas in `notes` had
  spilled it to 25 columns, corrupting As `coanalytes_required`/
  `normalisation_formula`/`moa_group`). All 43 rows now parse to 23 columns.

**Follow-up (Chesterton's fence + design, 2026-06-02):** searched prior sessions
and found ammonia was singled out as "external" not for toxicological reasons but
because its correction uniquely needs **temperature**, which "may not be measured
consistently" — a `temperature_default_C` fallback was proposed but never built.
Resolution adopted: make temperature **mandatory** instead of optional.
- `add_amspaf(require_temperature = TRUE)` (default): hard error if any sample
  reports NH3-N but has no water `temperature` row (`.assert_temperature_present`).
  No silent ammonia drop. `require_temperature = FALSE` for ammonia-free datasets.
- New exported **`get_silo_air_temp(lat, lon, start, end)`** — wraps
  `weatherOz::get_data_drill()` (SILO Data Drill, CC-BY 4.0), returns daily mean
  air temp `(Tmax+Tmin)/2` in the shape `estimate_water_temp()` wants; disk-cached
  under `R_user_dir(..,"cache")/silo`. `weatherOz` added to Suggests.
- Air→water imputation kept as the existing per-site `lm` in
  `estimate_water_temp()` with its hard ≥5-pairs requirement (decision: no pooled
  fallback — a site with no water-temp pairs cannot proceed).
- The `correct_ammonia_ph_temp()` helper is now documented as a standalone
  `ssd_paf()` spot-check convenience only; `paf.R` header rewritten to state
  `ssd_paf()` applies no chemistry normalisation for any analyte.
- Tests: `test-silo.R` (mocked weatherOz: mean transform, caching, validation),
  `test-temperature-mandatory.R` (the gate). Full suite 164 pass / 0 fail.

**Update (2026-06-03):** `estimate_water_temp()` now optionally adds a
first-harmonic day-of-year term to the air→water regression and selects between
the air-only and air+season models by AICc (args `seasonal`, `seasonal_min_n`,
`seasonal_min_quarters`); covered by `test-watertemp`. This captures seasonal
hysteresis without the full Bayesian model.

Still genuinely optional: upgrade `estimate_water_temp()` to a hierarchical
Bayesian air→water model (partial pooling + predictive distributions) when the
uncertainty-propagation work lands — would also dissolve the no-pairs hard stop.

---

## 4. Prescreen — detection-frequency escape hatch — RESOLVED 2026-06-04

`prescreen_analytes()` drops analytes below a detection-frequency threshold.
Rare-but-potent toxicants (e.g. a pesticide detected in 2% of samples but at
ecotoxicologically significant concentrations) could be screened out.

**Implemented:** a potency-based escape hatch (`potency_keep = TRUE` by
default). An analyte that fails the frequency screen is still kept if any
detected concentration reaches `potency_frac` (default 1) times its
95%-species-protection guideline value (`dgv_95pct_ug_L` in the metadata).
Needs a numeric `value` column in µg/L (matching the DGV units); only analytes
carrying a DGV (toxicants) are eligible, so major ions are unaffected. Reported
via a new `potency_kept` column / cli message; pinned by `test-prescreen.R`.
Analytes with no DGV (a few no-SSD organics) still rely on manual `protect`.

---

## 5. Organics imputation under near-total BDL — RESOLVED (skip + document) 2026-06-05

**Decision: do not impute near-all-BDL organics; document the limitation
(option c).** In WMF leachate organics are almost entirely below detection,
with too few detected points to fit a defensible model and no reliable covariate
(DOC fixes bulk carbon, not a congener-specific ratio). Rather than manufacture
values we let the existing machinery exclude them and say so plainly:
- The **DOC-like hurdle** already skips the organics model unless ≥1 of
  {DOC,TOC,BOD,COD,cBOD} is present.
- `min_target_detect_freq` (default 0.05, §10) drops any organic detected in
  fewer than 5% of samples as an imputation target — which on real leachate
  panels removes essentially all of them.
Net effect = option (c): organics carrying no recoverable signal are left as
reported (non-)detections, not imputed. Documented in the imputation vignette
("When not to impute"). The DOC-scaled-ratio and censored-upper-bound ideas are
kept here as possible future enhancements if a dataset ever has enough detected
organics to justify them, but they are out of scope now.

---

## 6. Added Risk Approach (ARA) — reference statistic — RESOLVED (geom_mean) 2026-06-05

**Decision: keep the geometric mean as the default `ref_norm`; expose
median/arith_mean/p80/p90/p95 for frameworks that mandate otherwise.**
ARA computes `C_adj = max(C_norm − ref_norm, 0)` — the increment above the
background the local pool is already adapted to. Crommentuijn et al. (1997, RIVM
601501001) frame the background `Cb` as the *representative natural background*
(MPC = MPA + Cb): a single typical/central value, NOT an upper percentile, and
they do not prescribe a specific summary statistic — so there is no external
"standard" to defer to, only the framework's intent. Geometric mean is the right
default because (1) it matches that intent — a central tendency, so we don't
over-subtract and erode protection; (2) for log-normal aquatic chemistry it is
the maximum-likelihood central tendency and uses every reference observation
(robust on small reference sets); (3) it is consistent with the PICT reading —
the community integrates against *typical* exposure, not the upper tail. A high
percentile would assume adaptation to extremes, subtract more, and classify less
as impact — the less-protective, harder-to-defend choice. Already documented
with this rationale in the `chronic-amspaf-interpretation` vignette ("Reference
summary statistic"); `prepare_reference(summary=)` exposes the alternatives.
Lives in `R/reference.R` / `R/mspaf.R` (`.amspaf_adjust`).

---

## 7. Two SSD representations — document coexistence

The package uses an SSD in two distinct ways that currently coexist without a
clear narrative:
- **Analytic log-normal** (HC50 + sigma) for the CA mixture combination in
  `compute_ca_group_mspaf()`.
- **Model-averaged "multi" curve** (6 distributions via `ssdtools::ssd_hp`) for
  per-analyte diagnostic PAF in `ssd_paf()`.
**Why this is intrinsic, not accidental (clarified 2026-06-05):** summing toxic
units and reading the combined PAF off a single log-normal CDF
(`msPAF_CA = Φ(log10(ΣTU)/σ̄)`) *is* the standard msPAF-CA method (De Zwart &
Posthuma 2005, §2). CA treats similarly-acting chemicals as dilutions of one
"super-chemical", so it is intrinsically parametric — it needs each component as
HC50 + a slope. That is exactly why two representations exist: the diagnostic PAF
can use the rich model-averaged curve, but analytic CA *must* commit to a
parametric slope. Two consequences worth recording:
- **IA is shape-agnostic; CA is not.** Independent action across MoA groups is
  `msPAF_IA = 1 − ∏(1 − PAF_i)` and consumes any per-component PAF — so the
  across-group step could read off the model-averaged curve even though the
  within-group CA step cannot. Check we are consistent about which PAF feeds IA.
- **Non-log-normal SSDs.** The log-logistic has the same closed form (the
  framework converts via `β = √3·σ/π`), so other common-slope shapes are fine.
  For an arbitrary/model-averaged shape there is no closed form; the general
  alternative is **numerical CA** — Monte-Carlo the mixture by sampling species
  sensitivities from each component SSD (under a cross-chemical concordance
  assumption), sum TUs per virtual species, take the empirical affected fraction.
  Handles any shape but is heavier and needs the concordance assumption — a
  natural fit for the §1 uncertainty refactor (CA over draws).

**Decision (2026-06-05):**
- **Now (v1.x):** keep the analytic log-normal CA and **document the assumption
  plainly** — that the CA step represents each component as a log-normal
  (HC50 + mean slope) and that this is the standard De Zwart & Posthuma method,
  exact for log-normal/log-logistic SSDs and an approximation of the
  model-averaged curve otherwise. Done in the `mspaf.R` roxygen + a "How the
  mixture is combined" section of the chronic-AmsPAF vignette.
- **v2.0 (tracked: GitHub issue #1):** add numerical (Monte-Carlo) CA under full
  rank-concordance as a selectable `mspaf_method`, validate it against the closed
  form, and make it the default; keep `"lognormal"` as the fast option. Scope ≈ 1
  focused day: a
  ~40-line internal that, per CA group, samples N virtual species via shared
  inverse-CDF draws (one Uniform per species, reused across chemicals = full
  concordance), `c_{chem,s} = Q_chem(u_s)`, `TU_s = Σ C/c`, `msPAF = mean(TU_s ≥
  1)`; thread an `mspaf_method` arg through `compute_ca_group_mspaf()` /
  `add_amspaf()`; pull each analyte's quantile function from the fitted SSD.
  **Performance:** precompute each analyte's quantile function on a fixed
  `u`-grid *once* (the SSD doesn't change per sample), so each sample's CA is
  vectorised arithmetic over the grid — penalty is a small constant per sample
  (the one-time grid build, root-finding on the model-averaged mixture CDF, is
  the only real cost and is amortised). Validate with a test asserting MC ≈
  closed form on a log-normal analyte before flipping the default. Natural CA
  operator for the §1 uncertainty work (CA over posterior draws). Also then:
  quantify the log-normal-vs-model-averaged disagreement on real analytes.
  Lives in `R/mspaf.R` + `R/paf.R`.

---

## 8. Chronic chemistry aggregation — RESOLVED (Path B + standalone + single kernel) 2026-06-05

`time_weighted_aggregate()` (`R/chronic.R`) is value-agnostic (geomean for
chemistry, arithmetic mean for AmsPAF). Decisions:

- **Pipeline position = Path B (per-sample AmsPAF, then time-average the PAFs).**
  Argument: the community integrates the *effect* over time, not the
  concentration, and the SSD is non-linear — so averaging concentration first
  (Path A) then mapping through the concave SSD upper tail systematically
  *under-states* the time-averaged risk for pulsed exposures (Jensen's
  inequality; the storm-pulse worked example in the chronic vignette shows ~10×).
  Path A's only merits are fewer SSD evaluations and yielding a chronic
  *chemistry* profile as a by-product — neither outweighs the bias. Documented
  in `chronic-amspaf-interpretation` vignette (Path A vs B + Jensen table).
- **Standalone chemistry use = first-class.** The function is value-agnostic and
  independently useful for producing a chronic-averaged chemistry profile
  (`summary = "geom_mean"`) outside the AmsPAF path; the roxygen already frames
  this as a primary use, not a side effect.
- **Single weighting scheme, justified, kept.** Forward-step duration weighting
  (minimal correction for pulse-biased sampling) + exponential decay (standard
  memory kernel, one interpretable `tau_days`). Richer kernels would add
  parameters with no defensible way to fit them from routine data, so they are
  deliberately omitted — tune `tau_days`, don't swap the kernel. Rationale now
  stated in the roxygen.

---

## 9. Imputation generalisation refactor (after §1 design settles)

Generalise the imputation engine so it is domain-agnostic (a reusable Bayesian
multivariate-GAM + PCA imputation tool) with a leachate-specific preset layered
on top. End state: three clear user-facing parts to the package —
**LMF**, **imputation**, **AmsPAF**. Lives in `R/impute.R`. Coordinate with §1
so the draws/uncertainty interface is designed once.

---

## 9a. BUG — PCA train/predict score-scale mismatch (imputation) — FIXED (2026-06)

**Status: FIXED.** `.prepare_chem_pca()` previously returned `pc_scores` copied
from `nipals::nipals()$scores`, but `.compute_pca_scores()` (used at prediction)
produces a regression projection (`X_scaled %*% loadings` via `.nipals_score_row`).
These differ by a per-component eigenvalue factor
(`nipals$scores[,k] == projection[,k] / eig[k]`), so the brms smooths were
trained on one scale and predicted on another → silently wrong imputations.
Fix: `.prepare_chem_pca()` now builds `pc_scores` by calling
`.compute_pca_scores()` on the training data, guaranteeing train==predict by
construction. Pinned by `tests/testthat/test-impute-pca.R` (asserts equality and
that the scores are NOT the eigenvalue-shrunk nipals scores).

`impute_coanalytes()` was never affected — it both fits and predicts its GAM
with `.compute_pca_scores()`.

Also fixed in passing: `cum_var <- pca_fit$R2cum` was `NULL` under nipals 1.0
(per-component `$R2`, no `$R2cum`); now `cumsum(pca_fit$R2)` with an `R2cum`
fallback. Also pinned by `test-impute-pca.R`.

---

## 10. In-code TODOs (smaller, well-specified)

- **Probabilistic NO3-N hardness weighting — RESOLVED 2026-06-04** (`R/paf.R`).
  `ssd_paf("NO3-N", hardness_mg_L=)` now blends the soft/moderate/hard SSDs by
  hardness-class probability instead of snapping to one class:
  `p_soft = plnorm(30, log(h), s)`, `p_hard = 1 - plnorm(150, log(h), s)`,
  `p_mod = 1 - p_soft - p_hard`, `s = sqrt(log(1+cv^2))`, default `cv = 0.05`
  (`.no3_weights()`); `ExpectedPAF = Σ w·PAF_class`. Continuous across the
  class boundaries (at h=150 PAF is the mean of the mod/hard PAFs); pinned by
  `test-paf-values.R`. `ssd_hc50()` / the CA mixture step still use a single
  representative class (`.no3_class()`); revisit if CA needs the blend too.

- **Benchmark the three BDL/imputation configs — RESOLVED 2026-06-05.** The
  three configs are a `fit_imputation_model(impute_method=)` option
  (`"rescor_mi"` default, `"cens"`, `"cens_factor"`):
  1. `rescor = TRUE` + `mi()` + post-hoc DL cap (current = `"rescor_mi"`).
  2. `rescor = FALSE` + `cens("left")` (clean BDL; no coupling = `"cens"`).
  3. `rescor = FALSE` + `cens("left")` + shared `(1 | sample_id)` factor
     (= `"cens_factor"`).
  Definitive hold-out benchmark on B.S01 routine metals (As/Cr/Cu/Ni/Pb/Zn, 89
  samples, 19 masked detected cells censored at 2x truth, log10 scale; 2000
  iter / 1000 warmup, 3 chains **serial**, ndraws=1000, batch=15):
  | method | t_fit | t_impute | max R̂ | divergences | RMSE | MAE | bias | cov90 | cov50 | width90 |
  |--------|-------|----------|-------|-------------|------|-----|------|-------|-------|---------|
  | rescor_mi | 653 s | 6 s | 1.044 | 215 | **0.247** | **0.201** | **−0.063** | 1.00 | 0.53 | 1.31 |
  | cens | 186 s | 22 s | **1.007** | **14** | 0.533 | 0.379 | −0.366 | 0.89 | 0.53 | 1.60 |
  | cens_factor | 433 s | 30 s | 1.147 | 554 | 0.476 | 0.344 | −0.343 | 0.95 | 0.58 | 1.34 |
  **Verdict — keep `rescor_mi` as default.** Best accuracy by far and near-zero
  bias (−0.06 vs ≈−0.35): masking sets DL = 2×truth, so a model with no
  cross-metal signal lands at bias ≈ −0.30 (≈ DL/2); rescor_mi's near-zero bias
  shows the residual-correlation coupling genuinely recovers the level from the
  observed correlated metals. **Memory now holds** — `ndraws`/`batch_size` kept
  the draws-based imputation under 16 GB at this size (the earlier OOM is fixed).
  rescor_mi over-covers (cov90 1.00; cov50 well-calibrated at 0.53) — the safe
  direction. `cens` is the robust cheap fallback: cleanest convergence (R̂ 1.007,
  14 divergences), best-calibrated cov90 0.89, but biased low (no borrowing →
  returns ≈ DL/2). **`cens_factor` is experimental — NOT converged here (R̂ 1.147,
  554 divergences):** the shared latent factor adds posterior geometry the sampler
  can't handle at 2000 iter. Needs a **non-centred reparameterisation** (or many
  more iterations) before it is trustworthy; until then it does not deliver its
  "censoring + coupling for cheap" promise. Docs already describe the trade-off
  (imputation vignette). Caveats: single random mask / one seed (n=19, noisy);
  repeated CV across seeds still outstanding (folded into §1 / §10 follow-ups).
  (Scripts in gitignored `test data/`.)

- **Target detection-frequency filter — RESOLVED 2026-06-04.** On a dataset
  with ~100 mostly all-BDL organics (and several near-undetected metals),
  `fit_imputation_model` built a huge organics/metals model (107-way -> Stan
  parser failure; 18-way metals -> intractable). The existing `min_detect_freq`
  only filtered PCA *inputs*, not the imputation *targets*. Added
  `min_target_detect_freq` (default 0.05): metal/organic targets detected in
  fewer than that fraction of samples are dropped (with an info message). This
  is what keeps the model tractable on real leachate panels.

- **Update docs for the new imputation options — DONE 2026-06-05.** The
  `imputation` vignette now has a "Choosing the imputation method" section (the
  three methods + trade-offs), an updated detection-limit-cap note (applies to
  all methods), `min_target_detect_freq` and `ndraws`/`batch_size` mentions, and
  the method shown in the model print. The `chronic-amspaf-interpretation`
  vignette's BDL note was refreshed. README imputation text checked (still
  accurate). `fit_imputation_model` / `impute_chemistry` roxygen already updated.

- **BDL cap inspection** (`.check_bdl_imputed`, `R/impute.R`). The cap clips
  imputed BDL cells to the detection limit; for sites whose chemistry
  legitimately implies high concentrations the cap can fire often. Surface a
  summary of cap activations so the regime is auditable.
