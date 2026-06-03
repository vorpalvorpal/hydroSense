# leachatetools — project TODO

Living list of outstanding work, gathered from the conceptual review
(2026-06) and in-code `TODO` markers. Grouped by theme; each item notes
**why** it matters and **where** it lives. Roughly ordered most- to
least-impactful within each section.

------------------------------------------------------------------------

## 0. Expand test coverage and documentation across the whole package (cross-cutting)

Two recent bugs (PCA `R2cum`-NULL crash and the train/predict
score-scale mismatch, §9a) lived undetected because the imputation path
is only covered by brms smoke tests that skip when Stan is unavailable.
Treat that as a signal: **every module needs brms-independent unit tests
for its core numerics**, not just end-to-end smoke tests. Audit each
`R/*.R` file for: - untested internal helpers (especially anything doing
maths or reshaping); - functions whose only coverage is a skip-prone
smoke test; - public functions lacking worked-example documentation / a
vignette. Target: a `covr::package_coverage()` baseline, then raise it
deliberately. Pair each new test with a docstring/example so coverage
and docs grow together. Run `BRMS_SMOKE_TEST=1` in CI so the brms path
is actually exercised somewhere.

**Progress (2026-06-03) — substantially done.** - Brms-independent tests
added across the package: `test-paf-values`, `test-watertemp`,
`test-to-meq`, `test-lmf` (LMF end-to-end), `test-impute-helpers`
(`impute_coanalytes`, `.check_bdl_imputed`, `print`), plus existing
`test-impute-pca`, `test-silo`, `test-temperature-mandatory`,
`test-mspaf-ca`, `test-paf-as-envelope`. Suite now 242 pass / 0 fail / 2
skip. - Coverage baseline established with covr: **~56% (≈66% with the
brms path on)**; up from ~26%/41%. Big remaining gaps are `impute.R`
(Stan-gated) and `ssd_fit.R` (needs the gitignored guideline XLSX + warm
cache). - Docs: all 17 exports have <title/@return>/@examples and every
argument documented; 13/17 examples are runnable (remaining are
Stan-only + network-only). Added a bundled example dataset
([`leachate_demo()`](https://vorpalvorpal.github.io/leachatetools/reference/leachate_demo.md),
qs2 qdata format) and an LMF vignette; imputation + normalisation +
chronic-AmsPAF vignettes also present. - Writing the LMF test surfaced
and fixed three latent namespace bugs in `lmf.R` (unqualified dplyr
verbs, unqualified tidyr pivots, unregistered units in `to_meq`); R CMD
check now passes with 0 notes/warnings. - **Still open:** running
`BRMS_SMOKE_TEST=1` in CI. Attempted on the Ubuntu runner but reverted —
rstan cannot compile a model there (toolchain). The brms smoke tests
pass locally; a separate, reliable Stan-in-CI setup (e.g. cmdstanr, or a
dedicated job) is needed to exercise the brms path in CI.

**Also completed this iteration (infrastructure, not originally listed
here):** pkgdown website + GitHub Pages (now at
vorpalvorpal.github.io/leachatetools, after removing a stale
account-level CNAME), R-CMD-check + pkgdown CI workflows, README,
LICENSE, `.Rbuildignore`; migrated archived `qs` → `qs2`; declared
previously undeclared dependencies; fixed an unknown-analyte crash in
[`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md);
unified the codebase on qualified `pkg::` calls.

------------------------------------------------------------------------

## 1. Uncertainty propagation (cross-cutting — highest priority)

**Goal:** carry uncertainty end-to-end — through imputation,
normalisation, SSD lookup, CA/RA mixture combination, and temporal
aggregation — and expose it at the package boundary so downstream
consumers (beyond this package) can keep propagating it.

Current state: uncertainty exists in pieces but is collapsed too
early. - `impute_chemistry(return = "draws")` already emits posterior
draws, but the chronic/AmsPAF pipeline consumes point estimates
(`return = "point"`). -
[`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
can return CIs (`ci = TRUE`) but the AmsPAF path uses point PAFs. - The
analytic CA combination (`compute_ca_group_mspaf`) uses a single sigma
per group with no uncertainty on HC50/sigma themselves.

Design questions to settle before implementing: - **Representation:**
posterior draws (Monte Carlo, composable, heavy) vs. parametric moments
(mean + variance, light, needs a propagation rule at each step). Draws
are the honest default given imputation is already Bayesian. -
**Plumbing:** thread a `draw_id` dimension through
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/leachatetools/reference/time_weighted_aggregate.md)
and
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
so a single posterior draw flows coherently from chemistry → PAF → msPAF
→ time-average, then summarise at the very end. - **SSD parameter
uncertainty:** draw HC50/HC05 (hence sigma) from the fitted `ssdtools`
model rather than treating them as fixed. - **Boundary output:** return
per-endpoint posterior draws (or quantiles) so an external risk model
can ingest them.

This is a substantial refactor — do it deliberately, likely after the
imputation generalisation (§6) lands so the draws interface is stable.

------------------------------------------------------------------------

## 2. Mixture-toxicity method correctness (ecotox/stats)

- **CA mixture-slope (sigma) formula — RESOLVED & FIXED (2026-06).**
  Confirmed against the primary source (De Zwart & Posthuma 2005, eq. 6,
  `references/Zwart.Posthuma2005.pdf` p.2672). The canonical CA
  combination is `msPAF_CA = Φ(log10(ΣTU) / σ̄)` with `σ̄ = mean(σ)` — the
  **plain arithmetic average of the component SSD slopes** (their
  `β̄_TMoA`; `β = √3·σ/π`), NOT a TU-weighted mean and NOT the old
  `sqrt(sum(w²σ²))` variance form. HUs/TUs are summed (non-log) then the
  combined PAF is read off one log-normal CDF.
  [`compute_ca_group_mspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/compute_ca_group_mspaf.md)
  now uses `sigma_mix <- mean(sigma)`; pinned by
  `tests/testthat/test-mspaf-ca.R`.

- **TMoA / mode-of-action grouping — RESOLVED (hybrid, 2026-06-02).**
  The bundled `moa_group` column is now populated from each substance’s
  stated mechanism of action (per its ANZG/Warne guidance doc), mapped
  onto the De Zwart & Posthuma (2005, Table 3/4) TMoA vocabulary:

  - **Metals & inorganic ions — grouped by stated mechanism** (the
    primary source gives only a principle — CA within a shared *primary
    target receptor*, RA across — plus a single two-metal example, NOT a
    per-metal rule; the ANZG docs don’t use TMoA language but their
    mechanistic text points to shared modes). Gill-ionoregulatory metals
    Al/Cd/Cu/Ni/Pb/Zn → one CA group `ionoregulatory (gill)`; Hg/Se →
    `sulfhydryl binding`; As → `arsenic` (speciation-specific); Cr/B/Mn
    each solo. Evidence basis: Al/Cd/Pb/Mn/Hg/Se/As from the older Warne
    HTML factsheets; Cu/Ni/Zn/Cr/B from the NEWER ANZG PDF tech briefs
    (§2.1 Mechanism). Those briefs confirm Cu (Na+ regulation, gill
    primary target), Ni (ionoregulatory
    - respiratory + ROS, gill-dependent), Zn (Ca2+ uptake disruption at
      gill) → ionoregulatory group; and explicitly state Cr(III) has “no
      evidence of a specific mode of action” and B acts by passive
      membrane diffusion with species-specific effects → both correctly
      solo. Ni is the least-pure member (respiratory/ROS alongside
      ionoregulatory) but ionoregulatory is primary. `ammonia` and
      `nitrate` each their own group. Ammonia/nitrate checked against
      the NEW ANZG tech briefs (PDF, not Warne HTML): ammonia §2.1 =
      gill structural damage + uncoupling oxidative phosphorylation +
      CNS ATP depletion + osmoregulatory disruption (NOT metal-style
      gill ion-uptake competition); nitrate §2.1 = methaemoglobinaemia
      (O2-pigment conversion)
    - salinity-driven osmoregulatory effects, gills low-permeability to
      NO3. Both confirmed distinct from the ionoregulatory-metal
      mechanism → solo.
  - **Organics CA-grouped by shared mechanism:** `nonpolar narcosis`
    (PAHs, BTEX, HCB), `polar narcosis` (phenol, 2-chlorophenol),
    `uncoupler oxidative phosphorylation` (2,4,6-TCP,
    pentachlorophenol), `AChE inhibition: organophosphate` (parathion,
    azinphos-methyl, dimethoate, demeton-S-methyl),
    `neurotoxicant: cyclodiene` (aldrin, dieldrin, endrin,
    g-BHC/lindane), `neurotoxicant: DDT` (DDT, DDE, methoxychlor). Users
    can override the whole scheme by supplying their own metadata CSV.
    Doc comment in `R/mspaf.R` updated; classification lives in
    `inst/extdata/anzecc_analyte_metadata.csv`. Hg+Se are CA-combined as
    a shared sulfhydryl-binding TMoA (both docs cite sulfhydryl/enzyme
    disruption). Possible future refinement: Cr/B/Mn currently stay solo
    because their ANZG docs state no clear primary mode — revisit if
    better mechanism data emerge (e.g. whether Cr(III) belongs with the
    ionoregulatory group).

- **As speciation default — RESOLVED (2026-06-02).** `.SSD_NAME_MAP`
  keeps unspeciated `"As" → "As_V"`, but the rationale is now corrected
  and verified: this is the conservative **max-PAF envelope** of the two
  speciation SSDs, NOT a redox/oxygenation assumption. Empirically
  (checked 0.1–5000 µg/L) the As(V) SSD yields a higher PAF than the
  As(III) SSD at every realistic concentration — they cross only above
  ~10 mg/L dissolved As. So even though landfill leachate is reducing
  (As(III) dominant), As(III) oxidises only slowly in freshwater, and
  As(III) is more toxic to the single most-sensitive diatom, the
  whole-assemblage As(V) SSD sits lower (lower tail pulled down by
  corrected algal acute points) → As(V) is conservative regardless of
  speciation. Comment in `R/paf.R` + CSV note updated; pinned by
  `tests/testthat/test-paf-as-envelope.R` (fails loudly if a future
  SSD-data revision ever lets As(III) overtake As(V), at which point
  switch to a true runtime max-of-both envelope). Users can still pass
  As(III)/As(V) explicitly.

------------------------------------------------------------------------

## 3. Ammonia (NH3-N) normalisation — trace & de-risk double correction — RESOLVED 2026-06-02

Traced both paths and **found + fixed an inverted normalisation
formula**.

- **SSD basis confirmed:** the metadata declares NH3-N as total
  ammonia-N at the ANZG index condition pH 7.0 / 20 °C
  (`dgv_conditions`), toxicity driven by the un-ionised NH3 fraction
  `f = 1/(1+10^(pKa−pH))`, `pKa = 0.09018 + 2729.92/T(K)`.
- **Direction bug (fixed):** the metadata formula was
  `C * f_ref / f_sample`, which normalises a high-pH sample *downward* —
  backwards. A 900 µg/L sample at pH 8.5 holds 100 µg/L un-ionised NH3
  (28× the DGV-equivalent ~3.5 µg/L) yet the old formula scaled it to
  31.9 µg/L (appears harmless). Corrected to `C * f_sample / f_ref`; the
  same sample now normalises to ~25,400 µg/L. This understated ammonia
  risk at high pH and overstated it at low pH — substantial, so the
  corrected direction is pinned by tests in `test-normalise.R`.
- **Two-path asymmetry documented, not unified:**
  [`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
  is the manual path (caller pre-corrects);
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
  auto-applies the metadata formula. Added the exported helper
  **`correct_ammonia_ph_temp(conc, pH, temperature_C)`** so
  [`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
  callers can do the external step with the *identical* maths
  (`test-normalise.R` asserts helper == metadata formula at several pH/T
  points). Roxygen on the helper + the `R/paf.R` header now warn loudly
  NOT to pre-correct before
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
  (double-application).
- **Incidental fix:** repaired CSV encoding corruption (`<c2><b5>`→µ ×8,
  `<e2><80><94>`→— ×1) and a malformed As row (unquoted commas in
  `notes` had spilled it to 25 columns, corrupting As
  `coanalytes_required`/ `normalisation_formula`/`moa_group`). All 43
  rows now parse to 23 columns.

**Follow-up (Chesterton’s fence + design, 2026-06-02):** searched prior
sessions and found ammonia was singled out as “external” not for
toxicological reasons but because its correction uniquely needs
**temperature**, which “may not be measured consistently” — a
`temperature_default_C` fallback was proposed but never built.
Resolution adopted: make temperature **mandatory** instead of
optional. - `add_amspaf(require_temperature = TRUE)` (default): hard
error if any sample reports NH3-N but has no water `temperature` row
(`.assert_temperature_present`). No silent ammonia drop.
`require_temperature = FALSE` for ammonia-free datasets. - New exported
**`get_silo_air_temp(lat, lon, start, end)`** — wraps
[`weatherOz::get_data_drill()`](https://docs.ropensci.org/weatherOz/reference/get_data_drill.html)
(SILO Data Drill, CC-BY 4.0), returns daily mean air temp
`(Tmax+Tmin)/2` in the shape
[`estimate_water_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/estimate_water_temp.md)
wants; disk-cached under `R_user_dir(..,"cache")/silo`. `weatherOz`
added to Suggests. - Air→water imputation kept as the existing per-site
`lm` in
[`estimate_water_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/estimate_water_temp.md)
with its hard ≥5-pairs requirement (decision: no pooled fallback — a
site with no water-temp pairs cannot proceed). - The
[`correct_ammonia_ph_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/correct_ammonia_ph_temp.md)
helper is now documented as a standalone
[`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
spot-check convenience only; `paf.R` header rewritten to state
[`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
applies no chemistry normalisation for any analyte. - Tests:
`test-silo.R` (mocked weatherOz: mean transform, caching, validation),
`test-temperature-mandatory.R` (the gate). Full suite 164 pass / 0 fail.

**Update (2026-06-03):**
[`estimate_water_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/estimate_water_temp.md)
now optionally adds a first-harmonic day-of-year term to the air→water
regression and selects between the air-only and air+season models by
AICc (args `seasonal`, `seasonal_min_n`, `seasonal_min_quarters`);
covered by `test-watertemp`. This captures seasonal hysteresis without
the full Bayesian model.

Still genuinely optional: upgrade
[`estimate_water_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/estimate_water_temp.md)
to a hierarchical Bayesian air→water model (partial pooling + predictive
distributions) when the uncertainty-propagation work lands — would also
dissolve the no-pairs hard stop.

------------------------------------------------------------------------

## 4. Prescreen — detection-frequency escape hatch — RESOLVED 2026-06-04

[`prescreen_analytes()`](https://vorpalvorpal.github.io/leachatetools/reference/prescreen_analytes.md)
drops analytes below a detection-frequency threshold. Rare-but-potent
toxicants (e.g. a pesticide detected in 2% of samples but at
ecotoxicologically significant concentrations) could be screened out.

**Implemented:** a potency-based escape hatch (`potency_keep = TRUE` by
default). An analyte that fails the frequency screen is still kept if
any detected concentration reaches `potency_frac` (default 1) times its
95%-species-protection guideline value (`dgv_95pct_ug_L` in the
metadata). Needs a numeric `value` column in µg/L (matching the DGV
units); only analytes carrying a DGV (toxicants) are eligible, so major
ions are unaffected. Reported via a new `potency_kept` column / cli
message; pinned by `test-prescreen.R`. Analytes with no DGV (a few
no-SSD organics) still rely on manual `protect`.

------------------------------------------------------------------------

## 5. Organics imputation under near-total BDL (open problem)

In WMF leachate, organics are almost entirely below detection, and there
is no defensible way to impute them without at least DOC as a covariate.
Current hurdle: organics imputed only if a DOC-like variable is present.
Options to explore: - Treat organics as essentially “DOC-scaled” — model
detected organics as a ratio to DOC and propagate that ratio (with
uncertainty) to BDL samples. - Report organics as a censored upper bound
rather than a point estimate when detection frequency is below some
floor (honest about non-identifiability). - Skip organics imputation
entirely and document the limitation. Decision needed — currently
flagged, not resolved.

------------------------------------------------------------------------

## 6. Added Risk Approach (ARA) — measure of local adaptation

ARA computes `C_adj = max(C_norm − ref_norm, 0)` to ask “what change are
we imposing on the *local* species pool?” The reference `ref_norm` is
currently a geometric mean of reference-site chemistry. But local
species are selected/adapted to the chemistry they actually experience —
the right “reference” is whatever statistic best represents the
conditions the local pool is adapted to. Think through: geomean vs
median vs a high percentile (chronic exposure may track typical, but
adaptation may track extremes) vs a time-integrated exposure. Decide the
appropriate measure and justify it. Lives in `R/reference.R` /
`R/mspaf.R` (`.amspaf_adjust`).

------------------------------------------------------------------------

## 7. Two SSD representations — document coexistence

The package uses an SSD in two distinct ways that currently coexist
without a clear narrative: - **Analytic log-normal** (HC50 + sigma) for
the CA mixture combination in
[`compute_ca_group_mspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/compute_ca_group_mspaf.md). -
**Model-averaged “multi” curve** (6 distributions via
[`ssdtools::ssd_hp`](https://bcgov.github.io/ssdtools/reference/ssd_hp.html))
for per-analyte diagnostic PAF in
[`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md).
These can disagree (the log-normal is an approximation of the
model-averaged curve). Document why each is used where, quantify the
disagreement on real analytes, and decide whether the CA step should
draw sigma from the model-averaged fit for consistency. Lives in
`R/mspaf.R` + `R/paf.R`.

------------------------------------------------------------------------

## 8. Chronic chemistry aggregation — clarify role in the pipeline

[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/leachatetools/reference/time_weighted_aggregate.md)
(`R/chronic.R`) is value-agnostic (geomean for chemistry, arithmetic
mean for AmsPAF). Clarify and document: - Where it sits relative to
imputation and AmsPAF (Path B: per-sample AmsPAF then time-average). -
Whether users would call it standalone (e.g. to get a chronic-averaged
chemistry profile) outside the AmsPAF pipeline — if so, make that a
first-class, documented entry point. - Options for the temporal
weighting (exponential decay + forward-step duration) vs alternatives.

------------------------------------------------------------------------

## 9. Imputation generalisation refactor (after §1 design settles)

Generalise the imputation engine so it is domain-agnostic (a reusable
Bayesian multivariate-GAM + PCA imputation tool) with a
leachate-specific preset layered on top. End state: three clear
user-facing parts to the package — **LMF**, **imputation**, **AmsPAF**.
Lives in `R/impute.R`. Coordinate with §1 so the draws/uncertainty
interface is designed once.

------------------------------------------------------------------------

## 9a. BUG — PCA train/predict score-scale mismatch (imputation) — FIXED (2026-06)

**Status: FIXED.**
[`.prepare_chem_pca()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-prepare_chem_pca.md)
previously returned `pc_scores` copied from `nipals::nipals()$scores`,
but
[`.compute_pca_scores()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-compute_pca_scores.md)
(used at prediction) produces a regression projection
(`X_scaled %*% loadings` via `.nipals_score_row`). These differ by a
per-component eigenvalue factor
(`nipals$scores[,k] == projection[,k] / eig[k]`), so the brms smooths
were trained on one scale and predicted on another → silently wrong
imputations. Fix:
[`.prepare_chem_pca()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-prepare_chem_pca.md)
now builds `pc_scores` by calling
[`.compute_pca_scores()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-compute_pca_scores.md)
on the training data, guaranteeing train==predict by construction.
Pinned by `tests/testthat/test-impute-pca.R` (asserts equality and that
the scores are NOT the eigenvalue-shrunk nipals scores).

[`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md)
was never affected — it both fits and predicts its GAM with
[`.compute_pca_scores()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-compute_pca_scores.md).

Also fixed in passing: `cum_var <- pca_fit$R2cum` was `NULL` under
nipals 1.0 (per-component `$R2`, no `$R2cum`); now `cumsum(pca_fit$R2)`
with an `R2cum` fallback. Also pinned by `test-impute-pca.R`.

------------------------------------------------------------------------

## 10. In-code TODOs (smaller, well-specified)

- **Probabilistic NO3-N hardness weighting — RESOLVED 2026-06-04**
  (`R/paf.R`). `ssd_paf("NO3-N", hardness_mg_L=)` now blends the
  soft/moderate/hard SSDs by hardness-class probability instead of
  snapping to one class: `p_soft = plnorm(30, log(h), s)`,
  `p_hard = 1 - plnorm(150, log(h), s)`, `p_mod = 1 - p_soft - p_hard`,
  `s = sqrt(log(1+cv^2))`, default `cv = 0.05` (`.no3_weights()`);
  `ExpectedPAF = Σ w·PAF_class`. Continuous across the class boundaries
  (at h=150 PAF is the mean of the mod/hard PAFs); pinned by
  `test-paf-values.R`.
  [`ssd_hc50()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_hc50.md)
  / the CA mixture step still use a single representative class
  ([`.no3_class()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-no3_class.md));
  revisit if CA needs the blend too.

- **Benchmark the three BDL/imputation configs** (`R/impute.R` docstring
  of `fit_imputation_model`). On real hold-out data (mask 10% of
  detected cells, compare RMSE / coverage):

  1.  `rescor = TRUE` + `mi()` + post-hoc DL cap (current).
  2.  `rescor = FALSE` + `cens("left")` (clean BDL; loses cross-analyte
      coupling).
  3.  `rescor = FALSE` + `cens("left")` + shared `(1 | sample_id)`
      latent factor.

- **BDL cap inspection** (`.check_bdl_imputed`, `R/impute.R`). The cap
  clips imputed BDL cells to the detection limit; for sites whose
  chemistry legitimately implies high concentrations the cap can fire
  often. Surface a summary of cap activations so the regime is
  auditable.
