# Calibrating AmsPAF against macroinvertebrate community data

A planning document for an empirical validation and calibration layer on top
of the Adjusted multi-substance Potentially Affected Fraction index (AmsPAF).
The goal is to ground-truth and/or correct the model-derived AmsPAF score
against observed community impact at a paired reference and downstream site
with ~10 years of annual macroinvertebrate monitoring data.

This is exploratory by design. Multiple analytical variants are described,
each testing a different aspect of the question. The intent is to run several
in parallel during EDA and let the data determine which framing is most
informative — not to pre-commit to a single design.

---

## 1. Context and conceptual model

AmsPAF predicts the fraction of aquatic species potentially affected by a
toxicant mixture, based on laboratory species sensitivity distributions
adjusted for local geogenic background (ARA). It's a model-derived exposure
metric — it predicts harm rather than observing it.

The macroinvertebrate data is the inverse: an observed community state
without direct causal attribution to specific stressors. Linking the two
gives us either:

- **Validation**: does AmsPAF actually track observed community impact at
  this site? If so, the index works.
- **Calibration**: a function translating AmsPAF to expected observed
  impact, absorbing site-specific differences between lab-derived sensitivity
  and field reality (community tolerance, bioavailability effects not
  captured by chemistry normalisation, indirect effects, interaction with
  non-toxicant stressors).

### Conceptual model of the two sites

The upstream and downstream sites at this catchment have effectively
identical hydrology and meteorology. They differ in two ways:

- **Geomorphology**: roughly constant on the relevant time scale (years).
  Produces a baseline community difference between sites that doesn't
  vary year to year.
- **Leachate impact**: highly variable in time. Acts only on the downstream
  site.

Both site communities evolve in time, pushed by:

- **Shared forcings** (hydrology, meteorology, regional climate) move both
  communities in *parallel* directions in community space.
- **Geomorphology** creates a constant offset between sites in a direction
  that depends on habitat differences (riffle vs pool species, substrate
  preferences etc.).
- **Leachate** moves only the downstream community.

### The leachate impact has two components

Leachate doesn't act on the community through a single mechanism. It has at
least two distinct effects:

- **Toxicity**: metal and ammonia exposure causing mortality, reproductive
  failure, behavioural disruption. This is what AmsPAF is designed to
  predict.
- **Water quality disruption**: changes in pH, DOC, alkalinity, conductivity,
  dissolved oxygen, hardness etc. — non-toxic mechanisms that nonetheless
  shape community composition. Acidification stresses sensitive taxa; DOC
  enrichment alters light penetration and biofilm composition; conductivity
  shifts affect osmoregulatory species; DO depression in highly organic
  leachate kills sensitive taxa directly.

The leachate mixing fraction (LMF) — your derived index of leachate
contribution to downstream water — is a candidate predictor for the WQ
disruption component, because it tracks general leachate presence rather
than specifically toxic exposure.

AmsPAF and LMF are correlated in time (both rise with leachate impact) but
not identical. AmsPAF can be high even when LMF is moderate (e.g., when
trace metals are elevated without much bulk ion signature, or vice versa).
This correlation creates both opportunity (we can try to separate the two
components) and difficulty (multicollinearity).

### Full decomposition

$$D_{\mathrm{down}}(t) = D_{\mathrm{up}}(t) + \Delta D_{\mathrm{geomorph}} + \Delta D_{\mathrm{tox}}(t) + \Delta D_{\mathrm{wq}}(t) + \mathrm{noise}$$

Taking the site-pair difference cancels the shared D_up(t):

$$\Delta D(t) \equiv D_{\mathrm{down}}(t) - D_{\mathrm{up}}(t) = \Delta D_{\mathrm{geomorph}} + \Delta D_{\mathrm{tox}}(t) + \Delta D_{\mathrm{wq}}(t) + \mathrm{noise}$$

The shared hydrological forcing is removed by the difference design. The
remaining task is to separate ΔD_geomorph (constant), ΔD_tox(t) (predicted
by AmsPAF), and ΔD_wq(t) (predicted by LMF) — and to relate the toxicity
component to AmsPAF.

### Chronic vs per-sample AmsPAF

A further complication: AmsPAF as defined per chemistry sample is an
instantaneous toxicity estimate, while the macroinvertebrate community
integrates exposure over weeks to months. For the calibration to be
meaningful, AmsPAF (and LMF) must be chronic-integrated over a
biologically relevant window before comparing to macroinvertebrate data.

This produces two distinct AmsPAF outputs:

- **Per-sample AmsPAF**: one value per chemistry sample. Used for routine
  monitoring, sample-level reporting, and triggering management responses
  to specific exceedances.
- **Chronic AmsPAF**: a time-integrated value evaluated on any given day
  from the preceding ~365 days of chemistry. Used for the macroinvertebrate
  calibration and trend visualisation.

Both are useful; neither replaces the other. The chronic integration
framework is detailed in §3.6.

---

## 2. Two distinct framings

### Framing A: Predict community impact from AmsPAF

Treat AmsPAF (and/or LMF) as predictor variables. Use a community metric
(SIGNAL2, Shannon, taxon richness) as the response. Fit a regression. Use
it to predict expected community impact at future samples.

**Use case**: communicating ecological consequences of contamination levels
to management or regulators.

**Strength**: Output is biologically meaningful in familiar units.

**Weakness**: Doesn't directly test the index. A correlation between
AmsPAF and SIGNAL2 might exist for reasons unrelated to msPAF's specific
mechanism.

### Framing B: Ground-truth AmsPAF against observed species loss

Treat AmsPAF as a prediction to be tested. msPAF claims "X% of species are
potentially affected by this mixture". The corresponding field observation
is "X% of upstream species are missing or substantially impacted downstream"
— call this PAFobs.

Compare AmsPAF to PAFobs directly. They're in the same units (fraction of
species affected) and should be approximately equal if msPAF is correctly
calibrated to the resident community.

**Use case**: validating and correcting the index methodology itself.

**Strength**: Directly tests msPAF's quantitative claim. Slope and
intercept of the fit have biological meaning (does msPAF over- or
under-predict?).

**Weakness**: Requires careful handling of imperfect detection in the
macroinvertebrate sampling, and disentangling geomorphological from
leachate-driven community difference, and disentangling toxicity-driven
from WQ-driven leachate impact.

**Framing B is the preferred approach.** It's a stronger test, and the
resulting calibration is more interpretable. Framing A can be done as a
secondary analysis.

---

## 3. Design considerations

### 3.1 Site-pair difference removes parallel forcings

As established in §1, the difference ΔD(t) cancels shared hydrological
forcing because both sites move in parallel under it. The difference
contains geomorphology (constant) plus leachate effects (variable) plus
noise.

The orthogonality between hydrology-driven and leachate-driven community
change is biologically plausible: a taxon's flow preference is largely
uncorrelated with its metal tolerance. Strong-flow years and
strong-leachate years push communities in different directions in
multi-dimensional community space.

This is testable. If years with extreme hydrology *don't* show elevated
ΔD after the geomorphological offset is removed, the orthogonality
assumption holds. If they do, hydrology and leachate share some community
direction and the simple decomposition is incomplete.

It doesn't control for site-differential confounders that vary in time:

- Riparian condition: if vegetation works happen at one site only
- Sediment / habitat disturbance: bank erosion or channel works
- Reference site drift: slow community change at reference (recovery from
  past impact, succession, climate-driven shift) introduces structure in
  the difference unrelated to leachate

**Action**: review the 10-year record for years with known differential
disturbances at either site (riparian works, major channel disturbance,
sampling protocol changes). Flag these years for sensitivity analysis.

### 3.2 Separating geomorphology from leachate impact: regression intercept

ΔD_geomorph is constant in time but unknown in magnitude. We need to
estimate it to extract the time-varying leachate components.

**Approach: regression intercept.** Fit ΔD ~ AmsPAF (and/or LMF). The
intercept estimates ΔD_geomorph; the slope(s) estimate the leachate
effect(s). Uses all the data, avoids judgement calls about which years
constitute "baseline", and produces an interpretable intercept that becomes
a diagnostic for whether the geomorphological correction is working.

Alternatives considered and rejected:

- "Low-leachate baseline period": requires pre-specifying what "low
  leachate" means; chronic leachate presence may contaminate all years.
- "Low-AmsPAF years": circular when validating AmsPAF.

Diagnostic interpretation of the intercept α:

- **α ≈ 0**: geomorphology is negligible or already absorbed; the two
  sites would have effectively identical communities without leachate.
- **α > 0**: a stable baseline community difference exists between sites
  that isn't predicted by leachate exposure. Most likely interpretation:
  geomorphological habitat difference.
- **α < 0**: the downstream site has *more* taxa than upstream after
  controlling for predictors. Unexpected; probably indicates a problem
  with the PAFobs definition or the site characterisation.

### 3.3 Disentangling toxicity from WQ disruption: multiple model variants

The two leachate components (toxicity, WQ disruption) create a design
question: should the calibration test AmsPAF's specific contribution to
community impact (controlling for LMF-predicted WQ effects), or its total
correlation with community impact (without distinguishing mechanisms)?

There are arguments for both. Rather than committing to one, fit multiple
variants and use the comparison as a diagnostic. The key variants:

**Variant 1: Single-predictor (AmsPAF only).**

$$\mathrm{logit}(\mathrm{PAF}_{\mathrm{obs}}) = \alpha + \beta \cdot \mathrm{logit}(\mathrm{AmsPAF}) + \varepsilon$$

Simplest. β absorbs whatever community-impact variance correlates with
AmsPAF, including any LMF-correlated effects. Tests "does AmsPAF correlate
with community impact?" Doesn't distinguish mechanism.

Operational interpretation: "for a typical sample with this AmsPAF at
this site, expect this much community impact". Useful for management
reporting.

**Variant 2: Two-predictor (AmsPAF + LMF).**

$$\mathrm{logit}(\mathrm{PAF}_{\mathrm{obs}}) = \alpha + \beta_{\mathrm{tox}} \cdot \mathrm{logit}(\mathrm{AmsPAF}) + \beta_{\mathrm{wq}} \cdot \mathrm{logit}(\mathrm{LMF}) + \varepsilon$$

Sharper test of AmsPAF specifically. β_tox is "AmsPAF's effect on
community impact after controlling for LMF" — i.e., the toxicity-specific
contribution. β_wq captures the WQ-disruption contribution via LMF.

Operational interpretation: "AmsPAF predicts the toxicity-driven
component; LMF predicts the WQ-driven component; both contribute to
community impact". Mechanistically clean but requires understanding both
predictors.

**Variant 3: LMF only.**

$$\mathrm{logit}(\mathrm{PAF}_{\mathrm{obs}}) = \alpha + \gamma \cdot \mathrm{logit}(\mathrm{LMF}) + \varepsilon$$

The null model for AmsPAF: does LMF alone explain community impact? If yes,
AmsPAF's contribution beyond LMF (β_tox in Variant 2) is the relevant
question. If LMF explains community impact poorly, AmsPAF in Variant 1 is
doing the meaningful work.

**Variant 4: Orthogonal decomposition.**

Construct two orthogonal predictors from AmsPAF and LMF:

- A "shared leachate" component (the principal axis of variation common
  to both)
- A "toxicity-specific" component (AmsPAF residualised on LMF)

Fit ΔD against both. The toxicity-specific component slope is the cleanest
test of toxicity-specific community impact. More rigorous than Variant 2
but harder to explain and consumes a degree of freedom on the
orthogonalisation step.

### Diagnostic interpretation from the comparison

The comparison between variants is itself the most informative output:

- **β (Variant 1) ≈ β_tox (Variant 2)**: AmsPAF and LMF don't share much
  community-impact variance. The two leachate mechanisms are largely
  independent and AmsPAF is doing toxicity-specific work.

- **β (Variant 1) ≫ β_tox (Variant 2)**: most of AmsPAF's apparent
  community-impact correlation in Variant 1 was actually WQ-driven (LMF
  doing the work via the AmsPAF-LMF correlation). AmsPAF's specific
  toxicity contribution is smaller than Variant 1 suggests.

- **β_tox ≈ 0 while β > 0 and γ > 0**: AmsPAF has no detectable
  independent contribution to community impact beyond what LMF predicts.
  This is a meaningful negative result for the AmsPAF mechanism.

- **β_tox ≈ β and β_wq ≈ 0**: LMF doesn't add explanatory power. The
  simpler Variant 1 is preferred.

- **All three of β, β_tox, β_wq, γ similar magnitude**: AmsPAF and LMF
  are highly correlated and effectively interchangeable predictors. The
  separation can't be resolved with this dataset.

### Practical constraint: n = 10 and multicollinearity

With only 10 data points, adding LMF as a second predictor takes you
from 2 parameters to 3 parameters. Residual degrees of freedom drop from
8 to 7. Coefficient SEs widen, especially because AmsPAF and LMF will be
correlated.

This is workable but tight. Confidence intervals on β_tox in Variant 2
will be wider than on β in Variant 1. The honest framing is "with this
sample size, we can fit a two-predictor model but the individual
coefficients have substantial uncertainty; the comparison between the
variants is the diagnostic, not a precise decomposition".

If the AmsPAF-LMF correlation is very high (say r > 0.9), multicollinearity
will make Variant 2 unstable and the decomposition won't work. Variant 4
(orthogonal decomposition) is the formal escape valve but with n = 10 it's
also limited.

### 3.4 Community dissimilarity metric: Baselga partitioning

Standard Bray-Curtis dissimilarity is symmetric and conflates two types
of community change:

- **Species replacement**: taxa lost from one site are replaced by
  different taxa at the other site. The communities have similar richness
  but different membership.
- **Species loss / nestedness**: one community is a depleted subset of
  the other. Same taxa, fewer of them.

These represent biologically different processes. msPAF predicts species
*loss* (a fraction of the species pool affected by toxicity). The
biologically corresponding observation is nestedness, not turnover.

Baselga (2010, *Global Ecology and Biogeography* 19(1); 2012, *Global
Ecology and Biogeography* 21(12); 2013, *Methods in Ecology and Evolution*
4(6)) provides a formal decomposition:

**For presence/absence (Jaccard partitioning)**:

$$J_{\mathrm{total}} = J_{\mathrm{turnover}} + J_{\mathrm{nestedness}}$$

where given a = shared taxa, b = taxa at site A only, c = taxa at site B only:

$$J_{\mathrm{total}} = \frac{b + c}{a + b + c}$$

$$J_{\mathrm{turnover}} = \frac{2 \min(b, c)}{a + 2 \min(b, c)}$$

$$J_{\mathrm{nestedness}} = J_{\mathrm{total}} - J_{\mathrm{turnover}}$$

**For abundance data (Bray-Curtis partitioning)**:

$$BC_{\mathrm{total}} = BC_{\mathrm{balanced}} + BC_{\mathrm{gradient}}$$

where BC_balanced is the analogue of turnover (abundance replacement)
and BC_gradient is the analogue of nestedness (net abundance loss).

**For this analysis**:

- **J_nestedness** (or BC_gradient) is the primary PAFobs metric — it
  matches msPAF's species-loss framing
- **J_turnover** (or BC_balanced) is a diagnostic — it captures community
  shift unrelated to species loss
- The conceptual model predicts that J_turnover should be roughly
  constant in time (it reflects geomorphological community difference,
  which is stable) and J_nestedness should vary with leachate impact (it
  reflects toxic community depletion)
- If J_turnover varies substantially in time, either geomorphology is
  changing, or there is a time-varying community-replacement effect not
  predicted by the model (perhaps the WQ component of leachate is more
  about replacement than loss?)

If J_turnover is stable in time and J_nestedness varies with AmsPAF/LMF,
the framework is working as predicted. If both vary, more thought is
needed.

R package `betapart` (Baselga & Orme 2012, *Methods in Ecology and
Evolution* 3(5)) implements the partitioning.

### 3.5 Diversity metric choice for Framing A

The right metric depends on the framing:

**For Framing A (predict community impact)**: SIGNAL2 (Chessman 2003,
*Marine and Freshwater Research* 54(2)) is the Australian standard,
designed for pollution monitoring, weights taxa by pollution sensitivity.
EPT richness (Ephemeroptera + Plecoptera + Trichoptera) is a simpler
alternative used in NSW DPI monitoring.

**For Framing B (ground-truth)**: J_nestedness (or BC_gradient) from
Baselga partitioning, as described in §3.4. Don't use SIGNAL2 — it
doesn't have the right structure (weighted average sensitivity score,
not a count of affected taxa).

Don't pick the metric that gives the best correlation with AmsPAF; that's
circular. Pick on biological grounds first, then fit.

### 3.6 Chronic exposure integration

#### The mismatch problem

PAF and SSDs are derived from chronic toxicity endpoints — exposures over
days to weeks. The chemistry data is sporadic point-in-time measurements:
a slowly-varying baseline component (a reasonable proxy for chronic
exposure) interrupted by pulses with recession back to baseline. Sampling
is biased toward pulses because routine monitoring is supplemented by
event-triggered sampling when leachate is visibly elevated.

The macroinvertebrate community integrates exposure over its generation
time (weeks to months for most families). To relate chemistry to
community, the chemistry needs to be integrated over time as well, with
careful handling of:

- Time integration matching biological response timescales
- Sampling bias toward pulses (raw mean over-weights pulse events)
- Differentiation between same-mean low-variance vs same-mean
  high-variance (pulse) exposure

The per-sample AmsPAF computed from each chemistry sample is not what
the macroinvertebrate community responds to. A chronic-integrated AmsPAF
is needed for calibration purposes. The per-sample AmsPAF remains the
appropriate output for chemistry-sample-level reporting; the
chronic-integrated AmsPAF is computed separately for calibration.

#### Recommended integration framework

For each calendar day t_0 within the analysis window:

**(a) Window**: include all samples from t_0 - 365 days to t_0 (one year
of look-back). Extension beyond 365 days adds little under exponential
decay weighting; 365 days is a reasonable balance.

**(b) Temporal weighting**: exponential decay with characteristic time
τ. Weight for a sample at time t_i:

$$w_{\text{temporal},i} = \exp\left(-\frac{t_0 - t_i}{\tau}\right)$$

τ = 90 days is the suggested default. This gives roughly:

- Days 0–90 before t_0: ~63% of total weight
- Days 90–180: ~23%
- Days 180–270: ~9%
- Days 270–365: ~3%

The biological rationale is that organisms have a "memory" of past
exposure that decays as cells are renewed, generations turn over, and
damaged individuals are replaced. τ ≈ 90 days roughly corresponds to a
macroinvertebrate generation cycle.

Sensitivity-test τ ∈ {30, 60, 90, 120, 180}. The literature doesn't
strongly motivate a specific value; treat τ as a smoothing parameter
and report robustness.

**(c) Sample-duration weighting (sampling bias correction)**: the time
each sample represents must be accounted for, otherwise pulse-triggered
sampling artificially inflates the integrated mean.

Use a **forward-step rule**: each sample i represents the period from
its own time t_i to the time of the next sample t_{i+1}. The duration
weight is:

$$\Delta t_i = t_{i+1} - t_i \text{ (or } t_0 - t_i \text{ for the most recent sample)}$$

The rationale: when a pulse begins, sampling typically intensifies
*after* the pulse is observed, not before. So a baseline sample's
"representativeness" extends forward until the next sample is taken,
which is usually the start of intensified pulse sampling. Symmetric
trapezoidal weighting (half before, half after) would incorrectly assign
some of the baseline period to pulse-era weight.

**Window-edge handling**: at t_0 - 365, fetch the sample immediately
before the window start (from the broader historical record, not
artificially constrained). This anchors the start of the window so the
first in-window sample doesn't have undefined duration. Years of data
exist; just look up the previous sample.

**(d) Combined weight per sample**:

$$w_i = w_{\text{temporal},i} \times \Delta t_i$$

**(e) Per-analyte summary statistic — time-weighted geometric mean**:

$$\bar{C}_{\text{geom}} = \exp\left( \frac{\sum_i w_i \log(C_i)}{\sum_i w_i} \right)$$

The geometric mean is the natural average for log-normally distributed
concentrations (which most environmental concentrations approximate) and
is less sensitive to extreme pulses than arithmetic mean while still
moving in response to elevated periods. It matches the chronic-exposure
framing of PAF better than arithmetic mean.

**(f) BDL handling**:

The principled approach is to treat BDL values as **left-censored
observations in the Bayesian multivariate GAM imputation step** (which
operates upstream of the chronic integration in the AmsPAF pipeline).
`brms` supports this natively via `cens()`:

```r
chemistry_data <- chemistry_data |>
  mutate(
    censored = if_else(quantified, "none", "left"),
    value_for_likelihood = if_else(quantified, log(value), log(DL))
  )

bf_pb <- bf(value_for_likelihood | cens(censored) ~ s(pH) + s(log_ec) + ... )
```

The `cens()` indicator tells `brms` "this observation is left-censored at
the value given" — the true log-concentration is somewhere below log(DL).
For a multivariate model with `set_rescor(TRUE)`, the residual covariance
captures cross-analyte correlations, which means BDL values for one
analyte are informed by quantified values of correlated analytes.

This replaces the conventional DL/2 substitution. The posterior for a
censored value gives a distribution over [0, DL] informed by:

- The censoring constraint itself (must be below DL)
- The model's smooth functions of drivers (pH, EC, DOC etc.)
- Cross-analyte correlations through the residual covariance
- Any quantified observations of the same analyte at other samples
- Priors

No DL/2 substitution is needed. The reference/downstream asymmetry from
earlier drafts becomes unnecessary because the imputation handles
censoring properly for both sites within a single framework.

**Analyte inclusion criterion**: an analyte participates in msPAF only if
quantified above DL in at least k% of samples (suggest k = 5–10) over
the analysis window. Below this threshold, the analyte is excluded from
msPAF entirely.

The rationale: an analyte never quantified contributes no information
about presence at our sites — only a prior-driven posterior. The
multivariate model handles such analytes correctly (produces small
posterior values) but the computational and interpretational cost isn't
worth it. Explicit exclusion is honest documentation of analytical
limitations rather than a hidden imputation.

For Blaxland this matters most for the **pesticide class**. If most/all
pesticides have never been detected at our sites, they should be
excluded from msPAF. The reported msPAF then represents metal + ammonia
(+ H2S where measured) toxicity only, with a note that pesticide
contamination wasn't detected over the analysis window. If future
sampling produces pesticide detections, the inclusion criterion is
re-evaluated and pesticides re-enter the index naturally.

For analytes that are mostly-but-not-entirely-BDL (say, quantified in
10–50% of samples), the Bayesian framework is most valuable. The
quantified observations anchor μ and σ; the censored observations
provide the constraint; cross-correlations with measured analytes help
fill in the BDL gaps. This is real value-add over DL/2 substitution.

**Fallback for the interim period**: if the Bayesian multivariate model
isn't yet implemented, fall back to DL/2 for downstream and 0 for
reference, accepting the convention and its known limitations. Migrate
to the censored-likelihood approach when the Bayesian framework comes
online.

**(g) Pipeline order**:

1. Bayesian multivariate GAM imputation fills in missing analytes at each
   raw sample (the imputation step already in the AmsPAF pipeline). This
   means every sample has full analyte coverage for chronic integration.
2. Compute time-weighted geometric mean per analyte over 365-day window
   ending at t_0.
3. Apply ARA shift using equivalently-integrated reference concentrations.
4. Evaluate SSDs at ARA-adjusted chronic concentrations to get per-analyte
   chronic PAF.
5. Combine per-analyte PAFs into chronic AmsPAF via CA + IA as for the
   per-sample version.

**(h) Apply same window and weighting to reference**: chronic-integrated
reference concentrations match the chronic-integrated downstream
concentrations in time and weighting. ARA shift then becomes a chronic-
chronic subtraction, consistent in time.

**(i) Same treatment for LMF**: LMF (leachate mixing fraction) needs the
same chronic integration for the multi-predictor variants to be coherent.
LMF is a sample-level value; compute its time-weighted geometric mean
over the same 365-day window using the same τ.

#### Output structure

The macroinvertebrate calibration uses chronic AmsPAF and chronic LMF
evaluated on the day of biological sampling. With ~10 years of annual
biological samples, this yields ~10 paired (chronic AmsPAF, chronic LMF,
PAFobs) observations for the calibration fit.

The chronic-integrated AmsPAF computed daily across the record is also
useful for trend visualisation and as a smoothed version of the
per-sample series. Both per-sample and chronic AmsPAF should be retained
as separate outputs.

#### Sensitivity tests

- **τ** ∈ {30, 60, 90, 120, 180}: how much does the calibration shift?
  Stability across this range indicates the integration is robust.
- **Window length** ∈ {180, 365, 730}: with exponential decay, longer
  windows add little but check the boundary.
- **Summary statistic**: time-weighted geometric mean vs arithmetic mean
  vs 90th percentile. The differences here indicate how much pulse vs
  baseline weighting matters.
- **BDL handling**: Bayesian censored vs DL/2 substitution vs DL/10
  substitution. If the calibration shifts substantially, the Bayesian
  framework's information gain over DL/2 is real. If little difference,
  the simpler substitution is operationally adequate.
- **Inclusion threshold k**: try k ∈ {5%, 10%, 20%} for the minimum
  detection frequency. This determines which marginal analytes are
  included or excluded from msPAF.

#### Optional refinement: pulse-aware duration weighting

The forward-step duration weighting assumes the next sample represents
the start of a new exposure period (typically the start of a pulse). If
pulses are short and recovery is fast, this is fine. If pulses are long
and recovery happens within sampling gaps, the forward-step approach
may underestimate the duration of recovery to baseline.

A more sophisticated approach: detect pulses empirically (sample value
> 3× preceding baseline median, or similar) and use forward-step weighting
between baseline samples, but trapezoidal weighting between baseline and
pulse samples. This is harder to implement consistently and the gain is
probably small. Recommend deferring until sensitivity testing on the
simple forward-step approach indicates it's worth the additional
complexity.

#### Note on recruitment dynamics

The upstream site is near-pristine and well-connected to downstream
ecologically. Sensitive taxa lost from downstream during exposure events
recolonise from upstream once conditions improve, provided the impact
period is shorter than the recolonisation time.

At the 180–365 day integration timescale, recruitment is generally not
limiting — biota will have had time to attempt recolonisation. The
chronic-integrated AmsPAF therefore captures something close to "exposure
intensity that prevented recolonisation", which is exactly what the
species-loss interpretation of msPAF requires. This is a strength of the
chronic integration framework at this site.

Short pulses followed by long recovery periods will have low
chronic-integrated AmsPAF (geometric mean stays near baseline),
correctly predicting low community impact at the next biological sample.
Sustained moderate elevation will have moderate chronic-integrated
AmsPAF, correctly predicting persistent community impact. The framework
handles these two cases appropriately without explicit recovery modelling.

### 3.7 Imperfect detection

Single-event sampling at a site captures a subset of the true community.
Two equally-clean sites with the same 10-family pool can record different
7-family subsets just by sampling stochasticity, producing apparent
nestedness/loss of 30–40% even with zero actual impact. This is the
biggest threat to the lost-taxa calculation. Approaches in order of
complexity:

**(a) Restrict to common taxa.** Use only families detected at the
reference site in ≥ k of n sampling events. Common families have high
detection probability; partitioning on this restricted set is much less
inflated by stochasticity. With n=10, k=5 is a reasonable starting point.

**(b) Multi-year occupancy comparison.** For each family, count years
detected at each site. Lost-occupancy for a family = (years at reference -
years at downstream) / (years at reference). Aggregate across families.
Gives one number per 10-year window — loses year-by-year resolution.

**(c) Multi-species occupancy model.** Formal latent-variable model
separating true occupancy from detection probability per taxon (MacKenzie
et al. 2002, *Ecology* 83(8); Royle & Dorazio 2008, *Hierarchical Modeling
and Inference in Ecology*). With repeated sampling (across years), can
disentangle true occupancy from detection. R packages: `unmarked`
(frequentist), `spOccupancy` (Bayesian). Substantial modelling exercise
but the formally correct approach.

**Recommendation**: start with (a). It's a one-line filter. If the
restricted PAFobs in clearly clean years is still > 20%, escalate to (c).

### 3.8 Abundance vs presence/absence

Abundance data adds information when reliable:

- Distinguishes "barely there" from "common" taxa
- Detects reductions short of extirpation
- Better signal-to-noise for common taxa

It doesn't fix imperfect detection — both abundance-weighted and
presence/absence metrics record absence as zero abundance without knowing
whether that's true absence or non-detection.

For Baselga partitioning, presence/absence (Jaccard) and abundance
(Bray-Curtis) versions both exist. Use BC if abundance data is reliable;
use Jaccard otherwise.

**Caveat on abundance quality**: AUSRIVAS protocols often record relative
abundance categories (rare/common/abundant) rather than exact counts. If
your data is categorical, treat it as ordinal and avoid arithmetic
operations on it (Bray-Curtis on category numbers is mathematically
suspect).

**Action**: check the abundance format in your data before deciding
between Jaccard and Bray-Curtis partitioning.

### 3.9 Taxonomic resolution

Family-level data is much cleaner than species-level for invertebrates:

- Identification reliability is much better at family
- Cross-year taxonomic consistency is high at family
- Species-level differences include many spurious (cryptic species,
  identification drift between years/taxonomists)

Less sensitive to subtle community changes, but more robust for
calibration. Use family-level unless there's a specific reason to go finer.

---

## 4. Detailed analysis plan

### Stage 1: Data assembly and inspection

**1.1 Compile macroinvertebrate data into a tidy format.**

Per sampling event:

- Site (reference vs downstream)
- Year (and ideally date)
- Family (or other taxonomic unit)
- Abundance (count, category, or presence/absence — record what you have)

Resulting in one row per family per sampling event.

**1.2 Compute chronic AmsPAF and chronic LMF.**

For each macroinvertebrate sampling date, compute chronic AmsPAF and
chronic LMF following §3.6: 365-day window, exponential decay with
τ = 90 days, forward-step duration weighting, time-weighted geometric
mean per analyte, then through the AmsPAF pipeline (ARA shift + SSD
evaluation + CA/IA combination).

Compute the same chronic values at the reference site for consistency
with the ARA framework.

Compute the empirical correlation between chronic AmsPAF and chronic
LMF across years. If r > 0.9, multicollinearity will limit the
multi-predictor variants substantially.

**1.3 Inspect sampling effort.**

Check that sampling protocols (method, season, replicates per visit)
are consistent across years and matched between sites. Note any
discontinuities or known anomalies. Flag years with differential
disturbance (riparian works, channel modification, major flood, fire,
sampling protocol change).

**1.4 Define the site-pool of taxa.**

For each site, the cumulative family list across all 10 years. Identify
"core" families (detected in ≥ 5 of 10 years at reference). These are the
robust set for analysis. Record the per-family detection frequency at the
reference site — this is the empirical detection probability and tells
you how much stochasticity to expect.

### Stage 2: Compute Baselga partitioning

For each year, compute J_nestedness, J_turnover, and J_total (and
BC equivalents if abundance is reliable) using the `betapart` package.
Restrict to the core families from Stage 1.4.

This gives three time series:

- J_nestedness(t): community species loss (target for AmsPAF / LMF
  calibration)
- J_turnover(t): species replacement (diagnostic; should be roughly
  constant if geomorphology is constant)
- J_total(t): overall community dissimilarity (the sum)

### Stage 3: Diagnostic plots, no fitting

Look at the data before fitting anything.

**3.1 Time series.** Plot J_nestedness(t), J_turnover(t), and J_total(t)
over time. Plot AmsPAF and LMF (each aggregation) on the same time axis.

Diagnostic checks:

- Is J_turnover roughly constant? If yes, the conceptual model is
  supported — geomorphology is stable. If no, investigate before
  proceeding.
- Does J_nestedness vary in time? If not, there's no signal to calibrate.
- Does J_nestedness visually track AmsPAF, LMF, or both? If only one,
  that's the dominant driver. If both, the multi-predictor variants are
  needed to disentangle.

**3.2 Scatter: PAFobs vs each predictor (each aggregation).** With n=10
paired points. Visual inspection for monotonic relationship before
fitting. Do this for AmsPAF and LMF separately.

**3.3 Robustness across aggregations.** Compute Spearman ρ between
J_nestedness and each AmsPAF/LMF aggregation. If similar across
aggregations, signal is robust. If only one specific aggregation
correlates, view with suspicion.

**3.4 AmsPAF-LMF correlation.** Plot AmsPAF vs LMF directly. How
correlated are they? If highly correlated (r > 0.9), the multi-predictor
variants will struggle with multicollinearity. If moderately correlated
(0.5 < r < 0.9), the decomposition is feasible. If weakly correlated
(r < 0.5), the two predictors capture different aspects of leachate impact
and the multi-predictor variant is most informative.

**3.5 Orthogonality check.** Identify years with extreme hydrology
(major floods, severe droughts). Do these years show elevated
J_nestedness *after* controlling for AmsPAF and LMF? If yes, hydrology
is not fully orthogonal to leachate in its community effect, and the
simple decomposition is incomplete. This would be a warning rather than
a showstopper.

**3.6 Sensitivity to taxon set.** Compute J_nestedness using (a) all
detected families, (b) core families (≥ 50% detection at reference).
Compare. If they agree on the AmsPAF-PAFobs relationship, detection
robustness isn't a problem. If they differ substantially, the detection
issue matters and may require escalating to occupancy modelling.

**3.7 Sampling protocol check.** Compare years with and without known
disturbance flags from Stage 1.3. Are flagged years outliers in the
predictor-PAFobs relationships? If so, consider excluding them from the
main fit.

### Stage 4: Fit multiple model variants

Fit each of the variants from §3.3 if Stage 3 supports the design. Use
logit-logit on bounded [0,1] variables:

**Variant 1** (AmsPAF only):

```
logit(J_nestedness) ~ logit(AmsPAF)
```

**Variant 2** (AmsPAF + LMF):

```
logit(J_nestedness) ~ logit(AmsPAF) + logit(LMF)
```

**Variant 3** (LMF only):

```
logit(J_nestedness) ~ logit(LMF)
```

**Variant 4** (orthogonal decomposition):

```
shared    <- prcomp(cbind(logit(AmsPAF), logit(LMF)))$x[, 1]
tox_only  <- residuals(lm(logit(AmsPAF) ~ logit(LMF)))
logit(J_nestedness) ~ shared + tox_only
```

For each variant report:

- Coefficients with 95% CIs
- R² and adjusted R²
- Leave-one-out cross-validation prediction RMSE
- Residual diagnostics

### Stage 5: Compare variants

The comparison is the primary diagnostic, per §3.3:

- If β (V1) ≈ β_tox (V2): AmsPAF and LMF are doing independent work
- If β (V1) ≫ β_tox (V2): much of V1's apparent AmsPAF signal was
  actually LMF-correlated
- If β_tox (V2) ≈ 0 while β (V1) > 0 and γ (V3) > 0: AmsPAF doesn't
  contribute beyond LMF
- If β_tox (V2) ≈ β (V1) and β_wq (V2) ≈ 0: LMF doesn't add explanatory
  power; V1 is the parsimonious choice
- If all coefficients are similar magnitude with overlapping CIs: AmsPAF
  and LMF are too correlated to separate with this data

Report all variants. Select the operational one based on:

- Which variant has the lowest LOO-CV RMSE
- Which variant has parameters with non-overlapping CIs from zero
- Which variant matches your operational interpretation (mechanism-specific
  vs total community impact prediction)

### Stage 6: Sensitivity tests

Whatever variant is selected, sanity-check robustness:

- Exclude disturbance-flagged years from Stage 1.3
- Vary the chronic integration parameters:
  - τ ∈ {30, 60, 90, 120, 180} days
  - Window length ∈ {180, 365, 730} days
  - Summary statistic: time-weighted geometric mean vs arithmetic mean
    vs 90th percentile
  - BDL handling: Bayesian censored vs DL/2 vs DL/10
  - Analyte inclusion threshold: k ∈ {5%, 10%, 20%} minimum detection
    frequency
- Try different core taxon thresholds (k = 4, 5, 6)
- Try Jaccard vs Bray-Curtis partitioning

If the fit and conclusions are robust across all sensitivity tests, the
calibration is defensible. If conclusions shift substantially with any
of these choices, the underlying signal is fragile and should be reported
with that caveat.

The τ sensitivity is particularly diagnostic. If results are stable
across τ ∈ {30, ..., 180}, the chronic integration is robust. If only
one specific τ gives a clean calibration, the framework is fitting
something specific to that timescale and may not generalise.

### Stage 7: Reporting

Don't replace AmsPAF with calibrated impact. Report both AmsPAF and the
calibrated community impact prediction, including which variant was used
and why:

> "AmsPAF = X% (toxicity prediction); calibrated species loss at this
> site historically: J_nestedness ≈ Y% (95% CI: a to b%), based on
> [Variant chosen], with a baseline geomorphological offset of Z%."

The calibration is a translator with explicit uncertainty, not a
replacement.

### Stage 8: Annual refit

Each new year of macroinvertebrate data:

- Compute J_nestedness for the new year
- Plot the new (AmsPAF, LMF, J_nestedness) point against the existing
  fitted line(s)
- If within prediction interval, calibration is stable
- If outside, investigate the year (anomaly? methodological change?)
  before refitting
- Update parameter estimates and uncertainties annually
- Check whether α (the estimated geomorphological offset) is stable
  across refits. If it drifts, the underlying assumption that
  geomorphology is constant is being violated.
- Re-evaluate variant choice as data accumulates; with n = 15 or 20
  years the multi-predictor variants become more reliably fittable

---

## 5. Specific implementation outline

```r
# Pseudocode for the EDA — translate to actual data structures
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(betapart)

# ---- Stage -1: Pre-screen analytes for inclusion in msPAF ----
#
# An analyte participates in msPAF only if quantified above DL in at least
# k% of samples at the analysis window. This excludes analytes that have
# never been detected (e.g. pesticides at a non-industrial landfill).
# Explicit exclusion is honest documentation of analytical limitations
# rather than a hidden imputation of small prior-driven values.

K_INCLUSION <- 0.05   # 5% minimum detection frequency

included_analytes <- chemistry_per_sample |>
  group_by(analyte) |>
  summarise(
    n_samples     = n(),
    n_quantified  = sum(quantified, na.rm = TRUE),
    detect_freq   = n_quantified / n_samples,
    .groups       = "drop"
  ) |>
  filter(detect_freq >= K_INCLUSION) |>
  pull(analyte)

# Analytes excluded by this criterion are documented but not imputed.
# The msPAF reported is over included_analytes only.
chemistry_per_sample <- chemistry_per_sample |>
  filter(analyte %in% included_analytes)

# ---- Stage 0: Chronic integration of AmsPAF and LMF ----
#
# For each biological sampling date, compute time-weighted geometric mean
# concentrations over a 365-day window with exponential decay
# (tau = 90 days), then feed through the AmsPAF pipeline.
#
# Inputs assumed:
#   chemistry_per_sample - one row per sample x analyte, with quantified flag
#                          and detection limit. BDL values have been
#                          imputed by the upstream Bayesian multivariate
#                          GAM (treated as left-censored observations).
#                          Missing analytes have also been imputed there.
#   sample_dates         - one row per sample with date and site

TAU_DAYS    <- 90
WINDOW_DAYS <- 365

compute_chronic_concentration <- function(analyte_data, focal_date,
                                          site, tau = TAU_DAYS,
                                          window = WINDOW_DAYS) {
  # Subset to this site and window. Include the sample immediately
  # preceding the window start as the anchor for forward-step duration
  # weighting.
  site_data <- analyte_data |>
    filter(site == !!site) |>
    arrange(date)

  window_start <- focal_date - window
  in_window    <- site_data$date >= window_start & site_data$date <= focal_date

  # Anchor sample: most recent sample before window_start
  anchor_idx <- which(site_data$date < window_start) |> max()
  if (is.finite(anchor_idx)) {
    use_idx <- c(anchor_idx, which(in_window))
  } else {
    use_idx <- which(in_window)
  }
  d <- site_data[use_idx, ]
  if (nrow(d) == 0) return(NA_real_)

  # BDL handling: the upstream Bayesian multivariate GAM imputation step
  # treats BDL values as left-censored observations and produces posterior
  # estimates for them. By the time data arrives at this function, BDL
  # values have already been replaced with posterior estimates (e.g.,
  # posterior means or draws).
  #
  # Fallback for interim use before the Bayesian step is implemented:
  # reference site BDL -> 0, downstream site BDL -> DL/2.
  d <- d |>
    mutate(C = if_else(quantified, value, value))   # value already imputed

  # Forward-step duration: t_{i+1} - t_i, with last sample extending to
  # focal_date
  d$next_date <- c(d$date[-1], focal_date)
  d$dt       <- as.numeric(d$next_date - d$date)

  # Trim the anchor sample's duration so it only covers from window_start
  # forward (not its full forward-step duration which may extend back
  # beyond the window)
  if (is.finite(anchor_idx)) {
    d$dt[1] <- as.numeric(d$next_date[1] - window_start)
  }

  # Drop any samples now with non-positive duration (edge case)
  d <- d |> filter(dt > 0)
  if (nrow(d) == 0) return(NA_real_)

  # Exponential temporal weight at the midpoint of the sample's interval
  midpoint <- d$date + d$dt / 2
  d$w_time <- exp(-as.numeric(focal_date - midpoint) / tau)

  # Combined weight: temporal x duration
  d$w <- d$w_time * d$dt

  # Time-weighted geometric mean. Guard against C = 0 (entire window at
  # reference BDL) with a small epsilon.
  eps <- 1e-9
  exp(sum(d$w * log(d$C + eps)) / sum(d$w))
}

# Apply per analyte, per focal (biological sampling) date, per site
chronic_concentrations <- biological_sampling_dates |>
  crossing(analyte = unique(chemistry_per_sample$analyte)) |>
  crossing(site    = c("reference", "downstream")) |>
  rowwise() |>
  mutate(
    C_chronic = compute_chronic_concentration(
      chemistry_per_sample |> filter(analyte == !!analyte),
      focal_date = date,
      site       = site
    )
  ) |>
  ungroup()

# Feed chronic concentrations through the AmsPAF pipeline:
# - apply ARA shift using chronic reference concentrations
# - evaluate SSDs at ARA-adjusted chronic concentrations to get
#   per-analyte chronic PAF
# - combine via CA + IA to get chronic AmsPAF
#
# Same procedure for LMF: chronic LMF is the time-weighted geometric mean
# of per-sample LMF over the 365-day window.
#
# Output (one row per biological sampling date):
#   chronic_amspaf
#   chronic_lmf

# ---- Stage 1: Data preparation ----

# Define core families (detected at reference in >= 5 of 10 years)
core_families <- macroinverts |>
  filter(site == "reference") |>
  group_by(family) |>
  summarise(years_detected = n_distinct(year)) |>
  filter(years_detected >= 5) |>
  pull(family)

# Check chronic AmsPAF-LMF correlation (using chronic-integrated values)
cor(chronic_predictors$chronic_amspaf, chronic_predictors$chronic_lmf,
    method = "spearman")

# ---- Stage 2: Baselga partitioning per year ----

compute_partition_year <- function(df_year, use_abundance = FALSE) {
  mat <- df_year |>
    filter(family %in% core_families) |>
    pivot_wider(
      id_cols     = site,
      names_from  = family,
      values_from = if (use_abundance) "abundance" else "presence",
      values_fill = 0
    ) |>
    column_to_rownames("site") |>
    as.matrix()

  if (nrow(mat) != 2) return(NULL)

  if (use_abundance) {
    bp <- bray.part(mat)
    tibble(
      nestedness = bp$bray.gra[1],
      turnover   = bp$bray.bal[1],
      total      = bp$bray[1]
    )
  } else {
    bp <- beta.pair(mat, index.family = "jaccard")
    tibble(
      nestedness = bp$beta.jne[1],
      turnover   = bp$beta.jtu[1],
      total      = bp$beta.jac[1]
    )
  }
}

partition_by_year <- macroinverts |>
  group_by(year) |>
  group_modify(\(.x, .y) compute_partition_year(.x, use_abundance = FALSE)) |>
  ungroup()

# ---- Stage 3: Join with chronic AmsPAF and chronic LMF ----

calibration_df <- partition_by_year |>
  left_join(chronic_predictors, by = "year")

# ---- Stage 4: Diagnostic plots ----

# Time series — does J_turnover look stable?
ggplot(calibration_df, aes(year)) +
  geom_line(aes(y = nestedness, colour = "Nestedness")) +
  geom_line(aes(y = turnover,   colour = "Turnover")) +
  geom_line(aes(y = total,      colour = "Total")) +
  labs(y = "Jaccard component", colour = NULL)

# Chronic AmsPAF-LMF correlation
ggplot(calibration_df, aes(chronic_amspaf, chronic_lmf)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Chronic AmsPAF (%)", y = "Chronic LMF (%)")

# PAFobs vs each chronic predictor
calibration_df |>
  select(year, nestedness, chronic_amspaf, chronic_lmf) |>
  pivot_longer(c(chronic_amspaf, chronic_lmf),
               names_to = "predictor", values_to = "value") |>
  ggplot(aes(value, nestedness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ predictor, scales = "free_x") +
  labs(y = "J_nestedness")

# ---- Stage 5: Fit multiple variants ----

clip_for_logit <- function(p, eps = 0.001) {
  pmax(pmin(p, 1 - eps), eps)
}

calibration_df <- calibration_df |>
  mutate(
    logit_pafobs = qlogis(clip_for_logit(nestedness)),
    logit_amspaf = qlogis(clip_for_logit(chronic_amspaf / 100)),
    logit_lmf    = qlogis(clip_for_logit(chronic_lmf / 100))
  )

# Variant 1: AmsPAF only
v1 <- lm(logit_pafobs ~ logit_amspaf, data = calibration_df)

# Variant 2: AmsPAF + LMF
v2 <- lm(logit_pafobs ~ logit_amspaf + logit_lmf, data = calibration_df)

# Variant 3: LMF only
v3 <- lm(logit_pafobs ~ logit_lmf, data = calibration_df)

# Variant 4: orthogonal decomposition
pca       <- prcomp(cbind(calibration_df$logit_amspaf,
                          calibration_df$logit_lmf))
shared    <- pca$x[, 1]
tox_only  <- residuals(lm(logit_amspaf ~ logit_lmf, data = calibration_df))
v4 <- lm(logit_pafobs ~ shared + tox_only, data = calibration_df)

# Compare variants
list(v1 = v1, v2 = v2, v3 = v3, v4 = v4) |>
  lapply(\(m) list(
    coef_ci   = confint(m),
    r_squared = summary(m)$r.squared,
    adj_r2    = summary(m)$adj.r.squared
  ))

# LOO-CV for each
loo_rmse <- function(formula, data) {
  preds <- sapply(seq_len(nrow(data)), \(i) {
    fit_i <- lm(formula, data = data[-i, ])
    plogis(predict(fit_i, newdata = data[i, ]))
  })
  obs <- data$nestedness
  sqrt(mean((preds - obs)^2))
}

loo_rmse(logit_pafobs ~ logit_amspaf, calibration_df)              # V1
loo_rmse(logit_pafobs ~ logit_amspaf + logit_lmf, calibration_df)  # V2
loo_rmse(logit_pafobs ~ logit_lmf, calibration_df)                 # V3
```

---

## 6. Pros and cons of the overall approach

### Pros

- Directly tests msPAF's quantitative claim against field data
- Site-pair difference design cleanly removes shared hydrology
- Baselga partitioning separates species loss (the msPAF target) from
  species replacement (the geomorphology target), giving an additional
  diagnostic
- Regression intercept absorbs the geomorphological offset without
  requiring an arbitrary "baseline period" definition
- Multi-variant approach lets the data determine whether AmsPAF and LMF
  capture independent or overlapping aspects of leachate impact
- Output parameters have direct biological interpretation
- Easy to update annually as new data accumulates

### Cons

- n = 10 paired years is statistically tight. Multi-predictor variants
  in particular have wide CIs and are vulnerable to single influential
  years
- The calibration is site-specific. Doesn't transfer to other monitoring
  locations without separate data
- Imperfect detection adds noise that needs explicit handling
- Slow drift at the reference site (succession, climate-driven changes)
  introduces structure in ΔD unrelated to leachate, indistinguishable
  from time-varying geomorphology
- The orthogonality assumption between hydrology and leachate effects
  may not hold perfectly; if it fails, the difference design over-corrects
  or under-corrects for hydrology
- AmsPAF and LMF may be too correlated to separate with this data, in
  which case the multi-predictor variants don't add information beyond
  the single-predictor ones
- One year of unusual conditions materially affects the fit

### What success would look like

- J_turnover roughly stable in time (geomorphology stable, as expected)
- J_nestedness varies in time and is predicted by at least one variant
- One or more variants show non-zero coefficients with credible intervals
  excluding zero
- Cross-validation prediction RMSE smaller than meaningful differences
  in J_nestedness (say < 0.1 in proportion units)
- Robust across multiple time aggregations and taxon-set choices
- The comparison between variants is informative (clear winner in CV
  performance, or clear story about mechanism)

### What would suggest abandoning the calibration

- J_turnover varies substantially in time (geomorphology isn't constant,
  or partitioning isn't separating things cleanly)
- Scatter plot shows no monotonic relationship between J_nestedness and
  either AmsPAF or LMF
- Cross-validation RMSE comparable to the standard deviation of
  J_nestedness itself for all variants
- Pathological fits (slopes > 5 or < 0.1, very wide CIs in every variant)
- Strong dependence on a single year that can't be defensibly excluded

---

## 7. Open questions

These need attention but can be deferred to implementation time:

1. **What's the actual format of the macroinvertebrate data?** Counts,
   category abundance, presence/absence only? This determines whether
   Jaccard or Bray-Curtis partitioning is appropriate.

2. **What's the empirical detection frequency at the reference site for
   core families?** If typically ≥ 0.8 per family per year, simple
   Baselga partitioning on core families is probably fine. If lower,
   formal occupancy modelling may be needed.

3. **Are sampling protocols consistent across years?** If protocols
   changed mid-record, the calibration may need to treat pre/post-change
   data separately, which kills sample size.

4. **Are there known differential disturbances at either site that need
   flagging?** Riparian works, channel modification, major hydrological
   events.

5. **How should the time-aggregation be chosen?** Empirical (which
   correlates best with PAFobs) is circular. Mechanistic (which timescale
   matches biological integration) requires assumptions. Pre-registration
   of the aggregation choice before looking at the data would be ideal.

6. **Is J_turnover actually stable in time?** This is testable by Stage
   3.1. If not, the conceptual model needs revision and the
   regression-intercept approach to estimating geomorphology is
   incomplete.

7. **Does the orthogonality assumption hold?** Stage 3.5 tests it. If
   extreme-hydrology years show elevated J_nestedness independent of
   AmsPAF and LMF, hydrology and leachate aren't fully orthogonal at
   the community level.

8. **What's the AmsPAF-LMF correlation in this dataset?** Stage 1.2 and
   3.4. If very high (r > 0.9), the multi-predictor variants struggle
   with multicollinearity and the decomposition can't be resolved.

9. **Should the reference site itself be inspected for AmsPAF
   exceedances?** Even at "clean" sites, the ARA shift may have produced
   non-zero AmsPAF. If so, the design should compare actual reference
   AmsPAF to downstream AmsPAF, not assume reference = 0.

10. **Does LMF behave well as a community-impact predictor?** It was
    designed as a chemistry indicator, not an ecological one. Its
    relationship to community impact may be non-monotonic (modest LMF
    causing organic enrichment that *boosts* productivity, while high
    LMF crashes communities). Worth checking in the scatter plots.

11. **What's the prediction interval interpretation in regulatory or
    internal reporting?** Calibrated AmsPAF with a wide CI may be less
    useful than uncalibrated AmsPAF if the CI swamps decision thresholds.
    Worth thinking through how the output will actually be used.

12. **What's the right τ for the exponential decay?** 90 days is the
    default suggestion based on macroinvertebrate generation times, but
    the literature doesn't strongly motivate a specific value.
    Sensitivity testing across τ ∈ {30, ..., 180} is required. If
    results depend strongly on τ, the chronic integration is fitting
    something specific and may not generalise well to future samples.

13. **Should pulse-aware duration weighting be implemented?** The
    forward-step rule is the default. A more sophisticated approach
    (detect pulses empirically, use forward-step between baselines but
    trapezoidal across pulse onsets) might be worth implementing if
    sensitivity testing shows the simple rule biases results. Defer
    until evidence warrants.

14. **Is per-sample AmsPAF still useful?** Yes, for chemistry-sample-
    level reporting and triggering management responses. The chronic
    integration is for the macroinvertebrate calibration only. Make
    sure both outputs are retained in the pipeline.

15. **What's the right inclusion threshold k?** 5% minimum detection
    frequency is a reasonable default but the choice is judgement.
    Higher k (10–20%) excludes more marginal analytes; lower k (1–5%)
    includes more. The Bayesian framework handles low-detection analytes
    appropriately but they contribute little signal. Sensitivity-test
    k ∈ {1%, 5%, 10%, 20%} and see how msPAF changes.

16. **Are pesticides genuinely all-BDL at this site?** Worth confirming
    by inspecting `chemistry_per_sample` before running the pipeline.
    If any pesticide has been detected even once, the inclusion criterion
    keeps that pesticide in the Bayesian framework and its BDL values
    are imputed via cross-analyte correlations. If genuinely all-BDL
    across all pesticides over the analysis window, the class is
    excluded and msPAF is reported as "metals + ammonia + H2S, with
    no detected pesticide contribution".

17. **How should excluded analytes be documented?** The output should
    clearly state which analytes were considered, which were included
    in msPAF, and which were excluded by the detection-frequency
    criterion. This makes the index honest about its scope.

---

## 8. Reference list for when you get back to this

Core methodology:

- Baselga A (2010) Partitioning the turnover and nestedness components
  of beta diversity. *Global Ecology and Biogeography* 19(1): 134-143.
  Original presence/absence partitioning.

- Baselga A (2012) The relationship between species replacement,
  dissimilarity derived from nestedness, and nestedness. *Global Ecology
  and Biogeography* 21(12): 1223-1232. Clarification and refinement.

- Baselga A (2013) Separating the two components of abundance-based
  dissimilarity: balanced changes in abundance vs. abundance gradients.
  *Methods in Ecology and Evolution* 4(6): 552-557. Bray-Curtis
  partitioning into balanced + gradient components.

- Baselga A, Orme CDL (2012) betapart: an R package for the study of beta
  diversity. *Methods in Ecology and Evolution* 3(5): 808-812.
  R package implementation.

- Chessman BC (2003) New sensitivity grades for Australian river
  macroinvertebrates. *Marine and Freshwater Research* 54(2): 95-103.
  SIGNAL2 methodology (for Framing A secondary analysis).

- De Zwart D, Posthuma L (2005) Complex mixture toxicity for single and
  multiple species: proposed methodologies. *Environmental Toxicology and
  Chemistry* 24(10): 2665-2676. msPAF framework.

- Posthuma L, Suter II GW, Traas TP (eds.) (2001) *Species Sensitivity
  Distributions in Ecotoxicology*. CRC Press. Background on SSDs and
  community-level effects estimation.

Imperfect detection and occupancy (if needed):

- MacKenzie DI et al. (2002) Estimating site occupancy rates when
  detection probabilities are less than one. *Ecology* 83(8): 2248-2255.

- Royle JA, Dorazio RM (2008) *Hierarchical Modeling and Inference in
  Ecology*. Academic Press.

- Chao A, Jost L (2012) Coverage-based rarefaction and extrapolation:
  standardizing samples by completeness rather than size. *Ecology*
  93(12): 2533-2547.

R packages:

- `betapart` — Baselga beta-diversity partitioning (primary)
- `vegan` — community ecology including Bray-Curtis, ordination
- `unmarked` — frequentist occupancy modelling (if escalating)
- `spOccupancy` — Bayesian multi-species occupancy (if escalating)
- `iNEXT` — rarefaction and extrapolation

Community tolerance / PICT:

- Blanck H, Wallin G, Waengberg SA (1988) Species-dependent variation
  in algal sensitivity to chemical compounds. *Ecotoxicology and
  Environmental Safety* 8(4): 339-351. Foundation of pollution-induced
  community tolerance theory.

---

## 9. Decision tree for when you return

```
Start
├── Is the macroinvertebrate data in usable format?
│   └── No  → format conversion first
│
├── Are sampling protocols consistent across years?
│   └── No  → split data, calibration may not be feasible
│
├── Stage -1: Pre-screen analytes
│   └── Include only analytes quantified above DL in >= k% of samples
│       (k = 5% default). Exclude all-BDL or near-all-BDL analytes
│       (typically pesticides at non-industrial landfills) and document.
│
├── Stage 0: Chronic integration of AmsPAF and LMF
│   └── 365-day window, exp decay tau=90, forward-step duration
│       weighting, time-weighted geometric mean. BDL handled by upstream
│       Bayesian multivariate GAM via left-censored likelihood.
│       Output: chronic AmsPAF and chronic LMF per biological sampling date
│
├── Stage 1: Data assembly (include chronic predictors alongside biology)
│   └── Output: tidy macroinverts + chronic AmsPAF + chronic LMF
│
├── Stage 2: Baselga partitioning per year
│   └── Output: J_nestedness, J_turnover, J_total time series
│
├── Stage 3: Diagnostic plots
│   ├── J_turnover roughly constant in time?
│   │   ├── Yes → conceptual model supported
│   │   └── No  → investigate; geomorphology may not be constant
│   │
│   ├── J_nestedness varies in time?
│   │   ├── Yes → proceed
│   │   └── No  → calibration not feasible; investigate why
│   │
│   ├── Chronic AmsPAF-LMF correlation?
│   │   ├── r > 0.9         → multi-predictor variants struggle;
│   │   │                     rely on single-predictor or PCA decomposition
│   │   ├── 0.5 < r < 0.9   → multi-predictor variants viable
│   │   └── r < 0.5         → multi-predictor variants most informative
│   │
│   ├── Orthogonality check (extreme hydrology years)?
│   │   ├── Holds → proceed with simple decomposition
│   │   └── Fails → may need to model hydrology effect explicitly
│   │
│   └── Sensitivity to PAFobs definition?
│       ├── Robust → use simpler definition
│       └── Not robust → think harder about target metric
│
├── Stage 4: Fit V1, V2, V3, V4 in parallel
│
├── Stage 5: Compare variants
│   ├── V1 ≈ V2 (β_tox)         → AmsPAF and LMF independent;
│   │                              either single-predictor variant works
│   ├── V1 ≫ V2 (β_tox)         → much of V1 signal was actually LMF;
│   │                              prefer V2 for mechanism clarity
│   ├── V2 β_tox ≈ 0, V3 γ > 0  → AmsPAF doesn't contribute beyond LMF;
│   │                              negative result for AmsPAF mechanism
│   ├── V2 β_wq ≈ 0             → LMF doesn't add; prefer V1 (parsimony)
│   └── All similar with wide CIs → can't separate; report ambiguity
│
├── Stage 6: Sensitivity tests for chosen variant
│
├── Stage 7: Interpret α and β coefficients
│   ├── α ≈ 0, β > 0    → clean calibration, geomorphology negligible
│   ├── α > 0, β > 0    → stable geomorphological offset of magnitude α;
│   │                     AmsPAF/LMF tracks leachate-driven impact via β
│   ├── α > 0, β ≈ 0    → predictors don't track variation; calibration fails
│   ├── α < 0           → unexpected; investigate
│   └── pathological    → abandon calibration; report AmsPAF as-is with
│                         a note about validation failure
│
└── Stage 8: Document and integrate into reporting
```

---

*This is a planning document. When you return to this work, start by
reviewing Open Questions in §7 against your current data state, then
follow the analysis plan in §4. Don't fit a model before inspecting the
plots in Stage 3.*
