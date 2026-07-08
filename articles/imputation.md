# Imputing missing and below-detection chemistry

## The problem

Routine monitoring chemistry is rarely complete. Two kinds of gaps
dominate:

- **Below detection limit (BDL).** The analyte was measured but the
  result fell under the laboratory’s reporting limit. All we know is
  that the true value lies somewhere in `[0, DL]`. Substituting `DL`,
  `DL/2`, or `0` — the usual quick fixes — biases every downstream
  summary, and the direction of the bias depends on which fix you
  picked.
- **Missing.** The analyte was not measured at that sample at all, so
  there is no value and no detection limit.

Both gaps matter for mixture toxicity:
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
sums the potentially affected fraction across analytes, so a sample
where Cu was never measured and Zn came back BDL is not automatically
“clean” — it is *under-determined*. The imputation tools fill these gaps
with statistically defensible estimates that borrow strength from the
analytes that *were* measured.

## The approach in brief

hydroSense imputes chemistry with a **Bayesian multivariate model**
fitted once on your monitoring history and then applied to new samples.
The design has three ideas at its core:

1.  **A shared chemistry summary as predictors.** All the routinely
    measured water-quality variables — major ions, pH, EC, NH₃-N, DOC,
    nutrients, redox indicators — are compressed into a handful of
    principal components with
    [`nipals::nipals()`](https://kwstat.github.io/nipals/reference/nipals.html),
    which tolerates within-sample missing cells natively. These PC
    scores become the predictors, so a missing metal is estimated from
    the *whole* measured chemical context, not one or two surrogate
    analytes.

2.  **Cross-analyte coupling (optional).** Metals *can* be modelled
    jointly, so an *observed* metal at a sample informs the *missing*
    ones — the idea being that in a leachate- or AMD-impacted aquifer,
    conservative tracers move ahead of redox-controlled metals, leaving
    Cu/Pb/Zn/Mn co-elevated beyond what the major ions predict. In
    practice, on the B.S01 panel this borrowing did **not** help (the
    cross-metal correlation is too weak and ragged, and borrowing noise
    hurts sparse analytes), so the **default method does no borrowing**
    and predicts each metal from the shared chemistry PCA alone. Opt
    into coupling with `impute_method = "factor"`; see *Choosing the
    imputation method*.

3.  **Hurdles, so silence is not mistaken for absence.** A sample is
    only imputed for a group if it carries at least one of that group’s
    *hurdle* analytes — one metal for the metals model, one of {DOC,
    TOC, BOD, COD, cBOD} for the organics model. A sample with no metals
    recorded is left alone, because a leachate metal pulse may simply
    not have reached that location yet, and inventing values there would
    manufacture an impact that the data do not support.

The engine itself is **domain-agnostic**: it imputes an arbitrary set of
analyte *groups* (each a joint model with its own optional presence
hurdle), described with
[`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md).
The metals/organics structure above is just the leachate preset,
[`leachate_impute_groups()`](https://vorpalvorpal.github.io/hydroSense/reference/leachate_impute_groups.md),
supplied as the default — swap in your own groups to impute a different
chemistry (see *Imputing in another domain*).

The work is split across three functions:

| Function | Role | Engine |
|----|----|----|
| [`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md) | Fit one model per group once on training data | `brms` (Stan) |
| [`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md) / [`leachate_impute_groups()`](https://vorpalvorpal.github.io/hydroSense/reference/leachate_impute_groups.md) | Describe the groups to impute (or use the leachate preset) | — |
| [`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md) | Fill BDL and missing group analytes on new data | `brms` posterior |
| [`impute_coanalytes()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_coanalytes.md) | Fill missing normalisation co-analytes (DOC, Ca, Mg, hardness) | [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) |

## The data contract

Every function takes the same long-format frame used throughout the
package: one row per sample × analyte, with `sample_id`, `site_id`,
`datetime`, `analyte`, `value`, and `detected`. A BDL result is a row
where `detected` is `FALSE` and `value` holds the detection limit; a
missing result is simply the absence of a row.

``` r

chem <- tibble::tribble(
  ~sample_id, ~analyte, ~value, ~detected,
  "s1",       "pH",       7.2,   TRUE,
  "s1",       "EC",      410,    TRUE,
  "s1",       "Cu",        3.1,  TRUE,
  "s1",       "Zn",        0.5,  FALSE,   # BDL: value = detection limit
  "s1",       "DOC",       2.4,  TRUE
  # Ni simply has no row for s1 -> "missing"
)
chem
#> # A tibble: 5 × 4
#>   sample_id analyte value detected
#>   <chr>     <chr>   <dbl> <lgl>   
#> 1 s1        pH        7.2 TRUE    
#> 2 s1        EC      410   TRUE    
#> 3 s1        Cu        3.1 TRUE    
#> 4 s1        Zn        0.5 FALSE   
#> 5 s1        DOC       2.4 TRUE
```

`pH` and `EC` are **required variables**: a sample missing either is
dropped before fitting, because they anchor the chemistry PCA. When a
required variable is itself BDL, its detection-limit value is used as a
conservative stand-in and the sample is kept.

## Fitting the model

[`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md)
assigns analytes to the groups it was given (the leachate preset by
default): the `metals` group claims anything in the built-in metals
list, the catch-all `organics` group takes the remaining non-excluded
concentration analytes, and microbiological counts and qualitative
descriptors are excluded outright (`exclude =`). Fitting runs MCMC
through Stan, so it takes minutes, not seconds — do it once and reuse
the result.

``` r

model <- fit_imputation_model(
  monitoring_long,
  required_vars = c("pH", "EC"),
  save_dir      = "models"   # optional: persist the fit as a .qs file
)

model
#> <imputation_model>  fitted 2026-06-03 | 248 samples | 14 PCA vars | 3 PCs (78% var)
#>   method:   marginal
#>   metals:   8 analytes
#>   organics: 2 analytes
```

The PCA axis count is chosen to reach a target cumulative variance
(`min_var_explained`, default 75%) up to `max_pcs` (default 6), with a
floor of two axes. The console reports which variables entered the PCA,
how many axes were retained, and the analytes in each group, so you can
confirm the model saw what you expected.

A target analyte is only modelled if it is *detected* in at least
`min_target_detect_freq` of samples (default 5%). Near- or
all-below-detection analytes carry no signal to impute from and would
otherwise inflate the model — on a panel with a hundred
mostly-undetected organics this filter is what keeps the fit tractable.

## Imputing in another domain

The metals/organics split is a *preset*, not a hard-coded assumption. To
impute a different chemistry, pass your own list of
[`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md)
objects. Each group is fitted as its own joint model (analytes in the
same group borrow residual correlation from one another; separate groups
do not), and each may carry a presence hurdle. A group with
`targets = NULL` is the catch-all — it claims every remaining modellable
analyte not already taken by an earlier group.

``` r

model <- fit_imputation_model(
  monitoring_long,
  groups = list(
    impute_group("trace_metals", targets = c("As", "Cd", "Pb"),
                 hurdle = c("As", "Cd", "Pb")),
    impute_group("nutrients",    targets = NULL,            # catch-all
                 hurdle = c("NO3-N", "P-total"))
  ),
  exclude = c("Coliforms", "Turbidity")   # never modelled
)
```

[`leachate_impute_groups()`](https://vorpalvorpal.github.io/hydroSense/reference/leachate_impute_groups.md)
returns the default preset, so you can inspect or extend it rather than
starting from scratch.

## Choosing the imputation method

`fit_imputation_model(impute_method =)` controls how below-detection
cells and cross-analyte coupling are handled. The methods split on one
design question: **should an observed metal at a sample help predict the
missing/BDL ones at that same sample?** All share the same PCA
predictors and per-analyte spline mean; they differ in whether (and how)
they add cross-analyte coupling on top.

| `impute_method` | Engine | BDL handling | Cross-analyte coupling | Uncertainty |
|----|----|----|----|----|
| `"marginal"` *(default)* | mgcv | `cens` (truncated draws) | **none** | proper posterior predictive (t) |
| `"rescor_mi"` | brms/Stan | `mi()` + DL cap | full residual-correlation matrix | over-wide; hard to sample |
| `"cens"` | brms/Stan | `cens("left")` | none | calibrated, slow |
| `"cens_factor"` | brms/Stan | `cens("left")` | shared per-sample latent factor | calibrated, cheap |
| `"factor"` | mgcv + cmdstanr | `cens` (bounded latents) | low-rank latent factor (inferred at prediction) | calibrated |

### The benchmark, and why `"marginal"` is the default

A hold-out benchmark on real B.S01 data — masking ~10% of detected
routine-metal cells, pooled over three mask seeds (57 held-out cells),
scored on the log₁₀ scale — gives:

| `impute_method` |     RMSE |      MAE |  bias |    90% coverage | 90% width |  fit time |
|-----------------|---------:|---------:|------:|----------------:|----------:|----------:|
| `"marginal"`    |     0.30 | **0.19** | −0.14 |        **0.91** |  **0.77** | **0.3 s** |
| `"rescor_mi"`   | **0.25** |     0.21 | −0.05 | 1.00 (too wide) |      1.29 |    ~160 s |
| `"factor"`      |     0.34 |     0.23 | −0.21 |            0.98 |      1.00 |     ~60 s |

The headline is a **negative result about borrowing**: on B.S01, letting
a sample’s observed metals condition its missing ones does **not** help.
The `"factor"` method — a principled low-rank latent-factor model that
does exactly that borrowing — is the *worst* performer. A cell-level
breakdown shows why: the loss is concentrated in the **sparse analytes**
(Pb, Cr), whose factor loadings are unidentified, so conditioning on the
sample’s other metals through a spurious loading yanks the prediction
far down (on the worst Pb cell, ~30× too low). For the well-measured
metals it merely equals `"marginal"`. This matches the weak, ragged
cross-metal correlation seen on this panel: there is little real signal
to borrow, and borrowing noise actively hurts where data is thin.

`"marginal"` — which deliberately does **no** borrowing (a per-analyte
censored GAM) — is therefore the robust choice: it has the best MAE,
well-calibrated 90% coverage (0.91), tight informative intervals, needs
neither brms nor cmdstanr, and is ~400× faster than `"rescor_mi"`.
`"rescor_mi"`’s slim RMSE edge (0.25 vs 0.30) comes with uninformative
100%-coverage intervals (width 1.29 vs 0.77) and a posterior that
saturates tree depth and trips E-BFMI.

### `"marginal"`: assumptions and trade-offs

Per target analyte, it fits `log(conc) ~ s(PC1) + … + s(PCk)` on the
*detected* observations (an `mgcv` GAM), then predicts each BDL/missing
cell from a proper **posterior predictive**: GAM parameter uncertainty
(`β ~ N(coef, Vp)`, unconditional so smoothing-parameter uncertainty is
included) **plus** residual-variance uncertainty (`σ²` drawn from its
scaled-inverse-χ² posterior, making the predictive Student-t), with BDL
cells’ draws truncated at the detection limit. Its assumptions:

- **The chemistry PCA carries the signal.** A metal is predicted only
  from the water-quality context (major ions, pH, EC, nutrients, redox
  indicators) summarised by the PC scores — *not* from the other metals.
  If your panel has strong, well-observed cross-metal coupling that the
  PCA does not already capture, `"marginal"` leaves that on the table
  (use `"factor"`).
- **Log-normal, smoothly-varying concentrations.** The spline mean
  assumes the log concentration is a smooth function of the PC scores;
  the residual is Gaussian on the log scale.
- **Detected-only mean.** The GAM mean is fit on quantified values; on a
  heavily censored analyte the mean can sit slightly high (a small
  negative bias in the masked-detect benchmark, −0.14), because the fit
  never sees the low tail. In genuine BDL cells — where the truth really
  is sub-limit — the DL truncation is the correct conservatism.

### When to use the others

- **`"factor"`** — the principled cross-analyte model. Use it when you
  have evidence of genuine, reasonably-observed metal-to-metal coupling
  that the WQ PCA misses; it will beat `"marginal"` only when there is
  real signal to borrow. On weak/ragged panels it mis-conditions sparse
  analytes (see above).
- **`"rescor_mi"` / `"cens"` / `"cens_factor"`** — the brms methods,
  retained for comparison and continuity. `"cens_factor"` is the
  best-calibrated of the three; all are far slower than `"marginal"` and
  need a Stan toolchain.

**How much does the choice matter downstream?** Usually little. Running
the full B.S01 daily-msPAF pipeline (`mspaf_daily`) under different
methods, everything else fixed, the central estimates differ by a
day-to-day median of ~0.005 percentage points: on a panel site the
metals are usually measured, so imputation fills only a handful of cells
and the gap-uncertainty bracket dwarfs the imputation effect. The method
choice matters most for the **width** of the reported bands — where
`"marginal"`’s calibrated, tighter intervals are the honest default.

``` r

# Default — fast, calibrated, no Stan toolchain needed:
model <- fit_imputation_model(monitoring_long)                       # "marginal"

# Opt in to cross-analyte borrowing (only worth it with genuine coupling):
model <- fit_imputation_model(monitoring_long, impute_method = "factor")
```

## Imputing new chemistry

[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)
applies the fitted model, returning the input frame with BDL and missing
metals/organics cells replaced by posterior estimates and two new
columns:

- `imputed` — `TRUE` for any cell that was filled.
- `imputed_kind` — `"observed"`, `"censored_left"` (was BDL), or
  `"missing"` (had no row).

``` r

filled <- impute_chemistry(monitoring_long, model)

# Inspect what was filled, and how
dplyr::count(filled, analyte, imputed_kind)
```

By default the result is a **posterior mean per cell**
(`return = "point"`). For uncertainty propagation, request the full
posterior with `return = "draws"` to get one row per sample × analyte ×
draw — feed those draws through
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
to obtain a posterior distribution of mixture PAF rather than a point
estimate. Draws can be memory-heavy (especially for `"rescor_mi"`);
`ndraws =` subsamples the posterior and `batch_size =` predicts the
samples in batches, both of which bound peak memory.

### The detection-limit cap

An imputed below-detection value must never exceed its detection limit.
The default `"marginal"` method (and `"factor"`) **truncate** BDL draws
at the limit by construction, so no cap is needed and the check is
skipped. The brms methods (`"rescor_mi"`, `"cens"`, `"cens_factor"`)
emit an unconstrained prediction, so for those
[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)
**caps** any imputed BDL cell at its original detection limit
(`bdl_cap = TRUE`, the default) and warns when the cap fires, with a
per-analyte breakdown (how many cells, and the worst exceedance ratio)
so the activations are auditable rather than just an aggregate count.
Frequent cap firing signals tension between the modelled chemistry and
the reported limits — and is characteristic of `"rescor_mi"`, whose
`mi()` predictions are not held below the bound.

For a full per-cell audit, the result carries a `"bdl_cap_summary"`
attribute — one row per capped (`sample_id`, `analyte`) cell with its
detection limit, the uncapped estimate, and the exceedance ratio.
Retrieve it with
[`bdl_cap_summary()`](https://vorpalvorpal.github.io/hydroSense/reference/bdl_cap_summary.md)
on the frame as returned (plain attributes are dropped by most dplyr
verbs, so read it before further wrangling):

``` r

filled <- impute_chemistry(monitoring_long, model)
bdl_cap_summary(filled)   # NULL (with a message) if nothing was capped
```

## Imputing normalisation co-analytes

The bioavailability normalisations in
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
need co-analytes — DOC, Ca, Mg, and hardness — to convert field
concentrations to the SSD index condition (see the *Analyte
normalisation* article). When one of these is missing,
[`impute_coanalytes()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_coanalytes.md)
fills it with a fast, deterministic log-Gaussian GAM
([`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html)) on the *same* PC
scores the metals model used, so the co-analyte estimate is conditioned
on the same chemical context.

``` r

filled <- filled |>
  impute_coanalytes(model)        # fills DOC / Ca / Mg / hardness where absent
```

Two deliberate boundaries:

- It only fills analytes that are **entirely absent** for a sample;
  existing BDL observations are left untouched.
- Imputed co-analyte values are **never** fed back into the
  metals/organics model — that model was already fitted on measured
  values. Co-analyte imputation is a downstream convenience for
  normalisation, not part of the Bayesian fit.

## Where it fits in the pipeline

The imputation steps slot in before mixture toxicity and temporal
aggregation:

``` r

model  <- fit_imputation_model(monitoring_long)

result <- monitoring_long |>
  impute_chemistry(model) |>     # BDL + missing metals/organics
  impute_coanalytes(model) |>    # missing DOC / Ca / Mg / hardness
  add_mspaf(reference = prep_ref) |>
  time_weighted_aggregate()
```

Fit once on the longest, richest history you have; impute repeatedly as
new data arrive.

### One call: `mspaf_pipeline()`

The daily version of that whole flow — fit an imputation model, fit the
reference/impact model, then compute the daily msPAF for a target site —
is bundled into
[`mspaf_pipeline()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_pipeline.md),
with imputation **on by default**:

``` r

out <- mspaf_pipeline(
  target     = subset(monitoring_long, site_id == "downstream"),
  reference  = subset(monitoring_long, site_id == "reference"),
  hydro      = rainfall,
  daily_args = list(require_temperature = FALSE, ndraws = 50L, seed = 1L)
)
```

A single imputation model is fitted (on the reference chemistry by
default; `impute_on = "target"` to fit on the target instead) and
threaded into both
[`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)
and
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md).
Extra arguments go through `reference_args` and `daily_args`; the fitted
reference and imputation models are attached to the result as the
`"reference_model"` and `"imputation_model"` attributes.

Because it imputes by default, the orchestrator **moves the msPAF
result** relative to the non-imputed calculation: imputing
below-detection cells (replacing detection-limit values with modelled
sub-limit estimates) and fabricating entirely-absent analytes both
change the mixture, and in practice the net effect can be sizeable. That
is a missingness (MAR/MNAR) judgement, not a neutral default — pass
`impute = FALSE` to reproduce the non-imputed behaviour of chaining
[`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)
and
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
directly. See *When not to impute* below.

## When not to impute

Imputation borrows strength across analytes — it does not conjure
information that the data never contained. A few situations call for
restraint:

- **Sparse or idiosyncratic chemistry.** If a sample shares little
  measured context with the training data, its PC scores are poorly
  determined and the posterior is wide. The hurdles guard the worst
  cases, but a sample that barely clears a hurdle still deserves a look.
- **Near-all-below-detection organics.** In landfill/WMF leachate most
  organic analytes are below detection in almost every sample, so there
  is no signal to borrow from and no reliable covariate to anchor them —
  DOC fixes the bulk carbon, not a congener-specific ratio. The package
  deliberately does **not** impute these: the DOC-like hurdle skips the
  organics model without a carbon covariate, and
  `min_target_detect_freq` (default 5%) drops any analyte detected too
  rarely to model. Such organics are left as reported non-detections
  rather than given manufactured values; treat them as a censored “≤
  detection limit” reading in any downstream summary.
- **Frequent BDL-cap warnings for an analyte.** This signals tension
  between the modelled chemistry and the reported detection limits;
  treat those cells as flagged rather than final.
- **Reporting raw detections.** Imputed values are estimates. Keep the
  `imputed` / `imputed_kind` columns with your results so that a reader
  can always separate measured from modelled, and present the posterior
  draws where the downstream decision is sensitive to the gaps.

## Installation note

**brms** is a hard dependency of hydroSense (it is listed under
`Imports`), so it is installed automatically with the package. The
Bayesian step additionally needs a working **Stan** engine, which is
*not* pulled in automatically. Set it up once:

``` r

install.packages("brms")
```

If the install alone is not enough, follow the short Stan toolchain
setup at <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>.
[`impute_coanalytes()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_coanalytes.md)
uses `mgcv` (already a hard dependency) and does **not** require Stan.
