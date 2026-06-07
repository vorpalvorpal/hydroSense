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
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
sums the potentially affected fraction across analytes, so a sample
where Cu was never measured and Zn came back BDL is not automatically
“clean” — it is *under-determined*. The imputation tools fill these gaps
with statistically defensible estimates that borrow strength from the
analytes that *were* measured.

## The approach in brief

leachatetools imputes chemistry with a **Bayesian multivariate model**
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

2.  **Cross-analyte coupling.** Metals are modelled jointly. In a
    leachate- or AMD-impacted aquifer, conservative tracers move ahead
    of redox-controlled metals, so a post-pulse sample can look
    near-baseline on its major ions while Cu/Pb/Zn/Mn stay co-elevated.
    That co-elevation is residual correlation with no predictor driving
    it — capturing it lets an *observed* metal at a sample inform the
    *missing* ones. By default this is a full residual correlation
    matrix (`rescor = TRUE`); see *Choosing the imputation method* for
    the alternatives.

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
[`impute_group()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_group.md).
The metals/organics structure above is just the leachate preset,
[`leachate_impute_groups()`](https://vorpalvorpal.github.io/leachatetools/reference/leachate_impute_groups.md),
supplied as the default — swap in your own groups to impute a different
chemistry (see *Imputing in another domain*).

The work is split across three functions:

| Function | Role | Engine |
|----|----|----|
| [`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md) | Fit one model per group once on training data | `brms` (Stan) |
| [`impute_group()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_group.md) / [`leachate_impute_groups()`](https://vorpalvorpal.github.io/leachatetools/reference/leachate_impute_groups.md) | Describe the groups to impute (or use the leachate preset) | — |
| [`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md) | Fill BDL and missing group analytes on new data | `brms` posterior |
| [`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md) | Fill missing normalisation co-analytes (DOC, Ca, Mg, hardness) | [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) |

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

[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)
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
#>   method:   rescor_mi
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
[`impute_group()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_group.md)
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

[`leachate_impute_groups()`](https://vorpalvorpal.github.io/leachatetools/reference/leachate_impute_groups.md)
returns the default preset, so you can inspect or extend it rather than
starting from scratch.

## Choosing the imputation method

`fit_imputation_model(impute_method =)` controls how below-detection
cells and cross-analyte coupling are handled. All three share the same
PCA predictors; `"rescor_mi"` and `"cens"` fit a wide multivariate
model, while `"cens_factor"` fits a single long-format model with a
shared per-sample factor:

| `impute_method` | BDL handling | Coupling | Notes |
|----|----|----|----|
| `"rescor_mi"` *(default)* | `mi()` (imputed) + DL cap | full residual-correlation matrix | most accurate recovery; posterior **draws are memory-heavy** |
| `"cens"` | `cens("left")` (proper censoring) | none | fastest; cleanest convergence; no borrowing across metals |
| `"cens_factor"` | `cens("left")` | shared per-sample latent factor | proper censoring **and** coupling; calibrated, cheap uncertainty |

The trade-off is accuracy versus the cost of uncertainty. `brms` cannot
combine a full residual correlation (`rescor = TRUE`) with proper
left-censoring, so the default pairs residual correlation with `mi()`
imputation and a post-hoc cap. A hold-out benchmark on real data found
`"rescor_mi"` recovers masked metals most accurately (the full
residual-correlation matrix is the strongest coupling), which is why it
is the default — but its posterior **draws** are expensive and can
exhaust memory.

`"cens"` and `"cens_factor"` give well-calibrated intervals cheaply.
They differ in coupling: `"cens"` treats the metals independently, while
`"cens_factor"` adds a single per-sample latent factor shared across all
analytes — fitted in long format, so the factor is well-identified (each
sample contributes several analyte observations) and genuinely lets an
observed metal inform the unobserved/BDL ones at that sample. In the
B.S01 hold-out this coupling places `"cens_factor"` between `"cens"` and
`"rescor_mi"` on accuracy (RMSE 0.43 vs 0.53 and 0.25) while converging
best of the three (R̂ 1.006). Prefer a cens method when you need
uncertainty (`return = "draws"`) or hit memory limits — `"cens_factor"`
when you also want cross-analyte borrowing, `"cens"` when you want the
simplest, fastest fit; prefer the default `"rescor_mi"` for the most
accurate point estimates.

``` r

# Proper censoring + cross-analyte coupling, with cheap calibrated uncertainty:
model <- fit_imputation_model(monitoring_long, impute_method = "cens_factor")
```

## Imputing new chemistry

[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
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
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
to obtain a posterior distribution of mixture PAF rather than a point
estimate. Draws can be memory-heavy (especially for `"rescor_mi"`);
`ndraws =` subsamples the posterior and `batch_size =` predicts the
samples in batches, both of which bound peak memory.

### The detection-limit cap

An imputed below-detection value must never exceed its detection limit.
The posterior prediction is not itself constrained below the DL, so
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
**caps** any imputed BDL cell at its original detection limit
(`bdl_cap = TRUE`, the default) and warns when the cap fires, with a
per-analyte breakdown (how many cells, and the worst exceedance ratio)
so the activations are auditable rather than just an aggregate count.
Frequent cap firing signals tension between the modelled chemistry and
the reported limits.

For a full per-cell audit, the result carries a `"bdl_cap_summary"`
attribute — one row per capped (`sample_id`, `analyte`) cell with its
detection limit, the uncapped estimate, and the exceedance ratio.
Retrieve it with
[`bdl_cap_summary()`](https://vorpalvorpal.github.io/leachatetools/reference/bdl_cap_summary.md)
on the frame as returned (plain attributes are dropped by most dplyr
verbs, so read it before further wrangling):

``` r

filled <- impute_chemistry(monitoring_long, model)
bdl_cap_summary(filled)   # NULL (with a message) if nothing was capped
```

## Imputing normalisation co-analytes

The bioavailability normalisations in
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
need co-analytes — DOC, Ca, Mg, and hardness — to convert field
concentrations to the SSD index condition (see the *Analyte
normalisation* article). When one of these is missing,
[`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md)
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
  add_amspaf(reference = prep_ref) |>
  time_weighted_aggregate()
```

Fit once on the longest, richest history you have; impute repeatedly as
new data arrive.

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

The Bayesian step is the only part of leachatetools that needs **brms**,
so brms is an optional dependency rather than a hard requirement — the
LMF and AmsPAF tools install and run without it. To use
[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)
and
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md),
install brms and a working Stan engine once:

``` r

install.packages("brms")
```

If the install alone is not enough, follow the short Stan toolchain
setup at <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>.
[`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md)
uses `mgcv` (already a hard dependency) and does **not** require brms or
Stan.
