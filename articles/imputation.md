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

2.  **Residual correlation across analytes (`rescor = TRUE`).** Metals
    are modelled jointly. In a leachate- or AMD-impacted aquifer,
    conservative tracers move ahead of redox-controlled metals, so a
    post-pulse sample can look near-baseline on its major ions while
    Cu/Pb/Zn/Mn stay co-elevated. That co-elevation is pure residual
    correlation with no predictor driving it — modelling the residual
    correlation matrix lets an *observed* metal at a sample inform the
    *missing* ones.

3.  **Hurdles, so silence is not mistaken for absence.** A sample is
    only imputed for a group if it carries at least one analyte from
    that group — one metal for the metals model, one of {DOC, TOC, BOD,
    COD, cBOD} for the organics model. A sample with no metals recorded
    is left alone, because a leachate metal pulse may simply not have
    reached that location yet, and inventing values there would
    manufacture an impact that the data do not support.

The work is split across three functions:

| Function | Role | Engine |
|----|----|----|
| [`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md) | Fit the metals/organics model(s) once on training data | `brms` (Stan) |
| [`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md) | Fill BDL and missing metals/organics on new data | `brms` posterior |
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
discovers the analyte groups automatically: anything matching the
built-in metals list becomes the metals model, the remaining
non-excluded concentration analytes become the organics model, and
microbiological counts and qualitative descriptors are excluded
outright. Fitting runs MCMC through Stan, so it takes minutes, not
seconds — do it once and reuse the result.

``` r

model <- fit_imputation_model(
  monitoring_long,
  required_vars = c("pH", "EC"),
  save_dir      = "models"   # optional: persist the fit as a .qs file
)

model
#> <imputation_model>  fitted 2026-06-03 | 248 samples | 14 PCA vars | 3 PCs (78% var)
#>   metals:   8 analytes
#>   organics: 2 analytes
```

The PCA axis count is chosen to reach a target cumulative variance
(`min_var_explained`, default 75%) up to `max_pcs` (default 6), with a
floor of two axes. The console reports which variables entered the PCA,
how many axes were retained, and the analytes in each group, so you can
confirm the model saw what you expected.

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
estimate.

### A note on the BDL cap

Combining `rescor = TRUE` with proper left-censoring is not possible in
`brms`, so BDL cells are modelled as missing (`mi()`) and then
**capped** at their original detection limit when the posterior mean
lands above it (`bdl_cap = TRUE`, the default). The cap keeps BDL
imputations physically consistent with what the laboratory reported.
Where the surrounding chemistry genuinely points to a high concentration
the cap can fire often;
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
warns when it does, naming the affected analytes, so those samples can
be inspected rather than trusted blindly.

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
