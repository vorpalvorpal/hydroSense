# Detecting leachate with the Leachate Mixing Fraction (LMF)

## What the LMF answers

[`add_lmf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_lmf.md)
asks a single question of each water sample:

> *What fraction of this sample’s chemistry can be explained as a
> mixture of local reference water and landfill leachate?*

The answer — the **leachate mixing fraction (LMF)** — is reported as a
percentage: `0` means the sample is chemically indistinguishable from
local reference water, `100` means it matches the leachate end-member.
It is a *source-detection* index, deliberately separate from toxicity: a
sample can carry a clear leachate signature (high LMF) while still being
well below guideline values, and vice versa.

LMF is built on the classic **two-component end-member mixing** model
used in catchment hydrochemistry (Christophersen & Hooper, 1992).
Leachate-impacted freshwater is treated as a blend of two sources, and
the major-ion signature is used to estimate the blend ratio.

## The method in brief

For each sample,
[`add_lmf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_lmf.md):

1.  **Converts to milliequivalents.** Concentrations are converted to
    meq/L with \[to_meq()\] so that ions combine on a charge-equivalent
    basis. Transforming species are collapsed to conserved totals — the
    nitrogen species (NH₃-N, NO₃-N, NO₂-N) to *total N*, and the
    carbonate species to *total alkalinity* — so that redox or
    speciation shifts between source and sample do not break the mixing
    assumption.

2.  **Builds two end-members.** The **reference** end-member is the mean
    and standard deviation of each ion in local background water. The
    **leachate** end-member is anchored to chloride (a conservative
    tracer) and expressed as Cl-scaled ratios.

3.  **Computes a per-ion mixing fraction** for each measured ion *i*:
    ``` math
     f_i = \frac{x_i - R_i}{L_i - R_i} 
    ```
    where $`x_i`$ is the sample, $`R_i`$ the reference, and $`L_i`$ the
    leachate end-member.

4.  **Combines the ions** into one estimate by an inverse-variance
    weighted mean, so ions measured precisely and with a large
    leachate–reference gradient count for more. An **informativeness**
    screen (the reference scatter relative to the leachate–reference
    gap) admits a sample only if enough high-information ions are
    present.

5.  **Down-weights disagreeing ions** with a robust (Huber) reweighting
    pass, so a single ion behaving non-conservatively (sulfate
    reduction, an extra source) does not drag the estimate. A χ²
    diagnostic flags overall ion disagreement.

## A worked example

The package ships a small synthetic dataset, `leachate_demo`, with three
sites: an impacted `downstream` site, a clean `reference` site, and the
`leachate` end-member. (It is fictional — generated from a mixing model
so the example is reproducible.)

``` r

demo <- leachate_demo()
unique(demo$site_id)
#> [1] "downstream" "reference"  "leachate"
head(subset(demo, site_id == "downstream"), 4)
#>   sample_id    site_id   datetime analyte    value detected units.analyte
#> 1     DS-01 downstream 2024-01-15      Na 149.4775     TRUE          mg/L
#> 2     DS-01 downstream 2024-01-15       K  29.9655     TRUE          mg/L
#> 3     DS-01 downstream 2024-01-15      Ca  75.4548     TRUE          mg/L
#> 4     DS-01 downstream 2024-01-15      Mg  25.9333     TRUE          mg/L
#>   valence.analyte atomic_mass.analyte
#> 1               1               22.99
#> 2               1               39.10
#> 3               2               40.08
#> 4               2               24.31
```

Supplying the leachate and reference chemistry directly uses the
caller-provided end-member path (no site-matching infrastructure
required):

``` r

out <- add_lmf(
  df             = subset(demo, site_id == "downstream"),
  leachate_data  = subset(demo, site_id == "leachate"),
  reference_data = subset(demo, site_id == "reference")
)

lmf <- subset(out, analyte == "LMF",
              c("sample_id", "value", "sigma_lmf", "n_ions_used", "chi2_per_df"))
lmf
#> # A tibble: 6 × 5
#>   sample_id value sigma_lmf n_ions_used chi2_per_df
#>   <chr>     <dbl>     <dbl>       <int>       <dbl>
#> 1 DS-01      14.9     0.472           8       0.163
#> 2 DS-02      14.7     0.467           8       0.208
#> 3 DS-03      14.7     0.465           8       0.155
#> 4 DS-04      15.1     0.476           8       0.129
#> 5 DS-05      14.6     0.463           8       0.201
#> 6 DS-06      14.9     0.471           8       0.268
```

The downstream samples return an LMF of about 15% — consistent with how
the demo data were constructed (a 15% leachate blend). `sigma_lmf` is
the uncertainty on the estimate (percentage points), `n_ions_used` the
number of ions that contributed, and `chi2_per_df` the ion-agreement
diagnostic (values near 1 indicate the ions tell a consistent mixing
story).

## Interpreting the result

LMF is a continuous index: the raw percentage is the primary output. The
following tier breaks are a suggested starting point for interpretation:

| Tier | LMF (%) | Suggested interpretation                            |
|------|---------|-----------------------------------------------------|
| 1    | ≤ 1     | Background — indistinguishable from local reference |
| 2    | 1–5     | Trace impact — detectable signature, monitor trend  |
| 3    | 5–20    | Significant impact — clear signature, investigate   |
| 4    | \> 20   | Severe impact — urgent investigation                |

> **Note:** These breaks are **starting points**, not regulatory
> thresholds. Classification decisions should be made by the reporting
> system or analyst, not imposed automatically. The tier 1/2 boundary in
> particular should be reviewed against the LMF distribution of your own
> reference samples after deployment — ideally it should sit near the
> 95th percentile of LMF on known-clean water.

## When LMF returns `NA`

A sample is not scored (its `lmf_reason` explains why) when it lacks
enough high-information ions, when the estimate is too imprecise
(`sigma_lmf` above the `max_sigma_lsi` threshold), or when no ion has a
usable leachate–reference gradient. These are honest non-answers: the
major-ion panel measured at that sample simply cannot distinguish the
two sources.

## References

- Christophersen N, Hooper RP (1992). Multivariate analysis of stream
  water chemical data. *Water Resources Research* 28(1):99–107.
- Christensen JB et al. (2001). Biogeochemistry of landfill leachate
  plumes. *Applied Geochemistry* 16(7–8):659–718.
- Kjeldsen P et al. (2002). Present and long-term composition of MSW
  landfill leachate. *Critical Reviews in Environmental Science and
  Technology* 32(4).
