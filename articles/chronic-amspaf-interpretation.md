# Interpreting chronic AmsPAF outputs

This vignette covers how to interpret the values returned by
[`add_amspaf()`](https://www.kedumba.com.au/leachatetools/reference/add_amspaf.md)
and
[`time_weighted_aggregate()`](https://www.kedumba.com.au/leachatetools/reference/time_weighted_aggregate.md)
in the context of an environmental assessment. The package deliberately
stops short of providing regulatory tier breaks — those depend on the
assessment context and are a consumer decision. This document gives you
the conceptual tools to make those decisions.

## msPAF is not a single-substance trigger value

The “5%” threshold in ANZG/ANZECC guidelines (the 95% species protection
level) is derived from a **single-substance** species sensitivity
distribution (SSD). It says: “if a water body contains only this
substance at this concentration, 5% of species would be affected.”

AmsPAF outputs a number on the same 0–100% scale but with a different
denominator: “given the *mixture* of substances present at observed
concentrations, what fraction of species would be affected by *the
entire mixture*?”

The implications:

1.  A 5% AmsPAF and a 5% single-substance PAF do not correspond to the
    same “level of impact.” A mixture of 10 substances each at their 1%
    PAF contributes a combined AmsPAF substantially larger than any
    single substance’s 1%.

2.  The denominator depends on which substances enter the calculation. A
    site with Cu and Zn measured will have a different AmsPAF than the
    same site with Cu, Zn, and Pb measured — even if Pb is well below
    its single-substance trigger.

3.  **You cannot directly compare an AmsPAF value to an ANZG
    single-substance protection level.** They are different metrics on
    the same scale.

What you *can* do:

- Compare AmsPAF across sites to identify relative impact.
- Compare AmsPAF over time at one site to detect trends.
- Calibrate AmsPAF against biological response data (e.g.
  macroinvertebrate community state) to derive site-relevant
  interpretive thresholds.

## Chronic vs per-sample AmsPAF

[`add_amspaf()`](https://www.kedumba.com.au/leachatetools/reference/add_amspaf.md)
operates on a chemistry data frame regardless of whether each row is a
single field sample or a chronic-aggregated value — it treats every row
as a snapshot of exposure and computes the AmsPAF that snapshot implies.

For chronic interpretation (i.e. comparing AmsPAF against community
state that integrates exposure over weeks to months), there are two
pathways:

**Path A — concentration averaging.** Time-aggregate chemistry first,
then compute AmsPAF on the aggregated chemistry:

``` r

chr_chem <- time_weighted_aggregate(imp, focal_dates, summary = "geom_mean")
prep_ref <- prepare_reference(chr_ref)
chr_amspaf <- add_amspaf(chr_chem, reference = prep_ref)
```

This gives “AmsPAF at the geometric-mean chronic chemistry.” For sites
with stable chemistry it is acceptable, but for sites with pulsed
exposures (storm events, mobilisation pulses) it can substantially
under-state the toxic burden — see *Pulsed exposures* below.

**Path B — response averaging (recommended).** Compute per-sample AmsPAF
first, then time-aggregate the AmsPAF values:

``` r

prep_ref       <- prepare_reference(ref_chem)
ps_amspaf      <- add_amspaf(imp, reference = prep_ref)
chr_amspaf     <- time_weighted_aggregate(
  dplyr::filter(ps_amspaf, analyte == "AmsPAF"),
  focal_dates    = focal_dates,
  summary        = "arith_mean"
)
```

This time-averages the toxic-response signal directly. It is more
biologically defensible (the community integrates the effect over time,
not the concentration) and is the default in the bundled
`run_amspaf_pipeline.R` script.

### Pulsed exposures and Jensen’s inequality

For a log-normal SSD with HC50 = 100 µg/L and σ_log10 = 0.5, consider
100 days of exposure with 99 days at 1 µg/L (baseline) and 1 day at 1000
µg/L (storm pulse):

| Aggregation                   | Concentration | PAF    |
|-------------------------------|---------------|--------|
| Geometric mean concentration  | 1.99 µg/L     | 0.03 % |
| Arithmetic mean concentration | 11.0 µg/L     | 2.7 %  |
| Mean of per-sample PAFs       | —             | 0.91 % |
| Max concentration             | 1000 µg/L     | 97.7 % |

PAF at the geometric mean concentration under-states the time-averaged
PAF by an order of magnitude. This is Jensen’s inequality acting on the
concave upper tail of the SSD: a brief excursion to high concentration
contributes much more to the integrated response than its concentration-
space weight suggests.

For chronic assessment of pulsed exposures, Path B (mean of per-sample
AmsPAF values) is the correct aggregation.

## Choosing `tau_days`

[`time_weighted_aggregate()`](https://www.kedumba.com.au/leachatetools/reference/time_weighted_aggregate.md)
uses exponential-decay temporal weighting with parameter `tau_days`
(default 90). The effective half-life is `tau_days * log(2)` ≈ 62 days
at the default.

Choose `tau_days` to match the response time of the target biology:

| Target | Typical response timescale | Suggested `tau_days` |
|----|----|----|
| Periphyton / algae | days to weeks | 14–30 |
| Macroinvertebrates (general) | weeks to months | 60–120 |
| Fish (sublethal effects) | months to years | 180–365 |
| Long-lived benthic communities | seasonal to multi-year | 365+ |

These are starting points, not prescriptions. For a calibration exercise
the right move is to test a range of `tau_days` values and pick the one
that maximises explanatory power against your biological response data —
within a range that is mechanistically defensible.

## Reference summary statistic

The Added Risk Approach (ARA) subtracts a background “reference”
concentration before evaluating the SSD, so the AmsPAF reflects only the
*increment above* what the local community is already adapted to. The
summary statistic used to define that background is chosen by the
`summary` argument to
[`prepare_reference()`](https://www.kedumba.com.au/leachatetools/reference/prepare_reference.md).

The package default is `"geom_mean"`. Rationale:

- For log-normal concentration distributions (the typical case for
  aquatic chemistry), the geometric mean is the maximum-likelihood
  central tendency.
- It uses every reference observation, making it more robust to small
  reference datasets than a fixed quantile.
- It is consistent with the PICT (pollution-induced community tolerance)
  interpretation: the community has been integrated against the typical
  exposure, not the upper tail.

Alternatives:

- `"median"` — most robust to outliers; good when reference data may be
  contaminated.
- `"arith_mean"` — uses every observation but is sensitive to outliers.
- `"p80"`, `"p90"`, `"p95"` — quantiles. Use when the assessment
  framework specifies an upper-percentile reference (e.g., some ANZG
  compliance protocols).

Set `bootstrap_ci = TRUE` to compute a 95% bootstrap CI on the chosen
summary statistic per analyte. Useful when the reference dataset is
small (n \< 20 per analyte) and you want to propagate that uncertainty
visually into downstream interpretation.

## Why no tier breaks?

Previous versions of this package shipped four tier breaks (1%, 5%, 10%,
20% AmsPAF) mapped to ANZG species protection levels. These have been
removed.

Reasons:

1.  As discussed above, AmsPAF and single-substance protection levels
    are not interchangeable.
2.  The “right” threshold for community impact is empirical — it depends
    on biological calibration (which the package does not provide) and
    on the management context (which the user knows).
3.  Pre-baked tier breaks invite false confidence. Better to surface the
    continuous AmsPAF value and let consumers attach interpretive
    thresholds appropriate to their assessment.

If you need tier breaks for reporting, consider:

- Calibrating against macroinvertebrate community data and choosing
  breakpoints empirically (e.g., the AmsPAF value above which Baselga
  nestedness consistently exceeds a reference threshold).
- Adopting external standards (e.g., NEMP guidance documents) and
  documenting them explicitly as your interpretive choice.
- Using a discrete classification (background / elevated / impacted)
  with bands chosen from reference-site AmsPAF distribution percentiles.

## Known caveats

**Mercury and methylation.** The Hg SSD bundled with this package is for
inorganic Hg. Methylmercury (MeHg), which forms in low-redox aquatic
sediments and is the form most toxic to vertebrates, is not separately
handled — ANZG does not publish a methylmercury SSD. At sites where
reductive conditions are expected (leachate plumes, anaerobic sediment),
the inorganic Hg PAF will under-state the actual biological risk by
several orders of magnitude. This is a limitation in upstream guideline
data, not in the package.

**Mode-of-action grouping.** The current metadata places all metals into
a single `ionoregulatory` group (Concentration Addition within the
group). This is the conservative choice — it produces a higher AmsPAF
than independent action across the metals — but it is mechanistically
imprecise. Different metals act on different molecular targets (Cu/Zn on
gill Na⁺/K⁺-ATPase, Pb on calcium signalling, Cd on metallothionein, Hg
on thiol groups, As on phosphate metabolism). A more accurate treatment
would use multiple MOA groups with CA within and IA across. This
requires evidence-based assignment from the ANZG/ANZECC technical
briefs, which is tracked as a future task.

**Biological validation.** The SSDs are validated against published ANZG
DGVs (see HANDOVER.md §4) but the AmsPAF as a community-state predictor
is not yet calibrated. Until macroinvertebrate calibration is complete,
AmsPAF outputs are best-effort predictions of chronic toxicity; they
have not been empirically tied to observed community impact.

**Below-detection handling in imputation.** The
[`impute_chemistry()`](https://www.kedumba.com.au/leachatetools/reference/impute_chemistry.md)
function uses `mi()` (missing imputation) rather than `cens("left")`
(left-censoring) because brms does not currently support censored
likelihoods with multivariate residual correlation (`rescor = TRUE`).
The post-hoc cap (`bdl_cap = TRUE`) clips imputed BDL values to the
original detection limit when the model predicts above DL. For sites
where the chemistry context legitimately suggests high concentrations
the cap may fire frequently; results in that regime should be inspected.
See the roxygen documentation of
[`fit_imputation_model()`](https://www.kedumba.com.au/leachatetools/reference/fit_imputation_model.md)
for details on the trade-off.
