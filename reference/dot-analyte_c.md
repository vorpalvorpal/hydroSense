# Per-analyte transform scale `c` = SSD HC5

Returns the analyte's 5% hazard concentration (HC5) from the fitted SSD,
used as the additive-\>proportional crossover `c` in
[`.g_transform()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-g_transform.md).
HC5 is a self-scaling toxicological anchor in the same normalised
concentration space as the impact `I` (normalisation maps measured
concentration onto the SSD scale), and it lies within the species data
range so it is well-determined (unlike HC1, a lower-tail extrapolation).
Below HC5 the transform is additive, but that region carries PAF \< 5%
and so cannot inflate reported risk.

## Usage

``` r
.analyte_c(fit)
```

## Arguments

- fit:

  A fitted `ssdtools` SSD object (e.g. an element of
  `derive_ssd_params()$fit`).

## Value

A single finite positive HC5 concentration.

## Details

If HC5 cannot give a good fit to a particular analyte/site, HC1
(`proportion = 0.01`) is the alternative to try – a tighter crossover,
but a less stable lower-tail extrapolation.
