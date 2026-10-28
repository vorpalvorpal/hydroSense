# Apportion msPAF to mode-of-action groups and analytes

Rebuilds the concentration-addition / independent-action arithmetic of
[`compute_mspaf_per_sample()`](https://vorpalvorpal.github.io/hydroSense/reference/compute_mspaf_per_sample.md)
from the per-analyte breakdown and splits each sample's msPAF into
additive contributions. Per sample, group PAF `p_g` is clamped to
`1 - 1e-12`, `L_g = -log(1 - p_g)`, the group takes
`msPAF * L_g / sum(L)` and each analyte its `TU / sum(TU)` share of the
group (an equal split if the group's TU sum to zero), so the
contributions sum to the msPAF.

## Usage

``` r
.mspaf_contributions(pafs, ssd_params)
```

## Arguments

- pafs:

  Flat `analyte_pafs` breakdown (point mode: `sample_id`, `site_id`,
  `analyte`, `C_adj`, `PAF`, `moa_group`).

- ssd_params:

  Tibble from
  [`derive_ssd_params()`](https://vorpalvorpal.github.io/hydroSense/reference/derive_ssd_params.md)
  (`hc50`, `sigma`).

## Value

Tibble
`(sample_id, site_id, analyte, group, c_adj, tu, paf, group_mspaf, contribution)`,
percentages in `paf`, `group_mspaf` and `contribution`.
