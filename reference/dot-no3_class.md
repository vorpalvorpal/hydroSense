# Resolve the appropriate NO3-N hardness-class analyte name.

Resolve the appropriate NO3-N hardness-class analyte name.

## Usage

``` r
.no3_class(hardness_mg_L, hardness_cv = 0.05)
```

## Arguments

- hardness_mg_L:

  Numeric. Measured hardness in mg/L as CaCO3.

- hardness_cv:

  Numeric. Accepted for signature parity with
  [`ssd_paf()`](https://vorpalvorpal.github.io/hydroSense/reference/ssd_paf.md);
  not used here (this helper applies hard cutoffs). The probabilistic
  blend lives in
  [`.no3_weights()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-no3_weights.md).

## Value

Character: "NO3-N_soft", "NO3-N_mod", or "NO3-N_hard".
