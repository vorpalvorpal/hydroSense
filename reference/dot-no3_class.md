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

  Numeric. Coefficient of variation of hardness measurement (default
  0.05 = 5%). Reserved for future probabilistic weighting — currently
  unused (hard cutoffs applied).

## Value

Character: "NO3-N_soft", "NO3-N_mod", or "NO3-N_hard".
