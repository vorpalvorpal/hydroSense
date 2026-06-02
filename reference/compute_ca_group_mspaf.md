# Compute msPAF for a single Concentration Addition group

Implements the multispecies Concentration Addition model of De Zwart &
Posthuma (2005, Environ Toxicol Chem 24(10):2665-2676), eq. 6: hazard
units (here `TU = C_adj / HC50`) are summed across the components
sharing a toxic mode of action, and the combined proportion affected is
the log-normal CDF evaluated at `log10(ΣTU)` with a mixture slope equal
to the **arithmetic mean of the component slopes** (their `β̄_TMoA`; in
the normal-CDF form the slope is the standard deviation `σ̄ = mean(σ)`):

## Usage

``` r
compute_ca_group_mspaf(group_data)
```

## Arguments

- group_data:

  Tibble with columns `C_adj`, `hc50`, `sigma`, `moa_group`; one row per
  analyte in the CA group.

## Value

msPAF as a proportion (not percentage).

## Details

\$\$msPAF\_{CA} = \Phi\\\left( \frac{\log\_{10}(\sum TU)}{\bar{\sigma}}
\right)\$\$

Note this is a plain (unweighted) average of the per-analyte sigmas, per
the primary source — NOT a TU-weighted or variance-style combination. A
single-component group therefore reduces exactly to that component's own
SSD.
