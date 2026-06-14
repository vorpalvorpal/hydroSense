# Draw correlated residual trajectories via the DK simulation smoother

Produces `nsim` draws from the Durbin–Koopman (2002) simulation smoother
using **externally-supplied** standard-normal noise matrices. Because
the noise is passed in by the caller (rather than drawn internally), the
caller can correlate the innovations *across analytes* while keeping
every per-analyte marginal distribution unchanged.

## Usage

``` r
.kalman_draw_coupled(setup, eta_std, a1_z, eps_std)
```

## Arguments

- setup:

  Return value of
  [`.kalman_sim_smoother_setup()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-kalman_sim_smoother_setup.md).

- eta_std:

  `[n_grid × nsim]` matrix of iid N(0,1) process innovations.

- a1_z:

  Length-`nsim` vector of N(0,1) initial-state noise.

- eps_std:

  `[m × nsim]` matrix of iid N(0,1) observation noise (one row per
  anchor position).

## Value

Numeric matrix `[n_grid × nsim]` of residual draws in original
(un-standardised) units.

## Details

The DK identity (Algorithm 1 of Durbin & Koopman 2002):
\$\$\tilde{\alpha} = \hat{\alpha} + (\alpha^+ - L y^+)\$\$ where
\\\alpha^+\\ is a prior simulation driven by the supplied innovations
and \\L y^+ = \hat{\alpha}^+\\ is the KFS smoother applied to the
synthetic observations \\y^+\\.

All arithmetic is performed in standardised units; the result is
un-standardised by `setup$resid_scale` before returning.

## References

Durbin J, Koopman SJ (2002). A simple and efficient simulation smoother
for state space time series analysis. *Biometrika* **89**(3), 603–616.
<https://doi.org/10.1093/biomet/89.3.603>
