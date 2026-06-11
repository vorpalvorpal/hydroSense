# Estimate OU parameters (gamma by MoM, theta by MLE) with a fallback ladder

`gamma` (the marginal/stationary variance, the well-identified moment)
is the sample variance of the anchor residuals, times `scale`. `theta`
(the mean-reversion rate / inverse correlation length, weakly identified
on sparse data) is fitted by 1-D MLE of the exact irregular-spacing OU
likelihood, bounded so the correlation length stays in the range 0.5x to
10x the median anchor spacing.

## Usage

``` r
.estimate_ou_kalman_params(anchor_dates, anchor_S, n_fit_min = 8L, scale = 1)
```

## Arguments

- anchor_dates:

  Date (or coercible) vector of anchor dates.

- anchor_S:

  Numeric residual S at each anchor (same length).

- n_fit_min:

  Minimum anchors to attempt the theta MLE (default 8).

- scale:

  Positive multiplier on `gamma` (the envelope-width knob; correlation
  length / theta is unchanged). Default 1.

## Value

list(theta, gamma, tier).

## Details

Tiers: `"degenerate"` (\< 2 finite anchors or ~constant S -\> no
uncertainty), `"fallback"` (\>= 2 but \< `n_fit_min` anchors -\> theta
pinned at the long- correlation bound, approximating a Brownian bridge
without diffuse init), `"mle"` (\>= `n_fit_min` -\> theta by MLE).
