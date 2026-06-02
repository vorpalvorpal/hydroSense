# PCA variables that must NOT be log-transformed before the chemistry PCA

Every other PCA variable is concentration-like — strictly positive and
strongly right-skewed, spanning orders of magnitude — so it is
`log10`-transformed before centring/scaling. Without that, the PCA is
dominated by a handful of high-magnitude major ions (e.g. Cl, SO4, TDS)
and the leading axes mostly track absolute ionic strength rather than
the multiplicative covariance structure that drives metal/organic
behaviour.

## Usage

``` r
.PCA_NO_LOG_VARS
```

## Format

An object of class `character` of length 4.

## Details

The exclusions are the variables for which a log is meaningless or
undefined:

- `pH` — already a logarithmic scale (−log10 of H+ activity).

- `temperature` — interval scale (°C); zero/negative values are valid.

- `ORP` — redox potential (mV); routinely negative.

- `DO` — dissolved oxygen (mg/L); legitimately ~0 in anoxic leachate
  plumes and only spans a narrow, near-linear range.
