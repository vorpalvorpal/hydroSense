# Probabilistic NO3-N hardness-class weights.

The three nitrate SSDs (soft / moderate / hard) have abrupt class
boundaries at 30 and 150 mg/L CaCO3. Rather than snap a sample to one
class, treat the true hardness as log-normal around the measured value
(CV = `hardness_cv`) and weight the classes by the probability mass
falling in each band. This smooths the PAF across boundaries without
biologically interpolating between the distinct SSDs.

## Usage

``` r
.no3_weights(hardness_mg_L, hardness_cv = 0.05)
```

## Arguments

- hardness_mg_L:

  Numeric. Measured hardness in mg/L as CaCO3.

- hardness_cv:

  Numeric. CV of the hardness measurement (default 0.05).

## Value

Named numeric vector `c(soft, mod, hard)` summing to 1, or all-`NA` if
`hardness_mg_L` is missing/non-positive.
