# Duration-weighted 90th percentile for a set of (focal date x series) windows

Expands each window to its rows, weights them exactly as the per-window
definition does, and takes the first value (in value order, ties in row
order) whose cumulative weight share reaches 0.9.

## Usage

``` r
.twa_p90(
  qi,
  qg,
  gb,
  r_key,
  k_r,
  v_r,
  pos,
  mid,
  dlt,
  f,
  tau_days,
  a_k,
  hi,
  has_a,
  w_a,
  w_h
)
```
