# Select API memory windows per analyte using AIC

Evaluates every `(window_short, window_long)` candidate pair and picks
the pair with the lowest AIC. Among pairs within 2 ΔAIC units of the
best, the most parsimonious (smallest `window_short + window_long`) is
chosen.

## Usage

``` r
.select_api_windows(
  obs,
  hydro,
  hydro_type,
  api_windows_short,
  api_windows_long,
  eps = 1e-09
)
```

## Details

Returns a list `(window_short, window_long, best_aic, null_aic)`. If all
candidate fits fail, returns the first candidate pair with
`best_aic = Inf`.
