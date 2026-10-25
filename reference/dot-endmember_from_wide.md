# Build an end-member from a wide, collapsed LMF panel

The statistics half of
[`build_endmember_from_override`](https://vorpalvorpal.github.io/hydroSense/reference/build_endmember_from_override.md),
split out so that `add_lmf(calibration = "trailing")` can slice one
prepared panel per sample instead of re-converting the whole table each
time.

## Usage

``` r
.endmember_from_wide(
  override_wide,
  type,
  min_leachate_samples,
  window_start,
  window_end
)
```

## Arguments

- override_wide:

  Wide per-sample tibble from
  [`collapse_species`](https://vorpalvorpal.github.io/hydroSense/reference/collapse_species.md).

- type, min_leachate_samples:

  As in
  [`build_endmember_from_override`](https://vorpalvorpal.github.io/hydroSense/reference/build_endmember_from_override.md).

- window_start, window_end:

  Informational span of the data.

## Value

As
[`build_endmember_from_override`](https://vorpalvorpal.github.io/hydroSense/reference/build_endmember_from_override.md).
