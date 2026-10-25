# Build an end-member from a caller-supplied override dataframe

Internal helper called by
[`add_lmf`](https://vorpalvorpal.github.io/hydroSense/reference/add_lmf.md)
when `reference_data` or `leachate_data` is supplied. Processes the
override data through the same meq conversion, BDL substitution, and
species-collapsing machinery as the standard builders, but uses all
available samples without date filtering.

## Usage

``` r
build_endmember_from_override(override_df, type, min_leachate_samples = 10L)
```

## Arguments

- override_df:

  Long-format dataframe. All samples are used; samples missing the LMF
  panel are dropped.

- type:

  Character scalar, either `"reference"` or `"leachate"`, controlling
  which output list format is produced.

- min_leachate_samples:

  Minimum number of leachate samples with F measured required to include
  F in the panel. Ignored when `type = "reference"`. Default `10L`.

## Value

For `type = "reference"`: a list with `stats` (tibble of ion, n_ref, R,
sigma_R), `window_start`, `window_end`, `n_samples`. For
`type = "leachate"`: a list with `L_values` (tibble of ion, mean_ratio,
L), `cl_anchor`, `n_samples`, `f_included`, `window_start`,
`window_end`.
