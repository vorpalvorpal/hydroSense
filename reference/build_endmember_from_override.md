# Build an end-member from a caller-supplied override dataframe

Internal helper called by
[`add_lmf`](https://vorpalvorpal.github.io/leachatetools/reference/add_lmf.md)
when `reference_data` or `leachate_data` is supplied. Processes the
override data through the same meq conversion, BDL substitution, and
species-collapsing machinery as the standard builders, but uses all
available samples without date filtering.

## Usage

``` r
build_endmember_from_override(override_df, type)
```

## Arguments

- override_df:

  Long-format dataframe in `data_df()` structure. All samples are used;
  samples missing the LMF panel are dropped.

- type:

  Character scalar, either `"reference"` or `"leachate"`, controlling
  which output list format is produced.

## Value

For `type = "reference"`: a list with the same structure as
[`build_reference_endmember`](https://vorpalvorpal.github.io/leachatetools/reference/build_reference_endmember.md):
`stats` (tibble of ion, n_ref, R, sigma_R), `window_start`,
`window_end`, `n_samples`. For `type = "leachate"`: a list with the same
structure as
[`build_leachate_endmember`](https://vorpalvorpal.github.io/leachatetools/reference/build_leachate_endmember.md):
`L_values` (tibble of ion, mean_ratio, L), `cl_anchor`, `n_samples`,
`f_included`, `window_start`, `window_end`.
