# Record what computed a derived product

LMF and msPAF are derived products: they depend on this package's code
as much as on the data. `hydrosense_provenance()` returns a small named
list identifying the code and the time of a computation, for the caller
to store beside the derived values (hydroSense itself writes nothing
anywhere).

## Usage

``` r
hydrosense_provenance()
```

## Value

A named list:

- package_version:

  `character(1)`, the installed hydroSense version.

- git_sha:

  `character(1)` commit SHA, or `NA_character_`.

- r_version:

  `character(1)`, the R version.

- computed_at:

  `POSIXct(1)` in UTC, the time of the call.

## Details

`git_sha` comes from the `RemoteSha` field recorded when the package was
installed from GitHub. When the package is loaded from a source checkout
(for example with `devtools::load_all()`) it is the checkout's `HEAD`,
read with `git rev-parse`. Otherwise (a CRAN or local-tarball install,
or no usable `git`) it is `NA`.

## Examples

``` r
hydrosense_provenance()
#> $package_version
#> [1] "0.2.0"
#> 
#> $git_sha
#> [1] NA
#> 
#> $r_version
#> [1] "4.6.1"
#> 
#> $computed_at
#> [1] "2026-10-25 12:47:32 UTC"
#> 
```
