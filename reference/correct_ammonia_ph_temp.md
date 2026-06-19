# Correct total ammonia-N to the ANZG reference pH and temperature

The freshwater ammonia default guideline values (DGVs) and the bundled
ammonia SSD are expressed as **total ammonia-N at the ANZG index
condition of pH 7.0 and 20 °C**. Ammonia toxicity is driven by the
*un-ionised* fraction (NH3), which rises with both pH and temperature. A
measured total ammonia-N must therefore be converted to the equivalent
reference-condition concentration — the one that holds the same
un-ionised NH3 — before it is compared against a DGV or passed to
[`ssd_paf()`](https://vorpalvorpal.github.io/hydroSense/reference/ssd_paf.md)
as `"NH3-N"`.

## Usage

``` r
correct_ammonia_ph_temp(
  conc,
  conc_units = NULL,
  pH,
  temperature_C,
  ref_pH = 7,
  ref_temperature_C = 20
)
```

## Arguments

- conc:

  Measured total ammonia-N. Numeric vector or `units` object. Bare
  numeric requires `conc_units`.

- conc_units:

  Character. Unit of `conc` when it is bare numeric, e.g. `"ug/L"` or
  `"mg/L"`. Ignored when `conc` is a `units` object.

- pH:

  Sample pH. Recycled against `conc`.

- temperature_C:

  Sample temperature (°C). Recycled against `conc`.

- ref_pH, ref_temperature_C:

  Reference index condition. Defaults to the ANZG ammonia DGV basis (pH
  7.0, 20 °C); change only if your SSD/DGV uses a different reference.

## Value

Total ammonia-N normalised to the reference condition (µg/L), the length
of the recycled inputs. Pass this to
`ssd_paf("NH3-N", conc = ., conc_units = "ug/L")`.

## Details

The conversion is `C_ref = C * f_sample / f_ref`, where
`f = 1 / (1 + 10^(pKa - pH))` is the un-ionised fraction and
`pKa = 0.09018 + 2729.92 / T(K)` (Emerson et al. 1975), with
`T(K) = temperature_C + 273.15`. A sample at high pH/temperature is
normalised **upward** (more of its ammonia is in the toxic NH3 form than
at the reference), and a sample at low pH/temperature **downward**.

## Do not double-correct

[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
applies this same correction **automatically** from the metadata
`normalisation_formula` (it reads the per-sample `pH` and `temperature`
columns). Use this helper only for the manual
[`ssd_paf()`](https://vorpalvorpal.github.io/hydroSense/reference/ssd_paf.md)
path. Do **not** pre-correct with this helper and then pass the result
to
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md),
or ammonia will be corrected twice.

## References

Emerson K, Russo RC, Lund RE, Thurston RV (1975). Aqueous ammonia
equilibrium calculations: effect of pH and temperature. Journal of the
Fisheries Research Board of Canada 32(12):2379–2383.

## Examples

``` r
# 900 µg/L total ammonia-N measured at pH 8.5, 20 °C is far more toxic than
# the same number at the pH 7.0 reference, so it normalises sharply upward:
correct_ammonia_ph_temp(900, conc_units = "ug/L", pH = 8.5, temperature_C = 20)
#> [1] 25394.81
```
