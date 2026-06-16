# Fit a season-blind predictive model of the site impact

Models the anthropogenic increment `I = C_norm - ref_norm` (the ARA
"added risk", i.e.
[`ara_summary()`](https://vorpalvorpal.github.io/leachatetools/reference/ara_summary.md)'s
`C_excess`) at a target site as a function of hydrology and a persistent
latent state – **never** of day-of-year. Used by
[`amspaf_daily()`](https://vorpalvorpal.github.io/leachatetools/reference/amspaf_daily.md)
(`interpolation = "model"`) to fill the gaps between grab samples with a
chemistry-grounded impact estimate instead of a forward-filled
concentration.

## Usage

``` r
fit_target_model(
  target,
  reference_model,
  hydro = NULL,
  hydro_type = "rainfall",
  imputation_model = NULL,
  conc_units = NULL,
  analyte_metadata = NULL,
  method = c("multi", "anzecc"),
  guideline_dir = getOption("leachatetools.guideline_dir"),
  transform = c("pseudo_log", "additive"),
  analyte_c = NULL,
  api_windows_short = c(3L, 7L, 14L),
  api_windows_long = c(30L, 60L, 90L, 180L),
  auto_select = TRUE,
  min_obs_model = 12L,
  pool = TRUE,
  eps = 1e-09
)
```

## Arguments

- target:

  Long-format target chemistry. Required columns: `sample_id`,
  `datetime`, `analyte`, `value`, `detected`. Toxicants must be in ug/L;
  supply via a `units.analyte` column or `conc_units`. Co-analyte rows
  (pH, DOC, hardness, temperature) should be present for normalisation.

- reference_model:

  A `reference_model` from
  [`fit_reference_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_reference_model.md).
  Supplies `ref_norm(t)` and (by default) the catchment hydrology
  series.

- hydro:

  Optional target-specific daily hydrology data frame (`date`, `value`);
  when `NULL` (default) the reference model's hydrology is reused
  (shared-catchment assumption). A target-local stage/discharge gauge
  can capture breach-mobilising flow that catchment rainfall misses.

- hydro_type:

  `"rainfall"`, `"stage"`, or `"discharge"`; used only when `hydro` is
  supplied. Default `"rainfall"`.

- imputation_model:

  Optional `imputation_model` from
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)
  (fit on the target's chemistry). When supplied, missing analytes are
  imputed in raw concentration space before the impact is computed,
  adding more anchor days (tier 2). Requires **brms**.

- conc_units:

  Unit string for target chemistry when no `units.analyte` column is
  present.

- analyte_metadata:

  Analyte metadata, or `NULL` for the bundled CSV.

- method:

  SSD method (`"multi"` or `"anzecc"`) used to derive each analyte's HC5
  transform scale when `analyte_c` is not supplied.

- guideline_dir:

  Path to the ANZG guideline data folder (for the SSD fits); falls back
  to `getOption("leachatetools.guideline_dir")`.

- transform:

  `"pseudo_log"` (default) or `"additive"`. Controls the
  variance-stabilising transform applied to the impact residual before
  smoothing. `"pseudo_log"` uses `g = asinh(I / c)` with per-analyte
  scale `c = HC5` (issue \#15), which compresses the dynamic range and
  prevents event spikes from inflating the baseline draw spread.
  `"additive"` keeps `g = I` (pre-#15 behaviour): the smoother operates
  in the original additive impact space. Ignored when `analyte_c` is
  supplied directly.

- analyte_c:

  Optional named numeric vector of per-analyte transform scales `c` (SSD
  HC5; issue \#15). When `NULL` (default) it is computed from the fitted
  SSDs. The impact residual is smoothed on the variance-stabilising
  scale `g = asinh(I / c)`; an analyte with `NA`/absent `c` keeps the
  additive model.

- api_windows_short, api_windows_long:

  Candidate short/long antecedent memory windows (days) for `f(hydro)`,
  selected by AIC.

- auto_select:

  Logical; AIC window selection per analyte (default `TRUE`).

- min_obs_model:

  Integer; minimum impact anchors required to attempt the `f(hydro)`
  GAM. Below this, the analyte uses the bridge tier. Default `12L`.

- pool:

  Logical (default `TRUE`). When `TRUE`, the per-analyte hydro responses
  are **partially pooled**: a single factor-smooth GAM (`bs = "fs"`) is
  fitted across all sufficiently-sampled analytes at one common
  AIC-selected window, shrinking each analyte's response toward a shared
  shape. This *regularises* noisy, low-signal analytes (it does not add
  hydrological coverage – co-sampled analytes already share the same
  regimes), and falls back to independent fits if it fails or doesn't
  beat an analyte-intercept null. Set `pool = FALSE` to force
  independent per-analyte fits (appropriate only when all analytes are
  densely sampled).

- eps:

  Small positive guard. Default `1e-9`.

## Value

An object of class `target_model`:

- `$models`:

  Named per-analyte list: `impact_fit` (gam or `NULL`), `window_short`,
  `window_long`, `tier` (`"model"` or `"bridge"`), `n_obs`, and
  `anchors` (tibble `date`, `I`, `S`, `hydro_short`, `hydro_long`).

- `$reference_model`:

  The supplied reference model.

- `$hydro`,`$hydro_type`:

  Hydrology series used for `f(hydro)`.

- `$fit_date`:

  Date fitted.

## Details

**Why season-blind.** Site impacts are driven by management failure (a
leachate breach), not the calendar. Heavy summer rain *enables* an
existing breach to escape; it does not *cause* one. A perfect-management
year has zero impact under identical rainfall. Day-of-year is therefore
a confounded covariate for impact and is excluded; hydrology enters only
as a modulator of how an already-present impact expresses itself
(first-flush mobilisation, sorption, antecedent memory).

**Model (per analyte).** \$\$I(t) = \beta \cdot f(\mathrm{hydro}\_t) +
S(t)\$\$ where `f(hydro)` is a thin-plate GAM on the short- and
long-window antecedent indices (no cyclic seasonal term), fitted only
when it beats an intercept-only null by AIC and at least `min_obs_model`
anchors are available; otherwise `beta = 0` and the analyte falls to a
pure state-interpolation ("bridge") tier. The residual state `S(t)` is
interpolated between observation anchors by a hydrology-weighted bridge
that pinches to the observed residual at each anchor and leans toward
the hydrologically more similar bracketing anchor in between.

This is a season-blind, hydrology-modulated variant of WRTDS-Kalman
(Zhang & Hirsch 2019); `f(hydro)` follows concentration–discharge theory
(Godsey, Kirchner & Clow 2009).

## References

Zhang Q, Hirsch RM (2019) Water Resources Research 55(11):9705–9723.
Godsey SE, Kirchner JW, Clow DW (2009) Hydrological Processes
23:1844–1864.

## See also

[`amspaf_daily()`](https://vorpalvorpal.github.io/leachatetools/reference/amspaf_daily.md),
[`fit_reference_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_reference_model.md),
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md),
[`ara_summary()`](https://vorpalvorpal.github.io/leachatetools/reference/ara_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ref_model <- fit_reference_model(reference_chem, latitude = -33.8,
                                 longitude = 151.2, conc_units = "ug/L")
tgt_model <- fit_target_model(target_chem, ref_model, conc_units = "ug/L")
tgt_model
} # }
```
