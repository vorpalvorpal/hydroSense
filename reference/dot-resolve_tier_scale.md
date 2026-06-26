# Resolve a per-analyte OU variance scale from `ou_scale` and a tier

`ou_scale` may be a single number (applied to every analyte) or a named
numeric keyed by impact tier, e.g. `c(model = 1, bridge = 2.5)`. A tier
absent from a named vector falls back to `1`.

## Usage

``` r
.resolve_tier_scale(ou_scale, tier)
```
