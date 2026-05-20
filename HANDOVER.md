# leachatetools — Handover Notes

_Last updated: 2026-05-20_

---

## 1. What this package does

`leachatetools` is an R package providing three core capabilities for
assessing water quality in leachate-impacted freshwater systems:

| Function family | File | Exports |
|---|---|---|
| SSD-based PAF lookup | `R/paf.R` + `R/ssd_fit.R` | `ssd_paf()`, `ssd_hc50()` |
| Multi-substance PAF (msPAF / AmsPAF) | `R/mspaf.R` | `add_amspaf()` |
| Leachate mixing fraction | `R/lmf.R` | `add_lmf()` |

---

## 2. SSD fitting — what was done and why

### 2.1 Data sources

All ANZECC 2000 / ANZG analytes come from one of two sources:

**Warne (2000)** — the original derivation data.  Raw species concentrations
are in `inst/extdata/anzecc_warne2000_observations.csv`.  The vast majority of
analytes (organics, most metals) use this source.

**ANZG XLSX technical briefs** — modern updated datasets from ANZG 2021–2026.
These replace the ANZECC 2000 values where available.  The XLSX files are NOT
bundled in the package (they are large and proprietary); they must be present
in a local "guideline data/" folder at fit time.  Set the path via:
```r
options(leachatetools.guideline_dir = "/path/to/guideline data")
```
XLSX analytes: NH3-N, NO3-N (3 classes), B, Cr, Cu, Ni, Zn.

**Ni special case**: The nickel XLSX (col 31) holds raw test concentrations,
_not_ the MLR-normalised values that ANZG used to derive the DGVs.  The 26
normalised negligible-effect values from Stauber et al. (2021) Table 3
(index conditions: pH 7.5, Ca 6 mg/L, Mg 4 mg/L, DOC 0.5 mg/L) are
pre-tabulated in `inst/extdata/ni_mlr_normalised_table3.csv`.

### 2.2 Distributions

**method = "multi"** (default for `ssd_paf()`):
Fits all 6 BCANZ-recommended distributions and model-averages the result:
`lnorm`, `llogis`, `burrIII3`, `lgumbel`, `gamma`, `weibull`.
This is current best practice and produces good results for all validated
analytes regardless of the original derivation method.

**method = "anzecc"**:
Uses the per-analyte distribution recorded in the `fit_dist` column of
`inst/extdata/anzecc_analyte_metadata.csv`.  This attempts to replicate the
original ANZG derivation as closely as possible.  Useful for back-checking
against published DGVs.

Key distribution aliases in the metadata:
- `burrIII3` — Burr Type III; Burrlioz 2.0 default; used for most ANZECC 2000
  analytes and most modern ANZG briefs.
- `lgumbel` — log-Gumbel = Fréchet = **inverse Weibull** (same distributional
  family, different names in different software). Burrlioz calls it
  "inverse Weibull"; ssdtools calls it "lgumbel". Used for NO3-N soft and
  hard water.
- `multi` — expand to the full 6-distribution set (same as method="multi").
  Used for NH3-N and Cr where ANZG used shinyssdtools with all BCANZ dists.

### 2.3 NO3-N — three separate analytes

Nitrate toxicity is physically modulated by water hardness (calcium/magnesium
ions compete with nitrate at ion-exchange sites).  The same species can be
150× less sensitive to NO3-N in hard water vs soft water — this is a
mechanistic effect, not a species-assemblage difference.  Because hardness
changes the actual effect on any individual organism, there is no valid
biological interpolation between hardness classes; ANZG derived three
completely separate SSDs.

The three analyte names used throughout the package:

| Name | Hardness range | Distribution | HC5 (DGV 95%) |
|---|---|---|---|
| `NO3-N_soft` | < 30 mg/L CaCO₃ | lgumbel | 1,100 µg/L |
| `NO3-N_mod` | 30–150 mg/L CaCO₃ | burrIII3 | 2,600 µg/L |
| `NO3-N_hard` | > 150 mg/L CaCO₃ | lgumbel | 29,000 µg/L |

Callers can either supply the explicit class name or pass `analyte = "NO3-N"`
with `hardness_mg_L` to let `ssd_paf()` / `ssd_lookup()` resolve the class.

**Known problem — class-boundary discontinuity**: When measured hardness is
near a boundary (e.g. 28 vs 31 mg/L), measurement noise causes abrupt jumps
in modelled PAF with no change in actual conditions.  A probabilistic
smoothing approach has been designed but not yet implemented (see §4.1).

### 2.4 ACR adjustment (MR analytes)

Analytes where ANZECC derived triggers from acute LC50/EC50 data using the
Maximum Residue (MR) method have a trigger_divisor (= ACR) > 1.  The fitting
code divides all concentrations by ACR before fitting, so the resulting SSD
is on a chronic-equivalent scale.  `ssd_paf(conc_ug_L = ...)` therefore
accepts the measured chronic-equivalent concentration directly — no further
ACR adjustment is needed by the caller.

---

## 3. Validation state

Validated by `dashboard/R/make_ssd_validation.R` against published ANZG/ANZECC
DGVs.  The ratio `calc/ref` should be close to 1.0; deviations reflect
ssdtools vs Burrlioz 2.0 algorithm differences.

**Summary (ssdtools 2.6, method = "anzecc" / per-analyte dist):**

| Analyte | HC5 ratio | Status |
|---|---|---|
| NH3-N | 1.003 | ✓ |
| NO3-N_soft | 0.980 | ✓ |
| NO3-N_mod | 1.004 | ✓ |
| NO3-N_hard | 0.992 | ✓ |
| B | 1.001 | ✓ |
| Cu | 0.998 | ✓ |
| Ni | 1.016 | ✓ |
| Zn | 0.998 | ✓ |
| Cr | 0.856 | ⚠ accepted — see below |
| All other Warne 2000 analytes | mostly 0.8–1.3 | ✓ |
| Chlorpyrifos | 3.133 | ✗ parked |
| Diazinon | 0.466 | ✗ parked |

**Cr 15% gap (HC5 ratio = 0.856)**: The data extraction is confirmed correct
(13 species, all matching ANZG Table A.1 exactly).  The gap arises because
ANZG used shinyssdtools v0.4.0 while we use ssdtools 2.6; the model-averaging
weights differ between versions.  Accepted as a version artefact; the SSD
is otherwise sound.

**Parked analytes**: Chlorpyrifos and Diazinon have large HC5 discrepancies
that cannot be explained by simple version differences.  They are excluded
from the SSD pipeline and treated as "no model available".  See
`dashboard/data-raw/analyte_metadata_parked.csv` for data and TODO notes.

---

## 4. Pending tasks

### 4.1 NO3-N probabilistic hardness weighting (HIGH PRIORITY)

Design is complete; implementation is a TODO stub in `paf.R .no3_class()`.

The idea: treat measured hardness as log-normally distributed with the
measurement CV (default 5%) and compute the probability of belonging to
each hardness class.  Weight the three PAF values by class probability.

```r
sdlog  <- sqrt(log(1 + cv^2))
p_soft <- plnorm(30,  log(h), sdlog)
p_hard <- 1 - plnorm(150, log(h), sdlog)
p_mod  <- 1 - p_soft - p_hard
# For ssd_paf():
expected_paf <- p_soft * paf_soft + p_mod * paf_mod + p_hard * paf_hard
```

This requires fitting / loading all three NO3-N models on every call when
near a boundary, so the cache should be warm-started for all three at session
init if NO3-N is in the analyte list.

### 4.2 Replace DGV back-calculation in msPAF (HIGH PRIORITY)

`R/mspaf.R` currently back-calculates a log-normal SSD from four published
ANZG DGV values (HC1, HC5, HC10, HC20) via OLS regression on normal
z-scores.  This was a necessary hack before proper SSD fitting existed.

It should be replaced with direct calls to:
- `ssd_hc50(analyte, ...)` for the Concentration Addition denominator
- `ssd_paf(analyte, conc_ug_L, ...)` for the direct PAF calculation

This is a significant refactor.  The DGV-back-calculation approach will
silently give wrong answers for analytes where the true SSD is not
log-normal (burrIII3 analytes diverge from log-normal in the tails).

### 4.3 Generalise index function inputs

`add_amspaf()` and `add_lmf()` currently expect specific data-frame column
names that match the dashboard's internal data structures (`analyteDF`,
`guidelineDF`, etc.).  These should be refactored to accept generic,
documented inputs so they can be called outside the dashboard context.
Document the expected schema clearly in the roxygen docs.

### 4.4 Investigate parked analytes

See `data-raw/analyte_metadata_parked.csv`.  Specifically:
- **Diazinon**: ssdtools HC5 is ~2.1× the published value.  Try fitting with
  burrIII3 only (no model averaging) — Burrlioz does not model-average.
  Check whether ACR adjustment is applied consistently.
- **Chlorpyrifos**: HC5 ratio ~3.1×.  Check whether the published 0.01 µg/L
  is derived from the same 12 chronic values or a different subset.
  Consider whether the ANZG 2019 technical brief has additional data.

---

## 5. Model caching

Models are fitted lazily on first use and cached in two layers:

1. **In-memory** (per R session): `leachatetools:::.paf_mem_cache`
2. **On-disk** (`tools::R_user_dir("leachatetools", "cache")`):
   files named `<stem>__<method>.qs`, e.g. `NH3_N__multi.qs`.

The `method` parameter ("multi" or "anzecc") is part of the cache key, so
both variants can coexist.  Delete the cache directory to force a full
refit after updating source data.

For ANZG_XLSX analytes, `guideline_dir` must be set at fit time.  Once
cached, the model can be loaded without access to the XLSX files.

---

## 6. Key files

```
dashboard/
├── leachatetools/                   ← R package (this package)
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   ├── R/
│   │   ├── paf.R                   ← ssd_paf(), ssd_hc50(), public API
│   │   ├── ssd_fit.R               ← internal fitting infrastructure
│   │   ├── mspaf.R                 ← add_amspaf() (needs §4.2 refactor)
│   │   └── lmf.R                   ← add_lmf()
│   └── inst/extdata/
│       ├── anzecc_analyte_metadata.csv
│       ├── anzecc_warne2000_observations.csv
│       └── ni_mlr_normalised_table3.csv
│
├── data-raw/
│   ├── anzecc_analyte_metadata.csv  ← master metadata (sync to inst/extdata after edits)
│   ├── analyte_metadata_parked.csv  ← Chlorpyrifos, Diazinon — see §4.4
│   └── ni_mlr_normalised_table3.csv ← 26 MLR-normalised Ni values
│
├── R/
│   ├── fit_ssd_models.R            ← batch fitter (dev tool; writes data/ssd_models/)
│   ├── make_ssd_validation.R       ← validation against published DGVs
│   └── ssd_lookup.R                ← dashboard SSD wrapper (delegates to pre-fitted .qs)
│
└── data/
    ├── ssd_models/                 ← pre-fitted .qs files (built by fit_ssd_models.R)
    └── ssd_validation.csv          ← current validation table
```

**Important**: `inst/extdata/anzecc_analyte_metadata.csv` is a copy of
`data-raw/anzecc_analyte_metadata.csv`.  If you edit the master in `data-raw/`,
run:
```r
file.copy("data-raw/anzecc_analyte_metadata.csv",
          "leachatetools/inst/extdata/anzecc_analyte_metadata.csv",
          overwrite = TRUE)
```

---

