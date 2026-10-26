## Numerical stability of msPAF under ulp-level input changes
## ============================================================================
##
## CONTRACT (work in progress; tests follow).
##
## S1. "fsc" factor-smooth basis: fs with a canonical null-space basis.
## S2. Pooled tau AIC is continuous in tau.
## S3. Marginal imputation GAMs: REML, PC-count gate m = min(P, floor(n_det/3) - 1).
## S4. End-to-end perturbation stability of mspaf_pipeline().
## S5. Invariants.
##
## All data are synthetic (helper-perf-fixtures.R and frames built here).
