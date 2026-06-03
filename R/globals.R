# Quiet R CMD check's "no visible binding for global variable" notes that come
# from non-standard evaluation in dplyr/tidyr verbs (bare column names) and the
# rlang `.data` / `.env` pronouns. These are not real bindings -- they are
# resolved inside the data mask at run time -- so we declare them here.
#
# Also declared: the dashboard-provided functions (data_df, feature_df,
# feature_sfc, get_reference_site) that add_lmf() calls only on its
# non-override path. They are supplied by the host environment at run time and
# are intentionally not defined in this package.
utils::globalVariables(c(
  # column names used in NSE (dplyr/tidyr); the .data / .env pronouns are
  # imported from rlang in leachatetools-package.R instead.
  "analyte", "atomic_mass.analyte", "Cl_", "Conc", "datetime", "detected",
  "f", "f_pct", "gradient", "high_info", "informativeness", "ion", "L",
  "label", "n_values", "name", "R", "sample_id", "sigma_R", "site_id",
  "species_id", "Species", "ssd_available", "total_alk_", "total_N_",
  "valence.analyte", "value", "value_ug_L", "var_f", "weight", "x",
  # host-environment (dashboard) functions used only on add_lmf()'s non-override
  # path; not part of this package
  "data_df", "feature_df", "feature_sfc", "get_reference_site"
))
