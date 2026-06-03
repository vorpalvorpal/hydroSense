# Quiet R CMD check's "no visible binding for global variable" / "no visible
# global function definition" notes that come from non-standard evaluation in
# dplyr/tidyr verbs (bare column names) and from host-environment functions.
# These are not real bindings -- column names are resolved inside the data mask
# at run time -- so we declare them here. The rlang `.data` / `.env` pronouns
# are imported from rlang in leachatetools-package.R instead.
#
# The list is the complete set reported by
#   codetools::checkUsagePackage("leachatetools", all = FALSE)
# (the same analysis R CMD check runs). Regenerate it if that check ever flags
# a new name.
utils::globalVariables(c(
  # NSE column names (dplyr/tidyr)
  "analyte", "atomic_mass.analyte", "Cl_", "Conc", "datetime", "detected",
  "dw_flag", "f", "f_pct", "gradient", "hi_flag", "high_info",
  "informativeness", "ion", "L", "label", "mean_ratio", "n_ref", "n_values",
  "name", "R", "reference", "row_str", "sample_id", "sigma_meas", "sigma_R",
  "site_id", "species_id", "Species", "ssd_available", "total_alk_",
  "total_N_", "total_N_mgl", "uuid", "valence.analyte", "value", "value_ug_L",
  "var_f", "weight", "wt_orig_pct", "wt_rob_pct", "x",
  # host-environment (dashboard) functions used only on add_lmf()'s non-override
  # path; supplied at run time, not defined in this package
  "data_df", "feature_df", "feature_sfc", "get_reference_site"
))
