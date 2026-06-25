# Quiet R CMD check's "no visible binding for global variable" / "no visible
# global function definition" notes that come from non-standard evaluation in
# dplyr/tidyr verbs (bare column names) and from host-environment functions.
# These are not real bindings -- column names are resolved inside the data mask
# at run time -- so we declare them here. The rlang `.data` / `.env` pronouns
# are imported from rlang in hydroSense-package.R instead.
#
# The list is the complete set reported by
#   codetools::checkUsagePackage("hydroSense", all = FALSE)
# (the same analysis R CMD check runs). Regenerate it if that check ever flags
# a new name.
## ── Unit-resolution helpers ──────────────────────────────────────────────────
##
## Every unit-bearing parameter in the public API accepts either a `units`
## object (auto-converted to the target unit) or a bare numeric plus a
## companion `*_units` character argument that declares what the numeric is.
## Bare numeric without a companion → hard error (no silent unit assumptions).
##
## .resolve_to()             — scalar or vector, any target unit
## .convert_df_tox_to_ugL() — convert SSD-eligible rows in a long-format
##                            data frame to µg/L before SSD lookup

#' Resolve a bare numeric or units object to a target unit
#' @param x Numeric or units object.
#' @param target_unit Character udunits2 unit string (desired output).
#' @param units_str Companion unit string; required when `x` is bare numeric.
#' @param arg_name Name of the calling argument, used in error messages.
#' @return Bare numeric in `target_unit`.
#' @keywords internal
.resolve_to <- function(x, target_unit, units_str = NULL, arg_name = "x") {
  if (inherits(x, "units")) {
    return(as.numeric(units::set_units(x, target_unit, mode = "standard")))
  }
  if (is.null(units_str)) {
    cli::cli_abort(
      c("{.arg {arg_name}} is bare numeric \u2014 units are unknown.",
        "i" = "Pass a companion {.arg {arg_name}_units} string, e.g. \\
               {.code {arg_name}_units = \"mg/L\"}.",
        "i" = "Or wrap directly: \\
               {.code units::set_units({x[1L]}, \"{target_unit}\")}."
      ),
      call = rlang::caller_env()
    )
  }
  # Fast path: when the source and target unit strings are identical the
  # conversion factor is exactly 1, so the udunits round-trip is a no-op
  # (`x * 1 == x` bit-for-bit). Skip it. `units::set_units()` string parsing is
  # the single dominant cost in the SSD/PAF hot loop, where concentrations are
  # already in the target unit (e.g. `ssd_paf(..., conc_units = "ug/L")` runs on
  # a frame `mspaf_daily()` has pre-converted to ug/L) \u2014 there the round-trip
  # parses "ug/L" -> "ug/L" once per analyte per draw for nothing.
  if (identical(units_str, target_unit)) {
    return(as.numeric(x))
  }
  as.numeric(
    units::set_units(
      units::set_units(x, units_str, mode = "standard"),
      target_unit,
      mode = "standard"
    )
  )
}

#' Convert SSD-eligible analyte rows in a long-format data frame to µg/L
#'
#' Uses `units.analyte` per-row when the column is present; otherwise applies
#' `conc_units` uniformly to all SSD-eligible rows.  Both paths error loudly
#' rather than silently assuming a unit.
#'
#' @param df Long-format data frame with `analyte` and `value` columns.
#' @param ssd_analytes Character vector of SSD-eligible analyte names.
#' @param conc_units Character unit string; required when `units.analyte` is
#'   absent from `df`.
#' @param call_arg Name of the data-frame argument in the caller, for errors.
#' @return `df` with SSD-eligible `value` rows converted to µg/L.
#' @keywords internal
.convert_df_tox_to_ugL <- function(df, ssd_analytes, conc_units = NULL,
                                   call_arg = "df") {
  if (nrow(df) == 0L) {
    return(df)
  }
  tox_mask <- df$analyte %in% ssd_analytes
  if (!any(tox_mask)) {
    return(df)
  }

  if ("units.analyte" %in% names(df)) {
    u <- df$units.analyte[tox_mask]
    if (anyNA(u) || any(!nzchar(trimws(u)))) {
      cli::cli_abort(
        c("Missing {.field units.analyte} for an SSD-eligible analyte row.",
          "i" = "Every toxicant row must carry a non-empty \\
                 {.field units.analyte} value."
        )
      )
    }
    # Conversion to ug/L is linear with zero offset for concentration units, so
    # it depends only on the unit string, not the value. Resolve the factor once
    # per distinct unit (a handful of udunits calls) and apply it as a vectorised
    # multiply — bit-identical to the per-row set_units round-trip
    # (`v * slope == slope * v` in IEEE) but O(distinct units) parses instead of
    # O(rows x draws). This per-row pmap was ~58% of draws-mode runtime.
    uniq <- unique(u)
    facs <- vapply(uniq, function(uu) {
      as.numeric(units::set_units(
        units::set_units(1, uu, mode = "standard"), "ug/L", mode = "standard"
      ))
    }, numeric(1))
    df$value[tox_mask] <- df$value[tox_mask] * facs[match(u, uniq)]
    return(df)
  }

  if (is.null(conc_units)) {
    cli::cli_abort(
      c("{.arg {call_arg}} has no {.field units.analyte} column and \\
         {.arg conc_units} was not supplied.",
        "i" = "Add a {.field units.analyte} column (one unit string per row),",
        "i" = "or supply {.arg conc_units} to apply a uniform unit to all \\
               SSD-eligible rows, e.g. {.code conc_units = \"mg/L\"}."
      ),
      call = rlang::caller_env()
    )
  }

  df$value[tox_mask] <- as.numeric(
    units::set_units(
      units::set_units(df$value[tox_mask], conc_units, mode = "standard"),
      "ug/L",
      mode = "standard"
    )
  )
  df
}

utils::globalVariables(c(
  # NSE column names (dplyr/tidyr)
  "analyte", "atomic_mass.analyte", "Cl_", "Conc", "datetime", "detected",
  "draw_id", "n_draws", "value_lower", "value_upper",
  "dw_flag", "f", "f_pct", "gradient", "hi_flag", "high_info",
  "informativeness", "ion", "L", "label", "mean_ratio", "n_ref", "n_values",
  "name", "R", "reference", "row_str", "sample_id", "sigma_meas", "sigma_R",
  "site_id", "species_id", "Species", "ssd_available", "total_alk_",
  "total_N_", "total_N_mgl", "uuid", "valence.analyte", "value", "value_ug_L",
  "var_f", "weight", "wt_orig_pct", "wt_rob_pct", "x",
  # reference_model / temporal ARA (reference_model.R, mspaf.R)
  ".date", "ara_diag", "best_aic", "C_adj", "C_excess", "C_norm",
  "coanalytes_required", "date", "doy", "fit_date", "floor_fired",
  "focal_date", "hydro_long", "hydro_short", "log_pred", "n_obs",
  "normalisation_formula", "null_aic", "rainfall_mm", "ref_norm",
  "ref_source", "ref_tier", "tau_long", "tau_short", "value_norm",
  # mspaf_daily (mspaf_daily.R)
  ".measured", "mspaf", "days_since_last_sample", "n_measured_analytes",
  # target_model (target_model.R)
  "I", "S", "impact", "impact_tier",
  # host-environment (dashboard) functions used only on add_lmf()'s non-override
  # path; supplied at run time, not defined in this package
  "data_df", "feature_df", "feature_sfc", "get_reference_site"
))
