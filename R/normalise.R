# Internal helpers for chemistry normalisation
#
# Normalisation converts a measured field concentration to the ANZG index
# condition before SSD lookup, e.g. Cu at DOC=0.5 mg/L or Zn at
# hardness=30 mg/L. The formulas are stored as R-expression strings in the
# `normalisation_formula` column of `inst/extdata/anzecc_analyte_metadata.csv`.
#
# Bioavailability/index-condition formulas are populated for Cu, Ni, Zn, Cd,
# Pb (hardness/DOC/pH bioavailability corrections) and NH3-N (pH/temperature
# un-ionised fraction). Analytes with an empty `normalisation_formula` cell
# fall through `.apply_normalisation()` unchanged (identity), which is correct
# for analytes whose SSD is not chemistry-dependent.
#
# Direction convention: the formula maps a measured field concentration TO the
# reference/index condition. For NH3-N this is `C * f_sample / f_ref` — a
# high-pH sample (more un-ionised NH3, more toxic) normalises UPWARD. The
# exported `correct_ammonia_ph_temp()` applies the identical conversion for the
# manual `ssd_paf()` path; the two must stay in sync.

# Internal cache: parsed expressions are expensive-ish to parse repeatedly and
# are pure functions of the formula string. Cache in a private env keyed by
# the formula string.
.normalise_cache <- new.env(parent = emptyenv())

#' Parse a normalisation formula string into a quoted expression
#'
#' @param formula_str Character string containing an R expression. The
#'   expression may reference `C` (the raw concentration) and any co-analyte
#'   names listed in `coanalytes_required`. Empty string or `NA` → identity
#'   (returns `C` unchanged).
#' @return A language object (quoted expression), or `NULL` for identity.
#' @keywords internal
.parse_normalisation_formula <- function(formula_str) {
  if (is.na(formula_str) || !nzchar(trimws(formula_str))) {
    return(NULL)  # NULL = identity: return C as-is
  }

  key <- formula_str
  if (exists(key, envir = .normalise_cache, inherits = FALSE)) {
    return(get(key, envir = .normalise_cache, inherits = FALSE))
  }

  parsed <- tryCatch(
    parse(text = formula_str, keep.source = FALSE)[[1L]],
    error = function(e) {
      cli::cli_abort(
        c(
          "Failed to parse normalisation formula.",
          "x" = "Formula: {.code {formula_str}}",
          "x" = "Error: {conditionMessage(e)}"
        )
      )
    }
  )

  assign(key, parsed, envir = .normalise_cache)
  parsed
}

#' Apply a parsed normalisation formula to a concentration value
#'
#' @param parsed_expr Parsed expression from `.parse_normalisation_formula()`,
#'   or `NULL` for identity.
#' @param C Numeric concentration (µg/L or relevant unit).
#' @param coanalytes Named numeric vector of co-analyte values available at
#'   this sample (e.g. `c(DOC = 0.5, pH = 7.2)`). May be empty.
#' @return Normalised concentration. Returns `C` when `parsed_expr` is `NULL`.
#'   Returns `NA_real_` if evaluation fails (e.g. required co-analyte absent
#'   from `coanalytes`).
#' @keywords internal
.apply_normalisation <- function(parsed_expr, C, coanalytes = numeric(0)) {
  if (is.null(parsed_expr)) {
    return(C)
  }

  env <- list2env(as.list(coanalytes), parent = baseenv())
  env$C <- C

  tryCatch(
    eval(parsed_expr, envir = env),
    error = function(e) NA_real_
  )
}

# ── User-facing ammonia pH/temperature correction ─────────────────────────────

#' Un-ionised ammonia fraction (NH3) of total ammonia-N
#'
#' @param pH Sample pH.
#' @param temperature_C Sample temperature (°C).
#' @return Fraction (0–1) of total ammonia-N present as un-ionised NH3, using
#'   `f = 1 / (1 + 10^(pKa - pH))` with `pKa = 0.09018 + 2729.92 / T(K)`
#'   (Emerson et al. 1975), `T(K) = temperature_C + 273.15`.
#' @keywords internal
.nh3_unionised_fraction <- function(pH, temperature_C) {
  pKa <- 0.09018 + 2729.92 / (temperature_C + 273.15)
  1 / (1 + 10^(pKa - pH))
}

#' Correct total ammonia-N to the ANZG reference pH and temperature
#'
#' The freshwater ammonia default guideline values (DGVs) and the bundled
#' ammonia SSD are expressed as **total ammonia-N at the ANZG index condition of
#' pH 7.0 and 20 °C**. Ammonia toxicity is driven by the *un-ionised* fraction
#' (NH3), which rises with both pH and temperature. A measured total ammonia-N
#' must therefore be converted to the equivalent reference-condition
#' concentration — the one that holds the same un-ionised NH3 — before it is
#' compared against a DGV or passed to [ssd_paf()] as `"NH3-N"`.
#'
#' The conversion is `C_ref = C * f_sample / f_ref`, where
#' `f = 1 / (1 + 10^(pKa - pH))` is the un-ionised fraction and
#' `pKa = 0.09018 + 2729.92 / T(K)` (Emerson et al. 1975), with
#' `T(K) = temperature_C + 273.15`. A sample at high pH/temperature is
#' normalised **upward** (more of its ammonia is in the toxic NH3 form than at
#' the reference), and a sample at low pH/temperature **downward**.
#'
#' @param conc_ug_L Measured total ammonia-N (µg/L). Numeric vector.
#' @param pH Sample pH. Recycled against `conc_ug_L`.
#' @param temperature_C Sample temperature (°C). Recycled against `conc_ug_L`.
#' @param ref_pH,ref_temperature_C Reference index condition. Defaults to the
#'   ANZG ammonia DGV basis (pH 7.0, 20 °C); change only if your SSD/DGV uses a
#'   different reference.
#' @return Total ammonia-N normalised to the reference condition (µg/L), the
#'   length of the recycled inputs. Pass this to `ssd_paf("NH3-N", ...)`.
#'
#' @section Do not double-correct:
#' [add_amspaf()] applies this same correction **automatically** from the
#' metadata `normalisation_formula` (it reads the per-sample `pH` and
#' `temperature` columns). Use this helper only for the manual [ssd_paf()]
#' path. Do **not** pre-correct with this helper and then pass the result to
#' [add_amspaf()], or ammonia will be corrected twice.
#'
#' @examples
#' # 900 µg/L total ammonia-N measured at pH 8.5, 20 °C is far more toxic than
#' # the same number at the pH 7.0 reference, so it normalises sharply upward:
#' correct_ammonia_ph_temp(900, pH = 8.5, temperature_C = 20)
#' @references
#' Emerson K, Russo RC, Lund RE, Thurston RV (1975). Aqueous ammonia
#' equilibrium calculations: effect of pH and temperature. Journal of the
#' Fisheries Research Board of Canada 32(12):2379–2383.
#' @export
correct_ammonia_ph_temp <- function(conc_ug_L, pH, temperature_C,
                                    ref_pH = 7.0, ref_temperature_C = 20) {
  if (!is.numeric(conc_ug_L) || !is.numeric(pH) || !is.numeric(temperature_C)) {
    cli::cli_abort("{.arg conc_ug_L}, {.arg pH} and {.arg temperature_C} must be numeric.")
  }
  f_sample <- .nh3_unionised_fraction(pH, temperature_C)
  f_ref    <- .nh3_unionised_fraction(ref_pH, ref_temperature_C)
  conc_ug_L * f_sample / f_ref
}
