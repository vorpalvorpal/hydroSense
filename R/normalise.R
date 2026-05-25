# Internal helpers for chemistry normalisation
#
# Normalisation converts a measured field concentration to the ANZG index
# condition before SSD lookup, e.g. Cu at DOC=0.5 mg/L or Zn at
# hardness=30 mg/L. The formulas are stored as R-expression strings in the
# `normalisation_formula` column of `inst/extdata/anzecc_analyte_metadata.csv`.
#
# Day-1 state: all formula cells are empty, so `.apply_normalisation()` returns
# the raw concentration unchanged (identity). Populate `normalisation_formula`
# and `coanalytes_required` for Cu, Ni, Zn, NH3-N etc. once the ANZG equations
# are extracted from the technical briefs.

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
