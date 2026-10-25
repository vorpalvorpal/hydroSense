# -- Analyte-name aliases and ion chemistry ------------------------------------
#
# hydroSense computes with its own analyte names ("hardness", "temperature",
# "CO3-CaCO3", ...). sampleTidy stores the same quantities under different
# names ("Hardness-total-CaCO3", "Temperature", "Alkalinity-total-CaCO3", ...)
# and carries a character `units` column rather than `units.analyte`. The
# helpers here let both spellings through the public entry points without the
# caller renaming anything: inputs are canonicalised on the way in, and the
# caller's own rows (with their original names) are what come back out.

#' Analyte-name aliases
#'
#' The mapping from the analyte names used by sampleTidy (and other upstream
#' stores) to the names hydroSense computes with. [add_lmf()], [to_meq()],
#' [add_mspaf()] and [mspaf_pipeline()] apply it to their inputs, so data
#' stored under either spelling is accepted, and return the caller's rows
#' with their original names.
#'
#' Two upstream quantities have no hydroSense counterpart of their own and map
#' to themselves: `"Alkalinity-total-CaCO3"` (the sum of `"CO3-CaCO3"` and
#' `"HCO3-CaCO3"`) and `"NO2+NO3-N"` (the sum of `"NO3-N"` and `"NO2-N"`).
#' [add_lmf()] uses such a combined value in place of its two species, and
#' ignores it when the species themselves are present, so nothing is counted
#' twice.
#'
#' @return A data frame with character columns `alias` (an upstream name) and
#'   `canonical` (the hydroSense name), one row per alias. Names that already
#'   match map to themselves.
#'
#' @examples
#' analyte_aliases()
#' @export
analyte_aliases <- function() {
  so4 <- "SO4\u00b2\u207b"
  tibble::tibble(
    alias = c(
      "Na", "K", "Ca", "Mg", "Cl", so4, "F",
      "NH3-N", "NO3-N", "NO2-N", "CO3-CaCO3", "HCO3-CaCO3",
      "NO2+NO3-N", "Alkalinity-total-CaCO3",
      "pH", "EC", "DOC",
      "Hardness-total-CaCO3", "hardness",
      "Temperature", "temperature"
    ),
    canonical = c(
      "Na", "K", "Ca", "Mg", "Cl", so4, "F",
      "NH3-N", "NO3-N", "NO2-N", "CO3-CaCO3", "HCO3-CaCO3",
      "NO2+NO3-N", "Alkalinity-total-CaCO3",
      "pH", "EC", "DOC",
      "hardness", "hardness",
      "temperature", "temperature"
    )
  )
}

## Rename aliased analytes to their canonical hydroSense names.
.canonicalise_analytes <- function(df) {
  al  <- analyte_aliases()
  hit <- match(df$analyte, al$alias)
  df$analyte <- dplyr::if_else(is.na(hit), df$analyte, al$canonical[hit])
  df
}

## Unit strings as udunits parses them: the micro sign (U+00B5 or U+03BC)
## becomes "u".
.normalise_units <- function(x) {
  gsub("[\u00b5\u03bc]", "u", as.character(x))
}

## The name of the per-row unit column of `df`: legacy `units.analyte` wins
## over sampleTidy's `units`; NA when neither is present.
.units_col <- function(df) {
  if ("units.analyte" %in% names(df)) {
    "units.analyte"
  } else if ("units" %in% names(df)) {
    "units"
  } else {
    NA_character_
  }
}

## Canonicalise an msPAF-bound chemistry frame. Row-preserving: analytes are
## renamed and, when the frame has only sampleTidy's `units`, a `units.analyte`
## column is added from it. Idempotent.
.canonicalise_mspaf_input <- function(df) {
  df <- .canonicalise_analytes(df)
  if (!"units.analyte" %in% names(df) && "units" %in% names(df)) {
    df$units.analyte <- .normalise_units(df$units)
  }
  df
}

## Canonicalise a long-format frame for the LMF panel. Renames aliases and
## resolves the two combined quantities: Alkalinity-total-CaCO3 stands in for
## CO3-CaCO3 + HCO3-CaCO3 (both are expressed as CaCO3, so meq/L is identical
## whichever of the two it is booked against), and NO2+NO3-N for NO3-N +
## NO2-N (both as N). Per sample: when ALL the species a combined row stands
## for are also present, the combined row is dropped (it would double count);
## otherwise the combined value is used and any partial species are dropped.
.canonicalise_lmf_input <- function(df) {
  df <- .canonicalise_analytes(df)
  combos <- list(
    list(combined = "Alkalinity-total-CaCO3",
         parts = c("CO3-CaCO3", "HCO3-CaCO3"), as = "HCO3-CaCO3"),
    list(combined = "NO2+NO3-N",
         parts = c("NO3-N", "NO2-N"), as = "NO3-N")
  )
  for (cb in combos) {
    is_comb <- df$analyte == cb$combined
    if (!any(is_comb)) next
    ## Samples carrying every species the combined row stands for.
    n_parts <- tapply(df$analyte %in% cb$parts & !duplicated(df[c("sample_id", "analyte")]),
                      df$sample_id, sum)
    full <- df$sample_id %in% names(n_parts)[n_parts == length(cb$parts)]
    comb_ids <- unique(df$sample_id[is_comb])
    drop_part <- df$analyte %in% cb$parts & df$sample_id %in% comb_ids & !full
    df <- df[!(is_comb & full) & !drop_part, , drop = FALSE]
    df$analyte[df$analyte == cb$combined] <- cb$as
  }
  df
}

## Valence and molar mass of the LMF panel ions, keyed by canonical name.
## Nitrogen species are expressed as N and the carbonate species as CaCO3,
## matching how the analytes are reported. Molar masses in g/mol.
.ion_properties <- function() {
  tibble::tibble(
    analyte = c("Na", "K", "Ca", "Mg", "Cl", "SO4\u00b2\u207b", "F",
                "NH3-N", "NO3-N", "NO2-N", "CO3-CaCO3", "HCO3-CaCO3"),
    valence = c(1, 1, 2, 2, 1, 2, 1, 1, 1, 1, 2, 2),
    atomic_mass = c(22.98977, 39.0983, 40.078, 24.305, 35.453, 96.06, 18.998403,
                    14.0067, 14.0067, 14.0067, 100.0869, 100.0869)
  )
}
