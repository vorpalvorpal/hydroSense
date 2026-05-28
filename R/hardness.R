# ── derive_hardness() ─────────────────────────────────────────────────────────

#' Three-way reconciliation between Ca, Mg, and hardness
#'
#' Water hardness (total, expressed as CaCO3-equivalent mg/L) is stoichiometric
#' in Ca and Mg:
#'
#' \deqn{hardness = 2.497 \cdot Ca + 4.118 \cdot Mg}
#'
#' where Ca and Mg are in mg/L.  This helper fills in any missing member of
#' \{Ca, Mg, hardness\} when the other two are available, and warns if all
#' three are present but inconsistent.
#'
#' **Per-sample logic:**
#'
#' \describe{
#'   \item{All three present}{Check consistency.  If
#'     \eqn{|\hat{hardness} - hardness| / hardness > tolerance}, emit a
#'     per-sample warning and keep the user-supplied values.}
#'   \item{Exactly two present}{Compute the third exactly from stoichiometry
#'     and append as a new row with \code{detected = TRUE}.}
#'   \item{One or none present}{Leave alone — fill via imputation if needed.}
#' }
#'
#' **Recommended pipeline use:**
#'
#' Call twice — once **before** imputation to fill samples where the third
#' member can be derived from raw measurements, and again **after**
#' `impute_coanalytes()` to fill hardness for samples whose Ca and Mg were
#' just imputed.  Idempotent if all three are already consistent.
#'
#' @param df Long-format chemistry data frame with columns
#'   `sample_id`, `analyte`, `value`, `detected`.  Required columns
#'   `site_id` and `datetime` are propagated to new rows from the existing
#'   per-sample metadata.
#' @param tolerance Relative tolerance (proportion of measured hardness) for
#'   consistency check when all three are present.  Default `0.05` (5%).
#' @param verbose Logical.  If `TRUE`, prints a summary of how many rows were
#'   derived and how many inconsistencies were detected.  Default `TRUE`.
#'
#' @return The input `df` with derived rows appended for `Ca`, `Mg`, or
#'   `hardness` wherever exactly two of the three were available.  Derived
#'   rows are tagged with `imputed = TRUE` and `imputed_kind = "derived"`
#'   if those columns exist on the input (otherwise they are added).
#'
#' @examples
#' \dontrun{
#' # Fill missing hardness where Ca and Mg are both measured
#' chem2 <- derive_hardness(chem)
#'
#' # Call twice in a typical pipeline: pre- and post-imputation
#' chem      <- derive_hardness(chem)
#' chem_imp  <- impute_chemistry(chem, model)
#' chem_imp2 <- impute_coanalytes(chem_imp, model)
#' chem_imp3 <- derive_hardness(chem_imp2)
#' }
#'
#' @export
derive_hardness <- function(df, tolerance = 0.05, verbose = TRUE) {
  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "analyte", "value", "detected"))
  checkmate::assert_number(tolerance, lower = 0, upper = 1)
  checkmate::assert_flag(verbose)

  # Conversion factors: hardness (mg/L as CaCO3) = 2.497*Ca + 4.118*Mg (mg/L)
  k_Ca <- 2.497
  k_Mg <- 4.118

  # Pivot the three analytes wide per sample
  trio <- df |>
    dplyr::filter(.data$analyte %in% c("Ca", "Mg", "hardness"),
                  .data$detected) |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select("sample_id", "analyte", "value") |>
    tidyr::pivot_wider(names_from = "analyte", values_from = "value")

  # Ensure all three columns exist even if no observations
  for (col in c("Ca", "Mg", "hardness")) {
    if (!col %in% names(trio)) trio[[col]] <- NA_real_
  }

  # Classify per sample
  trio <- dplyr::mutate(
    trio,
    .case = dplyr::case_when(
      !is.na(.data$Ca) & !is.na(.data$Mg) & !is.na(.data$hardness) ~ "all_three",
      !is.na(.data$Ca) & !is.na(.data$Mg)                          ~ "derive_hardness",
      !is.na(.data$Ca) & !is.na(.data$hardness)                    ~ "derive_Mg",
      !is.na(.data$Mg) & !is.na(.data$hardness)                    ~ "derive_Ca",
      TRUE                                                          ~ "skip"
    )
  )

  # Per-sample metadata for constructing new rows (site_id, datetime, ...)
  meta_cols <- intersect(c("site_id", "datetime", "focal_date"), names(df))
  if (length(meta_cols) > 0L) {
    sample_meta <- df |>
      dplyr::group_by(.data$sample_id) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::select("sample_id", dplyr::all_of(meta_cols))
  } else {
    sample_meta <- dplyr::distinct(df, .data$sample_id)
  }

  # ── Consistency check on all-three samples ────────────────────────────────
  all3 <- dplyr::filter(trio, .data$.case == "all_three") |>
    dplyr::mutate(
      hardness_calc = k_Ca * .data$Ca + k_Mg * .data$Mg,
      rel_err       = abs(.data$hardness_calc - .data$hardness) / .data$hardness,
      inconsistent  = .data$rel_err > tolerance
    )

  n_inconsistent <- sum(all3$inconsistent, na.rm = TRUE)
  if (n_inconsistent > 0L) {
    bad <- dplyr::filter(all3, .data$inconsistent) |>
      dplyr::mutate(msg = sprintf(
        "  %s: hardness=%.2f, 2.497*Ca + 4.118*Mg = %.2f (rel err %.1f%%)",
        .data$sample_id, .data$hardness, .data$hardness_calc,
        100 * .data$rel_err
      ))
    cli::cli_warn(c(
      "!" = "{n_inconsistent} sample{?s} have Ca/Mg/hardness disagreement \\
             > {100*tolerance}%:",
      bad$msg,
      "i" = "User-supplied values kept as-is.  Common causes: hardness \\
             reported as Ca-only hardness, units mismatch, analytical bias \\
             on one component."
    ))
  }

  # ── Construct derived rows for the three derivation cases ────────────────
  new_rows <- list()

  derive_one <- function(case_name, analyte_nm, value_expr) {
    src <- dplyr::filter(trio, .data$.case == case_name)
    if (nrow(src) == 0L) return(NULL)
    src |>
      dplyr::transmute(
        .data$sample_id,
        analyte      = analyte_nm,
        value        = value_expr(src),
        detected     = TRUE,
        imputed      = TRUE,
        imputed_kind = "derived"
      ) |>
      dplyr::left_join(sample_meta, by = "sample_id")
  }

  new_rows$hardness <- derive_one("derive_hardness", "hardness",
    function(s) k_Ca * s$Ca + k_Mg * s$Mg)

  new_rows$Mg <- derive_one("derive_Mg", "Mg",
    function(s) (s$hardness - k_Ca * s$Ca) / k_Mg)

  new_rows$Ca <- derive_one("derive_Ca", "Ca",
    function(s) (s$hardness - k_Mg * s$Mg) / k_Ca)

  # Drop negative-Ca / negative-Mg derivations (occur when hardness is
  # inconsistent with the supplied Ca or Mg)
  for (nm in c("Mg", "Ca")) {
    if (!is.null(new_rows[[nm]])) {
      neg <- new_rows[[nm]]$value < 0
      if (any(neg)) {
        bad_ids <- new_rows[[nm]]$sample_id[neg]
        cli::cli_warn(c(
          "!" = "{sum(neg)} sample{?s} produced negative derived {nm} \\
                 (hardness inconsistent with supplied {.val {setdiff(c('Ca','Mg'), nm)}}). \\
                 Dropping these rows: {.val {bad_ids}}."
        ))
        new_rows[[nm]] <- new_rows[[nm]][!neg, , drop = FALSE]
      }
    }
  }

  derived <- dplyr::bind_rows(new_rows)
  n_derived <- nrow(derived)

  if (verbose && (n_derived > 0L || n_inconsistent > 0L)) {
    counts <- if (n_derived > 0L) table(derived$analyte) else integer(0)
    cli::cli_inform(c(
      "i" = "derive_hardness: {n_derived} row{?s} derived \\
             ({paste(sprintf('%s=%d', names(counts), as.integer(counts)), \\
                     collapse=', ')}); \\
             {n_inconsistent} sample{?s} flagged inconsistent."
    ))
  }

  # ── Ensure imputed/imputed_kind columns exist on the original df ─────────
  if (n_derived > 0L) {
    if (!"imputed" %in% names(df)) {
      df <- dplyr::mutate(df, imputed = FALSE, imputed_kind = "observed")
    }
    dplyr::bind_rows(df, derived) |>
      dplyr::arrange(.data$sample_id, .data$analyte)
  } else {
    df
  }
}
