#' Prepare reference chemistry for AmsPAF background subtraction
#'
#' Applies chemistry normalisation (if formulas are populated in the metadata)
#' and computes per-analyte quantile concentrations from reference-site data.
#' The resulting object is passed as the `reference` argument to
#' `add_amspaf()`.
#'
#' This is a pure function — it has no side effects and no internal cache.
#' In the chronic pipeline, call it once after computing chronic chemistry for
#' the reference feature(s), then pass the object into `add_amspaf()` for
#' every focal date. In the per-sample pipeline, call it once on the raw
#' reference chemistry.
#'
#' @param reference_data Long-format chemistry data frame for the reference
#'   (background) site(s). Same schema as the input to `add_amspaf()`:
#'   `name.analyte`, `value`, `quantified`. If `quantified == FALSE` for a row,
#'   it is excluded from the quantile calculation (BDL observations are treated
#'   as absent at the reference site).
#' @param analyte_metadata Data frame of analyte metadata, or `NULL` to load
#'   the bundled `inst/extdata/anzecc_analyte_metadata.csv`. Must contain
#'   columns `analyte`, `coanalytes_required`, and `normalisation_formula`.
#' @param percentile Quantile (0–1) of the reference distribution used as the
#'   background anchor for Added Risk Approach (ARA) subtraction. Default
#'   `0.80` (80th percentile — the concentration the reference site does not
#'   exceed 80 % of the time).
#'
#' @return A list of class `"prepared_reference"` with elements:
#'   - `$normalised_quantiles`: tibble with columns `name.analyte`,
#'     `ref_norm` (normalised concentration at `percentile`).
#'   - `$dropped`: character vector of analytes excluded due to no quantified
#'     reference observations.
#'   - `$percentile`: the `percentile` value used.
#'
#' @examples
#' \dontrun{
#' # Per-sample pipeline
#' prep_ref <- prepare_reference(ref_df)
#' out <- add_amspaf(sample_df, reference = prep_ref)
#'
#' # Chronic pipeline: integrate reference chemistry first
#' chr_ref  <- compute_chronic_chemistry(ref_df, focal_dates = bio_dates)
#' prep_chr <- prepare_reference(chr_ref)
#' out_chr  <- add_amspaf(chr_chem, reference = prep_chr)
#' }
#'
#' @export
prepare_reference <- function(
    reference_data,
    analyte_metadata = NULL,
    percentile       = 0.80
) {
  checkmate::assert_data_frame(reference_data)
  checkmate::assert_names(names(reference_data),
    must.include = c("name.analyte", "value", "quantified"))
  checkmate::assert_number(percentile, lower = 0, upper = 1)

  meta <- .load_analyte_metadata(analyte_metadata)

  # Keep only quantified reference observations
  ref_q <- dplyr::filter(reference_data, .data$quantified)

  if (nrow(ref_q) == 0L) {
    cli::cli_warn(
      "No quantified reference observations found; ARA subtraction will be zero \\
       for all analytes."
    )
    return(structure(
      list(
        normalised_quantiles = tibble::tibble(
          name.analyte = character(0),
          ref_norm     = numeric(0)
        ),
        dropped    = character(0),
        percentile = percentile
      ),
      class = "prepared_reference"
    ))
  }

  # Apply normalisation per analyte row
  nq <- ref_q |>
    dplyr::left_join(
      meta |>
        dplyr::select("analyte", "coanalytes_required", "normalisation_formula"),
      by = c("name.analyte" = "analyte")
    ) |>
    dplyr::mutate(
      value_norm = purrr::pmap_dbl(
        list(
          formula_str = .data$normalisation_formula,
          C           = .data$value,
          sample_uid  = .data$uuid.sample %||% NA_character_
        ),
        function(formula_str, C, sample_uid) {
          parsed <- .parse_normalisation_formula(formula_str %||% "")
          if (is.null(parsed)) return(C)  # identity
          # For reference data, co-analytes are looked up from same sample
          # (simplified: extract from the full reference_data by uuid.sample)
          coanalytes <- .extract_coanalytes_for_sample(
            reference_data, sample_uid,
            meta$coanalytes_required[meta$analyte == .data$name.analyte]
          )
          .apply_normalisation(parsed, C, coanalytes)
        }
      )
    )

  # Per-analyte quantile of normalised concentrations
  qnt <- nq |>
    dplyr::group_by(.data$name.analyte) |>
    dplyr::summarise(
      ref_norm = quantile(.data$value_norm, probs = percentile, na.rm = TRUE),
      n_obs    = sum(!is.na(.data$value_norm)),
      .groups  = "drop"
    )

  # Analytes that ended up with no usable reference observations
  all_analytes <- unique(reference_data$name.analyte)
  dropped      <- setdiff(all_analytes, qnt$name.analyte[qnt$n_obs > 0L])

  structure(
    list(
      normalised_quantiles = dplyr::select(qnt, "name.analyte", "ref_norm"),
      dropped    = dropped,
      percentile = percentile
    ),
    class = "prepared_reference"
  )
}

#' @export
print.prepared_reference <- function(x, ...) {
  cat(sprintf(
    "<prepared_reference>  %d analytes | %gth percentile | %d dropped\n",
    nrow(x$normalised_quantiles),
    x$percentile * 100,
    length(x$dropped)
  ))
  invisible(x)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Load and lightly validate the analyte metadata table
#'
#' Reads the bundled CSV when `meta` is NULL. Returns a tibble. Caches the
#' result in an internal environment to avoid re-reading on every
#' `prepare_reference()` / `add_amspaf()` call within a session.
#' @keywords internal
.meta_cache_env <- new.env(parent = emptyenv())

.load_analyte_metadata <- function(meta = NULL) {
  if (!is.null(meta)) {
    checkmate::assert_data_frame(meta)
    checkmate::assert_names(
      names(meta),
      must.include = c("analyte", "coanalytes_required", "normalisation_formula")
    )
    return(meta)
  }

  if (exists("meta", envir = .meta_cache_env, inherits = FALSE)) {
    return(get("meta", envir = .meta_cache_env, inherits = FALSE))
  }

  path <- system.file("extdata", "anzecc_analyte_metadata.csv",
                      package = "leachatetools")
  if (!nzchar(path)) {
    cli::cli_abort(
      "Cannot find {.file inst/extdata/anzecc_analyte_metadata.csv} inside the \\
       installed leachatetools package. Re-install or supply {.arg analyte_metadata} \\
       explicitly."
    )
  }

  m <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(
      coanalytes_required  = dplyr::if_else(
        is.na(.data$coanalytes_required), "", .data$coanalytes_required),
      normalisation_formula = dplyr::if_else(
        is.na(.data$normalisation_formula), "", .data$normalisation_formula)
    )

  assign("meta", m, envir = .meta_cache_env)
  m
}

#' Extract co-analyte values for a given sample from a long-format df
#'
#' @param df Long-format chemistry df
#' @param sample_uid uuid.sample identifier (may be NA for non-sample data)
#' @param coanalytes_str Comma-separated string of required co-analyte names
#' @return Named numeric vector; may be empty
#' @keywords internal
.extract_coanalytes_for_sample <- function(df, sample_uid, coanalytes_str) {
  if (is.na(sample_uid) || !nzchar(coanalytes_str %||% "")) {
    return(numeric(0))
  }
  required <- trimws(strsplit(coanalytes_str, ",")[[1L]])
  required <- required[nzchar(required)]
  if (length(required) == 0L) return(numeric(0))

  co <- dplyr::filter(df,
    .data$uuid.sample == .env$sample_uid,
    .data$name.analyte %in% .env$required,
    .data$quantified
  )
  if (nrow(co) == 0L) return(numeric(0))
  vals <- co$value
  names(vals) <- co$name.analyte
  vals[required[required %in% names(vals)]]
}

# Null-coalescing operator (base R doesn't have one before 4.4)
`%||%` <- function(x, y) if (is.null(x)) y else x
