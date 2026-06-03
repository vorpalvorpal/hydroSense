#' Pre-screen analytes by detection frequency
#'
#' Computes detection frequency (proportion of samples with `detected == TRUE`)
#' per analyte and returns the names that meet a minimum threshold. Use this
#' before `impute_chemistry()` to drop analytes that were almost never detected —
#' imputing such analytes from near-zero priors adds noise without ecological
#' signal.
#'
#' Analytes listed in `coanalytes_required` in the bundled metadata (e.g. pH,
#' DOC, Ca, Mg, hardness, temperature) are **automatically protected** from
#' exclusion regardless of detection frequency, because they are needed for
#' chemistry normalisation in [add_amspaf()].  Additional analytes can be
#' protected via the `protect` argument (typical use: pass the
#' `required_vars` you intend to use in [fit_imputation_model()] here so
#' those vars survive prescreen).  Protected analytes that fall below the
#' threshold are reported separately so the caller is aware.
#'
#' **Potency escape hatch.** Frequency alone can screen out a rare-but-potent
#' toxicant (e.g. a pesticide detected in 2 % of samples but at
#' ecotoxicologically significant concentrations). With `potency_keep = TRUE`
#' (default), an analyte that fails the frequency threshold is still kept if any
#' detected concentration reaches `potency_frac` times its 95 %-species-
#' protection guideline value (`dgv_95pct_ug_L` in the metadata). This needs a
#' numeric `value` column in the same units as the guideline (µg/L); only
#' analytes that carry a guideline value (toxicants) can be rescued this way, so
#' major ions and analytes with no guideline are unaffected.
#'
#' @param df Long-format chemistry data frame.  Required columns:
#'   `sample_id` (character), `analyte` (character), `detected` (logical).
#'   When `group_by_feature = TRUE`, also requires `site_id` (character).
#' @param k Minimum detection frequency (proportion, 0–1). Analytes with
#'   `n_detected / n_samples < k` are excluded unless protected. Default
#'   `0.05` (5 %).
#' @param protect Optional character vector of additional analyte names to
#'   protect from prescreen exclusion, on top of the metadata-derived
#'   co-analytes.  Default `NULL`.  Pass the `required_vars` from your
#'   downstream imputation step here.
#' @param potency_keep Logical. Enable the potency escape hatch (keep a
#'   frequency-failing analyte whose concentration reaches its guideline value).
#'   Default `TRUE`. Set `FALSE` for a frequency-only prescreen.
#' @param potency_frac Numeric (>= 0). Fraction of the 95 % guideline value
#'   (`dgv_95pct_ug_L`) a detected concentration must reach to rescue an
#'   analyte. Default `1` (must reach the guideline). Lower it (e.g. `0.5`) to
#'   be more precautionary.
#' @param group_by_feature Logical. If `TRUE`, detection frequency is computed
#'   per `site_id` and an analyte is included only if it passes in *all*
#'   features (worst-case feature, pooled counts). If `FALSE` (default),
#'   frequency is pooled across all samples.
#' @param analyte_metadata Data frame of analyte metadata, or `NULL` to load
#'   the bundled `inst/extdata/anzecc_analyte_metadata.csv`. Used to identify
#'   co-analytes that must be protected from exclusion.
#' @param return Either `"vector"` (default) to return a character vector of
#'   included analyte names, or `"table"` to return a tibble with one row per
#'   analyte showing detection statistics and inclusion flag.
#'
#' @return When `return = "vector"`: a character vector of passing analyte
#'   names. The vector carries an attribute `"excluded"` listing analytes that
#'   did not pass the threshold (non-protected only), so callers can record
#'   what was dropped.
#'
#'   When `return = "table"`: a tibble with columns `analyte`,
#'   `n_samples`, `n_detected`, `detect_freq`, `limiting_site`,
#'   `protected` (logical), `potency_kept` (logical; rescued by the potency
#'   escape hatch), `included` (logical).  When
#'   `group_by_feature = TRUE`, `n_samples` / `n_detected` / `detect_freq`
#'   describe the *limiting* site (the worst-case feature that drove the
#'   inclusion decision) and `limiting_site` names it; when
#'   `group_by_feature = FALSE`, counts are pooled across all samples and
#'   `limiting_site` is `NA`.
#'
#' @examples
#' included <- prescreen_analytes(leachate_demo(), k = 0.05)
#' included                    # analytes retained
#' attr(included, "excluded")  # see what was dropped (and why)
#'
#' @export
prescreen_analytes <- function(
    df,
    k = 0.05,
    protect          = NULL,
    potency_keep     = TRUE,
    potency_frac     = 1,
    group_by_feature = FALSE,
    analyte_metadata = NULL,
    return = c("vector", "table")
) {
  return <- match.arg(return)

  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "analyte", "detected"))
  checkmate::assert_number(k, lower = 0, upper = 1)
  checkmate::assert_character(protect, null.ok = TRUE, any.missing = FALSE)
  checkmate::assert_flag(potency_keep)
  checkmate::assert_number(potency_frac, lower = 0)
  checkmate::assert_flag(group_by_feature)

  if (group_by_feature) {
    checkmate::assert_names(names(df), must.include = "site_id",
      .var.name = "df (when group_by_feature = TRUE)")
  }

  if (!is.logical(df$detected)) {
    cli::cli_abort(
      "{.arg df$detected} must be logical, not {.cls {class(df$detected)}}."
    )
  }

  # ── Identify protected analytes ───────────────────────────────────────────
  # Co-analytes required by normalisation formulas are auto-protected.
  meta <- .load_analyte_metadata(analyte_metadata)
  co_req_all <- meta$coanalytes_required[!is.na(meta$coanalytes_required) &
                                           nzchar(meta$coanalytes_required)]
  coanalytes_from_meta <- unique(unlist(
    lapply(co_req_all, function(x) trimws(strsplit(x, ",")[[1L]]))
  ))
  coanalytes_from_meta <- coanalytes_from_meta[nzchar(coanalytes_from_meta)]

  # Warn if any user-supplied `protect` name does not appear in the data —
  # a common sign of a typo or a stale analyte name that would silently
  # protect nothing.
  if (!is.null(protect)) {
    protect_missing <- setdiff(protect, unique(df$analyte))
    if (length(protect_missing) > 0L) {
      cli::cli_warn(c(
        "!" = "{length(protect_missing)} {.arg protect} analyte{?s} not found \\
               in {.arg df$analyte} \u2014 protecting nothing: {.val {protect_missing}}.",
        "i" = "Check for typos or stale analyte names."
      ))
    }
  }

  protected_analytes <- unique(c(protect, coanalytes_from_meta))

  # ── Compute detection frequency ───────────────────────────────────────────
  # n_samples counts DISTINCT samples (not rows) so duplicated analyte rows in
  # a single sample don't inflate the denominator.
  if (group_by_feature) {
    per_feat <- df |>
      dplyr::group_by(.data$analyte, .data$site_id) |>
      dplyr::summarise(
        n_samples  = dplyr::n_distinct(.data$sample_id),
        n_detected = dplyr::n_distinct(.data$sample_id[.data$detected]),
        .groups    = "drop"
      ) |>
      dplyr::mutate(feat_freq = .data$n_detected / .data$n_samples)

    # An analyte's inclusion is decided by its WORST-CASE (limiting) site —
    # the site with the lowest detection frequency.  Report that site's own
    # counts (not pooled sums across sites), so n_samples / n_detected /
    # detect_freq are mutually consistent and describe the site that actually
    # drove the decision.  `limiting_site` names it.
    tbl <- per_feat |>
      dplyr::group_by(.data$analyte) |>
      dplyr::arrange(.data$feat_freq, .by_group = TRUE) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::transmute(
        .data$analyte,
        n_samples     = .data$n_samples,
        n_detected    = .data$n_detected,
        detect_freq   = .data$feat_freq,  # worst-case feature
        limiting_site = .data$site_id
      )
  } else {
    tbl <- df |>
      dplyr::group_by(.data$analyte) |>
      dplyr::summarise(
        n_samples     = dplyr::n_distinct(.data$sample_id),
        n_detected    = dplyr::n_distinct(.data$sample_id[.data$detected]),
        detect_freq   = .data$n_detected / .data$n_samples,
        limiting_site = NA_character_,  # schema parity with group_by_feature
        .groups       = "drop"
      )
  }

  # ── Potency-based escape hatch ────────────────────────────────────────────
  # A rare analyte (below k) is still kept if any detected concentration
  # reaches `potency_frac` x its 95%-species-protection guideline value
  # (`dgv_95pct_ug_L` in the metadata). This protects rare-but-potent
  # toxicants from being screened out by frequency alone. It needs a numeric
  # `value` column in ug/L (matching the DGV units); only analytes that carry
  # a DGV (i.e. toxicants) can be rescued this way.
  do_potency <- potency_keep
  if (potency_keep && !"value" %in% names(df)) {
    cli::cli_warn(c(
      "!" = "The potency escape hatch needs a {.field value} column; none \\
             found, so using a frequency-only prescreen.",
      "i" = "Pass {.code potency_keep = FALSE} to silence this."
    ))
    do_potency <- FALSE
  }

  if (do_potency) {
    max_conc <- df |>
      dplyr::filter(.data$detected, !is.na(.data$value)) |>
      dplyr::group_by(.data$analyte) |>
      dplyr::summarise(max_conc = max(.data$value), .groups = "drop")
    dgv <- dplyr::select(meta, "analyte", dgv = "dgv_95pct_ug_L")
    tbl <- tbl |>
      dplyr::left_join(max_conc, by = "analyte") |>
      dplyr::left_join(dgv, by = "analyte")
  } else {
    tbl <- dplyr::mutate(tbl, max_conc = NA_real_, dgv = NA_real_)
  }

  tbl <- dplyr::mutate(
    tbl,
    protected    = .data$analyte %in% .env$protected_analytes,
    potency_kept = do_potency & !.data$protected & .data$detect_freq < k &
                   !is.na(.data$dgv) & !is.na(.data$max_conc) &
                   .data$max_conc >= potency_frac * .data$dgv,
    included     = .data$detect_freq >= k | .data$protected | .data$potency_kept
  )

  n_excluded  <- sum(!tbl$included)
  n_potency_kept  <- sum(tbl$potency_kept)
  n_protected_low <- sum(tbl$protected & tbl$detect_freq < k)

  if (n_excluded > 0) {
    cli::cli_inform(c(
      "i" = "prescreen_analytes: {n_excluded} analyte{?s} excluded \\
             (below k = {k} detection frequency): \\
             {.val {tbl$analyte[!tbl$included]}}."
    ))
  }
  if (n_protected_low > 0) {
    cli::cli_inform(c(
      "i" = "prescreen_analytes: {n_protected_low} protected analyte{?s} \\
             below k = {k} threshold but kept (required for normalisation): \\
             {.val {tbl$analyte[tbl$protected & tbl$detect_freq < k]}}."
    ))
  }
  if (n_potency_kept > 0) {
    cli::cli_inform(c(
      "i" = "prescreen_analytes: {n_potency_kept} rare analyte{?s} below k = {k} \\
             but kept (a detected concentration reaches {potency_frac}x the \\
             95% guideline value): {.val {tbl$analyte[tbl$potency_kept]}}."
    ))
  }

  # Internal-only columns used to decide the potency rescue.
  tbl <- dplyr::select(tbl, -dplyr::any_of(c("max_conc", "dgv")))

  if (return == "table") {
    return(tbl)
  }

  included_names <- tbl$analyte[tbl$included]
  excluded_names <- tbl$analyte[!tbl$included]
  structure(included_names, excluded = excluded_names)
}
