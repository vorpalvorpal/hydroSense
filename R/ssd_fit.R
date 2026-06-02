# ssd_fit.R — Model-fitting infrastructure for paf.R
#
# This file provides a single internal function, .fit_for_analyte(), which
# paf.R calls on cache-miss.  It handles both data-source families:
#
#   ANZG_XLSX  — modern ANZG technical-brief data, read directly from the
#                guideline XLSX files.  Requires guideline_dir to point to the
#                folder containing the ANZG "guideline data/" XLSX files.
#                Set via: options(leachatetools.guideline_dir = "/path/to/...")
#
#   Warne2000  — ANZECC 2000 analytes; raw data from the package-bundled
#                inst/extdata/anzecc_warne2000_observations.csv.
#
# Special cases:
#   Ni         — reads MLR-normalised values from ni_mlr_normalised_table3.csv
#                (bundled in inst/extdata/), not from the XLSX directly.
#   NO3-N_soft/mod/hard — three separate entries read from different sheets of
#                the same nitrate XLSX.

# ── XLSX reading helpers ───────────────────────────────────────────────────────

.read_col_by_index <- function(path, sheet, skip, species_col, conc_col,
                                data_row_offset = 1, units_factor = 1,
                                media_col = NULL, freshwater_only = FALSE) {
  raw  <- readxl::read_excel(path, sheet = sheet, col_names = FALSE, skip = skip)
  data <- raw[seq(data_row_offset + 1, nrow(raw)), ]

  get_col <- function(df, idx) suppressWarnings(as.character(unlist(df[, idx])))
  species <- get_col(data, species_col)
  conc    <- suppressWarnings(as.numeric(get_col(data, conc_col))) * units_factor
  df <- tibble::tibble(Conc = conc, Species = species)

  if (freshwater_only && !is.null(media_col)) {
    media <- get_col(data, media_col)
    df    <- df[grepl("fresh", media, ignore.case = TRUE), ]
  }
  df |> dplyr::filter(!is.na(Conc), Conc > 0, !is.na(Species), Species != "NA")
}

.read_min_per_species <- function(path, sheet, skip, species_col, conc_col,
                                   data_row_offset = 1, units_factor = 1,
                                   media_col = NULL, freshwater_only = FALSE) {
  df <- .read_col_by_index(path, sheet, skip, species_col, conc_col,
                            data_row_offset, units_factor, media_col, freshwater_only)
  df |>
    dplyr::group_by(Species) |>
    dplyr::summarise(Conc = min(Conc), .groups = "drop")
}

# ── Per-analyte XLSX reader registry ──────────────────────────────────────────
#
# Each entry is a one-argument function(guideline_dir) that returns a
# tibble(Conc, Species).  Ni is handled separately (reads from extdata CSV).

.XLSX_READERS <- list(

  `NH3-N` = function(gdir) {
    .read_col_by_index(
      file.path(gdir, "ammonia-fresh-dgvs-data-entry.xlsx"),
      sheet = 1, skip = 3, species_col = 4, conc_col = 55, units_factor = 1000)
  },

  `NO3-N_soft` = function(gdir) {
    .read_col_by_index(
      file.path(gdir, "nitrate-fresh-dgvs-data-entry.xlsx"),
      sheet = "Nitrate - soft water", skip = 9, species_col = 5,
      conc_col = 49, data_row_offset = 2, units_factor = 1000)
  },

  `NO3-N_mod` = function(gdir) {
    .read_col_by_index(
      file.path(gdir, "nitrate-fresh-dgvs-data-entry.xlsx"),
      sheet = "Nitrate - moderately hard water", skip = 9, species_col = 5,
      conc_col = 49, data_row_offset = 2, units_factor = 1000)
  },

  `NO3-N_hard` = function(gdir) {
    .read_col_by_index(
      file.path(gdir, "nitrate-fresh-dgvs-data-entry.xlsx"),
      sheet = "Nitrate - hard water", skip = 9, species_col = 5,
      conc_col = 49, data_row_offset = 2, units_factor = 1000)
  },

  `B` = function(gdir) {
    .read_col_by_index(
      file.path(gdir, "boron_fresh_dgv_data-entry_final.xlsx"),
      sheet = 1, skip = 6, species_col = 5, conc_col = 49, units_factor = 1000)
  },

  `Cr` = function(gdir) {
    .read_col_by_index(
      file.path(gdir, "chromium-III-fresh-dgvs-data-entry.xlsx"),
      sheet = 1, skip = 6, species_col = 5, conc_col = 48,
      data_row_offset = 2, units_factor = 1, media_col = 4,
      freshwater_only = TRUE)
  },

  `Cu` = function(gdir) {
    .read_col_by_index(
      file.path(gdir, "copper-fresh-dgvs-data-entry.xlsx"),
      sheet = "Accepted Chronic Data", skip = 5,
      species_col = 7, conc_col = 57, units_factor = 1)
  },

  # Ni: reads from bundled pre-normalised CSV rather than the XLSX.
  # The XLSX col 31 holds raw measured concentrations; the 26 MLR-normalised
  # negligible-effect values from Stauber et al. 2021 Table 3 are stored in
  # inst/extdata/ni_mlr_normalised_table3.csv.
  `Ni` = function(gdir) {
    csv_path <- system.file("extdata", "ni_mlr_normalised_table3.csv",
                             package = "leachatetools")
    readr::read_csv(
      csv_path,
      col_types = readr::cols(Conc_ug_L = readr::col_double(),
                              .default  = readr::col_character()),
      show_col_types = FALSE
    ) |>
      dplyr::transmute(Conc = Conc_ug_L, Species = Species)
  },

  `Zn` = function(gdir) {
    .read_min_per_species(
      file.path(gdir, "zinc-fresh-dgvs-data-entry.xlsx"),
      sheet = "ForWordDoc", skip = 3, species_col = 2, conc_col = 11)
  }
)

# ── Main fitting function ──────────────────────────────────────────────────────

#' Internal: load data and fit an SSD for one analyte.
#'
#' Called by `.load_or_fit()` in paf.R on a cache miss.
#'
#' @param analyte     Character. Canonical analyte name.
#' @param stem        Character. Safe file stem (from .SSD_NAME_MAP).
#' @param meta        One-row data.frame from anzecc_analyte_metadata.csv.
#' @param dists       Character vector. Distribution names for ssd_fit_dists().
#' @param guideline_dir Character. Path to the "guideline data/" folder.
#'   Only required for ANZG_XLSX analytes; may be NULL for Warne2000 analytes.
#'
#' @return A fitted ssdtools object with provenance attributes, or NULL on error.
#' @keywords internal
.fit_for_analyte <- function(analyte, stem, meta, dists, guideline_dir) {

  # ── Load raw data ──────────────────────────────────────────────────────────

  df <- if (meta$data_source == "ANZG_XLSX") {

    # Priority order:
    #   1. XLSX files (freshest data; requires guideline_dir to be set)
    #   2. Bundled anzg_xlsx_observations.csv (ships with the package; no
    #      external files needed; used automatically when guideline_dir is unset
    #      or the XLSX file is missing)
    #
    # The bundled CSV was extracted from the same XLSX files and contains
    # identical data.  It exists so the package works out-of-the-box without
    # requiring users to download ANZG guideline XLSX files.

    use_xlsx <- !is.null(guideline_dir) && nzchar(guideline_dir)

    if (use_xlsx) {
      reader <- .XLSX_READERS[[analyte]]
      if (is.null(reader)) {
        warning("No XLSX reader registered for '", analyte, "' \u2014 falling back to bundled CSV.")
        use_xlsx <- FALSE
      }
    }

    if (use_xlsx) {
      result <- tryCatch(
        .XLSX_READERS[[analyte]](guideline_dir),
        error = function(e) {
          cli::cli_warn(c(
            "!" = "Error reading XLSX data for {.val {analyte}}: {conditionMessage(e)}",
            "i" = "Falling back to bundled {.file anzg_xlsx_observations.csv}."
          ))
          NULL
        }
      )
      if (is.null(result)) use_xlsx <- FALSE else result
    }

    if (!use_xlsx) {
      # Fallback: read from the bundled CSV that ships with the package
      csv_path <- system.file("extdata", "anzg_xlsx_observations.csv",
                               package = "leachatetools")
      if (!nzchar(csv_path)) {
        stop("Cannot find bundled anzg_xlsx_observations.csv inside the ",
             "installed leachatetools package. Re-install the package.")
      }
      obs_all <- readr::read_csv(
        csv_path,
        col_types = readr::cols(value_ug_L = readr::col_double(),
                                .default   = readr::col_character()),
        show_col_types = FALSE
      )
      obs_all |>
        dplyr::filter(analyte == !!analyte,
                      !is.na(value_ug_L), value_ug_L > 0) |>
        dplyr::transmute(Conc = value_ug_L, Species = species_id)
    }

  } else {
    # Warne2000 (or any other CSV-backed source)
    obs_path <- system.file("extdata", "anzecc_warne2000_observations.csv",
                             package = "leachatetools")
    obs_all <- readr::read_csv(
      obs_path,
      col_types = readr::cols(value_ug_L = readr::col_double(),
                              .default   = readr::col_character()),
      show_col_types = FALSE
    )

    obs_key <- meta$observations_analyte
    obs_all |>
      dplyr::filter(analyte == obs_key,
                    !is.na(value_ug_L), value_ug_L > 0) |>
      dplyr::transmute(Conc = value_ug_L, Species = species_id)
  }

  df <- df |> dplyr::filter(!is.na(Conc), Conc > 0)

  # ── Check minimum n ────────────────────────────────────────────────────────

  n_min <- if (!is.na(meta$n_min_override)) as.integer(meta$n_min_override) else 6L
  if (nrow(df) < n_min) {
    warning("Analyte '", analyte, "': only ", nrow(df), " observations (min ",
            n_min, ") \u2014 cannot fit model.")
    return(NULL)
  }

  # ── ACR adjustment (MR analytes) ──────────────────────────────────────────
  # For analytes derived from acute LC50/EC50 data, divide by ACR so the
  # fitted SSD is on a chronic-equivalent concentration scale.

  is_acute  <- isTRUE(meta$anzecc_data_type %in% c("acute_LC50", "acute_EC50"))
  acr_val   <- suppressWarnings(as.numeric(meta$trigger_divisor))
  acute_flag <- is_acute && !is.na(acr_val) && acr_val > 1

  if (acute_flag) {
    df <- df |> dplyr::mutate(Conc = Conc / acr_val)
  }

  # ── Fit ────────────────────────────────────────────────────────────────────

  fit <- tryCatch(
    ssdtools::ssd_fit_dists(df, left = "Conc", dists = dists),
    error = function(e) {
      stop("ssd_fit_dists failed for '", analyte, "': ", conditionMessage(e))
    }
  )

  # ── Attach provenance attributes ──────────────────────────────────────────

  attr(fit, "analyte")          <- analyte
  attr(fit, "data_source")      <- meta$data_source
  attr(fit, "dists_used")       <- paste(dists, collapse = "+")
  attr(fit, "n_obs")            <- nrow(df)
  attr(fit, "acute_data")       <- is_acute
  attr(fit, "acr")              <- suppressWarnings(as.numeric(meta$acr))
  attr(fit, "trigger_divisor")  <- acr_val
  attr(fit, "acr_applied")      <- acute_flag
  attr(fit, "anzecc_hc_pct")    <- suppressWarnings(as.integer(meta$anzecc_hc_percentile))
  attr(fit, "source_note")      <- as.character(meta$notes)

  fit
}
