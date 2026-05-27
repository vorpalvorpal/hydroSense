# run_amspaf_pipeline.R ──────────────────────────────────────────────────────
#
# Standalone helper script: reads chemistry from a monitoring DuckDB (same
# schema as leachatetools/test data/monitoring.duckdb), runs the full
# chronic AmsPAF pipeline, and writes results back to the DB.
#
# This script is NOT part of the leachatetools package — it is intended to be
# dropped next to a DB copy and run as needed.  It calls the leachatetools
# package from wherever it is installed (or loaded via devtools::load_all).
#
# ── Quick start ──────────────────────────────────────────────────────────────
#
#   # Install leachatetools once (if not already installed):
#   # remotes::install_github("vorpalvorpal/leachatetools")
#
#   source("run_amspaf_pipeline.R")
#
#   run_amspaf_pipeline(
#     db_path = "monitoring.duckdb"
#   )
#
# ── Full signature ────────────────────────────────────────────────────────────
#
#   run_amspaf_pipeline(
#     db_path,                      # path to DuckDB file (required)
#     leachatetools_dir  = NULL,    # path for devtools::load_all(); NULL = use installed pkg
#     focal_features     = NULL,    # character vector of feature names to analyse
#                                   #   NULL = all non-reference surface-water features
#     reference_features = NULL,    # character vector of reference feature names
#                                   #   NULL = features with reference=TRUE in DB
#     date_range         = NULL,    # Date[2]: c(start, end).  NULL = full history
#     focal_dates        = NULL,    # Date vector for chronic anchors.
#                                   #   NULL = first day of each month in date_range
#     tau_days           = 90,      # exponential-decay half-life (days)
#     window_days        = 365,     # look-back window for chronic integration (days)
#     percentile         = 0.80,    # ARA reference quantile
#     impute             = TRUE,    # Bayesian imputation (brms)?  FALSE = raw values only
#     impute_targets     = NULL,    # analytes to impute.  NULL = default set (see body)
#     impute_drivers     = c("pH", "EC"),
#     impute_surrogates  = list(DOC = c("TOC", "BOD", "COD")),
#     impute_iter        = 2000,
#     impute_warmup      = 1000,
#     impute_chains      = 4,
#     impute_cores       = parallel::detectCores(),
#     min_detect_freq    = 0.05,    # prescreen threshold
#     fill_temperature   = TRUE,    # interpolate missing Temperature for existing samples?
#     write_back         = TRUE,    # write AmsPAF rows back to DB?
#     replace_existing   = FALSE    # if TRUE, delete existing computed AmsPAF before writing
#   )
#
# ── Return value ──────────────────────────────────────────────────────────────
# A named list (invisibly):
#   $amspaf    tibble  — chronic AmsPAF per focal_date × feature
#   $imputed   tibble  — imputed chemistry (NULL if impute=FALSE)
#   $chronic   tibble  — chronic-integrated chemistry
#   $n_written integer — number of analysis rows written to DB (0 if write_back=FALSE)
#
# ── Notes ─────────────────────────────────────────────────────────────────────
# • All concentrations in the DB are in mg/L.  The script multiplies by 1000
#   to obtain µg/L for AmsPAF (which uses µg/L SSDs), with the exception of
#   analytes in NO_CONVERT (pH, Temperature, hardness, EC, etc.).
# • Temperature is renamed from "Temperature" → "temperature" to match the
#   NH3-N co-analyte name used in the analyte metadata.
# • NO3-N is classified into NO3-N_soft / NO3-N_mod / NO3-N_hard using the
#   co-sampled Hardness-total-CaCO3 value AFTER imputation.
# • Reference features are not imputed (raw field measurements only).
# • AmsPAF rows written back use a synthetic sample row with purpose =
#   "chronic_amspaf_computed".  Existing rows for the same feature × focal_date
#   are left untouched unless replace_existing = TRUE.
#
# ── DB schema assumed ─────────────────────────────────────────────────────────
#   analysis  (uuid, uuid_sample, uuid_lab, value, quantified, rl_low, rl_high,
#              purpose, comments)
#   sample    (uuid, uuid_feature, uuid_project, date, datetime, organisation,
#              purpose, comments)
#   lab_method(uuid, uuid_analyte, name, method, organisation, api,
#              uuid_project, uuid_feature, comments)
#   analyte   (uuid, name, units, ...)
#   feature   (uuid, name, reference, matrix, ...)
#
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(duckdb)
  library(DBI)
  library(lubridate)
})

# ── Analytes that must NOT be converted from mg/L to µg/L ────────────────────
.NO_CONVERT <- c(
  "pH", "Temperature", "Hardness-total-CaCO3",
  "ORP", "DO", "EC", "SAR", "Stage",
  "Water Height-10min", "Precipitation-10min", "Precipitation-1hr",
  "Precipitation-24hr", "Temperature-10min",
  "Wind Speed-10min", "Wind Direction-10min",
  # lower-case variants (already-renamed)
  "temperature"
)

# Analytes stored in the DB as derived rows — must never enter the pipeline
.DERIVED_ANALYTES <- c("AmsPAF", "LMF")

# Fixed UUIDs for computed analytes / methods (DB-specific, but stable)
.UUID_AMSPAF_ANALYTE  <- "b3b1f259-a05d-4dde-b153-3b91f03b51c5"
.UUID_AMSPAF_LABMETHOD <- "00000000-0000-0000-0000-000000000004"
.UUID_TEMP_LABMETHOD   <- "8ac3819d-92e4-4d66-a6f2-213d998b9354"  # "Internal" field temp

# ─────────────────────────────────────────────────────────────────────────────
# Main function
# ─────────────────────────────────────────────────────────────────────────────

run_amspaf_pipeline <- function(
    db_path,
    leachatetools_dir  = NULL,
    focal_features     = NULL,
    reference_features = NULL,
    date_range         = NULL,
    focal_dates        = NULL,
    tau_days           = 90,
    window_days        = 365,
    percentile         = 0.80,
    impute             = TRUE,
    impute_targets     = NULL,
    impute_drivers     = c("pH", "EC"),
    impute_surrogates  = list(DOC = c("TOC", "BOD", "COD")),
    impute_iter        = 2000,
    impute_warmup      = 1000,
    impute_chains      = 4,
    impute_cores       = parallel::detectCores(),
    min_detect_freq    = 0.05,
    fill_temperature   = TRUE,
    write_back         = TRUE,
    replace_existing   = FALSE
) {

  # ── 0. Load leachatetools ──────────────────────────────────────────────────
  if (!is.null(leachatetools_dir)) {
    message("Loading leachatetools from ", leachatetools_dir)
    devtools::load_all(leachatetools_dir, quiet = TRUE)
  } else {
    if (!requireNamespace("leachatetools", quietly = TRUE))
      stop("leachatetools package not found. ",
           "Install with: remotes::install_github('vorpalvorpal/leachatetools') ",
           "or pass leachatetools_dir = '<path>' to load_all().")
    library(leachatetools)
  }

  stopifnot(file.exists(db_path))

  # ── 1. Connect and resolve features ────────────────────────────────────────
  message("\n=== 1. Connecting to DB and resolving features ===")
  con <- dbConnect(duckdb(), db_path, read_only = FALSE)
  on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

  all_features <- dbGetQuery(con, "SELECT uuid, name, reference, matrix FROM feature")

  # Focal features: supplied, or all non-reference freshwater surface features
  if (is.null(focal_features)) {
    focal_features <- all_features |>
      filter(!reference, matrix == "freshwater") |>
      pull(name)
    message("  Auto-detected focal features: ", paste(focal_features, collapse = ", "))
  }

  # Reference features: supplied, or all features with reference == TRUE
  if (is.null(reference_features)) {
    reference_features <- all_features |>
      filter(reference) |>
      pull(name)
    message("  Auto-detected reference features: ", paste(reference_features, collapse = ", "))
  }

  if (length(focal_features) == 0)
    stop("No focal features found. Supply focal_features or ensure the DB has ",
         "non-reference freshwater features.")
  if (length(reference_features) == 0)
    stop("No reference features found. Supply reference_features or ensure the DB has ",
         "features with reference=TRUE.")

  all_names <- c(focal_features, reference_features)
  feature_map <- all_features |>
    filter(name %in% all_names) |>
    select(uuid, name)

  missing_f <- setdiff(all_names, feature_map$name)
  if (length(missing_f) > 0)
    stop("Features not found in DB: ", paste(missing_f, collapse = ", "))

  uuid_in <- paste(sprintf("'%s'", feature_map$uuid), collapse = ", ")

  # ── 2. Fetch raw chemistry ─────────────────────────────────────────────────
  message("\n=== 2. Fetching chemistry ===")

  date_filter <- if (!is.null(date_range)) {
    if (length(date_range) != 2 || !inherits(date_range, "Date"))
      stop("date_range must be a Date[2] vector: c(start_date, end_date)")
    sprintf(" AND s.date BETWEEN '%s' AND '%s'", date_range[1], date_range[2])
  } else ""

  raw <- dbGetQuery(con, sprintf("
    SELECT
      an.uuid_sample  AS sample_id,
      s.uuid_feature  AS site_uuid,
      s.datetime      AS datetime,
      a.name          AS analyte,
      an.value        AS value,
      an.quantified   AS detected,
      an.rl_low       AS rl_low
    FROM analysis an
    JOIN sample    s  ON an.uuid_sample  = s.uuid
    JOIN lab_method lm ON an.uuid_lab    = lm.uuid
    JOIN analyte   a  ON lm.uuid_analyte = a.uuid
    WHERE s.uuid_feature IN (%s)%s
  ", uuid_in, date_filter))

  message("  Raw rows fetched: ", nrow(raw))

  # Join feature names
  raw <- raw |>
    left_join(feature_map, by = c("site_uuid" = "uuid")) |>
    rename(site_id = name) |>
    filter(!analyte %in% .DERIVED_ANALYTES)

  # ── 3. Unit conversion & rename ────────────────────────────────────────────
  message("\n=== 3. Unit conversion ===")

  chem <- raw |>
    mutate(
      value    = if_else(!analyte %in% .NO_CONVERT, value * 1000, value),
      datetime = as.Date(datetime),
      detected = as.logical(detected),
      analyte  = if_else(analyte == "Temperature", "temperature", analyte)
    )

  message("  Rows after conversion: ", nrow(chem))

  # ── 4. Fill missing temperature (optional) ─────────────────────────────────
  if (fill_temperature) {
    message("\n=== 4. Temperature gap-filling ===")
    n_written_temp <- .fill_temperature_gaps(chem, con, feature_map, write_back)
    if (n_written_temp > 0) {
      # Re-fetch to include the newly inserted temperature rows
      new_temp_raw <- dbGetQuery(con, sprintf("
        SELECT
          an.uuid_sample AS sample_id,
          s.uuid_feature AS site_uuid,
          s.datetime     AS datetime,
          a.name         AS analyte,
          an.value       AS value,
          an.quantified  AS detected,
          an.rl_low      AS rl_low
        FROM analysis an
        JOIN sample     s  ON an.uuid_sample  = s.uuid
        JOIN lab_method lm ON an.uuid_lab     = lm.uuid
        JOIN analyte    a  ON lm.uuid_analyte = a.uuid
        WHERE lm.uuid = '%s'
          AND s.uuid_feature IN (%s)%s
      ", .UUID_TEMP_LABMETHOD, uuid_in, date_filter))
      new_temp <- new_temp_raw |>
        left_join(feature_map, by = c("site_uuid" = "uuid")) |>
        rename(site_id = name) |>
        mutate(
          datetime = as.Date(datetime),
          detected = as.logical(detected),
          analyte  = "temperature"
        )
      # Add to chem (replacing any existing temperature rows for those samples)
      chem <- chem |>
        filter(!(sample_id %in% new_temp$sample_id & analyte == "temperature")) |>
        bind_rows(new_temp)
      message("  Temperature rows after gap-fill: ",
              sum(chem$analyte == "temperature"), " total")
    }
  }

  # ── 5. NO3-N hardness classification ──────────────────────────────────────
  message("\n=== 5. NO3-N hardness classification ===")
  chem <- .classify_no3(chem)

  # ── 6. Split focal vs reference ───────────────────────────────────────────
  focal_chem <- filter(chem, site_id %in% focal_features)
  ref_chem   <- filter(chem, site_id %in% reference_features)

  if (nrow(focal_chem) == 0) stop("No chemistry rows for focal features after filtering.")
  if (nrow(ref_chem)   == 0) stop("No chemistry rows for reference features after filtering.")

  message("  Focal rows: ", nrow(focal_chem), " | Reference rows: ", nrow(ref_chem))

  # ── 7. Prescreen ──────────────────────────────────────────────────────────
  message("\n=== 6. Prescreen analytes (k = ", min_detect_freq, ") ===")
  included <- prescreen_analytes(focal_chem, k = min_detect_freq)
  message("  Analytes passing prescreen: ", length(included), ": ",
          paste(sort(included), collapse = ", "))

  # Drivers and co-analytes must be retained even if below detection threshold
  keep_always <- c(impute_drivers, unlist(impute_surrogates),
                   "temperature", "Hardness-total-CaCO3",
                   "NO3-N_soft", "NO3-N_mod", "NO3-N_hard")
  focal_screened <- filter(focal_chem, analyte %in% union(included, keep_always))

  # ── 8. Imputation (optional) ───────────────────────────────────────────────
  imp_result <- NULL
  focal_final <- NULL

  if (impute) {
    message("\n=== 7. Bayesian imputation (impute_chemistry) ===")

    # Default imputation targets: SSD-eligible metals + NH3-N + temperature
    # Exclude NO3-N variants (classified post-imputation) and hardness (driver)
    if (is.null(impute_targets)) {
      impute_targets <- intersect(
        included,
        c("Al", "As", "B", "Cd", "Cr", "Cu", "Hg", "Mn", "Ni", "Pb", "Se", "Zn",
          "NH3-N", "NH3.N", "temperature")
      )
    }
    message("  Imputing: ", paste(sort(impute_targets), collapse = ", "))

    imp_input <- filter(focal_screened,
                        analyte %in% c(impute_targets, impute_drivers,
                                       unlist(impute_surrogates)))

    t0 <- proc.time()
    imp_result <- impute_chemistry(
      imp_input,
      drivers          = impute_drivers,
      driver_surrogates = impute_surrogates,
      iter             = impute_iter,
      warmup           = impute_warmup,
      chains           = impute_chains,
      cores            = impute_cores
    )
    t_imp <- round((proc.time() - t0)[3], 1)
    n_imp <- sum(imp_result$imputed, na.rm = TRUE)
    message("  Imputation done in ", t_imp, " s | Imputed rows: ", n_imp,
            sprintf(" (%.0f%%)", 100 * n_imp / nrow(imp_result)))

    # Identify which samples made it through imputation (had all required drivers)
    imputed_sample_ids <- unique(imp_result$sample_id)

    # Classify NO3-N for those samples using measured (pre-imputation) hardness
    # Use original focal_screened so hardness values are not lost
    no3_for_imputed <- .classify_no3(
      filter(focal_screened,
             analyte %in% c("NO3-N", "Hardness-total-CaCO3"),
             sample_id %in% imputed_sample_ids)
    ) |>
      filter(grepl("^NO3-N_", analyte))

    # Combine: imputed metals/NH3-N/temperature + measured NO3-N variants
    focal_final <- bind_rows(
      filter(imp_result, !grepl("^NO3-N", analyte)),
      no3_for_imputed
    )

  } else {
    message("\n=== 7. Imputation skipped (impute = FALSE) ===")
    focal_final <- focal_screened
  }

  # ── 9. Reference chemistry ─────────────────────────────────────────────────
  # Reference is always raw measured data — no imputation
  ref_final <- ref_chem |>
    filter(analyte %in% union(included, keep_always))

  # ── 10. Focal dates ────────────────────────────────────────────────────────
  message("\n=== 8. Focal dates ===")
  if (is.null(focal_dates)) {
    date_min <- min(focal_final$datetime, na.rm = TRUE)
    date_max <- max(focal_final$datetime, na.rm = TRUE)
    # First day of each month from first sample month to current month
    focal_dates <- seq(
      as.Date(format(date_min, "%Y-%m-01")),
      as.Date(format(date_max, "%Y-%m-01")),
      by = "month"
    )
  }
  message("  Focal dates: ", length(focal_dates), "  (",
          format(min(focal_dates), "%Y-%m-%d"), " – ",
          format(max(focal_dates), "%Y-%m-%d"), ")")

  # ── 11. Chronic chemistry ──────────────────────────────────────────────────
  message("\n=== 9. Chronic chemistry (tau=", tau_days, "d, window=", window_days, "d) ===")

  t0 <- proc.time()
  chr_focal <- suppressWarnings(
    compute_chronic_chemistry(focal_final, focal_dates,
                              tau_days = tau_days, window_days = window_days))
  message("  Focal chronic: ", round((proc.time() - t0)[3], 1), " s")

  t0 <- proc.time()
  chr_ref <- suppressWarnings(
    compute_chronic_chemistry(ref_final, focal_dates,
                              tau_days = tau_days, window_days = window_days))
  message("  Reference chronic: ", round((proc.time() - t0)[3], 1), " s")

  # Keep only focal dates that have at least one in-window sample
  covered <- chr_focal |>
    group_by(focal_date) |>
    summarise(n = max(n_samples_in_window), .groups = "drop") |>
    filter(n >= 1)
  chr_focal <- semi_join(chr_focal, covered, by = "focal_date")
  message("  Focal dates with data: ", nrow(covered), " of ", length(focal_dates))

  # ── 12. Reference preparation ──────────────────────────────────────────────
  message("\n=== 10. prepare_reference (", percentile * 100, "th pct ARA) ===")
  suppressMessages(prep_ref <- prepare_reference(chr_ref, percentile = percentile))
  message("  Reference analytes: ", nrow(prep_ref$normalised_quantiles),
          " | Dropped: ", length(prep_ref$dropped))

  # ── 13. AmsPAF ─────────────────────────────────────────────────────────────
  message("\n=== 11. add_amspaf() ===")
  t0 <- proc.time()
  paf_out <- suppressMessages(add_amspaf(chr_focal, reference = prep_ref, min_analytes = 3))
  message("  add_amspaf: ", round((proc.time() - t0)[3], 1), " s")

  amspaf <- paf_out |>
    filter(analyte == "AmsPAF") |>
    select(focal_date, site_id,
           value, n_analytes_used,
           any_of(c("n_analytes_imputed", "dominant_analyte", "max_paf"))) |>
    arrange(site_id, focal_date)

  message("\n  Results:")
  message("    Focal-date × feature rows: ", nrow(amspaf))
  message("    AmsPAF range: ",
          round(min(amspaf$value, na.rm = TRUE), 3), " – ",
          round(max(amspaf$value, na.rm = TRUE), 2), " %")
  if ("dominant_analyte" %in% names(amspaf)) {
    dom <- paste(sort(unique(na.omit(amspaf$dominant_analyte))), collapse = ", ")
    message("    Dominant analytes: ", dom)
  }

  # ── 14. Write back to DB ───────────────────────────────────────────────────
  n_written <- 0L
  if (write_back) {
    message("\n=== 12. Writing AmsPAF rows to DB ===")
    n_written <- .write_amspaf_to_db(amspaf, con, feature_map, replace_existing)
    message("  Rows written: ", n_written)
  }

  message("\n=== Done ===")
  invisible(list(
    amspaf   = amspaf,
    imputed  = imp_result,
    chronic  = chr_focal,
    n_written = n_written
  ))
}


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

#' Classify NO3-N into hardness variants
#' @param df Long-format chemistry df containing "NO3-N" and optionally
#'   "Hardness-total-CaCO3" rows
#' @return df with NO3-N rows replaced by NO3-N_soft/mod/hard, hardness removed
.classify_no3 <- function(df) {
  hardness <- df |>
    filter(analyte == "Hardness-total-CaCO3", detected) |>
    select(sample_id, hardness_mgL = value)

  no3 <- df |>
    filter(analyte == "NO3-N") |>
    left_join(hardness, by = "sample_id") |>
    mutate(analyte = case_when(
      !is.na(hardness_mgL) & hardness_mgL < 30    ~ "NO3-N_soft",
      !is.na(hardness_mgL) & hardness_mgL <= 150   ~ "NO3-N_mod",
      !is.na(hardness_mgL) & hardness_mgL > 150    ~ "NO3-N_hard",
      TRUE                                          ~ "NO3-N_mod"
    )) |>
    select(-hardness_mgL)

  n_by_type <- count(no3, analyte)
  if (nrow(n_by_type) > 0) {
    msg <- paste(sprintf("%s: %d", n_by_type$analyte, n_by_type$n), collapse = "  ")
    message("  NO3-N classified → ", msg)
  }

  df |>
    filter(!analyte %in% c("NO3-N", "Hardness-total-CaCO3")) |>
    bind_rows(no3)
}


#' Interpolate and write missing Temperature rows for existing samples
#'
#' For each sample that has no Temperature analysis row, estimates temperature
#' by linear interpolation from the nearest Temperature measurements on the
#' same feature.  Falls back to 20 °C if no measurements exist within a
#' configurable window.
#'
#' @param chem  Chemistry df (already unit-converted; "temperature" rows present)
#' @param con   Open DuckDB connection (read-write)
#' @param feature_map  Tibble with uuid/name for all relevant features
#' @param write_back  If FALSE, only report gaps without writing
#' @param fallback_temp Default temperature (°C) when interpolation fails
#' @param max_gap_days Interpolate only when nearest measurement is within this
#'   many days; beyond this, use fallback_temp
#' @return Number of analysis rows written (0 if write_back = FALSE)
.fill_temperature_gaps <- function(chem, con, feature_map,
                                   write_back     = TRUE,
                                   fallback_temp  = 20,
                                   max_gap_days   = 180) {

  # Samples that already have temperature
  samples_with_temp <- chem |>
    filter(analyte == "temperature") |>
    pull(sample_id) |>
    unique()

  # All (sample_id, datetime, site_id) combos in the focal chemistry
  all_samples <- chem |>
    distinct(sample_id, datetime, site_id) |>
    filter(!sample_id %in% samples_with_temp)

  if (nrow(all_samples) == 0) {
    message("  All samples already have temperature — no gap-filling needed.")
    return(0L)
  }
  message("  Samples lacking temperature: ", nrow(all_samples))

  # Existing temperature readings per feature
  temp_known <- chem |>
    filter(analyte == "temperature", detected, !is.na(value)) |>
    distinct(site_id, datetime, value) |>
    arrange(site_id, datetime)

  # Estimate temperature for each gap via linear interpolation
  new_rows <- all_samples |>
    group_by(site_id) |>
    group_modify(function(gaps, key) {
      feat_temps <- filter(temp_known, site_id == key$site_id) |>
        arrange(datetime)
      if (nrow(feat_temps) == 0) {
        message("    ", key$site_id,
                ": no Temperature measurements found — using fallback ", fallback_temp, " °C")
        return(mutate(gaps, temp_est = fallback_temp))
      }

      # approx() does linear interpolation with constant extrapolation
      # Deduplicate by date first (average same-day readings)
      feat_temps <- feat_temps |>
        group_by(datetime) |>
        summarise(value = mean(value), .groups = "drop") |>
        arrange(datetime)
      t_num <- as.numeric(feat_temps$datetime)
      v     <- feat_temps$value
      g_num <- as.numeric(gaps$datetime)

      # Use approx with rule=2 (extrapolate using boundary values)
      interp <- approx(t_num, v, xout = g_num, rule = 2)$y

      # Cap extrapolation: replace with fallback if nearest known measurement
      # is more than max_gap_days away
      nearest_gap <- vapply(g_num, function(gd) {
        min(abs(t_num - gd))
      }, numeric(1))
      interp[nearest_gap > max_gap_days] <- fallback_temp

      gaps$temp_est <- round(interp, 1)
      gaps
    }) |>
    ungroup()

  n_fallback <- sum(new_rows$temp_est == fallback_temp)
  if (n_fallback > 0)
    message("  Using fallback temperature (", fallback_temp,
            " °C) for ", n_fallback, " samples (nearest measurement > ",
            max_gap_days, " d)")

  message("  Interpolated temperature for ", nrow(new_rows), " samples")

  if (!write_back) return(0L)

  # Build analysis rows to insert
  # Generate UUIDs for new rows
  if (!requireNamespace("uuid", quietly = TRUE))
    stop("Package 'uuid' is required for writing temperature rows. ",
         "Install with: install.packages('uuid')")

  insert_rows <- new_rows |>
    mutate(
      uuid      = vapply(seq_len(n()), function(i) uuid::UUIDgenerate(), character(1)),
      uuid_lab  = .UUID_TEMP_LABMETHOD,
      value     = temp_est,
      quantified = TRUE,
      rl_low    = NA_real_,
      rl_high   = NA_real_,
      purpose   = "temperature_interpolated",
      comments  = sprintf("Interpolated by run_amspaf_pipeline.R on %s", Sys.Date())
    ) |>
    select(uuid, uuid_sample = sample_id, uuid_lab, value, quantified,
           rl_low, rl_high, purpose, comments)

  # Check which columns actually exist in the analysis table
  analysis_cols <- dbListFields(con, "analysis")
  insert_rows <- select(insert_rows, any_of(analysis_cols))

  dbAppendTable(con, "analysis", insert_rows)
  message("  Wrote ", nrow(insert_rows), " temperature rows to analysis table")
  nrow(insert_rows)
}


#' Write chronic AmsPAF results back to the DB
#'
#' Creates a synthetic sample row for each (focal_date × feature) and an
#' analysis row containing the AmsPAF value.  Existing computed rows for the
#' same feature × focal_date are deleted first if replace_existing = TRUE.
#'
#' @param amspaf  Tibble with columns focal_date, site_id, value, ...
#' @param con     Open DuckDB connection (read-write)
#' @param feature_map  Tibble with uuid/name
#' @param replace_existing  Delete existing computed AmsPAF before writing?
#' @return Number of analysis rows written
.write_amspaf_to_db <- function(amspaf, con, feature_map, replace_existing) {

  if (!requireNamespace("uuid", quietly = TRUE))
    stop("Package 'uuid' is required for writing AmsPAF rows. ",
         "Install with: install.packages('uuid')")

  # Join feature UUIDs
  amspaf_w_uuid <- amspaf |>
    left_join(feature_map, by = c("site_id" = "name")) |>
    rename(feature_uuid = uuid)

  if (any(is.na(amspaf_w_uuid$feature_uuid)))
    warning("Some site_ids could not be matched to feature UUIDs — those rows skipped.")
  amspaf_w_uuid <- filter(amspaf_w_uuid, !is.na(feature_uuid))
  if (nrow(amspaf_w_uuid) == 0) return(0L)

  # Check for existing computed AmsPAF rows (by joining sample → analysis → lab_method)
  existing <- dbGetQuery(con, sprintf("
    SELECT s.uuid AS uuid_sample, s.uuid_feature, s.date AS focal_date
    FROM sample s
    JOIN analysis an ON an.uuid_sample = s.uuid
    WHERE s.purpose = 'chronic_amspaf_computed'
      AND an.uuid_lab = '%s'
      AND s.uuid_feature IN (%s)
  ", .UUID_AMSPAF_LABMETHOD,
     paste(sprintf("'%s'", unique(amspaf_w_uuid$feature_uuid)), collapse = ", ")))

  if (nrow(existing) > 0) {
    if (replace_existing) {
      message("  Deleting ", nrow(existing), " existing computed AmsPAF rows …")
      sample_uuids <- paste(sprintf("'%s'", existing$uuid_sample), collapse = ", ")
      dbExecute(con, sprintf("DELETE FROM analysis WHERE uuid_sample IN (%s)", sample_uuids))
      dbExecute(con, sprintf("DELETE FROM sample   WHERE uuid         IN (%s)", sample_uuids))
    } else {
      # Skip rows already in DB
      existing_keys <- existing |>
        mutate(focal_date = as.Date(focal_date)) |>
        select(feature_uuid = uuid_feature, focal_date)
      amspaf_w_uuid <- amspaf_w_uuid |>
        anti_join(existing_keys, by = c("feature_uuid", "focal_date"))
      if (nrow(amspaf_w_uuid) == 0) {
        message("  All AmsPAF rows already exist (replace_existing = FALSE) — nothing written.")
        return(0L)
      }
      message("  Skipping ", nrow(existing), " existing rows; writing ",
              nrow(amspaf_w_uuid), " new rows …")
    }
  }

  # Build sample rows (one per focal_date × feature)
  sample_rows <- amspaf_w_uuid |>
    mutate(
      uuid         = vapply(seq_len(n()), function(i) uuid::UUIDgenerate(), character(1)),
      date         = focal_date,
      datetime     = as.POSIXct(paste(focal_date, "00:00:00")),
      organisation = "leachatetools",
      purpose      = "chronic_amspaf_computed",
      comments     = sprintf("Chronic AmsPAF computed by run_amspaf_pipeline.R on %s", Sys.Date())
    ) |>
    select(uuid, uuid_feature = feature_uuid, date, datetime,
           organisation, purpose, comments)

  # Check which columns exist in sample table
  sample_cols <- dbListFields(con, "sample")
  sample_rows_insert <- select(sample_rows, any_of(sample_cols))
  dbAppendTable(con, "sample", sample_rows_insert)

  # Build analysis rows (one per focal_date × feature)
  analysis_rows <- amspaf_w_uuid |>
    mutate(
      uuid_sample = sample_rows$uuid,
      uuid        = vapply(seq_len(n()), function(i) uuid::UUIDgenerate(), character(1)),
      uuid_lab    = .UUID_AMSPAF_LABMETHOD,
      quantified  = TRUE,
      purpose     = "chronic_amspaf_computed",
      comments    = sprintf(
        "tau=%dd window=%dd pct=%g%% analytes_used=%s%s",
        # tau/window/pct are not in the amspaf tibble so we record them in comments
        NA_integer_, NA_integer_, NA_real_,
        if ("n_analytes_used" %in% names(amspaf_w_uuid)) as.character(amspaf_w_uuid$n_analytes_used) else "?",
        if ("dominant_analyte" %in% names(amspaf_w_uuid))
          paste0(" dominant=", amspaf_w_uuid$dominant_analyte) else ""
      )
    ) |>
    select(uuid, uuid_sample, uuid_lab, value,
           quantified, purpose, comments)

  analysis_cols <- dbListFields(con, "analysis")
  analysis_rows_insert <- select(analysis_rows, any_of(analysis_cols))
  dbAppendTable(con, "analysis", analysis_rows_insert)

  nrow(analysis_rows_insert)
}


# ─────────────────────────────────────────────────────────────────────────────
# Example / quick-run block (edit and uncomment to run directly)
# ─────────────────────────────────────────────────────────────────────────────

# source("run_amspaf_pipeline.R")
#
# result <- run_amspaf_pipeline(
#   db_path           = "monitoring.duckdb",
#   # leachatetools_dir = "/path/to/leachatetools",  # uncomment if not installed
#   focal_features    = c("B.S01"),                  # NULL = all non-reference features
#   reference_features = c("B.S03"),                 # NULL = auto-detect from DB
#   date_range        = c(as.Date("2017-01-01"), Sys.Date()),
#   tau_days          = 90,
#   window_days       = 365,
#   percentile        = 0.80,
#   impute            = TRUE,
#   impute_iter       = 2000,
#   impute_chains     = 4,
#   write_back        = TRUE,
#   replace_existing  = FALSE
# )
#
# # Inspect results
# print(result$amspaf)
# ggplot2::ggplot(result$amspaf, ggplot2::aes(focal_date, value)) +
#   ggplot2::geom_line() +
#   ggplot2::facet_wrap(~site_id)
