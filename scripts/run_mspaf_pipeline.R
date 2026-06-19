# run_mspaf_pipeline.R ──────────────────────────────────────────────────────
#
# Standalone helper script: reads chemistry from a monitoring DuckDB (same
# schema as hydroSense/test data/monitoring.duckdb), runs the full
# chronic msPAF pipeline, and writes results back to the DB.
#
# This script is NOT part of the hydroSense package — it is intended to be
# dropped next to a DB copy and run as needed.  It calls the hydroSense
# package from wherever it is installed (or loaded via devtools::load_all).
#
# ── Quick start ──────────────────────────────────────────────────────────────
#
#   # Install hydroSense once (if not already installed):
#   # remotes::install_github("vorpalvorpal/hydroSense")
#
#   source("run_mspaf_pipeline.R")
#
#   run_mspaf_pipeline(
#     db_path = "monitoring.duckdb"
#   )
#
# ── Full signature ────────────────────────────────────────────────────────────
#
#   run_mspaf_pipeline(
#     db_path,                         # path to DuckDB file (required)
#     hydroSense_dir  = NULL,       # path for devtools::load_all(); NULL = installed pkg
#     focal_features     = NULL,       # character vector of feature names to analyse
#                                      #   NULL = all non-reference surface-water features
#     reference_features = NULL,       # character vector of reference feature names
#                                      #   NULL = features with reference=TRUE in DB
#     date_range         = NULL,       # Date[2]: c(start, end).  NULL = full history
#     focal_dates        = NULL,       # Date vector for chronic anchors.
#                                      #   NULL = first day of each month in date_range
#     tau_days           = 90,         # exponential-decay half-life (days)
#     window_days        = 365,        # look-back window for chronic integration (days)
#     ref_summary        = "geom_mean", # ARA reference summary statistic
#                                       # ("geom_mean", "median", "p80", ...)
#     impute             = TRUE,       # Bayesian imputation (brms)?  FALSE = raw values only
#     required_vars      = c("pH", "EC"),
#     impute_iter        = 2000,
#     impute_warmup      = 1000,
#     impute_chains      = 4,
#     impute_cores       = parallel::detectCores(),
#     model_dir          = NULL,       # directory for .qs model files
#                                      #   NULL = <db_dir>/models/
#     refit_model        = FALSE,      # if TRUE, always refit even if a cached model exists
#     min_detect_freq    = 0.05,       # prescreen threshold
#     fill_temperature   = TRUE,       # interpolate missing Temperature?
#     write_back         = TRUE,       # write msPAF rows back to DB?
#     replace_existing   = FALSE       # if TRUE, delete existing computed msPAF before writing
#   )
#
# ── Return value ──────────────────────────────────────────────────────────────
# A named list (invisibly):
#   $mspaf      tibble  — chronic msPAF per focal_date × feature
#   $persample   tibble  — per-sample msPAF + diagnostics
#                          (n_analytes_used, dominant_analyte, max_paf, ...)
#   $imputed     tibble  — imputed chemistry (NULL if impute=FALSE)
#   $n_written   integer — number of analysis rows written to DB (0 if write_back=FALSE)
#
# ── Notes ─────────────────────────────────────────────────────────────────────
# • All concentrations in the DB are in mg/L.  The script multiplies by 1000
#   to obtain µg/L for msPAF (which uses µg/L SSDs), with the exception of
#   analytes in NO_CONVERT (pH, Temperature, hardness, EC, etc.).
# • Temperature is renamed from "Temperature" → "temperature" to match the
#   NH3-N co-analyte name used in the analyte metadata.
# • Duplicate (sample_id, analyte) rows (e.g. field + lab EC/pH) are resolved
#   by preferring the row where lab_method.method = "field".
# • NO3-N is classified into NO3-N_soft / NO3-N_mod / NO3-N_hard using the
#   co-sampled hardness value AFTER imputation (so that hardness
#   is available to the WQ PCA during model fitting).
# • Reference features are not imputed (raw field measurements only).
# • Fitted imputation models are saved as .qs files and their paths recorded
#   in the imputation_models table in the DB.  On subsequent runs the cached
#   model is reused unless refit_model = TRUE or the file is missing.
# • msPAF rows written back use a synthetic sample row with purpose =
#   "chronic_mspaf_computed".  Existing rows for the same feature × focal_date
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
.DERIVED_ANALYTES <- c("msPAF", "LMF")

# Fixed UUIDs for computed analytes / methods (DB-specific, but stable)
.UUID_MSPAF_ANALYTE   <- "b3b1f259-a05d-4dde-b153-3b91f03b51c5"
.UUID_MSPAF_LABMETHOD <- "00000000-0000-0000-0000-000000000004"
.UUID_TEMP_LABMETHOD   <- "8ac3819d-92e4-4d66-a6f2-213d998b9354"  # "Internal" field temp

# ─────────────────────────────────────────────────────────────────────────────
# Main function
# ─────────────────────────────────────────────────────────────────────────────

run_mspaf_pipeline <- function(
    db_path,
    hydroSense_dir  = NULL,
    focal_features     = NULL,
    reference_features = NULL,
    date_range         = NULL,
    focal_dates        = NULL,
    tau_days           = 90,
    window_days        = 365,
    ref_summary        = "geom_mean",
    impute             = TRUE,
    required_vars      = c("pH", "EC"),
    impute_iter        = 2000,
    impute_warmup      = 1000,
    impute_chains      = 4,
    impute_cores       = parallel::detectCores(),
    model_dir          = NULL,
    refit_model        = FALSE,
    min_detect_freq    = 0.05,
    fill_temperature   = TRUE,
    write_back         = TRUE,
    replace_existing   = FALSE
) {

  # ── 0. Load hydroSense ──────────────────────────────────────────────────
  if (!is.null(hydroSense_dir)) {
    message("Loading hydroSense from ", hydroSense_dir)
    devtools::load_all(hydroSense_dir, quiet = TRUE)
  } else {
    if (!requireNamespace("hydroSense", quietly = TRUE))
      stop("hydroSense package not found. ",
           "Install with: remotes::install_github('vorpalvorpal/hydroSense') ",
           "or pass hydroSense_dir = '<path>' to load_all().")
    library(hydroSense)
  }

  stopifnot(file.exists(db_path))

  # Resolve model_dir now so it is available before the DB connection opens
  if (is.null(model_dir)) {
    model_dir <- file.path(dirname(normalizePath(db_path)), "models")
  }

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
      an.rl_low       AS rl_low,
      lm.method       AS method
    FROM analysis an
    JOIN sample     s  ON an.uuid_sample  = s.uuid
    JOIN lab_method lm ON an.uuid_lab     = lm.uuid
    JOIN analyte    a  ON lm.uuid_analyte = a.uuid
    WHERE s.uuid_feature IN (%s)%s
  ", uuid_in, date_filter))

  message("  Raw rows fetched: ", nrow(raw))

  # Deduplicate: for (sample_id, analyte) pairs with more than one row,
  # prefer the row where method == "field" (case-insensitive, e.g. field EC/pH);
  # fall back to whichever row comes first.
  n_before_dedup <- nrow(raw)
  raw <- raw |>
    mutate(.field_first = tolower(method) == "field") |>
    arrange(sample_id, analyte, desc(.field_first)) |>
    group_by(sample_id, analyte) |>
    slice(1L) |>
    ungroup() |>
    select(-.field_first)
  n_dedup <- n_before_dedup - nrow(raw)
  if (n_dedup > 0L)
    message("  Deduplication: removed ", n_dedup,
            " duplicate row(s) — preferred 'field' method where available")

  # Join feature names; filter out derived analytes
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
      # Translate DB analyte names to the names the package uses internally.
      analyte  = dplyr::case_when(
        analyte == "Temperature"           ~ "temperature",
        analyte == "Hardness-total-CaCO3"  ~ "hardness",
        TRUE                                ~ analyte
      )
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
          an.rl_low      AS rl_low,
          lm.method      AS method
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

  # ── 5. Split focal vs reference ───────────────────────────────────────────
  # Note: NO3-N classification is deferred to step 8 (after imputation) so
  # that hardness remains available to the WQ PCA during fitting.
  focal_chem <- filter(chem, site_id %in% focal_features)
  ref_chem   <- filter(chem, site_id %in% reference_features)

  if (nrow(focal_chem) == 0) stop("No chemistry rows for focal features after filtering.")
  if (nrow(ref_chem)   == 0) stop("No chemistry rows for reference features after filtering.")

  message("  Focal rows: ", nrow(focal_chem), " | Reference rows: ", nrow(ref_chem))

  # ── 5b. Pre-imputation hardness reconciliation ────────────────────────────
  # Fill hardness from Ca+Mg (or Ca/Mg from hardness + the other) wherever
  # exactly two of the three are measured.  Warn if all three present but
  # inconsistent.  Cheap, no model required — derived from stoichiometry.
  message("\n=== 5b. Hardness reconciliation (pre-imputation) ===")
  focal_chem <- derive_hardness(focal_chem, verbose = TRUE)
  ref_chem   <- derive_hardness(ref_chem,   verbose = TRUE)

  # ── 6. Prescreen ──────────────────────────────────────────────────────────
  message("\n=== 6. Prescreen analytes (k = ", min_detect_freq, ") ===")
  included <- prescreen_analytes(focal_chem, k = min_detect_freq,
                                  protect = required_vars)
  message("  Analytes passing prescreen: ", length(included), ": ",
          paste(sort(included), collapse = ", "))

  # Always retain required_vars, temperature, hardness (for PCA), and raw NO3-N
  # (for classification post-imputation).  Classified NO3-N variants (NO3-N_soft
  # etc.) do not exist yet at this stage — they are produced in step 8.
  keep_always <- c(required_vars, "pH", "EC", "NH3-N",
                   "temperature", "hardness", "NO3-N")
  focal_screened <- filter(focal_chem, analyte %in% union(included, keep_always))

  # ── 7. Imputation (optional) ───────────────────────────────────────────────
  imp_result <- NULL
  focal_imp  <- NULL

  if (impute) {
    message("\n=== 7. Bayesian imputation ===")
    t0 <- proc.time()

    model <- .load_or_fit_model(
      con            = con,
      focal_features = focal_features,
      df             = focal_screened,
      required_vars  = required_vars,
      model_dir      = model_dir,
      refit_model    = refit_model,
      iter           = impute_iter,
      warmup         = impute_warmup,
      chains         = impute_chains,
      cores          = impute_cores
    )

    imp_result <- impute_chemistry(focal_screened, model)
    t_imp <- round((proc.time() - t0)[3], 1)
    n_imp <- sum(imp_result$imputed, na.rm = TRUE)
    message("  Metals/organics imputation done in ", t_imp, " s | Rows: ", n_imp,
            sprintf(" (%.0f%%)", 100 * n_imp / nrow(imp_result)))

    # Co-analyte imputation (DOC, Ca, Mg, hardness) — separate GAM step,
    # uses same PCA but never feeds back into metals model.
    message("  Imputing co-analytes (DOC, Ca, Mg, hardness) …")
    focal_imp <- impute_coanalytes(imp_result, model)
    n_coa <- sum(focal_imp$imputed & focal_imp$imputed_kind == "missing" &
                   focal_imp$analyte %in% hydroSense:::.COANALYTE_TARGETS,
                 na.rm = TRUE)
    message("  Co-analyte rows imputed: ", n_coa)

    # Post-imputation hardness reconciliation — fill hardness wherever
    # Ca and Mg were just imputed but hardness itself wasn't.
    message("  Reconciling hardness from imputed Ca/Mg …")
    focal_imp <- derive_hardness(focal_imp, verbose = TRUE)

  } else {
    message("\n=== 7. Imputation skipped (impute = FALSE) ===")
    focal_imp <- focal_screened
  }

  # ── 8. NO3-N hardness classification (after imputation) ───────────────────
  # Hardness rows were preserved through imputation as WQ block variables.
  # Classify NO3-N now; hardness rows are kept for downstream normalisation.
  message("\n=== 8. NO3-N hardness classification ===")
  focal_final <- .classify_no3(focal_imp)

  # ── 9. Reference chemistry ─────────────────────────────────────────────────
  # Reference is always raw measured data — no imputation.
  ref_screened <- ref_chem |>
    filter(analyte %in% union(included, keep_always))
  ref_final <- .classify_no3(ref_screened)

  # ── 10. Focal dates ────────────────────────────────────────────────────────
  message("\n=== 10. Focal dates ===")
  if (is.null(focal_dates)) {
    date_min <- min(focal_final$datetime, na.rm = TRUE)
    date_max <- max(focal_final$datetime, na.rm = TRUE)
    # First day of each month from first sample month to last sample month
    focal_dates <- seq(
      as.Date(format(date_min, "%Y-%m-01")),
      as.Date(format(date_max, "%Y-%m-01")),
      by = "month"
    )
  }
  message("  Focal dates: ", length(focal_dates), "  (",
          format(min(focal_dates), "%Y-%m-%d"), " – ",
          format(max(focal_dates), "%Y-%m-%d"), ")")

  # ── 11. Reference preparation ──────────────────────────────────────────────
  # Build the reference summary ONCE from per-sample reference chemistry.
  # The summary is a single value per analyte representing what the local
  # community is integrated against — it does not vary by focal_date.
  message("\n=== 11. prepare_reference (summary = ", ref_summary, ") ===")
  suppressMessages(
    prep_ref <- prepare_reference(ref_final, summary = ref_summary)
  )
  message("  Reference analytes: ", nrow(prep_ref$ref_table),
          " | Dropped: ", length(prep_ref$dropped))

  # ── 12. Per-sample msPAF ──────────────────────────────────────────────────
  # Compute msPAF on each individual focal sample's imputed chemistry.
  # This produces a per-sample msPAF value that we'll time-aggregate next.
  message("\n=== 12. Per-sample msPAF ===")
  t0 <- proc.time()
  paf_persample <- suppressMessages(
    add_mspaf(focal_final, reference = prep_ref, min_analytes = 3)
  ) |>
    filter(analyte == "msPAF")
  message("  add_mspaf: ", round((proc.time() - t0)[3], 1), " s | ",
          "Per-sample msPAF rows: ", nrow(paf_persample))

  # ── 13. Chronic msPAF (Path B: time-aggregate per-sample msPAFs) ─────────
  # Time-weighted ARITHMETIC mean of per-sample msPAF values.
  # Arithmetic mean is appropriate here because msPAF is a bounded
  # percentage representing fraction of species affected; biology integrates
  # the toxic-response signal linearly over time.
  message("\n=== 13. Chronic msPAF (tau=", tau_days, "d, window=", window_days, "d) ===")
  t0 <- proc.time()
  chr_mspaf <- suppressWarnings(
    time_weighted_aggregate(
      paf_persample,
      focal_dates    = focal_dates,
      tau_days       = tau_days,
      window_days    = window_days,
      summary        = "arith_mean"
    )
  )
  message("  time_weighted_aggregate: ", round((proc.time() - t0)[3], 1), " s")

  mspaf <- chr_mspaf |>
    select(focal_date, site_id, value,
           n_samples_in_window,
           any_of("n_imputed_in_window")) |>
    arrange(site_id, focal_date)

  message("\n  Results:")
  message("    Focal-date × feature rows: ", nrow(mspaf))
  message("    msPAF range: ",
          round(min(mspaf$value, na.rm = TRUE), 3), " – ",
          round(max(mspaf$value, na.rm = TRUE), 2), " %")
  # Per-sample diagnostics (if user wants to inspect)
  dom <- sort(unique(na.omit(paf_persample$dominant_analyte)))
  if (length(dom) > 0L)
    message("    Per-sample dominant analytes: ", paste(dom, collapse = ", "))

  # ── 14. Write back to DB ───────────────────────────────────────────────────
  n_written <- 0L
  if (write_back) {
    message("\n=== 14. Writing msPAF rows to DB ===")
    n_written <- .write_mspaf_to_db(mspaf, con, feature_map, replace_existing)
    message("  Rows written: ", n_written)
  }

  message("\n=== Done ===")
  invisible(list(
    mspaf        = mspaf,           # chronic msPAF per (focal_date, site)
    persample     = paf_persample,    # per-sample msPAF + diagnostics
    imputed       = imp_result,
    n_written     = n_written
  ))
}


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

# Null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Load a cached imputation model from the DB, or fit a new one
#'
#' Creates the `imputation_models` table in the DB if it does not exist.
#' Looks for a cached model matching `focal_features`; loads and returns it if
#' the `.qs` file is still on disk.  Otherwise fits via `fit_imputation_model()`,
#' saves the result, and inserts the path into the DB.
#'
#' @param con Open DuckDB connection (read-write)
#' @param focal_features Character vector of focal feature names (used as key)
#' @param df Long-format chemistry df passed to `fit_imputation_model()`
#' @param required_vars Character vector of required variable names
#' @param model_dir Directory for .qs model files
#' @param refit_model Logical: if TRUE, bypass cache and always refit
#' @param ... Passed to `fit_imputation_model()` (iter, warmup, chains, cores)
#' @return An object of class `"imputation_model"`
.load_or_fit_model <- function(
    con, focal_features, df, required_vars, model_dir,
    refit_model, ...
) {
  # Ensure table exists
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS imputation_models (
      uuid           TEXT PRIMARY KEY,
      focal_features TEXT,
      model_path     TEXT,
      fit_date       DATE,
      fit_settings   TEXT,
      comments       TEXT
    )
  ")

  # Cache key: sorted, comma-joined feature names
  features_key <- paste(sort(focal_features), collapse = ",")

  # Try to use cached model (most-recently fitted for this feature set)
  if (!refit_model) {
    cached <- dbGetQuery(con, sprintf(
      "SELECT * FROM imputation_models WHERE focal_features = '%s' ORDER BY fit_date DESC LIMIT 1",
      gsub("'", "''", features_key)   # escape single quotes in feature names
    ))
    if (nrow(cached) > 0 && nzchar(cached$model_path[1]) &&
        file.exists(cached$model_path[1])) {
      message("  Loading cached imputation model from ", cached$model_path[1],
              "  (fitted ", cached$fit_date[1], ")")
      if (!requireNamespace("qs", quietly = TRUE))
        stop("Package 'qs' required to load cached model. ",
             "Install with: install.packages('qs')")
      return(qs::qread(cached$model_path[1]))
    } else if (nrow(cached) > 0) {
      message("  Cached model record found but file missing — refitting: ",
              cached$model_path[1])
    }
  }

  # Fit new model
  message("  Fitting new imputation model (this may take several minutes) …")
  extra <- list(...)
  model <- fit_imputation_model(
    df,
    required_vars = required_vars,
    save_dir      = model_dir,
    ...
  )

  # Insert record into DB
  if (!requireNamespace("uuid", quietly = TRUE))
    stop("Package 'uuid' is required for recording model path. ",
         "Install with: install.packages('uuid')")

  save_path    <- attr(model, "save_path") %||% ""
  fit_settings <- sprintf(
    "required_vars=%s;n_pca_vars=%d;iter=%d;chains=%d",
    paste(model$required_vars, collapse = ","),
    length(model$pca_vars),
    extra$iter   %||% 2000L,
    extra$chains %||% 4L
  )

  dbExecute(con, sprintf(
    paste0(
      "INSERT INTO imputation_models ",
      "(uuid, focal_features, model_path, fit_date, fit_settings, comments) ",
      "VALUES ('%s', '%s', '%s', '%s', '%s', '%s')"
    ),
    uuid::UUIDgenerate(),
    gsub("'", "''", features_key),
    gsub("'", "''", save_path),
    format(Sys.Date(), "%Y-%m-%d"),
    gsub("'", "''", fit_settings),
    paste0("Fitted by run_mspaf_pipeline.R on ", Sys.Date())
  ))
  if (nzchar(save_path))
    message("  Model path recorded in imputation_models table: ", save_path)
  else
    message("  Warning: model path is empty (save_dir may not have been writable).")

  model
}


#' Classify NO3-N into hardness variants
#'
#' @param df Long-format chemistry df containing "NO3-N" and optionally
#'   "hardness" rows.
#' @return df with NO3-N rows replaced by NO3-N_soft/mod/hard; raw NO3-N
#'   rows are removed from the result.  Hardness rows are preserved
#'   downstream so that add_mspaf() can use them for Cd/Pb/Zn normalisation.
.classify_no3 <- function(df) {
  hardness <- df |>
    filter(analyte == "hardness", detected) |>
    select(sample_id, hardness_mgL = value)

  no3 <- df |>
    filter(analyte == "NO3-N") |>
    left_join(hardness, by = "sample_id") |>
    mutate(analyte = case_when(
      !is.na(hardness_mgL) & hardness_mgL <  30  ~ "NO3-N_soft",
      !is.na(hardness_mgL) & hardness_mgL <= 150 ~ "NO3-N_mod",
      !is.na(hardness_mgL) & hardness_mgL >  150 ~ "NO3-N_hard",
      TRUE                                        ~ "NO3-N_mod"
    )) |>
    select(-hardness_mgL)

  n_by_type <- count(no3, analyte)
  if (nrow(n_by_type) > 0) {
    msg <- paste(sprintf("%s: %d", n_by_type$analyte, n_by_type$n), collapse = "  ")
    message("  NO3-N classified → ", msg)
  }

  # Remove raw NO3-N rows; keep hardness for downstream normalisation.
  df |>
    filter(analyte != "NO3-N") |>
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
                                   write_back    = TRUE,
                                   fallback_temp = 20,
                                   max_gap_days  = 180) {

  # Samples that already have temperature
  samples_with_temp <- chem |>
    filter(analyte == "temperature") |>
    pull(sample_id) |>
    unique()

  # All (sample_id, datetime, site_id) combos in the chemistry
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
                ": no Temperature measurements found — using fallback ",
                fallback_temp, " °C")
        return(mutate(gaps, temp_est = fallback_temp))
      }

      # Deduplicate by date first (average same-day readings)
      feat_temps <- feat_temps |>
        group_by(datetime) |>
        summarise(value = mean(value), .groups = "drop") |>
        arrange(datetime)

      t_num <- as.numeric(feat_temps$datetime)
      v     <- feat_temps$value
      g_num <- as.numeric(gaps$datetime)

      # Linear interpolation; rule=2 extrapolates using boundary values
      interp <- approx(t_num, v, xout = g_num, rule = 2)$y

      # Replace far-extrapolated values with fallback
      nearest_gap <- vapply(g_num, function(gd) min(abs(t_num - gd)), numeric(1))
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

  if (!requireNamespace("uuid", quietly = TRUE))
    stop("Package 'uuid' is required for writing temperature rows. ",
         "Install with: install.packages('uuid')")

  insert_rows <- new_rows |>
    mutate(
      uuid       = vapply(seq_len(n()), function(i) uuid::UUIDgenerate(), character(1)),
      uuid_lab   = .UUID_TEMP_LABMETHOD,
      value      = temp_est,
      quantified = TRUE,
      rl_low     = NA_real_,
      rl_high    = NA_real_,
      purpose    = "temperature_interpolated",
      comments   = sprintf("Interpolated by run_mspaf_pipeline.R on %s", Sys.Date())
    ) |>
    select(uuid, uuid_sample = sample_id, uuid_lab, value, quantified,
           rl_low, rl_high, purpose, comments)

  # Restrict to columns that exist in the DB table
  analysis_cols <- dbListFields(con, "analysis")
  insert_rows   <- select(insert_rows, any_of(analysis_cols))

  dbAppendTable(con, "analysis", insert_rows)
  message("  Wrote ", nrow(insert_rows), " temperature rows to analysis table")
  nrow(insert_rows)
}


#' Write chronic msPAF results back to the DB
#'
#' Creates a synthetic sample row for each (focal_date × feature) and an
#' analysis row containing the msPAF value.  Existing computed rows for the
#' same feature × focal_date are deleted first if replace_existing = TRUE.
#'
#' @param mspaf  Tibble with columns focal_date, site_id, value, ...
#' @param con     Open DuckDB connection (read-write)
#' @param feature_map  Tibble with uuid/name
#' @param replace_existing  Delete existing computed msPAF before writing?
#' @return Number of analysis rows written
.write_mspaf_to_db <- function(mspaf, con, feature_map, replace_existing) {

  if (!requireNamespace("uuid", quietly = TRUE))
    stop("Package 'uuid' is required for writing msPAF rows. ",
         "Install with: install.packages('uuid')")

  # Join feature UUIDs
  mspaf_w_uuid <- mspaf |>
    left_join(feature_map, by = c("site_id" = "name")) |>
    rename(feature_uuid = uuid)

  if (any(is.na(mspaf_w_uuid$feature_uuid)))
    warning("Some site_ids could not be matched to feature UUIDs — those rows skipped.")
  mspaf_w_uuid <- filter(mspaf_w_uuid, !is.na(feature_uuid))
  if (nrow(mspaf_w_uuid) == 0) return(0L)

  # Check for existing computed msPAF rows
  existing <- dbGetQuery(con, sprintf("
    SELECT s.uuid AS uuid_sample, s.uuid_feature, s.date AS focal_date
    FROM sample s
    JOIN analysis an ON an.uuid_sample = s.uuid
    WHERE s.purpose = 'chronic_mspaf_computed'
      AND an.uuid_lab = '%s'
      AND s.uuid_feature IN (%s)
  ", .UUID_MSPAF_LABMETHOD,
     paste(sprintf("'%s'", unique(mspaf_w_uuid$feature_uuid)), collapse = ", ")))

  if (nrow(existing) > 0) {
    if (replace_existing) {
      message("  Deleting ", nrow(existing), " existing computed msPAF rows …")
      sample_uuids <- paste(sprintf("'%s'", existing$uuid_sample), collapse = ", ")
      dbExecute(con, sprintf("DELETE FROM analysis WHERE uuid_sample IN (%s)", sample_uuids))
      dbExecute(con, sprintf("DELETE FROM sample   WHERE uuid         IN (%s)", sample_uuids))
    } else {
      # Skip rows already in DB
      existing_keys <- existing |>
        mutate(focal_date = as.Date(focal_date)) |>
        select(feature_uuid = uuid_feature, focal_date)
      mspaf_w_uuid <- mspaf_w_uuid |>
        anti_join(existing_keys, by = c("feature_uuid", "focal_date"))
      if (nrow(mspaf_w_uuid) == 0) {
        message("  All msPAF rows already exist (replace_existing = FALSE) — nothing written.")
        return(0L)
      }
      message("  Skipping ", nrow(existing), " existing rows; writing ",
              nrow(mspaf_w_uuid), " new rows …")
    }
  }

  # Build sample rows (one per focal_date × feature)
  sample_rows <- mspaf_w_uuid |>
    mutate(
      uuid         = vapply(seq_len(n()), function(i) uuid::UUIDgenerate(), character(1)),
      date         = focal_date,
      datetime     = as.POSIXct(paste(focal_date, "00:00:00")),
      organisation = "hydroSense",
      purpose      = "chronic_mspaf_computed",
      comments     = sprintf("Chronic msPAF computed by run_mspaf_pipeline.R on %s",
                             Sys.Date())
    ) |>
    select(uuid, uuid_feature = feature_uuid, date, datetime,
           organisation, purpose, comments)

  sample_cols        <- dbListFields(con, "sample")
  sample_rows_insert <- select(sample_rows, any_of(sample_cols))
  dbAppendTable(con, "sample", sample_rows_insert)

  # Build analysis rows (one per focal_date × feature)
  analysis_rows <- mspaf_w_uuid |>
    mutate(
      uuid_sample = sample_rows$uuid,
      uuid        = vapply(seq_len(n()), function(i) uuid::UUIDgenerate(), character(1)),
      uuid_lab    = .UUID_MSPAF_LABMETHOD,
      quantified  = TRUE,
      purpose     = "chronic_mspaf_computed",
      comments    = sprintf(
        "n_samples_in_window=%s%s",
        if ("n_samples_in_window" %in% names(mspaf_w_uuid))
          as.character(mspaf_w_uuid$n_samples_in_window) else "?",
        if ("n_imputed_in_window" %in% names(mspaf_w_uuid))
          paste0(" n_imputed=", mspaf_w_uuid$n_imputed_in_window) else ""
      )
    ) |>
    select(uuid, uuid_sample, uuid_lab, value, quantified, purpose, comments)

  analysis_cols        <- dbListFields(con, "analysis")
  analysis_rows_insert <- select(analysis_rows, any_of(analysis_cols))
  dbAppendTable(con, "analysis", analysis_rows_insert)

  nrow(analysis_rows_insert)
}


# ─────────────────────────────────────────────────────────────────────────────
# Example / quick-run block (edit and uncomment to run directly)
# ─────────────────────────────────────────────────────────────────────────────

# source("run_mspaf_pipeline.R")
#
# result <- run_mspaf_pipeline(
#   db_path            = "monitoring.duckdb",
#   # hydroSense_dir = "/path/to/hydroSense",  # uncomment if not installed
#   focal_features     = c("B.S01"),                 # NULL = all non-reference features
#   reference_features = c("B.S03"),                 # NULL = auto-detect from DB
#   date_range         = c(as.Date("2017-01-01"), Sys.Date()),
#   tau_days           = 90,
#   window_days        = 365,
#   ref_summary        = "geom_mean",
#   impute             = TRUE,
#   impute_iter        = 2000,
#   impute_chains      = 4,
#   refit_model        = FALSE,   # set TRUE to discard cached model and refit
#   write_back         = TRUE,
#   replace_existing   = FALSE
# )
#
# # Inspect results
# print(result$mspaf)
# ggplot2::ggplot(result$mspaf, ggplot2::aes(focal_date, value)) +
#   ggplot2::geom_line() +
#   ggplot2::facet_wrap(~site_id)
