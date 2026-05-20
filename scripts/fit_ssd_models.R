#!/usr/bin/env Rscript
# Fit SSD models for all ANZECC/ANZG analytes and save as .qs files.
#
# Data sources:
#   leachatetools/data-raw/anzecc_warne2000_observations.csv  — Warne (2000) raw values
#   leachatetools/data-raw/anzecc_analyte_metadata.csv        — per-analyte fit parameters
#   leachatetools/guideline data/*.xlsx                        — ANZG modern technical briefs
#
# Design:
#   - Per-analyte distribution: the 'fit_dist' column in anzecc_analyte_metadata.csv
#     controls which ssdtools distribution(s) to use for each analyte.  A single
#     name (e.g. "burrIII3", "lnorm") uses that distribution alone; a
#     comma-separated list (e.g. "lnorm,llogis,gamma") fits all listed dists and
#     model-averages.  If fit_dist is empty/NA the global --dists default is used.
#   - Warne 2000 analytes use "burrIII3" (Burr Type III — what Warne 2000 used).
#   - ANZG XLSX analytes use "lnorm" (modern default).
#   - For MR analytes: all acute values are divided by ACR before fitting,
#     so the resulting SSD is on a chronic-equivalent concentration scale.
#     ssd_hp(fit, conc) then gives PAF directly for msPAF without further ACR
#     adjustment.
#   - For analytes where ANZECC adopted HC1 (not HC5) due to bioaccumulation
#     or acute toxicity concerns (Se, DDT, BaP, Azinphos-methyl), this is
#     recorded in the metadata; the SSD itself is fitted normally.
#
# Usage (run from the dashboard/ root):
#   Rscript leachatetools/scripts/fit_ssd_models.R [--dists lnorm|multi]
#   (--dists sets the fallback default when fit_dist is not specified per-analyte)
#
# Output: leachatetools/data-raw/ssd_models/<safe_name>.qs  (one per analyte)

suppressPackageStartupMessages({
  library(ssdtools)
  library(readxl)
  library(dplyr)
  library(readr)
  library(qs)
})

# ── Config ────────────────────────────────────────────────────────────────────

args   <- commandArgs(trailingOnly = TRUE)
DISTS  <- if ("--dists" %in% args) {
  args[which(args == "--dists") + 1]
} else "lnorm"   # "lnorm" to match ANZECC; "multi" for model-averaged

PKG_DIR   <- "leachatetools"
RAW_DIR   <- file.path(PKG_DIR, "data-raw")
GUIDE_DIR <- file.path(PKG_DIR, "guideline data")
OUT_DIR   <- file.path(RAW_DIR, "ssd_models")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DIST_SETS <- list(
  lnorm    = "lnorm",
  burrIII3 = "burrIII3",
  multi    = c("lnorm", "llogis", "burrIII3", "lgumbel", "gamma", "weibull")
)
dists_default <- DIST_SETS[[DISTS]]
if (is.null(dists_default)) stop("--dists must be 'lnorm', 'burrIII3', or 'multi'")

message("=== fit_ssd_models.R  |  dists = ", DISTS, " ===\n")

# ── Helpers ───────────────────────────────────────────────────────────────────

safe_name <- function(x) gsub("[^A-Za-z0-9]", "_", x)

fit_and_save <- function(analyte, df, meta, source_note) {
  n_min <- if (!is.na(meta$n_min_override)) meta$n_min_override else 6L
  fname <- file.path(OUT_DIR, paste0(safe_name(analyte), ".qs"))

  if (nrow(df) < n_min) {
    message("  SKIP (n=", nrow(df), " < ", n_min, "): ", analyte)
    return(invisible(NULL))
  }

  is_acute   <- meta$anzecc_data_type %in% c("acute_LC50", "acute_EC50")
  acute_flag <- isTRUE(is_acute) && !is.na(meta$trigger_divisor) &&
                meta$trigger_divisor > 1

  # Resolve per-analyte distribution: metadata fit_dist overrides --dists default.
  # fit_dist may be:
  #   - a DIST_SETS alias ("lnorm", "burrIII3", "multi") → expanded via lookup
  #   - a comma-separated list of dist names ("lnorm,llogis,gamma") → split directly
  #   - empty/NA → use global --dists default
  fd_raw <- meta$fit_dist
  per_analyte_dists <- if (!is.na(fd_raw) && nzchar(trimws(fd_raw))) {
    fd_trimmed <- trimws(fd_raw)
    if (!is.null(DIST_SETS[[fd_trimmed]])) {
      DIST_SETS[[fd_trimmed]]                        # expand alias
    } else {
      trimws(strsplit(fd_trimmed, ",")[[1]])          # comma-separated list
    }
  } else {
    dists_default          # global fallback from --dists argument
  }
  dists_label <- paste(per_analyte_dists, collapse = "+")

  message("  Fitting: ", analyte,
          " (n=", nrow(df), ", dists=", dists_label,
          if (acute_flag) paste0(", ACR-adjusted /", meta$trigger_divisor) else "",
          ")")

  # For MR analytes: divide values by ACR so SSD is on chronic-equivalent scale
  if (acute_flag) {
    df <- df |> mutate(Conc = Conc / meta$trigger_divisor)
  }

  fit <- tryCatch(
    ssd_fit_dists(df, left = "Conc", dists = per_analyte_dists),
    error = function(e) { message("    ERROR: ", conditionMessage(e)); NULL }
  )
  if (is.null(fit)) return(invisible(NULL))

  attr(fit, "source_note")       <- source_note
  attr(fit, "acute_data")        <- is_acute
  attr(fit, "acr")               <- meta$acr
  attr(fit, "anzecc_hc_pct")     <- meta$anzecc_hc_percentile
  attr(fit, "trigger_divisor")   <- meta$trigger_divisor
  attr(fit, "acr_applied")       <- acute_flag
  attr(fit, "dists_used")        <- dists_label
  attr(fit, "fit_dist_source")   <- if (!is.na(fd_raw) && nzchar(trimws(fd_raw)))
                                       "metadata" else "cli_default"

  qsave(fit, fname)
  message("    Saved: ", fname)
  invisible(fit)
}

# ── Load metadata ─────────────────────────────────────────────────────────────

meta_all <- read_csv(
  file.path(RAW_DIR, "anzecc_analyte_metadata.csv"),
  col_types = cols(
    acr                  = col_double(),
    anzecc_hc_percentile = col_integer(),
    trigger_divisor      = col_double(),
    n_min_override       = col_integer(),
    fit_dist             = col_character(),   # per-analyte ssdtools distribution(s)
    .default             = col_character()
  ),
  show_col_types = FALSE
) |>
  mutate(ssd_available = as.logical(ssd_available))

obs_all <- read_csv(
  file.path(RAW_DIR, "anzecc_warne2000_observations.csv"),
  col_types = cols(value_ug_L = col_double(), .default = col_character()),
  show_col_types = FALSE
)

# ── 1. Modern ANZG XLSX analytes ──────────────────────────────────────────────

message("=== Modern ANZG XLSX analytes ===")

read_col_by_index <- function(path, sheet, skip, species_col, conc_col,
                               data_row_offset = 1, units_factor = 1,
                               media_col = NULL, freshwater_only = FALSE) {
  raw  <- read_excel(path, sheet = sheet, col_names = FALSE, skip = skip)
  data <- raw[seq(data_row_offset + 1, nrow(raw)), ]
  get_col <- function(df, idx) suppressWarnings(as.character(unlist(df[, idx])))
  species <- get_col(data, species_col)
  conc    <- suppressWarnings(as.numeric(get_col(data, conc_col))) * units_factor
  df <- tibble(Conc = conc, Species = species)
  if (freshwater_only && !is.null(media_col)) {
    media <- get_col(data, media_col)
    df    <- df[grepl("fresh", media, ignore.case = TRUE), ]
  }
  df |> filter(!is.na(Conc), Conc > 0, !is.na(Species), Species != "NA")
}

read_min_per_species <- function(path, sheet, skip, species_col, conc_col,
                                  data_row_offset = 1, units_factor = 1,
                                  media_col = NULL, freshwater_only = FALSE) {
  df <- read_col_by_index(path, sheet, skip, species_col, conc_col,
                           data_row_offset, units_factor, media_col, freshwater_only)
  df |> group_by(Species) |> summarise(Conc = min(Conc), .groups = "drop")
}

xlsx_fits <- list(

  list(analyte = "NH3-N",
       fn = function() read_col_by_index(
         file.path(GUIDE_DIR, "ammonia-fresh-dgvs-data-entry.xlsx"),
         sheet = 1, skip = 3, species_col = 4, conc_col = 55, units_factor = 1000),
       note = "ANZG 2026 ammonia; col 55 'Lowest value per species' (mg NH3-N/L at pH7 20°C -> µg/L)"),

  # NO3-N is split into 3 separate SSDs because hardness physically modifies
  # toxicity to the same species (not just a species-assemblage effect).
  # Each hardness class uses col 49 'LOWEST VALUE FOR SPECIES (mg NO3-N/L)'
  # from its own XLSX sheet, converted to µg/L.  ANZG used Burrlioz 2.0:
  # inverse Weibull for soft and hard water → lgumbel in ssdtools (log-Gumbel
  # = Fréchet = inverse Weibull; same distributional family, different name).
  # Burr Type III for moderately hard water → burrIII3.

  list(analyte = "NO3-N_soft",
       fn = function() read_col_by_index(
         file.path(GUIDE_DIR, "nitrate-fresh-dgvs-data-entry.xlsx"),
         sheet = "Nitrate - soft water", skip = 9, species_col = 5,
         conc_col = 49, data_row_offset = 2, units_factor = 1000),
       note = "ANZG 2025 nitrate soft water (< 30 mg/L CaCO3); 14 species; col 49 'Lowest value per species' (mg/L -> µg/L)"),

  list(analyte = "NO3-N_mod",
       fn = function() read_col_by_index(
         file.path(GUIDE_DIR, "nitrate-fresh-dgvs-data-entry.xlsx"),
         sheet = "Nitrate - moderately hard water", skip = 9, species_col = 5,
         conc_col = 49, data_row_offset = 2, units_factor = 1000),
       note = "ANZG 2025 nitrate moderately hard water (30-150 mg/L CaCO3); 11 species; col 49 'Lowest value per species' (mg/L -> µg/L)"),

  list(analyte = "NO3-N_hard",
       fn = function() read_col_by_index(
         file.path(GUIDE_DIR, "nitrate-fresh-dgvs-data-entry.xlsx"),
         sheet = "Nitrate - hard water", skip = 9, species_col = 5,
         conc_col = 49, data_row_offset = 2, units_factor = 1000),
       note = "ANZG 2025 nitrate hard water (> 150 mg/L CaCO3); 12 species; col 49 'Lowest value per species' (mg/L -> µg/L)"),

  list(analyte = "B",
       fn = function() read_col_by_index(
         file.path(GUIDE_DIR, "boron_fresh_dgv_data-entry_final.xlsx"),
         sheet = 1, skip = 6, species_col = 5, conc_col = 49, units_factor = 1000),
       note = "ANZG 2021 boron; col 49 (values in mg/L despite label -> µg/L)"),

  list(analyte = "Cr",
       fn = function() read_col_by_index(
         file.path(GUIDE_DIR, "chromium-III-fresh-dgvs-data-entry.xlsx"),
         sheet = 1, skip = 6, species_col = 5, conc_col = 48,
         data_row_offset = 2, units_factor = 1, media_col = 4, freshwater_only = TRUE),
       note = "ANZG 2026 Cr(III); col 48 'Lowest value per species (µg/L)'; freshwater only"),

  list(analyte = "Cu",
       fn = function() read_col_by_index(
         file.path(GUIDE_DIR, "copper-fresh-dgvs-data-entry.xlsx"),
         sheet = "Accepted Chronic Data", skip = 5,
         species_col = 7, conc_col = 57, units_factor = 1),
       note = "ANZG 2023 Cu draft; col 57 'Best value per species (µg/L)'"),

  list(analyte = "Ni",
       fn = function() {
         # Ni requires MLR bioavailability normalisation (Peters et al. 2021) to
         # index condition (pH 7.5, Ca 6 mg/L, Mg 4 mg/L, DOC 0.5 mg/L) before
         # fitting the SSD — the raw XLSX col 31 holds measured concentrations, not
         # normalised values.  The 26 normalised negligible-effect values used by
         # ANZG (Stauber et al. 2021, Table 3 of the technical brief) are
         # pre-tabulated in data-raw/ni_mlr_normalised_table3.csv.
         read_csv(
           file.path(RAW_DIR, "ni_mlr_normalised_table3.csv"),
           col_types = cols(Conc_ug_L = col_double(), .default = col_character()),
           show_col_types = FALSE
         ) |>
           transmute(Conc = Conc_ug_L, Species = Species)
       },
       note = "ANZG 2024 Ni; 26 MLR-normalised negligible-effect values (Stauber et al. 2021 Table 3, index: pH 7.5 Ca 6 Mg 4 DOC 0.5 mg/L); Burrlioz 2.0"),

  list(analyte = "Zn",
       fn = function() read_min_per_species(
         file.path(GUIDE_DIR, "zinc-fresh-dgvs-data-entry.xlsx"),
         sheet = "ForWordDoc", skip = 3, species_col = 2, conc_col = 11),
       note = "ANZG 2024 Zn draft; per-species minimum of col 11 'Predicted EC10 (µg/L at reference)'")
)

for (item in xlsx_fits) {
  analyte <- item$analyte
  message("Processing: ", analyte)
  meta <- meta_all |> filter(analyte == !!analyte)
  tryCatch({
    df  <- item$fn()
    df  <- df |> select(Conc, Species) |> filter(!is.na(Conc), Conc > 0)
    fit_and_save(analyte, df, meta, item$note)
  }, error = function(e) message("  ERROR (", analyte, "): ", conditionMessage(e)))
}

# ── 2. ANZECC 2000 analytes (Warne 2000) ─────────────────────────────────────

message("\n=== ANZECC 2000 analytes (Warne 2000) ===")

warne_analytes <- meta_all |>
  filter(data_source == "Warne2000", ssd_available == TRUE,
         !is.na(observations_analyte))

for (i in seq_len(nrow(warne_analytes))) {
  meta    <- warne_analytes[i, ]
  analyte <- meta$analyte
  obs_key <- meta$observations_analyte

  message("Processing: ", analyte, " (obs key: '", obs_key, "')")

  df <- obs_all |>
    filter(analyte == obs_key, !is.na(value_ug_L), value_ug_L > 0) |>
    transmute(Conc = value_ug_L, Species = species_id)

  note <- paste0(
    "Warne (2000) p.", meta$source_page, ". ",
    meta$anzecc_reliability, " (", meta$anzecc_data_type, "). ",
    if (!is.na(meta$acr)) paste0("ACR=", meta$acr, ". ") else "",
    meta$notes
  )

  fit_and_save(analyte, df, meta, note)
}

message("\n=== Done. Models saved to: ", OUT_DIR, " ===")
message(paste(" ", sort(list.files(OUT_DIR)), collapse = "\n"))
