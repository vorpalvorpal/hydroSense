#!/usr/bin/env Rscript
# Build SSD validation table: our calculated DGVs vs published reference values.
#
# Reference DGV sources (priority order):
#   NH3-N  : ANZG 2026, "PCx at pH 7" sheet, pH 7 / 20°C row (mg N/L → ×1000 µg/L)
#   Cu     : ANZG 2023, "DGVs" sheet, DOC = 0.5 mg/L row (µg/L)
#   B      : ANZG 2021 (= ANZECC 2000; µg B/L)
#   Cr     : ANZG 2026 Cr(III) (µg/L)
#   Ni     : ANZG 2024 (µg Ni/L)
#   NO3-N  : ANZG 2025 soft water (mg NO3-N/L → ×1000 µg/L)
#   Zn     : ANZG 2024 (µg Zn/L)
#   Others : anzecc_analyte_metadata.csv dgv_Npct_ug_L columns (ANZECC 2000 guidelineDF)
#
# Calculation formula (applied after fitting):
#   If data were pre-divided by ACR before fitting (acr_applied = TRUE):
#       calc_N = HC_N   (already on chronic-equivalent scale)
#   If NOT pre-divided but trigger_divisor > 1 (e.g. BaP safety factor):
#       calc_N = HC_N / trigger_divisor
#   Otherwise:
#       calc_N = HC_N
#
# Output: leachatetools/data-raw/ssd_validation.csv
#
# Usage (run from the dashboard/ root): Rscript leachatetools/scripts/make_ssd_validation.R

suppressPackageStartupMessages({
  library(ssdtools)
  library(dplyr)
  library(readr)
  library(qs)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

PKG_DIR <- "leachatetools"
RAW_DIR <- file.path(PKG_DIR, "data-raw")
OUT_DIR <- file.path(PKG_DIR, "data-raw")
SSD_DIR <- file.path(PKG_DIR, "data-raw", "ssd_models")

safe_name <- function(x) gsub("[^A-Za-z0-9]", "_", x)
sig3      <- function(x) signif(x, 3)

# ── Metadata ──────────────────────────────────────────────────────────────────

meta_all <- read_csv(
  file.path(RAW_DIR, "anzecc_analyte_metadata.csv"),
  col_types = cols(
    acr                  = col_double(),
    anzecc_hc_percentile = col_integer(),
    trigger_divisor      = col_double(),
    n_min_override       = col_integer(),
    fit_dist             = col_character(),
    dgv_99pct_ug_L       = col_double(),
    dgv_95pct_ug_L       = col_double(),
    dgv_90pct_ug_L       = col_double(),
    dgv_80pct_ug_L       = col_double(),
    .default             = col_character()
  ),
  show_col_types = FALSE
) |>
  mutate(ssd_available = as.logical(ssd_available))

# ── Reference DGV table ───────────────────────────────────────────────────────
# ANZG analytes: hardcoded from technical brief XLSX files or confirmed by user.
# All concentrations in µg/L.

# NH3-N: ANZG 2026 "PCx at pH 7" sheet, temp = 20°C row
#   Columns: PC99=0.320, PC95=0.970, PC90=1.648, PC80=3.047 (mg total-NH3-N/L → ×1000 µg/L)
nh3_ref <- tibble(
  analyte    = "NH3-N",
  ref_99     = 320,   ref_95 = 970,  ref_90 = 1648, ref_80 = 3047,
  dgv_source = "ANZG2026 PCx pH7 20°C (µg total-NH3-N/L)"
)

# Cu: ANZG 2023 "DGVs" sheet, DOC = 0.5 mg/L row (µg/L)
cu_ref <- tibble(
  analyte    = "Cu",
  ref_99     = 0.20,  ref_95 = 0.47, ref_90 = 0.73, ref_80 = 1.30,
  dgv_source = "ANZG2023 DGVs sheet DOC=0.5 mg/L (µg/L)"
)

# B: ANZG 2021 — identical to ANZECC 2000 values
b_ref <- tibble(
  analyte    = "B",
  ref_99     = 340,   ref_95 = 940,  ref_90 = 1500, ref_80 = 2500,
  dgv_source = "ANZG2021 (= ANZECC2000; µg B/L)"
)

# Cr(III): ANZG 2026
cr_ref <- tibble(
  analyte    = "Cr",
  ref_99     = 0.95,  ref_95 = 6.7,  ref_90 = 16,   ref_80 = 39,
  dgv_source = "ANZG2026 Cr(III) (µg/L)"
)

# Ni: ANZG 2024 (µg Ni/L)
ni_ref <- tibble(
  analyte    = "Ni",
  ref_99     = 0.31,  ref_95 = 2.0,  ref_90 = 4.6,  ref_80 = 10,
  dgv_source = "ANZG2024 Ni (µg Ni/L)"
)

# NO3-N: three separate SSDs, one per hardness class (ANZG 2025)
no3_soft_ref <- tibble(
  analyte    = "NO3-N_soft",
  ref_99     = 640,   ref_95 = 1100,  ref_90 = 1500,  ref_80 = 2300,
  dgv_source = "ANZG2025 NO3-N soft (<30 mg/L CaCO3) (µg NO3-N/L)"
)
no3_mod_ref <- tibble(
  analyte    = "NO3-N_mod",
  ref_99     = 1000,  ref_95 = 2600,  ref_90 = 4200,  ref_80 = 7100,
  dgv_source = "ANZG2025 NO3-N moderate (30-150 mg/L CaCO3) (µg NO3-N/L)"
)
no3_hard_ref <- tibble(
  analyte    = "NO3-N_hard",
  ref_99     = 18000, ref_95 = 29000, ref_90 = 38000, ref_80 = 56000,
  dgv_source = "ANZG2025 NO3-N hard (>150 mg/L CaCO3) (µg NO3-N/L)"
)

# Zn: ANZG 2024 (µg Zn/L)
zn_ref <- tibble(
  analyte    = "Zn",
  ref_99     = 1.5,   ref_95 = 4.1,  ref_90 = 6.8,  ref_80 = 12,
  dgv_source = "ANZG2024 Zn (µg Zn/L)"
)

# All remaining: ANZECC 2000 from guidelineDF (via metadata)
other_refs <- meta_all |>
  filter(!analyte %in% c("NH3-N", "Cu", "B", "Cr", "Ni", "NO3-N_soft", "NO3-N_mod", "NO3-N_hard", "Zn")) |>
  transmute(
    analyte,
    ref_99     = dgv_99pct_ug_L,
    ref_95     = dgv_95pct_ug_L,
    ref_90     = dgv_90pct_ug_L,
    ref_80     = dgv_80pct_ug_L,
    dgv_source = "ANZECC2000 guidelineDF"
  )

refs <- bind_rows(nh3_ref, cu_ref, b_ref, cr_ref, ni_ref,
                  no3_soft_ref, no3_mod_ref, no3_hard_ref,
                  zn_ref, other_refs)

# ── HC values from fitted SSDs ─────────────────────────────────────────────────

qs_files <- list.files(SSD_DIR, pattern = "\\.qs$", full.names = TRUE)
message("Found ", length(qs_files), " .qs model files in ", SSD_DIR)

process_fit <- function(f) {
  fit <- tryCatch(qread(f), error = function(e) {
    message("  ERROR reading ", basename(f), ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)

  analyte_sn <- sub("\\.qs$", "", basename(f))
  meta_row   <- meta_all |> filter(safe_name(analyte) == analyte_sn)
  if (nrow(meta_row) == 0) {
    message("  No metadata match for file: ", basename(f)); return(NULL) }
  meta_row <- meta_row[1L, ]

  # Skip analytes where ssd_available = FALSE (old .qs files from prior runs)
  if (!isTRUE(meta_row$ssd_available)) {
    message("  SKIP (ssd_available=FALSE): ", meta_row$analyte); return(NULL) }


  hc_tbl <- tryCatch(
    # ssdtools >=2.0 uses `proportion` (0.01, 0.05, …); older used `percent` (1, 5, …)
    ssd_hc(fit, proportion = c(0.01, 0.05, 0.10, 0.20), ci = FALSE),
    error = function(e) {
      message("  ssd_hc error for ", meta_row$analyte, ": ", conditionMessage(e)); NULL })
  if (is.null(hc_tbl)) return(NULL)

  # One model-averaged row per proportion; name as hc1/hc5/hc10/hc20
  hc <- setNames(hc_tbl$est, paste0("hc", round(hc_tbl$proportion * 100)))

  acr_applied <- isTRUE(attr(fit, "acr_applied"))
  td          <- meta_row$trigger_divisor

  # Post-fit divisor: only when data was NOT pre-divided (e.g. BaP trigger_divisor=10)
  post_div <- if (!acr_applied && !is.na(td) && td > 1) td else 1L

  message("  ", meta_row$analyte,
          " | dists=", attr(fit, "dists_used") %||% "?",
          " | acr_applied=", acr_applied,
          " | post_div=", post_div)

  tibble(
    analyte         = meta_row$analyte,
    reliability     = meta_row$anzecc_reliability,
    data_type       = meta_row$anzecc_data_type,
    acr_applied     = acr_applied,
    trigger_divisor = td,
    hc_pct_anzecc   = meta_row$anzecc_hc_percentile,  # 1 or 5
    dists_used      = attr(fit, "dists_used") %||% NA_character_,
    calc_99         = hc["hc1"]  / post_div,
    calc_95         = hc["hc5"]  / post_div,
    calc_90         = hc["hc10"] / post_div,
    calc_80         = hc["hc20"] / post_div
  )
}

hc_df <- bind_rows(lapply(qs_files, process_fit))
message("\nComputed HC values for ", nrow(hc_df), " analytes")

# ── Join references, compute ratios ───────────────────────────────────────────

val <- hc_df |>
  left_join(refs, by = "analyte") |>
  mutate(
    ratio_99 = round(calc_99 / ref_99, 3),
    ratio_95 = round(calc_95 / ref_95, 3),
    ratio_90 = round(calc_90 / ref_90, 3),
    ratio_80 = round(calc_80 / ref_80, 3)
  )

# ── Flag ──────────────────────────────────────────────────────────────────────
# Primary comparison level:
#   anzecc_hc_percentile == 1  → ANZECC adopted HC1 as basis → compare at 99% level
#   anzecc_hc_percentile == 5  → compare at 95% level (standard)
#
# Flag meanings:
#   ok                       — primary ratio within 0.5–2×
#   ANZG_vs_ANZECC2000       — ANZG XLSX analyte; reference is ANZECC 2000 (expected diff)
#   lnorm_not_yet_burrIII3   — ratio outside 0.5–2× but analyte uses lnorm; will improve
#                              once models are refit with burrIII3 (all Warne2000 analytes)
#   large_discrepancy        — ratio 3–5× even accounting for fitting method; needs review
#   DATA_ISSUE               — ratio >5× or <0.2×; likely a data extraction problem
#   partial_reference        — reference DGV only available at one protection level;
#                              ratio shown is at the level where ref is available
#   no_reference             — no published reference DGV in any column

val <- val |>
  mutate(
    # Primary ratio: at the protection level ANZECC/ANZG used as their basis
    primary_ratio = if_else(!is.na(hc_pct_anzecc) & hc_pct_anzecc == 1L,
                             ratio_99, ratio_95),

    # Fallback ratio: for analytes where the primary level has no reference but
    # another level does (e.g. BaP: only ref_90 is available)
    fallback_ratio = case_when(
      !is.na(primary_ratio)            ~ NA_real_,   # no fallback needed
      !is.na(ratio_90)                 ~ ratio_90,
      !is.na(ratio_99)                 ~ ratio_99,
      !is.na(ratio_80)                 ~ ratio_80,
      TRUE                             ~ NA_real_
    ),

    # Effective ratio used for flagging
    eff_ratio = coalesce(primary_ratio, fallback_ratio),

    # Warne2000 analytes currently fitted with lnorm instead of burrIII3
    is_lnorm_anzecc = dists_used == "lnorm" & !grepl("ANZG", dgv_source),

    # Magnitude bounds beyond which fitting-method differences can NOT explain the gap:
    #   HC5-primary analytes: lnorm vs BurrIII at HC5 typically within 0.2–5×
    #   HC1-primary analytes: lnorm vs BurrIII at HC1 can reach 0.05–15× due to tail
    beyond_fitting_method = case_when(
      !is.na(hc_pct_anzecc) & hc_pct_anzecc == 1L ~
        (eff_ratio > 15 | eff_ratio < 0.05),
      TRUE ~                                          # HC5-based (or unknown)
        (eff_ratio > 5  | eff_ratio < 0.20)
    ),

    flag = case_when(
      grepl("no ANZG DGV", dgv_source)              ~ "ANZG_vs_ANZECC2000_by_design",
      is.na(eff_ratio)                               ~ "no_reference",
      !is.na(fallback_ratio)                         ~ "partial_reference",
      # Data issues: too extreme for any fitting-method explanation
      beyond_fitting_method                          ~ "DATA_ISSUE",
      # Fitting-method gap: lnorm where burrIII3 is not yet used
      (eff_ratio > 2 | eff_ratio < 0.50) & is_lnorm_anzecc ~
                                                       "lnorm_not_yet_burrIII3",
      # Residual discrepancy with the right distribution (or ANZG analytes)
      eff_ratio > 3  | eff_ratio < 0.33             ~ "large_discrepancy",
      eff_ratio > 2  | eff_ratio < 0.50             ~ "fitting_method_diff",
      TRUE                                           ~ "ok"
    )
  )

# ── Format and write ──────────────────────────────────────────────────────────

val_out <- val |>
  mutate(across(c(calc_99, calc_95, calc_90, calc_80,
                  ref_99,  ref_95,  ref_90,  ref_80),  sig3)) |>
  select(
    analyte, reliability, data_type, acr_applied, dists_used, dgv_source,
    calc_99, ref_99,  ratio_99,
    calc_95, ref_95,  ratio_95,
    calc_90, ref_90,  ratio_90,
    calc_80, ref_80,  ratio_80,
    flag
  ) |>
  arrange(flag, reliability, analyte)

out_path <- file.path(OUT_DIR, "ssd_validation.csv")
write_csv(val_out, out_path, na = "")

message("\nWrote: ", out_path, "  (", nrow(val_out), " rows)")
message("\nFlag summary:")
print(as.data.frame(count(val_out, flag)), row.names = FALSE)

message("\nRows needing attention:")
attention <- val_out |>
  filter(!flag %in% c("ok", "ANZG_vs_ANZECC2000_by_design")) |>
  select(analyte, reliability, flag, ratio_95, ratio_99)
if (nrow(attention) == 0) {
  message("  (none)")
} else {
  print(as.data.frame(attention), row.names = FALSE)
}
