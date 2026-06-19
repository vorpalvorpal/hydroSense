#!/usr/bin/env Rscript
# Generate golden snapshots of the CURRENT add_mspaf / mspaf_daily output, for
# the issue #30 equivalence harness. Run ONCE on the pre-change code; the saved
# fixture is the reference the vectorised engine must reproduce (numeric rel 1e-9,
# categorical/integer exact). Self-contained: bundled SSD data, no DB/guideline_dir.
#
#   Rscript dev/gen_mspaf_golden.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(30L)

FIX <- "tests/testthat/fixtures"
dir.create(FIX, showWarnings = FALSE, recursive = TRUE)

## ── add_mspaf input: multi-sample, with edge cases ──────────────────────────
co <- c(pH = 7.5, DOC = 2.0, Ca = 6.0, Mg = 4.0, hardness = 30.0)
mk_sample <- function(sid, cu, zn, ni, drop_doc = FALSE, cu_bdl = FALSE) {
  an <- c("Cu", "Zn", "Ni", names(co))
  va <- c(cu, zn, ni, unname(co))
  det <- rep(TRUE, length(an)); if (cu_bdl) det[an == "Cu"] <- FALSE
  keep <- if (drop_doc) an != "DOC" else rep(TRUE, length(an))   # missing co-analyte
  tibble::tibble(sample_id = sid, site_id = "f1",
                 datetime = as.Date("2024-01-01") + match(sid, paste0("s", 1:9)) - 1L,
                 analyte = an[keep], value = va[keep], detected = det[keep])
}
pt_df <- dplyr::bind_rows(
  mk_sample("s1", 5,   10, 0.3),
  mk_sample("s2", 50,  20, 0.1),
  mk_sample("s3", 0.5, 2,  0.05, drop_doc = TRUE),     # Cu dropped (no DOC)
  mk_sample("s4", 8,   15, 0.2,  cu_bdl   = TRUE)      # Cu BDL -> 0
)

## draw-bearing version: metals get draw_id 1..N (lognormal jitter), co exact.
N <- 6L
mk_draws <- function(pt) {
  metals <- pt |> dplyr::filter(.data$analyte %in% c("Cu", "Zn", "Ni"))
  cofix  <- pt |> dplyr::filter(!.data$analyte %in% c("Cu", "Zn", "Ni")) |>
    dplyr::mutate(draw_id = NA_integer_)
  md <- metals |> dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::group_modify(\(.x, .y) {
      tibble::tibble(site_id = .x$site_id[1], datetime = .x$datetime[1],
                     value = .x$value[1] * exp(stats::rnorm(N, 0, 0.2)),
                     detected = .x$detected[1], draw_id = seq_len(N))
    }) |> dplyr::ungroup()
  dplyr::bind_rows(md, cofix)
}
dr_df <- mk_draws(pt_df)

## A small reference chemistry frame (ARA on)
ref_df <- dplyr::bind_rows(
  mk_sample("r1", 1, 3, 0.02), mk_sample("r2", 2, 4, 0.03), mk_sample("r3", 1.5, 2, 0.01)
) |> dplyr::mutate(site_id = "ref")

normalise_out <- function(out) {
  rows <- dplyr::filter(out, .data$analyte == "msPAF")
  scal_cols <- intersect(c("sample_id", "draw_id", "value", "n_analytes_used",
                           "n_analytes_imputed", "dominant_analyte", "max_paf"),
                         names(rows))
  scal <- dplyr::arrange(rows[, scal_cols],
                         dplyr::across(dplyr::any_of(c("sample_id", "draw_id"))))
  brk <- if ("analyte_pafs" %in% names(rows)) {
    tidyr::unnest(dplyr::select(rows, dplyr::any_of(c("sample_id", "draw_id")),
                               "analyte_pafs"), cols = "analyte_pafs") |>
      dplyr::arrange(dplyr::across(dplyr::any_of(c("sample_id", "draw_id", "analyte"))))
  } else NULL
  list(scalars = scal, breakdown = brk, ara_summary = attr(out, "ara_summary"))
}

prep_ref <- suppressMessages(prepare_reference(ref_df, conc_units = "ug/L"))
run <- function(df, reference) suppressMessages(
  add_mspaf(df, reference = reference, conc_units = "ug/L",
             return = if ("draw_id" %in% names(df)) "draws" else "summary"))

golden <- list(
  pt_noara  = normalise_out(run(pt_df, NULL)),
  pt_ara    = normalise_out(run(pt_df, prep_ref)),
  dr_noara  = normalise_out(run(dr_df, NULL)),
  dr_ara    = normalise_out(run(dr_df, prep_ref)),
  inputs    = list(pt_df = pt_df, dr_df = dr_df, ref_df = ref_df, N = N)
)
qs2::qs_save(golden, file.path(FIX, "mspaf_golden.qs2"))
cat("WROTE", file.path(FIX, "mspaf_golden.qs2"),
    "- scalars rows:", nrow(golden$dr_ara$scalars), "\n")
