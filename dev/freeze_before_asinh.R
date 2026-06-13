## Capture the engine invariant (add_amspaf, untouched by #15) + grab dates on
## B.S01, under the CURRENT additive model, for before/after comparison.
suppressMessages({library(dplyr); devtools::load_all(".", quiet=TRUE)})
options(leachatetools.guideline_dir = "guideline data")
cc <- qs2::qs_read("test data/bs01_v3_cache.qs2")
grabs <- cc$target_chem

## INVARIANT 1: per-sample add_amspaf on the grab chemistry (point mode). #15
## does not touch R/mspaf.R, so this must be bit-identical afterwards.
eng_tot <- suppressMessages(add_amspaf(grabs, reference = NULL)) |>
  dplyr::filter(.data$analyte == "AmsPAF") |>
  dplyr::select(sample_id, value) |> dplyr::arrange(sample_id)

grab_dates <- sort(unique(as.Date(grabs$datetime)))

qs2::qs_save(list(
  engine_amspaf_total = eng_tot,
  grab_dates          = grab_dates,
  note = "B.S01 under ADDITIVE model, frozen before #15 asinh transform",
  generated = Sys.time()
), "dev/before_asinh/engine_invariant.qs2")
cat("WROTE dev/before_asinh/engine_invariant.qs2  (n grabs:", nrow(eng_tot),
    " grab dates:", length(grab_dates), ")\n")
