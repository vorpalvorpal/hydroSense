## data-raw/ssd_paf_lookup.R — generate inst/extdata/ssd_paf_lookup.qs2
## Run once: Rscript data-raw/ssd_paf_lookup.R
## Output: inst/extdata/ssd_paf_lookup.qs2 (shipped with the package)
##
## Algorithm (per approved plan, issue #36):
##   For each (method, analyte) pair with a fitdists object:
##   1. Forward-scan the full concentration range to locate the effective
##      support of the SSD CDF (avoid ssd_hc() which underflows at extremes).
##   2. Build a monotone-cubic spline (monoH.FC) over an adaptive grid;
##      double knot count until max|err| < 1e-8 on 8 000 random check points.
##   3. Store the knot grid and PAF values (not the closure) so the file is
##      self-contained and fast to deserialise.

suppressMessages(devtools::load_all(".", quiet = TRUE))

## ── helpers ─────────────────────────────────────────────────────────────────

## Deduplicate ssd_hp() output and realign to the requested concentration
## vector (ssd_hp() may drop or reorder rows for extreme values).
.align_hp <- function(fit, conc) {
  raw <- ssdtools::ssd_hp(fit, conc = conc, ci = FALSE, proportion = TRUE)
  raw <- raw[!duplicated(raw$conc), ]
  raw$est[match(conc, raw$conc)]
}

## Locate the effective support of the SSD CDF by a forward scan.
## Returns a named list(log10_lo, log10_hi).
.scan_range <- function(fit) {
  scan_conc <- 10^seq(log10(1e-6), log10(1e9), length.out = 4000L)
  scan_paf  <- .align_hp(fit, scan_conc)

  lo_idx   <- max(1L, which(scan_paf >= 1e-9)[1L] - 1L)
  hi_idx   <- min(length(scan_conc), tail(which(scan_paf <= 1 - 1e-9), 1L) + 1L)

  list(
    log10_lo = log10(scan_conc[lo_idx]),
    log10_hi = log10(scan_conc[hi_idx])
  )
}

## Build a validated monotone-cubic spline for (fit, range).
## Doubles knot count from M = 1025 until max|err| < 1e-8 or M > 16 384.
## Returns a named list(M, pg, max_err) where pg is the PAF grid vector.
.build_spline <- function(fit, log10_lo, log10_hi, analyte, method) {
  M <- 1025L
  repeat {
    lg    <- seq(log10_lo, log10_hi, length.out = M)
    pg    <- .align_hp(fit, 10^lg)

    spfun <- stats::splinefun(lg, pg, method = "monoH.FC")

    ## Validate on 8 000 random concentrations within the support range.
    ## Seed fixed so the check is reproducible across re-runs.
    set.seed(99L)
    q_check <- 10^stats::runif(8000L, log10_lo, log10_hi)
    tr      <- .align_hp(fit, q_check)

    ## Clamp spline output to [0, 1] before computing error (splines can
    ## overshoot slightly at boundaries; the clamp is applied at query time).
    max_err <- max(
      abs(pmin(pmax(spfun(log10(q_check)), 0), 1) - tr),
      na.rm = TRUE
    )

    if (max_err < 1e-8) break
    if (M >= 16384L) {
      stop(sprintf(
        "Cannot reach 1e-8 accuracy budget for: %s (%s) — max_err = %.2e",
        analyte, method, max_err
      ))
    }
    M <- M * 2L
  }

  list(M = M, pg = pg, max_err = max_err)
}

## ── main loop ────────────────────────────────────────────────────────────────

tables <- list()

for (mth in c("multi", "anzecc")) {
  meta <- leachatetools:::.load_analyte_metadata(NULL)
  sp   <- suppressMessages(
    leachatetools:::derive_ssd_params(meta, method = mth, guideline_dir = NULL)
  )

  for (k in seq_len(nrow(sp))) {
    fit     <- sp$fit[[k]]
    analyte <- sp$analyte[k]

    if (!inherits(fit, "fitdists")) {
      message(sprintf("  [SKIP] %s / %s — no fitdists", mth, analyte))
      next
    }

    key <- paste(mth, analyte, sep = "/")
    cat(sprintf("[%s / %s] building...\n", mth, analyte))

    ## Step 1: locate effective support
    rng <- .scan_range(fit)

    ## Step 2: build and validate adaptive spline
    spl <- .build_spline(fit, rng$log10_lo, rng$log10_hi, analyte, mth)

    cat(sprintf("  M=%d  max_err=%.2e\n", spl$M, spl$max_err))

    ## Step 3: store the knot grid (not a closure — self-contained, fast I/O)
    tables[[key]] <- list(
      log10_lo         = rng$log10_lo,
      log10_hi         = rng$log10_hi,
      n                = spl$M,
      paf              = spl$pg,          # length-M PAF values at knot grid
      ssdtools_version = as.character(utils::packageVersion("ssdtools")),
      generated_on     = Sys.time()
    )
  }
}

## ── serialise ────────────────────────────────────────────────────────────────

dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
out_path <- file.path("inst", "extdata", "ssd_paf_lookup.qs2")

## qs_save stores a generic R object (named list); qd_save is for qdata frames.
qs2::qs_save(tables, out_path)

cat(sprintf(
  "\nWrote %d entries to %s (%.1f KB)\n",
  length(tables), out_path, file.size(out_path) / 1024
))
