## ============================================================================
## Shared builders for the pipeline-inputs TDD specs
## ============================================================================
##
## Used by test-lmf-trailing.R, test-lmf-explicit-inputs.R,
## test-ion-chemistry.R, test-analyte-aliases.R and
## test-mspaf-pipeline-plain.R.
##
## The LMF builders reuse the two-component mixing construction of
## helper-lmf-review.R (.rev_ions()), but give every sample its own datetime so
## that the trailing calibration window can be exercised. Datetimes are POSIXct
## UTC at midnight, as sampleTidy delivers them.

## POSIXct UTC at midnight for a character/Date vector of days.
.pi_utc <- function(x) as.POSIXct(as.character(as.Date(x)), tz = "UTC")

## Long-format LMF-panel samples, one per element of `dates`, at mixing
## fraction `f` (0 = reference, 1 = leachate). `scale` multiplies every
## concentration (used to make "drifted" end-member data that would visibly
## move an LMF if it leaked into a calibration window). Package-native format:
## carries units.analyte, valence.analyte and atomic_mass.analyte.
.pi_samples <- function(dates, f, prefix, site, noise_sd = 0.02, scale = 1,
                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ions  <- .rev_ions()
  dates <- .pi_utc(dates)
  dplyr::bind_rows(lapply(seq_along(dates), function(i) {
    noise <- if (noise_sd > 0) stats::rnorm(nrow(ions), 0, noise_sd) else 0
    conc  <- ((1 - f) * ions$R + f * ions$L) * scale * (1 + noise)
    tibble::tibble(
      sample_id           = paste0(prefix, "_", format(dates[i], "%Y%m%d"), "_", i),
      site_id             = site,
      analyte             = ions$analyte,
      value               = pmax(conc, 1e-6),
      detected            = TRUE,
      datetime            = dates[i],
      units.analyte       = "mg/L",
      valence.analyte     = ions$valence,
      atomic_mass.analyte = ions$mass
    )
  }))
}

## Monthly date sequence (first of month), inclusive.
.pi_monthly <- function(from, to) {
  seq(as.Date(from), as.Date(to), by = "month")
}

## Convert package-native long data to sampleTidy's shape: only sample_id,
## site_id, analyte, value, detected, datetime and a character `units` column.
## No valence / atomic-mass columns.
.pi_as_sampletidy <- function(df) {
  tibble::tibble(
    sample_id = df$sample_id,
    site_id   = df$site_id,
    analyte   = df$analyte,
    value     = df$value,
    detected  = df$detected,
    datetime  = df$datetime,
    units     = as.character(df$units.analyte)
  )
}

## The LMF output rows of an add_lmf() result, ordered by sample_id.
.pi_lmf <- function(out) {
  lmf <- out[out$analyte == "LMF", , drop = FALSE]
  lmf[order(lmf$sample_id), , drop = FALSE]
}

## The LMF output row for one sample_id.
.pi_lmf_for <- function(out, sid) {
  lmf <- out[out$analyte == "LMF" & out$sample_id == sid, , drop = FALSE]
  testthat::expect_equal(nrow(lmf), 1L)
  lmf
}

## Stub every old dashboard global in the global environment so that, if
## add_lmf() ever falls back to them, it fails with a recognisable sentinel
## rather than "could not find function". Removed on exit of the caller.
.pi_stub_dashboard_globals <- function(env = parent.frame()) {
  nms <- c("feature_df", "data_df", "feature_sfc", "get_reference_site")
  had <- vapply(nms, exists, logical(1), envir = globalenv(), inherits = FALSE)
  old <- lapply(nms[had], get, envir = globalenv())
  names(old) <- nms[had]
  for (nm in nms) {
    assign(nm, function(...) stop("SENTINEL: dashboard global called"),
           envir = globalenv())
  }
  withr::defer({
    suppressWarnings(rm(list = nms, envir = globalenv()))
    for (nm in names(old)) assign(nm, old[[nm]], envir = globalenv())
  }, envir = env)
  invisible(nms)
}

## ---------------------------------------------------------------------------
## msPAF builders
## ---------------------------------------------------------------------------

## Long-format msPAF chemistry for `sites`, n samples each, every 30 days from
## `start`. Package-native names: toxicants Cu, Zn, Ni, NH3-N (ug/L, via
## units.analyte); drivers pH, EC, DOC, Ca, Mg, hardness, temperature
## (units.analyte NA: co-analytes are used in their natural units).
.pi_mspaf_chem <- function(sites, n = 6L, seed = 1L, mult = 1,
                           start = "2024-01-01") {
  set.seed(seed)
  an    <- c("Cu", "Zn", "Ni", "NH3-N",
             "pH", "EC", "DOC", "Ca", "Mg", "hardness", "temperature")
  dates <- .pi_utc(as.Date(start) + (seq_len(n) - 1L) * 30L)
  dplyr::bind_rows(lapply(sites, function(site) {
    g <- tidyr::expand_grid(i = seq_len(n), analyte = an)
    k <- nrow(g)
    g$sample_id <- paste0(site, "_", g$i)
    g$site_id   <- site
    g$datetime  <- dates[g$i]
    g$value <- dplyr::case_when(
      g$analyte == "pH"          ~ stats::runif(k, 6.8, 7.6),
      g$analyte == "EC"          ~ stats::runif(k, 150, 400),
      g$analyte == "DOC"         ~ stats::runif(k, 1, 4),
      g$analyte == "Ca"          ~ stats::runif(k, 5, 10),
      g$analyte == "Mg"          ~ stats::runif(k, 3, 6),
      g$analyte == "hardness"    ~ stats::runif(k, 30, 60),
      g$analyte == "temperature" ~ stats::runif(k, 12, 20),
      g$analyte == "NH3-N"       ~ mult * exp(stats::rnorm(k, log(500), 0.4)),
      TRUE                       ~ mult * exp(stats::rnorm(k, log(3), 0.4))
    )
    g$detected      <- TRUE
    g$units.analyte <- ifelse(g$analyte %in% c("Cu", "Zn", "Ni", "NH3-N"),
                              "ug/L", NA_character_)
    g$i <- NULL
    g
  }))
}

## Rename package-native msPAF driver names to sampleTidy's stored names.
.pi_mspaf_sampletidy_names <- function(df) {
  map <- c(hardness = "Hardness-total-CaCO3", temperature = "Temperature")
  hit <- df$analyte %in% names(map)
  df$analyte[hit] <- unname(map[df$analyte[hit]])
  df
}

## Daily rainfall series covering `from`..`to` (fit_reference_model() needs a
## hydrology series or SILO coordinates).
.pi_hydro <- function(from = "2023-10-01", to = "2024-12-31", seed = 3L) {
  set.seed(seed)
  d <- seq(as.Date(from), as.Date(to), by = "day")
  tibble::tibble(date = d, value = pmax(0, stats::rnorm(length(d), 2, 4)))
}
