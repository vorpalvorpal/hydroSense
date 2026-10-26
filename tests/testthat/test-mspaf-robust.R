## Robustness of the msPAF engine to two defects seen on real monitoring data
## ============================================================================
##
## CONTRACT.
##
## D1. A point with no SSD-eligible analyte must not abort the site.
##   A point whose rows are all non-toxicants (e.g. only "Standing water level"
##   and "temperature") has nothing to assess.  add_mspaf() must:
##     * on a frame mixing it with an assessable point, return the assessable
##       point's msPAF row (and its analyte_pafs / ara_summary diagnostics)
##       IDENTICAL to running that point alone, and no msPAF row for the
##       toxicant-free point;
##     * on the toxicant-free point alone, return its input rows with ZERO msPAF
##       rows and the usual msPAF output columns and types (the same schema as a
##       point that falls below min_analytes), without an error;
##     * in both cases, say which point(s) had no SSD-eligible analytes.  This
##       is signalled as a cli MESSAGE (cli_inform), because the package already
##       reports dropped data that way (.summarise_mspaf_diagnostics(),
##       .summarise_ara_coverage()); cli_warn is kept for "msPAF cannot be
##       computed at all".  The message must name the point id (site_id).
##   compute_mspaf_per_sample() on such a point returns the zero-row point-mode
##   schema, i.e. .mspaf_empty_row(with_sample_id = TRUE) without draw_id.
##   One level up, mspaf_daily() and mspaf_pipeline() (chronic, impute = FALSE)
##   over two normal points plus one toxicant-free point must not error, must
##   give the normal points results IDENTICAL to a run without the
##   toxicant-free point, and must give the toxicant-free point no rows.
##
## D2. Results must not depend on the row order of the reference input.
##   fit_imputation_model() must give the same fitted object whatever order its
##   training rows arrive in (seeded shuffles, and rows sorted under en_US
##   collation): its PCA (NIPALS scores and loadings, pc_scores, medians,
##   active variables) and the imputed predictions it produces must be
##   IDENTICAL (expect_identical).  Exact equality is the bar because it is
##   attainable: once the training rows are canonicalised (keyed order) the
##   computation is the same sequence of floating-point operations, and any
##   residual difference is the very defect being tested (today the PCA
##   differs at ~1e-14 and the predictions by up to ~0.25%).
##   The full mspaf_pipeline() chronic golden case with a row-shuffled
##   reference must equal the unshuffled case to 1e-10 relative (numeric
##   columns, same NA pattern; everything else exact).  The unshuffled case is
##   represented by its golden (_goldens/mspaf_pipeline_chronic.rds), which
##   test-perf-goldens.R already pins to the unshuffled run at 1e-10; comparing
##   against the golden saves a second ~30 s pipeline run.  Goldens are never
##   regenerated.
##
## All data are synthetic (helper-perf-fixtures.R and the frames built here).

library(testthat)

# ── D1 fixtures ───────────────────────────────────────────────────────────────

.mr_date <- as.Date("2024-03-01")

## Point "MB01": three SSD toxicants in ug/L (the SSD unit), one sample.
.mr_point_tox <- function() {
  tibble::tibble(
    sample_id = "MB01_s1", site_id = "MB01", datetime = .mr_date,
    analyte = c("B", "Cr", "Mn"), value = c(500, 5, 800),
    detected = TRUE, units.analyte = "ug/L"
  )
}

## Point "MB02_WL": no SSD-eligible analyte at all.
.mr_point_free <- function() {
  tibble::tibble(
    sample_id = "MB02_WL_s1", site_id = "MB02_WL", datetime = .mr_date,
    analyte = c("Standing water level", "temperature"), value = c(3.2, 17),
    detected = TRUE, units.analyte = c("m", NA_character_)
  )
}

.mr_add <- function(df) {
  add_mspaf(df, reference = NULL, require_temperature = FALSE)
}
.mr_mspaf_rows <- function(x) x[x$analyte == "msPAF", , drop = FALSE]

## The caller's own input columns of an add_mspaf() result, without the
## diagnostic attributes that `[` carries along.
.mr_input_cols <- function(x, rows = TRUE) {
  x <- x[rows, names(.mr_point_free()), drop = FALSE]
  attr(x, "analyte_pafs") <- NULL
  attr(x, "ara_summary")  <- NULL
  x
}

## The msPAF output columns and their types (point mode: no draw_id).
.mr_mspaf_schema <- function() {
  s <- hydroSense:::.mspaf_empty_row(with_sample_id = FALSE)
  s <- s[, setdiff(names(s), "draw_id"), drop = FALSE]
  vapply(s, function(x) class(x)[1L], character(1))
}


describe("D1: add_mspaf() with a point that has no SSD-eligible analyte", {

  it("returns the assessable point's row unchanged and none for the free point", {
    alone <- suppressMessages(.mr_add(.mr_point_tox()))
    both  <- suppressMessages(.mr_add(dplyr::bind_rows(.mr_point_tox(),
                                                        .mr_point_free())))
    expect_identical(nrow(.mr_mspaf_rows(alone)), 1L)
    expect_identical(.mr_mspaf_rows(both), .mr_mspaf_rows(alone))
    expect_false("MB02_WL" %in% .mr_mspaf_rows(both)$site_id)
    expect_identical(attr(both, "analyte_pafs"), attr(alone, "analyte_pafs"))
    expect_identical(attr(both, "ara_summary"), attr(alone, "ara_summary"))
    # The free point's own input rows are handed back untouched.
    expect_identical(.mr_input_cols(both, both$site_id == "MB02_WL"),
                     .mr_point_free())
  })

  it("names the toxicant-free point in a message (mixed frame)", {
    expect_message(
      .mr_add(dplyr::bind_rows(.mr_point_tox(), .mr_point_free())),
      "MB02_WL", fixed = TRUE
    )
  })

  it("returns zero msPAF rows with the usual schema for the free point alone", {
    out <- suppressMessages(.mr_add(.mr_point_free()))
    expect_identical(nrow(.mr_mspaf_rows(out)), 0L)
    # Input rows come back as supplied.
    expect_identical(.mr_input_cols(out), .mr_point_free())
    # The msPAF output columns are present with their normal types: the same
    # schema add_mspaf() returns for a point that falls below min_analytes.
    schema <- .mr_mspaf_schema()
    expect_true(all(names(schema) %in% names(out)))
    got <- vapply(out[names(schema)], function(x) class(x)[1L], character(1))
    expect_identical(got, schema)
    below_min <- suppressMessages(.mr_add(.mr_point_tox()[1L, ]))
    expect_identical(names(out), names(below_min))
    expect_identical(lapply(out, class), lapply(below_min, class))
  })

  it("names the toxicant-free point in a message (free point alone)", {
    expect_message(.mr_add(.mr_point_free()), "MB02_WL", fixed = TRUE)
  })

  it("compute_mspaf_per_sample() returns the zero-row point-mode schema", {
    ssd <- hydroSense:::derive_ssd_params(
      hydroSense:::.load_analyte_metadata(NULL), "multi",
      getOption("hydroSense.guideline_dir")
    )
    out <- hydroSense:::compute_mspaf_per_sample(
      sample_data   = .mr_point_free(),
      ref_table     = tibble::tibble(analyte = character(0),
                                     ref_norm = numeric(0)),
      ssd_params    = ssd,
      min_analytes  = 3L,
      method        = "multi",
      guideline_dir = getOption("hydroSense.guideline_dir"),
      ara_enabled   = FALSE
    )
    empty <- hydroSense:::.mspaf_empty_row(with_sample_id = TRUE)
    empty <- empty[, setdiff(names(empty), "draw_id"), drop = FALSE]
    expect_identical(nrow(out), 0L)
    expect_identical(names(out), names(empty))
    expect_identical(lapply(out, class), lapply(empty, class))
  })
})


# ── D1 one level up: mspaf_daily() and mspaf_pipeline() ───────────────────────

## Two normal points (P1, P2 of the small structural target) ...
.mr_normal_target <- function() {
  tgt <- .pf_small_target()
  tgt[tgt$site_id %in% c("P1", "P2"), , drop = FALSE]
}

## ... plus a toxicant-free point PX sampled monthly through the same period.
.mr_free_target <- function() {
  d   <- seq(as.Date("2022-01-10"), by = "30 days", length.out = 10L)
  sid <- sprintf("PX_%03d", seq_along(d))
  dplyr::bind_rows(
    tibble::tibble(sample_id = sid, site_id = "PX", datetime = d,
                   analyte = "Standing water level",
                   value = 3 + sin(seq_along(d)), detected = TRUE,
                   units.analyte = "m"),
    tibble::tibble(sample_id = sid, site_id = "PX", datetime = d,
                   analyte = "temperature",
                   value = 15 + seq_along(d) / 3, detected = TRUE,
                   units.analyte = NA_character_)
  )
}

## External temperature overhanging a two-month daily grid (kept short so the
## four runs below stay within a few seconds each).
.mr_temperature <- function() {
  tmp <- .pf_scenario()$temperature
  tmp[tmp$datetime >= as.Date("2021-12-15") &
        tmp$datetime <= as.Date("2022-03-10"), ]
}

.mr_daily <- function(target) {
  sc <- .pf_scenario()
  suppressWarnings(suppressMessages(mspaf_daily(
    target, temperature = .mr_temperature(), reference = sc$reference,
    reference_model = .pf_reference_model(), interpolation = "model",
    start = as.Date("2022-01-01"), end = as.Date("2022-02-28")
  )))
}

.mr_pipeline <- function(target) {
  sc <- .pf_scenario()
  suppressWarnings(suppressMessages(mspaf_pipeline(
    target, reference = sc$reference, hydro = sc$rain, impute = FALSE,
    reference_args = list(
      min_obs_model = 20L,
      api_tau_bounds_short = c(3, 10),
      api_tau_bounds_long  = c(30, 90)
    ),
    daily_args = list(
      reference = sc$reference, temperature = .mr_temperature(),
      start = as.Date("2022-01-01"), end = as.Date("2022-02-28")
    ),
    chronic = TRUE,
    focal_dates = as.Date(c("2022-01-20", "2022-02-10", "2022-02-28")),
    tau = 90, window = 365
  )))
}

## Chronic frame without the (model-object) attributes.
.mr_strip_pipeline <- function(x) {
  attributes(x)[c("daily", "reference_model", "imputation_model")] <- NULL
  x
}


describe("D1: a toxicant-free point inside mspaf_daily() / mspaf_pipeline()", {

  it("mspaf_daily(): normal points unchanged, no rows for the free point", {
    base <- .mr_daily(.mr_normal_target())
    expect_setequal(unique(base$site_id), c("P1", "P2"))
    with_free <- .mr_daily(dplyr::bind_rows(.mr_normal_target(),
                                            .mr_free_target()))
    expect_false("PX" %in% with_free$site_id)
    expect_false("PX" %in% attr(with_free, "analyte_pafs")$site_id)
    # Identical frame, including the analyte_pafs / ara_summary attributes.
    expect_identical(with_free, base)
  })

  it("mspaf_pipeline() chronic: normal points unchanged, no rows for the free point", {
    base <- .mr_pipeline(.mr_normal_target())
    expect_setequal(unique(base$site_id), c("P1", "P2"))
    with_free <- .mr_pipeline(dplyr::bind_rows(.mr_normal_target(),
                                               .mr_free_target()))
    expect_false("PX" %in% with_free$site_id)
    expect_false("PX" %in% attr(with_free, "daily")$site_id)
    expect_identical(.mr_strip_pipeline(with_free), .mr_strip_pipeline(base))
    expect_identical(attr(with_free, "daily"), attr(base, "daily"))
  })
})


# ── D2: row-order invariance of the imputation model ──────────────────────────

.mr_fit_imp <- function(ref) {
  suppressWarnings(suppressMessages(fit_imputation_model(ref)))
}
.mr_predict <- function(model) {
  suppressWarnings(suppressMessages(
    impute_chemistry(.pf_scenario()$target, model)
  ))
}

## Baseline fit on the reference rows in their fixture order (memoised).
.mr_imp_base <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      m <- .mr_fit_imp(.pf_scenario()$reference)
      cache <<- list(model = m, pred = .mr_predict(m))
    }
    cache
  }
})

.mr_shuffle <- function(df, seed) {
  df[withr::with_seed(seed, sample.int(nrow(df))), , drop = FALSE]
}

## Fitted object and predictions must be identical to the baseline.
.mr_expect_same_imp <- function(model, pred) {
  base <- .mr_imp_base()
  expect_identical(model$pca$fit$scores,   base$model$pca$fit$scores)
  expect_identical(model$pca$fit$loadings, base$model$pca$fit$loadings)
  expect_identical(model$pca$pc_scores,    base$model$pca$pc_scores)
  expect_identical(model$pca,              base$model$pca)
  expect_identical(pred,                   base$pred)
}


describe("D2: fit_imputation_model() is invariant to training-row order", {

  # Seeds 1-3 each change the fitted PCA and the predictions on today's code
  # (predictions by up to 0.02%, 0.08% and 0.25%).
  for (seed in 1:3) {
    local({
      seed <- seed
      it(paste("gives an identical fit on shuffled rows, seed", seed), {
        ref <- .mr_shuffle(.pf_scenario()$reference, seed)
        m <- .mr_fit_imp(ref)
        .mr_expect_same_imp(m, .mr_predict(m))
      })
    })
  }

  it("gives an identical fit on rows sorted (and fitted) under en_US collation", {
    old <- Sys.getlocale("LC_COLLATE")
    withr::defer(Sys.setlocale("LC_COLLATE", old))
    ok <- suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))
    skip_if_not(nzchar(ok), "en_US.UTF-8 collation is not available")
    ref <- .pf_scenario()$reference
    # en_US interleaves case ("Ca", "Cu", "DOC", ..., "pH", "temperature",
    # "Zn"), unlike C collation, so analytes first appear in another order.
    ref <- ref[order(ref$sample_id, ref$analyte), , drop = FALSE]
    m <- .mr_fit_imp(ref)
    p <- .mr_predict(m)
    Sys.setlocale("LC_COLLATE", old)
    .mr_expect_same_imp(m, p)
  })
})


# ── D2: the chronic pipeline golden case on a shuffled reference ──────────────

## Element-wise relative comparison: 1e-10 relative, 1e-12 absolute floor near
## zero, identical NA pattern (the test-perf-goldens.R contract).
.mr_expect_num <- function(actual, golden, label) {
  expect_identical(length(actual), length(golden), label = label)
  expect_identical(is.na(actual), is.na(golden),
                   label = paste(label, "NA pattern"))
  ok <- !is.na(golden) & !is.na(actual)
  diff  <- abs(actual[ok] - golden[ok])
  scale <- pmax(abs(golden[ok]), abs(actual[ok]))
  bad <- diff > 1e-10 * scale & diff > 1e-12
  expect(
    !any(bad),
    sprintf("%s: %d value(s) differ by more than 1e-10 relative; worst relative difference %.3g",
            label, sum(bad), max(diff[bad] / scale[bad]))
  )
}

.mr_expect_frame <- function(actual, golden, label) {
  expect_identical(class(actual), class(golden), label = paste(label, "class"))
  expect_identical(names(actual), names(golden), label = paste(label, "names"))
  expect_identical(nrow(actual), nrow(golden), label = paste(label, "nrow"))
  if (!identical(names(actual), names(golden)) ||
      !identical(nrow(actual), nrow(golden))) {
    return(invisible())
  }
  for (nm in names(golden)) {
    a <- actual[[nm]]
    g <- golden[[nm]]
    lab <- paste0(label, "$", nm)
    if (is.double(g) && !inherits(g, c("Date", "POSIXt"))) {
      .mr_expect_num(as.numeric(a), as.numeric(g), lab)
    } else {
      expect_identical(a, g, label = lab)
    }
  }
}


describe("D2: mspaf_pipeline() chronic golden case on a shuffled reference", {
  .pf_skip_unless_golden_platform()

  it("matches the unshuffled (golden) result to 1e-10 relative", {
    golden <- readRDS(.pf_golden_path("mspaf_pipeline_chronic"))
    attr(golden, "golden_meta") <- NULL
    # Seed 1 moves daily msPAF by ~1.2% and chronic msPAF by ~0.3% today.
    ref <- .mr_shuffle(.pf_scenario()$reference, 1L)
    actual <- .pf_pipeline_case(reference = ref)
    expect_identical(names(actual), names(golden))
    for (part in names(golden)) {
      .mr_expect_frame(actual[[part]], golden[[part]], part)
    }
  })
})
