## ============================================================================
## Analyte-name aliases: sampleTidy's stored names are recognised
## ============================================================================
##
## Contract (TDD):
##   * analyte_aliases() is exported and returns a data frame with character
##     columns `alias` and `canonical`, one row per alias (no duplicate
##     aliases). It covers at least sampleTidy's stored names for the LMF ions
##     and the msPAF drivers listed in .st_names below.
##   * Where the package already has a name for the quantity, `canonical` is
##     that name (identity for names that already match). For the two
##     sampleTidy quantities with no package name today
##     (Alkalinity-total-CaCO3, NO2+NO3-N) the canonical value is left to the
##     implementation, but must be a non-empty string.
##   * add_lmf() and add_mspaf() give the same results for sampleTidy-named
##     input as for package-named input, and return the caller's rows with
##     their original analyte names.
##   * Alkalinity-total-CaCO3 stands in for CO3-CaCO3 + HCO3-CaCO3 (it is
##     their sum) and NO2+NO3-N for NO3-N + NO2-N; neither is double counted
##     when the split species are also present.

library(testthat)
library(hydroSense)

.st_names <- c(
  "Na", "K", "Ca", "Mg", "Cl", "SO4²⁻", "F", "NH3-N", "NO3-N",
  "NO2-N", "NO2+NO3-N", "Alkalinity-total-CaCO3", "pH", "DOC",
  "Hardness-total-CaCO3", "Temperature"
)

.known_canonical <- c(
  "Na" = "Na", "K" = "K", "Ca" = "Ca", "Mg" = "Mg", "Cl" = "Cl",
  "SO4²⁻" = "SO4²⁻", "F" = "F", "NH3-N" = "NH3-N",
  "NO3-N" = "NO3-N", "NO2-N" = "NO2-N", "pH" = "pH", "DOC" = "DOC",
  "Hardness-total-CaCO3" = "hardness", "Temperature" = "temperature"
)


## ---------------------------------------------------------------------------
## analyte_aliases()
## ---------------------------------------------------------------------------

test_that("analyte_aliases() is exported", {
  expect_true("analyte_aliases" %in% getNamespaceExports("hydroSense"))
})

test_that("analyte_aliases() returns an alias -> canonical data frame", {
  al <- analyte_aliases()
  expect_s3_class(al, "data.frame")
  expect_true(all(c("alias", "canonical") %in% names(al)))
  expect_type(al$alias, "character")
  expect_type(al$canonical, "character")
  expect_false(anyNA(al$alias))
  expect_false(anyNA(al$canonical))
  expect_false(anyDuplicated(al$alias) > 0)
})

test_that("analyte_aliases() covers sampleTidy's names with the right canonical names", {
  al <- analyte_aliases()
  expect_true(all(.st_names %in% al$alias),
              info = paste("missing:", paste(setdiff(.st_names, al$alias), collapse = ", ")))
  got <- al$canonical[match(names(.known_canonical), al$alias)]
  expect_identical(got, unname(.known_canonical))

  for (nm in c("Alkalinity-total-CaCO3", "NO2+NO3-N")) {
    can <- al$canonical[al$alias == nm]
    expect_length(can, 1L)
    expect_true(nzchar(can))
  }
})

test_that("the SO4 alias is the Unicode-superscript spelling sampleTidy stores", {
  al <- analyte_aliases()
  so4 <- "SO4²⁻"
  expect_identical(utf8ToInt(so4), c(83L, 79L, 52L, 178L, 8315L))
  expect_true(so4 %in% al$alias)
})


## ---------------------------------------------------------------------------
## add_lmf(): combined alkalinity and combined oxidised N
## ---------------------------------------------------------------------------

.al_em <- function() {
  set.seed(101)
  list(
    ref = .pi_samples(rep("2022-06-01", 12), 0, "ref", "refsite"),
    lea = .pi_samples(rep("2022-06-01", 12), 1, "lea", "leasite")
  )
}

.al_query <- function() {
  fvec <- c(0.2, 0.5, 0.8)
  set.seed(7)
  dplyr::bind_rows(lapply(seq_along(fvec), function(i) {
    .pi_samples("2022-06-01", fvec[i], paste0("q", i), "query")
  }))
}

## Replace two species rows per sample by one combined row (value = sum).
.al_combine <- function(df, parts, combined, keep_parts = FALSE) {
  sums <- df[df$analyte %in% parts, ] |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(value = sum(value), .groups = "drop")
  tmpl <- df[df$analyte == parts[1], ]
  tmpl$value   <- sums$value[match(tmpl$sample_id, sums$sample_id)]
  tmpl$analyte <- combined
  base <- if (keep_parts) df else df[!df$analyte %in% parts, ]
  dplyr::bind_rows(base, tmpl)
}

.al_lmf_values <- function(q, em) {
  out <- add_lmf(q, leachate_data = em$lea, reference_data = em$ref)
  .pi_lmf(out)
}

test_that("Alkalinity-total-CaCO3 gives the same LMF as the CO3 + HCO3 split", {
  em  <- .al_em()
  q   <- .al_query()
  ref <- .pi_lmf(add_lmf(q, leachate_data = em$lea, reference_data = em$ref))

  alk  <- c("CO3-CaCO3", "HCO3-CaCO3")
  em_t <- list(ref = .al_combine(em$ref, alk, "Alkalinity-total-CaCO3"),
               lea = .al_combine(em$lea, alk, "Alkalinity-total-CaCO3"))
  q_t  <- .al_combine(q, alk, "Alkalinity-total-CaCO3")

  tot <- .al_lmf_values(q_t, em_t)
  expect_equal(tot$value, ref$value, tolerance = 1e-6)
  expect_equal(tot$n_ions_used, ref$n_ions_used)
})

test_that("Alkalinity-total-CaCO3 is not double counted alongside the split", {
  em  <- .al_em()
  q   <- .al_query()
  ref <- .al_lmf_values(q, em)

  alk  <- c("CO3-CaCO3", "HCO3-CaCO3")
  em_b <- list(ref = .al_combine(em$ref, alk, "Alkalinity-total-CaCO3", keep_parts = TRUE),
               lea = .al_combine(em$lea, alk, "Alkalinity-total-CaCO3", keep_parts = TRUE))
  q_b  <- .al_combine(q, alk, "Alkalinity-total-CaCO3", keep_parts = TRUE)

  both <- .al_lmf_values(q_b, em_b)
  expect_equal(both$value, ref$value, tolerance = 1e-6)
})

test_that("NO2+NO3-N gives the same LMF as separate NO3-N and NO2-N", {
  em  <- .al_em()
  q   <- .al_query()
  ref <- .al_lmf_values(q, em)

  nox  <- c("NO3-N", "NO2-N")
  em_c <- list(ref = .al_combine(em$ref, nox, "NO2+NO3-N"),
               lea = .al_combine(em$lea, nox, "NO2+NO3-N"))
  q_c  <- .al_combine(q, nox, "NO2+NO3-N")

  comb <- .al_lmf_values(q_c, em_c)
  expect_equal(comb$value, ref$value, tolerance = 1e-6)
  expect_equal(comb$n_ions_used, ref$n_ions_used)
})

test_that("NO2+NO3-N is not double counted alongside NO3-N and NO2-N", {
  em  <- .al_em()
  q   <- .al_query()
  ref <- .al_lmf_values(q, em)

  nox  <- c("NO3-N", "NO2-N")
  em_b <- list(ref = .al_combine(em$ref, nox, "NO2+NO3-N", keep_parts = TRUE),
               lea = .al_combine(em$lea, nox, "NO2+NO3-N", keep_parts = TRUE))
  q_b  <- .al_combine(q, nox, "NO2+NO3-N", keep_parts = TRUE)

  both <- .al_lmf_values(q_b, em_b)
  expect_equal(both$value, ref$value, tolerance = 1e-6)
})

test_that("add_lmf on fully sampleTidy-shaped, sampleTidy-named input matches package input", {
  em <- .al_em()
  q  <- .al_query()
  native <- .al_lmf_values(q, em)

  to_st <- function(df) {
    df <- .al_combine(df, c("CO3-CaCO3", "HCO3-CaCO3"), "Alkalinity-total-CaCO3")
    df <- .al_combine(df, c("NO3-N", "NO2-N"), "NO2+NO3-N")
    .pi_as_sampletidy(df)
  }
  q_st <- to_st(q)
  out  <- add_lmf(q_st, leachate_data = to_st(em$lea), reference_data = to_st(em$ref))
  st   <- .pi_lmf(out)
  expect_equal(st$value, native$value, tolerance = 1e-6)

  ## Caller's analyte names are returned untouched.
  expect_setequal(unique(out$analyte[out$analyte != "LMF"]), unique(q_st$analyte))
})


## ---------------------------------------------------------------------------
## add_mspaf(): Hardness-total-CaCO3 and Temperature drive normalisation
## ---------------------------------------------------------------------------

test_that("add_mspaf gives the same msPAF for sampleTidy-named drivers", {
  tgt <- .pi_mspaf_chem("ds", n = 6L, seed = 2L, mult = 5)
  ref <- .pi_mspaf_chem("ref", n = 6L, seed = 1L)

  pkg <- add_mspaf(tgt, reference = ref)
  pkg_m <- pkg[pkg$analyte == "msPAF", ]
  expect_equal(nrow(pkg_m), 6L)   # fixture sanity: package names work today

  tgt_st <- .pi_mspaf_sampletidy_names(tgt)
  ref_st <- .pi_mspaf_sampletidy_names(ref)
  expect_true(all(c("Hardness-total-CaCO3", "Temperature") %in% tgt_st$analyte))
  expect_false(any(c("hardness", "temperature") %in% tgt_st$analyte))

  st   <- add_mspaf(tgt_st, reference = ref_st)
  st_m <- st[st$analyte == "msPAF", ]

  pkg_m <- pkg_m[order(pkg_m$sample_id), ]
  st_m  <- st_m[order(st_m$sample_id), ]
  expect_identical(st_m$sample_id, pkg_m$sample_id)
  expect_equal(st_m$value, pkg_m$value, tolerance = 1e-8)
  expect_identical(st_m$n_analytes_used, pkg_m$n_analytes_used)

  ## Caller's rows come back with their original (sampleTidy) names.
  expect_setequal(unique(st$analyte[st$analyte != "msPAF"]), unique(tgt_st$analyte))
})

test_that("add_mspaf treats Temperature as the water temperature NH3-N needs", {
  tgt <- .pi_mspaf_sampletidy_names(.pi_mspaf_chem("ds", n = 3L, seed = 2L))
  ## require_temperature = TRUE (default): NH3-N with a Temperature row is fine.
  expect_no_error(add_mspaf(tgt))
  ## ...and dropping Temperature still triggers the mandatory-temperature error.
  expect_error(add_mspaf(tgt[tgt$analyte != "Temperature", ]), "temperature")
})
