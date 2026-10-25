## ============================================================================
## brms is optional (Suggests), and the default path never touches it
## ============================================================================
##
## Contract (TDD):
##   * DESCRIPTION lists brms under Suggests, not Imports.
##   * Asking for a brms-based imputation method ("rescor_mi", "cens",
##     "cens_factor") when brms is unavailable errors with a message that
##     mentions brms.
##   * The default impute_method ("marginal") neither calls the brms guard
##     nor loads the brms namespace.
##
## Mocking requireNamespace(): testthat cannot mock base functions directly.
## Its documented workaround is a package-level binding — the package defines
## `requireNamespace <- NULL` at top level (R skips non-function bindings
## when looking up a function call, so this changes nothing at run time) and
## the test then mocks that binding with local_mocked_bindings().

library(testthat)
library(hydroSense)

.desc_deps <- function(field) {
  path <- system.file("DESCRIPTION", package = "hydroSense")
  dcf  <- read.dcf(path, fields = field)[1, 1]
  if (is.na(dcf)) return(character(0))
  deps <- trimws(strsplit(dcf, ",")[[1]])
  deps <- sub("\\s*\\(.*\\)\\s*$", "", deps)
  deps[nzchar(deps)]
}

## requireNamespace() stand-in that reports brms as missing and defers
## everything else to base.
.no_brms <- function(package, ..., quietly = TRUE) {
  if (identical(package, "brms")) return(FALSE)
  base::requireNamespace(package, ..., quietly = quietly)
}

test_that("brms is in Suggests, not Imports", {
  expect_false("brms" %in% .desc_deps("Imports"))
  expect_true("brms" %in% .desc_deps("Suggests"))
})

test_that("NAMESPACE imports nothing from brms", {
  imports <- getNamespaceImports(asNamespace("hydroSense"))
  expect_false("brms" %in% names(imports))
})

for (meth in c("rescor_mi", "cens", "cens_factor")) {
  test_that(paste0("impute_method = '", meth, "' errors mentioning brms when brms is unavailable"), {
    local_mocked_bindings(requireNamespace = .no_brms, .package = "hydroSense")
    df <- .imp_chem(n = 30, seed = 1)
    expect_error(
      fit_imputation_model(df, impute_method = meth),
      "brms"
    )
  })
}

test_that("the brms guard itself names brms when brms is unavailable", {
  local_mocked_bindings(requireNamespace = .no_brms, .package = "hydroSense")
  expect_error(hydroSense:::.require_brms(), "brms")
})

test_that("the marginal default runs with brms unavailable", {
  local_mocked_bindings(requireNamespace = .no_brms, .package = "hydroSense")
  df <- .imp_chem(n = 40, seed = 1)
  model <- expect_no_error(suppressMessages(fit_imputation_model(
    df,
    groups = list(impute_group("metals", targets = c("Cu", "Zn"),
                               hurdle = c("Cu", "Zn")))
  )))
  expect_identical(model$impute_method, "marginal")
})

test_that("loading hydroSense and fitting the marginal default does not load brms", {
  skip_on_cran()
  skip_if_not_installed("callr")

  df   <- .imp_chem(n = 40, seed = 1)
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  from_source <- file.exists(file.path(root, "DESCRIPTION")) &&
    requireNamespace("pkgload", quietly = TRUE)

  loaded <- callr::r(
    function(df, root, from_source) {
      if (from_source) {
        suppressMessages(pkgload::load_all(root, quiet = TRUE, helpers = FALSE,
                                           export_all = FALSE))
      } else {
        suppressMessages(library(hydroSense))
      }
      before <- isNamespaceLoaded("brms")
      suppressMessages(hydroSense::fit_imputation_model(
        df,
        groups = list(hydroSense::impute_group("metals", targets = c("Cu", "Zn"),
                                               hurdle = c("Cu", "Zn")))
      ))
      c(on_load = before, after_fit = isNamespaceLoaded("brms"))
    },
    args = list(df = df, root = root, from_source = from_source)
  )
  expect_false(loaded[["on_load"]])
  expect_false(loaded[["after_fit"]])
})
