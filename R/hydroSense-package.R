#' @keywords internal
"_PACKAGE"

# All dplyr/tidyr/tibble calls across the package are namespace-qualified
# (dplyr::, tidyr::, tibble::), so those namespaces are not imported here -- the
# packages only need to be listed under Imports in DESCRIPTION. The rlang
# data-mask pronouns .data / .env are used unqualified inside dplyr verbs, so
# they are imported. The `:=` walrus (dynamic column names in dplyr::summarise,
# e.g. .summarise_bracket()) is likewise used unqualified and imported here.
#' @importFrom rlang .data .env :=
NULL

# testthat can only mock a binding that exists in the package namespace. brms
# is an optional dependency (Suggests), and the specs that exercise its absence
# mock requireNamespace() through this placeholder, as documented in
# ?testthat::local_mocked_bindings. A non-function binding is skipped when R
# looks up a function call, so at run time requireNamespace() still resolves to
# base::requireNamespace().
requireNamespace <- NULL
