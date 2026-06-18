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
