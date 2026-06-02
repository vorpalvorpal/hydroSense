#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#'
#' @import dplyr
#' @importFrom tibble tibble
#'
#' The leachate-mixing-fraction code in `lmf.R` is written against unqualified
#' dplyr verbs (`filter()`, `mutate()`, `group_by()`, ...) and `tibble()`, so
#' those generics are imported into the package namespace here. The rest of the
#' package qualifies its dplyr calls explicitly (`dplyr::...`); importing the
#' namespace is harmless for that code and makes the `lmf.R` verbs resolve to
#' dplyr rather than falling through to `stats::filter()` etc.
#'
## usethis namespace: end
NULL
