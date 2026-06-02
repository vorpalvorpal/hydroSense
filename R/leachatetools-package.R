#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#'
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_wider pivot_longer
#'
#' The leachate-mixing-fraction code in `lmf.R` is written against unqualified
#' dplyr verbs (`filter()`, `mutate()`, `group_by()`, ...), `tibble()`, and the
#' tidyr pivots, so those generics are imported into the package namespace here.
#' The rest of the package qualifies its dplyr/tidyr calls explicitly
#' (`dplyr::...`, `tidyr::...`); importing the namespace is harmless for that
#' code and makes the `lmf.R` verbs resolve correctly rather than failing to be
#' found or falling through to `stats::filter()` etc.
#'
## usethis namespace: end
NULL
