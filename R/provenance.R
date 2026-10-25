#' Record what computed a derived product
#'
#' LMF and msPAF are derived products: they depend on this package's code as
#' much as on the data. `hydrosense_provenance()` returns a small named list
#' identifying the code and the time of a computation, for the caller to store
#' beside the derived values (hydroSense itself writes nothing anywhere).
#'
#' `git_sha` comes from the `RemoteSha` field recorded when the package was
#' installed from GitHub. When the package is loaded from a source checkout
#' (for example with `devtools::load_all()`) it is the checkout's `HEAD`,
#' read with `git rev-parse`. Otherwise (a CRAN or local-tarball install, or no
#' usable `git`) it is `NA`.
#'
#' @return A named list:
#'   \describe{
#'     \item{package_version}{`character(1)`, the installed hydroSense version.}
#'     \item{git_sha}{`character(1)` commit SHA, or `NA_character_`.}
#'     \item{r_version}{`character(1)`, the R version.}
#'     \item{computed_at}{`POSIXct(1)` in UTC, the time of the call.}
#'   }
#'
#' @examples
#' hydrosense_provenance()
#' @export
hydrosense_provenance <- function() {
  list(
    package_version = as.character(utils::packageVersion("hydroSense")),
    git_sha         = .hydrosense_git_sha(),
    r_version       = as.character(getRversion()),
    computed_at     = as.POSIXct(Sys.time(), tz = "UTC")
  )
}

## The commit hydroSense was installed from or is loaded from; NA if unknown.
.hydrosense_git_sha <- function() {
  is_sha <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) &&
      grepl("^[0-9a-f]{7,40}$", x)
  }

  sha <- tryCatch(
    utils::packageDescription("hydroSense", fields = "RemoteSha"),
    error = function(e) NA_character_,
    warning = function(w) NA_character_
  )
  if (is_sha(sha)) return(sha)

  path <- system.file(package = "hydroSense")
  ## An installed package keeps its code in R/*.rdb; a source tree (as loaded
  ## by devtools::load_all()) has the .R files.
  from_source <- nzchar(path) &&
    length(list.files(file.path(path, "R"), pattern = "\\.R$")) > 0L
  if (!from_source) return(NA_character_)

  sha <- tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(path), "rev-parse", "HEAD"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) NA_character_
  )
  if (length(sha) == 1L && is_sha(sha)) sha else NA_character_
}
