#' Synthetic leachate-impacted water-quality data (demo)
#'
#' Returns a small, fictional long-format water-quality dataset used in the
#' examples and vignettes. It is **not** real monitoring data — concentrations
#' are generated from a simple two-component mixing model so the worked
#' examples are reproducible and self-contained.
#'
#' The data are bundled in the qs2 *qdata* format under `inst/extdata` (qs2 is
#' already a dependency) and read on demand by this accessor, rather than being
#' lazy-loaded R data.
#'
#' Three sites are included:
#' \describe{
#'   \item{`downstream`}{A leachate-impacted site, ~15% leachate mixing
#'     fraction relative to the reference.}
#'   \item{`reference`}{Clean background / upstream chemistry.}
#'   \item{`leachate`}{The leachate end-member.}
#' }
#'
#' Each site has six bi-monthly samples through 2024. The analyte panel carries
#' the major ions used by [add_lmf()] / [to_meq()] (`Na`, `K`, `Ca`, `Mg`,
#' `Cl`, sulfate, `F`, `NO3-N`, `NO2-N`, `CO3-CaCO3`, `HCO3-CaCO3`), the
#' toxicants assessed by [ssd_paf()] / [add_amspaf()] (`Cu`, `Zn`, `NH3-N`),
#' and the co-analytes needed for the ammonia and bioavailability
#' normalisations (`pH`, `temperature`, `DOC`, `hardness`).
#'
#' @return A [tibble][tibble::tibble] with one row per sample x analyte and the
#'   columns:
#' \describe{
#'   \item{sample_id}{Character. Unique sample identifier.}
#'   \item{site_id}{Character. One of `"downstream"`, `"reference"`,
#'     `"leachate"`.}
#'   \item{datetime}{Date. Sampling date.}
#'   \item{analyte}{Character. Analyte name.}
#'   \item{value}{Numeric. Concentration in `units.analyte`.}
#'   \item{detected}{Logical. `FALSE` marks a below-detection-limit result.}
#'   \item{units.analyte}{Character. Concentration units: major ions in
#'     `"mg/L"`; `Cu`/`Zn`/`NH3-N` in `"ug/L"` (the SSD scale); `pH` unitless;
#'     `temperature` in `"degC"`; `hardness` in `"mg/L CaCO3"`.}
#'   \item{valence.analyte}{Numeric. Ionic charge (for the meq conversion in
#'     [to_meq()]); `NA` for non-ionic analytes.}
#'   \item{atomic_mass.analyte}{Numeric. Molar/atomic mass in g/mol (for the
#'     meq conversion); `NA` for non-ionic analytes.}
#' }
#'
#' @examples
#' demo <- leachate_demo()
#' table(demo$site_id, demo$analyte == "Cu")
#' head(subset(demo, site_id == "downstream"), 4)
#'
#' @export
leachate_demo <- function() {
  path <- system.file("extdata", "leachate_demo.qs2",
                      package = "leachatetools")
  if (!nzchar(path))
    cli::cli_abort("Bundled {.file leachate_demo.qs2} not found.")
  qs2::qd_read(path)
}
