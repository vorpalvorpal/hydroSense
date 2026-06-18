#!/usr/bin/env Rscript
# Benchmark the #50 gap-uncertainty bracket overhead.
#
# The informative envelope reuses the ignorable simulation-smoother draws,
# freezing the residual at its posterior mean on in-gap days. The plan claims it
# is "~free given the ignorable draws": the marginal cost over the single
# (ignorable) envelope should be the one extra add_amspaf() pass on the
# frozen-residual synthetic frame, not a second full draw simulation. This
# benchmark records the wall-clock cost of each `gap_uncertainty` mode at a
# realistic daily horizon so regressions in that claim are visible.
# Outside R CMD check (lives under bench/).
#
# Usage: Rscript bench/bench-bracket-overhead.R [label]
suppressMessages({
  library(bench)
  devtools::load_all(".", quiet = TRUE)
})

label <- (commandArgs(trailingOnly = TRUE)[1]) %||% "after"
set.seed(50L)

make_chem <- function(site, dates, mult = 1) {
  analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
  purrr::map_dfr(dates, function(d) {
    tibble::tibble(
      sample_id = paste0(site, format(d, "%Y%m%d")), site_id = site,
      datetime = d, analyte = analytes,
      value = c(
        exp(rnorm(1, log(0.5), 0.3)) * mult, exp(rnorm(1, log(5), 0.4)) * mult,
        exp(rnorm(1, log(0.3), 0.3)) * mult, runif(1, 6.5, 8), runif(1, 1, 5),
        runif(1, 20, 60), runif(1, 4, 12), runif(1, 2, 8)
      ),
      detected = TRUE,
      units.analyte = ifelse(analyte %in% c("Cu", "Zn", "Ni"), "ug/L", NA)
    )
  })
}

dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 26L)
hydro <- tibble::tibble(
  date = seq(as.Date("2020-07-01"), by = "day", length.out = 760L),
  value = pmax(0, rnorm(760, 2, 4))
)
ref <- make_chem("reference", dates)
tgt <- make_chem("target", dates, mult = 5)
rm <- suppressMessages(fit_reference_model(
  ref,
  hydro = hydro, conc_units = "ug/L", min_obs_model = 10L,
  api_tau_bounds_short = c(7, 7), api_tau_bounds_long = c(30, 30)
))

run <- function(mode) {
  suppressWarnings(suppressMessages(amspaf_daily(
    tgt,
    reference_model = rm, interpolation = "model",
    require_temperature = FALSE, conc_units = "ug/L",
    gap_uncertainty = mode, ndraws = 100L, seed = 1L, return = "summary"
  )))
}

results <- bench::mark(
  ignorable = run("ignorable"),
  informative = run("informative"),
  bracket = run("bracket"),
  iterations = 5, check = FALSE
)

print(results[, c("expression", "median", "mem_alloc")])
saveRDS(results, sprintf("bench/results-bracket-overhead-%s.rds", label))
