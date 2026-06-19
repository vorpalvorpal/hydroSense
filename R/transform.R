## Variance-stabilising transform for the daily impact residual (issue #15).
##
## The ARA impact I = C_norm - ref_norm is smoothed in a transformed space so
## that the latent process variance is *proportional* to concentration rather
## than constant. A single homoscedastic additive variance over-disperses the
## baseline: estimated from anchor residuals spanning baseline-to-event, it
## applies the same absolute spread at trace levels as at peaks, which the
## convex SSD then magnifies into a spuriously high baseline msPAF (issue #39,
## diagnosed on B.S01). Environmental concentrations carry multiplicative
## (lognormal) error (Helsel & Hirsch, Statistical Methods in Water Resources,
## USGS TWRI 4-A3), so a variance-stabilising transform is the principled fix.
##
## Transform: g(I) = asinh(I / c), inverse g_inv(g) = c * sinh(g).
##   * asinh (not log(I + c)) because the ARA impact is SIGNED -- C < ref is a
##     real precipitation/dilution case -- and asinh is defined on all reals,
##     odd, and smooth through zero.
##   * slope g'(I) = 1 / sqrt(I^2 + c^2): for |I| << c the transform is additive
##     (g ~ I/c), for |I| >> c it is logarithmic (g ~ sign(I) log(2|I|/c)).
##   * c sets the additive->proportional crossover. The fix does not come from
##     the per-point slope alone but from re-estimating the process variance in
##     this compressed g-space: event spikes no longer inflate the baseline
##     variance, so baseline draws stay bounded near the c scale.
##
## ref enters ONLY via the ARA difference I = C_norm - ref_norm, exactly as in
## the additive model; it is absent from the transform itself (the scale c is a
## toxicological constant from the SSD, not a property of the reference site).
##
## Shared seam: these pure helpers are also the intended transform for the
## chronic OU path (#23); keep them dependency-free.
##
## Argument naming: the impact is `I` and the scale is `c` in the maths/plan;
## the formals are `impact` / `scale_c` to keep snake_case and to avoid
## shadowing `base::c()`.

#' Variance-stabilising transform of the ARA impact
#'
#' Maps the signed impact `impact` (`I` = `C_norm - ref_norm`, normalised
#' concentration units) onto the `asinh(I / c)` scale, on which the daily
#' residual smoother's process variance is proportional to concentration rather
#' than constant. See [.g_inverse()] for the inverse and [.analyte_c()] for the
#' per-analyte scale `c`.
#'
#' @param impact Numeric vector of impacts `C_norm - ref_norm` (may be negative,
#'   `NA`, or `Inf`).
#' @param scale_c Single positive scale `c` (the additive->proportional
#'   crossover; the per-analyte SSD HC5, see [.analyte_c()]).
#' @return Numeric vector `asinh(impact / scale_c)`, same length as `impact`.
#'   `NA`/`Inf` propagate; `+/-Inf` map to `+/-Inf`.
#' @keywords internal
.g_transform <- function(impact, scale_c) {
  .assert_pos_scale(scale_c)
  asinh(impact / scale_c)
}

#' Inverse variance-stabilising transform
#'
#' Inverts [.g_transform()]: `g_inverse(g_transform(I, c), c) == I`.
#'
#' @param g Numeric vector on the transformed scale.
#' @param scale_c Single positive scale `c` (must match the value used in
#'   [.g_transform()]).
#' @return Numeric vector `scale_c * sinh(g)`.
#' @keywords internal
.g_inverse <- function(g, scale_c) {
  .assert_pos_scale(scale_c)
  scale_c * sinh(g)
}

#' Per-analyte transform scale `c` = SSD HC5
#'
#' Returns the analyte's 5% hazard concentration (HC5) from the fitted SSD, used
#' as the additive->proportional crossover `c` in [.g_transform()]. HC5 is a
#' self-scaling toxicological anchor in the same normalised concentration space
#' as the impact `I` (normalisation maps measured concentration onto the SSD
#' scale), and it lies within the species data range so it is well-determined
#' (unlike HC1, a lower-tail extrapolation). Below HC5 the transform is
#' additive, but that region carries PAF < 5% and so cannot inflate reported
#' risk.
#'
#' @details If HC5 cannot give a good fit to a particular analyte/site, HC1
#'   (`proportion = 0.01`) is the alternative to try -- a tighter crossover, but
#'   a less stable lower-tail extrapolation.
#'
#' @param fit A fitted `ssdtools` SSD object (e.g. an element of
#'   `derive_ssd_params()$fit`).
#' @return A single finite positive HC5 concentration.
#' @keywords internal
.analyte_c <- function(fit) {
  ## HC5 = SSD 5% hazard concentration; same call the package uses elsewhere
  ## (R/paf.R for the CA sigma). ci = FALSE -> deterministic point estimate.
  hc5 <- tryCatch(
    ssdtools::ssd_hc(fit, proportion = 0.05, ci = FALSE)$est,
    error = function(e) NA_real_
  )
  if (length(hc5) != 1L || !is.finite(hc5) || hc5 <= 0) {
    cli::cli_abort(c(
      "Could not derive a positive HC5 transform scale from the SSD fit.",
      "i" = "Expected a single finite positive {.fn ssdtools::ssd_hc} estimate."
    ))
  }
  hc5
}

#' Validate a transform scale
#' @param scale_c Candidate scale `c`.
#' @keywords internal
.assert_pos_scale <- function(scale_c) {
  if (!is.numeric(scale_c) || length(scale_c) != 1L ||
        is.na(scale_c) || scale_c <= 0) {
    cli::cli_abort("{.arg scale_c} must be a single positive number.")
  }
  invisible(scale_c)
}

#' Map an observation variance onto the g scale (delta method)
#'
#' The grab measurement error (S6) is specified as a variance in
#' concentration/impact space; the residual smoother now works on the
#' `g = asinh(x / c)` scale (issue #15), so the observation noise must be
#' transformed. By the delta method `Var(g) = Var(x) * g'(x)^2`, and with
#' `g'(x) = 1 / sqrt(x^2 + c^2)` this is `Var(x) / (x^2 + c^2)`. With
#' multiplicative grab error `Var(x) = (cv*x)^2` this plateaus at `cv^2` for
#' `|x| >> c` (proportional error) and vanishes at baseline.
#'
#' @param var_x Observation variance in concentration/impact space (numeric).
#' @param x The level at which the error is evaluated (impact `I` for the impact
#'   tier, concentration for the WQ tier).
#' @param scale_c Single positive transform scale `c`.
#' @return The observation variance on the `g` scale, `var_x / (x^2 + c^2)`.
#' @keywords internal
.s6_var_to_g <- function(var_x, x, scale_c) {
  .assert_pos_scale(scale_c)
  var_x / (x^2 + scale_c^2)
}
