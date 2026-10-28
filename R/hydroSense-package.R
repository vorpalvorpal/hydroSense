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

# mgcv must be LOADED whenever hydroSense is.  The package stores fitted
# mgcv::gam objects (reference_model$models[[.]]$gamm_fit, target-model
# wq_fit / impact_fit, imputation-model GAMs) and calls the stats generics on
# them -- stats::predict(), stats::AIC(), stats::vcov() and so on.  Those
# dispatch to predict.gam / logLik.gam / vcov.gam only when the mgcv namespace
# is loaded, because only then are mgcv's S3 methods registered.  A gam object
# also inherits from "glm" and "lm", so without mgcv loaded the generics fall
# through to the glm/lm methods: predict() errors ("lm object does not have a
# proper 'qr' component") and a fallback catches the error; AIC() and vcov()
# silently return glm answers.  A `mgcv::` call loads mgcv only when it runs,
# and a fresh R session -- a future::multisession worker, or a session that
# reads a saved model -- can predict on a shipped or saved model before any
# such call.  Importing one mgcv function makes loading hydroSense load mgcv.
#' @importFrom mgcv gam
NULL

# testthat can only mock a binding that exists in the package namespace. brms
# is an optional dependency (Suggests), and the specs that exercise its absence
# mock requireNamespace() through this placeholder, as documented in
# ?testthat::local_mocked_bindings. A non-function binding is skipped when R
# looks up a function call, so at run time requireNamespace() still resolves to
# base::requireNamespace().
requireNamespace <- NULL
