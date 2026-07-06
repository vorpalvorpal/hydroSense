# ── Analyte group constants ───────────────────────────────────────────────────

#' WQ analytes used in the PCA pre-processing step
#'
#' All candidates; the actual variables used at any given site are the
#' intersection of this list with analytes that pass prescreen in the training
#' data.  ORP and DO are included because they are valuable at sites where they
#' are measured, but at most sites they will be filtered out by prescreen.
#'
#' Note: the sulfate entry encodes its superscript charge (the SO4
#' two-minus symbol) with Unicode escape sequences in the source string, so
#' that matching works regardless of how the source file is parsed.
#' @keywords internal
.WQ_BLOCK_CANDIDATES <- c(
  # Field parameters / redox
  "temperature", "ORP", "DO",
  # Major cations
  "Ca", "Mg", "Na", "K",
  # Major anions & alkalinity
  "Cl", "SO4\u00b2\u207b", "Alkalinity-total-CaCO3", "F",
  # Hardness: total water hardness in mg/L as CaCO3 equivalents.
  # Callers must convert their data to this convention before passing in.
  "hardness",
  # Carbonate / alkalinity species
  "HCO3-CaCO3", "CO3-CaCO3", "OH-CaCO3",
  # Dissolved solids & suspended solids
  "TDS", "TSS",
  # Ionic balance totals
  "Anions-total", "Cations-total",
  # Organic carbon / oxygen demand (also used for organics hurdle)
  "DOC", "TOC", "BOD", "COD", "cBOD",
  # Nitrogen (excluding NH3-N which is a required driver)
  "NO2-N", "NO3-N", "NO2+NO3-N", "TKN-N", "N-total",
  # Sulfur
  "S",
  # Phosphorus
  "P-total", "P-reactive",
  # Redox-active metals as redox *predictors*. Dissolved Fe and Mn are the most
  # widely-measured redox indicators in leachate WQ (high dissolved Fe/Mn ==
  # reducing conditions, which mobilise some metals and precipitate others as
  # sulfides). Listing them here puts them in the WQ predictor block; the
  # `setdiff(all_analytes, union(pca_vars, exclude))` guard in
  # fit_imputation_model() consequently removes them from the imputation
  # *targets* (they predict the other metals instead of being imputed) — the
  # same dual role NH3-N already plays. They remain assessed by add_mspaf()
  # whenever measured.
  "Fe", "Mn"
)

#' Analytes that must never enter the imputation model as response variables
#'
#' The default `exclude` set for [fit_imputation_model()] (a leachate-preset
#' default; override via the `exclude` argument for other domains).  These are
#' typically non-concentration measurements (counts, qualitative, physical)
#' for which a log-normal concentration model is inappropriate, so they are
#' excluded from every imputation group.
#' @keywords internal
.IMPUTE_EXCLUDED <- c(
  # Microbiological counts (colony-forming units — not concentrations)
  "Coliforms",
  "Escherichia coli",
  "Faecal Coliforms",
  "Heterotrophic Plate Count (22\u00b0C)",
  "Heterotrophic Plate Count (36\u00b0C)",
  "E. coli",
  # Qualitative / physical descriptors
  "Appearance",
  "Colour",
  "Turbidity",
  "Stage"
)

#' All metal-type analytes (the metals group in [leachate_impute_groups()])
#'
#' A leachate-preset default: these define the `metals` [impute_group()]'s
#' targets and presence hurdle. For other domains, build your own groups with
#' [impute_group()] instead.
#' @keywords internal
.METAL_ANALYTES <- c(
  "Al", "As", "B", "Ba", "Be", "Cd", "Co", "Cr", "Cr-6", "Cu",
  "Fe", "Hg", "Mn", "Mo", "Ni", "Pb", "Sb", "Se", "Sn", "Sr", "V", "Zn"
)

#' Analytes that satisfy the organic-carbon hurdle for the organics group
#'
#' A leachate-preset default: the presence hurdle for the `organics`
#' [impute_group()] in [leachate_impute_groups()].
#' @keywords internal
.DOC_LIKE_ANALYTES <- c("DOC", "TOC", "BOD", "COD", "cBOD")

#' Co-analytes required by ANZECC/ANZG metal normalisation formulas
#'
#' These are imputed separately (after metals/organics imputation) via
#' [impute_coanalytes()] so that [add_mspaf()] has values to normalise
#' against.  pH and EC are excluded — they are always present (required vars).
#' @keywords internal
.COANALYTE_TARGETS <- c("DOC", "Ca", "Mg", "hardness")

#' PCA variables that must NOT be log-transformed before the chemistry PCA
#'
#' Every other PCA variable is concentration-like — strictly positive and
#' strongly right-skewed, spanning orders of magnitude — so it is
#' `log10`-transformed before centring/scaling.  Without that, the PCA is
#' dominated by a handful of high-magnitude major ions (e.g. Cl, SO4, TDS) and
#' the leading axes mostly track absolute ionic strength rather than the
#' multiplicative covariance structure that drives metal/organic behaviour.
#'
#' The exclusions are the variables for which a log is meaningless or undefined:
#'   - `pH` — already a logarithmic scale (−log10 of H+ activity).
#'   - `temperature` — interval scale (°C); zero/negative values are valid.
#'   - `ORP` — redox potential (mV); routinely negative.
#'   - `DO` — dissolved oxygen (mg/L); legitimately ~0 in anoxic leachate
#'     plumes and only spans a narrow, near-linear range.
#' @keywords internal
.PCA_NO_LOG_VARS <- c("pH", "temperature", "ORP", "DO")


# ── brms availability guard ───────────────────────────────────────────────────

#' Stop with a friendly, actionable message if brms is not installed
#'
#' The Bayesian imputation step ([fit_imputation_model()] /
#' [impute_chemistry()]) is the only part of the package that needs
#' \pkg{brms}, so brms is an optional ("Suggests") dependency rather than a
#' hard requirement.  This keeps the package quick to install for users who
#' only need the LMF or msPAF tools.  When someone actually calls an
#' imputation function without brms installed, this guard explains — in plain
#' language — what to install and why.
#' @keywords internal
.require_brms <- function() {
  if (requireNamespace("brms", quietly = TRUE)) {
    return(invisible(TRUE))
  }
  cli::cli_abort(c(
    "The chemistry imputation step needs the {.pkg brms} package, which \\
     isn't installed yet.",
    "i" = "{.pkg brms} fits the Bayesian model that fills in missing and \\
           below-detection-limit results. It's optional, so it isn't \\
           installed automatically \u2014 only this imputation step uses it.",
    " " = "",
    "*" = "To install it, run this once at the R console:",
    " " = "{.code install.packages(\"brms\")}",
    " " = "",
    "i" = "{.pkg brms} also needs a working Stan engine (a C++ compiler). If \\
           the line above isn't enough, follow the short setup guide at \\
           {.url https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}.",
    "i" = "Once {.pkg brms} is installed, re-run this function \u2014 no other \\
           changes are needed."
  ))
}

#' Choose the brms/Stan backend for imputation fits
#'
#' Prefers **cmdstanr** when both the package and a working CmdStan install are
#' present, otherwise falls back to **rstan**.  cmdstanr caches compiled model
#' binaries (a warm refit skips compilation: ~3x faster in the package
#' benchmark) and samples faster, while being statistically equivalent to rstan
#' (posterior-mean correlation 0.99999, max standardized difference ~0.12 SD on
#' the imputation model; see `dev/bench-backend.R`).  Override with
#' `options(hydroSense.brms_backend = "rstan")` (or `"cmdstanr"`).  An explicit
#' `backend` passed through `...` to [fit_imputation_model()] always wins.
#' @keywords internal
.brms_backend <- function() {
  opt <- getOption("hydroSense.brms_backend", NULL)
  if (!is.null(opt)) {
    return(match.arg(opt, c("cmdstanr", "rstan")))
  }
  has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(tryCatch(suppressMessages(cmdstanr::cmdstan_version()),
                      error = function(e) NULL))
  if (has_cmdstan) "cmdstanr" else "rstan"
}


# ── Imputation group specification ────────────────────────────────────────────

#' Declare an imputation group
#'
#' The imputation engine ([fit_imputation_model()]) is domain-agnostic: it
#' imputes one or more **groups** of target analytes from a shared
#' PCA-compressed chemistry context, with cross-target residual correlation
#' within each group.  `impute_group()` describes a single group — which
#' analytes it models and which (if any) presence hurdle gates it.  Pass a list
#' of these as the `groups` argument of [fit_imputation_model()].
#'
#' Each group is fitted as its own joint `brms` model, so analytes in different
#' groups do not borrow residual correlation from one another.  Group together
#' analytes you expect to co-vary (e.g. metals that move together in a plume).
#'
#' @param name Group label (a non-empty string).  Used as the name of the
#'   group's slot in the fitted model's `$groups` list and in console messages.
#' @param targets Character vector of analyte names to model jointly in this
#'   group, or `NULL` to mark this as the **catch-all** group, which claims
#'   every remaining modellable analyte (those not excluded, not used as PCA
#'   predictors, and not already claimed by an earlier group).  At most one
#'   catch-all group is allowed in a `groups` list.
#' @param hurdle Character vector of analyte names defining a *presence hurdle*,
#'   or `NULL` for no hurdle.  When a hurdle is set, a sample is only imputed
#'   for this group if it carries at least one of these analytes (detected or
#'   below-detection).  This stops silence being mistaken for absence — e.g. a
#'   sample with no metals recorded is left alone rather than given invented
#'   metal values.
#'
#' @return An object of class `"impute_group"`: a list with elements `name`,
#'   `targets`, and `hurdle`.
#'
#' @seealso [fit_imputation_model()], [leachate_impute_groups()]
#' @examples
#' # A metals group hurdled on metal presence, plus an everything-else group
#' # hurdled on dissolved-organic-carbon presence:
#' groups <- list(
#'   impute_group("metals", targets = c("Cu", "Zn", "Ni", "Pb"),
#'                hurdle = c("Cu", "Zn", "Ni", "Pb")),
#'   impute_group("organics", targets = NULL, hurdle = c("DOC", "TOC"))
#' )
#' @export
impute_group <- function(name, targets = NULL, hurdle = NULL) {
  if (!checkmate::test_string(name, min.chars = 1L)) {
    cli::cli_abort(c(
      "{.arg name} must be a single non-empty string.",
      "x" = "You supplied {.obj_type_friendly {name}}."
    ))
  }
  if (!is.null(targets)) {
    checkmate::assert_character(targets, min.len = 1L, any.missing = FALSE,
                               min.chars = 1L, .var.name = "targets")
    targets <- unique(targets)
  }
  if (!is.null(hurdle)) {
    checkmate::assert_character(hurdle, min.len = 1L, any.missing = FALSE,
                               min.chars = 1L, .var.name = "hurdle")
    hurdle <- unique(hurdle)
  }
  structure(
    list(name = name, targets = targets, hurdle = hurdle),
    class = "impute_group"
  )
}

#' @export
print.impute_group <- function(x, ...) {
  cat(sprintf("<impute_group: %s>\n", x$name))
  cat(sprintf("  targets: %s\n",
              if (is.null(x$targets)) "<catch-all>"
              else paste(x$targets, collapse = ", ")))
  cat(sprintf("  hurdle:  %s\n",
              if (is.null(x$hurdle)) "<none>"
              else paste(x$hurdle, collapse = ", ")))
  invisible(x)
}

#' Leachate imputation-group preset
#'
#' Returns the default [impute_group()] specification used by
#' [fit_imputation_model()] for landfill-leachate monitoring chemistry: a
#' **metals** group (hurdled on metal presence) and a catch-all **organics**
#' group (hurdled on dissolved-organic-carbon presence).  This is the leachate
#' domain layer on top of the otherwise domain-agnostic engine — pass your own
#' list of [impute_group()] objects to model a different chemistry.
#'
#' The metals set is [.METAL_ANALYTES] and the organic-carbon hurdle set is
#' [.DOC_LIKE_ANALYTES]. The metals presence hurdle excludes the redox
#' indicators Fe and Mn (routine analytes kept as PCA predictors, not trace-
#' metal contamination signals), so a sample reporting only Fe/Mn does not
#' trigger fabrication of the other trace metals.
#'
#' @return A list of two `"impute_group"` objects (`metals`, `organics`).
#' @seealso [impute_group()], [fit_imputation_model()]
#' @examples
#' leachate_impute_groups()
#' @export
leachate_impute_groups <- function() {
  list(
    impute_group("metals",   targets = .METAL_ANALYTES,
                 hurdle = setdiff(.METAL_ANALYTES, c("Fe", "Mn"))),
    impute_group("organics", targets = NULL,
                 hurdle = .DOC_LIKE_ANALYTES)
  )
}

#' Validate and normalise a list of imputation groups
#'
#' Checks that `groups` is a non-empty list of [impute_group()] objects with
#' unique names and at most one catch-all (`targets = NULL`) entry.
#' @param groups A list of `impute_group` objects.
#' @return The validated `groups` list (invisibly the same object).
#' @keywords internal
.validate_impute_groups <- function(groups) {
  if (!is.list(groups) || length(groups) == 0L) {
    cli::cli_abort("{.arg groups} must be a non-empty list of {.fn impute_group} objects.")
  }
  ok <- vapply(groups, inherits, logical(1L), what = "impute_group")
  if (!all(ok)) {
    cli::cli_abort(c(
      "Every element of {.arg groups} must be an {.fn impute_group} object.",
      "x" = "Element{?s} {.val {which(!ok)}} {?is/are} not."
    ))
  }
  nms <- vapply(groups, function(g) g$name, character(1L))
  if (anyDuplicated(nms)) {
    cli::cli_abort(c(
      "Group names in {.arg groups} must be unique.",
      "x" = "Duplicated: {.val {unique(nms[duplicated(nms)])}}."
    ))
  }
  n_catchall <- sum(vapply(groups, function(g) is.null(g$targets), logical(1L)))
  if (n_catchall > 1L) {
    cli::cli_abort(c(
      "At most one catch-all group ({.code targets = NULL}) is allowed.",
      "x" = "You supplied {n_catchall} catch-all groups."
    ))
  }
  invisible(groups)
}

#' Route a candidate analyte pool into imputation groups
#'
#' Explicit-target groups claim their analytes first, in declaration order
#' (so an analyte listed in two groups goes to the earlier one); the single
#' catch-all group (`targets = NULL`) then takes whatever remains.
#' @param candidate_pool Character vector of modellable analyte names.
#' @param groups A validated list of [impute_group()] objects.
#' @return A named list (one entry per group, named by `name`) of the analyte
#'   names assigned to each group.
#' @keywords internal
.route_groups <- function(candidate_pool, groups) {
  group_targets <- vector("list", length(groups))
  names(group_targets) <- vapply(groups, function(g) g$name, character(1L))
  pool         <- candidate_pool
  catchall_idx <- NA_integer_
  for (i in seq_along(groups)) {
    g <- groups[[i]]
    if (is.null(g$targets)) { catchall_idx <- i; next }
    claimed            <- intersect(pool, g$targets)
    group_targets[[i]] <- claimed
    pool               <- setdiff(pool, claimed)
  }
  if (!is.na(catchall_idx)) group_targets[[catchall_idx]] <- pool
  group_targets
}


# ── fit_imputation_model() ────────────────────────────────────────────────────

#' Fit the Bayesian multivariate imputation model(s)
#'
#' Fits a `brms` multivariate GAM for each **imputation group** (see
#' [impute_group()]), using a PCA-compressed water-quality (WQ) block as
#' additional environmental predictors.  The engine itself is domain-agnostic:
#' the leachate-specific groups (a `metals` group and a catch-all `organics`
#' group) are supplied by the default `groups = leachate_impute_groups()` and
#' can be swapped for any other chemistry by passing your own list of
#' [impute_group()] objects.  The returned model object is passed to
#' `impute_chemistry()` for prediction on new data.
#'
#' **Model structure**
#'
#' For each group, the mean structure is:
#' ```
#' s(PC1) + s(PC2) + ... + s(PCk)
#' ```
#' where `PC*` are the leading principal components of the unified chemistry
#' PCA (see *Chemistry PCA* below).  A group's target analytes are modelled
#' jointly with `rescor = TRUE`, so observed analytes at a given sample
#' condition the posterior of the missing ones through the residual correlation
#' matrix.  Separate groups are fitted as separate models and do not share
#' residual correlation.
#'
#' **Why `rescor = TRUE`** — the PCA captures the *instantaneous* chemical
#' covariance structure (what is measured together at a single moment), but
#' it cannot capture the *temporal-lag* covariance characteristic of
#' AMD/leachate-impacted aquifers, where conservative tracers move ahead of
#' redox-controlled metals.  At a post-pulse sample the PCA scores have
#' returned toward baseline but Cu/Pb/Zn/Mn remain elevated together — that
#' co-elevation is pure residual correlation with no predictor signal driving
#' it.  `rescor = TRUE` is the right machinery for this and is what makes
#' multivariate imputation borrow strength across analytes.
#'
#' **Costs of `rescor = TRUE`** — brms cannot combine `set_rescor(TRUE)` with
#' `cens("left")`, so this implementation uses `mi()` for BDL values and
#' applies a post-hoc cap (see [impute_chemistry()]).  The cap clips imputed
#' BDL cells to the original detection limit when the model predicts above
#' DL.  For sites where the chemistry context legitimately suggests high
#' concentrations the cap can fire frequently; results in that regime should
#' be inspected.  Three alternative configurations are worth benchmarking on
#' real hold-out data if predictive performance becomes a concern:
#'   - `rescor = TRUE` + `mi()` (current; expected to win on plume-affected
#'     groundwater because cross-analyte residual coupling captures plume
#'     dynamics that the predictor PCA misses).
#'   - `rescor = FALSE` + `cens("left")` (statistically clean for BDL; loses
#'     cross-analyte residual coupling).
#'   - `rescor = FALSE` + `cens("left")` + shared `(1 | sample_id)` (proper
#'     BDL handling with rank-1 latent-factor coupling across analytes).
#' Benchmark methodology: mask 10% of detected cells, fit each configuration,
#' compare hold-out RMSE / coverage.
#'
#' **Chemistry PCA**
#'
#' All `pca_vars` — major ions, pH, EC, NH3-N, DOC, nutrients, redox
#' indicators — that are present in `df` and pass a detection-frequency check
#' are submitted to `nipals::nipals()`, which handles within-sample missing
#' cells natively without prior imputation.  Using a single unified PCA (rather
#' than separate driver + WQ-block sets) eliminates predictor collinearity and
#' ensures normalisation co-analytes (DOC, Ca, Mg) influence the imputed metal
#' concentrations.  Principal components are added until cumulative variance
#' explained reaches `min_var_explained` or `max_pcs` is reached.  A minimum
#' of `min(2, available components)` PCs is used — i.e. the floor of two
#' only applies when at least two components exist.
#'
#' **Hurdles (applied at prediction time by `impute_chemistry()`)**
#'
#' Each group may carry a *presence hurdle* (see [impute_group()]): a sample is
#' only imputed for that group if it carries at least one of the hurdle
#' analytes (detected or BDL).  Under the leachate preset the metals group is
#' hurdled on metal presence and the organics group on DOC-like presence.
#'
#' **BDL required variables**
#'
#' When a `required_vars` analyte (pH or EC) is below the detection limit for a
#' sample, the stored detection-limit value is used as-is (conservative upper
#' bound).  A message is issued but the sample is retained.
#'
#' @param df Long-format chemistry data frame with columns `sample_id`,
#'   `site_id`, `datetime`, `analyte`, `value`, `detected`.
#' @param pca_vars Analyte names to include in the unified chemistry PCA (used
#'   as predictors for the brms model via PC scores).  Default: `c("pH", "EC",
#'   "NH3-N")` plus all `.WQ_BLOCK_CANDIDATES`.  Normalisation co-analytes
#'   (DOC, Ca, Mg) are included in the default set.
#' @param required_vars Analyte names that must be present in a sample for it
#'   to be retained in training and prediction.  Default `c("pH", "EC")`.
#'   Samples missing any of these are dropped entirely.
#' @param groups A list of [impute_group()] objects describing which analytes
#'   to impute and how each group is hurdled.  Default `NULL` uses
#'   [leachate_impute_groups()] (a `metals` group plus a catch-all `organics`
#'   group).  Supply your own list to impute a different chemistry.
#' @param exclude Analyte names that must never be modelled as response
#'   variables in any group (e.g. counts, qualitative descriptors).  Default
#'   `NULL` uses [.IMPUTE_EXCLUDED].
#' @param no_log_vars Analyte names that must **not** be log-transformed before
#'   the chemistry PCA (interval-scale or already-logarithmic variables such as
#'   pH and temperature).  Default `NULL` uses [.PCA_NO_LOG_VARS].
#' @param min_target_detect_freq Minimum detection frequency (fraction of
#'   samples in which the analyte is *detected*) for a metal/organic to be
#'   included as an imputation target. Targets below this are dropped (they have
#'   too few detections to model and would otherwise inflate the model on
#'   near-all-BDL panels). Default `0.05`. Combined with `min_target_detect_n`
#'   (both gates must pass).
#' @param min_target_detect_n Minimum **absolute** number of distinct samples in
#'   which a metal/organic is detected for it to be included as an imputation
#'   target. The fraction gate above already implies a count floor of
#'   `min_target_detect_freq * n_samples`, but that scales with dataset size and
#'   collapses on small datasets; this absolute floor guarantees enough anchors
#'   to constrain the fit regardless of dataset size. Default `4L` (non-binding
#'   on typical panels, where the fraction gate dominates).
#' @param min_detect_freq Minimum detection frequency for a PCA variable to be
#'   retained.  Default `0.05`.  Required vars are always retained regardless.
#' @param min_samples Minimum training samples after required-var filtering.
#' @param min_var_explained Target cumulative variance for PCA axis selection.
#'   Default `0.75`.
#' @param max_pcs Maximum PCA axes to use.  Default `6L`.
#' @param family brms response family.  Must be `"gaussian"` (concentrations
#'   are log-transformed before fitting; residual correlations require
#'   Gaussian family).
#' @param impute_method How below-detection (BDL) values and cross-analyte
#'   coupling are handled. One of:
#'   \describe{
#'     \item{`"rescor_mi"`}{(default) Residual correlation across analytes
#'       (`rescor = TRUE`) with BDL/missing treated as imputable (`mi()`); the
#'       imputed BDL cells are capped at the detection limit post-hoc by
#'       [impute_chemistry()] (brms cannot combine `rescor` with `cens()`). Most
#'       accurate recovery (best hold-out RMSE by a wide margin), but the `mi()` +
#'       correlation geometry is funnel-prone, so `adapt_delta = 0.95` and an
#'       `lkj(2)` prior on the residual correlation are set by default (override
#'       via `control` / `prior` in `...`). Even so the geometry stays hard
#'       (tree-depth saturation, low E-BFMI, worst-case R̂ ≈ 1.6 on a hard mask):
#'       trust the **point estimate**, but check `brms::rhat()` before relying on
#'       the **draws** — for well-calibrated uncertainty prefer `"cens_factor"`.}
#'     \item{`"cens"`}{Proper left-censoring of BDL at the detection limit
#'       (`cens("left")`), no residual correlation -- clean BDL handling but no
#'       cross-analyte coupling.}
#'     \item{`"cens_factor"`}{Proper left-censoring **with** cross-analyte
#'       coupling. Fitted as a single long-format model with a shared per-sample
#'       latent factor `(1 | sample_id)` common to all analytes, so an observed
#'       metal informs the unobserved/BDL ones at that sample. The factor is
#'       well-identified (each sample contributes several analyte observations);
#'       `adapt_delta = 0.95` is set by default to clear the factor's mild
#'       funnel (override via `control` in `...`).}
#'     \item{`"factor"`}{**(Route C)** Low-rank left-censored latent factor
#'       model (`Sigma = Lambda Lambda' + Psi`, rank `k = 2` by default),
#'       resolving findings 1-3 at the source: BDL cells are censored in the
#'       likelihood (no post-hoc cap), and the latent factor is inferred from a
#'       sample's *observed* metals at prediction time, so measured metals
#'       genuinely condition the missing/BDL ones. Fitted in two stages: a
#'       per-analyte `mgcv::gam` mean on the PC scores, then a Stan factor
#'       model on the residuals (needs \pkg{cmdstanr}). See
#'       `dev/plan-route-c.md`.}
#'   }
#'   See `vignette("imputation")` and the package benchmark for guidance on
#'   which to prefer.
#' @param iter,warmup,chains,cores brms MCMC settings.
#' @param save_dir If non-NULL, save the returned model object as a `.qs` file
#'   in this directory using `qs2::qs_save()`.
#' @param ... Additional arguments passed to `brms::brm()`.  The Stan
#'   **`backend`** defaults to `"cmdstanr"` when the \pkg{cmdstanr} package and a
#'   CmdStan install are both available (cached compiled binaries + faster
#'   sampling, statistically equivalent to rstan), otherwise `"rstan"`; pass
#'   `backend = ...` here or set `options(hydroSense.brms_backend = ...)` to
#'   override.
#'
#' @return A named list of class `"imputation_model"`:
#'   - `$pca`: PCA fit + metadata (loadings, training medians, n_pcs,
#'     `no_log_vars`, …)
#'   - `$groups`: a named list (one entry per fitted group, named by the group's
#'     `name`); each entry has `$fit` (brmsfit), `$analytes`, `$safe_names`,
#'     `$name`, `$hurdle`.  Empty if no group had any modellable analytes.
#'   - `$group_specs`: the input list of [impute_group()] objects
#'   - `$required_vars`, `$pca_vars`, `$exclude`, `$impute_method`
#'   - `$fit_date`, `$n_samples`: metadata
#'   If `save_dir` is supplied, the path to the saved file is returned as
#'   `attr(result, "save_path")`.
#'
#' @seealso [impute_chemistry()], [impute_group()], [leachate_impute_groups()]
#' @examples
#' \dontrun{
#' # Requires a Stan toolchain (brms). Fit once, then reuse for imputation.
#' model <- fit_imputation_model(monitoring_long)
#' draws <- impute_chemistry(monitoring_long, model, return = "draws")
#'
#' # A different domain: two custom groups instead of the leachate preset.
#' model2 <- fit_imputation_model(
#'   monitoring_long,
#'   groups = list(
#'     impute_group("trace_metals", targets = c("As", "Cd", "Pb"),
#'                  hurdle = c("As", "Cd", "Pb")),
#'     impute_group("nutrients", targets = NULL, hurdle = c("NO3-N", "P-total"))
#'   )
#' )
#' }
#' @export
fit_imputation_model <- function(
    df,
    pca_vars          = NULL,           # default built in body
    required_vars     = c("pH", "EC"),
    groups            = NULL,
    exclude           = NULL,
    no_log_vars       = NULL,
    min_detect_freq   = 0.05,
    min_target_detect_freq = 0.05,
    min_target_detect_n = 4L,
    min_samples       = 10L,
    min_var_explained = 0.75,
    max_pcs           = 6L,
    family            = "gaussian",
    impute_method     = c("rescor_mi", "cens", "cens_factor", "factor"),
    iter              = 2000,
    warmup            = 1000,
    chains            = 4,
    cores             = parallel::detectCores(),
    save_dir          = NULL,
    ...
) {
  .require_brms()
  impute_method <- match.arg(impute_method)

  if (is.null(pca_vars))    pca_vars    <- c("pH", "EC", "NH3-N",
                                             .WQ_BLOCK_CANDIDATES)
  if (is.null(groups))      groups      <- leachate_impute_groups()
  if (is.null(exclude))     exclude     <- .IMPUTE_EXCLUDED
  if (is.null(no_log_vars)) no_log_vars <- .PCA_NO_LOG_VARS
  .validate_impute_groups(groups)
  # Guard declared group targets against a make.names() collision up front —
  # independent of whether they clear the detection-frequency gates below, so
  # a mistyped group spec fails fast rather than silently dropping analytes.
  declared_targets <- unique(unlist(lapply(groups, function(g) g$targets), use.names = FALSE))
  if (length(declared_targets) > 0L) .assert_safe_analyte_names(declared_targets)

  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected"))
  checkmate::assert_character(required_vars, min.len = 1L, any.missing = FALSE)
  checkmate::assert_character(exclude, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_character(no_log_vars, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_count(min_samples)
  checkmate::assert_number(min_target_detect_freq, lower = 0, upper = 1)
  checkmate::assert_count(min_target_detect_n)

  # ── 1. BDL required-variable handling ────────────────────────────────────
  # For required vars (pH, EC) where a value is BDL, use the stored DL value.
  # These are genuine low-level conditions; the sample is retained.
  n_bdl_req <- df |>
    dplyr::filter(.data$analyte %in% .env$required_vars, !.data$detected) |>
    nrow()
  if (n_bdl_req > 0L) {
    cli::cli_inform(c(
      "i" = "{n_bdl_req} BDL row{?s} for required variable{?s} \u2014 using \\
             detection-limit value{?s} as conservative estimate."
    ))
  }
  # `value` already holds the DL for BDL rows; no transformation needed.

  # ── 2. Drop samples missing any required variable entirely ────────────────
  samples_with_required <- df |>
    dplyr::filter(.data$analyte %in% .env$required_vars) |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::summarise(n_req = dplyr::n_distinct(.data$analyte), .groups = "drop") |>
    dplyr::filter(.data$n_req == length(required_vars)) |>
    dplyr::pull(.data$sample_id)

  n_dropped <- dplyr::n_distinct(df$sample_id) - length(samples_with_required)
  if (n_dropped > 0L) {
    cli::cli_inform(c(
      "!" = "{n_dropped} sample{?s} dropped: missing one or more required \\
             variable{?s} ({.val {required_vars}}) entirely."
    ))
    df <- dplyr::filter(df, .data$sample_id %in% samples_with_required)
  }

  if (length(samples_with_required) < min_samples) {
    cli::cli_abort(c(
      "Only {length(samples_with_required)} sample{?s} remain after \\
       required-var filtering \u2014 fewer than {.arg min_samples} = {min_samples}.",
      "i" = "Lower {.arg min_samples}, add more data, or choose different \\
             {.arg required_vars}."
    ))
  }

  # ── 3. Filter pca_vars by presence frequency ──────────────────────────────
  # Required vars always pass regardless of presence frequency.
  # NB: this is a PRESENCE frequency (fraction of samples that have a row for
  # the analyte at all), not a detection frequency — a PCA variable is useful
  # as a predictor whether or not it was above the detection limit, so BDL
  # rows count towards retention here. (Contrast prescreen_analytes(), which
  # screens toxicants on true detection frequency.)
  n_samples_total <- dplyr::n_distinct(df$sample_id)
  pca_vars_present <- df |>
    dplyr::filter(.data$analyte %in% .env$pca_vars) |>
    dplyr::group_by(.data$analyte) |>
    dplyr::summarise(
      presence_freq = dplyr::n_distinct(.data$sample_id) / n_samples_total,
      .groups = "drop"
    ) |>
    dplyr::filter(.data$presence_freq >= min_detect_freq) |>
    dplyr::pull(.data$analyte)

  pca_vars_present <- union(
    intersect(required_vars, unique(df$analyte)),
    pca_vars_present
  )

  cli::cli_inform(c(
    "i" = "Chemistry PCA: {length(pca_vars_present)} variable{?s} available: \\
           {.val {sort(pca_vars_present)}}."
  ))

  # ── 4. Assign analytes to groups ───────────────────────────────────────────
  # Candidate pool = everything not used as a PCA predictor and not excluded
  # (microbiological counts, qualitative/physical descriptors — not amenable to
  # a log-normal model).
  all_analytes <- unique(df$analyte)
  candidate_pool <- setdiff(all_analytes, union(pca_vars, exclude))

  excl_present <- intersect(all_analytes, exclude)
  if (length(excl_present) > 0L) {
    cli::cli_inform(c(
      "i" = "{length(excl_present)} analyte{?s} excluded from imputation \\
             (non-concentration data): {.val {sort(excl_present)}}."
    ))
  }

  # Route candidates into groups: explicit-target groups claim first (in the
  # order declared), then the single catch-all group (targets = NULL) takes
  # whatever remains in the pool.
  group_targets <- .route_groups(candidate_pool, groups)

  # Drop target analytes detected too rarely to model. A brms regression needs
  # enough *detected* observations; near-/all-BDL analytes (e.g. ~100 organics
  # in a leachate panel) carry no signal and otherwise explode the model size.
  # A target must clear BOTH gates: a detection *fraction* and an *absolute*
  # detection count. The fraction (`det_freq = n_detect / n_samples_total`)
  # already implies a count floor of `min_target_detect_freq * n_samples_total`,
  # but that scales with dataset size and collapses on small datasets; the
  # absolute `min_target_detect_n` guarantees a minimum number of anchors to
  # constrain the fit regardless of dataset size. (#59 item 4 / #61.)
  all_targets <- unlist(group_targets, use.names = FALSE)
  det_tbl <- df |>
    dplyr::filter(.data$analyte %in% .env$all_targets) |>
    dplyr::group_by(.data$analyte) |>
    dplyr::summarise(
      n_detect = dplyr::n_distinct(.data$sample_id[.data$detected]),
      det_freq = dplyr::n_distinct(.data$sample_id[.data$detected]) / n_samples_total,
      .groups  = "drop"
    )
  pass_freq  <- det_tbl$det_freq  >= min_target_detect_freq
  pass_n     <- det_tbl$n_detect  >= min_target_detect_n
  keep_targets  <- det_tbl$analyte[pass_freq & pass_n]
  dropped       <- setdiff(all_targets, keep_targets)
  group_targets <- lapply(group_targets, intersect, keep_targets)
  if (length(dropped) > 0L) {
    # Attribute each drop to the gate(s) it failed, for an auditable message.
    drop_n_only <- det_tbl$analyte[pass_freq & !pass_n]
    why <- if (length(drop_n_only))
      sprintf(" (%s below min_target_detect_n = %d)",
              paste(sort(drop_n_only), collapse = ", "), min_target_detect_n)
    else ""
    cli::cli_inform(c(
      "i" = "Dropping {length(dropped)} target{?s} below \\
             min_target_detect_freq = {min_target_detect_freq} / \\
             min_target_detect_n = {min_target_detect_n}: \\
             {.val {sort(dropped)}}{why}."
    ))
  }

  for (nm in names(group_targets)) {
    cli::cli_inform(c(
      "i" = "{nm} group: {length(group_targets[[nm]])} analyte{?s}: \\
             {.val {sort(group_targets[[nm]])}}."
    ))
  }

  if (sum(lengths(group_targets)) == 0L) {
    cli::cli_warn(
      "No target analytes found outside the PCA and excluded sets. \\
       Returning model with no fitted groups (imputation will be a no-op)."
    )
    return(structure(
      list(
        pca           = NULL,
        groups        = list(),
        group_specs   = groups,
        required_vars = required_vars,
        pca_vars      = pca_vars,
        exclude       = exclude,
        impute_method = impute_method,
        fit_date      = Sys.Date(),
        n_samples     = length(samples_with_required)
      ),
      class = "imputation_model"
    ))
  }

  # ── 5. Fit unified chemistry PCA ──────────────────────────────────────────
  pca_obj <- .prepare_chem_pca(
    df, wq_vars        = pca_vars_present,
    min_var_explained  = min_var_explained,
    max_pcs            = max_pcs,
    no_log_vars        = no_log_vars
  )

  cli::cli_inform(c(
    "i" = "Chemistry PCA: {pca_obj$n_pcs} axis/axes explain \\
           {round(100 * pca_obj$var_explained, 1)}% of variance."
  ))

  # ── 6. Fit one model per non-empty group ───────────────────────────────────
  fitted_groups <- list()
  for (g in groups) {
    tgts <- group_targets[[g$name]]
    if (length(tgts) == 0L) next
    cli::cli_inform(c("i" = "Fitting {g$name} model \u2026"))
    fit <- .fit_group_model(
      df              = df,
      target_analytes = tgts,
      pca_obj         = pca_obj,
      family          = family,
      iter            = iter,
      warmup          = warmup,
      chains          = chains,
      cores           = cores,
      impute_method   = impute_method,
      group_name      = g$name,
      ...
    )
    fit$name   <- g$name
    fit$hurdle <- g$hurdle
    fitted_groups[[g$name]] <- fit
  }

  # ── 7. Assemble result ──────────────────────────────────────────────────────
  result <- structure(
    list(
      pca           = pca_obj,
      groups        = fitted_groups,
      group_specs   = groups,
      required_vars = required_vars,
      pca_vars      = pca_vars,
      exclude       = exclude,
      impute_method = impute_method,
      fit_date      = Sys.Date(),
      n_samples     = length(samples_with_required)
    ),
    class = "imputation_model"
  )

  # ── 9. Save if requested ───────────────────────────────────────────────────
  if (!is.null(save_dir)) {
    if (!requireNamespace("qs2", quietly = TRUE))
      cli::cli_abort("Package {.pkg qs2} is required for saving models. \\
                      Install with: {.code install.packages('qs2')}")
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    fname     <- sprintf("imputation_model_%s.qs", format(Sys.Date(), "%Y%m%d"))
    save_path <- file.path(save_dir, fname)
    qs2::qs_save(result, save_path)
    cli::cli_inform(c("v" = "Model saved to {.path {save_path}}"))
    attr(result, "save_path") <- save_path
  }

  result
}

#' @export
print.imputation_model <- function(x, ...) {
  if (is.null(x$pca)) {
    cat(sprintf(
      "<imputation_model>  fitted %s | %d samples | %d PCA vars | no fitted groups\n",
      x$fit_date, x$n_samples, length(x$pca_vars)
    ))
    return(invisible(x))
  }
  cat(sprintf(
    "<imputation_model>  fitted %s | %d samples | %d PCA vars | %d PCs (%.0f%% var)\n",
    x$fit_date, x$n_samples, length(x$pca_vars),
    x$pca$n_pcs, 100 * x$pca$var_explained
  ))
  if (!is.null(x$impute_method))
    cat(sprintf("  method:   %s\n", x$impute_method))
  groups <- x$groups %||% list()
  if (length(groups)) {
    width <- max(nchar(names(groups)))
    for (nm in names(groups)) {
      cat(sprintf("  %-*s %d analytes\n",
                  width + 1L, paste0(nm, ":"), length(groups[[nm]]$analytes)))
    }
  }
  invisible(x)
}


# ── impute_chemistry() ────────────────────────────────────────────────────────

#' Impute missing and BDL chemistry using a fitted imputation model
#'
#' Applies the models fitted by [fit_imputation_model()] to `df`, returning
#' posterior estimates for below-detection-limit (BDL) and missing observations
#' in every fitted group.
#'
#' **Completing the panel**
#'
#' For each eligible sample (see *Hurdles*), every target analyte the group
#' models is filled: BDL cells are replaced with their posterior estimate, and
#' analytes that are **entirely absent** for a sample gain a new
#' model-anchored row (`imputed_kind = "missing"`, `detected = TRUE`). This is
#' what lets a well-sampled analyte lift a sparsely-sampled one — e.g. a sample
#' with Zn but no Cu gains a Cu row predicted from the fitted Cu–Zn
#' relationship. Fabricated rows carry the originating sample's `site_id`,
#' `datetime`, and any other sample-level columns.
#'
#' Fabricated rows are model predictions, not measurements. With
#' `return = "draws"` the imputation-model uncertainty propagates through the
#' draw carrier; with `return = "point"` the single anchor enters any
#' downstream model as if observed, which can understate that model's
#' uncertainty. Prefer `return = "draws"` when the imputed values feed a
#' further model (e.g. the reference/target GAMs).
#'
#' **Hurdles**
#'
#' Each fitted group may carry a presence hurdle (see [impute_group()]).  When
#' `apply_hurdles = TRUE`, imputed values for a group are only returned for
#' samples carrying at least one of that group's hurdle analytes (detected or
#' BDL) — e.g. under the leachate preset, a sample with no metals recorded is
#' not given imputed metals, because a leachate metal pulse may simply not have
#' arrived at that location yet.  Samples failing a hurdle pass through
#' unchanged (non-imputed values preserved; BDL values remain flagged as BDL).
#'
#' @param df Long-format chemistry data frame (same schema as used for fitting).
#' @param model Fitted model from [fit_imputation_model()].
#' @param apply_hurdles Logical.  Apply each group's presence hurdle?  Default
#'   `TRUE`.  When `FALSE`, every sample is eligible for every group.
#' @param bdl_cap Logical.  Cap imputed BDL values at the original detection
#'   limit?  Default `TRUE`.  Applied to all methods: an imputed below-detection
#'   value should never exceed its detection limit.
#' @param return `"point"` (default) for posterior mean per cell; `"draws"` for
#'   one row per (sample × analyte × draw).
#' @param ndraws Integer or `NULL`.  Use only this many posterior draws for
#'   prediction (subsampled).  `NULL` (default) uses all draws.  Lowering it
#'   reduces memory/time, at some cost to interval precision.
#' @param batch_size Integer or `NULL`.  Predict eligible samples in batches of
#'   this many rows to bound peak memory (important for `"rescor_mi"`, whose
#'   `mi()` prediction is memory-heavy).  `NULL` (default) predicts all at once.
#'
#' @return `df` with BDL and missing cells in every fitted group
#'   replaced by posterior mean estimates, plus columns:
#'   - `imputed` (logical) — `TRUE` for filled cells
#'   - `imputed_kind` — `"observed"`, `"censored_left"`, or `"missing"`
#'
#'   When `bdl_cap = TRUE` and any imputed BDL cell exceeded its detection
#'   limit, a per-cell audit of the cap activations is attached as the
#'   `"bdl_cap_summary"` attribute; retrieve it with [bdl_cap_summary()].
#'
#' @seealso [fit_imputation_model()], [bdl_cap_summary()]
#' @examples
#' \dontrun{
#' model <- fit_imputation_model(monitoring_long)
#' # Point estimates (default), or full posterior draws with return = "draws":
#' imputed <- impute_chemistry(monitoring_long, model, return = "point")
#' }
#' @export
impute_chemistry <- function(
    df,
    model,
    apply_hurdles  = TRUE,
    bdl_cap        = TRUE,
    return         = c("point", "draws"),
    ndraws         = NULL,
    batch_size     = NULL
) {
  return <- match.arg(return)
  .require_brms()
  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected"))
  if (!inherits(model, "imputation_model"))
    cli::cli_abort("{.arg model} must be an object returned by {.fn fit_imputation_model}.")

  groups <- model$groups %||% list()
  if (length(groups) == 0L) {
    cli::cli_warn("Model has no fitted groups; returning {.arg df} unchanged.")
    result <- df
    if (!"imputed" %in% names(result))
      result <- dplyr::mutate(result, imputed = FALSE, imputed_kind = "observed")
    return(result)
  }

  # ── Compute WQ PC scores for new data ─────────────────────────────────────
  pca_scores <- .compute_pca_scores(df, model$pca)

  # ── Collect BDL detection limits across every group's targets ─────────────
  all_targets <- unlist(lapply(groups, function(g) g$analytes), use.names = FALSE)
  dl_tbl <- df |>
    dplyr::filter(.data$analyte %in% .env$all_targets, !.data$detected) |>
    dplyr::select("sample_id", "analyte", detection_limit = "value")

  # ── Impute each group ──────────────────────────────────────────────────────
  result <- df  # start with original; overlay imputed values below

  for (g in groups) {
    eligible <- if (apply_hurdles && !is.null(g$hurdle)) {
      df |>
        dplyr::filter(.data$analyte %in% .env$g$hurdle) |>
        dplyr::pull(.data$sample_id) |>
        unique()
    } else {
      unique(df$sample_id)
    }
    n_skip <- dplyr::n_distinct(df$sample_id) - length(eligible)
    if (n_skip > 0L)
      cli::cli_inform(c(
        "i" = "{g$name} hurdle: skipping {n_skip} sample{?s} (no \\
               {g$name}-group analyte present)."
      ))

    result <- .predict_and_merge(
      df           = result,
      group        = g,
      pca_scores   = pca_scores,
      eligible_ids = eligible,
      return       = return,
      ndraws       = ndraws,
      batch_size   = batch_size
    )
  }

  # Tag non-target rows (drivers, WQ vars, etc.) that were never imputed
  if (!"imputed" %in% names(result)) {
    result <- dplyr::mutate(result,
      imputed      = FALSE,
      imputed_kind = "observed"
    )
  } else {
    result <- dplyr::mutate(result,
      imputed      = dplyr::coalesce(.data$imputed, FALSE),
      imputed_kind = dplyr::coalesce(.data$imputed_kind, "observed")
    )
  }

  # Cap imputed BDL values at their detection limit, for every method. An
  # imputed below-detection value must not exceed its DL. (For rescor_mi this
  # is the only enforcement of the censoring bound; for the cens methods the
  # bound is enforced in the likelihood during fitting, but the emitted
  # prediction is unconstrained, so the cap is still needed.)
  .check_bdl_imputed(result, dl_tbl, bdl_cap)
}


# ── impute_coanalytes() ───────────────────────────────────────────────────────

#' Impute missing normalisation co-analytes from the fitted chemistry PCA
#'
#' Fits a univariate log-Gaussian GAM (`mgcv::gam`) for each target
#' co-analyte using the PC scores already computed by
#' [fit_imputation_model()].  Only samples where the co-analyte is entirely
#' absent are filled; BDL observations are left unchanged.
#'
#' This step belongs **after** [impute_chemistry()] and **before**
#' [time_weighted_aggregate()].  Imputed co-analyte values are never fed
#' back into the metals/organics model — the brms model ran on measured values
#' only and is already done.
#'
#' Using the chemistry PCA as the sole predictor set is appropriate because:
#' (a) the PCA already captures DOC/Ca/Mg variation in its axes; (b) a
#' univariate GAM on PC scores is unbiased and fast (no Stan required); (c)
#' the same PCA is used for the metals model so the co-analyte predictions
#' are conditioned on the same chemical environment summary.
#'
#' @param df Long-format chemistry data frame (same schema as
#'   [impute_chemistry()], with `imputed`/`imputed_kind` columns if
#'   [impute_chemistry()] has already been called).
#' @param model Fitted model from [fit_imputation_model()] (provides the PCA
#'   object and the list of `pca_vars`).
#' @param targets Co-analyte names to impute when missing.  Default
#'   `.COANALYTE_TARGETS` (`"DOC"`, `"Ca"`, `"Mg"`, `"hardness"`).  Only
#'   targets present in `model$pca_vars` are processed; others are
#'   skipped with a warning.
#' @param min_obs Minimum number of quantified observations required to fit a
#'   GAM for a target.  Targets with fewer observations are skipped.
#'   Default `10L`.
#' @param return `"point"` (default) for the posterior mean of the GAM
#'   prediction — identical to pre-draws behaviour.  `"draws"` for
#'   full posterior-predictive draws: each missing co-analyte cell emits
#'   `N` rows keyed by `draw_id 1..N`, reflecting both parameter uncertainty
#'   (`beta ~ N(coef(gam), Vp)`) and residual Gaussian noise (`gam$sig2`).
#'   Observed co-analyte cells keep `draw_id = NA`.
#' @param ndraws Number of draws to generate.  Required when
#'   `return = "draws"` and `df` contains no existing draws.  When `df`
#'   already carries draws (from [impute_chemistry()]), `N` is inferred
#'   from the existing draw domain; `ndraws` must be `NULL` or equal to
#'   that count.
#' @param seed Optional integer seed passed to [set.seed()] before the
#'   sampling calls, for reproducibility.
#'
#' @return `df` with missing co-analyte rows filled in, tagged with
#'   `imputed = TRUE` and `imputed_kind = "missing"`.  In `"draws"` mode
#'   each imputed cell is replicated `N` times with `draw_id 1..N`;
#'   observed cells keep `draw_id = NA`.  In `"point"` mode the output
#'   schema is unchanged from the pre-draws behaviour.
#'
#' @seealso [fit_imputation_model()], [impute_chemistry()], [summarise_draws()]
#' @examples
#' \dontrun{
#' # Deterministic GAM-based imputation (default, point mode)
#' impute_coanalytes(monitoring_long, model)
#'
#' # Posterior-predictive draws when df already carries metals draws
#' impute_coanalytes(metals_draws, model, return = "draws")
#' }
#' @export
impute_coanalytes <- function(
    df,
    model,
    targets  = NULL,
    min_obs  = 10L,
    return   = c("point", "draws"),
    ndraws   = NULL,
    seed     = NULL
) {
  return <- match.arg(return)
  if (is.null(targets)) targets <- .COANALYTE_TARGETS
  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected"))
  if (!inherits(model, "imputation_model"))
    cli::cli_abort(
      "{.arg model} must be an object returned by {.fn fit_imputation_model}."
    )
  if (is.null(model$pca))
    cli::cli_abort(
      "Model has no fitted PCA \u2014 did {.fn fit_imputation_model} find any \\
       target analytes?"
    )

  # ── Resolve draw count N ──────────────────────────────────────────────────
  domain <- .draw_domain(df)
  if (return == "draws") {
    if (length(domain) > 0L) {
      if (!is.null(ndraws) && as.integer(ndraws) != length(domain)) {
        cli::cli_abort(c(
          "{.arg ndraws} = {ndraws} conflicts with the input frame's draw \\
           domain (N = {length(domain)}).",
          "i" = "Omit {.arg ndraws} to reuse the existing draw count."
        ))
      }
      N <- length(domain)
    } else {
      if (is.null(ndraws))
        cli::cli_abort(c(
          "{.arg ndraws} is required when {.arg return = \"draws\"} and the \\
           input frame carries no draws.",
          "i" = "E.g. {.code ndraws = 500L}."
        ))
      N <- as.integer(ndraws)
    }
    if (!is.null(seed)) set.seed(seed)
  }

  # ── Compute PC scores for all samples ─────────────────────────────────────
  # .compute_pca_scores() uses mean() per cell, so scores are deterministic
  # even when df carries draws — the draw_id column does not corrupt scoring.
  pca_scores <- .compute_pca_scores(df, model$pca)
  pc_cols    <- paste0("PC", seq_len(model$pca$n_pcs))

  # Skip targets not represented in the PCA (they can't be predicted)
  targets_ok <- intersect(targets, model$pca_vars)
  skipped    <- setdiff(targets, model$pca_vars)
  if (length(skipped) > 0L)
    cli::cli_warn(c(
      "!" = "{length(skipped)} co-analyte target{?s} not in fitted \\
             {.arg pca_vars} \u2014 skipping: {.val {skipped}}."
    ))

  # Per-sample metadata for constructing new rows
  sample_meta <- df |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select("sample_id", "site_id", "datetime")

  result <- df

  for (tgt in targets_ok) {

    # "Present" means at least one row exists (detected or BDL)
    present_ids <- unique(result$sample_id[result$analyte == tgt])
    missing_ids <- setdiff(unique(result$sample_id), present_ids)

    if (length(missing_ids) == 0L) next   # already complete

    # Count quantified observations available for GAM fitting
    n_obs <- dplyr::n_distinct(
      result$sample_id[result$analyte == tgt & result$detected]
    )
    if (n_obs < min_obs) {
      cli::cli_warn(c(
        "!" = "Co-analyte {.val {tgt}}: only {n_obs} quantified sample{?s} \\
               (< {.arg min_obs} = {min_obs}) \u2014 skipping."
      ))
      next
    }

    cli::cli_inform(c(
      "i" = "Co-analyte {.val {tgt}}: imputing {length(missing_ids)} \\
             missing sample{?s} via GAM on {model$pca$n_pcs} PC score{?s}."
    ))

    # ── Fit GAM on quantified observations ──────────────────────────────────
    # Floor tied to this co-analyte's own scale (finding 6), not one constant
    # shared across analytes of wildly different magnitude.
    tgt_floor <- .scale_aware_log_floor(
      result$value[result$analyte == tgt & result$detected]
    )
    obs_vals <- result |>
      dplyr::filter(.data$analyte == tgt, .data$detected) |>
      dplyr::group_by(.data$sample_id) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::transmute(.data$sample_id,
                        log_tgt = log(pmax(.data$value, tgt_floor)))

    gam_data <- pca_scores |>
      dplyr::filter(.data$sample_id %in% present_ids) |>
      dplyr::left_join(obs_vals, by = "sample_id") |>
      dplyr::filter(!is.na(.data$log_tgt))

    gam_formula <- stats::as.formula(
      paste("log_tgt ~",
            paste(paste0("s(", pc_cols, ")"), collapse = " + "))
    )

    gam_fit <- tryCatch(
      mgcv::gam(gam_formula, data = gam_data, family = stats::gaussian()),
      error = function(e) {
        cli::cli_warn(c(
          "!" = "GAM fit failed for {.val {tgt}}: {conditionMessage(e)}.",
          "i" = "Skipping imputation for this co-analyte."
        ))
        NULL
      }
    )
    if (is.null(gam_fit)) next

    # ── Predict for missing samples ──────────────────────────────────────────
    pred_data <- dplyr::filter(pca_scores, .data$sample_id %in% missing_ids)

    if (return == "point") {
      pred_vals <- exp(
        as.numeric(stats::predict(gam_fit, newdata = pred_data, type = "response"))
      )

      new_rows <- tibble::tibble(
        sample_id    = pred_data$sample_id,
        analyte      = tgt,
        value        = pred_vals,
        detected     = TRUE,
        imputed      = TRUE,
        imputed_kind = "missing"
      ) |>
        dplyr::left_join(sample_meta, by = "sample_id")

    } else {
      # Posterior-predictive draws: beta ~ N(coef, Vp) + residual Gaussian noise.
      # Xp is n_miss x n_coef; t(beta_draws) is n_coef x N → eta_mat is n_miss x N.
      Xp         <- stats::predict(gam_fit, newdata = pred_data, type = "lpmatrix")
      beta_draws <- mgcv::rmvn(N, stats::coef(gam_fit), gam_fit$Vp)  # N x n_coef
      eta_mat    <- Xp %*% t(beta_draws)                        # n_miss x N
      eps_mat    <- matrix(
        stats::rnorm(length(eta_mat), 0, sqrt(gam_fit$sig2)),
        nrow = nrow(eta_mat)
      )
      value_mat  <- exp(eta_mat + eps_mat)                       # n_miss x N

      # as.numeric reads column-major: all n_miss rows of draw 1, then draw 2,
      # etc. rep(..., times=N) and rep(..., each=n_miss) match that layout.
      n_miss <- nrow(pred_data)
      new_rows <- tibble::tibble(
        sample_id    = rep(pred_data$sample_id, times = N),
        draw_id      = rep(seq_len(N), each = n_miss),
        analyte      = tgt,
        value        = as.numeric(value_mat),
        detected     = TRUE,
        imputed      = TRUE,
        imputed_kind = "missing"
      ) |>
        dplyr::left_join(sample_meta, by = "sample_id")
    }

    result <- dplyr::bind_rows(result, new_rows)
  }

  # ── Ensure imputed/imputed_kind columns are populated on all rows ──────────
  if (!"imputed" %in% names(result)) {
    result <- dplyr::mutate(result, imputed = FALSE, imputed_kind = "observed")
  } else {
    result <- dplyr::mutate(result,
      imputed      = dplyr::coalesce(.data$imputed, FALSE),
      imputed_kind = dplyr::coalesce(.data$imputed_kind, "observed")
    )
  }

  dplyr::arrange(result, .data$sample_id, .data$analyte)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Scale-aware log floor for a single numeric vector
#'
#' Half the smallest observed positive value, so a genuine zero maps just
#' below the column's own real values rather than to a single absolute
#' constant shared across analytes of wildly different magnitude.  Falls back
#' to `eps` when the vector has no observed positive values at all.
#' @param x Numeric vector.
#' @param eps Absolute fallback floor (default `1e-9`).
#' @keywords internal
.scale_aware_log_floor <- function(x, eps = 1e-9) {
  pos <- x[!is.na(x) & x > 0]
  if (length(pos) == 0L) eps else min(pos) / 2
}

#' Log10-transform the concentration-like columns of a PCA matrix
#'
#' Applies `log10()` to every column whose name is not in [.PCA_NO_LOG_VARS],
#' leaving pH / temperature / ORP / DO on their native scale.  `NA` cells are
#' preserved (so NIPALS can still handle within-sample missingness).  A genuine
#' zero is floored at half the column's own smallest observed positive value
#' ([.scale_aware_log_floor()]) rather than one absolute constant shared across
#' analytes of wildly different scale.  Both the training PCA
#' (`.prepare_chem_pca()`) and the scoring projection (`.compute_pca_scores()`)
#' call this so the transform is identical on both paths; scoring passes the
#' training-derived `floors` so a cell transforms identically at fit and
#' predict regardless of how many rows the scoring call sees.
#' @param mat Numeric matrix with named columns (samples × variables).
#' @param no_log_vars Column names to leave on their native scale.  Default
#'   [.PCA_NO_LOG_VARS].
#' @param eps Absolute fallback floor for columns with no positive values
#'   (default `1e-9`).
#' @param floors Optional named numeric vector of per-column floors (as
#'   produced by a prior call's `attr(., "floors")`). When `NULL`, floors are
#'   computed from `mat` itself. Supply the training floors at scoring time.
#' @keywords internal
.log_transform_pca <- function(mat, no_log_vars = .PCA_NO_LOG_VARS, eps = 1e-9,
                                floors = NULL) {
  log_cols <- setdiff(colnames(mat), no_log_vars)
  if (length(log_cols) > 0L) {
    sub <- mat[, log_cols, drop = FALSE]
    if (is.null(floors)) {
      floors <- apply(sub, 2, .scale_aware_log_floor, eps = eps)
    } else {
      floors <- floors[log_cols]
      floors[is.na(floors)] <- eps
    }
    floor_mat <- matrix(floors, nrow = nrow(sub), ncol = ncol(sub), byrow = TRUE)
    mat[, log_cols] <- log10(pmax(sub, floor_mat))
  }
  attr(mat, "floors") <- floors
  mat
}

#' Pivot long chemistry rows to a wide per-sample predictor frame
#'
#' Shared by [.prepare_chem_pca()] (training) and [.compute_pca_scores()]
#' (scoring) so predictor construction is identical on both paths.
#' Below-detection cells (`detected == FALSE`) for concentration-like analytes
#' are set to `NA` before collapsing, so a BDL predictor cell is treated as
#' missing (scored by NIPALS from the sample's observed predictors) rather than
#' substituted at its detection limit or DL/2 — this package avoids substituting
#' magic numbers for non-detects. The interval-scale `no_log_vars` (pH,
#' temperature, ORP, DO) are kept as-is when BDL: "below DL" is a concentration
#' idea and they are essentially never BDL. Duplicate `(sample, analyte)`
#' rows collapse via geometric mean for concentration-like analytes (matching
#' how targets are logged before collapsing) and arithmetic mean for
#' `no_log_vars`. An all-`NA` collapse yields `NaN` from `mean()`/`exp(NaN)`;
#' coerced back to `NA` so it reaches NIPALS scoring as missing, not a bogus
#' numeric value (bug B3). Callers without a `detected` column (e.g. the
#' hydro-layer WQ frames in `R/target_model.R`, which carry no BDL concept)
#' are treated as fully detected — no halving.
#' @keywords internal
.pivot_chem_wide <- function(df, wq_vars, no_log_vars) {
  if (!"detected" %in% names(df)) df$detected <- TRUE
  df |>
    dplyr::filter(.data$analyte %in% .env$wq_vars) |>
    dplyr::select("sample_id", "analyte", "value", "detected") |>
    # A below-detection concentration cell carries no trustworthy quantitative
    # value for the predictor PCA. Rather than substitute a magic number (its
    # detection limit, or DL/2 — the substitution hacks this package exists to
    # avoid), drop it to NA and let NIPALS score the sample from its observed
    # predictors (its native missing-cell handling). The interval-scale
    # no_log_vars (pH, temperature, ORP, DO) are kept as-is when BDL: "below
    # DL" is a concentration idea, and they are essentially never BDL anyway.
    # (Quantitatively recovering the "it's low" signal is Route C's job, via a
    # censored predictor treatment; the honest interim is to omit, not invent.)
    dplyr::mutate(
      value = dplyr::if_else(
        .data$detected | .data$analyte %in% .env$no_log_vars,
        .data$value, NA_real_
      )
    ) |>
    dplyr::summarise(
      value = if (dplyr::first(.data$analyte) %in% .env$no_log_vars) {
        mean(.data$value, na.rm = TRUE)
      } else {
        exp(mean(log(pmax(.data$value, .Machine$double.eps)), na.rm = TRUE))
      },
      .by = c("sample_id", "analyte")
    ) |>
    dplyr::mutate(value = dplyr::if_else(is.nan(.data$value), NA_real_, .data$value)) |>
    tidyr::pivot_wider(names_from = "analyte", values_from = "value")
}

#' Fit the unified chemistry PCA on training data
#'
#' This PCA spans the full unified chemistry predictor set (`pca_vars` in
#' `fit_imputation_model()`) — major ions, pH, EC, NH3-N, DOC, nutrients and
#' redox indicators.  Concentration-like variables are `log10`-transformed (see
#' [.log_transform_pca()]) before centring/scaling.  PC score columns are named
#' `PC1`, `PC2`, ….
#' @keywords internal
.prepare_chem_pca <- function(df, wq_vars, min_var_explained = 0.75, max_pcs = 4L,
                              no_log_vars = .PCA_NO_LOG_VARS) {
  # Pivot chemistry vars to wide (one row per sample); missing cells → NA.
  # Below-detection cells are halved (DL/2, parity with the LMF path) and
  # duplicate (sample, analyte) rows collapse on the log scale for
  # concentration-like analytes (geometric mean), matching how targets are
  # collapsed elsewhere. See .pivot_chem_wide().
  wq_wide <- .pivot_chem_wide(df, wq_vars, no_log_vars)

  # Ensure every sample_id is present (even those with no WQ vars)
  all_samples <- tibble::tibble(sample_id = unique(df$sample_id))
  wq_wide     <- dplyr::left_join(all_samples, wq_wide, by = "sample_id")

  wq_matrix   <- as.matrix(dplyr::select(wq_wide, -"sample_id"))

  # Training medians — kept as fallback for columns entirely absent in
  # scoring data.  Per-cell NAs within a column are handled by nipals.
  # Stored on the RAW scale: scoring fills absent columns with these medians
  # and then re-applies the same log transform below.
  train_medians <- apply(wq_matrix, 2, stats::median, na.rm = TRUE)

  # Log10-transform concentration-like variables (everything except
  # pH / temperature / ORP / DO).  Done before centring/scaling so the PCA
  # reflects multiplicative chemical variation rather than being dominated by
  # the highest-magnitude major ions.  NAs are preserved for NIPALS.  Floors
  # are computed from the training matrix here and stored below so scoring
  # re-uses them instead of recomputing from (possibly single-row) new data.
  wq_matrix  <- .log_transform_pca(wq_matrix, no_log_vars = no_log_vars)
  log_floors <- attr(wq_matrix, "floors")

  # Remove zero-variance or all-NA columns
  col_sds   <- apply(wq_matrix, 2, stats::sd, na.rm = TRUE)
  keep_cols <- !is.na(col_sds) & col_sds > 0
  if (!any(keep_cols))
    cli::cli_abort("All PCA variables have zero variance \u2014 cannot fit PCA.")
  if (!all(keep_cols)) {
    dropped <- colnames(wq_matrix)[!keep_cols]
    cli::cli_inform(c("!" = "Chemistry PCA: dropping zero-variance variable{?s}: \\
                              {.val {dropped}}"))
    wq_matrix <- wq_matrix[, keep_cols, drop = FALSE]
  }

  # NIPALS PCA — handles missing cells without prior imputation
  ncomp   <- min(max_pcs, ncol(wq_matrix), nrow(wq_matrix) - 1L)
  pca_fit <- nipals::nipals(wq_matrix, ncomp = ncomp, center = TRUE, scale = TRUE)
  # nipals (>= 1.0) returns per-component proportions in `$R2`; older/other
  # builds expose a cumulative `$R2cum`. Accept either, deriving the cumulative
  # curve from `$R2` when needed.
  cum_var <- if (!is.null(pca_fit$R2cum)) {
    pca_fit$R2cum
  } else {
    cumsum(pca_fit$R2)
  }

  # Determine number of PCs
  n_needed <- which(cum_var >= min_var_explained)[1L]
  if (is.na(n_needed) || n_needed > max_pcs) {
    n_pcs <- min(max_pcs, length(cum_var))
    cli::cli_warn(c(
      "!" = "Chemistry PCA: first {n_pcs} axis/axes explain only \\
             {round(100 * cum_var[n_pcs], 1)}% of variance \\
             (target: {min_var_explained * 100}%).",
      "i" = "Consider adding more {.arg pca_vars} or loosening \\
             {.arg min_var_explained}."
    ))
  } else {
    n_pcs <- max(2L, n_needed)
  }
  n_pcs <- min(n_pcs, length(cum_var))

  pca_obj <- list(
    fit           = pca_fit,
    medians       = train_medians,         # all WQ vars (before zero-var removal)
    active_vars   = colnames(wq_matrix),   # after zero-var removal
    no_log_vars   = no_log_vars,           # linear-scale vars (train/predict parity)
    log_floors    = log_floors,            # training-derived per-column log floors
    n_pcs         = n_pcs,
    var_explained = cum_var[n_pcs]
  )

  # Training scores MUST be produced by the same projection used at prediction
  # time (`.compute_pca_scores()`), otherwise the brms model is trained on one
  # score scale and predicted on another.  `nipals$scores` are NOT equal to the
  # regression projection used at scoring — they differ by a per-component
  # factor (the component eigenvalue) — so copying them here would silently
  # break imputation.  Deriving `pc_scores` via `.compute_pca_scores()` on the
  # training data guarantees train/predict consistency by construction.
  pca_obj$pc_scores <- .compute_pca_scores(df, pca_obj)

  pca_obj
}

#' Project new data onto stored NIPALS PCA axes
#'
#' Handles within-row missing values via NIPALS regression scoring: each
#' component score is estimated from observed variables only, then the residual
#' is deflated before the next component.  Columns entirely absent in `df` (not
#' measured at all, not just BDL) are filled with training medians.
#' @keywords internal
.compute_pca_scores <- function(df, pca_obj) {
  wq_vars     <- pca_obj$active_vars  # columns used in training
  no_log_vars <- pca_obj$no_log_vars %||% .PCA_NO_LOG_VARS

  # Same BDL half-DL substitution and log-scale duplicate collapse as
  # training (.pivot_chem_wide()), for train/predict parity.
  wq_wide <- .pivot_chem_wide(df, wq_vars, no_log_vars)

  all_samples <- tibble::tibble(sample_id = unique(df$sample_id))
  wq_wide     <- dplyr::left_join(all_samples, wq_wide, by = "sample_id")

  # Add entirely-absent training columns using training medians.
  # Per-cell NAs within a column are handled by the NIPALS scoring below.
  for (v in wq_vars) {
    if (!v %in% names(wq_wide)) wq_wide[[v]] <- pca_obj$medians[[v]]
  }

  wq_mat <- as.matrix(dplyr::select(wq_wide, dplyr::all_of(wq_vars)))

  # Fill columns that are entirely NA (not measured in this dataset at all)
  for (j in seq_len(ncol(wq_mat))) {
    if (all(is.na(wq_mat[, j]))) {
      col_nm       <- colnames(wq_mat)[j]
      wq_mat[, j]  <- pca_obj$medians[[col_nm]]
    }
  }

  # Apply the same log10 transform used at training time.  Done AFTER the
  # raw-scale median fills above so filled values are transformed identically
  # to how the training medians were, keeping centre/scale parameters valid.
  # Use the training-time no_log_vars AND the training-derived per-column
  # floors for train/predict parity (older pca_obj without these fields falls
  # back to the package default / recomputing from this call's data).
  wq_mat <- .log_transform_pca(
    wq_mat, no_log_vars = no_log_vars, floors = pca_obj$log_floors
  )

  # Centre and scale using training parameters stored in the nipals object
  wq_scaled <- sweep(wq_mat,    2, pca_obj$fit$center, "-")
  wq_scaled <- sweep(wq_scaled, 2, pca_obj$fit$scale,  "/")

  # NIPALS regression scoring: per-row, per-component
  loadings   <- pca_obj$fit$loadings[, seq_len(pca_obj$n_pcs), drop = FALSE]
  scores_mat <- t(apply(wq_scaled, 1, .nipals_score_row,
                        loadings = loadings, n_pcs = pca_obj$n_pcs))
  colnames(scores_mat) <- paste0("PC", seq_len(pca_obj$n_pcs))

  tibble::as_tibble(scores_mat) |>
    dplyr::mutate(sample_id = wq_wide$sample_id)
}

#' NIPALS regression scoring for one centred/scaled observation
#'
#' Computes PC scores for a single row by projecting onto each loading vector
#' using only observed (non-NA) elements, then deflating the residual before
#' moving to the next component.  For a fully observed row this is identical to
#' the standard `x %*% loadings` projection.  For a row with missing values it
#' correctly down-weights the loading vectors to the observed subspace — without
#' the bias that zero/median imputation introduces.
#'
#' @param x Numeric vector (centred + scaled); `NA` marks missing variables.
#' @param loadings p × K loading matrix (unit-normalised columns from nipals).
#' @param n_pcs Number of components to score.
#' @keywords internal
.nipals_score_row <- function(x, loadings, n_pcs) {
  scores  <- numeric(n_pcs)
  x_resid <- x
  for (k in seq_len(n_pcs)) {
    lk  <- loadings[, k]
    obs <- !is.na(x_resid)
    if (any(obs)) {
      # Regression onto observed sub-vector of the loading; renormalise by
      # sum(lk[obs]^2) because the loading is unit-normalised over *all* p vars
      scores[k] <- sum(x_resid[obs] * lk[obs]) / sum(lk[obs]^2)
    }
    # Deflate: remove this component's contribution from observed elements
    x_resid[obs] <- x_resid[obs] - scores[k] * lk[obs]
  }
  scores
}

#' Guard against make.names() collisions in a set of analyte names
#'
#' The fitting/prediction path maps original analyte names to syntactically
#' valid R names via `make.names()` and relies on that mapping being a
#' bijection. If two distinct analytes collide under `make.names()` (e.g.
#' `"Cr-6"` and `"Cr.6"` both become `"Cr.6"`), a value-equality inverse
#' lookup can silently return the wrong original name. This guard errors
#' early and informatively instead.
#'
#' @param targets Character vector of original analyte names.
#' @keywords internal
.assert_safe_analyte_names <- function(targets) {
  safe <- make.names(targets)
  dupe_safe <- unique(safe[duplicated(safe)])
  if (length(dupe_safe) > 0) {
    colliding <- targets[safe %in% dupe_safe]
    cli::cli_abort(c(
      "Safe name collision: distinct analyte names map to the same
       {.code make.names()} value, so the safe-name mapping is not unique.",
      "x" = "Colliding original name{?s}: {.val {colliding}}",
      "i" = "Colliding safe name{?s}: {.val {dupe_safe}}"
    ))
  }
  invisible(TRUE)
}

#' Fit a single brms group model (metals or organics)
#' @keywords internal
.fit_group_model <- function(df, target_analytes, pca_obj,
                              family, iter, warmup, chains, cores,
                              impute_method = "rescor_mi",
                              group_name = "group", ...) {
  .assert_safe_analyte_names(target_analytes)
  eps_log <- 1e-9

  safe_analytes <- stats::setNames(make.names(target_analytes), target_analytes)
  # safe_analytes: names = safe R names, values = original names
  safe_vec <- unname(safe_analytes)
  pc_cols  <- paste0("PC", seq_len(pca_obj$n_pcs))
  rhs      <- paste(paste0("s(", pc_cols, ")"), collapse = " + ")
  pc_wide  <- dplyr::select(pca_obj$pc_scores, "sample_id", dplyr::all_of(pc_cols))

  # Per-analyte log floor (half the analyte's own smallest positive value)
  # rather than one constant shared across analytes of wildly different scale
  # (finding 6). Computed once at fit time and reused at prediction time
  # (.predict_and_merge()) for train/predict parity.
  log_floors <- df |>
    dplyr::filter(.data$analyte %in% .env$target_analytes) |>
    dplyr::summarise(
      .floor = .scale_aware_log_floor(.data$value, eps = eps_log),
      .by = "analyte"
    ) |>
    tibble::deframe()

  # One row per (sample, analyte) with a safe column name and log value.
  base <- df |>
    dplyr::filter(.data$analyte %in% .env$target_analytes) |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      safe = unname(safe_analytes[.data$analyte]),
      lv   = log(pmax(.data$value, log_floors[.data$analyte]))
    )

  if (impute_method == "factor") {
    # Route C: low-rank censored factor model (dev/plan-route-c.md). Its
    # two-stage fit (per-analyte GAM mean + Stage-2 Stan factor model) is
    # structurally different from the brms branches below, so it is handled
    # entirely by .fit_group_model_factor() and returns directly.
    return(.fit_group_model_factor(
      target_analytes = target_analytes,
      safe_analytes   = safe_analytes,
      base            = base,
      pc_wide         = pc_wide,
      pc_cols         = pc_cols,
      log_floors      = log_floors,
      iter            = iter,
      warmup          = warmup,
      chains          = chains,
      cores           = cores,
      group_name      = group_name,
      ...
    ))
  }

  if (impute_method == "rescor_mi") {
    # Residual correlation + mi() for BDL/missing. BDL and missing are NA and
    # imputed; the post-hoc DL cap is applied in impute_chemistry().
    target_wide <- base |>
      dplyr::mutate(log_value = dplyr::if_else(.data$detected, .data$lv, NA_real_)) |>
      dplyr::select("sample_id", "safe", "log_value") |>
      tidyr::pivot_wider(names_from = "safe", values_from = "log_value")
    model_df <- dplyr::left_join(pc_wide, target_wide, by = "sample_id")
    bf_list <- purrr::map(safe_vec, function(s)
      brms::bf(stats::as.formula(paste0(s, " | mi() ~ ", rhs))))
    brms_formula <- do.call(brms::mvbf, c(bf_list, list(rescor = TRUE)))

  } else if (impute_method == "cens") {
    # Wide + subset(): left-censor BDL at its detection limit, rescor = FALSE,
    # responses independent (no cross-analyte coupling). subset() lets each
    # analyte use only the samples that measured it.
    resp_wide <- base |>
      dplyr::select("sample_id", "safe", "lv") |>
      tidyr::pivot_wider(names_from = "safe", values_from = "lv")
    cens_wide <- base |>
      dplyr::mutate(cf = dplyr::if_else(.data$detected, "none", "left")) |>
      dplyr::select("sample_id", "safe", "cf") |>
      tidyr::pivot_wider(names_from = "safe", values_from = "cf",
                         names_glue = "cens_{safe}")
    model_df <- pc_wide |>
      dplyr::left_join(resp_wide, by = "sample_id") |>
      dplyr::left_join(cens_wide, by = "sample_id")
    for (s in safe_vec) {
      model_df[[paste0("sub_", s)]]  <- !is.na(model_df[[s]])
      model_df[[s]]                  <- dplyr::coalesce(model_df[[s]], 0)
      model_df[[paste0("cens_", s)]] <- dplyr::coalesce(model_df[[paste0("cens_", s)]], "none")
    }
    bf_list <- purrr::map(safe_vec, function(s)
      brms::bf(stats::as.formula(sprintf(
        "%s | cens(cens_%s) + subset(sub_%s) ~ %s", s, s, s, rhs))))
    brms_formula <- do.call(brms::mvbf, c(bf_list, list(rescor = FALSE)))

  } else {
    # cens_factor: LONG-format shared-latent-factor model. One univariate
    # censored model with `(1 | sample_id)` shared across analytes provides
    # genuine cross-analyte coupling that is well-identified — each sample
    # contributes several analyte observations to pin down its single latent
    # value, so an observed metal informs the unobserved/BDL ones at that sample.
    # (The previous wide + `(1 |q| sample_id)` form silently dropped the
    # cross-response correlation under subset(), giving the cost of extra
    # parameters with none of the coupling.) `analyte` enters as a factor: a
    # per-analyte mean (`0 + safe`), per-analyte PC smooths (`by = safe`), and a
    # per-analyte residual SD (`sigma ~ 0 + safe`).
    long <- base |>
      dplyr::transmute(
        sample_id = .data$sample_id,
        safe      = factor(.data$safe, levels = safe_vec),
        lv        = .data$lv,
        cf        = dplyr::if_else(.data$detected, "none", "left")
      ) |>
      dplyr::left_join(pc_wide, by = "sample_id")
    long <- long[stats::complete.cases(long[, pc_cols, drop = FALSE]), , drop = FALSE]
    rhs_by <- paste(sprintf("s(%s, by = safe)", pc_cols), collapse = " + ")
    brms_formula <- brms::bf(
      stats::as.formula(paste0("lv | cens(cf) ~ 0 + safe + ", rhs_by,
                               " + (1 | sample_id)")),
      stats::as.formula("sigma ~ 0 + safe")
    )
    model_df <- long
  }

  n_units <- dplyr::n_distinct(model_df$sample_id)
  cli::cli_inform(c(
    "i" = "brms {group_name} ({impute_method}): {length(target_analytes)} \\
           analyte{?s} \u00d7 {n_units} sample{?s}. This may take several minutes."
  ))

  brm_args <- list(
    formula = brms_formula,
    data    = model_df,
    family  = family,
    iter    = iter,
    warmup  = warmup,
    chains  = chains,
    cores   = cores,
    ...
  )
  # Default to cmdstanr (cached binaries + faster sampling, statistically
  # equivalent to rstan) when available; an explicit backend in `...` wins.
  if (is.null(brm_args$backend)) brm_args$backend <- .brms_backend()
  # The cens_factor hierarchical model benefits from a higher adapt_delta to
  # clear the residual divergences from the per-sample factor's mild funnel.
  # Default it here; a user-supplied `control` (via ...) takes precedence.
  if (impute_method == "cens_factor" && is.null(brm_args$control)) {
    brm_args$control <- list(adapt_delta = 0.95)
  }
  # rescor_mi: the mi() + full residual-correlation geometry is funnel-prone.
  # adapt_delta = 0.95 + an lkj(2) prior on the residual correlation cut
  # divergences and were retained because rescor_mi still wins recovery by a wide
  # margin (3-seed B.S01 hold-out: RMSE 0.30 vs 0.75 cens_factor / 0.92 cens,
  # near-zero bias). HOWEVER the geometry remains hard: the sampler saturates the
  # NUTS tree depth (72-99% of transitions) and trips E-BFMI, and worst-case R-hat
  # reached ~1.6 on a hard mask. The point estimate is reliable; the DRAWS are
  # not guaranteed converged — production runs wanting trustworthy uncertainty
  # should check rhat() and consider raising max_treedepth / iter, or use
  # impute_method = "cens_factor" (best-calibrated intervals; see
  # vignette("imputation") and issue #59). All overridable via `...`.
  if (impute_method == "rescor_mi") {
    if (is.null(brm_args$control)) brm_args$control <- list(adapt_delta = 0.95)
    if (is.null(brm_args$prior))
      brm_args$prior <- brms::set_prior("lkj(2)", class = "rescor")
  }

  fit <- tryCatch(
    do.call(brms::brm, brm_args),
    error = function(e) {
      cli::cli_abort(c(
        "brms model fitting failed.",
        "x" = "{conditionMessage(e)}",
        "i" = "If this is a Stan compilation error, check your Stan toolchain \\
               ({.url https://paul-buerkner.github.io/brms/})."
      ))
    }
  )

  list(
    fit             = fit,
    analytes        = target_analytes,   # original names
    safe_names      = safe_analytes,     # names=safe, values=original
    pc_cols         = pc_cols,
    wide_sample_ids = unique(model_df$sample_id),
    impute_method   = impute_method,
    log_floors      = log_floors         # per-analyte log floor, train/predict parity
  )
}

#' Long-format prediction for the cens_factor shared-latent-factor model
#'
#' Builds one row per (eligible sample × analyte), predicts the latent log
#' concentration from the univariate model (the per-sample `(1 | sample_id)`
#' factor couples analytes), and returns the same `pm_long` shape the wide path
#' produces so the merge step is identical.
#' @keywords internal
.predict_factor_long <- function(group, pc_wide, return, ndraws, batch_size) {
  safe_analytes <- group$safe_names               # names = safe, values = orig
  safe_levels   <- unname(safe_analytes)
  orig_of       <- stats::setNames(names(safe_analytes), unname(safe_analytes))

  nd <- tidyr::expand_grid(sample_id = pc_wide$sample_id, safe = safe_levels) |>
    dplyr::left_join(pc_wide, by = "sample_id")
  nd$safe <- factor(nd$safe, levels = safe_levels)
  nd$cf   <- "none"   # placeholder; ignored by epred/predict
  nd$lv   <- 0        # response placeholder

  predict_rows <- function(ndi) {
    analyte <- unname(orig_of[as.character(ndi$safe)])
    if (return == "point") {
      ep <- brms::posterior_epred(group$fit, newdata = ndi, allow_new_levels = TRUE,
                                  sample_new_levels = "gaussian", ndraws = ndraws)
      tibble::tibble(sample_id = ndi$sample_id, analyte = analyte,
                     .post_mean = exp(colMeans(ep)))
    } else {
      pp  <- brms::posterior_predict(group$fit, newdata = ndi, allow_new_levels = TRUE,
                                     sample_new_levels = "gaussian", ndraws = ndraws)
      ndr <- nrow(pp)
      tibble::tibble(
        sample_id   = rep(ndi$sample_id, each = ndr),
        analyte     = rep(analyte, each = ndr),
        draw_id     = rep(seq_len(ndr), times = nrow(ndi)),
        .post_value = exp(as.vector(pp)))
    }
  }

  n_new <- nrow(nd)
  # batch_size is in samples; each sample spans length(safe_levels) rows.
  bs <- if (is.null(batch_size)) n_new
        else max(length(safe_levels), as.integer(batch_size) * length(safe_levels))
  batches <- split(seq_len(n_new), ceiling(seq_len(n_new) / bs))
  tryCatch(
    purrr::map_dfr(batches, function(idx) predict_rows(nd[idx, , drop = FALSE])),
    error = function(e) cli::cli_abort(c(
      "brms prediction failed during imputation (cens_factor).",
      "x" = "{conditionMessage(e)}"
    ))
  )
}

#' Predict and merge imputed values for one analyte group
#' @keywords internal
.predict_and_merge <- function(df, group, pca_scores, eligible_ids, return,
                                ndraws = NULL, batch_size = NULL) {
  eps_log    <- 1e-9
  log_floors <- group$log_floors  # per-analyte, from fit time (train/predict parity)

  target_analytes <- group$analytes
  safe_analytes   <- group$safe_names   # names=safe, values=original
  pc_cols         <- group$pc_cols

  # ── Build wide prediction df for eligible samples ─────────────────────────
  df_eligible <- dplyr::filter(df, .data$sample_id %in% .env$eligible_ids)
  if (nrow(df_eligible) == 0L) return(df)

  # PC scores for eligible samples (needed by every method).
  pc_wide <- dplyr::filter(pca_scores, .data$sample_id %in% .env$eligible_ids) |>
    dplyr::select("sample_id", dplyr::all_of(pc_cols))

  if (!is.null(group$impute_method) && group$impute_method == "cens_factor") {
    # ── Long-format shared-latent-factor prediction ──────────────────────────
    # Predict every (eligible sample × analyte) cell from the univariate model;
    # the merge below overlays only the BDL/missing cells.
    pm_long <- .predict_factor_long(group, pc_wide, return, ndraws, batch_size)

  } else {
    # ── Wide newdata + per-method prediction (rescor_mi / cens) ──────────────
    # Targets wide (log for detected, NA for BDL/missing → to be imputed)
    target_wide <- df_eligible |>
      dplyr::filter(.data$analyte %in% .env$target_analytes) |>
      dplyr::select("sample_id", "analyte", "value", "detected") |>
      dplyr::mutate(
        log_value = dplyr::if_else(
          .data$detected,
          log(pmax(.data$value, dplyr::coalesce(log_floors[.data$analyte], eps_log))),
          NA_real_
        )
      ) |>
      dplyr::summarise(
        log_value = mean(.data$log_value, na.rm = TRUE),
        .by        = c("sample_id", "analyte")
      ) |>
      tidyr::pivot_wider(names_from = "analyte", values_from = "log_value") |>
      dplyr::rename(dplyr::any_of(
        stats::setNames(names(safe_analytes), unname(safe_analytes))
      ))

    # Ensure all training analyte columns are present (add NA if absent)
    for (s in names(safe_analytes)) {
      if (!s %in% names(target_wide)) target_wide[[s]] <- NA_real_
    }

    wide_new <- pc_wide |>
      dplyr::left_join(target_wide, by = "sample_id")

    # cens models reference cens_<safe> / sub_<safe> columns in the formula. For
    # prediction set them so every target cell is predicted (cens "none", subset
    # TRUE); the response placeholder is ignored.
    if (!is.null(group$impute_method) && group$impute_method != "rescor_mi") {
      for (s in unname(safe_analytes)) {
        if (!s %in% names(wide_new)) wide_new[[s]] <- 0
        wide_new[[s]]                 <- dplyr::coalesce(wide_new[[s]], 0)
        wide_new[[paste0("cens_", s)]] <- "none"
        wide_new[[paste0("sub_", s)]]  <- TRUE
      }
    }

    # ── Posterior predictions (batched + optional draw subsampling) ─────────
    # Predicting all eligible samples at once can exhaust memory for rescor_mi
    # (mi() materialises every imputed cell as a latent parameter). Predict in
    # row-batches, and optionally use only `ndraws` posterior draws, to bound
    # peak memory.
    cens_method <- !is.null(group$impute_method) &&
      group$impute_method != "rescor_mi"

    predict_chunk <- function(wn) {
      if (cens_method) {
        # subset() models must be predicted one response at a time.
        purrr::map_dfr(unname(safe_analytes), function(s) {
          orig <- names(safe_analytes)[safe_analytes == s]
          if (return == "point") {
            ep <- brms::posterior_epred(group$fit, newdata = wn, resp = s,
                                        allow_new_levels = TRUE, ndraws = ndraws)
            tibble::tibble(sample_id = wn$sample_id, analyte = orig,
                           .post_mean = exp(colMeans(ep)))
          } else {
            pp <- brms::posterior_predict(group$fit, newdata = wn, resp = s,
                                          allow_new_levels = TRUE, ndraws = ndraws)
            nd <- nrow(pp)
            tibble::tibble(
              sample_id   = rep(wn$sample_id, each = nd),
              analyte     = orig,
              draw_id     = rep(seq_len(nd), times = ncol(pp)),
              .post_value = exp(as.vector(pp)))
          }
        })
      } else if (return == "point") {
        ep <- brms::posterior_epred(group$fit, newdata = wn,
                                    allow_new_levels = TRUE, ndraws = ndraws)
        .reshape_posterior_means(ep, wn$sample_id, safe_analytes)
      } else {
        pp <- brms::posterior_predict(group$fit, newdata = wn,
                                      allow_new_levels = TRUE, ndraws = ndraws)
        .reshape_posterior_draws(pp, wn$sample_id, safe_analytes)
      }
    }

    n_new <- nrow(wide_new)
    bs    <- if (is.null(batch_size)) n_new else max(1L, as.integer(batch_size))
    batches <- split(seq_len(n_new), ceiling(seq_len(n_new) / bs))
    pm_long <- tryCatch(
      purrr::map_dfr(batches, function(idx) predict_chunk(wide_new[idx, , drop = FALSE])),
      error = function(e) cli::cli_abort(c(
        "brms prediction failed during imputation.",
        "x" = "{conditionMessage(e)}"
      ))
    )
  }

  # ── Tag imputation kind for each (sample, analyte) ─────────────────────────
  impute_kind <- df_eligible |>
    dplyr::filter(.data$analyte %in% .env$target_analytes) |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      sample_id,
      analyte,
      .imputed_kind = dplyr::if_else(.data$detected, "observed", "censored_left")
    )

  # Samples missing the analyte entirely → "missing"
  all_combos <- tidyr::expand_grid(
    sample_id = unique(df_eligible$sample_id),
    analyte   = target_analytes
  )
  impute_kind <- dplyr::left_join(all_combos, impute_kind,
                                  by = c("sample_id", "analyte")) |>
    dplyr::mutate(
      .imputed_kind = dplyr::if_else(is.na(.data$.imputed_kind),
                                     "missing", .data$.imputed_kind)
    )

  # ── Merge posterior values back into df ───────────────────────────────────
  # For eligible samples: replace BDL / missing target values with posterior
  # means; add imputed/imputed_kind columns.  Non-target rows are unchanged.

  val_col <- if (return == "point") ".post_mean" else ".post_value"

  if (!"imputed" %in% names(df)) {
    df <- dplyr::mutate(df, imputed = FALSE, imputed_kind = "observed")
  }

  target_rows_eligible <- df |>
    dplyr::filter(
      .data$sample_id %in% .env$eligible_ids,
      .data$analyte   %in% .env$target_analytes
    ) |>
    dplyr::left_join(impute_kind, by = c("sample_id", "analyte")) |>
    dplyr::mutate(
      imputed      = .data$.imputed_kind != "observed",
      imputed_kind = .data$.imputed_kind
    ) |>
    dplyr::select(-".imputed_kind")

  # Overlay posterior means onto imputed rows
  imputed_rows   <- dplyr::filter(target_rows_eligible, .data$imputed)
  observed_rows  <- dplyr::filter(target_rows_eligible, !.data$imputed)

  join_cols <- if (return == "point") {
    c("sample_id", "analyte", val_col)
  } else {
    c("sample_id", "analyte", "draw_id", val_col)
  }

  imputed_filled <- dplyr::left_join(
    imputed_rows,
    dplyr::select(pm_long, dplyr::all_of(join_cols)),
    by = c("sample_id", "analyte")
  ) |>
    dplyr::mutate(
      value    = dplyr::coalesce(.data[[val_col]], .data$value),
      detected = TRUE
    ) |>
    dplyr::select(-dplyr::all_of(val_col))

  # ── Fabricate rows for entirely-absent target cells ───────────────────────
  # The overlay above can only fill cells that already have a row. An eligible
  # sample missing a target analyte *entirely* is classified "missing" in
  # `impute_kind`, and its prediction already sits in `pm_long` — emit a new row
  # for it so the group's model completes the panel (e.g. a Zn-only sample gains
  # a model-anchored Cu row). Hurdle-failing samples never reach `impute_kind`
  # (they are not in `eligible_ids`), so they are never fabricated.
  missing_combos <- impute_kind |>
    dplyr::filter(.data$.imputed_kind == "missing") |>
    dplyr::select("sample_id", "analyte")

  new_rows <- NULL
  if (nrow(missing_combos) > 0L) {
    # Per-sample carrier metadata: every column except the per-cell ones, which
    # we set explicitly below (and `draw_id` / `val_col`, supplied by pm_long).
    per_cell_cols <- c(
      "analyte", "value", "detected", "imputed",
      "imputed_kind", "draw_id", val_col
    )
    carrier <- df |>
      dplyr::filter(.data$sample_id %in% missing_combos$sample_id) |>
      dplyr::group_by(.data$sample_id) |>
      dplyr::slice(1L) |>
      dplyr::ungroup() |>
      dplyr::select(-dplyr::any_of(per_cell_cols))

    new_rows <- missing_combos |>
      dplyr::left_join(carrier, by = "sample_id") |>
      dplyr::left_join(
        dplyr::select(pm_long, dplyr::all_of(join_cols)),
        by = c("sample_id", "analyte")
      ) |>
      # Drop any combo without a prediction (e.g. a sample lacking PC scores):
      # fabricating an NA-valued row would feed garbage to downstream models.
      dplyr::filter(!is.na(.data[[val_col]])) |>
      dplyr::mutate(
        value        = .data[[val_col]],
        detected     = TRUE,
        imputed      = TRUE,
        imputed_kind = "missing"
      ) |>
      dplyr::select(-dplyr::all_of(val_col))
  }

  # Rows for non-eligible samples (failed hurdle): keep unchanged
  non_eligible_target_rows <- df |>
    dplyr::filter(
      !(.data$sample_id %in% .env$eligible_ids),
      .data$analyte %in% .env$target_analytes
    )

  # Non-target rows
  non_target_rows <- df |>
    dplyr::filter(!(.data$analyte %in% .env$target_analytes))

  dplyr::bind_rows(
    non_target_rows,
    observed_rows,
    imputed_filled,
    new_rows,
    non_eligible_target_rows
  ) |>
    dplyr::arrange(.data$sample_id, .data$analyte)
}


# ── Posterior reshape helpers (unchanged from original) ──────────────────────

#' @keywords internal
.reshape_posterior_means <- function(epred_draws, sample_ids, safe_analytes) {
  if (is.matrix(epred_draws)) {
    arr_list <- stats::setNames(list(epred_draws), unname(safe_analytes)[1L])
  } else {
    resp_nms <- dimnames(epred_draws)[[3L]]
    arr_list <- stats::setNames(
      lapply(seq_along(resp_nms), function(i) epred_draws[, , i]),
      resp_nms
    )
  }

  purrr::map2_dfr(arr_list, names(arr_list), function(mat, safe_nm) {
    orig_nm <- names(safe_analytes)[safe_analytes == safe_nm]
    if (length(orig_nm) == 0L) orig_nm <- safe_nm
    tibble::tibble(
      sample_id  = sample_ids,
      analyte    = orig_nm,
      .post_mean = exp(colMeans(mat))
    )
  })
}

#' @keywords internal
.reshape_posterior_draws <- function(post_draws, sample_ids, safe_analytes) {
  n_draws  <- dim(post_draws)[1L]
  resp_nms <- dimnames(post_draws)[[3L]]

  purrr::map2_dfr(seq_along(resp_nms), resp_nms, function(ri, safe_nm) {
    orig_nm <- names(safe_analytes)[safe_analytes == safe_nm]
    if (length(orig_nm) == 0L) orig_nm <- safe_nm
    mat <- post_draws[, , ri]
    tibble::tibble(
      sample_id   = rep(sample_ids, each = n_draws),
      analyte     = orig_nm,
      draw_id     = rep(seq_len(n_draws), times = length(sample_ids)),
      .post_value = exp(as.vector(t(mat)))
    )
  })
}


# ── BDL cap check + audit summary ─────────────────────────────────────────────

#' Inspect detection-limit cap activations from `impute_chemistry()`
#'
#' An imputed below-detection (BDL) cell must not exceed its detection limit
#' (DL).  The posterior prediction is not itself constrained below the DL, so
#' [impute_chemistry()] caps any imputed BDL cell whose estimate came out above
#' the limit (`bdl_cap = TRUE`).  Frequent capping signals tension between the
#' modelled chemistry and the reported limits, so the cells that triggered the
#' cap are worth auditing rather than trusting blindly.
#'
#' `impute_chemistry()` attaches a per-cell audit summary to its result as the
#' `"bdl_cap_summary"` attribute; this accessor returns it.  Because plain
#' attributes are dropped by most \pkg{dplyr} verbs, call this on the frame
#' **as returned by `impute_chemistry()`**, before further wrangling.
#'
#' @param x A data frame returned by [impute_chemistry()].
#' @return A tibble with one row per (`sample_id`, `analyte`) cell that exceeded
#'   its detection limit, with columns `detection_limit`, `n_rows` (rows over
#'   the DL — one per draw when `return = "draws"`), `max_imputed`, `max_ratio`
#'   (`max_imputed / detection_limit`) and `capped` (whether the cap was
#'   applied).  The tibble carries the overall **`fire_rate`** (fraction of
#'   capable BDL cells that were clipped) and **`n_bdl_cells`** as attributes,
#'   also reported via a message.  Returns `NULL` invisibly when no cell
#'   exceeded its DL.
#' @seealso [impute_chemistry()]
#' @export
bdl_cap_summary <- function(x) {
  s <- attr(x, "bdl_cap_summary", exact = TRUE)
  if (is.null(s)) {
    cli::cli_inform(c("v" = "No detection-limit cap activations recorded."))
    return(invisible(NULL))
  }
  fr <- attr(s, "fire_rate", exact = TRUE)
  nb <- attr(s, "n_bdl_cells", exact = TRUE)
  if (!is.null(fr) && !is.null(nb)) {
    cli::cli_inform(c(
      "i" = "BDL-cap fire-rate: {nrow(s)}/{nb} capable BDL cell{?s} \\
             ({round(100 * fr)}%) exceeded the DL{if (fr >= 0.5) ' - high; consider impute_method = \"cens_factor\"' else ''}."
    ))
  }
  s
}

#' @keywords internal
.check_bdl_imputed <- function(result, dl_tbl, cap = TRUE) {
  if (nrow(dl_tbl) == 0L) return(result)

  # dl_tbl is one row per (sample_id, analyte); join on the unique key so a
  # multi-draw `result` (return = "draws") is never duplicated by the join.
  dl_join <- dplyr::distinct(dl_tbl, .data$sample_id, .data$analyte,
                             .keep_all = TRUE)
  joined <- dplyr::left_join(result, dl_join, by = c("sample_id", "analyte"))
  if (!"detection_limit" %in% names(joined)) return(result)

  # BDL cells that *could* be capped (censored-left with a known DL) — the
  # denominator for the fire-rate. Counted per distinct (sample, analyte) cell
  # so draws-mode (many rows per cell) is not double-counted.
  bdl_at_risk <- joined$imputed_kind == "censored_left" &
                 !is.na(joined$detection_limit)
  n_bdl_cells <- dplyr::n_distinct(
    joined$sample_id[bdl_at_risk], joined$analyte[bdl_at_risk])

  exceed <- bdl_at_risk & joined$value > joined$detection_limit

  if (!any(exceed)) return(result)

  # Per-(sample, analyte) audit summary, computed on the *pre-cap* values so it
  # records the magnitude of each exceedance even when capping is applied.
  summary_tbl <- joined[exceed, , drop = FALSE] |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::summarise(
      detection_limit = dplyr::first(.data$detection_limit),
      n_rows          = dplyr::n(),
      max_imputed     = max(.data$value),
      max_ratio       = max(.data$value / .data$detection_limit),
      .groups = "drop"
    ) |>
    dplyr::mutate(capped = cap) |>
    dplyr::arrange(dplyr::desc(.data$max_ratio))

  # Per-analyte rollup for an auditable (not just aggregate) warning.
  by_analyte <- summary_tbl |>
    dplyr::group_by(.data$analyte) |>
    dplyr::summarise(cells     = dplyr::n(),
                     max_ratio = max(.data$max_ratio),
                     .groups   = "drop") |>
    dplyr::arrange(dplyr::desc(.data$max_ratio))

  analyte_lines <- stats::setNames(
    sprintf("%s: %d cell%s, up to %.1f\u00d7 DL",
            by_analyte$analyte, by_analyte$cells,
            ifelse(by_analyte$cells == 1L, "", "s"), by_analyte$max_ratio),
    rep("*", nrow(by_analyte))
  )

  n_cells   <- nrow(summary_tbl)
  fire_rate <- n_cells / max(n_bdl_cells, 1L)
  attr(summary_tbl, "fire_rate")   <- fire_rate
  attr(summary_tbl, "n_bdl_cells") <- n_bdl_cells
  # A high fire-rate is diagnostic: the model is over-predicting BDL cells. This
  # is characteristic of `rescor_mi` (mi() + post-hoc cap rather than respecting
  # the censoring bound); proper left-censoring (`cens` / `cens_factor`) barely
  # fires the cap. (#59 item 4 / #63.)
  high <- fire_rate >= 0.5
  cli::cli_warn(c(
    "!" = "{n_cells} of {n_bdl_cells} imputed below-detection cell{?s} \\
           ({round(100 * fire_rate)}%) exceeded the detection limit.",
    analyte_lines,
    "i" = "The posterior is not constrained below the DL, so some BDL cells \\
           came out above it.",
    if (high) c("!" = "A high cap fire-rate signals systematic over-prediction \\
                       of BDL cells; consider {.code impute_method = \"cens_factor\"} \\
                       (proper left-censoring) on BDL-heavy chemistry."),
    if (cap) c("i" = "Capped at DL ({.code bdl_cap = TRUE}).")
    else     c("i" = "NOT capped ({.code bdl_cap = FALSE})."),
    "i" = "Per-cell detail: {.run bdl_cap_summary(x)} on the returned frame."
  ))

  if (cap) joined$value[exceed] <- joined$detection_limit[exceed]

  out <- dplyr::select(joined, -"detection_limit")
  attr(out, "bdl_cap_summary") <- summary_tbl
  out
}
