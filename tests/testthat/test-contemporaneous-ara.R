## Contemporaneous ARA through mspaf_pipeline(), and tier-1 direct matching
## ============================================================================
##
## CONTRACT.
##
## PART 1 -- `ara` argument on mspaf_pipeline().
##
## Today the only ARA route through mspaf_pipeline() is the STATIC one:
## `daily_args$reference = <raw reference frame>`, which add_mspaf() turns into
## a per-analyte geometric mean (prepare_reference()).  mspaf_daily() and
## add_mspaf() already accept a fitted `reference_model` as `reference`
## (contemporaneous ARA: ref_norm(t) is the reference model's resolution for
## day t, see .resolve_ref_norm()), but the pipeline fits its reference model
## internally after daily_args is assembled, so a caller cannot ask for it
## without fitting the model twice.
##
## P1. New argument `ara = c("static", "contemporaneous", "none")`.
##     * `ara` NOT supplied: behaviour is exactly today's -- static ARA iff
##       `daily_args$reference` is given, no ARA otherwise.
##     * `ara = "static"`: requires `daily_args$reference` (else classed error
##       `hydroSense_ara_static_needs_reference`); output equals today's
##       static run (ara not supplied, same daily_args$reference).
##     * `ara = "none"`: no ARA; equals a run with no daily_args$reference.
##     * `ara = "contemporaneous"`: the pipeline's OWN fitted reference_model
##       is passed to mspaf_daily() as `reference`.  fit_reference_model() is
##       called exactly once.  Together with `daily_args$reference` it is a
##       classed error `hydroSense_ara_conflict`.  (The tests run it with
##       also_total = TRUE, see P3, so one pipeline run serves P1-P4.)
##     * an unknown `ara` value is an error.
## P2. ara = "contemporaneous" equals a manual replication (tol 1e-10):
##     ref <- fit_reference_model(reference, hydro = hydro, <reference_args>);
##     daily <- mspaf_daily(target, reference_model = ref, reference = ref,
##                          interpolation = "model", <daily_args>);
##     chronic <- time_weighted_aggregate() of the daily mspaf line exactly as
##     .chronic_from_daily() does in point mode.  Both the daily `mspaf` series
##     (attr "daily") and the chronic `value` must match.  The daily
##     ara_summary carries non-NA ref_tier (a temporal reference was used).
## P3. also_total = TRUE with ara = "contemporaneous": `value` is unchanged
##     and `value_total` equals `value` of the ara = "none" run.
## P4. Property: when the target point IS the reference site, contemporaneous
##     ARA removes most of the signal, so its chronic msPAF has a markedly lower
##     median than static ARA's (median(contemp) < 0.5 * median(static)).
##     Pointwise <= is NOT required.
##
## PART 2 -- tier-1 direct matching in .resolve_ref_norm_instant() never fires.
##
## Verified mechanism (R/reference_model.R, origin/main 8fe6de3):
##   * fit_reference_model() stores, for each analyte, `obs$hydro_short` =
##     .compute_hydro_features(hydro, obs$date, max(api_tau_bounds_short), ...),
##     i.e. the API at tau = max(api_tau_bounds_short) (30 d by default);
##   * .resolve_ref_norm_instant() compares it against the target's
##     `event_api` = .compute_api(..., tau = match_window_days) (5 d default);
##   * the gate width `match_hydro_tol` = 0.5 * IQR of the API at
##     tau = api_tau_bounds_short[1] (1 d by default) over the reference dates.
##   .compute_api() is an un-normalised store (S_t = k S_{t-1} + P_t), whose
##   level scales roughly with tau, so even for the SAME date and the SAME
##   rainfall series, API(30) - API(5) is typically tens of mm while the gate
##   is ~1-4 mm: tier 1 only fires after long dry spells when both stores are
##   near zero.  On the synthetic R1 fixture below with default tau bounds,
##   0 of 143 same-site same-day (sample, analyte) cells match (6 of 143 with
##   the narrowed bounds c(3, 10)); with match_hydro_tol = Inf, 138 do, and
##   every one has ref_norm == C_norm exactly.
##
## T1. The tier-1 rain-index comparison uses the same API definition (same
##     tau) on both sides, with a tolerance on that same scale.  Hence when the
##     target IS the reference site and a detected reference observation exists
##     on the target date, tier 1 fires: ref_tier == "matched" and ref_norm ==
##     that observation's normalised value, which for the same sample is the
##     target's own C_norm (tol 1e-10).  Checked for default tau bounds and the
##     narrowed bounds, on a fixture with varied rain (gaps and NA days).
## T2. A target date with very different antecedent rain from every reference
##     observation within the match window (a 600 mm storm between a reference
##     grab on day D and a target sample on D + 3) does NOT match (it falls to
##     the model/static tiers), while the same-day sample on D does match.
## T3. The tier counts exposed by ara_summary() show tier-1 ("matched") cells
##     in the same-site case: > 0, and in fact every detected same-day cell.
##
## All data are synthetic (helper-perf-fixtures.R and the storm fixture below).

## ---------------------------------------------------------------------------
## Part 1 fixtures
## ---------------------------------------------------------------------------

.ca_sc <- function() .pf_scenario()

## Target: the reference site itself (R1).  One point keeps the seven
## pipeline runs inside the runtime budget; the equality tests (P1-P3) do not
## depend on the target being the reference site, P4 does.
.ca_target <- function() .ca_sc()$reference

.ca_tmp <- function() {
  tmp <- .ca_sc()$temperature
  tmp[tmp$datetime >= as.Date("2022-04-15") &
        tmp$datetime <= as.Date("2022-09-15"), ]
}

.ca_reference_args <- function() {
  list(
    min_obs_model = 20L,
    api_tau_bounds_short = c(3, 10),
    api_tau_bounds_long  = c(30, 90)
  )
}

.ca_daily_base <- function() {
  list(
    temperature = .ca_tmp(),
    start = as.Date("2022-05-01"), end = as.Date("2022-08-31")
  )
}

.ca_focal <- as.Date(c("2022-06-15", "2022-07-31", "2022-08-31"))

## Call `fn` with `args` bound to symbols, so an error's deparsed call stays
## short (do.call() would inline every data frame into the message).
.ca_call <- function(fn, args) {
  env <- list2env(args, parent = environment())
  cl <- as.call(c(as.name(fn), stats::setNames(
    lapply(names(args), as.name), names(args)
  )))
  eval(cl, env)
}

## Run mspaf_pipeline(); `...` may carry `ara`, `also_total`, etc.
.ca_pipeline <- function(daily_args = list(), ..., keep_attrs = FALSE) {
  sc <- .ca_sc()
  args <- list(
    target = .ca_target(), reference = sc$reference, hydro = sc$rain, impute = FALSE,
    reference_args = .ca_reference_args(),
    daily_args = utils::modifyList(.ca_daily_base(), daily_args),
    chronic = TRUE, focal_dates = .ca_focal, tau = 90, window = 365
  )
  args <- utils::modifyList(args, list(...), keep.null = TRUE)
  out <- suppressWarnings(suppressMessages(.ca_call("mspaf_pipeline", args)))
  if (!keep_attrs) {
    attributes(out)[c("daily", "reference_model", "imputation_model")] <- NULL
  }
  out
}

## Memoised runs shared across tests.
.ca_cache <- new.env()
.ca_get <- function(key, f) {
  if (is.null(.ca_cache[[key]])) .ca_cache[[key]] <- f()
  .ca_cache[[key]]
}

## ara = "contemporaneous" (with also_total = TRUE, so one run serves P1-P4:
## `value` must be unaffected by also_total, and P2 pins it to the manual
## replication), counting fit_reference_model() calls.
.ca_p_contemp <- function() .ca_get("p_contemp", function() {
  n_fit <- 0L
  real_fit <- hydroSense:::fit_reference_model
  testthat::local_mocked_bindings(
    fit_reference_model = function(...) {
      n_fit <<- n_fit + 1L
      real_fit(...)
    },
    .package = "hydroSense"
  )
  out <- .ca_pipeline(ara = "contemporaneous", also_total = TRUE,
                      keep_attrs = TRUE)
  list(out = out, n_fit = n_fit)
})
.ca_p_none     <- function() .ca_get("p_none", function()
  .ca_pipeline(ara = "none"))
.ca_p_none_def <- function() .ca_get("p_none_def", function()
  .ca_pipeline())
.ca_p_static   <- function() .ca_get("p_static", function()
  .ca_pipeline(list(reference = .ca_sc()$reference), ara = "static",
               keep_attrs = TRUE))
.ca_p_static_def <- function() .ca_get("p_static_def", function()
  .ca_pipeline(list(reference = .ca_sc()$reference), keep_attrs = TRUE))

## Manual replication of the contemporaneous pipeline (P2).
.ca_manual <- function() .ca_get("manual", function() {
  ## .pf_reference_model() is fit_reference_model(sc$reference,
  ## hydro = sc$rain, <.ca_reference_args()>) with imputation_model = NULL:
  ## exactly the pipeline's own call under impute = FALSE (shared, memoised).
  ref <- .pf_reference_model()
  daily <- suppressWarnings(suppressMessages(.ca_call(
    "mspaf_daily",
    c(list(df = .ca_target(), reference_model = ref, reference = ref,
           interpolation = "model"),
      .ca_daily_base())
  )))
  long <- tibble::tibble(
    sample_id = paste0("daily_", daily$site_id, "_", daily$date),
    site_id   = daily$site_id,
    datetime  = daily$date,
    analyte   = "msPAF",
    detected  = TRUE,
    value     = daily$mspaf
  )
  chronic <- time_weighted_aggregate(
    long, focal_dates = .ca_focal, tau = 90, tau_units = "d",
    window = 365, window_units = "d", summary = "arith_mean",
    return = "summary"
  )
  list(daily = daily, chronic = chronic)
})

## Key a chronic frame by (focal_date, site_id) and keep `cols`.
.ca_key <- function(x, cols = "value") {
  x <- as.data.frame(x)[, c("focal_date", "site_id", cols), drop = FALSE]
  x <- x[order(x$site_id, x$focal_date, method = "radix"), , drop = FALSE]
  rownames(x) <- NULL
  x
}

.ca_daily_key <- function(d) {
  d <- as.data.frame(d)[, c("date", "site_id", "mspaf"), drop = FALSE]
  d <- d[order(d$site_id, d$date, method = "radix"), , drop = FALSE]
  rownames(d) <- NULL
  d
}

## ---------------------------------------------------------------------------
## Part 1 tests
## ---------------------------------------------------------------------------

describe("P1: mspaf_pipeline(ara =) argument", {
  skip_on_cran()

  it("has an `ara` argument with choices static / contemporaneous / none", {
    fmls <- formals(mspaf_pipeline)
    expect_true("ara" %in% names(fmls))
    expect_setequal(eval(fmls$ara), c("static", "contemporaneous", "none"))
  })

  it("ara = 'contemporaneous' fits the reference model exactly once", {
    expect_identical(.ca_p_contemp()$n_fit, 1L)
  })

  it("ara = 'contemporaneous' passes the fitted reference_model as reference", {
    out <- .ca_p_contemp()$out
    daily <- attr(out, "daily")
    ara <- attr(daily, "ara_summary")
    expect_false(is.null(ara))
    expect_true("ref_tier" %in% names(ara))
    expect_true(all(!is.na(ara$ref_tier)))
    expect_true(all(ara$ref_tier %in%
                      c("matched", "model", "model_integrated", "static")))
  })

  it("ara = 'none' equals a run with no daily_args$reference", {
    expect_equal(.ca_key(.ca_p_none()), .ca_key(.ca_p_none_def()),
                 tolerance = 1e-10)
  })

  it("ara = 'static' equals today's static run (ara not supplied)", {
    s <- .ca_p_static(); d <- .ca_p_static_def()
    expect_equal(.ca_key(s), .ca_key(d), tolerance = 1e-10)
    expect_equal(.ca_daily_key(attr(s, "daily")),
                 .ca_daily_key(attr(d, "daily")), tolerance = 1e-10)
  })

  it("ara = 'static' without daily_args$reference is a classed error", {
    expect_error(
      suppressWarnings(suppressMessages(.ca_pipeline(ara = "static"))),
      class = "hydroSense_ara_static_needs_reference"
    )
  })

  it("ara = 'contemporaneous' with daily_args$reference is a classed error", {
    expect_error(
      suppressWarnings(suppressMessages(.ca_pipeline(
        list(reference = .ca_sc()$reference), ara = "contemporaneous"
      ))),
      class = "hydroSense_ara_conflict"
    )
  })

  it("an unknown ara value is an error", {
    expect_error(suppressWarnings(suppressMessages(
      .ca_pipeline(ara = "dynamic")
    )))
  })
})

describe("P2: ara = 'contemporaneous' equals the manual replication", {
  skip_on_cran()

  it("daily mspaf equals mspaf_daily(reference = fitted model)", {
    out <- .ca_p_contemp()$out
    man <- .ca_manual()
    expect_equal(.ca_daily_key(attr(out, "daily")),
                 .ca_daily_key(man$daily), tolerance = 1e-10)
  })

  it("chronic value equals time_weighted_aggregate() of that daily series", {
    out <- .ca_p_contemp()$out
    man <- .ca_manual()
    expect_equal(.ca_key(out), .ca_key(man$chronic), tolerance = 1e-10)
  })

  it("differs from both the static and the no-ARA runs", {
    c_val <- .ca_key(.ca_p_contemp()$out)$value
    expect_false(isTRUE(all.equal(c_val, .ca_key(.ca_p_static())$value)))
    expect_false(isTRUE(all.equal(c_val, .ca_key(.ca_p_none())$value)))
  })
})

describe("P3: also_total with ara = 'contemporaneous'", {
  skip_on_cran()

  it("value equals the manual replication and value_total the ara = 'none' value", {
    tot <- .ca_p_contemp()$out
    expect_true("value_total" %in% names(tot))
    expect_equal(.ca_key(tot), .ca_key(.ca_manual()$chronic), tolerance = 1e-10)
    none <- .ca_key(.ca_p_none())
    tk <- .ca_key(tot, c("value", "value_total"))
    expect_equal(tk$value_total, none$value, tolerance = 1e-10)
  })
})

describe("P4: at the reference site itself, contemporaneous ARA << static", {
  skip_on_cran()

  it("chronic msPAF median is well below static's at the reference site", {
    con <- .ca_key(.ca_p_contemp()$out)
    sta <- .ca_key(.ca_p_static())
    con <- con$value[con$site_id == "R1"]
    sta <- sta$value[sta$site_id == "R1"]
    expect_length(con, length(.ca_focal))
    expect_length(sta, length(.ca_focal))
    expect_lt(stats::median(con), 0.5 * stats::median(sta))
  })

  it("daily msPAF median is well below static's at the reference site", {
    con <- attr(.ca_p_contemp()$out, "daily")
    sta <- attr(.ca_p_static(), "daily")
    con <- con$mspaf[con$site_id == "R1"]
    sta <- sta$mspaf[sta$site_id == "R1"]
    expect_gt(length(con), 0L)
    expect_lt(stats::median(con), 0.5 * stats::median(sta))
  })
})

## ---------------------------------------------------------------------------
## Part 2 fixtures
## ---------------------------------------------------------------------------

## Same-site same-day: R1 chemistry assessed against a reference model fitted
## on R1 itself.  `bounds = "default"` uses fit_reference_model()'s default tau
## bounds (the configuration seen on real data); "narrow" is the shared
## .pf_reference_model() fixture.
.ca_ref_model <- function(bounds = c("default", "narrow")) {
  bounds <- match.arg(bounds)
  .ca_get(paste0("rm_", bounds), function() {
    if (bounds == "narrow") return(.pf_reference_model())
    sc <- .ca_sc()
    suppressWarnings(suppressMessages(fit_reference_model(
      sc$reference, hydro = sc$rain, min_obs_model = 20L
    )))
  })
}

## ara_summary() of add_mspaf(df, reference = model), joined to `detected`.
.ca_ara <- function(df, model, min_analytes = 3L) {
  out <- suppressWarnings(suppressMessages(
    add_mspaf(df, reference = model, require_temperature = FALSE,
              min_analytes = min_analytes)
  ))
  a <- ara_summary(out)
  det <- dplyr::distinct(
    dplyr::select(df, "sample_id", "analyte", "detected")
  )
  dplyr::left_join(a, det, by = c("sample_id", "analyte"))
}

## Storm fixture: a reference grab on day D in an ordinary spell, then 600 mm
## of rain on D + 1 and D + 2; a target sample on D + 3 (inside the +-5 d
## window, no other reference grab within 5 d of it) and one on D itself.
.ca_storm_D <- as.Date("2020-06-01") + 370L  # a reference grab date

.ca_storm <- function() .ca_get("storm", function() {
  set.seed(501L)
  d <- seq(as.Date("2020-01-01"), as.Date("2022-12-31"), by = "day")
  wet <- stats::runif(length(d)) < 0.35
  v <- ifelse(wet, stats::rgamma(length(d), shape = 0.8, scale = 9), 0)
  D <- .ca_storm_D
  v[d == D + 1] <- 300
  v[d == D + 2] <- 300
  rain <- tibble::tibble(date = d, value = v)

  ref_dates <- seq(as.Date("2020-06-01"), as.Date("2022-06-01"), by = "10 days")
  stopifnot(D %in% ref_dates)
  ref <- .pf_point_chem("RS", ref_dates, mult = 1, seed = 502L, ni_every = 1L,
                        bdl_frac = 0, rain = rain)
  ## min_obs_model above n keeps every analyte at the static tier (no GAM
  ## fits, fast); tier 1 is independent of the model tier.
  model <- suppressWarnings(suppressMessages(fit_reference_model(
    ref, hydro = rain, min_obs_model = 1000L
  )))

  on_D <- ref[as.Date(ref$datetime) == D, , drop = FALSE]
  on_D$sample_id <- "T_D"
  after <- on_D
  after$sample_id <- "T_D_plus_3"
  after$datetime <- D + 3
  list(rain = rain, ref = ref, model = model,
       target = dplyr::bind_rows(on_D, after))
})

## ---------------------------------------------------------------------------
## Part 2 tests
## ---------------------------------------------------------------------------

describe("T1: same-site same-day reference observations match at tier 1", {
  skip_on_cran()

  for (b in c("default", "narrow")) {
    it(paste0("tier 1 fires and ref_norm == the observed value (", b,
              " tau bounds)"), {
      sc <- .ca_sc()
      a <- .ca_ara(sc$reference, .ca_ref_model(b))
      det <- a[a$detected %in% TRUE, , drop = FALSE]
      expect_gt(nrow(det), 50L)
      not_matched <- det[det$ref_tier != "matched", , drop = FALSE]
      expect_identical(
        nrow(not_matched), 0L,
        label = paste0("detected same-day cells not at tier 1 (of ",
                       nrow(det), ")")
      )
      expect_equal(det$ref_norm, det$C_norm, tolerance = 1e-10)
    })
  }
})

describe("T2: very different antecedent rain does not match", {
  skip_on_cran()

  it("a post-storm target date falls past tier 1; the same-day date matches", {
    st <- .ca_storm()
    a <- .ca_ara(st$target, st$model, min_analytes = 1L)
    expect_true(all(c("T_D", "T_D_plus_3") %in% a$sample_id))

    after <- a[a$sample_id == "T_D_plus_3", , drop = FALSE]
    expect_gt(nrow(after), 0L)
    expect_false(any(after$ref_tier == "matched"))
    expect_true(all(after$ref_tier %in% c("model", "static")))

    on_D <- a[a$sample_id == "T_D", , drop = FALSE]
    expect_gt(nrow(on_D), 0L)
    expect_true(all(on_D$ref_tier == "matched"))
    expect_equal(on_D$ref_norm, on_D$C_norm, tolerance = 1e-10)
  })
})

describe("T3: ara_summary() reports tier-1 matches in the same-site case", {
  skip_on_cran()

  it("counts 'matched' cells > 0 (default tau bounds)", {
    a <- .ca_ara(.ca_sc()$reference, .ca_ref_model("default"))
    n_matched <- sum(a$ref_tier == "matched", na.rm = TRUE)
    expect_gt(n_matched, 0L)
    expect_gte(n_matched, sum(a$detected %in% TRUE))
  })
})
