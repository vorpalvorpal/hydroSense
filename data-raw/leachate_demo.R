# Generates `leachate_demo`: a small synthetic long-format water-quality
# dataset used in examples and the LMF vignette. NOT real monitoring data.
#
# Three sites:
#   "downstream" — leachate-impacted (~15% leachate mixing fraction)
#   "reference"  — clean background / upstream
#   "leachate"   — the leachate end-member
#
# Carries the major-ion panel (for add_lmf / to_meq), the toxicants Cu/Zn/NH3-N
# (for ssd_paf / add_amspaf), and the co-analytes pH/temperature/DOC/hardness
# needed for the ammonia and bioavailability normalisations.
#
# Concentration units: major ions in mg/L; Cu/Zn/NH3-N in ug/L (the SSD scale);
# pH unitless; temperature in degC; DOC in mg/L; hardness in mg/L CaCO3.

set.seed(2024)

# analyte | units | valence | atomic_mass | reference mean | leachate mean
spec <- tibble::tribble(
  ~analyte,        ~units,        ~valence, ~mass,   ~ref,    ~leach,
  # major ions (mg/L) -- LMF panel
  "Na",            "mg/L",        1,        22.99,   12,      900,
  "K",             "mg/L",        1,        39.10,   3,       180,
  "Ca",            "mg/L",        2,        40.08,   22,      380,
  "Mg",            "mg/L",        2,        24.31,   6,       140,
  "Cl",            "mg/L",        1,        35.45,   18,      1800,
  "SO4²⁻","mg/L", 2,        96.06,   12,      280,
  "F",             "mg/L",        1,        19.00,   0.2,     4.5,
  "NO3-N",         "mg/L",        1,        14.01,   0.8,     6,
  "NO2-N",         "mg/L",        1,        14.01,   0.03,    0.9,
  "CO3-CaCO3",     "mg/L",        2,        100.09,  4,       90,
  "HCO3-CaCO3",    "mg/L",        2,        100.09,  55,      750,
  # toxicants (ug/L) -- have SSDs; NH3-N is also the N-panel driver
  "NH3-N",         "ug/L",        1,        14.01,   50,      180000,
  "Cu",            "ug/L",        NA,       NA,      1.0,     45,
  "Zn",            "ug/L",        NA,       NA,      4,       220,
  # co-analytes / drivers (not LMF ions, no meq conversion)
  "pH",            "pH",          NA,       NA,      7.6,     7.9,
  "temperature",   "degC",        NA,       NA,      16,      20,
  "DOC",           "mg/L",        NA,       NA,      2,       120,
  "hardness",      "mg/L CaCO3",  NA,       NA,      80,      1500
)

# Site sample schedule: bi-monthly through 2024 for downstream/reference;
# leachate sampled the same dates (dates are immaterial for the end-member).
dates <- as.Date(c("2024-01-15", "2024-03-15", "2024-05-15",
                   "2024-07-15", "2024-09-15", "2024-11-15"))

# Mixing fraction of leachate into the downstream site.
f_mix <- 0.15

make_site <- function(site, prefix, frac = NULL, ref_only = FALSE) {
  rows <- lapply(seq_along(dates), function(i) {
    # Target mean per analyte for this site.
    mu <- if (!is.null(frac)) {
      spec$ref + frac * (spec$leach - spec$ref)
    } else if (ref_only) {
      spec$ref
    } else {
      spec$leach
    }
    # pH / temperature are not mixed linearly -- set them directly with jitter.
    is_ph   <- spec$analyte == "pH"
    is_temp <- spec$analyte == "temperature"
    val <- mu * (1 + stats::rnorm(nrow(spec), 0, 0.03))
    val[is_ph]   <- mu[is_ph]   + stats::rnorm(sum(is_ph), 0, 0.05)
    val[is_temp] <- mu[is_temp] + stats::rnorm(sum(is_temp), 0, 1.0)
    tibble::tibble(
      sample_id           = sprintf("%s-%02d", prefix, i),
      site_id             = site,
      datetime            = dates[i],
      analyte             = spec$analyte,
      value               = round(pmax(val, 1e-6), 4),
      detected            = TRUE,
      units.analyte       = spec$units,
      valence.analyte     = spec$valence,
      atomic_mass.analyte = spec$mass
    )
  })
  dplyr::bind_rows(rows)
}

leachate_demo <- dplyr::bind_rows(
  make_site("downstream", "DS", frac = f_mix),
  make_site("reference",  "REF", ref_only = TRUE),
  make_site("leachate",   "LEACH")
)

# Sanity: drop the synthetic Na+ value occasionally below detection to show the
# `detected` flag in use (one downstream F measurement).
leachate_demo$detected[leachate_demo$site_id == "downstream" &
                       leachate_demo$analyte == "F" &
                       leachate_demo$sample_id == "DS-03"] <- FALSE

usethis::use_data(leachate_demo, overwrite = TRUE)
