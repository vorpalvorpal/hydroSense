## Tests for get_silo_air_temp() — the SILO Data Drill air-temperature helper.
## The network/API call to weatherOz::get_data_drill() is mocked so these run
## offline and deterministically.

library(testthat)
library(hydroSense)

skip_if_not_installed("weatherOz")

# A fake SILO Data Drill response: daily max/min air temp for a short range.
fake_drill <- function(longitude, latitude, start_date, end_date, values, api_key) {
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  tibble::tibble(
    longitude = longitude,
    latitude  = latitude,
    date      = dates,
    max_temp  = seq(20, 20 + length(dates) - 1),       # 20, 21, 22, ...
    min_temp  = seq(10, 10 + length(dates) - 1)        # 10, 11, 12, ...
  )
}

test_that("get_silo_air_temp returns mean air temp in estimate_water_temp shape", {
  testthat::local_mocked_bindings(
    get_data_drill = fake_drill,
    .package = "weatherOz"
  )

  out <- get_silo_air_temp(
    latitude   = -33.87, longitude = 151.21,
    start_date = "2020-01-01", end_date = "2020-01-05",
    api_key    = "test@example.com", cache = FALSE
  )

  expect_s3_class(out, "tbl_df")
  expect_named(out, c("datetime", "air_temp_mean_C"))
  expect_s3_class(out$datetime, "Date")
  expect_equal(nrow(out), 5L)
  # mean of (max=20, min=10) = 15; (21,11)=16; ...
  expect_equal(out$air_temp_mean_C, c(15, 16, 17, 18, 19))
})

test_that("get_silo_air_temp caches and reuses results without re-hitting the API", {
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  call_count <- 0L
  counting_drill <- function(...) {
    call_count <<- call_count + 1L
    fake_drill(...)
  }
  testthat::local_mocked_bindings(
    get_data_drill = counting_drill,
    .package = "weatherOz"
  )

  args <- list(latitude = -33.87, longitude = 151.21,
               start_date = "2020-01-01", end_date = "2020-01-03",
               api_key = "test@example.com", cache = TRUE)

  first  <- do.call(get_silo_air_temp, args)
  second <- do.call(get_silo_air_temp, args)   # should hit cache

  expect_equal(call_count, 1L)                 # API called only once
  expect_identical(first, second)

  # refresh = TRUE forces a re-fetch
  do.call(get_silo_air_temp, c(args, refresh = TRUE))
  expect_equal(call_count, 2L)
})

test_that("get_silo_air_temp validates coordinates and dates", {
  expect_error(
    get_silo_air_temp(latitude = 10, longitude = 151,           # lat outside grid
                      start_date = "2020-01-01", end_date = "2020-01-02"),
    "latitude"
  )
  expect_error(
    get_silo_air_temp(latitude = -33, longitude = 200,          # lon outside grid
                      start_date = "2020-01-01", end_date = "2020-01-02"),
    "longitude"
  )
  testthat::local_mocked_bindings(get_data_drill = fake_drill, .package = "weatherOz")
  expect_error(
    get_silo_air_temp(latitude = -33, longitude = 151,
                      start_date = "2020-01-05", end_date = "2020-01-01",  # end < start
                      api_key = "x@y.com", cache = FALSE),
    "before"
  )
})
