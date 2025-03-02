library(testthat)

# For testing, we define simple functions for initialization, transition, and
# likelihood.
# These functions produce deterministic outputs:
# - init_fn returns a vector of zeros.
# - transition_fn increments each particle by 1.
# - likelihood_fn returns uniform weights (ones) so that after normalization,
# weights are equal.

init_fn <- function(n, ...) rep(0, n)
transition_fn <- function(particles, t, ...) particles + 1
likelihood_fn <- function(y, particles, t, ...) rep(1, length(particles))

# A simple observation vector for testing (5 time steps)
y <- rep(0, 5)

test_that("particle_filter returns correct structure", {
  result <- particle_filter(y, n = 10, init_fn, transition_fn, likelihood_fn,
                            algorithm = "SIS")
  expect_true(is.list(result))
  expect_true("state_est" %in% names(result))
  expect_true("ess" %in% names(result))
  expect_true("algorithm" %in% names(result))
  expect_true("particles_history" %in% names(result))
})

test_that("particle_filter returns no particles_history when requested", {
  result <- particle_filter(y, n = 10, init_fn, transition_fn, likelihood_fn,
                            algorithm = "SIS", return_particles = FALSE)
  expect_false("particles_history" %in% names(result))
})

test_that("particle_filter validates n input", {
  expect_error(
    particle_filter(y, n = -10, init_fn, transition_fn, likelihood_fn,
                    algorithm = "SIS"),
    "n must be a positive integer"
  )
})

test_that("particle_filter state estimates update correctly", {
  # With the simple functions:
  # - At time t = 1: particles are 0, so state_est[1] = 0.
  # - At each subsequent time, particles are incremented by 1.
  # Thus, the state estimate should follow: 0, 1, 2, 3, 4.
  result <- particle_filter(y, n = 10, init_fn, transition_fn, likelihood_fn,
                            algorithm = "SIS")
  expect_equal(result$state_est, 0:4)
})

test_that("particle_filter ESS remains consistent", {
  # Since likelihood_fn returns uniform weights,
  # after normalization each weight equals 1/n, so ESS = n at every time step.
  result <- particle_filter(y, n = 10, init_fn, transition_fn, likelihood_fn,
                            algorithm = "SIS")
  expect_equal(result$ess, rep(10, 5))
})

test_that("particle_filter using SISR resets weights every step", {
  result <- particle_filter(y, n = 10, init_fn, transition_fn, likelihood_fn,
                            algorithm = "SISR")
  # For SISR, resampling is performed at every time step, so the ESS should be
  # reset to n.
  expect_equal(result$ess, rep(10, 5))
  expect_equal(result$algorithm, "SISR")
})

test_that("particle_filter using SISAR with forced resampling resets weights", {
  # Use a high threshold to force resampling (since ESS is lower than threshold).
  result <- particle_filter(y, n = 10, init_fn, transition_fn, likelihood_fn,
                            algorithm = "SISAR", threshold = 15)
  expect_equal(result$ess, rep(10, 5))
  expect_equal(result$algorithm, "SISAR")
})
