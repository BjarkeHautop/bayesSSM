test_that("particle_filter_core returns errors on wrong input", {
  init_fn <- function(num_particles) rep(0, num_particles)
  transition_fn <- function(particles) particles + 1
  weight_fn <- function(y, particles) rep(1, length(particles))

  wrong_init_fn <- function(num_particles) rep(0, num_particles + 1)
  wrong_transition_fn <- function(particles) c(particles, 1)
  wrong_weight_fn <- function(y, particles) {
    rep(1, length(particles) + 1)
  }

  wrong_init_fn_matrix <- function(num_particles) {
    matrix(rep(0, num_particles * 2), ncol = 5)
  }
  # A simple observation vector for testing (5 time steps)
  y <- rep(0, 5)
  wrong_y <- "hi"

  expect_error(
    .particle_filter_core(y,
                          num_particles = 0, init_fn, transition_fn,
                          weight_fn, resample_algorithm = "SIS"
    ),
    "Assertion on 'num_particles' failed"
  )

  expect_error(
    .particle_filter_core(y,
                          num_particles = 10, wrong_init_fn, transition_fn,
                          weight_fn, resample_resample_algorithm = "SISAR"
    ),
    "init_fn must return num_particles"
  )

  expect_error(
    .particle_filter_core(y,
                          num_particles = 10, wrong_init_fn_matrix,
                          transition_fn, weight_fn, resample_algorithm = "SIS"
    ),
    "init_fn must return num_particles rows"
  )

  expect_error(
    .particle_filter_core(y,
                          num_particles = 10, init_fn,
                          wrong_transition_fn, weight_fn, resample_algorithm = "SIS"
    ),
    "transition_fn must return num_particles"
  )

  expect_error(
    .particle_filter_core(y,
                          num_particles = 10, init_fn,
                          transition_fn, wrong_weight_fn, resample_algorithm = "SIS"
    ),
    "weight_fn must return num_particles"
  )

  expect_error(
    .particle_filter_core(wrong_y,
                          num_particles = 10, init_fn,
                          transition_fn, weight_fn, resample_algorithm = "SIS"
    ),
    "Assertion on 'y' failed"
  )

  wrong_len_obs_times <- 1:4

  expect_error(
    .particle_filter_core(
      y,
      num_particles = 10, init_fn, transition_fn, weight_fn,
      obs_times = wrong_len_obs_times
    ),
    "Assertion on 'obs_times' failed"
  )

  not_numeric_obs_times <- "hi"

  expect_error(
    .particle_filter_core(
      y,
      num_particles = 10, init_fn, transition_fn, weight_fn,
      obs_times = not_numeric_obs_times
    ),
    "Assertion on 'obs_times' failed"
  )

  non_integer_obs_times <- c(1.5, 2.5, 3.5, 4.5, 5.5)
  expect_error(
    .particle_filter_core(
      y,
      num_particles = 10, init_fn, transition_fn, weight_fn,
      obs_times = non_integer_obs_times
    ),
    "Assertion on 'obs_times' failed"
  )

  non_inc_obs_times <- c(1, 2, 3, 5, 4)
  expect_error(
    .particle_filter_core(
      y,
      num_particles = 10, init_fn, transition_fn, weight_fn,
      obs_times = non_inc_obs_times
    ),
    "Assertion on 'obs_times' failed"
  )
})
