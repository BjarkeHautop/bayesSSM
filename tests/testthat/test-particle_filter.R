test_that("particle_filter returns errors on wrong input", {
  init_fn <- function(particles) rep(0, particles)
  transition_fn <- function(particles) particles + 1
  log_likelihood_fn <- function(y, particles) rep(1, length(particles))

  wrong_init_fn <- function(particles) rep(0, particles + 1)
  wrong_transition_fn <- function(particles) c(particles, 1)
  wrong_log_likelihood_fn <- function(y, particles) {
    rep(1, length(particles) + 1)
  }

  wrong_init_fn_matrix <- function(particles) {
    matrix(rep(0, particles * 2), ncol = 5)
  }
  # A simple observation vector for testing (5 time steps)
  y <- rep(0, 5)

  expect_error(
    particle_filter(y,
      num_particles = 0, init_fn, transition_fn,
      log_likelihood_fn, algorithm = "SIS"
    ),
    "particles must be a positive integer"
  )

  expect_error(
    particle_filter(y,
      num_particles = 10, wrong_init_fn, transition_fn,
      log_likelihood_fn, algorithm = "SIS"
    ),
    "init_fn must return a vector of length num_particles"
  )

  expect_error(
    particle_filter(
      y, num_particles = 10, wrong_init_fn_matrix, transition_fn,
      log_likelihood_fn, algorithm = "SIS"
    ),
    "init_fn must return a matrix with num_particles rows"
  )

  expect_error(
    particle_filter(
      y, num_particles = 10, init_fn, wrong_transition_fn,
      log_likelihood_fn, algorithm = "SIS"
    ),
    "transistion_fn must return the same dimensions as"
  )

  expect_error(
    particle_filter(
      y, num_particles = 10, init_fn, transition_fn,
      wrong_log_likelihood_fn, algorithm = "SIS"
    ),
    "log_likelihood_fn must return dimensions matching num_particles"
  )
})


test_that("particle_filter returns correct structure", {
  init_fn <- function(particles) rep(0, particles)
  transition_fn <- function(particles) particles + 1
  log_likelihood_fn <- function(y, particles) rep(1, length(particles))

  # A simple observation vector for testing (5 time steps)
  y <- rep(0, 5)

  result <- particle_filter(y,
    num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
    algorithm = "SIS"
  )
  expect_true(is.list(result))
  expect_true("state_est" %in% names(result))
  expect_true("ess" %in% names(result))
  expect_true("algorithm" %in% names(result))
  expect_true("particles_history" %in% names(result))
})

test_that("particle_filter returns no particles_history when requested", {
  init_fn <- function(particles) rep(0, particles)
  transition_fn <- function(particles) particles + 1
  log_likelihood_fn <- function(y, particles) rep(1, length(particles))

  # A simple observation vector for testing (5 time steps)
  y <- rep(0, 5)

  result <- particle_filter(y,
    num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
    algorithm = "SIS", return_particles = FALSE
  )
  expect_false("particles_history" %in% names(result))
})

test_that("particle_filter works non-trivial setup", {
  # Works with more complicated setup
  init_fn_ssm <- function(particles) {
    rnorm(particles, mean = 0, sd = 1)
  }

  transition_fn_ssm <- function(particles, phi, sigma_x) {
    # X_t = phi*X_{t-1} + sin(X_{t-1}) + sigma_x*V_t,  V_t ~ N(0,1)
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }

  log_likelihood_fn_ssm <- function(y, particles, sigma_y) {
    dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
  }

  simulate_ssm <- function(num_steps, phi, sigma_x, sigma_y) {
    x <- numeric(num_steps)
    y <- numeric(num_steps)
    x[1] <- rnorm(1, mean = 0, sd = sigma_x)
    y[1] <- rnorm(1, mean = x[1], sd = sigma_y)
    for (t in 2:num_steps) {
      x[t] <- phi * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma_x)
      y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma_y)
    }
    list(x = x, y = y)
  }
  my_data <- simulate_ssm(50, phi = 0.8, sigma_x = 1, sigma_y = 0.5)

  result <- particle_filter(
    y = my_data$y,
    num_particles = 100,
    init_fn = init_fn_ssm,
    transition_fn = transition_fn_ssm,
    log_likelihood_fn = log_likelihood_fn_ssm,
    phi = 0.8,
    sigma_x = 1,
    sigma_y = 0.5,
    algorithm = "SISAR",
    resample_fn = "systematic"
  )
  rmse <- sqrt(mean((result$state_est - my_data$x)^2))
  expect_lt(rmse, 1)
})


# Check works with matrix
test_that("Multi-dim particle filter works", {
  init_fn <- function(particles) matrix(rnorm(particles * 2), ncol = 2)
  transition_fn <- function(particles) {
    particles + rnorm(nrow(particles) * 2)
  }
  log_likelihood_fn <- function(y, particles) rep(1, nrow(particles))

  # A vector of observation
  y <- rep(0, 5)

  result <- particle_filter(y,
    num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
    algorithm = "SIS"
  )
  expect_true(is.list(result))
  expect_true("state_est" %in% names(result))
  expect_true("ess" %in% names(result))
  expect_true("algorithm" %in% names(result))
  expect_true("particles_history" %in% names(result))
})
