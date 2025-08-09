test_that("bootstrap_filter returns errors on wrong input", {
  init_fn <- function(num_particles) rep(0, num_particles)
  transition_fn <- function(particles) particles + 1
  log_likelihood_fn <- function(y, particles) rep(1, length(particles))

  wrong_init_fn <- function(num_particles) rep(0, num_particles + 1)
  wrong_transition_fn <- function(particles) c(particles, 1)
  wrong_log_likelihood_fn <- function(y, particles) {
    rep(1, length(particles) + 1)
  }

  wrong_init_fn_matrix <- function(num_particles) {
    matrix(rep(0, num_particles * 2), ncol = 5)
  }
  # A simple observation vector for testing (5 time steps)
  y <- rep(0, 5)
  wrong_y <- "hi"

  expect_error(
    bootstrap_filter(y,
      num_particles = 0, init_fn, transition_fn,
      log_likelihood_fn, resample_algorithm = "SISR"
    ),
    "Assertion on 'num_particles' failed"
  )

  expect_error(
    bootstrap_filter(y,
      num_particles = 10, wrong_init_fn, transition_fn,
      log_likelihood_fn, resample_algorithm = "SISR"
    ),
    "init_fn must return num_particles"
  )

  expect_error(
    bootstrap_filter(
      y,
      num_particles = 10, wrong_init_fn_matrix, transition_fn,
      log_likelihood_fn, resample_algorithm = "SISR"
    ),
    "init_fn must return num_particles rows"
  )

  expect_error(
    bootstrap_filter(
      y,
      num_particles = 10, init_fn, wrong_transition_fn,
      log_likelihood_fn, resample_algorithm = "SISR"
    ),
    "transition_fn must return num_particles"
  )

  expect_error(
    bootstrap_filter(
      y,
      num_particles = 10, init_fn, transition_fn,
      wrong_log_likelihood_fn, resample_algorithm = "SISR"
    ),
    "weight_fn must return num_particles"
  )

  expect_error(
    bootstrap_filter(
      wrong_y,
      num_particles = 10, init_fn, transition_fn,
      log_likelihood_fn
    ),
    "Assertion on 'y' failed"
  )

  wrong_len_obs_times <- 1:4

  expect_error(
    bootstrap_filter(
      y,
      num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
      obs_times = wrong_len_obs_times
    ),
    "Assertion on 'obs_times' failed"
  )

  not_numeric_obs_times <- "hi"

  expect_error(
    bootstrap_filter(
      y,
      num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
      obs_times = not_numeric_obs_times
    ),
    "Assertion on 'obs_times' failed"
  )

  non_integer_obs_times <- c(1.5, 2.5, 3.5, 4.5, 5.5)
  expect_error(
    bootstrap_filter(
      y,
      num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
      obs_times = non_integer_obs_times
    ),
    "Assertion on 'obs_times' failed"
  )

  non_inc_obs_times <- c(1, 2, 3, 5, 4)
  expect_error(
    bootstrap_filter(
      y,
      num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
      obs_times = non_inc_obs_times
    ),
    "Assertion on 'obs_times' failed"
  )
})


test_that("bootstrap_filter returns correct structure", {
  init_fn <- function(num_particles) rep(0, num_particles)
  transition_fn <- function(particles) particles + 1
  log_likelihood_fn <- function(y, particles) rep(1, length(particles))

  # A simple observation vector for testing (5 time steps)
  y <- rep(0, 5)

  result <- bootstrap_filter(y,
    num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
    resample_algorithm = "SIS"
  )
  expect_true(is.list(result))
  expect_true("state_est" %in% names(result))
  expect_true("ess" %in% names(result))
  expect_true("resample_algorithm" %in% names(result))
  expect_true("particles_history" %in% names(result))
})

test_that("bootstrap_filter returns no particles_history when requested", {
  init_fn <- function(num_particles) rep(0, num_particles)
  transition_fn <- function(particles) particles + 1
  log_likelihood_fn <- function(y, particles) rep(1, length(particles))

  # A simple observation vector for testing (5 time steps)
  y <- rep(0, 5)

  result <- bootstrap_filter(y,
    num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
    resample_algorithm = "SIS", return_particles = FALSE
  )
  expect_false("particles_history" %in% names(result))
})

test_that("bootstrap_filter works non-trivial setup", {
  set.seed(1405)
  # Works with more complicated setup
  init_fn_ssm <- function(num_particles) {
    rnorm(num_particles, mean = 0, sd = 1)
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
    init_state <- rnorm(1, mean = 0, sd = sigma_x)
    x <- numeric(num_steps)
    y <- numeric(num_steps)
    x[1] <- phi * init_state + sin(init_state) +
      rnorm(1, mean = 0, sd = sigma_x)
    y[1] <- x[1] + rnorm(1, mean = 0, sd = sigma_y)
    for (t in 2:num_steps) {
      x[t] <- phi * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma_x)
      y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma_y)
    }
    x <- c(init_state, x)
    list(x = x, y = y)
  }
  my_data <- simulate_ssm(50, phi = 0.8, sigma_x = 1, sigma_y = 0.5)

  result <- bootstrap_filter(
    y = my_data$y,
    num_particles = 100,
    init_fn = init_fn_ssm,
    transition_fn = transition_fn_ssm,
    log_likelihood_fn = log_likelihood_fn_ssm,
    phi = 0.8,
    sigma_x = 1,
    sigma_y = 0.5,
    resample_algorithm = "SISAR",
    resample_fn = "systematic"
  )
  # Expect length of result$state_est to be equal to length of x
  expect_equal(length(result$state_est), length(my_data$x))
  rmse <- sqrt(mean((result$state_est - my_data$x)^2))
  plot(
    result$state_est,
    type = "l",
    col = "blue",
    main = "Particle Filter State Estimate vs True State",
    xlab = "Time",
    ylab = "State Estimate"
  )
  lines(my_data$x, col = "red")
  expect_lt(rmse, 0.5)
})


# Check works with matrix
test_that("Multi-dim particle filter works", {
  init_fn <- function(num_particles) matrix(rnorm(num_particles * 2), ncol = 2)
  transition_fn <- function(particles) {
    particles + rnorm(nrow(particles) * 2)
  }
  log_likelihood_fn <- function(y, particles) rep(1, nrow(particles))

  # A vector of observation
  y <- rep(0, 5)

  result <- bootstrap_filter(y,
    num_particles = 10, init_fn, transition_fn, log_likelihood_fn,
    resample_algorithm = "SIS"
  )
  expect_true(is.list(result))
  expect_true("state_est" %in% names(result))
  expect_true("ess" %in% names(result))
  expect_true("resample_algorithm" %in% names(result))
  expect_true("particles_history" %in% names(result))
})
