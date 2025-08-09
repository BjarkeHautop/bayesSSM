##### Tests for .pilot_run ######

test_that(".pilot_run works with bootstrap_filter", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    rnorm(num_particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, t, phi, sigma_x, ...) {
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }

  log_likelihood_fn <- function(y, particles, t, sigma_y, ...) {
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

  result <- .pilot_run(
    pf_wrapper = bootstrap_filter,
    y = my_data$y,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    phi = 0.8,
    sigma_x = 1,
    sigma_y = 0.5,
    resample_algorithm = "SISAR",
    resample_fn = "stratified"
  )

  expect_lt(result$target_n, 1000)
})

test_that(".pilot_run works with auxiliary_filter", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    rnorm(num_particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, t, phi, sigma_x, ...) {
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }

  log_likelihood_fn <- function(y, particles, t, sigma_y, ...) {
    dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
  }

  aux_log_likelihood_fn <- function(y, particles, t, sigma_y, ...) {
    # Can use same structure or something simpler
    dnorm(y, mean = particles, sd = sigma_y * 1.5, log = TRUE)
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

  result <- .pilot_run(
    pf_wrapper = auxiliary_filter,
    y = my_data$y,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    aux_log_likelihood_fn = aux_log_likelihood_fn,
    phi = 0.8,
    sigma_x = 1,
    sigma_y = 0.5,
    resample_algorithm = "SISR",
    resample_fn = "systematic"
  )

  expect_lt(result$target_n, 1000)
})

test_that(".pilot_run works with resample_move_filter", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    rnorm(num_particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, t, phi, sigma_x, ...) {
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }

  log_likelihood_fn <- function(y, particles, t, sigma_y, ...) {
    dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
  }

  move_fn <- function(particles, t, ...) {
    particles + rnorm(length(particles), mean = 0, sd = 0.1)
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

  result <- .pilot_run(
    pf_wrapper = resample_move_filter,
    y = my_data$y,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    move_fn = move_fn,
    phi = 0.8,
    sigma_x = 1,
    sigma_y = 0.5,
    resample_fn = "multinomial"
  )

  expect_lt(result$target_n, 1000)
})


##### Tests for .run_pilot_chain ######

test_that(".run_pilot_chain works", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    rnorm(num_particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, t, phi, ...) {
    phi * particles + rnorm(length(particles), mean = 0, sd = 1)
  }

  log_likelihood_fn <- function(y, particles, t, ...) {
    dnorm(y, mean = particles, sd = 1, log = TRUE)
  }

  simulate_ssm <- function(num_steps, phi) {
    x <- numeric(num_steps)
    y <- numeric(num_steps)
    x[1] <- rnorm(1, mean = 0, sd = 1)
    y[1] <- rnorm(1, mean = x[1], sd = 1)
    for (t in 2:num_steps) {
      x[t] <- phi * x[t - 1] + rnorm(1, mean = 0, sd = 1)
      y[t] <- x[t] + rnorm(1, mean = 0, sd = 1)
    }
    list(x = x, y = y)
  }

  my_data <- simulate_ssm(50, phi = 0.8)

  log_prior_phi <- function(phi) {
    dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }

  log_priors <- list(
    phi = log_prior_phi
  )

  result <- .run_pilot_chain(
    pf_wrapper = bootstrap_filter,
    y = my_data$y,
    pilot_m = 100,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    log_priors = log_priors,
    proposal_sd = 0.1,
    resample_algorithm = "SISAR",
    resample_fn = "stratified"
  )


  expect_lt(result$target_n, 500)
})

test_that(".run_pilot_chain works with auxiliary_filter", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    rnorm(num_particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, t, phi, ...) {
    phi * particles + rnorm(length(particles), mean = 0, sd = 1)
  }

  log_likelihood_fn <- function(y, particles, t, ...) {
    dnorm(y, mean = particles, sd = 1, log = TRUE)
  }

  aux_log_likelihood_fn <- function(y, particles, t, ...) {
    dnorm(y, mean = particles, sd = 1.5, log = TRUE)
  }

  simulate_ssm <- function(num_steps, phi) {
    x <- numeric(num_steps)
    y <- numeric(num_steps)
    x[1] <- rnorm(1, mean = 0, sd = 1)
    y[1] <- rnorm(1, mean = x[1], sd = 1)
    for (t in 2:num_steps) {
      x[t] <- phi * x[t - 1] + rnorm(1, mean = 0, sd = 1)
      y[t] <- x[t] + rnorm(1, mean = 0, sd = 1)
    }
    list(x = x, y = y)
  }

  my_data <- simulate_ssm(50, phi = 0.8)

  log_prior_phi <- function(phi) {
    dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }

  log_priors <- list(
    phi = log_prior_phi
  )

  result <- .run_pilot_chain(
    pf_wrapper = auxiliary_filter,
    y = my_data$y,
    pilot_m = 100,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    aux_log_likelihood_fn = aux_log_likelihood_fn,
    log_priors = log_priors,
    proposal_sd = 0.1,
    resample_algorithm = "SISR",
    resample_fn = "systematic"
  )

  expect_lt(result$target_n, 500)
})

test_that(".run_pilot_chain works with resample_move_filter", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    rnorm(num_particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, t, phi, ...) {
    phi * particles + rnorm(length(particles), mean = 0, sd = 1)
  }

  log_likelihood_fn <- function(y, particles, t, ...) {
    dnorm(y, mean = particles, sd = 1, log = TRUE)
  }

  move_fn <- function(particles, t, ...) {
    particles + rnorm(length(particles), mean = 0, sd = 0.1)
  }

  simulate_ssm <- function(num_steps, phi) {
    x <- numeric(num_steps)
    y <- numeric(num_steps)
    x[1] <- rnorm(1, mean = 0, sd = 1)
    y[1] <- rnorm(1, mean = x[1], sd = 1)
    for (t in 2:num_steps) {
      x[t] <- phi * x[t - 1] + rnorm(1, mean = 0, sd = 1)
      y[t] <- x[t] + rnorm(1, mean = 0, sd = 1)
    }
    list(x = x, y = y)
  }

  my_data <- simulate_ssm(50, phi = 0.8)

  log_prior_phi <- function(phi) {
    dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }

  log_priors <- list(
    phi = log_prior_phi
  )

  result <- .run_pilot_chain(
    pf_wrapper = resample_move_filter,
    y = my_data$y,
    pilot_m = 100,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    move_fn = move_fn,
    log_priors = log_priors,
    proposal_sd = 0.1,
    resample_fn = "multinomial"
  )

  expect_lt(result$target_n, 500)
})

test_that(".run_pilot_chain works and handles inputs correctly", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    rnorm(num_particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, t, phi, ...) {
    phi * particles + rnorm(length(particles), mean = 0, sd = 1)
  }

  log_likelihood_fn <- function(y, particles, t, ...) {
    dnorm(y, mean = particles, sd = 1, log = TRUE)
  }

  simulate_ssm <- function(num_steps, phi) {
    x <- numeric(num_steps)
    y <- numeric(num_steps)
    x[1] <- rnorm(1, mean = 0, sd = 1)
    y[1] <- rnorm(1, mean = x[1], sd = 1)
    for (t in 2:num_steps) {
      x[t] <- phi * x[t - 1] + rnorm(1, mean = 0, sd = 1)
      y[t] <- x[t] + rnorm(1, mean = 0, sd = 1)
    }
    list(x = x, y = y)
  }

  my_data <- simulate_ssm(50, phi = 0.8)

  log_prior_phi <- function(phi) {
    dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }

  log_priors <- list(phi = log_prior_phi)

  # Error: init_param outside prior support
  log_prior_phi <- function(phi) dunif(phi, min = 0, max = 1, log = TRUE)
  log_priors <- list(phi = log_prior_phi)

  expect_error(
    .run_pilot_chain(
      pf_wrapper = bootstrap_filter,
      y = my_data$y,
      pilot_m = 100,
      pilot_n = 100,
      pilot_reps = 10,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      proposal_sd = 0.1,
      pilot_init_params = c(phi = 1.5),
      resample_algorithm = "SISR",
      resample_fn = "systematic"
    ),
    "Initial parameter values are invalid:"
  )

  # Error: param_transform missing parameter key (phi)
  expect_error(
    .run_pilot_chain(
      pf_wrapper = bootstrap_filter,
      y = my_data$y,
      pilot_m = 100,
      pilot_n = 100,
      pilot_reps = 10,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      proposal_sd = 0.1,
      pilot_init_params = c(phi = 0.5),
      resample_algorithm = "SISR",
      resample_fn = "systematic",
      param_transform = list(sigma_x = "log")
    ),
    "param_transform must include every parameter"
  )

  # Error: param_transform not a list
  expect_error(
    .run_pilot_chain(
      pf_wrapper = bootstrap_filter,
      y = my_data$y,
      pilot_m = 100,
      pilot_n = 100,
      pilot_reps = 10,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      proposal_sd = 0.1,
      pilot_init_params = c(phi = 0.5),
      resample_algorithm = "SISR",
      resample_fn = "systematic",
      param_transform = "log"
    ),
    "param_transform must be a list"
  )

  # Warning: unsupported param_transform value
  expect_warning(
    .run_pilot_chain(
      pf_wrapper = bootstrap_filter,
      y = my_data$y,
      pilot_m = 100,
      pilot_n = 100,
      pilot_reps = 10,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      proposal_sd = 0.1,
      pilot_init_params = c(phi = 0.5),
      resample_algorithm = "SISR",
      resample_fn = "systematic",
      param_transform = list(phi = "arctan")
    ),
    "Only 'log', 'logit', 'identity' supported."
  )

  # Verbose output test
  expect_message(
    .run_pilot_chain(
      pf_wrapper = bootstrap_filter,
      y = my_data$y,
      pilot_m = 100,
      pilot_n = 100,
      pilot_reps = 10,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      proposal_sd = 0.1,
      pilot_init_params = c(phi = 0.5),
      resample_algorithm = "SISR",
      resample_fn = "systematic",
      verbose = TRUE
    ),
    "Pilot chain posterior variance:"
  )

  # Verbose output with parameter transformation
  expect_message(
    .run_pilot_chain(
      pf_wrapper = bootstrap_filter,
      y = my_data$y,
      pilot_m = 100,
      pilot_n = 100,
      pilot_reps = 10,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      proposal_sd = 0.1,
      pilot_init_params = c(phi = 0.5),
      resample_algorithm = "SISR",
      resample_fn = "systematic",
      param_transform = list(phi = "log"),
      verbose = TRUE
    ),
    "Pilot chain posterior variance \\(transformed space\\):"
  )

  # Add additional tests for multi-dimensional param_transform here if needed...
})


# More complicated example

test_that("More complicated example", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    rnorm(num_particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, phi, sigma_x, ...) {
    # X_t = phi*X_{t-1} + sin(X_{t-1}) + sigma_x*V_t,  V_t ~ N(0,1)
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }

  log_likelihood_fn <- function(y, particles, sigma_y, ...) {
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
  phi <- 0.8
  sigma_x <- 1
  sigma_y <- 0.5
  my_data <- simulate_ssm(50, phi = 0.8, sigma_x = 1, sigma_y = 0.5)

  log_prior_phi <- function(phi) {
    dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }

  log_prior_sigma_x <- function(sigma) {
    dexp(sigma, rate = 1, log = TRUE)
  }

  log_prior_sigma_y <- function(sigma) {
    dexp(sigma, rate = 1, log = TRUE)
  }

  log_priors <- list(
    phi = log_prior_phi,
    sigma_x = log_prior_sigma_x,
    sigma_y = log_prior_sigma_y
  )

  result <- .run_pilot_chain(
    pf_wrapper = bootstrap_filter,
    y = my_data$y,
    pilot_m = 1000,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    log_priors = log_priors,
    proposal_sd = c(0.1, 0.1, 0.1),
    pilot_init_params = c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
    resample_algorithm = "SISAR",
    resample_fn = "stratified"
  )

  means <- unname(result$pilot_theta_mean)
  expect_equal(means[1], phi, tolerance = 0.5)
  expect_equal(means[2], sigma_x, tolerance = 0.5)
  expect_equal(means[3], sigma_y, tolerance = 0.5)
})

test_that("More complicated example with transformation", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    rnorm(num_particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, phi, sigma_x, ...) {
    # X_t = phi*X_{t-1} + sin(X_{t-1}) + sigma_x*V_t,  V_t ~ N(0,1)
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }

  log_likelihood_fn <- function(y, particles, t, sigma_y, ...) {
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
  phi <- 0.8
  sigma_x <- 1
  sigma_y <- 0.5
  my_data <- simulate_ssm(50, phi = 0.8, sigma_x = 1, sigma_y = 0.5)

  log_prior_phi <- function(phi) {
    dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }

  log_prior_sigma_x <- function(sigma) {
    dexp(sigma, rate = 1, log = TRUE)
  }

  log_prior_sigma_y <- function(sigma) {
    dexp(sigma, rate = 1, log = TRUE)
  }

  log_priors <- list(
    phi = log_prior_phi,
    sigma_x = log_prior_sigma_x,
    sigma_y = log_prior_sigma_y
  )

  result <- .run_pilot_chain(
    pf_wrapper = bootstrap_filter,
    y = my_data$y,
    pilot_m = 1000,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    log_priors = log_priors,
    proposal_sd = c(0.1, 0.1, 0.1),
    pilot_init_params = c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
    resample_algorithm = "SISAR",
    resample_fn = "stratified",
    param_transform = list(
      phi = "identity",
      sigma_x = "log",
      sigma_y = "log"
    )
  )
  means <- unname(result$pilot_theta_mean)
  expect_equal(means[1], phi, tolerance = 0.5)
  expect_equal(means[2], sigma_x, tolerance = 0.5)
  expect_equal(means[3], sigma_y, tolerance = 0.5)
})

# Multi-dimensional example
test_that("Multi dimensional works", {
  set.seed(1405)

  init_fn <- function(num_particles, ...) {
    matrix(rnorm(num_particles * 2), ncol = 2)
  }
  transition_fn <- function(particles, phi, ...) {
    particles + rnorm(ncol(particles), mean = phi)
  }
  log_likelihood_fn <- function(y, particles, ...) {
    dnorm(y[1], mean = particles[, 1], sd = 1, log = TRUE) +
      dnorm(y[2], mean = particles[, 2], sd = 1, log = TRUE)
  }


  init_state <- matrix(rnorm(2), ncol = 2)
  num_steps <- 50
  x <- matrix(0, nrow = num_steps, ncol = 2)
  y <- matrix(0, nrow = num_steps, ncol = 2)
  phi <- 1
  x[1, ] <- init_state + rnorm(2, mean = phi)
  y[1, ] <- rnorm(2, mean = x[1, ], sd = 1)

  for (t in 2:num_steps) {
    x[t, ] <- x[t - 1, ] + rnorm(2, mean = phi)
    y[t, ] <- rnorm(2, mean = x[t, ], sd = 1)
  }
  x <- rbind(init_state, x)

  log_prior_phi <- function(phi) {
    dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }

  log_priors <- list(
    phi = log_prior_phi
  )

  result <- .run_pilot_chain(
    pf_wrapper = bootstrap_filter,
    y = y,
    pilot_m = 1000,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    log_priors = log_priors,
    proposal_sd = c(0.1),
    pilot_init_params = c(phi = 0.8),
    resample_algorithm = "SISAR",
    resample_fn = "stratified"
  )
  phi_est <- unname(result$pilot_theta_mean)

  expect_equal(phi_est, phi, tolerance = 0.1)
})
