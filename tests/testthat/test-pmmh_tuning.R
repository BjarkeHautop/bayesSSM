test_that(".pilot_run works non-trivial setup", {
  set.seed(1405)
  init_fn <- function(particles, ...) {
    rnorm(particles, mean = 0, sd = 1)
  }

  transition_fn <- function(particles, t, phi, sigma_x, ...) {
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
  my_data <- simulate_ssm(50, phi = 0.8, sigma_x = 1, sigma_y = 0.5)
  result <- .pilot_run(
    y = my_data$y,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    phi = 0.8,
    sigma_x = 1,
    sigma_y = 0.5,
    algorithm = "SISAR",
    resample_fn = "systematic"
  )
  expect_lt(result$target_n, 1000)
})

# .run_pilot_chain works
test_that(".run_pilot_chain works", {
  set.seed(1405)
  init_fn <- function(particles, ...) {
    rnorm(particles, mean = 0, sd = 1)
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
    y = my_data$y,
    pilot_m = 100,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    log_priors = log_priors,
    proposal_sd = c(0.1),
    algorithm = "SISR",
    resample_fn = "systematic"
  )
  expect_lt(result$target_n, 500)

  # Expect error if init_param outside of log_priors domain
  log_prior_phi <- function(phi) {
    dunif(phi, min = 0, max = 1, log = TRUE)
  }
  log_priors <- list(
    phi = log_prior_phi
  )
  expect_error(.run_pilot_chain(
    y = my_data$y,
    pilot_m = 100,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    log_priors = log_priors,
    proposal_sd = c(0.1),
    pilot_init_params = c(phi = 1.5),
    algorithm = "SISR",
    resample_fn = "systematic"
  ), "Invalid initial parameters:")

  # Expect error if param_transform does not have phi
  expect_error(.run_pilot_chain(
    y = my_data$y,
    pilot_m = 100,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    log_priors = log_priors,
    proposal_sd = c(0.1),
    pilot_init_params = c(phi = 0.5),
    algorithm = "SISR",
    resample_fn = "systematic",
    param_transform = list(
      sigma_x = "log"
    )
  ), "param_transform must include an entry")

  # Expect error if param_transform not a list
  expect_error(.run_pilot_chain(
    y = my_data$y,
    pilot_m = 100,
    pilot_n = 100,
    pilot_reps = 10,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    log_priors = log_priors,
    proposal_sd = c(0.1),
    pilot_init_params = c(phi = 0.5),
    algorithm = "SISR",
    resample_fn = "systematic",
    param_transform = "log"
  ), "param_transform must be a list")

  # Check verbose works
  expect_message(
    .run_pilot_chain(
      y = my_data$y,
      pilot_m = 100,
      pilot_n = 100,
      pilot_reps = 10,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      proposal_sd = c(0.1),
      pilot_init_params = c(phi = 0.5),
      algorithm = "SISR",
      resample_fn = "systematic",
      verbose = TRUE
    ),
    "Pilot chain posterior variance:"
  )

  # Check works with transformation
  expect_message(
    .run_pilot_chain(
      y = my_data$y,
      pilot_m = 100,
      pilot_n = 100,
      pilot_reps = 10,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      proposal_sd = c(0.1),
      pilot_init_params = c(phi = 0.5),
      algorithm = "SISR",
      resample_fn = "systematic",
      param_transform = list(
        phi = "log"
      ),
      verbose = TRUE
    ),
    "Pilot chain posterior variance \\(on transformed space\\):"
  )
})


# More complicated example

test_that("More complicated example", {
  set.seed(1405)

  init_fn <- function(particles, ...) {
    rnorm(particles, mean = 0, sd = 1)
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
    algorithm = "SISAR",
    resample_fn = "stratified"
  )
  means <- unname(result$pilot_theta_mean)
  expect_equal(means[1], phi, tolerance = 0.5)
  expect_equal(means[2], sigma_x, tolerance = 0.5)
  expect_equal(means[3], sigma_y, tolerance = 0.5)
})

test_that("More complicated example with transformation", {
  set.seed(1405)

  init_fn <- function(particles, ...) {
    rnorm(particles, mean = 0, sd = 1)
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
    algorithm = "SISAR",
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

  init_fn <- function(particles, ...) matrix(rnorm(particles * 2), ncol = 2)
  transition_fn <- function(particles, phi, ...) {
    particles + rnorm(nrow(particles) * 2, mean = phi)
  }
  log_likelihood_fn <- function(y, particles, ...) rep(1, nrow(particles))

  # A vector of observation
  y <- rep(0, 5)

  log_prior_phi <- function(phi) {
    stats::dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }

  log_priors <- list(
    phi = log_prior_phi
  )

  result <- .run_pilot_chain(
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
    algorithm = "SISAR",
    resample_fn = "stratified"
  )
  phi_est <- unname(result$pilot_theta_mean)

  expect_equal(phi_est, 0, tolerance = 0.1)
})
