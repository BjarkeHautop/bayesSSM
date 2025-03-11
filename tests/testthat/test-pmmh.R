#--------------------------
# Tests for default_tune_control
#--------------------------

test_that("default_tune_control returns a list with correct defaults", {
  result <- default_tune_control()

  # Check result type and names
  expect_type(result, "list")
  expect_named(result, c(
    "pilot_proposal_sd", "pilot_n", "pilot_m", "target_var",
    "pilot_burn_in", "pilot_reps", "pilot_algorithm", "pilot_resample_fn"
  ))

  # Check default values
  expect_equal(result$pilot_proposal_sd, 1)
  expect_equal(result$pilot_n, 100)
  expect_equal(result$pilot_m, 2000)
  expect_equal(result$target_var, 1)
  expect_equal(result$pilot_burn_in, 1000)
  expect_equal(result$pilot_reps, 100)
  expect_equal(result$pilot_algorithm, "SISAR")
  expect_equal(result$pilot_resample_fn, "stratified")
})

test_that("default_tune_control handles valid inputs", {
  result <- default_tune_control(
    pilot_proposal_sd = 0.5, pilot_n = 500, pilot_m = 5000,
    pilot_target_var = 2, pilot_burn_in = 2000, pilot_reps = 5,
    pilot_algorithm = "SISR", pilot_resample_fn = "systematic"
  )

  # Check valid inputs
  expect_equal(result$pilot_proposal_sd, 0.5)
  expect_equal(result$pilot_n, 500)
  expect_equal(result$pilot_m, 5000)
  expect_equal(result$target_var, 2)
  expect_equal(result$pilot_burn_in, 2000)
  expect_equal(result$pilot_reps, 5)
  expect_equal(result$pilot_algorithm, "SISR")
  expect_equal(result$pilot_resample_fn, "systematic")
})

test_that("default_tune_control errors on invalid inputs", {
  expect_error(
    default_tune_control(pilot_proposal_sd = -0.1),
    "pilot_proposal_sd must be a positive numeric value."
  )
  expect_error(
    default_tune_control(pilot_n = 0),
    "pilot_n must be a positive numeric value."
  )
  expect_error(
    default_tune_control(pilot_m = -10),
    "pilot_m must be a positive numeric value."
  )
  expect_error(
    default_tune_control(pilot_target_var = "a"),
    "pilot_target_var must be a positive numeric value."
  )
  expect_error(
    default_tune_control(pilot_burn_in = -1),
    "pilot_burn_in must be a positive numeric value."
  )
  expect_error(
    default_tune_control(pilot_algorithm = "InvalidAlg"),
    "'arg' should be one of"
  )
  expect_error(
    default_tune_control(pilot_resample_fn = "InvalidFn"),
    "'arg' should be one of"
  )
})

#--------------------------
# Tests for pmmh
#--------------------------

# -----------------------------
# Input Validation Tests for pmmh
# -----------------------------

test_that("pmmh checks input types", {
  init_fn_ssm <- function(particles) rnorm(particles, mean = 0, sd = 1)
  transition_fn_ssm <- function(particles, phi, sigma_x) {
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }
  log_likelihood_fn_ssm <- function(y, particles, sigma_y) {
    dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
  }
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

  valid_init_params <- list(phi = 0.8, sigma_x = 1, sigma_y = 0.5)
  valid_log_priors <- list(
    phi = function(phi) 0,
    sigma_x = function(sigma_x) 0,
    sigma_y = function(sigma_y) 0
  )


  # y must be numeric
  expect_error(
    pmmh(
      y = "not numeric", m = 10, init_fn_ssm = init_fn_ssm,
      transition_fn_ssm = transition_fn_ssm,
      log_likelihood_fn_ssm = log_likelihood_fn_ssm,
      log_priors = log_priors, init_params = valid_init_params, burn_in = 2
    ),
    "y must be a numeric vector"
  )

  # m must be a positive integer
  expect_error(
    pmmh(
      y = rnorm(10), m = -5, init_fn_ssm = init_fn_ssm,
      transition_fn_ssm = transition_fn_ssm,
      log_likelihood_fn_ssm = log_likelihood_fn_ssm,
      log_priors = log_priors, init_params = valid_init_params, burn_in = 2
    ),
    "m must be a positive integer"
  )

  # burn_in must be positive
  expect_error(
    pmmh(
      y = rnorm(10), m = 10, burn_in = -1, init_fn_ssm = init_fn_ssm,
      transition_fn_ssm = transition_fn_ssm,
      log_likelihood_fn_ssm = log_likelihood_fn_ssm,
      log_priors = log_priors, init_params = valid_init_params
    ),
    "burn_in must be a positive integer"
  )

  # burn_in must be smaller than m
  expect_error(
    pmmh(
      y = rnorm(10), m = 10, burn_in = 10, init_fn_ssm = init_fn_ssm,
      transition_fn_ssm = transition_fn_ssm,
      log_likelihood_fn_ssm = log_likelihood_fn_ssm,
      log_priors = log_priors, init_params = valid_init_params
    ),
    "burn_in must be smaller than"
  )

  # num_chains must be a positive integer
  expect_error(
    pmmh(
      y = rnorm(10), m = 10, burn_in = 2, num_chains = 0,
      init_fn_ssm = init_fn_ssm,
      transition_fn_ssm = transition_fn_ssm,
      log_likelihood_fn_ssm = log_likelihood_fn_ssm,
      log_priors = log_priors, init_params = valid_init_params
    ),
    "num_chains must be a positive integer"
  )

  # log-likelihood must take y as an argument
  expect_error(
    pmmh(
      y = rnorm(10), m = 10, burn_in = 2, num_chains = 1,
      init_fn_ssm = init_fn_ssm,
      transition_fn_ssm = transition_fn_ssm,
      log_likelihood_fn_ssm = function(particles, sigma_y) particles,
      log_priors = log_priors, init_params = valid_init_params
    ),
    "log_likelihood_fn_ssm does not contain 'y'"
  )
})

# -----------------------------
# Function Argument Tests for pmmh
# -----------------------------

test_that("pmmh checks function arguments", {
  valid_init_params <- list(phi = 0.8, sigma_x = 1, sigma_y = 0.5)
  valid_log_priors <- list(
    phi = function(phi) 0,
    sigma_x = function(sigma_x) 0,
    sigma_y = function(sigma_y) 0
  )

  mock_init_fn <- function(particles, phi, sigma_x) particles
  mock_transition_fn <- function(particles, phi, sigma_x) particles
  mock_log_likelihood_fn <- function(y, particles, sigma_y) particles


  # Check if functions accept 'particles'
  expect_error(
    pmmh(
      y = numeric(50),
      m = 10,
      init_fn_ssm = function(phi, sigma_x) 0,
      transition_fn_ssm = mock_transition_fn,
      log_likelihood_fn_ssm = mock_log_likelihood_fn,
      log_priors = valid_log_priors,
      init_params = valid_init_params,
      burn_in = 1
    ),
    "init_fn_ssm does not contain 'particles' as an argument"
  )

  expect_error(
    pmmh(
      y = numeric(50), m = 10, init_fn_ssm = mock_init_fn,
      transition_fn_ssm = function(phi, sigma_x) 0,
      log_likelihood_fn_ssm = mock_log_likelihood_fn,
      log_priors = valid_log_priors,
      init_params = valid_init_params, burn_in = 1
    ),
    "transition_fn_ssm does not contain 'particles' as an argument"
  )

  expect_error(
    pmmh(
      y = numeric(50), m = 10, init_fn_ssm = mock_init_fn,
      transition_fn_ssm = mock_transition_fn,
      log_likelihood_fn_ssm = function(y, sigma_y) 0,
      log_priors = valid_log_priors,
      init_params = valid_init_params, burn_in = 1
    ),
    "log_likelihood_fn_ssm does not contain 'particles' as an argument"
  )
})

# -----------------------------
# Parameter Matching Tests for pmmh
# -----------------------------

test_that("pmmh checks that parameters match init_params and log_priors", {
  valid_init_params <- list(phi = 0.8, sigma_x = 1, sigma_y = 0.5)
  valid_log_priors <- list(
    phi = function(phi) 0,
    sigma_x = function(sigma_x) 0,
    sigma_y = function(sigma_y) 0
  )

  mock_init_fn <- function(particles, phi, sigma_x) particles
  mock_transition_fn <- function(particles, phi, sigma_x) particles
  mock_log_likelihood_fn <- function(y, particles, sigma_y) particles

  invalid_init_params <- list(phi = 0.8, sigma_x = 1)
  expect_error(
    pmmh(
      y = numeric(50), m = 10, init_fn_ssm = mock_init_fn,
      transition_fn_ssm = mock_transition_fn,
      log_likelihood_fn_ssm = mock_log_likelihood_fn,
      log_priors = valid_log_priors,
      init_params = invalid_init_params, burn_in = 1
    ),
    "Parameters in functions do not match the names in init_params"
  )

  invalid_log_priors <- list(phi = function(phi) 0)
  expect_error(
    pmmh(
      y = numeric(50), m = 10, init_fn_ssm = mock_init_fn,
      transition_fn_ssm = mock_transition_fn,
      log_likelihood_fn_ssm = mock_log_likelihood_fn,
      log_priors = invalid_log_priors,
      init_params = valid_init_params, burn_in = 1
    ),
    "Parameters in functions do not match the names in log_priors"
  )
})

test_that("pmmh works with valid arguments", {
  set.seed(1405)
  init_fn_ssm <- function(particles) {
    stats::rnorm(particles, mean = 0, sd = 1)
  }
  transition_fn_ssm <- function(particles, phi, sigma_x) {
    phi * particles + sin(particles) +
      stats::rnorm(length(particles), mean = 0, sd = sigma_x)
  }
  log_likelihood_fn_ssm <- function(y, particles, sigma_y) {
    stats::dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
  }
  log_prior_phi <- function(phi) {
    stats::dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }
  log_prior_sigma_x <- function(sigma) {
    stats::dexp(sigma, rate = 1, log = TRUE)
  }
  log_prior_sigma_y <- function(sigma) {
    stats::dexp(sigma, rate = 1, log = TRUE)
  }
  log_priors <- list(
    phi = log_prior_phi,
    sigma_x = log_prior_sigma_x,
    sigma_y = log_prior_sigma_y
  )

  # Generate data
  t_val <- 20
  x <- numeric(t_val)
  y <- numeric(t_val)
  x[1] <- rnorm(1, mean = 0, sd = 1)
  y[1] <- rnorm(1, mean = x[1], sd = 0.5)
  for (t in 2:t_val) {
    x[t] <- 0.8 * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = 1)
    y[t] <- x[t] + rnorm(1, mean = 0, sd = 0.5)
  }

  expect_error(
    {
      suppressWarnings({
        pmmh_result <- pmmh(
          y = y,
          m = 1000,
          init_fn_ssm = init_fn_ssm,
          transition_fn_ssm = transition_fn_ssm,
          log_likelihood_fn_ssm = log_likelihood_fn_ssm,
          log_priors = log_priors,
          init_params = c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
          burn_in = 100,
          num_chains = 2,
          param_transform = list(
            phi = "identity",
            sigma_x = "log",
            sigma_y = "log"
          ),
          seed = 1405
        )
      })
    },
    regexp = NA
  ) # Expects that no errors are thrown
})
