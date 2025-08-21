#--------------------------
# Tests for default_tune_control
#--------------------------

test_that("default_tune_control returns a list with correct defaults", {
  result <- default_tune_control()

  # Check result type and names
  expect_type(result, "list")
  expect_named(result, c(
    "pilot_proposal_sd", "pilot_n", "pilot_m", "pilot_target_var",
    "pilot_burn_in", "pilot_reps", "pilot_resample_algorithm",
    "pilot_resample_fn"
  ))

  # Check default values
  expect_equal(result$pilot_proposal_sd, 0.5)
  expect_equal(result$pilot_n, 100)
  expect_equal(result$pilot_m, 2000)
  expect_equal(result$pilot_target_var, 1)
  expect_equal(result$pilot_burn_in, 500)
  expect_equal(result$pilot_reps, 100)
  expect_equal(result$pilot_resample_algorithm, "SISAR")
  expect_equal(result$pilot_resample_fn, "stratified")
})

test_that("default_tune_control handles valid inputs", {
  result <- default_tune_control(
    pilot_proposal_sd = 0.5, pilot_n = 500, pilot_m = 5000,
    pilot_target_var = 2, pilot_burn_in = 2000, pilot_reps = 5,
    pilot_resample_algorithm = "SISR", pilot_resample_fn = "systematic"
  )

  # Check valid inputs
  expect_equal(result$pilot_proposal_sd, 0.5)
  expect_equal(result$pilot_n, 500)
  expect_equal(result$pilot_m, 5000)
  expect_equal(result$pilot_target_var, 2)
  expect_equal(result$pilot_burn_in, 2000)
  expect_equal(result$pilot_reps, 5)
  expect_equal(result$pilot_resample_algorithm, "SISR")
  expect_equal(result$pilot_resample_fn, "systematic")
})

test_that("default_tune_control errors on invalid inputs", {
  expect_error(
    default_tune_control(pilot_proposal_sd = -0.1),
    "Assertion on 'pilot_proposal_sd' failed"
  )
  expect_error(
    default_tune_control(pilot_n = 0),
    "Assertion on 'pilot_n' failed"
  )
  expect_error(
    default_tune_control(pilot_m = -10),
    "Assertion on 'pilot_m' failed"
  )
  expect_error(
    default_tune_control(pilot_target_var = "a"),
    "Assertion on 'pilot_target_var' failed"
  )
  expect_error(
    default_tune_control(pilot_burn_in = -1),
    "Assertion on 'pilot_burn_in' failed"
  )
  expect_error(
    default_tune_control(pilot_resample_algorithm = "InvalidAlg"),
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
  init_fn <- function(num_particles) rnorm(num_particles, mean = 0, sd = 1)
  transition_fn <- function(particles, phi, sigma_x) {
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }
  log_likelihood_fn <- function(y, particles, sigma_y) {
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

  valid_init_params <- list(c(phi = 0.8, sigma_x = 1, sigma_y = 0.5))
  valid_log_priors <- list(
    phi = function(phi) 0,
    sigma_x = function(sigma_x) 0,
    sigma_y = function(sigma_y) 0
  )

  wrong_init_params <- c(phi = 0.8, sigma_x = 1, sigma_y = 0.5)

  wrong_init_params_length <- list(
    c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
    c(phi = 0.8, sigma_x = 1, sigma_y = 0.5, phi = 0.5)
  )

  wrong_init_params_names <- list(
    c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
    c(phi = 0.8, sigma_x = 1, sigma_w = 0.5)
  )



  # y must be numeric
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = "not numeric", m = 10, init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors, pilot_init_params = valid_init_params,
      burn_in = 2, num_chains = 1
    ),
    "Assertion on 'y' failed"
  )

  # m must be a positive integer
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = rnorm(10), m = -5, init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors, pilot_init_params = valid_init_params,
      burn_in = 2, num_chains = 1
    ),
    "Assertion on 'm' failed"
  )

  # burn_in must be positive
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = rnorm(10), m = 10, burn_in = -1, init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors, pilot_init_params = valid_init_params,
      num_chains = 1
    ),
    "Assertion on 'burn_in' failed"
  )

  # burn_in must be smaller than m
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = rnorm(10), m = 10, burn_in = 10, init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors, pilot_init_params = valid_init_params,
      num_chains = 1
    ),
    "Assertion on 'burn_in' failed"
  )

  # num_chains must be a positive integer
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = rnorm(10), m = 10, burn_in = 2, num_chains = 0,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors, pilot_init_params = valid_init_params
    ),
    "Assertion on 'num_chains' failed"
  )

  # log-likelihood must take y as an argument
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = rnorm(10), m = 10, burn_in = 2, num_chains = 1,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = function(particles, sigma_y) particles,
      log_priors = log_priors, pilot_init_params = valid_init_params
    ),
    "log_likelihood_fn does not contain 'y' as an argument"
  )

  # Verify init_params
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = rnorm(10), m = 10, burn_in = 2, num_chains = 2,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors, pilot_init_params = wrong_init_params
    ),
    "Assertion on 'pilot_init_params' failed"
  )

  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = rnorm(10), m = 10, burn_in = 2, num_chains = 2,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors, pilot_init_params = wrong_init_params_length
    ),
    "Assertion on 'pilot_init_params' failed"
  )

  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = rnorm(10), m = 10, burn_in = 2, num_chains = 2,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors, pilot_init_params = wrong_init_params_names
    ),
    "Assertion on 'pilot_init_params' failed"
  )

  # Pilot_init_param not same length as num_chains
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = rnorm(10), m = 10, burn_in = 2, num_chains = 1,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors, pilot_init_params = wrong_init_params_names
    ),
    "Assertion on 'pilot_init_params' failed"
  )
})

# -----------------------------
# Function Argument Tests for pmmh
# -----------------------------

test_that("pmmh checks function arguments", {
  valid_init_params <- list(c(phi = 0.8, sigma_x = 1, sigma_y = 0.5))
  valid_log_priors <- list(
    phi = function(phi) 0,
    sigma_x = function(sigma_x) 0,
    sigma_y = function(sigma_y) 0
  )

  mock_init_fn <- function(num_particles, phi, sigma_x) num_particles
  mock_transition_fn <- function(particles, phi, sigma_x) particles
  mock_log_likelihood_fn <- function(y, particles, sigma_y) particles


  # Check if functions accept 'particles'
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = numeric(50),
      m = 10,
      init_fn = function(phi, sigma_x) 0,
      transition_fn = mock_transition_fn,
      log_likelihood_fn = mock_log_likelihood_fn,
      log_priors = valid_log_priors,
      pilot_init_params = valid_init_params,
      burn_in = 1,
      num_chains = 1
    ),
    "init_fn does not contain 'particles' or 'num_particles' as an argument"
  )

  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = numeric(50), m = 10, init_fn = mock_init_fn,
      transition_fn = function(phi, sigma_x) 0,
      log_likelihood_fn = mock_log_likelihood_fn,
      log_priors = valid_log_priors,
      pilot_init_params = valid_init_params, burn_in = 1,
      num_chains = 1
    ),
    "transition_fn does not contain 'particles' as an argument"
  )

  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = numeric(50), m = 10, init_fn = mock_init_fn,
      transition_fn = mock_transition_fn,
      log_likelihood_fn = function(y, sigma_y) 0,
      log_priors = valid_log_priors,
      pilot_init_params = valid_init_params, burn_in = 1,
      num_chains = 1
    ),
    "log_likelihood_fn does not contain 'particles' as an argument"
  )
})

# -----------------------------
# Parameter Matching Tests for pmmh
# -----------------------------

test_that("pmmh checks that parameters match init_params and log_priors", {
  valid_init_params <- list(c(phi = 0.8, sigma_x = 1, sigma_y = 0.5))
  valid_log_priors <- list(
    phi = function(phi) 0,
    sigma_x = function(sigma_x) 0,
    sigma_y = function(sigma_y) 0
  )

  mock_init_fn <- function(num_particles, phi, sigma_x) num_particles
  mock_transition_fn <- function(particles, phi, sigma_x) particles
  mock_log_likelihood_fn <- function(y, particles, sigma_y) particles

  invalid_init_params <- list(c(phi = 0.8, sigma_x = 1, sigmay = 0.5))
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = numeric(50), m = 10, init_fn = mock_init_fn,
      transition_fn = mock_transition_fn,
      log_likelihood_fn = mock_log_likelihood_fn,
      log_priors = valid_log_priors,
      pilot_init_params = invalid_init_params, burn_in = 1,
      num_chains = 1
    ),
    "Parameters in functions do not match the names in pilot_init_params"
  )

  invalid_log_priors <- list(phi = function(phi) 0)
  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = numeric(50), m = 10, init_fn = mock_init_fn,
      transition_fn = mock_transition_fn,
      log_likelihood_fn = mock_log_likelihood_fn,
      log_priors = invalid_log_priors,
      pilot_init_params = valid_init_params, burn_in = 1,
      num_chains = 1
    ),
    "Parameters in functions do not match the names in log_priors"
  )
})

test_that("pmmh works with valid arguments", {
  set.seed(1405)
  init_fn <- function(num_particles) {
    rnorm(num_particles, mean = 0, sd = 1)
  }
  transition_fn <- function(particles, phi, sigma_x) {
    phi * particles + sin(particles) +
      rnorm(length(particles), mean = 0, sd = sigma_x)
  }
  log_likelihood_fn <- function(y, particles, sigma_y) {
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

  # Generate data
  t_val <- 10
  init_state <- rnorm(1, mean = 0, sd = 1)
  x <- numeric(t_val)
  y <- numeric(t_val)
  x[1] <- 0.8 * init_state + sin(init_state) +
    rnorm(1, mean = 0, sd = 1)
  y[1] <- x[1] + rnorm(1, mean = 0, sd = 0.5)
  for (t in 2:t_val) {
    x[t] <- 0.8 * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = 1)
    y[t] <- x[t] + rnorm(1, mean = 0, sd = 0.5)
  }
  x <- c(init_state, x)

  expect_error(
    {
      suppressWarnings({
        pmmh_result <- pmmh(
          pf_wrapper = bootstrap_filter,
          y = y,
          m = 500,
          init_fn = init_fn,
          transition_fn = transition_fn,
          log_likelihood_fn = log_likelihood_fn,
          log_priors = log_priors,
          pilot_init_params = list(
            c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
            c(phi = 0.5, sigma_x = 0.5, sigma_y = 1)
          ),
          burn_in = 100,
          num_chains = 2,
          param_transform = list(
            phi = "identity",
            sigma_x = "log",
            sigma_y = "log"
          ),
          seed = 1405,
          resample_algorithm = "SISR",
          resample_fn = "systematic"
        )
      })
    },
    regexp = NA
  ) # Expects that no errors are thrown

  # Swapping order in param_transform should not affect the result
  expect_error(
    {
      suppressWarnings({
        pmmh_result_swapped <- pmmh(
          pf_wrapper = bootstrap_filter,
          y = y,
          m = 500,
          init_fn = init_fn,
          transition_fn = transition_fn,
          log_likelihood_fn = log_likelihood_fn,
          log_priors = log_priors,
          pilot_init_params = list(
            c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
            c(phi = 0.5, sigma_x = 0.5, sigma_y = 1)
          ),
          burn_in = 100,
          num_chains = 2,
          param_transform = list(
            sigma_x = "log",
            phi = "identity",
            sigma_y = "log"
          ),
          seed = 1405
        )
      })
    },
    regexp = NA
  ) # Expects that no errors are thrown

  # Check if the results are the same
  expect_equal(pmmh_result$theta_chain, pmmh_result_swapped$theta_chain)

  # Using several cores
  expect_error(
    {
      suppressWarnings({
        pmmh_result_2_cores <- pmmh(
          pf_wrapper = bootstrap_filter,
          y = y,
          m = 500,
          init_fn = init_fn,
          transition_fn = transition_fn,
          log_likelihood_fn = log_likelihood_fn,
          log_priors = log_priors,
          pilot_init_params = list(
            c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
            c(phi = 0.5, sigma_x = 0.5, sigma_y = 1)
          ),
          burn_in = 100,
          num_chains = 2,
          param_transform = list(
            phi = "identity",
            sigma_x = "log",
            sigma_y = "log"
          ),
          seed = 1405,
          num_cores = 2
        )
      })
    },
    regexp = NA
  ) # Expects that no errors are thrown

  # Verify pmmh_result and pmmh_result_2_cores are exactly the same
  expect_equal(
    pmmh_result$theta_chain[1:10, ],
    pmmh_result_2_cores$theta_chain[1:10, ]
  )

  expect_error(
    pmmh(
      pf_wrapper = bootstrap_filter,
      y = y,
      m = 500,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      pilot_init_params = list(
        c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
        c(phi = 0.5, sigma_x = 0.5, sigma_y = 1)
      ),
      burn_in = 100,
      num_chains = 2,
      param_transform = list(
        phi = "identity",
        sigma_x = "log",
        sigma_y = "log"
      ),
      seed = 1405,
      num_cores = 0
    ),
    "Assertion on 'num_cores' failed"
  )
})

test_that("pmmh works with resample_move_filter", {
  run_long_tests <- FALSE  # Change to TRUE to run this test

  if (run_long_tests) {
    test_that("pmmh works with resample_move_filter", {
      set.seed(1405)

      time <- 20
      mu <- 1
      num_particles <- 20

      x <- numeric(time + 1)
      y <- numeric(time)
      x[1] <- rnorm(1, 0, 1)
      for (t in 1:time) {
        x[t + 1] <- x[t] + rnorm(1, mu)
        y[t] <- rnorm(1, x[t + 1], 0.1)
      }

      init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
      transition_fn <- function(particles, mu) {
        particles + rnorm(length(particles), mean = mu)
      }
      log_likelihood_fn <- function(y, particles) {
        dnorm(y, mean = particles, sd = 0.1, log = TRUE)
      }

      move_fn <- function(particle, y) {
        proposal <- particle + rnorm(1, 0, 0.1)
        log_p_curr <- log_likelihood_fn(y = y, particles = particle)
        log_p_prop <- log_likelihood_fn(y = y, particles = proposal)
        if (log(runif(1)) < (log_p_prop - log_p_curr)) {
          proposal
        } else {
          particle
        }
      }

      log_prior_mu <- function(mu) {
        dnorm(mu, mean = 0, sd = 10, log = TRUE)
      }

      log_priors <- list(mu = log_prior_mu)

      expect_error(
        {
          suppressWarnings({
            pmmh_result <- pmmh(
              pf_wrapper = resample_move_filter,
              y = y,
              m = 200,
              init_fn = init_fn,
              transition_fn = transition_fn,
              log_likelihood_fn = log_likelihood_fn,
              log_priors = log_priors,
              pilot_init_params = list(
                c(mu = 1.5)
              ),
              burn_in = 50,
              num_chains = 1,
              param_transform = list(
                mu = "identity"
              ),
              seed = 1405,
              verbose = FALSE,
              move_fn = move_fn,
              tune_control = default_tune_control(
                pilot_m = 200, pilot_burn_in = 10
              )
            )
          })
        },
        regexp = NA
      )

      mu_chains <- as.data.frame(pmmh_result$theta_chain)
      mu_est <- mean(mu_chains$mu)

      expect_equal(mu_est, mu, tolerance = 0.1)
    })
  }

  expect_equal(run_long_tests, FALSE)
})






# Multi-dimensional example
test_that("Multi dimensional works", {
  set.seed(1405)

  init_fn <- function(num_particles) {
    matrix(rnorm(num_particles * 2), ncol = 2)
  }
  transition_fn <- function(particles, phi) {
    particles + rnorm(nrow(particles) * 2, mean = phi)
  }
  log_likelihood_fn <- function(y, particles) rep(1, nrow(particles))

  # A vector of observation
  y <- rep(0, 20)

  log_prior_phi <- function(phi) {
    dnorm(phi, mean = 0, sd = 1, log = TRUE)
  }

  log_priors <- list(
    phi = log_prior_phi
  )

  expect_error(
    {
      suppressWarnings({
        pmmh_result <- pmmh(
          pf_wrapper = bootstrap_filter,
          y = y,
          m = 500,
          init_fn = init_fn,
          transition_fn = transition_fn,
          log_likelihood_fn = log_likelihood_fn,
          log_priors = log_priors,
          pilot_init_params = list(c(phi = 0.8), c(phi = 0.5)),
          burn_in = 100,
          num_chains = 2,
          param_transform = list(
            phi = "identity"
          ),
          seed = 1405
        )
      })
    },
    regexp = NA
  ) # Expects that no errors are thrown
  phi_chains <- as.data.frame(pmmh_result$theta_chain)
  phi_est <- mean(phi_chains$phi)

  expect_equal(phi_est, 0, tolerance = 0.1)
})
