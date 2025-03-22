test_that("check_params_match stops if log_likelihood_fn lacks 'y'", {
  log_likelihood_fn <- function(particles) {
    sum(particles)
  }

  init_fn_ssm <- function(particles) {
    particles * 2
  }

  transition_fn_ssm <- function(particles) {
    particles + 1
  }

  init_params <- list(param1 = 0.5)
  log_priors <- list(param1 = function(x) dnorm(x, 0, 1, log = TRUE))

  expect_error(
    .check_params_match(
      init_fn_ssm, transition_fn_ssm, log_likelihood_fn,
      init_params, log_priors
    ),
    "log_likelihood_fn does not contain 'y' as an argument"
  )
})
