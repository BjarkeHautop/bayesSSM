test_that("check_params_match stops if log_likelihood_fn lacks 'y'", {
  init_fn_ssm <- function(num_particles) {
    num_particles * 2
  }

  log_likelihood_fn <- function(particles) {
    sum(particles)
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

test_that("logit transformation works correctly", {
  theta <- 0.5 # value in (0, 1)

  expect_equal(
    .transform_params(
      theta,
      c("logit")
    ),
    0
  )

  expect_equal(
    .transform_params(
      theta,
      c("logit")
    ),
    log(theta / (1 - theta))
  )

  expect_equal(
    .back_transform_params(
      log(theta / (1 - theta)),
      c("logit")
    ),
    theta
  )

  expect_equal(
    .compute_log_jacobian(
      theta,
      c("logit")
    ),
    -log(theta * (1 - theta))
  )
})
