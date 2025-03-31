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

test_that("logit works correctly", {
  expect_equal(.transform_params(c(0.5, 0.5), c("logit", "logit")), c(0, 0))
  x <- 2.19722457733622
  expect_equal(
    .transform_params(c(0.1, 0.9), c("logit", "logit")),
    c(-x, x)
  )
  expect_equal(
    .transform_params(c(0.9, 0.1), c("logit", "logit")),
    c(x, -x)
  )

  expect_equal(
    .back_transform_params(c(0, 0), c("logit", "logit")),
    c(0.5, 0.5)
  )
  expect_equal(
    .back_transform_params(c(-x, x), c("logit", "logit")),
    c(0.1, 0.9)
  )
  expect_equal(
    .back_transform_params(c(x, -x), c("logit", "logit")),
    c(0.9, 0.1)
  )

  expect_equal(
    .compute_log_jacobian(c(0.5, 0.5), c("logit", "logit")), -1.386294 * 2,
    tolerance = 1e-6
  )
  expect_equal(
    .compute_log_jacobian(c(0.1, 0.9), c("logit", "logit")), -2.407946 * 2,
    tolerance = 1e-6
  )
  expect_equal(
    .compute_log_jacobian(c(0.9, 0.1), c("logit", "logit")), -2.407946 * 2,
    tolerance = 1e-6
  )
})
