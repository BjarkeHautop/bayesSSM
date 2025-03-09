init_fn <- function(n, ...) rep(0, n)
transition_fn <- function(particles, ...) particles + 1
likelihood_fn <- function(y, particles, ...) rep(1, length(particles))

# A simple observation vector for testing (5 time steps)
y <- rep(0, 5)

test_that("particle_filter returns correct structure", {
  result <- particle_filter(y,
    n = 10, init_fn, transition_fn, likelihood_fn,
    algorithm = "SIS"
  )
  expect_true(is.list(result))
  expect_true("state_est" %in% names(result))
  expect_true("ess" %in% names(result))
  expect_true("algorithm" %in% names(result))
  expect_true("particles_history" %in% names(result))
})

test_that("particle_filter returns no particles_history when requested", {
  result <- particle_filter(y,
    n = 10, init_fn, transition_fn, likelihood_fn,
    algorithm = "SIS", return_particles = FALSE
  )
  expect_false("particles_history" %in% names(result))
})

test_that("particle_filter validates n input", {
  expect_error(
    particle_filter(y,
      n = -10, init_fn, transition_fn, likelihood_fn,
      algorithm = "SIS"
    ),
    "n must be a positive integer"
  )
})

# Works with more complicated setup
init_fn_ssm <- function(n, ...) {
  rnorm(n, mean = 0, sd = 1)
}

transition_fn_ssm <- function(particles, phi, sigma_x, ...) {
  # X_t = phi*X_{t-1} + sin(X_{t-1}) + sigma_x*V_t,  V_t ~ N(0,1)
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}

log_likelihood_fn_ssm <- function(y, particles, sigma_y, ...) {
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

test_that("particle_filter works non-trivial setup", {
  result <- particle_filter(
    y = my_data$y,
    n = 100,
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
