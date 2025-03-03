init_fn_ssm <- function(n, ...) {
  rnorm(n, mean = 0, sd = 1)
}

transition_fn_ssm <- function(particles, t, phi, sigma_x, ...) {
  # X_t = phi*X_{t-1} + sin(X_{t-1}) + sigma_x*V_t,  V_t ~ N(0,1)
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}

log_likelihood_fn_ssm <- function(y, particles, t, sigma_y, ...) {
  # Y_t ~ N(X_t, sigma_y^2)
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

test_that(".pilot_run works non-trivial setup", {
  result <- .pilot_run(
    y = my_data$y,
    pilot_n = 100,
    pilot_reps = 10,
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
