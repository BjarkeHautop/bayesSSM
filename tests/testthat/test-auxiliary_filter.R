test_that("APF outperforms BPF under informative observations", {
  set.seed(1405)

  time <- 50
  mu <- 1
  sigma <- 0.1
  num_particles <- 20

  x <- numeric(time + 1)
  y <- numeric(time)
  x[1] <- rnorm(1, 0, 1)
  for (t in 1:time) {
    x[t + 1] <- x[t] + rnorm(1, mu)
    y[t] <- rnorm(1, x[t + 1], sigma)
  }

  init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
  transition_fn <- function(particles, mu) {
    particles + rnorm(length(particles), mean = mu)
  }
  log_likelihood_fn <- function(y, particles, sigma) {
    dnorm(y, mean = particles, sd = sigma, log = TRUE)
  }
  aux_log_likelihood_fn <- function(y, particles, mu, sigma) {
    forecast <- particles + mu
    dnorm(y, mean = forecast, sd = sigma, log = TRUE)
  }

  bpf <- bootstrap_filter(
    y = y,
    num_particles = num_particles,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    mu = mu,
    sigma = sigma
  )

  apf <- auxiliary_filter(
    y = y,
    num_particles = num_particles,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    aux_log_likelihood_fn = aux_log_likelihood_fn,
    mu = mu,
    sigma = sigma
  )

  mse_bpf <- mean((bpf$state_est - x)^2)
  mse_apf <- mean((apf$state_est - x)^2)

  expect_lt(mse_apf, mse_bpf)
})
