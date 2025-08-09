test_that("RMPF outperforms BPF under strong particle degeneracy", {
  set.seed(1405)

  time <- 50
  mu <- 1
  sigma <- 0.05 # highly informative observations
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

  # Simple MH move function: random walk with fixed proposal
  move_fn <- function(particle, y, sigma) {
    proposal <- particle + rnorm(1, 0, 0.1)
    log_p_curr <- log_likelihood_fn(y = y, particles = particle, sigma = sigma)
    log_p_prop <- log_likelihood_fn(y = y, particles = proposal, sigma = sigma)
    if (log(runif(1)) < (log_p_prop - log_p_curr)) {
      proposal
    } else {
      particle
    }
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

  rmpf <- resample_move_filter(
    y = y,
    num_particles = num_particles,
    init_fn = init_fn,
    transition_fn = transition_fn,
    log_likelihood_fn = log_likelihood_fn,
    move_fn = move_fn,
    mu = mu,
    sigma = sigma
  )

  mse_bpf <- mean((bpf$state_est - x)^2)
  mse_rmpf <- mean((rmpf$state_est - x)^2)

  expect_lt(mse_rmpf, mse_bpf)
})
