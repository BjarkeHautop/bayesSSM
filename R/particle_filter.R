#' Particle Filter
#'
#' This function implements a bootstrap particle filter for estimating the
#' hidden states in a state space model using sequential Monte Carlo methods.
#' Three filtering variants are supported:
#' \enumerate{
#'   \item \strong{SIS:} Sequential Importance Sampling (without resampling).
#'   \item \strong{SISR:} Sequential Importance Sampling with resampling at
#'   every time step.
#'   \item \strong{SISAR:} SIS with adaptive resampling based on the Effective
#'   Sample Size (ESS). Resampling is triggered when the ESS falls below a
#'   given threshold (default \code{particles / 2}).
#' }
#' It is recommended to use either SISR or SISAR to avoid weight degeneracy.
#'
#' @param y A numeric vector or matrix of observations. Each row represents an
#' observation at a time step. If observations not equally spaced, use the
#' \code{obs_times} argument to specify the time points at which
#' observations are available.
#' @param num_particles A positive integer specifying the number of particles.
#' @param init_fn A function that initializes the particle states. It should
#' take `num_particles` as an argument for initializing the particles and return
#' a vector or matrix of initial particle states. It can include any
#' model-specific parameters as named arguments.
#' @param transition_fn A function describing the state transition model. It
#' should take `particles` as an argument and return the propagated particles.
#' The function can optionally depend on time by including a time step argument
#' `t`. It can include any model-specific parameters as named arguments.
#' @param log_likelihood_fn A function that computes the log likelihoods for the
#' particles. It should take a `y` argument for the observations,
#' the current particles, and return a numeric vector of log likelihood
#' values. The function can optionally depend on time by including a time step
#' argument `t`. It can include any model-specific parameters as named
#' arguments.
#' @param obs_times A numeric vector indicating the time points at which
#' observations in \code{y} are available. Must be of the same length as the
#' number of rows in \code{y}. If not specified, it is assumed that observations
#' are available at consecutive time steps, i.e., \code{obs_times = 1:nrow(y)}.
#' @param algorithm A character string specifying the particle filtering
#' algorithm to use. Must be one of \code{"SISAR"}, \code{"SISR"}, or
#' \code{"SIS"}. Defaults to \code{"SISAR"}.
#' @param resample_fn A character string specifying the resampling method.
#' Must be one of \code{"stratified"}, \code{"systematic"}, or
#' \code{"multinomial"}. Defaults to \code{"stratified"}.
#' @param threshold A numeric value specifying the ESS threshold for triggering
#' resampling in the \code{"SISAR"} algorithm. If not provided, it defaults to
#' \code{num_particles / 2}.
#' @param return_particles A logical value indicating whether to return the full
#' particle history. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to \code{init_fn},
#' \code{transition_fn}, and \code{log_likelihood_fn}. I.e., parameter values if
#' the functions requires them.
#'
#' @return A list containing:
#'   \describe{
#'     \item{state_est}{A numeric vector of estimated states over
#'     time, computed as the weighted average of particles.}
#'     \item{ess}{A numeric vector of the Effective Sample Size (ESS) at each
#'     time step.}
#'     \item{loglike}{The accumulated log-likelihood of the observations given
#'     the model.}
#'     \item{loglike_history}{A numeric vector of the log-likelihood at each
#'     time step.}
#'     \item{algorithm}{A character string indicating the filtering algorithm
#'     used.}
#'     \item{particles_history}{(Optional) A list of particle state matrices
#'     over time (one per time step), returned if \code{return_particles} is
#'     \code{TRUE}.}
#'     \item{weights_history}{(Optional) A list of particle weight vectors over
#'     time (one per time step), returned if \code{return_particles} is
#'     \code{TRUE}.}
#'   }
#'
#' @details
#' The particle filter is a sequential Monte Carlo method that approximates the
#' posterior distribution of the state in a state space model. The three
#' supported algorithms differ in their approach to resampling:
#' \enumerate{
#'   \item \strong{SIS:} Particles are propagated and weighted without any
#'    resampling, which may lead to weight degeneracy over time.
#'   \item \strong{SISR:} Resampling is performed at every time step to combat
#'   weight degeneracy.
#'   \item \strong{SISAR:} Resampling is performed adaptively; particles are
#'   resampled only when the Effective Sample Size (ESS) falls below a
#'   specified threshold (defaulting to \code{particles / 2}).
#' }
#' The Effective Sample Size (ESS) in context of particle filters is defined as
#' \deqn{ESS = \left(\sum_{i=1}^{\text{n}} w_i^2\right)^{-1},}
#' where \eqn{n} is the number of particles and \eqn{w_i} are the
#' normalized weights of the particles.
#'
#' The default resampling method is stratified resampling, as Douc et al., 2005
#' showed that it always gives a lower variance compared to
#' multinomial resampling.
#'
#' @references Douc, R., Cappé, O., & Moulines, E. (2005). Comparison of
#' Resampling Schemes for Particle Filtering.
#' Accessible at: https://arxiv.org/abs/cs/0507025
#'
#' @importFrom stats dnorm rnorm
#' @importFrom lifecycle deprecate_warn
#' @export
#'
#' @examples
#' init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
#' transition_fn <- function(particles) particles + rnorm(length(particles))
#' log_likelihood_fn <- function(y, particles) {
#'   dnorm(y, mean = particles, sd = 1, log = TRUE)
#' }
#'
#' # Generate data
#' y <- cumsum(rnorm(50))
#' num_particles <- 100
#'
#' # Run the particle filter using default settings.
#' result <- particle_filter(
#'   y = y,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn
#' )
#' plot(result$state_est, type = "l", col = "blue", main = "State Estimates",
#'   ylim = range(c(result$state_est, y)))
#' points(y, col = "red", pch = 20)
#'
#' # With parameters
#' init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
#' transition_fn <- function(particles, mu) {
#'   particles + rnorm(length(particles), mean = mu)
#' }
#' log_likelihood_fn <- function(y, particles, sigma) {
#'   dnorm(y, mean = particles, sd = sigma, log = TRUE)
#' }
#'
#' # Generate data
#' y <- cumsum(rnorm(50))
#' num_particles <- 100
#'
#' # Run the particle filter using default settings.
#' result <- particle_filter(
#'   y = y,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn,
#'   mu = 1,
#'   sigma = 1
#' )
#' plot(result$state_est, type = "l", col = "blue", main = "State Estimates",
#'   ylim = range(c(result$state_est, y)))
#' points(y, col = "red", pch = 20)
#'
#' # With observations gaps
#' init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
#' transition_fn <- function(particles, mu) {
#'   particles + rnorm(length(particles), mean = mu)
#' }
#' log_likelihood_fn <- function(y, particles, sigma) {
#'   dnorm(y, mean = particles, sd = sigma, log = TRUE)
#' }
#'
#' # Generate data using DGP
#' simulate_ssm <- function(num_steps, mu, sigma) {
#'   x <- numeric(num_steps)
#'   y <- numeric(num_steps)
#'   x[1] <- rnorm(1, mean = 0, sd = sigma)
#'   y[1] <- rnorm(1, mean = x[1], sd = sigma)
#'   for (t in 2:num_steps) {
#'     x[t] <- mu * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma)
#'     y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma)
#'   }
#'   y
#' }
#'
#' data <- simulate_ssm(10, mu = 1, sigma = 1)
#' # Suppose we have data for t=1,2,3,5,6,7,8,9,10 (i.e., missing at t=4)
#'
#' obs_times <- c(1, 2, 3, 5, 6, 7, 8, 9, 10)
#' data_obs <- data[obs_times]
#'
#' num_particles <- 100
#' # Run the particle filter
#' # Specify observation times in the particle filter using obs_times
#' result <- particle_filter(
#'   y = data_obs,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn,
#'   obs_times = obs_times,
#'   mu = 1,
#'   sigma = 1,
#' )
#' plot(result$state_est, type = "l", col = "blue", main = "State Estimates",
#'   ylim = range(c(result$state_est, data)))
#' points(data, col = "red", pch = 20)
particle_filter <- function(
    y, num_particles, init_fn, transition_fn, log_likelihood_fn,
    obs_times = NULL,
    algorithm = c("SISAR", "SISR", "SIS"),
    resample_fn = c("stratified", "systematic", "multinomial"),
    threshold = NULL, return_particles = TRUE, ...) {
  # Validate num_particles: must be a positive integer
  if (!is.numeric(num_particles) || num_particles <= 0) {
    stop("num_particles must be a positive integer")
  }

  # Add ... as arg to functions if not present
  has_dots <- function(fun) {
    "..." %in% names(formals(fun))
  }

  if (!has_dots(init_fn)) {
    formals(init_fn) <- c(formals(init_fn), alist(... = ))
  }
  if (!has_dots(transition_fn)) {
    formals(transition_fn) <- c(formals(transition_fn), alist(... = ))
  }
  if (!has_dots(log_likelihood_fn)) {
    formals(log_likelihood_fn) <- c(
      formals(log_likelihood_fn),
      alist(... = )
    )
  }

  # Add t ar arg to transition_fn and log_likelihood_fn if not present
  if (!"t" %in% names(formals(transition_fn))) {
    formals(transition_fn) <- c(
      formals(transition_fn),
      alist(t = NULL)
    )
  }
  if (!"t" %in% names(formals(log_likelihood_fn))) {
    formals(log_likelihood_fn) <- c(
      formals(log_likelihood_fn),
      alist(t = NULL)
    )
  }

  # ---------------------------
  # Input validation
  # ---------------------------

  # Check y matrix or vector
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }

  # Ensure y is a matrix.
  if (is.vector(y)) {
    y <- matrix(y, ncol = 1)
  }

  if (is.null(obs_times)) {
    obs_times <- seq_len(nrow(y))
  }

  # Check numeric
  if (!is.numeric(obs_times)) {
    stop("obs_times must be numeric")
  }

  if (length(obs_times) != nrow(y)) {
    stop("obs_times must match the number of observations (rows in y)")
  }

  # Check obs_times integers
  if (!all(obs_times == round(obs_times))) {
    stop("obs_times must be integers")
  }

  # Check obs_times increasing
  if (any(diff(obs_times) < 1)) {
    stop("obs_times must be strictly increasing")
  }

  # Each row of y is an observation, so set num_steps to the number of rows.
  num_steps <- nrow(y)

  # Match provided algorithm and resampling method to valid options
  algorithm <- match.arg(algorithm)
  resample_fn <- match.arg(resample_fn)

  resample_func <- switch(
    resample_fn,
    multinomial = .resample_multinomial,
    stratified = .resample_stratified,
    systematic = .resample_systematic
  )

  # Initialization: obtain initial particles and ensure they are in matrix form.

  # Deprecated argument 'particles' in init_fn.
  init_formals <- names(formals(init_fn))
  if ("particles" %in% init_formals && !"num_particles" %in% init_formals) {
    lifecycle::deprecate_warn(
      when = "1.0.0",
      what = "init_fn(particles)",
      with = "init_fn(num_particles)"
    )
    particles <- init_fn(particles = num_particles, ...)
  } else {
    particles <- init_fn(num_particles = num_particles, ...)
  }

  if (is.null(dim(particles))) {
    # If particles come as a vector, treat it as 1-dimensional.
    if (length(particles) != num_particles) {
      stop("init_fn must return a vector of length num_particles")
    }
    particles <- matrix(particles, nrow = num_particles, byrow = TRUE)
  }

  # Check that init_fn returns exactly num_particles rows
  if (nrow(particles) != num_particles) {
    stop("init_fn must return a matrix with num_particles rows")
  }

  # Save state dimension d
  d <- ncol(particles)

  # Decide whether the state is 1-dim or multi-dim.
  one_dim <- (d == 1)

  # Initialize state estimates:
  if (one_dim) {
    state_est <- numeric(num_steps)
  } else {
    state_est <- matrix(NA, nrow = num_steps, ncol = d)
  }
  ess_vec <- numeric(num_steps)
  loglike <- 0 # log-likelihood accumulator
  loglike_history <- numeric(num_steps)  # log-likelihood at each time step

  # To store history, use lists (each element is an num_particles x d matrix)
  if (return_particles) {
    particles_history <- vector("list", num_steps)
    weights_history <- vector("list", num_steps)
  }

  # Helper function: log-sum-exp trick for numerical stability
  logsumexp <- function(lw) {
    max_lw <- max(lw)
    max_lw + log(sum(exp(lw - max_lw)))
  }

  # ---------------------------
  # Time step 1 initialization
  # ---------------------------

  # Simulate from initial (t=1) to first observation time, if needed
  initial_time <- 1L
  gap1 <- obs_times[1] - initial_time
  if (gap1 > 0) {
    t_curr <- initial_time
    for (step in seq_len(gap1)) {
      t_curr <- t_curr + 1L
      particles <- transition_fn(particles = particles, t = t_curr, ...)
    }
  }
  # Weight at first observation
  log_weights <- log_likelihood_fn(
    y = y[1, ],
    particles = particles,
    t = obs_times[1],
    ...
  )

  # Check that log_likelihood_fn returns a log-likelihood of correct dimensions.
  if (!is.numeric(log_weights) || length(log_weights) != num_particles) {
    stop("log_likelihood_fn must return dimensions matching num_particles")
  }

  log_l_1 <- logsumexp(log_weights) - log(num_particles)
  loglike <- log_l_1
  loglike_history[1] <- log_l_1

  log_normalizer <- logsumexp(log_weights)
  log_weights <- log_weights - log_normalizer
  weights <- exp(log_weights)
  # Compute the weighted average for the first time step
  if (one_dim) {
    state_est[1] <- sum(particles * weights)
  } else {
    state_est[1, ] <- colSums(particles * weights)
  }
  ess_vec[1] <- 1 / sum(weights^2)
  if (return_particles) {
    particles_history[[1]] <- particles
    weights_history[[1]] <- weights
  }

  # Set adaptive resampling threshold if needed
  if (algorithm == "SISAR" && is.null(threshold)) {
    threshold <- num_particles / 2
  }

  # ----------------------------
  # Main loop over time steps
  # ----------------------------
  prev_time <- obs_times[1]

  for (i in 2:num_steps) {
    # Evolve particles from prev_time to obs_times[i]
    gap <- obs_times[i] - prev_time
    t_curr <- prev_time
    # Evolve the particles for the gap (if gap > 1, do multiple transitions)
    for (step in 1:gap) {
      t_curr <- t_curr + 1
      particles <- transition_fn(
        particles = particles,
        t = t_curr,
        ...
      )
      if (is.null(dim(particles))) {
        if (length(particles) != num_particles) {
          stop("transition_fn must return a vector of length num_particles")
        }
        particles <- matrix(
          particles,
          nrow = num_particles,
          byrow = TRUE
        )
      } else {
        if (nrow(particles) != num_particles) {
          stop(
            paste0(
              "transition_fn must return the same number of rows ",
              "as num_particles"
            )
          )
        }
      }
    }
    prev_time <- obs_times[i]

    # Now evaluate the log-likelihood at the current observation time.
    log_likelihood <- log_likelihood_fn(
      y = y[i, ],
      particles = particles,
      t = obs_times[i],
      ...
    )

    # Check if the likelihood is extremely low.
    if (all(log_likelihood < -1e8)) {
      loglike <- -Inf
      loglike_history[i] <- -Inf
      result <- list(
        state_est = state_est,
        ess = ess_vec,
        loglike = loglike,
        loglike_history = loglike_history,
        algorithm = algorithm
      )
      if (return_particles) {
        result$particles_history <- particles_history
        result$weights_history <- weights_history
      }
      return(result)
    }

    log_weights    <- log(weights) + log_likelihood
    log_normalizer <- logsumexp(log_weights)
    loglike        <- loglike + log_normalizer
    log_weights    <- log_weights - log_normalizer
    weights        <- exp(log_weights)

    ess_current <- 1 / sum(weights^2)
    ess_vec[i] <- ess_current

    # Resampling step for SISR and SISAR
    if (algorithm == "SISR") {
      particles <- resample_func(particles, weights)
      weights <- rep(1 / num_particles, num_particles)
      ess_vec[i] <- num_particles
    } else if (algorithm == "SISAR" && ess_current < threshold) {
      particles <- resample_func(particles, weights)
      weights <- rep(1 / num_particles, num_particles)
      ess_vec[i] <- num_particles
    }

    # Compute state estimates
    if (one_dim) {
      state_est[i] <- sum(particles * weights)
    } else {
      state_est[i, ] <- colSums(particles * weights)
    }

    if (return_particles) {
      particles_history[[i]] <- particles
      weights_history[[i]] <- weights
    }
  }

  result <- list(
    state_est = state_est,
    ess = ess_vec,
    loglike = loglike,
    loglike_history = loglike_history,
    algorithm = algorithm
  )
  if (return_particles) {
    result$particles_history <- particles_history
    result$weights_history <- weights_history
  }
  result
}
