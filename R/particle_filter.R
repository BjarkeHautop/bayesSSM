#' Particle Filter Implementation with Log-Sum-Exp Trick
#'
#' This function implements a particle filter for estimating the hidden states
#' in a state space model using sequential Monte Carlo methods. Three filtering
#' variants are supported:
#' \enumerate{
#'   \item \strong{SIS:} Sequential Importance Sampling (without resampling).
#'   \item \strong{SISR:} Sequential Importance Sampling with resampling at
#'   every time step.
#'   \item \strong{SISAR:} SIS with adaptive resampling based on the Effective
#'   Sample Size (ESS). Resampling is triggered when the ESS falls below a
#'   given threshold (default \code{n / 2}).
#' }
#'
#' @param y A numeric vector of observations.
#' @param n A positive integer specifying the number of particles.
#' @param init_fn A function that initializes the particle states. It should
#' accept an integer (the number of particles) as its first argument and return
#' a vector or matrix of initial particle states.
#' @param transition_fn A function describing the state transition model. It
#' should take the current particles and the current time step as arguments and
#' return the propagated particles.
#' @param log_likelihood_fn A function that computes the log likelihoods for the
#' particles. It should accept an observation, the current particles, and the
#' current time step as arguments and return a numeric vector of log likelihood
#' values.
#' @param algorithm A character string specifying the particle filtering
#' algorithm. Options are:
#'   \describe{
#'     \item{"SIS"}{Sequential Importance Sampling (without resampling).}
#'     \item{"SISR"}{Sequential Importance Sampling with resampling at every
#'     time step.}
#'     \item{"SISAR"}{SIS with adaptive resampling based on ESS.}
#'   }
#'   The default is \code{"SISAR"}.
#' @param resample_fn A character string specifying the resampling method
#' (used with \code{"SISR"} or \code{"SISAR"}). Options include:
#'   \describe{
#'     \item{"stratified"}{(default) Stratified resampling.}
#'     \item{"systematic"}{Systematic resampling.}
#'     \item{"multinomial"}{Multinomial resampling.}
#'   }
#' @param threshold A numeric value specifying the ESS threshold for triggering
#' resampling in the \code{"SISAR"} algorithm. If not provided, it defaults to
#' \code{n / 2}.
#' @param return_particles A logical value indicating whether to return the full
#' particle history. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to \code{init_fn},
#' \code{transition_fn}, and \code{log_likelihood_fn}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{state_est}{A numeric vector of estimated states over time,
#'     computed as the weighted average of particles.}
#'     \item{ess}{A numeric vector of the Effective Sample Size (ESS) at each
#'     time step.}
#'     \item{loglike}{The accumulated log-likelihood of the observations given
#'     the model.}
#'     \item{algorithm}{A character string indicating the filtering algorithm
#'     used.}
#'     \item{particles_history}{(Optional) A matrix of particle states over
#'     time (one row per time step), returned if \code{return_particles}
#'     is \code{TRUE}.}
#'     \item{weights_history}{(Optional) A matrix of particle weights over time
#'     (one row per time step), returned if \code{return_particles} is
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
#'   specified threshold (defaulting to \code{n / 2}).
#' }
#' The Effective Sample Size (ESS) in context of particle filters is defined as
#' \deqn{ESS = \left(\sum_{i=1}^n w_i^2\right)^{-1},}
#' where \eqn{w_i} are the normalized weights of the particles.
#'
#' The default resampling method is stratified resampling, as Douc et al., 2005
#' showed that it always gives a lower variance compared to
#' multinomial resampling.
#'
#' @references Douc, R., Cappé, O., & Moulines, E. (2005). Comparison of
#' Resampling Schemes for Particle Filtering.
#' Accessible at: https://arxiv.org/abs/cs/0507025
#'
#' @export
#'
#' @examples
#' init_fn <- function(n) rnorm(n, 0, 1)
#' transition_fn <- function(particles, t) particles + rnorm(length(particles))
#' log_likelihood_fn <- function(y, particles, t) {
#'   dnorm(y, mean = particles, sd = 1, log = TRUE)
#' }
#'
#' # Generate synthetic observations.
#' y <- cumsum(rnorm(50))
#' n <- 100
#'
#' # Run the particle filter using default settings.
#' result <- particle_filter(y, n, init_fn, transition_fn, log_likelihood_fn)
#' plot(result$state_est, type = "l", col = "blue", main = "State Estimates")
particle_filter <- function(
    y, n, init_fn, transition_fn, log_likelihood_fn,
    algorithm = c("SISAR", "SISR", "SIS"),
    resample_fn = c("stratified", "systematic", "multinomial"),
    threshold = NULL, return_particles = TRUE, ...) {
  # Validate input: ensure n is a positive integer
  if (!is.numeric(n) || n <= 0) {
    stop("n must be a positive integer")
  }

  # Match provided algorithm and resampling method to valid options
  algorithm <- match.arg(algorithm)
  resample_fn <- match.arg(resample_fn)

  resample_func <- switch(resample_fn,
    multinomial = .resample_multinomial,
    stratified = .resample_stratified,
    systematic = .resample_systematic
  )

  num_steps <- length(y)
  state_est <- numeric(num_steps)
  ess_vec <- numeric(num_steps)
  loglike <- 0 # log-likelihood accumulator
  if (return_particles) {
    particles_history <- matrix(NA, nrow = num_steps, ncol = n)
    weights_history <- matrix(NA, nrow = num_steps, ncol = n)
  }

  # Helper function: log-sum-exp trick for numerical stability
  logsumexp <- function(lw) {
    max_lw <- max(lw)
    max_lw + log(sum(exp(lw - max_lw)))
  }

  # Initialization at time step 1
  particles <- init_fn(n, ...)
  log_weights <- log_likelihood_fn(y[1], particles, t = 1, ...)

  # Compute incremental log likelihood l_1 and normalize weights
  log_l_t <- logsumexp(log_weights) - log(n)
  loglike <- log_l_t
  log_normalizer <- logsumexp(log_weights)
  log_weights <- log_weights - log_normalizer
  weights <- exp(log_weights)

  state_est[1] <- sum(particles * weights)
  ess_vec[1] <- 1 / sum(weights^2)
  if (return_particles) {
    particles_history[1, ] <- particles
    weights_history[1, ] <- weights
  }

  # Set adaptive resampling threshold if needed
  if (algorithm == "SISAR" && is.null(threshold)) {
    threshold <- n / 2
  }

  # Main loop over remaining time steps
  for (t in 2:num_steps) {
    # Propagate particles through the state transition model
    particles <- transition_fn(particles, t, ...)

    # Compute log likelihoods for the current observation
    log_likelihoods <- log_likelihood_fn(y[t], particles, t, ...)

    # Update log weights: previous log(weights) plus new log likelihoods
    log_weights <- log(weights) + log_likelihoods

    # Compute the incremental log likelihood L_t
    log_l_t <- logsumexp(log_weights) - log(n)
    loglike <- loglike + log_l_t

    # Normalize weights using the log-sum-exp trick
    log_normalizer <- logsumexp(log_weights)
    log_weights <- log_weights - log_normalizer
    weights <- exp(log_weights)

    # Compute Effective Sample Size (ESS)
    ess_current <- 1 / sum(weights^2)
    ess_vec[t] <- ess_current

    # Resampling based on chosen algorithm
    if (algorithm == "SISR") {
      particles <- resample_func(particles, weights)
      weights <- rep(1 / n, n)
      ess_vec[t] <- n # Reset ESS after resampling
    } else if (algorithm == "SISAR" && ess_current < threshold) {
      particles <- resample_func(particles, weights)
      weights <- rep(1 / n, n)
      ess_vec[t] <- n # Reset ESS after resampling
    }

    # Update the state estimate as the weighted average of particles
    state_est[t] <- sum(particles * weights)
    if (return_particles) {
      particles_history[t, ] <- particles
      weights_history[t, ] <- weights
    }
  }

  # Return results as a list
  result <- list(
    state_est = state_est,
    ess = ess_vec,
    loglike = loglike,
    algorithm = algorithm
  )
  if (return_particles) {
    result$particles_history <- particles_history
    result$weights_history <- weights_history
  }
  result
}
