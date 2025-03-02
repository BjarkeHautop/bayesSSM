#' Particle Filter Implementation
#'
#' This function implements a particle filter for estimating the hidden states
#' in a state space model using sequential Monte Carlo methods. It supports
#' three variants: Sequential Importance Sampling (SIS), Sequential Importance
#' Sampling with Resampling (SISR), and SIS with Adaptive Resampling (SISAR).
#' The latter adapts resampling based on the Effective Sample Size (ESS).
#'
#' @param y A numeric vector of observations.
#' @param n An integer specifying the number of particles. Must be positive.
#' @param init_fn A function that initializes the particle states. It should
#'   accept an integer (the number of particles) as its first argument and return
#'   a vector or matrix of initial particle states.
#' @param transition_fn A function describing the state transition model.
#'   It should take the current particles and the current time step as
#'   arguments, and return the propagated particles.
#' @param likelihood_fn A function that computes the likelihood weights for the
#'   particles. It should accept an observation, the current particles, and the
#'   current time step as arguments.
#' @param algorithm A character string specifying the particle filtering
#'   algorithm to use. Options are:
#'   \describe{
#'     \item{"SIS"}{Sequential Importance Sampling (without resampling).}
#'     \item{"SISR"}{Sequential Importance Sampling with resampling at every
#'     time step.}
#'     \item{"SISAR"}{SIS with adaptive resampling based on Effective Sample
#'     Size (ESS).}
#'   }
#'   The default is \code{"SISAR"}.
#' @param resample_fn A character string specifying the resampling method
#'   (used with \code{"SISR"} or \code{"SISAR"}).
#'   Options include:
#'   \describe{
#'     \item{"stratified"}{(default) Stratified resampling.}
#'     \item{"systematic"}{Systematic resampling.}
#'     \item{"multinomial"}{Multinomial resampling.}
#'   }
#' @param threshold A numeric value specifying the ESS threshold for triggering
#'   resampling in the \code{"SISR"} or \code{"SISAR"} algorithm.
#'   If not provided, it defaults to \code{n / 2}.
#' @param return_particles A logical value indicating whether to return the full
#'   particle history.
#'   Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to \code{init_fn},
#'   \code{transition_fn}, and \code{likelihood_fn}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{state_est}{A numeric vector of estimated states over time,
#'     computed as the weighted average of particles.}
#'     \item{ess}{A numeric vector containing the Effective Sample Size (ESS) at
#'     each time step.}
#'     \item{algorithm}{A character string indicating the filtering algorithm
#'     used.}
#'      \item{weights_history}{(Optional) A matrix of particle weights over time
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
#' \deqn{ESS = \frac{1}{\sum_{i=1}^n w_i^2},}
#' where \eqn{w_i} are the normalized weights of the particles.
#'
#' The default resampling method is stratified resampling, as Douc et al., 2005
#' showed that it always gives a lower variance compared to
#' multinomial resampling.
#'
#' @export
#'
#' @examples
#' # Define example functions for initialization, transition, and likelihood.
#' init_fn <- function(n) rnorm(n, 0, 1)
#' transition_fn <- function(particles, t) particles + rnorm(length(particles))
#' likelihood_fn <- function(y, particles, t) dnorm(y, mean = particles, sd = 1)
#'
#' # Generate some synthetic observations.
#' y <- cumsum(rnorm(50))
#'
#' # Run the particle filter using default settings.
#' result <- particle_filter(y, n = 100, init_fn, transition_fn, likelihood_fn)
#' plot(result$state_est, type = "l", col = "blue", main = "State Estimates")
particle_filter <- function(
    y, n, init_fn, transition_fn, likelihood_fn,
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

  # Select the appropriate resampling function
  resample_func <- switch(resample_fn,
                          multinomial = .resample_multinomial,
                          stratified = .resample_stratified,
                          systematic = .resample_systematic)

  num_steps <- length(y)
  state_est <- numeric(num_steps)
  ess_vec <- numeric(num_steps)
  if (return_particles) {
    particles_history <- matrix(NA, nrow = num_steps, ncol = n)
    weights_history <- matrix(NA, nrow = num_steps, ncol = n)
  }

  # Initialization at time step 1
  particles <- init_fn(n, ...)
  weights <- likelihood_fn(y[1], particles, t = 1, ...)
  weights <- weights / sum(weights)
  state_est[1] <- sum(particles * weights)
  ess_vec[1] <- 1 / sum(weights^2)
  if (return_particles) {
    particles_history[1, ] <- particles
    weights_history[1, ] <- weights
  }

  # For adaptive resampling, set threshold if not provided
  if (algorithm == "SISAR" && is.null(threshold)) {
    threshold <- n / 2
  }

  # Main loop over remaining time steps
  for (time_step in 2:num_steps) {
    # Propagate particles through the state transition model
    particles <- transition_fn(particles, time_step, ...)

    # Compute likelihood weights for the current observation
    likelihoods <- likelihood_fn(y[time_step], particles, time_step, ...)
    weights <- weights * likelihoods
    weights <- weights / sum(weights)

    # Compute Effective Sample Size (ess)
    ess_current <- 1 / sum(weights^2)
    ess_vec[time_step] <- ess_current

    # Resampling based on chosen algorithm
    if (algorithm == "SISR") {
      particles <- resample_func(particles, weights)
      weights <- rep(1 / n, n)
      ess_vec[time_step] <- n  # Reset ess after resampling
    } else if (algorithm == "SISAR" && ess_current < threshold) {
      particles <- resample_func(particles, weights)
      weights <- rep(1 / n, n)
      ess_vec[time_step] <- n  # Reset ess after resampling
    }

    # Update the state estimate as the weighted average of particles
    state_est[time_step] <- sum(particles * weights)

    if (return_particles) {
      particles_history[time_step, ] <- particles
      weights_history[time_step, ] <- weights
    }
  }

  # Return results as a list
  result <- list(state_est = state_est, ess = ess_vec, algorithm = algorithm)
  if (return_particles) {
    result$particles_history <- particles_history
    result$weights_history <- weights_history
  }
  return(result)
}
