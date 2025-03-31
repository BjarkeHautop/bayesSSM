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
#' @param y A numeric vector or matrix of observations.
#' @param num_particles A positive integer specifying the number of particles.
#' @param init_fn A function that initializes the particle states. It should
#' take the current particles as its first argument and return
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
#' \code{particles / 2}.
#' @param return_particles A logical value indicating whether to return the full
#' particle history. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to \code{init_fn},
#' \code{transition_fn}, and \code{log_likelihood_fn}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{state_est}{A numeric vector of estimated states over
#'     time, computed as the weighted average of particles.}
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
#' @export
#'
#' @examples
#' init_fn <- function(particles) rnorm(particles, 0, 1)
#' transition_fn <- function(particles) particles + rnorm(length(particles))
#' log_likelihood_fn <- function(y, particles) {
#'   stats::dnorm(y, mean = particles, sd = 1, log = TRUE)
#' }
#'
#' # Generate data.
#' y <- cumsum(rnorm(50))
#' num_particles <- 100
#'
#' # Run the particle filter using default settings.
#' result <- result <- particle_filter(
#'   y = y,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn
#' )
#' plot(result$state_est, type = "l", col = "blue", main = "State Estimates")
particle_filter <- function(
    y, num_particles, init_fn, transition_fn, log_likelihood_fn,
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

  # Ensure y is a matrix.
  if (is.vector(y)) {
    y <- matrix(y, ncol = 1)
  }

  # Each row of y is an observation, so set num_steps to the number of rows.
  num_steps <- nrow(y)

  # Match provided algorithm and resampling method to valid options
  algorithm <- match.arg(algorithm)
  resample_fn <- match.arg(resample_fn)

  resample_func <- switch(resample_fn,
    multinomial = .resample_multinomial,
    stratified = .resample_stratified,
    systematic = .resample_systematic
  )

  # Initialization: obtain initial particles and ensure they are in matrix form.
  particles <- init_fn(particles = num_particles, ...)
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
  # Evaluate log-likelihood function at the first observation.
  log_weights <- log_likelihood_fn(y = y[1, ], particles = particles, ...)

  # Check that log_likelihood_fn returns a log-likelihood of correct dimensions.
  if (!is.numeric(log_weights) || length(log_weights) != num_particles) {
    stop("log_likelihood_fn must return dimensions matching num_particles")
  }

  log_l_1 <- logsumexp(log_weights) - log(num_particles)
  loglike <- log_l_1
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
  for (i in 2:num_steps) {
    # Transition: update particles
    particles_new <- transition_fn(particles = particles, ...)
    if (is.null(dim(particles_new))) {
      # Check that transition_fn preserves correct output dimensions
      if (length(particles_new) != num_particles || ncol(particles_new) != d) {
        stop(paste0(
          "transistion_fn must return the same dimensions as ",
          "the initial particles"
        ))
      }
      particles_new <- matrix(particles_new, nrow = num_particles, byrow = TRUE)
    } else {
      if (nrow(particles_new) != num_particles || ncol(particles_new) != d) {
        stop(paste0(
          "transistion_fn must return the same dimensions as ",
          "the initial particles"
        ))
      }
    }

    particles <- particles_new

    # Evaluate log-likelihood function for the i-th observation.
    log_likelihood <- log_likelihood_fn(y = y[i, ], particles = particles, ...)

    if (!is.numeric(log_likelihood) ||
        length(log_likelihood) != num_particles) {
      stop(paste0(
        "log_likelihood_fn must return a numeric vector of ",
        "length num_particles"
      ))
    }

    # If all likelihoods are 0 return loglike as -Inf
    if (all(log_likelihood < -1e8)) {
      loglike <- -Inf
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
      return(result)
    }

    log_weights <- log(weights) + log_likelihood
    log_l_i <- logsumexp(log_weights) - log(num_particles)
    loglike <- loglike + log_l_i

    log_normalizer <- logsumexp(log_weights)
    log_weights <- log_weights - log_normalizer
    weights <- exp(log_weights)
    if (is.nan(sum(weights))) {
      stop(paste0(
        "Weights are NaN because of zero likelihoods. ",
        "Recheck model and consider modifying arguments in tune_control()."
      ))
    }

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
    algorithm = algorithm
  )
  if (return_particles) {
    result$particles_history <- particles_history
    result$weights_history <- weights_history
  }
  result
}
