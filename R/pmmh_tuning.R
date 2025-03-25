#' Pilot Run for Particle Filter Tuning
#'
#' This internal function repeatedly evaluates the particle filter in order to
#' estimate the variance of the log-likelihoods and to compute a recommended
#' target number of particles for the Particle Marginal Metropolis Hastings
#' (PMMH) algorithm.
#'
#' @param y A numeric vector of observations.
#' @param pilot_n An integer specifying the initial number of particles to use.
#' @param pilot_reps An integer specifying the number of repetitions for the
#' pilot run.
#' @param ... Additional arguments passed to `init_fn`,
#' `transition_fn`, and `log_likelihood_fn`.
#'
#' @return A list containing:
#' \describe{
#'   \item{variance_estimate}{The estimated variance of the log-likelihoods
#'   from the pilot run.}
#'   \item{target_N}{The number of particles used in PMMH algorithm.}
#'   \item{pilot_loglikes}{A numeric vector of log-likelihood values computed
#'   during the run.}
#' }
#'
#' @details The function performs \code{pilot_reps} evaluations of the particle
#' filter using the provided parameter vector \code{theta}. It then estimates
#' the variance of the log-likelihoods and scales the initial particle number
#' by this variance. The final number of particles is taken as the ceiling of
#' the scaled value with a minimum of 100 and a maximum of 1000.
#'
#' @keywords internal
.pilot_run <- function(y, pilot_n, pilot_reps, init_fn,
                       transition_fn, log_likelihood_fn,
                       algorithm = c("SISAR", "SISR", "SIS"),
                       resample_fn = c(
                         "stratified", "systematic",
                         "multinomial"
                       ),
                       ...) {
  pilot_loglikes <- numeric(pilot_reps)
  for (i in seq_len(pilot_reps)) {
    pf_result <- particle_filter(
      y = y,
      num_particles = pilot_n,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      algorithm = "SISAR",
      resample_fn = "stratified",
      return_particles = FALSE,
      ...
    )
    pilot_loglikes[i] <- pf_result$loglike
  }
  variance_estimate <- stats::var(pilot_loglikes)
  target_n <- ceiling(pilot_n * variance_estimate)
  target_n <- max(target_n, 100) # Ensure a minimum number of particles
  target_n <- min(target_n, 1000) # Limit to 1000 particles

  list(
    variance_estimate = variance_estimate,
    target_n = target_n,
    pilot_loglikes = pilot_loglikes
  )
}

#' Run Pilot Chain for Posterior Estimation
#'
#' @param y A numeric vector of observations.
#' @param pilot_m An integer specifying the number of iterations for the pilot
#' chain.
#' @param pilot_n An integer specifying the number of particles for the particle
#' filter.
#' @param pilot_reps An integer specifying the number of repetitions for the
#' pilot run.
#' @param init_fn A function to initialize the state space model.
#' @param transition_fn A function that defines the state transition
#' dynamics in the state space model.
#' @param log_likelihood_fn A function that computes the log-likelihood of
#' the observations given the model parameters.
#' @param log_priors A list of functions representing the log-priors for each
#' model parameter.
#' @param proposal_sd A numeric vector specifying the standard deviations for
#' the random walk proposal distribution for each parameter.
#' @param algorithm A character string specifying the particle filter algorithm
#' to use. One of "SISAR", "SISR", or "SIS". Default is "SISAR".
#' @param resample_fn A character string specifying the resampling method to
#' use. One of "stratified", "systematic", or "multinomial". Default is
#' "stratified".
#' @param param_transform A character vector specifying the parameter
#' transformations when proposing parameters using a random walk.
#' Currently only supports "log" for log-transformation and
#' "identity" for no transformation. Default is `NULL`, which correspond to
#' no transformation ("identity).
#' @param init_params A numeric vector of initial parameter values. If `NULL`,
#' the default is a vector of ones. Default is `NULL`.
#' @param ... Additional arguments passed to the particle filter function.
#'
#' @return A list containing:
#' \item{pilot_theta_mean}{A numeric vector of the posterior mean of the
#' parameters.}
#' \item{pilot_theta_cov}{A matrix of the posterior covariance (or variance if
#' only one parameter).}
#' \item{target_N}{The estimated target number of particles for the PMMH
#' algorithm.}
#' \item{pilot_theta_chain}{A matrix containing the chain of parameter values
#' throughout the pilot run.}
#' \item{pilot_loglike_chain}{A vector containing the log-likelihood values
#' associated with each iteration of the pilot chain.}
#'
#' @details
#' This function runs a pilot chain to estimate the posterior mean and
#' covariance of the model parameters using a particle filter. The chain is run
#' for `pilot_m` iterations, with each iteration proposing new parameters and
#' evaluating their likelihood and prior. The chain is then used to estimate
#' the posterior mean and covariance, which are used to tune the number of
#' particles for the Particle Marginal Metropolis Hastings (PMMH) algorithm.
#'
#' @keywords internal
.run_pilot_chain <- function(y, pilot_m, pilot_n, pilot_reps, init_fn,
                             transition_fn, log_likelihood_fn,
                             log_priors, proposal_sd,
                             algorithm = c("SISAR", "SISR", "SIS"),
                             resample_fn = c(
                               "stratified", "systematic",
                               "multinomial"
                             ),
                             param_transform = NULL,
                             init_params = NULL,
                             verbose = FALSE,
                             ...) {
  num_params <- length(log_priors)
  pilot_theta_chain <- matrix(NA, nrow = pilot_m, ncol = num_params)
  colnames(pilot_theta_chain) <- names(log_priors)
  pilot_loglike_chain <- numeric(pilot_m)

  if (is.null(init_params)) {
    init_params <- rep(1, num_params)
    names(init_params) <- names(log_priors)
  }

  # Validate initial parameters using user-supplied log-priors.
  log_prior_init <- sapply(seq_along(init_params), function(i) {
    log_priors[[i]](init_params[i])
  })
  if (any(!is.finite(log_prior_init))) {
    stop(paste0(
      "Invalid initial parameters: at least one initial value is",
      "outside the support of its prior. Modify them in the argument",
      "init_params."
    ))
  }

  current_theta <- init_params
  pilot_theta_chain[1, ] <- current_theta

  # Run particle filter with current parameters.
  current_theta_list <- as.list(current_theta)
  pf_result <- do.call(particle_filter, c(
    list(
      y = y, n = pilot_n, init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      algorithm = algorithm,
      resample_fn = "stratified"
    ),
    current_theta_list,
    list(...)
  ))
  current_loglike <- pf_result$loglike
  pilot_loglike_chain[1] <- current_loglike

  # Create default transformation if none provided.
  if (is.null(param_transform)) {
    param_transform <- rep("identity", num_params)
    names(param_transform) <- names(log_priors)
  } else if (is.list(param_transform)) {
    # Ensure every parameter in log_priors has a corresponding transform.
    if (!all(names(log_priors) %in% names(param_transform))) {
      stop(paste0(
        "param_transform must include an entry for every",
        "parameter in log_priors."
      ))
    }
    # Reorder and convert the list to a vector.
    param_transform <- unlist(param_transform[names(log_priors)])

    # Validate transformations and replace any invalid entries.
    invalid <- !(param_transform %in% c("log", "identity"))
    if (any(invalid)) {
      warning(paste0(
        "Only 'log' and 'identity' transformations are supported.",
        "Using 'identity' for invalid entries."
      ))
      param_transform[invalid] <- "identity"
    }
  } else {
    stop("param_transform must be a list.")
  }

  # Reorder param_transform to match the order of log_priors
  param_transform <- as.list(unlist(param_transform[names(log_priors)]))



  for (m in 2:pilot_m) {
    valid_theta <- FALSE
    while (!valid_theta) {
      # Transform the current theta.
      current_theta_trans <- .transform_params(
        theta = current_theta,
        transform = param_transform
      )
      # Propose in the transformed space.
      proposed_theta_trans <- current_theta_trans +
        stats::rnorm(length(current_theta_trans), mean = 0, sd = proposal_sd)
      # Back-transform to original scale.
      proposed_theta <- .back_transform_params(
        theta_trans = proposed_theta_trans,
        transform = param_transform
      )

      # Compute the Jacobian adjustments.
      log_jacobian_proposed <- .compute_jacobian(
        theta = proposed_theta,
        transform = param_transform
      )
      log_jacobian_current <- .compute_jacobian(
        theta = current_theta,
        transform = param_transform
      )

      # Compute the log-priors for the proposed parameters.
      log_prior_proposed <- mapply(
        function(fn, theta) fn(theta),
        log_priors, proposed_theta
      )

      if (all(is.finite(log_prior_proposed))) {
        valid_theta <- TRUE
      }
    }

    # Compute log-priors for current parameters.
    log_prior_current <- mapply(
      function(fn, theta) fn(theta),
      log_priors, current_theta
    )

    # Run particle filter with the proposed parameters.
    proposed_theta_list <- as.list(proposed_theta)
    pf_prop <- do.call(particle_filter, c(
      list(
        y = y, n = pilot_n, init_fn = init_fn,
        transition_fn = transition_fn,
        log_likelihood_fn = log_likelihood_fn,
        algorithm = algorithm,
        resample_fn = resample_fn
      ),
      proposed_theta_list,
      list(...)
    ))
    proposed_loglike <- pf_prop$loglike

    # Acceptance ratio with Jacobian adjustments.
    log_accept_num <- sum(log_prior_proposed) + proposed_loglike +
      log_jacobian_proposed
    log_accept_denom <- sum(log_prior_current) + current_loglike +
      log_jacobian_current
    log_accept_ratio <- log_accept_num - log_accept_denom

    if (log(stats::runif(1)) < log_accept_ratio) {
      current_theta <- proposed_theta
      current_loglike <- proposed_loglike
    }

    pilot_theta_chain[m, ] <- current_theta
    pilot_loglike_chain[m] <- current_loglike
  }

  burn_in <- floor(pilot_m / 2)
  pilot_theta_post <- pilot_theta_chain[(burn_in + 1):pilot_m, , drop = FALSE]
  pilot_theta_mean <- colMeans(pilot_theta_post)

  if (ncol(pilot_theta_post) > 1) {
    pilot_theta_cov <- stats::cov(pilot_theta_post)
  } else {
    pilot_theta_cov <- stats::var(pilot_theta_post)
  }

  if (verbose) {
    message("Pilot chain posterior mean:")
    print(pilot_theta_mean)
    if (ncol(pilot_theta_post) > 1) {
      if (any(param_transform == "log")) {
        message("Pilot chain posterior covariance (on transformed space):")
      } else {
        message("Pilot chain posterior covariance:")
      }
      print(pilot_theta_cov)
    } else {
      if (any(param_transform == "log")) {
        message("Pilot chain posterior variance (on transformed space):")
      } else {
        message("Pilot chain posterior variance:")
      }
      print(pilot_theta_cov)
    }
  }

  # Run the pilot run using the estimated posterior mean.
  pilot_result <- do.call(.pilot_run, c(
    list(
      y = y, pilot_n = pilot_n, pilot_reps = pilot_reps,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      algorithm = algorithm,
      resample_fn = resample_fn
    ),
    as.list(pilot_theta_mean),
    list(...)
  ))
  target_n <- pilot_result$target_n
  message("Using ", target_n, " particles for PMMH:")

  list(
    pilot_theta_mean = pilot_theta_mean,
    pilot_theta_cov = pilot_theta_cov,
    target_n = target_n,
    pilot_theta_chain = pilot_theta_chain,
    pilot_loglike_chain = pilot_loglike_chain
  )
}
