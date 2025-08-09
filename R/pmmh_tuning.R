#' Pilot Run for Particle Filter Tuning
#'
#' This internal function repeatedly evaluates the particle filter in order to
#' estimate the variance of the log-likelihoods and to compute a recommended
#' target number of particles for the Particle Marginal Metropolis Hastings
#' (PMMH) algorithm.
#'
#' @inheritParams particle_filter
#' @param pilot_n An integer specifying the initial number of particles to use.
#' @param pilot_reps An integer specifying the number of repetitions for the
#' pilot run.
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
#' the scaled value with a minimum of 50 and a maximum of 1000.
#'
#' @keywords internal
.pilot_run <- function(
  pf_wrapper,
  y, pilot_n, pilot_reps,
  init_fn, transition_fn, log_likelihood_fn,
  obs_times = NULL,
  resample_fn = NULL,
  ... # any extra args for the filter
) {
  pilot_loglikes <- numeric(pilot_reps)

  for (i in seq_len(pilot_reps)) {
    pf_result <- pf_wrapper(
      y = y,
      num_particles = pilot_n,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      obs_times = obs_times,
      resample_fn = resample_fn,
      return_particles = FALSE,
      ...
    )
    pilot_loglikes[i] <- pf_result$loglike
  }

  variance_estimate <- var(pilot_loglikes)
  target_n <- ceiling(pilot_n * variance_estimate)
  target_n <- max(target_n, 50)
  target_n <- min(target_n, 1000)

  list(
    variance_estimate = variance_estimate,
    target_n = target_n,
    pilot_loglikes = pilot_loglikes
  )
}



#' Run Pilot Chain for Posterior Estimation
#'
#' @inheritParams particle_filter
#' @param pilot_m An integer specifying the number of iterations for the pilot
#' chain.
#' @param pilot_n An integer specifying the number of particles for the particle
#' filter.
#' @param pilot_reps An integer specifying the number of repetitions for the
#' pilot run.
#' @param log_priors A list of functions representing the log-priors for each
#' model parameter.
#' @param proposal_sd A numeric vector specifying the standard deviations for
#' the random walk proposal distribution for each parameter.
#' @param param_transform A character vector specifying the parameter
#' transformations when proposing parameters using a random walk.
#' Currently only supports "log" for log-transformation, "logit" for logit
#' transformation, and "identity" for no transformation. Default is `NULL`,
#' which correspond to no transformation ("identity).
#' @param pilot_init_params A numeric vector of initial parameter values.
#' If `NULL`, it will default to a vector of ones. Default is `NULL`.
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
#' @importFrom stats rnorm runif var cov
#' @keywords internal
.run_pilot_chain <- function(
    pf_wrapper, # particle filter wrapper to use
    y,
    pilot_m, pilot_n, pilot_reps,
    init_fn, transition_fn, log_likelihood_fn,
    log_priors, proposal_sd,
    obs_times = NULL,
    param_transform = NULL,
    pilot_init_params = NULL,
    verbose = FALSE,
    ...) {
  num_params <- length(log_priors)
  pilot_theta_chain <- matrix(NA, nrow = pilot_m, ncol = num_params)
  colnames(pilot_theta_chain) <- names(log_priors)
  pilot_loglike_chain <- numeric(pilot_m)

  # Default initial parameters
  if (is.null(pilot_init_params)) {
    pilot_init_params <- rep(1, num_params)
    names(pilot_init_params) <- names(log_priors)
  }

  # Validate initial parameters
  log_prior_init <- sapply(seq_along(pilot_init_params), function(i) {
    log_priors[[i]](pilot_init_params[i])
  })
  if (any(!is.finite(log_prior_init))) {
    stop(paste0(
      "Initial parameter values are invalid: some lie outside the prior ",
      "support. Please provide valid starting values via pilot_init_params."
    ))
  }


  current_theta <- pilot_init_params
  pilot_theta_chain[1, ] <- current_theta

  # Initial log-likelihood
  current_theta_list <- as.list(current_theta)
  pf_result <- do.call(pf_wrapper, c(
    list(
      y = y,
      num_particles = pilot_n,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      obs_times = obs_times,
      return_particles = FALSE
    ),
    current_theta_list,
    list(...)
  ))
  current_loglike <- pf_result$loglike
  pilot_loglike_chain[1] <- current_loglike

  # Set default transforms
  if (is.null(param_transform)) {
    param_transform <- rep("identity", num_params)
    names(param_transform) <- names(log_priors)
  } else if (is.list(param_transform)) {
    if (!all(names(log_priors) %in% names(param_transform))) {
      stop("param_transform must include every parameter in log_priors.")
    }
    param_transform <- unlist(param_transform[names(log_priors)])
    invalid <- !(param_transform %in% c("log", "logit", "identity"))
    if (any(invalid)) {
      warning(paste0(
        "Only 'log', 'logit', 'identity' supported.\n",
        "Invalid replaced with 'identity'."
      ))
      param_transform[invalid] <- "identity"
    }
  } else {
    stop("param_transform must be a list.")
  }

  param_transform <- as.list(param_transform[names(log_priors)])

  # Pilot MCMC loop
  for (m in 2:pilot_m) {
    valid_theta <- FALSE
    while (!valid_theta) {
      # Transform and propose
      current_theta_trans <- .transform_params(current_theta, param_transform)
      proposed_theta_trans <- current_theta_trans +
        rnorm(length(current_theta_trans), mean = 0, sd = proposal_sd)
      proposed_theta <- .back_transform_params(
        proposed_theta_trans,
        param_transform
      )

      # Check priors
      log_prior_proposed <- mapply(
        function(fn, val) fn(val), log_priors, proposed_theta
      )
      if (all(is.finite(log_prior_proposed))) valid_theta <- TRUE
    }

    # Current priors
    log_prior_current <- mapply(
      function(fn, val) fn(val), log_priors, current_theta
    )

    # Proposed log-likelihood
    proposed_theta_list <- as.list(proposed_theta)
    pf_prop <- do.call(pf_wrapper, c(
      list(
        y = y,
        num_particles = pilot_n,
        init_fn = init_fn,
        transition_fn = transition_fn,
        log_likelihood_fn = log_likelihood_fn,
        obs_times = obs_times,
        return_particles = FALSE
      ),
      proposed_theta_list,
      list(...)
    ))
    proposed_loglike <- pf_prop$loglike

    # Jacobians
    log_jacobian_proposed <- .compute_log_jacobian(
      proposed_theta,
      param_transform
    )
    log_jacobian_current <- .compute_log_jacobian(
      current_theta,
      param_transform
    )

    # MH acceptance
    log_accept_num <- sum(log_prior_proposed) + proposed_loglike +
      log_jacobian_proposed
    log_accept_denom <- sum(log_prior_current) + current_loglike +
      log_jacobian_current
    log_accept_ratio <- log_accept_num - log_accept_denom
    if (is.na(log_accept_ratio)) log_accept_ratio <- -Inf

    if (log(runif(1)) < log_accept_ratio) {
      current_theta <- proposed_theta
      current_loglike <- proposed_loglike
    }

    pilot_theta_chain[m, ] <- current_theta
    pilot_loglike_chain[m] <- current_loglike
  }

  # Posterior summaries
  burn_in <- floor(pilot_m / 2)
  pilot_theta_post <- pilot_theta_chain[(burn_in + 1):pilot_m, , drop = FALSE]
  pilot_theta_mean <- colMeans(pilot_theta_post)
  pilot_theta_cov <- if (ncol(pilot_theta_post) > 1) {
    cov(pilot_theta_post)
  } else {
    var(pilot_theta_post)
  }

  if (verbose) {
    message("Pilot chain posterior mean:")
    print(pilot_theta_mean)
    if (ncol(pilot_theta_post) > 1) {
      msg <- if (any(param_transform != "identity")) {
        "Pilot chain posterior covariance (transformed space):"
      } else {
        "Pilot chain posterior covariance:"
      }
      message(msg)
      print(pilot_theta_cov)
    } else {
      msg <- if (any(param_transform != "identity")) {
        "Pilot chain posterior variance (transformed space):"
      } else {
        "Pilot chain posterior variance:"
      }
      message(msg)
      print(pilot_theta_cov)
    }
  }

  # Pilot run for target_n
  pilot_result <- do.call(.pilot_run, c(
    list(
      pf_wrapper = pf_wrapper,
      y = y,
      pilot_n = pilot_n,
      pilot_reps = pilot_reps,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      obs_times = obs_times
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
