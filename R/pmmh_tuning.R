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
#' @param ... Additional arguments passed to `init_fn_ssm`,
#' `transition_fn_ssm`, and `log_likelihood_fn_ssm`.
#'
#' @return A list containing:
#' \describe{
#'   \item{variance_estimate}{The estimated variance of the log-likelihoods
#'   from the pilot run.}
#'   \item{target_N}{The recommended target number of particles (ensuring a
#'   minimum of 100).}
#'   \item{pilot_loglikes}{A numeric vector of log-likelihood values computed
#'   during the run.}
#' }
#'
#' @details The function performs \code{pilot_reps} evaluations of the particle
#' filter using the provided parameter vector \code{theta}. It then estimates
#' the variance of the log-likelihoods and scales the initial particle number
#' by this variance. The final number of particles is taken as the ceiling of
#' the scaled value with a minimum of 100.
#'
#' @keywords internal
.pilot_run <- function(y, pilot_n, pilot_reps, init_fn_ssm,
                       transition_fn_ssm, log_likelihood_fn_ssm,
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
      n = pilot_n,
      init_fn = init_fn_ssm,
      transition_fn = transition_fn_ssm,
      log_likelihood_fn = log_likelihood_fn_ssm,
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
  target_n <- min(target_n, 2000) # Limit to 2000 particles

  list(
    variance_estimate = variance_estimate,
    target_n = target_n,
    pilot_loglikes = pilot_loglikes
  )
}

#' Run Pilot Chain for Posterior Estimation
#'
#'
#' @param y A numeric vector of observations.
#' @param pilot_m An integer specifying the number of iterations for the pilot
#' chain.
#' @param pilot_n An integer specifying the number of particles for the particle
#' filter.
#' @param pilot_reps An integer specifying the number of repetitions for the
#' pilot run.
#' @param init_fn_ssm A function to initialize the state space model.
#' @param transition_fn_ssm A function that defines the state transition
#' dynamics in the state space model.
#' @param log_likelihood_fn_ssm A function that computes the log-likelihood of
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
.run_pilot_chain <- function(y, pilot_m, pilot_n, pilot_reps, init_fn_ssm,
                             transition_fn_ssm, log_likelihood_fn_ssm,
                             log_priors, proposal_sd,
                             algorithm = c("SISAR", "SISR", "SIS"),
                             resample_fn = c(
                               "stratified", "systematic",
                               "multinomial"
                             ),
                             param_transform = NULL,
                             init_params = NULL,
                             ...) {
  cat("Running pilot chain to estimate posterior mean and covariance...\n")

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
    stop("Invalid initial parameters: at least one initial value is outside the
         support of its prior. Modify them in the argument init_params.")
  }

  current_theta <- init_params
  pilot_theta_chain[1, ] <- current_theta

  # Prepare extra arguments from current_theta as a list.
  current_theta_list <- as.list(current_theta)

  # Run particle filter with current parameters using do.call
  pf_result <- do.call(particle_filter, c(
    list(
      y = y, n = pilot_n, init_fn = init_fn_ssm,
      transition_fn = transition_fn_ssm,
      log_likelihood_fn = log_likelihood_fn_ssm,
      algorithm = algorithm,
      resample_fn = "stratified"
    ),
    current_theta_list,
    list(...)
  ))

  current_loglike <- pf_result$loglike
  pilot_loglike_chain[1] <- current_loglike

  log_jacobian_proposed <- 0
  log_jacobian_current <- 0

  # Check if the param_transform argument is valid
  if (!is.null(param_transform)) {
    invalid_transform <- which(!(param_transform %in% c("log", "identity")))
    if (length(invalid_transform) > 0) {
      warning("Currently only supports 'log' and 'identity' transformation.
              Using 'identity' (no transformation) instead.")
      param_transform[invalid_transform] <- "identity"
    }
  }

  for (m in 2:pilot_m) {
    # Proposal step with transformation (if provided)
    valid_theta <- FALSE
    while (!valid_theta) {
      proposed_theta <- sapply(seq_along(current_theta), function(i) {
        if (!is.null(param_transform) && param_transform[i] == "log") {
          # Propose on the log scale
          log_theta <- log(current_theta[i]) + stats::rnorm(1,
            mean = 0,
            sd = proposal_sd[i]
          )
          exp(log_theta) # Transform back to original scale
        } else {
          # Regular proposal (identity transformation)
          current_theta[i] + stats::rnorm(1, mean = 0, sd = proposal_sd[i])
        }
      })

      # Calculate the log-jacobian adjustment
      log_jacobian_proposed <- sum(sapply(
        seq_along(proposed_theta),
        function(i) {
          if (!is.null(param_transform) && param_transform[i] == "log") {
            -log(proposed_theta[i]) # Log transformation adjustment
          } else {
            0 # No transformation, no adjustment
          }
        }
      ))

      # Compute the log-priors for the proposed parameters
      log_prior_proposed <- mapply(
        function(fn, theta) fn(theta), log_priors,
        proposed_theta
      )
      log_prior_proposed_total <- sum(log_prior_proposed)

      # Check if proposed theta is valid (log-prior is finite)
      if (all(is.finite(log_prior_proposed))) {
        valid_theta <- TRUE
      }
    }

    # Compute the log-priors for current and proposed parameters.
    log_prior_current <- mapply(
      function(fn, theta) fn(theta), log_priors,
      current_theta
    )
    log_prior_proposed <- mapply(
      function(fn, theta) fn(theta), log_priors,
      proposed_theta
    )

    log_prior_current_total <- sum(log_prior_current)
    log_prior_proposed_total <- sum(log_prior_proposed)

    # Prepare extra arguments from proposed_theta as a list.
    proposed_theta_list <- as.list(proposed_theta)

    # Run particle filter with the proposed parameters.
    pf_prop <- do.call(particle_filter, c(
      list(
        y = y, n = pilot_n, init_fn = init_fn_ssm,
        transition_fn = transition_fn_ssm,
        log_likelihood_fn = log_likelihood_fn_ssm,
        algorithm = algorithm,
        resample_fn = resample_fn
      ),
      proposed_theta_list,
      list(...)
    ))

    proposed_loglike <- pf_prop$loglike
    log_num <- log_prior_proposed_total + proposed_loglike +
      log_jacobian_proposed
    log_denom <- log_prior_current_total + current_loglike +
      log_jacobian_current
    log_accept_ratio <- log_num - log_denom

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

  cat("Pilot chain posterior mean:\n")
  print(pilot_theta_mean)
  if (ncol(pilot_theta_post) > 1) {
    if (!is.null(param_transform) && any(param_transform == "log")) {
      cat("Pilot chain posterior covariance (on transformed space):\n")
    } else {
      cat("Pilot chain posterior covariance:\n")
    }
    print(pilot_theta_cov)
  } else {
    if (!is.null(param_transform) && any(param_transform == "log")) {
      cat("Pilot chain posterior variance (on transformed space):\n")
    } else {
      cat("Pilot chain posterior variance:\n")
    }
    print(pilot_theta_cov)
  }

  # Run the pilot run using the estimated posterior mean.
  pilot_result <- do.call(.pilot_run, c(
    list(
      y = y, pilot_n = pilot_n, pilot_reps = pilot_reps,
      init_fn_ssm = init_fn_ssm,
      transition_fn_ssm = transition_fn_ssm,
      log_likelihood_fn_ssm = log_likelihood_fn_ssm,
      algorithm = algorithm,
      resample_fn = resample_fn
    ),
    as.list(pilot_theta_mean),
    list(...)
  ))

  target_n <- pilot_result$target_n
  cat("Estimated target number of particles for PMMH:", target_n, "\n")

  list(
    pilot_theta_mean = pilot_theta_mean,
    pilot_theta_cov = pilot_theta_cov,
    target_n = target_n,
    pilot_theta_chain = pilot_theta_chain,
    pilot_loglike_chain = pilot_loglike_chain
  )
}
