#' Pilot Run for Particle Filter Tuning
#'
#' This internal function repeatedly evaluates the particle filter in order to
#' estimate the variance of the log-likelihoods and to compute a recommended
#' target number of particles for the Particle Marginal Metropolis Hastings
#' (PMMH) algorithm.
#'
#' @param y A numeric vector or time series of observations.
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
                      resample_fn = c("stratified", "systematic",
                                      "multinomial"),
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
  variance_estimate <- var(pilot_loglikes)
  target_N <- ceiling(pilot_n * variance_estimate)
  target_N <- max(target_N, 100)  # Ensure a minimum number of particles

  list(
    variance_estimate = variance_estimate,
    target_N = target_N,
    pilot_loglikes = pilot_loglikes
  )
}

.run_pilot_chain_new <- function(y, pilot_m, pilot_n, pilot_reps, init_fn_ssm,
                                 transition_fn_ssm, log_likelihood_fn_ssm,
                                 log_priors, proposal_sd,
                                 algorithm = c("SISAR", "SISR", "SIS"),
                                 resample_fn = c("stratified", "systematic",
                                                 "multinomial"),
                                 proposal_fn = NULL,
                                 init_params = NULL,
                                 ...) {
  if (is.null(proposal_fn)) {
    proposal_fn <- function(theta, proposal_sd) {
      proposed_theta <- theta + rnorm(length(theta), mean = 0, sd = proposal_sd)
      list(proposed_theta = proposed_theta, log_jacobian = 0)
    }
  }

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
    stop("Invalid initial parameters: at least one initial value is outside the support of its prior.
         Modify them in the argument init_params.")
  }

  current_theta <- init_params
  pilot_theta_chain[1, ] <- current_theta

  # Prepare extra arguments from current_theta as a list.
  current_theta_list <- as.list(current_theta)

  # Run particle filter with current parameters using do.call
  pf_result <- do.call(particle_filter, c(
    list(y = y, n = pilot_n, init_fn = init_fn_ssm,
         transition_fn = transition_fn_ssm,
         log_likelihood_fn = log_likelihood_fn_ssm,
         algorithm = algorithm,
         resample_fn = "stratified"),
    current_theta_list,
    list(...)
  ))

  current_loglike <- pf_result$loglike
  pilot_loglike_chain[1] <- current_loglike

  log_jacobian_proposed <- 0
  log_jacobian_current <- 0

  for (m in 2:pilot_m) {
    prop <- proposal_fn(current_theta, proposal_sd)
    proposed_theta <- prop$proposed_theta
    log_jacobian_proposed <- prop$log_jacobian

    # Compute the log-priors for current and proposed parameters.
    log_prior_current <- mapply(function(fn, theta) fn(theta), log_priors, current_theta)
    log_prior_proposed <- mapply(function(fn, theta) fn(theta), log_priors, proposed_theta)

    log_prior_current_total <- sum(log_prior_current)
    log_prior_proposed_total <- sum(log_prior_proposed)

    # Prepare extra arguments from proposed_theta as a list.
    proposed_theta_list <- as.list(proposed_theta)

    # Run particle filter with the proposed parameters.
    pf_prop <- do.call(particle_filter, c(
      list(y = y, n = pilot_n, init_fn = init_fn_ssm,
           transition_fn = transition_fn_ssm,
           log_likelihood_fn = log_likelihood_fn_ssm,
           algorithm = algorithm,
           resample_fn = "stratified"),
      proposed_theta_list,
      list(...)
    ))

    proposed_loglike <- pf_prop$loglike

    log_accept_ratio <- (log_prior_proposed_total + proposed_loglike + log_jacobian_proposed) -
      (log_prior_current_total + current_loglike + log_jacobian_current)

    if (log(runif(1)) < log_accept_ratio) {
      current_theta <- proposed_theta
      current_loglike <- proposed_loglike
    }

    pilot_theta_chain[m, ] <- current_theta
    pilot_loglike_chain[m] <- current_loglike
  }

  burn_in <- floor(pilot_m / 2)

  pilot_theta_post <- pilot_theta_chain[(burn_in + 1):pilot_m]

  if (is.matrix(pilot_theta_post)) {
    pilot_theta_mean <- colMeans(pilot_theta_post)
    pilot_theta_cov <- cov(pilot_theta_post)
  } else {
    pilot_theta_mean <- mean(pilot_theta_post)
    pilot_theta_cov <- var(pilot_theta_post)
  }

  cat("Pilot chain posterior mean:\n")
  print(pilot_theta_mean)
  if (is.matrix(pilot_theta_post)) {
    cat("Pilot chain posterior covariance:\n")
    print(pilot_theta_cov)
  } else {
    cat("Pilot chain posterior variance:\n")
    print(pilot_theta_cov)
  }

  # Prepare the parameter list from the pilot mean to pass to .pilot_run.
  pilot_theta_mean_list <- as.list(pilot_theta_mean)

  # Run the pilot run using the estimated posterior mean.
  pilot_result <- do.call(.pilot_run, c(
    list(y = y, pilot_n = pilot_n, pilot_reps = pilot_reps,
         init_fn_ssm = init_fn_ssm,
         transition_fn_ssm = transition_fn_ssm,
         log_likelihood_fn_ssm = log_likelihood_fn_ssm,
         algorithm = algorithm,
         resample_fn = resample_fn),
    pilot_theta_mean_list,
    list(...)
  ))

  target_N <- pilot_result$target_N
  cat("Estimated target number of particles for PMMH:", target_N, "\n")

  list(
    pilot_theta_mean = pilot_theta_mean,
    pilot_theta_cov = pilot_theta_cov,
    target_N = target_N,
    pilot_theta_chain = pilot_theta_chain,
    pilot_loglike_chain = pilot_loglike_chain
  )
}


#' Run Pilot Chain for PMMH Tuning
#'
#' This internal function executes a pilot chain using the Particle Marginal
#' Metropolis Hastings (PMMH) algorithm to estimate the posterior mean and
#' covariance of the model parameters. It then uses these estimates
#' to compute a recommended target number of particles for the full PMMH run.
#'
#' @param y A numeric vector or time series of observations.
#' @param pilot_M An integer specifying the number of iterations for the pilot
#' chain.
#' @param pilot_N An integer specifying the number of particles to use in each
#' particle filter evaluation.
#' @param init_theta A numeric vector (length 3) of initial parameter values.
#' @param prop_sd A numeric vector (length 3) of proposal standard deviations
#' for the parameters.
#'
#' @return A list containing:
#' \describe{
#'   \item{pilot_theta_mean}{A numeric vector representing the posterior mean
#'   of the parameters from the pilot chain.}
#'   \item{pilot_theta_cov}{The covariance matrix of the parameters estimated
#'   from the pilot chain.}
#'   \item{target_N}{The estimated target number of particles for the PMMH
#'   algorithm.}
#'   \item{pilot_theta_chain}{A matrix containing the chain of parameter values
#'   (each row corresponds to an iteration).}
#'   \item{pilot_loglike_chain}{A numeric vector containing the log-likelihood
#'   values along the chain.}
#' }
#'
#' @details The function starts by initializing the pilot chain with
#' \code{init_theta} and runs a total of \code{pilot_M} iterations. At each
#' iteration, it proposes new parameter values using random normal perturbations
#' (with log-transform adjustments for scale parameters), computes the
#' acceptance probability, and either accepts or rejects the new parameters.
#' The first half of the chain is discarded as burn-in. Finally, the posterior
#' mean from the pilot chain is used in a call to \code{pilot_run} to estimate
#' a target number of particles.
#'
#' @keywords internal
.run_pilot_chain <- function(y, pilot_M, pilot_N, init_theta, prop_sd,
                             init_fn_ssm, transition_fn_ssm,
                             log_likelihood_fn_ssm, log_prior_fn,
                             algorithm = c("SISAR", "SISR", "SIS"),
                             resample_fn = c("stratified", "systematic",
                                             "multinomial"),
                             pilot_reps = 10,
                             proposal_fn = NULL,
                             ...) {
  algorithm <- match.arg(algorithm)
  resample_fn <- match.arg(resample_fn)

  # Default proposal: additive Gaussian random walk on the current scale.
  # Returns a list with proposed_theta and (optionally) a log_jacobian
  # adjustment.
  if (is.null(proposal_fn)) {
    proposal_fn <- function(theta, prop_sd) {
      proposed_theta <- theta + rnorm(length(theta), mean = 0, sd = prop_sd)
      list(proposed_theta = proposed_theta, log_jacobian = 0)
    }
  }

  cat("Running pilot chain to estimate posterior mean and covariance...\n")
  num_params <- length(init_theta)
  pilot_theta_chain <- matrix(NA, nrow = pilot_M, ncol = num_params)
  colnames(pilot_theta_chain) <- names(init_theta)
  pilot_loglike_chain <- numeric(pilot_M)

  current_theta <- init_theta
  pilot_theta_chain[1, ] <- current_theta

  # Initial particle filter call using current_theta (passed as named arguments)
  current_pf_call <- c(
    list(y = y,
         n = pilot_N,
         init_fn = init_fn_ssm,
         transition_fn = transition_fn_ssm,
         log_likelihood_fn = log_likelihood_fn_ssm,
         algorithm = algorithm,
         resample_fn = resample_fn),
    as.list(current_theta),
    list(...)
  )
  pf_result <- do.call(particle_filter, current_pf_call)
  current_loglike <- pf_result$loglike
  pilot_loglike_chain[1] <- current_loglike

  for (m in 2:pilot_M) {
    # Propose new parameters using the proposal function
    prop <- proposal_fn(current_theta, prop_sd)
    proposed_theta <- prop$proposed_theta
    log_jacobian_ratio <- prop$log_jacobian  # typically 0 for symmetric proposals

    # Compute log priors (user-supplied function that handles arbitrary θ)
    log_prior_current <- log_prior_fn(current_theta)
    log_prior_proposed <- log_prior_fn(proposed_theta)

    # Particle filter call for proposed parameters
    pf_prop_call <- c(
      list(y = y,
           n = pilot_N,
           init_fn = init_fn_ssm,
           transition_fn = transition_fn_ssm,
           log_likelihood_fn = log_likelihood_fn_ssm,
           algorithm = algorithm,
           resample_fn = resample_fn),
      as.list(proposed_theta),
      list(...)
    )
    pf_prop <- do.call(particle_filter, pf_prop_call)
    proposed_loglike <- pf_prop$loglike

    # Calculate acceptance ratio on the log-scale
    log_num <- log_prior_proposed + proposed_loglike
    log_denom <- log_prior_current + current_loglike
    log_accept_ratio <- log_num - log_denom + log_jacobian_ratio

    if (log(runif(1)) < log_accept_ratio) {
      current_theta <- proposed_theta
      current_loglike <- proposed_loglike
    }

    pilot_theta_chain[m, ] <- current_theta
    pilot_loglike_chain[m] <- current_loglike
  }

  # Discard the first half of the chain as burn-in
  burn_in <- floor(pilot_M / 2)
  pilot_theta_post <- pilot_theta_chain[(burn_in + 1):pilot_M, , drop = FALSE]

  if (is.matrix(pilot_theta_post)) {
    pilot_theta_mean <- colMeans(pilot_theta_post)
    pilot_theta_cov <- cov(pilot_theta_post)
  } else {
    pilot_theta_mean <- mean(pilot_theta_post)
    pilot_theta_cov <- var(pilot_theta_post)
  }

  cat("Pilot chain posterior mean:\n")
  print(pilot_theta_mean)
  if (is.matrix(pilot_theta_post)) {
    cat("Pilot chain posterior covariance:\n")
    print(pilot_theta_cov)
  } else {
    cat("Pilot chain posterior variance:\n")
    print(pilot_theta_cov)
  }

  # Estimate the target number of particles using the pilot run.
  # Pass the posterior mean θ to .pilot_run via do.call.
  pilot_run_call <- c(
    list(y = y,
         pilot_n = pilot_N,   # pilot_N is used as the number of particles in .pilot_run
         pilot_reps = pilot_reps,
         init_fn_ssm = init_fn_ssm,
         transition_fn_ssm = transition_fn_ssm,
         log_likelihood_fn_ssm = log_likelihood_fn_ssm),
    as.list(pilot_theta_mean),
    list(...)
  )
  pilot_result <- do.call(.pilot_run, pilot_run_call)
  target_N <- pilot_result$target_N
  cat("Estimated target number of particles for PMMH:", target_N, "\n")

  list(
    pilot_theta_mean = pilot_theta_mean,
    pilot_theta_cov = pilot_theta_cov,
    target_N = target_N,
    pilot_theta_chain = pilot_theta_chain,
    pilot_loglike_chain = pilot_loglike_chain
  )
}

