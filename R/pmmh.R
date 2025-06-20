#' Create Tuning Control Parameters
#'
#' This function creates a list of tuning parameters used by the
#' \code{\link{pmmh}} function. The tuning choices are inspired by Pitt et al.
#' [2012] and Dahlin and Schön [2019].
#'
#' @param pilot_proposal_sd Standard deviation for pilot proposals. Default is
#' 0.5.
#' @param pilot_n Number of pilot particles for particle filter. Default is 100.
#' @param pilot_m Number of iterations for MCMC. Default is 2000.
#' @param pilot_target_var The target variance for the posterior log-likelihood
#' evaluated at estimated posterior mean. Default is 1.
#' @param pilot_burn_in Number of burn-in iterations for MCMC. Default is 500.
#' @param pilot_reps Number of times a particle filter is run. Default is 100.
#' @param pilot_algorithm The algorithm used for the pilot particle filter.
#' Default is "SISAR".
#' @param pilot_resample_fn The resampling function used for the pilot particle
#' filter. Default is "stratified".
#'
#' @return A list of tuning control parameters.
#'
#' @references M. K. Pitt, R. d. S. Silva, P. Giordani, and R. Kohn.
#' On some properties of Markov chain Monte Carlo simulation methods based on
#' the particle filter. Journal of Econometrics, 171(2):134–151, 2012.
#' doi: https://doi.org/10.1016/j.jeconom.2012.06.004
#'
#' J. Dahlin and T. B. Schön. Getting started with particle
#' Metropolis-Hastings for inference in nonlinear dynamical models. Journal
#' of Statistical Software, 88(2):1–41, 2019. doi: 10.18637/jss.v088.c02
#'
#' @export
default_tune_control <- function(
    pilot_proposal_sd = 0.5, pilot_n = 100, pilot_m = 2000,
    pilot_target_var = 1, pilot_burn_in = 500, pilot_reps = 100,
    pilot_algorithm = c("SISAR", "SISR", "SIS"),
    pilot_resample_fn = c("stratified", "systematic", "multinomial")) {
  if (!is.numeric(pilot_proposal_sd) || pilot_proposal_sd <= 0) {
    stop("pilot_proposal_sd must be a positive numeric value.")
  }
  if (!is.numeric(pilot_n) || pilot_n <= 0) {
    stop("pilot_n must be a positive numeric value.")
  }
  if (!is.numeric(pilot_m) || pilot_m <= 0) {
    stop("pilot_m must be a positive numeric value.")
  }
  if (!is.numeric(pilot_target_var) || pilot_target_var <= 0) {
    stop("pilot_target_var must be a positive numeric value.")
  }
  if (!is.numeric(pilot_burn_in) || pilot_burn_in <= 0) {
    stop("pilot_burn_in must be a positive numeric value.")
  }
  pilot_algorithm <- match.arg(pilot_algorithm)
  pilot_resample_fn <- match.arg(pilot_resample_fn)
  list(
    pilot_proposal_sd = pilot_proposal_sd,
    pilot_n = pilot_n,
    pilot_m = pilot_m,
    pilot_target_var = pilot_target_var,
    pilot_burn_in = pilot_burn_in,
    pilot_reps = pilot_reps,
    pilot_algorithm = pilot_algorithm,
    pilot_resample_fn = pilot_resample_fn
  )
}

#' Particle Marginal Metropolis-Hastings (PMMH) for State-Space Models
#'
#' This function implements a Particle Marginal Metropolis-Hastings (PMMH)
#' algorithm to perform Bayesian inference in state-space models. It first
#' runs a pilot chain to tune the proposal distribution and the number of
#' particles for the particle filter, and then runs the main PMMH chain.
#'
#' @inheritParams particle_filter
#' @param m An integer specifying the total number of MCMC iterations.
#' @param log_priors A list of functions for computing the log-prior of each
#' parameter.
#' @param pilot_init_params A list of initial parameter values. Should be a list
#' of length \code{num_chains} where each element is a named vector of initial
#' parameter values.
#' @param burn_in An integer indicating the number of initial MCMC iterations
#' to discard as burn-in.
#' @param num_chains An integer specifying the number of PMMH chains to run.
#' @param param_transform An optional character vector that specifies the
#' transformation applied to each parameter before proposing. The proposal is
#' made using a multivariate normal distribution on the transformed scale.
#' Parameters are then mapped back to their original scale before evaluation.
#' Currently supports \code{"log"}, \code{"logit"}, and \code{"identity"}.
#' If \code{NULL}, the \code{"identity"} transformation is used for all
#' parameters.
#' @param tune_control A list of pilot tuning controls
#' (e.g., \code{pilot_m}, \code{pilot_reps}).
#' See \code{\link{default_tune_control}}.
#' @param verbose A logical value indicating whether to print information about
#' pilot_run tuning. Defaults to \code{FALSE}.
#' @param return_latent_state_est A logical value indicating whether to return
#' the latent state estimates for each time step. Defaults to \code{FALSE}.
#' @param seed An optional integer to set the seed for reproducibility.
#' @param num_cores An integer specifying the number of cores to use for
#' parallel processing. Defaults to 1. Each chain is assigned to its own core,
#' so the number of cores cannot exceed the number of chains
#' (\code{num_chains}). The progress information given to user is limited if
#' using more than one core.
#'
#' @details The PMMH algorithm is essentially a Metropolis Hastings algorithm
#' where instead of using the intractable marginal likelihood
#' \eqn{p(y_{1:T}\mid \theta)} it instead uses the estimated likelihood using
#' a particle filter (see also \code{\link{particle_filter}}). Values are
#' proposed using a multivariate normal distribution in the transformed space.
#' The proposal covariance and the number of particles is chosen based on a
#' pilot run. The minimum number of particles is chosen as 50 and maximum as
#' 1000.
#'
#' @return A list containing:

#' \describe{
#'   \item{\code{theta_chain}}{A dataframe of post burn-in parameter samples.}
#'   \item{\code{latent_state_chain}}{If \code{return_latent_state_est} is
#'   \code{TRUE}, a list of matrices containing the latent state estimates
#'   for each time step.}
#'   \item{\code{diagnostics}}{Diagnostics containing ESS and Rhat
#'   for each parameter (see \code{\link{ess}} and \code{\link{rhat}} for
#'   documentation).}
#' }
#'
#' @references Andrieu et al. (2010). Particle Markov chain Monte Carlo methods.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology),
#' 72(3):269–342. doi: 10.1111/j.1467-9868.2009.00736.x
#'
#' @importFrom stats rnorm dnorm runif dexp
#' @importFrom MASS mvrnorm
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples
#' init_fn <- function(num_particles) {
#'   rnorm(num_particles, mean = 0, sd = 1)
#' }
#' transition_fn <- function(particles, phi, sigma_x) {
#'   phi * particles + sin(particles) +
#'     rnorm(length(particles), mean = 0, sd = sigma_x)
#' }
#' log_likelihood_fn <- function(y, particles, sigma_y) {
#'   dnorm(y, mean = cos(particles), sd = sigma_y, log = TRUE)
#' }
#' log_prior_phi <- function(phi) {
#'   dnorm(phi, mean = 0, sd = 1, log = TRUE)
#' }
#' log_prior_sigma_x <- function(sigma) {
#'   dexp(sigma, rate = 1, log = TRUE)
#' }
#' log_prior_sigma_y <- function(sigma) {
#'   dexp(sigma, rate = 1, log = TRUE)
#' }
#' log_priors <- list(
#'   phi = log_prior_phi,
#'   sigma_x = log_prior_sigma_x,
#'   sigma_y = log_prior_sigma_y
#' )
#' # Generate data
#' t_val <- 10
#' x <- numeric(t_val)
#' y <- numeric(t_val)
#' phi <- 0.8
#' sigma_x <- 1
#' sigma_y <- 0.5
#'
#' init_state <- rnorm(1, mean = 0, sd = 1)
#' x[1] <- phi * init_state + sin(init_state) + rnorm(1, mean = 0, sd = sigma_x)
#' y[1] <- x[1] + rnorm(1, mean = 0, sd = sigma_y)
#' for (t in 2:t_val) {
#'   x[t] <- phi * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma_x)
#'   y[t] <- cos(x[t]) + rnorm(1, mean = 0, sd = sigma_y)
#' }
#' x <- c(init_state, x)
#'
#' # Should use much higher MCMC iterations in practice (m)
#' pmmh_result <- pmmh(
#'   y = y,
#'   m = 1000,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn,
#'   log_priors = log_priors,
#'   pilot_init_params = list(
#'     c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
#'     c(phi = 1, sigma_x = 0.5, sigma_y = 1)
#'   ),
#'   burn_in = 100,
#'   num_chains = 2,
#'   param_transform = list(
#'     phi = "identity",
#'     sigma_x = "log",
#'     sigma_y = "log"
#'   ),
#'   tune_control = default_tune_control(pilot_m = 500, pilot_burn_in = 100)
#' )
#' # Convergence warning is expected with such low MCMC iterations.
#'
#' # Suppose we have data for t=1,2,3,5,6,7,8,9,10 (i.e., missing at t=4)
#'
#' obs_times <- c(1, 2, 3, 5, 6, 7, 8, 9, 10)
#' y <- y[obs_times]
#'
#' # Specify observation times in the pmmh using obs_times
#' pmmh_result <- pmmh(
#'   y = y,
#'   m = 1000,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn,
#'   log_priors = log_priors,
#'   pilot_init_params = list(
#'     c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
#'     c(phi = 1, sigma_x = 0.5, sigma_y = 1)
#'   ),
#'   burn_in = 100,
#'   num_chains = 2,
#'   obs_times = obs_times,
#'   param_transform = list(
#'     phi = "identity",
#'     sigma_x = "log",
#'     sigma_y = "log"
#'   ),
#'   tune_control = default_tune_control(pilot_m = 500, pilot_burn_in = 100)
#' )
pmmh <- function(y, m, init_fn, transition_fn, log_likelihood_fn,
                 log_priors, pilot_init_params, burn_in, num_chains = 4,
                 obs_times = NULL,
                 algorithm = c("SISAR", "SISR", "SIS"),
                 resample_fn = c("stratified", "systematic", "multinomial"),
                 param_transform = NULL,
                 tune_control = default_tune_control(),
                 verbose = FALSE,
                 return_latent_state_est = FALSE,
                 seed = NULL,
                 num_cores = 1) {
  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    seed <- sample.int(.Machine$integer.max, 1) # Random seed if not provided
  }
  # ---------------------------
  # Input validation
  # ---------------------------
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.numeric(m) || m <= 0) stop("m must be a positive integer")
  if (!is.numeric(burn_in) || burn_in < 0) {
    stop("burn_in must be a positive integer")
  }
  if (burn_in >= m) {
    stop("burn_in must be smaller than the number of MCMC iterations (m)")
  }
  if (!is.numeric(num_chains) || num_chains <= 0) {
    stop("num_chains must be a positive integer")
  }
  if (!is.numeric(num_cores) || num_cores < 1) {
    stop("num_cores must be a positive integer")
  }
  if (num_cores > num_chains) {
    message("num_cores exceeds num_chains; setting num_cores to num_chains")
    num_cores <- num_chains
  }

  # Check pilot_init_params is a list
  if (!is.list(pilot_init_params)) {
    stop("pilot_init_params must be a list.")
  }

  # Check pilot_init_params all have same length
  lengths <- vapply(pilot_init_params, length, FUN.VALUE = integer(1))
  all_same <- all(lengths == lengths[1])
  if (!all_same) {
    stop("pilot_init_params must be a list of vectors of the same length.")
  }

  # Check pilot_init_params list of length num_chains
  if (length(lengths) != num_chains) {
    stop("pilot_init_params must be a list of length num_chains.")
  }

  # Check pilot_init_params all have same names
  if (!all(sapply(pilot_init_params, function(x) {
    all(names(x) == names(pilot_init_params[[1]]))
  }))) {
    stop("pilot_init_params must have the same parameter names.")
  }

  .check_params_match(
    init_fn, transition_fn, log_likelihood_fn,
    pilot_init_params[[1]], log_priors
  )

  algorithm <- match.arg(algorithm)
  resample_fn <- match.arg(resample_fn)

  num_params <- lengths[1]
  if (num_params == 0) {
    stop("pilot_init_params must contain at least one parameter.")
  }

  # Create default transformation if none provided.
  if (is.null(param_transform)) {
    param_transform <- rep("identity", num_params)
    names(param_transform) <- names(log_priors)
  } else if (is.list(param_transform)) {
    # Ensure every parameter in log_priors has a corresponding transform.
    if (!all(names(log_priors) %in% names(param_transform))) {
      stop(paste0(
        "param_transform must include an entry for every ",
        "parameter in log_priors."
      ))
    }

    # Validate transformations and replace any invalid entries.
    invalid <- !(param_transform %in% c("log", "logit", "identity"))
    if (any(invalid)) {
      warning(paste0(
        "Only 'log', 'logit', and 'identity' transformations are supported.",
        " Using 'identity' for invalid entries."
      ))
      param_transform[invalid] <- "identity"
    }
  } else {
    stop("param_transform must be a list.")
  }

  # Reorder param_transform to match the order of log_priors
  param_transform <- as.list(unlist(param_transform[names(log_priors)]))

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

  tune_control$pilot_proposal_sd <- rep(
    tune_control$pilot_proposal_sd,
    length.out = num_params
  )

  # ---------------------------
  # Define inner function for a single chain run
  # ---------------------------
  chain_result <- function(chain_index, seed) {
    set.seed(seed)
    message("Running chain ", chain_index, "...")

    # ---------------------------
    # Step 1: Run the pilot (particle) chain for tuning
    # ---------------------------
    message("Running pilot chain for tuning...")
    pilot_chain <- .run_pilot_chain(
      y = y,
      pilot_m = tune_control$pilot_m,
      pilot_n = tune_control$pilot_n,
      pilot_reps = tune_control$pilot_reps,
      init_fn = init_fn,
      transition_fn = transition_fn,
      log_likelihood_fn = log_likelihood_fn,
      log_priors = log_priors,
      proposal_sd = tune_control$pilot_proposal_sd,
      obs_times = obs_times,
      pilot_init_params = pilot_init_params[[chain_index]],
      algorithm = tune_control$pilot_algorithm,
      resample_fn = tune_control$pilot_resample_fn,
      param_transform = param_transform,
      verbose = verbose
    )

    init_theta <- pilot_chain$pilot_theta_mean
    proposal_cov <- pilot_chain$pilot_theta_cov
    target_n <- pilot_chain$target_n

    # Precompute the transformed proposal covariance:
    scale_vec <- sapply(seq_along(init_theta), function(j) {
      if (param_transform[j] == "log") {
        1 / init_theta[j] # dz/dtheta = 1 / theta
      } else if (param_transform[j] == "logit") {
        # dz/dtheta = 1 / (theta * (1 - theta))
        1 / (init_theta[j] * (1 - init_theta[j]))
      } else {
        1
      }
    })
    # Cov(z)=J Cov(theta) J^T, where J=diag(dz/dtheta)
    proposal_cov_trans <- diag(scale_vec) %*% proposal_cov %*% diag(scale_vec)

    # ---------------------------
    # Step 2: Run the PMMH chain using the tuned settings
    # ---------------------------

    message("Running Particle MCMC chain with tuned settings...")

    current_theta <- init_theta
    theta_chain <- matrix(NA, nrow = m, ncol = num_params)
    colnames(theta_chain) <- names(current_theta)
    state_est_chain <- vector("list", m)

    # Evaluate the particle filter at the initial parameter value.
    pf_result <- do.call(particle_filter, c(
      list(
        y = y,
        n = target_n,
        init_fn = init_fn,
        transition_fn = transition_fn,
        log_likelihood_fn = log_likelihood_fn,
        obs_times = obs_times,
        algorithm = algorithm,
        resample_fn = resample_fn,
        return_particles = FALSE
      ),
      as.list(current_theta)
    ))
    current_loglike <- pf_result$loglike
    current_state_est <- pf_result$state_est

    theta_chain[1, ] <- current_theta
    state_est_chain[[1]] <- current_state_est

    for (i in 2:m) {
      # --- Propose in the transformed space ---
      current_theta_trans <- .transform_params(current_theta, param_transform)
      proposed_theta_trans <- mvrnorm(
        n = 1, mu = current_theta_trans,
        Sigma = proposal_cov_trans
      )
      proposed_theta <- .back_transform_params(
        proposed_theta_trans,
        param_transform
      )

      # --- Check the validity of proposed parameters ---
      log_prior_proposed <- sapply(seq_along(proposed_theta), function(j) {
        log_priors[[j]](proposed_theta[j])
      })
      if (any(!is.finite(log_prior_proposed))) {
        theta_chain[i, ] <- current_theta
        state_est_chain[[i]] <- current_state_est
        next
      }

      # Run the particle filter for the proposed parameters.
      pf_proposed <- do.call(particle_filter, c(
        list(
          y = y,
          n = target_n,
          init_fn = init_fn,
          transition_fn = transition_fn,
          log_likelihood_fn = log_likelihood_fn,
          obs_times = obs_times,
          algorithm = algorithm,
          resample_fn = resample_fn,
          return_particles = FALSE
        ),
        as.list(proposed_theta)
      ))
      proposed_loglike <- pf_proposed$loglike

      # --- Compute the Jacobian adjustments ---
      log_jacobian_proposed <- .compute_log_jacobian(
        theta = proposed_theta,
        transform = param_transform
      )

      log_jacobian_current <- .compute_log_jacobian(
        theta = current_theta,
        transform = param_transform
      )

      # --- Compute the acceptance ratio ---
      log_prior_current <- sapply(seq_along(current_theta), function(j) {
        log_priors[[j]](current_theta[j])
      })
      log_accept_num <- (proposed_loglike + sum(log_prior_proposed) +
                           log_jacobian_proposed)
      log_accept_denom <- (current_loglike + sum(log_prior_current) +
                             log_jacobian_current)
      log_accept_ratio <- log_accept_num - log_accept_denom

      # If it’s NA/NaN, force it to -Inf.
      if (is.na(log_accept_ratio)) {
        log_accept_ratio <- -Inf
      }

      if (log(runif(1)) < log_accept_ratio) {
        current_theta <- proposed_theta
        current_loglike <- proposed_loglike
        current_state_est <- pf_proposed$state_est
      }

      theta_chain[i, ] <- current_theta
      state_est_chain[[i]] <- current_state_est
    }
    list(theta_chain = theta_chain,
         state_est_chain = state_est_chain)
  }

  # ---------------------------
  # Run chains (in parallel if more than one core is requested)
  # ---------------------------
  # Generate seeds
  seeds <- sample.int(.Machine$integer.max, num_chains)
  if (num_cores > 1) {
    tryCatch({
      # Execute the future parallel code
      chain_results <- future.apply::future_lapply(
        1:num_chains,
        function(i) chain_result(i, seeds[i]),
        future.seed = NULL
      )
    }, error = function(e) {
      message("An error occurred: ", e$message)
    }, finally = {
      # Ensure that the plan is reset to sequential
      future::plan(future::sequential)
    })
  } else {
    chain_results <- lapply(
      1:num_chains,
      function(i) chain_result(i, seeds[i])
    )
  }

  # Unpack the results from each chain
  theta_chains <- lapply(chain_results, function(res) res$theta_chain)
  state_est_chains <- lapply(chain_results, function(res) res$state_est_chain)

  # ---------------------------
  # Step 3: Post-processing - discard burn-in and compute latent state estimate.
  # ---------------------------
  theta_chain_post <- lapply(
    theta_chains, function(chain) chain[(burn_in + 1):m, , drop = FALSE]
  )
  state_est_chain_post <- lapply(
    state_est_chains, function(chain) chain[(burn_in + 1):m]
  )

  # ---------------------------
  # Step 4: Compute diagnostics (ESS and Rhat) for each parameter.
  # ---------------------------

  separate_parameters <- function(theta_chain_post) {
    param_names <- colnames(theta_chain_post[[1]])
    result <- list()

    for (param in param_names) {
      param_combined <- do.call(
        cbind, lapply(seq_along(theta_chain_post), function(i) {
          chain_data <- theta_chain_post[[i]]
          param_data <- chain_data[, param, drop = FALSE]
          colnames(param_data) <- paste(param, "chain", i, sep = "_")
          param_data
        })
      )
      result[[param]] <- param_combined
    }

    result
  }

  theta_chain_per_param <- separate_parameters(theta_chain_post)
  param_ess <- list()
  param_rhat <- list()

  num_chains <- length(theta_chain_post)
  ess_message_shown <- FALSE

  for (param in names(theta_chain_per_param)) {
    param_chain <- theta_chain_per_param[[param]]

    if (num_chains > 1) {
      param_ess[[param]] <- ess(param_chain)
    } else {
      param_ess[[param]] <- NA
      if (!ess_message_shown) {
        message(paste0("ESS cannot be computed with only one chain ",
                       "Run at least 2 chains."))
        ess_message_shown <- TRUE
      }
    }

    param_rhat[[param]] <- rhat(param_chain)
  }

  theta_chain_post <- lapply(theta_chain_post, as.data.frame)
  theta_chain_post <- bind_rows(theta_chain_post, .id = "chain")

  result <- list(
    theta_chain = theta_chain_post,
    diagnostics = list(ess = param_ess, rhat = param_rhat)
  )

  if (return_latent_state_est) {
    result$latent_state_chain <- state_est_chain_post
  }

  class(result) <- "pmmh_output"

  print(result)

  # If any ESS<400 print a warning
  if (any(sapply(param_ess, function(x) !is.na(x) && x < 400))) {
    warning(paste0(
      "Some ESS values are below 400, indicating poor mixing. ",
      "Consider running the chains for more iterations."
    ))
  }

  # If any Rhat>1.01 print a warning
  if (any(sapply(param_rhat, function(x) x > 1.01 && !is.na(x)))) {
    warning(paste0(
      "\nSome Rhat values are above 1.01, indicating that the chains ",
      "have not converged. \nConsider running the chains for more iterations ",
      "and/or increase burn_in."
    ))
  }

  result
}
