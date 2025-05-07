#' Helper function to validate input of user-defined functions and priors
#'
#' @param init_fn A function to initialize the state-space model.
#' @param transition_fn A function that defines the state transition of the
#' state-space model.
#' @param log_likelihood_fn A function that calculates the log-likelihood
#' for the state-space model given latent states.
#' @param log_priors A list of functions for computing the log-prior of each
#' parameter.
#' @param pilot_init_params A vector of initial parameter values.
#'
#' @returns NULL
#'
#' @keywords internal
.check_params_match <- function(
    init_fn, transition_fn, log_likelihood_fn, pilot_init_params,
    log_priors) {
  get_fn_params <- function(fn) {
    names(formals(fn))
  }
  # Check if 'particles' or 'num_particles' is in init_fn. 'particles' will
  # be deprecated in the future.
  check_init_fn <- function(fn, fn_name) {
    fn_args <- get_fn_params(fn)
    if (!("particles" %in% fn_args || "num_particles" %in% fn_args)) {
      stop(
        paste(
          fn_name,
          "does not contain 'particles' or 'num_particles' as an argument"
        )
      )
    }
  }

  check_init_fn(init_fn, "init_fn")

  # Check if 'particles' is in function
  check_particles <- function(fn, fn_name) {
    fn_args <- get_fn_params(fn)
    if (!"particles" %in% fn_args) {
      stop(paste(fn_name, "does not contain 'particles' as an argument"))
    }
  }

  check_particles(transition_fn, "transition_fn")
  check_particles(log_likelihood_fn, "log_likelihood_fn")

  # Check if 'y' is in log_likelihood_fn
  if (!"y" %in% get_fn_params(log_likelihood_fn)) {
    stop("log_likelihood_fn does not contain 'y' as an argument")
  }


  # Combine parameters from all three functions
  fn_params <- unique(c(
    get_fn_params(init_fn),
    get_fn_params(transition_fn),
    get_fn_params(log_likelihood_fn)
  ))
  # Drop 'num_particles', 'particles', 'y', 't', and '...' from the parameters
  drop_names <- c("num_particles", "particles", "y", "t", "...")
  fn_params <- fn_params[!(fn_params %in% drop_names)]

  # Check if the parameters match init_params
  if (!all(fn_params %in% names(pilot_init_params))) {
    stop("Parameters in functions do not match the names in pilot_init_params")
  }

  # Check if the parameters match log_priors
  if (!all(fn_params %in% names(log_priors))) {
    stop("Parameters in functions do not match the names in log_priors")
  }
}

# ---------------------------
# Helper functions for parameter transformation
# ---------------------------

#' Internal function to transform parameters
#'
#' @param theta parameter vector
#' @param transform transformation type for each parameter
#'
#' @returns transformed parameter vector
#'
#' @keywords internal
.transform_params <- function(theta, transform) {
  sapply(seq_along(theta), function(j) {
    if (transform[j] == "log") {
      log(theta[j])
    } else if (transform[j] == "invlogit") {
      1 / (1 + exp(-theta[j]))
    } else {
      theta[j]
    }
  })
}

#' Internal function to back-transform parameters
#'
#' @param theta_trans transformed parameter vector
#' @param transform transformation type for each parameter
#'
#' @returns original parameter vector
#'
#' @keywords internal
.back_transform_params <- function(theta_trans, transform) {
  sapply(seq_along(theta_trans), function(j) {
    if (transform[j] == "log") {
      exp(theta_trans[j])
    } else if (transform[j] == "invlogit") {
      log(theta_trans[j] / (1 - theta_trans[j])) # back-transforming from logit
    } else {
      theta_trans[j]
    }
  })
}

#' Internal function to compute the Jacobian of the transformation
#'
#' @param theta parameter vector (on original scale)
#' @param transform transformation type for each parameter
#'
#' @returns log-Jacobian of the transformation
#'
#' @keywords internal
.compute_log_jacobian <- function(theta, transform) {
  sum(sapply(seq_along(theta), function(j) {
    if (transform[j] == "log") {
      log(theta[j]) # log|dx/dz| = log(x)
    } else if (transform[j] == "invlogit") {
      val <- 1 / (1 + exp(-theta[j])) # inverse logit
      log(val * (1 - val)) # log|dx/dz| = log(x * (1 - x))
    } else {
      0 # no transformation
    }
  }))
}
