.check_params_match <- function(
    init_fn_ssm, transition_fn_ssm, log_likelihood_fn_ssm, init_params,
    log_priors) {

  # Helper function to get parameter names excluding 'particles' and 'y'
  get_fn_params <- function(fn) {
    fn_args <- names(formals(fn))
    return(fn_args)
  }

  # Check if 'particles' is in init_fn_ssm, transition_fn_ssm, and
  # log_likelihood_fn_ssm
  check_particles <- function(fn, fn_name) {
    fn_args <- get_fn_params(fn)
    if (!"particles" %in% fn_args) {
      stop(paste(fn_name, "does not contain 'particles' as an argument"))
    }
  }

  # Check if 'y' is in log_likelihood_fn_ssm
  if (!"y" %in% get_fn_params(log_likelihood_fn_ssm)) {
    stop("log_likelihood_fn_ssm does not contain 'y' as an argument")
  }

  # Check if 'particles' is in all functions
  check_particles(init_fn_ssm, "init_fn_ssm")
  check_particles(transition_fn_ssm, "transition_fn_ssm")
  check_particles(log_likelihood_fn_ssm, "log_likelihood_fn_ssm")

  # Combine parameters from all three functions
  # (ignoring 'particles' and 'y' in the check)
  fn_params <- unique(c(get_fn_params(init_fn_ssm),
                        get_fn_params(transition_fn_ssm),
                        get_fn_params(log_likelihood_fn_ssm)))
  # Drop 'particles' and 'y'
  fn_params <- fn_params[!(fn_params %in% c("particles", "y"))]

  # Check if the parameters match init_params
  if (!all(fn_params %in% names(init_params))) {
    stop("Parameters in functions do not match the names in init_params")
  }

  # Check if the parameters match log_priors
  if (!all(fn_params %in% names(log_priors))) {
    stop("Parameters in functions do not match the names in log_priors")
  }
}
