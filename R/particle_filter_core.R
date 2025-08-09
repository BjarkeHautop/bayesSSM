#' Core Particle Filter Function
#'
#' This function implements the underlying logic used for particle filters in
#' a state space model using sequential Monte Carlo methods.
#'
#' @inheritSection particle_filter_model_specification Model Specification
#'
#' @inheritParams particle_filter_common_params
#' @inherit particle_filter_returns return
#'
#' @param weight_fn A function that computes the log weights for the particles
#' given the observations and the current particles. It should take
#' `y`, `particles`, and `t` as arguments. The function can include any
#' model-specific parameters as named arguments.
#'
#' @importFrom stats dnorm rnorm
#' @importFrom checkmate assert_count assert_numeric assert_integerish
#' @keywords internal
.particle_filter_core <- function(
    y, num_particles, init_fn, transition_fn,
    weight_fn, # function(y, particles, t, ...) returning log weights
    aux_weight_fn = NULL, # for APF only
    move_fn = NULL, # for RMPF only
    obs_times = NULL,
    algorithm = c("BPF", "APF", "RMPF"),
    resample_algorithm = c("SIS", "SISR", "SISAR"),
    resample_fn = c("stratified", "systematic", "multinomial"),
    threshold = NULL,
    return_particles = TRUE,
    ...) {
  # Validate inputs
  assert_count(num_particles, positive = TRUE)

  algorithm <- match.arg(algorithm)
  resample_algorithm <- match.arg(resample_algorithm)
  resample_fn <- switch(match.arg(resample_fn),
    multinomial = .resample_multinomial,
    stratified = .resample_stratified,
    systematic = .resample_systematic
  )

  # Auto-set threshold
  if (is.null(threshold)) {
    threshold <- switch(resample_algorithm,
      SIS = Inf,
      SISR = num_particles,
      SISAR = num_particles / 2
    )
  }

  init_fn <- .ensure_dots(init_fn)
  transition_fn <- .ensure_dots(transition_fn)
  weight_fn <- .ensure_dots(weight_fn)
  if (!"t" %in% names(formals(transition_fn))) {
    formals(transition_fn) <- c(formals(transition_fn), alist(t = NULL))
  }
  if (!"t" %in% names(formals(weight_fn))) {
    formals(weight_fn) <- c(formals(weight_fn), alist(t = NULL))
  }

  if (!is.null(aux_weight_fn)) {
    aux_weight_fn <- .ensure_dots(aux_weight_fn)
    if (!"t" %in% names(formals(aux_weight_fn))) {
      formals(aux_weight_fn) <- c(formals(aux_weight_fn), alist(t = NULL))
    }
  }

  assert_numeric(y, any.missing = FALSE)
  if (is.vector(y)) y <- matrix(y, ncol = 1)
  if (is.null(obs_times)) obs_times <- seq_len(nrow(y))
  num_obs <- nrow(y)
  assert_integerish(obs_times, len = num_obs, lower = 1, sorted = TRUE)

  # Initialize particles
  particles <- init_fn(num_particles = num_particles, ...)
  if (is.null(dim(particles))) {
    if (length(particles) != num_particles) {
      stop("init_fn must return num_particles")
    }
    particles <- matrix(particles, nrow = num_particles)
  }
  if (nrow(particles) != num_particles) {
    stop("init_fn must return num_particles rows")
  }

  d <- ncol(particles)
  one_dim <- (d == 1)

  out_steps <- num_obs + 1
  state_est <- if (one_dim) {
    numeric(out_steps)
  } else {
    matrix(NA, nrow = out_steps, ncol = d)
  }
  ess_vec <- numeric(out_steps)
  loglike_history <- numeric(num_obs)
  loglike <- 0

  if (return_particles) {
    particles_history <- vector("list", out_steps)
    weights_history <- vector("list", out_steps)
  }

  # Initialization at t = 0
  weights <- rep(1 / num_particles, num_particles)
  ess_vec[1] <- 1 / sum(weights^2)
  if (one_dim) {
    state_est[1] <- sum(particles * weights)
  } else {
    state_est[1, ] <- colSums(particles * weights)
  }
  if (return_particles) {
    particles_history[[1]] <- particles
    weights_history[[1]] <- weights
  }

  if (algorithm == "SISAR" && is.null(threshold)) {
    threshold <- num_particles / 2
  }

  prev_t <- 0L
  for (i in seq_len(num_obs)) {
    gap <- obs_times[i] - prev_t
    for (step in seq_len(gap)) {
      tnow <- prev_t + step
      particles <- transition_fn(particles = particles, t = tnow, ...)
      if (is.null(dim(particles))) {
        if (length(particles) != num_particles) {
          stop("transition_fn must return num_particles")
        }
        particles <- matrix(particles, nrow = num_particles)
      } else if (nrow(particles) != num_particles) {
        stop("transition_fn must return num_particles rows")
      }
    }
    prev_t <- obs_times[i]

    # APF branch
    if (algorithm == "APF") {
      if (is.null(aux_weight_fn)) stop("APF requires aux_weight_fn")
      aux_log_weights <- aux_weight_fn(
        y = y[i, ],
        particles = particles,
        t = prev_t,
        ...
      )
      if (length(aux_log_weights) != num_particles) {
        stop("aux_weight_fn must return num_particles")
      }

      max_aux <- max(aux_log_weights)
      aux_weights <- exp(aux_log_weights - max_aux)
      aux_weights <- aux_weights / sum(aux_weights)
      ancestors <- resample_fn(seq_len(num_particles), aux_weights)
      old_particles <- particles
      particles <- old_particles[ancestors, , drop = FALSE]

      particles <- transition_fn(particles = particles, t = prev_t, ...)
      if (is.null(dim(particles))) {
        if (length(particles) != num_particles) {
          stop("transition_fn must return num_particles")
        }
        particles <- matrix(particles, nrow = num_particles)
      } else if (nrow(particles) != num_particles) {
        stop("transition_fn must return num_particles rows")
      }

      log_weights <- weight_fn(
        y = y[i, ],
        particles = particles,
        t = prev_t,
        ...
      )
      log_weights <- log_weights - aux_log_weights[ancestors]
    } else {
      log_weights <- weight_fn(
        y = y[i, ],
        particles = particles,
        t = prev_t,
        ...
      )
    }

    if (length(log_weights) != num_particles) {
      stop("weight_fn must return num_particles")
    }

    if (all(log_weights < -1e8)) {
      loglike <- -Inf
      loglike_history[i] <- -Inf
      result <- list(
        state_est = state_est, ess = ess_vec,
        loglike = loglike, loglike_history = loglike_history,
        algorithm = algorithm
      )
      if (return_particles) {
        result$particles_history <- particles_history
        result$weights_history <- weights_history
      }
      return(result)
    }

    max_logw <- max(log_weights)
    unnorm_weights <- exp(log_weights - max_logw)
    weight_sum <- sum(unnorm_weights)
    weights <- unnorm_weights / weight_sum
    loglike <- loglike + (max_logw + log(weight_sum) - log(num_particles))
    loglike_history[i] <- loglike

    ess <- 1 / sum(weights^2)
    ess_vec[i + 1] <- ess

    should_resample <- switch(resample_algorithm,
      SIS = FALSE,
      SISR = TRUE,
      SISAR = ess < threshold
    )

    if (algorithm == "RMPF" || should_resample) {
      particles <- resample_fn(particles, weights)
      weights <- rep(1 / num_particles, num_particles)
      ess_vec[i + 1] <- num_particles
    }

    if (algorithm == "RMPF") {
      if (is.null(move_fn)) stop("RMPF requires a move_fn")
      for (j in seq_len(num_particles)) {
        particles[j, ] <- move_fn(
          particle = particles[j, ],
          y = y[i, ], t = prev_t, ...
        )
      }
    }

    # Store output
    if (one_dim) {
      state_est[i + 1] <- sum(particles * weights)
    } else {
      state_est[i + 1, ] <- colSums(particles * weights)
    }
    if (return_particles) {
      particles_history[[i + 1]] <- particles
      weights_history[[i + 1]] <- weights
    }
  }

  result <- list(
    state_est = state_est,
    ess = ess_vec,
    loglike = loglike,
    loglike_history = loglike_history,
    algorithm = algorithm,
    resample_algorithm = resample_algorithm
  )

  if (return_particles) {
    result$particles_history <- do.call(
      rbind, lapply(particles_history, as.numeric)
    )
    result$weights_history <- do.call(
      rbind, lapply(weights_history, as.numeric)
    )
  }

  result
}
