#' Bootstrap Particle Filter (BPF)
#'
#' Implements a bootstrap particle filter for sequential Bayesian inference in
#' state space models using sequential Monte Carlo methods.
#'
#'
#' @inheritParams particle_filter_common_params
#'
#' @inherit particle_filter_returns return
#'
#' @inheritSection particle_filter_model_specification Model Specification
#'
#' @section The Effective Sample Size (ESS) is defined as:
#' \deqn{ESS = \left(\sum_{i=1}^{n} w_i^2\right)^{-1},}
#' where \eqn{w_i} are the normalized weights of the particles.
#'
#' Default resampling method is stratified resampling, which has lower variance
#' than multinomial resampling (Douc et al., 2005).
#'
#' @references
#' Gordon, N. J., Salmond, D. J., & Smith, A. F. M. (1993). Novel approach to
#' nonlinear/non-Gaussian Bayesian state estimation. IEE Proceedings F (Radar
#' and Signal Processing), 140(2), 107–113.
#' \doi{doi:10.1049/ip-f-2.1993.0015}
#'
#' Douc, R., Cappé, O., & Moulines, E. (2005). Comparison of
#' Resampling Schemes for Particle Filtering.
#' Accessible at: \url{https://arxiv.org/abs/cs/0507025}
#'
#' @export
#'
#' @examples
#' init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
#' transition_fn <- function(particles) particles + rnorm(length(particles))
#' log_likelihood_fn <- function(y, particles) {
#'   dnorm(y, mean = particles, sd = 1, log = TRUE)
#' }
#'
#' y <- cumsum(rnorm(50)) # dummy data
#' num_particles <- 100
#'
#' # Run the particle filter using default settings.
#' result <- bootstrap_filter(
#'   y = y,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn
#' )
#' plot(result$state_est,
#'   type = "l", col = "blue", main = "State Estimates",
#'   ylim = range(c(result$state_est, y))
#' )
#' points(y, col = "red", pch = 20)
#'
#' # With parameters
#' init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
#' transition_fn <- function(particles, mu) {
#'   particles + rnorm(length(particles), mean = mu)
#' }
#' log_likelihood_fn <- function(y, particles, sigma) {
#'   dnorm(y, mean = particles, sd = sigma, log = TRUE)
#' }
#'
#' y <- cumsum(rnorm(50)) # dummy data
#' num_particles <- 100
#'
#' # Run the bootstrap particle filter using default settings.
#' result <- bootstrap_filter(
#'   y = y,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn,
#'   mu = 1,
#'   sigma = 1
#' )
#' plot(result$state_est,
#'   type = "l", col = "blue", main = "State Estimates",
#'   ylim = range(c(result$state_est, y))
#' )
#' points(y, col = "red", pch = 20)
#'
#' # With observations gaps
#' init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
#' transition_fn <- function(particles, mu) {
#'   particles + rnorm(length(particles), mean = mu)
#' }
#' log_likelihood_fn <- function(y, particles, sigma) {
#'   dnorm(y, mean = particles, sd = sigma, log = TRUE)
#' }
#'
#' # Generate data using DGP
#' simulate_ssm <- function(num_steps, mu, sigma) {
#'   x <- numeric(num_steps)
#'   y <- numeric(num_steps)
#'   x[1] <- rnorm(1, mean = 0, sd = sigma)
#'   y[1] <- rnorm(1, mean = x[1], sd = sigma)
#'   for (t in 2:num_steps) {
#'     x[t] <- mu * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma)
#'     y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma)
#'   }
#'   y
#' }
#'
#' data <- simulate_ssm(10, mu = 1, sigma = 1)
#' # Suppose we have data for t=1,2,3,5,6,7,8,9,10 (i.e., missing at t=4)
#'
#' obs_times <- c(1, 2, 3, 5, 6, 7, 8, 9, 10)
#' data_obs <- data[obs_times]
#'
#' num_particles <- 100
#' # Specify observation times in the bootstrap particle filter using obs_times
#' result <- bootstrap_filter(
#'   y = data_obs,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn,
#'   obs_times = obs_times,
#'   mu = 1,
#'   sigma = 1,
#' )
#' plot(result$state_est,
#'   type = "l", col = "blue", main = "State Estimates",
#'   ylim = range(c(result$state_est, data))
#' )
#' points(data_obs, col = "red", pch = 20)
bootstrap_filter <- function(
    y, num_particles, init_fn, transition_fn,
    log_likelihood_fn, obs_times = NULL,
    resample_algorithm = c("SISAR", "SISR", "SIS"),
    resample_fn = c("stratified", "systematic", "multinomial"),
    threshold = NULL,
    return_particles = TRUE,
    ...) {
  resample_algorithm <- match.arg(resample_algorithm)
  resample_fn <- match.arg(resample_fn)

  # Ensure log_likelihood_fn can accept additional parameters
  log_likelihood_fn <- .ensure_dots(log_likelihood_fn)
  # Add t arg if in log_likelihood_fn
  if (!"t" %in% names(formals(log_likelihood_fn))) {
    formals(log_likelihood_fn) <- c(
      formals(log_likelihood_fn), alist(t = NULL)
    )
  }

  # Weight function just calls log_likelihood_fn
  weight_fn <- function(y, particles, t, ...) {
    log_likelihood_fn(y = y, particles = particles, t = t, ...)
  }

  # Call core filter engine
  .particle_filter_core(
    y = y,
    num_particles = num_particles,
    init_fn = init_fn,
    transition_fn = transition_fn,
    weight_fn = weight_fn,
    aux_weight_fn = NULL,
    obs_times = obs_times,
    algorithm = "BPF",
    resample_algorithm = resample_algorithm,
    resample_fn = resample_fn,
    threshold = threshold,
    return_particles = return_particles,
    ...
  )
}
