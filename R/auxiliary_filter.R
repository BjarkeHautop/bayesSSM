#' Auxiliary Particle Filter (APF)
#'
#' The Auxiliary Particle Filter enhances the bootstrap filter by using a
#' look ahead step: particles are weighted based on the predicted likelihood of
#' the next observation before resampling. This anticipatory weighting reduces
#' particle degeneracy and improves filter efficiency.
#'
#' @inheritParams particle_filter_common_params
#' @param aux_log_likelihood_fn A function that computes the log-likelihood of
#' the next observation given the current particles. It should accept
#' arguments `y`, `particles`, optionally `t`, and any additional model-specific
#' parameters via \code{...}. It returns a numeric vector of log-likelihoods.
#'
#' @inherit particle_filter_returns return
#'
#' @inheritSection particle_filter_model_specification Model Specification
#'
#' @section The Auxiliary Particle Filter (APF):
#'
#' The Auxiliary Particle Filter (APF) was introduced by Pitt and Shephard
#' (1999) to improve upon the standard bootstrap filter by incorporating a
#' look ahead step. Before resampling at time \eqn{t}, particles are weighted by
#' an auxiliary weight proportional to an estimate of the likelihood of the next
#' observation, guiding resampling to favor particles likely to contribute to
#' future predictions.
#'
#' Specifically, if \eqn{w_{t-1}^i} are the normalized weights and
#' \eqn{x_{t-1}^i} are the particles at time \eqn{t-1}, then auxiliary weights
#' are computed as
#' \deqn{
#'   \tilde{w}_t^i \propto w_{t-1}^i \, p(y_t | \mu_t^i),
#' }
#' where \eqn{\mu_t^i} is a predictive summary (e.g., the expected next state)
#' of the particle \eqn{x_{t-1}^i}. Resampling is performed using
#' \eqn{\tilde{w}_t^i} instead of \eqn{w_{t-1}^i}.
#' This can reduce the variance of the importance weights at time \eqn{t} and
#' help mitigate particle degeneracy, especially if the auxiliary weights are
#' chosen well.
#'
#' Default resampling method is stratified resampling, which has lower variance
#' than multinomial resampling (Douc et al., 2005).
#'
#' @references
#' Pitt, M. K., & Shephard, N. (1999). Filtering via simulation: Auxiliary
#' particle filters. Journal of the American Statistical Association, 94(446),
#' 590–599. \doi{doi:10.1080/01621459.1999.10474153}
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
#' aux_log_likelihood_fn <- function(y, particles) {
#'   # Predict next state (mean stays same) and compute log p(y | x)
#'   mean_forecast <- particles # since E[x'] = x in this model
#'   dnorm(y, mean = mean_forecast, sd = 1, log = TRUE)
#' }
#'
#' y <- cumsum(rnorm(50)) # dummy data
#' num_particles <- 100
#'
#' result <- auxiliary_filter(
#'   y = y,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn,
#'   aux_log_likelihood_fn = aux_log_likelihood_fn
#' )
#' plot(result$state_est,
#'   type = "l", col = "blue", main = "APF: State Estimates",
#'   ylim = range(c(result$state_est, y))
#' )
#' points(y, col = "red", pch = 20)
#'
#' # ---- With parameters ----
#' init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
#' transition_fn <- function(particles, mu) {
#'   particles + rnorm(length(particles), mean = mu)
#' }
#' log_likelihood_fn <- function(y, particles, sigma) {
#'   dnorm(y, mean = particles, sd = sigma, log = TRUE)
#' }
#' aux_log_likelihood_fn <- function(y, particles, mu, sigma) {
#'   # Forecast mean of x' given x, then evaluate p(y | forecast)
#'   forecast <- particles + mu
#'   dnorm(y, mean = forecast, sd = sigma, log = TRUE)
#' }
#'
#' y <- cumsum(rnorm(50)) # dummy data
#' num_particles <- 100
#'
#' result <- auxiliary_filter(
#'   y = y,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn,
#'   aux_log_likelihood_fn = aux_log_likelihood_fn,
#'   mu = 1,
#'   sigma = 1
#' )
#' plot(result$state_est,
#'   type = "l", col = "blue", main = "APF with Parameters",
#'   ylim = range(c(result$state_est, y))
#' )
#' points(y, col = "red", pch = 20)
#'
#' # ---- With observation gaps ----
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
#' obs_times <- c(1, 2, 3, 5, 6, 7, 8, 9, 10) # Missing at t = 4
#' data_obs <- data[obs_times]
#'
#' init_fn <- function(num_particles) rnorm(num_particles, 0, 1)
#' transition_fn <- function(particles, mu) {
#'   particles + rnorm(length(particles), mean = mu)
#' }
#' log_likelihood_fn <- function(y, particles, sigma) {
#'   dnorm(y, mean = particles, sd = sigma, log = TRUE)
#' }
#' aux_log_likelihood_fn <- function(y, particles, mu, sigma) {
#'   forecast <- particles + mu
#'   dnorm(y, mean = forecast, sd = sigma, log = TRUE)
#' }
#'
#' num_particles <- 100
#' result <- auxiliary_filter(
#'   y = data_obs,
#'   num_particles = num_particles,
#'   init_fn = init_fn,
#'   transition_fn = transition_fn,
#'   log_likelihood_fn = log_likelihood_fn,
#'   aux_log_likelihood_fn = aux_log_likelihood_fn,
#'   obs_times = obs_times,
#'   mu = 1,
#'   sigma = 1
#' )
#' plot(result$state_est,
#'   type = "l", col = "blue", main = "APF with Observation Gaps",
#'   ylim = range(c(result$state_est, data))
#' )
#' points(data_obs, col = "red", pch = 20)
auxiliary_filter <- function(
    y, num_particles, init_fn, transition_fn,
    log_likelihood_fn, aux_log_likelihood_fn,
    obs_times = NULL,
    resample_algorithm = c("SISAR", "SISR", "SIS"),
    resample_fn = c("stratified", "systematic", "multinomial"),
    threshold = NULL,
    return_particles = TRUE,
    ...) {
  resample_fn <- match.arg(resample_fn)
  resample_algorithm <- match.arg(resample_algorithm)

  # Ensure log-likelihood functions accept `t` and `...`
  log_likelihood_fn <- .ensure_dots(log_likelihood_fn)
  aux_log_likelihood_fn <- .ensure_dots(aux_log_likelihood_fn)

  if (!"t" %in% names(formals(log_likelihood_fn))) {
    formals(log_likelihood_fn) <- c(formals(log_likelihood_fn), alist(t = NULL))
  }
  if (!"t" %in% names(formals(aux_log_likelihood_fn))) {
    formals(aux_log_likelihood_fn) <- c(
      formals(aux_log_likelihood_fn),
      alist(t = NULL)
    )
  }

  # Define standard weight function (log-likelihood)
  weight_fn <- function(y, particles, t, ...) {
    log_likelihood_fn(y = y, particles = particles, t = t, ...)
  }

  # Define auxiliary weight function
  aux_weight_fn <- function(y, particles, t, ...) {
    aux_log_likelihood_fn(y = y, particles = particles, t = t, ...)
  }

  # Call the core filter with algorithm set to APF
  .particle_filter_core(
    y = y,
    num_particles = num_particles,
    init_fn = init_fn,
    transition_fn = transition_fn,
    weight_fn = weight_fn,
    aux_weight_fn = aux_weight_fn,
    obs_times = obs_times,
    algorithm = "APF",
    resample_fn = resample_fn,
    resample_algorithm = resample_algorithm,
    return_particles = return_particles,
    ...
  )
}
