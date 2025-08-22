#' Common Parameters for Particle Filters
#'
#' These parameters are shared by particle filter implementations such as the
#' bootstrap filter, auxiliary particle filter, and resample-move particle
#' filter.
#'
#' @param y A numeric vector or matrix of observations. Each row represents an
#' observation at a time step. If observations are not equally spaced, use the
#' \code{obs_times} argument.
#' @param num_particles A positive integer specifying the number of particles.
#' @param init_fn A function to initialize the particles. Should take
#' `num_particles` and return a matrix or vector of initial states. Additional
#' model parameters can be passed via \code{...}.
#' @param transition_fn A function for propagating particles. Should take
#' `particles` and optionally `t`. Additional model parameters via \code{...}.
#' @param log_likelihood_fn A function that returns the log-likelihood for each
#' particle given the current observation, particles, and optionally `t`.
#' Additional parameters via \code{...}.
#' @param obs_times A numeric vector specifying observation time points. Must
#' match the number of rows in \code{y}, or defaults to \code{1:nrow(y)}.
#' @param resample_algorithm A character string specifying the filtering
#' resample algorithm:
#' \code{"SIS"} for no resampling, \code{"SISR"} for resampling at every time
#' step, or \code{"SISAR"} for adaptive resampling when ESS
#' drops below \code{threshold}. Using \code{"SISR"} or \code{"SISAR"} to
#' avoid weight degeneracy is recommedended. Default is \code{"SISAR"}.
#' @param resample_fn A string indicating the resampling method:
#' \code{"stratified"}, \code{"systematic"}, or \code{"multinomial"}.
#' Default is \code{"stratified"}.
#' @param threshold A numeric value specifying the ESS threshold for
#' \code{"SISAR"}. Defaults to \code{num_particles / 2}.
#' @param return_particles Logical; if \code{TRUE}, returns the full particle
#' and weight histories.
#' @param ... Additional arguments passed to \code{init_fn},
#' \code{transition_fn}, and \code{log_likelihood_fn}.
#'
#' @name particle_filter_common_params
#' @keywords internal
NULL

#' Model Specification for Particle Filters
#'
#' @section Model Specification:
#' Particle filter implementations in this package assume a discrete-time
#' state-space model defined by:
#'
#' \itemize{
#'   \item A sequence of latent states \eqn{x_0, x_1, \ldots, x_T} evolving
#'   according to a Markov process.
#'   \item Observations \eqn{y_1, \ldots, y_T} that are conditionally independent
#'   given the corresponding latent states.
#' }
#'
#' The model is specified as:
#' \deqn{x_0 \sim \mu_\theta}
#' \deqn{x_t \sim f_\theta(x_t \mid x_{t-1}), \quad t = 1, \ldots, T}
#' \deqn{y_t \sim g_\theta(y_t \mid x_t), \quad t = 1, \ldots, T}
#'
#' where \eqn{\theta} denotes model parameters passed via \code{...}.
#'
#' The user provides the following functions:
#' \itemize{
#'   \item \code{init_fn}: draws from the initial distribution
#'   \eqn{\mu_\theta}.
#'   \item \code{transition_fn}: generates or evaluates the transition
#'   density \eqn{f_\theta}.
#'   \item \code{weight_fn}: evaluates the observation likelihood
#'   \eqn{g_\theta}.
#' }
#'
#' @name particle_filter_model_specification
#' @keywords internal
NULL




#' Shared Return Values for Particle Filters
#'
#' This block documents the common return value for particle filtering
#' functions.
#'
#' @return A list with components:
#' \describe{
#'   \item{state_est}{Estimated states over time (weighted mean of particles).}
#'   \item{ess}{Effective sample size at each time step.}
#'   \item{loglike}{Total log-likelihood.}
#'   \item{loglike_history}{Log-likelihood at each time step.}
#'   \item{algorithm}{The filtering algorithm used.}
#'   \item{particles_history}{Matrix of particle states over time
#'   (if \code{return_particles = TRUE}).}
#'   \item{weights_history}{Matrix of particle weights over time
#'   (if \code{return_particles = TRUE}).}
#' }
#'
#' @name particle_filter_returns
#' @keywords internal
NULL
