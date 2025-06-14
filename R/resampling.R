#' @title Internal Resampling Functions
#' @description Helper functions for resampling particles in a particle filter.
#' These functions implement multinomial, stratified, and systematic resampling.
#'
#' @importFrom stats runif
#' @importFrom Rcpp sourceCpp evalCpp
#' @useDynLib bayesSSM, .registration = TRUE
#'
#' @keywords internal
#' @noRd

# Multinomial resampling: samples indices with replacement based on weights.
.resample_multinomial <- function(particles, weights) {
  if (is.matrix(particles)) {
    n <- nrow(particles)
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    indices <- resample_multinomial_cpp(n, weights)
    particles[indices, , drop = FALSE]
  } else {
    n <- length(particles)
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    indices <- resample_multinomial_cpp(n, weights)
    particles[indices]
  }
}

# Stratified resampling: divides [0,1] into N strata and
# samples one point per stratum.
.resample_stratified <- function(particles, weights) {
  if (is.matrix(particles)) {
    n <- nrow(particles)
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    indices <- resample_stratified_cpp(n, weights)
    particles[indices, , drop = FALSE]
  } else {
    n <- length(particles)
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    indices <- resample_stratified_cpp(n, weights)
    particles[indices]
  }
}

# Systematic resampling: similar to stratified sampling but with a single
# random start
.resample_systematic <- function(particles, weights) {
  if (is.matrix(particles)) {
    n <- nrow(particles)
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    indices <- resample_systematic_cpp(n, weights)
    particles[indices, , drop = FALSE]
  } else {
    n <- length(particles)
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    indices <- resample_systematic_cpp(n, weights)
    particles[indices]
  }
}
