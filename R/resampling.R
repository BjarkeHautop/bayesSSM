#' @title Internal Resampling Functions
#' @description Helper functions for resampling particles in a particle filter.
#' These functions implement multinomial, stratified, and systematic resampling.
#'
#' @importFrom stats runif
#'
#' @keywords internal

# Multinomial resampling: samples indices with replacement based on weights.
.resample_multinomial <- function(particles, weights) {
  # Check all are non-negative
  if (any(weights < 0)) {
    stop("Weights must be non-negative")
  }

  total_weight <- sum(weights)
  if (total_weight == 0) {
    stop("Sum of weights must be greater than 0")
  }

  if (is.matrix(particles)) {
    n <- nrow(particles)
    # Throw error if particles dim doesn't match weights length
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    indices <- sample(1:n, size = n, replace = TRUE, prob = weights)
    particles[indices, , drop = FALSE]
  } else {
    n <- length(particles)
    # Throw error if particles dim doesn't match weights length
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    indices <- sample(1:n, size = n, replace = TRUE, prob = weights)
    particles[indices]
  }
}

# Stratified resampling: divides [0,1] into N strata and samples one point per
# stratum.
.resample_stratified <- function(particles, weights) {
  # Check all are non-negative
  if (any(weights < 0)) {
    stop("Weights must be non-negative")
  }

  total_weight <- sum(weights)
  if (total_weight == 0) {
    stop("Sum of weights must be greater than 0")
  }

  if (is.matrix(particles)) {
    n <- nrow(particles)
    # Throw error if particles dim doesn't match weights length
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    positions <- (runif(n) + seq_len(n) - 1) / n
    cumulative_sum <- cumsum(weights)
    indices <- numeric(n)
    i <- 1
    j <- 1
    while (i <= n) {
      if (positions[i] < cumulative_sum[j]) {
        indices[i] <- j
        i <- i + 1
      } else {
        j <- j + 1
      }
    }
    particles[indices, , drop = FALSE]
  } else {
    n <- length(particles)
    # Throw error if particles dim doesn't match weights length
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    positions <- (runif(n) + seq_len(n) - 1) / n
    cumulative_sum <- cumsum(weights)
    indices <- numeric(n)
    i <- 1
    j <- 1
    while (i <= n) {
      if (positions[i] < cumulative_sum[j]) {
        indices[i] <- j
        i <- i + 1
      } else {
        j <- j + 1
      }
    }
    particles[indices]
  }
}

# Systematic resampling: similar to stratified sampling but with a single
# random start
.resample_systematic <- function(particles, weights) {
  # Check all are non-negative
  if (any(weights < 0)) {
    stop("Weights must be non-negative")
  }

  total_weight <- sum(weights)
  if (total_weight == 0) {
    stop("Sum of weights must be greater than 0")
  }

  if (is.matrix(particles)) {
    n <- nrow(particles)
    # Throw error if particles dim doesn't match weights length
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    positions <- (runif(1) + seq_len(n) - 1) / n
    cumulative_sum <- cumsum(weights)
    indices <- numeric(n)
    i <- 1
    j <- 1
    while (i <= n) {
      if (positions[i] < cumulative_sum[j]) {
        indices[i] <- j
        i <- i + 1
      } else {
        j <- j + 1
      }
    }
    particles[indices, , drop = FALSE]
  } else {
    n <- length(particles)
    # Throw error if particles dim doesn't match weights length
    if (n != length(weights)) {
      stop("Number of particles must match the length of weights")
    }
    positions <- (runif(1) + seq_len(n) - 1) / n
    cumulative_sum <- cumsum(weights)
    indices <- numeric(n)
    i <- 1
    j <- 1
    while (i <= n) {
      if (positions[i] < cumulative_sum[j]) {
        indices[i] <- j
        i <- i + 1
      } else {
        j <- j + 1
      }
    }
    particles[indices]
  }
}
