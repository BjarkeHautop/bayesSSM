#' @title Internal Resampling Functions
#' @description Helper functions for resampling particles in a particle filter.
#' These functions implement multinomial, stratified, and systematic resampling.
#' @keywords internal

# Multinomial resampling: samples indices with replacement based on weights
.resample_multinomial <- function(particles, weights) {
  n <- length(weights)
  indices <- sample(seq_len(n), size = n, replace = TRUE, prob = weights)

  particles[indices]
}

# Stratified resampling: divides [0,1] into N strata and samples one point per
# stratum.
.resample_stratified <- function(particles, weights) {
  n <- length(weights)
  positions <- (stats::runif(1) + seq_len(n) - 1) / n
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

# Systematic resampling: similar to stratified sampling but with a single
# random start
.resample_systematic <- function(particles, weights) {
  n <- length(weights)
  u0 <- stats::runif(1, 0, 1 / n)
  positions <- u0 + (seq_len(n) - 1) / n
  cumulative_sum <- cumsum(weights)
  indices <- numeric(n)
  j <- 1
  for (i in seq_len(n)) {
    while (positions[i] > cumulative_sum[j]) {
      j <- j + 1
    }
    indices[i] <- j
  }

  particles[indices]
}
