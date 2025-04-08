#' Estimate effective sample size (ESS) of MCMC chains.
#'
#' @param chains A matrix of dimensions m (iterations) x k (chains).
#'
#' @returns The estimated effective sample size (ess) of the chains.
#'
#' @details Uses the formula for ESS proposed by Vehtari et al. (2021).
#'
#' @references Vehtari et al. (2021). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of MCMC.
#' Available at: https://doi.org/10.1214/20-BA1221
#'
#' @importFrom stats var acf
#' @export
#'
#' @examples
#' chains <- matrix(rnorm(3000), nrow = 1000, ncol = 3)
#' ess(chains)
ess <- function(chains) {
  if (!is.matrix(chains)) {
    stop(
      paste0(
        "Input 'chains' must be a matrix with dimensions M (iterations) ",
        "x K (chains)."
      )
    )
  }

  m <- nrow(chains)
  k <- ncol(chains)
  if (m < 2) {
    stop("Number of iterations must be at least 2.")
  }
  if (k < 2) {
    stop("Number of chains must be at least 2.")
  }

  # --- Compute within-chain and between-chain variances ---
  chain_means <- colMeans(chains)
  overall_mean <- mean(chain_means)

  # Between-chain variance, b
  b <- m / (k - 1) * sum((chain_means - overall_mean)^2)

  # Within-chain variances, W
  chain_vars <- apply(chains, 2, var)

  # If any chain_vars is zero give error
  if (any(chain_vars == 0)) {
    warning("One or more chains have zero variance.")
    return(NA)
  }

  w <- mean(chain_vars)

  # Marginal posterior variance estimator
  var_hat <- ((m - 1) / m) * w + (1 / m) * b

  # --- Compute autocorrelations ---
  # Create matrix to store acf values (lags 0 to m-1)
  acf_matrix <- matrix(NA, nrow = m, ncol = k)
  for (i in 1:k) {
    acf_obj <- acf(chains[, i], lag.max = m - 1, plot = FALSE)
    # acf_obj$acf is an array with dimensions (m, 1, 1); extract as vector.
    acf_matrix[, i] <- acf_obj$acf[, 1, 1]
  }

  hat_rho <- numeric(m)
  for (t in 0:(m - 1)) {
    term <- (1 / k) * sum(chain_vars * acf_matrix[t + 1, ])
    hat_rho[t + 1] <- 1 - (w - term) / var_hat
  }

  # --- Apply Geyer's initial monotone sequence method ---
  # Form pairs: for t = 1, 2, ..., define P_t = hat_rho[2*t] + hat_rho[2*t + 1]
  max_pairs <- floor((length(hat_rho) - 1) / 2)
  pairs <- numeric(max_pairs)
  for (t in 1:max_pairs) {
    idx1 <- 2 * t     # corresponds to lag (2*t - 1)
    idx2 <- 2 * t + 1 # corresponds to lag (2*t)
    # If the second index exceeds available lags, use only the first
    if (idx2 > length(hat_rho)) {
      pairs[t] <- hat_rho[idx1]
    } else {
      pairs[t] <- hat_rho[idx1] + hat_rho[idx2]
    }
  }

  # Enforce monotonicity on the pairs:
  # Each pair should not exceed the smallest preceding pair.
  if (length(pairs) >= 2) {
    for (t in 2:length(pairs)) {
      if (pairs[t] > pairs[t - 1]) {
        pairs[t] <- pairs[t - 1]
      }
    }
  }

  # Sum the pairs until the first negative pair is encountered.
  sum_rho <- 0
  for (t in seq_along(pairs)) {
    if (pairs[t] < 0) break
    sum_rho <- sum_rho + pairs[t]
  }

  tau <- 1 + 2 * sum_rho

  # Effective sample size
  ess <- (k * m) / tau

  ess
}
