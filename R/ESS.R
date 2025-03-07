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
#' @export
#'
#' @examples
#' chains <- matrix(rnorm(3000), nrow = 1000, ncol = 3)
#' ess(chains)
ess <- function(chains) {
  if (!is.matrix(chains)) {
    stop("Input 'chains' must be a matrix with dimensions M (iterations) x K
         (chains).")
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
  chain_vars <- apply(chains, 2, stats::var)

  # If any chain_vars is zero give error
  if (any(chain_vars == 0)) {
    warning("One or more chains have zero variance.")
    return(NA)
  }

  w <- mean(chain_vars)

  # Marginal posterior variance estimator
  var_hat <- ((m - 1) / m) * w + (1 / m) * b

  # --- Compute autocorrelation using the proposed formula ---
  # First, we calculate the autocorrelation function for each chain.
  # Note: stats::acf returns autocorrelation values.
  acf_matrix <- matrix(NA, nrow = m, ncol = k)
  for (i in 1:k) {
    acf_obj <- stats::acf(chains[, i], lag.max = m - 1, plot = FALSE)
    # acf_obj$acf is an array of dimensions (M, 1, 1); extract as a vector.
    acf_matrix[, i] <- acf_obj$acf[, 1, 1]
  }

  # Now compute the average autocovariance for each lag.
  # For chain k, the autocovariance at lag t is given by:
  #   gamma_t,k = acf_matrix[t+1, k] * chain_vars[k]
  # (so that for t=0, gamma_0,k equals chain_vars[k])
  avg_autocov <- numeric(m)
  for (t in 0:(m - 1)) {
    gamma_t_k <- acf_matrix[t + 1, ] * chain_vars
    avg_autocov[t + 1] <- mean(gamma_t_k)
  }

  hat_rho <- 1 - (w - avg_autocov) / var_hat

  # --- Apply Geyer's initial monotone sequence method ---
  # We form pairs P_t = hat_rho_{2t-1} + hat_rho_{2t} (with t starting at 1)
  # and accumulate until the first non-positive pair.
  sum_rho <- 0
  t_pair <- 1
  while (TRUE) {
    idx1 <- 2 * t_pair # corresponds to lag (2*t_pair - 1)
    idx2 <- 2 * t_pair + 1 # corresponds to lag (2*t_pair)

    if (idx1 > length(hat_rho)) break # no further lags available

    pair_sum <- if (idx2 > length(hat_rho)) {
      hat_rho[idx1]
    } else {
      hat_rho[idx1] + hat_rho[idx2]
    }

    if (pair_sum < 0) break # stop if the sum of the pair is negative

    sum_rho <- sum_rho + pair_sum
    t_pair <- t_pair + 1
  }

  tau <- 1 + 2 * sum_rho

  # Effective sample size
  ess <- (k * m) / tau

  ess
}
