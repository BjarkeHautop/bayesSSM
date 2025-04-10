#' Compute split Rhat statistic
#'
#' @param chains A matrix of dimensions m (iterations) x k (chains).
#'
#' @returns The split-Rhat statistic.
#'
#' @details Uses the formula for split-Rhat proposed by Gelman et al. (2013).
#'
#' @references Gelman et al. (2013). Bayesian Data Analysis, 3rd Edition.
#'
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#' chains <- matrix(rnorm(3000), nrow = 1000, ncol = 3)
#' rhat(chains)
rhat <- function(chains) {
  # Check that the input is a matrix
  if (!is.matrix(chains)) {
    stop("Input 'chains' must be a matrix with dimensions m (iterations) x k
         (chains).")
  }

  m <- nrow(chains)
  k <- ncol(chains)

  # Ensure even number of iterations
  if (m %% 2 == 1) {
    chains <- chains[-m, ]  # Drop last iteration
    m <- nrow(chains)       # Update m
  }

  # Split each chain into two parts
  chains_split <- matrix(NA, nrow = m %/% 2, ncol = 2 * k)
  for (i in 1:k) {
    # Split the k-th chain into two parts
    chains_split[, 2 * i - 1] <- chains[1:(m %/% 2), k]
    chains_split[, 2 * i] <- chains[(m %/% 2 + 1):m, k]
  }

  # Compute the mean of each new chain
  chain_means <- colMeans(chains_split)
  overall_mean <- mean(chain_means)

  # Compute the between-chain variance, B
  b <- m / (2 * k - 1) * sum((chain_means - overall_mean)^2)

  # Compute the within-chain variances
  chain_vars <- apply(chains_split, 2, var)

  # If any chain_vars is zero give warning and return NA
  if (any(chain_vars == 0)) {
    warning("One or more chains have zero variance.")
    return(NA)
  }

  w <- mean(chain_vars)

  # Compute the marginal posterior variance estimator
  var_hat <- ((m - 1) / m) * w + (1 / m) * b

  # Compute the potential scale reduction statistic, Rhat
  r_hat <- sqrt(var_hat / w)

  if (r_hat >= 0.99 && r_hat <= 1) {
    # Numerical issue can make it slightly less than 1, change it to 1 if so
    r_hat <- 1
  }
  r_hat
}
