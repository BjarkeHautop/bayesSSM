#' Estimate effective sample size (ESS) of MCMC chains.
#'
#' @param chains A matrix (iterations x chains) or a data.frame with a 'chain'
#' column and parameter columns.
#'
#' @returns The estimated effective sample size (ESS) if given a matrix, or a
#' named vector of ESS values if given a data frame.
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
#' # With a matrix:
#' chains <- matrix(rnorm(3000), nrow = 1000, ncol = 3)
#' ess(chains)
#'
#' # With a data frame:
#' chains_df <- data.frame(
#'   chain = rep(1:3, each = 1000),
#'   param1 = rnorm(3000),
#'   param2 = rnorm(3000)
#' )
#' ess(chains_df)
ess <- function(chains) {
  # Helper function to compute ESS from a matrix.
  compute_ess_matrix <- function(mat) {
    m <- nrow(mat)
    k <- ncol(mat)
    if (m < 2) {
      stop("Number of iterations must be at least 2.")
    }
    if (k < 2) {
      stop("Number of chains must be at least 2.")
    }

    # --- Compute within-chain and between-chain variances ---
    chain_means <- colMeans(mat)
    overall_mean <- mean(chain_means)

    # Between-chain variance, b
    b <- m / (k - 1) * sum((chain_means - overall_mean)^2)

    # Within-chain variances, W
    chain_vars <- apply(mat, 2, var)

    if (any(chain_vars == 0)) {
      warning("One or more chains have zero variance.")
      return(NA)
    }
    w <- mean(chain_vars)

    # Marginal posterior variance estimator
    var_hat <- ((m - 1) / m) * w + (1 / m) * b

    # --- Compute autocorrelations ---
    acf_matrix <- matrix(NA, nrow = m, ncol = k)
    for (i in 1:k) {
      acf_obj <- acf(mat[, i], lag.max = m - 1, plot = FALSE)
      acf_matrix[, i] <- acf_obj$acf[, 1, 1] # Extract as vector
    }

    hat_rho <- numeric(m)
    for (t in 0:(m - 1)) {
      term <- (1 / k) * sum(chain_vars * acf_matrix[t + 1, ])
      hat_rho[t + 1] <- 1 - (w - term) / var_hat
    }

    # --- Apply Geyer's initial monotone sequence method ---
    max_pairs <- floor((length(hat_rho) - 1) / 2)
    pairs <- numeric(max_pairs)
    for (t in 1:max_pairs) {
      idx1 <- 2 * t # corresponds to lag (2*t - 1)
      idx2 <- 2 * t + 1 # corresponds to lag (2*t)

      pairs[t] <- hat_rho[idx1] + hat_rho[idx2]
    }

    # Enforce monotonicity on the pairs:
    if (length(pairs) >= 2) {
      for (t in 2:length(pairs)) {
        if (pairs[t] > pairs[t - 1]) {
          pairs[t] <- pairs[t - 1]
        }
      }
    }

    # Sum pairs until the first negative pair
    sum_rho <- 0
    for (t in seq_along(pairs)) {
      if (pairs[t] < 0) break
      sum_rho <- sum_rho + pairs[t]
    }

    tau <- 1 + 2 * sum_rho
    ess <- (k * m) / tau

    ess
  }
  # Check if input is a matrix or a data frame
  if (!is.matrix(chains) && !is.data.frame(chains)) {
    stop("Input must be a matrix or a data frame with a 'chain' column.")
  }

  if (is.matrix(chains)) {
    return(compute_ess_matrix(chains))
  }

  if (is.data.frame(chains)) {
    if (!"chain" %in% names(chains)) {
      stop("Data frame must contain a 'chain' column.")
    }

    # Parameter columns are all except the 'chain' column.
    param_cols <- setdiff(names(chains), "chain")
    chain_ids <- unique(chains$chain)

    # Compute ESS for each parameter.
    ess_values <- sapply(param_cols, function(param) {
      param_data <- lapply(chain_ids, function(chain) {
        chains[[param]][chains$chain == chain]
      })

      # Ensure that all chains have the same number of iterations.
      chain_lengths <- sapply(param_data, length)
      if (length(unique(chain_lengths)) != 1) {
        stop(
          paste0(
            "Not all chains have the same number of iterations."
          )
        )
      }

      # Combine the data from each chain into a matrix.
      mat <- do.call(cbind, param_data)

      compute_ess_matrix(mat)
    })

    ess_values
  }
}
