#' Compute split Rhat statistic
#'
#' @param chains A matrix (iterations x chains) or a data.frame with a 'chain'
#' column and parameter columns.
#'
#' @return Rhat value (matrix input) or named vector of Rhat values.
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
#' # Example with matrix
#' chains <- matrix(rnorm(3000), nrow = 1000, ncol = 3)
#' rhat(chains)
#' #' # Example with data frame
#' chains_df <- data.frame(
#'   chain = rep(1:3, each = 1000),
#'   param1 = rnorm(3000),
#'   param2 = rnorm(3000)
#' )
#' rhat(chains_df)
rhat <- function(chains) {
  compute_rhat_matrix <- function(mat) {
    m <- nrow(mat)
    k <- ncol(mat)
    if (m < 2) {
      stop("Number of iterations must be at least 2.")
    }

    # Ensure even number of iterations
    if (m %% 2 == 1) {
      mat <- mat[-m, , drop = FALSE] # Drop last iteration
      m <- nrow(mat)
    }
    # Split each chain into two parts
    chains_split <- matrix(NA, nrow = m %/% 2, ncol = 2 * k)
    for (i in 1:k) {
      # Split the k-th chain into two parts
      chains_split[, 2 * i - 1] <- mat[1:(m %/% 2), i]
      chains_split[, 2 * i]     <- mat[(m %/% 2 + 1):m, i]
    }

    chain_means <- colMeans(chains_split)
    overall_mean <- mean(chain_means)

    b <- m / (2 * k - 1) * sum((chain_means - overall_mean)^2)
    chain_vars <- apply(chains_split, 2, var)

    if (any(chain_vars == 0)) {
      warning("One or more chains have zero variance.")
      return(NA)
    }

    w <- mean(chain_vars)
    var_hat <- ((m - 1) / m) * w + (1 / m) * b
    r_hat <- sqrt(var_hat / w)

    if (r_hat >= 0.99 && r_hat <= 1) {
      r_hat <- 1
    }
    r_hat
  }

  # Check if input is a matrix or a data frame
  if (!is.matrix(chains) && !is.data.frame(chains)) {
    stop("Input must be a matrix or a data frame with a 'chain' column.")
  }

  if (is.matrix(chains)) {
    return(compute_rhat_matrix(chains))
  }

  if (is.data.frame(chains)) {
    if (!"chain" %in% names(chains)) {
      stop("Data frame must contain a 'chain' column.")
    }

    param_cols <- setdiff(names(chains), "chain")
    chain_ids <- unique(chains$chain)

    rhat_values <- sapply(param_cols, function(param) {
      param_data <- lapply(chain_ids, function(chain) {
        chains[[param]][chains$chain == chain]
      })

      # Ensure that all chains have the same number of iterations.
      chain_lengths <- sapply(param_data, length)
      if (length(unique(chain_lengths)) != 1) {
        stop(
          paste0(
            "Not all chains have the same number of iterations for
            parameter: ", param
          )
        )
      }
      mat <- do.call(cbind, param_data)
      compute_rhat_matrix(mat)
    })

    rhat_values
  }
}
