#' Print method for PMMH output
#'
#' @param x An object of class `pmmh_output`.
#' @param ... Additional arguments.
#'
#' @returns The object `x` invisibly.
#'
#' @importFrom stats sd median quantile
#'
#' @export
#'
#' @examples
#' # Create dummy chains for two parameters across two chains
#' set.seed(1405)
#' chain1 <- data.frame(param1 = rnorm(100), param2 = rnorm(100))
#' chain2 <- data.frame(param1 = rnorm(100), param2 = rnorm(100))
#' dummy_output <- list(
#'   theta_chain = list(chain1, chain2),
#'   diagnostics = list(
#'     ess = c(param1 = 200, param2 = 190),
#'     rhat = c(param1 = 1.01, param2 = 1.00)
#'   )
#' )
#' class(dummy_output) <- "pmmh_output"
#' print(dummy_output)
print.pmmh_output <- function(x, ...) {
  # Extract parameter names from the first chain's columns
  param_names <- colnames(x$theta_chain[[1]])

  # Compute posterior summaries for each parameter by combining chains
  summary_stats <- t(sapply(param_names, function(param) {
    # Combine samples for the parameter across all chains
    samples <- unlist(lapply(x$theta_chain, function(chain) chain[, param]))

    # Compute summary statistics
    mean_val <- mean(samples)
    sd_val <- sd(samples)
    median_val <- median(samples)
    ci_lower <- quantile(samples, 0.025)
    ci_upper <- quantile(samples, 0.975)

    # Round these values to 2 decimal places
    c(
      Mean = round(mean_val, 2),
      SD = round(sd_val, 2),
      Median = round(median_val, 2),
      CI = round(ci_lower, 2),
      CI = round(ci_upper, 2)
    )
  }))

  # Extract ESS and Rhat diagnostics (assumed to be named vectors/lists)
  ess_values <- unlist(x$diagnostics$ess)
  rhat_values <- unlist(x$diagnostics$rhat)

  # Round ESS as integers and Rhat to 3 decimal places
  ess_values <- floor(ess_values)
  rhat_values <- round(rhat_values, 3)

  # Create a data frame with all the summary statistics and diagnostics
  diagnostics_df <- data.frame(
    Parameter = param_names,
    summary_stats,
    ESS = ess_values[param_names],
    Rhat = rhat_values[param_names],
    row.names = NULL,
    check.names = FALSE
  )

  # Print the summary table
  cat("PMMH Results Summary:\n")
  print(diagnostics_df, row.names = FALSE)

  # Return the object invisibly
  invisible(x)
}
