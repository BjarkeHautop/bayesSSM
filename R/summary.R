#' Summary method for PMMH output
#'
#' This function returns summary statistics for PMMH output objects, including
#' means, standard deviations, medians, credible intervals, and diagnostics.
#'
#' @param object An object of class `pmmh_output`.
#' @param ... Additional arguments.
#'
#' @returns A data frame containing summary statistics for each parameter.
#'
#' @importFrom stats sd median quantile
#'
#' @export
#'
#' @examples
#' # Create dummy chains for two parameters across two chains
#' chain1 <- data.frame(param1 = rnorm(100), param2 = rnorm(100), chain = 1)
#' chain2 <- data.frame(param1 = rnorm(100), param2 = rnorm(100), chain = 2)
#' dummy_output <- list(
#'   theta_chain = rbind(chain1, chain2),
#'   diagnostics = list(
#'     ess = c(param1 = 200, param2 = 190),
#'     rhat = c(param1 = 1.01, param2 = 1.00)
#'   )
#' )
#' class(dummy_output) <- "pmmh_output"
#' summary(dummy_output)
summary.pmmh_output <- function(object, ...) {
  # Extract parameter names, excluding 'chain' column
  param_names <- setdiff(colnames(object$theta_chain), "chain")

  # Compute summary statistics
  summary_stats <- sapply(param_names, function(param) {
    samples <- object$theta_chain[[param]]
    stats <- c(
      mean = mean(samples),
      sd = sd(samples),
      median = median(samples),
      quantile(samples, probs = c(0.025, 0.975))
    )
    names(stats)[4:5] <- c("2.5%", "97.5%") # clean names
    stats
  })


  # Convert to data frame
  summary_df <- as.data.frame(t(summary_stats))

  # Add ESS and Rhat diagnostics
  summary_df$ESS <- object$diagnostics$ess[param_names]
  summary_df$Rhat <- object$diagnostics$rhat[param_names]

  summary_df
}
