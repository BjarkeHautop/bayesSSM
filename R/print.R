#' Print method for PMMH output
#'
#' Displays a concise summary of parameter estimates from a PMMH output object,
#' including means, standard deviations, medians, 95\% credible intervals,
#' effective sample sizes (ESS), and Rhat. This provides
#' a quick overview of the posterior distribution and convergence diagnostics.
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
#' print(dummy_output)
print.pmmh_output <- function(x, ...) {
  param_names <- colnames(x$theta_chain)
  param_names <- param_names[!(param_names %in% c("chain"))]

  summary_stats <- t(sapply(param_names, function(param) {
    samples <- x$theta_chain[[param]]
    mean_val <- mean(samples)
    sd_val <- sd(samples)
    median_val <- median(samples)
    ci <- as.numeric(quantile(samples, c(0.025, 0.975)))

    c(
      Mean = round(mean_val, 2),
      SD = round(sd_val, 2),
      Median = round(median_val, 2),
      `2.5%` = round(ci[1], 2),
      `97.5%` = round(ci[2], 2)
    )
  }))

  ess_values <- floor(unlist(x$diagnostics$ess))
  rhat_values <- round(unlist(x$diagnostics$rhat), 3)

  diagnostics_df <- data.frame(
    Parameter = param_names,
    summary_stats,
    ESS = ess_values[param_names],
    Rhat = rhat_values[param_names],
    row.names = NULL,
    check.names = FALSE
  )

  cat("PMMH Results Summary:\n")
  print(diagnostics_df, row.names = FALSE)

  invisible(x)
}
