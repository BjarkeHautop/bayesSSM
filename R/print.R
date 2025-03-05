print.pmmh_output <- function(x, ...) {
  # Convert diagnostics to a data frame for readability
  diagnostics_df <- data.frame(
    Parameter = names(x$diagnostics),
    ESS = sapply(x$diagnostics, function(diag) diag$ess),
    Rhat = sapply(x$diagnostics, function(diag) diag$rhat)
  )

  # Print summary
  cat("PMMH Results Summary:\n")
  print(diagnostics_df, row.names = FALSE)

  # Return object invisibly to allow further use
  invisible(x)
}
