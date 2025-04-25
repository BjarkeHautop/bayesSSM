test_that("Summary works", {
  chain1 <- data.frame(param1 = rnorm(100), param2 = rnorm(100), chain = 1)
  chain2 <- data.frame(param1 = rnorm(100), param2 = rnorm(100), chain = 2)

  dummy_output <- list(
    theta_chain = rbind(chain1, chain2),
    diagnostics = list(
      ess = c(param1 = 200, param2 = 190),
      rhat = c(param1 = 1.01, param2 = 1.00)
    )
  )

  class(dummy_output) <- "pmmh_output"

  summary_out <- summary(dummy_output)

  # Check structure and content
  expect_s3_class(summary_out, "data.frame")
  expect_named(
    summary_out,
    c("mean", "sd", "median", "2.5%", "97.5%", "ESS", "Rhat")
  )
  expect_equal(rownames(summary_out), c("param1", "param2"))
  expect_equal(summary_out["param1", "ESS"], 200)
  expect_equal(summary_out["param1", "Rhat"], 1.01)
})
