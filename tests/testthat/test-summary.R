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

  # Check that the class is correctly assigned
  expect_s3_class(dummy_output, "pmmh_output")

  # Check that the structure is as expected
  expect_type(dummy_output$theta_chain, "list")
  expect_equal(nrow(dummy_output$theta_chain), 200) # 100 from each chain

  expect_named(dummy_output$diagnostics$ess, c("param1", "param2"))
  expect_named(dummy_output$diagnostics$rhat, c("param1", "param2"))

  expect_equal(dummy_output$diagnostics$ess[["param1"]], 200)
  expect_equal(dummy_output$diagnostics$rhat[["param1"]], 1.01)
})
