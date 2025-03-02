test_that("ESS independent", {
  n <- 3000
  chains <- matrix(rnorm(n), nrow = n / 3, ncol = 3)
  expect_equal(ESS(chains), n, tolerance = 0.05 * n)
})

test_that("ESS autocorrelated", {
  n <- 3000
  rho <- 0.9
  chains <- matrix(0, nrow = n / 3, ncol = 3)

  # Generate AR(1) process
  for (j in 1:3) {
    chains[1, j] <- rnorm(1)
    for (i in 2:(n / 3)) {
      chains[i, j] <- rho * chains[i - 1, j] + rnorm(1)
    }
  }

  # The ESS should be lower than `n` due to autocorrelation
  expect_lt(ESS(chains), n)
})
