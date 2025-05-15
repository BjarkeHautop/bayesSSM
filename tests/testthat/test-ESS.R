test_that("ESS independent", {
  n <- 3000
  chains <- matrix(rnorm(n), nrow = n / 3, ncol = 3)
  expect_equal(ess(chains), n, tolerance = 0.05 * n)
})

test_that("ess autocorrelated", {
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

  # The ess should be lower than `n` due to autocorrelation
  expect_lt(ess(chains), n)
})

test_that("dataframe input", {
  n <- 3000
  chains_df <- data.frame(
    chain = rep(1:3, each = n / 3),
    param1 = rnorm(n),
    param2 = rnorm(n)
  )
  expect_equal(ess(chains_df)[[1]], n, tolerance = 0.05 * n)
  expect_equal(ess(chains_df)[[2]], n, tolerance = 0.05 * n)
})

test_that("ess stops for non-valid input", {
  expect_error(
    ess(list(1, 2, 3)),
    "Input must be a matrix or a data frame with a 'chain' column."
  )
})

test_that("ess stops for too few iterations", {
  chains <- matrix(rnorm(3), nrow = 1, ncol = 3)
  expect_error(ess(chains), "Number of iterations must be at least 2")
})

test_that("ess stops for too few chains", {
  chains <- matrix(rnorm(6), nrow = 6, ncol = 1)
  expect_error(ess(chains), "Number of chains must be at least 2")
})

test_that("ess stops for zero variance chains", {
  chains <- matrix(rep(1, 9), nrow = 3, ncol = 3)
  expect_warning(ess(chains), "One or more chains have zero variance")
})

test_that("ess stops dataframe no 'chain' column ", {
  chains_df <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))
  expect_error(
    ess(chains_df),
    "Data frame must contain a 'chain' column."
  )
})

test_that("ess stops if every chain not same number of iterations", {
  chains_df <- data.frame(
    chain = c(1, 1, 1, 1, 1, 2, 2, 2),
    param1 = rnorm(8),
    param2 = rnorm(8)
  )
  expect_error(
    ess(chains_df),
    "Not all chains have the same number of iterations."
  )
})
