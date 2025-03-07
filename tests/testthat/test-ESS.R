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

test_that("ess stops for non-matrix input", {
  expect_error(ess(list(1, 2, 3)), "Input 'chains' must be a matrix")
  expect_error(ess(c(1, 2, 3)), "Input 'chains' must be a matrix")
  expect_error(ess(data.frame(a = c(1, 2, 3), b = c(4, 5, 6))),
               "Input 'chains' must be a matrix")
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
