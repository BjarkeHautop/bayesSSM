test_that("Statinoary gives less than 1.01", {
  m <- 3000
  chains <- matrix(rnorm(m), nrow = m / 3, ncol = 3)
  expect_lt(rhat(chains), 1.01)
})

test_that("Non-stationary gives Rhat > 1.01", {
  # Create a matrix with one chain that is clearly non-convergent
  m <- 100
  chains <- matrix(c(rep(1, m / 2), rep(100, m / 2)), nrow = m, ncol = 1)

  expect_warning(
    rhat(chains),
    "The split Rhat statistic is greater than 1.01"
  )
})
