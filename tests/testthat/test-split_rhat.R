test_that("Independent gives less than 1.01", {
  m <- 3000
  chains <- matrix(rnorm(m), nrow = m/3, ncol = 3)
  expect_lt(split_rhat(chains), 1.01)
})

test_that("split_rhat warns when Rhat > 1.01", {
  # Create a matrix with one chain that is clearly non-convergent
  m <- 100
  chains <- matrix(c(rep(1, m/2), rep(100, m/2)), nrow = m, ncol = 1)

  expect_warning(
    split_rhat(chains),
    "The split Rhat statistic is greater than 1.01"
  )
})
