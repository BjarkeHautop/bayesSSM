test_that("Statinoary gives less than 1.01", {
  m <- 3000
  chains <- matrix(rnorm(m), nrow = m / 3, ncol = 3)
  expect_lt(rhat(chains), 1.01)
})

test_that("Non-stationary gives Rhat > 1.01", {
  # Create a matrix with one chain that is clearly non-convergent
  m <- 100
  chains <- matrix(c(rnorm(m/2), rnorm(m/2) + 10),
                   nrow = m, ncol = 1)

  expect_gt(rhat(chains), 1.01)
})
