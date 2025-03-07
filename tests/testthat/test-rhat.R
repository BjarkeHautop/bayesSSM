test_that("Statinoary gives less than 1.01", {
  m <- 3000
  chains <- matrix(rnorm(m), nrow = m / 3, ncol = 3)
  expect_lt(rhat(chains), 1.01)
})

test_that("Non-stationary gives Rhat > 1.01", {
  # Create a matrix with one chain that is clearly non-convergent
  m <- 100
  chains <- matrix(c(rnorm(m / 2), rnorm(m / 2) + 10),
    nrow = m, ncol = 1
  )

  expect_gt(rhat(chains), 1.01)
})


test_that("rhat stops for non-matrix input", {
  expect_error(rhat(list(1, 2, 3)), "Input 'chains' must be a matrix")
  expect_error(rhat(c(1, 2, 3)), "Input 'chains' must be a matrix")
  expect_error(
    rhat(data.frame(a = c(1, 2, 3), b = c(4, 5, 6))),
    "Input 'chains' must be a matrix"
  )
})
