test_that("Stationary gives less than 1.01", {
  set.seed(1405)
  m <- 4000
  chains <- matrix(rnorm(m), nrow = m / 4, ncol = 4)
  expect_lt(rhat(chains), 1.01)

  set.seed(1405)
  m <- 4000
  chains_df <- data.frame(
    chain = rep(1:4, each = m / 4),
    param1 = rnorm(m),
    param2 = rnorm(m)
  )
  expect_lt(rhat(chains_df)[[1]], 1.01)
  expect_lt(rhat(chains_df)[[2]], 1.01)
})

test_that("Non-stationary gives large Rhat", {
  # Create a matrix with one chain that is clearly non-convergent
  set.seed(1405)
  m <- 100
  chains <- matrix(c(rnorm(m / 2), rnorm(m / 2) + 10),
    nrow = m, ncol = 1
  )

  expect_gt(rhat(chains), 2)
})


test_that("rhat stops for non-valid input", {
  expect_error(rhat(list(1, 2, 3)),
               "Input must be a matrix or a data frame with a 'chain' column.")
  # Data frame with no 'chain' column
  expect_error(rhat(data.frame(a = c(1, 2, 3), b = c(4, 5, 6))),
               "Data frame must contain a 'chain' column.")
})

test_that("rhat 0 variance", {
  chains <- matrix(rep(1, 16), nrow = 4, ncol = 4)
  expect_warning(rhat(chains), "One or more chains have zero variance")
})

test_that("rhat too few observations", {
  chains <- matrix(rep(1, 2), nrow = 1, ncol = 2)
  expect_error(rhat(chains), "Number of iterations must be at least 2.")
})

test_that("not same number of iterations", {
  m <- 8
  chains_df <- data.frame(
    chain = c(1, 1, 1, 1, 1, 2, 2, 2),
    param1 = rnorm(m),
    param2 = rnorm(m)
  )
  expect_error(
    rhat(chains_df),
    "Not all chains have the same number of iterations"
  )
})

test_that("rhat odd number of iterations", {
  set.seed(1405)
  m <- 4004
  chains <- matrix(rnorm(m), nrow = m / 4, ncol = 4)
  expect_lt(rhat(chains), 1.01)
})

test_that("rhat odd number of iterations with data frame", {
  set.seed(1405)
  m <- 4004
  chains_df <- data.frame(
    chain = rep(1:4, each = m / 4),
    param1 = rnorm(m),
    param2 = rnorm(m)
  )
  expect_lt(rhat(chains_df)[[1]], 1.01)
  expect_lt(rhat(chains_df)[[2]], 1.01)
})
