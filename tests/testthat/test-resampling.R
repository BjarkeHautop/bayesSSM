test_that("Throws error if particles dim doesn't match weight", {
  particles <- 1:3
  weights <- c(0.1, 0.2, 0.3, 0.2)
  expect_error(
    .resample_multinomial(particles, weights),
    "Number of particles must match the length of weights"
  )
  expect_error(
    .resample_stratified(particles, weights),
    "Number of particles must match the length of weights"
  )
  expect_error(
    .resample_systematic(particles, weights),
    "Number of particles must match the length of weights"
  )

  # Matrix case
  particles <- matrix(1:6, nrow = 3)
  weights <- c(0.1, 0.2, 0.3, 0.2)
  expect_error(
    .resample_multinomial(particles, weights),
    "Number of particles must match the length of weights"
  )
  expect_error(
    .resample_stratified(particles, weights),
    "Number of particles must match the length of weights"
  )
  expect_error(
    .resample_systematic(particles, weights),
    "Number of particles must match the length of weights"
  )
})


test_that("Throws error for non-valid weights", {
  particles <- 1:3
  weights <- rep(0, 3)
  expect_error(
    .resample_multinomial(particles, weights),
    "Sum of weights must be greater than 0"
  )
  expect_error(
    .resample_stratified(particles, weights),
    "Sum of weights must be greater than 0"
  )
  expect_error(
    .resample_systematic(particles, weights),
    "Sum of weights must be greater than 0"
  )

  weights <- c(-0.1, 0.5, 0.4)
  expect_error(
    .resample_multinomial(particles, weights),
    "Weights must be non-negative"
  )
  expect_error(
    .resample_stratified(particles, weights),
    "Weights must be non-negative"
  )
  expect_error(
    .resample_systematic(particles, weights),
    "Weights must be non-negative"
  )
})

test_that("Multinomial resampling produces valid output", {
  set.seed(123)
  particles <- 1:5
  weights <- c(0.1, 0.2, 0.3, 0.2, 0.2)

  resampled <- .resample_multinomial(particles, weights)

  # Check that all resampled values are from the original set
  expect_true(all(resampled %in% particles))

  # Check that the length is preserved
  expect_equal(length(resampled), length(particles))
})

test_that("Stratified resampling produces valid output", {
  set.seed(123)
  particles <- 1:5
  weights <- c(0.1, 0.2, 0.3, 0.2, 0.2)

  resampled <- .resample_stratified(particles, weights)

  expect_true(all(resampled %in% particles))
  expect_equal(length(resampled), length(particles))
})

test_that("Systematic resampling produces valid output", {
  set.seed(123)
  particles <- 1:5
  weights <- c(0.1, 0.2, 0.3, 0.2, 0.2)

  resampled <- .resample_systematic(particles, weights)

  expect_true(all(resampled %in% particles))
  expect_equal(length(resampled), length(particles))
})

test_that("Resampling handles uniform weights correctly", {
  set.seed(123)
  particles <- 1:10
  weights <- rep(1 / 10, 10)

  resampled_multinomial <- .resample_multinomial(particles, weights)
  resampled_stratified <- .resample_stratified(particles, weights)
  resampled_systematic <- .resample_systematic(particles, weights)

  expect_true(all(resampled_multinomial %in% particles))
  expect_true(all(resampled_stratified %in% particles))
  expect_true(all(resampled_systematic %in% particles))

  expect_equal(length(resampled_multinomial), length(particles))
  expect_equal(length(resampled_stratified), length(particles))
  expect_equal(length(resampled_systematic), length(particles))
})

test_that("Resampling handles extreme weights correctly", {
  set.seed(123)
  particles <- 1:5
  weights <- c(0, 0, 1, 0, 0) # All weight on one particle

  resampled_multinomial <- .resample_multinomial(particles, weights)
  resampled_stratified <- .resample_stratified(particles, weights)
  resampled_systematic <- .resample_systematic(particles, weights)

  expect_true(all(resampled_multinomial == 3))
  expect_true(all(resampled_stratified == 3))
  expect_true(all(resampled_systematic == 3))
})

test_that("Works with matrix particles", {
  set.seed(123)
  particles <- matrix(1:6, nrow = 3)
  weights <- rep(1 / 3, 3)

  resampled_multinomial <- .resample_multinomial(particles, weights)
  resampled_stratified <- .resample_stratified(particles, weights)
  resampled_systematic <- .resample_systematic(particles, weights)

  expect_true(all(resampled_multinomial %in% particles))
  expect_true(all(resampled_stratified %in% particles))
  expect_true(all(resampled_systematic %in% particles))

  expect_equal(nrow(resampled_multinomial), nrow(particles))
  expect_equal(nrow(resampled_stratified), nrow(particles))
  expect_equal(nrow(resampled_systematic), nrow(particles))
})
