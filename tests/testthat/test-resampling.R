##################### Tests for cpp code #####################
test_that("Validation of weights work", {
  expect_error(
    resample_multinomial_cpp(3, c(-1, 1, 2)),
    "Weights must be non-negative"
  )
  expect_error(
    resample_stratified_cpp(3, c(-1, 1, 2)),
    "Weights must be non-negative"
  )
  expect_error(
    resample_stratified_cpp(3, c(-1, 1, 2)),
    "Weights must be non-negative"
  )

  expect_error(
    resample_multinomial_cpp(3, c(0, 0, 0)),
    "Sum of weights must be greater than 0"
  )
  expect_error(
    resample_stratified_cpp(3, c(0, 0, 0)),
    "Sum of weights must be greater than 0"
  )
  expect_error(
    resample_systematic_cpp(3, c(0, 0, 0)),
    "Sum of weights must be greater than 0"
  )
})
test_that("Resampling functions return correct proportions", {
  set.seed(1405)
  weights <- c(0.1, 0.2, 0.3, 0.2, 0.2)

  # Repeat many times
  n <- 10000
  indices_multinomial <- replicate(n, resample_multinomial_cpp(5, weights))
  indices_stratified <- replicate(n, resample_stratified_cpp(5, weights))
  indices_systematic <- replicate(n, resample_systematic_cpp(5, weights))
  # Calculate proportions
  prop_multinomial <- as.numeric(table(indices_multinomial)) / (n * 5)
  prop_stratified <- as.numeric(table(indices_stratified)) / (n * 5)
  prop_systematic <- as.numeric(table(indices_systematic)) / (n * 5)

  # Check proportions
  expect_equal(prop_multinomial, weights, tolerance = 0.05)
  expect_equal(prop_stratified, weights, tolerance = 0.05)
  expect_equal(prop_systematic, weights, tolerance = 0.05)
})
test_that("Check stratified and systematic correctly uses cumulative weights", {
  weights <- c(0.1, 0.5, 0.1, 0.15, 0.15)
  # For stratified resampling index 2 should always be for
  # 2nd and 3rd sample
  indices_stratified <- replicate(100, resample_stratified_cpp(5, weights))
  expect_true(all(indices_stratified[2, ] == 2))
  expect_true(all(indices_stratified[3, ] == 2))
  # For systematic resampling index 2 should always be for
  # 2nd and 3rd sample. If index 1 was selected for sample 1 then index 3 should
  # be selected for sample 4. If index 2 was selected for sample 1 then
  # index 4 should be selected for sample 4.
  indices_systematic <- replicate(100, resample_systematic_cpp(5, weights))
  expect_true(all(indices_systematic[2, ] == 2))
  expect_true(all(indices_systematic[3, ] == 2))

  first <- indices_systematic[1, ]
  fourth <- indices_systematic[4, ]

  expect_true(all(fourth[first == 1] == 3))
  expect_true(all(fourth[first == 2] == 4))
})



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
