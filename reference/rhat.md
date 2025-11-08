# Compute split Rhat statistic

Compute split Rhat statistic

## Usage

``` r
rhat(chains)
```

## Arguments

- chains:

  A matrix (iterations x chains) or a data.frame with a 'chain' column
  and parameter columns.

## Value

Rhat value (matrix input) or named vector of Rhat values.

## Details

Uses the formula for split-Rhat proposed by Gelman et al. (2013).

## References

Gelman et al. (2013). Bayesian Data Analysis, 3rd Edition.

## Examples

``` r
# Example with matrix
chains <- matrix(rnorm(3000), nrow = 1000, ncol = 3)
rhat(chains)
#> [1] 1.000749
#' # Example with data frame
chains_df <- data.frame(
  chain = rep(1:3, each = 1000),
  param1 = rnorm(3000),
  param2 = rnorm(3000)
)
rhat(chains_df)
#> param1 param2 
#>      1      1 
```
