# Estimate effective sample size (ESS) of MCMC chains.

Estimate effective sample size (ESS) of MCMC chains.

## Usage

``` r
ess(chains)
```

## Arguments

- chains:

  A matrix (iterations x chains) or a data.frame with a 'chain' column
  and parameter columns.

## Value

The estimated effective sample size (ESS) if given a matrix, or a named
vector of ESS values if given a data frame.

## Details

Uses the formula for ESS proposed by Vehtari et al. (2021).

## References

Vehtari et al. (2021). Rank-normalization, folding, and localization: An
improved R-hat for assessing convergence of MCMC. Available at:
https://doi.org/10.1214/20-BA1221

## Examples

``` r
# With a matrix:
chains <- matrix(rnorm(3000), nrow = 1000, ncol = 3)
ess(chains)
#> [1] 2752.906

# With a data frame:
chains_df <- data.frame(
  chain = rep(1:3, each = 1000),
  param1 = rnorm(3000),
  param2 = rnorm(3000)
)
ess(chains_df)
#>   param1   param2 
#> 2910.196 3000.000 
```
