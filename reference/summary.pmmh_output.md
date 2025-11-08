# Summary method for PMMH output

This function returns summary statistics for PMMH output objects,
including means, standard deviations, medians, credible intervals, and
diagnostics.

## Usage

``` r
# S3 method for class 'pmmh_output'
summary(object, ...)
```

## Arguments

- object:

  An object of class \`pmmh_output\`.

- ...:

  Additional arguments.

## Value

A data frame containing summary statistics for each parameter.

## Examples

``` r
# Create dummy chains for two parameters across two chains
chain1 <- data.frame(param1 = rnorm(100), param2 = rnorm(100), chain = 1)
chain2 <- data.frame(param1 = rnorm(100), param2 = rnorm(100), chain = 2)
dummy_output <- list(
  theta_chain = rbind(chain1, chain2),
  diagnostics = list(
    ess = c(param1 = 200, param2 = 190),
    rhat = c(param1 = 1.01, param2 = 1.00)
  )
)
class(dummy_output) <- "pmmh_output"
summary(dummy_output)
#>               mean        sd      median      2.5%    97.5% ESS Rhat
#> param1  0.01499914 1.0339426 -0.07112074 -1.923310 2.193800 200 1.01
#> param2 -0.06911637 0.9336208 -0.09527880 -1.661481 1.571068 190 1.00
```
