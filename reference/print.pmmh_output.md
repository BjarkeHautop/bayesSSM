# Print method for PMMH output

Displays a concise summary of parameter estimates from a PMMH output
object, including means, standard deviations, medians, 95% credible
intervals, effective sample sizes (ESS), and Rhat. This provides a quick
overview of the posterior distribution and convergence diagnostics.

## Usage

``` r
# S3 method for class 'pmmh_output'
print(x, ...)
```

## Arguments

- x:

  An object of class \`pmmh_output\`.

- ...:

  Additional arguments.

## Value

The object \`x\` invisibly.

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
print(dummy_output)
#> PMMH Results Summary:
#>  Parameter Mean   SD Median  2.5% 97.5% ESS Rhat
#>     param1 0.10 1.03   0.09 -1.88  1.93 200 1.01
#>     param2 0.04 0.94   0.02 -1.71  2.08 190 1.00
```
