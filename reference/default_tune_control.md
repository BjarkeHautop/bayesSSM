# Create Tuning Control Parameters

This function creates a list of tuning parameters used by the
[`pmmh`](https://bjarkehautop.github.io/bayesSSM/reference/pmmh.md)
function. The tuning choices are inspired by Pitt et al. \[2012\] and
Dahlin and Schön \[2019\].

## Usage

``` r
default_tune_control(
  pilot_proposal_sd = 0.5,
  pilot_n = 100,
  pilot_m = 2000,
  pilot_target_var = 1,
  pilot_burn_in = 500,
  pilot_reps = 100,
  pilot_resample_algorithm = c("SISAR", "SISR", "SIS"),
  pilot_resample_fn = c("stratified", "systematic", "multinomial")
)
```

## Arguments

- pilot_proposal_sd:

  Standard deviation for pilot proposals. Default is 0.5.

- pilot_n:

  Number of pilot particles for particle filter. Default is 100.

- pilot_m:

  Number of iterations for MCMC. Default is 2000.

- pilot_target_var:

  The target variance for the posterior log-likelihood evaluated at
  estimated posterior mean. Default is 1.

- pilot_burn_in:

  Number of burn-in iterations for MCMC. Default is 500.

- pilot_reps:

  Number of times a particle filter is run. Default is 100.

- pilot_resample_algorithm:

  The resample_algorithm used for the pilot particle filter. Default is
  `"SISAR"`.

- pilot_resample_fn:

  The resampling function used for the pilot particle filter. Default is
  `"stratified"`.

## Value

A list of tuning control parameters.

## References

M. K. Pitt, R. d. S. Silva, P. Giordani, and R. Kohn. On some properties
of Markov chain Monte Carlo simulation methods based on the particle
filter. Journal of Econometrics, 171(2):134–151, 2012. doi:
https://doi.org/10.1016/j.jeconom.2012.06.004

J. Dahlin and T. B. Schön. Getting started with particle
Metropolis-Hastings for inference in nonlinear dynamical models. Journal
of Statistical Software, 88(2):1–41, 2019. doi: 10.18637/jss.v088.c02
