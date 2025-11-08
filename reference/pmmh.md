# Particle Marginal Metropolis-Hastings (PMMH) for State-Space Models

This function implements a Particle Marginal Metropolis-Hastings (PMMH)
resample_algorithm to perform Bayesian inference in state-space models.
It first runs a pilot chain to tune the proposal distribution and the
number of particles for the particle filter, and then runs the main PMMH
chain.

## Usage

``` r
pmmh(
  pf_wrapper,
  y,
  m,
  init_fn,
  transition_fn,
  log_likelihood_fn,
  log_priors,
  pilot_init_params,
  burn_in,
  num_chains = 4,
  obs_times = NULL,
  resample_algorithm = c("SISAR", "SISR", "SIS"),
  resample_fn = c("stratified", "systematic", "multinomial"),
  param_transform = NULL,
  tune_control = default_tune_control(),
  verbose = FALSE,
  return_latent_state_est = FALSE,
  seed = NULL,
  num_cores = 1,
  ...
)
```

## Arguments

- pf_wrapper:

  A particle filter wrapper function. See
  [`particle_filter`](https://bjarkehautop.github.io/bayesSSM/reference/particle_filter.md)
  for the particle filters implemented in this package.

- y:

  A numeric vector or matrix of observations. Each row represents an
  observation at a time step. If observations are not equally spaced,
  use the `obs_times` argument.

- m:

  An integer specifying the number of MCMC iterations for each chain.

- init_fn:

  A function to initialize the particles. Should take \`num_particles\`
  and return a matrix or vector of initial states. Additional model
  parameters can be passed via `...`.

- transition_fn:

  A function for propagating particles. Should take \`particles\` and
  optionally \`t\`. Additional model parameters via `...`.

- log_likelihood_fn:

  A function that returns the log-likelihood for each particle given the
  current observation, particles, and optionally \`t\`. Additional
  parameters via `...`.

- log_priors:

  A list of functions for computing the log-prior of each parameter.

- pilot_init_params:

  A list of initial parameter values. Should be a list of length
  `num_chains` where each element is a named vector of initial parameter
  values.

- burn_in:

  An integer indicating the number of initial MCMC iterations to discard
  as burn-in.

- num_chains:

  An integer specifying the number of PMMH chains to run.

- obs_times:

  A numeric vector specifying observation time points. Must match the
  number of rows in `y`, or defaults to `1:nrow(y)`.

- resample_algorithm:

  A character string specifying the resampling algorithm to use in the
  particle filter. Options are: \#'

  - SIS: Sequential Importance Sampling (without resampling).

  - SISR: Sequential Importance Sampling with resampling at every time
    step.

  - SISAR: SIS with adaptive resampling based on the Effective Sample
    Size (ESS). Resampling is triggered when the ESS falls below a given
    threshold (default `particles / 2`). Can be modified by specifying
    the `threshold` argument (in `...`), see also
    [`particle_filter`](https://bjarkehautop.github.io/bayesSSM/reference/particle_filter.md).

- resample_fn:

  A string indicating the resampling method: `"stratified"`,
  `"systematic"`, or `"multinomial"`. Default is `"stratified"`.

- param_transform:

  An optional character vector that specifies the transformation applied
  to each parameter before proposing. The proposal is made using a
  multivariate normal distribution on the transformed scale. Parameters
  are then mapped back to their original scale before evaluation.
  Currently supports `"log"`, `"logit"`, and `"identity"`. If `NULL`,
  the `"identity"` transformation is used for all parameters.

- tune_control:

  A list of pilot tuning controls (e.g., `pilot_m`, `pilot_reps`). See
  [`default_tune_control`](https://bjarkehautop.github.io/bayesSSM/reference/default_tune_control.md).

- verbose:

  A logical value indicating whether to print information about
  pilot_run tuning. Defaults to `FALSE`.

- return_latent_state_est:

  A logical value indicating whether to return the latent state
  estimates for each time step. Defaults to `FALSE`.

- seed:

  An optional integer to set the seed for reproducibility.

- num_cores:

  An integer specifying the number of cores to use for parallel
  processing. Defaults to 1. Each chain is assigned to its own core, so
  the number of cores cannot exceed the number of chains (`num_chains`).
  The progress information given to user is limited if using more than
  one core.

- ...:

  Additional arguments passed to `init_fn`, `transition_fn`, and
  `log_likelihood_fn`.

## Value

A list containing:

- `theta_chain`:

  A dataframe of post burn-in parameter samples.

- `latent_state_chain`:

  If `return_latent_state_est` is `TRUE`, a list of matrices containing
  the latent state estimates for each time step.

- `diagnostics`:

  Diagnostics containing ESS and Rhat for each parameter (see
  [`ess`](https://bjarkehautop.github.io/bayesSSM/reference/ess.md) and
  [`rhat`](https://bjarkehautop.github.io/bayesSSM/reference/rhat.md)
  for documentation).

## Details

The PMMH resample_algorithm is essentially a Metropolis Hastings
algorithm, where instead of using the intractable marginal likelihood
\\p(y\_{1:T}\mid \theta)\\ it instead uses the estimated likelihood
using a particle filter (see
[`particle_filter`](https://bjarkehautop.github.io/bayesSSM/reference/particle_filter.md)
for available particle filters). Values are proposed using a
multivariate normal distribution in the transformed space (specified
using \`param_transform\`). The proposal covariance and the number of
particles is chosen based on a pilot run. The number of particles is
chosen such that the variance of the log-likelihood estimate at the
estimated posterior mean is approximately 1 (with a minimum of 50
particles and a maximum of 1000).

## References

Andrieu et al. (2010). Particle Markov chain Monte Carlo methods.
Journal of the Royal Statistical Society: Series B (Statistical
Methodology), 72(3):269–342. doi: 10.1111/j.1467-9868.2009.00736.x

## Examples

``` r
init_fn <- function(num_particles) {
  rnorm(num_particles, mean = 0, sd = 1)
}
transition_fn <- function(particles, phi, sigma_x) {
  phi * particles + sin(particles) +
    rnorm(length(particles), mean = 0, sd = sigma_x)
}
log_likelihood_fn <- function(y, particles, sigma_y) {
  dnorm(y, mean = cos(particles), sd = sigma_y, log = TRUE)
}
log_prior_phi <- function(phi) {
  dnorm(phi, mean = 0, sd = 1, log = TRUE)
}
log_prior_sigma_x <- function(sigma) {
  dexp(sigma, rate = 1, log = TRUE)
}
log_prior_sigma_y <- function(sigma) {
  dexp(sigma, rate = 1, log = TRUE)
}
log_priors <- list(
  phi = log_prior_phi,
  sigma_x = log_prior_sigma_x,
  sigma_y = log_prior_sigma_y
)
# Generate data
t_val <- 10
x <- numeric(t_val)
y <- numeric(t_val)
phi <- 0.8
sigma_x <- 1
sigma_y <- 0.5

init_state <- rnorm(1, mean = 0, sd = 1)
x[1] <- phi * init_state + sin(init_state) + rnorm(1, mean = 0, sd = sigma_x)
y[1] <- x[1] + rnorm(1, mean = 0, sd = sigma_y)
for (t in 2:t_val) {
  x[t] <- phi * x[t - 1] + sin(x[t - 1]) + rnorm(1, mean = 0, sd = sigma_x)
  y[t] <- cos(x[t]) + rnorm(1, mean = 0, sd = sigma_y)
}
x <- c(init_state, x)

# Should use higher MCMC iterations in practice (m)
pmmh_result <- pmmh(
  pf_wrapper = bootstrap_filter,
  y = y,
  m = 1000,
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  log_priors = log_priors,
  pilot_init_params = list(
    c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
    c(phi = 1, sigma_x = 0.5, sigma_y = 1)
  ),
  burn_in = 100,
  num_chains = 2,
  param_transform = list(
    phi = "identity",
    sigma_x = "log",
    sigma_y = "log"
  ),
  tune_control = default_tune_control(pilot_m = 500, pilot_burn_in = 100)
)
#> Running chain 1...
#> Running pilot chain for tuning...
#> Using 50 particles for PMMH:
#> Running Particle MCMC chain with tuned settings...
#> Running chain 2...
#> Running pilot chain for tuning...
#> Using 50 particles for PMMH:
#> Running Particle MCMC chain with tuned settings...
#> PMMH Results Summary:
#>  Parameter Mean   SD Median  2.5% 97.5% ESS  Rhat
#>        phi 0.40 1.15   0.69 -2.30  2.61   8 1.103
#>    sigma_x 1.25 1.20   0.86  0.00  4.28  66 1.050
#>    sigma_y 1.08 0.33   1.02  0.61  1.86 116 1.001
#> Warning: Some ESS values are below 400, indicating poor mixing. Consider running the chains for more iterations.
#> Warning: 
#> Some Rhat values are above 1.01, indicating that the chains have not converged. 
#> Consider running the chains for more iterations and/or increase burn_in.
# Convergence warning is expected with such low MCMC iterations.

# Suppose we have data for t=1,2,3,5,6,7,8,9,10 (i.e., missing at t=4)

obs_times <- c(1, 2, 3, 5, 6, 7, 8, 9, 10)
y <- y[obs_times]

# Specify observation times in the pmmh using obs_times
pmmh_result <- pmmh(
  pf_wrapper = bootstrap_filter,
  y = y,
  m = 1000,
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  log_priors = log_priors,
  pilot_init_params = list(
    c(phi = 0.8, sigma_x = 1, sigma_y = 0.5),
    c(phi = 1, sigma_x = 0.5, sigma_y = 1)
  ),
  burn_in = 100,
  num_chains = 2,
  obs_times = obs_times,
  param_transform = list(
    phi = "identity",
    sigma_x = "log",
    sigma_y = "log"
  ),
  tune_control = default_tune_control(pilot_m = 500, pilot_burn_in = 100)
)
#> Running chain 1...
#> Running pilot chain for tuning...
#> Using 50 particles for PMMH:
#> Running Particle MCMC chain with tuned settings...
#> Running chain 2...
#> Running pilot chain for tuning...
#> Using 50 particles for PMMH:
#> Running Particle MCMC chain with tuned settings...
#> PMMH Results Summary:
#>  Parameter Mean   SD Median  2.5% 97.5% ESS  Rhat
#>        phi 0.43 1.04   0.73 -1.89  1.93  20 1.106
#>    sigma_x 1.20 1.02   0.88  0.08  3.98 125 1.020
#>    sigma_y 1.08 0.33   1.02  0.63  1.87 151 1.007
#> Warning: Some ESS values are below 400, indicating poor mixing. Consider running the chains for more iterations.
#> Warning: 
#> Some Rhat values are above 1.01, indicating that the chains have not converged. 
#> Consider running the chains for more iterations and/or increase burn_in.
```
