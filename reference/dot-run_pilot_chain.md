# Run Pilot Chain for Posterior Estimation

Run Pilot Chain for Posterior Estimation

## Usage

``` r
.run_pilot_chain(
  pf_wrapper,
  y,
  pilot_m,
  pilot_n,
  pilot_reps,
  init_fn,
  transition_fn,
  log_likelihood_fn,
  log_priors,
  proposal_sd,
  obs_times = NULL,
  param_transform = NULL,
  pilot_init_params = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- pilot_m:

  An integer specifying the number of iterations for the pilot chain.

- pilot_n:

  An integer specifying the number of particles for the particle filter.

- pilot_reps:

  An integer specifying the number of repetitions for the pilot run.

- log_priors:

  A list of functions representing the log-priors for each model
  parameter.

- proposal_sd:

  A numeric vector specifying the standard deviations for the random
  walk proposal distribution for each parameter.

- param_transform:

  A character vector specifying the parameter transformations when
  proposing parameters using a random walk. Currently only supports
  "log" for log-transformation, "logit" for logit transformation, and
  "identity" for no transformation. Default is \`NULL\`, which
  correspond to no transformation ("identity).

- pilot_init_params:

  A numeric vector of initial parameter values. If \`NULL\`, it will
  default to a vector of ones. Default is \`NULL\`.

- ...:

  Additional arguments passed to the particle filter function.

## Value

A list containing:

- pilot_theta_mean:

  A numeric vector of the posterior mean of the parameters.

- pilot_theta_cov:

  A matrix of the posterior covariance (or variance if only one
  parameter).

- target_N:

  The estimated target number of particles for the PMMH algorithm.

- pilot_theta_chain:

  A matrix containing the chain of parameter values throughout the pilot
  run.

- pilot_loglike_chain:

  A vector containing the log-likelihood values associated with each
  iteration of the pilot chain.

## Details

This function runs a pilot chain to estimate the posterior mean and
covariance of the model parameters using a particle filter. The chain is
run for \`pilot_m\` iterations, with each iteration proposing new
parameters and evaluating their likelihood and prior. The chain is then
used to estimate the posterior mean and covariance, which are used to
tune the number of particles for the Particle Marginal Metropolis
Hastings (PMMH) algorithm.
