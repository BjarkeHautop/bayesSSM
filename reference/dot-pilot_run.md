# Pilot Run for Particle Filter Tuning

This internal function repeatedly evaluates the particle filter in order
to estimate the variance of the log-likelihoods and to compute a
recommended target number of particles for the Particle Marginal
Metropolis Hastings (PMMH) algorithm.

## Usage

``` r
.pilot_run(
  pf_wrapper,
  y,
  pilot_n,
  pilot_reps,
  init_fn,
  transition_fn,
  log_likelihood_fn,
  obs_times = NULL,
  resample_fn = NULL,
  ...
)
```

## Arguments

- pilot_n:

  An integer specifying the initial number of particles to use.

- pilot_reps:

  An integer specifying the number of repetitions for the pilot run.

## Value

A list containing:

- variance_estimate:

  The estimated variance of the log-likelihoods from the pilot run.

- target_N:

  The number of particles used in PMMH algorithm.

- pilot_loglikes:

  A numeric vector of log-likelihood values computed during the run.

## Details

The function performs `pilot_reps` evaluations of the particle filter
using the provided parameter vector `theta`. It then estimates the
variance of the log-likelihoods and scales the initial particle number
by this variance. The final number of particles is taken as the ceiling
of the scaled value with a minimum of 50 and a maximum of 1000.
