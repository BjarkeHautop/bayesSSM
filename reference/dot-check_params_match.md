# Helper function to validate input of user-defined functions and priors

Helper function to validate input of user-defined functions and priors

## Usage

``` r
.check_params_match(
  init_fn,
  transition_fn,
  log_likelihood_fn,
  pilot_init_params,
  log_priors
)
```

## Arguments

- init_fn:

  A function to initialize the state-space model.

- transition_fn:

  A function that defines the state transition of the state-space model.

- log_likelihood_fn:

  A function that calculates the log-likelihood for the state-space
  model given latent states.

- pilot_init_params:

  A vector of initial parameter values.

- log_priors:

  A list of functions for computing the log-prior of each parameter.
