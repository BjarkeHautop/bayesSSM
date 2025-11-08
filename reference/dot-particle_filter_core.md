# Core Particle Filter Function

This function implements the underlying logic used for particle filters
in a state space model using sequential Monte Carlo methods.

## Usage

``` r
.particle_filter_core(
  y,
  num_particles,
  init_fn,
  transition_fn,
  weight_fn,
  aux_weight_fn = NULL,
  move_fn = NULL,
  obs_times = NULL,
  algorithm = c("BPF", "APF", "RMPF"),
  resample_algorithm = c("SIS", "SISR", "SISAR"),
  resample_fn = c("stratified", "systematic", "multinomial"),
  threshold = NULL,
  return_particles = TRUE,
  ...
)
```

## Arguments

- y:

  A numeric vector or matrix of observations. Each row represents an
  observation at a time step. If observations are not equally spaced,
  use the `obs_times` argument.

- num_particles:

  A positive integer specifying the number of particles.

- init_fn:

  A function to initialize the particles. Should take \`num_particles\`
  and return a matrix or vector of initial states. Additional model
  parameters can be passed via `...`.

- transition_fn:

  A function for propagating particles. Should take \`particles\` and
  optionally \`t\`. Additional model parameters via `...`.

- weight_fn:

  A function that computes the log weights for the particles given the
  observations and the current particles. It should take \`y\`,
  \`particles\`, and \`t\` as arguments. The function can include any
  model-specific parameters as named arguments.

- obs_times:

  A numeric vector specifying observation time points. Must match the
  number of rows in `y`, or defaults to `1:nrow(y)`.

- resample_algorithm:

  A character string specifying the filtering resample algorithm:
  `"SIS"` for no resampling, `"SISR"` for resampling at every time step,
  or `"SISAR"` for adaptive resampling when ESS drops below `threshold`.
  Using `"SISR"` or `"SISAR"` to avoid weight degeneracy is recommended.
  Default is `"SISAR"`.

- resample_fn:

  A string indicating the resampling method: `"stratified"`,
  `"systematic"`, or `"multinomial"`. Default is `"stratified"`.

- threshold:

  A numeric value specifying the ESS threshold for `"SISAR"`. Defaults
  to `num_particles / 2`.

- return_particles:

  Logical; if `TRUE`, returns the full particle and weight histories.

- ...:

  Additional arguments passed to `init_fn`, `transition_fn`, and
  `log_likelihood_fn`.

## Value

A list with components:

- state_est:

  Estimated states over time (weighted mean of particles).

- ess:

  Effective sample size at each time step.

- loglike:

  Total log-likelihood.

- loglike_history:

  Log-likelihood at each time step.

- algorithm:

  The filtering algorithm used.

- particles_history:

  Matrix of particle states over time (if `return_particles = TRUE`).

- weights_history:

  Matrix of particle weights over time (if `return_particles = TRUE`).

## Model Specification

Particle filter implementations in this package assume a discrete-time
state-space model defined by:

- A sequence of latent states \\x_0, x_1, \ldots, x_T\\ evolving
  according to a Markov process.

- Observations \\y_1, \ldots, y_T\\ that are conditionally independent
  given the corresponding latent states.

The model is specified as: \$\$x_0 \sim \mu\_\theta\$\$ \$\$x_t \sim
f\_\theta(x_t \mid x\_{t-1}), \quad t = 1, \ldots, T\$\$ \$\$y_t \sim
g\_\theta(y_t \mid x_t), \quad t = 1, \ldots, T\$\$

where \\\theta\\ denotes model parameters passed via `...`.

The user provides the following functions:

- `init_fn`: draws from the initial distribution \\\mu\_\theta\\.

- `transition_fn`: generates or evaluates the transition density
  \\f\_\theta\\.

- `weight_fn`: evaluates the observation likelihood \\g\_\theta\\.
