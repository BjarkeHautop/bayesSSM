# Common Parameters for Particle Filters

These parameters are shared by particle filter implementations such as
the bootstrap filter, auxiliary particle filter, and resample-move
particle filter.

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

- log_likelihood_fn:

  A function that returns the log-likelihood for each particle given the
  current observation, particles, and optionally \`t\`. Additional
  parameters via `...`.

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
