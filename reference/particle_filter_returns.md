# Shared Return Values for Particle Filters

This block documents the common return value for particle filtering
functions.

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
