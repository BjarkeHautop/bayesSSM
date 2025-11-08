# Model Specification for Particle Filters

Model Specification for Particle Filters

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
