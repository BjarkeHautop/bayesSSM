# bayesSSM: Bayesian Inference for State-Space Models

The bayesSSM package provides implementations of particle filtering,
Particle MCMC, and related methods for Bayesian inference in state-space
models. It includes tools for simulation, posterior inference, and
diagnostics.

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

## See also

Useful links:

- <https://github.com/BjarkeHautop/bayesSSM>

- <https://bjarkehautop.github.io/bayesSSM/>

- Report bugs at <https://github.com/BjarkeHautop/bayesSSM/issues>

## Author

**Maintainer**: Bjarke Hautop <bjarke.hautop@gmail.com> \[copyright
holder\]
