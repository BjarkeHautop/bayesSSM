# Particle filter functions

The package provides several particle filter implementations for
state-space models for estimating the intractable marginal likelihood
\\p(y\_{1:T}\mid \theta)\\:

- [`auxiliary_filter`](https://bjarkehautop.github.io/bayesSSM/reference/auxiliary_filter.md)

- [`bootstrap_filter`](https://bjarkehautop.github.io/bayesSSM/reference/bootstrap_filter.md)

- [`resample_move_filter`](https://bjarkehautop.github.io/bayesSSM/reference/resample_move_filter.md)

The simplest one is the
[`bootstrap_filter`](https://bjarkehautop.github.io/bayesSSM/reference/bootstrap_filter.md),
and is thus recommended as a starting point.
