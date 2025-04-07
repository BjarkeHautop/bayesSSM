
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesSSM <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/BjarkeHautop/bayesSSM/graph/badge.svg)](https://app.codecov.io/gh/BjarkeHautop/bayesSSM)
[![R-CMD-check](https://github.com/BjarkeHautop/bayesSSM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BjarkeHautop/bayesSSM/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

bayesSSM is an R package offering a set of tools for performing Bayesian
inference in state-space models (SSMs). It implements the Particle
Marginal Metropolis-Hastings (PMMH) for Bayesian inference.

## Why bayesSSM?

While there are several alternatives available for performing particle
MCMC, such as the [POMP](https://kingaa.github.io/pomp/) package, I
designed bayesSSM with ease of use in mind. It was developed as a
procrastination task during my Master’s thesis about Particle MCMC,
since I was implementing everything from scratch anyway. Everything is
written in R, so performance is not the best, but it is easy to use.

## Installation

You can install the development version of bayesSSM from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("BjarkeHautop/bayesSSM")
```

## Example

Consider the following SSM:

$$
\begin{aligned}
        X_1 &\sim N(0,1) \\
        X_t&=\phi X_{t-1}+\sin(X_{t-1})+\sigma_x V_t, \quad V_t \sim N(0,1) \\
        Y_t&=X_t+\sigma_y W_t, \quad W_t \sim N(0, \, 1).
\end{aligned}
$$

Let’s first simulate some data from this model with $\phi = 0.8$,
$\sigma_x = 1$, and $\sigma_y = 0.5$.

``` r
t_val <- 20
phi_val <- 0.8
sigma_x_val <- 1
sigma_y_val <- 0.5

x <- numeric(t_val)
y <- numeric(t_val)
x[1] <- rnorm(1, mean = 0, sd = sigma_x_val)
y[1] <- rnorm(1, mean = x[1], sd = sigma_y_val)
for (t in 2:t_val) {
  x[t] <- phi_val * x[t - 1] + sin(x[t - 1]) + rnorm(1,
    mean = 0,
    sd = sigma_x_val
  )
  y[t] <- x[t] + rnorm(1, mean = 0, sd = sigma_y_val)
}
```

We define the priors for our model as follows:

$$
\begin{aligned}
        \phi &\sim N(0,1), \\
        \sigma_x &\sim \text{Exp}(1), \\
        \sigma_y &\sim \text{Exp}(1).
\end{aligned}
$$

We can use `pmmh` to perform Bayesian inference on this model. To use
`pmmh` we need to define the functions for the SSM and the priors. The
functions `init_fn`, `transition_fn` should be functions that simulates
the latent states. They must contain the argument `particles`, which is
a vector of particles, and can contain any other arguments. The function
`log_likelihood_fn` should be a function that calculates the
log-likelihood of the observed data given the latent state variables. It
must contain the arguments `y` and `particles`.

The priors for the parameters must be defined as log-prior functions.
Every parameter from `init_fn`, `transition_fn`, and `log_likelihood_fn`
must have a corresponding log-prior function.

``` r
init_fn <- function(particles) {
  stats::rnorm(particles, mean = 0, sd = 1)
}
transition_fn <- function(particles, phi, sigma_x) {
  phi * particles + sin(particles) +
    stats::rnorm(length(particles), mean = 0, sd = sigma_x)
}
log_likelihood_fn <- function(y, particles, sigma_y) {
  stats::dnorm(y, mean = particles, sd = sigma_y, log = TRUE)
}
log_prior_phi <- function(phi) {
  stats::dnorm(phi, mean = 0, sd = 1, log = TRUE)
}
log_prior_sigma_x <- function(sigma) {
  stats::dexp(sigma, rate = 1, log = TRUE)
}
log_prior_sigma_y <- function(sigma) {
  stats::dexp(sigma, rate = 1, log = TRUE)
}

log_priors <- list(
  phi = log_prior_phi,
  sigma_x = log_prior_sigma_x,
  sigma_y = log_prior_sigma_y
)
```

Now we can run the PMMH algorithm using the `pmmh` function. We run 2
chains for 200 MCMC samples with a burn-in of 10. We also modify the
tuning to only use 200 pilot samples and a burn-in of 10. In practice
you would want to run it for a much larger number of samples.

``` r
library(bayesSSM)

result <- pmmh(
  y = y,
  m = 500, # number of MCMC samples
  init_fn = init_fn,
  transition_fn = transition_fn,
  log_likelihood_fn = log_likelihood_fn,
  log_priors = log_priors,
  init_params = c(phi = 0.5, sigma_x = 0.5, sigma_y = 0.5),
  burn_in = 50,
  num_chains = 2,
  seed = 1405,
  tune_control = default_tune_control(pilot_m = 200, pilot_burn_in = 10)
)
#> Running chain 1...
#> Running pilot chain for tuning...
#> Using 100 particles for PMMH:
#> Running particle MCMC chain with tuned settings...
#> Running chain 2...
#> Running pilot chain for tuning...
#> Using 100 particles for PMMH:
#> Running particle MCMC chain with tuned settings...
#> PMMH Results Summary:
#>  Parameter Mean   SD Median CI Lower.2.5% CI Upper.97.5% ESS  Rhat
#>        phi 0.78 0.11   0.78          0.53           1.00  27 1.054
#>    sigma_x 1.14 0.20   1.15          0.78           1.59   5 1.008
#>    sigma_y 0.36 0.22   0.29          0.06           0.91  17 1.143
#> Warning in pmmh(y = y, m = 500, init_fn = init_fn, transition_fn =
#> transition_fn, : Some ESS values are below 400, indicating poor mixing.
#> Consider running the chains for more iterations.
#> Warning in pmmh(y = y, m = 500, init_fn = init_fn, transition_fn = transition_fn, : 
#> Some Rhat values are above 1.01, indicating that the chains have not converged. Consider running the chains for more iterations and/or increase burn_in.
```

We get convergence warnings as expected due to the small number of
iterations.

## State-space Models

A state-space model (SSM) has the structure given in the following
directed acyclic graph (DAG):

![](man/figures/DAG_SSM.png)

The core function, `pmmh`, implements the Particle Marginal
Metropolis-Hastings, which is an algorithm that first generates a set of
$N$ particles to approximate the likelihood and then uses this
approximation in the acceptance probability. The implementation
automatically tunes the number of particles and the proposal
distribution for the parameters.
