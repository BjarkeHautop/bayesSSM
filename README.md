
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesSSM <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/BjarkeHautop/bayesSSM/graph/badge.svg)](https://app.codecov.io/gh/BjarkeHautop/bayesSSM)
[![R-CMD-check](https://github.com/BjarkeHautop/bayesSSM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BjarkeHautop/bayesSSM/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

bayesSSM is an R package offering a set of tools for performing Bayesian
inference in state-space models (SSMs). It implements the Particle
Marginal Metropolis-Hastings (PMMH) in the main function `pmmh`  
for Bayesian inference in SSMs.

## Why bayesSSM?

While there are several alternative packages available for performing
Particle MCMC bayesSSM is designed to be simple and easy to use. It was
developed as a procrastination task during my Master’s thesis about
Particle MCMC, since I was implementing everything from scratch anyway.
Everything is written in R, so performance is not the best.

## Installation

You can install the latest stable version of bayesSSM from CRAN with:

``` r
install.packages("bayesSSM")
#> Installing package into 'C:/Users/bjark/AppData/Local/Temp/RtmpgrSY3x/temp_libpath47d41614778d'
#> (as 'lib' is unspecified)
#> Warning: package 'bayesSSM' is not available for this version of R
#> 
#> A version of this package for your version of R might be available elsewhere,
#> see the ideas at
#> https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
```

or the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("BjarkeHautop/bayesSSM")
#> ℹ Loading metadata database✔ Loading metadata database ... done
#>  
#> → Will install 22 packages.
#> → Will update 1 package.
#> → Will download 21 CRAN packages (11.07 MB), cached: 1 (136.95 kB).
#> → Will download 1 package with unknown size.
#> + bayesSSM     0.4.6 → 0.4.6  [bld][cmp][dl] (GitHub: 416d836)
#> + cli                  3.6.4  [dl] (1.39 MB)
#> + digest               0.6.37 [dl] (223.14 kB)
#> + dplyr                1.1.4  [dl] (1.58 MB)
#> + fansi                1.0.6  [dl] (322.97 kB)
#> + future               1.40.0 [dl] (962.29 kB)
#> + future.apply         1.11.3 [dl] (160.32 kB)
#> + generics             0.1.3  [dl] (83.69 kB)
#> + globals              0.17.0 
#> + glue                 1.8.0  [dl] (183.78 kB)
#> + lifecycle            1.0.4  [dl] (140.93 kB)
#> + listenv              0.9.1  [dl] (109.21 kB)
#> + magrittr             2.0.3  [dl] (229.42 kB)
#> + parallelly           1.43.0 [dl] (601.63 kB)
#> + pillar               1.10.2 [dl] (672.16 kB)
#> + pkgconfig            2.0.3  [dl] (22.81 kB)
#> + R6                   2.6.1  [dl] (88.64 kB)
#> + rlang                1.1.6  [dl] (1.63 MB)
#> + tibble               3.2.1  [dl] (695.05 kB)
#> + tidyselect           1.2.1  [dl] (228.15 kB)
#> + utf8                 1.2.4  [dl] (150.80 kB)
#> + vctrs                0.6.5  [dl] (1.36 MB)
#> + withr                3.0.2  [dl] (231.37 kB)
#> ℹ Getting 21 pkgs (11.07 MB) and 1 pkg with unknown size, 1 (136.95 kB) cached
#> ✔ Cached copy of R6 2.6.1 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of cli 3.6.4 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of digest 0.6.37 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of dplyr 1.1.4 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of fansi 1.0.6 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of future.apply 1.11.3 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of future 1.40.0 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of generics 0.1.3 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of glue 1.8.0 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of lifecycle 1.0.4 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of listenv 0.9.1 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of magrittr 2.0.3 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of parallelly 1.43.0 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of pillar 1.10.2 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of pkgconfig 2.0.3 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of rlang 1.1.6 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of tibble 3.2.1 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of tidyselect 1.2.1 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of utf8 1.2.4 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of vctrs 0.6.5 (x86_64-w64-mingw32) is the latest build
#> ✔ Cached copy of withr 3.0.2 (i386+x86_64-w64-mingw32) is the latest build
#> ✔ Got bayesSSM 0.4.6 (source) (210.37 kB)
#> ✔ Installed R6 2.6.1  (588ms)
#> ✔ Installed cli 3.6.4  (848ms)
#> ✔ Installed digest 0.6.37  (1.1s)
#> ✔ Installed dplyr 1.1.4  (1.1s)
#> ✔ Installed fansi 1.0.6  (1.2s)
#> ✔ Installed generics 0.1.3  (1.1s)
#> ✔ Installed future.apply 1.11.3  (1.2s)
#> ✔ Installed listenv 0.9.1  (1.2s)
#> ✔ Installed magrittr 2.0.3  (1.3s)
#> ✔ Installed glue 1.8.0  (1.5s)
#> ✔ Installed lifecycle 1.0.4  (1.6s)
#> ✔ Installed pillar 1.10.2  (1.5s)
#> ✔ Installed pkgconfig 2.0.3  (1.5s)
#> ✔ Installed globals 0.17.0  (1.8s)
#> ✔ Installed future 1.40.0  (1.9s)
#> ✔ Installed parallelly 1.43.0  (1.7s)
#> ✔ Installed rlang 1.1.6  (869ms)
#> ✔ Installed tidyselect 1.2.1  (873ms)
#> ✔ Installed tibble 3.2.1  (936ms)
#> ✔ Installed utf8 1.2.4  (773ms)
#> ✔ Installed withr 3.0.2  (732ms)
#> ✔ Installed vctrs 0.6.5  (791ms)
#> ℹ Packaging bayesSSM 0.4.6
#> ✔ Packaged bayesSSM 0.4.6 (810ms)
#> ℹ Building bayesSSM 0.4.6
#> ✔ Built bayesSSM 0.4.6 (2.7s)
#> ✔ Installed bayesSSM 0.4.6 (github::BjarkeHautop/bayesSSM@416d836) (85ms)
#> ✔ 1 pkg + 24 deps: kept 1, upd 1, added 22, dld 1 (NA B) [15.4s]
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
        \phi &\sim \text{Uniform}(0,1), \\
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
  stats::dunif(phi, min = 0, max = 1, log = TRUE)
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
  pilot_init_params = list(
    c(phi = 0.4, sigma_x = 0.4, sigma_y = 0.4),
    c(phi = 0.8, sigma_x = 0.8, sigma_y = 0.8)
  ),
  burn_in = 50,
  num_chains = 2,
  seed = 1405,
  tune_control = default_tune_control(pilot_m = 200, pilot_burn_in = 10)
)
#> Running chain 1...
#> Running pilot chain for tuning...
#> Using 50 particles for PMMH:
#> Running particle MCMC chain with tuned settings...
#> Running chain 2...
#> Running pilot chain for tuning...
#> Using 50 particles for PMMH:
#> Running particle MCMC chain with tuned settings...
#> PMMH Results Summary:
#>  Parameter Mean   SD Median CI Lower.2.5% CI Upper.97.5% ESS  Rhat
#>        phi 0.82 0.09   0.83          0.65           0.95  65 1.009
#>    sigma_x 0.69 0.40   0.68          0.03           1.36  34 1.152
#>    sigma_y 0.82 0.37   0.86          0.21           1.49  18 1.264
#> Warning in pmmh(y = y, m = 500, init_fn = init_fn, transition_fn =
#> transition_fn, : Some ESS values are below 400, indicating poor mixing.
#> Consider running the chains for more iterations.
#> Warning in pmmh(y = y, m = 500, init_fn = init_fn, transition_fn = transition_fn, : 
#> Some Rhat values are above 1.01, indicating that the chains have not converged. 
#> Consider running the chains for more iterations and/or increase burn_in.
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
