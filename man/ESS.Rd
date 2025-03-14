% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ess.R
\name{ess}
\alias{ess}
\title{Estimate effective sample size (ESS) of MCMC chains.}
\usage{
ess(chains)
}
\arguments{
\item{chains}{A matrix of dimensions m (iterations) x k (chains).}
}
\value{
The estimated effective sample size (ess) of the chains.
}
\description{
Estimate effective sample size (ESS) of MCMC chains.
}
\details{
Uses the formula for ESS proposed by Vehtari et al. (2021).
}
\examples{
chains <- matrix(rnorm(3000), nrow = 1000, ncol = 3)
ess(chains)
}
\references{
Vehtari et al. (2021). Rank-normalization, folding, and
localization: An improved R-hat for assessing convergence of MCMC.
Available at: https://doi.org/10.1214/20-BA1221
}
