#' Particle filter functions
#'
#' The package provides several particle filter implementations for state-space
#' models for estimating the intractable marginal likelihood
#' \eqn{p(y_{1:T}\mid \theta)}:
#' \itemize{
#'   \item \code{\link{auxiliary_filter}}
#'   \item \code{\link{bootstrap_filter}}
#'   \item \code{\link{resample_move_filter}}
#' }
#' The simplest one is the \code{\link{bootstrap_filter}}, and is thus
#' recommended as a starting point.
#'
#' @name particle_filter
NULL
