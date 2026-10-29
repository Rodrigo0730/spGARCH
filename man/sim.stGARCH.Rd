\name{sim.stGARCH}
\alias{sim.stGARCH}
\title{Simulation of spatiotemporal GARCH models}
\description{
  Generates outcomes of a spatiotemporal GARCH (stGARCH) process on a spatial lattice over time,
  following Otto et al. (2025). Extends the stARCH model with a temporal GARCH
  component (beta) in the log-volatility equation.
}
\usage{
sim.stGARCH(n = nrow(W), TT, rho, gamma, delta, beta, alpha, W, control = list(), seed = NULL)
}
\arguments{
  \item{n}{number of spatial locations}
  \item{TT}{number of time periods (after burn-in)}
  \item{rho}{spatial autoregressive parameter}
  \item{gamma}{temporal dependence parameter}
  \item{delta}{spatial-temporal interaction parameter}
  \item{beta}{temporal GARCH parameter}
  \item{alpha}{unconditional log-volatility level}
  \item{W}{\code{n} times \code{n} spatial weight matrix}
  \item{control}{list of control arguments: \code{seed} (RNG seed) and \code{burnin} (default 50)}
  \item{seed}{RNG seed; overrides \code{control$seed}}
}
\value{
  An \code{n x TT} matrix of simulated spatiotemporal GARCH outcomes.
}
\seealso{
  \code{\link{sim.stARCH}}
}
