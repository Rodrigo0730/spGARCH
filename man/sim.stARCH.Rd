\name{sim.stARCH}
\alias{sim.stARCH}
\title{Simulation of spatiotemporal ARCH models}
\description{
  Generates outcomes of a spatiotemporal ARCH (stARCH) process on a spatial lattice over time,
  following Otto et al. (2025). The model incorporates both spatial and temporal
  dependencies in log-volatility.
}
\usage{
sim.stARCH(n = nrow(W), TT, rho, gamma, delta, alpha, W, control = list(), seed = NULL)
}
\arguments{
  \item{n}{number of spatial locations}
  \item{TT}{number of time periods (after burn-in)}
  \item{rho}{spatial autoregressive parameter}
  \item{gamma}{temporal dependence parameter}
  \item{delta}{spatial-temporal interaction parameter}
  \item{alpha}{unconditional log-volatility level}
  \item{W}{\code{n} times \code{n} spatial weight matrix}
  \item{control}{list of control arguments: \code{seed} (RNG seed) and \code{burnin} (default 50)}
  \item{seed}{RNG seed; overrides \code{control$seed}}
}
\value{
  An \code{n x TT} matrix of simulated spatiotemporal ARCH outcomes.
}
\seealso{
  \code{\link{sim.stGARCH}}
}
