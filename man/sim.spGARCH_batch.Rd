\name{sim.spGARCH_batch}
\alias{sim.spGARCH_batch}
\title{Parallel batch simulation of spatial GARCH models}
\description{
  Generates multiple independent realisations of a spatial GARCH process in parallel.
}
\usage{
sim.spGARCH_batch(n_sim = 1, ..., seed = NULL, workers = 1)
}
\arguments{
  \item{n_sim}{number of independent simulations}
  \item{...}{arguments passed to \code{\link{sim.spGARCH}}}
  \item{seed}{positive integer to initialize the random number generator (RNG)}
  \item{workers}{number of parallel workers; if 1, sequential execution is used}
}
\value{
  A list of length \code{n_sim} with the simulated spatial GARCH processes.
}
\seealso{
  \code{\link{sim.spGARCH}}
}
