\name{sim.spARCH_batch}
\alias{sim.spARCH_batch}
\title{Parallel batch simulation of spatial ARCH models}
\description{
  Generates multiple independent realisations of a spatial ARCH process in parallel.
}
\usage{
sim.spARCH_batch(n_sim = 1, ..., seed = NULL, workers = 1)
}
\arguments{
  \item{n_sim}{number of independent simulations}
  \item{...}{arguments passed to \code{\link{sim.spARCH}}}
  \item{seed}{positive integer to initialize the random number generator (RNG)}
  \item{workers}{number of parallel workers; if 1, sequential execution is used}
}
\value{
  A list of length \code{n_sim} with the simulated spatial ARCH processes.
}
\seealso{
  \code{\link{sim.spARCH}}
}
