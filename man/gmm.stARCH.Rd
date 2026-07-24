\name{gmm.stARCH}
\alias{gmm.stARCH}
\title{GMM estimation of spatiotemporal ARCH models}
\description{
  Estimates the parameters of a spatiotemporal ARCH (stARCH) model using the
  Generalised Method of Moments (GMM). The estimation uses Helmert-transformed
  data and spatial instrumental variables.
}
\usage{
gmm.stARCH(Y, W, p = 1, control = list())
}
\arguments{
  \item{Y}{an \code{n x T} matrix of spatiotemporal observations}
  \item{W}{\code{n} times \code{n} spatial weight matrix}
  \item{p}{autoregressive order (currently fixed to 1)}
  \item{control}{list of control arguments}
}
\value{
  A list with components \code{coefficients} (estimated \code{rho}, \code{gamma}, \code{delta}),
  \code{Y_double_star} (transformed data), \code{P_1} (quadratic moment matrix),
  \code{Q_N} (instrument matrix), and \code{opt_res} (optimisation details).
}
\seealso{
  \code{\link{sim.stARCH}}, \code{\link{sim.stGARCH}}
}
