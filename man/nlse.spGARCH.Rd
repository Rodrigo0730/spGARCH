\name{nlse.spGARCH}
\alias{nlse.spGARCH}
\title{Non-linear least squares estimation of spGARCH models}
\description{
  Estimates the parameters of an additive spatial GARCH model using non-linear least squares (NLSE).
  The additive form models log-volatility as a linear spatial process.
}
\usage{
nlse.spGARCH(y, W1, W2, start = c(1, 0.5, 0.5), control = list())
}
\arguments{
  \item{y}{numeric vector of observations}
  \item{W1}{\code{n} times \code{n} spatial weight matrix for the ARCH term}
  \item{W2}{\code{n} times \code{n} spatial weight matrix for the GARCH term}
  \item{start}{vector of starting values: \code{c(alpha, rho, lambda)}}
  \item{control}{list of control arguments for \code{\link[Rsolnp]{solnp}}}
}
\value{
  A list with components \code{coefficients}, \code{residuals}, \code{h} (fitted volatilities),
  \code{SSE} (sum of squared errors), and \code{opt_details}.
}
\seealso{
  \code{\link{qml.spGARCH}}
}
