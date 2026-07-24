\name{qml.spGARCH}
\alias{qml.spGARCH}
\title{Quasi-maximum-likelihood estimation of spatial GARCH models}
\description{
  Estimates the parameters of a spatial GARCH model via quasi-maximum-likelihood (QML).
}
\usage{
qml.spGARCH(formula, W1, W2, type = "spGARCH", data = NULL, b = 2,
            start = NULL, control = list())
}
\arguments{
  \item{formula}{an object of class \code{"formula"}: a symbolic description of the model to be fitted}
  \item{W1}{\code{n} times \code{n} spatial weight matrix for the ARCH term}
  \item{W2}{\code{n} times \code{n} spatial weight matrix for the GARCH term}
  \item{type}{type of spGARCH model: \code{"spGARCH"}, \code{"log-spGARCH"}, or \code{"e-spGARCH"}}
  \item{data}{an optional data frame containing the variables in the model}
  \item{b}{parameter b for logarithmic and exponential spGARCH (default 2)}
  \item{start}{vector of starting values for the optimisation}
  \item{control}{list of control arguments for \code{\link[Rsolnp]{solnp}}}
}
\value{
  An object of class \code{"spGARCH"} with the estimation results.
}
\seealso{
  \code{\link{qml.spARCH}}, \code{\link{sim.spGARCH}}
}
