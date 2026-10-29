\name{sim.stSV}
\alias{sim.stSV}
\title{Simulation of spatiotemporal stochastic volatility models}
\description{
  Generates outcomes of a spatiotemporal stochastic volatility (stSV) process on a spatial
  lattice over time. The latent log-volatility follows a spatiotemporal autoregression,
  either in the classic spatial autoregressive form or in the matrix exponential spatial
  specification (MESS).
}
\usage{
sim.stSV(n = nrow(W), TT, rho1, rho2, rho3, sigma2, W, mu = 0,
         control = list(), seed = NULL, model_type = c("classic", "MESS"))
}
\arguments{
  \item{n}{number of spatial locations}
  \item{TT}{number of time periods (after burn-in)}
  \item{rho1}{contemporaneous spatial dependence parameter}
  \item{rho2}{temporal persistence parameter}
  \item{rho3}{spatiotemporal (lagged neighbour) parameter}
  \item{sigma2}{variance of the volatility innovations}
  \item{W}{\code{n} times \code{n} spatial weight matrix}
  \item{mu}{log-volatility level, a scalar or a vector of length \code{n}}
  \item{control}{list of control arguments: \code{seed} (RNG seed) and \code{burnin} (default 50)}
  \item{seed}{RNG seed; overrides \code{control$seed}}
  \item{model_type}{\code{"classic"} or \code{"MESS"}}
}
\details{
  Observations are \eqn{Y_t = \exp(h_t / 2) \odot V_t} with \eqn{h_t = \mu + z_t},
  \eqn{V_t \sim N(0, I_n)} and \eqn{U_t \sim N(0, \sigma^2 I_n)}, starting from
  \eqn{z_0 = 0}.

  Classic model:
  \deqn{(I - \rho_1 W) z_t = (\rho_2 I + \rho_3 W) z_{t-1} + U_t.}
  This requires \eqn{|\rho_1| < 1}.

  MESS model:
  \deqn{z_t = e^{-\rho_1 W} \left( ((\rho_2 + 1) I - e^{\rho_3 W}) z_{t-1} + U_t \right).}
  To first order, \eqn{e^{-\rho_1 W} \approx I - \rho_1 W} and
  \eqn{(\rho_2 + 1) I - e^{\rho_3 W} \approx \rho_2 I - \rho_3 W}.

  In both models the spectral radius of the transition matrix must be below 1; otherwise
  the function stops with an error.
}
\value{
  An \code{n x TT} matrix of simulated outcomes.
}
\seealso{
  \code{\link{sim.stARCH}}, \code{\link{sim.stGARCH}}
}
\examples{
d <- 5
n <- d^2
coords <- expand.grid(x = 1:d, y = 1:d)
W <- outer(1:n, 1:n, function(i, j) {
  abs(coords$x[i] - coords$x[j]) + abs(coords$y[i] - coords$y[j]) == 1
}) * 1
W <- W / rowSums(W)
Y <- sim.stSV(TT = 100, rho1 = 0.3, rho2 = 0.5, rho3 = 0.1, sigma2 = 0.2,
              W = W, mu = -0.5, seed = 1)
dim(Y)
}
