# Checks of the C++ simulation backend against the model equations.

set.seed(1)

ring_W <- function(n) {
  W <- Matrix::sparseMatrix(i = c(1:n, 1:n), j = c(2:n, 1, n, 1:(n - 1)),
                            x = 1, dims = c(n, n))
  W / Matrix::rowSums(W)
}

empty_W <- function(n) {
  Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
}

test_that("spARCH backend: zero W gives h = alpha", {
  n <- 10
  alpha <- 1.5
  h <- .Simulation_CPP(rnorm(n), empty_W(n), alpha)
  expect_equal(h, rep(alpha, n), ignore_attr = TRUE)
})

test_that("spARCH backend satisfies h = alpha + rho W diag(Y) Y", {
  n <- 10
  rhoW <- as(0.3 * ring_W(n), "dgCMatrix")
  eps <- rnorm(n)
  alpha <- 1
  h <- .Simulation_CPP(eps, rhoW, alpha)
  Y <- sqrt(h) * eps
  expect_equal(as.numeric(h), alpha + as.numeric(rhoW %*% Y^2), tolerance = 1e-10)
})

test_that("log-spARCH backend satisfies ln h = alpha + rho W ln Y^2", {
  n <- 10
  rhoW <- as(0.5 * ring_W(n), "dgCMatrix")
  eps <- rnorm(n)
  alpha <- 0.2
  Y <- .Simulation_CPP_LOG(eps, rhoW, alpha, 1)
  log_h <- log(Y^2 / eps^2)
  expect_equal(as.numeric(log_h), alpha + as.numeric(rhoW %*% log(Y^2)), tolerance = 1e-10)
})

test_that("e-spARCH backend satisfies ln h = alpha + W b ln|eps|, also for tiny eps", {
  n <- 10
  W <- as(0.3 * ring_W(n), "dgCMatrix")
  eps <- runif(n, 1e-8, 1e-4)
  alpha <- 1
  b <- 2
  expected_Y <- sqrt(exp(alpha + as.numeric(W %*% (b * log(abs(eps)))))) * eps
  expect_equal(.Simulation_CPP_E(eps, W, alpha, b), expected_Y, ignore_attr = TRUE)
})

test_that("spGARCH backend: zero W2 reduces X = (I - W2)^-1 (alpha + W1 tau) to alpha + W1 tau", {
  n <- 10
  W1 <- as(0.3 * ring_W(n), "dgCMatrix")
  alpha <- 1.5
  tau <- runif(n)
  X <- .Simulation_spGARCH_CPP(alpha, tau, W1, empty_W(n))
  expect_equal(X, alpha + as.numeric(W1 %*% tau), ignore_attr = TRUE)
})
