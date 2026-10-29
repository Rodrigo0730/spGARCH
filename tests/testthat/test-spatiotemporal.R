# Tests for the spatiotemporal simulators sim.stARCH, sim.stGARCH and sim.stSV.

# ---- Setup: 5 x 5 rook lattice, row-standardised ---------------------------

d <- 5
n <- d^2
coords <- expand.grid(x = 1:d, y = 1:d)
W <- outer(1:n, 1:n, function(i, j) {
  abs(coords$x[i] - coords$x[j]) + abs(coords$y[i] - coords$y[j]) == 1
}) * 1
W <- W / rowSums(W)
TT <- 30

# ---- sim.stARCH -------------------------------------------------------------

test_that("sim.stARCH returns a finite n x TT matrix", {
  Y <- sim.stARCH(TT = TT, rho = 0.3, gamma = 0.2, delta = 0.1, alpha = 0.5,
                  W = W, seed = 1)
  expect_equal(dim(Y), c(n, TT))
  expect_true(all(is.finite(Y)))
})

test_that("sim.stARCH is reproducible, with the seed as argument or in control", {
  args <- list(TT = TT, rho = 0.3, gamma = 0.2, delta = 0.1, alpha = 0.5, W = W)
  Y1 <- do.call(sim.stARCH, c(args, seed = 1))
  Y2 <- do.call(sim.stARCH, c(args, list(control = list(seed = 1))))
  expect_equal(Y1, Y2)
})

# ---- sim.stGARCH ------------------------------------------------------------

test_that("sim.stGARCH returns a finite n x TT matrix and is reproducible", {
  args <- list(TT = TT, rho = 0.3, gamma = 0.2, delta = 0.1, beta = 0.1,
               alpha = 0.5, W = W, seed = 1)
  Y <- do.call(sim.stGARCH, args)
  expect_equal(dim(Y), c(n, TT))
  expect_true(all(is.finite(Y)))
  expect_equal(Y, do.call(sim.stGARCH, args))
})

# ---- sim.stSV ---------------------------------------------------------------

# Rebuilds the recursion from the same random draws (V first, then U).
reference_stSV <- function(TT, rho1, rho2, rho3, sigma2, mu, seed, burnin, type) {
  set.seed(seed)
  total <- TT + burnin
  V <- matrix(rnorm(n * total), n, total)
  U <- matrix(rnorm(n * total, sd = sqrt(sigma2)), n, total)
  if (type == "classic") {
    S_inv <- solve(diag(n) - rho1 * W)
    Tmat <- S_inv %*% (rho2 * diag(n) + rho3 * W)
    Bmat <- S_inv
  } else {
    E1 <- as.matrix(Matrix::expm(Matrix::Matrix(-rho1 * W)))
    E3 <- as.matrix(Matrix::expm(Matrix::Matrix(rho3 * W)))
    Tmat <- E1 %*% ((rho2 + 1) * diag(n) - E3)
    Bmat <- E1
  }
  z <- rep(0, n)
  Y <- matrix(0, n, total)
  for (t in seq_len(total)) {
    z <- as.numeric(Tmat %*% z + Bmat %*% U[, t])
    Y[, t] <- exp((mu + z) / 2) * V[, t]
  }
  Y[, (burnin + 1):total]
}

for (type in c("classic", "MESS")) {
  test_that(paste0("sim.stSV (", type, ") follows the model recursion"), {
    Y <- sim.stSV(TT = TT, rho1 = 0.3, rho2 = 0.5, rho3 = 0.1, sigma2 = 0.2,
                  W = W, mu = -0.5, seed = 3, model_type = type)
    expected <- reference_stSV(TT, 0.3, 0.5, 0.1, 0.2, mu = -0.5, seed = 3,
                               burnin = 50, type = type)
    expect_equal(dim(Y), c(n, TT))
    expect_equal(Y, expected, tolerance = 1e-10)
  })
}

test_that("sim.stSV honours control$burnin", {
  Y <- sim.stSV(TT = TT, rho1 = 0.3, rho2 = 0.5, rho3 = 0, sigma2 = 0.2, W = W,
                seed = 4, control = list(burnin = 10))
  expected <- reference_stSV(TT, 0.3, 0.5, 0, 0.2, mu = 0, seed = 4,
                             burnin = 10, type = "classic")
  expect_equal(Y, expected, tolerance = 1e-10)
})

test_that("sim.stSV rejects |rho1| >= 1 in the classic model", {
  expect_error(
    sim.stSV(TT = TT, rho1 = 1, rho2 = 0.5, rho3 = 0, sigma2 = 0.2, W = W),
    "rho1"
  )
})

test_that("sim.stSV rejects non-stationary parameters", {
  expect_error(
    sim.stSV(TT = TT, rho1 = 0.5, rho2 = 0.9, rho3 = 0.3, sigma2 = 0.2, W = W),
    "Unstable"
  )
  expect_error(
    sim.stSV(TT = TT, rho1 = -2, rho2 = 0.9, rho3 = -2, sigma2 = 0.2, W = W,
             model_type = "MESS"),
    "Unstable"
  )
})
