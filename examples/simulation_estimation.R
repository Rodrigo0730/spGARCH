## simulation -> estimation framework

#libraries
library(Matrix)
devtools::load_all(".") #no need to rerun after first run

#Lower triangular W, oriented process
n <- 5000
TT <- 1000
set.seed(1)
W <- matrix(0, n, n)
for (i in 2:n) {
  n_neigh <- min(i - 1, 3)
  cols <- sample(seq_len(i - 1), n_neigh)
  W[i, cols] <- 1
}
W <- W / pmax(rowSums(W), 1)
W <- as(W, "dgCMatrix")

cat(sprintf(
  "W: %dx%d, lower tri=%s, row sums ok=%s\n",
  nrow(W),
  ncol(W),
  all(W[upper.tri(W, diag = TRUE)] == 0),
  all(abs(rowSums(W) - 1) < 1e-10 | rowSums(W) == 0)
))

cat("\n--- Simulating spatial ARCH (spARCH) ---\n")
set.seed(42)
rho <- 0.3
alpha <- 1

#tell spARCH specifically that its lower triangular not to truncate
Y <- sim.spARCH(
  n = nrow(W),
  rho,
  alpha,
  W,
  b = 2,
  type = "spARCH",
  seed = NULL,
  control = list(triangular = TRUE)
)

cat(sprintf(
  "spARCH dim: %dx%d, all finite: %s, all real: %s\n",
  nrow(Y),
  ncol(Y),
  all(is.finite(Y)),
  all(is.numeric(Y))
))

fit <- tryCatch(
  qml.spARCH(Y ~ 0, W),
  error = function(e) {
    cat("ERROR:", e$message, "\n")
    traceback(3)
    NULL
  }
)

if (!is.null(fit)) {
  cat(sprintf("\nTrue:  alpha=%.2f rho=%.2f\n", alpha, rho))
  cat(sprintf(
    "Estimated:   alpha=%.3f rho=%.3f\n",
    fit$coefficients[1],
    fit$coefficients[2]
  ))
}

rm(Y, rho, alpha)


cat("\n--- Simulating spatiotemporal ARCH (stARCH) ---\n")
set.seed(42)
rho <- 0.3
gamma <- 0.2
delta <- 0.1
alpha <- 1

Y <- sim.stARCH(
  n = n,
  TT = TT,
  rho = rho,
  gamma = gamma,
  delta = delta,
  alpha = alpha,
  W = W,
  control = list(seed = 42, burnin = 100)
)
cat(sprintf(
  "stARCH dim: %dx%d, all finite: %s, all real: %s\n",
  nrow(Y),
  ncol(Y),
  all(is.finite(Y)),
  all(is.numeric(Y))
))

# ---- GMM estimation ----
cat("\n--- GMM estimation ---\n")

fit <- tryCatch(
  gmm.stARCH(Y, W),
  error = function(e) {
    cat("ERROR:", e$message, "\n")
    traceback(3)
    NULL
  }
)

if (!is.null(fit)) {
  cat(sprintf("\nTrue:  rho=%.2f gamma=%.2f delta=%.2f\n", rho, gamma, delta))
  cat(sprintf(
    "Est:   rho=%.3f gamma=%.3f delta=%.3f\n",
    fit$coefficients[1],
    fit$coefficients[2],
    fit$coefficients[3]
  ))
}


rm(Y, rho, gamma, delta, alpha)

cat("\n--- Simulating spatial GARCH (spGARCH) ---\n")
set.seed(42)
rho <- 0.2
lambda <- 0.1
alpha <- 1

W1 <- W
W2 <- matrix(0, n, n)
for (i in 2:n) {
  n_neigh <- min(i - 1, 3)
  cols <- sample(seq_len(i - 1), n_neigh)
  W2[i, cols] <- 1
}
W2 <- W2 / pmax(rowSums(W2), 1)
W2 <- as(W2, "dgCMatrix")

Y <- sim.spGARCH(
  n = n,
  rho = rho,
  lambda = lambda,
  alpha = alpha,
  W1 = W1,
  W2 = W2,
  type = "spGARCH",
  seed = NULL,
  control = list(triangular = TRUE)
)

cat(sprintf(
  "spGARCH length: %d, all finite: %s, all real: %s\n",
  length(Y),
  all(is.finite(Y)),
  all(is.numeric(Y))
))

# fit <- tryCatch(
#   qml.spGARCH(Y ~ 0, W, W, type = "spGARCH"),
#   error = function(e) {
#     cat("ERROR:", e$message, "\n")
#     traceback(3)
#     NULL
#   }
# )

fit <- tryCatch(
  nlse.spGARCH(y = Y, W1 = W1, W2 = W2),
  error = function(e) {
    cat("ERROR:", e$message, "\n")
    NULL
  }
)

if (!is.null(fit)) {
  cat(sprintf("\nTrue:  alpha=%.2f rho=%.2f lambda=%.2f\n", alpha, rho, lambda))
  cat(sprintf(
    "Est:   alpha=%.3f rho=%.3f lambda=%.3f\n",
    fit$coefficients[1],
    fit$coefficients[2],
    fit$coefficients[3]
  ))
}

rm(Y, rho, lambda, alpha)

cat("\nDone.\n")
