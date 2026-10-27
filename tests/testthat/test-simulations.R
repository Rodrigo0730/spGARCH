# Interface tests for sim.spARCH and sim.spGARCH: output type/length for every
# model type, reproducibility, and input validation.

# ---- Setup: 5-location path graph, row-standardised -----------------------

n <- 5
W <- matrix(0, n, n)
W[1, 2] <- 1
W[2, c(1, 3)] <- 1
W[3, c(2, 4)] <- 1
W[4, c(3, 5)] <- 1
W[5, 4] <- 1
W <- W / rowSums(W)
W_sparse <- as(W, "dgCMatrix")

W1 <- W_sparse
W2 <- W_sparse

control_args <- list(seed = 123, silent = TRUE, triangular = FALSE)

sim_arch <- function(type, ...) {
  sim.spARCH(n = n, rho = 0.5, alpha = 0.2, W = W_sparse, type = type,
             control = control_args, ...)
}

sim_garch <- function(type, ...) {
  sim.spGARCH(n = n, rho = 0.5, lambda = 0.3, alpha = 0.2, W1 = W1, W2 = W2,
              type = type, control = control_args, ...)
}

# ---- sim.spARCH: model types ----------------------------------------------

for (type in c("spARCH", "log-spARCH", "e-spARCH")) {
  test_that(paste0("sim.spARCH (", type, ") returns finite numeric vector of length n"), {
    res <- sim_arch(type)
    expect_type(res, "double")
    expect_length(res, n)
    expect_true(all(is.finite(res)))
  })
}

test_that("sim.spARCH (complex-spARCH) returns complex vector of length n", {
  res <- sim_arch("complex-spARCH")
  expect_type(res, "complex")
  expect_length(res, n)
  expect_true(all(is.finite(Re(res))))
})

test_that("sim.spARCH gives reproducible results with a fixed seed", {
  expect_equal(sim_arch("spARCH"), sim_arch("spARCH"))
})

test_that("sim.spARCH with triangular = TRUE works on a triangular W", {
  W_tri <- Matrix::sparseMatrix(i = 2:n, j = 1:(n - 1), x = 1, dims = c(n, n))
  expect_silent(
    res <- sim.spARCH(n = n, rho = 0.5, alpha = 1, W = W_tri,
                      control = list(seed = 123, silent = TRUE, triangular = TRUE))
  )
  expect_length(res, n)
  expect_true(all(is.finite(res)))
})

# ---- sim.spARCH: input validation -----------------------------------------

test_that("sim.spARCH rejects fewer than two locations", {
  expect_error(
    sim.spARCH(n = 1, rho = 0.5, alpha = 0.2, W = W_sparse[1, 1, drop = FALSE]),
    "n must be greater than 1"
  )
})

test_that("sim.spARCH errors for non-numeric n", {
  expect_error(
    sim.spARCH(n = "5", rho = 0.5, alpha = 0.2, W = W_sparse),
    "non-numeric argument n"
  )
})

test_that("sim.spARCH errors for non-quadratic W", {
  expect_error(
    sim.spARCH(n = n, rho = 0.5, alpha = 0.2, W = matrix(1, n, n + 1)),
    "must be quadratic"
  )
})

test_that("sim.spARCH errors for dimension mismatch between n and W", {
  W_wrong <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1, dims = c(3, 3))
  expect_error(
    sim.spARCH(n = n, rho = 0.5, alpha = 0.2, W = W_wrong),
    "must have dimension n x n"
  )
})

test_that("sim.spARCH errors for invalid type", {
  expect_error(sim_arch("invalid"), "invalid type:")
})

test_that("sim.spARCH rejects non-scalar alpha", {
  expect_error(
    sim.spARCH(n = n, rho = 0.5, alpha = rep(0.2, n), W = W_sparse,
               control = control_args),
    "scalar 'alpha'"
  )
})

test_that("sim.spARCH warns on unknown control names", {
  expect_warning(
    sim.spARCH(n = n, rho = 0.5, alpha = 0.2, W = W_sparse,
               control = list(silent = TRUE, foo = 1)),
    "unknown names in control: foo"
  )
})

test_that("sim.spARCH errors on NA in W", {
  W_na <- as.matrix(W_sparse)
  W_na[1, 2] <- NA_real_
  expect_error(sim.spARCH(n = n, rho = 0.5, alpha = 0.2, W = W_na))
})

test_that("sim.spARCH errors on NULL W", {
  expect_error(sim.spARCH(n = n, rho = 0.5, alpha = 0.2, W = NULL))
})

# ---- sim.spGARCH: model types ---------------------------------------------

for (type in c("spGARCH", "log-spGARCH", "e-spGARCH", "h-spGARCH")) {
  test_that(paste0("sim.spGARCH (", type, ") returns finite numeric vector of length n"), {
    res <- sim_garch(type)
    expect_type(res, "double")
    expect_length(res, n)
    expect_true(all(is.finite(res)))
  })
}

test_that("sim.spGARCH (complex-spGARCH) returns complex vector of length n", {
  res <- sim_garch("complex-spGARCH")
  expect_type(res, "complex")
  expect_length(res, n)
  expect_true(all(is.finite(Re(res))))
})

test_that("sim.spGARCH gives reproducible results with a fixed seed", {
  expect_equal(sim_garch("spGARCH"), sim_garch("spGARCH"))
})

# ---- sim.spGARCH: input validation ----------------------------------------

test_that("sim.spGARCH rejects fewer than two locations", {
  W_1 <- W_sparse[1, 1, drop = FALSE]
  expect_error(
    sim.spGARCH(n = 1, rho = 0.5, lambda = 0.3, alpha = 0.2, W1 = W_1, W2 = W_1),
    "n must be greater than 1"
  )
})

test_that("sim.spGARCH errors for non-numeric n", {
  expect_error(
    sim.spGARCH(n = "abc", rho = 0.5, lambda = 0.3, alpha = 0.2, W1 = W1, W2 = W2),
    "non-numeric argument n"
  )
})

test_that("sim.spGARCH errors for non-quadratic W1 or W2", {
  W_rect <- matrix(1, n, n + 1)
  expect_error(
    sim.spGARCH(n = n, rho = 0.5, lambda = 0.3, alpha = 0.2, W1 = W_rect, W2 = W2),
    "must be quadratic"
  )
  expect_error(
    sim.spGARCH(n = n, rho = 0.5, lambda = 0.3, alpha = 0.2, W1 = W1, W2 = W_rect),
    "must be quadratic"
  )
})

test_that("sim.spGARCH errors for dimension mismatch between n and W1 or W2", {
  W_wrong <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1, dims = c(3, 3))
  expect_error(
    sim.spGARCH(n = n, rho = 0.5, lambda = 0.3, alpha = 0.2, W1 = W_wrong, W2 = W2),
    "must have dimension n x n"
  )
  expect_error(
    sim.spGARCH(n = n, rho = 0.5, lambda = 0.3, alpha = 0.2, W1 = W1, W2 = W_wrong),
    "must have dimension n x n"
  )
})

test_that("sim.spGARCH errors for invalid type", {
  expect_error(sim_garch("spARCH"), "invalid type:")
})

test_that("sim.spGARCH rejects non-scalar alpha", {
  expect_error(
    sim.spGARCH(n = n, rho = 0.5, lambda = 0.3, alpha = rep(0.2, n),
                W1 = W1, W2 = W2, control = control_args),
    "scalar 'alpha'"
  )
})

test_that("sim.spGARCH warns on unknown control names", {
  expect_warning(
    sim.spGARCH(n = n, rho = 0.5, lambda = 0.3, alpha = 0.2, W1 = W1, W2 = W2,
                control = list(silent = TRUE, unknown = 1)),
    "unknown names in control: unknown"
  )
})

test_that("sim.spGARCH errors on NA in W1", {
  W_na <- as.matrix(W1)
  W_na[1, 2] <- NA_real_
  expect_error(
    sim.spGARCH(n = n, rho = 0.5, lambda = 0.3, alpha = 0.2, W1 = W_na, W2 = W2)
  )
})
