##################################################
#
# Simulation of spatial ARCH and GARCH models
#
##################################################

# spatial ARCH process

sim.spARCH <- function(
  n = nrow(W),
  rho,
  alpha,
  W,
  b = 2,
  type = "spARCH",
  seed = NULL,
  control = list()
) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)
  }

  con <- list(seed = sample(1:1e9, 1), silent = FALSE, triangular = FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) {
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  }
  if (!is.null(seed)) {
    con$seed <- seed
  }

  if (!is.element(type, c("spARCH", "log-spARCH", "e-spARCH", "complex-spARCH"))) {
    stop("invalid type: ", type)
  }

  R.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))

  if (length(n) > 1) {
    n <- length(n)
  }
  if (!is.numeric(n)) {
    stop("non-numeric argument n. \n")
  }
  if (n < 1) {
    stop("n must be greater than 1 \n")
  }
  if (dim(W)[1] != dim(W)[2]) {
    stop("Spatial weight matrix W must be quadratic \n")
  }
  if (n != dim(W)[1]) {
    stop("Spatial weight matrix W must have dimension n x n \n")
  }
  if (!inherits(W, "dgCMatrix")) {
    W <- methods::as(W, "dgCMatrix")
  }
  W <- rho * W

  set.seed(con$seed)

  if (type == "spARCH") {
    if (!con$triangular) {
      W_norm <- norm(W, "1")
      if (W_norm > 1e-12) {
        a <- 1 / sqrt(W_norm)
        eps <- truncnorm::rtruncnorm(n, a = -a, b = a)
      } else {
        eps <- rnorm(n)
      }
    } else {
      eps <- rnorm(n)
    }
  } else {
    eps <- rnorm(n)
  }

  if (!con$silent) {
    cat(
      "Current random seed: ",
      con$seed,
      " (To reproduce, use argument 'seed = ",
      con$seed,
      "') \n",
      sep = ""
    )
  }
  if (type == "log-spARCH") {
    return(.Simulation_CPP_LOG(eps, W, alpha, b))
  } else if (type == "e-spARCH") {
    return(.Simulation_CPP_E(eps, W, alpha, b))
  }
  h <- .Simulation_CPP(eps, W, alpha)
  if (type == "spARCH") {
    return(as.numeric(sqrt(h) * eps))
  }
  if (type == "complex-spARCH") {
    return(as.complex(sqrt(as.complex(h)) * eps))
  }
}

sim.spARCH_batch <- function(n_sim = 1, ..., seed = NULL, workers = 1) {
  if (workers > 1) {
    cl <- parallel::makeCluster(workers)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    foreach::registerDoSEQ()
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  `%dorng%` <- doRNG::`%dorng%`
  i <- NULL

  results <- foreach::foreach(i = seq_len(n_sim), .packages = "spGARCH") %dorng%
    {
      sim.spARCH(..., seed = NULL)
    }

  return(results)
}


sim.spGARCH <- function(
  n = nrow(W1),
  rho,
  lambda,
  alpha,
  W1,
  W2,
  b = 2,
  zeta = 0.5,
  theta = 0.5,
  type = "spGARCH",
  seed = NULL,
  control = list()
) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)
  }

  con <- list(seed = sample(1:1e9, 1), silent = FALSE, triangular = FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) {
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  }
  if (!is.null(seed)) {
    con$seed <- seed
  }

  R.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))

  if (length(n) > 1) {
    n <- length(n)
  }
  if (!is.numeric(n)) {
    stop("non-numeric argument n. \n")
  }
  if (
    !is.element(
      type,
      c("spGARCH", "e-spGARCH", "log-spGARCH", "h-spGARCH", "complex-spGARCH")
    )
  ) {
    stop("invalid type: ", type)
  }
  if (n < 1) {
    stop("n must be greater than 1 \n")
  }
  if (dim(W1)[1] != dim(W1)[2]) {
    stop("Spatial weight matrix W1 must be quadratic \n")
  }
  if (n != dim(W1)[1]) {
    stop("Spatial weight matrix W1 must have dimension n x n \n")
  }
  if (dim(W2)[1] != dim(W2)[2]) {
    stop("Spatial weight matrix W2 must be quadratic \n")
  }
  if (n != dim(W2)[1]) {
    stop("Spatial weight matrix W2 must have dimension n x n \n")
  }
  if (dim(W2)[1] != dim(W1)[1] | dim(W2)[2] != dim(W1)[2]) {
    stop("Spatial weight matrices W1 and W2 must have equal dimensions \n")
  }
  if (!inherits(W1, "dgCMatrix")) {
    W1 <- methods::as(W1, "dgCMatrix")
  }
  W1 <- rho * W1
  if (!inherits(W2, "dgCMatrix")) {
    W2 <- methods::as(W2, "dgCMatrix")
  }
  W2 <- lambda * W2

  set.seed(con$seed)

  if (length(alpha) == 1) {
    alpha <- as.numeric(alpha)
  } else if (length(alpha) == n) {
    stop(
      "spGARCH implementation expects scalar 'alpha' for the simulation interface; pass a scalar or adapt C++ API."
    )
  } else {
    stop("'alpha' must be a scalar or length n.")
  }

  functions <- .choose_functions(type)
  f_inv <- functions$f_inv
  tau_eps <- functions$tau_eps

  # simulation

  if (type %in% c("spGARCH", "complex-spGARCH")) {
    if (!con$triangular) {
      W1_norm <- norm(W1, "1")
      if (W1_norm > 1e-12) {
        a <- 1 / sqrt(sqrt(W1_norm))
        eps <- truncnorm::rtruncnorm(n, a = -a, b = a)
      } else {
        eps <- rnorm(n)
      }
    } else {
      eps <- rnorm(n)
    }
  } else {
    eps <- rnorm(n)
  }

  if (type %in% c("spGARCH", "complex-spGARCH")) {
    E <- Matrix::Diagonal(x = eps^2)
    M <- Matrix::Diagonal(n) - W1 %*% E - W2
    X <- as.numeric(solve(M, rep(alpha, n)))
  } else {
    tau <- tau_eps(eps, W1, W2, alpha, theta, zeta, b)
    X <- .Simulation_spGARCH_CPP(alpha, as.numeric(tau), W1, W2)
  }

  if (type == "complex-spGARCH") {
    Y <- sqrt(as.complex(as.vector(f_inv(X)))) * eps
  } else {
    Y <- sqrt(as.vector(f_inv(X))) * eps
  }

  return(as.vector(Y))
}

sim.spGARCH_batch <- function(n_sim = 1, ..., seed = NULL, workers = 1) {
  if (workers > 1) {
    cl <- parallel::makeCluster(workers)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    foreach::registerDoSEQ()
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  `%dorng%` <- doRNG::`%dorng%`
  i <- NULL

  results <- foreach::foreach(i = seq_len(n_sim), .packages = "spGARCH") %dorng%
    {
      sim.spGARCH(..., seed = NULL)
    }

  return(results)
}


##################################################
#
# Choosing Functions for Unified Approach
#
##################################################

.choose_functions <- function(model) {
  if (model == "spGARCH" | model == "complex-spGARCH") {
    f_inv <- function(x) {
      return(x)
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b) {
      # function tau in paper
      eps_2 <- eps^2
      n <- length(eps)
      return(
        diag(eps_2) %*%
          Matrix::solve(
            Matrix::Diagonal(n) - Matrix(W_1) %*% diag(eps_2) - Matrix(W_2)
          ) %*%
          alpha
      )
    }
    tau_y <- function(y) {
      return(y^2)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv) {
      # mapping from y to eps
      n <- length(y)
      X <- solve(diag(n) - Matrix(lambdaW_2)) %*%
        (alpha + Matrix(rhoW_1) %*% tau_y(y))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    # sourceCpp("Cpp_Functions.cpp")

    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) {
      n <- length(eps)
      aux <- solve(diag(n) - rhoW_1 %*% diag(as.vector(eps^2)) - lambdaW_2)
      result <- array(, dim = c(n, n))
      for (j in 1:n) {
        aux_rhoW_1 <- rhoW_1
        aux_rhoW_1[, -j] <- array(0, dim = c(n, n - 1))
        eps_j <- eps[j]
        result[, j] <- 2 * eps_j * aux %*% aux_rhoW_1 %*% aux %*% alpha
      }
      return(result)
    }
  } else if (model == "log-spGARCH") {
    f_inv <- function(x) {
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b) {
      return(log(eps^2))
    }
    tau_y <- function(y) {
      return(NULL)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv) {
      n <- length(y)
      X <- as.numeric(solve(
        Matrix::Diagonal(n) - lambdaW_2,
        alpha + as.numeric(rhoW_1 %*% log(y^2))
      ))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, W_1, W_2, rho, lambda) {
      tau_eps_prime <- function(eps) {
        return(2 / eps)
      }
      n <- length(eps)
      rhoW_1 <- rho * W_1
      lambdaW_2 <- lambda * W_2
      aux <- solve(Matrix::Diagonal(n) - rhoW_1 - lambdaW_2, rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i <- array(rep(h, n), dim = c(n, n))
      return(aux * tau_eps_prime(eps_j) * h_i)
    }
  } else if (model == "e-spGARCH") {
    f_inv <- function(x) {
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b) {
      return(theta * eps)
    }
    tau_y <- function(y) {
      return(log(abs(y)))
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv) {
      n <- length(y)
      M <- Matrix::Diagonal(n) - lambdaW_2 + (b / 2) * rhoW_1
      rhs <- rep(alpha, n) + as.numeric(rhoW_1 %*% (b * tau_y(y)))
      X <- as.numeric(solve(M, rhs))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, W_1, W_2, rho, lambda) {
      n <- length(eps)
      rhoW_1 <- rho * W_1
      lambdaW_2 <- lambda * W_2
      M <- Matrix::Diagonal(n) - lambdaW_2
      aux <- as.matrix(solve(M, rhoW_1))
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i <- array(rep(h, n), dim = c(n, n))
      return(b * aux / eps_j * h_i)
    }
  } else if (model == "h-spGARCH") {
    f_inv <- function(x) {
      return(exp(x))
    }
    tau_eps <- function(eps, W_1, W_2, alpha, theta, zeta, b) {
      # function tau in paper
      eps_2 <- log(eps^2)
      n <- length(eps)
      return(
        solve(diag(n) - Matrix(W_1) - Matrix(W_2)) %*%
          alpha +
          (diag(n) +
            solve(diag(n) - Matrix(W_1) - Matrix(W_2)) %*% Matrix(W_1)) %*%
            eps_2
      )
    }
    tau_y <- function(y) {
      log(y^2)
    }
    g <- function(y, alpha, rhoW_1, lambdaW_2, theta, zeta, b, tau_y, f_inv) {
      # mapping from y to eps
      n <- length(y)
      X <- solve(diag(n) - Matrix(lambdaW_2)) %*% (alpha + rhoW_1 %*% tau_y(y))
      eps <- y / sqrt(as.vector(f_inv(X)))
      return(list(eps = eps, h = f_inv(X)))
    }
    d_h_d_eps <- function(eps, h, alpha, theta, zeta, b, rhoW_1, lambdaW_2) {
      n <- length(eps)
      aux <- solve(diag(n) - Matrix(rhoW_1) - Matrix(lambdaW_2)) %*%
        Matrix(rhoW_1)
      eps_j <- t(array(rep(eps, n), dim = c(n, n)))
      h_i <- array(rep(h, n), dim = c(n, n))
      return(2 * aux / eps_j * h_i)
    }
  }

  return(list(
    f_inv = f_inv,
    tau_eps = tau_eps,
    tau_y = tau_y,
    g = g,
    d_h_d_eps = d_h_d_eps
  ))
}
