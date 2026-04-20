## clf_generate.R
## Vendored conditional linear family utilities for correlated binary data.
## Adapted from the archived binarySimCLF package implementation of the
## Qaqish (2003) CLF algorithm; attribution retained in inst/CREDITS.md.
## Original package archive: https://github.com/cran/binarySimCLF
## Distributed here under pgeeVar's GPL (>= 2) license with retained provenance.
## No external dependencies (pure base R).

# Exchangeable correlation matrix (n x n)
clf_xch <- function(n, rho) {
  if (n <= 0) stop("n must be at least 1")
  if (n == 1) return(matrix(1))
  toeplitz(c(1, rep(rho, n - 1)))
}

# AR(1) correlation matrix (n x n)
clf_ar1 <- function(n, rho) {
  if (n <= 0) stop("n must be at least 1")
  if (n <= 2) return(clf_xch(n, rho))
  toeplitz(c(1, rho^(1:(n - 1))))
}

# Correlation to covariance matrix for binary data
# r: correlation matrix, mu: marginal mean vector
clf_cor2var <- function(r, mu) {
  p <- length(mu)
  d <- sqrt(mu * (1 - mu))
  V <- d * r * rep(d, each = p)
  rownames(V) <- NULL
  V
}

# Compute sequential regression coefficients (CLF B matrix)
# V: covariance matrix (n x n)
# Returns B where B[1:(i-1), i] = V[1:(i-1), 1:(i-1)]^{-1} V[1:(i-1), i]
clf_allReg <- function(V) {
  n <- nrow(V)
  if (n < 2) return(V)
  B <- V
  for (i in 2:n) {
    i1 <- i - 1
    B[1:i1, i] <- solve(V[1:i1, 1:i1], V[1:i1, i])
  }
  B
}

# Get bounds for conditional mean of Y_i given Y_1,...,Y_{i-1}
clf_getBnds <- function(i, u, b) {
  y_max <- (b > 0)
  nuMax <- u[i] + crossprod(y_max - u[1:(i - 1)], b)
  nuMin <- u[i] + crossprod(!y_max - u[1:(i - 1)], b)
  list(nuMin = as.numeric(nuMin), nuMax = as.numeric(nuMax))
}

# Check that all conditional means stay in [0,1] for all binary histories
# u: marginal means, B: regression coefficient matrix from clf_allReg
clf_blrchk1 <- function(u, B) {
  n <- nrow(B)
  if (n < 2) return(TRUE)
  for (i in 2:n) {
    bounds <- clf_getBnds(i, u, B[1:(i - 1), i])
    if (bounds$nuMax > 1.0 || bounds$nuMin < 0.0) return(FALSE)
  }
  TRUE
}

# Wrapper: check CLF compatibility from covariance matrix
clf_blrchk <- function(u, V) {
  clf_blrchk1(u, clf_allReg(V))
}

# Generate one correlated binary vector using CLF algorithm
# u: marginal means, B: regression coefficient matrix
clf_mbsclf1 <- function(u, B) {
  n <- length(u)
  y <- numeric(n)
  y[1] <- as.numeric(runif(1) <= u[1])
  if (n < 2) return(list(succeed = TRUE, y = y))
  r <- y - u
  for (i in 2:n) {
    i1 <- i - 1
    ci <- u[i] + sum(r[1:i1] * B[1:i1, i])
    if (ci < 0 || ci > 1) {
      return(list(succeed = FALSE, y = NULL))
    }
    y[i] <- as.numeric(runif(1) <= ci)
    r[i] <- y[i] - u[i]
  }
  list(succeed = TRUE, y = y)
}

# Generate m correlated binary vectors
# Returns m x n matrix (each row is one vector)
clf_mbsclf <- function(m, u, B, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(B)
  y <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    result <- clf_mbsclf1(u, B)
    if (!result$succeed) return(result)
    y[i, ] <- result$y
  }
  list(succeed = TRUE, y = y)
}
