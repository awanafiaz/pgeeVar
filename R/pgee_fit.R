## pgee_fit.R
## Penalized GEE (PGEE) fitting via Firth-type bias reduction.
## Cleaned from thesis All_functions.R (geefirth, lines 460-550).
## That manuscript/thesis implementation descends from the public GPL-2
## geefirthr code associated with Momenul Haque Mondol and M. Shafiqur Rahman,
## then was modified further for the pgee-variance project.
## Changes from thesis:
##   - corpcor::pseudoinverse -> MASS::ginv
##   - Returns structured output with all intermediates for variance estimation
##   - Max iteration limit and convergence flag
##   - Supports both exchangeable and AR(1) working correlation
##   - Uses logistf for initialization when available, otherwise falls back
##     to glm / zero starts so the package can install without logistf
## Dependencies: MASS, clf_generate.R (for clf_xch, clf_ar1)

#' Fit a penalized generalized estimating equation model
#'
#' `pgee_fit()` fits a penalized generalized estimating equation (PGEE) model
#' for clustered binary outcomes using the low-level fitting routine currently
#' shipped in `pgeeVar`. The function is designed to retain the matrices and
#' residual objects needed by [compute_pgee_variances()], so its return value is
#' deliberately more detailed than a typical user-facing model fit.
#'
#' This is the low-level fitting engine underlying [pgee()]. Most end users
#' should call [pgee()] instead.
#'
#' @param y Binary response vector.
#' @param x Covariate matrix or data frame **without** an intercept column.
#' @param id Cluster identifier vector. Every cluster must contain at least two
#'   observations after missing values are removed.
#' @param corr_type Working correlation structure. Must be one of
#'   `"exchangeable"`, `"ar1"`, or `"independence"`.
#' @param tol Convergence tolerance for the maximum absolute scoring step.
#' @param max_iter Maximum number of scoring iterations.
#'
#' @details
#' Rows with missing values in `y`, `x`, or `id` are removed before fitting.
#'
#' When the optional `logistf` package is available, `pgee_fit()` uses Firth
#' logistic regression for starting values. Otherwise it falls back to
#' `glm(..., family = binomial())`, and finally to zero starts if needed.
#'
#' The returned object is intentionally verbose because the package's variance
#' estimators operate on retained cluster-level matrices such as `H_list`,
#' `V_list`, and `DT_list`.
#'
#' @return A list with fitted values and internal matrices. The most useful
#'   components for end users are:
#'   * `beta`: estimated regression coefficients, including the intercept.
#'   * `converged`: logical convergence flag.
#'   * `iter`: number of scoring iterations used.
#'   * `phi`: estimated dispersion parameter.
#'   * `alpha_hat`: estimated working-correlation parameter. For
#'     `corr_type = "independence"`, this value is always returned as `0`.
#'   * `N`, `p`, `n_i`, `n_star`: design summaries for the fitted data.
#'
#'   The object also retains advanced components used by
#'   [compute_pgee_variances()], including `X_list`, `y_list`, `mu_list`,
#'   `W_list`, `V_list`, `R_list`, `DT_list`, `r_list`, `H_list`, `I0`, and
#'   `Delta`.
#'
#' @references
#' Firth D (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika*, 80(1), 27-38.
#'
#' Mondol MH and Rahman MS (2019). Bias-reduced and separation-proof GEE with
#' small or sparse longitudinal binary data. *Statistics in Medicine*, 38(14),
#' 2544-2560.
#'
#' @seealso [pgee()], [dgp_generate()], [compute_pgee_variances()],
#'   [available_estimators()]
#'
#' @examples
#' set.seed(1)
#' sim <- dgp_generate(
#'   N = 40,
#'   gamma = 0.5,
#'   n_i = 4,
#'   beta = c(-1.2, 0.6, 0.3),
#'   rho = 0.2
#' )
#'
#' fit <- pgee_fit(
#'   y = sim$data$y,
#'   x = sim$data[c("X1", "obstime")],
#'   id = sim$data$id,
#'   corr_type = sim$corr_type
#' )
#'
#' fit$converged
#' round(fit$beta, 3)
#' fit$alpha_hat
#'
#' @export
pgee_fit <- function(y, x, id, corr_type = c("exchangeable", "ar1", "independence"),
                     tol = 1e-4, max_iter = 100) {
  corr_type <- match.arg(corr_type)

  # Remove clusters with any NA
  dat <- data.frame(y = y, x, id = id)
  dat <- dat[complete.cases(dat), ]

  y_full <- dat$y
  x_cols <- setdiff(names(dat), c("y", "id"))
  x_full <- dat[, x_cols, drop = FALSE]
  id_full <- dat$id

  # Check for single-observation clusters
  tab <- table(id_full)
  singles <- names(tab)[tab == 1]
  if (length(singles) > 0) {
    stop("Clusters with single observation: ", paste(singles, collapse = ", "))
  }

  # Split by cluster
  X_list <- split(data.frame(intercept = 1, x_full), id_full)
  y_list <- split(y_full, id_full)
  N <- length(y_list)  # number of clusters
  p <- ncol(X_list[[1]])  # number of parameters

  # Initial values from Firth logistic regression when available; otherwise
  # fall back to glm or zeros so the package remains installable without logistf.
  beta <- init_pgee_beta(y_full = y_full, x_full = x_full, p = p)

  converged <- FALSE
  iter <- 0

  while (iter < max_iter) {
    iter <- iter + 1

    # Compute mu, W for each cluster
    mu_list <- lapply(X_list, function(xi) {
      eta <- as.matrix(xi) %*% beta
      as.vector(1 / (1 + exp(-eta)))
    })
    var_mu_list <- lapply(mu_list, function(mu) mu * (1 - mu))
    W_list <- lapply(var_mu_list, function(vm) diag(as.vector(vm), nrow = length(vm)))

    # Standardized residuals for correlation estimation
    e_list <- mapply(function(y, mu, vm) (y - mu) / sqrt(vm),
                     y_list, mu_list, var_mu_list, SIMPLIFY = FALSE)

    # Estimate working correlation
    R_list <- estimate_corr(e_list, corr_type)

    # Dispersion parameter
    all_e2 <- unlist(lapply(e_list, function(e) e^2))
    n_star_iter <- sum(sapply(e_list, length))
    phi <- sum(all_e2) / (n_star_iter - p)

    # Working covariance V_i = phi * W_i^{1/2} R_i W_i^{1/2}
    V_list <- compute_V(W_list, R_list, phi)

    # D_i^T = X_i^T W_i (p x n_i)
    DT_list <- mapply(function(xi, Wi) t(as.matrix(xi)) %*% Wi,
                      X_list, W_list, SIMPLIFY = FALSE)

    # Expected Fisher information for Firth penalty
    # I(beta, alpha) = (1/phi) sum X_i^T W_i^{1/2} R_i^{-1} W_i^{1/2} X_i
    I_info <- matrix(0, p, p)
    for (i in 1:N) {
      xi <- as.matrix(X_list[[i]])
      W12 <- sqrt(W_list[[i]])
      Ri_inv <- MASS::ginv(R_list[[i]])
      I_info <- I_info + t(xi) %*% W12 %*% Ri_inv %*% W12 %*% xi
    }
    I_info <- I_info / phi

    # Compute Q_i = diag(0.5 - mu_i) for Firth penalty
    Q_list <- lapply(mu_list, function(mu) diag(0.5 - mu, nrow = length(mu)))

    # Compute Z_i: list of p diagonal matrices, one per covariate column
    # Z_i[[r]] = diag(X_i[, r])
    Z_list <- lapply(X_list, function(xi) {
      xi <- as.matrix(xi)
      lapply(1:ncol(xi), function(r) diag(xi[, r], nrow = nrow(xi)))
    })

    # Expected Fisher information derivatives: FF[[r]] = (2/phi) sum X_i^T W^{1/2} R^{-1} W^{1/2} Q_i Z_i[[r]] X_i
    FF <- vector("list", p)
    for (r in 1:p) {
      FF[[r]] <- matrix(0, p, p)
      for (i in 1:N) {
        xi <- as.matrix(X_list[[i]])
        W12 <- sqrt(W_list[[i]])
        Ri_inv <- MASS::ginv(R_list[[i]])
        FF[[r]] <- FF[[r]] + t(xi) %*% W12 %*% Ri_inv %*% W12 %*% Q_list[[i]] %*% Z_list[[i]][[r]] %*% xi
      }
      FF[[r]] <- 2 * FF[[r]] / phi
    }

    # Firth penalty: b_r = trace(I^{-1} FF[[r]])
    # The 0.5 is applied once in the Ustar formula below (matching thesis Pi.f + Ustari.f)
    I_inv <- MASS::ginv(I_info)
    penalty <- numeric(p)
    for (r in 1:p) {
      penalty[r] <- sum(diag(I_inv %*% FF[[r]]))
    }

    # Penalized score: U* = sum D_i^T V_i^{-1} (y_i - mu_i) + 0.5 * penalty
    U <- rep(0, p)
    for (i in 1:N) {
      Vi_inv <- MASS::ginv(V_list[[i]])
      U <- U + DT_list[[i]] %*% Vi_inv %*% (y_list[[i]] - mu_list[[i]])
    }
    Ustar <- U + 0.5 * penalty

    # Method of scoring update: beta <- beta + I^{-1} U*
    step <- as.vector(MASS::ginv(I_info) %*% Ustar)
    del <- max(abs(step))

    beta <- beta + step

    if (is.finite(del) && del < tol) {
      converged <- TRUE
      break
    }
  }

  # Final quantities at converged beta (recompute)
  mu_list <- lapply(X_list, function(xi) {
    eta <- as.matrix(xi) %*% beta
    as.vector(1 / (1 + exp(-eta)))
  })
  var_mu_list <- lapply(mu_list, function(mu) mu * (1 - mu))
  W_list <- lapply(var_mu_list, function(vm) diag(as.vector(vm), nrow = length(vm)))
  e_list <- mapply(function(y, mu, vm) (y - mu) / sqrt(vm),
                   y_list, mu_list, var_mu_list, SIMPLIFY = FALSE)
  R_list <- estimate_corr(e_list, corr_type)
  all_e2 <- unlist(lapply(e_list, function(e) e^2))
  n_star_final <- sum(sapply(e_list, length))
  phi <- sum(all_e2) / (n_star_final - p)
  V_list <- compute_V(W_list, R_list, phi)
  DT_list <- mapply(function(xi, Wi) t(as.matrix(xi)) %*% Wi,
                    X_list, W_list, SIMPLIFY = FALSE)

  # Recompute I0 and Delta at final beta
  I0 <- matrix(0, p, p)
  for (i in 1:N) {
    Vi_inv <- MASS::ginv(V_list[[i]])
    I0 <- I0 + DT_list[[i]] %*% Vi_inv %*% t(DT_list[[i]])
  }
  Delta <- MASS::ginv(I0)

  # Residuals
  r_list <- mapply(function(y, mu) y - mu, y_list, mu_list, SIMPLIFY = FALSE)

  # Hat matrices H_{ii} = D_i Delta D_i^T V_i^{-1}
  H_list <- vector("list", N)
  for (i in 1:N) {
    Vi_inv <- MASS::ginv(V_list[[i]])
    H_list[[i]] <- t(DT_list[[i]]) %*% Delta %*% DT_list[[i]] %*% Vi_inv
  }

  # Cluster sizes
  n_i <- sapply(y_list, length)
  n_star <- sum(n_i)

  # Estimated correlation parameter (scalar)
  alpha_hat <- R_list[[1]][1, 2]
  if (is.na(alpha_hat)) alpha_hat <- 0

  list(
    beta      = beta,
    converged = converged,
    iter      = iter,
    N         = N,
    p         = p,
    n_i       = n_i,
    n_star    = n_star,
    phi       = phi,
    alpha_hat = alpha_hat,
    X_list    = X_list,
    y_list    = y_list,
    mu_list   = mu_list,
    W_list    = W_list,
    V_list    = V_list,
    R_list    = R_list,
    DT_list   = DT_list,
    r_list    = r_list,
    H_list    = H_list,
    I0        = I0,
    Delta     = Delta
  )
}


## Helper: estimate working correlation from standardized residuals
estimate_corr <- function(e_list, corr_type) {
  N <- length(e_list)
  n_i <- sapply(e_list, length)

  if (corr_type == "independence") {
    R_list <- lapply(n_i, function(ni) diag(ni))
    return(R_list)
  }

  if (corr_type == "exchangeable") {
    # Average pairwise products within clusters
    alpha_vals <- sapply(e_list, function(e) {
      ni <- length(e)
      if (ni < 2) return(0)
      s <- sum(outer(e, e)) - sum(e^2)
      s / (ni * (ni - 1))
    })
    alpha <- mean(alpha_vals)
    R_list <- lapply(n_i, function(ni) clf_xch(ni, alpha))
  } else {
    # AR(1): average lag-1 products
    alpha_vals <- sapply(e_list, function(e) {
      ni <- length(e)
      if (ni < 2) return(0)
      sum(e[-ni] * e[-1]) / (ni - 1)
    })
    alpha <- mean(alpha_vals)
    R_list <- lapply(n_i, function(ni) clf_ar1(ni, alpha))
  }
  R_list
}


## Helper: compute working covariance V_i = phi * W_i^{1/2} R_i W_i^{1/2}
compute_V <- function(W_list, R_list, phi) {
  mapply(function(W, R) {
    W12 <- sqrt(W)
    phi * W12 %*% R %*% W12
  }, W_list, R_list, SIMPLIFY = FALSE)
}


init_pgee_beta <- function(y_full, x_full, p) {
  init_dat <- data.frame(y = y_full, x_full)

  if (requireNamespace("logistf", quietly = TRUE)) {
    beta <- tryCatch(
      as.vector(logistf::logistf(
        y ~ .,
        family = binomial,
        data = init_dat,
        control = logistf::logistf.control(maxit = 50)
      )$coef),
      error = function(e) rep(0, p)
    )
    if (length(beta) == p && all(is.finite(beta))) return(beta)
  }

  beta <- tryCatch(
    as.vector(stats::coef(stats::glm(y ~ ., family = stats::binomial(), data = init_dat))),
    warning = function(w) rep(0, p),
    error = function(e) rep(0, p)
  )
  if (length(beta) != p) return(rep(0, p))
  beta[!is.finite(beta)] <- 0
  beta
}
