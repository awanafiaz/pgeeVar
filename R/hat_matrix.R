## hat_matrix.R
## Hat matrix computation and matrix square root utilities.
## H_{ii} is already computed in pgee_fit.R. This file provides:
##   - Cross-cluster hat blocks H_{ij} (needed only for V_FZ)
##   - Eigendecomposition-based matrix power: (I - H_{ii})^{-c}
## No external dependencies beyond base R and MASS.

# Compute cross-cluster hat block H_{ij} = D_i Delta D_j^T V_j^{-1}
# Only needed for V_FZ. Call with compute_cross = TRUE in run_scenario.
# fit: output from pgee_fit()
# Returns: N x N list-of-lists, H_cross[[i]][[j]] for j != i
compute_H_cross <- function(fit) {
  N <- fit$N
  # Precompute V_j^{-1} for all j
  Vj_inv_list <- lapply(fit$V_list, MASS::ginv)
  H_cross <- vector("list", N)
  for (i in 1:N) {
    H_cross[[i]] <- vector("list", N)
    for (j in 1:N) {
      if (i == j) next
      # H_{ij} = D_i Delta D_j^T V_j^{-1}
      # t(DT_i) is n_i x p, Delta is p x p, DT_j is p x n_j, Vj_inv is n_j x n_j
      H_cross[[i]][[j]] <- t(fit$DT_list[[i]]) %*% fit$Delta %*%
        fit$DT_list[[j]] %*% Vj_inv_list[[j]]
    }
  }
  H_cross
}


# Compute (I - H_{ii})^{-c} via eigendecomposition
# H_ii: n_i x n_i hat matrix block (generally asymmetric)
# c: power parameter (1 for MD, 0.5 for KC)
# Returns: n_i x n_i matrix
#
# H_{ii} is similar to a symmetric PSD matrix via:
#   H_{ii} = V_i^{1/2} tilde{H}_{ii} V_i^{-1/2}
# so its eigenvalues are real and in [0, 1). We eigendecompose H_{ii}
# directly; for small imaginary parts (floating point noise), we take
# the real part. For genuinely complex eigenvalues, we fall back to ginv.
mat_power_IminusH <- function(H_ii, c) {
  n <- nrow(H_ii)

  eig <- eigen(H_ii)

  # Check for complex eigenvalues (signals numerical degeneracy)
  if (any(abs(Im(eig$values)) > 1e-8)) {
    IminusH <- diag(n) - H_ii
    if (c == 1) return(MASS::ginv(IminusH))
    if (c == 0.5) {
      # Compute (I-H)^{-1/2} via symmetrized pseudoinverse.
      # G = ginv(I-H) is not symmetric, so eigen(G) has complex eigenvectors
      # and P D^{1/2} P^T != P D^{1/2} P^{-1}. Symmetrize first.
      G <- MASS::ginv(IminusH)
      G_sym <- (G + t(G)) / 2
      eig_G <- eigen(G_sym, symmetric = TRUE)
      vals_G <- pmax(eig_G$values, 0)
      G_half <- eig_G$vectors %*% diag(sqrt(vals_G), nrow = n) %*% t(eig_G$vectors)
      return(G_half)
    }
    warning("mat_power_IminusH: complex eigenvalues with c=", c,
            "; falling back to ginv(I-H)")
    return(MASS::ginv(IminusH))
  }

  vals <- Re(eig$values)
  vecs <- Re(eig$vectors)

  # Clamp eigenvalues of H to [0, 1-eps] for numerical safety
  vals <- pmin(pmax(vals, 0), 1 - 1e-10)

  D_new <- diag((1 - vals)^(-c), nrow = n)
  # Use ginv instead of solve for numerical stability with near-singular
  # eigenvector matrices (common with small clusters)
  result <- vecs %*% D_new %*% MASS::ginv(vecs)

  Re(result)
}


# Symmetric positive definite matrix square root via eigendecomposition
mat_sqrt_sym <- function(A) {
  eig <- eigen(A, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  eig$vectors %*% diag(sqrt(vals), nrow = length(vals)) %*% t(eig$vectors)
}
