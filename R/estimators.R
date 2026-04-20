## estimators.R
## Variance estimators for PGEE.
## Each function takes a common argument list (from pgee_fit output)
## and returns a p x p variance-covariance matrix.
##
## Dependencies: MASS (for ginv), hat_matrix.R (for mat_power_IminusH)
##
## Notation (from derivations.qmd / notation-inventory.md):
##   N       = number of clusters
##   p       = number of parameters
##   n_i     = cluster size for cluster i
##   n_star  = sum(n_i)
##   Delta   = I_0^{-1}
##   H_ii    = D_i Delta D_i^T V_i^{-1}  (n_i x n_i)
##   r_i     = y_i - mu_hat_i  (n_i x 1)
##   d_i     = D_i^T V_i^{-1} r_i  (p x 1, score contribution)
##   D_i^T   = X_i^T W_i  (p x n_i, stored as DT_list)
##   W_i     = diag(mu_ij(1-mu_ij))
##   V_i     = phi W_i^{1/2} R_i W_i^{1/2}

## Check if all clusters have the same size (required for pooling estimators)
is_balanced <- function(n_i) length(unique(n_i)) == 1

# Covariance estimators should be symmetric. We only enforce that explicitly for
# the leverage-adjusted estimators, where matrix-power and inversion paths can
# introduce small numerical asymmetries. The remaining estimators are left
# untouched because they are symmetric by construction.
symmetrize_covariance <- function(V) {
  0.5 * (V + t(V))
}

#' Internal estimator groups used for targeted dispatch.
#'
#' @keywords internal
#' @noRd
score_level_estimators <- function() {
  c("V_LZ", "V_DF", "V_FG", "V_MBN")
}

## Common argument extractor: unpack fit object into components
## fit: output from pgee_fit()
## Returns named list of all needed quantities
#' Internal fit unpacker for variance estimators.
#'
#' @keywords internal
#' @noRd
unpack_fit <- function(fit, include_scores = TRUE) {
  N <- fit$N
  p <- fit$p
  n_i <- fit$n_i
  n_star <- fit$n_star

  d_list <- if (isTRUE(include_scores)) vector("list", N) else NULL
  Vi_inv_list <- vector("list", N)
  for (i in 1:N) {
    Vi_inv_list[[i]] <- MASS::ginv(fit$V_list[[i]])
    if (isTRUE(include_scores)) {
      # d_i = D_i^T V_i^{-1} r_i
      d_list[[i]] <- as.vector(fit$DT_list[[i]] %*% Vi_inv_list[[i]] %*% fit$r_list[[i]])
    }
  }

  list(
    N = N, p = p, n_i = n_i, n_star = n_star,
    Delta = fit$Delta, I0 = fit$I0,
    DT_list = fit$DT_list, V_list = fit$V_list,
    W_list = fit$W_list, R_list = fit$R_list,
    r_list = fit$r_list, H_list = fit$H_list,
    mu_list = fit$mu_list, y_list = fit$y_list,
    d_list = d_list, Vi_inv_list = Vi_inv_list,
    phi = fit$phi
  )
}

#' Internal targeted estimator dispatcher.
#'
#' @keywords internal
#' @noRd
compute_requested_estimators <- function(fit, estimator_names) {
  estimator_names <- as.character(estimator_names)
  requested_names <- unique(estimator_names)
  needs_scores <- any(requested_names %in% score_level_estimators())
  u <- unpack_fit(fit, include_scores = needs_scores)

  H_cross <- NULL
  if ("V_FZ" %in% requested_names) {
    H_cross <- compute_H_cross(fit)
  }

  results <- stats::setNames(vector("list", length(requested_names)), requested_names)

  for (estimator_name in requested_names) {
    results[[estimator_name]] <- switch(
      estimator_name,
      V_LZ = v_lz(u),
      V_DF = v_df(u),
      V_KC = v_kc(u),
      V_MD = v_md(u),
      V_FG = v_fg(u),
      V_MBN = v_morel(u),
      V_Pan = v_pan(u),
      V_GST = v_gsk(u),
      V_WL = v_wl(u),
      V_WB = v_wb(u),
      V_RS = v_rs(u),
      V_FW = v_fw(u),
      V_FZ = v_fz(u, H_cross),
      V_AR = v_ar(u),
      stop(
        "Unknown estimator: ", estimator_name,
        ". Use available_estimator_names() or available_estimators() to see valid names."
      )
    )
  }

  results[estimator_names]
}


###############################################################################
## 1. V_LZ (Liang-Zeger, uncorrected sandwich)  [Eq A.1]
###############################################################################
v_lz <- function(u) {
  mid <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    mid <- mid + u$d_list[[i]] %*% t(u$d_list[[i]])
  }
  u$Delta %*% mid %*% u$Delta
}


###############################################################################
## 2. V_DF (MacKinnon-White: N/(N-p) * V_LZ)  [Eq A.2]
###############################################################################
v_df <- function(u) {
  (u$N / (u$N - u$p)) * v_lz(u)
}


###############################################################################
## 3. V_KC (Kauermann-Carroll: (I-H)^{-1/2} on residuals)  [Eq A.3]
###############################################################################
v_kc <- function(u) {
  mid <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    IH_half_inv <- mat_power_IminusH(u$H_list[[i]], c = 0.5)
    IHT_half_inv <- t(IH_half_inv)
    adj_r <- IH_half_inv %*% u$r_list[[i]]
    adj_rT <- t(u$r_list[[i]]) %*% IHT_half_inv
    mid <- mid + u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% adj_r %*%
      adj_rT %*% u$Vi_inv_list[[i]] %*% t(u$DT_list[[i]])
  }
  symmetrize_covariance(u$Delta %*% mid %*% u$Delta)
}


###############################################################################
## 4. V_MD (Mancl-DeRouen: (I-H)^{-1} on residuals)  [Eq A.4]
###############################################################################
v_md <- function(u) {
  mid <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    IH_inv <- mat_power_IminusH(u$H_list[[i]], c = 1)
    IHT_inv <- t(IH_inv)
    adj_r <- IH_inv %*% u$r_list[[i]]
    adj_rT <- t(u$r_list[[i]]) %*% IHT_inv
    mid <- mid + u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% adj_r %*%
      adj_rT %*% u$Vi_inv_list[[i]] %*% t(u$DT_list[[i]])
  }
  symmetrize_covariance(u$Delta %*% mid %*% u$Delta)
}


###############################################################################
## 5. V_FG (Fay-Graubard: diagonal clipping at d=0.75)  [Eq A.5]
###############################################################################
v_fg <- function(u, b_clip = 0.75) {
  mid <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    # N_i = D_i^T V_i^{-1} D_i Delta  (p x p)
    Ni <- u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% t(u$DT_list[[i]]) %*% u$Delta
    # F_i = diag{ (1 - min(b, [N_i]_ss))^{-1/2} }
    diag_Ni <- diag(Ni)
    fi_diag <- (1 - pmin(b_clip, diag_Ni))^(-0.5)
    Fi <- diag(fi_diag, nrow = u$p)
    # FG correction at score level
    di <- u$d_list[[i]]
    mid <- mid + Fi %*% (di %*% t(di)) %*% Fi
  }
  u$Delta %*% mid %*% u$Delta
}


###############################################################################
## 6. V_MBN (FPC + additive trace inflation)  [Eqs 2.7-2.9]
###############################################################################
v_morel <- function(u) {
  # I_1 = (n*-1)/(n*-p) * N/(N-1) * sum (d_i - d_bar)(d_i - d_bar)^T
  d_bar <- Reduce("+", u$d_list) / u$N
  I1 <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    diff <- u$d_list[[i]] - d_bar
    I1 <- I1 + diff %*% t(diff)
  }
  fpc <- (u$n_star - 1) / (u$n_star - u$p)
  bessel <- u$N / (u$N - 1)
  I1 <- fpc * bessel * I1

  # Sandwich
  Vs <- u$Delta %*% I1 %*% u$Delta

  # Additive inflation
  kappa <- max(1, sum(diag(u$Delta %*% I1)) / u$p)
  delta_n <- min(0.5, u$p / (u$N - u$p))
  Vs + kappa * delta_n * u$Delta
}


###############################################################################
## 7. V_Pan (pooled covariance)  [Eq A.6]
###############################################################################
v_pan <- function(u) {
  # Pooling requires balanced clusters (R_u is n_i x n_i)
  if (!is_balanced(u$n_i)) return(matrix(NA, u$p, u$p))

  ni <- u$n_i[1]
  # R_u = (1/N) sum W_i^{-1/2} r_i r_i^T W_i^{-1/2}
  R_u <- matrix(0, ni, ni)
  for (i in 1:u$N) {
    W_half_inv <- diag(1 / sqrt(diag(u$W_list[[i]])), nrow = ni)
    rrt <- u$r_list[[i]] %*% t(u$r_list[[i]])
    R_u <- R_u + W_half_inv %*% rrt %*% W_half_inv
  }
  R_u <- R_u / u$N

  # Middle term: sum D_i^T V_i^{-1} W_i^{1/2} R_u W_i^{1/2} V_i^{-1} D_i
  mid <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    W_half <- diag(sqrt(diag(u$W_list[[i]])), nrow = ni)
    G_i <- W_half %*% R_u %*% W_half
    mid <- mid + u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% G_i %*%
      u$Vi_inv_list[[i]] %*% t(u$DT_list[[i]])
  }
  u$Delta %*% mid %*% u$Delta
}


###############################################################################
## 8. V_GST (Gosho: Pan + 1/(N-p) DF correction)  [Eq A.7]
###############################################################################
v_gsk <- function(u) {
  if (!is_balanced(u$n_i)) return(matrix(NA, u$p, u$p))

  ni <- u$n_i[1]
  # Same as Pan but with 1/(N-p) instead of 1/N in pooled R_u
  R_u <- matrix(0, ni, ni)
  for (i in 1:u$N) {
    W_half_inv <- diag(1 / sqrt(diag(u$W_list[[i]])), nrow = ni)
    rrt <- u$r_list[[i]] %*% t(u$r_list[[i]])
    R_u <- R_u + W_half_inv %*% rrt %*% W_half_inv
  }
  R_u <- R_u / (u$N - u$p)

  mid <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    W_half <- diag(sqrt(diag(u$W_list[[i]])), nrow = ni)
    G_i <- W_half %*% R_u %*% W_half
    mid <- mid + u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% G_i %*%
      u$Vi_inv_list[[i]] %*% t(u$DT_list[[i]])
  }
  u$Delta %*% mid %*% u$Delta
}


###############################################################################
## 9. V_WL (Wang-Long: Pan pooling + MD correction)  [Eq A.8]
###############################################################################
v_wl <- function(u) {
  if (!is_balanced(u$n_i)) return(matrix(NA, u$p, u$p))

  ni <- u$n_i[1]
  # R_u^{MD} = (1/N) sum W^{-1/2} (I-H)^{-1} r r^T (I-H^T)^{-1} W^{-1/2}
  R_u <- matrix(0, ni, ni)
  for (i in 1:u$N) {
    W_half_inv <- diag(1 / sqrt(diag(u$W_list[[i]])), nrow = ni)
    IH_inv <- mat_power_IminusH(u$H_list[[i]], c = 1)
    IHT_inv <- t(IH_inv)
    adj_r <- IH_inv %*% u$r_list[[i]]
    adj_rT <- t(u$r_list[[i]]) %*% IHT_inv
    R_u <- R_u + W_half_inv %*% adj_r %*% adj_rT %*% W_half_inv
  }
  R_u <- R_u / u$N

  mid <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    W_half <- diag(sqrt(diag(u$W_list[[i]])), nrow = ni)
    G_i <- W_half %*% R_u %*% W_half
    mid <- mid + u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% G_i %*%
      u$Vi_inv_list[[i]] %*% t(u$DT_list[[i]])
  }
  symmetrize_covariance(u$Delta %*% mid %*% u$Delta)
}


###############################################################################
## 10. V_WB (Westgate-Burchett: Pan pooling + (I-H)^{-c}, c=1/2)  [Eq A.9]
###############################################################################
v_wb <- function(u, c_param = 0.5) {
  if (!is_balanced(u$n_i)) return(matrix(NA, u$p, u$p))

  ni <- u$n_i[1]
  # R_u^{(c)} = (1/N) sum W^{-1/2} (I-H)^{-c} r r^T (I-H^T)^{-c} W^{-1/2}
  R_u <- matrix(0, ni, ni)
  for (i in 1:u$N) {
    W_half_inv <- diag(1 / sqrt(diag(u$W_list[[i]])), nrow = ni)
    IH_c_inv <- mat_power_IminusH(u$H_list[[i]], c = c_param)
    IHT_c_inv <- t(IH_c_inv)
    adj_r <- IH_c_inv %*% u$r_list[[i]]
    adj_rT <- t(u$r_list[[i]]) %*% IHT_c_inv
    R_u <- R_u + W_half_inv %*% adj_r %*% adj_rT %*% W_half_inv
  }
  R_u <- R_u / u$N

  mid <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    W_half <- diag(sqrt(diag(u$W_list[[i]])), nrow = ni)
    G_i <- W_half %*% R_u %*% W_half
    mid <- mid + u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% G_i %*%
      u$Vi_inv_list[[i]] %*% t(u$DT_list[[i]])
  }
  symmetrize_covariance(u$Delta %*% mid %*% u$Delta)
}


###############################################################################
## 11. V_RS (Rogers-Stoner: Pan + determinant-based Morel inflation)  [Eq A.11]
###############################################################################
v_rs <- function(u) {
  if (!is_balanced(u$n_i)) return(matrix(NA, u$p, u$p))

  V_pan_mat <- v_pan(u)

  ni <- u$n_i[1]
  # Compute Pan middle term for the inflation factor
  R_u <- matrix(0, ni, ni)
  for (i in 1:u$N) {
    W_half_inv <- diag(1 / sqrt(diag(u$W_list[[i]])), nrow = ni)
    rrt <- u$r_list[[i]] %*% t(u$r_list[[i]])
    R_u <- R_u + W_half_inv %*% rrt %*% W_half_inv
  }
  R_u <- R_u / u$N

  mid_pan <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    W_half <- diag(sqrt(diag(u$W_list[[i]])), nrow = ni)
    G_i <- W_half %*% R_u %*% W_half
    mid_pan <- mid_pan + u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% G_i %*%
      u$Vi_inv_list[[i]] %*% t(u$DT_list[[i]])
  }

  # d_det = max(1, |det(Delta I_1)|^{1/p})  [Rogers-Stoner Eq 5]
  det_mat <- u$Delta %*% mid_pan
  d_det <- max(1, abs(det(det_mat))^(1 / u$p))

  delta_n <- min(0.5, u$p / (u$N - u$p))

  V_pan_mat + d_det * delta_n * u$Delta
}


###############################################################################
## 12. V_FW (Ford-Westgate: (V_KC + V_MD)/2)  [Eq A.12]
###############################################################################
v_fw <- function(u) {
  0.5 * (v_kc(u) + v_md(u))
}


###############################################################################
## 13. V_FZ (Fan-Zhang-Zhang: MD + cross-cluster subtraction)  [Eq A.13]
## Requires H_cross from compute_H_cross(fit)
###############################################################################
v_fz <- function(u, H_cross) {
  mid <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    IH_inv <- mat_power_IminusH(u$H_list[[i]], c = 1)
    IHT_inv <- t(IH_inv)

    # Cross-cluster contamination: sum_{j != i} H_{ij} r_j r_j^T H_{ij}^T
    cross <- matrix(0, u$n_i[i], u$n_i[i])
    for (j in 1:u$N) {
      if (j == i) next
      H_ij <- H_cross[[i]][[j]]
      cross <- cross + H_ij %*% (u$r_list[[j]] %*% t(u$r_list[[j]])) %*% t(H_ij)
    }

    # Corrected inner term: r_i r_i^T - cross
    inner <- u$r_list[[i]] %*% t(u$r_list[[i]]) - cross

    adj_inner <- IH_inv %*% inner %*% IHT_inv
    mid <- mid + u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% adj_inner %*%
      u$Vi_inv_list[[i]] %*% t(u$DT_list[[i]])
  }
  symmetrize_covariance(u$Delta %*% mid %*% u$Delta)
}


###############################################################################
## 14. V_AR (Afiaz-Rahman: (I-H)^{-1} on scores + full Morel FPC)  [Eq 5.2]
###############################################################################
v_ar <- function(u) {
  # f_i = D_i^T V_i^{-1} (I - H_ii)^{-1} r_i
  f_list <- vector("list", u$N)
  for (i in 1:u$N) {
    IH_inv <- mat_power_IminusH(u$H_list[[i]], c = 1)
    f_list[[i]] <- as.vector(u$DT_list[[i]] %*% u$Vi_inv_list[[i]] %*% IH_inv %*% u$r_list[[i]])
  }
  f_bar <- Reduce("+", f_list) / u$N

  # M_AR = FPC * Bessel * sum (f_i - f_bar)(f_i - f_bar)^T
  M_ar <- matrix(0, u$p, u$p)
  for (i in 1:u$N) {
    diff <- f_list[[i]] - f_bar
    M_ar <- M_ar + diff %*% t(diff)
  }
  fpc <- (u$n_star - 1) / (u$n_star - u$p)
  bessel <- u$N / (u$N - 1)
  M_ar <- fpc * bessel * M_ar

  u$Delta %*% M_ar %*% u$Delta
}


###############################################################################
## Master function: compute all estimators
###############################################################################
# fit: output from pgee_fit()
# compute_cross: if TRUE, compute cross-cluster hat blocks for V_FZ (O(N^2))
# Returns: named list of p x p matrices
#
# Active set (14): 13 literature estimators plus V_AR from this manuscript.
# Literature: V_LZ, V_DF, V_KC, V_MD, V_FG, V_MBN, V_Pan, V_GST, V_WL, V_WB,
#             V_RS, V_FW, V_FZ
# Ours:       V_AR
# Pooling estimators (V_Pan, V_GST, V_WL, V_WB, V_RS) return NA under
# unbalanced clusters by design.
compute_all_estimators <- function(fit, compute_cross = TRUE) {
  estimator_names <- estimator_notation_names()

  if (!isTRUE(compute_cross)) {
    estimator_names <- setdiff(estimator_names, "V_FZ")
  }

  compute_requested_estimators(fit, estimator_names)
}
