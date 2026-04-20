## dgp.R
## Data generating process for PGEE variance estimation simulation.
## Generates longitudinal binary data with cluster-level and time covariates.
## Dependencies: clf_generate.R (for clf_xch, clf_ar1, clf_cor2var, clf_allReg,
##               clf_blrchk1, clf_mbsclf)
##
## Supports p=2 (intercept + X1), p=3 (+ time), p=4 (+ continuous x3).

#' Generate correlated binary longitudinal data
#'
#' `dgp_generate()` creates one simulated binary longitudinal dataset for the
#' PGEE settings studied in the `pgeeVar` manuscript. The function uses the
#' vendored Qaqish-style conditional linear family (CLF) machinery to generate
#' correlated binary outcomes under exchangeable or AR(1) dependence.
#'
#' This is the low-level simulation engine underlying
#' [simulate_correlated_binary_data()]. The returned object includes both the
#' observation-level dataset and the scenario metadata used to generate it, so
#' it can be passed directly into [pgee_fit()] or reused in simulation scripts.
#'
#' @param N Number of clusters.
#' @param gamma Proportion of clusters assigned the binary cluster-level
#'   covariate `X1 = 1`.
#' @param n_i Cluster-size specification. Use a single integer for a balanced
#'   design, or a length-2 vector `c(min, max)` to draw cluster sizes uniformly
#'   between those limits.
#' @param beta Regression coefficients, including the intercept. Supported
#'   lengths are 2, 3, and 4:
#'   * length 2: intercept and `X1`
#'   * length 3: intercept, `X1`, and `obstime`
#'   * length 4: intercept, `X1`, `obstime`, and a cluster-level continuous
#'     covariate `X3`
#' @param rho True within-cluster correlation parameter used by the CLF
#'   generator.
#' @param corr_type Correlation structure used for data generation. Must be
#'   either `"exchangeable"` or `"ar1"`.
#'
#' @details
#' The returned dataset always includes an `id` column and a binary response
#' `y`. The covariate columns depend on the length of `beta`.
#'
#' If the CLF construction fails for any cluster, `clf_ok` is set to `FALSE`
#' and the affected response values are returned as `NA`. The downstream
#' [pgee_fit()] pipeline is still safe in that case because the fitting routine
#' removes incomplete rows with `complete.cases()` before estimation.
#'
#' @return A list with the following components:
#'   * `data`: observation-level data frame containing `id`, `y`, and the
#'     generated covariates.
#'   * `covnames`: names of the non-intercept covariates to pass into
#'     [pgee_fit()].
#'   * `N`, `N_min`, `N1`, `n_i`, `n_star`, `p`: basic design summaries.
#'   * `clf_ok`: logical flag indicating whether CLF generation succeeded for all
#'     clusters.
#'   * `beta`, `rho`, `corr_type`: the generating parameters used for the
#'     scenario.
#'
#' @references
#' Qaqish BF (2003). A family of multivariate binary distributions for
#' simulating correlated binary variables with specified marginal means and
#' correlations. *Biometrika*, 90(2), 455-463.
#'
#' The vendored helper layer also follows the archived `binarySimCLF` package
#' implementation of By and Qaqish (2009), CRAN package source.
#'
#' @seealso [simulate_correlated_binary_data()], [pgee_fit()],
#'   [compute_pgee_variances()]
#'
#' @examples
#' set.seed(1)
#' sim <- dgp_generate(
#'   N = 12,
#'   gamma = 0.5,
#'   n_i = 4,
#'   beta = c(-1.2, 0.6, 0.3),
#'   rho = 0.2
#' )
#'
#' sim$clf_ok
#' head(sim$data)
#' sim$covnames
#'
#' @export
dgp_generate <- function(N, gamma, n_i, beta, rho,
                          corr_type = c("exchangeable", "ar1")) {
  corr_type <- match.arg(corr_type)
  p <- length(beta)

  # Determine cluster sizes
  if (length(n_i) == 1) {
    cl_sizes <- rep(n_i, N)
  } else {
    cl_sizes <- round(runif(N, min = n_i[1], max = n_i[2]))
  }
  n_star <- sum(cl_sizes)

  # Generate cluster-level binary covariate (treatment)
  N1 <- max(1, round(gamma * N))
  x1_cluster <- c(rep(1, N1), rep(0, N - N1))
  x1_cluster <- sample(x1_cluster)
  N_min <- min(sum(x1_cluster), sum(1 - x1_cluster))

  # Build observation-level data
  id <- rep(1:N, times = cl_sizes)
  X1 <- rep(x1_cluster, times = cl_sizes)

  if (p == 2) {
    # p=2: intercept + X1 only
    Xmat <- cbind(intercept = 1, X1 = X1)
    covariate_names <- c("X1")
  } else if (p == 3) {
    # p=3: intercept + X1 + time
    obstime <- unlist(lapply(cl_sizes, function(ni) (1:ni) * 0.2))
    Xmat <- cbind(intercept = 1, X1 = X1, obstime = obstime)
    covariate_names <- c("X1", "obstime")
  } else if (p == 4) {
    # p=4: intercept + X1 + time + x3 (cluster-level continuous)
    obstime <- unlist(lapply(cl_sizes, function(ni) (1:ni) * 0.2))
    x3_cluster <- rnorm(N)
    X3 <- rep(x3_cluster, times = cl_sizes)
    Xmat <- cbind(intercept = 1, X1 = X1, obstime = obstime, X3 = X3)
    covariate_names <- c("X1", "obstime", "X3")
  } else {
    stop("p must be 2, 3, or 4")
  }

  # Linear predictor and marginal means
  eta <- as.vector(Xmat %*% beta)
  mu <- 1 / (1 + exp(-eta))

  # Split mu by cluster
  mu_list <- split(mu, id)

  # Generate correlated binary outcomes via CLF
  clf_ok <- TRUE
  y_all <- numeric(n_star)
  idx <- 1

  for (i in 1:N) {
    ni <- cl_sizes[i]
    mu_i <- mu_list[[i]]

    if (corr_type == "exchangeable") {
      R_i <- clf_xch(ni, rho)
    } else {
      R_i <- clf_ar1(ni, rho)
    }

    V_i <- clf_cor2var(R_i, mu_i)
    B_i <- clf_allReg(V_i)
    if (!clf_blrchk1(mu_i, B_i)) {
      clf_ok <- FALSE
      y_all[idx:(idx + ni - 1)] <- NA
      idx <- idx + ni
      next
    }

    result <- clf_mbsclf(m = 1, u = mu_i, B = B_i)
    if (!result$succeed) {
      clf_ok <- FALSE
      y_all[idx:(idx + ni - 1)] <- NA
    } else {
      y_all[idx:(idx + ni - 1)] <- result$y[1, ]
    }
    idx <- idx + ni
  }

  # Build output data frame
  dat <- data.frame(id = id, y = y_all, Xmat)
  # Remove intercept column from the covariate frame passed to pgee_fit
  # (pgee_fit adds its own intercept)

  list(
    data      = dat,
    covnames  = covariate_names,
    N         = N,
    N_min     = N_min,
    N1        = N1,
    n_i       = cl_sizes,
    n_star    = n_star,
    p         = p,
    clf_ok    = clf_ok,
    beta      = beta,
    rho       = rho,
    corr_type = corr_type
  )
}
