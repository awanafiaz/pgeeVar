test_that("dgp_generate returns the expected low-level scenario structure", {
  set.seed(1)
  sim <- dgp_generate(
    N = 12,
    gamma = 0.5,
    n_i = 4,
    beta = c(-1.2, 0.6, 0.3),
    rho = 0.2
  )

  expect_true(all(c(
    "data", "covnames", "N", "N_min", "N1", "n_i", "n_star",
    "p", "clf_ok", "beta", "rho", "corr_type"
  ) %in% names(sim)))
  expect_s3_class(sim$data, "data.frame")
  expect_true(all(c("id", "y", "intercept", "X1", "obstime") %in% names(sim$data)))
  expect_equal(sim$covnames, c("X1", "obstime"))
  expect_equal(sim$N, 12)
  expect_equal(sim$p, 3)
})

test_that("pgee_fit returns the expected low-level fit structure", {
  set.seed(1)
  sim <- dgp_generate(
    N = 20,
    gamma = 0.5,
    n_i = 4,
    beta = c(-1.2, 0.6, 0.3),
    rho = 0.2
  )

  fit <- pgee_fit(
    y = sim$data$y,
    x = sim$data[c("X1", "obstime")],
    id = sim$data$id,
    corr_type = sim$corr_type
  )

  expect_true(all(c(
    "beta", "converged", "iter", "N", "p", "n_i", "n_star", "phi",
    "alpha_hat", "X_list", "y_list", "mu_list", "W_list", "V_list",
    "R_list", "DT_list", "r_list", "H_list", "I0", "Delta"
  ) %in% names(fit)))
  expect_equal(length(fit$beta), 3L)
  expect_true(all(is.finite(fit$beta)))
  expect_true(is.logical(fit$converged))
  expect_true(is.numeric(fit$alpha_hat))
})
