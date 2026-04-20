test_that("pgee accepts all supported id modes and subset forms", {
  fit_bare <- make_balanced_fit(id_mode = "bare")
  fit_char <- make_balanced_fit(id_mode = "character")
  fit_vec <- make_balanced_fit(id_mode = "vector")
  sim <- make_balanced_sim()
  keep_vec <- sim$X1 == 1
  fit_subset_expr <- pgee(y ~ X1 + obstime, data = sim, id = id, subset = X1 == 1)
  fit_subset_vec <- pgee(y ~ X1 + obstime, data = sim, id = id, subset = keep_vec)
  fit_subset_data <- pgee(y ~ X1 + obstime, data = sim[sim$X1 == 1, , drop = FALSE], id = id)

  expect_same_coef(fit_bare, fit_char)
  expect_same_coef(fit_bare, fit_vec)
  expect_same_coef(fit_subset_expr, fit_subset_vec)
  expect_same_coef(fit_subset_expr, fit_subset_data)
})

test_that("pgee validates id, formula, variance, and convergence edges", {
  sim <- make_balanced_sim()
  fit_one_iter <- NULL
  fit_alias <- make_balanced_fit(variance = "V_AR")
  expect_warning(
    fit_one_iter <- pgee(y ~ X1 + obstime, data = sim, id = id, max_iter = 1),
    "did not converge"
  )

  expect_false(fit_one_iter$converged)
  expect_identical(fit_alias$default_vcov, "AR")
  expect_error(
    pgee(y ~ X1 + obstime, data = sim, id = does_not_exist),
    "`id` must be a column in `data` or a vector of length nrow\\(data\\)"
  )
  expect_error(
    pgee(y ~ X1 + obstime, data = sim, id = "does_not_exist"),
    "`id` must be a column in `data` or a vector of length nrow\\(data\\)"
  )
  expect_error(
    pgee(y ~ 0 + X1 + obstime, data = sim, id = id),
    "requires an intercept"
  )
  expect_error(
    pgee(y ~ 1, data = sim, id = id),
    "requires at least one covariate"
  )
  expect_error(
    pgee(y ~ X1 + obstime + id, data = sim, id = id),
    "must not also appear on the right-hand side"
  )
  expect_error(
    pgee(y ~ X1 + obstime, data = sim, id = id, variance = "V_ZZ"),
    "available_estimator_names"
  )

  sim_bad <- sim[c(1, 5:8), , drop = FALSE]
  expect_error(
    pgee(y ~ X1 + obstime, data = sim_bad, id = id),
    "Each cluster must contain at least two observations"
  )
})

test_that("pgee handles na.action choices and subset plus na.action interaction", {
  sim_na <- make_balanced_sim()
  sim_na$obstime[c(1, 5)] <- NA

  fit_omit_fun <- pgee(y ~ X1 + obstime, data = sim_na, id = id, na.action = stats::na.omit)
  fit_omit_chr <- pgee(y ~ X1 + obstime, data = sim_na, id = id, na.action = "omit")

  expect_same_coef(fit_omit_fun, fit_omit_chr)
  expect_equal(nobs(fit_omit_fun), nobs(fit_omit_chr))
  expect_error(
    pgee(y ~ X1 + obstime, data = sim_na, id = id, na.action = stats::na.fail),
    "missing values in object"
  )
  expect_error(
    pgee(y ~ X1 + obstime, data = sim_na, id = id, na.action = "fail"),
    "missing values in object"
  )

  sim_subset_na <- make_balanced_sim()
  subset_rows <- which(sim_subset_na$X1 == 1)
  sim_subset_na$obstime[subset_rows[c(1, 5)]] <- NA
  sim_subset_na$y[subset_rows[9]] <- NA

  fit_subset_na <- pgee(
    y ~ X1 + obstime,
    data = sim_subset_na,
    id = id,
    subset = X1 == 1,
    na.action = stats::na.omit
  )

  expected_n <- sum(
    sim_subset_na$X1 == 1 &
      complete.cases(sim_subset_na[, c("y", "X1", "obstime", "id")]),
    na.rm = TRUE
  )
  expect_equal(nobs(fit_subset_na), expected_n)
})
