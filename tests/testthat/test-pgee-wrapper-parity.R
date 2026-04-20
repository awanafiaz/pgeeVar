# Scenarios in this file are chosen to be stable under both the optional
# logistf initialization path and the package's fallback initialization path.

test_that("wrapper and low-level engine agree on the same processed data", {
  data <- make_balanced_sim()
  fit_wrapper <- pgee(y ~ X1 + obstime, data = data, id = id)
  fit_engine <- make_engine_fit_from_balanced(data = data)
  strip_dimnames <- function(x) {
    dimnames(x) <- NULL
    x
  }

  expect_equal(unname(stats::coef(fit_wrapper)), unname(fit_engine$beta), tolerance = 1e-10)
  expect_equal(fit_wrapper$alpha_hat, fit_engine$alpha_hat, tolerance = 1e-10)
  expect_equal(fit_wrapper$phi, fit_engine$phi, tolerance = 1e-10)
  expect_identical(fit_wrapper$converged, fit_engine$converged)

  vars_wrapper <- compute_pgee_variances(fit_wrapper)
  vars_engine <- compute_pgee_variances(fit_engine)

  expect_equal(rownames(vars_wrapper$AR), names(coef(fit_wrapper)))
  expect_equal(colnames(vars_wrapper$AR), names(coef(fit_wrapper)))
  expect_equal(strip_dimnames(vars_wrapper$AR), vars_engine$AR, tolerance = 1e-10)
  expect_equal(strip_dimnames(vars_wrapper$MD), vars_engine$MD, tolerance = 1e-10)
  expect_equal(strip_dimnames(vars_wrapper$FZ), vars_engine$FZ, tolerance = 1e-10)
  expect_equal(
    lapply(vars_wrapper, strip_dimnames),
    vars_engine,
    tolerance = 1e-10
  )
})
