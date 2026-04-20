test_that("compute_pgee_variances returns the full package set and preserves request order", {
  engine_fit <- make_engine_fit_from_balanced()
  wrapper_fit <- make_balanced_fit()
  all_names <- available_estimator_names()

  vars_default <- compute_pgee_variances(engine_fit)
  vars_all <- compute_pgee_variances(engine_fit, "all")
  vars_subset <- compute_pgee_variances(engine_fit, c("AR", "MD", "V_AR"))
  vars_wrapper <- compute_pgee_variances(wrapper_fit, c("MD", "AR"))

  expect_equal(names(vars_default), all_names)
  expect_equal(vars_default, vars_all, tolerance = 1e-10)
  expect_equal(names(vars_subset), c("AR", "MD", "AR"))
  expect_equal(vars_subset[[1]], vars_default$AR, tolerance = 1e-10)
  expect_equal(vars_subset[[2]], vars_default$MD, tolerance = 1e-10)
  expect_equal(vars_subset[[3]], vars_default$AR, tolerance = 1e-10)
  expect_identical(vars_wrapper, wrapper_fit$vcovs[c("MD", "AR")])
  expect_equal(rownames(vars_wrapper[[1]]), names(coef(wrapper_fit)))
  expect_equal(colnames(vars_wrapper[[1]]), names(coef(wrapper_fit)))
  expect_error(
    compute_pgee_variances(engine_fit, "V_ZZ"),
    "available_estimator_names"
  )
  expect_error(
    compute_pgee_variances(list(), "AR"),
    "output of pgee_fit\\(\\) or pgee\\(\\)"
  )
})

test_that("variance identities and unbalanced behavior are locked down", {
  engine_fit <- make_engine_fit_from_balanced()
  vars <- compute_pgee_variances(engine_fit)
  finite_names <- names(vars)[vapply(vars, function(mat) all(is.finite(mat)), logical(1))]

  expect_equal(
    vars$DF,
    (engine_fit$N / (engine_fit$N - engine_fit$p)) * vars$LZ,
    tolerance = 1e-10
  )
  expect_equal(
    vars$FW,
    (vars$KC + vars$MD) / 2,
    tolerance = 1e-10
  )

  for (name in finite_names) {
    expect_equal(vars[[name]], t(vars[[name]]), tolerance = 1e-12)
  }

  unbalanced_fit <- make_unbalanced_fit()
  vars_ub <- compute_pgee_variances(unbalanced_fit)
  pooled_names <- c("Pan", "GST", "WL", "WB", "RS")
  non_pooling_names <- c("LZ", "DF", "KC", "MD", "FG", "MBN", "FW", "FZ", "AR")
  expected_dim <- rep(length(coef(unbalanced_fit)), 2)

  for (name in pooled_names) {
    expect_all_na_matrix(vars_ub[[name]], expected_dim = expected_dim)
  }

  expect_true(all(vapply(
    vars_ub[non_pooling_names],
    function(mat) all(is.finite(mat)),
    logical(1)
  )))
})
