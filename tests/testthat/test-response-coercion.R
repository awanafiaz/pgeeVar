test_that("numeric, logical, and factor responses produce the same fit", {
  fit_numeric <- make_balanced_fit(response = "numeric")
  fit_logical <- make_balanced_fit(response = "logical")
  fit_factor <- make_balanced_fit(response = "factor")

  expect_equal(unname(coef(fit_numeric)), unname(coef(fit_logical)), tolerance = 1e-10)
  expect_equal(unname(coef(fit_numeric)), unname(coef(fit_factor)), tolerance = 1e-10)
})

test_that("invalid response encodings are rejected", {
  sim_bad_numeric <- make_balanced_sim()
  sim_bad_numeric$y[1] <- 2

  sim_bad_factor <- make_balanced_sim()
  sim_bad_factor$y <- factor(rep(c("a", "b", "c"), length.out = nrow(sim_bad_factor)))

  expect_error(
    pgee(y ~ X1 + obstime, data = sim_bad_numeric, id = id),
    "coded as 0/1"
  )
  expect_error(
    pgee(y ~ X1 + obstime, data = sim_bad_factor, id = id),
    "two-level factor"
  )
})
