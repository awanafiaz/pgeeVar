test_that("simulate_correlated_binary_data returns analysis-ready output and scenario metadata", {
  sim_data <- make_balanced_sim(output = "data")
  sim_full <- make_balanced_sim(seed = 1, output = "full")
  scenario_names <- c("N", "gamma", "n_i", "beta", "rho", "corstr", "clf_ok", "N_min", "n_star", "event_rate")

  expect_s3_class(sim_data, "data.frame")
  expect_false("intercept" %in% names(sim_data))
  expect_equal(unique(sim_data$id), seq_len(length(unique(sim_data$id))))
  expect_equal(sort(names(attr(sim_data, "pgee_scenario"))), sort(scenario_names))

  expect_true(is.list(sim_full))
  expect_s3_class(sim_full$data, "data.frame")
  expect_false("intercept" %in% names(sim_full$data))
})

test_that("simulation wrapper is deterministic and handles RNG state as documented", {
  sim_a <- make_balanced_sim(seed = 99)
  sim_b <- make_balanced_sim(seed = 99)

  expect_equal(sim_a, sim_b)

  set.seed(123)
  expected_next <- runif(5)
  set.seed(123)
  simulate_correlated_binary_data(
    N = 20,
    gamma = 0.5,
    n_i = 4,
    beta = c(-1.2, 0.6, 0.3),
    rho = 0.2,
    seed = 42
  )
  observed_next <- runif(5)
  expect_equal(observed_next, expected_next, tolerance = 1e-10)

  set.seed(123)
  before_state <- .Random.seed
  simulate_correlated_binary_data(
    N = 20,
    gamma = 0.5,
    n_i = 4,
    beta = c(-1.2, 0.6, 0.3),
    rho = 0.2,
    seed = NULL
  )
  after_state <- .Random.seed
  expect_false(identical(before_state, after_state))
})

test_that("simulation wrapper validates user inputs", {
  expect_error(
    simulate_correlated_binary_data(
      N = 10, gamma = 0.5, n_i = 4, beta = c(-1.2, 0.6, 0.3), rho = 0.2, seed = 1.5
    ),
    "single finite integer"
  )
  expect_error(
    simulate_correlated_binary_data(
      N = 1, gamma = 0.5, n_i = 4, beta = c(-1.2, 0.6, 0.3), rho = 0.2
    ),
    "`N`"
  )
  expect_error(
    simulate_correlated_binary_data(
      N = 10, gamma = 1, n_i = 4, beta = c(-1.2, 0.6, 0.3), rho = 0.2
    ),
    "`gamma`"
  )
  expect_error(
    simulate_correlated_binary_data(
      N = 10, gamma = 0.5, n_i = 1, beta = c(-1.2, 0.6, 0.3), rho = 0.2
    ),
    "`n_i`"
  )
  expect_error(
    simulate_correlated_binary_data(
      N = 10, gamma = 0.5, n_i = 4, beta = c(-1.2), rho = 0.2
    ),
    "`beta`"
  )
  expect_error(
    simulate_correlated_binary_data(
      N = 10, gamma = 0.5, n_i = 4, beta = c(-1.2, 0.6, 0.3), rho = 1
    ),
    "`rho`"
  )
  expect_error(
    simulate_correlated_binary_data(
      N = 10, gamma = 0.5, n_i = 4, beta = c(-1.2, 0.6, 0.3), rho = 0.2, corstr = "bad"
    )
  )
  expect_error(
    simulate_correlated_binary_data(
      N = 10, gamma = 0.5, n_i = 4, beta = c(-1.2, 0.6, 0.3), rho = 0.2, on_clf_failure = "ignore"
    )
  )
  expect_error(
    simulate_correlated_binary_data(
      N = 10, gamma = 0.5, n_i = 4, beta = c(-1.2, 0.6, 0.3), rho = 0.2, on_clf_failure = "foo"
    )
  )
})
