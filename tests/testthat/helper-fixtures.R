options(warn = 2)

make_balanced_sim <- function(seed = 42, output = "data") {
  simulate_correlated_binary_data(
    N = 30,
    gamma = 0.5,
    n_i = 4,
    beta = c(-1.2, 0.6, 0.3),
    rho = 0.2,
    seed = seed,
    output = output
  )
}

make_unbalanced_sim <- function(seed = 7, output = "data") {
  simulate_correlated_binary_data(
    N = 25,
    gamma = 0.5,
    n_i = c(3, 5),
    beta = c(-1.2, 0.6, 0.3),
    rho = 0.2,
    seed = seed,
    output = output
  )
}

coerce_response_for_test <- function(data, response = c("numeric", "logical", "factor")) {
  response <- match.arg(response)

  if (identical(response, "logical")) {
    data$y <- as.logical(data$y)
  } else if (identical(response, "factor")) {
    data$y <- factor(
      ifelse(data$y == 1, "yes", "no"),
      levels = c("no", "yes")
    )
  }

  data
}

make_balanced_fit <- function(
  data = NULL,
  id_mode = c("bare", "character", "vector"),
  response = c("numeric", "logical", "factor"),
  corstr = "exchangeable",
  variance = "AR",
  na.action = stats::na.omit,
  max_iter = 100
) {
  id_mode <- match.arg(id_mode)
  response <- match.arg(response)

  if (is.null(data)) {
    data <- make_balanced_sim()
  }

  data <- coerce_response_for_test(data, response)

  if (identical(id_mode, "bare")) {
    return(
      pgee(
        y ~ X1 + obstime,
        data = data,
        id = id,
        corstr = corstr,
        variance = variance,
        na.action = na.action,
        max_iter = max_iter
      )
    )
  }

  if (identical(id_mode, "character")) {
    return(
      pgee(
        y ~ X1 + obstime,
        data = data,
        id = "id",
        corstr = corstr,
        variance = variance,
        na.action = na.action,
        max_iter = max_iter
      )
    )
  }

  pgee(
    y ~ X1 + obstime,
    data = data,
    id = data$id,
    corstr = corstr,
    variance = variance,
    na.action = na.action,
    max_iter = max_iter
  )
}

make_unbalanced_fit <- function(
  data = NULL,
  corstr = "exchangeable",
  variance = "AR"
) {
  if (is.null(data)) {
    data <- make_unbalanced_sim()
  }

  pgee(
    y ~ X1 + obstime,
    data = data,
    id = id,
    corstr = corstr,
    variance = variance
  )
}

make_engine_fit_from_balanced <- function(data = NULL, corstr = "exchangeable") {
  if (is.null(data)) {
    data <- make_balanced_sim()
  }

  pgee_fit(
    y = data$y,
    x = data[c("X1", "obstime")],
    id = data$id,
    corr_type = corstr
  )
}

expect_all_na_matrix <- function(x, expected_dim = NULL) {
  testthat::expect_true(is.matrix(x))
  if (is.null(expected_dim)) {
    testthat::expect_equal(dim(x), rep(nrow(x), 2))
  } else {
    testthat::expect_equal(dim(x), expected_dim)
  }
  testthat::expect_true(all(is.na(x)))
}

expect_same_coef <- function(object, expected, tolerance = 1e-10) {
  testthat::expect_equal(
    unname(stats::coef(object)),
    unname(stats::coef(expected)),
    tolerance = tolerance
  )
}
