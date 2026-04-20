resolve_na_action <- function(na.action) {
  if (is.character(na.action)) {
    if (length(na.action) != 1L) {
      stop("`na.action` must be `na.omit`, `na.fail`, \"omit\", or \"fail\".")
    }

    return(
      switch(
        na.action,
        omit = stats::na.omit,
        fail = stats::na.fail,
        stop("`na.action` must be `na.omit`, `na.fail`, \"omit\", or \"fail\".")
      )
    )
  }

  if (identical(na.action, stats::na.omit)) {
    return(stats::na.omit)
  }

  if (identical(na.action, stats::na.fail)) {
    return(stats::na.fail)
  }

  stop("`na.action` must be `na.omit`, `na.fail`, \"omit\", or \"fail\".")
}

resolve_id_spec <- function(id_expr, data, env) {
  if (is.symbol(id_expr)) {
    id_name <- as.character(id_expr)
    if (id_name %in% names(data)) {
      return(list(values = data[[id_name]], name = id_name))
    }
  }

  id_value <- tryCatch(
    eval(id_expr, envir = env),
    error = function(e) {
      stop("`id` must be a column in `data` or a vector of length nrow(data).")
    }
  )

  if (is.character(id_value) && length(id_value) == 1L && id_value %in% names(data)) {
    return(list(values = data[[id_value]], name = id_value))
  }

  if (length(id_value) != nrow(data)) {
    stop("`id` must be a column in `data` or a vector of length nrow(data).")
  }

  list(values = id_value, name = NULL)
}

coerce_binary_response <- function(y) {
  if (is.logical(y)) {
    return(as.integer(y))
  }

  if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      stop("The response must be binary: numeric 0/1, logical, or a two-level factor.")
    }
    return(as.integer(y == levels(y)[2L]))
  }

  if (is.numeric(y)) {
    if (!all(y %in% c(0, 1))) {
      stop("The response must be coded as 0/1.")
    }
    return(as.numeric(y))
  }

  stop("The response must be binary: numeric 0/1, logical, or a two-level factor.")
}

extract_xlevels <- function(model_frame, terms_obj) {
  rhs_vars <- all.vars(stats::delete.response(terms_obj))
  xlevels <- lapply(rhs_vars, function(var_name) {
    if (!var_name %in% names(model_frame)) {
      return(NULL)
    }

    column <- model_frame[[var_name]]
    if (is.factor(column)) levels(column) else NULL
  })
  names(xlevels) <- rhs_vars
  xlevels[!vapply(xlevels, is.null, logical(1))]
}

resolve_vcov_type <- function(object, type) {
  actual_type <- if (is.null(type)) object$default_vcov else type
  normalized <- normalize_estimator_names(actual_type)

  if (length(normalized) != 1L) {
    stop("`type` must be a single estimator name.")
  }

  normalized
}

resolve_parm_indices <- function(parm, coefficient_names) {
  if (is.null(parm)) {
    return(seq_along(coefficient_names))
  }

  if (is.character(parm)) {
    idx <- match(parm, coefficient_names)
    if (anyNA(idx)) {
      stop("Unknown coefficient name(s): ", paste(parm[is.na(idx)], collapse = ", "))
    }
    return(idx)
  }

  if (is.logical(parm)) {
    if (length(parm) != length(coefficient_names)) {
      stop("Logical `parm` must have length equal to the number of coefficients.")
    }
    return(which(parm))
  }

  if (is.numeric(parm)) {
    idx <- as.integer(parm)
    if (any(idx < 1L | idx > length(coefficient_names))) {
      stop("Numeric `parm` indices must be between 1 and ", length(coefficient_names), ".")
    }
    return(idx)
  }

  stop("`parm` must be NULL, character, logical, or numeric.")
}

pooling_estimator_names <- local({
  cached_names <- NULL

  function() {
    if (is.null(cached_names)) {
      catalog <- available_estimators()
      cached_names <<- catalog$name[catalog$class == "pooled"]
    }
    cached_names
  }
})

name_vcov_matrix <- function(vcov_mat, coefficient_names) {
  if (is.null(dimnames(vcov_mat))) {
    dimnames(vcov_mat) <- list(coefficient_names, coefficient_names)
  }
  vcov_mat
}

strip_simulated_data <- function(sim) {
  keep_cols <- c("id", "y", sim$covnames)
  sim$data[, keep_cols, drop = FALSE]
}

#' Simulate correlated binary longitudinal data
#'
#' `simulate_correlated_binary_data()` is the user-facing wrapper around the
#' lower-level [dgp_generate()] engine. It validates inputs, optionally manages
#' the RNG seed, and returns either a plain analysis-ready data frame or the
#' full scenario object.
#'
#' @param N Number of clusters. Must be an integer greater than or equal to 2.
#' @param gamma Proportion of clusters assigned to `X1 = 1`. Must lie strictly
#'   between 0 and 1.
#' @param n_i Cluster-size specification. Use a single integer for a balanced
#'   design or a length-2 vector `c(min, max)` for a discrete uniform range.
#'   Every realized cluster size must be at least 2.
#' @param beta Regression coefficients, including the intercept. Supported
#'   lengths are 2, 3, and 4.
#' @param rho Correlation parameter used by the CLF generator. Must satisfy
#'   `abs(rho) < 1`.
#' @param corstr Correlation structure. Must be `"exchangeable"` or `"ar1"`.
#' @param seed Optional integer seed. If `NULL`, the current RNG state is left
#'   unchanged. If supplied, the wrapper temporarily sets the seed and restores
#'   the prior RNG state on exit.
#' @param on_clf_failure How to handle CLF-generation failures. `"warn"` returns
#'   the generated object with a warning; `"error"` stops immediately.
#' @param output Either `"data"` for a plain data frame or `"full"` for the
#'   full scenario list.
#'
#' @details
#' The `"data"` output mode strips the internal intercept column used by the
#' engine and returns only `id`, `y`, and the user-facing covariates. The same
#' stripped frame is stored in `$data` when `output = "full"`, so the wrapper's
#' output intentionally differs from the lower-level [dgp_generate()] engine
#' object.
#'
#' When `output = "data"`, the returned data frame carries a
#' `"pgee_scenario"` attribute containing the generating metadata, including the
#' realized event rate.
#'
#' @return
#' If `output = "data"`, a data frame with columns `id`, `y`, and the simulated
#' covariates, plus a `"pgee_scenario"` attribute.
#'
#' If `output = "full"`, a list containing the full low-level scenario object,
#' with its `$data` component replaced by the stripped user-facing data frame.
#'
#' @seealso [pgee()], [dgp_generate()], [pgee_fit()]
#'
#' @examples
#' sim_data <- simulate_correlated_binary_data(
#'   N = 12,
#'   gamma = 0.5,
#'   n_i = 4,
#'   beta = c(-1.2, 0.6, 0.3),
#'   rho = 0.2,
#'   seed = 1
#' )
#'
#' head(sim_data)
#' attr(sim_data, "pgee_scenario")$clf_ok
#'
#' sim_full <- simulate_correlated_binary_data(
#'   N = 12,
#'   gamma = 0.5,
#'   n_i = 4,
#'   beta = c(-1.2, 0.6, 0.3),
#'   rho = 0.2,
#'   seed = 1,
#'   output = "full"
#' )
#'
#' sim_full$clf_ok
#' names(sim_full$data)
#'
#' @export
simulate_correlated_binary_data <- function(
  N,
  gamma,
  n_i,
  beta,
  rho,
  corstr = c("exchangeable", "ar1"),
  seed = NULL,
  on_clf_failure = c("warn", "error"),
  output = c("data", "full")
) {
  corstr <- match.arg(corstr)
  on_clf_failure <- match.arg(on_clf_failure)
  output <- match.arg(output)

  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 2 || N != as.integer(N)) {
    stop("`N` must be a single integer greater than or equal to 2.")
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma <= 0 || gamma >= 1) {
    stop("`gamma` must be a single number strictly between 0 and 1.")
  }
  if (!is.numeric(n_i) || !length(n_i) %in% c(1L, 2L) || any(!is.finite(n_i)) || any(n_i < 2)) {
    stop("`n_i` must be a numeric vector of length 1 or 2 with values at least 2.")
  }
  if (!is.numeric(beta) || !length(beta) %in% c(2L, 3L, 4L) || any(!is.finite(beta))) {
    stop("`beta` must be a numeric vector of length 2, 3, or 4.")
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop("`rho` must be a single finite number with absolute value less than 1.")
  }
  if (!is.null(seed)) {
    seed_int <- suppressWarnings(as.integer(seed))
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || is.na(seed_int) || seed != seed_int) {
      stop("`seed` must be NULL or a single finite integer value.")
    }
  }

  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }

    set.seed(seed_int)

    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  sim <- dgp_generate(
    N = as.integer(N),
    gamma = gamma,
    n_i = n_i,
    beta = beta,
    rho = rho,
    corr_type = corstr
  )

  data_out <- strip_simulated_data(sim)
  event_rate <- mean(data_out$y, na.rm = TRUE)
  scenario <- list(
    N = sim$N,
    gamma = gamma,
    n_i = sim$n_i,
    beta = sim$beta,
    rho = sim$rho,
    corstr = sim$corr_type,
    clf_ok = sim$clf_ok,
    N_min = sim$N_min,
    n_star = sim$n_star,
    event_rate = event_rate
  )
  attr(data_out, "pgee_scenario") <- scenario

  if (!sim$clf_ok) {
    message <- "CLF generation failed for at least one cluster; affected responses were set to NA."
    if (identical(on_clf_failure, "error")) {
      stop(message)
    }
    warning(message)
  }

  if (is.finite(event_rate) && event_rate %in% c(0, 1)) {
    warning("The realized event rate is degenerate (all outcomes are 0 or all are 1).")
  }

  if (identical(output, "data")) {
    return(data_out)
  }

  sim$data <- data_out
  sim
}

#' Fit a formula-based penalized generalized estimating equation model
#'
#' `pgee()` is the user-facing formula wrapper around the lower-level
#' [pgee_fit()] engine. It uses standard model-frame and model-matrix machinery,
#' caches all 14 supported variance estimators at fit time, and returns an
#' object of class `"pgee"` with a first-pass set of S3 methods for inference,
#' prediction, and extraction.
#'
#' @param formula A model formula with a binary response on the left-hand side.
#'   The formula must include an intercept and at least one covariate beyond the
#'   intercept.
#' @param data A data frame containing the variables in `formula` and, if `id`
#'   is given as a column name, the clustering variable.
#' @param id The clustering variable. This may be given as a bare column name, a
#'   character column name, or a vector of length `nrow(data)`.
#' @param corstr Working correlation structure. Must be one of
#'   `"exchangeable"`, `"ar1"`, or `"independence"`.
#' @param variance The default variance estimator name to use in
#'   [vcov.pgee()], [summary.pgee()], [confint.pgee()], and
#'   [predict.pgee()] when no explicit type is supplied. All 14 supported
#'   estimators are still precomputed and stored on the fitted object. Canonical
#'   package names are the short forms returned by [available_estimator_names()];
#'   legacy `V_*` aliases are still accepted.
#' @param subset Optional subset expression, evaluated in `data`.
#' @param na.action Missing-data action. Supported values are `na.omit`,
#'   `na.fail`, `"omit"`, and `"fail"`.
#' @param tol Convergence tolerance passed to [pgee_fit()].
#' @param max_iter Maximum number of scoring iterations passed to [pgee_fit()].
#'
#' @details
#' `pgee()` removes the intercept column from the model matrix before calling the
#' low-level engine because [pgee_fit()] adds its own intercept internally.
#'
#' The fitted object stores both the low-level engine output and the full set of
#' 14 cached covariance matrices so later calls to [vcov.pgee()],
#' [summary.pgee()], or [confint.pgee()] do not need to recompute them.
#'
#' The response may be supplied as numeric 0/1, logical, or a two-level factor.
#' Two-level factors are converted using the second level as the event level.
#'
#' @return An object of class `"pgee"`, containing named coefficients, cached
#'   covariance matrices, fitted values, residuals, model-frame components, and
#'   the underlying low-level engine fit.
#'
#' @seealso [simulate_correlated_binary_data()], [pgee_fit()],
#'   [compute_pgee_variances()], [available_estimators()]
#'
#' @examples
#' sim_data <- simulate_correlated_binary_data(
#'   N = 20,
#'   gamma = 0.5,
#'   n_i = 4,
#'   beta = c(-1.2, 0.6, 0.3),
#'   rho = 0.2,
#'   seed = 1
#' )
#'
#' fit <- pgee(y ~ X1 + obstime, data = sim_data, id = id)
#'
#' coef(fit)
#' summary(fit)
#' vcov(fit, type = "MD")
#' head(fitted(fit))
#'
#' @export
pgee <- function(
  formula,
  data,
  id,
  corstr = c("exchangeable", "ar1", "independence"),
  variance = "AR",
  subset = NULL,
  na.action = stats::na.omit,
  tol = 1e-4,
  max_iter = 100
) {
  corstr <- match.arg(corstr)

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a valid model formula.")
  }

  data <- as.data.frame(data)
  if (nrow(data) == 0L) {
    stop("`data` must contain at least one row.")
  }

  variance <- normalize_estimator_names(variance)
  if (length(variance) != 1L) {
    stop("`variance` must be a single supported estimator name.")
  }

  na_fun <- resolve_na_action(na.action)
  terms_obj <- stats::terms(formula, data = data)

  if (attr(terms_obj, "response") != 1L) {
    stop("`formula` must include a response variable.")
  }
  if (attr(terms_obj, "intercept") != 1L) {
    stop("`pgee()` currently requires an intercept in the formula.")
  }
  if (length(attr(terms_obj, "offset")) > 0L) {
    stop("Offsets are not supported in `pgee()` at this stage.")
  }

  id_info <- resolve_id_spec(substitute(id), data, parent.frame())
  rhs_vars <- all.vars(stats::delete.response(terms_obj))
  if (!is.null(id_info$name) && id_info$name %in% rhs_vars) {
    stop("The clustering variable supplied through `id` must not also appear on the right-hand side of `formula`.")
  }

  temp_rowid <- ".pgee_rowid_internal"
  temp_id <- ".pgee_id_internal"
  temp_subset <- ".pgee_subset_internal"
  data_with_rowid <- data
  data_with_rowid[[temp_rowid]] <- seq_len(nrow(data_with_rowid))
  subset_expr <- substitute(subset)
  has_subset <- !missing(subset) && !identical(subset_expr, quote(NULL))

  response_label <- paste(deparse(formula[[2L]]), collapse = "")
  formula_with_rowid <- stats::reformulate(
    termlabels = c(attr(terms_obj, "term.labels"), temp_rowid),
    response = response_label,
    intercept = attr(terms_obj, "intercept"),
    env = environment(formula)
  )

  if (!has_subset) {
    model_frame_raw <- stats::model.frame(
      formula_with_rowid,
      data = data_with_rowid,
      na.action = stats::na.pass,
      drop.unused.levels = TRUE
    )
  } else {
    subset_value <- eval(subset_expr, envir = data_with_rowid, enclos = parent.frame())
    # Pass the evaluated subset through a sentinel column so model.frame()
    # resolves it in the data context while keeping lm()/glm()-style behavior.
    data_with_rowid[[temp_subset]] <- subset_value
    model_frame_raw <- stats::model.frame(
      formula_with_rowid,
      data = data_with_rowid,
      subset = .pgee_subset_internal,
      na.action = stats::na.pass,
      drop.unused.levels = TRUE
    )
  }

  row_ids <- model_frame_raw[[temp_rowid]]
  model_frame_raw[[temp_rowid]] <- NULL

  id_subset <- id_info$values[row_ids]
  if (anyNA(id_subset)) {
    stop("Missing values are not allowed in `id`.")
  }

  model_frame_with_id <- model_frame_raw
  model_frame_with_id[[temp_id]] <- id_subset
  model_frame_used <- na_fun(model_frame_with_id)
  id_used <- model_frame_used[[temp_id]]
  model_frame_used[[temp_id]] <- NULL

  model_terms <- stats::terms(formula, data = model_frame_used)
  design_matrix <- stats::model.matrix(stats::delete.response(model_terms), model_frame_used)
  coefficient_names <- colnames(design_matrix)

  if (ncol(design_matrix) <= 1L) {
    stop("`pgee()` requires at least one covariate in addition to the intercept.")
  }

  response <- coerce_binary_response(stats::model.response(model_frame_used))
  if (length(unique(id_used)) < 2L) {
    stop("`pgee()` requires at least two clusters after subsetting and missing-value handling.")
  }

  cluster_sizes <- table(id_used)
  if (any(cluster_sizes < 2L)) {
    bad_clusters <- names(cluster_sizes)[cluster_sizes < 2L]
    stop(
      "Each cluster must contain at least two observations. ",
      "The following cluster(s) have only one observation after preprocessing: ",
      paste(bad_clusters, collapse = ", "), "."
    )
  }

  engine_fit <- pgee_fit(
    y = response,
    x = as.data.frame(design_matrix[, -1L, drop = FALSE]),
    id = id_used,
    corr_type = corstr,
    tol = tol,
    max_iter = max_iter
  )

  if (!engine_fit$converged) {
    warning("`pgee()` did not converge after ", max_iter, " iterations.")
  }

  vcovs <- compute_pgee_variances(engine_fit, "all")
  vcovs <- lapply(vcovs, name_vcov_matrix, coefficient_names = coefficient_names)
  coefficients <- engine_fit$beta
  names(coefficients) <- coefficient_names
  linear_predictors <- as.vector(design_matrix %*% engine_fit$beta)
  fitted_values <- stats::plogis(linear_predictors)
  residuals_response <- response - fitted_values
  xlevels <- extract_xlevels(model_frame_used, model_terms)

  structure(
    list(
      coefficients = coefficients,
      vcovs = vcovs,
      default_vcov = variance,
      fitted_values = fitted_values,
      linear_predictors = linear_predictors,
      residuals = residuals_response,
      response = response,
      call = match.call(),
      formula = formula,
      terms = model_terms,
      xlevels = xlevels,
      model_frame = model_frame_used,
      model_matrix = design_matrix,
      id = id_used,
      corstr = corstr,
      alpha_hat = engine_fit$alpha_hat,
      phi = engine_fit$phi,
      converged = engine_fit$converged,
      iter = engine_fit$iter,
      N = engine_fit$N,
      p = engine_fit$p,
      n_i = engine_fit$n_i,
      n_star = engine_fit$n_star,
      df.residual = engine_fit$N - engine_fit$p,
      engine_fit = engine_fit
    ),
    class = "pgee"
  )
}

#' Methods for `pgee` objects
#'
#' These methods provide the first-pass user interface for fitted objects
#' returned by [pgee()].
#'
#' @param object A fitted `"pgee"` object.
#' @param x A `"summary.pgee"` object.
#' @param formula A fitted `"pgee"` object for [model.frame.pgee()].
#' @param type Variance-estimator name for methods that need a covariance
#'   matrix. If `NULL`, the object's stored default variance is used. Canonical
#'   package names are the short forms returned by [available_estimator_names()];
#'   legacy `V_*` aliases are still accepted.
#' @param parm Coefficients to extract in [confint.pgee()]. May be `NULL`,
#'   character names, logical positions, or numeric indices.
#' @param level Confidence level for [confint.pgee()].
#' @param newdata Optional new data for [predict.pgee()]. If `NULL`, the method
#'   uses the training data.
#' @param se.fit Logical; if `TRUE`, [predict.pgee()] also returns standard
#'   errors.
#' @param vcov_type Variance-estimator name to use for prediction standard
#'   errors. If `NULL`, the object's stored default variance is used. Canonical
#'   package names are the short forms returned by [available_estimator_names()];
#'   legacy `V_*` aliases are still accepted.
#' @param ... Additional arguments, currently unused.
#'
#' @name pgee_methods
NULL

#' @rdname pgee_methods
#' @export
print.pgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat(
    "\nWorking correlation:", x$corstr,
    "\nAlpha estimate:", format(x$alpha_hat, digits = 4),
    "\nClusters:", x$N,
    "\nConverged:", if (isTRUE(x$converged)) "yes" else "no",
    "after", x$iter, "iteration(s)",
    "\nDefault variance:", x$default_vcov, "\n"
  )
  invisible(x)
}

#' @rdname pgee_methods
#' @return `coef.pgee()` returns the named coefficient vector.
#' @export
coef.pgee <- function(object, ...) {
  object$coefficients
}

#' @rdname pgee_methods
#' @return `vcov.pgee()` returns a single covariance matrix.
#' @export
vcov.pgee <- function(object, type = NULL, ...) {
  actual_type <- resolve_vcov_type(object, type)
  vcov_mat <- object$vcovs[[actual_type]]
  vcov_mat <- name_vcov_matrix(vcov_mat, names(object$coefficients))

  if (anyNA(vcov_mat) && actual_type %in% pooling_estimator_names() &&
      length(unique(object$n_i)) != 1L) {
    warning(
      "Estimator ", actual_type,
      " is not directly computable for unbalanced cluster sizes; returning an all-NA matrix."
    )
  }

  vcov_mat
}

#' @rdname pgee_methods
#' @return `summary.pgee()` returns an object of class `"summary.pgee"`.
#' @export
summary.pgee <- function(object, type = NULL, ...) {
  actual_type <- resolve_vcov_type(object, type)
  vcov_mat <- stats::vcov(object, type = actual_type)
  df_resid <- object$df.residual

  if (df_resid <= 0) {
    stop("Residual degrees of freedom `N - p` must be positive for t-based inference.")
  }

  std_error <- sqrt(diag(vcov_mat))
  t_value <- object$coefficients / std_error
  p_value <- 2 * stats::pt(-abs(t_value), df = df_resid)
  coef_table <- cbind(
    Estimate = object$coefficients,
    `Std. Error` = std_error,
    `t value` = t_value,
    `Pr(>|t|)` = p_value
  )

  structure(
    list(
      call = object$call,
      coefficients = coef_table,
      variance_type = actual_type,
      df = df_resid,
      corstr = object$corstr,
      alpha_hat = object$alpha_hat,
      phi = object$phi,
      converged = object$converged,
      iter = object$iter,
      N = object$N,
      n_i_summary = c(
        min = min(object$n_i),
        median = stats::median(object$n_i),
        max = max(object$n_i)
      ),
      event_rate = mean(object$response),
      default_vcov = object$default_vcov
    ),
    class = "summary.pgee"
  )
}

#' @rdname pgee_methods
#' @export
print.summary.pgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat(
    "\nVariance estimator:", x$variance_type,
    "\nResidual df:", x$df,
    "\nWorking correlation:", x$corstr,
    "\nAlpha estimate:", format(x$alpha_hat, digits = 4),
    "\nDispersion estimate:", format(x$phi, digits = 4),
    "\nClusters:", x$N,
    "\nCluster sizes (min/median/max):",
    paste(format(x$n_i_summary, digits = 4), collapse = " / "),
    "\nEvent rate:", format(x$event_rate, digits = 4),
    "\nConverged:", if (isTRUE(x$converged)) "yes" else "no",
    "after", x$iter, "iteration(s)\n\n"
  )
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}

#' @rdname pgee_methods
#' @return `confint.pgee()` returns a two-column confidence-interval matrix.
#' @export
confint.pgee <- function(object, parm = NULL, level = 0.95, type = NULL, ...) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.")
  }

  actual_type <- resolve_vcov_type(object, type)
  vcov_mat <- stats::vcov(object, type = actual_type)
  df_resid <- object$df.residual

  if (df_resid <= 0) {
    stop("Residual degrees of freedom `N - p` must be positive for t-based inference.")
  }

  idx <- resolve_parm_indices(parm, names(object$coefficients))
  std_error <- sqrt(diag(vcov_mat))[idx]
  estimates <- object$coefficients[idx]
  crit <- stats::qt(1 - (1 - level) / 2, df = df_resid)
  ci <- cbind(
    `Lower` = estimates - crit * std_error,
    `Upper` = estimates + crit * std_error
  )
  rownames(ci) <- names(object$coefficients)[idx]
  ci
}

#' @rdname pgee_methods
#' @return `fitted.pgee()` returns the response-scale fitted values.
#' @export
fitted.pgee <- function(object, ...) {
  object$fitted_values
}

#' @rdname pgee_methods
#' @return `residuals.pgee()` returns response or Pearson residuals.
#' @export
residuals.pgee <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)

  if (identical(type, "response")) {
    return(object$residuals)
  }

  variance_mu <- pmax(object$fitted_values * (1 - object$fitted_values), .Machine$double.eps)
  object$residuals / sqrt(variance_mu)
}

#' @rdname pgee_methods
#' @return `predict.pgee()` returns fitted values, or a list containing fitted
#'   values and standard errors when `se.fit = TRUE`.
#' @export
predict.pgee <- function(
  object,
  newdata = NULL,
  type = c("response", "link"),
  se.fit = FALSE,
  vcov_type = NULL,
  ...
) {
  type <- match.arg(type)
  actual_type <- resolve_vcov_type(object, vcov_type)

  if (is.null(newdata)) {
    design_matrix <- object$model_matrix
    eta <- object$linear_predictors
  } else {
    newdata <- as.data.frame(newdata)
    new_terms <- stats::delete.response(object$terms)
    new_model_frame <- stats::model.frame(
      new_terms,
      data = newdata,
      na.action = stats::na.pass,
      xlev = object$xlevels
    )
    design_matrix <- stats::model.matrix(new_terms, new_model_frame)
    eta <- as.vector(design_matrix %*% object$coefficients)
  }

  fit <- if (identical(type, "link")) eta else stats::plogis(eta)

  if (!isTRUE(se.fit)) {
    return(fit)
  }

  vcov_mat <- stats::vcov(object, type = actual_type)
  link_variance <- rowSums((design_matrix %*% vcov_mat) * design_matrix)
  link_se <- sqrt(pmax(link_variance, 0))
  se_values <- if (identical(type, "link")) {
    link_se
  } else {
    abs(fit * (1 - fit)) * link_se
  }

  list(
    fit = fit,
    se.fit = se_values,
    type = type,
    vcov_type = actual_type
  )
}

#' @rdname pgee_methods
#' @return `model.frame.pgee()` returns the stored training model frame.
#' @export
model.frame.pgee <- function(formula, ...) {
  formula$model_frame
}

#' @rdname pgee_methods
#' @return `model.matrix.pgee()` returns the stored training design matrix,
#'   including the intercept column.
#' @export
model.matrix.pgee <- function(object, ...) {
  object$model_matrix
}

#' @rdname pgee_methods
#' @return `formula.pgee()` returns the stored model formula.
#' @export
formula.pgee <- function(x, ...) {
  x$formula
}

#' @rdname pgee_methods
#' @return `nobs.pgee()` returns the total number of observations used in the fit.
#' @export
nobs.pgee <- function(object, ...) {
  object$n_star
}
