estimator_catalog <- function() {
  data.frame(
    name = c(
      "LZ", "DF", "KC", "MD", "FG", "MBN",
      "Pan", "GST", "WL", "WB", "RS",
      "FW", "FZ", "AR"
    ),
    notation = c(
      "V_LZ", "V_DF", "V_KC", "V_MD", "V_FG", "V_MBN",
      "V_Pan", "V_GST", "V_WL", "V_WB", "V_RS",
      "V_FW", "V_FZ", "V_AR"
    ),
    class = c(
      "sandwich", "sandwich", "leverage", "leverage", "leverage", "hybrid",
      "pooled", "pooled", "pooled", "pooled", "pooled",
      "hybrid", "hybrid", "score"
    ),
    requires_balance = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, TRUE, TRUE,
      FALSE, FALSE, FALSE
    ),
    reference_key = c(
      "liang_zeger_1986",
      "mackinnon_white_1985",
      "kauermann_carroll_2001",
      "mancl_derouen_2001",
      "fay_graubard_2001",
      "morel_bokossa_neerchal_2003",
      "pan_2001",
      "gosho_sato_takeuchi_2014",
      "wang_long_2011",
      "westgate_burchett_2016",
      "rogers_stoner_2015",
      "ford_westgate_2018",
      "fan_zhang_zhang_2013",
      "afiaz_rahman_2026"
    ),
    default_order = seq_len(14),
    description = c(
      "Uncorrected Liang-Zeger sandwich estimator.",
      "Degrees-of-freedom adjusted sandwich estimator.",
      "Kauermann-Carroll leverage correction.",
      "Mancl-DeRouen leverage correction.",
      "Fay-Graubard diagonal clipping correction.",
      "Morel-Bokossa-Neerchal hybrid correction.",
      "Pan pooled covariance estimator.",
      "Gosho pooled covariance estimator with DF adjustment.",
      "Wang-Long pooled covariance estimator with MD correction.",
      "Westgate-Burchett pooled leverage correction.",
      "Rogers-Stoner pooled hybrid correction.",
      "Ford-Westgate hybrid of KC and MD.",
      "Fan-Zhang-Zhang cross-cluster correction.",
      "Afiaz-Rahman score-based estimator."
    ),
    stringsAsFactors = FALSE
  )
}

estimator_notation_names <- function() {
  estimator_catalog()$notation
}

estimator_notation_from_name <- function(estimator_names) {
  catalog <- estimator_catalog()
  catalog$notation[match(estimator_names, catalog$name)]
}

normalize_estimator_names <- function(estimator_names, allow_all = FALSE) {
  if (length(estimator_names) == 0L) {
    return(character(0))
  }

  if (!is.character(estimator_names)) {
    estimator_names <- as.character(estimator_names)
  }

  if (isTRUE(allow_all) &&
      length(estimator_names) == 1L &&
      !is.na(estimator_names) &&
      identical(tolower(estimator_names), "all")) {
    return("all")
  }

  catalog <- estimator_catalog()
  alias_map <- stats::setNames(
    c(catalog$name, catalog$name),
    tolower(c(catalog$name, catalog$notation))
  )
  normalized <- unname(alias_map[tolower(estimator_names)])

  if (anyNA(normalized)) {
    bad <- estimator_names[is.na(normalized)]
    stop(
      "Unknown estimator(s): ", paste(bad, collapse = ", "),
      ". Use available_estimator_names() or available_estimators() to see valid names. ",
      "Legacy V_* aliases are still accepted."
    )
  }

  normalized
}

#' List the available variance estimators
#'
#' `available_estimators()` returns the estimator catalog shipped in
#' `pgeeVar`. It is the main reference for end users who want to see which
#' variance estimators the package currently supports and which names should be
#' passed to [compute_pgee_variances()].
#'
#' The catalog is meant to be read directly. In addition to the estimator name,
#' it records the broad estimator family, whether the method requires balanced
#' cluster sizes, and a short plain-language description.
#'
#' @return A data frame with one row per supported estimator and the following
#'   columns:
#'   * `name`: canonical estimator name used in user-facing function calls.
#'   * `notation`: manuscript/internal notation used in the package
#'     implementation.
#'   * `class`: broad estimator family such as `"sandwich"`, `"leverage"`,
#'     `"pooled"`, `"hybrid"`, or `"score"`.
#'   * `requires_balance`: logical flag indicating whether the estimator assumes
#'     equal cluster sizes.
#'   * `reference_key`: manuscript-style reference key for the estimator's main
#'     source.
#'   * `default_order`: display order used by the package.
#'   * `description`: short end-user description of what the estimator does.
#'
#' @references
#' The estimator catalog corresponds to the implementations studied in the
#' `pgee-variance` manuscript and traces to Liang and Zeger (1986),
#' *Biometrika*, 73(1), 13-22; MacKinnon and White (1985), *Journal of
#' Econometrics*, 29(3), 305-325; Kauermann and Carroll (2001), *Journal of the
#' American Statistical Association*, 96(456), 1387-1396; Mancl and DeRouen
#' (2001), *Biometrics*, 57(1), 126-134; Fay and Graubard (2001), *Biometrics*,
#' 57(4), 1198-1206; Morel, Bokossa, and Neerchal (2003), *Biometrical
#' Journal*, 45(4), 395-409; Pan (2001), *Biometrika*, 88(3), 901-906; Gosho,
#' Sato, and Takeuchi (2014), *Science Journal of Applied Mathematics and
#' Statistics*, 2(1), 20-25; Wang and Long (2011), *Statistics in Medicine*,
#' 30(11), 1278-1291; Westgate and Burchett (2016), *Statistics in Medicine*,
#' 35(21), 3733-3744; Rogers and Stoner (2015), *American Journal of Applied
#' Mathematics and Statistics*, 3(6), 243-251; Ford and Westgate (2018),
#' *Statistics in Medicine*, 37(28), 4318-4329; Fan, Zhang, and Zhang (2013),
#' *Journal of Biopharmaceutical Statistics*, 23(5), 1172-1187; and the current
#' Afiaz-Rahman manuscript estimator implementation.
#'
#' @seealso [available_estimator_names()], [compute_pgee_variances()]
#'
#' @examples
#' est <- available_estimators()
#' est[, c("name", "notation", "class", "requires_balance")]
#'
#' @export
available_estimators <- function() {
  out <- estimator_catalog()
  rownames(out) <- NULL
  out
}

#' Get the valid estimator names
#'
#' `available_estimator_names()` returns the canonical short names accepted by
#' the public estimator arguments in `pgeeVar`. Legacy `V_*` names are still
#' accepted as compatibility aliases, but the returned vector gives the names
#' shown in the package documentation and examples.
#'
#' @return A character vector of valid estimator names, in the package's default
#'   catalog order.
#'
#' @seealso [available_estimators()], [compute_pgee_variances()]
#'
#' @examples
#' available_estimator_names()
#'
#' @export
available_estimator_names <- function() {
  available_estimators()$name
}

#' Compute PGEE variance estimators
#'
#' `compute_pgee_variances()` computes one or more covariance-matrix estimators
#' from a fitted PGEE object produced by [pgee_fit()] or [pgee()]. Use
#' `estimators = "all"` to compute the full 14-estimator package set, or pass a
#' character vector of estimator names to request only the methods you want.
#'
#' This is currently an advanced API. The function returns the raw covariance
#' matrices rather than a formatted model summary.
#'
#' @param fit A fitted object returned by [pgee_fit()] or [pgee()].
#' @param estimators Either `"all"` for the full package set, or a character
#'   vector of estimator names returned by [available_estimator_names()].
#'   Legacy `V_*` aliases are accepted and normalized internally.
#'
#' @details
#' The output order matches the requested order in `estimators`.
#'
#' Estimators that require balanced cluster sizes are still listed in the
#' returned object for unbalanced designs, but their covariance matrices are
#' filled with `NA` values because the direct implementation is not computable in
#' that setting.
#'
#' If you supply an invalid estimator name, the error message points back to
#' [available_estimator_names()] and [available_estimators()] so you can inspect
#' the valid choices.
#'
#' @return A named list of covariance matrices. Each element is a `p x p` matrix,
#'   where `p` is the number of regression parameters in `fit`. For pooling
#'   estimators applied to unbalanced data, the returned matrix is still `p x p`
#'   but is filled with `NA` values because the direct implementation is not
#'   computable in that setting.
#'
#' @seealso [available_estimators()], [available_estimator_names()],
#'   [pgee_fit()]
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
#' fit <- pgee_fit(
#'   y = sim$data$y,
#'   x = sim$data[c("X1", "obstime")],
#'   id = sim$data$id,
#'   corr_type = sim$corr_type
#' )
#'
#' vars <- compute_pgee_variances(fit, c("MD", "AR"))
#' names(vars)
#' dim(vars$AR)
#' length(compute_pgee_variances(fit))
#'
#' @export
compute_pgee_variances <- function(fit, estimators = "all") {
  normalized_estimators <- normalize_estimator_names(estimators, allow_all = TRUE)

  if (length(normalized_estimators) == 1L && identical(normalized_estimators, "all")) {
    estimator_names <- available_estimator_names()
  } else {
    estimator_names <- normalized_estimators
  }

  if (inherits(fit, "pgee")) {
    coefficient_names <- names(fit$coefficients)
    result <- fit$vcovs[estimator_names]
    result[] <- lapply(result, function(vcov_mat) {
      if (is.null(dimnames(vcov_mat))) {
        dimnames(vcov_mat) <- list(coefficient_names, coefficient_names)
      }
      vcov_mat
    })
    return(result)
  }

  if (!is.list(fit) || is.null(fit$N) || is.null(fit$p)) {
    stop("fit must be the output of pgee_fit() or pgee().")
  }

  internal_names <- estimator_notation_from_name(estimator_names)
  result <- compute_requested_estimators(fit, internal_names)
  names(result) <- estimator_names
  result
}
