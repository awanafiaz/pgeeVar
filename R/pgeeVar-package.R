#' pgeeVar: Penalized GEE Variance Estimation Tools for Binary Longitudinal Data
#'
#' @description
#' `pgeeVar` provides the current numerical core for the `pgee-variance`
#' methods project. It supplies low-level simulation, fitting, and
#' variance-estimation tools for clustered binary longitudinal data.
#'
#' @details
#' The current public API is centered on three tasks:
#'
#' * generating correlated binary longitudinal data with
#'   [simulate_correlated_binary_data()]
#' * fitting penalized generalized estimating equation models with [pgee()]
#' * computing the supported small-sample covariance estimators with
#'   [compute_pgee_variances()]
#'
#' To inspect the supported estimator set, use [available_estimators()] or
#' [available_estimator_names()]. The package is intentionally focused on the
#' manuscript's 14-estimator variance-comparison workflow.
#'
#' @section Current package stage:
#' The package now exposes both user-facing wrappers and lower-level engines.
#' [simulate_correlated_binary_data()] and [pgee()] are the main entry points,
#' while [dgp_generate()] and [pgee_fit()] remain available for advanced users
#' who need the lower-level workflow.
#'
#' @name pgeeVar-package
#' @aliases pgeeVar
#' @keywords package
#' @importFrom stats binomial complete.cases nobs rnorm runif toeplitz
#' @importFrom utils globalVariables
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  globalVariables(".pgee_subset_internal")
}
