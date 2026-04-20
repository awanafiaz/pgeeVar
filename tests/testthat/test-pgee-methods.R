test_that("pgee methods return the expected structures and defaults", {
  fit <- make_balanced_fit()
  summary_fit <- summary(fit)
  vcov_default <- vcov(fit)
  vcov_md <- vcov(fit, type = "MD")
  vcov_null <- vcov(fit, type = NULL)
  ci_95 <- confint(fit)
  ci_99 <- confint(fit, level = 0.99)
  ci_kc <- confint(fit, parm = "X1", type = "KC")
  ci_vec <- confint(fit, parm = c("X1", "obstime"))
  pred_train <- predict(fit)
  pred_link <- predict(fit, type = "link")
  pred_train_rows <- predict(fit, newdata = model.frame(fit)[1:5, , drop = FALSE])
  pred_new <- predict(fit, newdata = make_balanced_sim()[1:5, , drop = FALSE])
  pred_se <- predict(fit, newdata = make_balanced_sim()[1:5, , drop = FALSE], se.fit = TRUE)
  pred_se_legacy <- predict(
    fit,
    newdata = make_balanced_sim()[1:5, , drop = FALSE],
    se.fit = TRUE,
    vcov_type = "V_AR"
  )

  expect_s3_class(fit, "pgee")
  expect_equal(colnames(model.matrix(fit)), c("(Intercept)", "X1", "obstime"))
  expect_equal(length(fitted(fit)), nobs(fit))
  expect_equal(length(residuals(fit, type = "response")), nobs(fit))
  expect_equal(length(residuals(fit, type = "pearson")), nobs(fit))
  expect_equal(nrow(model.frame(fit)), nobs(fit))
  expect_equal(nrow(model.matrix(fit)), nobs(fit))
  expect_equal(deparse(formula(fit)), deparse(y ~ X1 + obstime))

  expect_s3_class(summary_fit, "summary.pgee")
  expect_equal(
    colnames(summary_fit$coefficients),
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  )
  expect_equal(rownames(vcov_default), names(coef(fit)))
  expect_equal(colnames(vcov_default), names(coef(fit)))
  expect_equal(vcov_default, fit$vcovs[[fit$default_vcov]], tolerance = 1e-10)
  expect_equal(vcov_null, vcov_default, tolerance = 1e-10)
  expect_equal(vcov_md, fit$vcovs$MD, tolerance = 1e-10)
  expect_equal(vcov(fit, type = "V_MD"), vcov_md, tolerance = 1e-10)
  expect_true(is.finite(sqrt(diag(vcov_default))["X1"]))
  expect_error(vcov(fit, type = "V_ZZ"), "available_estimator_names")

  expect_equal(pred_train, fitted(fit), tolerance = 1e-10)
  expect_equal(pred_link, fit$linear_predictors, tolerance = 1e-10)
  expect_equal(pred_train_rows, fitted(fit)[1:5], tolerance = 1e-10)
  expect_length(pred_new, 5L)
  expect_identical(names(pred_se), c("fit", "se.fit", "type", "vcov_type"))
  expect_equal(pred_se$vcov_type, "AR")
  expect_equal(pred_se_legacy$vcov_type, "AR")
  expect_equal(dim(ci_kc), c(1L, 2L))
  expect_equal(dim(ci_vec), c(2L, 2L))
  expect_true(all((ci_99[, 2] - ci_99[, 1]) > (ci_95[, 2] - ci_95[, 1])))

  fit_ind <- make_balanced_fit(corstr = "independence")
  expect_equal(fit_ind$alpha_hat, 0)

  print_out <- capture.output(print(fit))
  summary_out <- capture.output(print(summary_fit))
  expect_true(any(grepl("Working correlation", print_out, fixed = TRUE)))
  expect_true(any(grepl("Default variance", print_out, fixed = TRUE)))
  expect_true(any(grepl("Residual df", summary_out, fixed = TRUE)))
})

test_that("pooling variance methods warn exactly once on unbalanced fits", {
  fit_unbalanced <- make_unbalanced_fit()
  warning_count <- 0L
  expected_dim <- rep(length(coef(fit_unbalanced)), 2)

  v_pan <- withCallingHandlers(
    vcov(fit_unbalanced, type = "Pan"),
    warning = function(w) {
      warning_count <<- warning_count + 1L
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(warning_count, 1L)
  expect_all_na_matrix(v_pan, expected_dim = expected_dim)

  warning_count_predict <- 0L
  pred_pan <- withCallingHandlers(
    predict(
      fit_unbalanced,
      newdata = model.frame(fit_unbalanced)[1:5, , drop = FALSE],
      se.fit = TRUE,
      vcov_type = "V_Pan"
    ),
    warning = function(w) {
      warning_count_predict <<- warning_count_predict + 1L
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(warning_count_predict, 1L)
  expect_equal(pred_pan$vcov_type, "Pan")
  expect_true(all(is.na(pred_pan$se.fit)))
  expect_length(pred_pan$fit, 5L)
})
