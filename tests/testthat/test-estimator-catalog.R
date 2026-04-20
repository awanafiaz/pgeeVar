test_that("available estimator catalog matches the fixed paper set", {
  catalog <- available_estimators()
  expected_names <- c(
    "LZ", "DF", "KC", "MD", "FG", "MBN",
    "Pan", "GST", "WL", "WB", "RS",
    "FW", "FZ", "AR"
  )
  expected_notation <- c(
    "V_LZ", "V_DF", "V_KC", "V_MD", "V_FG", "V_MBN",
    "V_Pan", "V_GST", "V_WL", "V_WB", "V_RS",
    "V_FW", "V_FZ", "V_AR"
  )
  expected_pooled <- c("Pan", "GST", "WL", "WB", "RS")

  expect_s3_class(catalog, "data.frame")
  expect_equal(nrow(catalog), 14L)
  expect_equal(
    names(catalog),
    c(
      "name", "notation", "class", "requires_balance", "reference_key",
      "default_order", "description"
    )
  )
  expect_equal(catalog$name, expected_names)
  expect_equal(catalog$notation, expected_notation)
  expect_equal(available_estimator_names(), expected_names)
  expect_equal(catalog$default_order, seq_len(14))
  expect_equal(catalog$name[catalog$class == "pooled"], expected_pooled)
  expect_false(anyNA(catalog$reference_key))
  expect_equal(length(unique(catalog$reference_key)), nrow(catalog))
  expect_setequal(
    pgeeVar:::pooling_estimator_names(),
    catalog$name[catalog$class == "pooled"]
  )
})
