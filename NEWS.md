# pgeeVar news

## pgeeVar 0.0.0.9000

Initial GitHub-oriented development release.

### Added

- `pgee()` as the main formula interface for penalized GEE fitting.
- `simulate_correlated_binary_data()` as the user-facing simulation wrapper.
- Fourteen covariance estimators via `compute_pgee_variances()` and `vcov(type = ...)`.
- `available_estimators()` and `available_estimator_names()` for estimator discovery.
- S3 methods for fitted `"pgee"` objects, including `summary()`, `vcov()`, `confint()`,
  `predict()`, `fitted()`, `residuals()`, `model.frame()`, `model.matrix()`,
  `formula()`, and `nobs()`.
- Exported help pages, a package overview page, runnable examples, and two vignettes.
- A `testthat` suite covering the current public package surface.

### Changed

- Canonical public estimator names are now the short forms (`AR`, `MD`, `Pan`, etc.).
  Legacy manuscript-style `V_*` names remain accepted as compatibility aliases.
- The estimator catalog now separates user-facing names from manuscript/internal
  notation through the `name` and `notation` columns.
- The archived `binarySimCLF` and `geefirthr` derivations retained in-package are
  now documented explicitly in package metadata and credits files.

### Notes

- This is still a development version intended for GitHub use.
- Citation metadata will be updated again once the companion manuscript has a
  public preprint identifier.
