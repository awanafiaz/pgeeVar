# pgeeVar

`pgeeVar` fits penalized generalized estimating equations for correlated binary
data and computes the finite-sample covariance estimators studied in the
`pgee-variance` project. The package is built around three practical tasks:

- simulate correlated binary longitudinal data;
- fit a penalized GEE model with a formula interface;
- compare one or more covariance estimators from the package catalog.

Most users will work with `simulate_correlated_binary_data()` and `pgee()`.
The lower-level `dgp_generate()` and `pgee_fit()` functions remain available for
simulation scripts that need direct access to the generated object or the raw
fit structure.

## Install from GitHub

```r
install.packages("pak")
pak::pak("awanafiaz/pgeeVar")
```

## Quick start

The example below simulates one small dataset, fits a penalized GEE model, and
extracts a few covariance estimators. The commented lines show the shape of the
returned objects for this fixed-seed example.

```r
library(pgeeVar)

example_data <- simulate_correlated_binary_data(
  N = 20,
  gamma = 0.5,
  n_i = 4,
  beta = c(-1.2, 0.6, 0.3),
  rho = 0.2,
  corstr = "exchangeable",
  seed = 1
)

fitted_model <- pgee(
  y ~ X1 + obstime,
  data = example_data,
  id = id,
  corstr = "exchangeable"
)

head(example_data, 4)
#>   id y X1 obstime
#> 1  1 1  1     0.2
#> 2  1 1  1     0.4
#> 3  1 1  1     0.6
#> 4  1 1  1     0.8

round(summary(fitted_model)$coefficients, 3)
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)   -1.388      0.941  -1.475    0.159
#> X1            -0.010      0.898  -0.012    0.991
#> obstime        0.547      0.965   0.567    0.578

round(sqrt(diag(vcov(fitted_model, type = "AR"))), 3)
#> (Intercept)          X1     obstime
#>       0.941       0.898       0.965

names(compute_pgee_variances(fitted_model, c("AR", "MBN")))
#> [1] "AR"  "MBN"
```

The default summary reports the coefficient table, the working-correlation
estimate, convergence information, and the default covariance choice. Use
`vcov()` when you want one estimator at a time, or
`compute_pgee_variances()` when you want a named subset or the full catalog.

## Estimator names

The public interface uses short estimator names such as `AR`, `MD`, `MBN`,
`Pan`, and `FZ`. The manuscript notation (`V_AR`, `V_MD`, and so on) is kept in
the catalog returned by `available_estimators()`. Legacy `V_*` names are still
accepted as aliases.

```r
available_estimators()[, c("name", "notation", "class", "requires_balance")]
available_estimator_names()
```

## Vignettes

The package includes two short vignettes:

- `analysis-workflow`: a first-pass walkthrough from simulation to fitting,
  covariance extraction, prediction, and low-level access;
- `simulation-design-workflow`: a compact simulation pattern for repeated fits,
  balanced versus unbalanced designs, and low-level scripting.

## Vendored CLF Credit

The `clf_*` functions in [`R/clf_generate.R`](R/clf_generate.R) are adapted from
the archived `binarySimCLF` package implementation of the conditional linear
family generator described by Qaqish (2003) and By and Qaqish (2009). Those
functions are retained in-package with attribution because the original package
is no longer available on CRAN.

See [`inst/CREDITS.md`](inst/CREDITS.md) and [`LICENSE.note`](LICENSE.note) for
the current provenance and licensing notes.
