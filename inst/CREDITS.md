# Credits

## Vendored and adapted code provenance

### `binarySimCLF`-derived CLF utilities

The conditional linear family helper functions in `R/clf_generate.R`

- `clf_xch()`
- `clf_ar1()`
- `clf_cor2var()`
- `clf_allReg()`
- `clf_getBnds()`
- `clf_blrchk1()`
- `clf_blrchk()`
- `clf_mbsclf1()`
- `clf_mbsclf()`

are adapted from the archived `binarySimCLF` package implementation of the
Qaqish (2003) conditional linear family generator.

Relevant upstream sources:

- Qaqish BF. 2003. A family of multivariate binary distributions for simulating
  correlated binary variables with specified marginal means and correlations.
- By K and Qaqish BF. 2009. `binarySimCLF`, archived CRAN package.

### `geefirthr` / manuscript-derived PGEE fitting code

The PGEE fitting code in `R/pgee_fit.R` was extracted from the manuscript/thesis
workflow and descends from the `geefirthr` implementation associated with
Momenul Haque Mondol and M. Shafiqur Rahman, then further modified during the
`pgee-variance` project for package integration, initialization fallback, and
variance-estimation support. Public `geefirthr` package metadata reports
`License: GPL-2`.

## Licensing note

`pgeeVar` is licensed under `GPL (>= 2)`. This package keeps explicit
attribution for adapted GPL-compatible upstream code in source headers and in
`LICENSE.note`.
