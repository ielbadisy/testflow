## Resubmission

This release fixes the r-devel errors CRAN reported against the currently
published version (0.8.2) on 2026-07-16:

  Running 'testthat.R' [...] failed / Error in
  `if (is.finite(resvar) && resvar < (mean(fitted)^2 + var(c(fitted))) *
  1e-30) perfect.fit[y] <- TRUE`: missing value where TRUE/FALSE needed

The error originated in `repeated_anova_test()`, which called `summary()`
on an `aovlist` fit (`aov(... + Error(id/within))`). That walks every
`Error()` stratum, including the degenerate intercept-only stratum, which
has a single fitted value. A recent r-devel change to
`stats:::summary.aov()`'s "perfect fit" check computes `var()` on that
stratum's fitted values, which is `NA` for a length-1 vector, turning the
guard's result into `NA` instead of `TRUE`/`FALSE`. This broke 5 tests
(`test-assumptions-repeated.R`, `test-assumptions-structure.R`,
`test-print.R`, `test-test-repeated.R`) and re-building the
`statistical-test-workflows.Rmd` vignette on r-devel (Debian, Fedora,
Windows), though not on release R.

`repeated_anova_test()` (R/repeated-methods.R) now summarizes each
`Error()` stratum independently and skips any that fail to summarize,
since the intercept-only stratum was never used for the reported test.
Verified with 421 passing tests locally (release R doesn't reproduce the
bug, so this can't be confirmed against r-devel until the next CRAN
pretest).

This version also folds in unrelated work already queued on the
development branch before the CRAN report arrived: closed remaining
`"not checked"` assumption gaps (sphericity, symmetry-of-deviations,
linearity, BH-adjusted correlation-matrix p-values; `car` moved from
`Suggests` to `Imports`), and removed the unverified Whitehead ordinal
sample-size method (`sample_size_ordinal(method = "whitehead")`) since its
formula couldn't be checked against a published worked example.

## R CMD check results

0 errors | 0 warnings | 1 note

The note is:

- unable to verify current time. This appears to be local environment-related
  during `R CMD check --as-cran`; no future timestamps were found in package
  files.

## test environments

- Local Ubuntu 24.04.3 LTS, R 4.5.1
