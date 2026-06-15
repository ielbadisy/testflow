## R CMD check results

0 errors | 0 warnings | 2 notes

The notes are:

- New submission.
- unable to verify current time. This appears to be local environment-related
  during `R CMD check --as-cran`; no future timestamps were found in package
  files.

## test environments

- Local Ubuntu 24.04.3 LTS, R 4.5.1

## Release summary

This release improves package documentation and CRAN readiness. It adds a
dedicated README example for `sumtab()`, documents the implemented effect-size
formulas in a vignette, adds focused Cohen's d formula validation tests, and
removes generated vignette PDFs from the source repository.
