## Update submission (1.0.0 -> 1.0.5)

Bug-fix releases; see NEWS.md for details. Main changes:

* `sample_size_precision()` (binary, `method = "wilson"` or `"exact"`):
  the returned n is now the smallest n from which the precision target
  holds for every larger n (achieved precision is not monotone in n with
  integer counts).
* `test_groups()` and other workflows no longer fail when an integer
  outcome has groups of odd and even size.
* `test_linear_regression()`, `test_logistic_regression()` and
  `test_factorial()` no longer fail for n > 5000 (Shapiro-Wilk limit).
* Dropped the 'dplyr', 'tidyr' and 'purrr' dependencies.
* `test_categorical()` falls back to a Monte Carlo p-value when Fisher's
  exact test runs out of workspace.

## Test environments

* Local: Ubuntu 24.04, R 4.5.1

## R CMD check results

0 errors | 0 warnings | 0 notes

The tarball was built with `R CMD build --compact-vignettes=both` to
compact the vignette PDFs.

## Reverse dependencies

None on CRAN.
