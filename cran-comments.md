## Update submission (1.0.5 -> 1.1.0)

* Dropped the 'tibble' and 'tidyselect' dependencies; result tables are
  built with 'basetable' (already imported).
* New `as.data.frame()` methods for testflow and sample-size results;
  `as_tibble()` keeps working.

## Test environments

* Local: Ubuntu 24.04, R 4.5.1

## R CMD check results

0 errors | 0 warnings | 0 notes

The tarball was built with `R CMD build --compact-vignettes=both` to
compact the vignette PDFs.

## Reverse dependencies

None on CRAN.
