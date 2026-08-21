## Resubmission

This is a resubmission. The previously published CRAN version was 0.9.0.

Since 0.9.0, this release:

- Fixes two scientific-validity bugs: `sample_size_bioequivalence(method =
  "iterative_tost")` computed TOST power via a normal approximation instead
  of the exact noncentral-t formula (Phillips 1990), under-sizing studies by
  up to ~20% in some designs; `test_agreement()`'s test of `kappa = 0` reused
  the confidence-interval standard error instead of the correct
  null-hypothesis standard error (Fleiss 1981), inflating z and
  anti-conservative p-values. Both verified against reference implementations
  (`PowerTOST::power.TOST()`, `irr::kappa2()`).
- Adds rare-prevalence support to `sample_size_precision(endpoint =
  "binary")` (new Wilson/exact search methods, rare-event diagnostics), and
  fixes a bug where two-sample binary precision planning silently ignored
  unequal allocation.
- Closes remaining "not checked" assumption-diagnostic gaps (sphericity,
  symmetry-of-deviations, linearity, BH-adjusted correlation matrix); `car`
  moved from `Suggests` to `Imports` accordingly.
- Fixes an R-devel error in `test_repeated()`/`test_repeated_long()` (see
  prior resubmission note below, now folded into this release).
- Removes the unverified `sample_size_ordinal(method = "whitehead")`.
- Refocuses the DESCRIPTION Title/Description on statistical-test workflows
  (drops a stale `ggplot2`-visualization framing) and adds a new consolidated
  pedagogical reference vignette.
- Package version bumped to 1.0.0 to mark the first release with the full
  statistical-testing + sample-size-planning surface validated against
  reference implementations across the board (see NEWS.md for the complete
  per-function detail).

## R CMD check results

0 errors | 0 warnings | 2 notes

The notes are:

- unable to verify current time. This appears to be local environment-related
  during `R CMD check --as-cran`; no future timestamps were found in package
  files.
- "New maintainer" — the Maintainer field's family-name casing changed from
  "EL BADISY" to "El Badisy" (title case); same person, same email address,
  no change in maintainership.

`R CMD build --compact-vignettes=both` was used to compact the 5 vignette
PDFs (largest: 648Kb to 346Kb), which otherwise triggered a "checking sizes
of PDF files under 'inst/doc'" WARNING under `--as-cran`.

Also fixed a "no visible binding for global variable '.data'" NOTE across
6 functions (`descriptives_numeric`, `iqr_outliers`, `make_plot`,
`plot.sample_size`, `repeated_core`, `test_agreement`) by adding
`@importFrom rlang .data` via a new `R/testflow-package.R`.

## test environments

- Local Ubuntu 24.04.3 LTS, R 4.5.1
