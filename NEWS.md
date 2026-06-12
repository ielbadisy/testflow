# testflow 0.7.0

## Workflow API

- Added formula-first support for core study designs:
  - `test_two_groups(outcome ~ group, data = data)`
  - `test_groups(outcome ~ group, data = data)`
  - `test_categorical(x ~ y, data = data)`
  - `test_correlation(y ~ x, data = data)`
  - `test_paired(after ~ before, data = data)`
  - `test_factorial(outcome ~ factor1 * factor2, data = data)`
- Preserved data-first and pipe-friendly calls for existing examples.

## Results

- Primary test results now expose the null hypothesis, confidence interval
  columns when available, and effect size fields.
- `as_tibble()` includes H0, confidence interval, effect size name, effect size
  estimate, effect size magnitude, and decision.

## Documentation

- Added rendered HTML vignettes in `inst/doc`.
- Expanded the README with a list of implemented workflows and tests considered.

## Quality

- Added focused tests for formula-first calls, H0, confidence intervals, and
  effect size exposure.
- Added GitHub Actions R CMD check workflow.
