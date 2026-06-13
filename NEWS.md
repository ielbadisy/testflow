# testflow 0.8.1

## Summary tables

- Added `sumtab()` for formula-driven descriptive summary tables with optional
  automatic p-value selection.

## Validation and documentation

- Added canonical references to workflow documentation and the main vignette.
- Added deeper validation checks against base R/reference calculations for
  one-sample tests, proportions, correlations, categorical association, and
  factorial ANOVA.
- Kept the hardened repeated-measures and post-hoc workflows from `0.8.0`.

# testflow 0.8.0

## Method hardening

- Replaced the repeated-categorical Friedman approximation with an explicit
  Cochran Q implementation for binary repeated measures.
- Added pairwise McNemar post-hoc comparisons for repeated categorical
  workflows.
- Added base R repeated-measures ANOVA extraction for normal repeated numeric
  workflows.
- Added paired post-hoc comparisons for repeated numeric workflows:
  paired t-tests after repeated-measures ANOVA and paired Wilcoxon tests after
  Friedman.
- Added method-specific post-hoc selection for multi-group comparisons:
  Tukey HSD after classical ANOVA, Welch pairwise t-tests after Welch ANOVA,
  and pairwise Wilcoxon tests after Kruskal-Wallis.
- Added validation tests comparing repeated ANOVA and Cochran Q results against
  explicit base R/reference formulas.

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
