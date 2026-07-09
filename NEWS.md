# testflow 0.8.2

## Statistical test correctness

An exhaustive correctness review of every `test_*()` workflow (run against
base-R references and the package's own documented formulas) found and
fixed the following:

- Fixed `test_two_groups()` reporting Cohen's d (and rank-biserial
  correlation) with the opposite sign from the primary t-test/Wilcoxon test
  it's attached to, whenever group labels weren't already alphabetically
  ordered in the source data. The shared `assert_two_groups()` helper now
  derives group order from `levels(as.factor(...))`, matching what
  `t.test`/`wilcox.test`'s formula interface uses internally, instead of
  first-appearance order in the raw data.
- Fixed `test_repeated()` computing the wrong F-statistic, error degrees of
  freedom, and p-value for the parametric (repeated-measures ANOVA) branch:
  the subject `id` column was never coerced to a factor before
  `aov(y ~ time + Error(id/time))`, so `Error()` built a spurious extra
  stratum instead of the correct `id:time` error term. This affected the
  **default path** (auto-generated subject ids are always integer), and the
  test suite's own reference computation reproduced the same bug rather than
  catching it.
- Fixed `test_categorical()` recommending Fisher's exact test while
  simultaneously showing an "Expected cell counts: acceptable" assumption
  panel for the same table, because the recommendation used an
  any-cell-below-threshold rule while the displayed panel used a different
  (Cochran 80%-of-cells) rule. Both now use the same `fisher_threshold`
  rule, so the panel explains the recommendation actually made.
- Fixed `test_correlation_matrix()` silently doing listwise deletion across
  *all* selected variables before computing pairwise correlations, which
  made `use = "pairwise.complete.obs"` a no-op, discarded usable rows, and
  contradicted its own "Pairwise complete observations are used" assumption
  message. Each pairwise correlation and test now uses its own
  pairwise-complete rows.
- Wired up `test_factorial()`'s `type` parameter, which defaulted to `2` but
  was documented as "a placeholder for future car integration" and always
  computed Type I (sequential) sums of squares regardless of its value. For
  unbalanced factorial designs, Type I vs. Type II vs. Type III sums of
  squares can give substantially different p-values; `type = 2`/`3` now
  dispatch to `car::Anova()` (Type III with the required sum-to-zero
  contrasts), and the reported eta squared is computed from the same sums of
  squares as the reported test rather than always from Type I.
- Fixed `test_multinomial()`'s pairwise post-hoc binomial tests not applying
  a multiple-comparison correction, unlike every other post-hoc helper in
  the package (BH by default). `p.adj`/`p.adjust.method` columns are now
  included.
- Relabeled `posthoc_groups()`'s Welch/Kruskal-Wallis follow-up tests, which
  were called "Games-Howell-style" and "Dunn-style" but are plain
  BH-adjusted pairwise Welch t-tests and pairwise Wilcoxon rank-sum tests
  respectively, not the named procedures (which use a studentized-range or
  pooled-rank-variance correction). Labels now describe what's actually run.
- Fixed `test_one_sample()` and `test_correlation()` labeling unchecked
  assumptions (Wilcoxon symmetry, Pearson linearity) as `"acceptable"` with
  no diagnostic actually performed. Both now use the honest `"warning"`/
  `"not checked"` status already used for the identical situation elsewhere
  in the package (e.g. `test_paired()`'s symmetry check).

## Sample size planning

- Fixed `sample_size_survival()` understating total sample size by a factor
  of 2 whenever `survival_a`/`survival_b` were supplied: the events-to-N
  conversion now correctly applies `N_total = 2 * D_total / (qA + qB)`
  instead of `D_total / (qA + qB)`. This affected superiority,
  non-inferiority, and equivalence planning alike.
- Fixed the sample-size power curve for paired and repeated-measures
  continuous superiority planning, which used a two-sample (parallel-group)
  power formula on a one-sample (paired-difference) result and understated
  plotted power at the target sample size.
- Fixed `sample_size_binary(objective = "equivalence")` silently ignoring
  the `allocation` ratio when `p1 == p2`, which produced incorrect per-group
  sizes for unequal allocation designs.
- Fixed `plot(x, type = "curve")` and `plot(x, type = "both")` erroring for
  every `sample_size` object with curve data: a local tibble column named
  `x` shadowed the `x` (sample-size object) argument during construction.
- Moved sample-size formulas and literature references out of the values
  returned by `sample_size()` and its helpers and into the function
  documentation (`?sample_size_continuous`, etc.), consistent with normal R
  package conventions. `summary()`, `print()`, and `as_tibble()` output no
  longer carries a `formula`/`reference` field.
- Added a "Sample size planning" vignette walking through the supported
  endpoints, designs, and objectives with worked examples.
- Added a "Sample size: complete API reference and formulas" vignette
  covering every `design`/`objective`/`method` combination across
  `sample_size_continuous()`, `sample_size_binary()`,
  `sample_size_survival()`, and `sample_size_ordinal()`, each paired with
  its formula and a runnable example.

## Outlier workflow

- Standardized `test_outliers()` on the same `test_*(formula, data, ...)`
  calling convention used by the other workflow helpers.
- Kept outlier visualization on the generic `plot()` S3 path via
  `plot.testflow()`.
- Hardened assumption printing so workflows with extra diagnostic columns no
  longer break `print()`.

## Output

- Structured printed `testflow` and summary output with `cli`-based colors,
  section headings, field labels, and assumption bullets.
- Interactive console colors are enabled by default for printed `testflow`
  objects and summaries. Set `options(testflow.cli_colors = FALSE)` to disable
  them, or `options(testflow.cli_colors = TRUE)` to force colors in scripts.

## Documentation and release prep

- Added a dedicated `sumtab()` example to the README.
- Added an effect-size formulas vignette documenting the implemented formulas
  for Cohen's d, eta squared, Cramer's V, rank-biserial correlation, Kendall's
  W, and related workflow summaries.
- Added validation tests that pin the documented Cohen's d formulas to the
  implementation.
- Removed generated vignette PDFs from the source repository.

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

First public release.

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
