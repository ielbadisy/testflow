# testflow (development version)

## Assumption checking

Closed remaining gaps where an assumption was labeled `"not checked"`
because no diagnostic had been implemented (as opposed to being genuinely
uncheckable, e.g. too few observations, which is unchanged). Every
assumption status now reflects either a real check or an honest "not
applicable"/"not enough data" label - never an unimplemented one.

- `car` is now a hard dependency (`Imports`, not `Suggests`): regression
  homoscedasticity (Breusch-Pagan), multicollinearity (VIF), and
  ANOVA/factorial variance homogeneity (Levene's test) no longer fall back
  to `"not checked"` when `car` isn't installed, and `test_factorial(type
  = 2/3)` no longer errors for the same reason.
- Added a real sphericity check for `test_repeated()`/`test_repeated_long()`
  with 3+ conditions, via Mauchly's test (`stats::mauchly.test()` on the
  wide-format multivariate model, `X = ~1`) - verified to match
  `mauchly.test()` and `rstatix::anova_test()` exactly. Reports "not
  applicable" for exactly two conditions, where sphericity isn't defined.
- Added a real symmetry-of-deviations check for `test_one_sample()` and
  `test_paired()` (used to justify the Wilcoxon signed-rank test), via the
  closed-form Cabilio & Masaro (1996) test - verified to match
  `lawstat::symmetry.test(option = "CM", boot = FALSE)` exactly.
- Added a real linearity check for `test_correlation(method = "pearson")`,
  via a quadratic-term F-test (`anova()` comparing `y ~ x` to `y ~ x +
  I(x^2)`) - equivalent to Ramsey's RESET test with a single quadratic
  regressor term, verified to match `lmtest::resettest(power = 2, type =
  "regressor")` exactly.
- `test_correlation_matrix()` now applies a Benjamini-Hochberg correction
  to its pairwise p-values by default (`p.adj` column in the correlation
  table), instead of reporting pairwise p-values with no correction
  applied.

## Sample size planning

- Removed `sample_size_ordinal(method = "whitehead")`. Unlike every other
  formula in the sample-size module, the proportional-odds formula had not
  been checked against a published worked example or an external
  reference implementation (only against itself), and the input
  convention for `probs` (which arm's category distribution to supply) is
  ambiguous across secondary sources describing Whitehead (1993).
  `sample_size_ordinal()` now implements only Noether's method, which was
  independently verified. It no longer takes `method`, `probs`, or
  `odds_ratio` arguments.

# testflow 0.9.0

## Display ergonomics: per-term results tables in console output

Second of two planned display-ergonomics changes (the first added
workflow-specific print titles and field labels). Regression coefficients,
Cox hazard ratios, the diagnostic-accuracy table, the ROC/Youden threshold,
the ICC(1)/(2)/(3) comparison table, correlation-method/pairwise-correlation
tables, factorial ANOVA tables, and every tidy post-hoc result (repeated
measures, repeated categorical, multinomial, paired categorical) were
computed but never shown in `print()`/`summary()` console output - only
reachable via `x$alternative_tests`/`x$posthoc`. They now print directly.

- `R/workflow-display.R`'s registry gained a `table`/`table_label` field
  per applicable workflow, resolved by `workflow_meta()` into a small
  data frame when the referenced field is already tidy (a `data.frame`/
  `tibble`) or a base `table`/`matrix` (auto-converted, with generic
  `Var1`/`Var2`/`Freq` column names relabeled using the same outcome/group
  labels already resolved for that workflow, rather than left meaningless).
- New `tf_table()` helper in `R/print.R` prints these tables, applying the
  same `format_stat()`/`format_p()` conventions already used elsewhere
  (p-value columns get "<0.001" instead of a misleading truncated "0.00";
  integer count columns are formatted without decimal places).
- `groups`/`repeated`'s omnibus-ANOVA post-hoc (`posthoc_groups()`,
  `R/posthoc.R`) wraps a raw `TukeyHSD`/`pairwise.htest` object rather than
  a tidy data frame; it's silently skipped rather than crashing or dumping
  a raw object, matching prior behavior. (`repeated`'s own separate
  paired-comparison post-hoc, from `paired_posthoc_numeric()`, *is* already
  tidy and renders normally - only `groups`' omnibus-ANOVA path is
  affected.) Formatting an arbitrary htest/TukeyHSD object nicely is a
  larger follow-up, not attempted here.
- Fixed a bug caught while verifying wide tables (e.g. the ICC comparison
  table, which combines a long text column with several numeric ones):
  base `data.frame` printing wraps onto a second, visually misaligned block
  of rows once the default 80-character console width is exceeded.
  `tf_table()` now temporarily widens `options(width = ...)` for the
  duration of the print call rather than truncating labels or letting rows
  wrap. The first version of this width calculation didn't account for `NA`
  values (which produce an `NA` width via plain `max()`/`sum()`, in turn
  making `options(width = NA)` error) - caught by the existing assumption-
  structure test suite hitting a real `car::Anova()` aliased-coefficients
  edge case with an all-`NA` ANOVA row, not by a synthetic test written
  after the fact.

Testing: devtools::test() 409/409 passing, including new assertions
pinning table rendering for six representative workflows (regression, Cox,
diagnostic, ICC's no-wrap behavior, agreement's relabeled columns, and the
groups-skips/repeated-renders post-hoc distinction). devtools::check(): 0
errors, 0 warnings, 1 environment-only NOTE, including vignette rebuild
(console-output examples updated automatically).

## Display ergonomics: workflow-specific print titles and field labels

- `print.testflow()` and `print.summary.testflow()` now show a title
  specific to the workflow that actually ran (e.g. "Linear Regression",
  "Cox Proportional Hazards Regression", "Intraclass Correlation") instead
  of the same generic "Statistical test workflow" for every single
  `test_*()` result. The now-redundant separate `Design:` field (which
  always duplicated the same information less legibly) is dropped.
- The `Outcome:`/`Group:` console field labels are now workflow-appropriate
  where `Group:` previously held something other than a comparison group:
  "Predictors" for `test_linear_regression()`/`test_logistic_regression()`/
  `test_cox()`, "Reference" for `test_diagnostic()`, "Rater 1"/"Rater 2" for
  `test_agreement()`, "Time" in place of "Outcome" for the two survival
  workflows, and so on.
- New internal `R/workflow-display.R` centralizes this per-workflow display
  metadata in one registry keyed by `x$workflow`, rather than threading new
  parameters through every `new_testflow()` call site across `R/test_*.R`
  (none of those files changed). `sample_size`-class objects are
  unaffected: their generic "Sample size planning" title already fits every
  call, since they share one purpose unlike the widely varying `test_*()`
  workflows.
- This is the first of two planned display-ergonomics changes; a follow-up
  will surface per-term/per-metric result tables (regression coefficients,
  Cox hazard ratios, the diagnostic-accuracy table, the ROC/Youden
  threshold, the ICC comparison table, and post-hoc results) directly in
  console output, which are currently computed but only reachable via
  `x$alternative_tests`/`x$posthoc`.

## New sample-size functions

- Added `sample_size_bioequivalence()`: two one-sided test (TOST)
  bioequivalence planning on the log scale, for `crossover` or `parallel`
  designs. `method = "iterative_tost"` (default) searches via
  `stats::uniroot()` for the smallest n achieving the exact TOST power,
  verified to recover the documented closed-form special case at GMR=1
  exactly and to match an independently computed `uniroot()` search
  off-center; `method = "normal_approx"` offers the closed-form
  approximation for comparison.
- Added `sample_size_precision()`: sample size for a target confidence-
  interval half-width instead of power against an effect size, covering one-
  and two-sample continuous and binary designs plus log-odds-ratio
  precision.
- Added `sample_size_cluster_adjust()`: design-effect adjustment
  (`DE = 1+(m-1)rho`, plus an unequal-cluster-size approximation) for
  converting an individually randomized sample size to a cluster-randomized
  one. Like `sample_size_adjust_dropout()`, it returns a plain integer, not
  a `sample_size` object.
- `sample_size_ordinal()` gained `method = c("noether", "whitehead")`;
  Whitehead's proportional-odds formula is included with an explicit
  caveat (in both the documentation and the returned assumptions) that,
  unlike every other formula in this release, it has not been validated
  against a published worked example and is sensitive to the anticipated
  category-probability distribution.
- `sample_size_survival()` gained optional `accrual_duration`/`follow_up`
  arguments: when supplied with `survival_a`/`survival_b`, event
  probabilities use a closed-form uniform-accrual adjustment instead of
  assuming every subject is followed for the full study duration. Verified
  to reduce to the existing flat conversion as `accrual_duration -> 0`, and
  to always require more subjects than instantaneous accrual for the same
  event target.
- This is Phase 4 (the last currently planned phase) of the teaching/
  clinical-research expansion that began with regression, survival
  analysis, and diagnostic/agreement statistics. A Bayesian assurance
  function (`ss_assurance` in the original formula spec) remains explicitly
  out of scope: it is simulation-based rather than closed-form and was
  flagged in the spec itself as "a later module, not v1."

## New workflows: diagnostic and agreement statistics

- Added `test_diagnostic()`: sensitivity, specificity, positive/negative
  predictive values, and positive/negative likelihood ratios from a 2x2
  table against a gold-standard reference, each with a confidence interval
  (exact/Clopper-Pearson for the proportions, the closed-form log-scale
  interval of Simel, Samsa & Matchar (1991) for the likelihood ratios). The
  primary test compares overall accuracy to the no-information rate, the
  same convention `caret::confusionMatrix()` uses.
- Added `test_roc()`: AUC computed from the Mann-Whitney U statistic (the
  same identity already used for the package's rank-biserial correlation),
  with the closed-form Hanley & McNeil (1982) confidence interval and the
  Youden's-J-optimal threshold, avoiding a `pROC` dependency.
- Added `test_agreement()`: Cohen's kappa with the Fleiss, Cohen & Everitt
  (1969) large-sample standard error (verified to match
  `psych::cohen.kappa()`'s confidence interval exactly), not the simpler
  approximation that understates variance under uneven category marginals.
- Added `test_icc()`: ICC(1,1), ICC(2,1), and ICC(3,1) via Shrout & Fleiss
  (1979) one-way/two-way ANOVA variance decomposition, with McGraw & Wong
  (1996) F-distribution-based confidence intervals (verified to match
  `irr::icc()` exactly for all three forms). ICC(2,1) (two-way random,
  absolute agreement) is reported as the primary/effect-size estimate,
  following the reliability-study recommendation of Koo & Li (2016).
- None of these four functions needed a new package dependency; every
  formula is implemented directly against a published closed-form reference
  and cross-checked against `psych`/`irr` during development (not shipped
  as dependencies).
- AUC, Cohen's kappa, and ICC get real, citable magnitude conventions
  (Hosmer/Lemeshow/Sturdivant 2013; Landis & Koch 1977; Koo & Li 2016
  respectively) - unlike the hazard ratio and concordance index added in the
  survival-analysis release, where no such convention exists in the
  literature and none was invented.
- Fixed a bug caught while building `test_icc()`: the ICC(2,1) F-test
  (testing whether ICC exceeds 0) was accidentally using the same
  Satterthwaite coefficients as the confidence interval, which are
  evaluated at different points (the null value 0 for the test, the
  estimated ICC for the interval) and must not be shared - this silently
  produced a wildly wrong p-value (0.48 instead of <0.001) for data with
  strong measured reliability.
- Fixed a second, general bug this surfaced: `safe_tidy_htest()` (used by
  most `test_*()` workflows) left stray name attributes on numeric columns
  inherited from the underlying htest object (e.g. `binom.test()$statistic`
  carrying the name `"number of successes"`), which could trip up
  `all.equal()`-based comparisons against independently computed reference
  values despite the numbers being identical.
- This is Phase 3 of a broader expansion (additional sample-size functions
  are planned next).

## New workflows: survival analysis

- Added `test_survival()`: Kaplan-Meier estimation and the log-rank test
  (`survival::survfit()`/`survival::survdiff()`) for a two-group time-to-
  event comparison, with a companion univariate Cox hazard ratio (and its
  confidence interval) as the effect size.
- Added `test_cox()`: Cox proportional hazards regression
  (`survival::coxph()`), with the overall likelihood-ratio test as the
  primary result, per-term hazard ratios, a proportional-hazards assumption
  check via the Schoenfeld residual test (`survival::cox.zph()`), and the
  concordance index as the effect size.
- `survival` is now a hard dependency (`Imports`, not `Suggests`): it ships
  with every standard R installation as a "Recommended" package, so
  `requireNamespace()`-guarding it added complexity without meaningfully
  widening compatibility. `Surv()` is re-exported from `testflow` so it's
  available after `library(testflow)` alone.
- Neither the hazard ratio nor the concordance index is assigned a
  magnitude label (negligible/small/moderate/large): unlike Cohen's d or
  R-squared, there is no widely agreed convention for either, and inventing
  one would overstate how standardized that judgment is. The numeric
  estimate is always reported.
- Fixed a related pre-existing bug this surfaced: `print.testflow()` and
  `make_report()` hid an effect size's estimate entirely whenever its
  magnitude label was `NA`, even when the estimate itself was available.
  Both now show the estimate unconditionally and append the magnitude label
  only when one is assigned.
- This is Phase 2 of a broader expansion (diagnostic/agreement statistics
  and additional sample-size functions are planned next).

## New workflows: regression

- Added `test_linear_regression()`: multiple linear regression via
  `stats::lm()`, with residual normality (Shapiro-Wilk), homoscedasticity
  (Breusch-Pagan), and multicollinearity (variance inflation factor)
  assumption checks, the overall model F-test as the primary result, the
  per-term coefficient table, and R squared / adjusted R squared as the
  effect size.
- Added `test_logistic_regression()`: logistic regression via
  `stats::glm(family = binomial)`, with multicollinearity and influential-
  observation (Cook's distance) assumption checks, the likelihood-ratio test
  against the intercept-only model as the primary result, coefficients on
  both the log-odds and odds-ratio scale, and McFadden's pseudo R squared as
  the effect size.
- Both are documented in `vignettes/statistical-test-workflows.Rmd` and
  `vignettes/effect-size-formulas.Rmd`, and every reported statistic is
  pinned in tests against `summary(lm())`/`summary(glm())` computed
  independently.
- This is Phase 1 of a broader expansion (survival analysis, diagnostic and
  agreement statistics, and additional sample-size functions are planned
  next).

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
