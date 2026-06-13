testflow
================

# testflow

`testflow` is organized around study design, not test names. Instead of
asking “which function runs a t-test?”, testflow asks: “what is the
design of my study question?”

The core grammar is:

``` r
testflow = test + interpret + plot
```

``` r
library(testflow)

cardio <- make_cardio_data()

x <- cardio |>
  test_two_groups(sbp_3m ~ sex)

x
#> Statistical test workflow
#> 
#> Outcome: sbp_3m 
#> Group: sex 
#> Design: two independent groups 
#> 
#> Assumptions:
#> - Normality by group: acceptable
#> - Homogeneity of variance: acceptable
#> - F-test variance comparison: acceptable
#> 
#> Recommended test:
#> Student independent t-test 
#> 
#> Result:
#> H0: the population mean or location of sbp_3m is equal across levels of sex. 
#> statistic = -1.91, df = 178.00, p = 0.058, 95% CI [-11.22, 0.18] 
#> 
#> Effect size:
#> Cohen's d = -0.29, small
#> 
#> Report:
#> The two independent groups workflow for sbp_3m did not show a statistically significant result using Student independent t-test, statistic = -1.91, df = 178.00, p = 0.058. The 95% confidence interval was [-11.22, 0.18]. The effect size was small (Cohen's d = -0.29). H0: the population mean or location of sbp_3m is equal across levels of sex.
plot(x)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
report(x)
#> [1] "The two independent groups workflow for sbp_3m did not show a statistically significant result using Student independent t-test, statistic = -1.91, df = 178.00, p = 0.058. The 95% confidence interval was [-11.22, 0.18]. The effect size was small (Cohen's d = -0.29). H0: the population mean or location of sbp_3m is equal across levels of sex."
as_tibble(x)
#> # A tibble: 1 × 15
#>   workflow design outcome group recommended_test null_hypothesis statistic    df
#>   <chr>    <chr>  <chr>   <chr> <chr>            <chr>               <dbl> <dbl>
#> 1 two_gro… two i… sbp_3m  sex   Student indepen… H0: the popula…     -1.91   178
#> # ℹ 7 more variables: p <dbl>, conf.low <dbl>, conf.high <dbl>,
#> #   effect_size <dbl>, effect_size_name <chr>, effect_size_magnitude <chr>,
#> #   decision <chr>
```

The public workflows are named for common study designs:

``` r
test_one_sample(cardio, sbp_3m, mu = 140)
test_two_groups(sbp_3m ~ sex, data = cardio)
test_paired(sbp_3m ~ sbp_baseline, data = cardio)
test_groups(sbp_3m ~ treatment, data = cardio)
test_factorial(sbp_3m ~ sex * treatment, data = cardio)
test_categorical(treatment ~ controlled_3m, data = cardio)
test_correlation(sbp_3m ~ age, data = cardio)
test_outliers(cardio, c(sbp_3m, ldl, crp))
```

## Implemented workflows and tests

Each workflow returns a `testflow` object with the recommended test, H0,
p-value, confidence interval when available, appropriate effect size,
report text, and plot.

| Workflow                | Formula-oriented call                                                                     | Tests considered                                                                                |
|-------------------------|-------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| One sample              | `test_one_sample(cardio, sbp_3m, mu = 140)`                                               | one-sample t-test, Wilcoxon signed-rank, sign test                                              |
| Two independent groups  | `test_two_groups(sbp_3m ~ sex, data = cardio)`                                            | Student t-test, Welch t-test, Wilcoxon rank-sum                                                 |
| Paired measurements     | `test_paired(sbp_3m ~ sbp_baseline, data = cardio)`                                       | paired t-test, Wilcoxon signed-rank, sign test                                                  |
| More than two groups    | `test_groups(sbp_3m ~ treatment, data = cardio)`                                          | one-way ANOVA + Tukey, Welch ANOVA + Welch pairwise t-tests, Kruskal-Wallis + pairwise Wilcoxon |
| Factorial design        | `test_factorial(sbp_3m ~ sex * treatment, data = cardio)`                                 | factorial ANOVA with main effects and interactions                                              |
| Repeated measurements   | `test_repeated(cardio, c(sbp_baseline, sbp_3m, sbp_6m), id = id)`                         | repeated-measures ANOVA + paired t-tests, Friedman + paired Wilcoxon                            |
| Categorical association | `test_categorical(treatment ~ controlled_3m, data = cardio)`                              | chi-square independence test, Fisher exact test                                                 |
| Paired categorical      | `test_paired_categorical(cardio, controlled_baseline, controlled_3m)`                     | McNemar test                                                                                    |
| Repeated categorical    | `test_repeated_categorical(cardio, c(controlled_baseline, controlled_3m, controlled_6m))` | Cochran Q test + pairwise McNemar tests                                                         |
| One proportion          | `test_proportion(cardio, controlled_3m, success = "yes", p = 0.5)`                        | exact binomial test, one-sample proportion test                                                 |
| Multinomial             | `test_multinomial(cardio, treatment)`                                                     | chi-square goodness-of-fit, pairwise binomial checks                                            |
| Correlation             | `test_correlation(sbp_3m ~ age, data = cardio)`                                           | Pearson, Spearman, Kendall                                                                      |
| Correlation matrix      | `test_correlation_matrix(cardio, c(age, sbp_3m, ldl))`                                    | matrix of Pearson/Spearman/Kendall correlations                                                 |
| Outliers                | `test_outliers(cardio, c(sbp_3m, ldl, crp))`                                              | IQR outliers, Mahalanobis distance                                                              |

## References

- Fisher, R. A. (1925). .
- Gosset, W. S. (1908). The probable error of a mean.
- Welch, B. L. (1947). Generalization of Student’s problem with unequal
  variances.
- Wilcoxon, F. (1945). Individual comparisons by ranking methods.
- Mann, H. B., & Whitney, D. R. (1947). On a test of whether one of two
  random variables is stochastically larger than the other.
- Levene, H. (1960). Robust tests for equality of variances.
- Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion
  variance analysis.
- Tukey, J. W. (1949). Comparing individual means in the analysis of
  variance.
- Dunn, O. J. (1964). Multiple comparisons using rank sums.
- Friedman, M. (1937). The use of ranks to avoid the assumption of
  normality implicit in the analysis of variance.
- Cochran, W. G. (1950). The comparison of percentages in matched
  samples.
- McNemar, Q. (1947). Note on the sampling error of the difference
  between correlated proportions or percentages.
- Pearson, K. (1895, 1900).
- Spearman, C. (1904). The proof and measurement of association between
  two things.
- Kendall, M. G. (1938). A new measure of rank correlation.
- Cramer, H. (1946). .
- Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or
  fiducial limits illustrated in the case of the binomial.
- Cohen, J. (1988). .
