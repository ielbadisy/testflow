# testflow

`testflow` is organized around study design, not test names.
Instead of asking "which function runs a t-test?", testflow asks:
"what is the design of my clinical question?"

The core grammar is:

```r
testflow = test + interpret + plot
```

```r
library(testflow)

cardio <- make_cardio_data()

x <- cardio |>
  test_two_groups(sbp_3m, sex)

x
plot(x)
report(x)
as_tibble(x)
```

The public workflows are named for common study designs:

```r
test_one_sample(cardio, sbp_3m, mu = 140)
test_two_groups(cardio, sbp_3m, sex)
test_paired(cardio, sbp_baseline, sbp_3m)
test_groups(cardio, sbp_3m, treatment)
test_categorical(cardio, treatment, controlled_3m)
test_correlation(cardio, age, sbp_3m)
test_outliers(cardio, c(sbp_3m, ldl, crp))
```
