# Per-workflow console-display metadata, keyed by x$workflow. Centralizing
# this here means print.testflow()/print.summary.testflow() get a
# workflow-specific title, field labels, and (where one exists) a per-term/
# per-metric results table without threading new parameters through every
# new_testflow() call site in R/test_*.R.
workflow_display <- list(
  two_groups = list(title = "Two Independent Groups"),
  paired = list(title = "Paired Measurements"),
  groups = list(
    title = "More Than Two Independent Groups",
    table = "posthoc", table_label = "Post-hoc comparisons"
  ),
  factorial = list(
    title = "Factorial ANOVA",
    table = "anova_table", table_label = "ANOVA table"
  ),
  repeated = list(
    title = "Repeated Measures",
    table = "posthoc", table_label = "Post-hoc comparisons"
  ),
  categorical = list(title = "Categorical Association"),
  paired_categorical = list(
    title = "Paired Categorical Measurements",
    table = "discordant_pairs", table_label = "Discordant pairs"
  ),
  repeated_categorical = list(
    title = "Repeated Categorical Measurements",
    table = "posthoc", table_label = "Pairwise McNemar"
  ),
  one_sample = list(title = "One-Sample Test"),
  proportion = list(title = "One-Sample Proportion"),
  multinomial = list(
    title = "Multinomial Goodness of Fit",
    table = "posthoc", table_label = "Pairwise category checks"
  ),
  correlation = list(
    title = "Correlation",
    table = "correlation_table", table_label = "Method comparison"
  ),
  correlation_matrix = list(
    title = "Correlation Matrix",
    table = "correlation_table", table_label = "Pairwise correlations"
  ),
  outliers = list(title = "Outlier Screening"),
  linear_regression = list(
    title = "Linear Regression",
    outcome_label = "Outcome", group_label = "Predictors",
    table = "coefficients", table_label = "Coefficients"
  ),
  logistic_regression = list(
    title = "Logistic Regression",
    outcome_label = "Outcome", group_label = "Predictors",
    table = "odds_ratios", table_label = "Odds ratios"
  ),
  survival = list(
    title = "Kaplan-Meier / Log-Rank",
    outcome_label = "Time", group_label = "Group",
    table = "cox", table_label = "Companion hazard ratio"
  ),
  cox = list(
    title = "Cox Proportional Hazards Regression",
    outcome_label = "Time", group_label = "Predictors",
    table = "hazard_ratios", table_label = "Hazard ratios"
  ),
  diagnostic = list(
    title = "Diagnostic Test Accuracy",
    outcome_label = "Test", group_label = "Reference",
    table = "diagnostic_table", table_label = "Diagnostic accuracy"
  ),
  roc = list(
    title = "ROC Curve",
    outcome_label = "Predictor", group_label = "Outcome",
    table = "optimal_threshold", table_label = "Optimal threshold (Youden's J)"
  ),
  agreement = list(
    title = "Inter-Rater Agreement",
    outcome_label = "Rater 1", group_label = "Rater 2",
    table = "agreement_table", table_label = "Agreement table"
  ),
  icc = list(
    title = "Intraclass Correlation",
    table = "icc_table", table_label = "ICC comparison"
  )
)

# Resolves display metadata for a testflow object: title, outcome/group
# field labels, and (when the registered field is already a data.frame/
# tibble, or a base table/matrix that can be tidied into one) a results
# table to print. Anything else the registered field could point at - most
# notably groups/repeated's post-hoc, which posthoc_groups() (R/posthoc.R)
# wraps as list(method=, result=, p.adjust.method=) around a raw TukeyHSD/
# pairwise.htest object rather than a tidy data frame - resolves table_data
# to NULL and is silently skipped by the caller, matching current behavior
# rather than crashing or dumping a raw object to the console.
workflow_meta <- function(x) {
  meta <- workflow_display[[x$workflow]]
  outcome_label <- meta$outcome_label %||% "Outcome"
  group_label <- meta$group_label %||% "Group"
  table_field <- meta$table
  table_data <- NULL
  if (!is.null(table_field)) {
    raw <- if (identical(table_field, "posthoc")) x$posthoc else x$alternative_tests[[table_field]]
    if (is.data.frame(raw)) {
      table_data <- raw
    } else if (is.table(raw) || is.matrix(raw)) {
      # as.data.frame(as.table(...)) always names its columns "Var1"/"Var2"/
      # "Freq", which is meaningless in printed output; relabel using the
      # same outcome/group labels already resolved for this workflow.
      table_data <- as.data.frame(as.table(raw))
      names(table_data) <- c(outcome_label, group_label, "n")
      table_data$n <- as.integer(table_data$n)
    }
  }
  list(
    title = meta$title %||% x$design %||% "Statistical test workflow",
    outcome_label = outcome_label,
    group_label = group_label,
    table_data = table_data,
    table_label = meta$table_label
  )
}
