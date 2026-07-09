# Per-workflow console-display metadata, keyed by x$workflow. Centralizing
# this here means print.testflow()/print.summary.testflow() get a
# workflow-specific title and field labels without threading new
# parameters through every new_testflow() call site in R/test_*.R.
workflow_display <- list(
  two_groups = list(title = "Two Independent Groups"),
  paired = list(title = "Paired Measurements"),
  groups = list(title = "More Than Two Independent Groups"),
  factorial = list(title = "Factorial ANOVA"),
  repeated = list(title = "Repeated Measures"),
  categorical = list(title = "Categorical Association"),
  paired_categorical = list(title = "Paired Categorical Measurements"),
  repeated_categorical = list(title = "Repeated Categorical Measurements"),
  one_sample = list(title = "One-Sample Test"),
  proportion = list(title = "One-Sample Proportion"),
  multinomial = list(title = "Multinomial Goodness of Fit"),
  correlation = list(title = "Correlation"),
  correlation_matrix = list(title = "Correlation Matrix"),
  outliers = list(title = "Outlier Screening"),
  linear_regression = list(
    title = "Linear Regression",
    outcome_label = "Outcome", group_label = "Predictors"
  ),
  logistic_regression = list(
    title = "Logistic Regression",
    outcome_label = "Outcome", group_label = "Predictors"
  ),
  survival = list(
    title = "Kaplan-Meier / Log-Rank",
    outcome_label = "Time", group_label = "Group"
  ),
  cox = list(
    title = "Cox Proportional Hazards Regression",
    outcome_label = "Time", group_label = "Predictors"
  ),
  diagnostic = list(
    title = "Diagnostic Test Accuracy",
    outcome_label = "Test", group_label = "Reference"
  ),
  roc = list(
    title = "ROC Curve",
    outcome_label = "Predictor", group_label = "Outcome"
  ),
  agreement = list(
    title = "Inter-Rater Agreement",
    outcome_label = "Rater 1", group_label = "Rater 2"
  ),
  icc = list(title = "Intraclass Correlation")
)

# Resolves display metadata (title, outcome/group field labels) for a
# testflow object, falling back to sensible generic defaults for any
# workflow not (yet) in the registry.
workflow_meta <- function(x) {
  meta <- workflow_display[[x$workflow]]
  list(
    title = meta$title %||% x$design %||% "Statistical test workflow",
    outcome_label = meta$outcome_label %||% "Outcome",
    group_label = meta$group_label %||% "Group"
  )
}
