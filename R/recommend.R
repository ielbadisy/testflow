recommend_two_groups <- function(normality, variance) {
  normal_ok <- all(normality$status == "acceptable")
  variance_ok <- variance$status[1] == "acceptable"
  if (normal_ok && variance_ok) "Student independent t-test" else if (normal_ok) "Welch t-test" else "Wilcoxon rank-sum test"
}

recommend_one_sample <- function(normality) {
  if (normality$status[1] == "acceptable") "One-sample t-test" else "Wilcoxon signed-rank test"
}

recommend_paired <- function(normality) {
  if (normality$status[1] == "acceptable") "Paired t-test" else "Wilcoxon signed-rank test"
}

recommend_groups <- function(normality, variance) {
  normal_ok <- all(normality$status == "acceptable")
  variance_ok <- variance$status[1] == "acceptable"
  if (normal_ok && variance_ok) "One-way ANOVA" else if (normal_ok) "Welch ANOVA" else "Kruskal-Wallis test"
}
