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

title_case_method <- function(x) {
  tools::toTitleCase(as.character(x))
}

# Specific, data-driven reason for an automatically chosen test: lists the
# Shapiro-Wilk (and, when used, Levene) p-values behind the decision.
recommend_rationale <- function(test, normality, variance = NULL, alpha = 0.05, subject = NULL, extra = NULL) {
  np <- if ("p_value" %in% names(normality)) normality$p_value else normality$p
  lab <- if ("group" %in% names(normality)) normality$group else sub("^Normality: ", "", normality$name)
  p_eq <- function(p) ifelse(is.na(p), "p not available", ifelse(grepl("^<", format_p(p)), paste0("p ", format_p(p)), paste0("p = ", format_p(p))))
  norm_txt <- paste0(
    "Shapiro-Wilk normality", if (!is.null(subject)) paste0(" of ", subject) else "", " ",
    paste0(p_eq(np), if (length(lab) > 1) paste0(" (", lab, ")") else "", collapse = ", ")
  )
  normal_ok <- all(normality$status == "acceptable")
  var_txt <- NULL
  variance_ok <- TRUE
  if (!is.null(variance)) {
    vp <- if ("p_value" %in% names(variance)) variance$p_value[1] else variance$p[1]
    var_txt <- paste0("Levene ", p_eq(vp))
    variance_ok <- variance$status[1] == "acceptable"
  }
  reason <- if (!normal_ok) {
    "at least one normality check failed or could not be run, so a rank-based test was used"
  } else if (!variance_ok) {
    "normality looks acceptable but variances look unequal, so a test that does not assume equal variances was used"
  } else if (!is.null(variance)) {
    "normality and equal variances look acceptable"
  } else {
    "normality looks acceptable"
  }
  paste0(
    test, " chosen: ", paste(c(norm_txt, var_txt), collapse = "; "),
    " (alpha = ", format_stat(alpha), "); ", reason,
    if (!is.null(extra)) paste0("; ", extra) else "", "."
  )
}
