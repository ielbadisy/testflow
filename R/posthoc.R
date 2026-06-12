posthoc_groups <- function(data, outcome, group, recommendation, alpha = 0.05) {
  formula <- stats::as.formula(paste(outcome, "~", group))
  if (recommendation == "One-way ANOVA") {
    fit <- stats::aov(formula, data = data)
    return(list(method = "Tukey HSD", result = stats::TukeyHSD(fit), p.adjust.method = "Tukey"))
  }

  if (recommendation == "Welch ANOVA") {
    pw <- stats::pairwise.t.test(
      data[[outcome]],
      data[[group]],
      pool.sd = FALSE,
      p.adjust.method = "BH"
    )
    return(list(method = "Games-Howell-style Welch pairwise t-tests", result = pw, p.adjust.method = "BH"))
  }

  pw <- stats::pairwise.wilcox.test(
    data[[outcome]],
    data[[group]],
    p.adjust.method = "BH",
    exact = FALSE
  )
  list(method = "Dunn-style pairwise Wilcoxon rank-sum tests", result = pw, p.adjust.method = "BH")
}

paired_posthoc_numeric <- function(wide, method = c("t", "wilcox"), p.adjust.method = "BH") {
  method <- match.arg(method)
  measures <- names(wide)
  pairs <- utils::combn(measures, 2, simplify = FALSE)
  out <- purrr::map_dfr(pairs, function(pair) {
    x <- wide[[pair[1]]]
    y <- wide[[pair[2]]]
    test <- if (method == "t") {
      stats::t.test(y, x, paired = TRUE)
    } else {
      stats::wilcox.test(y, x, paired = TRUE, exact = FALSE)
    }
    tibble::tibble(
      group1 = pair[1],
      group2 = pair[2],
      method = ifelse(method == "t", "paired t-test", "paired Wilcoxon signed-rank test"),
      statistic = unname(test$statistic[1]),
      p = test$p.value
    )
  })
  out$p.adj <- stats::p.adjust(out$p, method = p.adjust.method)
  out$p.adjust.method <- p.adjust.method
  out
}

pairwise_mcnemar <- function(mat, measure_names, p.adjust.method = "BH") {
  pairs <- utils::combn(seq_along(measure_names), 2, simplify = FALSE)
  out <- purrr::map_dfr(pairs, function(pair) {
    tab <- table(mat[, pair[1]], mat[, pair[2]])
    test <- stats::mcnemar.test(tab)
    tibble::tibble(
      group1 = measure_names[pair[1]],
      group2 = measure_names[pair[2]],
      method = "McNemar test",
      statistic = unname(test$statistic[1]),
      p = test$p.value
    )
  })
  out$p.adj <- stats::p.adjust(out$p, method = p.adjust.method)
  out$p.adjust.method <- p.adjust.method
  out
}
