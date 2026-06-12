#' @export
plot.testflow <- function(x, ...) {
  x$plot
}

plot_subtitle <- function(recommended, primary) {
  paste0(recommended, ", p = ", format_p(primary$p.value %||% NA_real_))
}

make_plot <- function(workflow, data, outcome = NULL, group = NULL, recommended = NULL, primary = NULL, effect = NULL, extra = list()) {
  caption <- if (!is.null(effect) && nrow(effect) > 0) paste0(effect$name[1], " = ", format_stat(effect$estimate[1]), ", ", effect$magnitude[1]) else NULL
  subtitle <- if (!is.null(recommended) && !is.null(primary)) plot_subtitle(recommended, primary) else NULL

  if (workflow == "one_sample") {
    mu <- extra$mu %||% 0
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = .data[[outcome]])) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 20, fill = "#4C78A8", color = "white") +
        ggplot2::geom_density(color = "#F58518", linewidth = 0.8) +
        ggplot2::geom_vline(xintercept = mu, linetype = "dashed") +
        ggplot2::labs(title = paste("One-sample workflow:", outcome), subtitle = subtitle, caption = caption, x = outcome, y = "Density") +
        ggplot2::theme_minimal()
    )
  }
  if (workflow == "two_groups" || workflow == "groups") {
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = .data[[group]], y = .data[[outcome]], color = .data[[group]])) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5) +
        ggplot2::geom_jitter(width = 0.12, alpha = 0.65) +
        ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +
        ggplot2::labs(title = paste("Group comparison:", outcome), subtitle = subtitle, caption = caption, x = group, y = outcome) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
    )
  }
  if (workflow == "paired") {
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = .data$time, y = .data$value, group = .data$row_id)) +
        ggplot2::geom_line(alpha = 0.35) +
        ggplot2::geom_point(size = 2, color = "#4C78A8") +
        ggplot2::stat_summary(ggplot2::aes(group = 1), fun = mean, geom = "line", linewidth = 1.1, color = "#F58518") +
        ggplot2::labs(title = "Paired measurement workflow", subtitle = subtitle, caption = caption, x = NULL, y = outcome) +
        ggplot2::theme_minimal()
    )
  }
  ggplot2::ggplot(data, ggplot2::aes(x = seq_len(nrow(data)))) +
    ggplot2::geom_point(ggplot2::aes(y = 1)) +
    ggplot2::labs(title = paste("testflow:", workflow), subtitle = subtitle, caption = caption, x = NULL, y = NULL) +
    ggplot2::theme_minimal()
}
