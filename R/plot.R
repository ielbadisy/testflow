#' Plot a testflow object
#' @param x A testflow object.
#' @param title Optional plot title override. Defaults to the stored title.
#' @param subtitle Optional plot subtitle override. Defaults to the stored
#' subtitle.
#' @param caption Optional plot caption override. Defaults to the stored
#' caption.
#' @param ... Unused.
#' @return A `ggplot` object stored in the `testflow` object, with optional
#' title, subtitle, and caption overrides applied. If the workflow was created
#' with `plot = FALSE`, returns `NULL`.
#' @export
plot.testflow <- function(x, title = NULL, subtitle = NULL, caption = NULL, ...) {
  customize_plot_labels(x$plot, title = title, subtitle = subtitle, caption = caption)
}

customize_plot_labels <- function(plot, title = NULL, subtitle = NULL, caption = NULL) {
  if (is.null(plot)) return(NULL)
  labels <- list()
  if (!is.null(title)) labels$title <- title
  if (!is.null(subtitle)) labels$subtitle <- subtitle
  if (!is.null(caption)) labels$caption <- caption
  if (!length(labels)) return(plot)
  plot + do.call(ggplot2::labs, labels)
}

# Sober plot palette: greys/black plus two colour-blind-safe accents (Okabe-Ito blue and vermillion).
TFCOL <- c(ink = "grey15", pt = "grey45", line = "grey75", a1 = "#0072B2", a2 = "#D55E00")
tfpal <- function(k) c(TFCOL[c("ink", "a1", "a2")], paste0("grey", round(seq(40, 70, length.out = max(k - 3, 0)))))[
  seq_len(k)]
tftheme <- function() ggplot2::theme_minimal() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
  plot.caption = ggplot2::element_text(hjust = 0), legend.position = "right")
meanci <- function(y) {
  n <- length(y); m <- mean(y)
  h <- if (n > 1) stats::qt(0.975, n - 1) * stats::sd(y) / sqrt(n) else 0
  data.frame(y = m, ymin = m - h, ymax = m + h)
}
pstr <- function(p) {
  s <- format_p(p)
  if (startsWith(s, "<")) paste("p <", substring(s, 2)) else paste("p =", s)
}
effstr <- function(e) {
  if (is.null(e) || !nrow(e)) return(NULL)
  ci <- if (!is.null(e$conf.low) && !is.na(e$conf.low[1]))
    paste0(" (95% CI ", format_stat(e$conf.low[1]), " to ", format_stat(e$conf.high[1]), ")")
  paste0(e$name[1], " = ", format_stat(e$estimate[1]), ci)
}
# Wrap between terms only, never inside "term p = x" (hard spaces until wrapped).
termp <- function(lab, tab) {
  tt <- gsub(" ", "\u00a0", paste0(tab$term, " ", vapply(tab$p.value, pstr, "")))
  gsub("\u00a0", " ", paste(strwrap(paste0(lab, ": ", paste(tt, collapse = "; ")), 45), collapse = "\n"))
}
# Within levels in design order: keep factor levels, sort numbers, else order of first appearance.
lvls <- function(x) if (is.factor(x)) levels(x) else if (is.numeric(x)) sort(unique(x)) else unique(x)
MEANCAP <- "Blue: mean and 95% CI."

plot_subtitle <- function(recommended, primary) paste0(recommended, ", ", pstr(primary$p.value %||% NA_real_))

make_plot <- function(workflow, data, outcome = NULL, group = NULL, recommended = NULL, primary = NULL, effect = NULL,
                      extra = list()) {
  eff <- effstr(effect)
  subtitle <- if (!is.null(recommended) && !is.null(primary)) plot_subtitle(recommended, primary) else NULL
  if (workflow == "one_sample") {
    mu <- extra$mu %||% 0
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = .data[[outcome]])) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 20, fill = TFCOL[["line"]],
                                colour = "white") +
        ggplot2::geom_density(color = TFCOL[["a1"]], linewidth = 0.8) +
        ggplot2::geom_vline(xintercept = mu, linetype = "dashed") +
        ggplot2::labs(subtitle = subtitle, caption = eff, x = outcome, y = "Density") +
        tftheme()
    )
  }
  if (workflow == "two_groups" || workflow == "groups") {
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = .data[[group]], y = .data[[outcome]])) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.45, colour = TFCOL[["line"]], fill = NA) +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1),
                            colour = TFCOL[["pt"]], alpha = 0.6, size = 1.3) +
        ggplot2::stat_summary(fun.data = meanci, geom = "pointrange", colour = TFCOL[["a1"]],
                              position = ggplot2::position_nudge(x = 0.32)) +
        ggplot2::labs(subtitle = subtitle, caption = paste(c(eff, MEANCAP), collapse = "\n"), x = group,
                      y = outcome) +
        tftheme()
    )
  }
  if (workflow == "paired") {
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = .data$time, y = .data$value)) +
        ggplot2::geom_line(ggplot2::aes(group = .data$row_id), colour = TFCOL[["line"]]) +
        ggplot2::geom_point(colour = TFCOL[["pt"]], size = 1.3) +
        ggplot2::stat_summary(ggplot2::aes(group = 1), fun = mean, geom = "line", linewidth = 1,
                              colour = TFCOL[["a1"]]) +
        ggplot2::stat_summary(fun.data = meanci, geom = "pointrange", colour = TFCOL[["a1"]]) +
        ggplot2::labs(subtitle = subtitle, caption = paste(c(eff, MEANCAP), collapse = "\n"),
                      x = extra$x %||% "Measurement", y = extra$y %||% "Value") +
        tftheme()
    )
  }
  ggplot2::ggplot(data, ggplot2::aes(x = seq_len(nrow(data)))) +
    ggplot2::geom_point(ggplot2::aes(y = 1)) +
    ggplot2::labs(title = paste("testflow:", workflow), subtitle = subtitle, caption = eff, x = NULL, y = NULL) +
    ggplot2::theme_minimal()
}

# Means with 95% CI of y over x, one colour/shape per level of z; faint subject lines when id is given.
profile_plot <- function(df, y, x, z = NULL, id = NULL, subtitle = NULL, caption = NULL, xlab = x, ylab = y) {
  df <- as.data.frame(df)
  df[[x]] <- factor(df[[x]], lvls(df[[x]]))
  if (is.null(z)) { df$.z <- factor(""); z <- ".z" } else df[[z]] <- factor(df[[z]], lvls(df[[z]]))
  k <- nlevels(df[[z]]); one <- k == 1
  col <- if (one) TFCOL[["a1"]] else tfpal(k)
  dg <- ggplot2::position_dodge(width = if (one) 0 else 0.3)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data[[z]], shape = .data[[z]],
                                        group = .data[[z]]))
  grey <- if (one) list(colour = TFCOL[[if (is.null(id)) "pt" else "line"]])
  p <- p + if (is.null(id)) {
    do.call(ggplot2::geom_point, c(list(position = ggplot2::position_jitterdodge(jitter.width = 0.12,
      jitter.height = 0, dodge.width = if (one) 0.01 else 0.3, seed = 1), alpha = 0.3, size = 1.1), grey))
  } else {
    do.call(ggplot2::geom_line, c(list(ggplot2::aes(group = .data[[id]]), alpha = if (one) 0.6 else 0.15), grey))
  }
  p + ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 0.9, position = dg) +
    ggplot2::stat_summary(fun.data = meanci, geom = "pointrange", position = dg) +
    ggplot2::scale_colour_manual(values = unname(col)) +
    ggplot2::scale_shape_manual(values = rep(c(16, 17, 15, 18, 1, 2), length.out = k)) +
    ggplot2::labs(subtitle = subtitle, caption = caption, x = xlab, y = ylab,
                  colour = if (one) NULL else z, shape = if (one) NULL else z) +
    tftheme() + if (one) ggplot2::theme(legend.position = "none")
}
