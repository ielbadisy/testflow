lt <- function(m) m[lower.tri(m, TRUE)]
ph_dat <- function() {
  set.seed(3)
  list(tooth = transform(ToothGrowth, grp = interaction(supp, dose)),
       unbal = data.frame(len = c(round(rnorm(7), 1), round(rnorm(15, .8), 1), rnorm(4, 1.5), rexp(9), NA),
                          grp = factor(c(rep(letters[1:4], c(7, 15, 4, 9)), "a"), levels = c(letters[1:4], "z"))))
}

test_that("Conover-Iman and DSCF match PMCMRplus with ties, unequal sizes and NA", {
  skip_if_not_installed("PMCMRplus")
  for (d in ph_dat()) {
    cc <- d[complete.cases(d), ]; y <- cc$len; g <- droplevels(cc$grp)
    for (a in c("none", "holm", "BH", "hommel")) {
      ct <- testgroups(len ~ grp, data = d, posthoc = "conover", p.adjust.method = a, plot = FALSE)$posthoc
      pc <- suppressWarnings(PMCMRplus::kwAllPairsConoverTest(y, g, p.adjust.method = a))
      expect_equal(ct$table$p.adj, lt(pc$p.value), tolerance = 1e-8)
      expect_equal(ct$result$p.value, pc$p.value, tolerance = 1e-8)
    }
    expect_equal(ct$table$statistic, lt(pc$statistic), tolerance = 1e-8)
    ds <- testgroups(len ~ grp, data = d, posthoc = "dscf", plot = FALSE)$posthoc
    pd <- PMCMRplus::dscfAllPairsTest(y, g)
    expect_equal(ds$table$p.adj, lt(pd$p.value), tolerance = 1e-8)
    expect_equal(abs(ds$table$statistic), abs(lt(pd$statistic)), tolerance = 1e-8)
    expect_equal(ds$table$estimate[1], median(outer(y[g == levels(g)[2]], y[g == levels(g)[1]], "-")))
    dn <- testgroups(len ~ grp, data = d, posthoc = "dunn", p.adjust.method = "hochberg", plot = FALSE)$posthoc
    pn <- suppressWarnings(PMCMRplus::kwAllPairsDunnTest(y, g, p.adjust.method = "hochberg"))
    expect_equal(dn$table$p.adj, lt(pn$p.value), tolerance = 1e-8)
    gh <- testgroups(len ~ grp, data = d, posthoc = "gameshowell", plot = FALSE)$posthoc
    pg <- PMCMRplus::gamesHowellTest(y, g)
    expect_equal(gh$table$p.adj, lt(pg$p.value), tolerance = 1e-6)
  }
})

test_that("rank post hoc tests agree with base R where base R has an equivalent", {
  d <- ph_dat()$unbal; d <- d[complete.cases(d), ]; d$grp <- droplevels(d$grp)
  ct <- testgroups(len ~ grp, data = d, posthoc = "conover", p.adjust.method = "none", plot = FALSE)$posthoc$table
  # The Conover-Iman t is the pooled-variance t test on the ranks scaled by the residual of the rank ANOVA.
  h <- unname(kruskal.test(len ~ grp, d)$statistic); r <- rank(d$len); N <- nrow(d); k <- nlevels(d$grp)
  rb <- tapply(r, d$grp, mean); n <- tapply(r, d$grp, length)
  s2 <- var(r) * (N - 1 - h) / (N - k)
  expect_equal(ct$statistic[1], unname((rb[2] - rb[1]) / sqrt(s2 * (1 / n[1] + 1 / n[2]))), tolerance = 1e-10)
  expect_equal(ct$estimate, unname(c(outer(rb, rb, "-")[lower.tri(diag(k))])), tolerance = 1e-10)
  tk <- testgroups(len ~ grp, data = d, posthoc = "tukey", plot = FALSE)$posthoc$table
  ref <- TukeyHSD(aov(len ~ grp, d))$grp
  expect_equal(tk$conf.low, unname(ref[, "lwr"]), tolerance = 1e-10)
})

test_that("small and singleton groups are handled", {
  d <- data.frame(y = c(1, 2, 3, 4, 5, 6, 7, 8, 9), g = rep(c("a", "b", "c"), c(4, 4, 1)))
  x <- testgroups(y ~ g, data = d, posthoc = "dunn", plot = FALSE)
  expect_equal(x$assumptions$status[x$assumptions$name == "Bartlett test"], "not checked")
  expect_true(is.na(x$alternative_tests$welch$p.value))
  expect_equal(nrow(x$posthoc$table), 3)
  gh <- testgroups(y ~ g, data = d, posthoc = "gameshowell", plot = FALSE)$posthoc$table
  expect_true(all(is.na(gh$p[gh$group2 == "c"])) && !is.na(gh$p[1]))
  expect_error(testgroups(y ~ g, data = d, posthoc = "welch", plot = FALSE), "at least 2 observations")
  expect_equal(nrow(testgroups(y ~ g, data = d, posthoc = "dscf", plot = FALSE)$posthoc$table), 3)
})

test_that("Friedman Nemenyi and Conover match PMCMRplus, with ties and incomplete subjects", {
  skip_if_not_installed("PMCMRplus")
  set.seed(4)
  w <- data.frame(id = 1:40, t1 = round(rnorm(40)), t2 = round(rnorm(40, .4)), t3 = round(rnorm(40, .8)),
                  t4 = round(rnorm(40, 1)))
  w$t2[3] <- NA
  m <- as.matrix(w[complete.cases(w), -1]); rownames(m) <- seq_len(nrow(m))
  # The incomplete subject makes aov() warn "Error() model is singular" for the ANOVA; post hoc tests use
  # the complete subjects.
  rp <- function(...) suppressWarnings(testrepeated(w, c(t1, t2, t3, t4), id = id, plot = FALSE, ...))$posthoc
  nm <- rp(posthoc = "nemenyi")
  pn <- PMCMRplus::frdAllPairsNemenyiTest(m)
  expect_equal(nm$p.adj, lt(pn$p.value), tolerance = 1e-8)
  expect_equal(abs(nm$statistic), lt(pn$statistic), tolerance = 1e-8)
  for (a in c("none", "holm", "BY")) {
    cv <- rp(posthoc = "conover", p.adjust.method = a)
    pv <- PMCMRplus::frdAllPairsConoverTest(m, p.adjust.method = a)
    expect_equal(cv$p.adj, lt(pv$p.value), tolerance = 1e-8)
  }
  expect_equal(abs(cv$statistic), abs(lt(pv$statistic)), tolerance = 1e-8)
  r <- t(apply(m, 1, rank)); expect_equal(cv$estimate[1], mean(r[, 2] - r[, 1]))
})

test_that("repeated posthoc options keep the default and accept p.adjust.method", {
  set.seed(5)
  w <- data.frame(id = 1:25, t1 = rexp(25), t2 = rexp(25, .7), t3 = rexp(25, .5))
  a <- testrepeated(w, c(t1, t2, t3), id = id, plot = FALSE)$posthoc
  expect_equal(unique(a$method), "paired Wilcoxon signed-rank test")
  t <- testrepeated(w, c(t1, t2, t3), id = id, posthoc = "t", p.adjust.method = "bonferroni", plot = FALSE)$posthoc
  expect_equal(t$p[1], t.test(w$t2, w$t1, paired = TRUE)$p.value)
  expect_equal(t$p.adj, p.adjust(t$p, "bonferroni"))
  expect_null(testrepeated(w, c(t1, t2, t3), id = id, posthoc = "none", plot = FALSE)$posthoc)
  l <- data.frame(id = rep(w$id, 3), time = rep(c("t1", "t2", "t3"), each = 25), v = c(w$t1, w$t2, w$t3))
  expect_equal(testrepeatedlong(l, v, time, id, posthoc = "nemenyi", plot = FALSE)$posthoc$p,
               testrepeated(w, c(t1, t2, t3), id = id, posthoc = "nemenyi", plot = FALSE)$posthoc$p)
})

test_that("post hoc plots draw for every method", {
  d <- ph_dat()$tooth
  for (m in c("tukey", "gameshowell", "welch", "wilcoxon", "dunn", "conover", "dscf")) {
    x <- testgroups(len ~ grp, data = d, posthoc = m)
    p <- plot(x, type = "posthoc")
    expect_s3_class(p, "ggplot"); expect_silent(ggplot2::ggplot_build(p))
  }
  expect_equal(plot(x, title = "T")$labels$title, "T")
  expect_error(plot(testgroups(len ~ grp, data = d, posthoc = "none"), type = "posthoc"), "No post hoc")
  w <- data.frame(id = 1:20, t1 = rexp(20), t2 = rexp(20), t3 = rexp(20))
  for (m in c("t", "nemenyi", "conover")) {
    p <- plot(testrepeated(w, c(t1, t2, t3), id = id, posthoc = m), type = "posthoc")
    expect_s3_class(p, "ggplot"); expect_silent(ggplot2::ggplot_build(p))
  }
})
