an <- function() {
  skip_if_not_installed("MASS")
  transform(MASS::anorexia, change = Postwt - Prewt)
}
rmdat <- function() {
  set.seed(7); n <- 24
  s <- rnorm(n, 0, 3); g <- rep(c("a", "b"), each = n / 2)
  M <- sapply(1:4, function(t) 10 + t * (g == "b") + s * t / 2 + rnorm(n))
  data.frame(id = seq_len(n), g = g, setNames(as.data.frame(M), paste0("m", 1:4)))
}

test_that("testrepeated reports Greenhouse-Geisser and Huynh-Feldt corrections as anova.mlm", {
  d <- rmdat(); M <- as.matrix(d[, paste0("m", 1:4)])
  x <- testrepeated(d, measures = m1:m4, id = id, plot = FALSE)
  a <- x$alternative_tests$anova
  ref <- stats::anova(lm(M ~ 1), X = ~1, test = "Spherical")
  expect_equal(a$statistic, ref$F[1], tolerance = 1e-8)
  expect_equal(a$p.gg, ref$`G-G Pr`[1], tolerance = 1e-8)
  expect_equal(a$p.hf, ref$`H-F Pr`[1], tolerance = 1e-8)
  adj <- summary(car::Anova(lm(M ~ 1), idata = data.frame(time = factor(1:4)), idesign = ~time),
                 multivariate = FALSE)$pval.adjustments
  expect_equal(c(a$gg.eps, a$hf.eps), unname(adj[1, c("GG eps", "HF eps")]), tolerance = 1e-8)
  expect_match(x$interpretation, "Greenhouse-Geisser epsilon")
})

test_that("mixed ANOVA within terms carry sphericity corrections (anova.mlm, afex)", {
  d <- rmdat(); M <- as.matrix(d[, paste0("m", 1:4)])
  m <- testrepeated(d, measures = m1:m4, id = id, between = g, plot = FALSE)$alternative_tests$mixed_anova
  ref <- stats::anova(lm(M ~ g, d), X = ~1, test = "Spherical")
  w <- m$term != "g"
  expect_true(all(is.na(m$gg.eps[!w])))
  expect_equal(m$p.gg[w], ref$`G-G Pr`[1:2], tolerance = 1e-8)
  expect_equal(m$p.hf[w], ref$`H-F Pr`[1:2], tolerance = 1e-8)
  skip_if_not_installed("afex")
  L <- reshape(d, direction = "long", varying = paste0("m", 1:4), v.names = "y", timevar = "time", idvar = "id")
  L$time <- factor(L$time)
  af <- suppressMessages(afex::aov_ez("id", "y", L, within = "time", between = "g"))
  adj <- summary(af)$pval.adjustments
  expect_equal(unname(m$gg.eps[w]), unname(adj[, "GG eps"]), tolerance = 1e-8)
  expect_equal(unname(m$hf.eps[w]), unname(adj[, "HF eps"]), tolerance = 1e-8)
})

test_that("testtwogroups gives Hedges g with CI as effectsize", {
  skip_if_not_installed("effectsize")
  d <- droplevels(subset(an(), Treat %in% c("FT", "Cont")))
  e <- testtwogroups(change ~ Treat, data = d, method = "welch", plot = FALSE)$effect_size
  expect_equal(e$name, c("Cohen's d", "Hedges' g"))
  g <- effectsize::hedges_g(change ~ Treat, data = d); cd <- effectsize::cohens_d(change ~ Treat, data = d)
  expect_equal(e$estimate[2], g$Hedges_g, tolerance = 1e-6)
  expect_equal(c(e$conf.low[2], e$conf.high[2]), c(g$CI_low, g$CI_high), tolerance = 1e-5)
  expect_equal(c(e$conf.low[1], e$conf.high[1]), c(cd$CI_low, cd$CI_high), tolerance = 1e-5)
  expect_equal(round(e$estimate[2], 3), -0.987)
  expect_equal(e$comparison[1], "Cont vs FT")
})

test_that("rank-biserial r has a DeLong CI and states which group U refers to", {
  skip_if_not_installed("pROC")
  d <- droplevels(subset(an(), Treat %in% c("FT", "Cont")))
  x <- testtwogroups(change ~ Treat, data = d, method = "wilcoxon", plot = FALSE)
  e <- x$effect_size
  a <- d$change[d$Treat == "Cont"]; b <- d$change[d$Treat == "FT"]
  r <- pROC::roc(controls = b, cases = a, direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  expect_equal(e$estimate[2], ci[2], tolerance = 1e-10)
  expect_equal(c(e$conf.low[2], e$conf.high[2]), ci[c(1, 3)], tolerance = 1e-10)
  expect_equal(c(e$conf.low[1], e$conf.high[1]), 2 * ci[c(1, 3)] - 1, tolerance = 1e-10)
  expect_equal(e$estimate[1], 2 * x$primary_test$statistic / (length(a) * length(b)) - 1)
  expect_match(x$interpretation, "W is the Mann-Whitney U of Cont")
})

test_that("testgroups lets the user pick the post hoc method (Tukey, Games-Howell, Dunn)", {
  d <- an()
  tk <- testgroups(change ~ Treat, data = d, posthoc = "tukey", plot = FALSE)$posthoc
  ref <- TukeyHSD(aov(change ~ Treat, d))$Treat
  expect_equal(tk$method, "Tukey HSD")
  expect_equal(tk$table$p.adj, unname(ref[, "p adj"]), tolerance = 1e-10)
  expect_equal(tk$table$estimate, unname(ref[, "diff"]), tolerance = 1e-10)
  expect_equal(paste(tk$table$group2, tk$table$group1, sep = "-"), rownames(ref))
  expect_null(testgroups(change ~ Treat, data = d, posthoc = FALSE, plot = FALSE)$posthoc)
  expect_equal(testgroups(change ~ Treat, data = d, posthoc = TRUE, plot = FALSE)$posthoc$method,
               "Pairwise Wilcoxon rank-sum tests (BH-adjusted)")
  skip_if_not_installed("PMCMRplus")
  gh <- testgroups(change ~ Treat, data = d, posthoc = "gameshowell", plot = FALSE)$posthoc
  pg <- PMCMRplus::gamesHowellTest(change ~ Treat, data = d)
  expect_equal(gh$table$p.adj, pg$p.value[lower.tri(pg$p.value, TRUE)], tolerance = 1e-6)
  expect_equal(gh$result$p.value, pg$p.value, tolerance = 1e-6)
  dn <- testgroups(change ~ Treat, data = d, posthoc = "dunn", p.adjust.method = "holm", plot = FALSE)$posthoc
  pd <- suppressWarnings(PMCMRplus::kwAllPairsDunnTest(change ~ Treat, data = d, p.adjust.method = "holm"))
  expect_equal(dn$result$p.value, pd$p.value, tolerance = 1e-10)
  expect_equal(abs(dn$table$statistic), pd$statistic[lower.tri(pd$statistic, TRUE)], tolerance = 1e-10)
})

test_that("testdunn matches kwAllPairsDunnTest with ties", {
  skip_if_not_installed("PMCMRplus")
  d <- transform(ToothGrowth, grp = interaction(supp, dose))
  x <- testdunn(len ~ grp, data = d, p.adjust.method = "none")
  ref <- suppressWarnings(PMCMRplus::kwAllPairsDunnTest(len ~ grp, data = d, p.adjust.method = "none"))
  expect_equal(x$p, ref$p.value[lower.tri(ref$p.value, TRUE)], tolerance = 1e-10)
  expect_equal(nrow(x), choose(6, 2))
  expect_equal(testdunn(d, len, grp)$p.adj, stats::p.adjust(x$p, "BH"))
})

test_that("Kruskal-Wallis effect sizes are labelled eta^2_H and epsilon^2_H", {
  d <- an()
  e <- testgroups(change ~ Treat, data = d, plot = FALSE)$effect_size
  h <- unname(kruskal.test(change ~ Treat, d)$statistic)
  expect_equal(e$name, c("Kruskal eta squared (H)", "Kruskal epsilon squared (H)"))
  expect_equal(e$estimate, c((h - 2) / 69, h / 71))
  skip_if_not_installed("effectsize")
  expect_equal(e$estimate[2], effectsize::rank_epsilon_squared(change ~ Treat, data = d)$rank_epsilon_squared,
               tolerance = 1e-8)
})

test_that("testpaired adds the matched-pairs rank-biserial r_c for the Wilcoxon test", {
  skip_if_not_installed("effectsize")
  d <- subset(MASS::anorexia, Treat == "CBT")
  e <- testpaired(Postwt ~ Prewt, data = d, plot = FALSE)$effect_size
  expect_equal(e$name, c("Cohen's dz", "Matched-pairs rank-biserial r"))
  expect_equal(e$estimate[2], effectsize::rank_biserial(d$Postwt, d$Prewt, paired = TRUE)$r_rank_biserial,
               tolerance = 1e-8)
  expect_equal(round(e$estimate[2], 3), 0.395)
  dz <- testpaired(after ~ before, data = data.frame(before = 1:30, after = 1:30 + rnorm(30, 1)), plot = FALSE)
  expect_equal(nrow(dz$effect_size), if (dz$recommended$test == "Paired t-test") 1 else 2)
})
