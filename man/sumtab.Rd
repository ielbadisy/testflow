\name{sumtab}
\alias{sumtab}
\title{Build a compact descriptive summary table}
\usage{
sumtab(
  formula,
  data,
  p_value = FALSE,
  overall = TRUE,
  digits = 1,
  p_digits = 3,
  alpha = 0.05,
  fisher_threshold = 5,
  na.rm = TRUE
)
}
\arguments{
\item{formula}{A one-sided formula such as \code{~ age + sex} or
\code{~ age + sex | treatment}.}

\item{data}{A data frame.}

\item{p_value}{Logical; add a p-value column when a grouping variable is
supplied.}

\item{overall}{Logical; include an overall summary column.}

\item{digits}{Number of digits for summary statistics.}

\item{p_digits}{Number of digits for formatted p-values.}

\item{alpha}{Significance level used by automatic test selection.}

\item{fisher_threshold}{Expected-count threshold for Fisher's exact test.}

\item{na.rm}{Logical; remove missing values before summaries and tests.}
}
\value{
A tibble with one row per numeric variable and one row per categorical level.
}
\description{
Builds a compact table of descriptive summaries using a formula interface.
Numeric variables are summarized as mean (SD); median [Q1, Q3]; n.
Categorical variables are summarized as n (percent).
}
\details{
When \code{p_value = TRUE} and a grouping variable is supplied,
\code{sumtab()} chooses the p-value test automatically. Numeric variables use
Student t-test, Welch t-test, or Wilcoxon rank-sum test for two groups, and
one-way ANOVA, Welch ANOVA, or Kruskal-Wallis test for more than two groups.
Categorical variables use a chi-square test unless expected counts fall below
\code{fisher_threshold}, in which case Fisher's exact test is used.
}
\examples{
dat <- make_cardio_data(80, seed = 1)
sumtab(~ age + sex | treatment, dat, p_value = TRUE)
}
