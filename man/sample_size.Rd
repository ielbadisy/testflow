\name{sample_size}
\alias{sample_size}
\alias{sample_size_continuous}
\alias{sample_size_binary}
\alias{sample_size_survival}
\alias{sample_size_ordinal}
\alias{sample_size_adjust_dropout}
\alias{plot.sample_size}
\alias{print.sample_size}
\alias{summary.sample_size}
\alias{as_tibble.sample_size}
\alias{report.sample_size}
\title{Sample size planning for common trial designs}
\usage{
sample_size(
  endpoint = c("continuous", "binary", "survival", "ordinal"),
  design = c("parallel", "paired", "repeated"),
  objective = c("superiority", "noninferiority", "equivalence"),
  ...
)

sample_size_continuous(
  design = c("parallel", "paired", "repeated"),
  objective = c("superiority", "noninferiority", "equivalence"),
  delta = NULL,
  sd = NULL,
  sd_diff = NULL,
  expected_difference = 0,
  margin = NULL,
  alpha = 0.05,
  power = 0.9,
  allocation = 1,
  dropout = 0,
  n_time = 2
)

sample_size_binary(
  design = c("parallel", "paired", "repeated"),
  objective = c("superiority", "noninferiority", "equivalence"),
  p1,
  p2,
  margin = NULL,
  discordant_or = NULL,
  discordance_rate = NULL,
  p10 = NULL,
  p01 = NULL,
  alpha = 0.05,
  power = 0.9,
  allocation = 1,
  dropout = 0,
  n_time = 2
)

sample_size_survival(
  design = c("parallel"),
  objective = c("superiority", "noninferiority", "equivalence"),
  hr,
  margin_hr = NULL,
  lower = NULL,
  upper = NULL,
  survival_a = NULL,
  survival_b = NULL,
  alpha = 0.05,
  power = 0.9,
  allocation = 1,
  dropout = 0,
  method = c("exponential", "ph_only")
)

sample_size_ordinal(
  design = c("parallel"),
  objective = c("superiority"),
  p_superiority,
  alpha = 0.05,
  power = 0.9,
  dropout = 0
)

sample_size_adjust_dropout(n, dropout = 0)
}
\arguments{
\item{endpoint}{Endpoint family: continuous, binary, survival, or ordinal.}
\item{design}{Study design: parallel, paired, or repeated. Repeated designs are routed to the paired formulas when \code{n_time = 2}.}
\item{objective}{Planning objective.}
\item{delta}{Effect size or margin on the continuous scale.}
\item{sd}{Between-subject standard deviation.}
\item{sd_diff}{Paired-difference standard deviation.}
\item{expected_difference}{Planned difference from the null boundary.}
\item{margin}{Non-inferiority or equivalence margin on the risk-difference or continuous scale.}
\item{p1, p2}{Anticipated binary response probabilities.}
\item{discordant_or}{Discordant odds ratio for paired binary planning.}
\item{discordance_rate}{Probability of a discordant pair.}
\item{p10, p01}{Alternative paired-binary inputs for the discordant cells.}
\item{hr}{Hazard ratio.}
\item{margin_hr}{Non-inferiority hazard-ratio margin.}
\item{lower, upper}{Equivalence bounds on the hazard-ratio scale.}
\item{survival_a, survival_b}{Survival probabilities at the planning horizon.}
\item{p_superiority}{Probability that a randomly selected subject from group A exceeds one from group B.}
\item{alpha}{Type I error rate.}
\item{power}{Target power.}
\item{allocation}{Allocation ratio \eqn{n_B / n_A}.}
\item{dropout}{Expected dropout proportion.}
\item{method}{Survival event-count method.}
\item{n}{Evaluable sample size before dropout adjustment.}
}
\value{
A \code{sample_size} object with the raw sample size, dropout-adjusted sample size, formula label, reference, report text, and a plot-ready summary.
}
\description{
Planning functions for common trial designs using formulas drawn from Julious \emph{Sample Sizes for Clinical Trials}.
}
\details{
The implementation currently covers:
\itemize{
\item Continuous endpoints: parallel and paired superiority, non-inferiority, and equivalence.
\item Binary endpoints: parallel risk-difference planning, paired discordant-pair planning, non-inferiority, and equivalence.
\item Survival endpoints: event-based superiority, non-inferiority, and equivalence.
\item Ordinal endpoints: Noether's superiority approximation.
}

The helper \code{sample_size_adjust_dropout()} applies \eqn{\lceil n / (1 - d) \rceil}.
}
\examples{
sample_size_continuous(
  design = "paired",
  objective = "superiority",
  delta = 5,
  sd_diff = 10,
  alpha = 0.05,
  power = 0.90
)

sample_size_binary(
  design = "paired",
  objective = "superiority",
  p10 = 0.2,
  p01 = 0.1,
  alpha = 0.05,
  power = 0.90
)
}
