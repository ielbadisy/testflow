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
  n_time = 2,
  correlation = 0.5
)

sample_size_binary(
  design = c("parallel", "paired", "repeated"),
  objective = c("superiority", "noninferiority", "equivalence"),
  p1,
  p2,
  margin = NULL,
  method = c("pooled", "anticipated"),
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
  method = c("exponential", "ph_only"),
  accrual_duration = NULL,
  follow_up = NULL
)

sample_size_ordinal(
  design = c("parallel"),
  objective = c("superiority"),
  p_superiority = NULL,
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
\item{method}{For \code{sample_size_binary()}, the parallel superiority variance estimator: \code{pooled} (null variance) or \code{anticipated} (alternative-hypothesis variance). For \code{sample_size_survival()}, the event-count method: \code{exponential} or \code{ph_only}.}
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
\item{n_time}{Number of repeated measures.}
\item{correlation}{Within-subject correlation for repeated measures when \code{design = "repeated"} and \code{n_time > 2}.}
\item{n}{Evaluable sample size before dropout adjustment.}
\item{accrual_duration}{Uniform accrual (enrollment) duration for \code{sample_size_survival()}. When supplied with \code{follow_up}, event probabilities account for staggered enrollment. Requires \code{survival_a}/\code{survival_b}.}
\item{follow_up}{Additional follow-up duration after accrual ends, for \code{sample_size_survival()}. Total study duration is \code{accrual_duration + follow_up}.}
\item{...}{Endpoint-specific arguments passed to the selected helper.}
}
\value{
A \code{sample_size} object with the raw sample size, dropout-adjusted sample size, method label, assumptions, report text, and a plot-ready summary. Formulas and literature references are documented in the individual function help pages, not returned in the object.
}
\description{
Planning functions for common trial designs across continuous, binary, survival, and ordinal endpoints.
}
\details{
The implementation currently covers:
\itemize{
\item Continuous endpoints: parallel and paired superiority, non-inferiority, and equivalence, plus repeated-measures approximation for 3+ time points.
\item Binary endpoints: parallel risk-difference planning, paired discordant-pair planning, non-inferiority, and equivalence.
\item Survival endpoints: event-based superiority, non-inferiority, and equivalence.
\item Ordinal endpoints: Noether's superiority approximation.
}

The helper \code{sample_size_adjust_dropout()} applies \eqn{\lceil n / (1 - d) \rceil}.

See the individual function help pages (e.g. \code{?sample_size_continuous}) for the underlying formulas and their literature references.
}
\references{
Julious SA. \emph{Sample Sizes for Clinical Trials}. Chapman & Hall/CRC; 2010.
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

sample_size_continuous(
  design = "repeated",
  n_time = 4,
  correlation = 0.5,
  objective = "superiority",
  delta = 5,
  sd_diff = 10,
  alpha = 0.05,
  power = 0.90
)
}
