% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LR.R
\name{setHypothesis}
\alias{setHypothesis}
\title{Set a hypothesis.}
\usage{
setHypothesis(sampleTibble, numberOfContributors, knownProfilesList, theta)
}
\arguments{
\item{sampleTibble}{A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.}

\item{numberOfContributors}{The total number of contributors to the mixture. Note: can be a vector of possible hypotheses, but elements should always be larger than or equal to the number of known profiles.}

\item{knownProfilesList}{A list of tibbles containing the alleles of the known contributors.}

\item{theta}{The inbreeding coefficient (Fst).}
}
\value{
A list containing the information relavent to the hypothesis.
}
\description{
Set a hypothesis given the sample, the number of total contributors, a list of the known contributors, and population parameters and data.
}
