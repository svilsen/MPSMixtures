% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.R
\name{potentialParentsMultiCore}
\alias{potentialParentsMultiCore}
\title{Multi-core potential parents}
\usage{
potentialParentsMultiCore(coverageTibble, stutterRatioModel,
  numberOfThreads = 4)
}
\arguments{
\item{coverageTibble}{A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.}

\item{stutterRatioModel}{A linear model fit; created by the \link{lm}-function.}

\item{numberOfThreads}{The maximum number of threads allowed.}
}
\value{
A list of the potential parents.
}
\description{
Simplified multi-core implementation of the potential parents function.
}
