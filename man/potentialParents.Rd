% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary.R
\name{potentialParents}
\alias{potentialParents}
\title{Potential parents list}
\usage{
potentialParents(coverageTibble, stutterRatioModel = NULL, trace = FALSE,
  simplifiedReturn = TRUE)
}
\arguments{
\item{coverageTibble}{A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.}

\item{stutterRatioModel}{A linear model fit; created by the \link{lm}-function.}

\item{trace}{If 'TRUE' adds a simple marker trace.}

\item{simplifiedReturn}{Should the returned list be simplified (TRUE/FALSE)?}
}
\value{
A list of the potential parents.
}
\description{
Creates a list of the potential parents for every region in the provided sample.
}
