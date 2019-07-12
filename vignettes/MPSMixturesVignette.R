## ----hidden, echo = FALSE------------------------------------------------
load_rdata <- function(fileName) {
    load(fileName)
    get(ls()[ls() != "fileName"])
}

## ----loadPackages, message = FALSE---------------------------------------
library("Biostrings")
library("tidyverse")
library("stringr")
library("Rsolnp")
library("MPSMixtures")
library("microbenchmark")

## ----populationParameters------------------------------------------------
theta = -1.0
data(examplePopulation)

## ----markerImbalances----------------------------------------------------
set.seed(123456)
markers = unique(examplePopulation$Marker)
numberOfMarkers = length(markers)

markerImbalances <- rgamma(numberOfMarkers, shape = 10, rate = 10)
markerImbalances <- markerImbalances / mean(markerImbalances)
markerImbalances

## ------------------------------------------------------------------------
stutterRatioTibble <- 
    tibble(BlockLengthMissingMotif = sample(1:40, 1000, replace = T),
           StutterRatio = 0.005 * BlockLengthMissingMotif + rnorm(1000, 0, 0.02)) %>% 
    mutate(StutterRatio = ifelse(StutterRatio < 0, 0, StutterRatio))

stutterRatioModel = lm(StutterRatio ~ 0 + I(BlockLengthMissingMotif - 1), 
                       data = stutterRatioTibble)

## ----genotypes-----------------------------------------------------------
numberOfContributors = 2
trueProfiles <- sampleGenotypesHWE(numberOfContributors, examplePopulation, c("V", "P"))

## ----alleleParameters----------------------------------------------------
nu = 1400
gamma = 8
phi = c(0.7, 0.3)

## ----noiseParameters-----------------------------------------------------
psi = 3
rho = 2

## ----sampledCoverage-----------------------------------------------------
sampleTibble = sampleCoverage(
    trueProfiles = trueProfiles, markerImbalances = markerImbalances,
    populationLadder = examplePopulation, stutterRatioModel = stutterRatioModel,
    alleleCoverageParameters = list(nu = nu, eta = gamma, phi = phi),
    noiseParameters = list(psi = psi, rho = rho, maxElements = 20), 
    p = NULL)

## ----potentialParentsList------------------------------------------------
potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel)

## ----bothKnown-----------------------------------------------------------
knownProfilesIndividual <- 
    estimateParametersOfKnownProfiles(sampleTibble, markerImbalances, trueProfiles, 
                                      potentialParentsList, 
                                      NULL, 1, 0.8, rep(1e-8, 4))

knownProfilesIndividual[c("Parameters", "LogLikelihoods", "Fitness")]

