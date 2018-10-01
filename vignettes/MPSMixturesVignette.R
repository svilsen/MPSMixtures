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

## ----sampledCoverage, cache = TRUE---------------------------------------
sampleTibble = sampleCoverage(
    trueProfiles = trueProfiles, markerImbalances = markerImbalances,
    populationLadder = examplePopulation, stutterRatioModel = stutterRatioModel,
    alleleCoverageParameters = list(nu = nu, eta = gamma, phi = phi),
    noiseParameters = list(psi = psi, rho = rho, maxElements = 20), 
    p = NULL)

## ----pp, cache = TRUE----------------------------------------------------
potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel)

## ----bothKnown-----------------------------------------------------------
knownProfilesIndividual <- 
    estimateParametersOfKnownProfiles(sampleTibble, markerImbalances, trueProfiles, 
                                      potentialParentsList, 
                                      NULL, 1, 0.8, rep(1e-8, 4))

knownProfilesIndividual[c("Parameters", "LogLikelihoods", "Fitness")]

## ----figQQ, echo = FALSE, fig.align = "center", fig.cap = "\\label{fig:qqtrue}Q-Q plot of the fitted model for both the allele and noise coverage models.", fig.width = 8, fig.height = 4, fig.pos="ht!"----
ggplotQQPlotProfiles(sampleTibble, knownProfilesIndividual)

## ----figEC, echo = FALSE, fig.align = "center", fig.cap = "\\label{fig:ectrue}Bar-plot of the coverage against the allele length for the each marker and profile, shown in coloumns and rows, respectively. Furthermore, the prediction interval for each observation is seen as the black bars, with the expected coverage at their center.", fig.width = 17, fig.height = 7, fig.pos="ht!"----
ggplotPredictionIntervals(sampleTibble, knownProfilesIndividual, contributorNames = c("V", "P"))

## ----perpProfile---------------------------------------------------------
knownPerpetrator <- setHypothesis(sampleTibble, numberOfContributors, trueProfiles["P"], theta)

## ----singleMajor, cache = TRUE-------------------------------------------
controlSinglePopulation <-  
    optimalUnknownProfileCombination.control(
        numberOfPopulations = 1, numberOfIterations = 350,
        populationSize = 2000, numberOfFittestIndividuals = 1,
        mutationDecayRate = 1, hillClimbingIterations = 0, 
        parentSelectionWindowSize = 15,
        allowParentSurvival = TRUE, trace = F, 
        levelsOfStutterRecursion = 1, 
        numberOfIterationsEqualMinMax = 10, 
        tolerance = rep(1e-8, 4)
    )

singlePopulationBenchmarkMajor <- microbenchmark(
    optimalSingleMajor <- 
        optimalUnknownProfileCombination(sampleTibble = sampleTibble, 
                                         markerImbalances = markerImbalances, 
                                         H = knownPerpetrator[[1]], 
                                         potentialParentsList = potentialParentsList, 
                                         control = controlSinglePopulation),
    times = 1)

## ----sampleParmsSingleMajor----------------------------------------------
optimalSingleMajor$U[[1]][["Parameters"]]

## ----multipleMajor, cache = TRUE-----------------------------------------
controlMultiplePopulation <- 
    optimalUnknownProfileCombination.control(
        numberOfPopulations = 32, numberOfIterations = 100,
        populationSize = 75, numberOfFittestIndividuals = 1,
        numberOfIterationsEqualMinMax = 10,
        mutationDecayRate = 1, hillClimbingIterations = 0, 
        parentSelectionWindowSize = 6,
        allowParentSurvival = TRUE, trace = F, 
        levelsOfStutterRecursion = 1, 
        tolerance = rep(1e-8, 4)
    )

multiplePopulationBenchmarkMajor <- microbenchmark(
    optimalMultipleMajor <- 
        optimalUnknownProfileCombination(sampleTibble = sampleTibble, 
                                         markerImbalances = markerImbalances, 
                                         H = knownPerpetrator[[1]], 
                                         potentialParentsList = potentialParentsList, 
                                         control = controlMultiplePopulation),
    times = 1)

## ----multipleMajorParameters---------------------------------------------
optimalMultipleMajor$U[[1]][c("LogLikelihoods", "Fitness", "Parameters")]

## ------------------------------------------------------------------------
knownVictim <- setHypothesis(sampleTibble, numberOfContributors, trueProfiles["V"], theta)

## ----multipleMinor, cache = TRUE-----------------------------------------
optimalMultipleMinor <- optimalUnknownProfileCombination(sampleTibble = sampleTibble,
                                                         markerImbalances = markerImbalances, 
                                                         H = knownVictim[[1]], 
                                                         potentialParentsList = potentialParentsList, 
                                                         control = controlMultiplePopulation)

## ----multipleMinorParameters---------------------------------------------
optimalMultipleMinor$U[[1]][c("Parameters", "LogLikelihoods", "Fitness")]

## ----plottingDifferences-------------------------------------------------
optimalMinor <- optimalMultipleMinor$U[[1]]
optimalVSTruePerpetratorProfile <- 
    ggplotComparingProfiles(sampleTibble, knownProfilesIndividual, 
                            optimalMinor, TRUE, c("True", "Optimal"))

## ----figDifference, echo = FALSE, fig.align = "center", fig.cap = "\\label{fig:minor}The coverage against allele, for both profiles, shown in coloumns, and each marker, shown in rows, where the two profiles differ. Profiles one and two correspond to the optimal profiles under the model and the true profiles, respectively. The red and blue colouring corresonds to the major and minor contributors, respectively. Lastely, the black bars are regions attributed to the noise distribution.", fig.width = 8, fig.height = 4, fig.pos="ht!"----
show(optimalVSTruePerpetratorProfile)

