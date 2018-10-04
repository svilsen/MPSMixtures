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

## ----figQQ, echo = FALSE, fig.align = "center", fig.cap = "\\label{fig:qqtrue}Q-Q plot of the fitted model for both the allele and noise coverage models.", fig.width = 8, fig.height = 4, fig.pos="ht!"----
ggplotQQPlotProfiles(sampleTibble, knownProfilesIndividual)

## ----figEC, echo = FALSE, fig.align = "center", fig.cap = "\\label{fig:ectrue}Bar-plot of the coverage against the allele length for the each marker and profile, shown in coloumns and rows, respectively. Furthermore, the prediction interval for each observation is seen as the black bars, with the expected coverage at their center.", fig.width = 17, fig.height = 7, fig.pos="ht!"----
ggplotPredictionIntervals(sampleTibble, knownProfilesIndividual, contributorNames = c("V", "P"))

## ----perpProfile---------------------------------------------------------
knownPerpetrator <- setHypothesis(sampleTibble, numberOfContributors, trueProfiles["P"], theta)

## ----singleMajor, eval = FALSE-------------------------------------------
#  controlSinglePopulation <-
#      optimalUnknownProfileCombination.control(
#          numberOfPopulations = 1, numberOfIterations = 350,
#          populationSize = 2000, numberOfFittestIndividuals = 1,
#          mutationDecayRate = 1, hillClimbingIterations = 0,
#          parentSelectionWindowSize = 15,
#          allowParentSurvival = TRUE, trace = F,
#          levelsOfStutterRecursion = 1,
#          numberOfIterationsEqualMinMax = 10,
#          tolerance = rep(1e-8, 4)
#      )
#  
#  singlePopulationBenchmarkMajor <- microbenchmark(
#      optimalSingleMajor <-
#          optimalUnknownProfileCombination(sampleTibble = sampleTibble,
#                                           markerImbalances = markerImbalances,
#                                           H = knownPerpetrator[[1]],
#                                           potentialParentsList = potentialParentsList,
#                                           control = controlSinglePopulation),
#      times = 1)

## ----singleMajor2, echo = FALSE------------------------------------------
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

optimalSingleMajor <- load_rdata(system.file('extdata', "singleMajor.RData", package = 'MPSMixtures'))
singlePopulationBenchmarkMajor <- load_rdata(system.file('extdata', "singleMajorBench.RData", package = 'MPSMixtures'))

## ----sampleParmsSingleMajor----------------------------------------------
optimalSingleMajor$U[[1]][["Parameters"]]

## ----multipleMajor, eval = FALSE-----------------------------------------
#  controlMultiplePopulation <-
#      optimalUnknownProfileCombination.control(
#          numberOfPopulations = 32, numberOfIterations = 100,
#          populationSize = 75, numberOfFittestIndividuals = 1,
#          numberOfIterationsEqualMinMax = 10,
#          mutationDecayRate = 1, hillClimbingIterations = 0,
#          parentSelectionWindowSize = 6,
#          allowParentSurvival = TRUE, trace = F,
#          levelsOfStutterRecursion = 1,
#          tolerance = rep(1e-8, 4)
#      )
#  
#  multiplePopulationBenchmarkMajor <- microbenchmark(
#      optimalMultipleMajor <-
#          optimalUnknownProfileCombination(sampleTibble = sampleTibble,
#                                           markerImbalances = markerImbalances,
#                                           H = knownPerpetrator[[1]],
#                                           potentialParentsList = potentialParentsList,
#                                           control = controlMultiplePopulation),
#      times = 1)
#  
#  save(optimalMultipleMajor, file = "../inst/extdata/multipleMajor.RData")
#  save(multiplePopulationBenchmarkMajor, file = "../inst/extdata/multipleMajorBench.RData")

## ----multipleMajor2, echo = FALSE----------------------------------------
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

optimalMultipleMajor <- load_rdata(system.file('extdata', "multipleMajor.RData", package = 'MPSMixtures'))
multiplePopulationBenchmarkMajor <- load_rdata(system.file('extdata', "multipleMajorBench.RData", package = 'MPSMixtures'))

## ----multipleMajorParameters---------------------------------------------
optimalMultipleMajor$U[[1]][c("LogLikelihoods", "Fitness", "Parameters")]

## ------------------------------------------------------------------------
knownVictim <- setHypothesis(sampleTibble, numberOfContributors, trueProfiles["V"], theta)

## ----multipleMinor, eval = FALSE-----------------------------------------
#  optimalMultipleMinor <- optimalUnknownProfileCombination(sampleTibble = sampleTibble,
#                                                           markerImbalances = markerImbalances,
#                                                           H = knownVictim[[1]],
#                                                           potentialParentsList = potentialParentsList,
#                                                           control = controlMultiplePopulation)

## ----multipleMinor2, echo = FALSE----------------------------------------
optimalMultipleMinor <- load_rdata(system.file('extdata', "multipleMinor.RData", package = 'MPSMixtures'))

## ----multipleMinorParameters---------------------------------------------
optimalMultipleMinor$U[[1]][c("Parameters", "LogLikelihoods", "Fitness")]

## ----plottingDifferences-------------------------------------------------
optimalMinor <- optimalMultipleMinor$U[[1]]
optimalVSTruePerpetratorProfile <- 
    ggplotComparingProfiles(sampleTibble, knownProfilesIndividual, 
                            optimalMinor, TRUE, c("True", "Optimal"))

## ----figDifference, echo = FALSE, fig.align = "center", fig.cap = "\\label{fig:minor}The coverage against allele, for both profiles, shown in coloumns, and each marker, shown in rows, where the two profiles differ. Profiles one and two correspond to the optimal profiles under the model and the true profiles, respectively. The red and blue colouring corresonds to the major and minor contributors, respectively. Lastely, the black bars are regions attributed to the noise distribution.", fig.width = 8, fig.height = 6, fig.pos="ht!"----
show(optimalVSTruePerpetratorProfile)

## ----multipleUnknown, eval = FALSE---------------------------------------
#  bothUnknown <- setHypothesis(sampleTibble, numberOfContributors, list(), theta)
#  
#  controlMultiplePopulation <-
#      optimalUnknownProfileCombination.control(
#          numberOfPopulations = 64,
#          numberOfIterations = 100,
#          populationSize = 300,
#          numberOfFittestIndividuals = 1,
#          numberOfIterationsEqualMinMax = 50,
#          hillClimbingIterations = 0,
#          mutationDecayRate = 1,
#          parentSelectionWindowSize = 12,
#          trace = FALSE,
#          tolerance = rep(1e-8, 4)
#      )
#  
#  multiplePopulationBenchmarkUnknown <- microbenchmark(
#      optimalMultipleUnknown <-
#          optimalUnknownProfileCombination(sampleTibble = sampleTibble,
#                                           markerImbalances = markerImbalances,
#                                           H = bothUnknown[[1]],
#                                           potentialParentsList = potentialParentsList,
#                                           control = controlMultiplePopulation),
#      times = 1)

## ---- echo = FALSE-------------------------------------------------------
bothUnknown <- setHypothesis(sampleTibble, numberOfContributors, list(), theta)

controlMultiplePopulation <- 
    optimalUnknownProfileCombination.control(
        numberOfPopulations = 64, 
        numberOfIterations = 100,
        populationSize = 300,
        numberOfFittestIndividuals = 1,
        numberOfIterationsEqualMinMax = 50,
        hillClimbingIterations = 0,
        mutationDecayRate = 1,
        parentSelectionWindowSize = 12,
        trace = FALSE, 
        tolerance = rep(1e-8, 4)
    )

optimalMultipleUnknown <- load_rdata(system.file('extdata', "multipleUnknown.RData", package = 'MPSMixtures'))
multiplePopulationBenchmarkUnknown <- load_rdata(system.file('extdata', "multipleUnknownBench.RData", package = 'MPSMixtures'))

## ----multipleUnknownParameteres------------------------------------------
optimalMultipleUnknown$U[[1]][c("Parameters", "LogLikelihoods", "Fitness")]

## ----hypotheses----------------------------------------------------------
theta = 0
knownProfilesHp <- trueProfiles
knownProfilesHd <- trueProfiles["V"]

Hp <- setHypothesis(sampleTibble, numberOfContributors, knownProfilesHp, theta)
Hd <- setHypothesis(sampleTibble, numberOfContributors, knownProfilesHd, theta)

## ----multipleLRMinor, eval = FALSE---------------------------------------
#  LRParallelPopulations <-
#      LR(sampleTibble = sampleTibble, Hp = Hp, Hd = Hd, markerImbalances = markerImbalances,
#         potentialParentsList = potentialParentsList, stutterRatioModel = NULL,
#         control = optimalUnknownProfileCombination.control(
#             numberOfPopulations = 32, numberOfIterations = 150,
#             populationSize = 75, numberOfFittestIndividuals = 1000,
#             hillClimbingIterations = 0, mutationDecayRate = 1,
#             parentSelectionWindowSize = 6, simplifiedReturn = FALSE,
#             allowParentSurvival = TRUE, trace = FALSE,
#             tolerance = rep(1e-8, 4)))

## ----multipleLRMinor2----------------------------------------------------
LRParallelPopulations <- load_rdata(system.file('extdata', "multipleLRMinor.RData", package = 'MPSMixtures'))

## ----multipleLRMinorTable------------------------------------------------
LRParallelPopulations$ComparisonTable

