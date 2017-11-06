#' @title Set a hypothesis.
#'
#' @description Set a hypothesis given the sample, the number of total contributors, a list of the known contributors, and population parameters and data.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param numberOfContributors The total number of contributors to the mixture. Note: can be a vector of possible hypotheses, but elements should always be larger than or equal to the number of known profiles.
#' @param knownProfilesList A list of tibbles containing the alleles of the known contributors.
#' @param theta The inbreeding coefficient (Fst).
#' @param alleleFrequencies The allele frequencies of the population.
#'
#' @return A list containing the information relavent to the hypothesis.
#' @export
setHypothesis <- function(sampleTibble, numberOfContributors, knownProfilesList, theta, alleleFrequencies) {
    knownGenotypeMatrix <- genotypeMatrix(sampleTibble, knownProfilesList)
    numberOfKnownContributors <- length(knownProfilesList)

    res <- lapply(numberOfContributors[which(numberOfContributors >= numberOfKnownContributors)], function(ii) {
        list(NumberOfContributors = ii, NumberOfKnownProfiles = numberOfKnownContributors,
             KnownProfiles = knownGenotypeMatrix, ThetaCorrection = theta, AlleleFrequencies = alleleFrequencies)
    })

    return(res)
}

# knownProfiles = H$KnownProfiles; coverage = sampleTibble$Coverage; markerImbalances = sampleTibble$MarkerImbalance; tolerance = control$tolerance; theta = H$ThetaCorrection; alleleFrequencies = H$AlleleFrequencies; numberOfPopulations = control$numberOfPopulations; populationSize = control$populationSize; numberOfIterations = control$numberOfIterations; numberOfIterationsEqualMinMax = control$numberOfIterationsEqualMinMax; numberOfFittestIndividuals = control$numberOfFittestIndividuals; parentSelectionWindowSize = control$parentSelectionWindowSize; allowParentSurvival = control$allowParentSurvival; mutationDegreesOfFreedom = control$mutationDegreesOfFreedom; mutationDecay = control$mutationDecay; fractionFittestIndividuals = control$fractionFittestIndividuals; hillClimbingDirections = control$hillClimbingDirections; hillClimbingIterations = control$hillClimbingIterations; simplifiedReturn = FALSE; seed = control$seed; trace = control$trace
# H = Hd[[i]]

#' @title Optimal unknown profile combinations under H_i
#'
#' @description Finds the set of unknown profile combinations contributing the most to the probability of the evidence under H_i. Note: this is mostly for use internally in the \link{LR}-function.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param H A hypothesis (see the \link{setHypothesis}-function for the general structure).
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample.
#' @param control An \link{LR.control} object.
#'
#' @return A list of the unknown profile combinations contributing the most to the probability of the evidence under the provided hypothesis. The size of the list is controlled by control$numberOfFittestIndividuals.
#' @export
optimalUnknownProfilesHi <- function(sampleTibble, H, potentialParentsList, allKnownProfiles, control) {
    numberOfMarkers = dim(sampleTibble %>% distinct(Marker))[1]
    numberOfAlleles = (sampleTibble %>% group_by(Marker) %>% summarise(Count = n()))$Count

    numberOfContributors = H$NumberOfContributors
    numberOfKnownContributors = H$NumberOfKnownProfiles

    if ((numberOfContributors - numberOfKnownContributors) == 0) {
        creatingIndividualObject <- .setupIndividual(numberOfMarkers, numberOfAlleles,
                                                    numberOfContributors, numberOfKnownContributors, H$KnownProfiles,
                                                    sampleTibble$Coverage, potentialParentsList, sampleTibble$MarkerImbalance,
                                                    control$tolerance, H$ThetaCorrection, H$AlleleFrequencies)

        optimalUnknownProfiles <- list(creatingIndividualObject)
    }
    else {
        crossoverProbability <- ifelse(is.null(control$crossoverProbability), 1 / (2 * (numberOfContributors - numberOfKnownContributors) * numberOfMarkers), control$crossoverProbability)
        mutationProbabilityLowerLimit <- ifelse(is.null(control$mutationProbabilityLowerLimit), 1 / (2 * (numberOfContributors - numberOfKnownContributors) * numberOfMarkers), control$mutationProbabilityLowerLimit)

        if (control$numberOfPopulations == 1) {
            optimalUnknownProfiles <- .runningSinglePopulationEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, H$KnownProfiles, allKnownProfiles,
                                                                                   sampleTibble$Coverage, potentialParentsList, sampleTibble$MarkerImbalance, control$tolerance, H$ThetaCorrection, H$AlleleFrequencies,
                                                                                   control$populationSize, control$numberOfIterations, control$numberOfIterationsEqualMinMax, control$numberOfFittestIndividuals,
                                                                                   control$parentSelectionWindowSize, control$allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, control$mutationDegreesOfFreedom,
                                                                                   control$mutationDecay, control$hillClimbingDirections, control$hillClimbingIterations,
                                                                                   control$seed, control$trace)
        }
        else {
            optimalUnknownProfiles <- .runningParallelPopulationEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, H$KnownProfiles, allKnownProfiles,
                                                                                     sampleTibble$Coverage, potentialParentsList, sampleTibble$MarkerImbalance, control$tolerance, H$ThetaCorrection, H$AlleleFrequencies,
                                                                                     control$numberOfPopulations, control$populationSize, control$numberOfIterations, control$numberOfInnerIterations,
                                                                                     control$numberOfIterationsEqualMinMax, control$fractionOfPopulationsMax, control$numberOfFittestIndividuals,
                                                                                     control$parentSelectionWindowSize, control$allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, control$mutationDegreesOfFreedom,
                                                                                     control$mutationDecay, control$hillClimbingDirections, control$hillClimbingIterations,
                                                                                     control$seed, control$trace, control$numberOfThreads)
        }
    }

    return(optimalUnknownProfiles)
}

#' @title LR control function
#'
#' @description A function setting all relavent parameters used internally in the \link{LR}-function.
#'
#' @param numberOfPopulations The number of sub-populations.
#' @param populationSize The size of the sub-populations.
#' @param numberOfIterations The maximum number of (outer) iterations.
#' @param numberOfInnerIterations The number of inner iterations.
#' @param numberOfIterationsEqualMax The number of iterations with the number of sub-populations with the same 'largest' individual divided 'numberOfPopulations' less than or equal to 'fractionOfPopulationsMax'.
#' @param fractionOfPopulationsMax Fraction of the number of sub-populations with the same maximum needed before convergence counter initiates.
#' @param numberOfFittestIndividuals The number of unique individuals stored and returned.
#' @param parentSelectionWindowSize The size of the parent selecetion window.
#' @param allowParentSurvival Should parents be allowed to survive from iteration to iteration (TRUE/FALSE).
#' @param crossoverProbability The cross-over probability.
#' @param mutationProbabilityLowerLimit The lower limit to the mutation probability.
#' @param mutationDegreesOfFreedom The degrees of freedom of the t-distribution used to create the mutation probabilities.
#' @param mutationDecay The decay of the mutation probability to the lower limit.
#' @param hillClimbingDirections The number of hill climbing directions (not in use at the moment).
#' @param hillClimbingIterations The number of iterations each child or parent is hill climbed.
#' @param tolerance Tolerance of internal log-likelihood maximisation.
#' @param seed A seed for the c++ implementaion.
#' @param trace Show trace (TRUE/FALSE)?
#' @param simplifiedReturn Should the returned list be simplified (TRUE/FALSE)?
#' @param numberOfThreads The maximum number of threads allowed.
#'
#' @details
#' If 'numberOfFittestIndividuals' is 'NULL' (default), then the number of fittest individuals stored across iterations is ceiling(0.1 * 'populationSize'). If the number of populations is larger than 1, then 'populationSize' refers to the sub-population size. Thus, it is recommended this parameter is set if running multiple populations.
#'
#' If 'fractionOfPopulationsMax' is 'NULL' (default), it is set as the maximum of 0.05 and 1 / 'numberOfPopulations'.
#'
#' If 'mutationDecay' is 'NULL' (default) and if the number of populations is 1, then it is set as seq(0, by = mutationDecayRate * 4 / (numberOfIterations), length.out = numberOfIterations). If the population size is larger than 1 the parameter 'numberOfIterations' is replaced by 'numberOfIterations * numberOfInnerIterations' i.e. the total number of iterations.
#'
#' If 'seed' is 'NULL', then it is sampled from uniformly from [1; 10^6].
#'
#' @return A list containing all relevant parameters.
LR.control <- function(numberOfPopulations = 4, populationSize = 10, numberOfIterations = 25, numberOfInnerIterations = 10,
                       numberOfIterationsEqualMinMax = 10, fractionOfPopulationsMax = NULL, numberOfFittestIndividuals = 10,
                       parentSelectionWindowSize = 5, allowParentSurvival = TRUE, crossoverProbability = NULL, mutationProbabilityLowerLimit = NULL, mutationDegreesOfFreedom = 100,
                       mutationDecayRate = 2, mutationDecay = NULL, fractionFittestIndividuals = 1, hillClimbingDirections = 1, hillClimbingIterations = 1,
                       tolerance = 1e-4, seed = NULL, trace = TRUE, simplifiedReturn = FALSE, numberOfThreads = 4) {

    numberOfFittestIndividuals <- if (is.null(numberOfFittestIndividuals)) ceiling(0.1 * populationSize) else numberOfFittestIndividuals
    fractionOfPopulationsMax <- if (is.null(fractionOfPopulationsMax)) max(c(0.05, 1 / numberOfPopulations)) else fractionOfPopulationsMax

    if (is.null(mutationDecay)) {
        if (numberOfPopulations == 1) {
            mutationDecay <- seq(0, by = mutationDecayRate * 4 / (numberOfIterations), length.out = numberOfIterations)
        }
        else {
            mutationDecay <- seq(0, by = mutationDecayRate * 4 / (numberOfIterations * numberOfInnerIterations), length.out = numberOfIterations * numberOfInnerIterations)
        }
    }

    seed <- if(is.null(seed)) sample(1e6, 1) else seed

    controlList <- list(numberOfPopulations = numberOfPopulations, populationSize = populationSize, numberOfIterations = numberOfIterations,
                        numberOfInnerIterations = numberOfInnerIterations, numberOfIterationsEqualMinMax = numberOfIterationsEqualMinMax,
                        fractionOfPopulationsMax = fractionOfPopulationsMax,
                        numberOfFittestIndividuals = numberOfFittestIndividuals, parentSelectionWindowSize = parentSelectionWindowSize, allowParentSurvival = allowParentSurvival, crossoverProbability = crossoverProbability,
                        mutationProbabilityLowerLimit = mutationProbabilityLowerLimit, mutationDegreesOfFreedom = mutationDegreesOfFreedom, mutationDecayRate = mutationDecayRate,
                        mutationDecay = mutationDecay, fractionFittestIndividuals = 1.0, hillClimbingDirections = hillClimbingDirections, hillClimbingIterations = hillClimbingIterations,
                        tolerance = tolerance, seed = seed, trace = trace, simplifiedReturn = simplifiedReturn, numberOfThreads = numberOfThreads)
    return(controlList)
}

# sampleTibble = coverageTibble; potentialParentsList = rrpp; control = LR.control(numberOfPopulations = 4, numberOfIterations = 25, populationSize = 10, numberOfFittestIndividuals = 100, hillClimbingIterations = 2, parentSelectionWindowSize = 2, simplifiedReturn = F, allowParentSurvival = T)

#' @title Likelihood ratio
#'
#' @description Calculates an approximate likelihood ratio given two competing hypotheses (in the future lists of competing hypotheses).
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param Hp The prosecutor hypothesis (see \link{setHypothesis}).
#' @param Hd The defence hypothesis.
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample. If NULL then a stutterRatioModel should be provided.
#' @param stutterRatioModel A linear model of class \link{lm} modelling the relationship between coverage and stutter. Only needed if the potential parents list is not provided.
#' @param control An \link{LR.control} object.
#'
#' @return A list of likelihood ratios comparing the two hypotheses (always calculated as Hp / Hd).
#' @export
LR <- function(sampleTibble, Hp, Hd, potentialParentsList, stutterRatioModel = NULL, control = LR.control()) {
    ## Set-up
    if (is.null(potentialParentsList)) {
        if (control$trace)
            cat("Building potential parents list.\n")

        potentialParentsTibble <- potentialParents(sampleTibble, stutterRatioModel, trace = control$trace)
        potentialParentsList <- lapply(potentialParentsTibble, function(ppp) lapply(ppp, function(pppp) {
            pppp %>% select(PotentialParent, StutterRatio) %>% as.matrix()
        }))
    }

    allKnownProfiles = unique(do.call("cbind", lapply(append(Hp, Hd), function(H) H$KnownProfiles)), MARGIN = 2)

    if (control$trace)
        cat("Running Hp-list.\n")

    optimalUnknownGenotypesHp <- vector("list", length(Hp))
    for (i in seq_along(Hp)) {
        optimalUnknownGenotypesHp[[i]] <- optimalUnknownProfilesHi(sampleTibble, Hp[[i]], potentialParentsList, allKnownProfiles, control)
    }

    if (control$trace)
        cat("Running Hd-list.\n")

    optimalUnknownGenotypesHd <- vector("list", length(Hd))
    for (i in seq_along(Hd)) {
        optimalUnknownGenotypesHd[[i]] <- optimalUnknownProfilesHi(sampleTibble, Hd[[i]], potentialParentsList, allKnownProfiles, control)
    }

    if (control$trace)
        cat("Optimising parameters and calculating LR's.\n")

    allPairwiseCombinations <- expand.grid(Hp = seq_along(optimalUnknownGenotypesHp), Hd = seq_along(optimalUnknownGenotypesHd))
    pairwiseComparisonResults <- vector("list", dim(allPairwiseCombinations)[1])
    for (i in 1:dim(allPairwiseCombinations)[1]) {
        if (control$trace)
            cat("  Combination:", paste0("Hp = ", allPairwiseCombinations[i, 1], ", Hd = ", allPairwiseCombinations[i, 2]), "\n")

        LHpNormaliser <- max(sapply(optimalUnknownGenotypesHp[[allPairwiseCombinations[i, 1]]], function(hh) hh$Fitness))
        LHdNormaliser <- max(sapply(optimalUnknownGenotypesHd[[allPairwiseCombinations[i, 2]]], function(hh) hh$Fitness))

        parametersLRHp <- optimiseParametersLargeLikelihood(sampleTibble, optimalUnknownGenotypesHp[[allPairwiseCombinations[i, 1]]],
                                                            Hp[[allPairwiseCombinations[i, 1]]]$NumberOfContributors, LHpNormaliser)
        parametersLRHd <- optimiseParametersLargeLikelihood(sampleTibble, optimalUnknownGenotypesHd[[allPairwiseCombinations[i, 2]]],
                                                            Hd[[allPairwiseCombinations[i, 2]]]$NumberOfContributors, LHdNormaliser)

        logLR = log(parametersLRHp$Likelihood) + LHpNormaliser - log(parametersLRHd$Likelihood) - LHdNormaliser

        resultsList <- list(LR = exp(logLR), Log10LR = logLR * log10(exp(1)), Hp = Hp[[allPairwiseCombinations[i, 1]]], Hd = Hd[[allPairwiseCombinations[i, 1]]])
        if (!control$simplifiedReturn) {
            resultsList$HpOptimalUnknownGenotypes = optimalUnknownGenotypesHp[[allPairwiseCombinations[i, 1]]]
            resultsList$HdOptimalUnknownGenotypes = optimalUnknownGenotypesHd[[allPairwiseCombinations[i, 2]]]

            resultsList$HpParameterEstimates <- parametersLRHp
            resultsList$HdParameterEstimates <- parametersLRHd
        }

        pairwiseComparisonResults[[i]] <- resultsList
    }


    resList <- list(AllPairwiseComparisonData =pairwiseComparisonResults)

    comparisonTable <- expand.grid(Hp = sapply(Hp, function(H) H$NumberOfContributors), Hd = sapply(Hd, function(H) H$NumberOfContributors))
    comparisonTable$Log10LR <- sapply(pairwiseComparisonResults, function(L) L$Log10LR)

    resList$ComparisonTable <- comparisonTable
    return(resList)
}
