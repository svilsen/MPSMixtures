#' @title Set a hypothesis.
#'
#' @description Set a hypothesis given the sample, the number of total contributors, a list of the known contributors, and population parameters and data.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param numberOfContributors The total number of contributors to the mixture. Note: can be a vector of possible hypotheses, but elements should always be larger than or equal to the number of known profiles.
#' @param knownProfilesList A list of tibbles containing the alleles of the known contributors.
#' @param theta The inbreeding coefficient (Fst).
#'
#' @return A list containing the information relavent to the hypothesis.
#' @export
setHypothesis <- function(sampleTibble, numberOfContributors, knownProfilesList, theta) {
    knownGenotypeMatrix <- genotypeMatrix(sampleTibble, knownProfilesList)
    numberOfKnownContributors <- length(knownProfilesList)

    res <- lapply(numberOfContributors[which(numberOfContributors >= numberOfKnownContributors)], function(ii) {
        list(NumberOfContributors = ii, NumberOfKnownProfiles = numberOfKnownContributors,
             KnownProfiles = knownGenotypeMatrix, ThetaCorrection = theta)
    })

    return(res)
}

# allKnownProfiles = H$KnownProfiles; coverage = sampleTibble$Coverage; markerImbalances = sampleTibble$MarkerImbalance; tolerance = control$tolerance; theta = H$ThetaCorrection; alleleFrequencies = sampleTibble$AlleleFrequencies; numberOfPopulations = control$numberOfPopulations; populationSize = control$populationSize; numberOfIterations = control$numberOfIterations; numberOfIterationsEqualMinMax = control$numberOfIterationsEqualMinMax; numberOfFittestIndividuals = control$numberOfFittestIndividuals; parentSelectionWindowSize = control$parentSelectionWindowSize; allowParentSurvival = control$allowParentSurvival; mutationDegreesOfFreedom = control$mutationDegreesOfFreedom; mutationDecay = control$mutationDecay; fractionFittestIndividuals = control$fractionFittestIndividuals; hillClimbingDirections = control$hillClimbingDirections; hillClimbingIterations = control$hillClimbingIterations; simplifiedReturn = FALSE; seed = control$seed; trace = control$trace
# H = Hd[[i]]

.optimalUnknownProfilesHi <- function(sampleTibble, H, potentialParentsList, allKnownProfiles, control) {
    numberOfMarkers = dim(sampleTibble %>% distinct(Marker))[1]
    numberOfAlleles = (sampleTibble %>% group_by(Marker) %>% summarise(Count = n()))$Count

    numberOfContributors = H$NumberOfContributors
    numberOfKnownContributors = H$NumberOfKnownProfiles

    if ((numberOfContributors - numberOfKnownContributors) == 0) {
        creatingIndividualObject <- .setupIndividual(numberOfMarkers, numberOfAlleles,
                                                    numberOfContributors, numberOfKnownContributors, H$KnownProfiles,
                                                    sampleTibble$Coverage, potentialParentsList, sampleTibble$MarkerImbalance,
                                                    control$tolerance, H$ThetaCorrection, sampleTibble$AlleleFrequencies, control$levelsOfStutterRecursion)

        optimalUnknownProfiles <- list(creatingIndividualObject)
    }
    else {
        crossoverProbability <- ifelse(is.null(control$crossoverProbability), 1 / (2 * (numberOfContributors - numberOfKnownContributors) * numberOfMarkers), control$crossoverProbability)
        mutationProbabilityLowerLimit <- ifelse(is.null(control$mutationProbabilityLowerLimit), 1 / (2 * (numberOfContributors - numberOfKnownContributors) * numberOfMarkers), control$mutationProbabilityLowerLimit)

        if (control$numberOfPopulations == 1) {
            optimalUnknownProfiles <- MPSMixtures:::.runningSinglePopulationEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, H$KnownProfiles, allKnownProfiles,
                                                                                   sampleTibble$Coverage, potentialParentsList, sampleTibble$MarkerImbalance, control$tolerance, H$ThetaCorrection, sampleTibble$AlleleFrequencies,
                                                                                   control$populationSize, control$numberOfIterations, control$numberOfIterationsEqualMinMax, control$numberOfFittestIndividuals,
                                                                                   control$parentSelectionWindowSize, control$allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, control$mutationDegreesOfFreedom,
                                                                                   control$mutationDecay, control$hillClimbingDirections, control$hillClimbingIterations,
                                                                                   control$seed, control$trace, control$levelsOfStutterRecursion)

            optimalUnknownProfiles <- optimalUnknownProfiles[order(sapply(optimalUnknownProfiles, function(oup) oup$Fitness), decreasing = TRUE)]
        }
        else {
            optimalUnknownProfiles <- MPSMixtures:::.runningParallelPopulationEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, H$KnownProfiles, allKnownProfiles,
                                                                                     sampleTibble$Coverage, potentialParentsList, sampleTibble$MarkerImbalance, control$tolerance, H$ThetaCorrection, sampleTibble$AlleleFrequencies,
                                                                                     control$numberOfPopulations, control$populationSize, control$numberOfIterations, control$numberOfInnerIterations,
                                                                                     control$numberOfIterationsEqualMinMax, control$fractionOfPopulationsMax, control$numberOfFittestIndividuals,
                                                                                     control$parentSelectionWindowSize, control$allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, control$mutationDegreesOfFreedom,
                                                                                     control$mutationDecay, control$hillClimbingDirections, control$hillClimbingIterations,
                                                                                     control$seed, control$trace, control$numberOfThreads, control$levelsOfStutterRecursion)
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
                       tolerance = 1e-6, seed = NULL, trace = TRUE, simplifiedReturn = FALSE, numberOfThreads = 4, levelsOfStutterRecursion = 2) {

    if (numberOfPopulations == 1) {
        numberOfFittestIndividuals <- if (is.null(numberOfFittestIndividuals)) ceiling(0.1 * populationSize) else min(numberOfFittestIndividuals, populationSize)
    }
    else {
        numberOfFittestIndividuals <- if (is.null(numberOfFittestIndividuals)) ceiling(0.1 * populationSize) else numberOfFittestIndividuals
    }

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
                        tolerance = tolerance, seed = seed, trace = trace, simplifiedReturn = simplifiedReturn, numberOfThreads = numberOfThreads, levelsOfStutterRecursion = levelsOfStutterRecursion)
    return(controlList)
}

# control = LR.control(numberOfPopulations = 12, numberOfIterations = 25, populationSize = 10, numberOfFittestIndividuals = 100, hillClimbingIterations = 2, parentSelectionWindowSize = 2, simplifiedReturn = F, allowParentSurvival = T, trace = T)

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
LR <- function(sampleTibble, Hp, Hd, potentialParentsList = NULL, stutterRatioModel = NULL, control = LR.control()) {
    ## Set-up
    if (is.null(potentialParentsList)) {
        if (control$trace)
            cat("Building potential parents list.\n")

        potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel, control$numberOfThreads)
    }

    allKnownProfiles = unique(do.call("cbind", lapply(append(Hp, Hd), function(H) H$KnownProfiles)), MARGIN = 2)

    if (control$trace)
        cat("Running Hp-list.\n")

    if ((length(Hp) != length(Hd)) & ((length(Hp) != 1) | (length(Hd) != 1)))  {
        stop("'Hp' must have the same length as 'Hd', or either 'Hp' or 'Hd' must have length '1'.")
    }

    optimalUnknownGenotypesHp <- vector("list", length(Hp))
    for (i in seq_along(Hp)) {
        optimalUnknownGenotypesHp[[i]] <- MPSMixtures:::.optimalUnknownProfilesHi(sampleTibble, Hp[[i]], potentialParentsList, allKnownProfiles, control)
    }

    if ((length(Hp) == 1) & (length(Hd) != 1)) {
        Hp = rep(Hp[[i]], length(Hd))
        optimalUnknownGenotypesHp = rep(optimalUnknownGenotypesHp[[1]], length(Hd))
    }

    if (control$trace)
        cat("Running Hd-list.\n")

    optimalUnknownGenotypesHd <- vector("list", length(Hd))
    for (i in seq_along(Hd)) {
        optimalUnknownGenotypesHd[[i]] <- MPSMixtures:::.optimalUnknownProfilesHi(sampleTibble, Hd[[i]], potentialParentsList, allKnownProfiles, control)
    }

    if ((length(Hd) == 1) & (length(Hp) != 1)) {
        Hd = rep(Hd[[i]], length(Hp))
        optimalUnknownGenotypesHd = rep(optimalUnknownGenotypesHd[[1]], length(Hp))
    }

    if (control$trace)
        cat("Optimising parameters and calculating LR's.\n")


    pairwiseComparisonResults <- vector("list", length(Hp))
    for (i in seq_along(pairwiseComparisonResults)) {
        if (control$trace)
            cat("  Comparison:", i, "\n")

        optimalUnknownGenotypesHp_i <- optimalUnknownGenotypesHp[[i]][sapply(optimalUnknownGenotypesHp[[i]], function(hh) !is.null(hh$Fitness))]
        optimalUnknownGenotypesHd_i <- optimalUnknownGenotypesHd[[i]][sapply(optimalUnknownGenotypesHd[[i]], function(hh) !is.null(hh$Fitness))]

        LHpNormaliser <- max(sapply(optimalUnknownGenotypesHp_i, function(hh) hh$Fitness))
        LHdNormaliser <- max(sapply(optimalUnknownGenotypesHd_i, function(hh) hh$Fitness))

        parametersLRHp <- optimiseParametersLargeLikelihood(sampleTibble, optimalUnknownGenotypesHp_i,
                                                            Hp[[i]]$NumberOfContributors, LHpNormaliser)
        parametersLRHd <- optimiseParametersLargeLikelihood(sampleTibble, optimalUnknownGenotypesHd_i,
                                                            Hd[[i]]$NumberOfContributors, LHdNormaliser)

        logLR = log(parametersLRHp$Likelihood) + LHpNormaliser - log(parametersLRHd$Likelihood) - LHdNormaliser

        resultsList <- list(LR = exp(logLR), Log10LR = logLR * log10(exp(1)), Hp = Hp[[i]], Hd = Hd[[i]])
        if (!control$simplifiedReturn) {
            resultsList$HpOptimalUnknownGenotypes = optimalUnknownGenotypesHp[[i]]
            resultsList$HdOptimalUnknownGenotypes = optimalUnknownGenotypesHd[[i]]

            resultsList$HpParameterEstimates <- parametersLRHp
            resultsList$HdParameterEstimates <- parametersLRHd
        }

        pairwiseComparisonResults[[i]] <- resultsList
    }


    resList <- list(AllPairwiseComparisonData = pairwiseComparisonResults)

    comparisonTable <- data.frame(Hp = sapply(Hp, function(H) H$NumberOfContributors), Hd = sapply(Hd, function(H) H$NumberOfContributors))
    comparisonTable$Log10LR <- sapply(pairwiseComparisonResults, function(L) L$Log10LR)

    resList$ComparisonTable <- comparisonTable
    return(resList)
}
