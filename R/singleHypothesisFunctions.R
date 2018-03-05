#' Estimate parameters
#'
#' @details Estimates the allele and noise parameters in cases where all contirbutors are assumed known.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param knownProfilesList A list of tibbles containing the alleles of the known contributors.
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample.
#' @param stutterRatioModel An \link{lm}-object modelling the stutter ratio given the BLMM. Note: this is only used if 'potentialParentsList' is NULL, but even then it can be NULL itself, implying no stuttering in the sample.
#' @param tolerance The tolerance used for termination in parameter estimation.
#'
#' @return A list of estimated parameters.
#' @export
estimateParametersOfKnownProfiles <- function(sampleTibble, markerImbalance, knownProfilesList, potentialParentsList, stutterRatioModel = NULL, levelsOfStutterRecursion = 2, convexMarkerImbalanceInterpolation = 0.8, tolerance = 1e-6) {
    if (is.null(potentialParentsList)) {
        potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel, control$numberOfThreads)
    }

    H <- setHypothesis(sampleTibble, length(knownProfilesList), knownProfilesList, 0)[[1]]

    numberOfMarkers = dim(sampleTibble %>% distinct(Marker))[1]
    numberOfAlleles = (sampleTibble %>% group_by(Marker) %>% summarise(Count = n()))$Count

    numberOfContributors = H$NumberOfContributors
    numberOfKnownContributors = H$NumberOfKnownProfiles

    creatingIndividualObject <- MPSMixtures:::.setupIndividual(numberOfMarkers, numberOfAlleles,
                                                 numberOfContributors, numberOfKnownContributors, H$KnownProfiles,
                                                 sampleTibble$Coverage, potentialParentsList, markerImbalance,
                                                 convexMarkerImbalanceInterpolation,
                                                 tolerance, H$ThetaCorrection, sampleTibble$AlleleFrequencies, levelsOfStutterRecursion)

    return(creatingIndividualObject)
}

#' Control function for the 'optimalUnknownProfileCombination' function
#'
#' @details A control function setting default parameters for the \link{optimalUnknownProfileCombination} function. It takes the same parameters as the \link{LR.control} function, with the expection of 'simplifiedReturn'.
#'
#' @return A list of default parameters.
#' @export
optimalUnknownProfileCombination.control <- function(numberOfPopulations = 4, populationSize = 10, numberOfIterations = 25, numberOfInnerIterations = 10,
                                                     numberOfIterationsEqualMinMax = 10, fractionOfPopulationsMax = NULL, numberOfFittestIndividuals = 10,
                                                     parentSelectionWindowSize = 5, allowParentSurvival = TRUE, crossoverProbability = NULL, mutationProbabilityLowerLimit = NULL, mutationDegreesOfFreedom = 100,
                                                     mutationDecayRate = 2, mutationDecay = NULL, fractionFittestIndividuals = 1, hillClimbingDirections = 1, hillClimbingIterations = 1,
                                                     convexMarkerImbalanceInterpolation = 0.8, tolerance = 1e-6, seed = NULL, trace = TRUE, numberOfThreads = 4, levelsOfStutterRecursion = 2) {
    controlList <- LR.control(numberOfPopulations, populationSize, numberOfIterations, numberOfInnerIterations, numberOfIterationsEqualMinMax, fractionOfPopulationsMax, numberOfFittestIndividuals,
                              parentSelectionWindowSize, allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, mutationDegreesOfFreedom, mutationDecayRate,
                              mutationDecay, fractionFittestIndividuals, hillClimbingDirections, hillClimbingIterations, convexMarkerImbalanceInterpolation, tolerance, seed, trace, FALSE, numberOfThreads, levelsOfStutterRecursion)
    return(controlList)
}

#' @title Optimal unknown profile combinations.
#'
#' @description Finds the set of unknown profile combinations, as well as parameters, maximising the log-likehood function.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param numberOfContributors The total number of contributors to the mixture.
#' @param knownProfilesList A list of tibbles containing the alleles of the known contributors.
#' @param theta The inbreeding coefficient (Fst).
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample.
#' @param stutterRatioModel An \link{lm}-object modelling the stutter ratio given the BLMM. Note: this is only used if 'potentialParentsList' is NULL, but even then it can be NULL itself, implying no stuttering in the sample.
#' @param control An \link{LR.control} object.
#'
#' @return A list of the unknown profile combinations contributing the most to the probability of the evidence under the provided hypothesis. The size of the list is controlled by 'control$numberOfFittestIndividuals'.
#' @export
optimalUnknownProfileCombination <- function(sampleTibble, markerImbalances, numberOfContributors, knownProfilesList, theta,
                                             potentialParentsList, stutterRatioModel = NULL, control = optimalUnknownProfileCombination.control()) {
    if (is.null(potentialParentsList)) {
        if (control$trace)
            cat("Building potential parents list.\n")

        potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel, control$numberOfThreads)
    }

    H <- setHypothesis(sampleTibble, numberOfContributors, knownProfilesList, theta)[[1]]
    optimalUnknownProfiles <- .optimalUnknownProfilesHi(sampleTibble, H, markerImbalances, potentialParentsList, H$KnownProfiles, control)

    return(optimalUnknownProfiles)
}



sumUnknownsKStepApproximation <- function(optimalUnkownProfileCombinationList, sampleTibble, H, k) {
    optimalCombinationIndex <- which.max(sapply(optimalUnkownProfileCombinationList, function(xx) xx$Fitness))
    optimalCombination <- optimalUnkownProfileCombinationList[[optimalCombinationIndex]]

    optimalCombinationLogLikelihoodAlleles <- optimalCombination$LogLikelihoodAlleleCoverage
    normalisingConstant <- optimalCombination$Fitness

    estimatedParameters <- optimalCombination$Parameters
    if (length(optimalUnkownProfileCombinationList) > 1) {
        estimatedParameters <- .optimiseParametersLargeLikelihood(sampleTibble, optimalUnkownProfileCombinationList,
                                                                  H$NumberOfContributors, normalisingConstant)
    }

    numberOfAlleles = (sampleTibble %>% group_by(Marker) %>% summarise(N = n()))$N
    numberOfUknownContributors = H$NumberOfContributors - H$NumberOfKnownProfiles

    encodedGenotypeSplit <- split(optimalCombination$EncodedUnknownProfiles, rep(1:length(numberOfAlleles), each = 2 * (numberOfUknownContributors)))
    ecmList <- vector("list", length(encodedGenotypeSplit))
    for (m in seq_along(encodedGenotypeSplit)) {
        ecmList[[m]] <- .
    }



    return(0)
}
