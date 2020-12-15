#' @title Set a hypothesis.
#'
#' @description Set a hypothesis given the sample, the number of total contributors, a list of the known contributors, and population parameters and data.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param numberOfContributors The total number of contributors to the mixture. Note: can be a vector of possible hypotheses, but elements should always be larger than or equal to the number of known profiles.
#' @param knownProfilesList A list of tibbles containing the alleles of the known contributors.
#' @param theta The inbreeding coefficient (Fst).
#'
#' @return A list containing the information relevant to the hypothesis.
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

.optimalUnknownProfilesHi <- function(sampleTibble, H, markerImbalances, noiseParameters, potentialParentsList,
                                      allKnownProfiles, dualEstimation, control) {
    numberOfMarkers = dim(sampleTibble %>% distinct(Marker))[1]
    numberOfAlleles = (sampleTibble %>% group_by(Marker) %>% summarise(Count = n(), .groups = "drop"))$Count

    numberOfContributors = H$NumberOfContributors
    numberOfKnownContributors = H$NumberOfKnownProfiles

    if ((numberOfContributors - numberOfKnownContributors) == 0) {
        creatingIndividualObject <- .setupIndividual(numberOfMarkers, numberOfAlleles,
                                                     numberOfContributors, numberOfKnownContributors, H$KnownProfiles,
                                                     sampleTibble$Coverage, potentialParentsList, markerImbalances, control$convexMarkerImbalanceInterpolation,
                                                     noiseParameters, control$tolerance, H$ThetaCorrection,
                                                     sampleTibble$AlleleFrequencies, control$levelsOfStutterRecursion,
                                                     dualEstimation)

        optimalUnknownProfiles <- list(U = list(creatingIndividualObject))
    }
    else {
        crossoverProbability <- ifelse(is.null(control$crossoverProbability), 1 / (2 * (numberOfContributors - numberOfKnownContributors) * numberOfMarkers), control$crossoverProbability)
        mutationProbabilityLowerLimit <- ifelse(is.null(control$mutationProbabilityLowerLimit), 1 / (2 * (numberOfContributors - numberOfKnownContributors) * numberOfMarkers), control$mutationProbabilityLowerLimit)
        mutationProbabilityUpperLimit <- ifelse(is.null(control$mutationProbabilityUpperLimit), 1 / (2 * (numberOfContributors - numberOfKnownContributors) * numberOfMarkers), control$mutationProbabilityUpperLimit)

        mutationDecay = control$mutationDecay
        if (control$numberOfPopulations == 1) {
            if (is.null(mutationDecay)) {
                mutationDecay <- c(seq(mutationProbabilityUpperLimit, mutationProbabilityLowerLimit, length.out = floor(control$mutationDecayRate / control$numberOfIterations)), rep(mutationProbabilityLowerLimit, times = control$numberOfIterations - floor(control$mutationDecayRate / control$numberOfIterations)))
            }

            optimalUnknownProfiles <- .runningSinglePopulationEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, H$KnownProfiles, allKnownProfiles,
                                                                                    sampleTibble$Coverage, potentialParentsList, markerImbalances, control$convexMarkerImbalanceInterpolation, noiseParameters,
                                                                                    control$tolerance, H$ThetaCorrection, sampleTibble$AlleleFrequencies,
                                                                                    control$populationSize, control$numberOfIterations, control$numberOfIterationsEqualMinMax, control$numberOfFittestIndividuals,
                                                                                    control$parentSelectionWindowSize, control$allowParentSurvival, control$fractionFittestIndividuals,
                                                                                    crossoverProbability, mutationProbabilityLowerLimit, control$mutationIterations,
                                                                                    control$mutationDegreesOfFreedom, mutationDecay, control$hillClimbingIterations,
                                                                                    control$seed, control$trace, control$levelsOfStutterRecursion, dualEstimation)

            optimalUnknownProfiles$U <- optimalUnknownProfiles$U[order(sapply(optimalUnknownProfiles$U, function(oup) oup$Fitness), decreasing = TRUE)]
        }
        else {
            if (is.null(mutationDecay)) {
                mutationDecay <- c(seq(mutationProbabilityUpperLimit, mutationProbabilityLowerLimit, length.out = floor(control$mutationDecayRate / (control$numberOfIterations * control$numberOfInnerIterations))), rep(mutationProbabilityLowerLimit, times = control$numberOfIterations * control$numberOfInnerIterations - floor(control$mutationDecayRate / (control$numberOfIterations * control$numberOfInnerIterations))))
            }

            optimalUnknownProfiles <- .runningParallelPopulationEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, H$KnownProfiles, allKnownProfiles,
                                                                                      sampleTibble$Coverage, potentialParentsList, markerImbalances, control$convexMarkerImbalanceInterpolation,
                                                                                      noiseParameters, control$tolerance, H$ThetaCorrection, sampleTibble$AlleleFrequencies,
                                                                                      control$numberOfPopulations, control$populationSize, control$numberOfIterations, control$numberOfInnerIterations,
                                                                                      control$numberOfIterationsEqualMinMax, control$fractionOfPopulationsMax, control$numberOfFittestIndividuals,
                                                                                      control$parentSelectionWindowSize, control$allowParentSurvival, control$fractionFittestIndividuals,
                                                                                      crossoverProbability, mutationProbabilityLowerLimit, control$mutationIterations, control$mutationDegreesOfFreedom,
                                                                                      mutationDecay, control$hillClimbingIterations, control$seed, control$trace, control$numberOfThreads,
                                                                                      control$levelsOfStutterRecursion, dualEstimation, control$traceLimit)
        }
    }

    return(optimalUnknownProfiles)
}


#' @title Likelihood ratio
#'
#' @description Calculates an approximate likelihood ratio given two competing hypotheses (in the future lists of competing hypotheses).
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param Hp The prosecutor hypothesis (see \link{setHypothesis}).
#' @param Hd The defence hypothesis.
#' @param markerImbalances A vector of prior marker imbalances.
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample. If NULL then a 'stutterRatioModel' should be provided.
#' @param noiseParameters A vector of the prior noise parameters.
#' @param stutterRatioModel A linear model of class \link{lm} modelling the relationship between coverage and stutter. Only needed if the potential parents list is not provided.
#' @param control An \link{optimalUnknownProfileCombination.control} object.
#'
#' @return A list of likelihood ratios comparing the two hypotheses (always calculated as Hp / Hd).
#' @export
LR <- function(sampleTibble, Hp, Hd, markerImbalances = NULL, potentialParentsList = NULL, noiseParameters = NULL, stutterRatioModel = NULL, control = optimalUnknownProfileCombination.control()) {
    ## Set-up
    if (is.null(potentialParentsList)) {
        if (control$trace)
            cat("Building potential parents list.\n")

        potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel, control$numberOfThreads)
    }

    numberOfMarkers = length(unique(sampleTibble$Marker))
    if (is.null(markerImbalances)) {
        markerImbalances = rep(1, numberOfMarkers)
    } else if (numberOfMarkers != length(markerImbalances)) {
        stop("The length of 'markerImbalances' is not equal to the number of unique markers in 'sampleTibble'.")
    }

    dualEstimation <- FALSE
    if (is.null(noiseParameters)) {
        dualEstimation <- TRUE
        noiseParameters <- c()
    }

    allKnownProfiles = unique(do.call("cbind", lapply(append(Hp, Hd), function(H) H$KnownProfiles)), MARGIN = 2)

    if (control$trace)
        cat("Running Hp-list.\n")

    optimalUnknownGenotypesHp <- vector("list", length(Hp))
    for (i in seq_along(Hp)) {
        optimalUnknownGenotypesHp[[i]] <- .optimalUnknownProfilesHi(sampleTibble, Hp[[i]], markerImbalances, noiseParameters, potentialParentsList, allKnownProfiles, dualEstimation, control)$U
    }

    if (control$trace)
        cat("Running Hd-list.\n")

    optimalUnknownGenotypesHd <- vector("list", length(Hd))
    for (i in seq_along(Hd)) {
        optimalUnknownGenotypesHd[[i]] <- .optimalUnknownProfilesHi(sampleTibble, Hd[[i]], markerImbalances, noiseParameters, potentialParentsList, allKnownProfiles, dualEstimation, control)$U
    }

    if (control$trace)
        cat("Optimising parameters and calculating LR's.\n")

    pairwiseComparisons <- expand.grid("Hp" = seq_along(Hp), "Hd" = seq_along(Hd))
    pairwiseComparisonResults <- vector("list", dim(pairwiseComparisons)[1])
    for (i in 1:dim(pairwiseComparisons)[1]) {
        if (control$trace)
            cat("  Comparison:", i, "\n")

        pairwiseComparisons_i <- as.numeric(pairwiseComparisons[i, ])
        optimalUnknownGenotypesHp_i <- optimalUnknownGenotypesHp[[pairwiseComparisons_i[1]]][sapply(optimalUnknownGenotypesHp[[pairwiseComparisons_i[1]]], function(hh) !is.null(hh$Fitness))]
        optimalUnknownGenotypesHd_i <- optimalUnknownGenotypesHd[[pairwiseComparisons_i[2]]][sapply(optimalUnknownGenotypesHd[[pairwiseComparisons_i[2]]], function(hh) !is.null(hh$Fitness))]

        LHpNormaliser <- max(sapply(optimalUnknownGenotypesHp_i, function(hh) hh$Fitness))
        LHdNormaliser <- max(sapply(optimalUnknownGenotypesHd_i, function(hh) hh$Fitness))

        parametersLRHp <- .optimiseParametersLargeLikelihood(sampleTibble, optimalUnknownGenotypesHp_i,
                                                             Hp[[pairwiseComparisons_i[1]]]$NumberOfContributors, LHpNormaliser)
        parametersLRHd <- .optimiseParametersLargeLikelihood(sampleTibble, optimalUnknownGenotypesHd_i,
                                                             Hd[[pairwiseComparisons_i[2]]]$NumberOfContributors, LHdNormaliser)

        logLR = log(parametersLRHp$Likelihood) + LHpNormaliser - log(parametersLRHd$Likelihood) - LHdNormaliser

        resultsList <- list(LR = exp(logLR), Log10LR = logLR * log10(exp(1L)), Hp = Hp[[pairwiseComparisons_i[1]]], Hd = Hd[[pairwiseComparisons_i[2]]])
        if (!control$simplifiedReturn) {
            resultsList$HpOptimalUnknownGenotypes = approximationSetUnknownGenotypeCombinations(optimalUnknownGenotypesHp[[pairwiseComparisons_i[1]]], method = "mh",
                                                                                                sampleTibble, Hp[[pairwiseComparisons_i[1]]],
                                                                                                potentialParentsList, stutterRatioModel,
                                                                                                control = approximationSetUnknownGenotypeCombinations.control(
                                                                                                    levelsOfStutterRecursion = control$levelsOfStutterRecursion,
                                                                                                    numberOfThreads = control$numberOfThreads, trace = control$trace,
                                                                                                    numberOfSimulationsMH = control$numberOfSimulationsMH))

            resultsList$HdOptimalUnknownGenotypes = approximationSetUnknownGenotypeCombinations(optimalUnknownGenotypesHd[[pairwiseComparisons_i[2]]], method = "mh",
                                                                                                sampleTibble, Hd[[pairwiseComparisons_i[2]]],
                                                                                                potentialParentsList, stutterRatioModel,
                                                                                                control = approximationSetUnknownGenotypeCombinations.control(
                                                                                                    levelsOfStutterRecursion = control$levelsOfStutterRecursion,
                                                                                                    numberOfThreads = control$numberOfThreads, trace = control$trace,
                                                                                                    numberOfSimulationsMH = control$numberOfSimulationsMH))

            resultsList$HpParameterEstimates <- parametersLRHp
            resultsList$HdParameterEstimates <- parametersLRHd
        }

        pairwiseComparisonResults[[i]] <- resultsList
    }


    resList <- list(AllPairwiseComparisonData = pairwiseComparisonResults)

    comparisonTable <- data.frame(Hp = sapply(Hp, function(H) H$NumberOfContributors)[pairwiseComparisons[, 1]],
                                  Hd = sapply(Hd, function(H) H$NumberOfContributors)[pairwiseComparisons[, 2]])

    comparisonTable$Log10LR <- sapply(pairwiseComparisonResults, function(L) L$Log10LR)

    resList$ComparisonTable <- comparisonTable
    return(resList)
}
