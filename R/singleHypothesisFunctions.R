#' Estimate parameters
#'
#' @details Estimates the allele and noise parameters in cases where all contributors are assumed known.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param markerImbalances A vector of the prior marker imbalances.
#' @param knownProfilesList A list of tibbles containing the alleles of the known contributors.
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample.
#' @param noiseParameters A vector of the prior noise parameters.
#' @param stutterRatioModel An \link{lm}-object modelling the stutter ratio given the BLMM. Note: this is only used if 'potentialParentsList' is NULL, but even then it can be NULL itself, implying no stuttering in the sample.
#' @param levelsOfStutterRecursion The number of layers used in the stutter recursion.
#' @param convexMarkerImbalanceInterpolation A fraction used to create a convex combination of the of MoM and the prior estimates of the marker imbalances.
#' @param tolerance The tolerance used for termination in parameter estimation.
#' @param numberOfThreads The number of threads passed to the \link{potentialParentsMultiCore}-function if applicable.
#'
#' @return A list of estimated parameters.
#' @export
estimateParametersOfKnownProfiles <- function(sampleTibble, markerImbalances, knownProfilesList, potentialParentsList, noiseParameters = NULL,
                                              stutterRatioModel = NULL, levelsOfStutterRecursion = 2, convexMarkerImbalanceInterpolation = 0.8,
                                              tolerance = 1e-6, numberOfThreads = 4) {
    if (is.null(potentialParentsList)) {
        potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel, numberOfThreads)
    }

    dualEstimation <- FALSE
    if (is.null(noiseParameters)) {
        dualEstimation <- TRUE
        noiseParameters <- c()
    }

    H <- setHypothesis(sampleTibble, length(knownProfilesList), knownProfilesList, 0)[[1]]

    numberOfMarkers = dim(sampleTibble %>% distinct(Marker))[1]
    numberOfAlleles = (sampleTibble %>% group_by(Marker) %>% summarise(Count = n(), .groups = "drop"))$Count

    numberOfContributors = H$NumberOfContributors
    numberOfKnownContributors = H$NumberOfKnownProfiles

    creatingIndividualObject <- .setupIndividual(numberOfMarkers, numberOfAlleles,
                                                 numberOfContributors, numberOfKnownContributors, H$KnownProfiles,
                                                 sampleTibble$Coverage, potentialParentsList, markerImbalances,
                                                 convexMarkerImbalanceInterpolation, noiseParameters,
                                                 tolerance, H$ThetaCorrection, sampleTibble$AlleleFrequencies,
                                                 levelsOfStutterRecursion, dualEstimation)

    return(creatingIndividualObject)
}


#' @title 'optimalUnknownProfileCombination' control function
#'
#' @description A function setting all relavent parameters used internally in the \link{LR}-function.
#'
#' @param numberOfPopulations The number of sub-populations.
#' @param populationSize The size of the sub-populations.
#' @param numberOfIterations The maximum number of (outer) iterations.
#' @param numberOfInnerIterations The number of inner iterations.
#' @param numberOfIterationsEqualMinMax The number of iterations with the number of sub-populations with the same 'largest' individual divided 'numberOfPopulations' less than or equal to 'fractionOfPopulationsMax'.
#' @param fractionOfPopulationsMax Fraction of the number of sub-populations with the same maximum needed before convergence counter initiates.
#' @param numberOfFittestIndividuals The number of unique individuals stored and returned.
#' @param parentSelectionWindowSize The size of the parent selecetion window.
#' @param allowParentSurvival Should parents be allowed to survive from iteration to iteration (TRUE/FALSE).
#' @param crossoverProbability The cross-over probability.
#' @param mutationProbabilityLowerLimit The lower limit on the probability of mutation.
#' @param mutationProbabilityUpperLimit The upper limit on the probability of mutation.
#' @param mutationIterations The number of iterations used in mutation stage (not currently in use).
#' @param mutationDegreesOfFreedom The degrees of freedom of the t-distribution used to create the mutation probabilities.
#' @param mutationDecayRate The rate of the decay of the mutation probability. Needed if 'mutationDecay' is 'NULL'.
#' @param mutationDecay The decay of the mutation probability to the lower limit.
#' @param fractionFittestIndividuals The fraction of the population forced to survive (not currently in use)
#' @param hillClimbingIterations The number of iterations each child or parent is hill climbed.
#' @param convexMarkerImbalanceInterpolation A fraction used to create a convex combination of the of MoM and the prior estimates of the marker imbalances.
#' @param tolerance Tolerance of internal log-likelihood maximisation (can be vector of upto size '4').
#' @param seed A seed for the c++ implementaion (deprecate).
#' @param trace Show trace (TRUE/FALSE)?
#' @param traceLimit Limits the trace of the parallel implementation.
#' @param simplifiedReturn Should the returned list be simplified (TRUE/FALSE)?
#' @param numberOfThreads The maximum number of threads allowed.
#' @param levelsOfStutterRecursion The number of layers used in the stutter recursion.
#' @param numberOfSimulationsMH The number of simulations used for each marker of the Metropolis-Hastings algorithm.
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
optimalUnknownProfileCombination.control <- function(numberOfPopulations = 4, populationSize = 10, numberOfIterations = 25, numberOfInnerIterations = 10,
                                                     numberOfIterationsEqualMinMax = 10, fractionOfPopulationsMax = NULL, numberOfFittestIndividuals = 10,
                                                     parentSelectionWindowSize = 5, allowParentSurvival = TRUE, crossoverProbability = NULL, mutationProbabilityLowerLimit = NULL, mutationProbabilityUpperLimit = 0.95,
                                                     mutationIterations = 2, mutationDegreesOfFreedom = 100,
                                                     mutationDecayRate = 1 / 2, mutationDecay = NULL, fractionFittestIndividuals = 1, hillClimbingIterations = 2,
                                                     convexMarkerImbalanceInterpolation = 0.8, tolerance = 1e-6, seed = NULL, trace = TRUE, simplifiedReturn = FALSE, numberOfThreads = 4, levelsOfStutterRecursion = 2,
                                                     traceLimit = 100, numberOfSimulationsMH = 10000) {

    if (numberOfPopulations == 1) {
        numberOfFittestIndividuals <- if (is.null(numberOfFittestIndividuals)) ceiling(0.1 * populationSize) else min(numberOfFittestIndividuals, populationSize)
    }
    else {
        numberOfFittestIndividuals <- if (is.null(numberOfFittestIndividuals)) ceiling(0.1 * populationSize) else numberOfFittestIndividuals
    }

    fractionOfPopulationsMax <- if (is.null(fractionOfPopulationsMax)) max(c(0.05, 1 / numberOfPopulations)) else fractionOfPopulationsMax

    seed <- if(is.null(seed)) sample(1e6, 1) else seed

    if (length(tolerance) == 1) {
        tolerance = c(tolerance, -rep(1e-4, 3))
    }
    else if (length(tolerance) == 2) {
        tolerance = c(tolerance, -1e-4, -1e-4)
    }
    else if (length(tolerance) == 3) {
        tolerance = c(tolerance, -1e-4)
    }

    if (length(tolerance) != 4) {
        tolerance = c(1e-8, -rep(1e-4, 3))
    }

    if (is.null(mutationDecayRate)) {
        mutationDecayRate = 0
    }

    controlList <- list(numberOfPopulations = numberOfPopulations, populationSize = populationSize, numberOfIterations = numberOfIterations,
                        numberOfInnerIterations = numberOfInnerIterations, numberOfIterationsEqualMinMax = numberOfIterationsEqualMinMax,
                        fractionOfPopulationsMax = fractionOfPopulationsMax,
                        numberOfFittestIndividuals = numberOfFittestIndividuals, parentSelectionWindowSize = parentSelectionWindowSize, allowParentSurvival = allowParentSurvival, crossoverProbability = crossoverProbability,
                        mutationProbabilityLowerLimit = mutationProbabilityLowerLimit, mutationProbabilityUpperLimit = mutationProbabilityUpperLimit,
                        mutationIterations = mutationIterations, mutationDegreesOfFreedom = mutationDegreesOfFreedom, mutationDecayRate = mutationDecayRate,
                        mutationDecay = mutationDecay, fractionFittestIndividuals = fractionFittestIndividuals, hillClimbingIterations = hillClimbingIterations,
                        convexMarkerImbalanceInterpolation = convexMarkerImbalanceInterpolation, tolerance = tolerance, seed = seed, trace = trace, simplifiedReturn = simplifiedReturn, numberOfThreads = numberOfThreads, levelsOfStutterRecursion = levelsOfStutterRecursion,
                        traceLimit = traceLimit, numberOfSimulationsMH = numberOfSimulationsMH)
    return(controlList)
}

#' @title Optimal unknown profile combinations.
#'
#' @description Finds the set of unknown profile combinations, as well as parameters, maximising the log-likehood function.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param markerImbalances The prior marker imbalances.
#' @param H An hypothesis.
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample.
#' @param noiseParameters A vector of the prior noise parameters.
#' @param stutterRatioModel An \link{lm}-object modelling the stutter ratio given the BLMM. Note: this is only used if 'potentialParentsList' is NULL, but even then it can be NULL itself, implying no stuttering in the sample.
#' @param HArguments If 'H' is 'NULL' a list of arguments passed to the 'setHypothesis' function. It must contain 'numberOfContributors', 'knownProfilesList', and 'theta'.
#' @param control An \link{optimalUnknownProfileCombination.control} object.
#'
#' @return A list of the unknown profile combinations contributing the most to the probability of the evidence under the provided hypothesis. The size of the list is controlled by 'control$numberOfFittestIndividuals'.
#' @export
optimalUnknownProfileCombination <- function(sampleTibble, markerImbalances, H, potentialParentsList,
                                             noiseParameters = NULL, stutterRatioModel = NULL, HArguments = NULL,
                                             control = optimalUnknownProfileCombination.control()) {
    if (is.null(potentialParentsList)) {
        if (control$trace)
            cat("Building potential parents list.\n")

        potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel, control$numberOfThreads)
    }

    if (is.null(H)) {
        H <- do.call(setHypothesis, c(list(sampleTibble = sampleTibble), HArguments))[[1]]
    }

    dualEstimation <- FALSE
    if (is.null(noiseParameters)) {
        dualEstimation <- TRUE
        noiseParameters <- c()
    }

    optimalUnknownProfiles <- .optimalUnknownProfilesHi(sampleTibble, H, markerImbalances, noiseParameters, potentialParentsList,
                                                        H$KnownProfiles, dualEstimation, control)
    return(optimalUnknownProfiles)
}

####################

.oneStepApproximation <- function(optimalUnkownProfileCombinationList, sampleTibble, H, numberOfUnknownContributors, numberOfAlleles, partialSumAlleles, estimatedParameters,
                                  optimalCombinationIndex, levelsOfStutterRecursion, numberOfThreads, returnSimplified, potentialParentsList,
                                  EAApproximationType, numberOfSimulations, suggestionPrior) {

    optimalCombination <- optimalUnkownProfileCombinationList[[optimalCombinationIndex]]
    oneStep <- lapply(seq_along(numberOfAlleles), function(m) {
        cat("m:", m, "\n")
        markerIndices <- (partialSumAlleles[m] + 1):(partialSumAlleles[m] + numberOfAlleles[m])
        res <- .oneStepApproximationCpp(optimalCombination$EncodedUnknownProfiles[((2 * (m - 1) * numberOfUnknownContributors) + 1):(2 * m * numberOfUnknownContributors)],
                                        estimatedParameters$SampleParameters, estimatedParameters$NoiseParameters, estimatedParameters$MixtureParameters,
                                        sampleTibble$Coverage[markerIndices], estimatedParameters$MarkerImbalanceParameters[m], potentialParentsList[m],
                                        as.matrix(H$KnownProfiles[markerIndices, ]), as.matrix(H$KnownProfiles[markerIndices, ]),
                                        sampleTibble$AlleleFrequencies[markerIndices], H$Theta, H$NumberOfContributors,
                                        1, numberOfAlleles[m], levelsOfStutterRecursion)

        sorting <- order(res$NormalisedProbabilities, decreasing = T)
        res <- lapply(res, function(xx) {
            if (length(xx) == 1)
                return(xx)
            else
                xx[sorting]
        })

        return(res)
    })

    res <- oneStep
    if (returnSimplified) {
        res <- sapply(oneStep, function(xx) xx$LogUnnormalisedProbabilitySum)
    }

    return(res)
}

.EAApproximation <- function(optimalUnkownProfileCombinationList, sampleTibble, H, numberOfUnknownContributors, numberOfAlleles, partialSumAlleles, estimatedParameters,
                             optimalCombinationIndex, levelsOfStutterRecursion, numberOfThreads, simplifiedReturn, potentialParentsList,
                             EAApproximationType, numberOfSimulations, suggestionPrior) {

    if (EAApproximationType == 1) {
        optimalUnkownProfileCombinationList_ = lapply(optimalUnkownProfileCombinationList, function(xx) as.matrix(xx$EncodedUnknownProfiles))
        res <- .EAApproximationCpp(optimalUnkownProfileCombinationList_,
                                   estimatedParameters$SampleParameters, estimatedParameters$NoiseParameters, estimatedParameters$MixtureParameters,
                                   sampleTibble$Coverage, estimatedParameters$MarkerImbalanceParameters, potentialParentsList,
                                   H$KnownProfiles, H$KnownProfiles, sampleTibble$AlleleFrequencies, H$Theta, H$NumberOfContributors,
                                   length(numberOfAlleles), numberOfAlleles, levelsOfStutterRecursion, TRUE)
    }
    else if (EAApproximationType == 2) {
        EAApproximation_m <- mclapply(seq_along(numberOfAlleles), function(m) {
            markerIndices <- (partialSumAlleles[m] + 1):(partialSumAlleles[m] + numberOfAlleles[m])
            optimalUnkownProfileCombinationList_m <- unique(lapply(optimalUnkownProfileCombinationList, function(xx) as.matrix(xx$EncodedUnknownProfiles[((2 * (m - 1) * numberOfUnknownContributors) + 1):(2 * m * numberOfUnknownContributors)])))

            res <- .EAApproximationCpp(optimalUnkownProfileCombinationList_m,
                                       estimatedParameters$SampleParameters, estimatedParameters$NoiseParameters, estimatedParameters$MixtureParameters,
                                       sampleTibble$Coverage[markerIndices], estimatedParameters$MarkerImbalanceParameters[m], potentialParentsList[m],
                                       as.matrix(H$KnownProfiles[markerIndices, ]), as.matrix(H$KnownProfiles[markerIndices, ]),
                                       sampleTibble$AlleleFrequencies[markerIndices], H$Theta, H$NumberOfContributors,
                                       1, numberOfAlleles[m], levelsOfStutterRecursion, FALSE)

            sorting <- order(res$NormalisedProbabilities, decreasing = T)
            res <- lapply(res, function(xx) {
                if (length(xx) == 1)
                    return(xx)
                else
                    xx[sorting]
            })

            return(res)
        }, mc.cores = numberOfThreads)

        res <- EAApproximation_m
        if (simplifiedReturn) {
            res <- sapply(EAApproximation_m, function(xx) xx$LogUnnormalisedProbabilitySum)
        }
    }

    return(res)
}


.samplePosteriorGenotypes <- function(optimalUnkownProfileCombinationList, sampleTibble, H,
                                      numberOfUnknownContributors, numberOfAlleles, partialSumAlleles, estimatedParameters,
                                      optimalCombinationIndex, levelsOfStutterRecursion, numberOfThreads, simplifiedReturn, potentialParentsList,
                                      EAApproximationType, numberOfSimulationsMH, suggestionMH) {
    suggestionBool = switch(tolower(suggestionMH),
                            "guided" = TRUE,
                            "random" = FALSE)

    optimalCombination <- optimalUnkownProfileCombinationList[[optimalCombinationIndex]]
    sampledPosterior <- mclapply(seq_along(numberOfAlleles), function(m) {
        markerIndices <- (partialSumAlleles[m] + 1):(partialSumAlleles[m] + numberOfAlleles[m])
        sampledGenotypesAll <- .samplePosteriorGenotypesGuidedCpp(optimalCombination$EncodedUnknownProfiles[((2 * (m - 1) * numberOfUnknownContributors) + 1):(2 * m * numberOfUnknownContributors)],
                                                                                estimatedParameters$SampleParameters, estimatedParameters$NoiseParameters,
                                                                                estimatedParameters$MixtureParameters, estimatedParameters$MarkerImbalanceParameters[m],
                                                                                sampleTibble$Coverage[markerIndices], potentialParentsList[m],
                                                                                as.matrix(H$KnownProfiles[markerIndices, ]), as.matrix(H$KnownProfiles[markerIndices, ]),
                                                                                sampleTibble$AlleleFrequencies[markerIndices], H$Theta, H$NumberOfContributors,
                                                                                numberOfAlleles[m], levelsOfStutterRecursion,
                                                                                numberOfSimulationsMH, suggestionBool, sample(1e6, 1))

        sampledGenotypes = sampledGenotypesAll[["SampledGenotypes"]]
        sampledLogLikelihoods = sampledGenotypesAll[["UnnormalisedLogLikelihood"]]

        uniqueGenotypes = unique(sampledGenotypes)
        countUniqueGenotypes = rep(NA, length(uniqueGenotypes))
        unnormalisedLogLikelihood = rep(NA, length(uniqueGenotypes))
        for (i in seq_along(uniqueGenotypes)) {
            identifiedDuplicates <- sapply(sampledGenotypes, function(xx) sum(xx != uniqueGenotypes[[i]]) == 0)
            unnormalisedLogLikelihood[i] = mean(sampledLogLikelihoods[identifiedDuplicates])
            countUniqueGenotypes[i] = sum(identifiedDuplicates)
        }

        posteriorProb = countUniqueGenotypes / sum(countUniqueGenotypes)
        sorted <- order(posteriorProb, decreasing = T)

        res <- list(GenotypeMatrix = uniqueGenotypes[sorted], LogUnnormalisedProbabilitySum = log(sum(exp(unnormalisedLogLikelihood))),
                    LogUnnormalisedProbability = unnormalisedLogLikelihood[sorted], NormalisedProbabilities = posteriorProb[sorted],
                    AcceptedProposals = sampledGenotypesAll[["AcceptedProposals"]],
                    AcceptRate = sampledGenotypesAll[["AcceptedProposals"]] / numberOfSimulationsMH)
        return(res)
    }, mc.cores = numberOfThreads)

    res <- sampledPosterior
    return(res)
}

#' Control function for the 'approximationSetUnknownGenotypeCombinations' function.
#'
#' @description A control function setting default parameters for the \link{approximationSetUnknownGenotypeCombinations} function.
#'
#' @param numberOfThreads The maximum number of threads allowed.
#' @param levelsOfStutterRecursion The number of layers used in the stutter recursion.
#' @param EAApproximationType The type of EA approximation used to find the normalising constant.
#' @param numberOfSimulationsMH The number of simulations used for each marker of the Metropolis-Hastings algorithm.
#' @param suggestionMH The proposal distribution used in the Metropolis-Hastings algorithm.
#' @param trace Show trace (TRUE/FALSE)?
#' @param simplifiedReturn Should the returned list be simplified (TRUE/FALSE)?
#'
#' @return A list of default parameters.
approximationSetUnknownGenotypeCombinations.control <- function(simplifiedReturn = FALSE, levelsOfStutterRecursion = 2, numberOfThreads = 4, trace = FALSE,
                                                                EAApproximationType = 2, numberOfSimulationsMH = 10000, suggestionMH = "random") {
    res = list(simplifiedReturn = simplifiedReturn, levelsOfStutterRecursion = levelsOfStutterRecursion, numberOfThreads = numberOfThreads,
               trace = trace, EAApproximationType = EAApproximationType, numberOfSimulationsMH = numberOfSimulationsMH, suggestionMH = suggestionMH)
    return(res)
}


#' Approximations for the set of combinations of unknown contributors
#'
#' @description Three approximation to the sum over the set of unknown genotype combinations using: (1) a simple one-step approximation, (2) the set of the fittest individuals, and (3) a Metropolis-Hastings sampler.
#'
#' @param optimalUnkownProfileCombinationList An object created by the 'optimalUnknownProfileCombination' function.
#' @param method The approximation method, either: 'EA', 'oneStep', or 'MH'.
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param H A hypothesis.
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample.
#' @param stutterRatioModel An \link{lm}-object modelling the stutter ratio given the BLMM. Note: this is only used if 'potentialParentsList' is NULL, but even then it can be NULL itself, implying no stuttering in the sample.
#' @param HArguments If 'H' is 'NULL' a list of arguments passed to the 'setHypothesis' function. It must contain 'numberOfContributors', 'knownProfilesList', and 'theta'.
#' @param control An \link{optimalUnknownProfileCombination.control} object.
#'
#' @return A list containing the unnormalised probabilitis for every marker, the corresponding genotypes, and approximations to the posterior probabilities for every unknown genotype combination for every marker.
#' @export
approximationSetUnknownGenotypeCombinations <- function(optimalUnkownProfileCombinationList, method = 'EA', sampleTibble,
                                                        H = NULL, potentialParentsList, stutterRatioModel = NULL, HArguments = NULL,
                                                        control = approximationSetUnknownGenotypeCombinations.control()) {
    method = tolower(method)
    approximationMethod <- switch (method,
                                   'ea' = .EAApproximation,
                                   'onestep' = .oneStepApproximation,
                                   'mh' = .samplePosteriorGenotypes
    )

    if (is.null(potentialParentsList)) {
        if (control$trace)
            cat("Building potential parents list.\n")

        potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel, control$numberOfThreads)
    }

    if (is.null(H)) {
        H <- do.call(setHypothesis, c(list(sampleTibble = sampleTibble), HArguments))[[1]]
    }

    optimalCombinationIndex <- which.max(sapply(optimalUnkownProfileCombinationList, function(xx) xx$Fitness))
    optimalCombination <- optimalUnkownProfileCombinationList[[optimalCombinationIndex]]

    estimatedParameters <- optimalCombination$Parameters
    if (length(optimalUnkownProfileCombinationList) > 1) {
        normalisingConstant <- optimalCombination$Fitness
        estimatedParameters <- suppressWarnings(.optimiseParametersLargeLikelihood(sampleTibble, optimalUnkownProfileCombinationList,
                                                                                   H$NumberOfContributors, normalisingConstant))
    }

    numberOfAlleles = (sampleTibble %>% group_by_(~Marker) %>% summarise(N = n()))$N
    partialSumAlleles = .partialSumEigen(numberOfAlleles)
    numberOfUnknownContributors = H$NumberOfContributors - H$NumberOfKnownProfiles

    res <- approximationMethod(optimalUnkownProfileCombinationList, sampleTibble, H, numberOfUnknownContributors, numberOfAlleles, partialSumAlleles,
                               estimatedParameters, optimalCombinationIndex, control$levelsOfStutterRecursion, control$numberOfThreads,
                               control$simplifiedReturn, potentialParentsList, control$EAApproximationType, control$numberOfSimulationsMH, control$suggestionMH)

    class(res) <- "setOfUnknownGenotypes"
    return(res)
}

#' Head of 'setOfUnknownGenotypes'-object
#'
#' @description Returns the most prevalent elements of 'setOfUnknownGenotypes'-object for each entry in the list.
#'
#' @param x A list of unknown genotypes, normalised, and unnormalised posterior probabilities.
#' @param ... Additional arguments: The number of elements, 'n', and the elements of the list to be returned, 'outputElements'.
#' @return A reduced 'setOfUnknownGenotypesList'-object.
head.setOfUnknownGenotypes <- function(x, ...) {
    setOfUnknownGenotypesList <- x
    res <- vector("list", length(setOfUnknownGenotypesList))
    argsList <- list(...)

    if (is.null(argsList$n)) {
        n = 6L
    }
    else {
        n = argsList$n
    }

    if (n < 0) {
        argsList$setOfUnknownGenotypesList = setOfUnknownGenotypesList
        argsList$n <- -n
        return(do.call("tail", argsList))
    }

    if (is.null(argsList$outputElements)) {
        outputElements <- 1:length(setOfUnknownGenotypesList[[1]])
        argsList$outputElements <- outputElements
    }
    else {
        outputElements <- which(names(setOfUnknownGenotypesList[[1]]) %in% argsList$outputElements)
    }

    if (length(argsList$outputElements) != length(outputElements)) {
        elementsNotFound <- argsList$outputElements[-which(argsList$outputElements %in% names(setOfUnknownGenotypesList[[1]]))]
        objList <- ifelse(length(elementsNotFound) == 1, paste0(elementsNotFound[length(elementsNotFound)], "."),
                          paste0(paste(elementsNotFound[-length(elementsNotFound)], collapse = ", " ), ", and ", elementsNotFound[length(elementsNotFound)], "."))

        warning(paste0("The following requested elements were not found in the 'setOfUnknownGenotypesList' object: ", objList))
    }

    for (m in seq_along(setOfUnknownGenotypesList)) {
        res[[m]] <- lapply(setOfUnknownGenotypesList[[m]][outputElements], function(xx) xx[1:min(length(xx), n)])
    }

    class(res) <- "setOfUnknownGenotypesList"
    return(res)
}

#' Tail of 'setOfUnknownGenotypes'-object
#'
#' @description Returns the least prevalent elements of 'setOfUnknownGenotypes'-object for each entry in the list.
#'
#' @param x A list of unknown genotypes, normalised, and unnormalised posterior probabilities.
#' @param ... Additional arguments: The number of elements, 'n', and the elements of the list to be returned, 'outputElements'.
#'
#' @return A reduced 'setOfUnknownGenotypesList'-object.
tail.setOfUnknownGenotypes <- function(x, ...) {
    setOfUnknownGenotypesList <- x
    res <- vector("list", length(setOfUnknownGenotypesList[[1]]))
    argsList <- list(...)

    if (is.null(argsList$n)) {
        n = 6L
    }
    else {
        n = argsList$n
    }

    if (n < 0) {
        argsList$setOfUnknownGenotypesList = setOfUnknownGenotypesList
        argsList$n <- -n
        return(do.call("head", argsList))
    }

    if (is.null(argsList$outputElements)) {
        outputElements <- 1:length(setOfUnknownGenotypesList[[1]])
    }
    else {
        outputElements <- which(names(setOfUnknownGenotypesList[[1]]) %in% argsList$outputElements)
    }

    if (length(argsList$outputElements) != length(outputElements)) {
        elementsNotFound <- argsList$outputElements[-which(argsList$outputElements %in% names(setOfUnknownGenotypesList[[1]]))]
        objList <- ifelse(length(elementsNotFound) == 1, paste0(elementsNotFound[length(elementsNotFound)], "."),
                          paste0(paste(elementsNotFound[-length(elementsNotFound)], collapse = ", " ), ", and ", elementsNotFound[length(elementsNotFound)], "."))

        warning(paste0("The following requested elements were not found in the 'setOfUnknownGenotypesList' object: ", objList))
    }

    for (m in seq_along(setOfUnknownGenotypesList)) {
        res[[m]] <- lapply(setOfUnknownGenotypesList[[m]][outputElements], function(xx) xx[max(1, length(xx) - n):length(xx)])
    }

    class(res) <- "setOfUnknownGenotypesList"
    return(res)
}

