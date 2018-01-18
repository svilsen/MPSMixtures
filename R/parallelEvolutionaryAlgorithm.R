# P1 = currentPopulationList[[1]]; P2 = currentPopulationList[[migrationForward[1]]]; migrant = migrant_i
.migration <- function(P1, P2, migrant, randomReplacement = FALSE) {
    if (randomReplacement) {
        replacedIndividual = sample(1:length(P2$Fitness), 1)
    }
    else {
        replacedIndividual = which.min(P2$Fitness)
    }

    P2$EncodedProfiles[, replacedIndividual] <- P1$EncodedProfiles[, migrant]
    P2$SampleParameters[, replacedIndividual] <- P1$SampleParameters[, migrant]
    P2$NoiseParameters[, replacedIndividual] <- P1$NoiseParameters[, migrant]
    P2$MixtureParameters[, replacedIndividual] <- P1$MixtureParameters[, migrant]
    P2$Fitness[replacedIndividual] <- P1$Fitness[migrant]

    return(P2)
}

.parallelPopulationEvolutionaryAlgorithmMigration <- function(currentPopulationList) {
    numberOfPopulations = length(currentPopulationList)
    migrationForward <- ((1:numberOfPopulations) %% numberOfPopulations) + 1
    migrationBackward <- ((1:numberOfPopulations - 3) %% numberOfPopulations) + 1

    migratedPopulation <- currentPopulationList
    for (i in seq_along(currentPopulationList)) {
        migrant_i = which.max(currentPopulationList[[i]]$Fitness)
        migratedPopulation[[migrationForward[i]]] <- MPSMixtures:::.migration(currentPopulationList[[i]], migratedPopulation[[migrationForward[i]]], migrant_i)
        migratedPopulation[[migrationBackward[i]]] <- MPSMixtures:::.migration(currentPopulationList[[i]], migratedPopulation[[migrationBackward[i]]], migrant_i)
    }

    return(migratedPopulation)
}

# control <- optimalUnknownProfileCombination.control( numberOfPopulations = 12, numberOfIterations = 100, populationSize = 10, numberOfFittestIndividuals = 100, hillClimbingIterations = 5, parentSelectionWindowSize = 2, allowParentSurvival = TRUE, trace = TRUE); knownProfilesList <- knownPerpetrator; H <- setHypothesis(sampleTibble, numberOfContributors, knownProfilesList, theta)[[1]]; numberOfMaxThreads = control$numberOfThreads; numberOfInnerIterations = control$numberOfInnerIterations; numberOfIterationsEqualMax = control$numberOfIterationsEqualMinMax; fractionOfPopulationsMax = control$fractionOfPopulationsMax; allKnownProfiles = H$KnownProfiles; coverage = sampleTibble$Coverage; markerImbalances = sampleTibble$MarkerImbalance; tolerance = control$tolerance; theta = H$ThetaCorrection; alleleFrequencies = sampleTibble$AlleleFrequencies; numberOfPopulations = control$numberOfPopulations; populationSize = control$populationSize; numberOfIterations = control$numberOfIterations; numberOfIterationsEqualMinMax = control$numberOfIterationsEqualMinMax; numberOfFittestIndividuals = control$numberOfFittestIndividuals; parentSelectionWindowSize = control$parentSelectionWindowSize; allowParentSurvival = control$allowParentSurvival; mutationDegreesOfFreedom = control$mutationDegreesOfFreedom; mutationDecay = control$mutationDecay; fractionFittestIndividuals = control$fractionFittestIndividuals; hillClimbingDirections = control$hillClimbingDirections; hillClimbingIterations = control$hillClimbingIterations; simplifiedReturn = FALSE; seed = control$seed; trace = control$trace; numberOfMarkers = dim(sampleTibble %>% distinct(Marker))[1]; numberOfAlleles = (sampleTibble %>% group_by(Marker) %>% summarise(Count = n()))$Count; numberOfContributors = H$NumberOfContributors; numberOfKnownContributors = H$NumberOfKnownProfiles; crossoverProbability <- ifelse(is.null(control$crossoverProbability), 1 / (2 * (numberOfContributors - numberOfKnownContributors) * numberOfMarkers), control$crossoverProbability); mutationProbabilityLowerLimit <- ifelse(is.null(control$mutationProbabilityLowerLimit), 1 / (2 * (numberOfContributors - numberOfKnownContributors) * numberOfMarkers), control$mutationProbabilityLowerLimit); knownProfiles <- H$KnownProfiles;
# knownProfiles = H$KnownProfiles; numberOfMaxThreads = 4; numberOfInnerIterations = 10; fractionOfPopulationsMax = control$fractionOfPopulationsMax; numberOfIterationsEqualMax = control$numberOfIterationsEqualMinMax; numberOfIterations = control$numberOfIterations

.runningParallelPopulationEvolutionaryAlgorithm <- function(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                                                            coverage, potentialParentsList, markerImbalances, tolerance, theta, alleleFrequencies,
                                                            numberOfPopulations, populationSize, numberOfIterations, numberOfInnerIterations,
                                                            numberOfIterationsEqualMax, fractionOfPopulationsMax, numberOfFittestIndividuals,
                                                            parentSelectionWindowSize, allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, mutationDegreesOfFreedom,
                                                            mutationDecay, hillClimbingDirections, hillClimbingIterations,
                                                            seed, trace, numberOfMaxThreads) {
    if (length(seed) == 1) {
        seed = sample(1:1e6, numberOfPopulations)
    }

    if (length(seed) != numberOfPopulations) {
        warning("All initial populations are equal.")
    }

    numberOfThreads <- min(numberOfPopulations, numberOfMaxThreads)
    mutationDecaySplit <- split(mutationDecay, rep(1:numberOfIterations, each = numberOfInnerIterations))
    currentPopulationList <- mclapply(1:numberOfPopulations, function(i) MPSMixtures:::.initialisingParallelEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                                                                                                                                  coverage, potentialParentsList, markerImbalances, tolerance, theta, alleleFrequencies, populationSize, seed[i]), mc.cores = numberOfThreads)

    j = 0
    k = 0
    converged = FALSE
    fittestIndividuals <- list()
    while (!converged) {
        ## Migration between subpopulations
        currentPopulationMigratedList <- MPSMixtures:::.parallelPopulationEvolutionaryAlgorithmMigration(currentPopulationList)

        ## Updating subpopulations
        newPopulation <- mclapply(1:numberOfPopulations, function(i) {
            PEA_i <- MPSMixtures:::.runningParallelEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                                                                         coverage, potentialParentsList, markerImbalances, tolerance, theta, alleleFrequencies,
                                                                         numberOfInnerIterations, numberOfInnerIterations, populationSize,
                                                                         parentSelectionWindowSize, allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, mutationDegreesOfFreedom,
                                                                         mutationDecaySplit[[j + 1]], hillClimbingDirections, hillClimbingIterations,
                                                                         sample(1:1e6, 1), FALSE,
                                                                         currentPopulationMigratedList[[i]]$EncodedProfiles, currentPopulationMigratedList[[i]]$SampleParameters,
                                                                         currentPopulationMigratedList[[i]]$NoiseParameters, currentPopulationMigratedList[[i]]$MixtureParameters,
                                                                         currentPopulationMigratedList[[i]]$Fitness)
            return(PEA_i)
        }, mc.cores = numberOfThreads)

        if (any(sapply(newPopulation, function(np) is.null(np$LastPopulation))))
            stop("Contact your local support staff: a thread returned a 'NULL' object indicating that the internal 'c++' code crashed.")

        currentPopulationList <- lapply(newPopulation, function(np) np$LastPopulation)
        currentFittestList <- unlist(lapply(newPopulation, function(np) np$FittestIndividuals), recursive = FALSE)

        ## Updating list of fittest individuals
        fittestIndividuals <- append(fittestIndividuals, currentFittestList)
        fittestIndividuals <- fittestIndividuals[!duplicated(lapply(fittestIndividuals, function(cfl) cfl$EncodedUnknownProfiles))]
        fittestIndividuals <- fittestIndividuals[order(unlist(lapply(fittestIndividuals, function(fi) fi$Fitness)), decreasing = TRUE)[1:(min(c(numberOfFittestIndividuals, length(fittestIndividuals))))]]

        ## Updating convergence condition
        populationFitness <- do.call("c", lapply(currentPopulationList, function(cpl) cpl$Fitness))
        uniqueSubpopulationMaxima <- unique.matrix(do.call("cbind", lapply(currentPopulationList, function(cpl) cpl$EncodedProfiles[, which.max(cpl$Fitness)])), MARGIN = 2)
        fractionOfUniqueSubpopulationMaxima = dim(uniqueSubpopulationMaxima)[2] / numberOfPopulations

        if (fractionOfUniqueSubpopulationMaxima <= fractionOfPopulationsMax) {
            k = k + 1
        } else {
            k = 0
        }

        j = j + 1
        converged = (k == numberOfIterationsEqualMax) | (j == numberOfIterations)

        if (trace)
            cat("\tOuter iteration:", j, "\n",
                "\t\tFitness:\n",
                "\t\t  Highest:", max(populationFitness), "\n",
                "\t\t  Average:", mean(populationFitness), "\n",
                "\t\t  Lowest:", min(populationFitness), "\n",
                "\t\tSubpopulations /w unique maxima:", fractionOfUniqueSubpopulationMaxima, "\n",
                "\t\tTermination counter:", k, "/", numberOfIterationsEqualMax, "\n")
    }

    return(fittestIndividuals)
}
