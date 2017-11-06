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
        migratedPopulation[[migrationForward[i]]] <- .migration(currentPopulationList[[i]], migratedPopulation[[migrationForward[i]]], migrant_i)
        migratedPopulation[[migrationBackward[i]]] <- .migration(currentPopulationList[[i]], migratedPopulation[[migrationBackward[i]]], migrant_i)
    }

    return(migratedPopulation)
}

# numberOfMaxThreads = control$numberOfMaxThreads; numberOfInnerIterations = control$numberOfInnerIterations;

.runningParallelPopulationEvolutionaryAlgorithm <- function(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                                                           coverage, potentialParentsList, markerImbalances, tolerance, theta, alleleFrequencies,
                                                           numberOfPopulations, populationSize, numberOfIterations, numberOfInnerIterations,
                                                           numberOfIterationsEqualMax, fractionOfPopulationsMax,
                                                           numberOfFittestIndividuals,
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
            PEA_i <- .runningParallelEvolutionaryAlgorithm(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                                                        coverage, potentialParentsList, markerImbalances, tolerance, theta, alleleFrequencies,
                                                        numberOfInnerIterations, numberOfInnerIterations, populationSize,
                                                        parentSelectionWindowSize, allowParentSurvival, crossoverProbability, mutationProbabilityLowerLimit, mutationDegreesOfFreedom,
                                                        mutationDecaySplit[[j + 1]], hillClimbingDirections, hillClimbingIterations, sample(1:1e6, 1), FALSE,
                                                        currentPopulationMigratedList[[i]]$EncodedProfiles, currentPopulationMigratedList[[i]]$SampleParameters,
                                                        currentPopulationMigratedList[[i]]$NoiseParameters, currentPopulationMigratedList[[i]]$MixtureParameters)
            return(PEA_i)
        }, mc.cores = numberOfThreads)

        currentPopulationList <- lapply(newPopulation, function(np) np$LastPopulation)
        currentFittestList <- unlist(lapply(newPopulation, function(np) np$FittestIndividuals), recursive = FALSE)

        ## Updating list of fittest individuals
        fittestIndividuals <- append(fittestIndividuals, currentFittestList)
        duplicateIndividuals <- fittestIndividuals[!duplicated(lapply(fittestIndividuals, function(cfl) cfl$EncodedUnknownProfiles))]
        fittestIndividuals <- fittestIndividuals[order(unlist(lapply(fittestIndividuals, function(fi) fi$Fitness)), decreasing = TRUE)[1:(min(c(numberOfFittestIndividuals, length(fittestIndividuals))))]]

        ## Updating convergence condition
        populationFitness <- do.call("c", lapply(currentPopulationList, function(cpl) cpl$Fitness))
        uniqueSubpopulationMaxima <- unique.matrix(do.call("cbind", lapply(currentPopulationList, function(cpl) cpl$EncodedProfiles[, which.max(cpl$Fitness)])), MARGIN = 2)
        numberOfUniqueSubpopulationMaxima = dim(uniqueSubpopulationMaxima)[2]

        if ((numberOfUniqueSubpopulationMaxima / numberOfPopulations) <= fractionOfPopulationsMax) {
            k = k + 1
        }
        else {
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
                "\t\tSubpopulations /w unique maxima:", numberOfUniqueSubpopulationMaxima, "\n",
                "\t\tTermination counter:", k, "/", numberOfIterationsEqualMax, "\n")
    }

    return(fittestIndividuals)
}
