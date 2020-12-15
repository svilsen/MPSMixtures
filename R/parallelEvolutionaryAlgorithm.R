
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
    P2$MarkerImbalanceParameters[, replacedIndividual] <- P1$MarkerImbalanceParameters[, migrant]
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

.runningParallelPopulationEvolutionaryAlgorithm <- function(numberOfMarkers, numberOfAlleles, numberOfContributors, numberOfKnownContributors, knownProfiles, allKnownProfiles,
                                                            coverage, potentialParentsList, markerImbalances, convexMarkerImbalanceInterpolation, noiseParameters, tolerance, theta, alleleFrequencies,
                                                            numberOfPopulations, populationSize, numberOfIterations, numberOfInnerIterations,
                                                            numberOfIterationsEqualMax, fractionOfPopulationsMax, numberOfFittestIndividuals,
                                                            parentSelectionWindowSize, allowParentSurvival, fractionFittestIndividuals, crossoverProbability,
                                                            mutationProbabilityLowerLimit, mutationIterations, mutationDegreesOfFreedom,
                                                            mutationDecay, hillClimbingIterations, seed, trace, numberOfMaxThreads,
                                                            levelsOfStutterRecursion, dualEstimation, traceLimit) {
    if (length(seed) == 1) {
        seed = sample(1:1e6, numberOfPopulations)
    }

    if (length(seed) != numberOfPopulations) {
        warning("All initial populations are equal.")
    }

    numberOfThreads <- min(numberOfPopulations, numberOfMaxThreads)
    mutationDecaySplit <- split(mutationDecay, rep(1:numberOfIterations, each = numberOfInnerIterations))
    currentPopulationList <- mclapply(1:numberOfPopulations, function(i) {
        .initialisingParallelEvolutionaryAlgorithm(
            numberOfMarkers = numberOfMarkers,
            numberOfAlleles = numberOfAlleles,
            numberOfContributors = numberOfContributors,
            numberOfKnownContributors = numberOfKnownContributors,
            knownProfiles = knownProfiles,
            allKnownProfiles = allKnownProfiles,
            coverage = coverage,
            potentialParents = potentialParentsList,
            markerImbalances = markerImbalances,
            convexMarkerImbalanceInterpolation = convexMarkerImbalanceInterpolation,
            noiseParameters = noiseParameters,
            tolerance = tolerance,
            theta = theta,
            alleleFrequencies = alleleFrequencies,
            populationSize = populationSize,
            seed = seed[i],
            levelsOfStutterRecursion = levelsOfStutterRecursion,
            dualEstimation = dualEstimation)
    }, mc.cores = numberOfThreads)

    j = 0
    k = 0
    converged = FALSE
    topFittestIndividuals <- list()
    while (!converged) {
        ## Migration between sub-populations
        currentPopulationMigratedList <- .parallelPopulationEvolutionaryAlgorithmMigration(currentPopulationList)

        ## Updating sub-populations
        newPopulation <- mclapply(1:numberOfPopulations, function(i) {
            PEA_i <- .runningParallelEvolutionaryAlgorithm(
                numberOfMarkers = numberOfMarkers,
                numberOfAlleles = numberOfAlleles,
                numberOfContributors = numberOfContributors,
                numberOfKnownContributors = numberOfKnownContributors,
                knownProfiles = knownProfiles,
                allKnownProfiles = allKnownProfiles,
                coverage = coverage,
                potentialParents = potentialParentsList,
                markerImbalances = markerImbalances,
                convexMarkerImbalanceInterpolation = convexMarkerImbalanceInterpolation,
                noiseParameters = noiseParameters,
                tolerance = tolerance,
                theta = theta,
                alleleFrequencies = alleleFrequencies,
                numberOfIterations = numberOfInnerIterations,
                numberOfIterationsEqualMinMax = numberOfInnerIterations,
                numberOfFittestIndividuals = populationSize,
                parentSelectionWindowSize = parentSelectionWindowSize,
                allowParentSurvival = allowParentSurvival,
                fractionEnsuredSurvival = fractionFittestIndividuals,
                crossoverProbability = crossoverProbability,
                mutationProbabilityLowerLimit = mutationProbabilityLowerLimit,
                mutationIterations = mutationIterations,
                mutationDegreesOfFreedom = mutationDegreesOfFreedom,
                mutationDecay = mutationDecaySplit[[j + 1]],
                hillClimbingIterations = hillClimbingIterations,
                seed = sample(1:1e6, 1),
                trace = FALSE,
                encodedPopulationList = currentPopulationMigratedList[[i]]$EncodedProfiles,
                sampleParametersList = currentPopulationMigratedList[[i]]$SampleParameters,
                noiseParametersList = currentPopulationMigratedList[[i]]$NoiseParameters,
                mixtureParametersList = currentPopulationMigratedList[[i]]$MixtureParameters,
                markerParametersList = currentPopulationMigratedList[[i]]$MarkerImbalanceParameters,
                fitnessList = currentPopulationMigratedList[[i]]$Fitness,
                levelsOfStutterRecursion = levelsOfStutterRecursion,
                dualEstimation = dualEstimation
            )
            return(PEA_i)
        }, mc.cores = numberOfThreads)

        if (any(sapply(newPopulation, function(np) is.null(np))))
            stop("Contact your local support staff: a thread returned a 'NULL' object indicating that the internal 'c++' code crashed.")

        currentPopulationList <- lapply(newPopulation, function(np) np$LastPopulation)
        currentFittestList <- unlist(lapply(newPopulation, function(np) np$FittestIndividuals), recursive = FALSE)

        ## Updating list of fittest individuals
        currentFittestList <- currentFittestList[!sapply(currentFittestList, is.null)]
        fittestIndividuals <- append(topFittestIndividuals, currentFittestList)

        duplicatedFittestIndividuals <- duplicated(lapply(fittestIndividuals, function(cfl) cfl$EncodedUnknownProfiles))
        fittestIndividuals <- fittestIndividuals[!duplicatedFittestIndividuals]

        numberOfKeptIndividuals <- max(1, min(numberOfFittestIndividuals, length(fittestIndividuals)))
        sortedFittestIndividuals <- order(unlist(lapply(fittestIndividuals, function(fi) fi$Fitness)), decreasing = TRUE)

        topFittestIndividuals <- fittestIndividuals[sortedFittestIndividuals[1:numberOfKeptIndividuals]]

        ## Updating convergence condition
        populationFitness <- do.call("c", lapply(currentPopulationList, function(cpl) cpl$Fitness))
        maxPopulationFitness <- do.call("c", lapply(currentPopulationList, function(cpl) max(cpl$Fitness)))

        if ((max(maxPopulationFitness) - min(maxPopulationFitness)) < tolerance[1]) {
            k = k + 1
        } else {
            k = 0
        }

        j = j + 1
        converged = (k == numberOfIterationsEqualMax) | (j == numberOfIterations)

        if (trace) {
            if ((j == 1) | converged | ((j %% traceLimit) == 0))
                cat("\tOuter iteration:", j, "\n",
                    "\t\tFitness:\n",
                    "\t\t  Highest:", max(populationFitness), "\n",
                    "\t\t  Average:", mean(populationFitness), "\n",
                    "\t\t  Lowest:", min(populationFitness), "\n",
                    "\t\tSubpopulations maxima difference:", (max(maxPopulationFitness) - min(maxPopulationFitness)), "\n", #fractionOfUniqueSubpopulationMaxima, "\n",
                    "\t\tTermination counter:", k, "/", numberOfIterationsEqualMax, "\n")
        }
    }

    resList <- list(U = topFittestIndividuals, UMaxSize = numberOfFittestIndividuals, NumberOfIterations = j,
                    NumberOfIterationsEqualMax = numberOfIterationsEqualMax, NumberOfPopulations = numberOfPopulations)
    return(resList)
}
