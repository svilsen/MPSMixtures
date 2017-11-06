if (FALSE) {
    ## EXAMPLE SET UP
    #
    library("STRMPS")
    library("tidyverse")
    library("Rsolnp")

    #
    loadRData <- function(fileName){
        load(fileName)
        get(ls()[ls() != "fileName"])
    }

    load("~/AAU/PhD/Articles/DropOut/R/Files/noiseModelTibblesExactAUTOSOMAL.RData")
    noiseDropOutsMixturesExactAUTOSOMAL <- lapply(seq_along(noiseModelTibbleExactAUTOSOMAL), function(nn) {
        tt <- noiseModelTibbleExactAUTOSOMAL[[nn]] %>% mutate(Threshold = 0)
        return(tt)
    })

    noiseDropOutsMixturesExactThresholdAUTOSOMAL <- lapply(seq_along(noiseModelTibbleCollapsedExactAUTOSOMAL), function(nn) {
        tt <- noiseModelTibbleCollapsedExactAUTOSOMAL[[nn]] %>% mutate(DropOut = DropOutThreshold) %>% select(-DropOutThreshold) %>%
            filter((DropOut - ThresholdedObservation) == 0)
        return(tt)
    })

    ## True profiles mixture data
    trueProfiles <- loadRData("~/AAU/PhD/Articles/DropOut/R/Files/trueProfilesMixturesRegion.RData")
    mixtureDescription <- loadRData("~/AAU/PhD/Articles/DropOut/R/Files/mixtureDescription.RData") %>%
        select(Sample, SampleInternal, Replicate, MixType, Major, Minor, MajorRelationship) %>%
        mutate(MajorRelationshipNumeric = sapply(strsplit(MajorRelationship, ":"), function(mr) 1 - min(as.numeric(mr)) / (max(as.numeric(mr)) + 1) ),
               MinorRelationshipNumeric = 1 - MajorRelationshipNumeric)

    ## BLMM stutter model
    stutterRatioModel <- loadRData("~/AAU/PhD/Articles/DropOut/R/Files/stutterRatioBLMMModel.RData")
    stutterRatioModelTibble <- tibble(Marker = sapply(strsplit(names(coef(stutterRatioModel)), "Marker"), function(x) x[2]),
                                      Intercept = 0, Slope = unname(coef(stutterRatioModel))) %>% arrange(Marker)

    load("~/AAU/PhD/Articles/DropOut/R/Files/estimatedParametersSimpleModelReferenceDataBase.RData")
    markerParametersTibble <- MarkerParametersTibbleFixedHomozygoteScalar

    unique(noiseDropOutsMixturesExactAUTOSOMAL[[20]]$MajorRelationship)
    coverageTibble <-  noiseDropOutsMixturesExactAUTOSOMAL[[20]] %>%
        mutate(MarkerNumeric = as.numeric(relevel(factor(Marker), markerParametersTibble$Marker[1]))) %>%
        dplyr::select(Marker, MarkerNumeric, Allele, MotifLength, Region, Coverage) %>%
        left_join(markerParametersTibble %>% select(-MarkerNumeric), by = "Marker")

    pp <- potentialParents(coverageTibble, stutterRatioModel, trace = TRUE)
    rrpp <- lapply(pp, function(ppp) lapply(ppp, function(pppp) {
        pppp %>% select(PotentialParent, StutterRatio) %>% as.matrix()
    }))

    numberOfAlleles <- (coverageTibble %>% group_by(Marker) %>% summarise(Count = n()))$Count
    alleleFrequencies = rep(1 / numberOfAlleles, times = numberOfAlleles)

    ## Second profile is major when n > 16 (15 and 16 are 1:1 mixture)
    Hp <- setHypothesis(coverageTibble, c(2), trueProfiles, 0, alleleFrequencies)
    Hd <- setHypothesis(coverageTibble, c(2), trueProfiles[2], 0, alleleFrequencies)

    allknowntest <- MPSMixtures:::.setupIndividual(length(numberOfAlleles), numberOfAlleles, 2, Hp[[1]]$NumberOfKnownProfiles, Hp[[1]]$KnownProfiles,
                                                   coverageTibble$Coverage, rrpp, coverageTibble$MarkerImbalance, 0.001, 0, alleleFrequencies)
    allknowntest$SampleParameters
    allknowntest$NoiseParameters
    allknowntest$MixtureParameters


    EC <- coverageTibble$MarkerImbalance * allknowntest$SampleParameters[1] * allknowntest$ExpectedContributionMatrix %*% allknowntest$MixtureParameters
    coverageTibble$Coverage - EC

    set.seed(123)
    LRSinglePopulation <- LR(sampleTibble = coverageTibble, Hp = Hp, Hd = Hd, potentialParentsList = rrpp,
                             control = LR.control(numberOfPopulations = 1, numberOfIterations = 200, populationSize = 1000,
                                                  numberOfFittestIndividuals = 100, hillClimbingIterations = 2,
                                                  parentSelectionWindowSize = 15, simplifiedReturn = F, allowParentSurvival = T))

    LRParallelPopulations <- LR(sampleTibble = coverageTibble, Hp = Hp, Hd = Hd, potentialParentsList = rrpp,
                                control = LR.control(numberOfPopulations = 8, numberOfIterations = 100, populationSize = 10,
                                                     numberOfFittestIndividuals = 100, hillClimbingIterations = 2,
                                                     parentSelectionWindowSize = 2, simplifiedReturn = F, allowParentSurvival = T))

    # load("LRTest.RData")
    # save(LRSinglePopulation, LRParallelPopulations, file = "LRTest.RData")

    LRSinglePopulation$AllPairwiseComparisonData[[1]]$HpParameterEstimates
    LRParallelPopulations$AllPairwiseComparisonData[[1]]$HpParameterEstimates

    LRSinglePopulation$AllPairwiseComparisonData[[1]]$HdParameterEstimates
    LRParallelPopulations$AllPairwiseComparisonData[[1]]$HdParameterEstimates

    LRSinglePopulation$ComparisonTable
    LRParallelPopulations$ComparisonTable

    # Adj...
    adjacencyMatrix <- function(n) {
        library("igraph")
        get.adjacency(graph.edgelist(matrix(c(rep(seq(1, n), times = 2), (1:n) %% n + 1, (1:n - 3) %% n + 1), ncol = 2), directed = TRUE), sparse = FALSE)
    }

    numberOfIterationsUntilHighFitnessIndividualPropagatesThroughPopulation <- function(M) {
        library("expm")
        N = nrow(M)
        I = diag(1, N)
        A = M + I

        k = 1
        while ((k <= N) && (!all(A != 0))) {
            k = k + 1
            A = (M + I) %^% k
        }

        return(k)
    }

    res <- matrix(c(seq(1, 14), ceiling((seq(1, 14) + 1) / 3), floor((seq(1, 14)) / 3) + 1, rep(0, 14)), ncol = 4)
    for (i in 1:14) {
        M = adjacencyMatrix(i)
        res[i, 4] <- numberOfIterationsUntilHighFitnessIndividualPropagatesThroughPopulation(M)
    }

}
