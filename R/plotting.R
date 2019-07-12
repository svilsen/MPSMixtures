#' Q-Q plots for profile list
#'
#' @description Q-Q plots for the allele coverage and noise Poisson-gamma coverage models.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param profileList A list containing the expected contribution matrix and paramters of the fitted profile.
#' @param component Which of the two model components should be constructed? 'Allele', 'Noise', or 'both'?
#'
#' @return ggplot2-object.
ggplotQQPlotProfiles <- function(sampleTibble, profileList, component = "allele", residualEnvelopes = NULL){
    Type = NULL

    SP <- profileList$Parameters$SampleParameters
    NP <- profileList$Parameters$NoiseParameters
    MP <- profileList$Parameters$MixtureParameters
    MIP <- profileList$Parameters$MarkerImbalanceParameters

    MIPTibble <- sampleTibble %>% distinct_(~Marker) %>% mutate_(MIP = ~MIP)

    C <- sampleTibble$Coverage
    MI <- (sampleTibble %>% left_join(MIPTibble, by = "Marker"))$MIP

    ECM <- profileList$ExpectedContributionMatrix
    NV <- profileList$NoiseVector

    MU <- c(SP[1] * MI * (ECM %*% MP) + NP[1] * NV)
    D <- SP[2] * (1.0 - NV) + NP[2] * NV

    DRTibble <- tibble(C = C, NV = NV, MU = MU, D = D) %>%
        mutate("EC" = apply(ECM, 1, sum)) %>%
        rowwise() %>%
        mutate("DR" = MPSMixtures:::.devianceResidualPoissonGammaDistribution(C, MU, ifelse(NV != 1, MU / D, D)),
               "Type" = if (NV == 1) "Noise coverage" else if (EC > 0) "Allele coverage" else NA) %>%
        filter(!is.na(Type)) %>% ungroup() %>% arrange(Type)

    noise_ones_inflation = NULL
    if (NP[3] > 2e-8) {
        noise_ones <- which(DRTibble$C == 1 & DRTibble$NV == 1)
        noise_ones_inflation <- sample(noise_ones, floor(length(noise_ones) * NP[3]))

        if (length(noise_ones_inflation) > 0) {
            DRTibble <- DRTibble[-noise_ones_inflation, ]
        }
    }

    qq_allele <- qqnorm((DRTibble %>% filter(Type == "Allele coverage"))$DR, plot.it = FALSE)
    qq_noise <- qqnorm((DRTibble %>% filter(Type == "Noise coverage"))$DR, plot.it = FALSE)

    DRTibble <- DRTibble %>% mutate(X = c(qq_allele$x, qq_noise$x), Y = c(qq_allele$y, qq_noise$y))
    if (tolower(component) == "allele") {
        p <- ggplot(DRTibbleAllele, aes_(x = ~X, y = ~Y)) +
            geom_point() + facet_wrap(~Type) + theme_bw() +
            xlab("Theoretical quantiles") + ylab("Sample quantiles") + ggtitle("Q-Q plot") +
            geom_smooth(method = "lm", se = FALSE)
    }
    else if (tolower(component) == "noise") {
        p <- ggplot(DRTibble %>% filter(Type == "Noise coverage"), aes_(x = ~X, y = ~Y)) +
            geom_point() + facet_wrap(~Type) + theme_bw() +
            xlab("Theoretical quantiles") + ylab("Sample quantiles") + ggtitle("Q-Q plot") +
            geom_smooth(method = "lm", se = FALSE)
    }
    else if (tolower(component) == "both") {
        p <- ggplot(DRTibble, aes_(x = ~X, y = ~Y)) +
            geom_point() + facet_wrap(~Type) + theme_bw() +
            xlab("Theoretical quantiles") + ylab("Sample quantiles") + ggtitle("Q-Q plot") +
            geom_smooth(method = "lm", se = FALSE)
    }
    else {
        p <- ggplot()
    }

    return(p)
}

#' Q-Q plot and envelope of a DNA profile
#'
#' @description Q-Q plot and envelopes for the allele coverage component of a DNA profile.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param profileList A list containing the expected contribution matrix and paramters of the fitted profile.
#' @param markerImbalances A vector of prior marker imbalances.
#' @param knownProfilesList A list of tibbles containing the alleles of the known contributors.
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample. If NULL then a stutterRatioModel should be provided.
#' @param stutterRatioModel A linear model of class \link{lm} modelling the relationship between coverage and stutter. Only needed if the potential parents list is not provided.
#' @param levelsOfStutterRecursion The number of layers used in the stutter recursion.
#' @param convexMarkerImbalanceInterpolation A fraction used to create a convex combination of the of MoM and the prior estimates of the marker imbalances.
#' @param tolerance The tolerance used for termination in parameter estimation.
#' @param numberOfThreads The number of threads passed to the \link{potentialParentsMultiCore}-function if applicable.
#' @param numberOfSimulations The number of simulations used to create the envelopes.
#' @param envelopes The envelope quantiles.
#' @param trace Show trace (TRUE/FALSE)?
#' @param traceLimit Limits the trace of the parallel implementation.
#'
#' @return ggplot2-object.
ggplotQQPlotProfilesEnvelopes <- function(sampleTibble, profileList, markerImbalances, knownProfilesList, potentialParentsList,
                                          stutterRatioModel = NULL, levelsOfStutterRecursion = 2, convexMarkerImbalanceInterpolation = 0.8,
                                          tolerance = rep(1e-8, 4), numberOfThreads = 4, numberOfSimulations = 1000, envelopes = c(0.005, 0.995),
                                          trace = FALSE, traceLimit = 100) {
    Type = NULL

    SP <- profileList$Parameters$SampleParameters
    NP <- profileList$Parameters$NoiseParameters
    MP <- profileList$Parameters$MixtureParameters
    MIP <- profileList$Parameters$MarkerImbalanceParameters

    MIPTibble <- sampleTibble %>% distinct_(~Marker) %>% mutate_(MIP = ~MIP)

    C <- sampleTibble$Coverage
    MI <- (sampleTibble %>% left_join(MIPTibble, by = "Marker"))$MIP

    ECM <- profileList$ExpectedContributionMatrix
    NV <- profileList$NoiseVector

    MU <- c(SP[1] * MI * (ECM %*% MP) + NP[1] * NV)
    D <- SP[2] * (1.0 - NV) + NP[2] * NV

    DRTibble <- tibble(C = C, NV = NV, MU = MU, D = D) %>%
        mutate("EC" = apply(ECM, 1, sum)) %>%
        rowwise() %>%
        mutate("DR" = MPSMixtures:::.devianceResidualPoissonGammaDistribution(C, MU, ifelse(NV != 1, MU / D, D)),
               "Type" = if (NV == 1) "Noise coverage" else if (EC > 0) "Allele coverage" else NA) %>%
        filter(!is.na(Type)) %>% ungroup() %>% arrange(Type)

    noise_ones_inflation = NULL
    if (NP[3] > 2e-8) {
        noise_ones <- which(DRTibble$C == 1 & DRTibble$NV == 1)
        noise_ones_inflation <- sample(noise_ones, floor(length(noise_ones) * NP[3]))

        if (length(noise_ones_inflation) > 0) {
            DRTibble <- DRTibble[-noise_ones_inflation, ]
        }
    }

    qq_allele <- qqnorm((DRTibble %>% filter(Type == "Allele coverage"))$DR, plot.it = FALSE)
    qq_noise <- qqnorm((DRTibble %>% filter(Type == "Noise coverage"))$DR, plot.it = FALSE)

    DRTibble <- DRTibble %>% mutate(X = c(qq_allele$x, qq_noise$x), Y = c(qq_allele$y, qq_noise$y))

    allele_coverage_indices <- (sampleTibble %>% mutate(NV = NV, Index = 1:n()) %>% filter(NV == 0))$Index
    DRTibbleAllele <- DRTibble %>% filter(Type == "Allele coverage")

    sampleTibble_i = sampleTibble
    resMatrix = matrix(0, nrow = dim(DRTibbleAllele)[1], ncol = numberOfSimulations)
    for (i in 1:numberOfSimulations) {
        if (trace & ((i == 1) | (i == numberOfSimulations) | ((i %% traceLimit) == 0)))
            cat(i, " ")

        sample_coverage_i = rnbinom(dim(DRTibbleAllele)[1], mu = DRTibbleAllele$MU, size = DRTibbleAllele$MU / DRTibbleAllele$D)
        sampleTibble_i$Coverage[allele_coverage_indices] <- sample_coverage_i

        known_i <- MPSMixtures:::estimateParametersOfKnownProfiles(sampleTibble_i, markerImbalances, knownProfilesList,
                                                                   potentialParentsList, stutterRatioModel,
                                                                   levelsOfStutterRecursion, convexMarkerImbalanceInterpolation,
                                                                   tolerance, numberOfThreads)

        SP_i <- known_i$Parameters$SampleParameters
        MP_i <- known_i$Parameters$MixtureParameters
        MI_i <- rep(known_i$Parameters$MarkerImbalanceParameters, (sampleTibble[allele_coverage_indices, ] %>% group_by(Marker) %>% summarise(Count = n()))$Count)

        ECM_i <- known_i$ExpectedContributionMatrix

        MU_i <- c(SP_i[1] * MI_i * (ECM_i %*% MP_i)[allele_coverage_indices, ])
        R_i <- sapply(seq_along(MU_i), function(j) MPSMixtures:::.devianceResidualPoissonGammaDistribution(sample_coverage_i[j], MU_i[j], MU_i[j] / SP_i[2]))

        resMatrix[, i] <- sort(R_i)
    }

    CI = apply(resMatrix, 1, quantile, prob = envelopes)
    reverse_order <- order((tibble("X" = qq_allele$x) %>% mutate("Count" = 1:n()) %>% arrange_(~X))$Count)

    DRTibbleAllele <- DRTibbleAllele %>% mutate("LowerEnvelope" = CI[1, reverse_order], "UpperEnvelope" = CI[2, reverse_order])
    p <- ggplot(DRTibbleAllele, aes_(x = ~X, y = ~DR)) +
        facet_wrap(~Type) +
        geom_ribbon(aes_(ymin = ~LowerEnvelope, ymax = ~UpperEnvelope),
                    alpha = 0.3, colour = "dodgerblue2", fill = "dodgerblue2") +
        geom_point() +
        theme_bw() + xlab("Theoretical quantiles") + ylab("Sample quantiles") + ggtitle("Q-Q plot") +
        geom_abline(intercept = 0, slope = 1, colour = "black", alpha = 0.9, size = 0.8)

    return(p)
}

#' P-P plot and envelope of a DNA profile
#'
#' @description P-P plot and envelopes for the allele coverage component of a DNA profile.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param profileList A list containing the expected contribution matrix and paramters of the fitted profile.
#' @param markerImbalances A vector of prior marker imbalances.
#' @param knownProfilesList A list of tibbles containing the alleles of the known contributors.
#' @param potentialParentsList A list containing a list of potential parents for each allele in the sample. If NULL then a stutterRatioModel should be provided.
#' @param stutterRatioModel A linear model of class \link{lm} modelling the relationship between coverage and stutter. Only needed if the potential parents list is not provided.
#' @param levelsOfStutterRecursion The number of layers used in the stutter recursion.
#' @param convexMarkerImbalanceInterpolation A fraction used to create a convex combination of the of MoM and the prior estimates of the marker imbalances.
#' @param tolerance The tolerance used for termination in parameter estimation.
#' @param numberOfThreads The number of threads passed to the \link{potentialParentsMultiCore}-function if applicable.
#' @param numberOfSimulations The number of simulations used to create the envelopes.
#' @param envelopes The envelope quantiles.
#' @param trace Show trace (TRUE/FALSE)?
#' @param traceLimit Limits the trace of the parallel implementation.
#'
#' @return ggplot2-object.
ggplotPPPlotProfilesEnvelopes <- function(sampleTibble, profileList, markerImbalances, knownProfilesList, potentialParentsList,
                                          stutterRatioModel = NULL, levelsOfStutterRecursion = 2, convexMarkerImbalanceInterpolation = 0.8,
                                          tolerance = rep(1e-8, 4), numberOfThreads = 4, numberOfSimulations = 1000, envelopes = c(0.005, 0.995),
                                          trace = FALSE, traceLimit = 100) {
    Type = NULL

    SP <- profileList$Parameters$SampleParameters
    NP <- profileList$Parameters$NoiseParameters
    MP <- profileList$Parameters$MixtureParameters
    MIP <- profileList$Parameters$MarkerImbalanceParameters

    MIPTibble <- sampleTibble %>% distinct_(~Marker) %>% mutate_(MIP = ~MIP)

    C <- sampleTibble$Coverage
    MI <- (sampleTibble %>% left_join(MIPTibble, by = "Marker"))$MIP

    ECM <- profileList$ExpectedContributionMatrix
    NV <- profileList$NoiseVector

    MU <- c(SP[1] * MI * (ECM %*% MP) + NP[1] * NV)
    D <- SP[2] * (1.0 - NV) + NP[2] * NV

    DRTibble <- tibble(C = C, NV = NV, MU = MU, D = D) %>%
        mutate("EC" = apply(ECM, 1, sum)) %>%
        rowwise() %>%
        mutate("DR" = MPSMixtures:::.devianceResidualPoissonGammaDistribution(C, MU, ifelse(NV != 1, MU / D, D)),
               "Type" = if (NV == 1) "Noise coverage" else if (EC > 0) "Allele coverage" else NA) %>%
        filter(!is.na(Type)) %>% ungroup() %>% arrange(Type)

    noise_ones_inflation = NULL
    if (NP[3] > 2e-8) {
        noise_ones <- which(DRTibble$C == 1 & DRTibble$NV == 1)
        noise_ones_inflation <- sample(noise_ones, floor(length(noise_ones) * NP[3]))

        if (length(noise_ones_inflation) > 0) {
            DRTibble <- DRTibble[-noise_ones_inflation, ]
        }
    }

    allele_coverage_indices <- (sampleTibble %>% mutate(NV = NV, Index = 1:n()) %>% filter(NV == 0))$Index
    DRTibbleAllele <- DRTibble %>%
        filter(Type == "Allele coverage")

    sampleTibble_i = sampleTibble
    resMatrix = matrix(0, nrow = dim(DRTibbleAllele)[1], ncol = numberOfSimulations)
    for (i in 1:numberOfSimulations) {
        if (trace & ((i == 1) | (i == numberOfSimulations) | ((i %% traceLimit) == 0)))
            cat(i, " ")

        sample_coverage_i = rnbinom(dim(DRTibbleAllele)[1], mu = DRTibbleAllele$MU, size = DRTibbleAllele$MU / DRTibbleAllele$D)
        sampleTibble_i$Coverage[allele_coverage_indices] <- sample_coverage_i

        known_i <- MPSMixtures:::estimateParametersOfKnownProfiles(sampleTibble_i, markerImbalances, knownProfilesList,
                                                                   potentialParentsList, stutterRatioModel,
                                                                   levelsOfStutterRecursion, convexMarkerImbalanceInterpolation,
                                                                   tolerance, numberOfThreads)

        SP_i <- known_i$Parameters$SampleParameters
        MP_i <- known_i$Parameters$MixtureParameters
        MI_i <- rep(known_i$Parameters$MarkerImbalanceParameters, (sampleTibble[allele_coverage_indices, ] %>% group_by(Marker) %>% summarise(Count = n()))$Count)

        ECM_i <- known_i$ExpectedContributionMatrix

        MU_i <- c(SP_i[1] * MI_i * (ECM_i %*% MP_i)[allele_coverage_indices, ])
        P_i <- pnbinom(sample_coverage_i, mu = MU_i, size = MU_i / SP_i[2])

        resMatrix[, i] <- sort(P_i)
    }

    CI = apply(resMatrix, 1, quantile, prob = envelopes)

    DRTibbleAllele <- DRTibbleAllele %>%
        mutate(Y = pnbinom(C, mu = MU, size = MU / D)) %>%
        ungroup() %>%
        arrange_(~Y) %>%
        mutate(X = (1 : n()) / n() - 0.5 / n(),
               "LowerEnvelope" = CI[1, ],
               "UpperEnvelope" = CI[2, ])

    p <- ggplot(DRTibbleAllele, aes_(x = ~X, y = ~Y)) +
        facet_wrap(~Type) +
        geom_ribbon(aes_(ymin = ~LowerEnvelope, ymax = ~UpperEnvelope),
                    alpha = 0.3, colour = "dodgerblue2", fill = "dodgerblue2") +
        geom_point() + theme_bw() +
        xlab("Uniform CDF") + ylab("Sample CDF") + ggtitle("P-P plot") +
        geom_abline(intercept = 0, slope = 1, colour = "black", alpha = 0.9, size = 0.8)

    return(p)
}

#' Prediction intervals.
#'
#' @description Prediction intervals of the Poisson-gamma model.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param profileList A list containing the expected contribution matrix and paramters of the fitted profile.
#' @param predictionInterval The upper and lower quantiles of the prediction interval.
#' @param numberOfSimulations The number of simulations used to create the prediction interval.
#' @param contributorNames The names of the contributors; must be 'NULL', or a vector of length '1' or the number of contributors.
#'
#' @return ggplot2-object.
ggplotPredictionIntervals <- function(sampleTibble, profileList, predictionInterval = c(0.005, 0.995),
                                      numberOfSimulations = 10000, contributorNames = NULL) {
    SP <- profileList$Parameters$SampleParameters
    NP <- profileList$Parameters$NoiseParameters
    MP <- profileList$Parameters$MixtureParameters

    MIP <- profileList$Parameters$MarkerImbalanceParameters
    MIPTibble <- sampleTibble %>% distinct_(~Marker) %>% mutate_(MIP = ~MIP)
    MI <- (sampleTibble %>% left_join(MIPTibble, by = "Marker"))$MIP

    numberOfProfiles <- length(MP)
    if (is.null(contributorNames)) {
        contributorNames <- paste("C", 1:numberOfProfiles, sep = "")
    } else if (length(contributorNames) == 1) {
        contributorNames <- paste(contributorNames, 1:numberOfProfiles, sep = "")
    } else if (length(contributorNames) > numberOfProfiles) {
        stop("The length of 'contributorNames' should be equal to '1' or the number of contributors.")
    }

    ECM <- profileList$ExpectedContributionMatrix
    NV <- profileList$NoiseVector

    MU <- c(SP[1] * MI * (ECM %*% MP) + NP[1] * NV)
    D <- SP[2] * (1.0 - NV) + NP[2] * NV

    C <- sampleTibble$Coverage

    simulations <- matrix(0, ncol = length(C), nrow = numberOfSimulations)
    for (N in 1:numberOfSimulations) {
        simulations[N, ] <- rnbinom(length(C), mu = MU, size = ifelse(NV == 1.0, D, MU / D))
    }

    predictedInterval <- apply(simulations, 2, function(si) quantile(si, probs = predictionInterval))

    plotTibbleCollected <- sampleTibble %>%
        mutate(ExpectedCoverage = MU, Dispersion = D,
               LowerPI = predictedInterval[1, ],
               UpperPI = predictedInterval[2, ])

    profileECMReweighted <- t(apply((t(profileList$ExpectedContributionMatrix) * profileList$Parameters$MixtureParameters), 2, function(ii) {
        res <- if (sum(ii) != 0) ii / sum(ii) else rep(0, length(ii))
        return(res)
    }))

    plotTibbleList <- vector("list", numberOfProfiles)
    for (j in 1:numberOfProfiles) {
        plotTibbleList[[j]] <- plotTibbleCollected %>%
            mutate_(ProfileCoverage = ~Coverage * profileECMReweighted[, j],
                    Contributor = ~contributorNames[j]) %>%
            filter_(~ProfileCoverage != 0)
    }

    plotTibble <- bind_rows(plotTibbleList) %>%
        group_by_(~Marker, ~Allele, ~Region) %>%
        summarise_(ProfileCoverage = ~sum(ProfileCoverage),
                   ExpectedCoverage = ~unique(ExpectedCoverage),
                   LowerPI = ~unique(LowerPI),
                   UpperPI = ~unique(UpperPI),
                   Contributor = ~paste(Contributor, collapse = "/")) %>%
        ungroup() %>%
        mutate_(Contributor = ~paste("Contributor:", Contributor))

    p <- ggplot(plotTibble, aes_(x = ~Allele)) +
        geom_bar(aes_(y = ~ProfileCoverage, fill = ~Contributor), stat = "identity") +
        geom_point(aes_(y = ~ExpectedCoverage)) +
        geom_errorbar(aes_(ymin = ~LowerPI, ymax = ~UpperPI), width = 0.3) +
        facet_grid(Contributor~Marker, scales = "free_x", space = "free_x") +
        scale_x_continuous(breaks = seq(min(plotTibble$Allele), max(plotTibble$Allele)),
                           labels = seq(min(plotTibble$Allele), max(plotTibble$Allele))) +
        theme_bw() + xlab("Allele length") + ylab("Coverage") +
        theme(legend.position = "none")

    return(p)
}

