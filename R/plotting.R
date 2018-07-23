#' Compare profile combinations
#'
#' @details Function used for visually comparing two profiles.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param profileList1 A list containing the expected contribution matrix and the mixture paramters of profile 1.
#' @param profileList2 A list containing the expected contribution matrix and the mixture paramters of profile 2.
#' @param plotDifferencesOnly Restrict the returned \link{ggplot}-object to markers where the two profiles mismatch.
#' @param profileNames Vector of profile names; should have length '2' or be 'NULL'.
#' @param contributorNames Vector of contributor names; should have the same length as the number of contributors or be 'NULL'.
#'
#' @return ggplot2-object.
ggplotComparingProfiles <- function(sampleTibble, profileList1, profileList2, plotDifferencesOnly = FALSE,
                                    profileNames = NULL, contributorNames = NULL) {
    if (is.null(profileNames)) {
        profileNames <- c(1, 2)
    }
    else if (length(profileNames) == 1) {
        profileNames <- paste(profileNames, 1:2, sep = "_")
    }
    else if (length(profileNames) > 2) {
        stop("'profileNames' must have length 1 or 2, or set to 'NULL'.")
    }

    numberOfProfiles1 <- dim(profileList1$ExpectedContributionMatrix)[2]
    numberOfProfiles2 <- dim(profileList2$ExpectedContributionMatrix)[2]

    matching <- match(data.frame(profileList1$ExpectedContributionMatrix), data.frame(profileList2$ExpectedContributionMatrix))
    matchingProfiles <- cbind(which(matching != 0), matching[which(matching != 0)])

    restructuringList1 <- profileList1[c("ExpectedContributionMatrix", "Parameters")]
    restructuringList2 <- profileList2[c("ExpectedContributionMatrix", "Parameters")]

    if (dim(matchingProfiles)[1] != 0) {
        for (i in 1:dim(matchingProfiles)[1]) {
            switching <- 1:numberOfProfiles2
            switching[matchingProfiles[i, 1]] <- matchingProfiles[i, 2]
            switching[matchingProfiles[i, 2]] <- matchingProfiles[i, 1]

            restructuringList2$ExpectedContributionMatrix <- restructuringList2$ExpectedContributionMatrix[, switching]
            restructuringList2$MixtureParameters <- restructuringList2$Parameters$MixtureParameters[switching]
        }
    }

    profile1ECMReweighted <- t(apply((t(restructuringList1$ExpectedContributionMatrix) * restructuringList1$Parameters$MixtureParameters), 2, function(ii) {
        res <- if (sum(ii) != 0) ii / sum(ii) else rep(0, length(ii))
        return(res)
    }))

    profile2ECMReweighted <- t(apply((t(restructuringList2$ExpectedContributionMatrix) * restructuringList2$Parameters$MixtureParameters), 2, function(ii) {
        res <- if (sum(ii) != 0) ii / sum(ii) else rep(0, length(ii))
        return(res)
    }))

    profiles <- list(profile1ECMReweighted, profile2ECMReweighted)
    numberOfProfiles <- c(numberOfProfiles1, numberOfProfiles2)
    partialSumProfiles <- .partialSumEigen(numberOfProfiles)

    plotTibble <- vector("list", sum(numberOfProfiles) + 2)
    for (i in 1:2) {
        for (j in 1:numberOfProfiles[i]) {
            k = partialSumProfiles[i] + j
            plotTibble[[k]] <- sampleTibble %>%
                mutate_(ProfileCoverage = ~Coverage * profiles[[i]][, j],
                       Profile = ~paste("Profile:", profileNames[i]),
                       Contributor = ~paste0("C", j)) %>%
                filter_(~ProfileCoverage != 0)
        }

        plotTibble[[sum(numberOfProfiles) + i]] <- sampleTibble %>%
            mutate_(ProfileCoverage = ~Coverage * as.numeric(rowSums(profiles[[i]]) == 0),
                    Profile = ~paste("Profile:", profileNames[i]),
                    Contributor = ~paste0("N")) %>%
            filter_(~ProfileCoverage != 0)
    }

    plotTibble <- bind_rows(plotTibble)
    gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
    }

    colours <- gg_color_hue(2)
    if (plotDifferencesOnly) {
        P1 <- P2 <- NULL
        profileFrame <- tibble(P1 = apply(profileList1$ExpectedContributionMatrix, 1, sum),
                               P2 = apply(profileList2$ExpectedContributionMatrix, 1, sum),
                               Difference = P1 - P2)
        splitMarker <- split(profileFrame, sampleTibble$Marker)
        differenceMarker <- sapply(splitMarker, function(xx) sum(abs(xx[, 3])) != 0)

        plotTibble <- plotTibble %>% filter_(~Marker %in% names(differenceMarker)[differenceMarker])
    }

    if (dim(plotTibble)[1] > 0) {
        p.Profile <- ggplot(plotTibble, aes_(x = ~Allele, y = ~ProfileCoverage, fill = ~Contributor)) + geom_bar(stat = "identity") +
            facet_grid(Marker ~ Profile, scales = "free_y") + scale_fill_manual(values = c(colours[1], colours[2], "#000000")) +
            theme_bw() + theme(legend.position = "top") + ylab("Coverage") + xlab("Allele length") +
            scale_x_continuous(breaks = seq(min(plotTibble$Allele), max(plotTibble$Allele)))
    }
    else {
        p.Profile <- ggplot()
    }

    return(p.Profile)
}

#' Q-Q plots for profile list
#'
#' @description Q-Q plots for the allele coverage and noise Poisson-gamma coverage models.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param profileList A list containing the expected contribution matrix and paramters of the fitted profile.
#'
#' @return ggplot2-object.
ggplotQQPlotProfiles <- function(sampleTibble, profileList){
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
        mutate("DR" = .devianceResidualPoissonGammaDistribution(C, MU, MU / D),
               "Type" = if (NV == 1) "Noise coverage" else if (EC > 0) "Allele coverage" else NA) %>%
        filter(!is.na(Type)) %>% ungroup() %>% arrange(Type)

    qq_allele <- qqnorm((DRTibble %>% filter(Type == "Allele coverage"))$DR, plot.it = FALSE)
    qq_noise <- qqnorm((DRTibble %>% filter(Type == "Noise coverage"))$DR, plot.it = FALSE)
    p <- ggplot(DRTibble %>% mutate(X = c(qq_allele$x, qq_noise$x), Y = c(qq_allele$y, qq_noise$y)),
           aes_(x = ~X, y = ~Y)) +
        geom_point() + facet_wrap(~Type) + theme_bw() +
        xlab("Theoretical quantiles") + ylab("Sample quantiles") + ggtitle("Q-Q plot") +
        geom_smooth(method = "lm", se = FALSE)

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
        simulations[N, ] <- rnbinom(length(C), mu = MU, size = D)
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

