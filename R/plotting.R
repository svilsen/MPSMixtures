#' Compare profile combinations
#'
#' @details Function used for visually comparing two profiles.
#'
#' @param sampleTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param profileList1 A list containing the expected contribution matrix and the mixture paramters of profile 1.
#' @param profileList2 A list containing the expected contribution matrix and the mixture paramters of profile 2.
#' @param plotDifferencesOnly Restrict the returned \link{ggplot}-object to markers where the two profiles mismatch.
#'
#' @return ggplot2-object.
ggplotComparingProfiles <- function(sampleTibble, profileList1, profileList2, plotDifferencesOnly = FALSE, profileNames = NULL) {
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
    partialSumProfiles <- partialSumEigen(numberOfProfiles)
    plotTibble <- vector("list", sum(numberOfProfiles) + 2)
    for (i in 1:2) {
        for (j in 1:numberOfProfiles[i]) {
            k = partialSumProfiles[i] + j
            plotTibble[[k]] <- sampleTibble %>%
                mutate(ProfileCoverage = Coverage * profiles[[i]][, j], Profile = paste("Profile:", profileNames[i]),
                       Contributor = paste0("C", j)) %>%
                filter(ProfileCoverage != 0)
        }

        plotTibble[[length(plotTibble) + i - 1]] <- sampleTibble %>%
            mutate(ProfileCoverage = Coverage * as.numeric(rowSums(profiles[[i]]) == 0), Profile = paste("Profile:", profileNames[i]),
                   Contributor = paste0("N")) %>%
            filter(ProfileCoverage != 0)
    }

    plotTibble <- bind_rows(plotTibble)

    gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
    }

    colours <- gg_color_hue(2)
    if (plotDifferencesOnly) {
        profileFrame <- tibble(P1 = apply(profileList1$ExpectedContributionMatrix, 1, sum),
                               P2 = apply(profileList2$ExpectedContributionMatrix, 1, sum),
                               Difference = P1 - P2)
        splitMarker <- split(profileFrame, sampleTibble$Marker)
        differenceMarker <- sapply(splitMarker, function(xx) sum(abs(xx[, 3])) != 0)

        plotTibble <- plotTibble %>% filter(Marker %in% names(differenceMarker)[differenceMarker])
    }

    if (dim(plotTibble)[1] > 0) {
        p.Profile <- ggplot(plotTibble, aes(x = Allele, y = ProfileCoverage, fill = Contributor)) + geom_bar(stat = "identity") +
            facet_grid(Marker ~ Profile, scales = "free_y") + scale_fill_manual(values = c(colours[1], colours[2], "#000000")) +
            theme_bw() + theme(legend.position = "top") + ylab("Coverage") +
            scale_x_continuous(breaks = seq(min(plotTibble$Allele), max(plotTibble$Allele)))
    }
    else {
        p.Profile <- ggplot()
    }

    return(p.Profile)
}
