.rtnbinom = function(n, tau, mu, theta) {
    fill = rep(NA, n)
    for(i in 1:n){
        val = tau
        while(val <= tau){
            val = rnbinom(1, size = theta, mu = mu)
        }
        fill[i] = val
    }
    fill
}

#' Create genotype matrix
#'
#' @description Creates genotype matrix, given a sample and the true profiles of the contributors.
#'
#' @param coverageTibble A tibble containing sample information.
#' @param trueProfiles A list of tibbles, an element for each contributor.
#'
#' @return A matrix with the allele counts of of each contributor.
genotypeMatrix <- function(coverageTibble, trueProfiles) {
    all(unique(bind_rows(trueProfiles)$Region) %in% coverageTibble$Region)

    markersInProfiles <- sort(as.character(unique(do.call(rbind, trueProfiles)$Marker)))

    allelesOnMarker <- coverageTibble %>% ungroup() %>% group_by_(~Marker) %>% arrange_(~Marker) %>% summarise(Count = n())
    partialSumAlleles <- cumsum(c(0, allelesOnMarker$Count))

    nMatrix <- matrix(0, ncol = length(trueProfiles), nrow = dim(coverageTibble)[1])
    for (m in seq_along(markersInProfiles)) {
        for (i in seq_along(trueProfiles)) {
            contributor_i <- trueProfiles[[i]] %>% filter_(~Marker == markersInProfiles[m]) %>% mutate(Genotype = 3 - n())
            coverageTibble_mi <- coverageTibble %>% filter_(~Marker == markersInProfiles[m])

            j <- partialSumAlleles[m] + which(coverageTibble_mi$Region %in% contributor_i$Region)
            nMatrix[j, i] <- contributor_i$Genotype
        }
    }

    return(nMatrix)
}

.cyclicRotation <- function(x, y) {
    (nchar(y) == nchar(x)) && (grepl(y, strrep(x, 2), fixed = TRUE))
}

.getEntireRepeatStructure <- function(s, motifLength = 4) {
    motifLength <- if (!is.integer(motifLength)) as.integer(motifLength) else motifLength

    sD <- DNAString(s)

    typesOfMotifs <- oligonucleotideFrequency(sD, motifLength)
    typesOfMotifs <- typesOfMotifs[typesOfMotifs > 0]
    motifs <- names(typesOfMotifs)

    positionOfMotifs <- lapply(as.list(motifs), function(y, text = s) unlist(gregexpr2(y, text = text)))

    allRepeats <- structure(vector("list", length(positionOfMotifs)), .Names = motifs)
    for (i in seq_along(positionOfMotifs)) {
        y = positionOfMotifs[[i]]
        rleValues <- rle(y-(motifLength*(0:(length(y) - 1))))
        end <- y[cumsum(rleValues$lengths)] + motifLength

        startEndFrame <- data.frame(Start = end - rleValues$length*motifLength, End = end, Repeats = rleValues$length)
        j = 1
        while (j <= dim(startEndFrame)[1]) {
            whichExtending <- which(startEndFrame$Start == startEndFrame$End[j])
            if (length(whichExtending) > 0) {
                startEndFrame$End[j] <- startEndFrame$End[whichExtending]
                startEndFrame$Repeats[j] <- startEndFrame$Repeats[j] + startEndFrame$Repeats[whichExtending]
                startEndFrame <- startEndFrame[-whichExtending, ]
            }

            j = j + 1
        }

        allRepeats[[i]] <- startEndFrame
    }

    if (length(allRepeats) == 0) {
        allRepeats <- list("NA" = data.frame(Start = NA, End = NA, Repeats = NA))
    }
    allRepeats <- enframe(allRepeats, name = "Motif") %>% unnest()

    reducedRepeats <- allRepeats
    i = 1
    while ((i <= dim(reducedRepeats)[1]) & !is.na(reducedRepeats$Repeats[1])) {
        whichWithin <- rep(FALSE, dim(reducedRepeats)[1])
        whichWithin[-i] <- IRanges(reducedRepeats$Start[-i], reducedRepeats$End[-i]) %within% IRanges(reducedRepeats$Start[i], reducedRepeats$End[i])

        reducedRepeats <- reducedRepeats[!whichWithin, ]
        i = i + 1
    }

    return(reducedRepeats)
}

#' Potential parents list
#'
#' @description Creates a list of the potential parents for every region in the provided sample.
#'
#' @param coverageTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param stutterRatioModel A linear model fit; created by the \link{lm}-function.
#' @param trace If 'TRUE' adds a simple marker trace.
#' @param simplifiedReturn Should the returned list be simplified (TRUE/FALSE)?
#'
#' @return A list of the potential parents.
potentialParents <- function(coverageTibble, stutterRatioModel = NULL, trace = FALSE, simplifiedReturn = TRUE) {
    gapOpeningPenalty = 6
    gapExtensionPenalty = 1

    markersInProfiles <- sort(as.character(unique(coverageTibble$Marker)))

    potentialParentsAll <- structure(vector("list", length(markersInProfiles)), .Names = markersInProfiles)
    for (m in markersInProfiles) {
        if (trace)
            cat("Marker:", m, ":: ", which(markersInProfiles == m), "/", length(markersInProfiles), "\n")

        coverageTibble_m <- coverageTibble %>% filter_(~Marker == m)
        entireParentRepeatStructure_m <- vector("list", length(coverageTibble_m$Region))

        potentialParentsAll_m <- structure(vector("list", dim(coverageTibble_m)[1]), .Names = coverageTibble_m$Region)
        for (i in 1:dim(coverageTibble_m)[1]) {
            coverageTibble_mi <- coverageTibble_m[i, ]
            res_m <- coverageTibble_m %>% mutate(PotentialParent = 0, StutterRatio = 0)

            whichPotentialParents <- which(coverageTibble_m$Allele == coverageTibble_mi$Allele + 1)
            alignedPotentialParents <- rep(FALSE, length(whichPotentialParents))
            stutterRatioPotentialParents <- rep(0, length(whichPotentialParents))

            if (!is.null(stutterRatioModel)) {
                for (j in seq_along(whichPotentialParents)) {
                    subMatrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -nchar(coverageTibble_m$Region[whichPotentialParents[j]]), baseOnly = T)
                    stutterAligned <- pairwiseAlignment(DNAStringSet(as.character(coverageTibble_mi$Region)), coverageTibble_m$Region[whichPotentialParents[j]], substitutionMatrix = subMatrix,
                                                        gapOpening = -gapOpeningPenalty, gapExtension = -gapExtensionPenalty)

                    alignedPotentialParents[j] <- (stutterAligned@score == (nchar(coverageTibble_m$Region[whichPotentialParents[j]]) - coverageTibble_mi$MotifLength - (gapOpeningPenalty + coverageTibble_mi$MotifLength*gapExtensionPenalty)))

                    if (alignedPotentialParents[j]) {
                        if (is.null(entireParentRepeatStructure_m[[whichPotentialParents[j]]])) {
                            entireParentRepeatStructure_m[[whichPotentialParents[j]]] <- .getEntireRepeatStructure(coverageTibble_m$Region[whichPotentialParents[j]], round(coverageTibble_m$MotifLength[whichPotentialParents[j]]))
                        }
                        entireParentRepeatStructure <- entireParentRepeatStructure_m[[whichPotentialParents[j]]]

                        missingRepeatUnitStartPosition <- which(unlist(strsplit(as.character(aligned(stutterAligned)), "")) == "-")[1]
                        entireParentRepeatStructure_k <- entireParentRepeatStructure[which((missingRepeatUnitStartPosition >= entireParentRepeatStructure$Start) & (missingRepeatUnitStartPosition < entireParentRepeatStructure$End)),]
                        endingMotif <- entireParentRepeatStructure_k$Motif[which(entireParentRepeatStructure_k$End == (missingRepeatUnitStartPosition + coverageTibble_mi$MotifLength))]

                        occurenceInParent <- entireParentRepeatStructure_k$Repeats
                        motifCycles <- sapply(entireParentRepeatStructure_k$Motif, function(m) .cyclicRotation(endingMotif, m))

                        BLMM <- max(occurenceInParent[motifCycles])
                        stutterRatioPotentialParents[j] <- unname(suppressWarnings(predict.lm(stutterRatioModel, newdata = tibble(Marker = m, BlockLengthMissingMotif = BLMM))))
                    }
                }

                whichPotentialParents <- whichPotentialParents[alignedPotentialParents]
                res_m$PotentialParent[whichPotentialParents] <- 1
                res_m$StutterRatio[whichPotentialParents] <- stutterRatioPotentialParents[alignedPotentialParents]
            }

            if (simplifiedReturn) {
                if (sum(res_m$PotentialParent) > 0) {
                    res_m <- res_m %>% mutate(Index = 1:n()) %>% filter_(~PotentialParent == 1) %>%
                        select_(~PotentialParent, ~Index, ~StutterRatio) %>% as.matrix()
                }
                else {
                    res_m <- matrix(c(0, -1, 0), nrow = 1)
                    colnames(res_m) <- c("PotentialParent", "Index", "StutterRatio")
                }

            }
            potentialParentsAll_m[[i]] <- res_m
        }

        potentialParentsAll[[m]] <- potentialParentsAll_m
    }

    return(potentialParentsAll)
}

#' Multi-core potential parents
#'
#' @description Simplified multi-core implementation of the potential parents function.
#'
#' @param coverageTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param stutterRatioModel A linear model fit; created by the \link{lm}-function.
#' @param numberOfThreads The maximum number of threads allowed.
#'
#' @return A list of the potential parents.
potentialParentsMultiCore <- function(coverageTibble, stutterRatioModel, numberOfThreads = 4) {
    markersInProfiles <- sort(as.character(unique(coverageTibble$Marker)))

    potentialParentsAll <- mclapply(markersInProfiles, function(mm) {
        coverageTibble_m <- coverageTibble %>% filter_(~Marker == mm)

        potentialParentsAll_m <- structure(vector("list", dim(coverageTibble_m)[1]), .Names = coverageTibble_m$Region)
        for (i in 1:dim(coverageTibble_m)[1]) {
            coverageTibble_mi <- coverageTibble_m[i, ]

            whichPotentialParents <- which(coverageTibble_m$Allele == coverageTibble_mi$Allele + 1)
            alignedPotentialParents <- rep(FALSE, length(whichPotentialParents))
            stutterRatioPotentialParents <- rep(0, length(whichPotentialParents))

            if (!is.null(stutterRatioModel)) {
                for (j in seq_along(whichPotentialParents)) {
                    child <- as.character(coverageTibble_mi$Region)
                    parent <- coverageTibble_m$Region[whichPotentialParents[j]]
                    motifLength <- round(coverageTibble_mi$MotifLength)

                    x = 1
                    stringsMatch = TRUE
                    while (stringsMatch) {
                        if ((stringr::str_sub(child, x, x) == stringr::str_sub(parent, x, x)) & (x <= nchar(child))) {
                            x = x + 1
                        }
                        else {
                            stringsMatch = FALSE
                        }
                    }

                    missingRepeatUnitStartPosition <- x - 1
                    if (missingRepeatUnitStartPosition == nchar(child)) {
                        alignedPotentialParents[j] = TRUE
                    } else if (stringr::str_sub(child, x - 1) == stringr::str_sub(parent, x - 1 + motifLength)) {
                        alignedPotentialParents[j] = TRUE
                    }


                    if (alignedPotentialParents[j]) {
                        endingMotif <- stringr::str_sub(parent, x, x - 1 + motifLength)

                        BLMM <- 1
                        moreInBlock <- TRUE
                        while (moreInBlock) {
                            adjacentMotif <- stringr::str_sub(parent, max(1, x - BLMM * motifLength), x - 1 + motifLength - BLMM * motifLength)
                            if (adjacentMotif == endingMotif) {
                                BLMM = BLMM + 1
                            }
                            else {
                                moreInBlock = FALSE
                            }
                        }

                        stutterRatioPotentialParents[j] <- unname(suppressWarnings(predict.lm(stutterRatioModel, newdata = data.frame(Marker = mm, BlockLengthMissingMotif = BLMM))))
                    }
                }

                whichPotentialParents <- whichPotentialParents[alignedPotentialParents]
                if (length(whichPotentialParents) > 0) {
                    res_m <- matrix(c(as.numeric(alignedPotentialParents[alignedPotentialParents]),
                                      whichPotentialParents,
                                      stutterRatioPotentialParents[alignedPotentialParents]), ncol = 3, byrow = F)
                }
                else {
                    res_m <- matrix(c(0, -1, 0), nrow = 1)
                }

                potentialParentsAll_m[[i]] <- res_m
            }
        }

        return(potentialParentsAll_m)
    }, mc.cores = numberOfThreads)

    return(potentialParentsAll)
}

.expectedContributionMatrix <- function(coverageTibble, genotypeMatrix_s, potentialParentsList_s) {
    coverageTibble_s <- coverageTibble %>%
        left_join(coverageTibble %>% distinct_(~Marker) %>% mutate(MarkerNumeric = 1:n()), by = "Marker")

    resMatrix <- do.call(cbind, lapply(1:dim(genotypeMatrix_s)[2], function(contributor) {
        coverageTibble_c <- coverageTibble_s %>% mutate(Contributor = genotypeMatrix_s[, contributor])

        PotentialParentsContribution <- list()
        ContributorGenotype <- list(coverageTibble_c$Contributor)
        for(p in seq_along(potentialParentsList_s)) {
            potentialParentsList_sp <- potentialParentsList_s[[p]]
            coverageTibble_cm <- coverageTibble_c %>% filter_(~MarkerNumeric == p)

            parentContribution = rep(0, length(potentialParentsList_sp))
            for (pp in seq_along(potentialParentsList_sp)) {
                potentialParentsListIndex <- potentialParentsList_sp[[pp]][, 2]
                potentialParentsListStutterRatio <- potentialParentsList_sp[[pp]][, 3]

                if (potentialParentsListIndex[1] != -1) {
                    for (ppp in seq_along(potentialParentsListIndex)) {
                        parentContribution[pp] <- parentContribution[pp] + coverageTibble_cm$Contributor[potentialParentsListIndex[ppp]] * potentialParentsListStutterRatio[ppp]
                    }

                }
            }

            PotentialParentsContribution[[p]] <- parentContribution
        }

        res <- unname(as.matrix(
            tibble(ContributorGenotype = do.call(c, ContributorGenotype),
                   PotentialParentContribution = do.call(c, PotentialParentsContribution)) %>%
                mutate_(ExpectedContribution = ~ContributorGenotype + PotentialParentContribution) %>% select_(~ExpectedContribution)))
        return(res)
    }))

    return(resMatrix)
}


.encodeProfile <- function(profile, numberOfContributors, numberOfMarkers, numberOfAlleles) {
    partialSumAlleles <- .partialSumEigen(numberOfAlleles)
    encodedProfile <- lapply(seq_along(numberOfAlleles), function(a) {
        profile_a <- as.matrix(profile[(partialSumAlleles[a] + 1):partialSumAlleles[a + 1], ])
        possibleValues <- (1:numberOfAlleles[a]) - 1
        encoded <- unlist(lapply(1:dim(profile_a)[2], function(j) {
            rep(possibleValues, times = profile_a[, j])
        }))
    })

    if(!all(unlist(lapply(encodedProfile, function(e) length(e) == 2*numberOfContributors)))) {
        stop("Encoding failed.")
    }

    encodedProfile <- matrix(unlist(encodedProfile), ncol = 1)
    return(encodedProfile)
}

.alleleCoverageParameters.control <- function(nu, eta, phi) {
    return(list(nu = ifelse(is.null(nu), 1000, nu), eta = ifelse(is.null(eta), 2, eta), phi = if(is.null(phi)) c(0.7, 0.3) else phi))
}

.noiseParameters.control <- function(psi, rho) {
    return(list(psi = ifelse(is.null(psi), 2, psi), rho = ifelse(is.null(rho), 2, rho)))
}


#' Generate an allelic ladder
#'
#' @description Generates an allelic ladder / population data-base given a set of markers, (STR) regions, frequencies, and the motifLengths.
#'
#' @param markers A vector of the typed markers in the population.
#' @param regions A vector, or a list of vectors, of possible STR regions in the population.
#' @param frequencies A vector or a list of vectors, of the allele frequencies in the population.
#' @param motifLength A vector of the motif lengths of the markers.
#'
#' @details The 'regions', 'frequencies', and marker list/vectors should have the same length.
#'
#' The regions should be given as either full or compressed STR regions. A compressed region is given as e.g. '[AATG]x[CGTT]y' corresponding to 'AATG' repeated 'x' times followed by 'CGTT' repeated 'y' times.
#' Note if the region has motifs which are not repeated e.g. '[AATG]xCTT', the non-repeated regions should also be encapsuled in square brackets, i.e. '[AATG]x[CTT]'.
#'
#' The 'motifLength' parameter needs length 1, length equal to the length of 'markers', or 'NULL'. If 'NULL' it is assumed to be 4.
#'
#' @return A \link{tibble}.
#' @example inst/examples/generateExamplePopulation.R
generateAllelicLadder <- function(markers, regions, frequencies, motifLength = NULL) {
    if (length(unlist(regions)) != length(unlist(frequencies)))
        stop("The length of 'regions' and 'frequencies' must be equal.")

    if (is.null(motifLength))
        motifLength = rep(4, length(markers))
    else if (length(motifLength) == 1)
        motifLength = rep(motifLength, length(markers))

    if (length(markers) != length(motifLength))
        stop("The length of 'motifLength' must be equal to the length of 'markers'.")

    if (is.list(regions) & is.list(frequencies)) {
        if (any(sapply(regions, length) != sapply(frequencies, length)))
            stop("The length of 'regions' and 'frequencies' must be equal.")

        if (length(markers) != length(regions))
            stop("The length of 'markers' must be equal to the length of 'frequencies' and 'regions'.")

        markers <- rep(markers, times = sapply(regions, length))
        motifLength <- rep(motifLength, times = sapply(regions, length))
        regions <- unlist(regions)
        frequencies <- unlist(frequencies)
    }
    else if (is.list(regions)) {
        if (length(regions) != length(markers))
            stop("The length of 'markers' must be equal to the length of 'frequencies' and 'regions'.")

        markers <- rep(markers, times = sapply(regions, length))
        motifLength <- rep(motifLength, times = sapply(regions, length))
        regions <- unlist(regions)
    }
    else if (is.list(frequencies)) {
        if (length(frequencies) != length(markers))
            stop("The length of 'markers' must be equal to the length of 'frequencies' and 'regions'.")

        markers <- rep(markers, times = sapply(frequencies, length))
        motifLength <- rep(motifLength, times = sapply(frequencies, length))
        frequencies <- unlist(frequencies)
    }

    decompressString <- function(s) {
        s_split <- str_split(s, "\\[", simplify = TRUE)
        decompressed <- sapply(s_split[, 2:dim(s_split)[2]], function(i) {
            ss <- unlist(str_split(i, "\\]"))
            return(paste(rep(ss[1], times = ifelse(is.na(as.numeric(ss[2])), 1, as.numeric(ss[2]))), collapse = ""))
        })

        return(paste(decompressed, collapse = ""))
    }

    compressedStrings <- str_detect(regions, "\\[")
    regions <- sapply(seq_along(compressedStrings), function(i) if (compressedStrings[i]) decompressString(regions[i]) else regions[i])

    populationTibble <- tibble(Marker = markers, MotifLength = motifLength,
                               Allele = str_length(regions) / motifLength,
                               Region = regions,
                               AlleleFrequencies = frequencies) %>%
        group_by_(~Marker, ~Allele, ~MotifLength, ~Region) %>%
        summarise_(AlleleFrequencies = ~sum(AlleleFrequencies))
    return(populationTibble)
}

#' Sample MPS coverage
#'
#' @description Sample coverage under the Poisson-gamma model of Vilsen et. al 2017
#'
#' @param trueProfiles A list of \link{tibble}'s containing the true profiles of the mixture.
#' @param markerImbalances A vector of the marker imbalances.
#' @param populationLadder An allelic ladder containing the possible STR regions of each marker and their allele frequencies.
#' @param stutterRatioModel A linear model of the relationship between BLMM and stutter ratio for each marker in the data.
#' @param alleleCoverageParameters A list containing the parameters of the allele coverage model.
#' @param noiseParameters A list containing the paramters of the noise model.
#' @param p A probability threshold for the created stutters.
#' @param numberOfThreads The number of threads passed to the \link{potentialParentsMultiCore} function.
#'
#' @return A \link{tibble} with the sampled coverage, the marker, and the marker imbalances.
sampleCoverage <- function(trueProfiles, markerImbalances, populationLadder, stutterRatioModel,
                           alleleCoverageParameters = list(), noiseParameters = list(),
                           p = NULL, numberOfThreads = 4) {
    alleleCoverageParameters <- .alleleCoverageParameters.control(nu = alleleCoverageParameters$nu, eta = alleleCoverageParameters$eta, phi = alleleCoverageParameters$phi)
    noiseParameters <- .noiseParameters.control(psi = noiseParameters$psi, rho = noiseParameters$rho)

    p <- ifelse(is.null(p), 0.025, p)

    markers <- sort(unique(populationLadder$Marker))

    numberOfContributors <- length(trueProfiles)
    if (length(alleleCoverageParameters$phi) != numberOfContributors)
        stop("The number elements of 'phi' should be equal to the length of the 'trueProfiles' list.")

    ## Create regions
    sampleTibble_m <- vector("list", length(markers))
    for (m in seq_along(markers)) {
        trueProfiles_m <- lapply(trueProfiles, function(tp) tp %>% filter_(~Marker == markers[m]))
        populationLadder_m <- populationLadder %>% filter_(~Marker == markers[m])

        collectedProfiles <- bind_rows(trueProfiles_m) %>% select_(~-Name) %>%
            distinct_(~Region, .keep_all = T) %>% mutate(IsAllele = 1, IsStutter = 0)

        collectedProfilesStutters <- bind_rows(lapply(1:nrow(collectedProfiles), function(i) {
            BLMMs_i <- .getEntireRepeatStructure(collectedProfiles$Region[i], collectedProfiles$MotifLength[i])
            possibleStutters_i <- sapply(1:nrow(BLMMs_i), function(j) {
                paste(str_sub(collectedProfiles$Region[i], start = 1, end = BLMMs_i[j, 2] - 1),
                      str_sub(collectedProfiles$Region[i], start =  BLMMs_i[j, 2] + collectedProfiles$MotifLength[i]),
                      sep = "")
            })

            BLMMs_i <- BLMMs_i[which(!duplicated(possibleStutters_i)), ]
            possibleStutters_i <- possibleStutters_i[which(!duplicated(possibleStutters_i))]

            res <- populationLadder_m[which(populationLadder_m$Region %in% possibleStutters_i), ]
            if (dim(res)[1] != 0) {
                Marker = BlockLengthMissingMotif = NULL
                res <- res %>% mutate("BlockLengthMissingMotif" = BLMMs_i$Repeats[which(possibleStutters_i %in% populationLadder_m$Region)],
                                      "IsStutter" = 1,
                                      "IsAllele" = predict.lm(stutterRatioModel, newdata = data.frame("Marker" = Marker, "BlockLengthMissingMotif" = BlockLengthMissingMotif)))
            }
        }))

        profileStutters <- bind_rows(collectedProfiles, collectedProfilesStutters) %>%
            group_by_(~Marker, ~MotifLength, ~Allele, ~Region, ~AlleleFrequencies) %>%
            summarise_(IsAllele = ~max(IsAllele),
                      IsStutter = ~sum(IsStutter),
                      BlockLengthMissingMotif = ~sum(BlockLengthMissingMotif, na.rm = TRUE)) %>%
            ungroup() %>%
            filter_(~IsAllele > p)

        sampleTibble_m[[m]] <- profileStutters
    }

    sampleTibble <- bind_rows(sampleTibble_m) %>% arrange_(~Marker, ~Allele, ~Region)

    ## Sample coverage
    profilesMatrix <- genotypeMatrix(sampleTibble, trueProfiles)
    potentialParentsList <- potentialParentsMultiCore(sampleTibble, stutterRatioModel, numberOfThreads)
    ECM <- .expectedContributionMatrix(sampleTibble, profilesMatrix, potentialParentsList)

    numberOfAlleles <- (sampleTibble %>% group_by_(~Marker) %>% summarise(Count = n()))$Count
    markerImbalancesRep <- rep(markerImbalances, times = numberOfAlleles)
    expectedAlleleCoverage <- alleleCoverageParameters$nu * markerImbalancesRep * (ECM %*% alleleCoverageParameters$phi)

    sampledAllele <- sampleTibble %>%
        mutate("Coverage" = rnbinom(n(), mu = expectedAlleleCoverage, size = expectedAlleleCoverage / alleleCoverageParameters$eta))

    sampledNoise <- populationLadder %>% filter_(~!(Region %in% sampleTibble$Region)) %>% group_by_(~Marker) %>%
        mutate(Coverage = rnbinom(n(), mu = noiseParameters$psi, size = noiseParameters$psi / noiseParameters$rho)) %>%
        filter_(~Coverage > 0)

    sampledCoverage <- bind_rows(sampledAllele, sampledNoise) %>%
        arrange_(~Marker, ~Allele, ~Region) %>%
        group_by_(~Marker, ~MotifLength, ~Allele, ~AlleleFrequencies, ~Region) %>%
        summarise_(Coverage = ~sum(Coverage)) %>%
        left_join(tibble(Marker = markers, MarkerImbalance = markerImbalances), by = "Marker") %>%
        ungroup()

    return(sampledCoverage)
}


#' Sample genotypes
#'
#' @description Sample genotypes, given the number of profiles and a population, assuming the population was in HWE.
#'
#' @param numberOfContributors The number of profiles which should be sampled.
#' @param populationLadder A \link{tibble} containing the markers, regions, and allele frequencies of the population.
#' @param contributorNames A vector of names. The length should be the number of profiles.
#'
#' @return A list of tibbles, one for each sampled profile.
sampleGenotypesHWE <- function(numberOfContributors, populationLadder, contributorNames = NULL) {
    if (!is.null(contributorNames)) {
        if (length(contributorNames) == 1)
            contributorNames = rep(contributorNames, numberOfContributors)

        if (numberOfContributors != length(contributorNames))
            stop("The length of 'contributorNames' must equal to the number of contributors, or 1.")
    }
    else {
        contributorNames <- paste("C", 1:numberOfContributors, sep = "")
    }

    markers <- unique(populationLadder$Marker)

    sampledProfiles <- vector("list", numberOfContributors)
    for(cc in 1:numberOfContributors) {
        sampledProfiles_mm <- vector("list", length(markers))
        for (mm in seq_along(markers)) {
            populationLadder_m <- populationLadder %>% filter_(~Marker == markers[mm])
            numberOfAlleles_m <- dim(populationLadder_m)[1]

            ## Set-up probabilities
            outerProduct <- outer(populationLadder_m$AlleleFrequencies, populationLadder_m$AlleleFrequencies, "*")

            homozygoteProbability <- diag(diag(outerProduct))
            heterozygoteProbability <- 2 * outerProduct
            heterozygoteProbability[lower.tri(heterozygoteProbability, diag = T)] <- 0

            partialSumStacked <- .partialSumEigen(matrix(homozygoteProbability + heterozygoteProbability, ncol = 1))

            ## Draw genotype
            randomDraw <- runif(1)
            whichCombination <- which.max(partialSumStacked >= randomDraw) - 1
            firstAllele <-  ifelse(whichCombination %% numberOfAlleles_m == 0, numberOfAlleles_m, whichCombination %% numberOfAlleles_m)
            secondAllele <- ceiling(whichCombination / numberOfAlleles_m)

            sampledProfiles_mm[[mm]] <- populationLadder_m[c(firstAllele, secondAllele), ] %>% mutate(Name = contributorNames[cc])
        }

        sampledProfiles[[cc]] <- bind_rows(sampledProfiles_mm) %>% group_by_(~Marker) %>%
            distinct_(~Region, .keep_all = TRUE) %>% arrange_(~Marker, ~Region)
    }

    sampledProfiles <- structure(sampledProfiles, .Names = contributorNames)
    return(sampledProfiles)
}


#' Collect coverage and population information
#'
#' @description Collects a coverage tibble and the population ladder, by the 'Marker' column.
#'
#' @param coverageTibble A \link{tibble} containing the coverage, the marker, and marker imbalance scalar of each allele in the sample.
#' @param populationLadder A \link{tibble} containing the markers, regions, and allele frequencies of the population.
#'
#' @return A joined tibble, having added any regions not seen in the 'coverageTibble' with a coverage set to 0, and a column of allele frequencies from the 'populationLadder'.
collectSamplePopulation <- function(coverageTibble, populationLadder) {
    markers <- unique(coverageTibble$Marker)
    res <- vector("list", length(markers))
    for (mm in seq_along(markers)) {
        populationLadder_m = populationLadder %>% filter_(~Marker == markers[mm])
        coverageTibble_m = coverageTibble %>% filter_(~Marker == markers[mm])

        populationLadderCoverageLimited <- (populationLadder_m %>% mutate(Coverage = 0))[-which(populationLadder_m$Region %in% coverageTibble_m$Region), ]
        columnExtention <- unique(c(which(names(coverageTibble_m) == "Marker"), which(!(names(coverageTibble_m) %in% names(populationLadderCoverageLimited)))))

        res[[mm]] <- bind_rows(coverageTibble_m, populationLadderCoverageLimited %>%
                             left_join(coverageTibble_m %>% select_(~columnExtention) %>%
                                           distinct_(~Marker, .keep_all = TRUE), by = "Marker")) %>%
            ungroup() %>% arrange_(~Marker, ~Allele, ~Region)
    }

    return(bind_rows(res))
}


#' An example population.
#'
#' An example populations generated by the \link{generateAllelicLadder}-function.
#'
#' @docType data
#' @keywords datasets
#' @name examplePopulation
#' @usage data(examplePopulation)
NULL
