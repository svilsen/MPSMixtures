genotypeMatrix <- function(coverageTibble, trueProfiles) {
    markersInProfiles <- sort(as.character(unique(do.call(rbind, trueProfiles)$Marker)))

    allelesOnMarker <- coverageTibble %>% ungroup() %>% group_by(Marker) %>% arrange(Marker) %>% summarise(Count = n())
    partialSumAlleles <- cumsum(c(0, allelesOnMarker$Count))

    nMatrix <- matrix(0, ncol = length(trueProfiles), nrow = dim(coverageTibble)[1])
    for (m in seq_along(markersInProfiles)) {
        for (i in seq_along(trueProfiles)) {
            contributor_i <- trueProfiles[[i]] %>% filter(Marker == markersInProfiles[m]) %>% mutate(Genotype = 3 - n())
            coverageTibble_mi <- coverageTibble %>% filter(Marker == markersInProfiles[m])

            j <- partialSumAlleles[m] + which(coverageTibble_mi$Region %in% contributor_i$Region)
            nMatrix[j, i] <- contributor_i$Genotype
        }
    }

    return(nMatrix)
}

# s = coverageTibble_m$Region[1]; motifLength = round(coverageTibble_m$MotifLength[1])
getEntireRepeatStructure <- function(s, motifLength = 4) {
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

# gapOpeningPenalty = 6; gapExtensionPenalty = 1; trace = T
potentialParents <- function(coverageTibble, stutterRatioModel = NULL, gapOpeningPenalty = 6, gapExtensionPenalty = 1, trace = FALSE) {
    markersInProfiles <- sort(as.character(unique(coverageTibble$Marker)))

    potentialParentsAll <- list()
    for (m in markersInProfiles) {
        if (trace)
            cat("Marker:", m, ":: ", which(markersInProfiles == m), "/", length(markersInProfiles), "\n")

        coverageTibble_m <- coverageTibble %>% filter(Marker == m)
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
                            entireParentRepeatStructure_m[[whichPotentialParents[j]]] <- getEntireRepeatStructure(coverageTibble_m$Region[whichPotentialParents[j]], round(coverageTibble_m$MotifLength[whichPotentialParents[j]]))
                        }
                        entireParentRepeatStructure <- entireParentRepeatStructure_m[[whichPotentialParents[j]]]

                        missingRepeatUnitStartPosition <- which(unlist(strsplit(as.character(aligned(stutterAligned)), "")) == "-")[1]
                        entireParentRepeatStructure_k <- entireParentRepeatStructure[which((missingRepeatUnitStartPosition >= entireParentRepeatStructure$Start) & (missingRepeatUnitStartPosition < entireParentRepeatStructure$End)),]
                        endingMotif <- entireParentRepeatStructure_k$Motif[which(entireParentRepeatStructure_k$End == (missingRepeatUnitStartPosition + coverageTibble_mi$MotifLength))]

                        occurenceInParent <- entireParentRepeatStructure_k$Repeats
                        motifCycles <- sapply(entireParentRepeatStructure_k$Motif, function(m) STRMPS:::.cyclicRotation(endingMotif, m))

                        BLMM <- max(occurenceInParent[motifCycles])
                        stutterRatioPotentialParents[j] <- unname(suppressWarnings(predict.lm(stutterRatioModel, newdata = tibble(Marker = m, BlockLengthMissingMotif = BLMM))))
                    }
                }

                whichPotentialParents <- whichPotentialParents[alignedPotentialParents]
                res_m$PotentialParent[whichPotentialParents] <- 1
                res_m$StutterRatio[whichPotentialParents] <- stutterRatioPotentialParents[alignedPotentialParents]
            }

            potentialParentsAll_m[[i]] <- res_m
        }

        potentialParentsAll[[m]] <- potentialParentsAll_m
    }

    return(potentialParentsAll)
}


encodeProfile <- function(profile, numberOfContributors, numberOfMarkers, numberOfAlleles) {
    partialSumAlleles <- partialSumEigen(numberOfAlleles)
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
