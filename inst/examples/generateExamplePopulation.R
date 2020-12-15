library("stringr")
library("htmltab")
library("dplyr")

markers <- c("TPOX", "D3S1358", "D5S818", "CSF1PO", "D7S820", "TH01", "VWA",
             "D13S317", "D16S539", "D18S51")
regions <- lapply(seq_along(markers), function(mm) {
    htmlTable <- htmltab(
        doc = paste("https://strbase.nist.gov/str_", markers[mm], ".htm", sep = ""),
        which = 3
    )

    names(htmlTable) <- str_replace_all(names(htmlTable), fixed(" "), "")

    regions_mm <- htmlTable$RepeatStructure[str_detect(htmlTable$RepeatStructure, "\\[")]

    regionsRestructured <- sapply(str_split(regions_mm, " "), function(s) {
        s <- s[str_length(s) > 0]
        s <- unlist(str_split(unlist(str_split(s, "\\[")), "\\]"))
        s <- s[str_length(s) > 0]
        s <- unlist(str_split(str_replace(s, ".*?([0-9]+).*", "~\\1~"), "~"))
        s <- s[str_length(s) > 0]

        s <- ifelse(!str_detect(s, ".*?([0-9]+).*"), paste("[", s, "]", sep = ""), s)
        return(paste(s, collapse = ""))
    })

    return(regionsRestructured)
})

frequencies <- lapply(regions, function(r) {
    gammaDraw <- rgamma(length(r), shape = 1, rate = 4)
    gammaDraw / sum(gammaDraw)
})

examplePopulation <- generateAllelicLadder(markers, regions, frequencies)
