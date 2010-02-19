findOverlaps1 <- function(pat, target)
{
    if (is(pat, "Alignments1")) pat <- as(pat, "GRangesList")
    if (!is(target, "GRanges")) stop("'target' must be a GRanges object")
    seqnames_target <- as.character(seqnames(target))
    strand_target <- as.character(strand(target))

    patgroup <- togroup(pat@partitioning)
    ulData <- pat@unlistData
    ourranges <- ulData@ranges
    seqnames_pat <- as.character(ulData@seqnames)
    strand_pat <- as.character(ulData@strand)

    uniqueSeqn <- sort(intersect(unique(seqnames_pat), unique(seqnames_target)))
    ans <- vector(mode = "list", length = length(uniqueSeqn) * 2L)
    i <- 1L
    for (seqn in uniqueSeqn) {
        patWant0 <- seqnames_pat == seqn
        tarWant0 <- seqnames_target == seqn
        for (strd in c("+", "-")) {
            patWant <- which(patWant0 & (strand_pat == strd))
            patRanges <- ourranges[patWant]
            patIdx <- patgroup[patWant]

            tarIdx <- which(tarWant0 & (strand_target == strd))
            tarRanges <- ranges(target[tarIdx])

            overlaps <- IRanges::findOverlaps(patRanges, tarRanges)
            queryIdx <- patIdx[queryHits(overlaps)]
            subjectIdx <- tarIdx[subjectHits(overlaps)]
            ans[[i]] <- cbind(query = queryIdx, subject = subjectIdx)
            i <- i + 1L
        }
    }
    DIM <- c(length(target), length(pat))
    m <- do.call(rbind, ans)
    new("RangesMatching",
        matchMatrix = m[order(m[ , 1L], m[ , 2L]), , drop=FALSE],
        DIM = DIM)
}
