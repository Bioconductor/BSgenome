findOverlaps1 <- function(pat, target)
{
    if (is(pat, "Alignments1")) {
        pat <- as(pat, "GRangesList")
    }
    if (!is(target, "GRanges")) {
        stop("'target' must be a GRanges object")
    }
    ans <- list()
    uniqueSeqn <-
      sort(intersect(unique(unlist(seqnames(pat), use.names=FALSE)),
                     unique(seqnames(target))))
    for (seqn in uniqueSeqn) {
        for (strd in c("+", "-")) {
            patWant <- seqnames(pat) == seqn & strand(pat) == strd
            patRanges <- unlist(ranges(pat[patWant]), use.names=FALSE)
            patIdx <-
              seqselect(togroup(pat@partitioning),
                        unlist(patWant, use.names=FALSE))

            tarWant <- seqnames(target) == seqn & strand(target) == strd
            tarIdx <- which(tarWant)
            tarRanges <- ranges(target[tarIdx])
            
            overlaps <- IRanges::findOverlaps(patRanges, tarRanges)
            queryIdx <- patIdx[queryHits(overlaps)]
            subjectIdx <- tarIdx[subjectHits(overlaps)]
            ans <- c(ans, list(cbind(query = queryIdx, subject = subjectIdx)))
        }
    }
    DIM <- c(length(target), length(pat))
    m <- do.call(rbind, ans)
    new("RangesMatching",
        matchMatrix = m[order(m[ , 1L], m[ , 2L]), , drop=FALSE],
        DIM = DIM)
}
