findOverlaps1 <- function(pat, target)
{
    if (is(pat, "Alignments1")) {
        pat <- as(pat, "GRangesList")
    }
    ans <- list()
    for (seqn in sort(unique(unlist(seqnames(pat))))) {
        for (strd in c("+", "-")) {
            patranges <- ranges(pat[seqnames(pat) == seqn & strand(pat) == strd])
            tarwant <- seqnames(target) == seqn & strand(target) == strd
            taridx <- which(tarwant)
            tarranges <- ranges(target[taridx])
            lapply(seq_len(length(patranges)), function(i)
               {
                   x <- IRanges::findOverlaps(patranges[[i]], tarranges)
                   subjectIdx <- taridx[subjectHits(x)]
                   ans <<- c(ans, list(cbind(query=rep(i, length(subjectIdx)),
                                             subject=subjectIdx)))
               })
        }
    }
    DIM <- c(length(target), length(pat))
    m <- do.call(rbind, ans)
    new("RangesMatching", matchMatrix=m[order(m[ , 1L], m[ , 2L]), ],
        DIM = DIM)
}
