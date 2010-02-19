.similarSeqnameConvention <- function(seqs1, seqs2) {
    funList <-
      list(isRoman = function(x) grepl("[iv]", tolower(x)),
           isArabic = function(x) grepl("[1-9]", tolower(x)),
           haschr = function(x) grepl("^chr", x),
           hasChr = function(x) grepl("^Chr", x),
           hasCHR = function(x) grepl("^CHR", x),
           haschrom = function(x) grepl("^chrom", x),
           hasChrom = function(x) grepl("^Chrom", x),
           hasCHROM = function(x) grepl("^CHROM", x))
    all(sapply(funList, function(f) any(f(seqs1)) == any(f(seqs2))))
}


setMethod("findOverlaps", c("GRangesList", "GRanges"),
    function(query, subject, maxgap = 0, multiple = TRUE,
             type = c("any", "start", "end", "within", "equal"))
    {
        if (is(query, "Alignments1"))
            query <- as(query, "GRangesList")

        DIM <- c(length(subject), length(query))
        if (DIM[1L] == 0L || DIM[2L] == 0L) {
            matchMatrix <-
              matrix(integer(), ncol = 2,
                     dimnames = list(NULL, c("query", "subject")))
        } else {
            queryGroup <- togroup(query@partitioning)
            queryUnlistedRanges <- query@unlistData@ranges
            querySeqnames <- as.character(query@unlistData@seqnames)
            uniqueQuerySeqnames <- unique(querySeqnames)
            queryStrand <- as.character(query@unlistData@strand)

            subjectSeqnames <- as.character(seqnames(subject))
            uniqueSubjectSeqnames <- unique(subjectSeqnames)
            subjectStrand <- as.character(strand(subject))

            if (!.similarSeqnameConvention(uniqueQuerySeqnames,
                    uniqueSubjectSeqnames))
                stop("'query' and 'subject' do not use a similiar naming ",
                     "convention for seqnames")

            uniqueSeqnames <-
              sort(intersect(uniqueQuerySeqnames, uniqueSubjectSeqnames))
            matrixList <- vector("list", length = length(uniqueSeqnames) * 2L)
            i <- 1L
            for (seqnm in uniqueSeqnames) {
                queryWant0 <- querySeqnames == seqnm
                tarWant0 <- subjectSeqnames == seqnm
                for (strd in c("+", "-")) {
                    queryWant <- which(queryWant0 & (queryStrand == strd))
                    queryRanges <- queryUnlistedRanges[queryWant]
                    queryIdx <- queryGroup[queryWant]

                    tarIdx <- which(tarWant0 & (subjectStrand == strd))
                    tarRanges <- ranges(subject[tarIdx])

                    overlaps <- callGeneric(queryRanges, tarRanges)
                    queryIdx <- queryIdx[queryHits(overlaps)]
                    subjectIdx <- tarIdx[subjectHits(overlaps)]
                    matrixList[[i]] <-
                      cbind(query = queryIdx, subject = subjectIdx)
                    i <- i + 1L
                }
            }
            matchMatrix <- do.call(rbind, matrixList)
            matchMatrix <-
              matchMatrix[order(matchMatrix[ , 1L], matchMatrix[ , 2L]), ,
                          drop=FALSE]
            matchMatrix <- matchMatrix[!duplicated(matchMatrix), , drop=FALSE]
        }
        new("RangesMatching", matchMatrix = matchMatrix, DIM = DIM)
    }
)
