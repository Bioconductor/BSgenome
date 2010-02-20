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


setMethod("findOverlaps", c("GRanges", "GRanges"),
    function(query, subject, maxgap = 0, multiple = TRUE,
             type = c("any", "start", "end", "within", "equal"))
    {
        DIM <- c(length(subject), length(query))
        if (DIM[1L] == 0L || DIM[2L] == 0L) {
            matchMatrix <-
              matrix(integer(), ncol = 2,
                     dimnames = list(NULL, c("query", "subject")))
        } else {
            querySeqnames <- as.character(seqnames(query))
            uniqueQuerySeqnames <- unique(querySeqnames)
            queryStrand <- as.character(strand(query))

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
                subjectWant0 <- subjectSeqnames == seqnm
                for (strd in c("+", "-")) {
                    queryIdx <- which(queryWant0 & (queryStrand == strd))
                    subjectIdx <- which(subjectWant0 & (subjectStrand == strd))

                    queryRanges <- ranges(query[queryIdx])
                    subjectRanges <- ranges(subject[subjectIdx])

                    overlaps <-
                      callGeneric(queryRanges, subjectRanges,
                                  maxgap = maxgap, multiple = multiple,
                                  type = type)

                    matrixList[[i]] <-
                      cbind(query = queryIdx[queryHits(overlaps)],
                            subject = subjectIdx[subjectHits(overlaps)])

                    i <- i + 1L
                }
            }
            matchMatrix <- do.call(rbind, matrixList)
            matchMatrix <-
              matchMatrix[order(matchMatrix[ , 1L], matchMatrix[ , 2L]), ,
                          drop=FALSE]
        }
        new("RangesMatching", matchMatrix = matchMatrix, DIM = DIM)
    }
)

setMethod("findOverlaps", c("GRangesList", "GRanges"),
    function(query, subject, maxgap = 0, multiple = TRUE,
             type = c("any", "start", "end", "within", "equal"))
    {
        ans <-
          callGeneric(unlist(query, use.names=FALSE), subject,
                      maxgap = maxgap, multiple = multiple, type = type)
        matchMatrix <- ans@matchMatrix
        matchMatrix[, 1L] <- togroup(query@partitioning)[matchMatrix[, 1L]]
        matchMatrix <- matchMatrix[!duplicated(matchMatrix), , drop=FALSE]
        DIM <- c(length(subject), length(query))
        initialize(ans, matchMatrix = matchMatrix, DIM = DIM)
    }
)

setMethod("findOverlaps", c("GRanges", "GRangesList"),
    function(query, subject, maxgap = 0, multiple = TRUE,
             type = c("any", "start", "end", "within", "equal"))
    {
        ans <-
          callGeneric(query, unlist(subject, use.names=FALSE),
                      maxgap = maxgap, multiple = multiple, type = type)
        matchMatrix <- ans@matchMatrix
        matchMatrix[, 2L] <- togroup(subject@partitioning)[matchMatrix[, 2L]]
        matchMatrix <- matchMatrix[!duplicated(matchMatrix), , drop=FALSE]
        DIM <- c(length(subject), length(query))
        initialize(ans, matchMatrix = matchMatrix, DIM = DIM)
    }
)

setMethod("findOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, maxgap = 0, multiple = TRUE,
             type = c("any", "start", "end", "within", "equal"))
    {
        ans <-
          callGeneric(unlist(query, use.names=FALSE),
                      unlist(subject, use.names=FALSE),
                      maxgap = maxgap, multiple = multiple, type = type)
        matchMatrix <- ans@matchMatrix
        matchMatrix[, 1L] <- togroup(query@partitioning)[matchMatrix[, 1L]]
        matchMatrix[, 2L] <- togroup(subject@partitioning)[matchMatrix[, 2L]]
        matchMatrix <- matchMatrix[!duplicated(matchMatrix), , drop=FALSE]
        DIM <- c(length(subject), length(query))
        initialize(ans, matchMatrix = matchMatrix, DIM = DIM)
    }
)
