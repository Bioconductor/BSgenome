.similarSeqnameConvention <- function(seqs1, seqs2) {
    funList <-
      list(isRoman = function(x) grepl("[ivIV]", x),
           isArabic = function(x) grepl("[1-9]", x),
           haschr = function(x) grepl("^chr", x),
           hasChr = function(x) grepl("^Chr", x),
           hasCHR = function(x) grepl("^CHR", x),
           haschrom = function(x) grepl("^chrom", x),
           hasChrom = function(x) grepl("^Chrom", x),
           hasCHROM = function(x) grepl("^CHROM", x))
    all(sapply(funList, function(f) any(f(seqs1)) == any(f(seqs2))))
}

.cleanMatchMatrix <- function(matchMatrix) {
    fastDiff <- IRanges:::diffWithInitialZero
    nr <- nrow(matchMatrix)
    nc <- ncol(matchMatrix)
    if (nr <= 1L) {
        matchMatrix
    } else {
        matchMatrix <-
          matchMatrix[IRanges:::orderTwoIntegers(matchMatrix[ , 1L, drop=TRUE],
                                                 matchMatrix[ , 2L, drop=TRUE]), ,
                      drop=FALSE]
        matchMatrix[fastDiff(matchMatrix[,1L,drop=TRUE]) != 0L |
                    fastDiff(matchMatrix[,2L,drop=TRUE]) != 0L, , drop=FALSE]
    }
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
            querySeqnames <- seqnames(query)
            querySplitRanges <- splitRanges(querySeqnames)
            uniqueQuerySeqnames <- names(querySplitRanges)

            subjectSeqnames <- seqnames(subject)
            subjectSplitRanges <- splitRanges(subjectSeqnames)
            uniqueSubjectSeqnames <- names(subjectSplitRanges)

            if (!.similarSeqnameConvention(uniqueQuerySeqnames,
                                           uniqueSubjectSeqnames))
                stop("'query' and 'subject' do not use a similiar naming ",
                     "convention for seqnames")

            commonSeqnames <-
              intersect(uniqueQuerySeqnames, uniqueSubjectSeqnames)

            queryStrand <- as.character(strand(query))
            queryRanges <- unname(ranges(query))

            subjectStrand <- as.character(strand(subject))
            subjectRanges <- unname(ranges(subject))

            matchMatrix <-
              do.call(rbind,
                      lapply(commonSeqnames, function(seqnm)
                      {
                          qIdxs <- querySplitRanges[[seqnm]]
                          sIdxs <- subjectSplitRanges[[seqnm]]
                          overlaps <-
                            findOverlaps(seqselect(queryRanges, qIdxs),
                                         seqselect(subjectRanges, sIdxs),
                                         maxgap = maxgap, multiple = multiple,
                                         type = type)
                          qHits <- queryHits(overlaps)
                          sHits <- subjectHits(overlaps)
                          matches <-
                            cbind(query = as.integer(qIdxs)[qHits],
                                  subject = as.integer(sIdxs)[sHits])
                          matches[seqselect(queryStrand, qIdxs)[qHits] ==
                                  seqselect(subjectStrand, sIdxs)[sHits], ,
                                  drop=FALSE]
                      }))
            matchMatrix <-
              matchMatrix[IRanges:::orderTwoIntegers(matchMatrix[ , 1L, drop=TRUE],
                                                     matchMatrix[ , 2L, drop=TRUE]), ,
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
        matchMatrix <- .cleanMatchMatrix(matchMatrix)
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
        matchMatrix <- .cleanMatchMatrix(matchMatrix)
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
        matchMatrix <- .cleanMatchMatrix(matchMatrix)
        DIM <- c(length(subject), length(query))
        initialize(ans, matchMatrix = matchMatrix, DIM = DIM)
    }
)
