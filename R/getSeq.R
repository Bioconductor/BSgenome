### =========================================================================
### getSeq()
### -------------------------------------------------------------------------

.normargStrand <- function(strand, length)
{
    if (is(strand, "Rle"))
        strand <- as.vector(strand)
    if (is.factor(strand))
        strand <- as.vector(strand)
    if (!is.character(strand))
        stop("invalid 'strand'")
    if (length(strand) > length)
        stop("too many elements in 'strand'")
    if (length(strand) < length) {
        if (length(strand) == 0L)
            stop("cannot recycle zero-length 'strand'")
        strand <- IRanges:::recycleVector(strand, length)
    }
    strand
}

### Assumes 'names' is a character vector.
.normGetSeqArgs <- function(names, start, end, width, strand)
{
    start <- IRanges:::.normargSEW(start, "start")
    end <- IRanges:::.normargSEW(end, "end")
    width <- IRanges:::.normargSEW(width, "width")
    l0 <- length(names)
    l1 <- length(start)
    l2 <- length(end)
    l3 <- length(width)
    max0123 <- max(l0, l1, l2, l3)
    ## Recycling will fail for vectors of length 0.
    names <- IRanges:::recycleVector(names, max0123)
    start <- IRanges:::recycleVector(start, max0123)
    end <- IRanges:::recycleVector(end, max0123)
    width <- IRanges:::recycleVector(width, max0123)
    strand <- .normargStrand(strand, max0123)
    list(names=names, start=start, end=end, width=width, strand=strand)
}

### Assumes 'seqlengths' is a named integer vector and 'names' is a character
### vector with values in 'names(seqlengths)'.
.newGRangesFromGetSeqArgs <- function(seqlengths, names,
                                      start, end, width, strand)
{
    refwidths <- seqlengths[names]
    ranges <- solveUserSEW(refwidths, start=start, end=end, width=width)
    GRanges(seqnames=names, ranges=ranges, strand=strand)
}

### Assumes 'x' is a RangesList object with names.
.newGRangesFromNamedRangesList <- function(x, strand)
{
    seqnames <- rep.int(names(x), elementLengths(x))
    ranges <- unlist(x, use.names=FALSE)
    strand <- .normargStrand(strand, length(ranges))
    GRanges(seqnames=seqnames, ranges=ranges, strand=strand)
}

### Assumes 'x' is a Ranges object with names.
.newGRangesFromNamedRanges <- function(x, strand)
{
    strand <- .normargStrand(strand, length(x))
    GRanges(seqnames=names(x), ranges=unname(x), strand=strand)
}

.newGRanges <- function(names, strand)
{
    if (is(names, "GRanges")) {
        if (!identical(strand, "+"))
            stop("'strand' cannot be specified ",
                 "when 'names' is a GRanges object")
        return(names)
    }
    if (is(names, "RangedData")) {
        if ("strand" %in% colnames(names)) {
            if (!identical(strand, "+"))
                stop("'strand' cannot be specified when 'names' ",
                     "is a RangedData object with a strand column")
        } else {
            strand <- Rle(strand(.normargStrand(strand, length(names))))
            values(names) <- DataFrame(strand=strand)
        }
        return(as(names, "GRanges"))
    }
    if (is(names, "RangesList") || is(names, "Ranges")) {
        if (is.null(names(names)))
            stop("when 'names' is a RangesList or Ranges object, ",
                 "it must be named with the sequence names")
        if (is(names, "RangesList"))
            ans <- .newGRangesFromNamedRangesList(names, strand)
        else
            ans <- .newGRangesFromNamedRanges(names, strand)
        return(ans)
    } 
    stop("invalid 'names'")
}

.dropUnusedGRangesSeqnamesLevels <- function(x)
{
    ## There should be a better way to drop unused levels!
    seqnames <- as.character(seqnames(x))
    GRanges(seqnames=seqnames,
            ranges=ranges(x),
            strand=strand(x),
            seqlengths=seqlengths(x)[unique(seqnames)])
}

### Assumes 'ranges' and 'strand' have the same length and that this
### length is >= 1.
.extractSeqsFromDNAString <- function(subject, ranges, strand)
{
    rglist <- split(unname(ranges), strand)
    if (elementLengths(rglist)[["*"]] != 0L)
        stop("'strand' elements must be \"+\" or \"-\"")
    plus_ranges <- rglist[["+"]]
    if (length(plus_ranges) == 0L) {
        plus_dnaset <- DNAStringSet()
    } else {
        plus_dnaset <- DNAStringSet(subject,
                                    start=start(plus_ranges),
                                    width=width(plus_ranges))
        plus_dnaset <- xvcopy(plus_dnaset)
    }
    minus_ranges <- rglist[["-"]]
    if (length(minus_ranges) == 0L) {
        minus_dnaset <- DNAStringSet()
    } else {
        minus_dnaset <- DNAStringSet(subject,
                                     start=start(minus_ranges),
                                     width=width(minus_ranges))
        minus_dnaset <- reverseComplement(minus_dnaset)
    }
    unsplit.list.of.XStringSet("DNAStringSet",
                               list(plus_dnaset, minus_dnaset),
                               as.factor(strand))
}

.extractSeqsFromBSgenome <- function(x, granges)
{
    ## Check that 'seqlengths(granges)' is compatible with 'x'.
    #if (!all(names(seqlengths(granges)) %in% names(seqlengths(x))))
    #    stop("sequence names in GRanges are incompatible with BSgenome object")
    #if (!all(is.na(seqlengths(granges)))
    # && !identical(seqlengths(granges),
    #               seqlengths(x)[names(seqlengths(granges))]))
    #    stop("sequence lengths in GRanges are incompatible ",
    #         "with BSgenome object")

    ## Split 'granges' by sequence names.
    grglist <- split(granges, seqnames(granges))
    ## Loop over the sequence names and extract the ranges.
    dnaset_list <- lapply(seq_len(length(grglist)),
        function(i)
        {
            grg <- grglist[[i]]
            if (length(grg) == 0L)
                return(DNAStringSet())
            subject <- x[[names(grglist)[i]]]
            masks(subject) <- NULL
            .extractSeqsFromDNAString(subject, ranges(grg), strand(grg))
        }
    )
    ## "unsplit" 'dnaset_list'.
    unsplit.list.of.XStringSet("DNAStringSet", dnaset_list,
                               as.factor(seqnames(granges)))
}

.extractSpecialSeq <- function(x, name, start, end, width, strand)
{
    nhits <- 0L
    for (mseqname in mseqnames(x)) {
        mseq <- x[[mseqname]]
        ii <- grep(name, names(mseq))
        nhits <- nhits + length(ii)
        if (length(ii) == 1L)
            subject <- mseq[[ii]]
    }
    if (nhits == 0L)
        stop("sequence ", name, " not found")
    if (nhits > 1L)
        stop("sequence ", name, " found more than once, ",
             "please use a non-ambiguous name")
    ans <- subseq(subject, start=start, width=width)
    if (strand == "+")
        return(xvcopy(ans))
    if (strand == "-")
        return(reverseComplement(ans))
    stop("'strand' elements must be \"+\" or \"-\"")
}

### Assumes 'x' is a list of DNAString objects and turns it into a
### DNAStringSet object.
### TODO: Find a better way to do this (current implementation is not very
### efficient). Also maybe make this a "DNAStringSet" method for list objects.
.listOfDNAStringToDNAStringSet <- function(x)
{
    if (length(x) == 0L)
        return(DNAStringSet())
    subject <- do.call(xscat, x)
    DNAStringSet(successiveViews(subject, elementLengths(x)))
}

.extractSpecialSeqsFromBSgenome <- function(x, granges)
{
    ans <- lapply(seq_len(length(granges)),
                  function(i)
                  {
                      name <- as.character(seqnames(granges))[i]
                      start <- start(granges)[i]
                      width <- width(granges)[i]
                      strand <- as.character(strand(granges))[i]
                      .extractSpecialSeq(x, name, start, NA, width, strand)
                  }
           )
    .listOfDNAStringToDNAStringSet(ans)
}

.extractSpecialSeqsFromBSgenome2 <- function(x, names,
                                             start, end, width, strand)
{
    ans <- lapply(seq_len(length(names)),
               function(i)
                   .extractSpecialSeq(x, names[i],
                                      start[i], end[i], width[i], strand[i])
           )
    .listOfDNAStringToDNAStringSet(ans)
}

setGeneric("getSeq", function(x, ...) standardGeneric("getSeq"))

setMethod("getSeq", "BSgenome",
    function(x, names, start=NA, end=NA, width=NA, strand="+",
             as.character=TRUE)
    {
        if (!isTRUEorFALSE(as.character))
            stop("'as.character' must be TRUE or FALSE")
        if (missing(names)) {
            names <- seqnames(x)
        } else {
            if (is(names, "Rle"))
                names <- as.vector(names)
            if (is.factor(names))
                names <- as.vector(names)
        }
        if (is.character(names)) {
            args <- .normGetSeqArgs(names, start, end, width, strand)
            is_not_special <- args$names %in% seqnames(x)
            is_special <- !is_not_special
            ans <- rep.int(DNAStringSet(""), length(args$names))
            ans[is_special] <- .extractSpecialSeqsFromBSgenome2(x,
                                   args$names[is_special],
                                   args$start[is_special],
                                   args$end[is_special],
                                   args$width[is_special],
                                   args$strand[is_special])
            granges <- .newGRangesFromGetSeqArgs(seqlengths(x),
                           args$names[is_not_special],
                           args$start[is_not_special],
                           args$end[is_not_special],
                           args$width[is_not_special],
                           args$strand[is_not_special])
            ans[is_not_special] <- .extractSeqsFromBSgenome(x, granges)
            if (length(ans) == 1L)
                ans <- ans[[1L]]  # turn 1st and unique element into DNAString
        } else {
            if (!identical(c(start, end, width), c(NA, NA, NA)))
                stop("'start', 'end' and 'width' can only be specified when ",
                     "'names' is either missing, a character vector/factor, ",
                     "a character-Rle, or a factor-Rle")
            granges <- .newGRanges(names, strand)
            seqnames <- as.character(seqnames(granges))
            is_not_special <- seqnames %in% seqnames(x)
            is_special <- !is_not_special
            if (any(is_special)) {
                ans <- rep.int(DNAStringSet(""), length(is_special))
                ans[is_special] <-
                    .extractSpecialSeqsFromBSgenome(x, granges[is_special])
                granges <- .dropUnusedGRangesSeqnamesLevels(
                               granges[is_not_special])
                ans[is_not_special] <- .extractSeqsFromBSgenome(x, granges)
            } else {
                ans <- .extractSeqsFromBSgenome(x, granges)
            }
        }
        if (as.character)
            ans <- as.character(ans)
        ans
    }
)

