### =========================================================================
### getSeq()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers called by "getSeq" method for BSgenome objects.
###

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extractFromBSgenomeSingleSequences()
###

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

### 'names' must be a character vector or a GRanges object.
### If 'names' is character vector, then 'start', 'end', 'width', and 'strand'
### are assumed to be already normalized. Otherwise, they are ignored.
.extractFromBSgenomeSingleSequences <- function(x, names,
                                                start, end, width, strand)
{
    if (is.character(names)) {
        refwidths <- seqlengths(x)[names]
        ranges <- solveUserSEW(refwidths, start=start, end=end, width=width)
        names <- GRanges(seqnames=names, ranges=ranges, strand=strand)
    }

    ## Check that 'seqlengths(names)' is compatible with 'x'.
    #if (!all(names(seqlengths(names)) %in% names(seqlengths(x))))
    #    stop("sequence names in GRanges are incompatible with BSgenome object")
    #if (!all(is.na(seqlengths(names)))
    # && !identical(seqlengths(names),
    #               seqlengths(x)[names(seqlengths(names))]))
    #    stop("sequence lengths in GRanges are incompatible ",
    #         "with BSgenome object")

    ## Split 'names' by sequence names.
    grglist <- split(names, seqnames(names))
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
                               as.factor(seqnames(names)))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extractFromBSgenomeMultipleSequences()
###

.getOneSeqFromBSgenomeMultipleSequences <- function(x, name,
                                                    start, end, width, strand)
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

### 'names' must be a character vector or a GRanges object.
### If 'names' is character vector, then 'start', 'end', 'width', and 'strand'
### are assumed to be already normalized. Otherwise, they are ignored.
.extractFromBSgenomeMultipleSequences <- function(x, names,
                                                  start, end, width, strand)
{
    if (is.character(names)) {
        ans <- lapply(seq_len(length(names)),
            function(i)
                .getOneSeqFromBSgenomeMultipleSequences(x, names[i],
                                start[i], end[i], width[i], strand[i])
        )
    } else {
        ans <- lapply(seq_len(length(names)),
            function(i)
            {
                name <- as.character(seqnames(names))[i]
                start <- start(names)[i]
                width <- width(names)[i]
                strand <- as.character(strand(names))[i]
                .getOneSeqFromBSgenomeMultipleSequences(x, name,
                                start, NA, width, strand)
            }
        )
    }
    .listOfDNAStringToDNAStringSet(ans)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "getSeq" generic and method for BSgenome objects.
###

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
            sseq_idx <- args$names %in% seqnames(x)
            sseq_args <- list(names=args$names[sseq_idx],
                              start=args$start[sseq_idx],
                              end=args$end[sseq_idx],
                              width=args$width[sseq_idx],
                              strand=args$strand[sseq_idx])
            mseq_idx <- !sseq_idx
            mseq_args <- list(names=args$names[mseq_idx],
                              start=args$start[mseq_idx],
                              end=args$end[mseq_idx],
                              width=args$width[mseq_idx],
                              strand=args$strand[mseq_idx])
        } else {
            if (!identical(c(start, end, width), c(NA, NA, NA)))
                stop("'start', 'end' and 'width' can only be specified when ",
                     "'names' is either missing, a character vector/factor, ",
                     "a character-Rle, or a factor-Rle")
            names <- .newGRanges(names, strand)
            seqnames <- as.character(seqnames(names))
            sseq_idx <- seqnames %in% seqnames(x)
            sseq_args <- list(names=.dropUnusedGRangesSeqnamesLevels(
                                    names[sseq_idx]))
            mseq_idx <- !sseq_idx
            mseq_args <- list(names=names[mseq_idx])
        }
        ans <- rep.int(DNAStringSet(""), length(sseq_idx))
        ans[sseq_idx] <- .extractFromBSgenomeSingleSequences(x,
                                   sseq_args$names,
                                   sseq_args$start,
                                   sseq_args$end,
                                   sseq_args$width,
                                   sseq_args$strand)
        ans[mseq_idx] <- .extractFromBSgenomeMultipleSequences(x,
                                   mseq_args$names,
                                   mseq_args$start,
                                   mseq_args$end,
                                   mseq_args$width,
                                   mseq_args$strand)
        if (as.character)
            ans <- as.character(ans)
        else if (length(ans) == 1L && is.character(names))
            ans <- ans[[1L]]  # turn 1st and unique element into DNAString
        ans
    }
)

