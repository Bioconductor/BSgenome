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
    Rle(strand(strand))
}

.mkGRangesFromRangesList <- function(x, strand)
{
    seqnames <- rep.int(names(x), elementLengths(x))
    ranges <- unlist(x, use.names=FALSE)
    strand <- .normargStrand(strand, length(ranges))
    GRanges(seqnames=seqnames, ranges=ranges, strand=strand)
}

.mkGRanges <- function(names, start, end, width, strand)
{
    if (is(names, "GRanges")) {
        if (!identical(c(start, end, width, strand), c(NA, NA, NA, "+")))
            stop("'start', 'end', 'width' and 'strand' cannot be ",
                 "specified when 'names' is a GRanges object")
        return(names)
    }
    if (is(names, "RangedData")) {
        if (!identical(c(start, end, width), c(NA, NA, NA)))
            stop("'start', 'end' and 'width' cannot be ",
                 "specified when 'names' is a RangedData object")
        if ("strand" %in% colnames(names)) {
            if (!identical(strand, "+"))
                stop("'strand' cannot be specified ",
                     "when 'names' is a stranded RangedData object")
        } else {
            strand <- .normargStrand(strand, length(names))
            values(names) <- DataFrame(strand=strand)
        }
        return(as(names, "GRanges"))
    }
    if (is(names, "RangesList")) {
        if (!identical(c(start, end, width), c(NA, NA, NA)))
            stop("'start', 'end' and 'width' cannot be ",
                 "specified when 'names' is a RangesList object")
        if (is.null(names(names)))
            stop("when 'names' is a RangesList object, it must be ",
                 "named with the sequence names")
        return(.mkGRangesFromRangesList(names, strand))
    } 
    if (is.factor(names))
        names <- as.vector(names)
    if (!is.character(names) || any(is.na(names)))
        stop("'names' must be a character vector (with no NAs)")
    ranges <- IRanges(start=start, end=end, width=width)
    strand <- .normargStrand(strand, length(ranges))
    GRanges(seqnames=names, ranges=ranges, strand)
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
    if (!all(names(seqlengths(granges)) %in% names(seqlengths(x))))
        stop("sequence names in GRanges are incompatible with BSgenome object")
    if (!all(is.na(seqlengths(granges)))
     && !identical(seqlengths(granges), seqlengths(x)[names(seqlengths(granges))]))
        stop("sequence lengths in GRanges are incompatible with BSgenome object")
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

.extractSpecialSeq <- function(x, name, start, width, strand)
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

.extractSpecialSeqsFromBSgenome <- function(x, granges)
{
    ans <- lapply(seq_len(length(granges)),
                  function(i)
                  {
                      name <- as.character(seqnames(granges))[i]
                      start <- start(granges)[i]
                      width <- width(granges)[i]
                      strand <- as.character(strand(granges))[i]
                      .extractSpecialSeq(x, name, start, width, strand)
                  }
           )
    ## Inefficient way to turn a list of DNAString objects into a DNAStringSet
    ## object. TODO: Find a better way (then maybe make this a "DNAStringSet"
    ## method for list objects).
    subject <- do.call(xscat, ans)
    DNAStringSet(successiveViews(subject, elementLengths(ans)))
}

setGeneric("getSeq", function(x, ...) standardGeneric("getSeq"))

setMethod("getSeq", "BSgenome",
    function(x, names, start=NA, end=NA, width=NA, strand="+",
             as.character=TRUE)
    {
        if (!isTRUEorFALSE(as.character))
            stop("'as.character' must be TRUE or FALSE")
        if (missing(names))
            names <- seqnames(x)
        granges <- .mkGRanges(names, start, end, width, strand)
        seqnames <- as.character(seqnames(granges))
        is_special <- !(seqnames %in% seqnames(x))
        if (any(is_special)) {
            ans <- rep.int(DNAStringSet(""), length(is_special))
            ans[is_special] <-
                .extractSpecialSeqsFromBSgenome(x, granges[is_special])
            granges <- .dropUnusedGRangesSeqnamesLevels(granges[!is_special])
            ans[!is_special] <- .extractSeqsFromBSgenome(x, granges)
        } else {
            ans <- .extractSeqsFromBSgenome(x, granges)
        }
        if (as.character)
            ans <- as.character(ans)
        ans
    }
)

