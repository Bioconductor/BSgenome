### =========================================================================
### The "getSeq" method for BSgenome objects.
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
    lvls <- levels(strand())
    if (!all(strand %in% lvls))
        stop("strand values must be in '", paste(lvls, collapse="' '"), "'")
    if (length(strand) > length)
        stop("too many elements in 'strand'")
    if (length(strand) < length) {
        if (length(strand) == 0L)
            stop("cannot recycle zero-length 'strand'")
        strand <- IRanges:::recycleVector(strand, length)
    }
    strand
}

### Make 'strand' star-free.
### If applied on the output of .normargStrand(), then returns a character
### vector guaranteed to have only "+" or "-" elements (no "*").
### If applied on the output of strand(GRanges), then returns a 'factor' Rle
### guaranteed to have only "+" or "-" elements (no "*").
.starfreeStrand <- function(strand)
{
    ustrand <- unique(strand)
    if (length(ustrand) == 1L && ustrand == "*") {
        ## Strand is "*" for all the elements. Replace it by "+".
        strand[] <- "+"
        return(strand)
    }
    if ("*" %in% ustrand)
        stop("cannot mix \"*\" with other strand values")
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
    strand <- .starfreeStrand(.normargStrand(strand, max0123))
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

### Error messages refer to 'x' as "names" because of the context in which
### .toGRanges() is called.
.toGRanges <- function(x, strand)
{
    if (is(x, "RangedData")) {
        msg <- c("Starting with BioC 2.14, passing a RangedData object ",
                 "to getSeq() is deprecated. Please use a GRanges object ",
                 "instead. Note that you can obtain a GRanges object by ",
                 "simply coercing the RangedData object with ",
                 "'as(..., \"RangedData\")'. However, we strongly recommend ",
                 "that you start migrating your code to operate on GRanges ",
                 "objects instead of RangedData objects. Please ask on the ",
                 "bioc-devel mailing list if you have questions or concerns ",
                 "about this (http://bioconductor.org/help/mailing-list/).")
        .Deprecated(msg=msg)
        if ("strand" %in% colnames(x)) {
            if (!identical(strand, "+"))
                stop("'strand' cannot be specified when 'names' ",
                     "is a RangedData object with a strand column")
        } else {
            strand <- Rle(strand(.normargStrand(strand, nrow(x))))
            values(x) <- DataFrame(strand=strand)
        }
        return(as(x, "GRanges"))
    }
    if (is(x, "RangesList") || is(x, "Ranges")) {
        if (is.null(names(x)))
            stop("when 'names' is a RangesList or Ranges object, ",
                 "it must be named with the sequence names")
        if (is(x, "RangesList"))
            ans <- .newGRangesFromNamedRangesList(x, strand)
        else
            ans <- .newGRangesFromNamedRanges(x, strand)
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

### 'names' must be a character vector or a GRanges object.
### If 'names' is character vector, then 'start', 'end', 'width', and 'strand'
### are assumed to be already normalized (and star-free for 'strand').
### Otherwise, they are ignored and 'strand(names)' is assumed to be star-free.
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
    grl <- split(names, seqnames(names))
    grl_seqlevels <- names(grl)

    ## Loop over the sequence names and extract the ranges.
    dnaset_list <- lapply(seq_len(length(grl)),
        function(i)
        {
            gr <- grl[[i]]
            if (length(gr) == 0L)
                return(DNAStringSet())
            seqlevel <- grl_seqlevels[i]
            is_circular <- isCircular(x)[[seqlevel]]
            idx <- match(seqlevel, x@user_seqnames)
            if (is.na(idx))
                stop("invalid sequence name: ", seqlevel)
            seqname <- names(x@user_seqnames)[[idx]]
            if (is.null(SNPlocs(x, seqname))) {
                ## Try to grab the subject sequence from the cache.
                subject <- try(get(seqname, envir=x@.seqs_cache,
                                   inherits=FALSE), silent=TRUE)
                if (is(subject, "try-error")) {
                    ## The subject sequence was not in the cache.
                    ans <- loadSubseqsFromStrandedSequence(
                                   x@single_sequences,
                                   seqname, ranges(gr), strand(gr),
                                   is_circular)
                    return(ans)
                }
                ## The subject sequence was in the cache.
                .linkToCachedObject(subject) <- .newLinkToCachedObject(
                                                    seqname,
                                                    x@.seqs_cache,
                                                    x@.link_counts)
            } else {
                subject <- x[[seqlevel]]
            }
            masks(subject) <- NULL
            loadSubseqsFromStrandedSequence(
                                   subject,
                                   seqlevel, ranges(gr), strand(gr),
                                   is_circular)
        }
    )
    ## "unsplit" 'dnaset_list'.
    unsplit_list_of_XVectorList("DNAStringSet", dnaset_list,
                                as.factor(seqnames(names)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extractFromBSgenomeMultipleSequences()
###

### Assumes that 'strand' is star-free.
.getOneSeqFromBSgenomeMultipleSequences <- function(x, name,
                                                    start, end, width, strand,
                                                    exact.match)
{
    nhits <- 0L
    for (mseqname in mseqnames(x)) {
        mseq <- x[[mseqname]]
        if (exact.match) {
            ii <- which(names(mseq) == name)
        } else {
            ii <- grep(name, names(mseq))
        }
        if (nhits == 0L && length(ii) == 1L)
            subject <- mseq[[ii]]
        nhits <- nhits + length(ii)
    }
    if (nhits == 0L)
        stop("sequence ", name, " not found")
    if (nhits > 1L)
        stop("sequence ", name, " found more than once, ",
             "please use a non-ambiguous name")
    ans <- subseq(subject, start=start, end=end, width=width)
    if (strand == "+")
        return(xvcopy(ans))
    else
        return(reverseComplement(ans))
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
### are assumed to be already normalized (and star-free for 'strand').
### Otherwise, they are ignored and 'strand(names)' is assumed to be star-free.
.extractFromBSgenomeMultipleSequences <- function(x, names,
                                                  start, end, width, strand)
{
    if (is.character(names)) {
        ans <- lapply(seq_len(length(names)),
            function(i)
                .getOneSeqFromBSgenomeMultipleSequences(x, names[i],
                                start[i], end[i], width[i], strand[i],
                                exact.match=FALSE)
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
                                start, NA, width, strand,
                                exact.match=TRUE)
            }
        )
    }
    .listOfDNAStringToDNAStringSet(ans)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "getSeq" generic and method for BSgenome objects.
###

setMethod("getSeq", "BSgenome",
    function(x, names, start=NA, end=NA, width=NA, strand="+",
             as.character=FALSE)
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
            if (is(names, "GRanges") || is(names, "GRangesList")) {
                if (!identical(strand, "+"))
                    stop("'strand' cannot be specified ",
                         "when 'names' is a GRanges or GRangesList object")
                if (is(names, "GRangesList")) {
                    unlisted_names <- unlist(names, use.names=FALSE)
                    unlisted_ans <- getSeq(x, unlisted_names,
                                           as.character=as.character)
                    return(relist(unlisted_ans, names))
                }
            } else {
                names <- .toGRanges(names, strand)
            }
            ## We don't need the result of merge(). By calling merge() here
            ## we're just making sure that 'x' and 'names' are based on
            ## compatible reference genomes (merge() will raise an error if
            ## they are not).
            merge(seqinfo(x), seqinfo(names))
            strand(names) <- .starfreeStrand(strand(names))
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
        if (is(names, "GRanges"))
            names(ans) <- names(names)
        if (as.character)
            return(as.character(ans))
        if (length(ans) == 1L && is.character(names))
            return(ans[[1L]])  # turn 1st and unique element into DNAString
        ans
    }
)

