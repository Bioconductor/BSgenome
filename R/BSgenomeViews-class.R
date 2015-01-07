### =========================================================================
### BSgenomeViews objects
### -------------------------------------------------------------------------
###
### The BSgenomeViews class is a container for storing a set of genomic
### positions on a BSgenome object, called the "subject".
###


### We cannot (and should not try to) extend Views here, for the same
### reasons that we didn't make GRanges a subclass of IRanges.
###
### TODO: A cleaner class design would be to have 2 abstractions: IViews and
### GViews. For IViews: the subject is a Vector and the ranges slot is an
### IRanges object. Note that this is how the current Views class is defined. 
### For GViews: the subject is a named List and the granges slot is a
### GRanges object. Both IViews and GViews would be direct subclasses of a
### more general Views class that contains List and has a subject slot.
### BSgenomeViews below then should become a subclass of GViews.

setClass("BSgenomeViews",
    contains="List",
    representation(
        subject="BSgenome",
        granges="GRanges",
        elementMetadata="DataFrame"
    ),
    prototype(
        elementType="DNAString"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("subject", "BSgenomeViews", function(x) x@subject)

setMethod("granges", "BSgenomeViews",
    function(x, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        ans <- x@granges
        if (use.mcols)
            mcols(ans) <- mcols(x)
        ans
    }
)

setMethod("length", "BSgenomeViews", function(x) length(granges(x)))
setMethod("names", "BSgenomeViews", function(x) names(granges(x)))
setMethod("seqnames", "BSgenomeViews", function(x) seqnames(granges(x)))
setMethod("start", "BSgenomeViews", function(x) start(granges(x)))
setMethod("end", "BSgenomeViews", function(x) end(granges(x)))
setMethod("width", "BSgenomeViews", function(x) width(granges(x)))
setMethod("ranges", "BSgenomeViews", function(x, use.mcols=FALSE) ranges(granges(x)))
setMethod("strand", "BSgenomeViews", function(x) strand(granges(x)))
setMethod("seqinfo", "BSgenomeViews", function(x) seqinfo(granges(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

BSgenomeViews <- function(genome, granges)
{
    genome <- getBSgenome(genome)
    if (!is(granges, "GenomicRanges"))
        stop("'granges' must be a GRanges object")
    ans_seqinfo <- seqinfo(genome)
    ## Calling merge() is the standard way to check that 'genome' and
    ## 'granges' are based on the same reference genome.
    merge(ans_seqinfo, seqinfo(granges))
    ans_granges <- granges(granges)
    seqlevels(ans_granges) <- seqlevels(ans_seqinfo)
    seqinfo(ans_granges) <- ans_seqinfo
    ans_mcols <- mcols(granges)
    if (is.null(ans_mcols))
        ans_mcols <- new("DataFrame", nrows=length(granges))
    new("BSgenomeViews", subject=genome, granges=ans_granges,
                         elementMetadata=ans_mcols)
}

### Provided for convenience. Need to do some ugly tweaks with the supplied
### args because of the weird signature of the Views() generic.
setMethod("Views", "BSgenome",
    function(subject, start=NULL, end=NULL, width=NULL, names=NULL)
    {
        if (!(is.null(end) && is.null(width) && is.null(names)))
            stop(wmsg("use call of the form 'Views(genome, gr)', ",
                      "where 'genome' is a BSgenome object ",
                      "and 'gr' a GRanges object, to create a ",
                      "BSgenomeViews object"))
        if (!is(start, "GenomicRanges"))
            stop("the location of the views on the genome must be ",
                 "specified as a GRanges object")
        BSgenomeViews(subject, start)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Displaying
###

.makeNakedMatFromBSgenomeViews <- function(x)
{
    lx <- length(x)
    nc <- ncol(mcols(x))
    ans_seqnames <- as.character(seqnames(x))
    ans_ranges <- showAsCell(ranges(x))
    ans_strand <- as.character(strand(x))
    ans_dna <- sapply(getSeq(subject(x), granges(x)),
                      Biostrings:::toSeqSnippet, 23L)
    if (lx != 0L)
        ans_dna <- paste0("[", ans_dna, "]")
    ans <- cbind(seqnames=as.character(seqnames(x)),
                 ranges=showAsCell(ranges(x)),
                 strand=as.character(strand(x)),
                 dna=ans_dna)
    if (nc > 0L) {
        tmp <- do.call(data.frame, lapply(mcols(x), showAsCell))
        ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
    }
    ans
}

showBSgenomeViews <- function(x, margin="",
                                 print.classinfo=FALSE,
                                 print.seqinfo=FALSE)
{
    lx <- length(x)
    nc <- ncol(mcols(x))
    cat(class(x), " object with ",
        lx, " view", ifelse(lx == 1L, "", "s"),
        " and ",
        nc, " metadata column", ifelse(nc == 1L, "", "s"),
        ":\n", sep="")
    out <- S4Vectors:::makePrettyMatrixForCompactPrinting(x,
               .makeNakedMatFromBSgenomeViews)
    if (print.classinfo) {
        .COL2CLASS <- c(
            seqnames="Rle",
            ranges="IRanges",
            strand="Rle",
            dna="DNAStringSet"
        )
        classinfo <-
            S4Vectors:::makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        ## A sanity check, but this should never happen!
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste0(margin, rownames(out))
    ## We set 'max' to 'length(out)' to avoid the getOption("max.print")
    ## limit that would typically be reached when 'showHeadLines' global
    ## option is set to Inf.
    print(out, quote=FALSE, right=TRUE, max=length(out))
    if (print.seqinfo) {
        cat(margin, "-------\n", sep="")
        cat(margin, "seqinfo: ", summary(seqinfo(x)), "\n", sep="")
    }
}

setMethod("show", "BSgenomeViews",
    function(object)
        showBSgenomeViews(object, margin="  ",
                          print.classinfo=TRUE, print.seqinfo=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

.extract_dna_from_BSgenomeViews <- function(x, use.mcols=FALSE)
{
    ## Doesn't work because getSeq() doesn't propagate the metadata cols of
    ## the GRanges object. Maybe it should?
    #getSeq(subject(x), granges(x, use.mcols=use.mcols))
    if (!isTRUEorFALSE(use.mcols))
        stop("'use.mcols' must be TRUE or FALSE")
    ans <- getSeq(subject(x), granges(x))
    if (use.mcols)
        mcols(ans) <- mcols(x)
    ans
}

.from_BSgenomeViews_to_DNAStringSet <- function(from)
    .extract_dna_from_BSgenomeViews(from, use.mcols=TRUE)

setAs("BSgenomeViews", "DNAStringSet", .from_BSgenomeViews_to_DNAStringSet)
setAs("BSgenomeViews", "XStringSet", .from_BSgenomeViews_to_DNAStringSet)

setMethod("as.character", "BSgenomeViews",
    function(x, ...)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric(x, ...)
    }
)

setMethod("as.data.frame", "BSgenomeViews",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        df1 <- as.data.frame(granges(x, use.mcols=TRUE),
                             row.names=NULL, optional=optional)
        df2 <- as.data.frame(.extract_dna_from_BSgenomeViews(x),
                             row.names=NULL, optional=optional)
        colnames(df2) <- "dna"
        cbind(df1, df2)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("extractROWS", "BSgenomeViews",
    function(x, i)
    {
        x@granges <- extractROWS(x@granges, i)
        x@elementMetadata <- extractROWS(x@elementMetadata, i)
        x
    }
)

### Extracting a view.
setMethod("getListElement", "BSgenomeViews",
    function(x, i, exact=TRUE)
    {
        i <- normalizeDoubleBracketSubscript(i, x, exact=exact,
                                             error.if.nomatch=TRUE)
        getSeq(subject(x), granges(x)[i])[[1L]]
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DNAStringSet methods
###

setMethod("seqtype", "BSgenomeViews",
    function(x) seqtype(new(x@elementType))
)

setMethod("nchar", "BSgenomeViews",
    function(x, type="chars", allowNA=FALSE) width(x)
)

setMethod("unlist", "BSgenomeViews",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

setMethod("alphabetFrequency", "BSgenomeViews",
    function(x, as.prob=FALSE, ...)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

setMethod("hasOnlyBaseLetters", "BSgenomeViews",
    function(x)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

setMethod("uniqueLetters", "BSgenomeViews",
    function(x)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

setMethod("letterFrequencyInSlidingView", "BSgenomeViews",
    function(x, view.width, letters, OR="|", as.prob=FALSE)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

setMethod("letterFrequency", "BSgenomeViews",
    function(x, letters, OR="|", as.prob=FALSE, ...)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

setMethod("oligonucleotideFrequency", "BSgenomeViews",
    function(x, width, step=1,
             as.prob=FALSE, as.array=FALSE,
             fast.moving.side="right", with.labels=TRUE, ...)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

setMethod("nucleotideFrequencyAt", "BSgenomeViews",
    function(x, at, as.prob=FALSE, as.array=TRUE,
             fast.moving.side="right", with.labels=TRUE, ...)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

setMethod("consensusMatrix", "BSgenomeViews",
    function(x, as.prob=FALSE, shift=0L, width=NULL, ...)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

setMethod("consensusString", "BSgenomeViews",
    function(x, ...)
    {
        x <- .extract_dna_from_BSgenomeViews(x)
        callGeneric()
    }
)

### TODO: Add more DNAStringSet methods...

