### =========================================================================
### OnDiskNamedSequences objects
### -------------------------------------------------------------------------


setClass("OnDiskNamedSequences")  # VIRTUAL class with no slots

### OnDiskNamedSequences API:
###   - length()
###   - names()
###   - [[
###   - seqlengthsFilepath()
###   - seqinfo API (seqinfo(), seqlengths(), seqlevels(), etc...)
###   - seqnames()

setGeneric("seqlengthsFilepath",
    function(x) standardGeneric("seqlengthsFilepath")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RdaNamedSequences objects
###
### The "dirpath" slot should contain 1 serialized XString object per
### sequence + a serialized named integer vector ('seqlengths.rda')
### containing the sequence names and lengths.
###

setClass("RdaNamedSequences",
    contains=c("RdaCollection", "OnDiskNamedSequences"),
    representation(
        seqlengths="RdaCollection"
    )
)

setMethod("[[", "RdaNamedSequences",
    function(x, i, j, ...)
    {
        ans <- callNextMethod()
        if (!is(ans, "XString"))
            stop("serialized object in file '", rdaPath(x, i), "' ",
                 "must be an XString object")
        updateObject(ans)
    }
)

setMethod("seqlengthsFilepath", "RdaNamedSequences",
    function(x) rdaPath(x@seqlengths, "seqlengths")
)

setMethod("seqlevels", "RdaNamedSequences", function(x) names(x))

setMethod("seqlengths", "RdaNamedSequences",
    function(x)
    {
        ans <- x@seqlengths[["seqlengths"]]
        if (!is.integer(ans) || is.null(names(ans)))
            stop("serialized object in file '", seqlengthsFilepath(x), "' ",
                 "must be a named integer vector")
        ans
    }
)

setMethod("seqinfo", "RdaNamedSequences",
    function(x)
    {
        x_seqlengths <- seqlengths(x)
        x_seqlevels <- seqlevels(x)
        Seqinfo(x_seqlevels, x_seqlengths)
    }
)

### Constructor.
RdaNamedSequences <- function(dirpath, seqnames)
{
    sequences <- RdaCollection(dirpath, seqnames)
    seqlengths <- RdaCollection(dirpath, "seqlengths")
    new("RdaNamedSequences", sequences, seqlengths=seqlengths)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FastaNamedSequences objects
###

setClass("FastaNamedSequences",
    contains="OnDiskNamedSequences",
    representation(
        ## FaFile object pointing to the .fa/.fa.fai (or .fa.rz/.fa.rz.fai
        ## if RAZip compressed) files containing the sequences/index.
        fafile="FaFile"
    )
)

setMethod("names", "FastaNamedSequences", function(x) seqlevels(x@fafile))

setMethod("length", "FastaNamedSequences", function(x) length(names(x)))

### We only support subetting by name.
### Returns a DNAString object.
setMethod("[[", "FastaNamedSequences",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        if (!is.character(i))
            stop("a FastaNamedSequences object can only be subsetted by name")
        if (length(i) < 1L)
            stop("attempt to select less than one element")
        if (length(i) > 1L)
            stop("attempt to select more than one element")
        fafile <- x@fafile
        fafile_seqlengths <- seqlengths(fafile)
        idx <- match(i, names(fafile_seqlengths))
        if (is.na(idx))
            stop("invalid sequence name: ", i)
        param <- GRanges(i, IRanges(1L, fafile_seqlengths[[idx]]))
        scanFa(fafile, param=param)[[1L]]
    }
)

setMethod("seqlengthsFilepath", "FastaNamedSequences",
    function(x) path(x@fafile)
)

setMethod("seqinfo", "FastaNamedSequences", function(x) seqinfo(x@fafile))

### Constructor.
FastaNamedSequences <- function(filepath)
{
    fafile <- FaFile(filepath)
    open(fafile)
    new("FastaNamedSequences", fafile=fafile)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### OnDiskNamedSequences methods
###

setMethod("seqnames", "OnDiskNamedSequences", function(x) seqlevels(x))

