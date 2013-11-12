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
### RdaSequences objects
###
### The "dirpath" slot should contain 1 serialized XString object per
### sequence + a serialized named integer vector ('seqlengths.rda')
### containing the sequence names and lengths.
###

setClass("RdaSequences",
    contains=c("RdaCollection", "OnDiskNamedSequences"),
    representation(
        seqlengths="RdaCollection"
    )
)

setMethod("[[", "RdaSequences",
    function(x, i, j, ...)
    {
        ans <- callNextMethod()
        if (!is(ans, "XString"))
            stop("serialized object in file '", rdaPath(x, i), "' ",
                 "must be an XString object")
        updateObject(ans)
    }
)

setMethod("seqlengthsFilepath", "RdaSequences",
    function(x) rdaPath(x@seqlengths, "seqlengths")
)

setMethod("seqlevels", "RdaSequences", function(x) names(x))

setMethod("seqlengths", "RdaSequences",
    function(x)
    {
        ans <- x@seqlengths[["seqlengths"]]
        if (!is.integer(ans) || is.null(names(ans)))
            stop("serialized object in file '", seqlengthsFilepath(x), "' ",
                 "must be a named integer vector")
        ans
    }
)

setMethod("seqinfo", "RdaSequences",
    function(x)
    {
        x_seqlengths <- seqlengths(x)
        x_seqlevels <- seqlevels(x)
        Seqinfo(x_seqlevels, x_seqlengths)
    }
)

### Constructor.
RdaSequences <- function(dirpath, seqnames)
{
    sequences <- RdaCollection(dirpath, seqnames)
    seqlengths <- RdaCollection(dirpath, "seqlengths")
    new("RdaSequences", sequences, seqlengths=seqlengths)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FaRzSequences objects
###

setClass("FaRzSequences",
    contains="OnDiskNamedSequences",
    representation(
        ## FaFile object pointing to the .fa/.fa.fai (or .fa.rz/.fa.rz.fai
        ## if RAZip compressed) files containing the sequences/index.
        fafile="FaFile"
    )
)

### We only support subetting by name.
### Returns a DNAString object.
setMethod("[[", "FaRzSequences",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        if (!is.character(i))
            stop("a FaRzSequences object can only be subsetted by name")
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

setMethod("seqlengthsFilepath", "FaRzSequences",
    function(x) path(x@fafile)
)

setMethod("seqinfo", "FaRzSequences", function(x) seqinfo(x@fafile))

setMethod("names", "FaRzSequences", function(x) seqlevels(x@fafile))

### Constructor.
FaRzSequences <- function(filepath)
{
    fafile <- FaFile(filepath)
    open(fafile)
    new("FaRzSequences", fafile=fafile)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### OnDiskNamedSequences methods
###

setMethod("seqnames", "OnDiskNamedSequences", function(x) seqlevels(x))

