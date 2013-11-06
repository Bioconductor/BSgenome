### =========================================================================
### OnDiskNamedSequences objects
### -------------------------------------------------------------------------


setClass("OnDiskNamedSequences")  # VIRTUAL class with no slots


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

setClass("RdaSequences", contains=c("RdaCollection", "OnDiskNamedSequences"))

setMethod("seqlengthsFilepath", "RdaSequences",
    function(x) rdaPath(x, "seqlengths")
)

setMethod("seqinfo", "RdaSequences",
    function(x)
    {
        x_seqlengths <- seqlengths(x)
        x_seqnames <- names(x_seqlengths)
        Seqinfo(x_seqnames, x_seqlengths)
    }
)

setMethod("[[", "RdaSequences",
    function(x, i, j, ...)
    {
        ans <- callNextMethod()
        if (i == "seqlengths") {
            if (!is.integer(ans) || is.null(names(ans)))
                stop("serialized object in file '", rdaPath(x, i), "' ",
                     "must be a named integer vector")
            return(ans)
        }
        if (!is(ans, "XString"))
            stop("serialized object in file '", rdaPath(x, i), "' ",
                 "must be an XString object")
        updateObject(ans)
    }
)

setMethod("seqlengths", "RdaSequences", function(x) x[["seqlengths"]])

### Constructor.
RdaSequences <- function(dirpath, seqnames)
{
    rdas <- RdaCollection(dirpath, c("seqlengths", seqnames))
    new("RdaSequences", rdas)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FaRzSequences objects
###

setClass("FaRzSequences",
    contains="OnDiskNamedSequences",
    representation(
        ## FaFile object pointing to the .fa.rz file (and .fa.rz.fai index)
        ## containing the sequences.
        fafile="FaFile"
    )
)

setMethod("seqlengthsFilepath", "FaRzSequences",
    function(x) path(x@fafile)
)

setMethod("seqinfo", "FaRzSequences", function(x) seqinfo(x@fafile))

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
        param <- GRanges(seqname, IRanges(1L, fafile_seqlengths[[idx]]))
        scanFa(fafile, param=param)[[1L]]
    }
)

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

