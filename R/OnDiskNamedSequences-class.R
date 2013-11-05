### =========================================================================
### OnDiskNamedSequences objects
### -------------------------------------------------------------------------


setClass("OnDiskNamedSequences")  # VIRTUAL class with no slots


setGeneric("loadOnDiskNamedSequence", signature="x",
    function(x, seqname) standardGeneric("loadOnDiskNamedSequence")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RdaSequences objects
###

setClass("RdaSequences",
    contains="OnDiskNamedSequences",
    representation(
        ## 'dirpath' should contain 1 serialized XString object per sequence
        ## + a serialized named integer vector ('seqlengths.rda') containing
        ## the sequence names and lengths.
        dirpath="character"
    )
)

.get_rda_filepath <- function(dirpath, objname)
{
    filename <- paste0(objname, ".rda")
    file.path(dirpath, filename)
}

.load_serialized_object <- function(dirpath, objname)
{
    filepath <- .get_rda_filepath(dirpath, objname)
    tempenv <- new.env(parent=emptyenv())
    loaded_names <- load(filepath, envir=tempenv)
    if (length(loaded_names) != 1L)
        stop("file '", filepath, "' contains 0 or more ",
             "than 1 serialized object")
    if (loaded_names != objname)
        stop("serialized object in file '", filepath, "' ", 
             "doesn't have the expected name (expected: ", objname,
             " -- current: ", loaded_names, ")")
    get(objname, envir=tempenv)
}

setMethod("seqlengths", "RdaSequences",
    function(x)
    {
        ans <- .load_serialized_object(x@dirpath, "seqlengths")
        if (!is.integer(ans) || is.null(names(ans)))
            stop("serialized object in file '", 
                 .get_rda_filepath(x@dirpath, "seqlengths"), "' ",
                 "must be a named integer vector")
        ans
    }
)

setMethod("seqinfo", "RdaSequences",
    function(x)
    {
        x_seqlengths <- seqlengths(x)
        x_seqnames <- names(x_seqlengths)
        Seqinfo(x_seqnames, x_seqlengths)
    }
)

### 'seqname' is assumed to be a single character string (NOT checked).
### Returns an XString object.
setMethod("loadOnDiskNamedSequence", "RdaSequences",
    function(x, seqname)
    {
        ans <- .load_serialized_object(x@dirpath, seqname)
        if (!is(ans, "XString"))
            stop("serialized object in file '",
                 .get_rda_filepath(x@dirpath, seqname), "' ",
                 "must be an XString object")
        updateObject(ans)
    }
)

### Constructor.
RdaSequences <- function(dirpath)
{
    new("RdaSequences", dirpath=dirpath)
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

setMethod("seqinfo", "FaRzSequences", function(x) seqinfo(x@fafile))

### 'seqname' is assumed to be a single character string (NOT checked).
### Returns a DNAString object.
setMethod("loadOnDiskNamedSequence", "FaRzSequences",
    function(x, seqname)
    {
        fafile <- x@fafile
        fafile_seqlengths <- seqlengths(fafile)
        i <- match(seqname, names(fafile_seqlengths))
        if (is.na(i))
            stop("invalid sequence name: ", seqname)
        param <- GRanges(seqname, IRanges(1L, fafile_seqlengths[[i]]))
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

### We only support subetting by name.
setMethod("[[", "OnDiskNamedSequences",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        if (!is.character(i))
            stop("an OnDiskNamedSequences can only be subsetted by name")
        if (length(i) < 1L)
            stop("attempt to select less than one element")
        if (length(i) > 1L)
            stop("attempt to select more than one element")
        loadOnDiskNamedSequence(x, i)
    }
)

