### =========================================================================
### OnDiskNamedSequences objects
### -------------------------------------------------------------------------


setClass("OnDiskNamedSequences")  # VIRTUAL class with no slots

### OnDiskNamedSequences API
### ------------------------
### Concrete subclasses need to implement the following:
###   - names()
###   - seqlengthsFilepath()
###   - seqinfo API (implementing seqinfo() is enough to make seqlengths(),
###                  seqlevels(), etc... work)
###   - getListElement() -- load full sequence as XString derivative
### Default methods are provided for the following:
###   - length()
###   - seqnames()
###   - show()
###   - loadSubseqsFromLinearSequence()

setGeneric("seqlengthsFilepath",
    function(x) standardGeneric("seqlengthsFilepath")
)

setGeneric("loadSubseqsFromLinearSequence",
    function(x, seqname, ranges)
        standardGeneric("loadSubseqsFromLinearSequence")
)

### Low-level utilities

.check_getListElement_index <- function(i, what)
{
    if (!is.character(i))
        stop(wmsg(what, " can only be subsetted by name"))
    if (length(i) < 1L)
        stop(wmsg("attempt to select less than one element"))
    if (length(i) > 1L)
        stop(wmsg("attempt to select more than one element"))
}

### Works on an XString derivative or any object 'x' for which seqlengths()
### is defined.
.get_seqlength <- function(x, seqname)
{
    if (!isSingleString(seqname))
        stop(wmsg("'seqname' must be a single string"))
    if (is(x, "XString"))
        return(length(x))
    x_seqlengths <- seqlengths(x)
    idx <- match(seqname, names(x_seqlengths))
    if (is.na(idx))
        stop(wmsg("invalid sequence name: ", seqname))
    x_seqlengths[[idx]]
}

### Default methods

setMethod("length", "OnDiskNamedSequences", function(x) length(names(x)))

setMethod("seqnames", "OnDiskNamedSequences", function(x) seqlevels(x))

setMethod("show", "OnDiskNamedSequences",
    function(object)
    {
        cat(class(object), " instance of length ", length(object),
            ":\n", sep="")
        GenomeInfoDb:::compactPrintNamedAtomicVector(seqlengths(object))
    }
)

### Load regions from a single sequence as an XStringSet derivative.

### 'seqname' is ignored.
setMethod("loadSubseqsFromLinearSequence", "XString",
    function(x, seqname, ranges)
        xvcopy(extractAt(x, ranges))  # ignores 'seqname'
)

setMethod("loadSubseqsFromLinearSequence", "OnDiskNamedSequences",
    function(x, seqname, ranges)
    {
        seq <- getListElement(x, seqname)
        loadSubseqsFromLinearSequence(seq, seqname, ranges)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RdsNamedSequences objects
###
### The 'dirpath' slot should point to a directory that contains one .rds
### file per XString derivative + a seqlengths.rds file that contains a
### serialized named integer vector with the sequence names and lengths.
###

setClass("RdsNamedSequences",
    contains=c("RdsCollection", "OnDiskNamedSequences"),
    prototype=prototype(
        elementType="XString"
    )
)

setMethod("seqlevels", "RdsNamedSequences", function(x) names(x))

setMethod("seqlengthsFilepath", "RdsNamedSequences",
    function(x) file.path(path(x), "seqlengths.rds")
)

.read_seqlengths_from_file <- function(x) readRDS(seqlengthsFilepath(x))

setMethod("seqlengths", "RdsNamedSequences",
    function(x)
    {
        noext_ends <- nchar(x@filenames) - nchar(".rds")
        noext_filenames <- substr(x@filenames, 1L, noext_ends)
        setNames(.read_seqlengths_from_file(x)[noext_filenames], names(x))
    }
)

setMethod("seqinfo", "RdsNamedSequences",
    function(x)
    {
        x_seqlengths <- seqlengths(x)
        Seqinfo(names(x_seqlengths), unname(x_seqlengths))
    }
)

setAs("RdsCollection", "RdsNamedSequences",
    function(from)
    {
        ans <- new2("RdsNamedSequences", from, elementType="XString",
                                         check=FALSE)
        seqlengths_from_file <- .read_seqlengths_from_file(ans)
        if (!is.integer(seqlengths_from_file))
            stop(wmsg("object serialized in ",
                      "file '", seqlengthsFilepath(ans), "' ",
                      "must be a named integer vector"))
        if (!all(names(from) %in% names(seqlengths_from_file)))
            stop(wmsg("the names on the RdsCollection object to coerce ",
                      "to RdsNamedSequences must be a subset of ",
                      "the names on the integer vector serialized ",
                      "in file '", seqlengthsFilepath(ans), "'"))
        ans
    }
)

### Constructor.
RdsNamedSequences <- function(path, seqnames)
{
    filenames <- paste0(seqnames, ".rds")
    as(RdsCollection(path, filenames), "RdsNamedSequences")
}

### Load a full sequence as an XString derivative.
setMethod("getListElement", "RdsNamedSequences",
    function(x, i, exact=TRUE)
    {
        .check_getListElement_index(i, "an RdsNamedSequences object")
        ans <- callNextMethod()
        if (!is(ans, "XString")) {
            filepath <- file.path(path(x), x@filenames[[i]])
            stop(wmsg("serialized object in file '", filepath, "' ",
                      "must be an XString derivative"))
        }
        updateObject(ans)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RdaNamedSequences objects
###
### June 2020: THE RdaNamedSequences CLASS IS SUPERSEDED BY THE
### RdsNamedSequences CLASS!
### TODO: Deprecate the RdaNamedSequences class.
###
### The "dirpath" slot should contain 1 serialized XString derivative
### per sequence + a serialized named integer vector ('seqlengths.rda')
### containing the sequence names and lengths.
###

setClass("RdaNamedSequences",
    contains=c("RdaCollection", "OnDiskNamedSequences"),
    representation(
        seqlengths="RdaCollection"
    )
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
            stop(wmsg("serialized object in file '", seqlengthsFilepath(x),
                      "' must be a named integer vector"))
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

### Load a full sequence as an XString derivative.
setMethod("getListElement", "RdaNamedSequences",
    function(x, i, exact=TRUE)
    {
        .check_getListElement_index(i, "an RdaNamedSequences object")
        ans <- x[[i]]
        if (!is(ans, "XString"))
            stop(wmsg("serialized object in file '", rdaPath(x, i), "' ",
                      "must be an XString derivative"))
        updateObject(ans)
    }
)


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

### Load a full sequence as a DNAString object.
setMethod("getListElement", "FastaNamedSequences",
    function(x, i, exact=TRUE)
    {
        .check_getListElement_index(i, "a FastaNamedSequences object")
        fafile <- x@fafile
        seqlength <- .get_seqlength(fafile, i)
        param <- GRanges(i, IRanges(1L, seqlength))
        scanFa(fafile, param=param)[[1L]]
    }
)

### Load regions from a single sequence as a DNAStringSet object.
### TODO: Try the following optimization: reduce 'ranges' before calling
### scanFa(), then extract the regions from the DNAStringSet returned by
### scanFa(). This way, a given nucleotide is loaded only once.
setMethod("loadSubseqsFromLinearSequence", "FastaNamedSequences",
    function(x, seqname, ranges)
    {
        if (!is(ranges, "IntegerRanges"))
            stop(wmsg("'ranges' must be an IntegerRanges object"))
        fafile <- x@fafile
        seqlength <- .get_seqlength(fafile, seqname)
        if (length(ranges) != 0L &&
            (min(start(ranges)) < 1L || max(end(ranges)) > seqlength))
            stop(wmsg("trying to load regions beyond the boundaries ",
                      "of non-circular sequence \"", seqname, "\""))
        param <- GRanges(seqname, ranges)
        scanFa(fafile, param=param)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TwobitNamedSequences objects
###

setClass("TwobitNamedSequences",
    contains="OnDiskNamedSequences",
    representation(
        ## TwoBitFile object pointing to the .2bit file containing the
        ## sequences.
        twobitfile="TwoBitFile"
    )
)

setMethod("names", "TwobitNamedSequences", function(x) seqlevels(x@twobitfile))

setMethod("seqlengthsFilepath", "TwobitNamedSequences",
    function(x) path(x@twobitfile)
)

setMethod("seqinfo", "TwobitNamedSequences", function(x) seqinfo(x@twobitfile))

### Constructor.
TwobitNamedSequences <- function(filepath)
{
    twobitfile <- TwoBitFile(filepath)
    new("TwobitNamedSequences", twobitfile=twobitfile)
}

### Load a full sequence as a DNAString object.
setMethod("getListElement", "TwobitNamedSequences",
    function(x, i, exact=TRUE)
    {
        .check_getListElement_index(i, "a TwobitNamedSequences object")
        twobitfile <- x@twobitfile
        seqlength <- .get_seqlength(twobitfile, i)
        which <- GRanges(i, IRanges(1L, seqlength))
        import(twobitfile, which=which)[[1L]]
    }
)

### Load regions from a single sequence as a DNAStringSet object.
setMethod("loadSubseqsFromLinearSequence", "TwobitNamedSequences",
    function(x, seqname, ranges)
    {
        if (!is(ranges, "IntegerRanges"))
            stop(wmsg("'ranges' must be an IntegerRanges object"))
        twobitfile <- x@twobitfile
        seqlength <- .get_seqlength(twobitfile, seqname)
        if (length(ranges) != 0L &&
            (min(start(ranges)) < 1L || max(end(ranges)) > seqlength))
            stop(wmsg("trying to load regions beyond the boundaries ",
                      "of non-circular sequence \"", seqname, "\""))
        which <- GRanges(seqname, ranges)
        import(twobitfile, which=which)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### loadSubseqsFromStrandedSequence()
###

.loadSubseqsFromCircularSequence <- function(x, seqname, ranges)
{
    if (!is(ranges, "IntegerRanges"))
        stop(wmsg("'ranges' must be an IntegerRanges object"))
    seqlength <- .get_seqlength(x, seqname)
    LRranges <- splitLRranges(ranges, seqlength, seqname)
    Lans <- loadSubseqsFromLinearSequence(x, seqname, LRranges$L)
    Rans <- loadSubseqsFromLinearSequence(x, seqname, LRranges$R)
    xscat(Lans, Rans)
}

loadSubseqsFromStrandedSequence <- function(x, seqname, ranges, strand,
                                            is_circular=NA)
{
    if (length(ranges) != length(strand))
        stop(wmsg("'ranges' and 'strand' must have the same length"))
    if (!is.logical(is_circular) || length(is_circular) != 1L)
        stop(wmsg("'is_circular' must be a single logical"))
    if (is_circular %in% c(NA, FALSE)) {
        loadFUN <- loadSubseqsFromLinearSequence
    } else {
        loadFUN <- .loadSubseqsFromCircularSequence
    }
    ans <- loadFUN(x, seqname, ranges)
    idx <- which(strand == "-")
    if (length(idx) != 0L)
        ans[idx] <- reverseComplement(ans[idx])
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### OnDiskNamedSequences() constructor
###

OnDiskNamedSequences <- function(dirpath, seqnames=NULL)
{
    filename <- "seqlengths.rds"
    filepath <- file.path(dirpath, filename)
    if (file.exists(filepath))
        return(RdsNamedSequences(dirpath, seqnames))

    filename <- "seqlengths.rda"
    filepath <- file.path(dirpath, filename)
    if (file.exists(filepath))
        return(RdaNamedSequences(dirpath, seqnames))

    filename <- "single_sequences.fa"
    filepath <- file.path(dirpath, filename)
    if (file.exists(filepath))
        return(FastaNamedSequences(filepath))

    filename <- "single_sequences.fa.rz"
    filepath <- file.path(dirpath, filename)
    if (file.exists(filepath))
        return(FastaNamedSequences(filepath))

    filename <- "single_sequences.2bit"
    filepath <- file.path(dirpath, filename)
    if (file.exists(filepath))
        return(TwobitNamedSequences(filepath))

    stop(wmsg("invalid directory content at ", dirpath))
}

