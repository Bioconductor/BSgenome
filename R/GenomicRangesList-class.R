### =========================================================================
### GenomicRangesList objects
### -------------------------------------------------------------------------
###
### Class definition

setClass("GenomicRangesList", contains = "CompressedList",
         prototype = prototype(elementType = "GenomicRanges",
                               unlistData = new("GenomicRanges")))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GenomicRangesList <- function(...)
{
    listData <- list(...)
    if (length(listData) == 1 && is.list(listData[[1L]]))
        listData <- listData[[1L]]
    if (!all(sapply(listData, is, "GenomicRanges")))
        stop("all elements in '...' must be GenomicRanges objects")
    IRanges:::newCompressedList("GenomicRangesList", listData)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("seqnames", "GenomicRangesList",
    function(x)
        new2("CompressedRleList",
             unlistData = x@unlistData@seqnames, partitioning = x@partitioning,
             check=FALSE))

setMethod("ranges", "GenomicRangesList",
    function(x, ...)
        new2("CompressedIRangesList",
             unlistData = x@unlistData@ranges, partitioning = x@partitioning,
             check=FALSE))

setMethod("strand", "GenomicRangesList",
    function(x)
        new2("CompressedRleList",
             unlistData = x@unlistData@strand, partitioning = x@partitioning,
             check=FALSE))

setMethod("values", "GenomicRangesList",
    function(x, ...)
        new2("CompressedSplitDataFrameList",
             unlistData = x@unlistData@values, partitioning = x@partitioning,
             check=FALSE))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RangesList methods.
###

setMethod("start", "GenomicRangesList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData = start(x@unlistData@ranges),
             partitioning = x@partitioning, check=FALSE))

setMethod("end", "GenomicRangesList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData = end(x@unlistData@ranges),
             partitioning = x@partitioning, check=FALSE))

setMethod("width", "GenomicRangesList",
    function(x)
        new2("CompressedIntegerList",
             unlistData = width(x@unlistData@ranges),
             partitioning = x@partitioning, check=FALSE))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SplitDataFrameList methods.
###

setMethod("[", "GenomicRangesList",
    function(x, i, j, ..., drop)
    {
        if (!missing(i))
            x <- x[i]
        if (!missing(j))
            values(x) <- values(x)[, j, drop=FALSE]
        x
    }
)

setMethod("ncol", "GenomicRangesList", function(x) ncol(x@unlistData@values))

setMethod("colnames", "GenomicRangesList",
    function(x, do.NULL = TRUE, prefix = "col") 
        colnames(x@unlistData@values, do.NULL = do.NULL, prefix = prefix))
setReplaceMethod("colnames", "GenomicRangesList",
    function(x, value)
    {
        colnames(x@unlistData@values) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

setMethod("show", "GenomicRangesList",
    function(object)
    {
        k <- length(object)
        cumsumN <- cumsum(elementLengths(object))
        N <- tail(cumsumN, 1)
        cat(class(object), " of length ", k, "\n", sep = "")
        if (k == 0L) {
            cat("<0 elements>\n")
        } else if ((k == 1L) || (N <= 20L)) {
            show(as.list(object))
        } else {
            sketch <- function(x) c(head(x, 3), "...", tail(x, 3))
            if (k >= 3 && cumsumN[3L] <= 20)
                showK <- 3
            else if (k >= 2 && cumsumN[2L] <= 20)
                showK <- 2
            else
                showK <- 1
            diffK <- k - showK
            show(as.list(object[seq_len(showK)]))
            if (diffK > 0)
                cat("...\n<", k - showK,
                    ifelse(diffK == 1, " more element>\n", " more elements>\n"),
                    sep="")
        }
    }
)
