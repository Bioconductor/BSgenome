### =========================================================================
### GenomicFeatureList objects
### -------------------------------------------------------------------------
###
### Class definition

setClass("GenomicFeatureList", contains = "CompressedList",
         prototype = prototype(elementType = "GenomicFeature",
                               unlistData = new("GenomicFeature")))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GenomicFeatureList <- function(...)
{
    listData <- list(...)
    if (length(listData) == 1 && is.list(listData[[1L]]))
        listData <- listData[[1L]]
    if (!all(sapply(listData, is, "GenomicFeature")))
        stop("all elements in '...' must be GenomicFeature objects")
    IRanges:::newCompressedList("GenomicFeatureList", listData)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("RangedData", "GenomicFeatureList",
    function(from)
    {
        ranges <- unlist(ranges(from), use.names=FALSE)
        values <- unlist(values(from), use.names=FALSE)
        nms <- rownames(from)
        rownames(values) <- NULL
        whichStrand <- which(colnames(values) == "strand")
        if (length(whichStrand) > 0)
            values <- values[-whichStrand]
        new("GenomicFeatureList",
            unlistData =
            GenomicFeature(seqnames = space(from),
                           ranges = ranges,
                           strand = Rle(strand(from)),
                           values),
            partitioning =
            PartitioningByEnd(end = seq_len(nrow(from)), names = nms))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("seqnames", "GenomicFeatureList",
    function(x)
        new2("CompressedRleList",
             unlistData = x@unlistData@seqnames, partitioning = x@partitioning,
             check=FALSE))

setMethod("ranges", "GenomicFeatureList",
    function(x, ...)
        new2("CompressedIRangesList",
             unlistData = x@unlistData@ranges, partitioning = x@partitioning,
             check=FALSE))

setMethod("strand", "GenomicFeatureList",
    function(x)
        new2("CompressedRleList",
             unlistData = x@unlistData@strand, partitioning = x@partitioning,
             check=FALSE))

setMethod("values", "GenomicFeatureList",
    function(x, ...)
        new2("CompressedSplitDataFrameList",
             unlistData = x@unlistData@values, partitioning = x@partitioning,
             check=FALSE))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RangesList methods.
###

setMethod("start", "GenomicFeatureList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData = start(x@unlistData@ranges),
             partitioning = x@partitioning, check=FALSE))

setMethod("end", "GenomicFeatureList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData = end(x@unlistData@ranges),
             partitioning = x@partitioning, check=FALSE))

setMethod("width", "GenomicFeatureList",
    function(x)
        new2("CompressedIntegerList",
             unlistData = width(x@unlistData@ranges),
             partitioning = x@partitioning, check=FALSE))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SplitDataFrameList methods.
###

setMethod("[", "GenomicFeatureList",
    function(x, i, j, ..., drop)
    {
        if (!missing(i))
            x <- x[i]
        if (!missing(j))
            values(x) <- values(x)[, j, drop=FALSE]
        x
    }
)

setMethod("ncol", "GenomicFeatureList", function(x) ncol(x@unlistData@values))

setMethod("colnames", "GenomicFeatureList",
    function(x, do.NULL = TRUE, prefix = "col") 
        colnames(x@unlistData@values, do.NULL = do.NULL, prefix = prefix))
setReplaceMethod("colnames", "GenomicFeatureList",
    function(x, value)
    {
        colnames(x@unlistData@values) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

setMethod("show", "GenomicFeatureList",
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
