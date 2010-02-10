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

setMethod("as.data.frame", "GenomicFeatureList",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (is.null(names(x)))
            feature <- rep(seq_len(length(x)), elementLengths(x))
        else
            feature <- rep(names(x), elementLengths(x))
        data.frame(feature = feature,
                   as.data.frame(unlist(x, use.names = FALSE),
                                 row.names = row.names),
                   stringsAsFactors = FALSE)
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

setReplaceMethod("seqnames", "GenomicFeatureList",
    function(x, value) 
    {
        if (!is(value, "AtomicList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an AtomicList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        if (!is(value, "Rle"))
            value <- Rle(as.character(value))
        else if (!is.character(runValue(value)))
            runValue(value) <- as.character(runValue)
        x@unlistData@seqnames <- value
        x
    }
)

setReplaceMethod("ranges", "GenomicFeatureList",
    function(x, value) 
    {
        if (!is(value, "RangesList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not a RangesList with the same ",
                 "elementLengths as 'x'")
        x@unlistData@ranges <- as(unlist(value, use.names = FALSE), "IRanges")
        x
    }
)

setReplaceMethod("strand", "GenomicFeatureList",
    function(x, value) 
    {
        if (!is(value, "AtomicList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an AtomicList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        if (!is(value, "Rle"))
            value <- Rle(strand(value))
        else if (!is.factor(runValue(value)) ||
                 !identical(levels(runValue(value)), levels(strand())))
            runValue(value) <- strand(runValue)
        x@unlistData@strand <- value
        x
    }
)

setReplaceMethod("values", "GenomicFeatureList",
    function(x, value) 
    {
        if (!is(value, "SplitDataFrameList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not a SplitDataFrameList with the ",
                 "same elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        x@unlistData@values <- value
        x
    }
)


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

setReplaceMethod("start", "GenomicFeatureList",
    function(x, value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        start(x@unlistData@ranges) <- value
        x
    }
)

setReplaceMethod("end", "GenomicFeatureList",
    function(x, value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        end(x@unlistData@ranges) <- value
        x
    }
)

setReplaceMethod("width", "GenomicFeatureList",
    function(x, value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        width(x@unlistData@ranges) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SplitDataFrameList methods.
###

setMethod("[", "GenomicFeatureList",
    function(x, i, j, ..., drop)
    {
        if (!missing(i))
            x <- callNextMethod(x = x, i = i)
        if (!missing(j))
            values(x) <- values(x)[, j, drop=FALSE]
        x
    }
)

setReplaceMethod("[", "GenomicFeatureList",
    function(x, i, j, ..., value)
    {
        if (!is(value, "GenomicFeatureList"))
            stop("replacement value must be a GenomicFeatureList object")
        if (!missing(j)) {
            subvalues <- values(x)[i,,drop=FALSE]
            subvalues[,j] <- values(value)
            values(value) <- subvalues
        }
        callNextMethod(x = x, i = i, value = value)
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
