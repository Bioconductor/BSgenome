### =========================================================================
### GenomicFeature objects
### -------------------------------------------------------------------------
###
### Class definition

setClass("GenomicFeature", contains = "Sequence",
         representation(seqnames = "Rle",
                        ranges = "IRanges",
                        strand = "Rle",
                        values = "DataFrame"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GenomicFeature.slots <- function(x)
{
    n <- length(x@seqnames)
    if ((length(x@ranges) != n) || (length(x@strand) != n) ||
        (nrow(x@values) != n))
        "slot lengths are not all equal"
    else
        NULL
}

.valid.GenomicFeature.seqnames <- function(x)
{
    if (!is.character(runValue(x@seqnames)))
        "slot 'seqnames' should be a 'character' Rle"
    else
        NULL
}

.valid.GenomicFeature.strand <- function(x)
{
    if (!is.factor(runValue(x@strand)) ||
        !identical(levels(runValue(x@strand)), levels(strand())))
        paste("slot 'strand' should be a 'factor' Rle with levels c(",
              paste('"', levels(strand()), '"', sep = "", collapse = ", "),
                    ")", sep = "")
    else
        NULL
}


.valid.GenomicFeature.values <- function(x)
{
    msg <- NULL
    if (any(c("seqnames", "ranges", "strand", "start", "end", "width",
              "feature") %in% colnames(x@values)))
        msg <-
          paste("slot 'values' cannot use \"seqnames\", \"ranges\",",
                "\"strand\", \"start\", \"end\", \"width\", or \"feature\"",
                "as column names")
    if (!is.null(rownames(x@values)))
        msg <- c(msg, "slot 'values' cannot contain row names")
    msg
}

.valid.GenomicFeature <- function(x)
{
    c(.valid.GenomicFeature.slots(x),
      .valid.GenomicFeature.seqnames(x),
      .valid.GenomicFeature.strand(x),
      .valid.GenomicFeature.values(x))
}

setValidity2("GenomicFeature", .valid.GenomicFeature)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GenomicFeature <-
function(seqnames = Rle(), ranges = IRanges(),
         strand = Rle(NA_character_, length(seqnames)), ...)
{
    if (!is(seqnames, "Rle"))
        seqnames <- Rle(seqnames)
    if (!is.character(runValue(seqnames)))
        runValue(seqnames) <- as.character(runValue(seqnames))

    if (!is(ranges, "IRanges"))
        ranges <- as(ranges, "IRanges")

    if (!is(strand, "Rle"))
        strand <- Rle(strand)
    if (!is.factor(runValue(strand)))
        runValue(strand) <- strand(runValue(strand))

    values <- DataFrame(...)
    if (ncol(values) == 0)
        values <- new("DataFrame", nrows = length(seqnames))
    if (!is.null(rownames(values))) {
        if (!is.null(names(ranges)))
            names(ranges) <- rownames(values)
        rownames(values) <- NULL
    }

    new("GenomicFeature", seqnames = seqnames, ranges = ranges, strand = strand,
        values = values)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setMethod("as.data.frame", "GenomicFeature",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        ranges <- ranges(x)
        if (missing(row.names))
            row.names <- names(x)
        if (!is.null(names(x)))
            names(x) <- NULL
        data.frame(seqnames = as.vector(seqnames(x)),
                   start = start(x),
                   end = end(x),
                   width = width(x),
                   strand = as.vector(strand(x)),
                   as.data.frame(values(x)),
                   row.names = row.names,
                   stringsAsFactors = FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("seqnames", "GenomicFeature", function(x) x@seqnames)
setMethod("ranges", "GenomicFeature", function(x, ...) x@ranges)
setMethod("strand", "GenomicFeature", function(x) x@strand)
setMethod("values", "GenomicFeature",
    function(x, ...)
    {
        ans <- x@values
        if (!is.null(names(x)))
            rownames(ans) <- names(x)
        ans
    }
)

setReplaceMethod("seqnames", "GenomicFeature",
    function(x, value) 
    {
        if (!is(value, "Rle"))
            value <- Rle(value)
        n <- length(x)
        k <- length(value)
        if (k != n) {
            if ((k == 0) || (k > n) || (n %% k != 0))
                stop(k, " elements in value to replace ", n, "elements")
            value <- rep(value, length.out = n)
        }
        initialize(x, seqnames = value)
    }
)
setReplaceMethod("ranges", "GenomicFeature",
    function(x, value) 
    {
        if (!is(value, "IRanges"))
            value <- as(value, "IRanges")
        n <- length(x)
        k <- length(value)
        if (k != n) {
            if ((k == 0) || (k > n) || (n %% k != 0))
                stop(k, " elements in value to replace ", n, "elements")
            value <- rep(value, length.out = n)
        }
        initialize(x, ranges = value)
    }
)
setReplaceMethod("strand", "GenomicFeature",
    function(x, value) 
    {
        if (!is(value, "Rle"))
            value <- Rle(value)
        if (!is.factor(runValue(value)))
            runValue(value) <- strand(runValue(value))
        n <- length(x)
        k <- length(value)
        if (k != n) {
            if ((k == 0) || (k > n) || (n %% k != 0))
                stop(k, " elements in value to replace ", n, "elements")
            value <- rep(value, length.out = n)
        }
        initialize(x, strand = value)
    }
)
setReplaceMethod("values", "GenomicFeature",
    function(x, value)
    {
        if (is.null(value))
            value <- new("DataFrame", nrows = length(x))
        else if (!is(value, "DataFrame"))
            value <- DataFrame(value)
        if (!is.null(rownames(value)))
            rownames(value) <- NULL
        n <- length(x)
        k <- nrow(value)
        if (k != n) {
            if ((k == 0) || (k > n) || (n %% k != 0))
                stop(k, " rows in value to replace ", n, "rows")
            value <- value[rep(seq_len(k), length.out = n), , drop=FALSE]
        }
        initialize(x, values = value)
    }
)

setMethod("names", "GenomicFeature", function(x) names(x@ranges))
setReplaceMethod("names", "GenomicFeature",
    function(x, value)
    {
        names(x@ranges) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ranges methods.
###

setMethod("start", "GenomicFeature", function(x, ...) start(x@ranges))
setMethod("end", "GenomicFeature", function(x, ...) end(x@ranges))
setMethod("width", "GenomicFeature", function(x) width(x@ranges))

setReplaceMethod("start", "GenomicFeature",
    function(x, value)
    {
        start(x@ranges) <- value
        x
    }
)

setReplaceMethod("end", "GenomicFeature",
    function(x, value)
    {
        end(x@ranges) <- value
        x
    }
)

setReplaceMethod("width", "GenomicFeature",
    function(x, value)
    {
        width(x@ranges) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DataTable methods.
###

setMethod("ncol", "GenomicFeature", function(x) ncol(x@values))

setMethod("colnames", "GenomicFeature",
    function(x, do.NULL = TRUE, prefix = "col") 
        colnames(x@values, do.NULL = do.NULL, prefix = prefix))
setReplaceMethod("colnames", "GenomicFeature",
    function(x, value)
    {
        colnames(x@values) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sequence methods.
###

setMethod("[", "GenomicFeature",
    function(x, i, j, ..., drop)
    {
        if (!missing(i)) {
            iInfo <- IRanges:::.bracket.Index(i, names(x), length(x))
            if (!is.null(iInfo[["msg"]]))
                stop(iInfo[["msg"]])
            i <- iInfo[["idx"]]
        }
        ans <-
          initialize(x, seqnames = x@seqnames[i], ranges = x@ranges[i],
                     strand = x@strand[i], values = x@values[i, j, drop=FALSE])
        nms <- names(ans)
        if (!is.null(nms)) {
            whichEmpty <- which(nms == "")
            nms[whichEmpty] <- as.character(whichEmpty)
            nms2 <- make.unique(nms)
            if (length(whichEmpty) > 0 || !identical(nms, nms2))
                names(ans) <- nms2
        }
        ans
    }
)

setReplaceMethod("[", "GenomicFeature",
    function(x, i, j, ..., value)
    {
        if (!is(value, "GenomicFeature"))
            stop("replacement value must be a GenomicFeature object")
        seqnames <- x@seqnames
        ranges <- x@ranges
        strand <- x@strand
        values <- x@values
        if (missing(i)) {
            seqnames[] <- value@seqnames
            ranges[] <- value@ranges
            strand[] <- value@strand
            values[,j] <- value@values
        } else {
            iInfo <- IRanges:::.bracket.Index(i, names(x), length(x))
            if (!is.null(iInfo[["msg"]]))
                stop(iInfo[["msg"]])
            i <- iInfo[["idx"]]
            seqnames[i] <- value@seqnames
            ranges[i] <- value@ranges
            strand[i] <- value@strand
            values[i,j] <- value@values
        }
        initialize(x, seqnames = seqnames, ranges = ranges,
                   strand = strand, values = values)
    }
)

setMethod("c", "GenomicFeature",
    function(x, ..., recursive = FALSE)
    {
        if (recursive)
            stop("'recursive' mode not supported")
        args <- list(x, ...)
        ans <-
          initialize(x,
                     seqnames = do.call(c, lapply(args, slot, "seqnames")),
                     ranges = do.call(c, lapply(args, slot, "ranges")),
                     strand = do.call(c, lapply(args, slot, "strand")),
                     values = do.call(rbind, lapply(args, slot, "values")))
        nms <- names(ans)
        if (!is.null(nms)) {
            whichEmpty <- which(nms == "")
            nms[whichEmpty] <- as.character(whichEmpty)
            nms2 <- make.unique(nms)
            if (length(whichEmpty) > 0 || !identical(nms, nms2))
                names(ans) <- nms2
        }
        ans
    }
)

setMethod("length", "GenomicFeature", function(x) length(x@seqnames))

setMethod("rev", "GenomicFeature",
    function(x)
    {
        if (length(x) == 0)
            x
        else
            initialize(x, seqnames = rev(x@seqnames), ranges = rev(x@ranges),
                       strand = rev(x@strand),
                       values = x@values[length(x):1, , drop=FALSE])
    }
)

setMethod("seqselect", "GenomicFeature",
    function(x, start = NULL, end = NULL, width = NULL)
    {
        ans <-
          initialize(x,
                     seqnames =
                     seqselect(x@seqnames, start = start, end = end, width = width),
                     ranges =
                     seqselect(x@ranges, start = start, end = end, width = width),
                     strand =
                     seqselect(x@strand, start = start, end = end, width = width),
                     values =
                     seqselect(x@values, start = start, end = end, width = width))
        nms <- names(ans)
        if (!is.null(nms)) {
            whichEmpty <- which(nms == "")
            nms[whichEmpty] <- as.character(whichEmpty)
            nms2 <- make.unique(nms)
            if (length(whichEmpty) > 0 || !identical(nms, nms2))
                names(ans) <- nms2
        }
        ans
    }
)

setReplaceMethod("seqselect", "GenomicFeature",
    function(x, start = NULL, end = NULL, width = NULL, value)
    {
        if (!is(value, "GenomicFeature"))
            stop("replacement value must be a GenomicFeature object")
        seqnames <- as.vector(x@seqnames)
        ranges <- x@ranges
        strand <- as.vector(x@strand)
        values <- x@values
        seqselect(seqnames, start = start, end = end, width = width) <- value@seqnames
        seqselect(ranges, start = start, end = end, width = width) <- value@ranges
        seqselect(strand, start = start, end = end, width = width) <- value@strand
        seqselect(values, start = start, end = end, width = width) <- value@values
        initialize(x, seqnames = Rle(seqnames), ranges = ranges, 
                   strand = Rle(strand), values = values)
    }
)

setMethod("split", "GenomicFeature",
    function(x, f, drop = FALSE)
    {
        if (missing(f)) {
            nms <- names(x)
            if (is.null(nms))
                f <- seq_len(length(x))
            else
                f <- factor(nms, levels = nms)
        }
        IRanges:::newCompressedList("GenomicFeatureList", x, splitFactor = f,
                                    drop = drop)
    }
)

setMethod("window", "GenomicFeature",
    function(x, start = NA, end = NA, width = NA,
             frequency = NULL, delta = NULL, ...)
    {
        initialize(x,
                   seqnames =
                   window(x@seqnames, start = start, end = end, width = width,
                          frequency = frequency, delta = delta),
                   ranges =
                   window(x@ranges, start = start, end = end, width = width,
                          frequency = frequency, delta = delta),
                   strand =
                   window(x@strand, start = start, end = end, width = width,
                          frequency = frequency, delta = delta),
                   values =
                   window(x@values, start = start, end = end, width = width,
                          frequency = frequency, delta = delta))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

setMethod("show", "GenomicFeature",
    function(object)
    {
        lo <- length(object)
        nc <- ncol(object)
        cat(class(object), " with ",
            lo, ifelse(lo == 1, " range and ", " ranges and "),
            nc, ifelse(nc == 1, " values column\n", " values columns\n"),
            sep = "")
        if (lo > 0) {
            k <- ifelse(lo <= 12L, lo, min(lo, 10L))
            subset  <- object[seq_len(k)]
            out <-
              cbind(seqnames = as.character(seqnames(subset)),
                    ranges = IRanges:::showAsCell(ranges(subset)),
                    strand = as.character(strand(subset)),
                    "|" = rep.int("|", k))
            if (nc > 0)
                out <-
                  cbind(out,
                        as.matrix(format.data.frame(do.call(data.frame,
                                                            lapply(values(subset),
                                                                    IRanges:::showAsCell)))))
            if (is.null(names(subset)))
                rownames(out) <- seq_len(k)
            else
                rownames(out) <- names(subset)
            classinfo <-
              matrix(c("<Rle>", "<IRanges>", "<Rle>", "|",
                       unlist(lapply(values(subset), function(x)
                                     paste("<", class(x), ">", sep = "")),
                              use.names = FALSE)), nrow = 1,
                     dimnames = list("", colnames(out)))
            out <- rbind(classinfo, out)
            print(out, quote = FALSE, right = TRUE)
            diffK <- lo - k
            if (diffK > 0)
                cat("...\n<", diffK,
                    ifelse(diffK == 1, " more range>\n", " more ranges>\n"),
                    sep="")
        }
    }
)
