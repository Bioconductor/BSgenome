gdreduce <- function(...) .Deprecated("gdReduce")

setGeneric("gdReduce",
           function(f, ..., init, right=FALSE, accumulate=FALSE,
                    gdArgs=list()) standardGeneric("gdReduce"),
           signature="...")

setMethod(gdReduce, "GenomeDataList",
          function(f, ..., init, right=FALSE, accumulate=FALSE,
                   gdArgs=list())
{
    args <- list(...)
    if (length(args) > 1L)
        stop("'...' must satisfy 'length(...) == 1' in\n  ",
             "'gdReduce,GenomeDataList-method'")
    gdlist <- args[[1]]
    gdnames <- lapply(gdlist, names)
    nms <- unique(unlist(gdnames))
    listData <-
        if (missing(init))
            lapply(nms, function(nm) {
                ok <- sapply(gdnames, match, x=nm, nomatch=0L) > 0L
                Reduce(f, lapply(gdlist[ok], "[[", nm), 
                       right=right, accumulate=accumulate)
            })
        else
            lapply(nms, function(nm) {
                ok <- sapply(gdnames, match, x=nm, nomatch=0L) > 0L
                Reduce(f, lapply(gdlist[ok], "[[", nm), init=init,
                       right=right, accumulate=accumulate)
            })
    names(listData) <- nms
    do.call(GenomeData, c(list(listData=listData), gdArgs))
})

setMethod(gdReduce, "GenomeData",
          function(f, ..., init, right=FALSE, accumulate=FALSE,
                   gdArgs=list())
{
    if (missing(init))
        gdReduce(f, GenomeDataList(list(...)), right=right,
                 accumulate=accumulate, gdArgs=gdArgs)
    else
        gdReduce(f, GenomeDataList(list(...)), init=init, right=right,
                 accumulate=accumulate, gdArgs=gdArgs)

})
