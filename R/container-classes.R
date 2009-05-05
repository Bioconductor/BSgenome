## A container for data in the form of a list of chromosomes.  Each
## sub-element can be anything

setClass("GenomeData", contains = "SimpleTypedList")

### FIXME: ideally there would be some sort of GenomeDescription
### object that encapsulates this information. For now, we store all
### the metadata fields in the SimpleTypedList metadata list. At least
### that way, all the metadata is in one place.
setMethod("providerVersion", "GenomeData",
          function(x) metadata(x)$providerVersion)
setMethod("organism", "GenomeData", function(x) metadata(x)$organism)
setMethod("provider", "GenomeData", function(x) metadata(x)$provider)

GenomeData <- function(listData = list(), providerVersion = NULL,
                       organism = NULL, provider = NULL, metadata = list(),
                       elementMetadata = NULL, ...)
{
  if (!is.list(metadata))
    stop("'metadata' must be a list")
  md <- list(organism = organism, provider = provider,
             providerVersion = providerVersion, ...)
  metadata[names(md)] <- md
  new("GenomeData", listData = listData, metadata = metadata,
      elementMetadata = elementMetadata)
}

if (FALSE)
{
    ## instead of copying, just use IRanges:::labeledLine() for now

### This stuff copy/pasted from IRanges, should go in utils package
ellipsize <- function(obj, width = getOption("width"), sep = " ",
                      ellipsis = "...")
{
  str <- encodeString(obj)
  ## get order selectSome() would print
  half <- seq_len(ceiling(length(obj) / 2))
  ind <- t(cbind(half, length(obj) - half + 1))[seq_along(obj)]
  nc <- cumsum(nchar(str[ind]) + nchar(sep)) - nchar(sep)
  last <- findInterval(width, nc)
  if (length(obj) > last) {
    ## make sure ellipsis fits
    while (last && (nc[last] + nchar(sep)*2^(last>1) + nchar(ellipsis)) > width)
      last <- last - 1
    if (last == 0) ## have to truncate the first element
      str <- paste(substring(str[1], 1, width - nchar(ellipsis)), ellipsis,
                   sep = "")
    else if (last == 1) ## can only show the first
      str <- c(str[1], "...")
    else str <- selectSome(str, last+1)
  }
  paste(str, collapse = sep)
}

labeledLine <- function(label, els, count = TRUE, sep = " ", ellipsis = "...") {
  if (count)
    label <- paste(label, "(", length(els), ")", sep = "")
  label <- paste(label, ": ", sep = "")
  width <- getOption("width") - nchar(label)
  line <- ellipsize(els, width, sep, ellipsis)
  paste(label, line, "\n", sep = "")
}

### End copy/paste

}


setMethod("show", "GenomeData", function(object) {
  cat("A GenomeData instance")
  if (!is.null(organism(object)))
    cat(" for", organism(object))
  if (!is.null(provider(object)) || !is.null(providerVersion(object))) {
    cat("\nbuild: ")
    if (!is.null(provider(object)))
      cat(provider(object), " ", sep = "")
    if (!is.null(providerVersion(object)))
      cat(providerVersion(object))
  }
  cat("\n")
  nms <- names(object)
  if (is.null(nms))
    nms <- seq_len(length(object))
  cat(IRanges:::labeledLine("chromosomes", nms))
})

setValidity("GenomeData",
            function(object) {
              org <- organism(object)
              prov <- provider(object)
              provVer <- providerVersion(object)
              if (!is.null(org) && !isSingleString(org))
                "organism must be a single string or NULL"
              else if (!is.null(prov) && !isSingleString(prov))
                "provider must be a single string or NULL"
              else if (!is.null(provVer) && !isSingleString(provVer))
                "providerVersion must be a single string or NULL"
              else NULL
            })

setClass("GenomeDataList", prototype = prototype(elementType = "GenomeData"),
         contains = "SimpleTypedList")

setValidity("GenomeDataList",
            function(object) {
              ## each element must be a "GenomeData"
                if (!identical(elementType(object), "GenomeData"))
                    return("The elementType(object) is not 'GenomeData'")
                TRUE
            })

GenomeDataList <- function(listData = list(), metadata = list(),
                           elementMetadata = NULL)
{
  new("GenomeDataList", listData = listData, metadata = metadata,
      elementMetadata = elementMetadata)
}

## An apply-type function for working with GenomeData and GenomeDataList objects:

gdApply <- function(...) 
{
    .Deprecated("gdapply",
                msg = "gdApply() is deprecated. Use gdapply() instead [all lowercase].")
    gdapply(...)
}


setGeneric("gdapply",
           function(X, FUN, ...) {
               standardGeneric("gdapply")
           })

setMethod("gdapply",
          signature(X = "GenomeDataList", FUN = "function"),
          function(X, FUN, ...) {
              NX <- length(X)
              new.elements <- vector(mode = "list", length = NX)
              names(new.elements) <- names(X)
              for (i in seq_len(NX))
              {
                  new.elements[[i]] <- gdapply(X[[i]], FUN, ...)
              }
              cls <- lapply(new.elements, class)
              ucl <- unique(unlist(cls))
              if (identical(ucl, "GenomeData"))
                  GenomeDataList(new.elements, metadata(X))
              else new.elements
          })

setMethod("gdapply",
          signature(X = "GenomeData", FUN = "function"),
          function(X, FUN, ...) {
              NX <- length(X)
              new.elements <- vector(mode = "list", length = NX)
              names(new.elements) <- names(X)
              for (i in seq_len(NX))
              {
                  new.elements[[i]] <- FUN(X[[i]], ...)
              }
              cls <- lapply(new.elements, class)
              ucl <- unique(unlist(cls))
              if (length(ucl) == 1)
                  GenomeData(new.elements, metadata = metadata(X))
              else new.elements
          })

setAs("GenomeData", "data.frame",
      function(from) {
          ans <- 
              do.call(rbind, 
                      sapply(names(from),
                             function(chr) {
                                 cbind(as(from[[chr]], "data.frame"), chromosome = chr)
                             }, simplify = FALSE))
          row.names(ans) <- NULL
          ans
      })

setAs("GenomeDataList", "data.frame",
      function(from) {
          ans <- 
              do.call(rbind, 
                      sapply(names(from),
                             function(sample) {
                                 cbind(as(from[[sample]], "data.frame"), sample = sample)
                             }, simplify = FALSE))
          row.names(ans) <- NULL
          ans
      })


