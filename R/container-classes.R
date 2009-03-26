## A container for data in the form of a list of chromosomes.  Each
## sub-element can be anything

setClass("GenomeData",
         representation(organism = "characterORNULL",
                        provider = "characterORNULL"),
         contains = "AnnotatedList")

### FIXME: ideally there would be some sort of GenomeDescription
### object that encapsulates this information. For historical reasons,
### the 'providerVersion' is stored in the metadata, while organism
### and provider are explicit slots.
setMethod("providerVersion", "GenomeData",
          function(x) metadata(x)$providerVersion)
setMethod("organism", "GenomeData", function(x) x@organism)
setMethod("provider", "GenomeData", function(x) x@provider)

GenomeData <- function(elements = list(), providerVersion = NULL,
                       organism = NULL, provider = NULL,
                       elementMetadata = NULL, ...)
{
  new("GenomeData", elements = elements,
      metadata = list(providerVersion = providerVersion),
      elementMetadata = elementMetadata, organism = organism,
      provider = provider)
}

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

## Name:         elements             NAMES      elementClass    elementLengths
## Class:            list   characterORNULL         character           integer
## Name:         compress
## Class:         logical

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
  cat(labeledLine("chromsomes", nms))
})

setValidity("GenomeData",
            function(object) {
              org <- organism(object)
              prov <- provider(object)
              if (!is.null(org) && !isSingleString(org))
                "organism must be a single string or NULL"
              else if (!is.null(prov) && !isSingleString(prov))
                "provider must be a single string or NULL"
              else NULL
            })

setClass("GenomeDataList", prototype = prototype(elementClass = "GenomeData"),
         contains = "AnnotatedList")

setValidity("GenomeDataList",
            function(object) {
              ## each element must be a "GenomeData"
                if (!identical(elementClass(object), "GenomeData"))
                    return("The elementClass(object) is not 'GenomeData'")
                TRUE
            })

GenomeDataList <- function(elements = list(), annotation = NULL,
                           elementMetadata = NULL)
{
  new("GenomeDataList", elements = elements, annotation = annotation,
      elementMetadata = elementMetadata)
}

## An apply-type function for working with GenomeData and GenomeDataList objects:


setGeneric("gdApply",
           function(X, FUN, ...) {
               standardGeneric("gdApply")
           })

setMethod("gdApply",
          signature(X = "GenomeDataList", FUN = "function"),
          function(X, FUN, ...) {
              NX <- length(X)
              new.elements <- vector(mode = "list", length = NX)
              names(new.elements) <- names(X)
              for (i in seq_len(NX))
              {
                  new.elements[[i]] <- gdApply(X[[i]], FUN, ...)
              }
              cls <- lapply(new.elements, class)
              ucl <- unique(unlist(cls))
              if (identical(ucl, "GenomeData"))
                  GenomeDataList(new.elements)
              else new.elements
          })

setMethod("gdApply",
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
              if (length(ucl) = = 1)
                  GenomeData(new.elements, organism = X@organism)
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


