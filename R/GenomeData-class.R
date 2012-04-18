## A container for data in the form of a list of chromosomes.  Each
## sub-element can be anything

setClass("GenomeData", contains = "SimpleList")

### FIXME: ideally there would be some sort of GenomeDescription
### object that encapsulates this information. For now, we store all
### the metadata fields in the SimpleList metadata list.
### At least that way, all the metadata is in one place.
setMethod("providerVersion", "GenomeData",
          function(x) metadata(x)$providerVersion)
setMethod("organism", "GenomeData", function(x) metadata(x)$organism)
setMethod("provider", "GenomeData", function(x) metadata(x)$provider)

GenomeData <- function(listData = list(),
                       providerVersion = metadata[["providerVersion"]],
                       organism = metadata[["organism"]],
                       provider = metadata[["provider"]],
                       metadata = list(),
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

## coersion to data.frame methods
setAs("GenomeData", "data.frame",
      function(from) {
          ind <- names(from)
          if (is.null(ind))
              ind <- seq(length(from))
          ans <-
              do.call(rbind,
                      sapply(ind,
                             function(chr) {
                               cbind(as(from[[chr]], "data.frame"),
                                     chromosome = chr)
                             }, simplify = FALSE))
          row.names(ans) <- NULL
          ans
      })

## seems common to store Ranges inside GenomeData, so tie into IRanges here
setAs("GenomeData", "RangesList", function(from) {
  ans <- do.call("RangesList", lapply(from, as, "Ranges"))
  metadata(ans) <- metadata(from)
  universe(ans) <- providerVersion(from)
  ans
})

setAs("GenomeData", "RangedData", function(from) {
  ans <- do.call("c", lapply(from, as, "RangedData"))
  names(ans) <- names(from)
  metadata(ans) <- metadata(from)
  universe(ans) <- providerVersion(from)
  ans
})

