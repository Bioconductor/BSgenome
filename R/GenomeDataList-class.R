## Conceptually a list of GenomeData objects

setClass("GenomeDataList", prototype = prototype(elementType = "GenomeData"),
         contains = "SimpleList")

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

setAs("GenomeDataList", "data.frame",
      function(from) {
          ind <- names(from)
          if (is.null(ind))
              ind <- seq(length(from))
          ans <- 
              do.call(rbind, 
                      sapply(ind,
                             function(sample) {
                                 cbind(as(from[[sample]], "data.frame"),
                                       sample = sample)
                             }, simplify = FALSE))
          row.names(ans) <- NULL
          ans
      })

setAs("GenomeDataList", "RangedDataList", function(from) {
  ans <- do.call(RangedDataList, lapply(from, as, "RangedData"))
  metadata(ans) <- metadata(from)
  elementMetadata(ans) <- elementMetadata(from)
  ans
})

