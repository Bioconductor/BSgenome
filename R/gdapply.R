## An apply-type function for working with GenomeData and GenomeDataList objects:
gdApply <- function(...) 
{
    msg <- "gdApply() is defunct. Use gdapply() instead [all lowercase]."
    .Defunct(msg=msg)
}


setGeneric("gdapply",
           function(X, FUN, ...) {
               standardGeneric("gdapply")
           })

setMethod("gdapply",
          signature(X = "GenomeDataList", FUN = "function"),
          function(X, FUN, ...) {
              new.elements <- lapply(X, gdapply, FUN, ...)
              cls <- lapply(new.elements, class)
              ucl <- unique(unlist(cls))
              if (identical(ucl, "GenomeData"))
                  GenomeDataList(new.elements, metadata(X))
              else new.elements
          })

setMethod("gdapply",
          signature(X = "GenomeData", FUN = "function"),
          function(X, FUN, ...) {
              new.elements <- lapply(X, FUN, ...)
              cls <- lapply(new.elements, class)
              ucl <- unique(unlist(cls))
              if (length(ucl) == 1)
                  GenomeData(new.elements, metadata = metadata(X))
              else new.elements
          })

