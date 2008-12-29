## A container for data in the form of a list of chromosomes.  Each
## sub-element can be anything

## AnnotatedList provides "annotation" slot, which should be
## provider-version (e.g. hg18). AnnotatedList also provides
## elementMetadata (formerly pData).  We add organism (e.g. Mus
## musculus) and provider (e.g. UCSC). From this, it should be possible
## to uniquely lookup a BSgenome genome. But it might be simpler in the
## long-run to have some unique "BSgenome ID" that we could stick
## into the 'annotation' slot and look everything up with that.
setClass("GenomeData",
         representation(organism = "characterORNULL",
                        provider = "characterORNULL"),
         contains = "AnnotatedList")

setMethod("providerVersion", "GenomeData", function(x) x@annotation)
setMethod("organism", "GenomeData", function(x) x@organism)
setMethod("provider", "GenomeData", function(x) x@provider)

## > showClass("GenomeData")
## Class "GenomeData"

## Slots:
  
## Name:         organism          provider        annotation   elementMetadata
## Class: characterORNULL   characterORNULL   characterORNULL  XDataFrameORNULL

## Name:         elements             NAMES      elementClass elementCumLengths
## Class:            list   characterORNULL         character           integer

## Name:  compressedIndices      compressible
## Class:           integer           logical

## Extends: 
## Class "AnnotatedList", directly
## Class "TypedList", by class "AnnotatedList", distance 2

isSingleString <- function(x)
{
  is.character(x) && length(x) == 1 && !is.na(x)
}

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
