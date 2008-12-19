


## A container for data in the form of a list of chromosomes.  Each
## sub-element can be anything

setClass("ChromosomeData",
         representation(genome="character", pData="data.frame"),
         contains = "TypedList")

## > showClass("ChromosomeData")
## Class “ChromosomeData”

## Slots:
                                                                      
## Name:           genome        elements           NAMES    elementClass
## Class:       character            list characterORNULL       character

## Extends: "TypedList"


setValidity("ChromosomeData",
            function(object) {
                validElements <- lapply(elements(object), function(x) is(x, elementClass(object)))
                if (!all(validElements))
                    return(sprintf("Not all elements are of class '%s'",
                                   elementClass(object)))
              
                if( !all.equal(names(elements(object)), rownames(object@pData)) )
                    return("names mismatch between elements and pData")
                TRUE
            })
            
            
setClass("SampleChromosomeData",
         contains = "TypedList")


setValidity("SampleChromosomeData",
            function(object) {
                if (!identical(elementClass(object), "ChromosomeData"))
                    return("The elementClass(object) is not 'ChromosomeData'")
                ## each element must be a "ChromosomeData"
                validElements <- lapply(elements(object), function(x) is(x, elementClass(object)))
                if (!all(validElements))
                    return(sprintf("Not all elements are of class '%s'",
                                   elementClass(object)))
                TRUE
            })


