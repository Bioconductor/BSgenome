### =========================================================================
### The "GenomeDescription" class
### -------------------------------------------------------------------------

setClass("GenomeDescription",
    representation(
        ## organism: "Homo sapiens", "Mus musculus", etc...
        organism="character",

        ## species: "Human", "Mouse", etc...
        species="character",

        ## provider: "UCSC", "BDGP", etc...
        provider="character",

        ## provider_version: "hg18", "mm8", "sacCer1", etc...
        provider_version="character",

        ## release_date: "Mar. 2006", "Feb. 2006", "Oct. 2003", etc...
        release_date="character",

        ## release_name: "NCBI Build 36.1", "NCBI Build 36",
        ## "SGD 1 Oct 2003 sequence", etc...
        release_name="character"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("organism", function(x) standardGeneric("organism"))
setMethod("organism", "GenomeDescription", function(x) x@organism)

setGeneric("species", function(x) standardGeneric("species"))
setMethod("species", "GenomeDescription", function(x) x@species)

setGeneric("provider", function(x) standardGeneric("provider"))
setMethod("provider", "GenomeDescription", function(x) x@provider)

setGeneric("providerVersion", function(x) standardGeneric("providerVersion"))
setMethod("providerVersion", "GenomeDescription", function(x) x@provider_version)

setGeneric("releaseDate", function(x) standardGeneric("releaseDate"))
setMethod("releaseDate", "GenomeDescription", function(x) x@release_date)

setGeneric("releaseName", function(x) standardGeneric("releaseName"))
setMethod("releaseName", "GenomeDescription", function(x) x@release_name)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

setValidity("GenomeDescription",
    function(object)
    {
        unlist(
            lapply(slotNames("GenomeDescription"),
                function(slotname)
                {
                    slotval <- slot(object, slotname)
                    if (isSingleStringOrNA(slotval))
                        return(NULL)
                    problem <- paste("slot '", slotname, "' must be a ",
                                     "single string (or NA)", sep="")
                    return(problem)
                }
            )
        )
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor-like functions
###

GenomeDescription <- function(organism, species,
                              provider, provider_version,
                              release_date, release_name)
{
    new("GenomeDescription",
        organism=organism,
        species=species,
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

setMethod("show", "GenomeDescription",
    function(object)
    {
        PREFIX <- "| "
        mystrwrap <- function(line)
            writeLines(strwrap(line, width=getOption("width")+1,
                               exdent=0, prefix=PREFIX))
        cat(PREFIX, "organism: ", organism(object), " (",  species(object), ")\n", sep="")
        cat(PREFIX, "provider: ", provider(object), "\n", sep="")
        cat(PREFIX, "provider version: ", providerVersion(object), "\n", sep="")
        cat(PREFIX, "release date: ", releaseDate(object), "\n", sep="")
        cat(PREFIX, "release name: ", releaseName(object), "\n", sep="")
    }
)

