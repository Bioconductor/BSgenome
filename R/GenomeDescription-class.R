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

setGeneric("bsgenomeName", function(x) standardGeneric("bsgenomeName"))
setMethod("bsgenomeName", "GenomeDescription",
    function(x)
    {
        part1 <- "BSgenome"
        tmp <- strsplit(organism(x), " ", fixed=TRUE)[[1L]]
        part2 <- paste(substr(tmp[1L], start=1L, stop=1L), tmp[2L], sep="")
        part3 <- provider(x)
        part4 <- providerVersion(x)
        paste(part1, part2, part3, part4, sep=".")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

setValidity("GenomeDescription",
    function(object)
    {
        .validSlot <- function(slotname)
        {
            slotval <- slot(object, slotname)
            if (isSingleStringOrNA(slotval))
                return(NULL)
            problem <- paste("slot '", slotname, "' must be a ",
                             "single string (or NA)", sep="")
            return(problem)
        }
        unlist(lapply(slotNames("GenomeDescription"), .validSlot))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor-like functions
###

GenomeDescription <- function(organism, species,
                              provider, provider_version,
                              release_date, release_name)
{
    if (identical(organism, "NA")) organism <- NA_character_
    if (identical(species, "NA")) species <- NA_character_
    if (identical(release_date, "NA")) release_date <- NA_character_
    if (identical(release_name, "NA")) release_name <- NA_character_
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

