### =========================================================================
### The "BSgenome" class
### -------------------------------------------------------------------------

setClass(
    "BSgenome",
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

        ## release_name: "NCBI Build 36.1", "NCBI Build 36", "SGD 1 Oct 2003 sequence", etc...
        release_name="character",

        ## source_url: permanent URL to the place where the FASTA files used
        ## to produce the sequences below can be found (and downloaded)
        source_url="character",

        ## seqnames: names of "single" sequences (e.g. chromosomes),
        ##           "single" sequences are stored as DNAString objects
        seqnames="character",

        ## mseqnames: names of "multiple" sequences (e.g. upstream),
        ##            "multiple" sequences are stored as BStringViews objects
        mseqnames="character",

        ## .data_env: env. where we define the data objects
        .data_env="environment",

        ## .cache_env: private data store
        .cache_env="environment"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###

setGeneric("organism", function(x) standardGeneric("organism"))
setMethod("organism", "BSgenome", function(x) x@organism)

setGeneric("species", function(x) standardGeneric("species"))
setMethod("species", "BSgenome", function(x) x@species)

setGeneric("provider", function(x) standardGeneric("provider"))
setMethod("provider", "BSgenome", function(x) x@provider)

setGeneric("providerVersion", function(x) standardGeneric("providerVersion"))
setMethod("providerVersion", "BSgenome", function(x) x@provider_version)

setGeneric("releaseDate", function(x) standardGeneric("releaseDate"))
setMethod("releaseDate", "BSgenome", function(x) x@release_date)

setGeneric("releaseName", function(x) standardGeneric("releaseName"))
setMethod("releaseName", "BSgenome", function(x) x@release_name)

setGeneric("sourceUrl", function(x) standardGeneric("sourceUrl"))
setMethod("sourceUrl", "BSgenome", function(x) x@source_url)

setGeneric("seqnames", function(x) standardGeneric("seqnames"))
setMethod("seqnames", "BSgenome", function(x) x@seqnames)

setGeneric("mseqnames", function(x) standardGeneric("mseqnames"))
setMethod("mseqnames", "BSgenome", function(x) x@mseqnames)

setMethod("names", "BSgenome", function(x) c(seqnames(x), mseqnames(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor-like functions and generics
###

### IMPORTANT: The "assignDataToNames" function below will NOT work on Windows
### if the "SaveImage" feature is "on" in the DESCRIPTION file of the "client"
### package. A "client" package will typically create a BSgenome object at
### load-time by having something like this in its R/zzz.R file:
###
###   Celegans <- new("BSgenome",
###       organism="Caenorhabditis elegans",
###       species="C. elegans",
###       provider="UCSC",
###       provider_version="ce2",
###       release_date="Mar. 2004",
###       release_name="WormBase v. WS120",
###       source_url="ftp://hgdownload.cse.ucsc.edu/goldenPath/ce2/bigZips/",
###       seqnames=c(
###           "chrI",
###           "chrII",
###           "chrIII",
###           "chrIV",
###           "chrV",
###           "chrX",
###           "chrM"
###       ),
###       mseqnames=c(
###           "upstream1000",
###           "upstream2000",
###           "upstream5000"
###       ),
###       package="BSgenome.Celegans.UCSC.ce2",
###       subdir="extdata"
###   )
###
### The "client" package NEEDS to have "SaveImage: no" in its DESCRIPTION file.
### The problem with "SaveImage: yes" is that "R CMD INSTALL" will execute
### this function BEFORE installing the "inst files", but this function needs
### the "inst files" to be installed BEFORE it can be called!
assignDataToNames <- function(x, package, subdir)
{
    names <- names(x)
    data_env <- x@.data_env
    cache_env <- x@.cache_env

    addCachedItem <- function(name, file)
    {
        ## Add a lazy-loading active binding thingie named 'name'
        ## containing the single object serialized in file 'file'.
        force(name)
        force(file)
        getter <- function()
        {
            if (!is.null(cache_env[[name]]))
              return(cache_env[[name]])
            found <- load(file, envir=cache_env)
            if (name != found)
              cache_env[[name]] <- get(found, envir=cache_env)
              #stop("bad data file: name of objects must match:\n",
              #     sQuote(name), " != ", sQuote(found))
            cache_env[[name]]
        }
        makeActiveBinding(name, getter, data_env)
    }

    for (name in names) {
        file <- system.file(subdir, paste(name, ".rda", sep=""), package=package)
        addCachedItem(name, file)
    }
}

setMethod("initialize", "BSgenome",
    function(.Object, organism, species, provider, provider_version,
                      release_date, release_name, source_url,
                      seqnames, mseqnames, package, subdir)
    {
        .Object@organism <- organism
        .Object@species <- species
        .Object@provider <- provider
        .Object@provider_version <- provider_version
        .Object@release_date <- release_date
        .Object@release_name <- release_name
        .Object@source_url <- source_url
        .Object@seqnames <- seqnames
        .Object@mseqnames <- mseqnames
        .Object@.data_env <- new.env(parent=emptyenv())
        .Object@.cache_env <- new.env(parent=emptyenv())
        assignDataToNames(.Object, package, subdir)
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

setMethod("show", "BSgenome",
    function(object)
    {
        showSequenceIndex <- function(names, indent)
        {
            index_width <- 81 - nchar(indent)
            col_width <- max(nchar(names))
            ncols <- index_width %/% (col_width + 2)
            col <- 1
            for (name in names) {
                if (col == 1) cat(indent)
                cat(format(name, width=col_width))
                if (col == ncols) {
                    cat("\n")
                    col <- 1
                } else {
                    cat("  ")
                    col <- col + 1
                }
            }
            if (col != 1) cat("\n")
        }
        cat(object@species, "genome:\n")
        cat("- organism: ", object@organism, "\n", sep="")
        cat("- provider: ", object@provider, "\n", sep="")
        cat("- provider version: ", object@provider_version, "\n", sep="")
        cat("- release date: ", object@release_date, "\n", sep="")
        cat("- release name: ", object@release_name, "\n", sep="")
        if (length(mseqnames(object)) != 0) {
            cat("- single sequences (DNAString objects, see '?seqnames'):\n")
            indent <- "    "
        } else
            indent <- "  "
        if (length(seqnames(object)) != 0)
            showSequenceIndex(seqnames(object), indent)
        else
            cat(indent, "NONE\n", sep="")
        if (length(mseqnames(object)) != 0) {
            cat("- multiple sequences (BStringViews objects, see '?mseqnames'):\n")
            showSequenceIndex(mseqnames(object), indent)
        }
        cat("  (use the '$' or '[[' operator to access a given sequence)\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("length", "BSgenome", function(x) length(names(x)))

setMethod("[[", "BSgenome",
    function(x, i, j, ...)
    {
        ## 'x' is guaranteed to be a "BSgenome" object (if it's not, then the
        ## method dispatch algo will not even call this method), so nargs() is
        ## guaranteed to be >= 1
        if (nargs() >= 3)
            stop("too many subscripts")
        subscripts <- list(...)
        if (!missing(i))
            subscripts$i <- i
        if (!missing(j))
            subscripts$j <- j
        ## At this point, 'subscripts' should be guaranteed
        ## to be of length <= 1
        if (length(subscripts) == 0)
            stop("no index specified")
        i <- subscripts[[1]]
        if (length(i) < 1)
            stop("attempt to select less than one element")
        if (length(i) > 1)
            stop("attempt to select more than one element")
        if (is.character(i)) {
            name <- try(match.arg(i, names(x)), silent=TRUE)
            if (is(name, "try-error"))
                stop("no such sequence")
        } else {
            if (!is.numeric(i) || is.na(i))
                stop("no such sequence")
            i <- as.integer(i)
            if (i < 1 || length(x) < i)
                stop("no such sequence")
            name <- names(x)[i]
        }
        get(name, envir=x@.data_env)
    }
)

setReplaceMethod("[[", "BSgenome",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a \"BSgenome\" object")
    }
)

setMethod("$", "BSgenome",
    function(x, name) x[[name]]
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other functions and generics
###

setGeneric("unload", function(x, what) standardGeneric("unload"))

setMethod("unload", "BSgenome",
    function(x, what)
    {
        if (missing(what))
            what <- ls(x@.cache_env) ## everything
        remove(list=what, envir=x@.cache_env)
    }
)

