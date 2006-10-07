# ===========================================================================
# The BSgenome class
# ---------------------------------------------------------------------------

setClass(
    "BSgenome",
    representation(
        organism="character",
        source_provider="character", # "UCSC", "BDGP", etc...
        source_release="character",  # replace the UCSCRelease slot
        source_url="character",      # replace the UCSCBaseUrl slot
        source_files="character",    # replace the UCSCFiles slot
        data_env="environment",      # env. where we define the data objects
        cache_env="environment"      # private data store
    )
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Constructor-like functions and generics

# IMPORTANT: This function will NOT work on Windows if the "SaveImage"
# feature is "on" in the DESCRIPTION file of the "client" package.
# A "client" package will typically create a BSgenome object at load-time
# by have something like this in its R/zzz.R file:
#   organism <- "Caenorhabditis elegans"
#   source_provider <- "UCSC"
#   source_release <- "ce2"
#   source_url <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/ce2/bigZips/"
#   source_files <- c(
#       "chrI",
#       "chrII",
#       "chrIII",
#       "chrIV",
#       "chrV",
#       "chrM",
#       "chrX",
#       "upstream1000",
#       "upstream2000",
#       "upstream5000"
#   )
#   Celegans <- new("BSgenome", organism, source_provider,
#                   source_release, source_url, source_files,
#                   "BSgenome.Celegans.UCSC.ce2", "extdata")
# The "client" package NEEDS to have "SaveImage: no" in its DESCRIPTION file.
# The problem with "SaveImage: yes" is that "R CMD INSTALL" will execute
# this function BEFORE installing the "inst files", but this function needs
# the "inst files" to be installed BEFORE it can be called!
assignDataToNames <- function(x, package, subdir)
{
    names <- x@source_files
    data_env <- x@data_env
    cache_env <- x@cache_env

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
    function(.Object, organism, source_provider, source_release, source_url,
                      source_files, package, subdir)
    {
        .Object@organism <- organism
        .Object@source_provider <- source_provider
        .Object@source_release <- source_release
        .Object@source_url <- source_url
        .Object@source_files <- source_files
        .Object@data_env <- new.env(parent=emptyenv())
        .Object@cache_env <- new.env(parent=emptyenv())
        assignDataToNames(.Object, package, subdir)
        .Object
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The 'names' and 'show' methods

setMethod("names", "BSgenome", function(x) x@source_files)

setMethod("show", "BSgenome",
    function(object)
    {
        cat(object@organism, "genome:\n")
        ans <- show(names(object))
        cat("(use the '$' or '[[' operator to access a given chromosome)\n")
        ans
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subsetting

setMethod("length", "BSgenome", function(x) length(names(x)))

setMethod("[[", "BSgenome",
    function(x, i, j, ...)
    {
        # 'x' is guaranteed to be a "BSgenome" object (if it's not, then the
        # method dispatch algo will not even call this method), so nargs() is
        # guaranteed to be >= 1
        if (nargs() >= 3)
            stop("too many subscripts")
        subscripts <- list(...)
        if (!missing(i))
            subscripts$i <- i
        if (!missing(j))
            subscripts$j <- j
        # At this point, 'subscripts' should be guaranteed
        # to be of length <= 1
        if (length(subscripts) == 0)
            stop("no index specified")
        i <- subscripts[[1]]
        if (length(i) < 1)
            stop("attempt to select less than one element")
        if (length(i) > 1)
            stop("attempt to select more than one element")
        if (is.character(i))
            name <- match.arg(i, names(x))
        else
            name <- names(x)[i]
        get(name, envir=x@data_env)
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Other functions and generics

setGeneric("unload", function(x, what) standardGeneric("unload"))

setMethod("unload", "BSgenome",
    function(x, what)
    {
        if (missing(what))
            what <- ls(x@cache_env) ## everything
        remove(list=what, envir=x@cache_env)
    }
)
