# ===========================================================================
# The BSgenome class
# ---------------------------------------------------------------------------

setClass(
    "BSgenome",
    representation(
        organism="character",
        provider="character",   # "UCSC", "BDGP", etc...
        release="character", 
        source_url="character", 
        seqnames="character",   # names of "single" sequences (e.g. chromosomes),
                                # "single" sequences are stored as DNAString objects
        mseqnames="character",  # names of "multiple" sequences (e.g. upstream),
                                # "multiple" sequences are stored as BStringViews objects
        data_env="environment", # env. where we define the data objects
        cache_env="environment" # private data store
    )
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The 'seqnames', 'mseqnames' and 'names' accessors

setGeneric("seqnames", function(x) standardGeneric("seqnames"))
setMethod("seqnames", "BSgenome", function(x) x@seqnames)

setGeneric("mseqnames", function(x) standardGeneric("mseqnames"))
setMethod("mseqnames", "BSgenome", function(x) x@mseqnames)

setMethod("names", "BSgenome", function(x) c(seqnames(x), mseqnames(x)))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Constructor-like functions and generics

# IMPORTANT: This function will NOT work on Windows if the "SaveImage"
# feature is "on" in the DESCRIPTION file of the "client" package.
# A "client" package will typically create a BSgenome object at load-time
# by have something like this in its R/zzz.R file:
#   organism <- "Caenorhabditis elegans"
#   provider <- "UCSC"
#   release <- "ce2"
#   source_url <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/ce2/bigZips/"
#   seqnames <- c(
#       "chrI",
#       "chrII",
#       "chrIII",
#       "chrIV",
#       "chrV",
#       "chrX",
#       "chrM"
#   }
#   mseqnames <- c(
#       "upstream1000",
#       "upstream2000",
#       "upstream5000"
#   )
#   Celegans <- new("BSgenome", organism, provider,
#                   release, source_url, seqnames, mseqnames,
#                   "BSgenome.Celegans.UCSC.ce2", "extdata")
# The "client" package NEEDS to have "SaveImage: no" in its DESCRIPTION file.
# The problem with "SaveImage: yes" is that "R CMD INSTALL" will execute
# this function BEFORE installing the "inst files", but this function needs
# the "inst files" to be installed BEFORE it can be called!
assignDataToNames <- function(x, package, subdir)
{
    names <- names(x)
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
    function(.Object, organism, provider, release, source_url,
                      seqnames, mseqnames, package, subdir)
    {
        .Object@organism <- organism
        .Object@provider <- provider
        .Object@release <- release
        .Object@source_url <- source_url
        .Object@seqnames <- seqnames
        .Object@mseqnames <- mseqnames
        .Object@data_env <- new.env(parent=emptyenv())
        .Object@cache_env <- new.env(parent=emptyenv())
        assignDataToNames(.Object, package, subdir)
        .Object
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The 'show' method

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
        cat(object@organism, "genome:\n")
        if (length(mseqnames(object)) != 0) {
            cat("\n  Single sequences (DNAString objects, see '?seqnames'):\n")
            indent <- "    "
        } else
            indent <- "  "
        if (length(seqnames(object)) != 0)
            showSequenceIndex(seqnames(object), indent)
        else
            cat(indent, "NONE\n", sep="")
        if (length(mseqnames(object)) != 0) {
            cat("\n  Multiple sequences (BStringViews objects, see '?mseqnames'):\n")
            showSequenceIndex(mseqnames(object), indent)
        }
        cat("\n  (use the '$' or '[[' operator to access a given sequence)\n")
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
