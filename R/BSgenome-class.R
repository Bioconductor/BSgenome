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

        ## for SNPs injection
        SNPlocs_pkgname="character",

        ## seqnames: names of "single" sequences (e.g. chromosomes)
        seqnames="character",

        ## mseqnames: names of "multiple" sequences (e.g. upstream)
        mseqnames="character",

        package="character",
        subdir="character",

        .activebindings_env="environment",
        .datacache_env="environment"
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

setGeneric("SNPlocs_pkgname", function(x) standardGeneric("SNPlocs_pkgname"))
setMethod("SNPlocs_pkgname", "BSgenome",
    function(x) if (length(x@SNPlocs_pkgname) == 0) NULL else x@SNPlocs_pkgname
)

setGeneric("seqnames", function(x) standardGeneric("seqnames"))
setMethod("seqnames", "BSgenome", function(x) x@seqnames)

setGeneric("mseqnames", function(x) standardGeneric("mseqnames"))
setMethod("mseqnames", "BSgenome", function(x) x@mseqnames)

setMethod("names", "BSgenome", function(x) c(seqnames(x), mseqnames(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor-like functions and generics
###

### IMPORTANT: The ".assignDataToNames" function below will NOT work on Windows
### if the "SaveImage" feature is "on" in the DESCRIPTION file of the "client"
### package. A "client" package will typically create a BSgenome object at
### load-time by having something like this in its R/zzz.R file:
###
###   Celegans <- BSgenome(
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
.assignDataToNames <- function(x, getSNPcount=NULL, getSNPlocs=NULL)
{
    names <- names(x)
    activebindings_env <- x@.activebindings_env
    datacache_env <- x@.datacache_env

    addCachedItem <- function(name, file)
    {
        ## Add a lazy-loading active binding thingie named 'name'
        ## containing the single object serialized in file 'file'.
        force(name)
        force(file)
        getter <- function()
        {
            if (exists(name, envir=datacache_env, inherits=FALSE))
                return(datacache_env[[name]])
            found <- load(file, envir=datacache_env)
            if (name != found) {
                datacache_env[[name]] <- get(found, envir=datacache_env)
            }

            ## Inject the SNPs, if any
            if (!is.null(getSNPcount) && name %in% names(getSNPcount())) {
                snps <- getSNPlocs(name)
                if (nrow(snps) != getSNPcount()[name])
                    warning("reported SNP count for ", name, " in package ",
                            SNPlocs_pkgname(x), " does not match the ",
                            "number of SNPs returned by ", SNPlocs_pkgname(x),
                            ":::getSNPlocs()")
                .inplaceReplaceLetterAtLoc(datacache_env[[name]], snps$loc, snps$alleles_as_ambig)
            }

            ## Load and put the (inactive) built-in masks, if any
            objname <- paste("masks.", name, sep="")
            filename <- system.file("data", paste(objname, ".rda", sep=""), package=x@package)
            if (file.exists(filename)) {
                load(filename)
                masks(datacache_env[[name]]) <- get(objname)
                active(masks(datacache_env[[name]])) <- FALSE
                remove(list=objname)
            }

            datacache_env[[name]]
        }
        makeActiveBinding(name, getter, activebindings_env)
    }

    for (name in names) {
        file <- system.file(x@subdir, paste(name, ".rda", sep=""), package=x@package)
        addCachedItem(name, file)
    }
}

BSgenome <- function(organism, species, provider, provider_version,
                     release_date, release_name, source_url,
                     seqnames, mseqnames, package, subdir)
{
    ans <- new("BSgenome",
        organism=organism,
        species=species,
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        source_url=source_url,
        seqnames=seqnames,
        mseqnames=mseqnames,
        package=package,
        subdir=subdir,
        .activebindings_env=new.env(parent=emptyenv()),
        .datacache_env=new.env(parent=emptyenv())
    )
    .assignDataToNames(ans)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The injectSNPs() function.
###

injectSNPs <- function(bsgenome, SNPlocs_pkgname)
{
    if (!is(bsgenome, "BSgenome"))
        stop("'bsgenome' must be a BSgenome object")
    if (!is.null(SNPlocs_pkgname(bsgenome)))
        stop("SNPs already injected in genome, injecting from more than 1 package is not supported")
    if (!is.character(SNPlocs_pkgname) || length(SNPlocs_pkgname) != 1 || is.na(SNPlocs_pkgname))
        stop("'SNPlocs_pkgname' must be a single string")
    library(SNPlocs_pkgname, character.only=TRUE)
    getSNPcount <- get("getSNPcount",
                      envir=as.environment(paste("package", SNPlocs_pkgname, sep=":")),
                      inherits=FALSE)
    if (!all(names(getSNPcount()) %in% seqnames(bsgenome)))
        stop("seqnames in package ", SNPlocs_pkgname, " are not compatible ",
             "with the seqnames of this BSgenome object")
    getSNPlocs <- get("getSNPlocs",
                      envir=as.environment(paste("package", SNPlocs_pkgname, sep=":")),
                      inherits=FALSE)
    bsgenome@SNPlocs_pkgname <- SNPlocs_pkgname
    bsgenome@.activebindings_env <- new.env(parent=emptyenv())
    bsgenome@.datacache_env <- new.env(parent=emptyenv())
    .assignDataToNames(bsgenome, getSNPcount, getSNPlocs)
    bsgenome
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

.SHOW_PREFIX <- "| "
.SHOW_SEQSECTION_PREFIX <- "|   "

setMethod("show", "BSgenome",
    function(object)
    {
        mystrwrap <- function(line)
            writeLines(strwrap(line, width=getOption("width")+1,
                               exdent=0, prefix=.SHOW_PREFIX))
        showSequenceIndex <- function(names, prefix)
        {
            index_width <- getOption("width") + 2 - nchar(prefix)
            col_width <- max(nchar(names))
            ncols <- index_width %/% (col_width + 2)
            col <- 1
            for (name in names) {
                if (col == 1) cat(prefix)
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
        cat(object@species, "genome\n")
        cat(.SHOW_PREFIX, "\n", sep="")
        cat(.SHOW_PREFIX, "organism: ", object@organism, "\n", sep="")
        cat(.SHOW_PREFIX, "provider: ", object@provider, "\n", sep="")
        cat(.SHOW_PREFIX, "provider version: ", object@provider_version, "\n", sep="")
        cat(.SHOW_PREFIX, "release date: ", object@release_date, "\n", sep="")
        cat(.SHOW_PREFIX, "release name: ", object@release_name, "\n", sep="")
        if (!is.null(SNPlocs_pkgname(object)))
            cat(.SHOW_PREFIX, "with SNPs injected from package: ", SNPlocs_pkgname(object), "\n", sep="")
        cat(.SHOW_PREFIX, "\n", sep="")
        if (length(mseqnames(object)) != 0)
            mystrwrap("single sequences (see '?seqnames'):")
        else
            mystrwrap("sequences (see '?seqnames'):")
        if (length(seqnames(object)) != 0)
            showSequenceIndex(seqnames(object), .SHOW_SEQSECTION_PREFIX)
        else
            cat(.SHOW_SEQSECTION_PREFIX, "NONE\n", sep="")
        cat(.SHOW_PREFIX, "\n", sep="")
        if (length(mseqnames(object)) != 0) {
            mystrwrap("multiple sequences (see '?mseqnames'):")
            showSequenceIndex(mseqnames(object), .SHOW_SEQSECTION_PREFIX)
            cat(.SHOW_PREFIX, "\n", sep="")
        }
        mystrwrap("(use the '$' or '[[' operator to access a given sequence)")
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
        get(name, envir=x@.activebindings_env)
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
            what <- ls(x@.datacache_env) ## everything
        remove(list=what, envir=x@.datacache_env)
    }
)


