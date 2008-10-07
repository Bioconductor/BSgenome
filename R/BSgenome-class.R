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

        ## release_name: "NCBI Build 36.1", "NCBI Build 36",
        ## "SGD 1 Oct 2003 sequence", etc...
        release_name="character",

        ## source_url: permanent URL to the place where the FASTA files used
        ## to produce the sequences below can be found (and downloaded)
        source_url="character",

        ## for SNPs injection
        SNPlocs_pkgname="character",

        ## seqnames: names of "single" sequences (e.g. chromosomes)
        seqnames="character",

        ## seqlengths: lengths of "single" sequences
        seqlengths="integer",

        ## mseqnames: names of "multiple" sequences (e.g. upstream)
        mseqnames="character",

        ## where to find the serialized objects containing the sequences
        seqs_pkgname="character",
        seqs_dir="character",

        ## where to find the serialized objects containing the masks
        nmask_per_seq="integer",
        masks_pkgname="character",
        masks_dir="character",

        .seqs_cache="environment",
        .link_counts="environment"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helper functions used for delayed-loading/caching/unloading the
### sequences.
###

.getObjFilepath <- function(objname, objdir)
{
    ## Should never happen.
    if (objdir == "")
        ## TODO: Put this kind of checking in a validity method for BSgenome
        ## objects (that's what validity methods are for).
        stop("internal anomaly: objdir is \"\"")
    filename <- paste(objname, ".rda", sep="")
    filepath <- file.path(objdir, filename)
    if (!file.exists(filepath))
        stop("file '", filepath, "' doesn't exist")
    filepath
}

.loadSingleObject <- function(objname, objdir, objpkgname)
{
    filepath <- .getObjFilepath(objname, objdir)
    tmp_env <- new.env(parent=emptyenv())
    loaded_names <- load(filepath, envir=tmp_env)
    ## Check that we get only 1 object
    if (length(loaded_names) != 1)
        stop("file '", filepath, "' contains 0 or more than 1 serialized object. ",
             "May be the ", objpkgname, " package is corrupted?")
    ## ... and that it has the expected name
    if (loaded_names != objname)
        stop("the serialized object in file '", filepath, "' ",
             "doesn't have the expected name. ",
             "May be the ", objpkgname, " package is corrupted?")
    get(objname, envir=tmp_env)
}

### Return a new link to a cached object.
### 'objname' is the name of the cached object.
### 'cache' is the caching environment.
### 'link_counts' is the environment where we keep track of the number of links
### for each cached object. When the number of links for a given cached object
### reaches 0, then it is removed from the cache.
.newLinkToCachedObject <- function(objname, cache, link_counts)
{
    ans <- new.env(parent=emptyenv())
    if (exists(objname, envir=link_counts, inherits=FALSE))
        link_count0 <- get(objname, envir=link_counts, inherits=FALSE) + 1L
    else
        link_count0 <- 1L
    reg.finalizer(ans,
        function(e)
        {
            link_count <- get(objname, envir=link_counts, inherits=FALSE) - 1L
            assign(objname, link_count, envir=link_counts)
            if (link_count == 0) {
                if (getOption("verbose"))
                    cat("uncaching ", objname, "\n", sep="")
                remove(list=objname, envir=cache)
            }
        }
    )
    assign(objname, link_count0, envir=link_counts)
    ans
}

setGeneric(".linkToCachedObject<-", signature="x",
    function(x, value) standardGeneric(".linkToCachedObject<-")
)

setReplaceMethod(".linkToCachedObject", "SequencePtr",
    function(x, value)
    {
        x@.link_to_cached_object <- value
        x
    }
)

setReplaceMethod(".linkToCachedObject", "XString",
    function(x, value)
    {
        .linkToCachedObject(x@xdata) <- value
        x
    }
)

setReplaceMethod(".linkToCachedObject", "MaskedXString",
    function(x, value)
    {
        .linkToCachedObject(x@unmasked) <- value
        x
    }
)

setReplaceMethod(".linkToCachedObject", "XStringSet",
    function(x, value)
    {
        .linkToCachedObject(x@super) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The names of the built-in masks.
###

BUILTIN_MASKNAMES <- c("AGAPS", "AMB", "RM", "TRF")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "length" and accessor methods.
###

setMethod("length", "BSgenome", function(x) length(names(x)))

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
    function(x)
    {
        if (length(x@SNPlocs_pkgname) == 0)
            return(NULL)
        library(x@SNPlocs_pkgname, character.only=TRUE)
        x@SNPlocs_pkgname
    }
)

setGeneric("SNPcount", function(x) standardGeneric("SNPcount"))
setMethod("SNPcount", "BSgenome",
    function(x)
    {
        pkg <- SNPlocs_pkgname(x)
        if (is.null(pkg))
            return(NULL)
        getSNPcount <- try(get("getSNPcount",
                               envir=as.environment(paste("package", pkg, sep=":")),
                               inherits=FALSE), silent=TRUE)
        if (!is.function(getSNPcount))
            stop("cannot use package ", SNPlocs_pkgname(x), " for SNP injection: ",
                 "it doesn't seem to define (and export) a function called ",
                 "'getSNPcount'")
        getSNPcount()
    }
)

setGeneric("SNPlocs", signature="x", function(x, seqname) standardGeneric("SNPlocs"))
setMethod("SNPlocs", "BSgenome",
    function(x, seqname)
    {
        pkg <- SNPlocs_pkgname(x)
        if (is.null(pkg))
            return(NULL)
        getSNPlocs <- try(get("getSNPlocs",
                              envir=as.environment(paste("package", pkg, sep=":")),
                              inherits=FALSE), silent=TRUE)
        if (!is.function(getSNPlocs))
            stop("cannot use package ", SNPlocs_pkgname(x), " for SNP injection: ",
                 "it doesn't seem to define (and export) a function called ",
                 "'getSNPlocs'")
        getSNPlocs(seqname)
    }
)

setGeneric("seqnames", function(x) standardGeneric("seqnames"))
setMethod("seqnames", "BSgenome",
    function(x) { if (length(x@seqnames) == 0) NULL else x@seqnames }
)

setGeneric("seqlengths", function(x) standardGeneric("seqlengths"))
setMethod("seqlengths", "BSgenome",
    function(x)
    {
        if (length(x@seqlengths) == 1 && is.na(x@seqlengths)) {
            objname <- "seqlengths"
            x@seqlengths <- .loadSingleObject(objname, x@seqs_dir, x@seqs_pkgname)
            if (!identical(names(x@seqlengths), x@seqnames)) {
                filepath <- .getObjFilepath(objname, x@seqs_dir)
                stop("sequence names found in file '", filepath, "' are not ",
                     "identical to the names returned by seqnames(). ",
                     "May be the ", x@seqs_pkgname, " package is corrupted?")
            }
        }
        x@seqlengths
    }
)

setGeneric("mseqnames", function(x) standardGeneric("mseqnames"))
setMethod("mseqnames", "BSgenome",
    function(x) { if (length(x@mseqnames) == 0) NULL else x@mseqnames }
)

setMethod("names", "BSgenome", function(x) c(seqnames(x), mseqnames(x)))

setGeneric("masknames", function(x) standardGeneric("masknames"))
setMethod("masknames", "BSgenome",
    function(x)
    {
        if (x@nmask_per_seq == 0)
            return(NULL)
        ## TODO: Put this kind of checking in a validity method for BSgenome
        ## objects (that's what validity methods are for).
        if (x@nmask_per_seq > length(BUILTIN_MASKNAMES))
            stop("internal anomaly: x@nmask_per_seq > ", length(BUILTIN_MASKNAMES))
        BUILTIN_MASKNAMES[seq_len(x@nmask_per_seq)]
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor-like functions and generics
###

BSgenome <- function(organism, species, provider, provider_version,
                     release_date, release_name, source_url,
                     seqnames, mseqnames, seqs_pkgname, seqs_dir,
                     nmask_per_seq, masks_pkgname, masks_dir)
{
    if (is.null(seqnames))
        seqnames <- character(0)
    if (is.null(mseqnames))
        mseqnames <- character(0)
    ans <- new("BSgenome",
        organism=organism,
        species=species,
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        source_url=source_url,
        seqnames=seqnames,
        seqlengths=as.integer(NA),
        mseqnames=mseqnames,
        seqs_pkgname=seqs_pkgname,
        seqs_dir=seqs_dir,
        nmask_per_seq=as.integer(nmask_per_seq),
        masks_pkgname=masks_pkgname,
        masks_dir=masks_dir,
        .seqs_cache=new.env(parent=emptyenv()),
        .link_counts=new.env(parent=emptyenv())
    )
    ans
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
### Subsetting.
###

.loadBSgenomeSequence <- function(name, bsgenome)
{
    seqs_pkgname <- bsgenome@seqs_pkgname
    seqs_dir <- bsgenome@seqs_dir
    nmask_per_seq <- length(masknames(bsgenome))
    masks_pkgname <- bsgenome@masks_pkgname
    masks_dir <- bsgenome@masks_dir
    ans <- .loadSingleObject(name, seqs_dir, seqs_pkgname)
    if (!is(ans, "XString"))
        return(ans)
    ## Check the length of the sequence
    if (length(ans) != seqlengths(bsgenome)[[name]]) {
        seq_filepath <- .getObjFilepath(name, seqs_dir)
        stop("sequence found in file '", seq_filepath, "' does ",
             "not have the length reported by seqlengths(). ",
             "May be the ", seqs_pkgname, " package is corrupted?")
    }
    ## Inject the SNPs, if any
    if (!is.null(SNPlocs_pkgname(bsgenome))) {
        snp_count <- SNPcount(bsgenome)
        if (name %in% names(snp_count)) {
            snps <- SNPlocs(bsgenome, name)
            if (nrow(snps) != snp_count[name])
                warning("reported SNP count for sequence ", name, " in package ",
                        SNPlocs_pkgname(bsgenome), " does not match the ",
                        "number of SNPs returned by ", SNPlocs_pkgname(bsgenome),
                        ":::getSNPlocs()")
            .inplaceReplaceLetterAt(ans, snps$loc, snps$alleles_as_ambig)
        } 
    }
    ## Load and set the built-in masks, if any
    if (nmask_per_seq > 0) {
        objname <- paste(name, ".masks", sep="")
        builtinmasks <- .loadSingleObject(objname, masks_dir, masks_pkgname)
        if (length(builtinmasks) < nmask_per_seq) {
            masks_filepath <- .getObjFilepath(objname, masks_dir)
            stop("expecting ", nmask_per_seq, " built-in masks per ",
                 "single sequence, found only ", length(builtinmasks),
                 " in file '", masks_filepath, "'. ",
                 "May be the ", masks_pkgname, " package is corrupted?")
        }
        if (length(builtinmasks) > nmask_per_seq)
            builtinmasks <- builtinmasks[seq_len(nmask_per_seq)]
        if (!identical(names(builtinmasks), masknames(bsgenome))) {
            masks_filepath <- .getObjFilepath(objname, masks_dir)
            stop("mask names found in file '", masks_filepath, "' are not ",
                 "identical to the names returned by masknames(). ",
                 "May be the ", masks_pkgname, " package is corrupted?")
        }
        masks(ans) <- builtinmasks
    }
    ans
}

.getBSgenomeSequence <- function(name, bsgenome)
{
    seqs_cache <- bsgenome@.seqs_cache
    if (!exists(name, envir=seqs_cache, inherits=FALSE)) {
        ans <- .loadBSgenomeSequence(name, bsgenome)
        assign(name, ans, envir=seqs_cache)
    }
    ans <- get(name, envir=seqs_cache)
    .linkToCachedObject(ans) <- .newLinkToCachedObject(name, seqs_cache,
                                                       bsgenome@.link_counts)
    ans
}

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
        .getBSgenomeSequence(name, x)
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

