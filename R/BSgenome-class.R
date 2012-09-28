### =========================================================================
### The "BSgenome" class
### -------------------------------------------------------------------------

setClass("BSgenome",
    contains="GenomeDescription",
    representation(
        ## source_url: permanent URL to the place where the FASTA files used
        ## to produce the sequences below can be found (and downloaded)
        source_url="character",

        ## mseqnames: names of "multiple" sequences (e.g. upstream)
        mseqnames="character",

        ## where to find the serialized objects containing the sequences
        seqs_pkgname="character",
        seqs_dirpath="character",

        ## where to find the serialized objects containing the masks
        nmask_per_seq="integer",
        masks_pkgname="character",
        masks_dirpath="character",

        ## for SNPs injection
        injectSNPs_handler="InjectSNPsHandler",

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
    if (length(loaded_names) != 1L)
        stop("file '", filepath, "' contains 0 or more than 1 serialized object. ",
             "May be the ", objpkgname, " package is corrupted?")
    ## ... and that it has the expected name
    if (loaded_names != objname)
        stop("the serialized object in file '", filepath, "' ",
             "doesn't have the expected name. ",
             "May be the ", objpkgname, " package is corrupted?")
    ans <- get(objname, envir=tmp_env)
    ## TODO: This is temporary code to make temporarily broken BSgenome
    ## data packages work despite the big internal renaming I made in
    ## IRanges >= 1.3.76. Remove it after all the BSgenome data packages
    ## have been reforged.
    if (is(ans, "XString") || is(ans, "XStringSet"))
        ans <- updateObject(ans)
    return(ans)
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

setReplaceMethod(".linkToCachedObject", "SharedVector",
    function(x, value)
    {
        x@.link_to_cached_object <- value
        x
    }
)

setReplaceMethod(".linkToCachedObject", "XString",
    function(x, value)
    {
        .linkToCachedObject(x@shared) <- value
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
        ### This means that an XStringSet object with a "pool" slot of
        ### length != 1 will be permanently cached.
        if (length(x@pool) == 1L)
            x@pool@.link_to_cached_object_list[[1L]] <- value
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

setGeneric("sourceUrl", function(x) standardGeneric("sourceUrl"))
setMethod("sourceUrl", "BSgenome", function(x) x@source_url)

setMethod("SNPlocs_pkgname", "BSgenome",
    function(x) SNPlocs_pkgname(x@injectSNPs_handler)
)

setMethod("SNPcount", "BSgenome",
    function(x) SNPcount(x@injectSNPs_handler)
)

setMethod("SNPlocs", "BSgenome",
    function(x, seqname) SNPlocs(x@injectSNPs_handler, seqname)
)

setGeneric("mseqnames", function(x) standardGeneric("mseqnames"))
setMethod("mseqnames", "BSgenome",
    function(x) { if (length(x@mseqnames) == 0L) NULL else x@mseqnames }
)

setMethod("names", "BSgenome", function(x) c(seqnames(x), mseqnames(x)))

setGeneric("masknames", function(x) standardGeneric("masknames"))
setMethod("masknames", "BSgenome",
    function(x)
    {
        if (x@nmask_per_seq == 0L)
            return(NULL)
        ## TODO: Put this kind of checking in a validity method for BSgenome
        ## objects (that's what validity methods are for).
        if (x@nmask_per_seq > length(BUILTIN_MASKNAMES))
            stop("internal anomaly: x@nmask_per_seq > ",
                 length(BUILTIN_MASKNAMES))
        BUILTIN_MASKNAMES[seq_len(x@nmask_per_seq)]
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor-like functions and generics
###

.makeSeqinfo <- function(seqnames, circ_seqs, seqs_pkgname, seqs_dirpath,
                         provider_version)
{
    objname <- "seqlengths"
    seqlengths <- .loadSingleObject(objname, seqs_dirpath, seqs_pkgname)
    if (!identical(names(seqlengths), seqnames)) {
        filepath <- .getObjFilepath(objname, seqs_dirpath)
        stop("sequence names found in file '", filepath, "' are not ",
             "identical to 'seqnames'. ",
             "May be the ", seqs_pkgname, " package is corrupted?")
    }
    if (identical(circ_seqs, NA))
        is_circ <- NA
    else
        is_circ <- seqnames %in% circ_seqs
    Seqinfo(seqnames=seqnames, seqlengths=seqlengths, isCircular=is_circ,
            genome=provider_version)
}

BSgenome <- function(organism, species, provider, provider_version,
                     release_date, release_name, source_url,
                     seqnames, circ_seqs=NA, mseqnames,
                     seqs_pkgname, seqs_dirpath,
                     nmask_per_seq, masks_pkgname, masks_dirpath)
{
    seqinfo <- .makeSeqinfo(seqnames, circ_seqs, seqs_pkgname, seqs_dirpath,
                            provider_version)
    if (is.null(mseqnames))
        mseqnames <- character(0)
    new("BSgenome",
        GenomeDescription(organism, species,
                          provider, provider_version,
                          release_date, release_name,
                          seqinfo),
        source_url=source_url,
        mseqnames=mseqnames,
        seqs_pkgname=seqs_pkgname,
        seqs_dirpath=seqs_dirpath,
        nmask_per_seq=as.integer(nmask_per_seq),
        masks_pkgname=masks_pkgname,
        masks_dirpath=masks_dirpath,
        .seqs_cache=new.env(parent=emptyenv()),
        .link_counts=new.env(parent=emptyenv())
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

.SHOW_BSGENOME_PREFIX <- "| "
.SHOW_SEQSECTION_PREFIX <- "|   "

setMethod("show", "BSgenome",
    function(object)
    {
        mystrwrap <- function(line)
            writeLines(strwrap(line, width=getOption("width")+1,
                               exdent=0L, prefix=.SHOW_BSGENOME_PREFIX))
        showSequenceIndex <- function(names, prefix)
        {
            index_width <- getOption("width") + 2L - nchar(prefix)
            col_width <- max(nchar(names))
            ncols <- index_width %/% (col_width + 2L)
            col <- 1L
            for (name in names) {
                if (col == 1L) cat(prefix)
                cat(format(name, width=col_width))
                if (col == ncols) {
                    cat("\n")
                    col <- 1L
                } else {
                    cat("  ")
                    col <- col + 1L
                }
            }
            if (col != 1L) cat("\n")
        }
        if (!is.na(object@species)) {
            cat(object@species, "genome\n")
            cat(.SHOW_BSGENOME_PREFIX, "\n", sep="")
        }
        showGenomeDescription(object, margin=.SHOW_BSGENOME_PREFIX)
        if (!is.null(SNPlocs_pkgname(object)))
            cat(.SHOW_BSGENOME_PREFIX, "with SNPs injected from package: ", SNPlocs_pkgname(object), "\n", sep="")
        cat(.SHOW_BSGENOME_PREFIX, "\n", sep="")
        if (length(mseqnames(object)) != 0L)
            mystrwrap("single sequences (see '?seqnames'):")
        else
            mystrwrap("sequences (see '?seqnames'):")
        if (length(seqnames(object)) != 0L)
            showSequenceIndex(seqnames(object), .SHOW_SEQSECTION_PREFIX)
        else
            cat(.SHOW_SEQSECTION_PREFIX, "NONE\n", sep="")
        cat(.SHOW_BSGENOME_PREFIX, "\n", sep="")
        if (length(mseqnames(object)) != 0L) {
            mystrwrap("multiple sequences (see '?mseqnames'):")
            showSequenceIndex(mseqnames(object), .SHOW_SEQSECTION_PREFIX)
            cat(.SHOW_BSGENOME_PREFIX, "\n", sep="")
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
    seqs_dirpath <- bsgenome@seqs_dirpath
    nmask_per_seq <- length(masknames(bsgenome))
    masks_pkgname <- bsgenome@masks_pkgname
    masks_dirpath <- bsgenome@masks_dirpath
    ans <- .loadSingleObject(name, seqs_dirpath, seqs_pkgname)
    if (!is(ans, "XString"))
        return(ans)
    ## Check the length of the sequence
    if (length(ans) != seqlengths(bsgenome)[[name]]) {
        seq_filepath <- .getObjFilepath(name, seqs_dirpath)
        stop("sequence found in file '", seq_filepath, "' does ",
             "not have the length reported by seqlengths(). ",
             "May be the ", seqs_pkgname, " package is corrupted?")
    }
    ## Inject the SNPs, if any
    snps <- SNPlocs(bsgenome, name)
    if (!is.null(snps)) 
        .inplaceReplaceLetterAt(ans, snps$loc, snps$alleles_as_ambig)
    ## Load and set the built-in masks, if any
    if (nmask_per_seq > 0L) {
        objname <- paste(name, ".masks", sep="")
        builtinmasks <- .loadSingleObject(objname, masks_dirpath, masks_pkgname)
        if (length(builtinmasks) < nmask_per_seq) {
            masks_filepath <- .getObjFilepath(objname, masks_dirpath)
            stop("expecting ", nmask_per_seq, " built-in masks per ",
                 "single sequence, found only ", length(builtinmasks),
                 " in file '", masks_filepath, "'. ",
                 "May be the ", masks_pkgname, " package is corrupted?")
        }
        if (length(builtinmasks) > nmask_per_seq)
            builtinmasks <- builtinmasks[seq_len(nmask_per_seq)]
        if (!identical(names(builtinmasks), masknames(bsgenome))) {
            masks_filepath <- .getObjFilepath(objname, masks_dirpath)
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
    ## Using the 'if (!exists()) assign(); get()' approach is NOT 100%
    ## reliable:
    ##
    ##   if (!exists(name, envir=seqs_cache, inherits=FALSE)) {
    ##       ...
    ##       assign(name, ans, envir=seqs_cache)
    ##   }
    ##   get(name, envir=seqs_cache, inherits=FALSE)
    ##
    ## because the symbol (name) can disappear from the cache between the
    ## moment we test for its presence and the moment we try to get it.
    ## It's not me being paranoid, we've seen this happen! One possible
    ## explanation for this is that the symbol was candidate for removal
    ## from the cache but that removal didn't happen yet because gc() had
    ## not yet been called (removal from the cache is implemented thru the
    ## finalizers registered on the objects that are copied from the cache
    ## and made available to the user). Then the call to get() would trigger
    ## garbbage collection and that in turn would trigger the removal of
    ## the symbol *before* get() had a chance to get to it. So now we use
    ## the 'try(get(...))' approach, which hopefully is 100% reliable!
    ans <- try(get(name, envir=seqs_cache, inherits=FALSE), silent=TRUE)
    if (is(ans, "try-error")) {
        ans <- .loadBSgenomeSequence(name, bsgenome)
        if (getOption("verbose"))
            cat("caching ", name, "\n", sep="")
        assign(name, ans, envir=seqs_cache)
    }
    .linkToCachedObject(ans) <- .newLinkToCachedObject(
                                    name,
                                    seqs_cache,
                                    bsgenome@.link_counts)
    ans
}

setMethod("[[", "BSgenome",
    function(x, i, j, ...)
    {
        ## 'x' is guaranteed to be a "BSgenome" object (if it's not, then the
        ## method dispatch algo will not even call this method), so nargs() is
        ## guaranteed to be >= 1
        if (nargs() >= 3L)
            stop("too many subscripts")
        subscripts <- list(...)
        if (!missing(i))
            subscripts$i <- i
        if (!missing(j))
            subscripts$j <- j
        ## At this point, 'subscripts' should be guaranteed
        ## to be of length <= 1
        if (length(subscripts) == 0L)
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
            if (i < 1L || length(x) < i)
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
