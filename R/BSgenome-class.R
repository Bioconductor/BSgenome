### =========================================================================
### The "BSgenome" class
### -------------------------------------------------------------------------

setClass("BSgenome",
    contains="GenomeDescription",
    representation(
        ## source_url: permanent URL to the place where the FASTA files used
        ## to produce the sequences below can be found (and downloaded)
        source_url="character",

        ## named vector representing the translation table from the original
        ## seqnames (as stored in self@seqinfo@seqnames, the 'seqinfo' slot
        ## being inherited from the GenomeDescription class) to the user
        ## seqnames
        user_seqnames="character",

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
### Accessor methods.
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
### seqinfo() accessor and related.
###

setMethod("seqinfo", "BSgenome",
    function(x)
    {
        ans <- x@seqinfo
        seqlevels(ans) <- x@user_seqnames
        ans
    }
)

### This is a restricted "seqinfo<-" method for BSgenome objects that
### only supports replacement of the sequence names, i.e., except for their
### sequence names, Seqinfo objects 'value' and 'seqinfo(x)' must be identical.
setReplaceMethod("seqinfo", "BSgenome",
    function(x, new2old=NULL, force=FALSE, value)
    {
        if (!is(value, "Seqinfo"))
            stop("the supplied 'seqinfo' must be a Seqinfo object")
        IN_THIS_CONTEXT <- paste0("when replacing the 'seqinfo' ",
                                  "of a BSgenome object")
        if (!identical(force, FALSE))
            stop("'force' not supported ", IN_THIS_CONTEXT)
        x_seqinfo <- seqinfo(x)
        if (is.null(new2old)) {
            ## Support no-op seqinfo(x) <- seqinfo(x).
            if (!identical(value, x_seqinfo))
                stop("'new2old' must be specified ", IN_THIS_CONTEXT)
            return(x)
        }
        if (length(value) != length(x_seqinfo))
            stop("the supplied 'seqinfo' must have the same length ",
                 "as the current 'seqinfo' ", IN_THIS_CONTEXT)
        if (!identical(new2old, seq_along(value)))
            stop("'new2old' must be NULL or equal to 'seq_along(value)' ",
                 IN_THIS_CONTEXT)
        new_seqnames <- seqnames(value)
        seqnames(x_seqinfo) <- new_seqnames
        if (!identical(value, x_seqinfo))
            stop("the supplied and current 'seqinfo' can differ only ",
                 "in their sequence names ", IN_THIS_CONTEXT)
        if (any(new_seqnames %in% mseqnames(x)))
            stop("the supplied 'seqnames' cannot match any of the ",
                 "multiple sequence names (as returned by 'mseqnames(x)')")
        x@user_seqnames[] <- new_seqnames  # using [] to preserve the names
        x
    }
)

setReplaceMethod("seqnames", "BSgenome",
    function(x, value)
    {
        x_seqinfo <- seqinfo(x)
        seqnames(x_seqinfo) <- value
        seqinfo(x, new2old=seq_along(x_seqinfo)) <- x_seqinfo
        x
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
    user_seqnames <- seqnames(seqinfo)
    names(user_seqnames) <- user_seqnames
    if (is.null(mseqnames))
        mseqnames <- character(0)
    new("BSgenome",
        GenomeDescription(organism, species,
                          provider, provider_version,
                          release_date, release_name,
                          seqinfo),
        source_url=source_url,
        user_seqnames=user_seqnames,
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

.loadBSgenomeSequence <- function(objname, bsgenome, user_seqname)
{
    seqs_pkgname <- bsgenome@seqs_pkgname
    seqs_dirpath <- bsgenome@seqs_dirpath
    nmask_per_seq <- length(masknames(bsgenome))
    masks_pkgname <- bsgenome@masks_pkgname
    masks_dirpath <- bsgenome@masks_dirpath
    ans <- .loadSingleObject(objname, seqs_dirpath, seqs_pkgname)
    if (!is(ans, "XString"))
        return(ans)
    ## Check the length of the sequence
    if (length(ans) != seqlengths(bsgenome)[[user_seqname]]) {
        seq_filepath <- .getObjFilepath(objname, seqs_dirpath)
        stop("sequence found in file '", seq_filepath, "' does ",
             "not have the length reported by seqlengths(). ",
             "May be the ", seqs_pkgname, " package is corrupted?")
    }
    ## Inject the SNPs, if any
    snps <- SNPlocs(bsgenome, objname)
    if (!is.null(snps)) 
        .inplaceReplaceLetterAt(ans, snps$loc, snps$alleles_as_ambig)
    ## Load and set the built-in masks, if any
    if (nmask_per_seq > 0L) {
        masks_objname <- paste(objname, ".masks", sep="")
        builtinmasks <- .loadSingleObject(masks_objname, masks_dirpath,
                                          masks_pkgname)
        if (length(builtinmasks) < nmask_per_seq) {
            masks_filepath <- .getObjFilepath(masks_objname, masks_dirpath)
            stop("expecting ", nmask_per_seq, " built-in masks per ",
                 "single sequence, found only ", length(builtinmasks),
                 " in file '", masks_filepath, "'. ",
                 "May be the ", masks_pkgname, " package is corrupted?")
        }
        if (length(builtinmasks) > nmask_per_seq)
            builtinmasks <- builtinmasks[seq_len(nmask_per_seq)]
        if (!identical(names(builtinmasks), masknames(bsgenome))) {
            masks_filepath <- .getObjFilepath(masks_objname, masks_dirpath)
            stop("mask names found in file '", masks_filepath, "' are not ",
                 "identical to the names returned by masknames(). ",
                 "May be the ", masks_pkgname, " package is corrupted?")
        }
        masks(ans) <- builtinmasks
    }
    ans
}

.getBSgenomeSequence <- function(objname, bsgenome, user_seqname)
{
    seqs_cache <- bsgenome@.seqs_cache
    ## Using the 'if (!exists()) assign(); get()' approach is NOT 100%
    ## reliable:
    ##
    ##   if (!exists(objname, envir=seqs_cache, inherits=FALSE)) {
    ##       ...
    ##       assign(objname, ans, envir=seqs_cache)
    ##   }
    ##   get(objname, envir=seqs_cache, inherits=FALSE)
    ##
    ## because the symbol (objname) can disappear from the cache between
    ## the moment we test for its presence and the moment we try to get it.
    ## It's not me being paranoid, we've seen this happen! One possible
    ## explanation for this is that the symbol was candidate for removal
    ## from the cache but that removal didn't happen yet because gc() had
    ## not yet been called (removal from the cache is implemented thru the
    ## finalizers registered on the objects that are copied from the cache
    ## and made available to the user). Then the call to get() would trigger
    ## garbbage collection and that in turn would trigger the removal of
    ## the symbol *before* get() had a chance to get to it. So now we use
    ## the 'try(get(...))' approach, which hopefully is 100% reliable!
    ans <- try(get(objname, envir=seqs_cache, inherits=FALSE), silent=TRUE)
    if (is(ans, "try-error")) {
        ans <- .loadBSgenomeSequence(objname, bsgenome, user_seqname)
        if (getOption("verbose"))
            cat("caching ", objname, "\n", sep="")
        assign(objname, ans, envir=seqs_cache)
    }
    .linkToCachedObject(ans) <- .newLinkToCachedObject(
                                    objname,
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
        ## Translate user-supplied sequence named back to original name.
        idx <- match(name, x@user_seqnames)
        if (is.na(idx)) {
            objname <- name
        } else {
            objname <- names(x@user_seqnames)[[idx]]
        }
        .getBSgenomeSequence(objname, x, name)
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
