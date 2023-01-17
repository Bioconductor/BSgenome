### =========================================================================
### BSgenome objects
### -------------------------------------------------------------------------


setClass("BSgenome",
    ## Expected metadata data are:
    ## - organism
    ## - common_name
    ## - provider
    ## - genome
    ## - release_date
    ## - source_url: permanent URL to the place where the 2bit and/or FASTA
    ##   files used to produce the "single" and "multiple" sequences can be
    ##   found (and downloaded).
    contains="Annotated",
    representation(
        ## Name of the BSgenome data package where the BSgenome object is
        ## defined.
        pkgname="character",

        ## On-disk "single" and "multiple" sequences.
        single_sequences="OnDiskNamedSequences",
        multiple_sequences="RdaCollection",

        ## Original seqinfo.
        seqinfo="Seqinfo",

        ## Named vector representing the translation table from the original
        ## seqnames (as stored in self@seqinfo@seqnames) to the user seqnames.
        user_seqnames="character",

        ## For SNPs injection.
        injectSNPs_handler="InjectSNPsHandler",

        ## Used for caching the single and multiple on-disk sequences.
        .seqs_cache="environment",
        .link_counts="environment"
    )
)

setClass("MaskedBSgenome",
    contains="BSgenome",
    representation(
        ## Name of the BSgenome data package where the MaskedBSgenome object
        ## is defined.
        masks_pkgname="character",

        ## Only the single sequences can be masked. They should all have the
        ## same set of masks.
        nmask_per_seq="integer",
        masks="RdaCollection"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helper functions used for delayed-loading/caching/unloading the
### sequences.
###

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
### Accessors
###

### Metadata.
setMethod("organism", "BSgenome",
    function(object) metadata(object)$organism
)

setMethod("commonName", "BSgenome",
    function(object) metadata(object)$common_name
)

setMethod("provider", "BSgenome", function(x) metadata(x)$provider)

setMethod("providerVersion", "BSgenome",
    function(x)
    {
        msg <- c("Using providerVersion() on a ", class(x), " object ",
                 "is deprecated. Please use 'metadata(x)$genome' instead.")
        .Deprecated(msg=c("  ", wmsg(msg)))
        metadata(x)$genome
    }
)

setMethod("releaseDate", "BSgenome", function(x) metadata(x)$release_date)

setGeneric("sourceUrl", function(x) standardGeneric("sourceUrl"))
setMethod("sourceUrl", "BSgenome", function(x) metadata(x)$source_url)

setMethod("length", "BSgenome", function(x) length(names(x)))

setGeneric("mseqnames", function(x) standardGeneric("mseqnames"))
setMethod("mseqnames", "BSgenome",
    function(x)
    {
        ans <- names(x@multiple_sequences)
        if (length(ans) == 0L)
            ans <- NULL
        ans
    }
)

setMethod("names", "BSgenome", function(x) c(seqnames(x), mseqnames(x)))

setGeneric("masknames", function(x) standardGeneric("masknames"))

setMethod("masknames", "BSgenome", function(x) NULL)

setMethod("masknames", "MaskedBSgenome",
    function(x)
    {
        ## TODO: Put this kind of checking in a validity method for
        ## MaskedBSgenome objects (that's what validity methods are for).
        if (x@nmask_per_seq > length(BUILTIN_MASKNAMES))
            stop("internal anomaly: x@nmask_per_seq > ",
                 length(BUILTIN_MASKNAMES))
        BUILTIN_MASKNAMES[seq_len(x@nmask_per_seq)]
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo() accessor and related.
###

setMethod("seqnames", "BSgenome", function(x) unname(x@user_seqnames))

setMethod("seqinfo", "BSgenome",
    function(x)
    {
        ans <- x@seqinfo
        seqlevels(ans) <- seqnames(x)
        ans
    }
)

### We implement a restricted seqinfo() setter for BSgenome object 'x' that
### supports altering **only** the seqlevels and/or genome of 'seqinfo(x)'.
### It does NOT allow subsetting 'seqinfo(x)' (by dropping/reordering some
### of its seqlevels), or altering its seqlengths or circularity flags!
### In other words, except for their seqnames() and genome(), Seqinfo
### objects 'new_seqinfo' and 'old_seqinfo' must be identical. This is all
### we need to make the seqlevelsStyle() setter work on a BSgenome object.
.check_new2old_and_new_seqinfo <-
    function(new2old, new_seqinfo, old_seqinfo, context="")
{
    if (length(new_seqinfo) != length(old_seqinfo))
        stop(wmsg("the supplied 'seqinfo' must have the same ",
                  "length as the current 'seqinfo'", context))
    if (!(is.null(new2old) || identical(new2old, seq_along(new_seqinfo))))
        stop(wmsg("'new2old' can only be set to NULL or ",
                  "'seq_along(seqinfo(x))'", context))
    seqnames(old_seqinfo) <- seqnames(new_seqinfo)
    genome(old_seqinfo) <- genome(new_seqinfo)
    if (!identical(new_seqinfo, old_seqinfo))
        stop(wmsg("seqlengths() and isCircular() of the supplied 'seqinfo' ",
                  "must be identical to seqlengths() and isCircular() of ",
                  "the current 'seqinfo'", context))
}

.set_BSgenome_seqinfo <-
    function(x, new2old=NULL,
                pruning.mode=c("error", "coarse", "fine", "tidy"),
                value)
{
    if (!is(value, "Seqinfo"))
        stop(wmsg("the supplied 'seqinfo' must be a Seqinfo object"))
    context <- paste0(" when replacing the 'seqinfo' of a ",
                      classNameForDisplay(x), " object")
    pruning.mode <- match.arg(pruning.mode)
    if (pruning.mode != "error")
        stop(wmsg("'pruning.mode' is not supported", context))
    .check_new2old_and_new_seqinfo(new2old, value, seqinfo(x), context)
    new_seqnames <- seqnames(value)
    if (any(new_seqnames %in% mseqnames(x)))
        stop(wmsg("the supplied 'seqnames' cannot match any of the ",
                  "multiple sequence names (as returned by 'mseqnames(x)')"))
    x@user_seqnames[] <- new_seqnames  # using [] to preserve the names
    genome(x@seqinfo) <- unname(genome(value))
    x
}

setReplaceMethod("seqinfo", "BSgenome", .set_BSgenome_seqinfo)

setReplaceMethod("seqnames", "BSgenome",
    function(x, value)
    {
        seqnames(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SNP related accessors
###

setMethod("SNPlocs_pkgname", "BSgenome",
    function(x) SNPlocs_pkgname(x@injectSNPs_handler)
)

setMethod("snpcount", "BSgenome",
    function(x) snpcount(x@injectSNPs_handler)
)

setMethod("snplocs", "BSgenome",
    function(x, seqname, ...) snplocs(x@injectSNPs_handler, seqname, ...)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BSgenome constructor
###

### 'seqnames' only used for sanity check.
.make_BSgenome_seqinfo <- function(single_sequences, circ_seqs,
                                   genome, seqnames)
{
    seqlengths <- seqlengths(single_sequences)
    if (length(seqnames) == 0L) {
        seqnames <- names(seqlengths)
    } else if (!identical(names(seqlengths), seqnames)) {
        stop("sequence names found in file '",
             seqlengthsFilepath(single_sequences),
             "' are not identical to 'seqnames'. ",
             "May be the data on disk is corrupted?")
    }
    if (identical(circ_seqs, NA)) {
        is_circ <- NA
    } else {
        is_circ <- seqnames %in% circ_seqs
    }
    Seqinfo(seqnames=seqnames,
            seqlengths=seqlengths,
            isCircular=is_circ,
            genome=genome)
}

### NOTES:
### - In BioC 2.14, the 'seqs_pkgname' BSgenome slot was renamed 'pkgname'
###   but the corresponding argument was not renamed for backward
###   compatibility with existing BSgenome packages.
### - In BioC 3.1, the 'species' argument was replaced with the 'common_name'
###   argument but the former was kept for backward compatibility with
###   existing BSgenome packages.
### - In BioC 3.12, the 'provider_version' argument was replaced with the
###   'genome' argument, and the 'release_name' became ignored and was only
###   kept for backward compatibility.
BSgenome <- function(organism, common_name, genome,
                     provider, provider_version,
                     release_date, release_name, source_url,
                     seqnames, circ_seqs=NA, mseqnames,
                     seqs_pkgname, seqs_dirpath,
                     species=NA_character_)
{
    single_sequences <- OnDiskNamedSequences(seqs_dirpath, seqnames=seqnames)
    if (missing(genome))
        genome <- provider_version
    seqinfo <- .make_BSgenome_seqinfo(single_sequences, circ_seqs,
                                      genome, seqnames)
    seqnames <- seqnames(seqinfo)
    if (missing(common_name))
        common_name <- species
    metadata <- list(organism=organism,
                     common_name=common_name,
                     provider=provider,
                     genome=genome,
                     release_date=release_date,
                     source_url=source_url)
    if (is.null(mseqnames))
        mseqnames <- character(0)
    multiple_sequences <- RdaCollection(seqs_dirpath, mseqnames)
    names(user_seqnames) <- user_seqnames <- seqnames

    new("BSgenome", metadata=metadata,
                    pkgname=seqs_pkgname,
                    single_sequences=single_sequences,
                    multiple_sequences=multiple_sequences,
                    seqinfo=seqinfo,
                    user_seqnames=user_seqnames,
                    .seqs_cache=new.env(parent=emptyenv()),
                    .link_counts=new.env(parent=emptyenv())
    )
}

MaskedBSgenome <- function(ref_bsgenome,
                           masks_pkgname, nmask_per_seq, masks_dirpath)
{
    masks <- RdaCollection(masks_dirpath,
                           paste0(seqnames(ref_bsgenome), ".masks"))
    new("MaskedBSgenome", ref_bsgenome,
        .seqs_cache=new.env(parent=emptyenv()),
        .link_counts=new.env(parent=emptyenv()),
        masks_pkgname=masks_pkgname,
        nmask_per_seq=as.integer(nmask_per_seq),
        masks=masks) 
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setMethod("as.list", "BSgenome",
    function(x)
        lapply(setNames(seqnames(x), seqnames(x)),
               function(seqname) x[[seqname]])
)

setAs("BSgenome", "GenomeDescription",
    function(from)
    {
        metadata <- metadata(from)
        GenomeDescription(organism=metadata$organism,
                          common_name=metadata$common_name,
                          provider=metadata$provider,
                          provider_version=metadata$genome,
                          release_date=metadata$release_date,
                          release_name=NA_character_,
                          seqinfo=seqinfo(from))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### bsgenomeName()
###

setMethod("bsgenomeName", "BSgenome",
    function(x) bsgenomeName(as(x, "GenomeDescription"))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method
###

.print_BSgenome_metadata <- function(metadata, margin="")
{
    cat(margin, "BSgenome object", sep="")
    common_name <- metadata$common_name
    if (!is.na(common_name))
        cat(" for ", common_name, sep="")
    cat("\n")
    cat(margin, "- organism: ", metadata$organism, "\n", sep="")
    cat(margin, "- provider: ", metadata$provider, "\n", sep="")
    cat(margin, "- genome: ", metadata$genome, "\n", sep="")
    release_date <- metadata$release_date
    if (!is.na(release_date))
        cat(margin, "- release date: ", release_date, "\n", sep="")
}

.show_BSgenome <- function(x, margin="")
{
    mystrwrap <- function(line)
        writeLines(strwrap(line, width=getOption("width")+1,
                           exdent=0L, prefix=margin))
    .print_BSgenome_metadata(metadata(x), margin=margin)
    if (!is.null(SNPlocs_pkgname(x)))
        cat(margin, "with SNPs injected from package: ",
            SNPlocs_pkgname(x), "\n", sep="")
    cat(margin, "- ", length(seqnames(x)), " sequence(s):\n", sep="")
    margin2 <- paste0(margin, "    ")
    if (length(seqnames(x)) == 0L) {
        cat(margin2, "NONE\n", sep="")
    } else {
        printAtomicVectorInAGrid(seqnames(x), prefix=margin2)
        cat(margin, "\n", sep="")
        mystrwrap(paste0("Tips: call 'seqnames()' on the object to get all ",
                         "the sequence names, call 'seqinfo()' to get the ",
                         "full sequence info, use the '$' or '[[' operator ",
                         "to access a given sequence, see '?BSgenome' for ",
                         "more information."))
    }
}

setMethod("show", "BSgenome",
    function(object) .show_BSgenome(object, margin="| ")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### List-element extraction (with [[).
###

.loadBSgenomeSequence <- function(x, seqname, user_seqname)
{
    idx <- match(user_seqname, x@user_seqnames)
    if (is.na(idx)) {  # multiple sequence
        if (substr(seqname, 1L, 8L) == "upstream") {
            msg <- c(
                "  Starting with BioC 3.0, the upstream sequences ",
                "are defunct.\n",
                "  However they can easily be extracted ",
                "from the full genome\n  sequences with something like ",
                "(for example for hg19):\n\n",
                "      library(TxDb.Hsapiens.UCSC.hg19.knownGene)\n",
                "      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene\n",
                "      gn <- sort(genes(txdb))\n",
                "      up1000 <- flank(gn, width=1000)\n",
                "      library(BSgenome.Hsapiens.UCSC.hg19)\n",
                "      genome <- BSgenome.Hsapiens.UCSC.hg19\n",
                "      up1000seqs <- getSeq(genome, up1000)\n\n",
                "  IMPORTANT: Make sure you use a TxDb package (or ",
                "TranscriptDb object)\n  that contains a gene model ",
                "based on the exact same reference genome\n",
                "  as the BSgenome object you pass to getSeq(). Note that ",
                "you can make\n  your own custom TranscriptDb object from ",
                "various annotation resources.\n",
                "  See the makeTxDbFromUCSC(), ",
                "makeTxDbFromBiomart(), and\n",
                "  makeTxDbFromGFF() functions in the ",
                "GenomicFeatures package."
            )
            .Defunct(msg=paste0(msg, collapse=""))
        }
        ans <- x@multiple_sequences[[seqname]]
        return(ans)
    }
    # single sequence
    ans <- getListElement(x@single_sequences, seqname)
    ## Check the length of the sequence
    if (length(ans) != seqlengths(x)[[user_seqname]]) {
        stop(user_seqname, " sequence does not have the expected length. ",
             "May be the data on disk is corrupted?")
    }
    ## Inject the SNPs, if any
    snps <- snplocs(x, seqname)
    if (!is.null(snps))
        .inplaceReplaceLetterAt(ans, snps$loc, snps$alleles_as_ambig)
    ## Load and set the built-in masks, if any
    nmask_per_seq <- length(masknames(x))
    if (nmask_per_seq > 0L) {
        masks_objname <- paste0(seqname, ".masks")
        builtinmasks <- x@masks[[masks_objname]]
        if (length(builtinmasks) < nmask_per_seq) {
            masks_path <- rdaPath(x@masks, masks_objname)
            stop("expecting ", nmask_per_seq, " built-in masks per ",
                 "single sequence, found only ", length(builtinmasks),
                 " in file '", masks_path, "'. ",
                 "May be the data on disk is corrupted?")
        }
        if (length(builtinmasks) > nmask_per_seq)
            builtinmasks <- builtinmasks[seq_len(nmask_per_seq)]
        if (!identical(names(builtinmasks), masknames(x))) {
            masks_path <- rdaPath(x@masks, masks_objname)
            stop("mask names found in file '", masks_path, "' are not ",
                 "identical to the names returned by masknames(). ",
                 "May be the data on disk is corrupted?")
        }
        masks(ans) <- builtinmasks
    }
    ans
}

.getBSgenomeSequence <- function(x, seqname, user_seqname)
{
    seqs_cache <- x@.seqs_cache
    ## Using the 'if (!exists()) assign(); get()' approach is NOT 100%
    ## reliable:
    ##
    ##   if (!exists(seqname, envir=seqs_cache, inherits=FALSE)) {
    ##       ...
    ##       assign(seqname, ans, envir=seqs_cache)
    ##   }
    ##   get(seqname, envir=seqs_cache, inherits=FALSE)
    ##
    ## because the symbol (seqname) can disappear from the cache between
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
    ans <- try(get(seqname, envir=seqs_cache, inherits=FALSE), silent=TRUE)
    if (is(ans, "try-error")) {
        ans <- .loadBSgenomeSequence(x, seqname, user_seqname)
        if (getOption("verbose"))
            cat("caching ", seqname, "\n", sep="")
        assign(seqname, ans, envir=seqs_cache)
    }
    .linkToCachedObject(ans) <- .newLinkToCachedObject(
                                    seqname,
                                    seqs_cache,
                                    x@.link_counts)
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
        i <- subscripts[[1L]]
        if (length(i) < 1L)
            stop("attempt to select less than one element")
        if (length(i) > 1L)
            stop("attempt to select more than one element")
        if (is.character(i)) {
            user_seqname <- try(match.arg(i, names(x)), silent=TRUE)
            if (is(user_seqname, "try-error"))
                stop("no such sequence")
        } else {
            if (!is.numeric(i) || is.na(i))
                stop("no such sequence")
            i <- as.integer(i)
            if (i < 1L || length(x) < i)
                stop("no such sequence")
            user_seqname <- names(x)[i]
        }
        ## Translate user-supplied sequence named back to original name.
        idx <- match(user_seqname, x@user_seqnames)
        if (is.na(idx)) {
            seqname <- user_seqname  # multiple sequence
        } else {
            seqname <- names(x@user_seqnames)[[idx]]  # single sequence
        }
        .getBSgenomeSequence(x, seqname, user_seqname)
    }
)

setMethod("$", "BSgenome",
    function(x, name) x[[name]]
)

