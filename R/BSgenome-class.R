### =========================================================================
### BSgenome objects
### -------------------------------------------------------------------------


setClass("BSgenome",
    contains="GenomeDescription",
    representation(
        ## Name of the BSgenome data package where the BSgenome object is
        ## defined.
        pkgname="character",

        ## On-disk "single" and "multiple" sequences.
        single_sequences="OnDiskNamedSequences",
        multiple_sequences="RdaCollection",

        ## Permanent URL to the place where the FASTA files that were used
        ## to produce the "single" and "multiple" sequences above can be
        ## found (and downloaded).
        source_url="character",

        ## Named vector representing the translation table from the original
        ## seqnames (as stored in self@seqinfo@seqnames, the 'seqinfo' slot
        ## being inherited from the GenomeDescription class) to the user
        ## seqnames.
        user_seqnames="character",

        ## For SNPs injection.
        injectSNPs_handler="InjectSNPsHandler",

        ## Used for caching the single and multiple on-disk sequences.
        .seqs_cache="environment",
        .link_counts="environment",

        ## TODO: Remove the 2 slots below once all the BSgenome data packages
        ## with masks have been split in 2 packages (naked BSgenome +
        ## MaskedBSgenome).
        nmask_per_seq="integer",
        masks="RdaCollection"
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
### BSgenome constructor
###

### 'seqnames' only used for sanity check.
.makeBSgenomeSeqinfo <- function(single_sequences,
                                 circ_seqs, provider_version,
                                 seqnames)
{
    seqlengths <- seqlengths(single_sequences)
    if (!identical(names(seqlengths), seqnames))
        stop("sequence names found in file '",
             seqlengthsFilepath(single_sequences),
             "' are not identical to 'seqnames'. ",
             "May be the data on disk is corrupted?")
    if (identical(circ_seqs, NA)) {
        is_circ <- NA
    } else {
        is_circ <- seqnames %in% circ_seqs
    }
    Seqinfo(seqnames=seqnames,
            seqlengths=seqlengths,
            isCircular=is_circ,
            genome=provider_version)
}

### NOTE: In BioC 2.14, the 'seqs_pkgname' and 'masks_pkgname' BSgenome slots
### have been replaced with the 'pkgname' slot but the 2 corresponding
### arguments have been kept for backward compatibility.
BSgenome <- function(organism, species, provider, provider_version,
                     release_date, release_name, source_url,
                     seqnames, circ_seqs=NA, mseqnames,
                     seqs_pkgname, seqs_dirpath,
                     nmask_per_seq=0, masks_pkgname=NA, masks_dirpath=NA)
{
    fa_filename <- "single_sequences.fa"
    fa_filepath <- file.path(seqs_dirpath, fa_filename)
    farz_filename <- paste0(fa_filename, ".rz")
    farz_filepath <- file.path(seqs_dirpath, farz_filename)
    if (file.exists(farz_filepath)) {
        single_sequences <- FastaNamedSequences(farz_filepath)
    } else if (file.exists(fa_filepath)) {
        single_sequences <- FastaNamedSequences(fa_filepath)
    } else {
        single_sequences <- RdaNamedSequences(seqs_dirpath, seqnames)
    }
    seqinfo <- .makeBSgenomeSeqinfo(single_sequences,
                                    circ_seqs, provider_version,
                                    seqnames)
    genome_description <- GenomeDescription(organism, species,
                                            provider, provider_version,
                                            release_date, release_name,
                                            seqinfo)
    if (is.null(mseqnames))
        mseqnames <- character(0)
    multiple_sequences <- RdaCollection(seqs_dirpath, mseqnames)
    user_seqnames <- seqnames(seqinfo)
    names(user_seqnames) <- user_seqnames
    nmask_per_seq <- as.integer(nmask_per_seq)
    masks_dirpath <- as.character(masks_dirpath)
    if (nmask_per_seq == 0L || length(seqnames) == 0L) {
        masks <- RdaCollection(masks_dirpath, character(0))
    } else {
        masks <- RdaCollection(masks_dirpath, paste0(seqnames, ".masks"))
    }

    new("BSgenome", genome_description,
        pkgname=seqs_pkgname,
        single_sequences=single_sequences,
        multiple_sequences=multiple_sequences,
        source_url=source_url,
        user_seqnames=user_seqnames,
        .seqs_cache=new.env(parent=emptyenv()),
        .link_counts=new.env(parent=emptyenv()),
        nmask_per_seq=nmask_per_seq,
        masks=masks
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
### List-element extraction (with [[).
###

.loadBSgenomeSequence <- function(x, seqname, user_seqname)
{
    idx <- match(user_seqname, x@user_seqnames)
    if (is.na(idx)) {  # multiple sequence
        if (substr(seqname, 1L, 8L) == "upstream") {
            msg <- c(
                "  Starting with BioC 2.14, upstream sequences ",
                "are deprecated.\n",
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
                "  See the makeTranscriptDbFromUCSC(), ",
                "makeTranscriptDbFromBiomart(), and\n",
                "  makeTranscriptDbFromGFF() functions in the ",
                "GenomicFeatures package."
            )
            .Deprecated(msg=paste0(msg, collapse=""))
        }
        ans <- x@multiple_sequences[[seqname]]
        return(ans)
    }
    # single sequence
    ans <- x@single_sequences[[seqname]]
    ## Check the length of the sequence
    if (length(ans) != seqlengths(x)[[user_seqname]]) {
        stop(user_seqname, " sequence does not have the expected length. ",
             "May be the data on disk is corrupted?")
    }
    ## Inject the SNPs, if any
    snps <- SNPlocs(x, seqname)
    if (!is.null(snps))
        .inplaceReplaceLetterAt(ans, snps$loc, snps$alleles_as_ambig)
    ## Load and set the built-in masks, if any
    nmask_per_seq <- length(masknames(x))
    if (nmask_per_seq > 0L) {
        masks_objname <- paste0(seqname, ".masks")
        builtinmasks <- x@masks[[masks_objname]]
        if (length(builtinmasks) < nmask_per_seq) {
            masks_filepath <- rdaPath(x@masks, masks_objname)
            stop("expecting ", nmask_per_seq, " built-in masks per ",
                 "single sequence, found only ", length(builtinmasks),
                 " in file '", masks_filepath, "'. ",
                 "May be the data on disk is corrupted?")
        }
        if (length(builtinmasks) > nmask_per_seq)
            builtinmasks <- builtinmasks[seq_len(nmask_per_seq)]
        if (!identical(names(builtinmasks), masknames(x))) {
            masks_filepath <- rdaPath(x@masks, masks_objname)
            stop("mask names found in file '", masks_filepath, "' are not ",
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
        i <- subscripts[[1]]
        if (length(i) < 1)
            stop("attempt to select less than one element")
        if (length(i) > 1)
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

setReplaceMethod("[[", "BSgenome",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a \"BSgenome\" object")
    }
)

setMethod("$", "BSgenome",
    function(x, name) x[[name]]
)

