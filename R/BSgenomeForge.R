### =========================================================================
### The BSgenomeForge functions
### -------------------------------------------------------------------------


.getMasksObjname <- function(seqnames)
{
    if (length(seqnames) == 0)
        return(character(0))
    paste(seqnames, ".masks", sep="")
}

.saveObject <- function(object, objname, destdir=".", verbose=TRUE)
{
    assign(objname, object)
    destfile <- file.path(destdir, paste(objname, ".rda", sep=""))
    if (verbose)
        cat("Saving '", objname, "' object to compressed data file '",
            destfile, "' ... ", sep="")
    ## Using compress="xz" (instead of compress="gzip") would produce a .rda
    ## file that is about 20% smaller on disk but it would also take almost 3x
    ## longer to load it later on with load(). Tested on hg19 chr1 with R-2.14
    ## (2011-09-20 r57033). This is why we stick to compress="gzip".
    save(list=objname, file=destfile, compress="gzip")
    if (verbose)
        cat("DONE\n")
    remove(list=objname)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "getSeqSrcpaths" function.
###

getSeqSrcpaths <- function(seqnames, prefix="", suffix=".fa", seqs_srcdir=".")
{
    if (!is.null(seqnames) && !is.character(seqnames))
        stop("'seqnames' must be a character vector (or NULL)")
    if (length(seqnames) == 0) {
        warning("'seqnames' is empty")
        return(character(0))
    }
    if (!isSingleString(prefix))
        stop("'prefix' must be a single string")
    if (!isSingleString(suffix))
        stop("'suffix' must be a single string")
    if (!isSingleString(seqs_srcdir))
        stop("'seqs_srcdir' must be a single string")
    srcfiles <- paste(prefix, seqnames, suffix, sep="")
    ans <- file.path(seqs_srcdir, srcfiles)
    is_OK <- file.exists(ans)
    if (!all(is_OK)) {
        missing_files <- paste(ans[!is_OK], collapse=", ")
        stop(missing_files, ": file(s) not found")
    }
    names(ans) <- seqnames
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The getSeqlengths() and forgeSeqlengthsFile() functions.
###

getSeqlengths <- function(seqnames, prefix="", suffix=".fa", seqs_srcdir=".")
{
    if (length(seqnames) == 0) {
        warning("'seqnames' is empty")
        return(integer(0))
    }
    srcpaths <- getSeqSrcpaths(seqnames, prefix=prefix, suffix=suffix,
                               seqs_srcdir=seqs_srcdir)
    sapply(seqnames, function(seqname)
           {
               srcpath <- srcpaths[[seqname]]
               ans <- fasta.info(srcpath)
               if (length(ans) == 0)
                   stop("In file '", srcpath, "': no sequence found")
               if (length(ans) > 1)
                   warning("In file '", srcpath, "': ", length(ans),
                           " sequences found, using first sequence only")
               if (names(ans)[1] != seqname)
                   warning("In file '", srcpath, "': sequence description \"",
                           names(ans), "\" doesn't match user-specified ",
                           "sequence name \"", seqname, "\"")
               ans[[1]]
           },
           USE.NAMES=TRUE
    )
}

forgeSeqlengthsFile <- function(seqnames, prefix="", suffix=".fa",
                                seqs_srcdir=".", seqs_destdir=".",
                                verbose=TRUE)
{
    seqlengths <- getSeqlengths(seqnames, prefix=prefix, suffix=suffix,
                                seqs_srcdir=seqs_srcdir)
    if (!isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    .saveObject(seqlengths, "seqlengths", destdir=seqs_destdir,
                verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeSeqFiles" function.
###

.forgeRdaSeqFile <- function(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                             is.single.seq=TRUE, verbose=TRUE)
{
    if (!isSingleString(name))
        stop("'name' must be a single string")
    srcpath <- getSeqSrcpaths(name, prefix=prefix, suffix=suffix,
                              seqs_srcdir=seqs_srcdir)
    if (verbose)
        cat("Loading FASTA file '", srcpath, "' in '", name, "' object ... ", sep="")
    seq <- readDNAStringSet(srcpath, "fasta")
    if (verbose)
        cat("DONE\n")
    if (is.single.seq) {
        if (length(seq) == 0)
            stop("file contains no DNA sequence")
        if (length(seq) > 1)
            warning("file contains ", length(seq), " sequences, ",
                    "using the first sequence only")
        seq <- seq[[1]] # now 'seq' is a DNAString object
    }
    .saveObject(seq, name, destdir=seqs_destdir, verbose=verbose)
}

.forgeFastaRzFile <- function(seqnames, prefix, suffix,
                              seqs_srcdir, seqs_destdir,
                              mode, verbose=TRUE)
{
    if (!is.character(seqnames))
        stop("'seqnames' must be a character vector")
    dest_filename <- "single_sequences.fa"
    dest_filepath <- file.path(seqs_destdir, dest_filename)
    for (seqname in seqnames) {
        srcpath <- getSeqSrcpaths(seqname, prefix=prefix, suffix=suffix,
                                  seqs_srcdir=seqs_srcdir)
        if (verbose)
            cat("Loading '", seqname, "' sequence from FASTA file '", srcpath, "' ... ", sep="")
        seq <- readDNAStringSet(srcpath, "fasta")
        if (verbose)
            cat("DONE\n")
        if (length(seq) == 0L)
            stop("file contains no DNA sequence")
        if (length(seq) > 1L) {
            warning("file contains ", length(seq), " sequences, ",
                    "using the first sequence only")
            seq <- seq[1L]
        }
        if (verbose)
            cat("Appending '", seqname, "' sequence to FASTA file '", dest_filepath, "' ... ", sep="")
        names(seq) <- seqname
        writeXStringSet(seq, dest_filepath, append=TRUE, format="fasta", width=50L)
        if (verbose)
            cat("DONE\n")
    }

    if (mode == "fa") {
        ## "fa" mode
        if (verbose)
            cat("Indexing FASTA file '", dest_filepath, "' ... ", sep="")
        indexFa(dest_filepath)
        if (verbose)
            cat("DONE\n")
    } else {
        ## "fa.rz" mode
        if (verbose)
            cat("Compressing FASTA file '", dest_filepath, "' ... ", sep="")
        farz_filepath <- sprintf("%s.rz", dest_filepath)
        razip(dest_filepath, dest=farz_filepath)
        unlink(dest_filepath)
        if (verbose)
            cat("DONE\n")
        if (verbose)
            cat("Indexing compressed FASTA file '", farz_filepath,
                "' ... ", sep="")
        indexFa(farz_filepath)
        if (verbose)
            cat("DONE\n")
    }
}

forgeSeqFiles <- function(seqnames, mseqnames=NULL, prefix="", suffix=".fa",
                          seqs_srcdir=".", seqs_destdir=".",
                          mode=c("rda", "fa", "fa.rz"), verbose=TRUE)
{
    if (length(seqnames) == 0) {
        warning("'seqnames' is empty")
    } else {
        ## just for the side effect of checking the arguments
        getSeqSrcpaths(seqnames, prefix=prefix, suffix=suffix,
                       seqs_srcdir=seqs_srcdir)
    }
    if (length(mseqnames) != 0) {
        if (!is.character(mseqnames))
            stop("'mseqnames' must be a character vector (or NULL)")
        ## just for the side effect of checking the arguments
        getSeqSrcpaths(mseqnames, prefix=prefix, suffix=suffix,
                       seqs_srcdir=seqs_srcdir)
    }
    if (!isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    mode <- match.arg(mode)
    if (mode == "rda") {  # "rda" mode
        for (name in seqnames) {
            .forgeRdaSeqFile(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                             is.single.seq=TRUE, verbose=verbose)
        }
    } else {  # "fa" and "fa.rz" modes
        .forgeFastaRzFile(seqnames, prefix, suffix, seqs_srcdir, seqs_destdir,
                          mode, verbose=verbose)
    }
    for (name in mseqnames) {
        .forgeRdaSeqFile(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                         is.single.seq=FALSE, verbose=verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeMasksFiles" function.
###

## AGAPS is the mask of "assembly gaps" (inter-contig Ns).
## If 'filetype' is:
##        NA: then AGAPS is an empty mask;
##     "gap": then AGAPS is extracted from a UCSC "gap" file;
##     "agp": then AGAPS is extracted from an NCBI "agp" file.
## AGAPS is active by default.
.forge.AGAPSmask <- function(seqname, mask_width, masks_srcdir,
                             filetype, filename,
                             fileprefix, filesuffix)
{
    if (!isSingleStringOrNA(filetype))
        stop("'filetype' must be a single string or NA")
    if (!isSingleStringOrNA(filename))
        stop("'filename' must be a single string or NA")
    if (!isSingleStringOrNA(fileprefix))
        stop("'fileprefix' must be a single string or NA")
    if (!isSingleStringOrNA(filesuffix))
        stop("'filesuffix' must be a single string or NA")
    if (is.na(filetype)) {
        ans <- Mask(mask_width)
        desc(ans) <- "assembly gaps (empty)"
    } else {
        if (is.na(filename))
            filename <- paste(fileprefix, seqname, filesuffix, sep="")
        filepath <- file.path(masks_srcdir, filename)
        if (filetype == "gap")
            ans <- read.gapMask(filepath, seqname=seqname, mask.width=mask_width)
        else
            ans <- read.agpMask(filepath, seqname=seqname, mask.width=mask_width)
    }
    active(ans) <- TRUE
    names(ans) <- "AGAPS"
    ans
}

## AMB is the mask of "intra-contig ambiguities".
## AMB is active by default.
.forge.AMBmask <- function(seq, AGAPSmask)
{
    active(AGAPSmask) <- TRUE
    masks(seq) <- AGAPSmask
    amb_letters <- names(IUPAC_CODE_MAP)[nchar(IUPAC_CODE_MAP) > 1]
    for (amb_letter in amb_letters)
        seq <- maskMotif(seq, amb_letter)
    ans <- collapse(masks(seq)[-1])
    desc(ans) <- "intra-contig ambiguities"
    if (isEmpty(ans))
        desc(ans) <- paste(desc(ans), "(empty)")
    active(ans) <- TRUE
    names(ans) <- "AMB"
    ans
}

## RM is the "RepeatMasker" mask (from the RepeatMasker .out file).
## RM is NOT active by default.
.forge.RMmask <- function(seqname, mask_width, masks_srcdir,
                          filename, fileprefix, filesuffix)
{
    if (!isSingleStringOrNA(filename))
        stop("'filename' must be a single string or NA")
    if (!isSingleStringOrNA(fileprefix))
        stop("'fileprefix' must be a single string or NA")
    if (!isSingleStringOrNA(filesuffix))
        stop("'filesuffix' must be a single string or NA")
    if (is.na(filename))
        filename <- paste(fileprefix, seqname, filesuffix, sep="")
    filepath <- file.path(masks_srcdir, filename)
    if (file.exists(filepath)) {
        ans <- read.rmMask(filepath, seqname=seqname, mask.width=mask_width)
        desc(ans) <- "RepeatMasker"
    } else {
        ans <- Mask(mask_width)
        desc(ans) <- "RepeatMasker (empty)"
    }
    active(ans) <- FALSE
    names(ans) <- "RM"
    ans
}

## TRF is the "Tandem Repeats Finder" mask (from the Tandem Repeats Finder
## .bed file).
## TRF is NOT active by default.
.forge.TRFmask <- function(seqname, mask_width, masks_srcdir,
                           filename, fileprefix, filesuffix)
{
    if (!isSingleStringOrNA(filename))
        stop("'filename' must be a single string or NA")
    if (!isSingleStringOrNA(fileprefix))
        stop("'fileprefix' must be a single string or NA")
    if (!isSingleStringOrNA(filesuffix))
        stop("'filesuffix' must be a single string or NA")
    if (is.na(filename))
        filename <- paste(fileprefix, seqname, filesuffix, sep="")
    filepath <- file.path(masks_srcdir, filename)
    if (file.exists(filepath)) {
        ans <- read.trfMask(filepath, seqname=seqname, mask.width=mask_width)
        desc(ans) <- "Tandem Repeats Finder [period<=12]"
    } else {
        ans <- Mask(mask_width)
        desc(ans) <- "Tandem Repeats Finder [period<=12] (empty)"
    }
    active(ans) <- FALSE
    names(ans) <- "TRF"
    ans
}

.forgeMasksFile <- function(seqname, nmask_per_seq,
                            seqs_destdir=".", mode=c("rda", "fa", "fa.rz"),
                            masks_srcdir=".", masks_destdir=".",
                            AGAPSfiles_type="gap", AGAPSfiles_name=NA,
                            AGAPSfiles_prefix="", AGAPSfiles_suffix="_gap.txt",
                            RMfiles_name=NA, RMfiles_prefix="", RMfiles_suffix=".fa.out",
                            TRFfiles_name=NA, TRFfiles_prefix="", TRFfiles_suffix=".bed",
                            verbose=TRUE)
{
    if (!isSingleString(seqname))
        stop("'seqname' must be a single string")
    if (!is.numeric(nmask_per_seq)
     || length(nmask_per_seq) != 1
     || !(nmask_per_seq %in% 0:4))
        stop("'nmask_per_seq' must be 0, 1, 2, 3 or 4")
    if (nmask_per_seq == 0)
        warning("forging an empty mask collection ('nmask_per_seq' is set to 0)")
    if (!isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    if (!isSingleString(masks_srcdir))
        stop("'masks_srcdir' must be a single string")
    if (!isSingleString(masks_destdir))
        stop("'masks_destdir' must be a single string")

    ## Load the sequence.
    mode <- match.arg(mode)
    if (mode == "rda") {  # "rda" mode
        seqfile <- file.path(seqs_destdir, paste(seqname, ".rda", sep=""))
        load(seqfile)
        seq <- get(seqname)
        remove(list=seqname)
    } else {  # "fa" and "fa.rz" modes
        fa_filename <- "single_sequences.fa"
        if (mode == "fa.rz")
            fa_filename <- paste0(fa_filename, ".rz")
        fa_filepath <- file.path(seqs_destdir, fa_filename)
        fafile <- FaFile(fa_filepath)
        ## TODO: Implement "[[" method for FaFile objects (in Rsamtools) and
        ## use it here (i.e. do 'seq <- fafile[[seqname]]' instead of the 2
        ## lines below).
        param <- GRanges(seqname, IRanges(1L, seqlengths(fafile)[[seqname]]))
        seq <- scanFa(fafile, param=param)[[1L]]
    }
    mask_width <- length(seq)

    ## Start with an empty mask collection (i.e. a MaskCollection of
    ## length 0).
    masks <- new("MaskCollection", width=mask_width)
    if (nmask_per_seq >= 1) {
        AGAPSmask <- .forge.AGAPSmask(seqname, mask_width, masks_srcdir,
                                      AGAPSfiles_type, AGAPSfiles_name,
                                      AGAPSfiles_prefix, AGAPSfiles_suffix)
        masks <- append(masks, AGAPSmask)
    }
    if (nmask_per_seq >= 2) {
        AMBmask <- .forge.AMBmask(seq, AGAPSmask)
        masks <- append(masks, AMBmask)
    }
    if (nmask_per_seq >= 3) {
        RMmask <- .forge.RMmask(seqname, mask_width, masks_srcdir,
                                RMfiles_name, RMfiles_prefix, RMfiles_suffix)
        masks <- append(masks, RMmask)
    }
    if (nmask_per_seq >= 4) {
        TRFmask <- .forge.TRFmask(seqname, mask_width, masks_srcdir,
                                  TRFfiles_name, TRFfiles_prefix, TRFfiles_suffix)
        masks <- append(masks, TRFmask)
    }
    objname <- .getMasksObjname(seqname)
    .saveObject(masks, objname, destdir=masks_destdir, verbose=verbose)
}

forgeMasksFiles <- function(seqnames, nmask_per_seq,
                            seqs_destdir=".", mode=c("rda", "fa", "fa.rz"),
                            masks_srcdir=".", masks_destdir=".",
                            AGAPSfiles_type="gap", AGAPSfiles_name=NA,
                            AGAPSfiles_prefix="", AGAPSfiles_suffix="_gap.txt",
                            RMfiles_name=NA, RMfiles_prefix="", RMfiles_suffix=".fa.out",
                            TRFfiles_name=NA, TRFfiles_prefix="", TRFfiles_suffix=".bed",
                            verbose=TRUE)
{
    if (length(seqnames) == 0)
        warning("'seqnames' is empty")
    for (seqname in seqnames) {
        .forgeMasksFile(seqname, nmask_per_seq,
                        seqs_destdir=seqs_destdir, mode=mode,
                        masks_srcdir=masks_srcdir, masks_destdir=masks_destdir,
                        AGAPSfiles_type=AGAPSfiles_type, AGAPSfiles_name=AGAPSfiles_name,
                        AGAPSfiles_prefix=AGAPSfiles_prefix, AGAPSfiles_suffix=AGAPSfiles_suffix,
                        RMfiles_name=RMfiles_name, RMfiles_prefix=RMfiles_prefix, RMfiles_suffix=RMfiles_suffix,
                        TRFfiles_name=TRFfiles_name, TRFfiles_prefix=TRFfiles_prefix, TRFfiles_suffix=TRFfiles_suffix,
                        verbose=verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The BSgenomeDataPkgSeed and MaskedBSgenomeDataPkgSeed classes and their
### low-level constructors.
###

setClass("BSgenomeDataPkgSeed",
    representation(
        Package="character",
        Title="character",
        Description="character",
        Version="character",
        Author="character",
        Maintainer="character",
        License="character",
        organism="character",
        species="character",
        provider="character",
        provider_version="character",
        release_date="character",
        release_name="character",
        source_url="character",
        organism_biocview="character",
        BSgenomeObjname="character",
        seqnames="character",         # a single string containing R source code
        circ_seqs="character",        # a single string containing R source code
        mseqnames="character",        # a single string containing R source code
        PkgDetails="character",
        SrcDataFiles1="character",
        SrcDataFiles2="character",
        PkgExamples="character",
        seqfiles_prefix="character",
        seqfiles_suffix="character"
    ),
    prototype(
        Author="The Bioconductor Dev Team",
        Maintainer="Bioconductor Package Maintainer <maintainer@bioconductor.org>",
        License="Artistic-2.0",
        source_url="-- information not available --",
        seqnames="NULL",              # equivalent to "character(0)"
        circ_seqs="NULL",             # equivalent to "character(0)"
        mseqnames="NULL",             # equivalent to "character(0)"
        PkgDetails="",
        SrcDataFiles1="-- information not available --",
        SrcDataFiles2="",
        PkgExamples="",
        seqfiles_prefix="",
        seqfiles_suffix=".fa"
    )
)   

setClass("MaskedBSgenomeDataPkgSeed",
    representation(
        Package="character",
        Title="character",
        Description="character",
        Version="character",
        Author="character",
        Maintainer="character",
        RefPkgname="character",
        License="character",
        source_url="character",
        organism_biocview="character",
        nmask_per_seq="integer",      # a single integer
        PkgDetails="character",
        SrcDataFiles1="character",
        SrcDataFiles2="character",
        PkgExamples="character",
        AGAPSfiles_type="character",
        AGAPSfiles_name="character",
        AGAPSfiles_prefix="character",
        AGAPSfiles_suffix="character",
        RMfiles_name="character",
        RMfiles_prefix="character",
        RMfiles_suffix="character",
        TRFfiles_name="character",
        TRFfiles_prefix="character",
        TRFfiles_suffix="character"
    ),
    prototype(
        Author="The Bioconductor Dev Team",
        Maintainer="Bioconductor Package Maintainer <maintainer@bioconductor.org>",
        License="Artistic-2.0",
        source_url="-- information not available --",
        nmask_per_seq=0L,
        PkgDetails="",
        SrcDataFiles1="-- information not available --",
        SrcDataFiles2="",
        PkgExamples="",
        AGAPSfiles_type="gap",
        AGAPSfiles_name=as.character(NA),
        AGAPSfiles_prefix="",
        AGAPSfiles_suffix="_gap.txt",
        RMfiles_name=as.character(NA),
        RMfiles_prefix="",
        RMfiles_suffix=".fa.out",
        TRFfiles_name=as.character(NA),
        TRFfiles_prefix="",
        TRFfiles_suffix=".bed"
    )
)   

### Generic transformation of a named list into an S4 object with automatic
### coercion of the list elements to the required types.
makeS4FromList <- function(Class, x)
{
    if (!is.list(x) || is.null(names(x)))
        stop("'x' must be a named list")
    explicit_slots <- getSlots(Class)[names(x)]
    if (any(is.na(explicit_slots))) {
        invalid_names <- setdiff(names(x), names(getSlots(Class)))
        stop("some names in 'x' are not valid ", Class, " slots (",
             paste(invalid_names, collapse=", "), ")")
    }
    y <- lapply(seq_len(length(x)),
                function(i) {
                    x_elt <- x[[i]]
                    if (is(x_elt, explicit_slots[i]))
                        return(x_elt)
                    as(x_elt, explicit_slots[i])
                })
    names(y) <- names(x)
    y$Class <- Class
    do.call(new, y)
}

BSgenomeDataPkgSeed <- function(x)
    makeS4FromList("BSgenomeDataPkgSeed", x)

MaskedBSgenomeDataPkgSeed <- function(x)
    makeS4FromList("MaskedBSgenomeDataPkgSeed", x)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeBSgenomeDataPkg" generic and methods.
###

setGeneric("forgeBSgenomeDataPkg", signature="x",
    function(x, seqs_srcdir=".", destdir=".",
                mode=c("rda", "fa", "fa.rz"), verbose=TRUE)
        standardGeneric("forgeBSgenomeDataPkg")
)

setMethod("forgeBSgenomeDataPkg", "BSgenomeDataPkgSeed",
    function(x, seqs_srcdir=".", destdir=".",
                mode=c("rda", "fa", "fa.rz"), verbose=TRUE)
    {
        require(Biobase) ||
            stop("the Biobase package is required")
        mode <- match.arg(mode)
        template_path <- system.file("pkgtemplates", "BSgenome_datapkg", package="BSgenome")
        BSgenome_version <- installed.packages()['BSgenome','Version']
        symvals <- list(
            PKGTITLE=x@Title,
            PKGDESCRIPTION=x@Description,
            PKGVERSION=x@Version,
            AUTHOR=x@Author,
            MAINTAINER=x@Maintainer,
            BSGENOMEVERSION=BSgenome_version,
            LIC=x@License,
            ORGANISM=x@organism,
            SPECIES=x@species,
            PROVIDER=x@provider,
            PROVIDERVERSION=x@provider_version,
            RELEASEDATE=x@release_date,
            RELEASENAME=x@release_name,
            SOURCEURL=x@source_url,
            ORGANISMBIOCVIEW=x@organism_biocview,
            BSGENOMEOBJNAME=x@BSgenomeObjname,
            SEQNAMES=x@seqnames,
            CIRCSEQS=x@circ_seqs,
            MSEQNAMES=x@mseqnames,
            PKGDETAILS=x@PkgDetails,
            SRCDATAFILES1=x@SrcDataFiles1,
            SRCDATAFILES2=x@SrcDataFiles2,
            PKGEXAMPLES=x@PkgExamples
        )
        ## Should never happen
        if (any(duplicated(names(symvals)))) {
            str(symvals)
            stop("'symvals' contains duplicated symbols")
        }
        ## All symvals should by single strings (non-NA)
        is_OK <- sapply(symvals, isSingleString)
        if (!all(is_OK)) {
            bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
            stop("values for symbols ", bad_syms, " are not single strings")
        }
        createPackage(x@Package, destdir, template_path, symvals)
        pkgdir <- file.path(destdir, x@Package)
        ## Just to avoid codetools "no visible binding" NOTEs
        .seqnames <- .mseqnames <- NULL
        ## Sourcing this file will set the values of the '.seqnames' and
        ## '.mseqnames' variables
        source(file.path(pkgdir, "R", "zzz.R"), local=TRUE)
        seqs_destdir <- file.path(pkgdir, "inst", "extdata")
        if (mode == "rda") {
            ## Forge the "seqlengths.rda" file
            forgeSeqlengthsFile(.seqnames,
                                prefix=x@seqfiles_prefix,
                                suffix=x@seqfiles_suffix,
                                seqs_srcdir=seqs_srcdir,
                                seqs_destdir=seqs_destdir,
                                verbose=verbose)
        }
        ## Forge the sequence files (either "rda", "fa", or "fa.rz")
        forgeSeqFiles(.seqnames, mseqnames=.mseqnames,
                      prefix=x@seqfiles_prefix,
                      suffix=x@seqfiles_suffix,
                      seqs_srcdir=seqs_srcdir,
                      seqs_destdir=seqs_destdir,
                      mode=mode,
                      verbose=verbose)
    }
)

setMethod("forgeBSgenomeDataPkg", "list",
    function(x, seqs_srcdir=".", destdir=".",
                mode=c("rda", "fa", "fa.rz"), verbose=TRUE)
    {
        y <- BSgenomeDataPkgSeed(x)
        forgeBSgenomeDataPkg(y,
            seqs_srcdir=seqs_srcdir, destdir=destdir,
            mode=mode, verbose=verbose)
    }
)

.removeCommentLines <- function(infile=stdin(), outfile=stdout())
{
    if (is.character(infile)) {
        infile <- file(infile, "r")
        on.exit(close(infile))
    }
    if (is.character(outfile)) {
        outfile <- file(outfile, "w")
        on.exit({close(infile); close(outfile)})
    }
    while (TRUE) {
        lines <- readLines(infile, n=25000L)
        if (length(lines) == 0L)
            return()
        keep_it <- substr(lines, 1L, 1L) != "#"
        writeLines(lines[keep_it], outfile)
    }
}

### A "comment aware" version of read.dcf().
read.dcf2 <- function(file, ...)
{
    clean_file <- tempfile()
    .removeCommentLines(file, clean_file)
    on.exit(file.remove(clean_file))
    read.dcf(clean_file, ...)
}

### Return a named character vector.
.readSeedFile <- function(file, verbose=TRUE)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if (file.exists(file)) {
        ## Using 'x["isdir"][[1]]' is safer than using 'x$isdir' or
        ## 'x[["isdir"]]' because it will fail if "isdir" is not a defined
        ## column
        isdir <- file.info(file)["isdir"][[1]]
        if (isdir)
            stop("'", file, "' is a directory, not a seed file")
    } else {
        file0 <- file
        seed_dir <- system.file("extdata", "GentlemanLab", package="BSgenome")
        file <- file.path(seed_dir, file)
        if (!file.exists(file)) {
            file <- paste(file, "-seed", sep="")
            if (!file.exists(file))
                stop("seed file '", file0, "' not found")
        }
        if (verbose)
            cat("Seed file '", file0, "' not found, using file '", file, "'\n",
                sep="")
    }
    ans <- read.dcf2(file)  # a character matrix
    if (nrow(ans) != 1)
        stop("seed file '", file, "' must have exactly 1 record")
    ans[1, , drop=TRUE]
}

setMethod("forgeBSgenomeDataPkg", "character",
    function(x, seqs_srcdir=".", destdir=".",
                mode=c("rda", "fa", "fa.rz"), verbose=TRUE)
    {
        y <- .readSeedFile(x, verbose=verbose)
        y <- as.list(y)
        if (missing(seqs_srcdir)) {
            seqs_srcdir <- y[["seqs_srcdir"]]
            if (is.null(seqs_srcdir))
                stop("'seqs_srcdir' argument is missing, and the ",
                     "'seqs_srcdir' field is missing in seed file")
        }
        y <- y[!(names(y) %in% "seqs_srcdir")]
        forgeBSgenomeDataPkg(y,
            seqs_srcdir=seqs_srcdir, destdir=destdir,
            mode=mode, verbose=verbose)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeMaskedBSgenomeDataPkg" generic and methods.
###

setGeneric("forgeMaskedBSgenomeDataPkg", signature="x",
    function(x, masks_srcdir=".", destdir=".", verbose=TRUE)
        standardGeneric("forgeMaskedBSgenomeDataPkg")
)

setMethod("forgeMaskedBSgenomeDataPkg", "MaskedBSgenomeDataPkgSeed",
    function(x, masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        require(Biobase) ||
            stop("the Biobase package is required")
        require(x@RefPkgname, character.only=TRUE) ||
            stop("the ", x@RefPkgname, " package is required")
        ref_envir <- as.environment(paste0("package:", x@RefPkgname))
        ref_bsgenome <- get(x@RefPkgname, envir=ref_envir)
        template_path <- system.file("pkgtemplates", "MaskedBSgenome_datapkg", package="BSgenome")
        BSgenome_version <- installed.packages()['BSgenome','Version']
        symvals <- list(
            PKGTITLE=x@Title,
            PKGDESCRIPTION=x@Description,
            PKGVERSION=x@Version,
            AUTHOR=x@Author,
            MAINTAINER=x@Maintainer,
            BSGENOMEVERSION=BSgenome_version,
            REFPKGNAME=x@RefPkgname,
            LIC=x@License,
            ORGANISM=organism(ref_bsgenome),
            SPECIES=species(ref_bsgenome),
            PROVIDER=provider(ref_bsgenome),
            PROVIDERVERSION=providerVersion(ref_bsgenome),
            RELEASEDATE=releaseDate(ref_bsgenome),
            RELEASENAME=releaseName(ref_bsgenome),
            SOURCEURL=x@source_url,
            ORGANISMBIOCVIEW=x@organism_biocview,
            NMASKPERSEQ=as.character(x@nmask_per_seq),
            PKGDETAILS=x@PkgDetails,
            SRCDATAFILES1=x@SrcDataFiles1,
            SRCDATAFILES2=x@SrcDataFiles2,
            PKGEXAMPLES=x@PkgExamples
        )
        ## Should never happen
        if (any(duplicated(names(symvals)))) {
            str(symvals)
            stop("'symvals' contains duplicated symbols")
        }
        ## All symvals should by single strings (non-NA)
        is_OK <- sapply(symvals, isSingleString)
        if (!all(is_OK)) {
            bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
            stop("values for symbols ", bad_syms, " are not single strings")
        }
        createPackage(x@Package, destdir, template_path, symvals)
        ## Forge the "*.masks.rda" files
        seqs_destdir <- dirname(seqlengthsFilepath(ref_bsgenome@single_sequences))
        pkgdir <- file.path(destdir, x@Package)
        masks_destdir <- file.path(pkgdir, "inst", "extdata")
        forgeMasksFiles(seqnames(ref_bsgenome), x@nmask_per_seq,
                        seqs_destdir=seqs_destdir, mode="fa.rz",
                        masks_srcdir=masks_srcdir, masks_destdir=masks_destdir,
                        AGAPSfiles_type=x@AGAPSfiles_type, AGAPSfiles_name=x@AGAPSfiles_name,
                        AGAPSfiles_prefix=x@AGAPSfiles_prefix, AGAPSfiles_suffix=x@AGAPSfiles_suffix,
                        RMfiles_name=x@RMfiles_name, RMfiles_prefix=x@RMfiles_prefix, RMfiles_suffix=x@RMfiles_suffix,
                        TRFfiles_name=x@TRFfiles_name, TRFfiles_prefix=x@TRFfiles_prefix, TRFfiles_suffix=x@TRFfiles_suffix,
                        verbose=verbose)
    }
)

setMethod("forgeMaskedBSgenomeDataPkg", "list",
    function(x, masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        y <- MaskedBSgenomeDataPkgSeed(x)
        forgeMaskedBSgenomeDataPkg(y,
            masks_srcdir=masks_srcdir, destdir=destdir,
            verbose=verbose)
    }
)

setMethod("forgeMaskedBSgenomeDataPkg", "character",
    function(x, masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        y <- .readSeedFile(x, verbose=verbose)
        y <- as.list(y)
        if (missing(masks_srcdir)) {
            masks_srcdir <- y[["masks_srcdir"]]
            if (is.null(masks_srcdir))
                stop("'masks_srcdir' argument is missing, and the ",
                     "'masks_srcdir' field is missing in seed file")
        }
        y <- y[!(names(y) %in% "masks_srcdir")]
        forgeMaskedBSgenomeDataPkg(y,
            masks_srcdir=masks_srcdir, destdir=destdir,
            verbose=verbose)
    }
)

