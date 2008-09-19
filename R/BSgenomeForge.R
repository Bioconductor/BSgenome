### =========================================================================
### The BSgenomeForge functions
### -------------------------------------------------------------------------


.isSingleString <- function(x)
{
    is.character(x) && length(x) == 1 && !is.na(x)
}

.isSingleStringOrNA <- function(x)
{
    is.vector(x) && is.atomic(x) && length(x) == 1 && (is.character(x) || is.na(x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeSeqFiles" function.
###

.forgeSeqFile <- function(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                          is.single.seq=TRUE, verbose=TRUE)
{
    if (!.isSingleString(name))
        stop("'name' must be a single string")
    srcfile <- paste(prefix, name, suffix, sep="")
    srcpath <- file.path(seqs_srcdir, srcfile)
    if (verbose)
        cat("Loading FASTA file '", srcpath, "' in '", name, "' object... ", sep="")
    seq <- read.DNAStringSet(srcpath, "fasta")
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
    assign(name, seq)
    dest <- file.path(seqs_destdir, paste(name, ".rda", sep=""))
    if (verbose)
        cat("Saving '", name, "' object to compressed data file '", dest, "'... ", sep="")
    save(list=name, file=dest, compress=TRUE)
    if (verbose)
        cat("DONE\n")
    remove(list=name)
}

forgeSeqFiles <- function(seqnames, mseqnames=NULL, prefix="", suffix=".fa",
                          seqs_srcdir=".", seqs_destdir=".", verbose=TRUE)
{
    if (!is.null(seqnames) && !is.character(seqnames))
        stop("'seqnames' must be a character vector (or NULL)")
    if (length(seqnames) == 0)
        warning("'seqnames' is empty")
    if (!is.null(mseqnames) && !is.character(mseqnames))
        stop("'mseqnames' must be a character vector (or NULL)")
    if (!.isSingleString(prefix))
        stop("'prefix' must be a single string")
    if (!.isSingleString(suffix))
        stop("'suffix' must be a single string")
    if (!.isSingleString(seqs_srcdir))
        stop("'seqs_srcdir' must be a single string")
    if (!.isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")

    srcfiles <- paste(prefix, c(seqnames, mseqnames), suffix, sep="")
    srcpaths <- file.path(seqs_srcdir, srcfiles)
    is_OK <- file.exists(srcpaths)
    if (!all(is_OK)) {
        missing_files <- paste(srcfiles[!is_OK], collapse=", ")
        stop(missing_files, ": files not found in directory ", seqs_srcdir)
    }
    for (name in seqnames) {
        .forgeSeqFile(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                      is.single.seq=TRUE, verbose=verbose)
    }
    for (name in mseqnames) {
        .forgeSeqFile(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                      is.single.seq=FALSE, verbose=verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeMaskFiles" function.
###

## mask1 is the mask of "assembly gaps" (inter-contig Ns).
## If 'mask1.srctype' is NA, then mask1 is an empty mask.
## If 'mask1.srctype' is "gap", then mask1 is extracted from UCSC "gap" files.
## If 'mask1.srctype' is "agp", then mask1 is extracted from NCBI "agp" files.
## mask1 is active by default.
.forgeMask1 <- function(seqname, mask_width, masks_srcdir,
                        mask1.srctype, mask1.prefix, mask1.suffix)
{
    if (!.isSingleStringOrNA(mask1.srctype))
        stop("'mask1.srctype' must be a single string or NA")
    if (!.isSingleStringOrNA(mask1.prefix))
        stop("'mask1.prefix' must be a single string or NA")
    if (!.isSingleStringOrNA(mask1.suffix))
        stop("'mask1.suffix' must be a single string or NA")
    if (is.na(mask1.srctype)) {
        ans <- Mask(mask_width)
        names(ans) <- "assembly gaps (empty)"
    } else {
        if (is.na(mask1.prefix))
            file <- "gap.txt"
        else
            file <- paste(mask1.prefix, seqname, mask1.suffix, sep="")
        file <- file.path(masks_srcdir, file)
        if (mask1.srctype == "gap")
            ans <- read.gapMask(file, mask_width, seqname=seqname)
        else
            ans <- read.agpMask(file, mask_width, seqname=seqname)
    }
    active(ans) <- TRUE
    ans
}

## mask2 is the mask of "intra-contig Ns".
## mask2 is active by default.
.forgeMask2 <- function(seq, mask1)
{
    active(mask1) <- TRUE
    masks(seq) <- mask1
    ans <- masks(maskMotif(seq, "N"))[2]
    names(ans) <- "intra-contig Ns"
    if (isEmpty(ans))
        names(ans) <- paste(names(ans), "(empty)")
    active(ans) <- TRUE
    ans
}

## mask3 is the "RepeatMasker" mask (from the RepeatMasker .out file).
## mask3 is NOT active by default.
.forgeMask3 <- function(seqname, mask_width, masks_srcdir)
{
    file <- file.path(masks_srcdir, paste(seqname, ".fa.out", sep=""))
    if (file.exists(file)) {
        ans <- read.rmMask(file, mask_width)
        names(ans) <- "RepeatMasker"
    } else {
        ans <- Mask(mask_width)
        names(ans) <- "RepeatMasker (empty)"
    }
    active(ans) <- FALSE
    ans
}

## mask4 is the "Tandem Repeats Finder" mask (from the Tandem Repeats Finder
## .bed file).
## mask4 is NOT active by default.
.forgeMask4 <- function(seqname, mask_width, masks_srcdir)
{
    file <- file.path(masks_srcdir, paste(seqname, ".bed", sep=""))
    if (file.exists(file)) {
        ans <- read.trfMask(file, mask_width)
        names(ans) <- "Tandem Repeats Finder [period<=12]"
    } else {
        ans <- Mask(mask_width)
        names(ans) <- "Tandem Repeats Finder [period<=12] (empty)"
    }
    active(ans) <- FALSE
    ans
}

.forgeMaskFile <- function(seqname, nmask_per_seq,
                           seqs_destdir=".", masks_srcdir=".", masks_destdir=".",
                           mask1.srctype="gap", mask1.prefix="", mask1.suffix="_gap.txt",
                           verbose=TRUE)
{
    if (!.isSingleString(seqname))
        stop("'seqname' must be a single string")
    if (!is.numeric(nmask_per_seq)
     || length(nmask_per_seq) != 1
     || !(nmask_per_seq %in% 0:4))
        stop("'nmask_per_seq' must be 0, 1, 2, 3 or 4")
    if (nmask_per_seq == 0)
        warning("forging an empty mask collection ('nmask_per_seq' is set to 0)")
    if (!.isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    if (!.isSingleString(masks_srcdir))
        stop("'masks_srcdir' must be a single string")
    if (!.isSingleString(masks_destdir))
        stop("'masks_destdir' must be a single string")

    ## Get the length of the sequence.
    seqfile <- file.path(seqs_destdir, paste(seqname, ".rda", sep=""))
    load(seqfile)
    seq <- get(seqname)
    mask_width <- length(seq)

    ## Start with an empty mask collection (i.e. a MaskCollection of
    ## length 0).
    masks <- new("MaskCollection", width=mask_width)
    if (nmask_per_seq >= 1) {
        mask1 <- .forgeMask1(seqname, mask_width, masks_srcdir,
                             mask1.srctype, mask1.prefix, mask1.suffix)
        masks <- append(masks, mask1)
    }
    if (nmask_per_seq >= 2) {
        mask2 <- .forgeMask2(seq, mask1)
        masks <- append(masks, mask2)
    }
    remove(seq, list=seqname)
    if (nmask_per_seq >= 3) {
        mask3 <- .forgeMask3(seqname, mask_width, masks_srcdir)
        masks <- append(masks, mask3)
    }
    if (nmask_per_seq >= 4) {
        mask4 <- .forgeMask4(seqname, mask_width, masks_srcdir)
        masks <- append(masks, mask4)
    }

    ## Save the masks.
    objname <- paste("masks.", seqname, sep="")
    assign(objname, masks)
    dest <- file.path(masks_destdir, paste(objname, ".rda", sep=""))
    if (verbose)
        cat("Saving '", objname, "' object to compressed data file '", dest, "'... ", sep="")
    save(list=objname, file=dest, compress=TRUE)
    if (verbose)
        cat("DONE\n")
    remove(list=objname)
}

forgeMaskFiles <- function(seqnames, nmask_per_seq,
                           seqs_destdir=".", masks_srcdir=".", masks_destdir=".",
                           mask1.srctype="gap", mask1.prefix="", mask1.suffix="_gap.txt",
                           verbose=TRUE)
{
    if (length(seqnames) == 0)
        warning("'seqnames' is empty")
    for (seqname in seqnames) {
        .forgeMaskFile(seqname, nmask_per_seq,
                       seqs_destdir=seqs_destdir,
                       masks_srcdir=masks_srcdir, masks_destdir=masks_destdir,
                       mask1.srctype=mask1.srctype,
                       mask1.prefix=mask1.prefix, mask1.suffix=mask1.suffix,
                       verbose=verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "BSgenomeDataPkgSeed" class.
###

setClass(
    "BSgenomeDataPkgSeed",
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
        seqfiles.prefix="character",
        seqfiles.suffix="character",
        seqnames="character",         # a single string containing R source code
        mseqnames="character",        # a single string containing R source code
        nmask_per_seq="integer",      # a single integer
        mask1.srctype="character",
        mask1.prefix="character",
        mask1.suffix="character",
        PkgDetails="character",
        SrcDataFiles="character",
        PkgExamples="character"
    ),
    prototype(
        Author="H. Pages",
        Maintainer="Biocore Team c/o BioC user list <bioconductor@stat.math.ethz.ch>",
        License="Artistic-2.0",
        source_url="-- information not available --",
        seqfiles.prefix="",
        seqfiles.suffix=".fa",
        seqnames="NULL",              # equivalent to "character(0)"
        mseqnames="NULL",
        nmask_per_seq=0L,
        mask1.srctype="gap",
        mask1.prefix="",
        mask1.suffix="_gap.txt",
        PkgDetails="",
        SrcDataFiles="-- information not available --",
        PkgExamples=""
    )
)   


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeBSgenomeDataPkg" generic and methods.
###

setGeneric("forgeBSgenomeDataPkg", signature="x",
    function(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)
        standardGeneric("forgeBSgenomeDataPkg")
)

setMethod("forgeBSgenomeDataPkg", "BSgenomeDataPkgSeed",
    function(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        template_path <- system.file("BSgenomeDataPkg-template", package="BSgenome")
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
            MSEQNAMES=x@mseqnames,
            NMASKPERSEQ=as.character(x@nmask_per_seq),
            PKGDETAILS=x@PkgDetails,
            SRCDATAFILES=x@SrcDataFiles,
            PKGEXAMPLES=x@PkgExamples
        )
        ## Should never happen
        if (any(duplicated(names(symvals)))) {
            str(symvals)
            stop("'symvals' contains duplicated symbols")
        }
        ## All symvals should by single strings (non-NA)
        is_OK <- sapply(symvals,
                        function(val)
                          { is.character(val) && length(val) == 1 && !is.na(val) }
                 )
        if (!all(is_OK)) {
            bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
            stop("values for symbols ", bad_syms, " are not single strings")
        }
        createPackage(x@Package, destdir, template_path, symvals)
        pkgdir <- file.path(destdir, x@Package)
        source(file.path(pkgdir, "R", "zzz.R"), local=TRUE)
        bsgenome <- get(x@BSgenomeObjname, inherits=FALSE)
        seqs_destdir <- file.path(pkgdir, "inst", "extdata")
        forgeSeqFiles(seqnames(bsgenome), mseqnames=mseqnames(bsgenome),
                      prefix=x@seqfiles.prefix, suffix=x@seqfiles.suffix,
                      seqs_srcdir=seqs_srcdir,
                      seqs_destdir=seqs_destdir,
                      verbose=verbose)
        if (x@nmask_per_seq > 0) {
            forgeMaskFiles(seqnames(bsgenome), x@nmask_per_seq,
                      seqs_destdir=seqs_destdir,
                      masks_srcdir=masks_srcdir, masks_destdir=file.path(pkgdir, "data"),
                      mask1.srctype=x@mask1.srctype,
                      mask1.prefix=x@mask1.prefix, mask1.suffix=x@mask1.suffix,
                      verbose=verbose)
        }
    }
)

setMethod("forgeBSgenomeDataPkg", "list",
    function(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        storage.mode(x$nmask_per_seq) <- "integer"
        x$Class <- "BSgenomeDataPkgSeed"
        y <- do.call("new", x)
        forgeBSgenomeDataPkg(y,
            seqs_srcdir=seqs_srcdir,
            masks_srcdir=masks_srcdir,
            destdir=destdir, verbose=verbose)
    }
)

removeCommentsFromFile <- function(infile, outfile)
{
    if (!.isSingleString(infile))
        stop("'infile' must be a single string")
    if (!.isSingleString(outfile))
        stop("'outfile' must be a single string")
    if (file.exists(outfile))
        stop("file '", outfile, "' already exists")
    infile <- file(infile, "r")
    #on.exit(close(infile))
    outfile <- file(outfile, "w")
    #on.exit(close(outfile)) # doesn't seem to work
    while (TRUE) {
        text <- readLines(infile, n=1)
        if (length(text) == 0)
            break
        if (substr(text, 1, 1) != "#")
            writeLines(text, outfile)
    }
    close(outfile)
    close(infile)
}

### Return a named character vector.
readSeedFile <- function(file)
{
    if (!.isSingleString(file))
        stop("'file' must be a single string")
    tmp_file <- file.path(tempdir(), "cleanseed9999.dcf")
    removeCommentsFromFile(file, tmp_file)
    ans <- read.dcf(tmp_file)  # a character matrix
    file.remove(tmp_file)
    if (nrow(ans) != 1)
        stop("seed file ", file, " must have exactly 1 record")
    ans[1, , drop=TRUE]
}

setMethod("forgeBSgenomeDataPkg", "character",
    function(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        if (!file.exists(x)) {
            x0 <- x
            seed_dir <- system.file("extdata", "GentlemanLab", package="BSgenome")
            x <- file.path(seed_dir, x)
            if (!file.exists(x)) {
                x <- paste(x, "-seed", sep="")
                if (!file.exists(x))
                    stop("seed file ", x0, " not found")
            }
            if (verbose)
                cat("Seed file ", x0, " not found, using file ", x, "\n", sep="")
        }
        y <- as.list(readSeedFile(x))
        forgeBSgenomeDataPkg(y,
            seqs_srcdir=seqs_srcdir,
            masks_srcdir=masks_srcdir,
            destdir=destdir, verbose=verbose)
    }
)

