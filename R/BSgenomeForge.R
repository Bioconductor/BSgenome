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

.forgeSeqFile <- function(name, prefix="", suffix=".fa",
                          srcdir=".", destdir=".", is.single.seq=TRUE, verbose=TRUE)
{
    if (!.isSingleString(name))
        stop("'name' must be a single string")
    if (!.isSingleString(prefix))
        stop("'prefix' must be a single string")
    if (!.isSingleString(suffix))
        stop("'suffix' must be a single string")
    if (!.isSingleString(srcdir))
        stop("'srcdir' must be a single string")
    if (!.isSingleString(destdir))
        stop("'destdir' must be a single string")
    srcfile <- paste(prefix, name, suffix, sep="")
    srcpath <- file.path(srcdir, srcfile)
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
    dest <- file.path(destdir, paste(name, ".rda", sep=""))
    if (verbose)
        cat("Saving '", name, "' object to compressed data file '", dest, "'... ", sep="")
    save(list=name, file=dest, compress=TRUE)
    if (verbose)
        cat("DONE\n")
    remove(list=name)
}

forgeSeqFiles <- function(seqnames, mseqnames=NULL, prefix="", suffix=".fa",
                          srcdir=".", destdir=".", verbose=TRUE)
{
    if (length(seqnames) == 0)
        warning("'seqnames' is empty")
    for (name in seqnames) {
        .forgeSeqFile(name, prefix=prefix, suffix=suffix,
                      srcdir=srcdir, destdir=destdir, is.single.seq=TRUE, verbose=verbose)
    }
    for (name in mseqnames) {
        .forgeSeqFile(name, prefix=prefix, suffix=suffix,
                      srcdir=srcdir, destdir=destdir, is.single.seq=FALSE, verbose=verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeMaskFiles" function.
###

## mask1 is the mask of "assembly gaps" (inter-contig Ns).
## If 'srctype1' is NA, then mask1 is an empty mask.
## If 'srctype1' is "gap", then mask1 is extracted from UCSC "gap" files.
## If 'srctype1' is "agp", then mask1 is extracted from NCBI "agp" files.
## mask1 is active by default.
.forgeMask1 <- function(seqname, mask_width, srcdir, srctype1, prefix1, suffix1)
{
    if (!.isSingleStringOrNA(srctype1))
        stop("'srctype1' must be a single string or NA")
    if (!.isSingleStringOrNA(prefix1))
        stop("'prefix1' must be a single string or NA")
    if (!.isSingleStringOrNA(suffix1))
        stop("'suffix1' must be a single string or NA")
    if (is.na(srctype1)) {
        ans <- Mask(mask_width)
        names(ans) <- "assembly gaps (empty)"
    } else {
        if (is.na(prefix1))
            file <- "gap.txt"
        else
            file <- paste(prefix1, seqname, suffix1, sep="")
        file <- file.path(srcdir, file)
        if (srctype1 == "gap")
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
.forgeMask3 <- function(seqname, mask_width, srcdir)
{
    file <- file.path(srcdir, paste(seqname, ".fa.out", sep=""))
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
.forgeMask4 <- function(seqname, mask_width, srcdir)
{
    file <- file.path(srcdir, paste(seqname, ".bed", sep=""))
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

.forgeMaskFile <- function(seqname, masks_per_seq,
                           seqdir=".", srcdir=".", destdir=".",
                           srctype1="gap", prefix1="", suffix1="_gap.txt", verbose=TRUE)
{
    if (!.isSingleString(seqname))
        stop("'seqname' must be a single string")
    if (!is.numeric(masks_per_seq)
     || length(masks_per_seq) != 1
     || !(masks_per_seq %in% 0:4))
        stop("'masks_per_seq' must be 0, 1, 2, 3 or 4")
    if (masks_per_seq == 0)
        warning("forging an empty mask collection ('masks_per_seq' is set to 0)")
    if (!.isSingleString(seqdir))
        stop("'seqdir' must be a single string")
    if (!.isSingleString(srcdir))
        stop("'srcdir' must be a single string")
    if (!.isSingleString(destdir))
        stop("'destdir' must be a single string")

    ## Get the length of the sequence.
    seqfile <- file.path(seqdir, paste(seqname, ".rda", sep=""))
    load(seqfile)
    seq <- get(seqname)
    mask_width <- length(seq)

    ## Start with an empty mask collection (i.e. a MaskCollection of
    ## length 0).
    masks <- new("MaskCollection", width=mask_width)
    if (masks_per_seq >= 1) {
        mask1 <- .forgeMask1(seqname, mask_width, srcdir, srctype1, prefix1, suffix1)
        masks <- append(masks, mask1)
    }
    if (masks_per_seq >= 2) {
        mask2 <- .forgeMask2(seq, mask1)
        masks <- append(masks, mask2)
    }
    remove(seq, list=seqname)
    if (masks_per_seq >= 3) {
        mask3 <- .forgeMask3(seqname, mask_width, srcdir)
        masks <- append(masks, mask3)
    }
    if (masks_per_seq >= 4) {
        mask4 <- .forgeMask4(seqname, mask_width, srcdir)
        masks <- append(masks, mask4)
    }

    ## Save the masks.
    objname <- paste("masks.", seqname, sep="")
    assign(objname, masks)
    dest <- file.path(destdir, paste(objname, ".rda", sep=""))
    if (verbose)
        cat("Saving '", objname, "' object to compressed data file '", dest, "'... ", sep="")
    save(list=objname, file=dest, compress=TRUE)
    if (verbose)
        cat("DONE\n")
    remove(list=objname)
}

forgeMaskFiles <- function(seqnames, masks_per_seq,
                           seqdir=".", srcdir=".", destdir=".",
                           srctype1="gap", prefix1="", suffix1="_gap.txt", verbose=TRUE)
{
    if (length(seqnames) == 0)
        warning("'seqnames' is empty")
    for (seqname in seqnames) {
        .forgeMaskFile(seqname, masks_per_seq,
                       seqdir=seqdir, srcdir=srcdir, destdir=destdir,
                       srctype1=srctype1, prefix1=prefix1, suffix1=suffix1, verbose=verbose)
    }
}

