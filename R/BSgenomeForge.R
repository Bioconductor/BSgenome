###
###

forgeSeqFiles <- function(srcdir, destdir, names, prefix="", suffix="", comments=NULL, single.seq=TRUE)
{
    if (length(names) == 0)
        return()
    for (i in 1:length(names)) {
        name <- names[i]
        srcfile <- paste(prefix, name, suffix, sep="")
        srcpath <- paste(srcdir, srcfile, sep="/")
        cat("Loading FASTA file '", srcpath, "' in '", name, "' object... ", sep="")
        seq <- read.DNAStringSet(srcpath, "fasta")
        cat("DONE\n")
        if (single.seq) {
            if (length(seq) != 1)
                stop("file should contain exactly one sequence, found ", length(seq))
            seq <- seq[[1]] # now 'seq' is a DNAString object
        }
        if (is.null(comments))
            c <- ""
        else
            c <- comments[i]
        comment(seq) <- paste(c, " (generated from FASTA file ", srcfile, ")", sep="")
        assign(name, seq)
        dest <- paste(destdir, "/", name, ".rda", sep="")
        cat("Saving '", name, "' object to compressed data file '", dest, "'... ", sep="")
        save(list=name, file=dest, compress=TRUE)
        cat("DONE\n")
        remove(list=name)
    }
}

## mask1 is the mask of "assembly gaps" (extracted from NCBI "agp" files or
## UCSC "gap" files).
## If 'agp_or_gap' is NA, mask N-blocks of width >= 4.
.makeMask1 <- function(agp_or_gap, srcdir, seqname, seq, prefix, suffix)
{
    if (is.na(agp_or_gap))
        return(masks(maskMotif(seq, "N", min.block.width=4)))
    file <- file.path(srcdir, paste(prefix, seqname, suffix, sep=""))
    if (agp_or_gap == "agp")
        read.agpMask(file, length(seq), seqname=seqname)
    else
        read.gapMask(file, length(seq), seqname=seqname)
}

## mask2 is the "rm" mask (from the RepeatMasker .out file).
.makeMask2 <- function(srcdir, seqname, mask_width)
{
    file <- file.path(srcdir, paste(seqname, ".fa.out", sep=""))
    if (file.exists(file)) {
        ans <- read.rmMask(file, mask_width)
        names(ans) <- "RepeatMasker"
    } else {
        ans <- Mask(mask_width, start=integer(0), width=integer(0))
        names(ans) <- "RepeatMasker (empty)"
    }
    ans
}

## mask3 is the "trf" mask (from the Tandem Repeats Finder .bed file).
.makeMask3 <- function(srcdir, seqname, mask_width)
{
    file <- file.path(srcdir, paste(seqname, ".bed", sep=""))
    if (file.exists(file)) {
        ans <- read.trfMask(file, mask_width)
        names(ans) <- "Tandem Repeats Finder [period<=12]"
    } else {
        ans <- Mask(mask_width, start=integer(0), width=integer(0))
        names(ans) <- "Tandem Repeats Finder [period<=12] (empty)"
    }
    ans
}

## 'masks_per_seq' must be 1, 2 or 3.
forgeMaskFiles <- function(srcdir, destdir, seqnames, seqdir,
                           masks_per_seq, agp_or_gap, prefix="", suffix="")
{
    for (seqname in seqnames) {
        ## Get the length of the sequence.
        seqfile <- file.path(seqdir, paste(seqname, ".rda", sep=""))
        load(seqfile)
        mask_width <- length(get(seqname))

        ## Start with an empty mask collection (i.e. a MaskCollection of
        ## length 0).
        masks <- new("MaskCollection", width=mask_width)
        if (masks_per_seq >= 1) {
            mask1 <- .makeMask1(agp_or_gap, srcdir, seqname, get(seqname),
                                prefix, suffix)
            masks <- append(masks, mask1)
        }
        remove(list=seqname)
        if (masks_per_seq >= 2) {
            mask2 <- .makeMask2(srcdir, seqname, mask_width)
            masks <- append(masks, mask2)
        }
        if (masks_per_seq >= 3) {
            mask3 <- .makeMask3(srcdir, seqname, mask_width)
            masks <- append(masks, mask3)
        }

        ## Save the masks.
        objname <- paste("masks.", seqname, sep="")
        assign(objname, masks)
        dest <- file.path(destdir, paste(objname, ".rda", sep=""))
        cat("Saving '", objname, "' object to compressed data file '", dest, "'... ", sep="")
        save(list=objname, file=dest, compress=TRUE)
        cat("DONE\n")
        remove(list=objname)
    }
}

