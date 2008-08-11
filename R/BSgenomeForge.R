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

## mask1 is the mask of "assembly gaps" (inter-contig Ns).
## If 'agp_or_gap' is "agp", then mask1 is extracted from NCBI "agp" files.
## If 'agp_or_gap' is "gap", then mask1 is extracted from UCSC "gap" files.
## If 'agp_or_gap' is NA, then mask1 is an empty mask.
## mask1 is active by default.
.makeMask1 <- function(agp_or_gap, srcdir, seqname, mask_width, prefix, suffix)
{
    if (is.na(agp_or_gap)) {
        ans <- Mask(mask_width)
        names(ans) <- "assembly gaps (empty)"
    } else {
        if (is.na(prefix))
            file <- file.path(srcdir, "gap.txt")
        else
            file <- file.path(srcdir, paste(prefix, seqname, suffix, sep=""))
        if (agp_or_gap == "agp")
            ans <- read.agpMask(file, mask_width, seqname=seqname)
        else
            ans <- read.gapMask(file, mask_width, seqname=seqname)
    }
    active(ans) <- TRUE
    ans
}

## mask2 is the mask of "intra-contig Ns".
## mask2 is active by default.
.makeMask2 <- function(seq, mask1)
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
.makeMask3 <- function(srcdir, seqname, mask_width)
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
.makeMask4 <- function(srcdir, seqname, mask_width)
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

## 'masks_per_seq' must be 1, 2, 3 or 4.
forgeMaskFiles <- function(srcdir, destdir, seqnames, seqdir,
                           masks_per_seq, agp_or_gap, prefix=NA, suffix=NA)
{
    for (seqname in seqnames) {
        ## Get the length of the sequence.
        seqfile <- file.path(seqdir, paste(seqname, ".rda", sep=""))
        load(seqfile)
        seq <- get(seqname)
        mask_width <- length(seq)

        ## Start with an empty mask collection (i.e. a MaskCollection of
        ## length 0).
        masks <- new("MaskCollection", width=mask_width)
        if (masks_per_seq >= 1) {
            mask1 <- .makeMask1(agp_or_gap, srcdir, seqname, mask_width,
                                prefix, suffix)
            masks <- append(masks, mask1)
        }
        if (masks_per_seq >= 2) {
            mask2 <- .makeMask2(seq, mask1)
            masks <- append(masks, mask2)
        }
        remove(seq, list=seqname)
        if (masks_per_seq >= 3) {
            mask3 <- .makeMask3(srcdir, seqname, mask_width)
            masks <- append(masks, mask3)
        }
        if (masks_per_seq >= 4) {
            mask4 <- .makeMask4(srcdir, seqname, mask_width)
            masks <- append(masks, mask4)
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

