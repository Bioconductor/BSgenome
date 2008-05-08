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

newEmptyMask <- function(width)
{
    nir1 <- Biostrings:::newEmptyNormalIRanges()
    new("MaskCollection", nir_list=list(nir1), width=width, active=TRUE)
}

forgeMaskFiles <- function(srcdir, destdir, names, seqdir)
{
    for (seqname in names) {
        ## Get the length of the sequence
        seqfile <- file.path(seqdir, paste(seqname, ".rda", sep=""))
        load(seqfile)
        seq_len <- length(get(seqname))
        remove(list=seqname)

        masks <- new("MaskCollection", width=seq_len)

        ## Make the "rm" mask (from RepeatMasker .out file)
        file1 <- file.path(srcdir, paste(seqname, ".fa.out", sep=""))
        if (file.exists(file1)) {
            mask1 <- read.rmMask(file1, seq_len)
            #names(mask1) <- "RepeatMasker"
        } else {
            mask1 <- newEmptyMask(seq_len)
            names(mask1) <- "rm (empty)"
        }
        masks <- append(masks, mask1)

        ## Make the "trf mask (Tandem Repeats Finder .bed file)
        file2 <- file.path(srcdir, paste(seqname, ".bed", sep=""))
        if (file.exists(file2)) {
            mask2 <- read.trfMask(file2, seq_len)
            #names(mask2) <- "Tandem Repeats Finder (period<=12)"
        } else {
            mask2 <- newEmptyMask(seq_len)
            names(mask2) <- "trf (empty)"
        }
        masks <- append(masks, mask2)

        ## Save the masks
        objname <- paste("masks.", seqname, sep="")
        assign(objname, masks)
        dest <- file.path(destdir, paste(objname, ".rda", sep=""))
        cat("Saving '", objname, "' object to compressed data file '", dest, "'... ", sep="")
        save(list=objname, file=dest, compress=TRUE)
        cat("DONE\n")
        remove(list=objname)
    }
}

