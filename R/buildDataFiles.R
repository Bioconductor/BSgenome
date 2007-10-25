buildDataFiles <- function(srcdir, destdir, names, prefix="", suffix="", comments=NULL, single.seq=TRUE)
{
    if (length(names) == 0)
        return()
    for (i in 1:length(names)) {
        name <- names[i]
        srcfile <- paste(prefix, name, suffix, sep="")
        srcpath <- paste(srcdir, srcfile, sep="/")
        cat("Loading FASTA file '", srcpath, "' in '", name, "' object... ", sep="")
        seq <- read.BStringViews(srcpath, "fasta", "DNAString")
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
