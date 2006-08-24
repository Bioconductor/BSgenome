buildDataFiles <- function(srcdir, destdir, files, comments)
{
    for (i in 1:length(files)) {
        varname <- files[i]
        srcfile <- paste(varname, ".fa", sep="")
        srcpath <- paste(srcdir, "/", srcfile, sep="")
        cat("Loading FASTA file '", srcpath, "' in '", varname, "' object... ", sep="")
        f <- file(srcpath)
        dnav <- BStringViews(f, "DNAString")
        close(f)
        cat("DONE\n")
        comment(dnav) <- paste(comments[i], " (generated from FASTA file ", srcfile, " from UCSC)", sep="")
        assign(varname, dnav)
        dest <- paste(destdir, "/", varname, ".rda", sep="")
        cat("Saving '", varname, "' object to compressed data file '", dest, "'... ", sep="")
        save(list=varname, file=dest, compress=TRUE)
        cat("DONE\n")
        remove(list=varname)
    }
}
