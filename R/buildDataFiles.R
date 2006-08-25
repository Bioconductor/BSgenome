buildDataFiles <- function(srcdir, destdir, files, prefix="", suffix="", comments=NULL)
{
    for (i in 1:length(files)) {
        varname <- files[i]
        srcfile <- paste(prefix, varname, suffix, sep="")
        srcpath <- paste(srcdir, srcfile, sep="/")
        cat("Loading FASTA file '", srcpath, "' in '", varname, "' object... ", sep="")
        f <- file(srcpath)
        bsv <- BStringViews(f, "DNAString")
        close(f)
        cat("DONE\n")
        if (is.null(comments))
            c <- ""
        else
            c <- comments[i]
        comment(bsv) <- paste(c, " (generated from FASTA file ", srcfile, ")", sep="")
        assign(varname, bsv)
        dest <- paste(destdir, "/", varname, ".rda", sep="")
        cat("Saving '", varname, "' object to compressed data file '", dest, "'... ", sep="")
        save(list=varname, file=dest, compress=TRUE)
        cat("DONE\n")
        remove(list=varname)
    }
}
