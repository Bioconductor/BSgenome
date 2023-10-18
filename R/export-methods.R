### =========================================================================
### Some "export" methods for BSgenome objects
### -------------------------------------------------------------------------
###
### Note: the export() generic is defined in the BiocIO package.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export as FASTA file
###

### This is exported but not yet documented.
### TODO: Document this in export-methods.Rd (add \usage entry).
writeBSgenomeToFasta <- function(x, filepath,
                                 compress=FALSE, compression_level=NA,
                                 verbose=TRUE)
{
    if (!is(x, "BSgenome"))
        stop("'x' must be a BSgenome object")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    append <- FALSE
    for (seqname in seqnames(x)) {
        dna <- x[[seqname]]
        masks(dna) <- NULL
        dna <- DNAStringSet(dna)
        names(dna) <- seqname
        if (verbose)
            cat("writing", seqname, "sequence to file ... ")
        writeXStringSet(dna, filepath=filepath, append=append,
                        compress=compress, compression_level=compression_level)
        if (verbose)
            cat("OK\n")
        append <- TRUE
    }
}

setMethod("export", c("BSgenome", "FastaFile"),
    function(object, con, format,
             compress=FALSE, compression_level=NA, verbose=TRUE)
    {
        writeBSgenomeToFasta(object, path(con),
                             compress=compress,
                             compression_level=compression_level,
                             verbose=verbose)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export as twoBit file
###

### This is exported but not yet documented.
### TODO: Document this in export-methods.Rd (add \usage entry).
writeBSgenomeToTwobit <- function(x, filepath, ...)
{
    i <- 0L
    object <- bsapply(
        new("BSParams", X=x, FUN=function(chr) {
            i <<- i + 1
            rtracklayer:::.DNAString_to_twoBit(chr, seqnames(x)[i])
        }))
    invisible(rtracklayer:::.TwoBits_export(object, filepath))
}

setMethod("export", c("BSgenome", "TwoBitFile"),
    function(object, con, format, ...)
    {
        if (!missing(format))
            rtracklayer:::checkArgFormat(con, format)
        writeBSgenomeToTwobit(object, rtracklayer:::twoBitPath(path(con)), ...)
    }
)

