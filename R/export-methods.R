### =========================================================================
### Some "export" methods for BSgenome objects
### -------------------------------------------------------------------------
###
### Note: the export() generic is defined in the rtracklayer package.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export as FASTA file.
###

setMethod("export", c("BSgenome", "FastaFile"),
    function(object, con, format, ...)
    {
        append <- FALSE
        for (seqname in seqnames(object)) {
            dna <- object[[seqname]]
            masks(dna) <- NULL
            dna <- DNAStringSet(dna)
            names(dna) <- seqname
            writeXStringSet(dna, path(con), append = append, ...)
            append <- TRUE
        }
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Export as twoBit file.
###

setMethod("export", c("BSgenome", "TwoBitFile"),
    function(object, con, format, ...)
    {
        if (!missing(format))
            rtracklayer:::checkArgFormat(con, format)
        i <- 0L
        object <- bsapply(
            new("BSParams", X=object, FUN=function(chr) {
                i <<- i + 1
                rtracklayer:::.DNAString_to_twoBit(chr, seqnames(object)[i])
            }, ...))
        con <- rtracklayer:::twoBitPath(path(con))
        invisible(rtracklayer:::.TwoBits_export(object, con))
    }
)

