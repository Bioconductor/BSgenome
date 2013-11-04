### =========================================================================
### RdaBSgenome objects
### -------------------------------------------------------------------------


setClass("RdaBSgenome",
    contains="BSgenome",
    representation(
        ## where to find the serialized objects containing the sequences
        seqs_pkgname="character",
        seqs_dirpath="character"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RdaBSgenome constructor
###

.makeRdaBSgenomeSeqinfo <- function(seqnames, circ_seqs,
                                    seqs_pkgname, seqs_dirpath,
                                    provider_version)
{
    objname <- "seqlengths"
    seqlengths <- loadSingleObject(objname, seqs_dirpath, seqs_pkgname)
    if (!identical(names(seqlengths), seqnames)) {
        filepath <- getObjFilepath(objname, seqs_dirpath)
        stop("sequence names found in file '", filepath, "' are not ",
             "identical to 'seqnames'. ",
             "May be the ", seqs_pkgname, " package is corrupted?")
    }
    if (identical(circ_seqs, NA)) {
        is_circ <- NA
    } else {
        is_circ <- seqnames %in% circ_seqs
    }
    Seqinfo(seqnames=seqnames, seqlengths=seqlengths, isCircular=is_circ,
            genome=provider_version)
}

RdaBSgenome <- function(organism, species, provider, provider_version,
                        release_date, release_name, source_url,
                        seqnames, circ_seqs=NA, mseqnames,
                        seqs_pkgname, seqs_dirpath,
                        nmask_per_seq, masks_pkgname, masks_dirpath)
{
    seqinfo <- .makeRdaBSgenomeSeqinfo(seqnames, circ_seqs,
                                       seqs_pkgname, seqs_dirpath,
                                       provider_version)
    user_seqnames <- seqnames(seqinfo)
    names(user_seqnames) <- user_seqnames
    if (is.null(mseqnames))
        mseqnames <- character(0)
    genome_description <- GenomeDescription(organism, species,
                                            provider, provider_version,
                                            release_date, release_name,
                                            seqinfo)
    new("RdaBSgenome", genome_description,
        source_url=source_url,
        user_seqnames=user_seqnames,
        mseqnames=mseqnames,
        seqs_pkgname=seqs_pkgname,
        seqs_dirpath=seqs_dirpath,
        nmask_per_seq=as.integer(nmask_per_seq),
        masks_pkgname=masks_pkgname,
        masks_dirpath=masks_dirpath,
        .seqs_cache=new.env(parent=emptyenv()),
        .link_counts=new.env(parent=emptyenv())
    )
}

### Temporary alias for backward compatibility.
BSgenome <- RdaBSgenome


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "loadBSgenomeNakedSequence" method
###

setMethod("loadBSgenomeNakedSequence", "RdaBSgenome",
    function(x, seqname)
        loadSingleObject(seqname, x@seqs_dirpath, x@seqs_pkgname)
)

