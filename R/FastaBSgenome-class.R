### =========================================================================
### FastaBSgenome objects
### -------------------------------------------------------------------------


setClass("FastaBSgenome",
    contains="BSgenome",
    representation(
        ## FaFile object containing the sequences
        fafile="FaFile"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FastaBSgenome constructor (NOT exported)
###

.makeFastaBSgenomeSeqinfo <- function(fafile, circ_seqs, provider_version)
{
    seqlengths <- seqlengths(fafile)
    seqnames <- names(seqlengths)
    if (identical(circ_seqs, NA)) {
        is_circ <- NA
    } else {
        is_circ <- seqnames %in% circ_seqs
    }
    Seqinfo(seqnames=seqnames,
            seqlengths=seqlengths,
            isCircular=is_circ,
            genome=provider_version)
}

FastaBSgenome <- function(filepath, circ_seqs=NA,
                          nmask_per_seq=0, masks_dirpath,
                          organism, species, provider, provider_version,
                          release_date, release_name, source_url,
                          pkgname=NA)
{
    fafile <- FaFile(filepath)
    open(fafile)
    seqinfo <- .makeFastaBSgenomeSeqinfo(fafile, circ_seqs, provider_version)
    user_seqnames <- seqnames(seqinfo)
    names(user_seqnames) <- user_seqnames
    genome_description <- GenomeDescription(organism, species,
                                            provider, provider_version,
                                            release_date, release_name,
                                            seqinfo)
    new("FastaBSgenome", genome_description,
        source_url=source_url,
        user_seqnames=user_seqnames,
        mseqnames=character(0),
        fafile=fafile,
        nmask_per_seq=as.integer(nmask_per_seq),
        masks_pkgname=as.character(pkgname),
        masks_dirpath=masks_dirpath,
        .seqs_cache=new.env(parent=emptyenv()),
        .link_counts=new.env(parent=emptyenv())
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BSgenome constructor
###
### Returns either an RdaBSgenome or a FastaBSgenome object.
###

BSgenome <- function(organism, species, provider, provider_version,
                     release_date, release_name, source_url,
                     seqnames, circ_seqs=NA, mseqnames,
                     seqs_pkgname, seqs_dirpath,
                     nmask_per_seq, masks_pkgname, masks_dirpath)
{
    fasta_filename <- paste0(provider_version, ".fa.rz")
    fasta_filepath <- file.path(seqs_dirpath, fasta_filename)
    if (file.exists(fasta_filepath)) {
        FastaBSgenome(fasta_filepath, circ_seqs,
                      nmask_per_seq, masks_dirpath,
                      organism, species, provider, provider_version,
                      release_date, release_name, source_url,
                      masks_pkgname)
    } else {
        RdaBSgenome(organism, species, provider, provider_version,
                    release_date, release_name, source_url,
                    seqnames, circ_seqs, mseqnames,
                    seqs_pkgname, seqs_dirpath,
                    nmask_per_seq, masks_pkgname, masks_dirpath)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "loadBSgenomeNakedSequence" method
###

.loadFastaSequence <- function(x, seqname)
{
    param <- GRanges(seqname, IRanges(1L, seqlengths(x)[[seqname]]))
    scanFa(x@fafile, param=param)[[1L]]
}

setMethod("loadBSgenomeNakedSequence", "FastaBSgenome",
    function(x, seqname) .loadFastaSequence(x, seqname)
)

