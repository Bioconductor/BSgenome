###
###
###

.seqnames <- @SEQNAMES@

.circ_seqs <- @CIRCSEQS@

.mseqnames <- @MSEQNAMES@

.nmask_per_seq <- @NMASKPERSEQ@

.onLoad <- function(libname, pkgname)
{
    extdata_dir <- system.file("extdata", package=pkgname, lib.loc=libname)
    bsgenome <- BSgenome(
        organism="@ORGANISM@",
        species="@SPECIES@",
        provider="@PROVIDER@",
        provider_version="@PROVIDERVERSION@",
        release_date="@RELEASEDATE@",
        release_name="@RELEASENAME@",
        source_url="@SOURCEURL@",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dir=extdata_dir,
        nmask_per_seq=.nmask_per_seq,
        masks_pkgname=pkgname,
        masks_dir=extdata_dir
    )
    objname <- "@BSGENOMEOBJNAME@"
    ns <- asNamespace(pkgname)
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)
}

