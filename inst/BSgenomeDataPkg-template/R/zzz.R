###
###

.pkgname <- "@PKGNAME@"

.seqnames <- @SEQNAMES@

.circ_seqs <- @CIRCSEQS@

.mseqnames <- @MSEQNAMES@

.nmask_per_seq <- @NMASKPERSEQ@

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
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
        seqs_dirpath=extdata_dirpath,
        nmask_per_seq=.nmask_per_seq,
        masks_pkgname=pkgname,
        masks_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "@BSGENOMEOBJNAME@"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

