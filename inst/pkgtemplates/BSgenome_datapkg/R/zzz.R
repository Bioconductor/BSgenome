###
###

.pkgname <- "@PKGNAME@"

.seqnames <- @SEQNAMES@

.circ_seqs <- @CIRCSEQS@

.mseqnames <- @MSEQNAMES@

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
        common_name="@COMMONNAME@",
        genome="@GENOME@",
        provider="@PROVIDER@",
        release_date="@RELEASEDATE@",
        source_url="@SOURCEURL@",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "@BSGENOMEOBJNAME@"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

