###
###

.pkgname <- "@PKGNAME@"

.nmask_per_seq <- @NMASKPERSEQ@

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export MaskedBSgenome object.
    bsgenome <- MaskedBSgenome(
        ref_bsgenome=@REFPKGNAME@,
        masks_pkgname=pkgname,
        nmask_per_seq=.nmask_per_seq,
        masks_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)
}

