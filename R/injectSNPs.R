### =========================================================================
### SNP injection
### -------------------------------------------------------------------------

.keep_snps_packages <- function(pkgs)
{
    idx1 <- which(substr(pkgs, 1, 8) == "SNPlocs.")
    idx2 <- which(substr(pkgs, 1, 12) == "XtraSNPlocs.")
    idx <- c(idx1, idx2)
    pkgs[idx]
}

installed.SNPs <- function()
{
    pkgs <- installed.packages()[ , "Package"]
    pkgs <- .keep_snps_packages(pkgs)
    names(pkgs) <- NULL
    return(pkgs)
}

available.SNPs <- function(type=getOption("pkgType"))
{
    contrib_url <- get_data_annotation_contrib_url(type)
    pkgs <- available.packages(contrib_url)[, "Package"]
    pkgs <- .keep_snps_packages(pkgs)
    names(pkgs) <- NULL
    return(pkgs)
}

setGeneric("injectSNPs", signature="x",
    function(x, snps) standardGeneric("injectSNPs")
)

### 'snps' can be a SNPlocs object or the name of a SNPlocs package.
setMethod("injectSNPs", "BSgenome",
    function(x, snps)
    {
        if (!is.null(SNPlocs_pkgname(x)))
            stop("SNPs were already injected in genome 'x'. ",
                 "Injecting from more than 1 package is not supported.")
        ans <- x
        ## We want the original sequence names, not the user sequence names,
        ## so we use 'seqnames(x@seqinfo)' instead of 'seqnames(x)'.
        ans@injectSNPs_handler <- InjectSNPsHandler(snps,
                                                    x@pkgname,
                                                    seqnames(x@seqinfo))
        ans@.seqs_cache <- new.env(parent=emptyenv())
        ans@.link_counts <- new.env(parent=emptyenv())
        ans
    }
)

