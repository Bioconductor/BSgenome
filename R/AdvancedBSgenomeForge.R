### Everything in this file has moved to BSgenomeForge!


### TODO: Move this to S4Vectors (or BiocBaseUtils).
load_package_gracefully <- function(package, ...)
{
    if (!requireNamespace(package, quietly=TRUE))
        stop("Could not load package ", package, ". Is it installed?\n\n  ",
             wmsg("Note that ", ..., " requires the ", package, " package. ",
                  "Please install it with:"),
             "\n\n    BiocManager::install(\"", package, "\")")
}

call_fun_in_BSgenomeForge <- function(fun, ...)
{
    load_package_gracefully("BSgenomeForge", "starting with BioC 3.19, ",
                            "calling ", fun, "()")
    msg <- c(fun, "() has moved from BSgenome to the BSgenomeForge package, ",
             "and is formally deprecated in BSgenome >= 1.75.1. Please call ",
             "BSgenomeForge::", fun, "() to get rid of this warning.")
    .Deprecated(msg=wmsg(msg))
    FUN <- base::get(fun, envir=asNamespace("BSgenomeForge"), inherits=FALSE)
    do.call(FUN, list(...))
}

forgeSeqlengthsRdsFile <- function(...)
{
    call_fun_in_BSgenomeForge("forgeSeqlengthsRdsFile", ...)
}

forgeSeqlengthsRdaFile <- function(...)
{
    call_fun_in_BSgenomeForge("forgeSeqlengthsRdaFile", ...)
}

forgeSeqFiles <- function(...)
{
    call_fun_in_BSgenomeForge("forgeSeqFiles", ...)
}

forgeMasksFiles <- function(...)
{
    call_fun_in_BSgenomeForge("forgeMasksFiles", ...)
}

forgeBSgenomeDataPkg <- function(...)
{
    call_fun_in_BSgenomeForge("forgeBSgenomeDataPkg", ...)
}

forgeMaskedBSgenomeDataPkg <- function(...)
{
    call_fun_in_BSgenomeForge("forgeMaskedBSgenomeDataPkg", ...)
}

