### Everything in this file has moved to BSgenomeForge!


call_fun_in_BSgenomeForge <- function(fun, ...)
{
    S4Vectors:::load_package_gracefully("BSgenomeForge", "in order to use ",
                                        fun, "() in BioC >= 3.19")
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

