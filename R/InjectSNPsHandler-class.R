### =========================================================================
### The "InjectSNPsHandler" class
### -------------------------------------------------------------------------

setClass("InjectSNPsHandler",
    representation(
        SNPlocs_pkgname="character",  # single string
        getSNPcount="function",
        loadLoc="function",
        loadAlleles="function",
        getSNPlocs="function",
        seqname_translation_table="character"  # named character vector
    )
)

.check.seqname_translation_table <- function(x, SNPlocs_seqnames,
                                             bsgenome_seqnames)
{
    if (!is.character(x) || any(is.na(x)))
        return("must be a character vector with no NAs")
    if (!all(x %in% SNPlocs_seqnames))
        return("must have all its elements in 'names(getSNPcount())'")
    if (is.null(names(x))
     || any(is.na(names(x)))
     || any(duplicated(names(x))))
        return("must have unique non-NA names")
    if (!all(names(x) %in% bsgenome_seqnames))
        return("has names incompatible with BSgenome seqnames")
    NULL
}

.get.seqname_translation_table <- function(SNPlocs_pkgname,
                                           bsgenome_pkgname,
                                           SNPlocs_seqnames,
                                           bsgenome_seqnames)
{
    library(SNPlocs_pkgname, character.only=TRUE)
    pkgenvir <- as.environment(paste("package", SNPlocs_pkgname, sep=":"))
    COMPATIBLE_BSGENOMES <- try(get("COMPATIBLE_BSGENOMES",
                                    envir=pkgenvir,
                                    inherits=FALSE), silent=TRUE)
    if (is(COMPATIBLE_BSGENOMES, "try-error"))
        return(character(0))
    if (!is.list(COMPATIBLE_BSGENOMES))
        stop("cannot use package ", SNPlocs_pkgname, " for SNP injection:\n",
             "  '", SNPlocs_pkgname, "::COMPATIBLE_BSGENOMES' is not a list")
    seqname_translation_table <- COMPATIBLE_BSGENOMES[[bsgenome_pkgname]]
    if (is.null(seqname_translation_table)) {
        warning(bsgenome_pkgname, " not in ",
                "'", SNPlocs_pkgname, "::COMPATIBLE_BSGENOMES'")
        return(character(0))
    }
    pb <- .check.seqname_translation_table(seqname_translation_table,
                                           SNPlocs_seqnames,
                                           bsgenome_seqnames)
    if (!is.null(pb))
        stop("cannot inject ", SNPlocs_pkgname, " in ",
             bsgenome_pkgname, ":\n",
             "  bad seqname translation table (it ", pb, ")")
    seqname_translation_table
}

.get_SNPlocs_FUN <- function(funname, SNPlocs_pkgname)
{
    pkgenvir <- as.environment(paste("package", SNPlocs_pkgname, sep=":"))
    FUN <- try(get(funname, envir=pkgenvir, inherits=FALSE), silent=TRUE)
    if (!is.function(FUN))
        stop("cannot use package ", SNPlocs_pkgname, " for SNP injection:\n",
             "  it doesn't seem to define (and export) a function called ",
             "'", funname, "'")
    FUN
}

### Calling this constructor has the side effect of loading the SNPlocs
### package!
InjectSNPsHandler <- function(SNPlocs_pkgname, bsgenome_pkgname,
                              bsgenome_seqnames)
{
    if (!isSingleString(SNPlocs_pkgname))
        stop("'SNPlocs_pkgname' must be a single string")
    library(SNPlocs_pkgname, character.only=TRUE)
    getSNPcount <- .get_SNPlocs_FUN("getSNPcount", SNPlocs_pkgname)
    loadLoc <- .get_SNPlocs_FUN(".loadLoc", SNPlocs_pkgname)
    loadAlleles <- .get_SNPlocs_FUN(".loadAlleles", SNPlocs_pkgname)
    getSNPlocs <- .get_SNPlocs_FUN("getSNPlocs", SNPlocs_pkgname)
    seqname_translation_table <- .get.seqname_translation_table(
                                           SNPlocs_pkgname,
                                           bsgenome_pkgname,
                                           names(getSNPcount()),
                                           bsgenome_seqnames)
    #if (length(seqname_translation_table) == 0L
    # && !all(bsgenome_seqnames %in% names(getSNPcount())))
    #    stop("cannot inject ", SNPlocs_pkgname, " in ",
    #         bsgenome_pkgname, ":\n",
    #         "  seqnames are incompatible and no seqname translation\n",
    #         "  table is provided")
    new("InjectSNPsHandler",
        SNPlocs_pkgname=SNPlocs_pkgname,
        getSNPcount=getSNPcount,
        loadLoc=loadLoc,
        loadAlleles=loadAlleles,
        getSNPlocs=getSNPlocs,
        seqname_translation_table=seqname_translation_table)
}

setGeneric("SNPlocs_pkgname", function(x) standardGeneric("SNPlocs_pkgname"))

setMethod("SNPlocs_pkgname", "InjectSNPsHandler",
    function(x)
    {
        if (length(x@SNPlocs_pkgname) == 0L)
            return(NULL)
        x@SNPlocs_pkgname
    }
)

setGeneric("SNPcount", function(x) standardGeneric("SNPcount"))

setMethod("SNPcount", "InjectSNPsHandler",
    function(x)
    {
        if (length(x@SNPlocs_pkgname) == 0L)
            return(NULL)
        ans <- x@getSNPcount()
        if (length(x@seqname_translation_table) == 0L)
            return(ans)
        ans <- ans[x@seqname_translation_table]
        names(ans) <- names(x@seqname_translation_table)
        ans
    }
)

setGeneric("loadLoc", signature="x",
    function(x, seqname) standardGeneric("loadLoc"))

setMethod("loadLoc", "InjectSNPsHandler",
    function(x, seqname)
    {
        if (length(x@SNPlocs_pkgname) == 0L)
            return(NULL)
        if (!isSingleString(seqname))
            stop("'seqname' must be a single string")
        if (!(seqname %in% names(SNPcount(x))))
            return(NULL)
        if (length(x@seqname_translation_table) != 0L)
            seqname <- x@seqname_translation_table[seqname]
        x@loadLoc(seqname)
    }
)

setGeneric("loadAlleles", signature="x",
    function(x, seqname) standardGeneric("loadAlleles"))

setMethod("loadAlleles", "InjectSNPsHandler",
    function(x, seqname)
    {
        if (length(x@SNPlocs_pkgname) == 0L)
            return(NULL)
        if (!isSingleString(seqname))
            stop("'seqname' must be a single string")
        if (!(seqname %in% names(SNPcount(x))))
            return(NULL)
        if (length(x@seqname_translation_table) != 0L)
            seqname <- x@seqname_translation_table[seqname]
        x@loadAlleles(seqname)
    }
)

setGeneric("SNPlocs", signature="x",
    function(x, seqname) standardGeneric("SNPlocs"))

setMethod("SNPlocs", "InjectSNPsHandler",
    function(x, seqname)
    {
        if (length(x@SNPlocs_pkgname) == 0L)
            return(NULL)
        if (!isSingleString(seqname))
            stop("'seqname' must be a single string")
        if (!(seqname %in% names(SNPcount(x))))
            return(NULL)
        if (length(x@seqname_translation_table) != 0L)
            seqname <- x@seqname_translation_table[seqname]
        x@getSNPlocs(seqname)
    }
)

