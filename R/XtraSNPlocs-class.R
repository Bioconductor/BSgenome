### =========================================================================
### XtraSNPlocs objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### XtraSnpTable class
###
### NOT exported. Used internally in XtraSNPlocs class.
###

setClass("XtraSnpTable", contains="OnDiskLongTable")

.XTRASNPTABLE_REAL_COLUMNS <- c(
    "snpClass",
    "alleles",
    "start", "width", "strand",
    "loctype"
)

.XTRASNPTABLE_COMPUTED_COLUMNS <- c("RefSNP_id", "seqnames", "end")

.XTRASNPTABLE_ALL_COLUMNS <- c(.XTRASNPTABLE_REAL_COLUMNS,
                               .XTRASNPTABLE_COMPUTED_COLUMNS)

.XtraSnpTable <- function(dirpath=".")
{
    odlt <- OnDiskLongTable(dirpath)
    stopifnot(identical(colnames(odlt), .XTRASNPTABLE_REAL_COLUMNS))
    new("XtraSnpTable", odlt)
}

.XtraSnpTable_check_user_supplied_columns <- function(columns)
{
    if (!is.character(columns) || any(is.na(columns)))
        stop(wmsg("'columns' must be a character vector with no NAs"))
    if (!all(columns %in% .XTRASNPTABLE_ALL_COLUMNS))
        stop(wmsg("valid columns: ",
                  paste(.XTRASNPTABLE_ALL_COLUMNS, collapse=", ")))
}

.XtraSnpTable_get_real_from_user_supplied_columns <- function(columns)
{
    real_columns <- setdiff(columns, .XTRASNPTABLE_COMPUTED_COLUMNS)
    if ("end" %in% columns)
        real_columns <- union(real_columns, c("start", "width"))
    real_columns
}

.XtraSnpTable_get_DF_for_seqname <- function(x, seqname, columns,
                                             drop.rs.prefix)
{
    real_columns <- .XtraSnpTable_get_real_from_user_supplied_columns(columns)
    x_breakpoints <- breakpoints(x)
    if (length(seqname) == 0L) {
        blockidx <- 0L
    } else {
        blockidx <- match(seqname, names(x_breakpoints), nomatch=0L)
    }
    DF0 <- data.frame(
        lapply(setNames(real_columns, real_columns),
            function(colname)
                getBlockFromOnDiskLongTable(x, colname, blockidx)
        ),
        stringsAsFactors=FALSE)
    if ("RefSNP_id" %in% columns) {
        if (nrow(DF0) == 0L) {
            ans_RefSNP_id <- integer(0)
        } else {
            x_rowids <- rowids(x)
            idx <- PartitioningByEnd(x_breakpoints)[[blockidx]]
            ans_RefSNP_id <- x_rowids[idx]
            if (!drop.rs.prefix)
                ans_RefSNP_id <- paste0("rs", ans_RefSNP_id)
        }
        DF0[["RefSNP_id"]] <- ans_RefSNP_id
    }
    if ("seqnames" %in% columns)
        DF0[["seqnames"]] <- factor(seqname)
    if ("strand" %in% columns) {
        strand0 <- DF0[["strand"]]
        ans_strand <- rep.int("+", length(strand0))
        minus_idx <- which(strand0 != as.raw(0L))
        ans_strand[minus_idx] <- "-"
        DF0[["strand"]] <- ans_strand
    }
    if ("end" %in% columns)
        DF0[["end"]] <- DF0[["start"]] + DF0[["width"]] - 1L
    if ("loctype" %in% columns)
        DF0[["loctype"]] <- as.integer(DF0[["loctype"]])
    DF0[columns]
}

.XtraSnpTable_get_DF_for_seqnames <- function(x, seqnames, columns,
                                              drop.rs.prefix)
{
    if (length(seqnames) <= 1L)
        return(.XtraSnpTable_get_DF_for_seqname(x, seqnames, columns,
                                                drop.rs.prefix))
    DF_list <- lapply(seqnames,
        function(seqname)
            .XtraSnpTable_get_DF_for_seqname(x, seqname, columns,
                                             drop.rs.prefix))
    do.call(rbind, DF_list)
}

.XtraSnpTable_get_DF_for_ids <- function(x, rowids, ids,
                                         ifnotfound, columns,
                                         x_seqlevels)
{
    real_columns <- .XtraSnpTable_get_real_from_user_supplied_columns(columns)
    DF0 <- getDataFromOnDiskLongTable(x, rowids, real_columns,
                                      with.batch_label=TRUE)
    if ("RefSNP_id" %in% columns)
        DF0[["RefSNP_id"]] <- ids
    if ("seqnames" %in% columns)
        DF0[["seqnames"]] <- factor(DF0[["batch_label"]], levels=x_seqlevels)
    if ("strand" %in% columns) {
        strand0 <- DF0[["strand"]]
        ans_strand <- rep.int("+", length(strand0))
        minus_idx <- which(strand0 != as.raw(0L))
        ans_strand[minus_idx] <- "-"
        DF0[["strand"]] <- ans_strand
    }
    if ("end" %in% columns)
        DF0[["end"]] <- DF0[["start"]] + DF0[["width"]] - 1L
    if ("loctype" %in% columns)
        DF0[["loctype"]] <- as.integer(DF0[["loctype"]])
    DF0[columns]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### XtraSNPlocs class
###

setClass("XtraSNPlocs",
    representation(
        ## Name of the XtraSNPlocs data package where the XtraSNPlocs
        ## object is defined.
        pkgname="character",

        ## XtraSnpTable object containing the SNP table.
        snp_table="XtraSnpTable",

        ## Provider of the SNPs (e.g. "dbSNP").
        provider="character",

        ## E.g. "dbSNP Human BUILD 141".
        provider_version="character",

        ## Official release date of the SNPs (e.g. "May 2014").
        release_date="character",

        ## Official release name of the SNPs (e.g. "dbSNP Human BUILD 141").
        release_name="character",

        ## URL to the place where the original SNP data was downloaded.
        download_url="character",

        ## Date the original SNP data was downloaded.
        download_date="character",

        ## Reference genome of the SNPs.
        reference_genome="GenomeDescription"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("provider", "XtraSNPlocs", function(x) x@provider)

setMethod("providerVersion", "XtraSNPlocs", function(x) x@provider_version)

setMethod("releaseDate", "XtraSNPlocs", function(x) x@release_date)

setMethod("releaseName", "XtraSNPlocs", function(x) x@release_name)

setMethod("referenceGenome", "XtraSNPlocs", function(x) x@reference_genome)

setMethod("organism", "XtraSNPlocs", function(x) organism(referenceGenome(x)))

setMethod("species", "XtraSNPlocs", function(x) species(referenceGenome(x)))

setMethod("seqinfo", "XtraSNPlocs", function(x) seqinfo(referenceGenome(x)))

setMethod("seqnames", "XtraSNPlocs", function(x) seqnames(referenceGenome(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Not intended to be used directly.
newXtraSNPlocs <- function(pkgname, snp_table_dirpath,
                           provider, provider_version,
                           release_date, release_name,
                           download_url, download_date,
                           reference_genome)
{
    snp_table <- .XtraSnpTable(snp_table_dirpath)
    new("XtraSNPlocs",
        pkgname=pkgname,
        snp_table=snp_table,
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        download_url=download_url,
        download_date=download_date,
        reference_genome=reference_genome)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpTable()
###

setGeneric("snpTable", function(x) standardGeneric("snpTable"))

setMethod("snpTable", "XtraSNPlocs", function(x) x@snp_table)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

setMethod("show", "XtraSNPlocs",
    function(object)
    {
        cat(class(object), " object for ", organism(object),
            " (", releaseName(object), ")\n", sep="")
        cat("| reference genome: ",
            providerVersion(referenceGenome(object)), "\n", sep="")
        cat("| nb of SNPs: ", nrow(snpTable(object)), "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpcount()
###

.get_snpcount_from_breakpoints <- function(breakpoints, seqnames,
                                           class, pkgname)
{
    ans <- integer(length(seqnames))
    names(ans) <- seqnames
    ## Most of the times, 'names(breakpoints)' will be identical to
    ## 'seqnames'. However, in the unlikely situation where some
    ## chromosomes have no SNPs, these chromosomes will be missing
    ## from 'names(breakpoints)'. The sanity check below accomodates
    ## for that.
    m <- match(names(breakpoints), seqnames)
    if (any(is.na(m)) || !isStrictlySorted(m))
        stop(wmsg("BSgenome internal error: ",
                  "the data in this ", class, " object is not in the ",
                  "expected format. Please contact the maintainer of the ",
                  pkgname, " package."))
    ans[m] <- diff(c(0L, breakpoints))
    ans
}

setMethod("snpcount", "XtraSNPlocs",
    function(x)
        .get_snpcount_from_breakpoints(breakpoints(snpTable(x)), seqlevels(x),
                                       class(x), x@pkgname)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpsBySeqname()
###

.normarg_columns <- function(columns)
{
    #IGNORED_COLUMNS <- c("seqnames", "start", "end", "width", "strand")
    #ignored_idx <- which(columns %in% IGNORED_COLUMNS)
    #if (length(ignored_idx) != 0L)
    #    warning(wmsg("ignored columns (because they're implicit): ",
    #                 paste(IGNORED_COLUMNS[ignored_idx], collapse=", ")))
    columns <- setdiff(columns, "width")
    union(columns, c("seqnames", "start", "end", "strand"))
}

.get_GRanges_by_seqname_from_XtraSNPlocs <- function(x, seqnames, columns,
                                                     drop.rs.prefix)
{
    columns <- .normarg_columns(columns)
    DF <- .XtraSnpTable_get_DF_for_seqnames(snpTable(x), seqnames, columns,
                                            drop.rs.prefix)
    makeGRangesFromDataFrame(DF, keep.extra.columns=TRUE, seqinfo=seqinfo(x))
}

setGeneric("snpsBySeqname", signature="x",
    function(x, seqnames, ...) standardGeneric("snpsBySeqname")
)

### Returns a GRanges object unless 'as.DataFrame=TRUE'.
setMethod("snpsBySeqname", "XtraSNPlocs",
    function(x, seqnames,
             columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
             drop.rs.prefix=FALSE,
             as.DataFrame=FALSE)
    {
        if (!is.character(seqnames)
         || any(is.na(seqnames))
         || any(duplicated(seqnames)))
            stop(wmsg("'seqnames' must be a character vector ",
                      "with no NAs and no duplicates"))
        if (!all(seqnames %in% seqlevels(x)))
            stop(wmsg("'seqnames' must be a subset of: ",
                      paste(seqlevels(x), collapse=", ")))
        .XtraSnpTable_check_user_supplied_columns(columns)
        if (!isTRUEorFALSE(drop.rs.prefix))
            stop(wmsg("'drop.rs.prefix' must be TRUE or FALSE"))
        if (!isTRUEorFALSE(as.DataFrame))
            stop(wmsg("'as.DataFrame' must be TRUE or FALSE"))
        if (as.DataFrame)
            return(.XtraSnpTable_get_DF_for_seqnames(snpTable(x),
                                                     seqnames, columns,
                                                     drop.rs.prefix))
        .get_GRanges_by_seqname_from_XtraSNPlocs(x, seqnames, columns,
                                                 drop.rs.prefix)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpsById()
###

.get_GRanges_by_id_from_XtraSNPlocs <- function(x, rowids, ids,
                                                ifnotfound, columns)
{
    columns <- .normarg_columns(columns)
    DF <- .XtraSnpTable_get_DF_for_ids(snpTable(x), rowids, ids,
                                       ifnotfound, columns,
                                       seqlevels(x))
    makeGRangesFromDataFrame(DF, keep.extra.columns=TRUE, seqinfo=seqinfo(x))
}

setGeneric("snpsById", signature="x",
    function(x, ids, ...) standardGeneric("snpsById")
)

.normarg_ids <- function(ids)
{
    if (!(is.character(ids) || is.numeric(ids)))
        stop(wmsg("'ids' must be a character or integer vector with no NAs"))
    if (S4Vectors:::anyMissing(ids))
        stop(wmsg("'ids' cannot contain NAs"))
    if (is.numeric(ids)) {
        if (!is.integer(ids))
            ids <- as.integer(ids)
        return(ids)
    }
    prefixes <- unique(substr(ids, 1L, 2L))
    if ("rs" %in% prefixes) {
        if (!setequal(prefixes, "rs"))
            stop(wmsg("'ids' cannot mix SNP ids that are prefixed ",
                      "with \"rs\" with SNP ids that are not"))
        ## Drop the "rs" prefix.
        ids <- substr(ids, 3L, nchar(ids))
    }
    ids <- suppressWarnings(as.integer(ids))
    if (S4Vectors:::anyMissing(ids))
        stop(wmsg("cannot extract the digital part of some SNP ids in 'ids'"))
    ids
}

### Returns a GRanges object unless 'as.DataFrame=TRUE'.
### If 'ifnotfound="error"' and if the function returns then the returned
### object is guaranteed to be parallel to 'ids'.
setMethod("snpsById", "XtraSNPlocs",
    function(x, ids,
             ifnotfound=c("error", "warning", "drop"),
             columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
             as.DataFrame=FALSE)
    {
        rowids <- .normarg_ids(ids)
        ifnotfound <- match.arg(ifnotfound)
        if (ifnotfound != "error")
            stop(wmsg("only 'ifnotfound=\"error\"' is supported for now"))
        .XtraSnpTable_check_user_supplied_columns(columns)
        if (!isTRUEorFALSE(as.DataFrame))
            stop(wmsg("'as.DataFrame' must be TRUE or FALSE"))
        if (as.DataFrame)
            return(.XtraSnpTable_get_DF_for_ids(snpTable(x), rowids, ids,
                                                ifnotfound, columns,
                                                seqlevels(x)))
        .get_GRanges_by_id_from_XtraSNPlocs(x, rowids, ids,
                                            ifnotfound, columns)
    }
)

