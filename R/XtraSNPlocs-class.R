### =========================================================================
### XtraSNPlocs objects
### -------------------------------------------------------------------------
###

setClass("XtraSNPlocs",
    representation(
        ## Name of the XtraSNPlocs data package where the XtraSNPlocs
        ## object is defined.
        pkgname="character",

        ## OnDiskLongTable object containing the SNP table.
        snp_table="OnDiskLongTable",

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
    snp_table <- OnDiskLongTable(snp_table_dirpath)
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
### snptable()
###

setGeneric("snptable", function(x) standardGeneric("snptable"))

setMethod("snptable", "XtraSNPlocs", function(x) x@snp_table)


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
        cat("| nb of SNPs: ", nrow(snptable(object)), "\n", sep="")
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
        .get_snpcount_from_breakpoints(breakpoints(snptable(x)), seqlevels(x),
                                       class(x), x@pkgname)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpsBySeqname()
###

.XTRASNPLOCS_VALID_COLUMNS <- c(
    "RefSNP_id",
    "snpClass",
    "seqnames", "start", "end", "width", "strand",
    "alleles",
    "loctype"
)

.get_real_columns_from_user_supplied_columns <- function(columns)
{
    if (!all(columns %in% .XTRASNPLOCS_VALID_COLUMNS))
        stop(wmsg("valid 'columns' are: ",
             paste(.XTRASNPLOCS_VALID_COLUMNS, collapse=", ")))
    real_columns <- setdiff(columns, c("RefSNP_id", "seqnames", "end"))
    if ("end" %in% columns)
        real_columns <- union(real_columns, c("start", "width"))
    real_columns
}

.get_df_for_one_seqname_from_XtraSNPlocs <- function(x, seqname, columns,
                                                     drop.rs.prefix)
{
    real_columns <- .get_real_columns_from_user_supplied_columns(columns)
    x_snptable <- snptable(x)
    x_breakpoints <- breakpoints(x_snptable)
    if (length(seqname) == 0L) {
        blockidx <- 0L
    } else {
        blockidx <- match(seqname, names(x_breakpoints), nomatch=0L)
    }
    df0 <- data.frame(
        lapply(setNames(real_columns, real_columns),
            function(colname)
                getBlockFromOnDiskLongTable(x_snptable, colname, blockidx)
        ),
        stringsAsFactors=FALSE)
    if ("RefSNP_id" %in% columns) {
        if (nrow(df0) == 0L) {
            ans_RefSNP_id <- integer(0)
        } else {
            x_rowids <- rowids(x_snptable)
            idx <- PartitioningByEnd(x_breakpoints)[[blockidx]]
            ans_RefSNP_id <- x_rowids[idx]
            if (!drop.rs.prefix)
                ans_RefSNP_id <- paste0("rs", ans_RefSNP_id)
        }
        df0[["RefSNP_id"]] <- ans_RefSNP_id
    }
    if ("seqnames" %in% columns)
        df0[["seqnames"]] <- factor(seqname)
    if ("strand" %in% columns) {
        strand0 <- df0[["strand"]]
        ans_strand <- rep.int("+", length(strand0))
        minus_idx <- which(strand0 != as.raw(0L))
        ans_strand[minus_idx] <- "-"
        df0[["strand"]] <- ans_strand
    }
    if ("end" %in% columns)
        df0[["end"]] <- df0[["start"]] + df0[["width"]] - 1L
    if ("loctype" %in% columns)
        df0[["loctype"]] <- as.integer(df0[["loctype"]])
    df0[columns]
}

.get_df_by_seqname_from_XtraSNPlocs <- function(x, seqnames, columns,
                                                drop.rs.prefix)
{
    if (length(seqnames) <= 1L)
        return(.get_df_for_one_seqname_from_XtraSNPlocs(x, seqnames, columns,
                                                         drop.rs.prefix))
    dfs <- lapply(seqnames,
               function(seqname)
                   .get_df_for_one_seqname_from_XtraSNPlocs(x, seqname,
                                                             columns,
                                                             drop.rs.prefix))
    do.call(rbind, dfs)
}

.get_GRanges_by_seqname_from_XtraSNPlocs <- function(x, seqnames, columns,
                                                     drop.rs.prefix)
{
    columns <- setdiff(columns, "width")
    columns <- union(columns, c("seqnames", "start", "end", "strand"))
    df <- .get_df_by_seqname_from_XtraSNPlocs(x, seqnames,
                                               columns, drop.rs.prefix)
    makeGRangesFromDataFrame(df, keep.extra.columns=TRUE, seqinfo=seqinfo(x))
}

setGeneric("snpsBySeqname", signature="x",
    function(x, seqnames, ...) standardGeneric("snpsBySeqname")
)

### Returns a GRanges object unless 'as.data.frame=TRUE'.
setMethod("snpsBySeqname", "XtraSNPlocs",
    function(x, seqnames,
             columns=c("seqnames", "start", "end", "strand"),
             drop.rs.prefix=FALSE,
             as.data.frame=FALSE)
    {
        if (!is.character(seqnames)
         || any(is.na(seqnames))
         || any(duplicated(seqnames)))
            stop("'seqnames' must be a character vector ",
                 "with no NAs and no duplicates")
        if (!all(seqnames %in% seqlevels(x)))
            stop(wmsg("'seqnames' must be a subset of: ",
                      paste(seqlevels(x), collapse=", ")))
        if (!is.character(columns) || any(is.na(columns)))
            stop("'columns' must be a character vector with no NAs")
        if (!isTRUEorFALSE(drop.rs.prefix))
            stop("'drop.rs.prefix' must be TRUE or FALSE")
        if (!isTRUEorFALSE(as.data.frame))
            stop("'as.data.frame' must be TRUE or FALSE")
        if (as.data.frame)
            return(.get_df_by_seqname_from_XtraSNPlocs(x, seqnames, columns,
                                                       drop.rs.prefix))
        .get_GRanges_by_seqname_from_XtraSNPlocs(x, seqnames, columns,
                                                 drop.rs.prefix)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpsById()
###

.get_df_by_id_from_XtraSNPlocs <- function(x, ids, ifnotfound, columns)
{
    real_columns <- .get_real_columns_from_user_supplied_columns(columns)
    x_snptable <- snptable(x)
    df0 <- getDataFromOnDiskLongTable(x_snptable, ids, real_columns,
                                      with.batch_label=TRUE,
                                      as.data.frame=TRUE)
    if ("RefSNP_id" %in% columns) {
        df0[["RefSNP_id"]] <- ids
    }
    if ("seqnames" %in% columns) {
        df0[["seqnames"]] <- factor(df0[["batch_label"]], levels=seqlevels(x))
    }
    if ("strand" %in% columns) {
        strand0 <- df0[["strand"]]
        ans_strand <- rep.int("+", length(strand0))
        minus_idx <- which(strand0 != as.raw(0L))
        ans_strand[minus_idx] <- "-"
        df0[["strand"]] <- ans_strand
    }
    if ("end" %in% columns)
        df0[["end"]] <- df0[["start"]] + df0[["width"]] - 1L
    if ("loctype" %in% columns)
        df0[["loctype"]] <- as.integer(df0[["loctype"]])
    df0[columns]
}

.get_GRanges_by_id_from_XtraSNPlocs <- function(x, ids, ifnotfound, columns)
{
}

setGeneric("snpsById", signature="x",
    function(x, ids, ...) standardGeneric("snpsById")
)

### Returns a GRanges object unless 'as.data.frame=TRUE'.
### If 'ifnotfound="error"' and if the function returns then the returned
### object is guaranteed to be parallel to 'ids'.
setMethod("snpsById", "XtraSNPlocs",
    function(x, ids,
             ifnotfound=c("error", "warning", "drop"),
             columns=c("seqnames", "start", "end", "strand"),
             as.data.frame=FALSE)
    {
        if (!(is.character(ids) || is.integer(ids)) || any(is.na(ids)))
            stop("'ids' must be a character or integer vector with no NAs")
        ifnotfound <- match.arg(ifnotfound)
        if (!is.character(columns) || any(is.na(columns)))
            stop("'columns' must be a character vector with no NAs")
        if (!isTRUEorFALSE(as.data.frame))
            stop("'as.data.frame' must be TRUE or FALSE")
        if (as.data.frame)
            return(.get_df_by_id_from_XtraSNPlocs(x, ids,
                                                  ifnotfound, columns))
        .get_GRanges_by_id_from_XtraSNPlocs(x, ids, ifnotfound, columns)
    }
)

