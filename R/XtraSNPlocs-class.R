### =========================================================================
### XtraSNPlocs objects
### -------------------------------------------------------------------------
###


setClassUnion("OnDiskLongTable_OR_old",
    c("OnDiskLongTable", "OnDiskLongTable_old")
)

setClass("XtraSNPlocs",
    representation(
        ## Name of the XtraSNPlocs data package where the XtraSNPlocs
        ## object is defined.
        pkgname="character",

        ## OnDiskLongTable_old object containing the SNP data.
        snp_table="OnDiskLongTable_OR_old",

        ## Provider of the SNPs (e.g. "dbSNP").
        provider="character",

        ## E.g. "dbSNP Human BUILD 141".
        provider_version="character",

        ## Official release date of the SNPs (e.g. "May 2014").
        release_date="character",

        ## Official release name of the SNPs (e.g. "dbSNP Human BUILD 141").
        release_name="character",

        ## URL to the place where the original SNP data was downloaded from.
        source_data_url="character",

        ## Date the original SNP data was downloaded.
        download_date="character",

        ## Reference genome of the SNPs.
        reference_genome="GenomeDescription"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers
###

.XTRASNPLOCS_COLUMNS <- c(
    "seqnames", "start", "end", "width", "strand",
    "RefSNP_id", "alleles", "snpClass", "loctype"
)

### Only some of the columns above are actually columns of the
### OnDiskLongTable_old component of an XtraSNPlocs object. We call them
### "physical" columns. The other columns are "computed" columns i.e. columns
### that we compute on demand.
.XTRASNPLOCS_PHYSICAL_COLUMNS <- c(
    "snpClass", "alleles", "start", "width", "strand", "loctype"
)

.XtraSNPlocs_check_user_supplied_columns <- function(columns)
{
    if (!is.character(columns) || any(is.na(columns)))
        stop(wmsg("'columns' must be a character vector with no NAs"))
    if (!all(columns %in% .XTRASNPLOCS_COLUMNS))
        stop(wmsg("valid columns: ",
                  paste(.XTRASNPLOCS_COLUMNS, collapse=", ")))
}

.XtraSNPlocs_get_physical_from_user_supplied_columns <- function(columns)
{
    physical_columns <- intersect(columns, .XTRASNPLOCS_PHYSICAL_COLUMNS)
    if ("end" %in% columns)
        physical_columns <- union(physical_columns, c("start", "width"))
    physical_columns
}

.XtraSNPlocs_get_DF_for_seqnames <- function(x, seqnames, columns,
                                             drop.rs.prefix)
{
    physical_columns <-
        .XtraSNPlocs_get_physical_from_user_supplied_columns(columns)
    with_RefSNP_id <- "RefSNP_id" %in% columns
    DF0 <- getBatchesFromOnDiskLongTable_old(x@snp_table, seqnames,
                                         physical_columns,
                                         with.batch_label=TRUE,
                                         with.rowids=with_RefSNP_id,
                                         as.data.frame=with_RefSNP_id)
    if (with_RefSNP_id) {
        ans_RefSNP_id <- rownames(DF0)
        rownames(DF0) <- NULL
        DF0 <- DataFrame(DF0)
        if (!drop.rs.prefix && length(ans_RefSNP_id) != 0L)
            ans_RefSNP_id <- paste0("rs", ans_RefSNP_id)
        DF0[["RefSNP_id"]] <- ans_RefSNP_id
    }
    if ("seqnames" %in% columns)
        DF0[["seqnames"]] <- factor(DF0[["batch_label"]], levels=seqlevels(x))
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

### 'rowidx' must be a list of 2 vectors parallel to each other, as returned
### by rowids2rowidx().
.XtraSNPlocs_get_DF_for_ids <- function(x, rowidx, columns)
{
    physical_columns <-
        .XtraSNPlocs_get_physical_from_user_supplied_columns(columns)
    DF0 <- getRowsByIndexFromOnDiskLongTable_old(x@snp_table, rowidx[[1L]],
                                                 physical_columns,
                                                 with.batch_label=TRUE)
    if ("RefSNP_id" %in% columns)
        DF0[["RefSNP_id"]] <- rowidx[[2L]]
    if ("seqnames" %in% columns)
        DF0[["seqnames"]] <- factor(DF0[["batch_label"]], levels=seqlevels(x))
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
### Getters
###

setMethod("dimnames", "XtraSNPlocs",
    function(x) list(NULL, .XTRASNPLOCS_COLUMNS)
)

setMethod("dim", "XtraSNPlocs",
    function(x)
    {
        x_nrow <- nrow(x@snp_table)
        x_ncol <- length(.XTRASNPLOCS_COLUMNS)
        c(x_nrow, x_ncol)
    }
)

setMethod("provider", "XtraSNPlocs", function(x) x@provider)

setMethod("providerVersion", "XtraSNPlocs", function(x) x@provider_version)

setMethod("releaseDate", "XtraSNPlocs", function(x) x@release_date)

setMethod("releaseName", "XtraSNPlocs", function(x) x@release_name)

setMethod("referenceGenome", "XtraSNPlocs", function(x) x@reference_genome)

setMethod("organism", "XtraSNPlocs",
    function(object) organism(referenceGenome(object))
)

setMethod("commonName", "XtraSNPlocs",
    function(object) commonName(referenceGenome(object))
)

setMethod("seqinfo", "XtraSNPlocs", function(x) seqinfo(referenceGenome(x)))

setMethod("seqnames", "XtraSNPlocs", function(x) seqnames(referenceGenome(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Not intended to be used directly.
### 'download_url' argument is for backward compatibility with XtraSNPlocs
### packages <= 0.99.12.
newXtraSNPlocs <- function(pkgname, snp_table_dirpath,
                           provider, provider_version,
                           release_date, release_name,
                           source_data_url, download_date,
                           reference_genome, download_url="")
{
    if (missing(source_data_url))
        source_data_url <- download_url
    snp_table <- OnDiskLongTable_old(snp_table_dirpath)
    stopifnot(identical(colnames(snp_table), .XTRASNPLOCS_PHYSICAL_COLUMNS))

    new("XtraSNPlocs",
        pkgname=pkgname,
        snp_table=snp_table,
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        source_data_url=source_data_url,
        download_date=download_date,
        reference_genome=reference_genome)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

setMethod("show", "XtraSNPlocs",
    function(object)
    {
        cat("# ", class(object), " object for ", organism(object),
            " (", releaseName(object), ")\n", sep="")
        cat("# reference genome: ",
            providerVersion(referenceGenome(object)), "\n", sep="")
        cat("# nb of SNPs: ", nrow(object), "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpcount()
###

.get_snpcount_from_blocksizes <- function(blocksizes, seqnames,
                                          class, pkgname)
{
    ans <- integer(length(seqnames))
    names(ans) <- seqnames
    ## Most of the times, 'names(blocksizes)' will be identical to
    ## 'seqnames'. However, in the unlikely situation where some
    ## chromosomes have no SNPs, these chromosomes will be missing
    ## from 'names(blocksizes)'. The sanity check below accomodates
    ## for that.
    m <- match(names(blocksizes), seqnames)
    if (any(is.na(m)) || !isStrictlySorted(m))
        stop(wmsg("BSgenome internal error: ",
                  "the data in this ", class, " object is not in the ",
                  "expected format. Please contact the maintainer of the ",
                  pkgname, " package."))
    ans[m] <- blocksizes
    ans
}

setMethod("snpcount", "XtraSNPlocs",
    function(x)
        .get_snpcount_from_blocksizes(blocksizes(x@snp_table), seqlevels(x),
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
    DF <- .XtraSNPlocs_get_DF_for_seqnames(x, seqnames, columns,
                                           drop.rs.prefix)
    makeGRangesFromDataFrame(DF, keep.extra.columns=TRUE, seqinfo=seqinfo(x))
}

### Returns a GRanges object unless 'as.DataFrame=TRUE'.
setMethod("snpsBySeqname", "XtraSNPlocs",
    function(x, seqnames,
             columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
             drop.rs.prefix=FALSE, as.DataFrame=FALSE)
    {
        if (!is.character(seqnames)
         || any(is.na(seqnames))
         || any(duplicated(seqnames)))
            stop(wmsg("'seqnames' must be a character vector ",
                      "with no NAs and no duplicates"))
        if (!all(seqnames %in% seqlevels(x)))
            stop(wmsg("'seqnames' must be a subset of: ",
                      paste(seqlevels(x), collapse=", ")))
        .XtraSNPlocs_check_user_supplied_columns(columns)
        if (!isTRUEorFALSE(drop.rs.prefix))
            stop(wmsg("'drop.rs.prefix' must be TRUE or FALSE"))
        if (!isTRUEorFALSE(as.DataFrame))
            stop(wmsg("'as.DataFrame' must be TRUE or FALSE"))
        if (as.DataFrame)
            return(.XtraSNPlocs_get_DF_for_seqnames(x, seqnames, columns,
                                                    drop.rs.prefix))
        .get_GRanges_by_seqname_from_XtraSNPlocs(x, seqnames, columns,
                                                 drop.rs.prefix)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpsByOverlaps()
###

.to_DataFrame <- function(x, columns)
{
    spatial_colnames <- c("seqnames", "start", "end", "width", "strand")
    spatial_colnames <- intersect(spatial_colnames, columns)
    spatial_cols <- vector(mode="list", length=length(spatial_colnames))
    names(spatial_cols) <- spatial_colnames
    if ("seqnames" %in% spatial_colnames)
        spatial_cols[["seqnames"]] <- seqnames(x)
    if ("start" %in% spatial_colnames)
        spatial_cols[["start"]] <- start(x)
    if ("end" %in% spatial_colnames)
        spatial_cols[["end"]] <- end(x)
    if ("width" %in% spatial_colnames)
        spatial_cols[["width"]] <- width(x)
    if ("strand" %in% spatial_colnames)
        spatial_cols[["strand"]] <- strand(x)
    cbind(DataFrame(spatial_cols), mcols(x))[ , columns, drop=FALSE]
}

### Returns a GRanges object unless 'as.DataFrame=TRUE'.
### Arguments passed thru ... are further arguments to be passed to
### subsetByOverlaps().
setMethod("snpsByOverlaps", "XtraSNPlocs",
    function(x, ranges,
             columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
             drop.rs.prefix=FALSE, as.DataFrame=FALSE,
             ...)
    {
        ranges <- normarg_ranges(ranges)
        dots <- list(...)
        if (isTRUE(dots$invert))
            stop(wmsg("snpsByOverlaps() does not support 'invert=TRUE'"))

        ## The only purpose of the line below is to check that 'x' and 'ranges'
        ## are based on the same reference genome (merge() will raise an error
        ## if they are not).
        merge(seqinfo(x), seqinfo(ranges))
        seqlevels(ranges, pruning.mode="coarse") <-
            intersect(seqlevels(x), seqlevelsInUse(ranges))
        if (!isTRUEorFALSE(as.DataFrame))
            stop(wmsg("'as.DataFrame' must be TRUE or FALSE"))
        snps_by_seqname <- snpsBySeqname(x, seqlevels(ranges),
                                            columns=columns,
                                            drop.rs.prefix=drop.rs.prefix)
        ans <- subsetByOverlaps(snps_by_seqname, ranges, ...)
        if (as.DataFrame)
            ans <- .to_DataFrame(ans, columns)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpsById()
###

.get_GRanges_by_id_from_XtraSNPlocs <- function(x, rowidx, columns)
{
    columns <- .normarg_columns(columns)
    DF <- .XtraSNPlocs_get_DF_for_ids(x, rowidx, columns)
    makeGRangesFromDataFrame(DF, keep.extra.columns=TRUE, seqinfo=seqinfo(x))
}

### Returns a GRanges object unless 'as.DataFrame=TRUE'.
### If 'ifnotfound="error"' and if the function returns then the returned
### object is guaranteed to be parallel to 'ids'.
setMethod("snpsById", "XtraSNPlocs",
    function(x, ids,
             columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
             ifnotfound=c("error", "warning", "drop"),
             as.DataFrame=FALSE)
    {
        user_rowids <- ids2rowids(ids)
        ifnotfound <- match.arg(ifnotfound)
        x_rowids_env <- get_rowids_env_old(x@snp_table)
        rowidx <- rowids2rowidx(user_rowids, ids, x_rowids_env, ifnotfound)
        .XtraSNPlocs_check_user_supplied_columns(columns)
        if (!isTRUEorFALSE(as.DataFrame))
            stop(wmsg("'as.DataFrame' must be TRUE or FALSE"))
        if (as.DataFrame)
            return(.XtraSNPlocs_get_DF_for_ids(x, rowidx, columns))
        .get_GRanges_by_id_from_XtraSNPlocs(x, rowidx, columns)
    }
)

