### =========================================================================
### XtraSNPlocs objects
### -------------------------------------------------------------------------
###


.XTRASNPLOCS_COLUMNS <- c(
    "seqnames", "start", "end", "width", "strand",
    "RefSNP_id", "alleles", "snpClass", "loctype"
)

### Only some of the columns above are actually columns of the OnDiskLongTable
### component of an XtraSNPlocs object. We call them "real" columns. The other
### columns are additional columns that we compute on demand.
.XTRASNPLOCS_REAL_COLUMNS <- c(
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

.XtraSNPlocs_get_real_from_user_supplied_columns <- function(columns)
{
    real_columns <- intersect(columns, .XTRASNPLOCS_REAL_COLUMNS)
    if ("end" %in% columns)
        real_columns <- union(real_columns, c("start", "width"))
    real_columns
}

.XtraSNPlocs_get_DF_for_seqnames <- function(x, seqnames, columns,
                                             drop.rs.prefix)
{
    real_columns <- .XtraSNPlocs_get_real_from_user_supplied_columns(columns)
    with_RefSNP_id <- "RefSNP_id" %in% columns
    DF0 <- getBatchesFromOnDiskLongTable(snpData(x), seqnames, real_columns,
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

### Return a list of 2 vectors parallel to each other. The 1st and 2nd vectors
### are respectively the row indices (integer vector) and the 'user_ids' vector
### (which can be character, numeric, or integer), with the not found entries
### removed from both of them.
### Note that, if 'ifnotfound="error"' then the 2 returned vectors are parallel
### to input vectors 'user_rowids' and 'user_ids'.
.XtraSNPlocs_rowids2rowidx <- function(user_rowids, user_ids, x_rowids,
                                       ifnotfound)
{
    if (is.null(x_rowids))
        stop(wmsg("BSgenome internal error: data contains no SNP ids"))
    rowidx <- match(user_rowids, x_rowids)
    notfound_idx <- which(is.na(rowidx))
    if (length(notfound_idx) != 0L) {
        if (length(notfound_idx) <= 10L) {
            ids_to_show <- user_ids[notfound_idx]
        } else {
            ids_to_show <- c(user_ids[notfound_idx[1:9]], "...")
        }
        ids_to_show <- paste0(ids_to_show, collapse=", ")
        if (ifnotfound == "error")
            stop(wmsg("SNP ids not found: ",
                      ids_to_show,
                      "\n\nUse 'ifnotfound=\"drop\"' to drop them."))
        if (ifnotfound == "warning")
            warning(wmsg("SNP ids not found: ", ids_to_show,
                         "\n\nThey were dropped."))
        rowidx <- rowidx[-notfound_idx]
        user_ids <- user_ids[-notfound_idx]
    }
    list(rowidx, user_ids)
}

### 'rowidx' must be a list of 2 vectors parallel to each other, as returned
### by .XtraSNPlocs_rowids2rowidx() above.
.XtraSNPlocs_get_DF_for_ids <- function(x, rowidx, columns)
{
    real_columns <- .XtraSNPlocs_get_real_from_user_supplied_columns(columns)
    DF0 <- getRowsByIndexFromOnDiskLongTable(snpData(x), rowidx[[1L]],
                                             real_columns,
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
### XtraSNPlocs class
###

setClass("XtraSNPlocs",
    representation(
        ## Name of the XtraSNPlocs data package where the XtraSNPlocs
        ## object is defined.
        pkgname="character",

        ## OnDiskLongTable object containing the SNP data.
        snp_data="OnDiskLongTable",

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
### Getters
###

### NOT exported.
setGeneric("snpData", function(x) standardGeneric("snpData"))
setMethod("snpData", "XtraSNPlocs", function(x) x@snp_data)

setMethod("colnames", "XtraSNPlocs",
    function(x, do.NULL=TRUE, prefix="col") .XTRASNPLOCS_COLUMNS
)

setMethod("nrow", "XtraSNPlocs", function(x) nrow(snpData(x)))

setMethod("ncol", "XtraSNPlocs", function(x) length(colnames(x)))

setMethod("dim", "XtraSNPlocs", function(x) c(nrow(x), ncol(x)))

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
newXtraSNPlocs <- function(pkgname, snp_data_dirpath,
                           provider, provider_version,
                           release_date, release_name,
                           download_url, download_date,
                           reference_genome)
{
    snp_data <- OnDiskLongTable(snp_data_dirpath)
    stopifnot(identical(colnames(snp_data), .XTRASNPLOCS_REAL_COLUMNS))

    new("XtraSNPlocs",
        pkgname=pkgname,
        snp_data=snp_data,
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        download_url=download_url,
        download_date=download_date,
        reference_genome=reference_genome)
}


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
        cat("| nb of SNPs: ", nrow(snpData(object)), "\n", sep="")
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
        .get_snpcount_from_blocksizes(blocksizes(snpData(x)), seqlevels(x),
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

setGeneric("snpsBySeqname", signature="x",
    function(x, seqnames, ...) standardGeneric("snpsBySeqname")
)

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

### Same args and signature as GenomicFeatures::transcriptsByOverlaps()
### EXCEPT for 'minoverlap' default value that we set to zero so we also
### get SNPs that are insertions.
setGeneric("snpsByOverlaps", signature="x",
    function(x, ranges, maxgap=0L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"), ...)
        standardGeneric("snpsByOverlaps")
)

### TODO: Avoid code duplication between .normarg_ranges() and
### GenomicAlignments:::.normarg_param().
.normarg_ranges <- function(ranges)
{
    if (isSingleString(ranges)) {
        tmp1 <- strsplit(ranges, ":", fixed=TRUE)[[1L]]
        if (length(tmp1) != 2L) 
            stop(wmsg("when a character string, 'ranges' must be ",
                      "of the form \"ch14:5201-5300\""))
        tmp2 <- as.integer(strsplit(tmp1[2L], "-", fixed=TRUE)[[1L]])
        if (length(tmp2) != 2L || any(is.na(tmp2)))
            stop(wmsg("when a character string, 'ranges' must be ", 
                      "of the form \"ch14:5201-5300\""))
        ranges <- GRanges(tmp1[1L], IRanges(tmp2[1L], tmp2[2L]))
        return(ranges)
    }
    if (!is(ranges, "GenomicRanges"))
        stop(wmsg("'ranges' ranges must be a GenomicRanges object ",
                  "or a character string of the form \"ch14:5201-5300\""))
    ranges
}

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
    function(x, ranges, maxgap=0L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
             drop.rs.prefix=FALSE, as.DataFrame=FALSE, ...)
    {
        ranges <- .normarg_ranges(ranges)
        ## The only purpose of the line below is to check that 'x' and 'ranges'
        ## are based on the same reference genome (merge() will raise an error
        ## if they are not).
        merge(seqinfo(x), seqinfo(ranges))
        seqlevels(ranges, force=TRUE) <- intersect(seqlevels(x),
                                                   seqlevelsInUse(ranges))
        if (!isTRUEorFALSE(as.DataFrame))
            stop(wmsg("'as.DataFrame' must be TRUE or FALSE"))
        snps_by_seqname <- snpsBySeqname(x, seqlevels(ranges),
                                            columns=columns,
                                            drop.rs.prefix=drop.rs.prefix)
        ans <- subsetByOverlaps(snps_by_seqname, ranges,
                                maxgap=maxgap, minoverlap=minoverlap,
                                type=type, ...)
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

setGeneric("snpsById", signature="x",
    function(x, ids, ...) standardGeneric("snpsById")
)

### Return an integer vector with no NAs parallel to 'ids'.
.ids2rowids <- function(ids)
{
    if (!(is.character(ids) || is.numeric(ids)))
        stop(wmsg("'ids' must be a character or integer vector with no NAs"))
    if (S4Vectors:::anyMissing(ids))
        stop(wmsg("'ids' cannot contain NAs"))
    if (is.character(ids)) {
        prefixes <- unique(substr(ids, 1L, 2L))
        if ("rs" %in% prefixes) {
            if (!setequal(prefixes, "rs"))
                stop(wmsg("'ids' cannot mix SNP ids that are prefixed ",
                          "with \"rs\" with SNP ids that are not"))
            ## Drop the "rs" prefix.
            ids <- substr(ids, 3L, nchar(ids))
        }
        ids <- suppressWarnings(as.numeric(ids))
        if (S4Vectors:::anyMissing(ids))
            stop(wmsg("cannot extract the digital part of ",
                      "some SNP ids in 'ids'"))
    }
    if (length(ids) != 0L && min(ids) < 0)
        stop(wmsg("'ids' contains unrealistic SNP ids"))
    if (!is.integer(ids)) {
        ids <- suppressWarnings(as.integer(ids))
        if (S4Vectors:::anyMissing(ids))
            stop(wmsg("'ids' contains SNP ids that are too big"))
    }
    ids
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
        user_rowids <- .ids2rowids(ids)
        ifnotfound <- match.arg(ifnotfound)
        x_rowids <- rowids(snpData(x))
        rowidx <- .XtraSNPlocs_rowids2rowidx(user_rowids, ids, x_rowids,
                                             ifnotfound)
        .XtraSNPlocs_check_user_supplied_columns(columns)
        if (!isTRUEorFALSE(as.DataFrame))
            stop(wmsg("'as.DataFrame' must be TRUE or FALSE"))
        if (as.DataFrame)
            return(.XtraSNPlocs_get_DF_for_ids(x, rowidx, columns))
        .get_GRanges_by_id_from_XtraSNPlocs(x, rowidx, columns)
    }
)

