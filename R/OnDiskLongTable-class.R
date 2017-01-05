### =========================================================================
### OnDiskLongTable objects
### -------------------------------------------------------------------------
###

setClassUnion("GRangesORNULL", c("GRanges", "NULL"))

setClass("OnDiskLongTable",
    representation(

        ## A single string.
        dirpath="character",

        ## Whether on-disk storage is column- (the default) or batch-oriented.
        #bybatch="logical",

        ## A named ordinary list of length the nb of cols where all list
        ## elements are vector-like objects of length 0. The names on the list
        ## are interpreted as the column names and thus cannot be NAs or empty
        ## strings.
        header="list",

        ## An integer vector of length the nb of batches. Contains the
        ## cumulated nb of rows in the batches. A batch cannot be empty so
        ## 'breakpoints' must contain *strictly* sorted positive integers.
        ## If the table has zero rows, then 'breakpoints' is empty. Otherwise,
        ## its last element is the nb of rows in the table.
        breakpoints="integer",

        ## [OPTIONAL] A *sorted* *unstranded* GRanges object with 1 range per
        ## batch. The object must be naked i.e. no names and no metadata
        ## columns.
        spatial_index="GRangesORNULL",

        ## Where to load and cache the row ids. The row ids are stored in a
        ## vector of *unique* integer values of length the nb of rows. This
        ## vector is typically big (OnDiskLongTable objects can have tens or
        ## hundreds of millions of rows), so can take a long time to load from
        ## disk to memory. Caching it allows fast translation from
        ## user-supplied row ids to row indices.
        .rowids_cache="environment"
    ),
    prototype(
        dirpath="."
        #bybatch=FALSE
    )
)

.get_OnDiskLongTable_nrow_from_breakpoints <- function(breakpoints)
{
    breakpoints_len <- length(breakpoints)
    if (breakpoints_len == 0L) 0L else breakpoints[[breakpoints_len]]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid_OnDiskLongTable_dirpath <- function(dirpath)
{
    if (!isSingleString(dirpath) || dirpath == "")
        return(wmsg("'dirpath' must be a single (non-empty) string"))
    if (!dir.exists(dirpath))
        return(wmsg("directory not found: ", dirpath))
    NULL
}

.valid_OnDiskLongTable_header <- function(header)
{
    if (!is.list(header))
        return(wmsg("'header' must be a list"))
    colnames <- names(header)
    if (is.null(colnames))
        return(wmsg("'header' must be a named list"))
    if (any(colnames %in% c(NA_character_, "")))
        return(wmsg("'names(header)' cannot contain NAs or empty strings"))
    if (!all(lengths(header) == 0L))
        return(wmsg("all list elements in 'header' must be ",
                    "vector-like objects of length 0"))
    NULL
}

.valid_OnDiskLongTable_breakpoints <- function(breakpoints)
{
    if (!is.integer(breakpoints)
     || S4Vectors:::anyMissing(breakpoints)
     || is.unsorted(breakpoints, strictly=TRUE)
     || length(breakpoints) != 0L && breakpoints[[1L]] <= 0L)
        return(wmsg("invalid breakpoints found for OnDiskLongTable object"))
    NULL
}

.valid_OnDiskLongTable_spatial_index <- function(spatial_index, breakpoints)
{
    if (is.null(spatial_index))
        return(NULL)
    if (!is(spatial_index, "GRanges"))
        return(wmsg("'spatial_index' must be a GRanges object or NULL"))
    if (length(spatial_index) != length(breakpoints))
        return(wmsg("'spatial_index' must have the length of 'breakpoints'"))
    if (is.unsorted(spatial_index))
        return(wmsg("'spatial_index' must be a sorted GRanges object"))
    if (!all(strand(spatial_index) == "*"))
        return(wmsg("'spatial_index' must be unstranded"))
    if (!is.null(names(spatial_index)))
        return(wmsg("'spatial_index' cannot have names"))
    if (ncol(mcols(spatial_index)) != 0L)
        return(wmsg("'spatial_index' cannot have metadata columns"))
    NULL
}

.valid_OnDiskLongTable <- function(x)
{
    c(.valid_OnDiskLongTable_dirpath(x@dirpath),
      .valid_OnDiskLongTable_header(x@header),
      .valid_OnDiskLongTable_breakpoints(x@breakpoints),
      .valid_OnDiskLongTable_spatial_index(x@spatial_index, x@breakpoints))
}

setValidity2("OnDiskLongTable", .valid_OnDiskLongTable)

.check_OnDiskLongTable_dirpath <- function(dirpath)
{
    errmsg <- .valid_OnDiskLongTable_dirpath(dirpath)
    if (!is.null(errmsg))
        stop(errmsg)
}

.check_OnDiskLongTable_rowids <- function(rowids)
{
    if (!is.integer(rowids))
        stop(wmsg("'rowids' must be an integer vector"))
    if (S4Vectors:::anyMissing(rowids) || anyDuplicated(rowids))
        stop(wmsg("'rowids' cannot contain NAs or duplicated values"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read/write data to/from disk
###

### 'opath' is the relative path (with respect to 'dirpath') to the file that
### is to be read/written. It must NOT have the ".rds" extension.
.make_filepath <- function(dirpath, opath)
    file.path(dirpath, paste0(opath, ".rds"))

.read_object <- function(dirpath, opath)
{
    filepath <- .make_filepath(dirpath, opath)
    suppressWarnings(readRDS(filepath))
}

.write_object <- function(object, dirpath, opath,
                          compress=TRUE, overwrite=FALSE)
{
    filepath <- .make_filepath(dirpath, opath)
    if (!overwrite && file.exists(filepath))
        stop(wmsg("file already exists: ", filepath))
    saveRDS(object, file=filepath, compress=compress)
}

### Make batch/col physical name from batch/col number.
.BATCH_FMT <- "b%05d"  # nb of batches must be <= 99999
.COL_FMT <- "c%03d"    # nb of cols must be <= 999
.batch_physname <- function(b) sprintf(.BATCH_FMT, b)
.col_physname <- function(c) sprintf(.COL_FMT, c)

.block_physname <- function(b, c, bybatch=FALSE)
{
    batch_physname <- .batch_physname(b)
    col_physname <- .col_physname(c)
    #if (bybatch) {
    #    block_physname <- file.path(batch_physname, col_physname)
    #} else {
        block_physname <- file.path(col_physname, batch_physname)
    #}
    block_physname
}

.read_OnDiskLongTable_block <- function(dirpath, b, c, bybatch=FALSE)
{
    block_physname <- .block_physname(b, c, bybatch=bybatch)
    .read_object(dirpath, block_physname)
}

.write_OnDiskLongTable_block <- function(block_data,
                                         dirpath, b, c, bybatch=FALSE,
                                         compress=TRUE)
{
    block_physname <- .block_physname(b, c, bybatch=bybatch)
    .write_object(block_data, dirpath, block_physname, compress=compress)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .write_zero_row_OnDiskLongTable()
###

.remove_file <- function(dirpath, filename)
{
    filepath <- file.path(dirpath, filename)
    if (file.exists(filepath)) {
        if (!file.remove(filepath))
            stop(wmsg("failed to remove file: ", filepath))
    }
}

.remove_dirs <- function(dirpaths)
{
    fail_idx <- which(as.logical(unlink(dirpaths, recursive=TRUE)))
    if (length(fail_idx) != 0L)
        stop(wmsg("failed to remove directories: ",
                  paste0(dirpaths[fail_idx], collapse=", ")))
}

.write_zero_row_OnDiskLongTable <- function(df, dirpath, spatial_index)
{
    ## Remove stuff.
    .remove_file(dirpath, "rowids.rds")
    col_dirpaths <- file.path(dirpath, .col_physname(seq_len(ncol(df))))
    .remove_dirs(col_dirpaths)

    ## Create stuff.
    header <- df[0, , drop=FALSE]
    .write_object(header, dirpath, "header", overwrite=TRUE)
    for (c in seq_len(ncol(df))) {
        col_dirpath <- col_dirpaths[[c]]
        if (!dir.create(col_dirpath, showWarnings=TRUE, mode="0775"))
            stop(wmsg("failed to create directory: ", col_dirpath))
    }
    .write_object(integer(0), dirpath, "breakpoints", overwrite=TRUE)
    if (!is.null(spatial_index)) {
        .write_object(spatial_index[0], dirpath, "spatial_index",
                      overwrite=TRUE)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .append_df_to_OnDiskLongTable()
###

.write_OnDiskLongTable_column <- function(df_col, dirpath,
                                          breakpoints, b_offset, c,
                                          compress=TRUE)
{
    batch_ranges <- as(PartitioningByEnd(breakpoints), "IRanges")
    for (b2 in seq_along(batch_ranges)) {
        b <- b_offset + b2
        block_data <- extractROWS(df_col, batch_ranges[b2])
        .write_OnDiskLongTable_block(block_data, dirpath, b, c,
                                     compress=compress)
    }
}

.append_df_to_OnDiskLongTable <- function(df, dirpath,
                                          batchsize, spatial_index,
                                          compress=TRUE)
{
    ## Check existence of rowids.rds file.
    filepath <- file.path(dirpath, "rowids.rds")
    if (file.exists(filepath))
        stop(wmsg("cannot append data to an ",
                  "OnDiskLongTable object with row ids"))

    ## Check 'df' colnames.
    old_colnames <- names(.read_object(dirpath, "header"))
    if (!identical(colnames(df), old_colnames))
        stop(wmsg("'colnames(df)' do not match the colnames found in ",
                  file.path(dirpath, "header.rds")))

    ## Load 'breakpoints' vector.
    breakpoints0 <- .read_object(dirpath, "breakpoints")

    ## Append 'df' columns.
    col_dirpaths <- file.path(dirpath, .col_physname(seq_len(ncol(df))))
    if (is.null(spatial_index)) {
        breakpoints <- end(breakInChunks(nrow(df), batchsize))
    } else {
        breakpoints <- cumsum(batchsize)
        spatial_index0 <- .read_object(dirpath, "spatial_index")
    }
    b_offset <- length(breakpoints0)
    for (c in seq_len(ncol(df))) {
        col_dirpath <- col_dirpaths[[c]]
        if (!dir.exists(col_dirpath))
            stop(wmsg("directory not found: ", col_dirpath))
        .write_OnDiskLongTable_column(df[[c]], dirpath,
                                      breakpoints, b_offset, c,
                                      compress=compress)
    }

    ## Update 'breakpoints' vector.
    nrow0 <- .get_OnDiskLongTable_nrow_from_breakpoints(breakpoints0)
    breakpoints <- c(breakpoints0, nrow0 + breakpoints)
    .write_object(breakpoints, dirpath, "breakpoints", overwrite=TRUE)

    if (!is.null(spatial_index)) {
        ## Update 'spatial_index'.
        ## We use suppressWarnings() to avoid the "The 2 combined objects have
        ## no sequence levels in common" warning that would typically happen
        ## when writeOnDiskLongTable() is used in a loop to append data from
        ## different chromosomes (unless at each iteration the user is
        ## cautiously building and passing a 'spatial_index' with seqlevels
        ## set to all the chromosomes).
        spatial_index <- suppressWarnings(c(spatial_index0, spatial_index))
        .write_object(spatial_index, dirpath, "spatial_index", overwrite=TRUE)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeOnDiskLongTable() and writeOnDiskLongTableRowids()
###

.normarg_batchsize <- function(batchsize, totalsize)
{
    if (!isSingleNumberOrNA(batchsize))
        stop(wmsg("'batchsize' must be a single integer or NA"))
    if (!is.integer(batchsize)) 
        batchsize <- as.integer(batchsize)
    if (is.na(batchsize))
        return(totalsize)
    if (batchsize < 0L) 
        stop(wmsg("'batchsize' cannot be negative"))
    if (batchsize == 0L && totalsize != 0L)
        stop(wmsg("'batchsize' can be 0 only if 'nrow(df)' is 0"))
    batchsize
}

### Ignore the row names on 'df'. If 'df' is a DataFrame, 'metadata(df)' and
### 'mcols(df)' are also ignored.
writeOnDiskLongTable <- function(df, dirpath=".",
                                 batchsize=NA, spatial_index=NULL,
                                 append=FALSE, compress=TRUE)
{
    if (!is.data.frame(df) && !is(df, "DataFrame")) 
        stop(wmsg("'df' must be a data.frame or DataFrame object"))
    .check_OnDiskLongTable_dirpath(dirpath)
    if (!isTRUEorFALSE(append))
        stop(wmsg("'append' must be TRUE or FALSE"))
    if (is.null(spatial_index)) {
        batchsize <- .normarg_batchsize(batchsize, nrow(df))
    } else {
        if (!identical(batchsize, NA))
            stop(wmsg("'batchsize' must be NA when 'spatial_index' ",
                      "is supplied"))
        if (!is(spatial_index, "GenomicRanges"))
            stop(wmsg("'spatial_index' must be a GenomicRanges object"))
        batchsize <- mcols(spatial_index)$batchsize
        if (is.null(batchsize))
            stop(wmsg("'spatial_index' must have a \"batchsize\" ",
                      "metadata column"))
        if (!is.numeric(batchsize))
            stop(wmsg("'spatial_index' metadata column \"batchsize\" ",
                      "must be numeric"))
        if (!is.integer(batchsize))
            batchsize <- as.integer(batchsize)
        if (!isTRUE(all(batchsize >= 1L)))
            stop(wmsg("'spatial_index' metadata column \"batchsize\" ",
                      "must contain positive numbers"))
        if (sum(batchsize) != nrow(df))
            stop(wmsg("sum of values in 'spatial_index' metadata ",
                      "column \"batchsize\" must equal 'nrow(df)'"))
        spatial_index <- unstrand(granges(spatial_index, use.names=FALSE))
    }
    if (!append)
        .write_zero_row_OnDiskLongTable(df, dirpath, spatial_index)
    .append_df_to_OnDiskLongTable(df, dirpath, batchsize, spatial_index,
                                  compress=TRUE)
}

### After row ids are saved, data cannot be appended to the OnDiskLongTable
### object.
writeOnDiskLongTableRowids <- function(rowids, dirpath=".", compress=TRUE)
{
    .check_OnDiskLongTable_rowids(rowids)
    .check_OnDiskLongTable_dirpath(dirpath)
    breakpoints <- .read_object(dirpath, "breakpoints")
    nrow <- .get_OnDiskLongTable_nrow_from_breakpoints(breakpoints)
    if (length(rowids) != nrow)
        stop(wmsg("length of 'rowids' is incompatible ",
                  "with the OnDiskLongTable object in ", dirpath))
    .write_object(rowids, dirpath, "rowids", compress=compress)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.friendly_read_object <- function(dirpath, opath)
{
    object <- try(.read_object(dirpath, opath), silent=TRUE)
    if (inherits(object, "try-error"))
        stop(wmsg("Cannot open ", .make_filepath(dirpath, opath), ". ",
                  "Please make sure that '", dirpath, "' is the path to a ",
                  "valid OnDiskLongTable directory structure, that is, one ",
                  "as written by writeOnDiskLongTable()."))
    object
}

OnDiskLongTable <- function(dirpath=".")
{
    .check_OnDiskLongTable_dirpath(dirpath)
    header <- .friendly_read_object(dirpath, "header")
    breakpoints <- .friendly_read_object(dirpath, "breakpoints")
    spatial_index <- try(.read_object(dirpath, "spatial_index"), silent=TRUE)
    if (inherits(spatial_index, "try-error")) {
        spatial_index <- NULL
    } else if (!is(spatial_index, "GRanges")) {
        stop(wmsg("invalid 'spatial_index'"))
    }
    .rowids_cache <- new.env(parent=emptyenv())
    new2("OnDiskLongTable", dirpath=dirpath,
                            header=header,
                            breakpoints=breakpoints,
                            spatial_index=spatial_index,
                            .rowids_cache=.rowids_cache)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Basic getters (from the matrix/data.frame API)
###

setMethod("dim", "OnDiskLongTable",
    function(x)
    {
        nrow <- .get_OnDiskLongTable_nrow_from_breakpoints(x@breakpoints)
        ncol <- length(x@header)
        c(nrow, ncol)
    }
)

setMethod("dimnames", "OnDiskLongTable",
    function(x) list(NULL, names(x@header))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method
###

setMethod("show", "OnDiskLongTable",
    function(object)
    {
        cat(class(object), " object with ",
            nrow(object), " row", if (nrow(object) >= 2L) "s" else "",
            " and ",
            ncol(object), " column", if (ncol(object) >= 2L) "s" else "",
            "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### OnDiskLongTable specific getters
###

### Generic defined in OnDiskLongTable_old-class.R
#setGeneric("breakpoints", function(x) standardGeneric("breakpoints"))

setMethod("breakpoints", "OnDiskLongTable", function(x) x@breakpoints)

setGeneric("batchsizes", function(x) standardGeneric("batchsizes"))

setMethod("batchsizes", "OnDiskLongTable",
    function(x)
    {
        x_breakpoints <- breakpoints(x)
        ans <- S4Vectors:::diffWithInitialZero(x_breakpoints)
        names(ans) <- names(x_breakpoints)
        ans
    }
)

setGeneric("spatialIndex", function(x) standardGeneric("spatialIndex"))

setMethod("spatialIndex", "OnDiskLongTable", function(x) x@spatial_index)

### Generic defined in OnDiskLongTable_old-class.R
#setGeneric("rowids", function(x) standardGeneric("rowids"))

### Return NULL or an integer vector with no NAs or duplicated values.
setMethod("rowids", "OnDiskLongTable",
    function(x)
    {
        objname <- "rowids"

        ## 1. Try to get the row ids from the cache.
        ans <- try(get(objname, envir=x@.rowids_cache, inherits=FALSE),
                   silent=TRUE)
        if (!inherits(ans, "try-error"))
            return(ans)

        filepath <- .make_filepath(x@dirpath, objname)
        if (file.exists(filepath)) {
            ## 2.a. Load the row ids from disk.
            ans <- readRDS(filepath)
            .check_OnDiskLongTable_rowids(ans)
            if (length(ans) != nrow(x))
                stop(wmsg("length(rowids) != nrow(x)"))
        } else {
            ## 2.b. No row ids.
            ans <- NULL
        }

        ## 3. Cache the row ids.
        assign(objname, ans, envir=x@.rowids_cache)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getRowsByIndexFromOnDiskLongTable()
###

.normarg_rowidx <- function(rowidx, x)
{
    if (!is.integer(rowidx))
        stop(wmsg("'rowidx' must be an integer vector"))
    if (S4Vectors:::anyMissingOrOutside(rowidx, 1L, nrow(x)))
        stop(wmsg("'rowidx' contains NAs or invalid row indices"))
    rowidx
}

.normarg_colidx <- function(colidx, x)
{
    if (is.character(colidx)) {
        colidx <- match(colidx, colnames(x))
        if (S4Vectors:::anyMissing(colidx))
            stop(wmsg("'colidx' contains invalid column names"))
        return(colidx)
    }
    if (!is.numeric(colidx))
        stop(wmsg("'colidx' must be an integer or character vector"))
    if (!is.integer(colidx))
        colidx <- as.integer(colidx)
    if (S4Vectors:::anyMissingOrOutside(colidx, 1L, ncol(x)))
        stop(wmsg("'colidx' contains invalid column indices"))
    colidx
}

.breakpoints2offsets <- function(breakpoints)
{
    breakpoints_len <- length(breakpoints)
    if (breakpoints_len == 0L)
        return(integer(0))
    c(0L, breakpoints[-breakpoints_len])
}

### Translation from row index to "row key".
###   breakpoints: integer vector of break points.
###   rowidx:      integer vector containing valid row indices.
### Return a list of 2 integer vectors parallel to 'rowidx'. The 1st vector
### contains batch indices. The 2nd vector contains row indices relative to
### the batch indicated by the corresponding element in the 1st vector.
### In addition, if 'breakpoints' has names (i.e. batch names) then the 1st
### vector also has names indicating the batch associated with each row index
### in 'rowidx'.
.rowidx2rowkeys <- function(breakpoints, rowidx)
{
    batchidx <- findInterval(rowidx - 1L, breakpoints) + 1L
    names(batchidx) <- names(breakpoints)[batchidx]
    rowrelidx <- rowidx - .breakpoints2offsets(unname(breakpoints))[batchidx]
    list(batchidx, rowrelidx)
}

### c: column index of length 1
### rowkeys: list of 2 integer vectors as returned by .rowidx2rowkeys()
.read_OnDiskLongTable_column <- function(x, c, rowkeys)
{
    if (length(rowkeys[[1L]]) == 0L)
        return(x@header[[c]])
    list_of_keys <- split(unname(rowkeys[[2L]]), rowkeys[[1L]])
    tmp <- lapply(seq_along(list_of_keys),
             function(i) {
               b <- as.integer(names(list_of_keys)[[i]])
               block_data <- .read_OnDiskLongTable_block(x@dirpath, b, c)
               block_data[list_of_keys[[i]]]
             })
    S4Vectors:::quick_unsplit(tmp, rowkeys[[1L]])
}

### rowidx: integer vector.
### colidx: integer or character vector.
### Return a DataFrame (or data.frame) with 1 row per row id in 'rowidx'.
### Note that we do NOT set the row names to 'rowidx' on the returned DataFrame
### because we want to support duplicates in 'rowidx'.
getRowsByIndexFromOnDiskLongTable <- function(x, rowidx, colidx,
                                              as.data.frame=FALSE)
{
    if (!is(x, "OnDiskLongTable"))
        stop(wmsg("'x' must be an OnDiskLongTable object"))
    rowidx <- .normarg_rowidx(rowidx, x)
    x_breakpoints <- breakpoints(x)
    rowkeys <- .rowidx2rowkeys(x_breakpoints, rowidx)
    colidx <- .normarg_colidx(colidx, x)
    if (!isTRUEorFALSE(as.data.frame))
        stop(wmsg("'as.data.frame' must be TRUE or FALSE"))
    names(colidx) <- colnames(x)[colidx]
    ans_listData <-
        lapply(colidx, function(c) .read_OnDiskLongTable_column(x, c, rowkeys))
    if (as.data.frame) {
        ans <- data.frame(ans_listData, stringsAsFactors=FALSE)
    } else {
        ans <- new("DataFrame", listData=ans_listData,
                                nrows=length(rowidx))
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getRowsByIdFromOnDiskLongTable()
###

.rowids2rowidx <- function(x, rowids)
{
    x_rowids <- rowids(x)
    if (is.null(x_rowids))
        stop(wmsg("'x' has no row ids: cannot use ",
                  "getRowsByIdFromOnDiskLongTable() on it"))
    if (!is.integer(rowids) || S4Vectors:::anyMissing(rowids))
        stop(wmsg("'rowids' must be an integer vector with no NAs"))
    rowidx <- match(rowids, x_rowids)
    if (S4Vectors:::anyMissing(rowidx))
        stop(wmsg("'rowids' contains invalid row ids"))
    rowidx
}

### rowids: integer vector.
### colidx: integer or character vector.
### Return a DataFrame (or data.frame) with 1 row per row id in 'rowids'.
### Note that we do NOT set the row names to 'rowids' on the returned DataFrame
### because we want to support duplicates in 'rowids'.
getRowsByIdFromOnDiskLongTable <- function(x, rowids, colidx,
                                              as.data.frame=FALSE)
{
    if (!is(x, "OnDiskLongTable"))
        stop(wmsg("'x' must be an OnDiskLongTable object"))
    rowidx <- .rowids2rowidx(x, rowids)
    getRowsByIndexFromOnDiskLongTable(x, rowidx, colidx,
                                         as.data.frame=as.data.frame)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getBatchesByOverlapsFromOnDiskLongTable()
###

### c: column index of length 1
.read_OnDiskLongTable_batch_column <- function(x, batchidx, c)
{
    if (length(batchidx) == 0L)
        return(x@header[[c]])
    tmp <- lapply(batchidx,
        function(b) .read_OnDiskLongTable_block(x@dirpath, b, c)
    )
    S4Vectors:::quick_unlist(tmp)
}

### seqnames: factor-Rle.
### rowids: integer vector or NULL.
.make_DataFrame_or_data_frame <- function(seqnames,
                                          listData=list(), rowids=NULL,
                                          as.data.frame=FALSE)
{
    if (as.data.frame) {
        seqnames <- S4Vectors:::decodeRle(seqnames)
        ans_listData <- c(list(seqnames=seqnames), listData)
        ans <- data.frame(ans_listData, row.names=rowids,
                          stringsAsFactors=FALSE)
        return(ans)
    }
    ans_listData <- c(list(seqnames=seqnames), listData)
    ## Unfortunately, DataFrame cannot store its row names as an integer
    ## vector.
    if (!is.null(rowids))
        rowids <- as.character(rowids)
    new("DataFrame", listData=ans_listData,
                     nrows=length(seqnames),
                     rownames=rowids)
}

### ranges: GenomicRanges object.
### colidx: integer or character vector.
getBatchesByOverlapsFromOnDiskLongTable <- function(x, ranges,
                                                    maxgap=0L, minoverlap=1L,
                                                    colidx,
                                                    with.rowids=FALSE,
                                                    as.data.frame=FALSE)
{
    if (!is(x, "OnDiskLongTable"))
        stop(wmsg("'x' must be an OnDiskLongTable object"))
    colidx <- .normarg_colidx(colidx, x)
    if (!isTRUEorFALSE(with.rowids))
        stop(wmsg("'with.rowids' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.data.frame))
        stop(wmsg("'as.data.frame' must be TRUE or FALSE"))
    x_spatial_index <- spatialIndex(x)
    if (is.null(x_spatial_index))
        stop(wmsg("'x' has no spatial_index"))
    hits <- findOverlaps(ranges, x_spatial_index,
                         maxgap=maxgap, minoverlap=minoverlap)
    batchidx <- sort(unique(subjectHits(hits)))
    names(colidx) <- colnames(x)[colidx]
    ans_listData <- lapply(colidx,
        function(c) .read_OnDiskLongTable_batch_column(x, batchidx, c)
    )
    x_batchsizes <- batchsizes(x)
    ans_batchsizes <- x_batchsizes[batchidx]
    ans_seqnames <- rep.int(seqnames(x_spatial_index)[batchidx],
                            ans_batchsizes)
    ans_rowids <- NULL
    if (with.rowids) {
        x_rowids <- rowids(x)
        if (!is.null(x_rowids)) {
            rowidx <- as.integer(successiveIRanges(x_batchsizes)[blockidx])
            ans_rowids <- x_rowids[rowidx]
        }
    }
    .make_DataFrame_or_data_frame(ans_seqnames, ans_listData, ans_rowids,
                                  as.data.frame)
}

