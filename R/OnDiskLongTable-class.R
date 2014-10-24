### =========================================================================
### OnDiskLongTable objects
### -------------------------------------------------------------------------
###

setClass("OnDiskLongTable",
    representation(
        dirpath="character",         # Single string.
        colnames="character",        # Vector of column names (no NAs, no
                                     # empty strings).
        breakpoints="integer",       # Values must be positive and strictly
                                     # sorted. Last value matches nb of rows
                                     # in the table.
        .rowids_cache="environment"  # Where to load and cache the row ids.
                                     # Row ids are unique integer values.
                                     # Used for fast translation from
                                     # user-supplied row ids to row indices.
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers
###

.check_dirpath <- function(dirpath)
{
    if (!isSingleString(dirpath) || dirpath == "")
        stop("'dirpath' must be a single (non-empty) string")
    if (!dir.exists(dirpath))
        stop("directory not found: ", dirpath)
}

.check_breakpoints <- function(breakpoints)
{
    if (!is.integer(breakpoints)
     || any(is.na(breakpoints))
     || is.unsorted(breakpoints, strictly=TRUE)
     || length(breakpoints) != 0L && breakpoints[[1L]] < 1L)
        stop("invalid breakpoints found for OnDiskLongTable object")
}

.check_rowids <- function(rowids)
{
    if (!is.integer(rowids))
        stop("'rowids' must be an integer vector")
    if (any(is.na(rowids)) || any(duplicated(rowids)))
        stop("'rowids' cannot contain NAs or duplicated values")
}

.get_nrow_from_breakpoints <- function(breakpoints)
{
    if (length(breakpoints) == 0L)
        return(0L)
    breakpoints[[length(breakpoints)]]
}

.colidx2dirname <- function(colidx)
    sprintf("col%03d", colidx)

.colidx2dirpath <- function(dirpath, colidx)
    file.path(dirpath, .colidx2dirname(colidx))

.blockidx2blockname <- function(blockidx)
{
    not_na_idx <- which(!is.na(blockidx))
    blockname <- character(length(blockidx))
    blockname[not_na_idx] <- sprintf("b%05d", blockidx[not_na_idx])
    blockname
}

.save_object <- function(dirpath, objname, value,
                         compress=TRUE, compression_level,
                         overwrite=FALSE)
{
    filename <- paste0(objname, ".rda")
    filepath <- file.path(dirpath, filename)
    if (!overwrite && file.exists(filepath))
        stop("file already exists: ", filepath)
    assign(objname, value)
    save(list=objname, file=filepath,
         compress=compress, compression_level=compression_level)
}

.load_object <- function(dirpath, objname)
{
    filename <- paste0(objname, ".rda")
    filepath <- file.path(dirpath, filename)
    if (!file.exists(filepath))
        stop("file not found: ", filepath)
    tmpenv <- new.env(parent=emptyenv())
    if (!identical(load(filepath, envir=tmpenv), objname))
        stop(filepath, " does not seem to belong ",
             "to an OnDiskLongTable object")
    get(objname, envir=tmpenv, inherits=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .write_zero_row_OnDiskLongTable()
###

.remove_file <- function(dirpath, filename)
{
    filepath <- file.path(dirpath, filename)
    if (file.exists(filepath)) {
        if (!file.remove(filepath))
            stop("failed to remove file: ", filepath)
    }
}

.remove_dirs <- function(dirpaths)
{
    fail_idx <- which(as.logical(unlink(dirpaths, recursive=TRUE)))
    if (length(fail_idx) != 0L)
        stop("failed to remove directories: ",
             paste0(dirpaths[fail_idx], collapse=", "))
}

.create_and_initialize_column_dir <- function(df_col, col_dirpath)
{
    if (!dir.create(col_dirpath, showWarnings=TRUE, mode="0775"))
        stop("failed to create directory: ", col_dirpath)
    ## A "special block" is written in newly created column dir. It contains
    ## serialized zero-length column.
    blockname <- .blockidx2blockname(0L)
    block <- df_col[integer(0)]
    .save_object(col_dirpath, blockname, block)
}

.write_zero_row_OnDiskLongTable <- function(df, dirpath)
{
    ## Remove stuff.
    .remove_file(dirpath, "rowids.rda")
    col_dirpaths <- .colidx2dirpath(dirpath, seq_len(ncol(df)))
    .remove_dirs(col_dirpaths)

    ## Create stuff.
    .save_object(dirpath, "colnames", colnames(df), overwrite=TRUE)
    for (j in seq_len(ncol(df)))
        .create_and_initialize_column_dir(df[[j]], col_dirpaths[[j]])
    .save_object(dirpath, "breakpoints", integer(0), overwrite=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .append_batch_to_OnDiskLongTable()
###

.write_column <- function(df_col, col_dirpath,
                          breakpoints, blockidx_offset,
                          compress, compression_level)
{
    if (!dir.exists(col_dirpath))
        stop("directory not found: ", col_dirpath)
    block_ranges <- PartitioningByEnd(breakpoints)
    for (b in seq_along(block_ranges)) {
        blockidx <- blockidx_offset + b
        blockname <- .blockidx2blockname(blockidx)
        block <- df_col[block_ranges[[b]]]
        .save_object(col_dirpath, blockname, block,
                     compress, compression_level)
    }
}

.append_batch_to_OnDiskLongTable <- function(df, dirpath,
                                             batch_label, blocksize,
                                             compress, compression_level)
{
    ## Check existence of rowids.rda file.
    filepath <- file.path(dirpath, "rowids.rda")
    if (file.exists(filepath))
        stop("cannot append data to an OnDiskLongTable object with row ids")

    ## Check 'df' colnames.
    old_colnames <- .load_object(dirpath, "colnames")
    if (!identical(colnames(df), old_colnames))
        stop("'colnames(df)' must be identical to colnames in ",
             file.path(dirpath, "colnames.rda"))

    ## Load 'breakpoints' vector.
    breakpoints1 <- .load_object(dirpath, "breakpoints")

    ## Append 'df' columns.
    col_dirpaths <- .colidx2dirpath(dirpath, seq_len(ncol(df)))
    breakpoints2 <- end(breakInChunks(nrow(df), blocksize))
    blockidx_offset <- length(breakpoints1)
    for (j in seq_len(ncol(df)))
        .write_column(df[[j]], col_dirpaths[[j]],
                      breakpoints2, blockidx_offset,
                      compress, compression_level)

    ## Update 'breakpoints' vector.
    if (!is.null(batch_label))
        names(breakpoints2) <- rep.int(batch_label, length(breakpoints2))
    nrow1 <- .get_nrow_from_breakpoints(breakpoints1)
    breakpoints <- c(breakpoints1, nrow1 + breakpoints2)
    .save_object(dirpath, "breakpoints", breakpoints, overwrite=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### saveAsOnDiskLongTable() and saveRowidsForOnDiskLongTable()
###

.normarg_blocksize <- function(blocksize, totalsize)
{
    if (!isSingleNumberOrNA(blocksize))
        stop("'blocksize' must be a single integer or NA")
    if (!is.integer(blocksize)) 
        blocksize <- as.integer(blocksize)
    if (is.na(blocksize))
        return(totalsize)
    if (blocksize < 0L) 
        stop("'blocksize' cannot be negative")
    if (blocksize == 0L && totalsize != 0L)
        stop("'blocksize' can be 0 only if 'nrow(df)' is 0")
    blocksize
}

### Ignore the row names on 'df'. If 'df' is a DataFrame, 'metadata(df)' and
### 'mcols(df)' are also ignored.
saveAsOnDiskLongTable <- function(df, dirpath=".", append=FALSE,
                                  batch_label=NULL, blocksize=NA,
                                  compress=TRUE, compression_level)
{
    if (!is.data.frame(df) && !is(df, "DataFrame")) 
        stop("'df' must be a data.frame or DataFrame object")
    .check_dirpath(dirpath)
    if (!isTRUEorFALSE(append))
        stop("'append' must be TRUE or FALSE")
    if (!(is.null(batch_label) ||
          isSingleString(batch_label) && batch_label != ""))
        stop("'batch_label' must be NULL or a single (non-empty) string")
    blocksize <- .normarg_blocksize(blocksize, nrow(df))

    if (!append)
        .write_zero_row_OnDiskLongTable(df, dirpath)
    .append_batch_to_OnDiskLongTable(df, dirpath,
                                     batch_label, blocksize,
                                     compress, compression_level)
}

### After row ids are saved, data cannot be appended to the OnDiskLongTable
### object.
saveRowidsForOnDiskLongTable <- function(rowids, dirpath=".",
                                         compress=TRUE, compression_level)
{
    .check_rowids(rowids)
    .check_dirpath(dirpath)
    breakpoints <- .load_object(dirpath, "breakpoints")
    nrow <- .get_nrow_from_breakpoints(breakpoints)
    if (length(rowids) != nrow)
        stop("length of 'rowids' is incompatible ",
             "with the OnDiskLongTable object in ", dirpath)
    .save_object(dirpath, "rowids", rowids,
                 compress, compression_level)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

OnDiskLongTable <- function(dirpath=".")
{
    .check_dirpath(dirpath)
    ans_colnames <- .load_object(dirpath, "colnames")
    ans_breakpoints <- .load_object(dirpath, "breakpoints")
    ans_rowids_cache <- new.env(parent=emptyenv())
    ans <- new("OnDiskLongTable", dirpath=dirpath,
                                  colnames=ans_colnames,
                                  breakpoints=ans_breakpoints,
                                  .rowids_cache=ans_rowids_cache)
    .check_breakpoints(ans_breakpoints)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "nrow", "rowids", "colnames", and "ncol" methods
###

setGeneric("breakpoints", function(x) standardGeneric("breakpoints"))

setMethod("breakpoints", "OnDiskLongTable", function(x) x@breakpoints)

setMethod("nrow", "OnDiskLongTable",
    function(x) .get_nrow_from_breakpoints(breakpoints(x))
)

setGeneric("rowids", function(x) standardGeneric("rowids"))

### Return NULL or an integer vector with no NAs or duplicated values.
setMethod("rowids", "OnDiskLongTable",
    function(x)
    {
        objname <- "rowids"

        ## 1. Try to get the row ids from the cache.
        ans <- try(get(objname, envir=x@.rowids_cache, inherits=FALSE),
                   silent=TRUE)
        if (!is(ans, "try-error"))
            return(ans)

        ## 2. Try to find them on disk.
        filename <- paste0(objname, ".rda")
        filepath <- file.path(x@dirpath, filename)
        if (!file.exists(filepath))
            return((assign(objname, NULL, envir=x@.rowids_cache)))

        ## 3. Load the row ids from disk.
        if (!identical(load(filepath, envir=x@.rowids_cache), objname))
            stop(filepath, " does not seem to belong ",
                 "to an OnDiskLongTable object")
        ans <- get(objname, envir=x@.rowids_cache, inherits=FALSE)
        .check_rowids(ans)
        if (length(ans) != nrow(x))
            stop("length(rowids) != nrow(x)")
        ans
    }
)

setMethod("colnames", "OnDiskLongTable",
    function(x, do.NULL=TRUE, prefix="col") x@colnames
)

setMethod("ncol", "OnDiskLongTable", function(x) length(colnames(x)))


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
### Load a single block
###

### colidx:    column index of length 1.
### blockname: a valid block name. Can also be a valid block index as a
###            single integer.
.load_block <- function(x, colidx, blockname)
{
    if (!is.character(blockname))
        blockname <- .blockidx2blockname(blockname)
    col_dirpath <- .colidx2dirpath(x@dirpath, colidx)
    .load_object(col_dirpath, blockname)
}

.normarg_colidx <- function(colidx, x)
{
    if (is.character(colidx)) {
        colidx <- match(colidx, colnames(x))
        if (any(is.na(colidx)))
            stop("'colidx' contains invalid column names")
        return(colidx)
    }
    if (!is.numeric(colidx))
        stop("'colidx' must be an integer or character vector")
    if (!is.integer(colidx))
        colidx <- as.integer(colidx)
    if (S4Vectors:::anyMissingOrOutside(colidx, 1L, ncol(x)))
        stop("'colidx' contains invalid column indices")
    colidx
}

### User-friendly wrapper to .load_block(). NOT exported.
### colidx:   column index of length 1. Can be a single column name.
### blockidx: block index of length 1. Can be a single block name.
getBlockFromOnDiskLongTable <- function(x, colidx, blockidx)
{
    if (!is(x, "OnDiskLongTable"))
        stop("'x' must be an OnDiskLongTable object")
    if (!(isSingleNumber(colidx) || isSingleString(colidx)))
        stop("'colidx' must be a single integer or string")
    if (!(isSingleNumber(blockidx) || isSingleString(blockidx)))
        stop("'blockidx' must be a single integer or string")
    colidx <- .normarg_colidx(colidx, x)
    .load_block(x, colidx, blockidx)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get data from OnDiskLongTable
###

.rowids2rowidx <- function(x, rowids)
{
    if (!is.integer(rowids))
        stop("'rowids' must be an integer vector")
    x_rowids <- rowids(x)
    if (!is.null(x_rowids)) {
        rowidx <- match(rowids, x_rowids)
        if (any(is.na(rowidx)))
            stop("'rowids' contains invalid row ids")
        return(rowidx)
    }
    ## When 'rowids(x)' is NULL, then row ids are implicit and considered
    ## to be 'seq_len(nrow(x))'.
    if (S4Vectors:::anyMissingOrOutside(rowids, 1L, nrow(x)))
        stop("'rowids' contains invalid row indices")
    rowids
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
### contains block indices. The 2nd vector contains row indices relative to
### the block indicated by the corresponding element in the 1st vector.
### In addition, if 'breakpoints' has names (i.e. "batch labels") then the
### 1st vector also has names indicating the label of the batch associated
### with each row index in 'rowidx'.
.rowidx2rowkeys <- function(breakpoints, rowidx)
{
    blockidx <- findInterval(rowidx - 1L, breakpoints) + 1L
    names(blockidx) <- names(breakpoints)[blockidx]
    rowrelidx <- rowidx - .breakpoints2offsets(unname(breakpoints))[blockidx]
    list(blockidx, rowrelidx)
}

### colidx: column index of length 1
.get_values_from_column <- function(x, colidx, rowkeys)
{
    if (length(rowkeys[[1L]]) == 0L)
        return(.load_block(x, colidx, 0L))
    list_of_keys <- split(unname(rowkeys[[2L]]), rowkeys[[1L]])
    tmp <- lapply(seq_along(list_of_keys),
             function(i) {
               blockidx <- as.integer(names(list_of_keys)[i])
               keys <- list_of_keys[[i]]
               block <- .load_block(x, colidx, blockidx)
               block[keys]
             })
    quickUnsplit(tmp, rowkeys[[1L]])
}

### rowids: integer vector.
### colidx: integer or character vector.
### Return a DataFrame or data.frame.
getDataFromOnDiskLongTable <- function(x, rowids, colidx,
                                       with.batch_label=FALSE,
                                       as.data.frame=FALSE)
{
    if (!is(x, "OnDiskLongTable"))
        stop("'x' must be an OnDiskLongTable object")
    rowidx <- .rowids2rowidx(x, rowids)
    x_breakpoints <- breakpoints(x)
    rowkeys <- .rowidx2rowkeys(x_breakpoints, rowidx)
    colidx <- .normarg_colidx(colidx, x)
    if (!isTRUEorFALSE(with.batch_label))
        stop("'with.batch_label' must be TRUE or FALSE")
    if (!isTRUEorFALSE(as.data.frame))
        stop("'as.data.frame' must be TRUE or FALSE")
    names(colidx) <- colnames(x)[colidx]
    ans_listData <-
        lapply(colidx, function(j) .get_values_from_column(x, j, rowkeys))
    if (with.batch_label)
        ans_listData[["batch_label"]] <- factor(names(rowkeys[[1L]]),
                                                levels=names(x_breakpoints))
    if (as.data.frame) {
        ans <- data.frame(ans_listData, stringsAsFactors=FALSE)
    } else {
        ans <- new("DataFrame", listData=ans_listData,
                                nrows=length(rowids),
                                rownames=as.character(rowids))
    }
    ans
}

