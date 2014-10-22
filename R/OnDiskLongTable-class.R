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
                                     # user-supplied row ids to row numbers.
    )
)

.get_nrow_from_breakpoints <- function(breakpoints)
{
    if (length(breakpoints) == 0L)
        return(0L)
    breakpoints[[length(breakpoints)]]
}

.colidx2dirname <- function(colidx) sprintf("col%03d", colidx)

.blockidx2blockname <- function(blockidx)
{
    not_na_idx <- which(!is.na(blockidx))
    blockname <- character(length(blockidx))
    blockname[not_na_idx] <- sprintf("b%05d", blockidx[not_na_idx])
    blockname
}

.save_object <- function(dirpath, objname, value,
                         compress=TRUE, compression_level)
{
    filename <- paste0(objname, ".rda")
    filepath <- file.path(dirpath, filename)
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
### Writing data in the "OnDiskLongTable format"
###

.save_zero_length_columns <- function(df, dirpath)
{
    col_dirnames <- .colidx2dirname(seq_len(ncol(df)))
    col_dirpaths <- file.path(dirpath, col_dirnames)
    fail_idx <- which(as.logical(unlink(col_dirpaths, recursive=TRUE)))
    if (length(fail_idx) != 0L)
        stop("failed to remove directories: ",
             paste0(col_dirpaths[fail_idx], collapse=", "))
    for (j in seq_len(ncol(df))) {
        col_dirpath <- col_dirpaths[j]
        if (!dir.create(col_dirpath, showWarnings=TRUE, mode="0775"))
            stop("failed to create directory: ", col_dirpath)
        ## Save block 0 (special block containing zero-length column).
        blockname <- .blockidx2blockname(0L)
        block <- df[[j]][integer(0)]
        .save_object(col_dirpath, blockname, block)
    }
}

### Return new break points.
.save_blocks <- function(df, dirpath, blocksize,
                         compress, compression_level,
                         prev_breakpoints=integer(0))
{
    col_dirnames <- .colidx2dirname(seq_len(ncol(df)))
    col_dirpaths <- file.path(dirpath, col_dirnames)
    block_ranges <- breakInChunks(nrow(df), blocksize)
    blockidx_offset <- length(prev_breakpoints)
    for (j in seq_len(ncol(df))) {
        col_dirpath <- col_dirpaths[j]
        if (!dir.exists(col_dirpath))
            stop("directory not found: ", col_dirpath)
        for (b in seq_along(block_ranges)) {
            blockidx <- blockidx_offset + b
            blockname <- .blockidx2blockname(blockidx)
            block <- df[[j]][block_ranges[[b]]]
            .save_object(col_dirpath, blockname, block,
                         compress, compression_level)
        }
    }
    prev_nrow <- .get_nrow_from_breakpoints(prev_breakpoints)
    c(prev_breakpoints, prev_nrow + end(block_ranges))
}

.save_new_OnDiskLongTable <- function(df, dirpath, blocksize,
                                      compress, compression_level)
{
    if (!dir.exists(dirpath))
        stop("directory not found: ", dirpath)

    ## Remove rowids.rda file.
    filepath <- file.path(dirpath, "rowids.rda")
    if (file.exists(filepath)) {
        if (!file.remove(filepath))
            stop("failed to remove file: ", filepath)
    }

    ## Save 'df' colnames.
    .save_object(dirpath, "colnames", colnames(df))

    ## Save 'df' columns.
    .save_zero_length_columns(df, dirpath)
    breakpoints <- .save_blocks(df, dirpath, blocksize,
                                compress, compression_level)

    ## Save 'breakpoints' vector.
    .save_object(dirpath, "breakpoints", breakpoints)
}

.append_to_OnDiskLongTable <- function(df, dirpath, blocksize,
                                       compress, compression_level)
{
    if (!dir.exists(dirpath))
        stop("directory not found: ", dirpath)

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
    prev_breakpoints <- .load_object(dirpath, "breakpoints")

    ## Append 'df' columns.
    breakpoints <- .save_blocks(df, dirpath, blocksize,
                                compress, compression_level,
                                prev_breakpoints)

    ## Update 'breakpoints' vector.
    .save_object(dirpath, "breakpoints", breakpoints)
}

### Ignore the row names on 'df'. If 'df' is a DataFrame, 'metadata(df)' and
### 'mcols(df)' are also ignored.
saveAsOnDiskLongTable <- function(df, dirpath=".", append=FALSE,
                                  blocksize=NA,
                                  compress=TRUE, compression_level)
{
    if (!is.data.frame(df) && !is(df, "DataFrame")) 
        stop("'df' must be a data.frame or DataFrame object")
    if (!isSingleString(dirpath)) 
        stop("'dirpath' must be a single string")
    if (!isTRUEorFALSE(append))
        stop("'append' must be TRUE or FALSE")
    if (!isSingleNumberOrNA(blocksize))
        stop("'blocksize' must be a single integer or NA")
    if (!is.integer(blocksize)) 
        blocksize <- as.integer(blocksize)
    if (is.na(blocksize))
        blocksize <- nrow(df)
    else if (blocksize < 0L) 
        stop("'blocksize' cannot be negative")
    else if (blocksize == 0L && nrow(df) != 0L)
        stop("'blocksize' can be 0 only if 'nrow(df)' is 0")
    if (!append)
        .save_new_OnDiskLongTable(df, dirpath, blocksize,
                                  compress, compression_level)
    else
        .append_to_OnDiskLongTable(df, dirpath, blocksize,
                                   compress, compression_level)
}

.check_rowids <- function(rowids)
{
    if (!is.integer(rowids))
        stop("'rowids' must be an integer vector")
    if (any(is.na(rowids)) || any(duplicated(rowids)))
        stop("'rowids' cannot contain NAs or duplicated values")
}

### After row ids are saved, data cannot be appended to the OnDiskLongTable
### object.
saveRowidsForOnDiskLongTable <- function(rowids, dirpath=".",
                                         compress=TRUE, compression_level)
{
    .check_rowids(rowids)
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
    if (!isSingleString(dirpath) || dirpath == "")
        stop("'dirpath' must be a single (non-empty) string")
    if (!dir.exists(dirpath))
        stop("directory not found: ", dirpath)
    ans_colnames <- .load_object(dirpath, "colnames")
    ans_breakpoints <- .load_object(dirpath, "breakpoints")
    ans_rowids_cache <- new.env(parent=emptyenv())
    ans <- new("OnDiskLongTable", dirpath=dirpath,
                                  colnames=ans_colnames,
                                  breakpoints=ans_breakpoints,
                                  .rowids_cache=ans_rowids_cache)
    if (any(is.na(ans_breakpoints))
     || is.unsorted(ans_breakpoints, strictly=TRUE)
     || (length(ans_breakpoints) != 0L && ans_breakpoints[[1L]] < 1L))
        stop("invalid breakpoints found for OnDiskLongTable object")
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "nrow", "rowids", "colnames", and "ncol" methods
###

setMethod("nrow", "OnDiskLongTable",
    function(x) .get_nrow_from_breakpoints(x@breakpoints)
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
    col_dirname <- .colidx2dirname(colidx)
    col_dirpath <- file.path(x@dirpath, col_dirname)
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
        stop("'colidx' contains invalid column numbers")
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
        stop("'rowids' contains invalid row numbers")
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
###   x:      OnDiskLongTable object.
###   rowidx: integer vector containing valid row numbers (NAs allowed).
### Return a named integer vector "parallel" to 'rowidx', i.e. each
### (name, value) pair in the returned vector represents the "row key" of
### the corresponding row index. The name is the "block name" where to find
### the row and the value is the row number *relative to the block*, that is,
### its position within the block.
.rowidx2rowkeys <- function(x, rowidx)
{
    blockidx <- findInterval(rowidx - 1L, x@breakpoints) + 1L
    rowkeys <- rowidx - .breakpoints2offsets(x@breakpoints)[blockidx]
    names(rowkeys) <- .blockidx2blockname(blockidx)
    rowkeys
}

### colidx: column index of length 1
.get_values_from_column <- function(x, colidx, rowkeys)
{
    if (length(rowkeys) == 0L)
        return(.load_block(x, colidx, 0L))
    list_of_keys <- split(unname(rowkeys), names(rowkeys))
    tmp <- lapply(seq_along(list_of_keys),
             function(i) {
               blockname <- names(list_of_keys)[i]
               keys <- list_of_keys[[i]]
               block <- .load_block(x, colidx, blockname)
               block[keys]
             })
    quickUnsplit(tmp, names(rowkeys))
}

### rowids: integer vector.
### colidx: integer or character vector.
### Return a DataFrame.
getDataFromOnDiskLongTable <- function(x, rowids, colidx)
{
    if (!is(x, "OnDiskLongTable"))
        stop("'x' must be an OnDiskLongTable object")
    rowidx <- .rowids2rowidx(x, rowids)
    rowkeys <- .rowidx2rowkeys(x, rowidx)
    colidx <- .normarg_colidx(colidx, x)
    names(colidx) <- colnames(x)[colidx]
    ans_listData <-
        lapply(colidx, function(j) .get_values_from_column(x, j, rowkeys))
    new("DataFrame", listData=ans_listData, nrows=length(rowids),
                     rownames=as.character(rowids))
}

