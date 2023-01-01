### =========================================================================
### OnDiskLongTable_old objects
### -------------------------------------------------------------------------
###


setClass("OnDiskLongTable_old",
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

.OnDiskLongTable_old_check_dirpath <- function(dirpath)
{
    if (!isSingleString(dirpath) || dirpath == "")
        stop(wmsg("'dirpath' must be a single (non-empty) string"))
    if (!dir.exists(dirpath))
        stop(wmsg("directory not found: ", dirpath))
}

.OnDiskLongTable_old_check_breakpoints <- function(breakpoints)
{
    if (!is.integer(breakpoints)
     || S4Vectors:::anyMissing(breakpoints)
     || is.unsorted(breakpoints, strictly=TRUE)
     || length(breakpoints) != 0L && breakpoints[[1L]] < 1L)
        stop(wmsg("invalid breakpoints found for OnDiskLongTable_old object"))
}

.OnDiskLongTable_old_check_rowids <- function(rowids)
{
    if (!is.integer(rowids))
        stop(wmsg("'rowids' must be an integer vector"))
    if (S4Vectors:::anyMissing(rowids) || anyDuplicated(rowids))
        stop(wmsg("'rowids' cannot contain NAs or duplicated values"))
}

.OnDiskLongTable_old_get_nrow_from_breakpoints <- function(breakpoints)
{
    if (length(breakpoints) == 0L)
        return(0L)
    breakpoints[[length(breakpoints)]]
}

.OnDiskLongTable_old_colidx2dirname <- function(colidx)
    sprintf("col%03d", colidx)

.OnDiskLongTable_old_colidx2dirpath <- function(dirpath, colidx)
    file.path(dirpath, .OnDiskLongTable_old_colidx2dirname(colidx))

.OnDiskLongTable_old_blockidx2blockname <- function(blockidx)
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
        stop(wmsg("file already exists: ", filepath))
    assign(objname, value)
    save(list=objname, file=filepath,
         compress=compress, compression_level=compression_level)
}

.load_object <- function(dirpath, objname)
{
    filename <- paste0(objname, ".rda")
    filepath <- file.path(dirpath, filename)
    if (!file.exists(filepath))
        stop(wmsg("file not found: ", filepath))
    tmpenv <- new.env(parent=emptyenv())
    if (!identical(load(filepath, envir=tmpenv), objname))
        stop(wmsg(filepath, " does not seem to belong ",
                  "to an OnDiskLongTable_old object"))
    get(objname, envir=tmpenv, inherits=FALSE)
}

### colidx:    column index of length 1.
### blockname: a valid block name. Can also be a valid block index as a
###            single integer.
.load_block <- function(x, colidx, blockname)
{
    if (!is.character(blockname))
        blockname <- .OnDiskLongTable_old_blockidx2blockname(blockname)
    col_dirpath <- .OnDiskLongTable_old_colidx2dirpath(x@dirpath, colidx)
    .load_object(col_dirpath, blockname)
}

### Return an environment (x@.rowids_cache) that contains the row ids as
### one big integer vector.
get_rowids_env_old <- function(x)
{
    stopifnot(is(x, "OnDiskLongTable_old"))

    objname <- "rowids"

    ## 1. Try to get the row ids from the cache.
    tmp <- try(get(objname, envir=x@.rowids_cache, inherits=FALSE),
               silent=TRUE)
    if (!is(tmp, "try-error"))
        return(x@.rowids_cache)

    ## 2. Try to find them on disk.
    filename <- paste0(objname, ".rda")
    filepath <- file.path(x@dirpath, filename)
    if (!file.exists(filepath))
        return(x@.rowids_cache)

    ## 3. Load the row ids from disk.
    if (!identical(load(filepath, envir=x@.rowids_cache), objname))
        stop(wmsg(filepath, " does not seem to belong ",
                  "to an OnDiskLongTable_old object"))
    tmp <- get(objname, envir=x@.rowids_cache, inherits=FALSE)
    .OnDiskLongTable_old_check_rowids(tmp)
    if (length(tmp) != nrow(x))
        stop(wmsg("length(rowids) != nrow(x)"))
    x@.rowids_cache
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .write_zero_row_OnDiskLongTable_old()
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

.create_and_initialize_column_dir <- function(df_col, col_dirpath)
{
    if (!dir.create(col_dirpath, showWarnings=TRUE, mode="0775"))
        stop(wmsg("failed to create directory: ", col_dirpath))
    ## A "special block" is written in newly created column dir. It contains
    ## serialized zero-length column.
    blockname <- .OnDiskLongTable_old_blockidx2blockname(0L)
    block <- df_col[integer(0)]
    .save_object(col_dirpath, blockname, block)
}

.write_zero_row_OnDiskLongTable_old <- function(df, dirpath)
{
    ## Remove stuff.
    .remove_file(dirpath, "rowids.rda")
    col_dirpaths <- .OnDiskLongTable_old_colidx2dirpath(dirpath,
                                                        seq_len(ncol(df)))
    .remove_dirs(col_dirpaths)

    ## Create stuff.
    .save_object(dirpath, "colnames", colnames(df), overwrite=TRUE)
    for (j in seq_len(ncol(df)))
        .create_and_initialize_column_dir(df[[j]], col_dirpaths[[j]])
    .save_object(dirpath, "breakpoints", integer(0), overwrite=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .append_batch_to_OnDiskLongTable_old()
###

.write_column <- function(df_col, col_dirpath,
                          breakpoints, blockidx_offset,
                          compress, compression_level)
{
    if (!dir.exists(col_dirpath))
        stop(wmsg("directory not found: ", col_dirpath))
    block_ranges <- PartitioningByEnd(breakpoints)
    for (b in seq_along(block_ranges)) {
        blockidx <- blockidx_offset + b
        blockname <- .OnDiskLongTable_old_blockidx2blockname(blockidx)
        block <- df_col[block_ranges[[b]]]
        .save_object(col_dirpath, blockname, block,
                     compress, compression_level)
    }
}

.append_batch_to_OnDiskLongTable_old <- function(df, dirpath,
                                             batch_label, blocksize,
                                             compress, compression_level)
{
    ## Check existence of rowids.rda file.
    filepath <- file.path(dirpath, "rowids.rda")
    if (file.exists(filepath))
        stop(wmsg("cannot append data to an ",
                  "OnDiskLongTable_old object with row ids"))

    ## Check 'df' colnames.
    old_colnames <- .load_object(dirpath, "colnames")
    if (!identical(colnames(df), old_colnames))
        stop(wmsg("'colnames(df)' must be identical to colnames in ",
                  file.path(dirpath, "colnames.rda")))

    ## Load 'breakpoints' vector.
    breakpoints1 <- .load_object(dirpath, "breakpoints")

    ## Append 'df' columns.
    col_dirpaths <- .OnDiskLongTable_old_colidx2dirpath(dirpath,
                                                        seq_len(ncol(df)))
    breakpoints2 <- end(breakInChunks(nrow(df), chunksize=blocksize))
    blockidx_offset <- length(breakpoints1)
    for (j in seq_len(ncol(df)))
        .write_column(df[[j]], col_dirpaths[[j]],
                      breakpoints2, blockidx_offset,
                      compress, compression_level)

    ## Update 'breakpoints' vector.
    if (!is.null(batch_label))
        names(breakpoints2) <- rep.int(batch_label, length(breakpoints2))
    nrow1 <- .OnDiskLongTable_old_get_nrow_from_breakpoints(breakpoints1)
    breakpoints <- c(breakpoints1, nrow1 + breakpoints2)
    .save_object(dirpath, "breakpoints", breakpoints, overwrite=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### saveAsOnDiskLongTable_old() and saveRowidsForOnDiskLongTable_old()
###

.normarg_blocksize <- function(blocksize, totalsize)
{
    if (!isSingleNumberOrNA(blocksize))
        stop(wmsg("'blocksize' must be a single integer or NA"))
    if (!is.integer(blocksize)) 
        blocksize <- as.integer(blocksize)
    if (is.na(blocksize))
        return(totalsize)
    if (blocksize < 0L) 
        stop(wmsg("'blocksize' cannot be negative"))
    if (blocksize == 0L && totalsize != 0L)
        stop(wmsg("'blocksize' can be 0 only if 'nrow(df)' is 0"))
    blocksize
}

### Ignore the row names on 'df'. If 'df' is a DataFrame, 'metadata(df)' and
### 'mcols(df)' are also ignored.
saveAsOnDiskLongTable_old <- function(df, dirpath=".", append=FALSE,
                                      batch_label=NULL, blocksize=NA,
                                      compress=TRUE, compression_level)
{
    if (!is.data.frame(df) && !is(df, "DataFrame")) 
        stop(wmsg("'df' must be a data.frame or DataFrame object"))
    .OnDiskLongTable_old_check_dirpath(dirpath)
    if (!isTRUEorFALSE(append))
        stop(wmsg("'append' must be TRUE or FALSE"))
    if (!(is.null(batch_label) ||
          isSingleString(batch_label) && batch_label != ""))
        stop(wmsg("'batch_label' must be NULL or a single (non-empty) string"))
    blocksize <- .normarg_blocksize(blocksize, nrow(df))

    if (!append)
        .write_zero_row_OnDiskLongTable_old(df, dirpath)
    .append_batch_to_OnDiskLongTable_old(df, dirpath,
                                         batch_label, blocksize,
                                     compress, compression_level)
}

### After row ids are saved, data cannot be appended to the OnDiskLongTable_old
### object.
saveRowidsForOnDiskLongTable_old <- function(rowids, dirpath=".",
                                             compress=TRUE, compression_level)
{
    .OnDiskLongTable_old_check_rowids(rowids)
    .OnDiskLongTable_old_check_dirpath(dirpath)
    breakpoints <- .load_object(dirpath, "breakpoints")
    nrow <- .OnDiskLongTable_old_get_nrow_from_breakpoints(breakpoints)
    if (length(rowids) != nrow)
        stop(wmsg("length of 'rowids' is incompatible ",
                  "with the OnDiskLongTable_old object in ", dirpath))
    .save_object(dirpath, "rowids", rowids,
                 compress, compression_level)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

OnDiskLongTable_old <- function(dirpath=".")
{
    .OnDiskLongTable_old_check_dirpath(dirpath)
    ans_colnames <- .load_object(dirpath, "colnames")
    ans_breakpoints <- .load_object(dirpath, "breakpoints")
    ans_rowids_cache <- new.env(parent=emptyenv())
    ans <- new("OnDiskLongTable_old", dirpath=dirpath,
                                      colnames=ans_colnames,
                                      breakpoints=ans_breakpoints,
                                      .rowids_cache=ans_rowids_cache)
    .OnDiskLongTable_old_check_breakpoints(ans_breakpoints)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setGeneric("breakpoints", function(x) standardGeneric("breakpoints"))

setMethod("breakpoints", "OnDiskLongTable_old", function(x) x@breakpoints)

setGeneric("blocksizes", function(x) standardGeneric("blocksizes"))

setMethod("blocksizes", "OnDiskLongTable_old",
    function(x)
    {
        x_breakpoints <- breakpoints(x)
        ans <- S4Vectors:::diffWithInitialZero(x_breakpoints)
        names(ans) <- names(x_breakpoints)
        ans
    }
)

setMethod("dim", "OnDiskLongTable_old",
    function(x)
    {
        x_nrow <- .OnDiskLongTable_old_get_nrow_from_breakpoints(breakpoints(x))
        x_ncol <- length(colnames(x))
        c(x_nrow, x_ncol)
    }
)

setMethod("dimnames", "OnDiskLongTable_old",
    function(x) list(NULL, x@colnames)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method
###

setMethod("show", "OnDiskLongTable_old",
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
### Get batches from OnDiskLongTable_old
###

.normarg_batch_labels <- function(batch_labels, x_batch_labels)
{
    if (is.null(x_batch_labels))
        stop(wmsg("'x' has no batch labels (i.e. 'breakpoints(x)' has ",
                  "no names): cannot use getBatchesFromOnDiskLongTable_old() ",
                  "on it"))
    if (!is.character(batch_labels)
     || S4Vectors:::anyMissing(batch_labels)
     || anyDuplicated(batch_labels))
        stop(wmsg("'batch_labels' must be a character vector ",
                  "with no NAs and no duplicates"))
    batch_labels
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

### colidx: column index of length 1
.get_blocks_from_column <- function(x, colidx, blockidx)
{
    if (length(blockidx) == 0L)
        blockidx <- 0L
    tmp <- lapply(blockidx, function(i) .load_block(x, colidx, i))
    S4Vectors:::quick_unlist(tmp)
}

### batch_labels: character vector of batch labels.
### colidx: integer or character vector.
### Return a DataFrame (or data.frame).
getBatchesFromOnDiskLongTable_old <- function(x, batch_labels, colidx,
                                              with.batch_label=FALSE,
                                              with.rowids=FALSE,
                                              as.data.frame=FALSE)
{
    if (!is(x, "OnDiskLongTable_old"))
        stop(wmsg("'x' must be an OnDiskLongTable_old object"))
    x_blocksizes <- blocksizes(x)
    x_batch_labels <- names(x_blocksizes)
    batch_labels <- .normarg_batch_labels(batch_labels, x_batch_labels)
    colidx <- .normarg_colidx(colidx, x)
    if (!isTRUEorFALSE(with.batch_label))
        stop(wmsg("'with.batch_label' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(with.rowids))
        stop(wmsg("'with.rowids' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.data.frame))
        stop(wmsg("'as.data.frame' must be TRUE or FALSE"))
    hits <- findMatches(batch_labels, x_batch_labels)
    blockidx <- subjectHits(hits)
    ans_blocksizes <- x_blocksizes[blockidx]
    names(colidx) <- colnames(x)[colidx]
    ans_listData <-
        lapply(colidx, function(j) .get_blocks_from_column(x, j, blockidx))
    if (with.batch_label) {
        batch_label_levels <- unique(x_batch_labels)
        ans_batch_label <- factor(batch_labels[queryHits(hits)],
                                  levels=batch_label_levels)
        ans_batch_label <- rep.int(ans_batch_label, ans_blocksizes)
        ans_listData[["batch_label"]] <- ans_batch_label
    }
    ans_rowids <- NULL
    if (with.rowids) {
        x_rowids_env <- get_rowids_env_old(x)
        if (length(ls(x_rowids_env)) != 0L) {
            rowidx <- IRanges:::unlist_as_integer(
                          successiveIRanges(x_blocksizes)[blockidx]
                      )
            ans_rowids <- extract_rowids(x_rowids_env, rowidx)
        }
    }
    if (as.data.frame) {
        ans <- data.frame(ans_listData, row.names=ans_rowids,
                          stringsAsFactors=FALSE)
    } else {
        ## Unfortunately, DataFrame cannot store its row names as an integer
        ## vector.
        if (!is.null(ans_rowids))
            ans_rowids <- as.character(ans_rowids)
        ans <- new("DFrame", listData=ans_listData,
                             nrows=sum(ans_blocksizes),
                             rownames=ans_rowids)
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getRowsByIndexFromOnDiskLongTable_old()
###

.normarg_rowidx <- function(rowidx, x)
{
    if (!is.integer(rowidx))
        stop(wmsg("'rowidx' must be an integer vector"))
    if (S4Vectors:::anyMissingOrOutside(rowidx, 1L, nrow(x)))
        stop(wmsg("'rowidx' contains NAs or invalid row indices"))
    rowidx
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
    S4Vectors:::quick_unsplit(tmp, rowkeys[[1L]])
}

### rowidx: integer vector.
### colidx: integer or character vector.
### Return a DataFrame (or data.frame) with 1 row per row id in 'rowidx'.
### Note that we do NOT set the row names to 'rowidx' on the returned DataFrame
### because we want to support duplicates in 'rowidx'.
getRowsByIndexFromOnDiskLongTable_old <- function(x, rowidx, colidx,
                                                  with.batch_label=FALSE,
                                                  as.data.frame=FALSE)
{
    if (!is(x, "OnDiskLongTable_old"))
        stop(wmsg("'x' must be an OnDiskLongTable_old object"))
    rowidx <- .normarg_rowidx(rowidx, x)
    x_breakpoints <- breakpoints(x)
    rowkeys <- .rowidx2rowkeys(x_breakpoints, rowidx)
    colidx <- .normarg_colidx(colidx, x)
    if (!isTRUEorFALSE(with.batch_label))
        stop(wmsg("'with.batch_label' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(as.data.frame))
        stop(wmsg("'as.data.frame' must be TRUE or FALSE"))
    names(colidx) <- colnames(x)[colidx]
    ans_listData <-
        lapply(colidx, function(j) .get_values_from_column(x, j, rowkeys))
    if (with.batch_label) {
        batch_label_levels <- unique(names(x_breakpoints))
        ans_listData[["batch_label"]] <- factor(names(rowkeys[[1L]]),
                                                levels=batch_label_levels)
    }
    if (as.data.frame) {
        ans <- data.frame(ans_listData, stringsAsFactors=FALSE)
    } else {
        ans <- S4Vectors:::new_DataFrame(ans_listData, nrows=length(rowidx))
    }
    ans
}

