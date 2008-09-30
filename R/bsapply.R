bsapply <- function(X, FUN, exclude = "", simplify = FALSE, ...)
{
    ##Some argument checking.
    if(!is(X, "BSgenome")) stop("'X' must be a BSgenome object")
    if(!is.function(FUN)) stop("'FUN' must be a function")
    if(!is.character(exclude)) stop("'exclude' must be a character vector")

    ##get the csomes:
    csomes <- seqnames(X)
    csomeLength <- length(csomes)

    ##Restrict the csomes based on our exclusion factor:
    pariahIndex <- sapply(exclude, grep, csomes)

    ##Added precaution in case some indices are found more than once...
    pariahIndex <- unique(pariahIndex)
    ##IF the grepping has found something then we need to change our csomes
    ##But let's not be silly and allow people to exclude EVERYTHING.
    if(length(pariahIndex) > 0 && length(pariahIndex) != csomeLength){
        csomes <- csomes[-pariahIndex]
    }
    
    ##Some stuff has to be done for each chromosome
    processSeqname <- function(seqname, ...) FUN(X[[seqname]], ...)

    ##Then apply the above function to each
    if(simplify == FALSE){
        ans <- lapply(csomes, processSeqname, ...)
        names(ans) <- csomes
    }else{
        ans <- sapply(csomes, processSeqname, ...)
    }
    ans
}

