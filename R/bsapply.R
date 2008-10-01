setClass("BSParams",
         representation=representation(
           X="BSgenome",
           FUN="function",
           exclude = "character",
           simplify="logical",
           maskList="ANY"),
         prototype=prototype(
           exclude="",
           simplify=FALSE,
           maskList=c()
           ))


bsapply <- function(BSParams, ...){##X, FUN, exclude = "", simplify = FALSE, maskList = c(), ...){

    ##Some argument checking in case someone changes the value inside the params object.
    if(!is(BSParams@X, "BSgenome")) stop("'X' must be a BSgenome object")
    if(!is.function(BSParams@FUN)) stop("'FUN' must be a function")
    if(!is.character(BSParams@exclude)) stop("'exclude' must be a character vector")

    ##Argument checking for maskList:
    if(length(BSParams@maskList)>0){ ##if there are masks
        for(i in seq_len(length(BSParams@maskList))){ #loop thru masks and deal with each in turn
            if( !( BSParams@maskList[i] %in% masknames(BSParams@X)) ){stop("'BSParams@maskList' must be vector of names corresponding to the default BSgenome masks.")}
            ##     if( !( BSParams@maskList[i] %in% masknames(BSParams@X) || is.function(BSParams@maskList[i]) || is(BSParams@maskList[i], "IRanges") ) ){stop("'masks' must be an integer corresponding to the default BSgenome masks that are desired, an IRanges object that be used to generate a mask or a function that returns an IRanges object so that the masks can be generated.")}
        }
    }
    
    
    ##get the seqnames:
    seqnames <- seqnames(BSParams@X)
    seqLength = length(seqnames)

    ##Restrict the seqnames based on our exclusion factor:
    pariahIndex = numeric()
    for(i in seq_len(length(BSParams@exclude))){
        ind = grep(BSParams@exclude[i],seqnames)
        ##Catenate the indices together as we go
        pariahIndex = c(pariahIndex, ind)
    }

    ##Added precaution in case some indices are found more than once...
    pariahIndex <- unique(pariahIndex)
    ##IF the grepping has found something then we need to change our seqnames
    ##But let's not be silly and allow people to exclude EVERYTHING.
    if(length(pariahIndex) > 0 && length(pariahIndex) != seqLength){
        seqnames <- seqnames[-pariahIndex]
    }
    

    ##Some stuff has to be done for each chromosome
    processSeqname <- function(seqname, ...){

        seq = BSParams@X[[seqname]]
        
        if(length(BSParams@maskList)>0){ ##if there are masks
            for(i in seq_len(length(BSParams@maskList))){ #loop thru masks and deal with each in turn
                if(BSParams@maskList[i] %in% masknames(BSParams@X)){#IF its one of the names, then change its active state to be different from the default...
                    if(active(masks(seq))[BSParams@maskList[i]] == FALSE){
                        active(masks(seq))[BSParams@maskList[i]] <- TRUE
                    }else{active(masks(seq))[BSParams@maskList[i]] <- FALSE}
                }
##                 if(is.function(BSParams@maskList[i])){#IF its a function, then lets call it to generate a GENOMIC RANGES object
##                     ##
##                 }
            }
        }


        

        ## This is where we finally get run the function that was passed in
        result <- BSParams@FUN(seq, ...)
        return(result)
    }

    
    ##Then apply the above function to each
    if(BSParams@simplify == FALSE){
        ans <- lapply(seqnames, processSeqname, ...)
        names(ans) <- seqnames
        return(ans)
    }else{
        ans <- sapply(seqnames, processSeqname, ...)
        return(ans)
    }
    ans
}

