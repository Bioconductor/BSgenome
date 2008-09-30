bsapply <- function(X, FUN, exclude = "", simplify = FALSE, maskList = c(), ...){

    ##Some argument checking.
    if(!is(X, "BSgenome")) stop("'X' must be a BSgenome object")
    if(!is.function(FUN)) stop("'FUN' must be a function")
    if(!is.character(exclude)) stop("'exclude' must be a character vector")

    #Argument checking for maskList:
    if(length(maskList)>0){ ##if there are masks
        for(i in seq_len(length(maskList))){ #loop thru masks and deal with each in turn
            if( !( maskList[i] %in% masknames(X)) ){stop("'maskList' must be vector of names corresponding to the default BSgenome masks.")}
            ##     if( !( maskList[i] %in% masknames(X) || is.function(maskList[i]) || is(maskList[i], "IRanges") ) ){stop("'masks' must be an integer corresponding to the default BSgenome masks that are desired, an IRanges object that be used to generate a mask or a function that returns an IRanges object so that the masks can be generated.")}
        }
    }

            

    
    
    ##get the seqnames:
    seqnames <- seqnames(X)
    seqLength = length(seqnames)

    ##Restrict the seqnames based on our exclusion factor:
    pariahIndex = numeric()
    for(i in seq_len(length(exclude))){
        ind = grep(exclude[i],seqnames)
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

        seq = X[[seqname]]
        
        if(length(maskList)>0){ ##if there are masks
            for(i in seq_len(length(maskList))){ #loop thru masks and deal with each in turn
                if(maskList[i] %in% masknames(X)){#IF its one of the names, then change its active state to be different from the default...
                    if(active(masks(seq))[maskList[i]] == FALSE){
                        active(masks(seq))[maskList[i]] <- TRUE
                    }else{active(masks(seq))[maskList[i]] <- FALSE}
                }
##                 if(is.function(maskList[i])){#IF its a function, then lets call it to generate a GENOMIC RANGES object
##                     ##
##                 }
            }
        }


        

        ## This is where we finally get run the function that was passed in
        result <- FUN(seq, ...)
        return(result)
    }

    
    ##Then apply the above function to each
    if(simplify == FALSE){
        ans <- lapply(seqnames, processSeqname, ...)
        names(ans) <- seqnames
        return(ans)
    }else{
        ans <- sapply(seqnames, processSeqname, ...)
        return(ans)
    }
    ans
}

