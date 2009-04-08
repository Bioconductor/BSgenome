setClass("BSParams",
         representation=representation(
           X="BSgenome",
           FUN="function",
           exclude = "character",
           simplify="logical",
           maskList="logical",
           motifList="character"),
         prototype=prototype(
           exclude="",
           simplify=FALSE,
           maskList=as.logical(vector()),
           motifList=as.character(vector())
           ))


bsapply <- function(BSParams, ...){##X, FUN, exclude = "", simplify = FALSE, maskList = c(), ...){

    ##Some argument checking in case someone changes the value inside the params object.
    if(!is(BSParams@X, "BSgenome")) stop("'X' must be a BSgenome object")
    if(!is.function(BSParams@FUN)) stop("'FUN' must be a function")
    if(!is.character(BSParams@exclude)) stop("'exclude' must be a character vector")

    ##Argument checking for maskList:
    if(length(BSParams@maskList)>0){ ##if there are masks
        for(i in seq_len(length(BSParams@maskList))){ #loop thru masks and deal with each in turn
            if( !( names(BSParams@maskList)[i] %in% masknames(BSParams@X)) ){stop("the names of 'BSParams@maskList' must be vector of names corresponding to the default BSgenome masks.")}
        }
    }
    
    ##get the seqnames:
    seqnames <- seqnames(BSParams@X)
    seqLength <- length(seqnames)

    ##Restrict the seqnames based on our exclusion factor:
    pariahIndex <- numeric()
    for(i in seq_len(length(BSParams@exclude))){
        ind <- grep(BSParams@exclude[i],seqnames)
        ##Catenate the indices together as we go
        pariahIndex <- c(pariahIndex, ind)
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

        seq <- BSParams@X[[seqname]]
        
        if(length(BSParams@maskList)>0){ ##if there are masks
            for(i in seq_len(length(BSParams@maskList))){ ##loop thru masks and deal with each in turn
                if(names(BSParams@maskList)[i] %in% masknames(BSParams@X)){#IF its one of the named masks, then set it accordingly
                        active(masks(seq))[names(BSParams@maskList)[i]] <- BSParams@maskList[[i]]
                }
            }
        }

        if(length(BSParams@motifList)>0){ ##if there are motifs
            for(i in seq_len(length(BSParams@motifList))){ ##loop thru motifs 
                seq <- maskMotif(seq, BSParams@motifList[i]) ##mask off each motif
            }
        }


        ## This is where we finally get to run the MAIN function (FUN) that was passed in
        result <- BSParams@FUN(seq, ...)
        return(result)
    }

    
    ##Then apply the above function to each
    if(BSParams@simplify == FALSE){
        ans <-
          GenomeData(lapply(structure(seqnames, names = seqnames), processSeqname, ...),
                     providerVersion = providerVersion(BSParams@X),
                     organism = organism(BSParams@X),
                     provider = provider(BSParams@X))
        return(ans)
    }else{
        ans <- sapply(seqnames, processSeqname, ...)
        return(ans)
    }
    ans
}

