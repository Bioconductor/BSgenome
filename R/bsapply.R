bsapply <- function(X, FUN, exclude = "", ...){

    ##Some argument checking.
    if(!is(X, "BSgenome")){stop("'X' must be a BSgenome object")}
    if(!is.function(FUN)){stop("'FUN' must be a function that takes a DNAString object as its 1st argument.")}
    if(!is.character(exclude) || length(exclude)>1){stop("'exclude' must be a single string.")}

    ##get the csomes:
    csomes <- seqnames(X)

    ##Restrict the csomes based on our exclusion factor:
    if(length(grep(exclude, csomes)) > 0 && length(grep(exclude, csomes)) != length(csomes)){
        csomes <- csomes[-grep(exclude, csomes)]
    }
    
    ##get a list of existing csome item the user has in memory 'right now'
    loadedCSomes <- ls(X@.datacache_env)
    
    ##Some stuff has to be done for each chromosome
    processSeqname <- function(seqname, ...){
        result <- FUN(X[[seqname]], ...)
        ## Here is where I need to clean up after each of these
        ## just check if the thing is in the list, and if not,
        ## then jettison
        if(seqname %in% loadedCSomes){
            ##output will be dropped as we gain confidence,
            ##or I can wrap it up so that it can be toggled off
            cat("Keeping Chromosome:", seqname, "\n")
        }else{
            cat("Unloading Chromosome:",seqname, "\n")
            unload(X, seqname)
        }
        return(result)
    }

    ##on.exit is called as an insurance policy. We just want to clear all
    ##the memory "just in case" something got interrupted.

    checkSeqs <- function(){
            ##cat("Checking for any large loaded sequences", "\n")
            for(i in seq_len(length(csomes))){
                ##1st we check to make sure that it was NOT loaded at the beginning
                if(!(csomes[i] %in% loadedCSomes)){
                    ##Then we should also check if it is loaded "right now"
                    if(csomes[i] %in% ls(X@.datacache_env)){
                        cat("Unloading unwanted Chromosome:",csomes[i], "\n")
                        unload(X, csomes[i])
                    }
                }}
    }    
    on.exit(checkSeqs())
    
    ##Then apply the above function to each
    list <- lapply(csomes, processSeqname, ...)
    names(list) <- csomes
    return(list)

}
