bsapply <- function(X, FUN, exclude = "rand", ...){

    #HP: Please check the arguments. E.g. what would happen if the 'X'
    #you get from the user was not a BSgenome object?

    #HP: Please use <- instead of =

    ##get the csomes:
    csomes = seqnames(X)

    #HP: The single line below is equivalent to the next 11 lines
    #csomes <- csomes[-grep(exclude, csomes)]

    ##Restrict the csomes based on our exclusion factor:
    pariahIndex = grep(exclude, csomes)
    keepIndex = logical()
    for(i in 1:length(csomes)){
        if(i %in% pariahIndex){
            keepIndex = c(keepIndex, FALSE)
        }else{
            keepIndex = c(keepIndex, TRUE)
        }
    }    
    csomes = csomes[keepIndex]
    
    #HP: You don't need this dirty trick. Just use X[[element]]
    #in your processElement() function below.

    ##Some dirty tricks so that I can get the name of this thing:
    pkg = X@package
    pkgName = ls(paste("package:",pkg,sep =""))

    ##get a list of existing csome item the user has in memory 'right now'
    loadedCSomes = ls(X@.datacache_env)
    
    #HP: Please use 'seqname' instead of 'element'
    ##Some stuff has to be done for each chromosome
    processElement = function(element, ...){
        str = paste(pkgName,"[[\"",element,"\"]]", sep ="")
        chr = eval(parse(text=str))
        #HP: or just do chr <- X[[element]]
        result = FUN(chr, ...)
        ## Here is where I need to clean up after each of these
        ## just check if the thing is in the list, and if not = jettison
        if(element %in% loadedCSomes){
            #HP: Or use cat()
            #cat("Refusing to unload Chromosome:", element, "\n")
            print(paste("Refusing to unload Chromosome:", element))
        }else{
            print(paste("Unloading Chromosome:",element ))
            unload(X, element)
        }
        return(result)
    }

    #HP: Before you start the main loop, call on.exit(). See ?on.exit
    #for the details. This guarantees that some code will be executed
    #before the function actually exit, even in case of user interrupt
    #or other kind of abnormal termination. This code must unload any
    #sequence that is currently loaded but was not loaded when the function
    #started.

    ##Then apply the above function to each
    list = lapply(csomes, processElement, ...)
    names(list) = csomes
    return(list)

}
