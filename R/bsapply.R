bsapply <- function(X, FUN, exclude = "rand", ...){

    ##get the csomes:
    csomes = seqnames(X)

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
    
    ##Some dirty tricks so that I can get the name of this thing:
    pkg = X@package
    pkgName = ls(paste("package:",pkg,sep =""))

    ##get a list of existing csome item the user has in memory 'right now'
    loadedCSomes = ls(X@.datacache_env)
    
    ##Some stuff has to be done for each chromosome
    processElement = function(element, ...){
        str = paste(pkgName,"[[\"",element,"\"]]", sep ="")
        chr = eval(parse(text=str))
        result = FUN(chr, ...)
        ## Here is where I need to clean up after each of these
        ## just check if the thing is in the list, and if not = jettison
        if(element %in% loadedCSomes){
            print(paste("Refusing to unload Chromosome:", element))
        }else{
            print(paste("Unloading Chromosome:",element ))
            unload(X, element)
        }
        return(result)
    }

    ##Then apply the above function to each
    list = lapply(csomes, processElement, ...)
    names(list) = csomes
    return(list)

}
