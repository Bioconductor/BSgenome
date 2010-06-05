setClass("BSParams",
         representation=representation(
           X="BSgenome",
           FUN="function",
           exclude = "character",
           simplify="logical",
           maskList="logical",
           motifList="character",
           userMask="RangesList",
           invertUserMask="logical"),
         prototype=prototype(
           exclude=character(0),
           simplify=FALSE,
           maskList=logical(0),
           motifList=character(0),
           invertUserMask=FALSE,
           userMask=RangesList()
           ),
         validity=function(object) {
             msg <- TRUE
             if (length(object@maskList) > 0) {
                 for(i in seq_len(length(object@maskList))) {
                     if(!(names(object@maskList)[i] %in% masknames(object@X))) {
                         msg <-
                           "the names of 'maskList' must be vector of names corresponding to the default BSgenome masks."
                         break
                     }
                 }
             }
             msg
         })

bsapply <-
function(BSParams, ...)
{
    if (!is(BSParams, "BSParams") ||
        !isTRUE(validObject(BSParams, test=TRUE)))
        stop("'X' must be a valid BSgenome object")

    processSeqname <- function(seqname, ...)
    {
        seq <- BSParams@X[[seqname]]

        if (length(BSParams@maskList) > 0) {
            for(i in seq_len(length(BSParams@maskList))){
                if(names(BSParams@maskList)[i] %in% masknames(BSParams@X)) {
                    active(masks(seq))[names(BSParams@maskList)[i]] <-
                      BSParams@maskList[[i]]
                }
            }
        }

        if (length(BSParams@motifList) > 0){
            for(i in seq_len(length(BSParams@motifList))) {
                seq <- maskMotif(seq, BSParams@motifList[i])
            }
        }

        seqMask <- userMask[[seqname]]
        if (length(seqMask) > 0) {
            seqMask <- Mask(length(seq), start(seqMask), end(seqMask))
            if (BSParams@invertUserMask)
                seqMask <- gaps(seqMask)
            if (!is.null(masks(seq)))
                seqMask <- append(masks(seq), seqMask)
            masks(seq) <- seqMask
        }

        BSParams@FUN(seq, ...)
    }

    if (length(seqnames(BSParams@X)) > 0)
        seqnames <- seqnames(BSParams@X)
    else
        seqnames <- mseqnames(BSParams@X)

    exclude <- unname(BSParams@exclude[nzchar(BSParams@exclude)])
    pariahIndex <- unique(unlist(lapply(exclude, grep, seqnames)))
    if (length(pariahIndex) > 0)
        seqnames <- seqnames[- pariahIndex]

    userMask <- reduce(BSParams@userMask)

    if (BSParams@simplify) {
        sapply(seqnames, processSeqname, ...)
    } else {
        GenomeData(lapply(structure(seqnames, names = seqnames), processSeqname, ...),
                   providerVersion = providerVersion(BSParams@X),
                   organism = organism(BSParams@X),
                   provider = provider(BSParams@X))
    }
}
