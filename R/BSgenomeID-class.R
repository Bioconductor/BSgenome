## Formal representation of a BSgenome ID (eg Celegans.UCSC.ce2)

setClass("BSgenomeID", contains = "character",
         validity = function(object) {
           if (length(object) != 1)
             "BSgenomeID must be a single string"
           else if (length(.BSgenomeID.tokens(object)) != 3)
             "BSgenomeID must have 3 tokens separated by '.'"
           else NULL
         })

.BSgenomeID.tokens <- function(x) strsplit(x, "\\.")[[1]]

setMethod("organism", "BSgenomeID", function(x) .BSgenomeID.tokens(x)[1])
setMethod("provider", "BSgenomeID", function(x) .BSgenomeID.tokens(x)[2])
setMethod("providerVersion", "BSgenomeID", function(x) .BSgenomeID.tokens(x)[3])

setMethod("genome", "AnnotatedList",
          function(x) {
            ann <- annotation(x)
            if (!is.null(ann))
              ann <- new("BSgenomeID", ann)
            ann
          })
