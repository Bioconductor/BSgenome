### The strand generic and the strand() method, for lack of a better home.

setGeneric("strand", function(x, ...) standardGeneric("strand"))

setMethod("strand", "missing", function(x) {
  factor(levels=c("-","+","*"))
})

setMethod(strand, "character", function(x, ...) {
  lvls <- c(NA, levels(strand()))
  if (!all(x %in% lvls | is.na(x)))
    stop("strand values must be in '", paste(lvls, collapse="' '"), "'")
  factor(x, levels=lvls)
})
