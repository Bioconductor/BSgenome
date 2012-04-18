## vmatchPattern/vcountPattern for BSgenome

setMethod("vmatchPattern", "BSgenome",
    function(pattern, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", exclude="", maskList=logical(0),
             userMask=RangesList(), invertUserMask=FALSE,
             asRangedData=FALSE)
    {
        matchFUN <- function(posPattern, negPattern, chr,
                             seqlengths,
                             max.mismatch = max.mismatch,
                             min.mismatch = min.mismatch,
                             with.indels = with.indels,
                             fixed = fixed,
                             algorithm = algorithm) {
            posMatches <-
              matchPattern(pattern = posPattern, subject = chr,
                           max.mismatch = max.mismatch,
                           min.mismatch = min.mismatch,
                           with.indels = with.indels,
                           fixed = fixed,
                           algorithm = algorithm)
            negMatches <-
              matchPattern(pattern = negPattern, subject = chr,
                           max.mismatch = max.mismatch,
                           min.mismatch = min.mismatch,
                           with.indels = with.indels,
                           fixed = fixed,
                           algorithm = algorithm)
            COUNTER <<- COUNTER + 1L
            seqnames <- names(seqlengths)
            GRanges(seqnames = 
                    Rle(factor(seqnames[COUNTER], levels = seqnames),
                        length(posMatches) + length(negMatches)),
                    ranges =
                    c(as(posMatches, "IRanges"), as(negMatches, "IRanges")),
                    strand =
                    Rle(strand(c("+", "-")),
                        c(length(posMatches), length(negMatches))),
                    seqlengths = seqlengths)
        }

        if (!is(pattern, "DNAString"))
            pattern <- DNAString(pattern)
        algorithm <- Biostrings:::normargAlgorithm(algorithm)
        if (Biostrings:::isCharacterAlgo(algorithm)) 
            stop("'subject' must be a single (non-empty) string ", 
                 "for this algorithm")
        pattern <- Biostrings:::normargPattern(pattern, DNAStringSet())
        max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
        min.mismatch <-
          Biostrings:::normargMinMismatch(min.mismatch, max.mismatch)
        with.indels <- Biostrings:::normargWithIndels(with.indels)
        fixed <- Biostrings:::normargFixed(fixed, DNAStringSet())
        if (!missing(asRangedData)) {
            msg <- c("the 'asRangedData' argument is deprecated ",
                     "and will be removed soon")
            .Deprecated(msg=msg)
        }
        if (asRangedData) {
            msg <- "using 'asRangedData=TRUE' is defunct"
            .Defunct(msg=msg)
        }

        posPattern <- pattern
        negPattern <- reverseComplement(posPattern)
        bsParams <-
          new("BSParams", X = subject, FUN = matchFUN, exclude = exclude,
              simplify = FALSE, maskList = logical(0), userMask = userMask,
              invertUserMask = invertUserMask)
        COUNTER <- 0L
        seqlengths <- seqlengths(subject)
        matches <-
          bsapply(bsParams, posPattern = posPattern, negPattern = negPattern,
                  seqlengths = seqlengths,
                  max.mismatch = max.mismatch, min.mismatch = min.mismatch,
                  with.indels = with.indels, fixed = fixed,
                  algorithm = algorithm)
        nms <- factor(names(matches), levels = names(seqlengths))
        nms <- nms[unlist(lapply(matches, length), use.names=FALSE) > 0]
        matches <- do.call(c, unname(as.list(matches)))
        runValue(seqnames(matches)) <- nms
        matches
    }
)

setMethod("vcountPattern", "BSgenome",
    function(pattern, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", exclude="", maskList=logical(0),
             userMask=RangesList(), invertUserMask=FALSE)
    {
        countFUN <- function(posPattern, negPattern, chr,
                             max.mismatch = max.mismatch,
                             min.mismatch = min.mismatch,
                             with.indels = with.indels,
                             fixed = fixed,
                             algorithm = algorithm) {
            data.frame(strand = strand(c("+", "-")),
                       count =
                       c(countPattern(pattern = posPattern, subject = chr,
                                      max.mismatch = max.mismatch,
                                      min.mismatch = min.mismatch,
                                      with.indels = with.indels,
                                      fixed = fixed,
                                      algorithm = algorithm),
                         countPattern(pattern = negPattern, subject = chr,
                                      max.mismatch = max.mismatch,
                                      min.mismatch = min.mismatch,
                                      with.indels = with.indels,
                                      fixed = fixed,
                                      algorithm = algorithm)))
        }

        if (!is(pattern, "DNAString"))
            pattern <- DNAString(pattern)
        algorithm <- Biostrings:::normargAlgorithm(algorithm)
        if (Biostrings:::isCharacterAlgo(algorithm)) 
            stop("'subject' must be a single (non-empty) string ", 
                 "for this algorithm")
        pattern <- Biostrings:::normargPattern(pattern, DNAStringSet())
        max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
        min.mismatch <- Biostrings:::normargMinMismatch(min.mismatch, max.mismatch)
        with.indels <- Biostrings:::normargWithIndels(with.indels)
        fixed <- Biostrings:::normargFixed(fixed, DNAStringSet())

        posPattern <- pattern
        negPattern <- reverseComplement(posPattern)
        bsParams <-
          new("BSParams", X = subject, FUN = countFUN, exclude = exclude,
              simplify = FALSE, maskList = logical(0), userMask = userMask,
              invertUserMask = invertUserMask)
        counts <-
          bsapply(bsParams, posPattern = posPattern, negPattern = negPattern,
                  max.mismatch = max.mismatch, min.mismatch = min.mismatch,
                  with.indels = with.indels, fixed = fixed,
                  algorithm = algorithm)
        cbind(data.frame(seqname =
                         rep(factor(names(counts), levels = names(counts)),
                             each = 2)),
              do.call(rbind, unname(as.list(counts))))
    }
)


## vmatchPDict/vcountPDict for BSgenome
setMethod("vmatchPDict", "BSgenome",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, fixed=TRUE,
             algorithm="auto", verbose=FALSE, exclude="",
             maskList=logical(0), asRangedData=FALSE)
    {
        matchFUN <- function(posPDict, negPDict, chr,
                             seqlengths,
                             max.mismatch = max.mismatch,
                             min.mismatch = min.mismatch,
                             fixed = fixed,
                             algorithm = algorithm,
                             verbose = verbose) {
            posMatches <-
              matchPDict(pdict = posPDict, subject = chr,
                         max.mismatch = max.mismatch,
                         min.mismatch = min.mismatch,
                         fixed = fixed,
                         algorithm = algorithm,
                         verbose = verbose)
            posCounts <- elementLengths(posMatches)
            negMatches <-
              matchPDict(pdict = negPDict, subject = chr,
                         max.mismatch = max.mismatch,
                         min.mismatch = min.mismatch,
                         fixed = fixed,
                         algorithm = algorithm,
                         verbose = verbose)
            negCounts <- elementLengths(negMatches)
            COUNTER <<- COUNTER + 1L
            seqnames <- names(seqlengths)
            GRanges(seqnames = 
                    Rle(factor(seqnames[COUNTER], levels = seqnames),
                        sum(posCounts) + sum(negCounts)),
                    ranges = c(unlist(posMatches), unlist(negMatches)),
                    strand =
                    Rle(strand(rep(c("+", "-"))),
                        c(sum(posCounts), sum(negCounts))),
                    index =
                    c(Rle(seq_len(length(posMatches)), posCounts),
                      Rle(seq_len(length(negMatches)), negCounts)),
                    seqlengths = seqlengths)
        }

        if (is(pdict, "PDict"))
            stop("'pdict' must be a DNAStringSet object ")
        if (!is(pdict, "DNAStringSet"))
            pdict <- DNAStringSet(pdict)
        max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
        min.mismatch <-
          Biostrings:::normargMinMismatch(min.mismatch, max.mismatch)
        fixed <- Biostrings:::normargFixed(fixed, DNAStringSet())
        algorithm <- Biostrings:::normargAlgorithm(algorithm)
        if (Biostrings:::isCharacterAlgo(algorithm)) 
            stop("'subject' must be a single (non-empty) string ", 
                 "for this algorithm")
        if (!isTRUEorFALSE(verbose))
            stop("'verbose' must be TRUE or FALSE")
        if (!missing(asRangedData)) {
            msg <- c("the 'asRangedData' argument is deprecated ",
                     "and will be removed soon")
            .Deprecated(msg=msg)
        }
        if (asRangedData) {
            msg <- "using 'asRangedData=TRUE' is defunct"
            .Defunct(msg=msg)
        }

        posPDict <- pdict
        negPDict <- reverseComplement(posPDict)
        bsParams <-
          new("BSParams", X = subject, FUN = matchFUN, exclude = exclude,
              simplify = FALSE, maskList = logical(0))
        COUNTER <- 0L
        seqlengths <- seqlengths(subject)
        matches <-
          bsapply(bsParams,
                  posPDict = PDict(posPDict, max.mismatch = max.mismatch),
                  negPDict = PDict(negPDict, max.mismatch = max.mismatch),
                  seqlengths = seqlengths,
                  max.mismatch = max.mismatch, min.mismatch = min.mismatch,
                  fixed = fixed, algorithm = algorithm, verbose = verbose)
        nms <- factor(names(matches), levels = names(seqlengths))
        nms <- nms[unlist(lapply(matches, length), use.names=FALSE) > 0]
        matches <- do.call(c, unname(as.list(matches)))
        runValue(seqnames(matches)) <- nms
        matches
    }
)


setMethod("vcountPDict", "BSgenome",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, fixed=TRUE,
             algorithm="auto", collapse=FALSE, weight=1L, verbose=FALSE,
             exclude="", maskList=logical(0))
    {
        countFUN <- function(posPDict, negPDict, chr,
                             max.mismatch = max.mismatch,
                             min.mismatch = min.mismatch,
                             fixed = fixed,
                             algorithm = algorithm,
                             verbose = verbose) {
            posCounts <-
              countPDict(pdict = posPDict, subject = chr,
                         max.mismatch = max.mismatch,
                         min.mismatch = min.mismatch,
                         fixed = fixed,
                         algorithm = algorithm,
                         verbose = verbose)
            negCounts <-
              countPDict(pdict = negPDict, subject = chr,
                         max.mismatch = max.mismatch,
                         min.mismatch = min.mismatch,
                         fixed = fixed,
                         algorithm = algorithm,
                         verbose = verbose)
            DataFrame(strand =
                      Rle(strand(c("+", "-")),
                          c(length(posCounts), length(negCounts))),
                      index =
                      c(seq_len(length(posCounts)),
                        seq_len(length(negCounts))),
                      count = c(Rle(posCounts), Rle(negCounts)))
        }

        if (is(pdict, "PDict"))
            stop("'pdict' must be a DNAStringSet object ")
        if (!is(pdict, "DNAStringSet"))
            pdict <- DNAStringSet(pdict)
        max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
        min.mismatch <-
          Biostrings:::normargMinMismatch(min.mismatch, max.mismatch)
        fixed <- Biostrings:::normargFixed(fixed, DNAStringSet())
        algorithm <- Biostrings:::normargAlgorithm(algorithm)
        if (Biostrings:::isCharacterAlgo(algorithm)) 
            stop("'subject' must be a single (non-empty) string ", 
                 "for this algorithm")
        if (!identical(collapse, FALSE))
            stop("'collapse' is not supported in BSgenome method")
        if (!identical(weight, 1L))
            stop("'weight' is not supported in BSgenome method")
        if (!isTRUEorFALSE(verbose))
            stop("'verbose' must be TRUE or FALSE")

        posPDict <- pdict
        negPDict <- reverseComplement(posPDict)
        bsParams <-
          new("BSParams", X = subject, FUN = countFUN, exclude = exclude,
              simplify = FALSE, maskList = logical(0))
        counts <-
          bsapply(bsParams,
                  posPDict = PDict(posPDict, max.mismatch = max.mismatch),
                  negPDict = PDict(negPDict, max.mismatch = max.mismatch),
                  max.mismatch = max.mismatch, min.mismatch = min.mismatch,
                  fixed = fixed, algorithm = algorithm, verbose = verbose)
        DataFrame(DataFrame(seqname =
                            Rle(factor(names(counts), levels = names(counts)),
                                unlist(lapply(counts, nrow),
                                       use.names = FALSE))),
                  do.call(rbind, unname(as.list(counts))))
    }
)

## matchPWM/countPWM for BSgenome
setMethod("matchPWM", "BSgenome",
    function(pwm, subject, min.score="80%", exclude="", maskList=logical(0),
             asRangedData=FALSE)
    {
        matchFUN <- function(posPWM, negPWM, chr, seqlengths, min.score) {
            posMatches <-
              matchPWM(pwm = posPWM, subject = chr, min.score = min.score)
            posScores <-
              PWMscoreStartingAt(pwm = posPWM, subject = chr,
                                 starting.at = start(posMatches))
            negMatches <-
              matchPWM(pwm = negPWM, subject = chr, min.score = min.score)
            negScores <-
              PWMscoreStartingAt(pwm = negPWM, subject = chr,
                                 starting.at = start(negMatches))
            COUNTER <<- COUNTER + 1L
            seqnames <- names(seqlengths)
            GRanges(seqnames = 
                    Rle(factor(seqnames[COUNTER], levels = seqnames),
                        length(posMatches) + length(negMatches)),
                    ranges =
                    c(as(posMatches, "IRanges"), as(negMatches, "IRanges")),
                    strand =
                    Rle(strand(c("+", "-")),
                        c(length(posMatches), length(negMatches))),
                    score = c(posScores, negScores),
                    string =
                    c(as.character(posMatches),
                      as.character(reverseComplement(DNAStringSet(negMatches)))),
                    seqlengths = seqlengths)
        }

        ## checking 'pwm'
        pwm <- Biostrings:::.normargPwm(pwm)
        ## checking 'min.score'
        min.score <- Biostrings:::.normargMinScore(min.score, pwm)
        if (!missing(asRangedData)) {
            msg <- c("the 'asRangedData' argument is deprecated ",
                     "and will be removed soon")
            .Deprecated(msg=msg)
        }
        if (asRangedData) {
            msg <- "using 'asRangedData=TRUE' is defunct"
            .Defunct(msg=msg)
        }

        posPWM <- pwm[DNA_BASES, , drop = FALSE]
        negPWM <- reverseComplement(posPWM)
        bsParams <-
          new("BSParams", X = subject, FUN = matchFUN, exclude = exclude,
              simplify = FALSE, maskList = logical(0))
        COUNTER <- 0L
        seqlengths <- seqlengths(subject)
        matches <-
          bsapply(bsParams, posPWM = posPWM, negPWM = negPWM,
                  seqlengths = seqlengths, min.score = min.score)
        nms <- factor(names(matches), levels = names(seqlengths))
        nms <- nms[unlist(lapply(matches, length), use.names=FALSE) > 0]
        matches <- do.call(c, unname(as.list(matches)))
        runValue(seqnames(matches)) <- nms
        ## create compact DNAStringSet
        string <- factor(elementMetadata(matches)[["string"]])
        elementMetadata(matches)[["string"]] <-
          DNAStringSet(levels(string))[as.integer(string)]
        matches
    }
)

setMethod("countPWM", "BSgenome",
    function(pwm, subject, min.score="80%", exclude="", maskList=logical(0))
    {
        countFUN <- function(posPWM, negPWM, chr, min.score) {
            data.frame(strand = strand(c("+", "-")),
                       count =
                       c(countPWM(pwm = posPWM, subject = chr,
                                  min.score = min.score),
                         countPWM(pwm = negPWM, subject = chr,
                                  min.score = min.score)))
        }

        ## checking 'pwm'
        pwm <- Biostrings:::.normargPwm(pwm)
        ## checking 'min.score'
        min.score <- Biostrings:::.normargMinScore(min.score, pwm)

        posPWM <- pwm[DNA_BASES, , drop = FALSE]
        negPWM <- reverseComplement(posPWM)
        bsParams <-
          new("BSParams", X = subject, FUN = countFUN, exclude = exclude,
              simplify = FALSE, maskList = logical(0))
        counts <-
          bsapply(bsParams, posPWM = posPWM, negPWM = negPWM,
                  min.score = min.score)
        cbind(data.frame(seqname =
                         rep(factor(names(counts), levels = names(counts)),
                             each = 2)),
              do.call(rbind, unname(as.list(counts))))
    }
)

setGeneric("vmatchLRPatterns", signature = "subject",
           function(Lpattern, Rpattern, max.gaplength, subject,
                    max.Lmismatch=0, max.Rmismatch=0,
                    with.Lindels=FALSE, with.Rindels=FALSE,
                    Lfixed=TRUE, Rfixed=TRUE, ...)
           standardGeneric("vmatchLRPatterns"))

setMethod("vmatchLRPatterns", "BSgenome",
          function(Lpattern, Rpattern, max.gaplength, subject,
                   max.Lmismatch=0, max.Rmismatch=0,
                   with.Lindels=FALSE, with.Rindels=FALSE,
                   Lfixed=TRUE, Rfixed=TRUE,
                   exclude = "", userMask = RangesList(),
                   invertUserMask = FALSE)
          {
            matchFUN <- function(posLPattern, negLPattern,
                                 posRPattern, negRPattern, chr,
                                 max.gaplength, subject,
                                 max.Lmismatch, max.Rmismatch,
                                 with.Lindels, with.Rindels,
                                 Lfixed, Rfixed) {
              posMatches <-
                matchLRPatterns(posLPattern, posRPattern,
                               max.gaplength, chr,
                               max.Lmismatch, max.Rmismatch,
                               with.Lindels, with.Rindels,
                               Lfixed, Rfixed)
              negMatches <- IRanges()
              if (append(posLPattern, posRPattern) !=
                  append(negRPattern, negLPattern))
               negMatches <-
                matchLRPatterns(negLPattern, negRPattern,
                               max.gaplength, chr,
                               max.Lmismatch, max.Rmismatch,
                               with.Lindels, with.Rindels,
                               Lfixed, Rfixed)
              strings <- c(as.character(posMatches), as.character(negMatches))
              strand <- rep(c("+", "-"),
                            c(length(posMatches), length(negMatches)))
              ord <- order(strand, strings)
              rngs <- c(as(posMatches, "IRanges"), as(negMatches, "IRanges"))
              RangedData(rngs[ord], string = Rle(strings[ord]),
                         strand = Rle(strand[ord]),
                         space = "1")
            }

            if (!is(Lpattern, "DNAString"))
              pattern <- DNAString(Lpattern)
            if (!is(Rpattern, "DNAString"))
              pattern <- DNAString(Rpattern)
            
            Lpattern <- Biostrings:::normargPattern(Lpattern, DNAStringSet())
            Rpattern <- Biostrings:::normargPattern(Rpattern, DNAStringSet())
            max.gaplength <- Biostrings:::normargMaxMismatch(max.gaplength)
            max.Lmismatch <- Biostrings:::normargMaxMismatch(max.Lmismatch)
            max.Rmismatch <- Biostrings:::normargMaxMismatch(max.Rmismatch)
            with.Lindels <- Biostrings:::normargWithIndels(with.Lindels)
            with.Rindels <- Biostrings:::normargWithIndels(with.Rindels)
            Lfixed <- Biostrings:::normargFixed(Lfixed, DNAStringSet())
            Rfixed <- Biostrings:::normargFixed(Rfixed, DNAStringSet())
            
            posLPattern <- Lpattern
            negLPattern <- reverseComplement(posLPattern)
            posRPattern <- Rpattern
            negRPattern <- reverseComplement(posRPattern)
            
            bsParams <-
              new("BSParams", X = subject, FUN = matchFUN, exclude = exclude,
                  simplify = TRUE, userMask = userMask,
                  invertUserMask = invertUserMask 
                  )
            matches <-
              bsapply(bsParams,
                      posLPattern = posLPattern, negLPattern = negLPattern,
                      posRPattern = posRPattern, negRPattern = negRPattern,
                      max.gaplength = max.gaplength,
                      max.Lmismatch = max.Lmismatch,
                      max.Rmismatch = max.Rmismatch, 
                      with.Lindels = with.Lindels, with.Rindels = with.Rindels,
                      Lfixed = Lfixed, Rfixed = Rfixed)
            nms <- names(matches)
            matches <- do.call(c, unname(matches))
            names(matches) <- nms
            matches
          }
          )
