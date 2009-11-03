## vmatchPattern/vcountPattern for BSgenome
setMethod("vmatchPattern", "BSgenome",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             exclude = "")
    {
        matchFUN <- function(posPattern, negPattern, chr,
                             algorithm = algorithm,
                             max.mismatch = max.mismatch,
                             min.mismatch = min.mismatch,
                             with.indels = with.indels,
                             fixed = fixed) {
            if (is(chr, "MaskedXString"))
                active(masks(chr)) <- FALSE
            posMatches <-
              matchPattern(pattern = posPattern, subject = chr,
                           algorithm = algorithm,
                           max.mismatch = max.mismatch,
                           min.mismatch = min.mismatch,
                           with.indels = with.indels,
                           fixed = fixed)
            negMatches <-
              matchPattern(pattern = negPattern, subject = chr,
                           algorithm = algorithm,
                           max.mismatch = max.mismatch,
                           min.mismatch = min.mismatch,
                           with.indels = with.indels,
                           fixed = fixed)
            list("ranges" =
                 IRangesList("+" = as(posMatches, "IRanges"),
                             "-" = as(negMatches, "IRanges")),
                 "string" =
                 CharacterList("+" = as.character(posMatches),
                               "-" = as.character(negMatches)))
        }

        if (!is(pattern, "DNAString"))
            pattern <- DNAString(pattern)
        algorithm <- Biostrings:::.normargAlgorithm(algorithm)
        if (Biostrings:::.is.character.algo(algorithm)) 
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
          new("BSParams", X = subject, FUN = matchFUN, exclude = exclude,
              simplify = FALSE)
        matches <-
         bsapply(bsParams, posPattern = posPattern, negPattern = negPattern,
                 algorithm = algorithm,
                 max.mismatch = max.mismatch, min.mismatch = min.mismatch,
                 with.indels = with.indels, fixed = fixed)
        ans <-
          do.call(c,
                  lapply(seq_len(length(matches)),
                         function(i) {
                             RangedData(ranges =
                                        unlist(matches[[i]][["ranges"]],
                                               use.names = FALSE),
                                        strand =
                                          strand(rep(c("+", "-"),
                                    elementLengths(matches[[i]][["ranges"]]))),
                                        space = names(matches)[i])
                         }))
        string <-
          factor(do.call(c,
                         unname(lapply(matches, function(x)
                                       unlist(x[["string"]],
                                              use.names = FALSE)))))
        ans$string <- DNAStringSet(levels(string))[as.integer(string)]
        ans
    }
)

setMethod("vcountPattern", "BSgenome",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             exclude = "")
    {
        countFUN <- function(posPattern, negPattern, chr,
                             algorithm = algorithm,
                             max.mismatch = max.mismatch,
                             min.mismatch = min.mismatch,
                             with.indels = with.indels,
                             fixed = fixed,
                             exclude = exclude) {
            if (is(chr, "MaskedXString"))
                active(masks(chr)) <- FALSE
            data.frame(strand = strand(c("+", "-")),
                       count =
                            c(countPattern(pattern = posPattern, subject = chr,
                                           algorithm = algorithm,
                                           max.mismatch = max.mismatch,
                                           min.mismatch = min.mismatch,
                                           with.indels = with.indels,
                                           fixed = fixed),
                              countPattern(pattern = negPattern, subject = chr,
                                           algorithm = algorithm,
                                           max.mismatch = max.mismatch,
                                           min.mismatch = min.mismatch,
                                           with.indels = with.indels,
                                           fixed = fixed)))
        }

        if (!is(pattern, "DNAString"))
            pattern <- DNAString(pattern)
        algorithm <- Biostrings:::.normargAlgorithm(algorithm)
        if (Biostrings:::.is.character.algo(algorithm)) 
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
              simplify = FALSE)
        counts <-
          bsapply(bsParams, posPattern = posPattern, negPattern = negPattern,
                  algorithm = algorithm,
                  max.mismatch = max.mismatch, min.mismatch = min.mismatch,
                  with.indels = with.indels, fixed = fixed)
        cbind(data.frame(seqname =
                         rep(factor(names(counts), levels = names(counts)),
                             each = 2)),
              do.call(rbind, unname(as.list(counts))))
  }
)


## matchPWM/countPWM for BSgenome
setMethod("matchPWM", "BSgenome",
    function(pwm, subject, min.score = "80%", exclude = "")
    {
        matchFUN <- function(posPWM, negPWM, chr, min.score) {
            if (is(chr, "MaskedXString"))
                active(masks(chr)) <- FALSE
            posMatches <-
              matchPWM(pwm = posPWM, subject = chr, min.score = min.score)
            negMatches <-
              matchPWM(pwm = negPWM, subject = chr, min.score = min.score)
            list("ranges" =
                 IRangesList("+" = as(posMatches, "IRanges"),
                             "-" = as(negMatches, "IRanges")),
                 "string" =
                 CharacterList("+" = as.character(posMatches),
                               "-" = as.character(negMatches)))
        }

        ## checking 'pwm'
        pwm <- Biostrings:::.normargPwm(pwm)
        ## checking 'min.score'
        min.score <- Biostrings:::.normargMinScore(min.score, pwm)

        posPWM <- pwm[DNA_BASES, , drop = FALSE]
        negPWM <- reverseComplement(posPWM)
        bsParams <-
          new("BSParams", X = subject, FUN = matchFUN, exclude = exclude,
              simplify = FALSE)
        matches <-
          bsapply(bsParams, posPWM = posPWM, negPWM = negPWM,
                  min.score = min.score)
        ans <-
          do.call(c,
                  lapply(seq_len(length(matches)),
                         function(i) {
                             RangedData(ranges =
                                        unlist(matches[[i]][["ranges"]],
                                               use.names = FALSE),
                                        strand =
                                          strand(rep(c("+", "-"),
                                    elementLengths(matches[[i]][["ranges"]]))),
                                        space = names(matches)[i])
                         }))
        string <-
          factor(do.call(c,
                         unname(lapply(matches, function(x)
                                       unlist(x[["string"]],
                                              use.names = FALSE)))))
        ans$string <- DNAStringSet(levels(string))[as.integer(string)]
        ans
    }
)

setMethod("countPWM", "BSgenome",
    function(pwm, subject, min.score = "80%", exclude = "")
    {
        countFUN <- function(posPWM, negPWM, chr, min.score) {
            if (is(chr, "MaskedXString"))
                active(masks(chr)) <- FALSE
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
              simplify = FALSE)
        counts <-
          bsapply(bsParams, posPWM = posPWM, negPWM = negPWM,
                  min.score = min.score)
        cbind(data.frame(seqname =
                         rep(factor(names(counts), levels = names(counts)),
                             each = 2)),
              do.call(rbind, unname(as.list(counts))))
    }
)
