\name{BSgenome-utils}

\alias{BSgenome-utils}

\alias{matchPWM,BSgenome-method}
\alias{countPWM,BSgenome-method}
\alias{vmatchPattern,BSgenome-method}
\alias{vcountPattern,BSgenome-method}


\title{BSgenome utilities}

\description{
  Utilities for BSgenome objects.
}

\usage{
  \S4method{matchPWM}{BSgenome}(pwm, subject, min.score = "80\%", exclude = "")
  \S4method{countPWM}{BSgenome}(pwm, subject, min.score = "80\%", exclude = "")
  \S4method{vmatchPattern}{BSgenome}(pattern, subject, algorithm="auto",
                                     max.mismatch=0, min.mismatch=0,
                                     with.indels=FALSE, fixed=TRUE,
                                     exclude = "")
  \S4method{vcountPattern}{BSgenome}(pattern, subject, algorithm="auto",
                                     max.mismatch=0, min.mismatch=0,
                                     with.indels=FALSE, fixed=TRUE,
                                     exclude = "")
}

\arguments{
  \item{pwm}{
    A numeric matrix with row names A, C, G and T representing a Position
    Weight Matrix.
  }
  \item{pattern}{
    The pattern string.
  }
  \item{subject}{
    A \link{BSgenome} object containing the subject sequences.
  }
  \item{min.score}{
    The minimum score for counting a match.
    Can be given as a character string containing a percentage (e.g.
    \code{"85\%"}) of the highest possible score or as a single number.
  }
  \item{algorithm}{
    One of the following: \code{"auto"}, \code{"naive-exact"},
    \code{"naive-inexact"}, \code{"boyer-moore"}, \code{"shift-or"}
    or \code{"indels"}.
  }
  \item{max.mismatch, min.mismatch}{
    The maximum and minimum number of mismatching letters allowed (see
    \code{?`\link[Biostrings]{lowlevel-matching}`} for the details).
    If non-zero, an inexact matching algorithm is used.
  }
  \item{with.indels}{
    If \code{TRUE} then indels are allowed. In that case, \code{min.mismatch}
    must be \code{0} and \code{max.mismatch} is interpreted as the maximum
    "edit distance" allowed between the pattern and a match.
    Note that in order to avoid pollution by redundant matches,
    only the "best local matches" are returned.
    Roughly speaking, a "best local match" is a match that is locally
    both the closest (to the pattern P) and the shortest.
    More precisely, a substring S' of the subject S is a "best local match" iff:
    \preformatted{
       (a) nedit(P, S') <= max.mismatch
       (b) for every substring S1 of S':
               nedit(P, S1) > nedit(P, S')
       (c) for every substring S2 of S that contains S':
               nedit(P, S2) <= nedit(P, S')
    }
    One nice property of "best local matches" is that their first and last
    letters are guaranteed to be aligned with letters in P (i.e. they match
    letters in P).
  }
  \item{fixed}{
    If \code{FALSE} then IUPAC extended letters are interpreted as ambiguities
    (see \code{?`\link[Biostrings]{lowlevel-matching}`} for the details).
  }
  \item{exclude}{
    A character vector with strings that will be used to filter out
    chromosomes whose names match these strings.
  }
}

\value{
  A \link[IRanges]{RangedData} object for \code{matchPWM} and
  \code{vmatchPattern} with two values columns: strand (factor) and
  string (DNAStringSet).

  A data.frame object for \code{countPWM} and \code{vcountPattern}
  with columns three columns: seqname (factor), strand (factor),
  and count (integer).
}

\author{P. Aboyoun}

\seealso{
  \code{\link[Biostrings]{matchPWM}},
  \code{\link[Biostrings]{matchPattern}},
  \code{\link[BSgenome]{bsapply}}
}

\examples{
  library(BSgenome.Celegans.UCSC.ce2)
  data(HNF4alpha)

  pwm <- PWM(HNF4alpha)
  matchPWM(pwm, Celegans)
  countPWM(pwm, Celegans)

  pattern <- consensusString(HNF4alpha)
  vmatchPattern(pattern, Celegans, fixed = "subject")
  vcountPattern(pattern, Celegans, fixed = "subject")
}

\keyword{methods}
\keyword{utilities}
