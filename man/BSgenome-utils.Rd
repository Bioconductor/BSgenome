\name{BSgenome-utils}

\alias{BSgenome-utils}

\alias{vmatchPattern,BSgenome-method}
\alias{vcountPattern,BSgenome-method}
\alias{vmatchPDict,BSgenome-method}
\alias{vcountPDict,BSgenome-method}
\alias{matchPWM,BSgenome-method}
\alias{countPWM,BSgenome-method}

\title{BSgenome utilities}

\description{
  Utilities for BSgenome objects.
}

\usage{
\S4method{vmatchPattern}{BSgenome}(pattern, subject, max.mismatch=0, min.mismatch=0,
              with.indels=FALSE, fixed=TRUE, algorithm="auto",
              exclude="", maskList=logical(0), userMask=IRangesList(),
              invertUserMask=FALSE)
\S4method{vcountPattern}{BSgenome}(pattern, subject, max.mismatch=0, min.mismatch=0,
              with.indels=FALSE, fixed=TRUE, algorithm="auto",
              exclude="", maskList=logical(0), userMask=IRangesList(),
              invertUserMask=FALSE)

\S4method{vmatchPDict}{BSgenome}(pdict, subject, max.mismatch=0, min.mismatch=0,
            fixed=TRUE, algorithm="auto", verbose=FALSE,
            exclude="", maskList=logical(0))
\S4method{vcountPDict}{BSgenome}(pdict, subject, max.mismatch=0, min.mismatch=0,
            fixed=TRUE, algorithm="auto", collapse=FALSE,
            weight=1L, verbose=FALSE, exclude="", maskList=logical(0))

\S4method{matchPWM}{BSgenome}(pwm, subject, min.score="80\%", exclude="", maskList=logical(0))
\S4method{countPWM}{BSgenome}(pwm, subject, min.score="80\%", exclude="", maskList=logical(0))
}

\arguments{
  \item{pattern}{
    A \link[Biostrings]{DNAString} object containing the pattern sequence.
  }
  \item{subject}{
    A \link{BSgenome} object containing the subject sequences.
  }
  \item{max.mismatch, min.mismatch}{
    The maximum and minimum number of mismatching letters allowed (see
    \code{?`\link[Biostrings]{lowlevel-matching}`} for the details).
    If non-zero, an inexact matching algorithm is used.
  }
  \item{with.indels}{
    If \code{TRUE} then indels are allowed. In that case, \code{min.mismatch}
    must be \code{0} and \code{max.mismatch} is interpreted as the maximum
    "edit distance" allowed between any pattern and any of its matches
    (see \code{?`\link[Biostrings]{matchPattern}`} for the details).
  }
  \item{fixed}{
    If \code{FALSE} then IUPAC extended letters are interpreted as ambiguities
    (see \code{?`\link[Biostrings]{lowlevel-matching}`} for the details).
  }
  \item{algorithm}{
    For \code{vmatchPattern} and \code{vcountPattern} one of the following:
    \code{"auto"}, \code{"naive-exact"}, \code{"naive-inexact"},
    \code{"boyer-moore"}, \code{"shift-or"}, or \code{"indels"}.

    For \code{vmatchPDict} and \code{vcountPDict} one of the following:
    \code{"auto"}, \code{"naive-exact"}, \code{"naive-inexact"},
    \code{"boyer-moore"}, or \code{"shift-or"}.
  }
  \item{exclude}{
    A character vector with strings that will be used to filter out
    chromosomes whose names match these strings.
  }
  \item{maskList}{
    A named logical vector of maskStates preferred when used with a
    BSGenome object.  When using the bsapply function, the masks will
    be set to the states in this vector.
  }
  \item{userMask}{
    An \linkS4class{IntegerRangesList}, containing a mask to be applied
    to each chromosome. See \code{\link{bsapply}}.
  }
  \item{invertUserMask}{
    Whether the \code{userMask} should be inverted.
  }
  \item{collapse, weight}{
    ignored arguments.
  }
  \item{pdict}{
    A \link[Biostrings]{PDict} or \link[Biostrings]{DNAStringSet} object
    containing the pattern sequences.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{pwm}{
    A numeric matrix with row names A, C, G and T representing a Position
    Weight Matrix.
  }
  \item{min.score}{
    The minimum score for counting a match.
    Can be given as a character string containing a percentage (e.g.
    \code{"85\%"}) of the highest possible score or as a single number.
  }
}

\value{
  A \link[GenomicRanges]{GRanges} object for \code{vmatchPattern}.
  \code{genome} and \code{seqinfo} information from "subject" are
  propagated to the return object.

  A data.frame object for \code{vcountPattern} and \code{countPWM}
  with three columns: "seqname" (factor), "strand" (factor), and
  "count" (integer).

  A \link[GenomicRanges]{GRanges} object for \code{vmatchPDict} with
  one metadata column: "index", which represents a mapping to a
  position in the original pattern dictionary. \code{genome} and
  \code{seqinfo} information from "subject" are propagated.

  A \link[S4Vectors]{DataFrame} object for \code{vcountPDict} with four
  columns: "seqname" ('factor' Rle), "strand" ('factor' Rle), 
  "index" (integer) and "count" ('integer' Rle). As with \code{vmatchPDict}
  the index column represents a mapping to a position in the original
  pattern dictionary.

  A \link[GenomicRanges]{GRanges} object for \code{matchPWM} with two
  metadata columns: "score" (numeric), and "string" (DNAStringSet).
  \code{genome} and \code{seqinfo} information from "subject" are
  included.

}

\author{P. Aboyoun}

\seealso{
  \code{\link[Biostrings]{matchPattern}},
  \code{\link[Biostrings]{matchPDict}},
  \code{\link[Biostrings]{matchPWM}},
  \code{\link{bsapply}}
}

\examples{
library(BSgenome.Celegans.UCSC.ce2)
data(HNF4alpha)

pattern <- consensusString(HNF4alpha)
vmatchPattern(pattern, Celegans, fixed="subject")
vcountPattern(pattern, Celegans, fixed="subject")

pdict <- PDict(HNF4alpha)
vmatchPDict(pdict, Celegans)
vcountPDict(pdict, Celegans)

pwm <- PWM(HNF4alpha)
matchPWM(pwm, Celegans)
countPWM(pwm, Celegans)
}

\keyword{methods}
\keyword{utilities}
