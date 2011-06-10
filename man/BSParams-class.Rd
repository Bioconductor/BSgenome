\name{BSParams-class}
\docType{class}
\alias{class:BSParams}
\alias{BSParams-class}
\alias{BSParams}

\title{Class "BSParams"}

\description{
  A parameter class for representing all parameters needed for running
  the \code{bsapply} method.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BSParams", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{X}:}{
      a BSgenome object that contains chromosomes that you wish to apply
      FUN on
    }
    \item{\code{FUN}:}{
      the function to apply to each chromosome in the BSgenome object 'X'
    }
    \item{\code{exclude}:}{
      this is a character vector with strings that will be used to filter
      out chromosomes whose names match these strings.
    }
    \item{\code{simplify}:}{
      TRUE/FALSE value to indicate whether or not the function should try
      to simplify the output for you.
    }
    \item{\code{maskList}:}{
      A named logical vector of maskStates preferred when used with a
      BSGenome object.  When using the bsapply function, the masks will
      be set to the states in this vector.
    }    
    \item{\code{motifList}:}{
      A character vector which should contain motifs that the user
      wishes to mask from the sequence.
    }
    \item{\code{userMask}:}{
      A \linkS4class{RangesList} object, where each element masks the
      corresponding chromosome in \code{X}. This allows the user to
      conveniently apply masks besides those included in \code{X}.
    }
    \item{\code{invertUserMask}:}{
      A \code{logical} indicating whether to invert each mask in
      \code{userMask}.
    }
  }
}
\section{Methods}{
  \describe{
    \item{\code{bsapply(p)}}{Performs the function FUN using the
  parameters contained within \code{BSParams}.}
  }
}

\author{Marc Carlson}

\seealso{
  \code{\link{bsapply}}
}

\keyword{classes}
