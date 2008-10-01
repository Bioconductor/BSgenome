\name{BSParams-class}
\docType{class}
\alias{BSParams-class}

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
      A list of mask names for the BSGenome object.  For each of these
      mask names, the active status will be toggled to be something other
      than the default state.
    }    
  }
}
% It seems quite likely that this class will become one of many used by
% bsapply, but for now we don't know what those will be yet so one
% class is fine for now
% \section{Extends}{
% Class \code{"BSBaseParams"}, directly.
% }
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
