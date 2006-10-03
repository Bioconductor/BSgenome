\name{BSgenome}
\docType{class}

% Classes
\alias{class:BSgenome}
\alias{BSgenome-class}

% Constructor-like functions and generics:
\alias{initialize,BSgenome-method}

% Standard generic methods:
\alias{show,BSgenome-method}
\alias{[[,BStringViews-method}
\alias{[[<-,BStringViews-method}
\alias{$,BSgenome-method}

% Other functions and generics:
\alias{unload}
\alias{unload,BSgenome-method}


\title{The BSgenome class}

\description{
  A container for the complete genome sequence of a given specie.
}

\details{
    [TODO: Put some details here]
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x} is a BSgenome object
  and \code{name} is the name of a chromosome (character-string).

  \describe{
    \item{}{
      \code{x[[name]]}:
      [TODO: Document me]
    }
    \item{}{
      \code{x$name}:
      [TODO: Document me]
    }
  }
}

\section{Other functions and generics}{
  In the code snippets below,
  \code{x} is a BSgenome object
  and \code{name} is the name of a chromosome (character-string).

  \describe{
    \item{}{
      \code{unload(x, name)}:
      [TODO: Document me]
    }
  }
}


\author{H. Pages}

\seealso{
  \link[Biostrings]{BStringViews}, \link[Biostrings]{DNAString}
}

\examples{
  ## [TODO: Put some examples here]
}

\keyword{methods}
\keyword{classes}
