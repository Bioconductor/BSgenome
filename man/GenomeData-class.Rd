\name{GenomeData-class}
\docType{class}

\alias{GenomeData-class}



% Accessors:
\alias{provider,GenomeData-method}
\alias{organism,GenomeData-method}
\alias{providerVersion,GenomeData-method}

% Constructor:
\alias{GenomeData}

% Show:
\alias{show,GenomeData-method}

% Coercion:
\alias{coerce,GenomeData,data.frame-method}

\title{Data on the genome}
\description{\code{GenomeData} formally represents genomic data as a
  list, with one element per chromosome in the genome.} 

\details{
  This class facilitates storing data on the genome by formalizing a set
  of metadata fields for storing the organism (e.g. Mmusculus), genome
  build provider (e.g. UCSC), and genome build version (e.g. mm9).

  The data is represented as a list, with one element per chromosome (or
  really any sequence, like a gene). There are no constraints as to the
  data type of the elements.

  Note that as an \code{\linkS4class{AnnotatedList}}, it is possible to
  store chromosome-level data (e.g. the lengths) in the
  \code{elementMetadata} slot. The \code{organism}, \code{provider} and
  \code{providerVersion} are all stored in the the \code{AnnotatedList}
  \code{metadata}, so they may be retrieved in list form by calling
  \code{metadata(x)}.
}

\section{Accesor methods}{
  In the code snippets below, \code{x} is a \code{GenomeData} object.

  \describe{
    \item{}{\code{organism(x)}:
      Get the single string indicating the
      organism, if specified, otherwise \code{NULL}.
    }
    \item{}{\code{provider(x)}:
      Get the single string indicating the
      genome build provider, if specified, otherwise \code{NULL}.
    }
    \item{}{\code{providerVersion(x)}:
      Get the single string indicating the
      genome build version, if specified, otherwise \code{NULL}.
    }
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomeData(elements = list(), providerVersion = NULL,
        organism = NULL, provider = NULL, metadata = list(),
        elementMetadata = NULL, ...)}:
      Creates a \code{GenomeData} with the elements from the
      \code{elements} parameter, a list. The other arguments correspond to
      the metadata fields, and, with the exception of \code{elementMetadata},
      should all be either single strings or \code{NULL} (unspecified).
      Additional global metadata elements may be
      passed in \code{metadata}, in list-form, and via \code{...}. The
      elements in \code{metadata} are always overridden by the explicit
      arguments, like \code{organism} and those in \code{...}. 
      \code{elementMetadata} should be an \code{\linkS4class{XDataFrame}}
      or \code{NULL}. 
    }
  }
}

\section{Coercion}{
  \describe{
    \item{}{\code{as(from, "data.frame")}: Coerces each subelement to a
      data frame, and binds them into a single data frame with an
      additional column indicating chromosome
    }
  }
}

\author{ Michael Lawrence }

\seealso{
  \linkS4class{GenomeDataList}, a container of this class and useful for
  storing data on multiple samples.

  \linkS4class{AnnotatedList}, the base of this class.
}

\examples{
gd <- GenomeData(list(chr1 = IRanges(1, 10), chrX = IRanges(2, 5)),
                 organism = "Mmusculus", provider = "UCSC",
                 providerVersion = "mm9")
organism(gd)
providerVersion(gd)
provider(gd)
gd[["chr1"]] # get data for chromsome 1
}

\keyword{methods}
\keyword{classes}
