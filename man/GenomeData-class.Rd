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

\title{Data on the genome}
\description{\code{GenomeData} formally represents genomic data as a
  list, with one element per chromosome in the genome.} 

\details{
  This class facilitates storing data on the genome by providing
  metadata slots for storing the organism (e.g. Mmusculus), genome build
  provider (e.g. UCSC), and genome build version (e.g. mm9).

  The data is represented as a list, with one element per chromosome (or
  really any sequence, like a gene). There are no constraints as to the
  data type of the elements.

  Note that as an \code{\linkS4class{AnnotatedList}}, it is possible to
  store chromosome-level data (e.g. the lengths) in the
  \code{elementMetadata} slot.
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
        organism = NULL, provider = NULL,
        elementMetadata = NULL, ...)}:
      Creates a \code{GenomeData} with the elements from the
      \code{elements} parameter, a list. The other arguments correspond to
      the metadata slots, and, with the exception of \code{elementMetadata},
      should all be either single strings or \code{NULL} (unspecified).
      \code{elementMetadata} should be an \code{\linkS4class{XDataFrame}}
      or \code{NULL}.
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
