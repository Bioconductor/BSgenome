\name{GenomeDataList-class}
\docType{class}

\alias{GenomeDataList-class}

% Constructor:
\alias{GenomeDataList}

% Coercion:
\alias{coerce,GenomeDataList,data.frame-method}

% Reduction
\alias{gdreduce,GenomeDataList-method}

\title{List of GenomeData objects}
\description{\code{GenomeDataList} is a list of
  \code{\linkS4class{GenomeData}} objects. It could be useful for
  storing data on multiple experiments or samples.} 

\details{
  This class inherits from \code{\linkS4class{SimpleList}} and requires
  that all of its elements to be instances of \code{GenomeData}.

  One should try to take advantage of the metadata storage facilities
  provided by \code{SimpleList}. The \code{elementMetadata} field,
  for example, could be used to store the experimental design, while
  the \code{metadata} field could store the experimental platform.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomeDataList(listData = list(), metadata = list(),
                           elementMetadata = NULL)}:
      Creates a \code{GenomeDataList} with the elements from the
      \code{listData} parameter, a list of \code{GenomeData}
      instances. The other arguments correspond to the optional metadata
      stored in \code{\linkS4class{SimpleList}}. 
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

\section{Reduction}{
  \describe{
	\item{}{\code{gdreduce(f, ..., init, right=FALSE, accumulate=FALSE, gdArgs=list())}:
	  
	  Currently this method works when a single \code{GenomeDataList} is
	  provided as \code{...}. It successively combines the
	  \code{\linkS4class{GenomeData}} elements in the
	  \code{GenomeDataList} using \code{f}. \code{f} is a function
	  accepting two objects returned by \code{"[["} applied to the
	  succesive \code{GenomeData} elements of \code{...}, returning a
	  single \code{GenomeData} object to be used in subsequent calls to
	  \code{f}. \code{init}, \code{right}, and \code{accumulate} are as
	  described for \code{\link{Reduce}}.  \code{gdArgs} can be used to
	  provide metadata information to the contructor used to create the
	  final \code{GenomeData} object.
  
	}
  }
}

\author{ Michael Lawrence }

\seealso{
  \linkS4class{GenomeData}, the type of elements stored in this class.
  \linkS4class{SimpleList}
}

\examples{
gd <- GenomeData(list(chr1 = IRanges(1, 10), chrX = IRanges(2, 5)),
                 organism = "Mmusculus", provider = "UCSC",
                 providerVersion = "mm9")
gdl <- GenomeDataList(list(gd), elementMetadata = DataFrame(induced = TRUE))
gdl[[1]] # get first element

gdr <- gdreduce(function(x, y) {
    ## "[[" returns IRanges instances, construct a synthetic version
    IRanges(c(start(x), start(y)), c(end(x), end(y)))
}, GenomeDataList(list(gd, gd[2])))
gdr[["chr1"]]
gdr[["chrX"]]
}

\keyword{methods}
\keyword{classes}
