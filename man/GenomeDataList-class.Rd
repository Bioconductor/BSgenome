\name{GenomeDataList-class}
\docType{class}

\alias{GenomeDataList-class}

% Constructor:
\alias{GenomeDataList}

\title{List of GenomeData objects}
\description{\code{GenomeDataList} is a list of
  \code{\linkS4class{GenomeData}} objects. It could be useful for
  storing data on multiple experiments or samples.} 

\details{
  This is simply an \code{\linkS4class{AnnotatedList}} that requires all of its
  elements to be instances of \code{GenomeData}.

  One should try to take advantage of the metadata storage facilities
  provided by \code{AnnotatedList}. The \code{elementMetadata} field,
  for example, could be used to store the experimental design, while the
  \code{annotation} field could identify the experimental platform.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomeDataList(elements = list(), annotation = NULL,
                           elementMetadata = NULL)}:
      Creates a \code{GenomeDataList} with the elements from the
      \code{elements} parameter, a list of \code{GenomeData}
      instances. The other arguments correspond to the optional metadata
      stored in \code{\linkS4class{AnnotatedList}}. 
    }
  }
}

\author{ Michael Lawrence }

\seealso{
  \linkS4class{GenomeData}, the type of elements stored in this class.
  \linkS4class{AnnotatedList}, the base of this class.
}

\examples{
gd <- GenomeData(list(chr1 = IRanges(1, 10), chrX = IRanges(2, 5)),
                 organism = "Mmusculus", provider = "UCSC",
                 providerVersion = "mm9")
gdl <- GenomeDataList(list(gd), elementMetadata = XDataFrame(induced = TRUE))
gdl[[1]] # get first element
}

\keyword{methods}
\keyword{classes}
