\name{GenomeDataList-class}
\docType{class}

\alias{GenomeDataList-class}

% Constructor:
\alias{GenomeDataList}

% Coercion:
\alias{coerce,GenomeDataList,data.frame-method}

\title{List of GenomeData objects}
\description{
  WARNING: Starting with BioC 3.1, GenomeData and GenomeDataList objects 
  are defunct. Note that the GenomeData/GenomeDataList containers
  predate the \link[GenomicRanges]{GRanges}/\link[GenomicRanges]{GRangesList} 
  containers and, most of the times, the latters can be used instead of
  the formers.
  Please let us know on the bioc-devel mailing list
  (\url{http://bioconductor.org/help/mailing-list/}) if you have
  a use case where you think there are significant benefits in
  using GenomeData/GenomeDataList over
  \link[GenomicRanges]{GRanges}/\link[GenomicRanges]{GRangesList},
  or if you have questions or concerns about this.

  \code{GenomeDataList} is a list of
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

\author{ Michael Lawrence }

\seealso{
  The \link[GenomicRanges]{GRanges} and \link[GenomicRanges]{GRangesList}
  classes defined and documented in the \pkg{GenomicRanges} package.

  \linkS4class{GenomeData}, the type of elements stored in this class.

  \linkS4class{SimpleList}
}

\examples{
\dontrun{
gd <- GenomeData(list(chr1 = IRanges(1, 10), chrX = IRanges(2, 5)),
                 organism = "Mmusculus", provider = "UCSC",
                 providerVersion = "mm9")
gdl <- GenomeDataList(list(gd), elementMetadata = DataFrame(induced = TRUE))
gdl[[1]] # get first element
}
}

\keyword{methods}
\keyword{classes}
