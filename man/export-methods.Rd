\name{export-methods}

\alias{export-methods}

\alias{writeBSgenomeToFasta}
\alias{export,BSgenome,FastaFile,ANY-method}

\alias{writeBSgenomeToTwobit}
\alias{export,BSgenome,TwoBitFile,ANY-method}


\title{Export a BSgenome object as a FASTA or twoBit file}

\description{
  \code{\link[BiocIO]{export}} methods for \link{BSgenome} objects.

  NOTE: The \code{\link[BiocIO]{export}} generic function and most of
  its methods are defined and documented in the \pkg{BiocIO} package.
  This man page only documents the 2 \code{\link[BiocIO]{export}}
  methods defined in the \pkg{BSgenome} package.
}

\usage{
\S4method{export}{BSgenome,FastaFile,ANY}(object, con, format, compress=FALSE, compression_level=NA, verbose=TRUE)
\S4method{export}{BSgenome,TwoBitFile,ANY}(object, con, format, ...)
}

\arguments{
  \item{object}{
    The \link{BSgenome} object to export.
  }
  \item{con}{
    A \link[rtracklayer]{FastaFile} or \link[rtracklayer]{TwoBitFile} object.

    Alternatively \code{con} can be a single string containing the
    path to a FASTA or twoBit file, in which case either the file extension
    or the \code{format} argument needs to be \code{"fasta"}, \code{"twoBit"},
    or \code{"2bit"}.
    Also note that in this case, the \code{\link[BiocIO]{export}} method
    that is called is either the method with signature
    \code{c("ANY", "character", "missing")} or the method with signature
    \code{c("ANY", "character", "character")}, both defined in the
    \pkg{BiocIO} package. If \code{object} is a \link{BSgenome} object
    and the file extension or the \code{format} argument is
    \code{"fasta"}, \code{"twoBit"}, or \code{"2bit"}, then the flow
    eventually reaches one of 2 methods documented here.
  }
  \item{format}{
    If not missing, should be \code{"fasta"}, \code{"twoBit"}, or
    \code{"2bit"} (case insensitive for \code{"twoBit"} and \code{"2bit"}).
  }
  \item{compress, compression_level}{
    Forwarded to \code{\link{writeXStringSet}}.
    See \code{?\link{writeXStringSet}} for the details.
  }
  \item{verbose}{
    Whether or not the function should display progress.
    \code{TRUE} by default.
  }
  \item{...}{
    Extra arguments. The method for \link[rtracklayer]{TwoBitFile} objects
    forwards them to \code{\link{bsapply}}.
  }
}

\seealso{
  \itemize{
    \item \link{BSgenome} objects.

    \item The \code{\link[BiocIO]{export}} generic function in the
          \pkg{BiocIO} package.

    \item \link[rtracklayer]{FastaFile} and \link[rtracklayer]{TwoBitFile}
          objects in the \pkg{rtracklayer} package.
  }
}

\author{Michael Lawrence}

\examples{
library(BSgenome.Celegans.UCSC.ce2)
genome <- BSgenome.Celegans.UCSC.ce2

## Export as FASTA file.
out1_file <- file.path(tempdir(), "Celegans.fasta")
export(genome, out1_file)

## Export as twoBit file.
out2_file <- file.path(tempdir(), "Celegans.2bit")
export(genome, out2_file)

## Sanity checks:
dna0 <- DNAStringSet(as.list(genome))

system.time(dna1 <- import(out1_file))
stopifnot(identical(names(dna0), names(dna1)) && all(dna0 == dna1))

system.time(dna2 <- import(out2_file))  # importing twoBit is 10-20x
                                        # faster than importing non
                                        # compressed FASTA 
stopifnot(identical(names(dna0), names(dna2)) && all(dna0 == dna2))
}

\keyword{methods}
\keyword{utilities}
