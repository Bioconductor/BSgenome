\name{GenomeDescription-class}
\docType{class}

% Classes
\alias{class:GenomeDescription}
\alias{GenomeDescription-class}
\alias{GenomeDescription}

% Accessor methods:
\alias{organism}
\alias{organism,GenomeDescription-method}
\alias{species}
\alias{species,GenomeDescription-method}
\alias{provider}
\alias{provider,GenomeDescription-method}
\alias{providerVersion}
\alias{providerVersion,GenomeDescription-method}
\alias{releaseDate}
\alias{releaseDate,GenomeDescription-method}
\alias{releaseName}
\alias{releaseName,GenomeDescription-method}
\alias{bsgenomeName}
\alias{bsgenomeName,GenomeDescription-method}
\alias{seqinfo,GenomeDescription-method}

% Standard generic methods:
\alias{show,GenomeDescription-method}


\title{GenomeDescription objects}

\description{
  A GenomeDescription object holds the meta information describing
  a given genome.
}

\details{
  In general the user will not need to manipulate directly a
  GenomeDescription instance but will manipulate instead a higher-level
  object that belongs to a class containing the GenomeDescription class.
  For example the top-level object defined in any BSgenome data package
  is a \link{BSgenome} object. But because the \link{BSgenome} class
  contains the GenomeDescription class, it is also a GenomeDescription
  object and can therefore be treated as such. In other words all the
  methods described below will work on it.
}

\section{Accessor methods}{
  In the code snippets below,
  \code{x} is a GenomeDescription object.

  \describe{
    \item{}{
      \code{organism(x)}:
      Return the target organism for this genome e.g. \code{"Homo sapiens"},
      \code{"Mus musculus"}, \code{"Caenorhabditis elegans"}, etc...
    }
    \item{}{
      \code{species(x)}:
      Return the target species for this genome e.g. \code{"Human"},
      \code{"Mouse"}, \code{"Worm"}, etc...
    }
    \item{}{
      \code{provider(x)}:
      Return the provider of this genome e.g. \code{"UCSC"}, \code{"BDGP"},
      \code{"FlyBase"}, etc...
    }
    \item{}{
      \code{providerVersion(x)}:
      Return the provider-side version of this genome.
      For example UCSC uses versions \code{"hg18"}, \code{"hg17"}, etc...
      for the different Builds of the Human genome.
    }
    \item{}{
      \code{releaseDate(x)}:
      Return the release date of this genome e.g. \code{"Mar. 2006"}.
    }
    \item{}{
      \code{releaseName(x)}:
      Return the release name of this genome, which is generally made of the
      name of the organization who assembled it plus its Build version.
      For example, UCSC uses \code{"hg18"} for the version of the
      Human genome corresponding to the Build 36.1 from NCBI hence
      the release name for this genome is \code{"NCBI Build 36.1"}.
    }
    \item{}{
      \code{bsgenomeName(x)}:
      Uses the meta information stored in \code{x} to make the name of
      the corresponding BSgenome data package (see \link{available.genomes}
      for details about the naming scheme used for those packages).
      Of course there is no guarantee that a package with that name
      actually exists.
    }
    \item{\code{seqinfo(x)}}{
      Gets information about the genome sequences. This information
      is returned in a \link[GenomicRanges]{Seqinfo} object.
      Each part of the information can be retrieved separately with
      \code{seqnames(x)}, \code{seqlengths(x)}, and \code{isCircular(x)},
      respectively, as described below.
    }
    \item{\code{seqnames(x)}}{
      Gets the names of the genome sequences.
      \code{seqnames(x)} is equivalent to \code{seqnames(seqinfo(x))}.
    }
    \item{\code{seqlengths(x)}}{
      Gets the lengths of the genome sequences.
      \code{seqlengths(x)} is equivalent to \code{seqlengths(seqinfo(x))}.
    }
    \item{\code{isCircular(x)}}{
      Returns the circularity flags of the genome sequences.
      \code{isCircular(x)} is equivalent to \code{isCircular(seqinfo(x))}.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{available.genomes}},
  \link[GenomicRanges]{Seqinfo-class},
  \link{BSgenome-class}
}

\examples{
library(BSgenome.Celegans.UCSC.ce2)
class(Celegans)
is(Celegans, "GenomeDescription")
provider(Celegans)
seqinfo(Celegans)
gendesc <- as(Celegans, "GenomeDescription")
class(gendesc)
gendesc
provider(gendesc)
seqinfo(gendesc)
bsgenomeName(gendesc)
}

\keyword{methods}
\keyword{classes}
