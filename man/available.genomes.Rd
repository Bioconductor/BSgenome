\name{available.genomes}

\alias{available.genomes}
\alias{installed.genomes}
\alias{getBSgenome}


\title{Find available/installed genomes}

\description{
  \code{available.genomes} gets the list of BSgenome data packages that
  are available in the Bioconductor repositories for your version of
  R/Bioconductor.

  \code{installed.genomes} gets the list of BSgenome data packages that
  are currently installed on your system.

  \code{getBSgenome} searchs the installed BSgenome data packages for the
  specified genome and returns it as a \link{BSgenome} object.
}

\usage{
available.genomes(splitNameParts=FALSE, type=getOption("pkgType"))

installed.genomes(splitNameParts=FALSE)

getBSgenome(genome, masked=FALSE)
}

\arguments{
  \item{splitNameParts}{
    Whether to split or not the package names in parts. In that case
    the result is returned in a data frame with 5 columns.
  }
  \item{type}{
    Character string indicating the type of package (\code{"source"},
    \code{"mac.binary"} or \code{"win.binary"}) to look for.
  }
  \item{genome}{
    A \link{BSgenome} object, or the full name of an installed BSgenome
    data package, or a short string specifying a genome assembly (a.k.a.
    provider version) that refers unambiguously to an installed BSgenome
    data package.
  }
  \item{masked}{
    \code{TRUE} or \code{FALSE}. Whether to search for the \emph{masked}
    BSgenome object (i.e. the object that contains the masked sequences)
    or not (the default).
  }
}

\details{
  A BSgenome data package contains the full genome sequences for a given
  organism.

  Its name typically has 4 parts (5 parts if it's a \emph{masked} BSgenome
  data package i.e. if it contains masked sequences) separated by a dot
  e.g. \code{BSgenome.Mmusculus.UCSC.mm10} or
  \code{BSgenome.Mmusculus.UCSC.mm10.masked}:
  \enumerate{
    \item The 1st part is always \code{BSgenome}.

    \item The 2nd part is the name of the organism in abbreviated form e.g.
          \code{Mmusculus}, \code{Hsapiens}, \code{Celegans},
          \code{Scerevisiae}, \code{Ecoli}, etc...

    \item The 3rd part is the name of the organisation who provided the
          genome sequences. We formally refer to it as the \emph{provider}
          of the genome. E.g. \code{UCSC}, \code{NCBI}, \code{TAIR}, etc...

    \item The 4th part is the release string or number used by this
          organisation for this particular genome assembly. We formally
          refer to it as the \emph{provider version} of the genome.
          E.g. \code{hg38}, \code{GRCh38}, \code{hg19}, \code{mm10},
          \code{susScr3}, etc...

    \item If the package contains masked sequences, its name has the
          \code{.masked} suffix added to it, which is typically the 5th part.
  }

  A BSgenome data package contains a single top-level object (a
  \link{BSgenome} object) named like the package itself that can be
  used to access the genome sequences.
}

\value{
  For \code{available.genomes} and \code{installed.genomes}: by default
  (i.e. if \code{splitNameParts=FALSE}), a character vector containing
  the names of the BSgenome data packages that are available (for
  \code{available.genomes}) or currently installed (for
  \code{installed.genomes}). If \code{splitNameParts=TRUE}, the list of
  packages is returned in a data frame with one row per package and the
  following columns: \code{pkgname} (character), \code{organism} (factor),
  \code{provider} (factor), \code{provider_version} (character), and
  \code{masked} (logical).

  For \code{getBSgenome}: the \link{BSgenome} object containing the sequences
  for the specified genome. Or an error if the object cannot be found in the
  BSgenome data packages currently installed.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item \link{BSgenome} objects.

    \item \code{\link[utils]{available.packages}}.
  }
}

\examples{
## ---------------------------------------------------------------------
## available.genomes() and installed.genomes()
## ---------------------------------------------------------------------

# What genomes are currently installed:
installed.genomes()

# What genomes are available:
available.genomes()

# Split the package names in parts:
av_gen <- available.genomes(splitNameParts=TRUE)
table(av_gen$organism)
table(av_gen$provider)

# Make your choice and install with:
library(BiocInstaller)
biocLite("BSgenome.Scerevisiae.UCSC.sacCer1")

# Have a coffee 8-)

# Load the package and display the index of sequences for this genome:
library(BSgenome.Scerevisiae.UCSC.sacCer1)
Scerevisiae  # same as BSgenome.Scerevisiae.UCSC.sacCer1

## ---------------------------------------------------------------------
## getBSgenome()
## ---------------------------------------------------------------------

## Specify the full name of an installed BSgenome data package:
genome <- getBSgenome("BSgenome.Celegans.UCSC.ce2")
genome

## Specify a genome assembly (a.k.a. provider version):
genome <- getBSgenome("hg38")
class(genome)  # BSgenome object
providerVersion(genome)
genome$chrM

genome <- getBSgenome("hg38", masked=TRUE)
class(genome)  # MaskedBSgenome object
providerVersion(genome)
genome$chr22
}

\keyword{manip}
