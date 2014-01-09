\name{available.genomes}

\alias{available.genomes}
\alias{installed.genomes}

\title{Find available/installed genomes}

\description{
  \code{available.genomes} gets the list of BSgenome data packages that
  are currently available on the Bioconductor repositories for your
  version of R/Bioconductor.
  \code{installed.genomes} gets the list of BSgenome data packages that
  are already installed on your machine.
}

\usage{
available.genomes(splitNameParts=FALSE, type=getOption("pkgType"))
installed.genomes(splitNameParts=FALSE)
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
}

\details{
  A BSgenome data package contains the full genome sequences for a given
  organism.

  Its name typically has 4 parts (5 parts if it contains masked sequences)
  separated by a dot e.g. \code{BSgenome.Celegans.UCSC.ce2}.
  The 1st part is always \code{BSgenome}, the 2nd part is the name of the
  organism (abbreviated), the 3rd part is the name of the organisation who
  provided the genome sequences, and the 4th part is the release string or
  number used by this organisation for this particular assembly.

  If the package contains masked sequences, its name has the \code{.masked}
  suffix added to it, which is typically the 5th part.

  A BSgenome data package contains a single top-level object (a
  \link[BSgenome]{BSgenome} object) named like the package itself that
  can be used to access the genome sequences.
}

\value{
  A character vector containing the names of the BSgenome data packages
  that are currently available (for \code{available.genomes}),
  or already installed (for \code{installed.genomes}).
}

\author{H. Pages}

\seealso{
  \link{BSgenome-class},
  \code{\link[utils:update.packages]{available.packages}}
}

\examples{
# What genomes are already installed:
installed.genomes()

# What genomes are available:
available.genomes()

# Split the package names in parts:
ag <- available.genomes(TRUE)
table(ag$organism)
table(ag$provider)

# Make your choice and install with:
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Scerevisiae.UCSC.sacCer1")

# Have a coffee ;-)

# Load the package and display the index of sequences for this genome:
library(BSgenome.Scerevisiae.UCSC.sacCer1)
BSgenome.Scerevisiae.UCSC.sacCer1
}

\keyword{manip}
