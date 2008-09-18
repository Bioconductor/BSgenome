\name{BSgenomeForge}

\alias{BSgenomeForge}
\alias{forgeSeqFiles}
\alias{forgeMaskFiles}

\title{The BSgenomeForge functions}

\description{
  The BSgenomeForge functions provide support during the process of
  forging a "Biostring-based genome data package" (aka "BSgenome data package")
  like the BSgenome.Hsapiens.UCSC.hg18 package.
  They are intended to be used only by users that want to make a new
  BSgenome data package, not by "regular" users of these packages.

  There are currently two BSgenomeForge functions: \code{forgeSeqFiles}
  and \code{forgeMaskFiles}. They must be used in this order.

  See the BSgenomeForge vignette in this package
  (\code{vignette("BSgenomeForge")}) for an extensive coverage
  of this topic.
}

\usage{
  forgeSeqFiles(seqnames, mseqnames=NULL, prefix="", suffix=".fa",
                srcdir=".", destdir=".", verbose=TRUE)
  forgeMaskFiles(seqnames, masks_per_seq,
                 seqdir=".", srcdir=".", destdir=".",
                 srctype1="gap", prefix1="", suffix1="_gap.txt", verbose=TRUE)
}

\arguments{
  \item{seqnames}{
    A character vector containing the names of the single sequences to
    forge for \code{forgeSeqFiles}, or for which masks must be forged for
    \code{forgeMaskFiles}.
  }
  \item{mseqnames}{
    A character vector containing the names of the multiple sequences to forge.
  }
  \item{prefix, suffix}{
    A single string specifying the prefix or suffix to add to each of the names
    given in \code{seqnames} and \code{mseqnames} in order to get the name of the
    corresponding source file.
    See the \code{srcdir} argument below for what the source files are.
  }
  \item{srcdir}{
    A single string indicating the path to the source folder i.e.
    to the folder containing the source files.

    For \code{forgeSeqFiles}, the source files must be FASTA files containing
    the original sequences to load and turn into Biostrings objects.

    For \code{forgeMaskFiles}, the source files must be UCSC "gap" or NCBI "agp"
    files for the mask of "assembly gaps", RepeatMasker .out files for the
    "RepeatMasker" mask, and Tandem Repeats Finder .bed files for
    the "Tandem Repeats Finder" mask.
  }
  \item{destdir}{
    A single string indicating the path to the destination folder i.e.
    to the folder where the serialized Biostrings objects derived from
    the source files should be saved.
    Both \code{forgeSeqFiles} and \code{forgeMaskFiles} will produce one
    .rda file per sequence.

    For \code{forgeSeqFiles}, each .rda file will contain a serialized
    \link[Biostrings]{DNAString} or \link[Biostrings]{DNAStringSet} object.
    The standard place for these files is the inst/extdata/ subfolder of the
    source tree of the BSgenome data packages that is under construction.

    For \code{forgeMaskFiles}, each .rda file will contain a serialized
    \link[IRanges]{MaskCOllection} object.
    The standard place for these files is the data/ subfolder of the
    source tree of the BSgenome data packages that is under construction.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{masks_per_seq}{
    [TODO: Document me]
  }
  \item{seqdir}{
    [TODO: Document me]
  }
  \item{srctype1}{
    [TODO: Document me]
  }
  \item{prefix1, suffix1}{
    [TODO: Document me]
  }
}

\author{H. Pages}

\examples{
  forgeSeqFiles("chrM", prefix="ce2", suffix=".fa",
                srcdir=system.file("extdata", package="BSgenome"),
                destdir=tempdir())
  load(file.path(tempdir(), "chrM.rda"))
  chrM
}

\keyword{manip}
