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
                seqs_srcdir=".", seqs_destdir=".", verbose=TRUE)

  forgeMaskFiles(seqnames, nmask_per_seq,
                 seqs_destdir=".", masks_srcdir=".", masks_destdir=".",
                 mask1.srctype="gap", mask1.prefix="", mask1.suffix="_gap.txt",
                 verbose=TRUE)
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
    See the \code{seqs_srcdir} argument below for what the source files are.
  }
  \item{seqs_srcdir, masks_srcdir}{
    A single string indicating the path to the source folder i.e.
    to the folder containing the source files (only read access is
    needed to this folder).

    For \code{forgeSeqFiles} (\code{seqs_srcdir} argument), the source files
    must be FASTA files containing the original sequences to load and turn
    into Biostrings objects.

    For \code{forgeMaskFiles} (\code{masks_srcdir} argument), the source files
    must be UCSC "gap" or NCBI "agp" files for the mask of "assembly gaps",
    RepeatMasker .out files for the "RepeatMasker" mask,
    and Tandem Repeats Finder .bed files for the "Tandem Repeats Finder" mask.
  }
  \item{seqs_destdir, masks_destdir}{
    A single string indicating the path to the destination folder i.e.
    to the folder where the serialized Biostrings objects derived from
    the source files will be saved (the folder must already exist, and
    read/write access is required to it).
    Both \code{forgeSeqFiles} and \code{forgeMaskFiles} will produce one
    .rda file per sequence.

    For \code{forgeSeqFiles} (\code{seqs_destdir} argument), each .rda file
    will contain a serialized \link[Biostrings]{DNAString}
    or \link[Biostrings]{DNAStringSet} object.
    The standard place for these files is the inst/extdata/ subfolder of the
    source tree of the BSgenome data packages that is under construction.

    For \code{forgeMaskFiles} (\code{masks_destdir} argument), each .rda file
    will contain a serialized \link[IRanges]{MaskCOllection} object.
    The standard place for these files is the data/ subfolder of the
    source tree of the BSgenome data packages that is under construction.

    In addition to the source files located in \code{masks_srcdir},
    \code{forgeMaskFiles} also needs read access to the files produced
    by \code{forgeSeqFiles}. This is specified thru the \code{seqs_destdir}
    argument of \code{forgeMaskFiles}.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{nmask_per_seq}{
    [TODO: Document me]
  }
  \item{mask1.srctype}{
    [TODO: Document me]
  }
  \item{mask1.prefix, mask1.suffix}{
    [TODO: Document me]
  }
}

\author{H. Pages}

\examples{
  forgeSeqFiles("chrM", prefix="ce2", suffix=".fa",
                seqs_srcdir=system.file("extdata", package="BSgenome"),
                seqs_destdir=tempdir())
  load(file.path(tempdir(), "chrM.rda"))
  chrM
}

\keyword{manip}
