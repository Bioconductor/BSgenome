\name{BSgenomeForge}

\alias{BSgenomeForge}

\alias{BSgenomeDataPkgSeed}
\alias{class:BSgenomeDataPkgSeed}
\alias{BSgenomeDataPkgSeed-class}

\alias{forgeBSgenomeDataPkg}
\alias{forgeBSgenomeDataPkg,BSgenomeDataPkgSeed-method}
\alias{forgeBSgenomeDataPkg,list-method}
\alias{forgeBSgenomeDataPkg,character-method}
\alias{forgeSeqlengthsFile}
\alias{forgeSeqFiles}
\alias{forgeMasksFiles}


\title{The BSgenomeForge functions}

\description{
  Set of functions for creating a BSgenome data package.

  These functions are intended for users who want to make a new
  BSgenome data package, not for regular users of these packages.

  See the BSgenomeForge vignette in this package
  (\code{vignette("BSgenomeForge")}) for an extensive coverage
  of this topic.
}

\usage{
  ## Top-level BSgenomeForge function:

  forgeBSgenomeDataPkg(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)

  ## Low-level BSgenomeForge functions:

  forgeSeqlengthsFile(seqnames, prefix="", suffix=".fa",
                      seqs_srcdir=".", seqs_destdir=".", verbose=TRUE)

  forgeSeqFiles(seqnames, mseqnames=NULL, prefix="", suffix=".fa",
                seqs_srcdir=".", seqs_destdir=".", verbose=TRUE)

  forgeMasksFiles(seqnames, nmask_per_seq,
                  seqs_destdir=".", masks_srcdir=".", masks_destdir=".",
                  AGAPSfiles.type="gap", AGAPSfiles.name=NA,
                  AGAPSfiles.prefix="", AGAPSfiles.suffix="_gap.txt",
                  verbose=TRUE)
}

\arguments{
  \item{x}{
    A BSgenomeDataPkgSeed object or the name of a BSgenome data package seed
    file.
  }
  \item{seqs_srcdir, masks_srcdir}{
    A single string indicating the path to the source folder i.e.
    to the folder containing the source files (only read access is
    needed to this folder).

    There are 2 groups of source files: (1) the source files containing the
    sequence data and (2) the source files containing the mask data.

    \code{seqs_srcdir} must point to the folder containing the first group.
    The source files in this folder must be FASTA files. During the forging
    process, these files will be turned into Biostrings objects.

    \code{masks_srcdir} must point to the folder containing the second group.
    The source files in this folder must be UCSC "gap" or NCBI "agp" files for
    the AGAPS masks (assembly gaps), RepeatMasker .out files for the RM masks
    (RepeatMasker), and Tandem Repeats Finder .bed files for the TRF masks
    (Tandem Repeats Finder).
  }
  \item{destdir}{
    The directory where the BSgenome data package will be created.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{seqnames}{
    A character vector containing the names of the single sequences to forge.
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
  \item{seqs_destdir, masks_destdir}{
    A single string indicating the path to the destination folder i.e.
    to the folder where the serialized Biostrings objects derived from
    the source files will be saved (the folder must already exist, and
    read/write access is required to it).

    \code{forgeSeqlengthsFile} will produce a single .rda file.
    Both \code{forgeSeqFiles} and \code{forgeMasksFiles} will produce one
    .rda file per sequence.

    For \code{forgeSeqFiles} (\code{seqs_destdir} argument), each .rda file
    will contain a serialized \link[Biostrings]{DNAString}
    or \link[Biostrings]{DNAStringSet} object.
    The standard place for these files is the inst/extdata/ subfolder of the
    source tree of the BSgenome data packages that is under construction.

    For \code{forgeMasksFiles} (\code{masks_destdir} argument), each .rda file
    will contain a serialized \link[IRanges]{MaskCOllection} object.
    The standard place for these files is the data/ subfolder of the
    source tree of the BSgenome data packages that is under construction.

    In addition to the source files located in \code{masks_srcdir},
    \code{forgeMasksFiles} also needs read access to the files produced
    by \code{forgeSeqFiles}. This is specified thru the \code{seqs_destdir}
    argument of \code{forgeMasksFiles}.
  }
  \item{nmask_per_seq}{
    [TODO: Document me]
  }
  \item{AGAPSfiles.type, AGAPSfiles.name, AGAPSfiles.prefix, AGAPSfiles.suffix}{
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
