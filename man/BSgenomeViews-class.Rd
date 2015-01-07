\name{BSgenomeViews-class}
\docType{class}

\alias{class:BSgenomeViews}
\alias{BSgenomeViews-class}

% constructors
\alias{BSgenomeViews}
\alias{Views,BSgenome-method}

% accessors
\alias{subject,BSgenomeViews-method}
\alias{granges,BSgenomeViews-method}
\alias{length,BSgenomeViews-method}
\alias{names,BSgenomeViews-method}
\alias{seqnames,BSgenomeViews-method}
\alias{start,BSgenomeViews-method}
\alias{end,BSgenomeViews-method}
\alias{width,BSgenomeViews-method}
\alias{strand,BSgenomeViews-method}
\alias{ranges,BSgenomeViews-method}
\alias{elementLengths,BSgenomeViews-method}
\alias{seqinfo,BSgenomeViews-method}

% coercion
\alias{coerce,BSgenomeViews,DNAStringSet-method}
\alias{coerce,BSgenomeViews,XStringSet-method}
\alias{as.character,BSgenomeViews-method}
\alias{as.data.frame,BSgenomeViews-method}

% subsetting
\alias{extractROWS,BSgenomeViews-method}
\alias{getListElement,BSgenomeViews-method}

% DNAStringSet methods
\alias{seqtype,BSgenomeViews-method}
\alias{nchar,BSgenomeViews-method}
\alias{unlist,BSgenomeViews-method}
\alias{alphabetFrequency,BSgenomeViews-method}
\alias{hasOnlyBaseLetters,BSgenomeViews-method}
\alias{uniqueLetters,BSgenomeViews-method}
\alias{letterFrequency,BSgenomeViews-method}
\alias{letterFrequencyInSlidingView,BSgenomeViews-method}
\alias{oligonucleotideFrequency,BSgenomeViews-method}
\alias{nucleotideFrequencyAt,BSgenomeViews-method}
\alias{consensusMatrix,BSgenomeViews-method}
\alias{consensusString,BSgenomeViews-method}

% displaying
\alias{show,BSgenomeViews-method}

\title{BSgenomeViews objects}

\description{
  The BSgenomeViews class is a container for storing a set of genomic
  positions on a \link{BSgenome} object, called the "subject" in this
  context.
}

\usage{
## Constructor
## ------------

BSgenomeViews(subject, granges)

## Accessors
## ---------

\S4method{subject}{BSgenomeViews}(x)
\S4method{granges}{BSgenomeViews}(x, use.mcols=FALSE)

\S4method{length}{BSgenomeViews}(x)
\S4method{names}{BSgenomeViews}(x)
\S4method{seqnames}{BSgenomeViews}(x)
\S4method{start}{BSgenomeViews}(x)
\S4method{end}{BSgenomeViews}(x)
\S4method{width}{BSgenomeViews}(x)
\S4method{strand}{BSgenomeViews}(x)
\S4method{ranges}{BSgenomeViews}(x, use.mcols=FALSE)
\S4method{elementLengths}{BSgenomeViews}(x)
\S4method{seqinfo}{BSgenomeViews}(x)

## DNAStringSet methods
## --------------------

\S4method{seqtype}{BSgenomeViews}(x)

\S4method{nchar}{BSgenomeViews}(x, type="chars", allowNA=FALSE)

\S4method{unlist}{BSgenomeViews}(x, recursive=TRUE, use.names=TRUE)

\S4method{alphabetFrequency}{BSgenomeViews}(x, as.prob=FALSE, ...)

\S4method{hasOnlyBaseLetters}{BSgenomeViews}(x)

\S4method{uniqueLetters}{BSgenomeViews}(x)

\S4method{letterFrequency}{BSgenomeViews}(x, letters, OR="|", as.prob=FALSE, ...)
\S4method{letterFrequencyInSlidingView}{BSgenomeViews}(x, view.width, letters, OR="|", as.prob=FALSE)

\S4method{oligonucleotideFrequency}{BSgenomeViews}(x, width, step=1,
                         as.prob=FALSE, as.array=FALSE,
                         fast.moving.side="right", with.labels=TRUE, ...)

\S4method{nucleotideFrequencyAt}{BSgenomeViews}(x, at, as.prob=FALSE, as.array=TRUE,
                      fast.moving.side="right", with.labels=TRUE, ...)

\S4method{consensusMatrix}{BSgenomeViews}(x, as.prob=FALSE, shift=0L, width=NULL, ...)

\S4method{consensusString}{BSgenomeViews}(x, ...)
}

\arguments{
  \item{subject}{
    A \link{BSgenome} object or the name of a reference genome specified
    in a way that is accepted by the \code{\link{getBSgenome}} function.
    In that case the corresponding BSgenome data package needs to be already
    installed (see \code{?\link{getBSgenome}} for the details).
  }
  \item{granges}{
    A \link[GenomicRanges]{GRanges} object containing ranges relative to
    the genomic sequences stored in \code{subject}.
  }
  \item{x}{
    A BSgenomeViews object.
  }
  \item{use.mcols}{
    \code{TRUE} or \code{FALSE} (the default).
    Whether the metadata columns on \code{x} (accessible with \code{mcols(x)})
    should be propagated to the returned object or not.
  }
  \item{type, allowNA, recursive, use.names}{
    Ignored.
  }
  \item{as.prob, view.width, letters, OR, width}{
    See \code{?\link[Biostrings]{alphabetFrequency}} and
    \code{?\link[Biostrings]{oligonucleotideFrequency}} in the
    \pkg{Biostrings} package.
  }
  \item{...}{
    Further arguments to be passed to or from other methods.
  }
  \item{step, as.array, fast.moving.side, with.labels, at}{
    See \code{?\link[Biostrings]{oligonucleotideFrequency}} in the
    \pkg{Biostrings} package.
  }
  \item{shift}{
    See \code{?\link[Biostrings]{consensusMatrix}} in the
    \pkg{Biostrings} package.
  }
}

\section{Constructors}{
  \describe{
    \item{}{
      \code{BSgenomeViews(subject, granges)}:
      Make a BSgenomeViews object by putting the views specified by
      \code{granges} on top of the genomic sequences stored in \code{subject}.
      See above for how argument \code{subject} and \code{granges} should be
      specified.
    }
    \item{}{
      \code{Views(subject, granges)}: Equivalent to
      \code{BSgenomeViews(subject, granges)}. Provided for convenience.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a BSgenomeViews object.
  \describe{
    \item{}{
      \code{subject(x)}: Return the \link{BSgenome} object containing the
      full genomic sequences on top of which the views in \code{x} are
      defined.
    }
    \item{}{
      \code{granges(x, use.mcols=FALSE)}: Return the genomic ranges of the
      views as a \link[GenomicRanges]{GRanges} object. These ranges are
      relative to the genomic sequences stored in \code{subject(x)}.
    }
    \item{}{
      \code{length(x)}: The number of views in \code{x}.
    }
    \item{}{
      \code{names(x)}: The names of the views in \code{x}.
    }
    \item{}{
      \code{seqnames(x)}, \code{start(x)}, \code{end(x)}, \code{width(x)},
      \code{strand(x)}: Equivalent to \code{seqnames(granges(x))},
      \code{start(granges(x))}, \code{end(granges(x))},
      \code{width(granges(x))}, \code{strand(granges(x))}, respectively.
    }
    \item{}{
      \code{ranges(x, use.mcols=FALSE)}: Equivalent to
      \code{ranges(granges(x, use.mcols), use.mcols)}.
    }
    \item{}{
      \code{elementLengths(x)}: Equivalent to \code{width(x)}.
    }
    \item{}{
      \code{seqinfo(x)}: Equivalent to \code{seqinfo(subject(x))} and to
      \code{seqinfo(granges(x))} (both are guaranteed to be the same).
      See \code{?\link[GenomeInfoDb]{seqinfo}} in the \pkg{GenomeInfoDb}
      package for more information.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a BSgenomeViews object.
  \describe{
    \item{}{
      \code{as(x, "DNAStringSet")}: Turn \code{x} into a
      \link[Biostrings]{DNAStringSet} object by extxracting the DNA sequence
      corresponding to each view. Alternatively \code{as(x, "XStringSet")}
      can be used for this, and is equivalent to \code{as(x, "DNAStringSet")}.
    }
    \item{}{
      \code{as.character(x)}: Equivalent to
      \code{as.character(as(x, "DNAStringSet"))}.
    }
    \item{}{
      \code{as.data.frame(x)}: Turn \code{x} into a data.frame.
    }
  }
}

\section{Subsetting}{
  \describe{
    \item{}{
      \code{x[i]}: Select the views specified by \code{i}.
    }
    \item{}{
      \code{x[[i]]}: Extract the one view specified by \code{i}.
    }
  }
}

\section{DNAStringSet methods}{
  For convenience, some methods defined for \link[Biostrings]{DNAStringSet}
  objects in the \pkg{Biostrings} package can be used directly on a
  BSgenomeViews object. In that case, everything happens like if the
  BSgenomeViews object \code{x} was turned into a
  \link[Biostrings]{DNAStringSet} object (with \code{as(x, "DNAStringSet")})
  before it's passed to the method for \link[Biostrings]{DNAStringSet} objects.

  At the moment, the list of such methods is:
  \code{\link[Biostrings]{seqtype}},
  \code{\link[Biostrings]{nchar,XStringSet-method}},
  \code{\link[Biostrings]{unlist,XStringSet-method}},
  \code{\link[Biostrings]{alphabetFrequency}},
  \code{\link[Biostrings]{hasOnlyBaseLetters}},
  \code{\link[Biostrings]{uniqueLetters}},
  \code{\link[Biostrings]{letterFrequency}},
  \code{\link[Biostrings]{letterFrequencyInSlidingView}},
  \code{\link[Biostrings]{oligonucleotideFrequency}},
  \code{\link[Biostrings]{nucleotideFrequencyAt}},
  \code{\link[Biostrings]{consensusMatrix}},
  and \code{\link[Biostrings]{consensusString}}.

  See the corresponding man page in the \pkg{Biostrings} package for a
  description of these methods.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item The \link{BSgenome} class.

    \item The \link[GenomicRanges]{GRanges} class in the \pkg{GenomicRanges}
          package.

    \item The \link[Biostrings]{DNAStringSet} class in the \pkg{Biostrings}
          package.

    \item The \link[GenomeInfoDb]{seqinfo}} and related getters in the
          \pkg{GenomeInfoDb} package for getting the sequence information
          stored in an object.

    \item \link[GenomicFeatures]{TxDb} objects in the \pkg{GenomicFeatures}
          package.
  }
}

\examples{
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
ex <- exons(txdb, columns=c("exon_id", "tx_name", "gene_id"))
v <- Views(genome, ex)
v

subject(v)
granges(v)
seqinfo(v)
as(v, "DNAStringSet")

v10 <- v[1:10]  # select the first 10 views
subject(v10)    # same as subject(v)
granges(v10)
seqinfo(v10)    # same as seqinfo(v)
as(v10, "DNAStringSet")
alphabetFrequency(v10)
alphabetFrequency(v10, collapse=TRUE)

## BSgenomeViews objects are list-like objects. That is, the
## BSgenomeViews class derives from List and typical list/List
## operations (e.g. [[, elementLengths(), unlist(), elementType(),
## etc...) work on these objects:
is(v, "List")  # TRUE
v[[2]]
head(elementLengths(v))  # elementLengths(v) is the same as width(v)
unlist(v) 
elementType(v)
}

\keyword{methods}
\keyword{classes}
