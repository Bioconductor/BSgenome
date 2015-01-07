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
\alias{ranges,BSgenomeViews-method}
\alias{strand,BSgenomeViews-method}
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
\alias{letterFrequencyInSlidingView,BSgenomeViews-method}
\alias{letterFrequency,BSgenomeViews-method}
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
## Constructors
## ------------

BSgenomeViews(genome, granges)

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
\S4method{ranges}{BSgenomeViews}(x, use.mcols=FALSE)
\S4method{strand}{BSgenomeViews}(x)
\S4method{seqinfo}{BSgenomeViews}(x)

## DNAStringSet methods
## --------------------

\S4method{seqtype}{BSgenomeViews}(x)

\S4method{nchar}{BSgenomeViews}(x, type="chars", allowNA=FALSE)

\S4method{unlist}{BSgenomeViews}(x, recursive=TRUE, use.names=TRUE)

\S4method{alphabetFrequency}{BSgenomeViews}(x, as.prob=FALSE, ...)

\S4method{hasOnlyBaseLetters}{BSgenomeViews}(x)

\S4method{uniqueLetters}{BSgenomeViews}(x)

\S4method{letterFrequencyInSlidingView}{BSgenomeViews}(x, view.width, letters, OR="|", as.prob=FALSE)

\S4method{letterFrequency}{BSgenomeViews}(x, letters, OR="|", as.prob=FALSE, ...)

\S4method{oligonucleotideFrequency}{BSgenomeViews}(x, width, step=1,
                         as.prob=FALSE, as.array=FALSE,
                         fast.moving.side="right", with.labels=TRUE, ...)

\S4method{nucleotideFrequencyAt}{BSgenomeViews}(x, at, as.prob=FALSE, as.array=TRUE,
                      fast.moving.side="right", with.labels=TRUE, ...)

\S4method{consensusMatrix}{BSgenomeViews}(x, as.prob=FALSE, shift=0L, width=NULL, ...)

\S4method{consensusString}{BSgenomeViews}(x, ...)
}

\arguments{
  \item{genome}{
    A \link{BSgenome} object or the name of a reference genome specified
    in a way that is accepted by the \code{\link{getBSgenome}} function.
    In that case the corresponding BSgenome data package needs to be already
    installed (see \code{?\link{getBSgenome}} for the details).
  }
  \item{granges}{
    A \link[GenomicRanges]{GRanges} object containing ranges relative to
    \code{genome}.
  }
  \item{x}{
    A BSgenomeViews object.
  }
  \item{use.mcols}{
  }
  \item{type, allowNA}{
  }
  \item{recursive, use.names}{
  }
  \item{as.prob}{
  }
  \item{...}{
    Further arguments to be passed to or from other methods.
  }
  \item{view.width, letters, OR}{
  }
  \item{width, step, as.array, fast.moving.side, with.labels}{
  }
  \item{at}{
  }
  \item{shift}{
  }
}

\section{Constructors}{
  \describe{
    \item{}{
      \code{BSgenomeViews(genome, granges)}:
    }
  }
}

\section{Accessors}{
}

\section{Coercion}{
  In the code snippets below, \code{x} is a BSgenomeViews object.
  \describe{
    \item{}{
      \code{as(x, "DNAStringSet")}:
      Turn \code{x} into a \link[Biostrings]{DNAStringSet} object.
      Alternatively \code{as(x, "XStringSet")} can be used for this, and is
      equivalent to \code{as(x, "DNAStringSet")}.
    }
    \item{}{
      \code{as.character(x)}:
    }
    \item{}{
      \code{as.data.frame(x)}:
    }
  }
}

\section{Subsetting}{
  \describe{
    \item{}{
      \code{x[i]}: Select the views specified by \code{i}.
    }
    \item{}{
      \code{x[[i]]}:
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
  \code{\link[Biostrings]{letterFrequencyInSlidingView}},
  \code{\link[Biostrings]{letterFrequency}},
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
  }
}

\examples{
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
ex <- exons(txdb, columns=c("tx_name", "gene_id"))
v <- Views(genome, ex)
v

subject(v)
granges(v)
v[[2]]
as(v, "DNAStringSet")
alphabetFrequency(v, collapse=TRUE)
}

\keyword{methods}
\keyword{classes}
