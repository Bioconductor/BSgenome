\name{BSgenome-class}
\docType{class}

% Classes
\alias{class:BSgenome}
\alias{BSgenome-class}
\alias{BSgenome}

% Accesor methods:
\alias{organism}
\alias{organism,BSgenome-method}
\alias{species}
\alias{species,BSgenome-method}
\alias{provider}
\alias{provider,BSgenome-method}
\alias{providerVersion}
\alias{providerVersion,BSgenome-method}
\alias{releaseDate}
\alias{releaseDate,BSgenome-method}
\alias{releaseName}
\alias{releaseName,BSgenome-method}
\alias{sourceUrl}
\alias{sourceUrl,BSgenome-method}
\alias{SNPlocs_pkg}
\alias{SNPlocs_pkg,BSgenome-method}
\alias{seqnames}
\alias{seqnames,BSgenome-method}
\alias{seqlengths}
\alias{seqlengths,BSgenome-method}
\alias{mseqnames}
\alias{mseqnames,BSgenome-method}
\alias{names,BSgenome-method}

% Constructor-like functions and generics:
\alias{BSgenome}

% Standard generic methods:
\alias{show,BSgenome-method}
\alias{length,BSgenome-method}
\alias{[[,BSgenome-method}
\alias{[[<-,BSgenome-method}
\alias{$,BSgenome-method}


\title{The BSgenome class}

\description{
  A container for the complete genome sequence of a given species.
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a BSgenome object
  and \code{name} is the name of a sequence (character-string).

  \describe{
    \item{}{
      \code{organism(x)}:
      Return the target organism for this genome e.g. \code{"Homo sapiens"},
      \code{"Mus musculus"}, \code{"Caenorhabditis elegans"}, etc...
    }
    \item{}{
      \code{species(x)}:
      Return the target species for this genome e.g. \code{"Human"},
      \code{"Mouse"}, \code{"C. elegans"}, etc...
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
      \code{sourceUrl(x)}:
      Return the source URL i.e. the permanent URL to the place where the
      FASTA files used to produce the sequences contained in \code{x} can
      be found (and downloaded).
    }
    \item{}{
      \code{SNPlocs_pkg(x)}:
      Return the name of the package from which the SNPs are injected,
      if any.
    }
    \item{}{
      \code{seqnames(x)}:
      Return the index of the single sequences contained in \code{x}.
      Each single sequence is stored in an \link[Biostrings]{XString}
      or \link[Biostrings]{MaskedXString} object and typically comes
      from a source file (FASTA) with a single record.
      The names returned by \code{seqnames(x)} usually reflect the names
      of those source files but a common prefix or suffix was eventually
      removed in order to keep them as short as possible.
    }
    \item{}{
      \code{seqlengths(x)}:
      Return the lengths of the single sequences contained in \code{x}.
      \code{names(seqlengths(x))} is guaranteed to be identical to
      \code{seqnames(x)}.
    }
    \item{}{
      \code{mseqnames(x)}:
      Return the index of the multiple sequences contained in \code{x}.
      Each multiple sequence is stored in an \link[Biostrings]{XStringSet}
      object and typically comes from a source file (FASTA) with multiple
      records.
      The names returned by \code{mseqnames(x)} usually reflect the names
      of those source files but a common prefix or suffix was eventually
      removed in order to keep them as short as possible.
    }
    \item{}{
      \code{names(x)}:
      Return the index of all sequences contained in \code{x}.
      This is the same as \code{c(seqnames(x), mseqnames(x))}.
    }
    \item{}{
      \code{length(x)}:
      Return the length of \code{x}, i.e., the number of all sequences
      that it contains. This is the same as \code{length(names(x))}.
    }
    \item{}{
      \code{x[[name]]}:
      Return sequence (single or multiple) named \code{name}.
      No sequence is actually loaded into memory until this is explicitely
      requested with a call to \code{x[[name]]} or \code{x$name}.
      When loaded, a sequence is kept in a cache. It will be automatically
      removed from the cache at garbage collection if it's not in use anymore
      i.e. if there are no reference to it (other than the reference stored
      in the cache). With \code{options(verbose=TRUE)}, a message is printed
      each time a sequence is removed from the cache. 
    }
    \item{}{
      \code{x$name}:
      Same as \code{x[[name]]} but \code{name} is not evaluated and
      therefore must be a literal character string or a name (possibly
      backtick quoted).
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{available.genomes}},
  \link[Biostrings]{XString-class},
  \link[Biostrings]{MaskedXString-class},
  \link[Biostrings]{XStringSet-class},
  \code{\link[IRanges]{subseq}},
  \code{\link{getSeq}},
  \code{\link[Biostrings]{matchPattern}},
  \code{\link{rm}},
  \code{\link{gc}}
}

\examples{
  ## Loading a BSgenome data package doesn't load its sequences
  ## into memory:
  library(BSgenome.Celegans.UCSC.ce2)

  ## Number of sequences in this genome:
  length(Celegans) 

  ## Display a summary of the sequences:
  Celegans

  ## Index of single sequences:
  seqnames(Celegans)

  ## Lengths (i.e. number of nucleotides) of the sequences:
  seqlengths(Celegans)

  ## Load chromosome I from disk to memory (hence takes some time)
  ## and keep a reference to it:
  chrI <- Celegans[["chrI"]]  # equivalent to Celegans$chrI

  chrI

  class(chrI)   # a DNAString instance
  length(chrI)  # with 15080483 nucleotides

  ## Multiple sequences:
  mseqnames(Celegans) 
  upstream1000 <- Celegans$upstream1000
  upstream1000
  class(upstream1000)  # a DNAStringSet instance
  ## Character vector containing the description lines of the first
  ## 4 sequences in the original FASTA file:
  names(upstream1000)[1:4]

  ## ---------------------------------------------------------------------
  ## PASS-BY-ADDRESS SEMANTIC, CACHING AND MEMORY USAGE
  ## ---------------------------------------------------------------------

  ## We want a message to be printed each time a sequence is removed
  ## from the cache:
  options(verbose=TRUE)

  gc()  # nothing seems to be removed from the cache
  rm(chrI, upstream1000)
  gc()  # chrI and upstream1000 are removed from the cache (they are
        # not in use anymore)

  options(verbose=FALSE)

  ## Get the current amount of data in memory (in Mb):
  mem0 <- gc()["Vcells", "(Mb)"]

  system.time(chrV <- Celegans[["chrV"]])  # read from disk
  
  gc()["Vcells", "(Mb)"] - mem0  # chrV occupies 20Mb in memory

  system.time(tmp <- Celegans[["chrV"]])  # much faster! (sequence
                                          # is in the cache)

  gc()["Vcells", "(Mb)"] - mem0  # we're still using 20Mb (sequences
                                 # have a pass-by-address semantic
                                 # i.e. the sequence data are not
                                 # duplicated)
  
  ## subseq() doesn't copy the sequence data either, hence it is very
  ## fast and memory efficient (but the returned object will hold a
  ## reference to chrV):
  y <- subseq(chrV, 10, 8000000) 
  gc()["Vcells", "(Mb)"] - mem0

  ## We must remove all references to chrV before it can be removed from
  ## the cache (so the 20Mb of memory used by this sequence are freed).
  options(verbose=TRUE)
  rm(chrV, tmp)
  gc()

  ## Remember that 'y' holds a reference to chrV too:
  rm(y)
  gc()

  options(verbose=FALSE)
  gc()["Vcells", "(Mb)"] - mem0
}

\keyword{methods}
\keyword{classes}
