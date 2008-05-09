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
\alias{SNPlocs_pkgname}
\alias{SNPlocs_pkgname,BSgenome-method}
\alias{seqnames}
\alias{seqnames,BSgenome-method}
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

% Other functions and generics:
\alias{unload}
\alias{unload,BSgenome-method}


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
      \code{SNPlocs_pkgname(x)}:
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
    }
    \item{}{
      \code{x$name}:
      Same as \code{x[[name]]} but \code{name} is not evaluated and
      therefore must be a literal character string or a name (possibly
      backtick quoted).
    }
  }
}

\section{Other functions and generics}{
  In the code snippets below,
  \code{x} is a BSgenome object
  and \code{name} is the name of a sequence (character-string).

  \describe{
    \item{}{
      \code{unload(x, name)}:
      Try to free the memory occupied by a loaded sequence by
      removing the 1st reference to this sequence. This 1st
      reference is a hidden reference that is created behind the scene
      by \code{x[[name]]} or \code{x$name}.
      See below for an example of how to make proper use of \code{unload()}.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{available.genomes}},
  \link[Biostrings]{XString-class},
  \link[Biostrings]{MaskedXString-class},
  \link[Biostrings]{XStringSet-class},
  \code{\link{getSeq}},
  \code{\link[Biostrings]{matchPattern}},
  \code{\link{rm}},
  \code{\link{gc}}
}

\examples{
  library(BSgenome.Celegans.UCSC.ce2)   # This doesn't load the chromosome 
                                        # sequences into memory.
  length(Celegans)                      # Number of sequences in this genome.
  Celegans                              # Displays a summary of the sequences
                                        # provided in this genome.
  seqnames(Celegans)                    # Index of single sequences.
  class(Celegans$chrI)                  # A DNAString instance.
  mseqnames(Celegans)                   # Index of multiple sequences.
  class(Celegans$upstream1000)          # A DNAStringSet instance.
  desc(Celegans$upstream1000)[1:4]      # Character vector containing the
                                        # description line found in the FASTA
                                        # file for the first 4 FASTA records.

  ## Some important considerations about memory usage:
  mem0 <- gc()["Vcells", "(Mb)"]        # Current amount of data in memory (in
                                        # Mb).
  Celegans[["chrV"]]                    # Loads chromosome V into memory (hence
                                        # takes a long time).
  gc()["Vcells", "(Mb)"] - mem0         # Chromosome V occupies 20Mb of memory.
  Celegans[["chrV"]]                    # Much faster (sequence is already in
                                        # memory, hence it's not loaded again).
  Celegans$chrV                         # Equivalent to Celegans[["chrV"]].
  class(Celegans$chrV)                  # Chromosome V (like any other
                                        # chromosome sequence) is a DNAString
                                        # object.
  nchar(Celegans$chrV)                  # It has 20922231 letters (nucleotides).
  x <- Celegans$chrV                    # Very fast because a BString object
                                        # doesn't contain the sequence, only a
                                        # pointer to the sequence, hence chrV
                                        # seq is not duplicated in memory. But
                                        # we now have 2 objects pointing to the
                                        # same place in memory.
  y <- substr(x, 10, 100)               # A 3rd object pointing to chrV seq.
  
  ## We must remove all references to chrV seq if we want the 20Mb of memory
  ## used by it to be freed (note that it can be hard to keep track of all the
  ## references to a given sequence).
  ## IMPORTANT: The 1st reference to this seq (Celegans$chrV) should be removed
  ## last. This is achieved with unload(). All other references are removed by
  ## just removing the referencing object.
  rm(x)
  rm(y)
  unload(Celegans, "chrV")
  gc()["Vcells", "(Mb)"]
}

\keyword{methods}
\keyword{classes}

