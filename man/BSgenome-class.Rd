\name{BSgenome-class}
\docType{class}

% Classes:
\alias{class:BSgenome}
\alias{BSgenome-class}
\alias{BSgenome}

% Accessors:
\alias{length,BSgenome-method}
\alias{sourceUrl}
\alias{sourceUrl,BSgenome-method}
\alias{mseqnames}
\alias{mseqnames,BSgenome-method}
\alias{names,BSgenome-method}
\alias{masknames}
\alias{masknames,BSgenome-method}
\alias{seqinfo,BSgenome-method}
\alias{seqinfo<-,BSgenome-method}
\alias{seqnames<-,BSgenome-method}

% Constructor-like functions and generics:
\alias{BSgenome}

% Other methods:
\alias{as.list,BSgenome-method}
\alias{show,BSgenome-method}
\alias{[[,BSgenome-method}
\alias{$,BSgenome-method}


\title{BSgenome objects}

\description{
  The BSgenome class is a container for storing the full genome sequences
  of a given organism. BSgenome objects are usually made in advance by
  a volunteer and made available to the Bioconductor community as
  "BSgenome data packages".
  See \code{?\link{available.genomes}} for how to get the list of
  "BSgenome data packages" curently available.
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is a BSgenome object.
  Note that, because the BSgenome class contains the
  \link[GenomeInfoDb]{GenomeDescription} class, then all the accessor methods
  for \link[GenomeInfoDb]{GenomeDescription} objects can also be used on
  \code{x}.

  \describe{
    \item{\code{sourceUrl(x)}}{
      Returns the source URL i.e. the permanent URL to the place where the
      FASTA files used to produce the sequences contained in \code{x} can
      be found (and downloaded).
    }
    \item{\code{seqnames(x)}, \code{seqnames(x) <- value}}{
      Gets or sets the names of the single sequences contained in \code{x}.
      Each single sequence is stored in a \link[Biostrings]{DNAString}
      or \link[Biostrings]{MaskedDNAString} object and typically comes
      from a source file (FASTA) with a single record.
      The names returned by \code{seqnames(x)} usually reflect the names
      of those source files but a common prefix or suffix was eventually
      removed in order to keep them as short as possible.
    }
    \item{\code{seqlengths(x)}}{
      Returns the lengths of the single sequences contained in \code{x}.

      See \code{?`\link[XVector]{length,XVector-method}`} and
      \code{?`\link[Biostrings]{length,MaskedXString-method}`} for
      the definition of the length of a \link[Biostrings]{DNAString}
      or \link[Biostrings]{MaskedDNAString} object.
      Note that the length of a masked sequence
      (\link[Biostrings]{MaskedXString} object) is not
      affected by the current set of active masks but the \code{nchar}
      method for \link[Biostrings]{MaskedXString} objects is.

      \code{names(seqlengths(x))} is guaranteed to be identical to
      \code{seqnames(x)}.
    }
    \item{\code{mseqnames(x)}}{
      Returns the index of the multiple sequences contained in \code{x}.
      Each multiple sequence is stored in a
      \link[Biostrings]{DNAStringSet} object and typically comes from
      a source file (FASTA) with multiple records.
      The names returned by \code{mseqnames(x)} usually reflect the names
      of those source files but a common prefix or suffix was eventually
      removed in order to keep them as short as possible.
    }
    \item{\code{names(x)}}{
      Returns the index of all sequences contained in \code{x}.
      This is the same as \code{c(seqnames(x), mseqnames(x))}.
    }
    \item{\code{length(x)}}{
      Returns the length of \code{x}, i.e., the total number of sequences
      in it (single and multiple sequences). This is the same as
      \code{length(names(x))}.
    }
    \item{\code{x[[name]]}}{
      Returns the sequence (single or multiple) in \code{x} named \code{name}
      (\code{name}  must be a single string).
      No sequence is actually loaded into memory until this is explicitely
      requested with a call to \code{x[[name]]} or \code{x$name}.
      When loaded, a sequence is kept in a cache. It will be automatically
      removed from the cache at garbage collection if it's not in use anymore
      i.e. if there are no reference to it (other than the reference stored
      in the cache). With \code{options(verbose=TRUE)}, a message is printed
      each time a sequence is removed from the cache. 
    }
    \item{\code{x$name}}{
      Same as \code{x[[name]]} but \code{name} is not evaluated and
      therefore must be a literal character string or a name (possibly
      backtick quoted).
    }
    \item{\code{masknames(x)}}{
      The names of the built-in masks that are defined for all the single
      sequences. There can be up to 4 built-in masks per sequence. These will
      always be (in this order):
      (1) the mask of assembly gaps, aka "the AGAPS mask";
      
      (2) the mask of intra-contig ambiguities, aka "the AMB mask";
      
      (3) the mask of repeat regions that were determined by the RepeatMasker
      software, aka "the RM mask";
      
      (4) the mask of repeat regions that were determined by the Tandem Repeats
          Finder software (where only repeats with period less than or equal to
          12 were kept), aka "the TRF mask".
	  
      All the single sequences in a given package are guaranteed to have the
      same collection of built-in masks (same number of masks and in the same
      order).

      \code{masknames(x)} gives the names of the masks in this collection.
      Therefore the value returned by \code{masknames(x)} is a character vector
      made of the first N elements of \code{c("AGAPS", "AMB", "RM", "TRF")},
      where N depends only on the BSgenome data package being looked at
      (0 <= N <= 4).
      The man page for most BSgenome data packages should provide the exact
      list and permanent URLs of the source data files that were used to
      extract the built-in masks.
      For example, if you've installed the BSgenome.Hsapiens.UCSC.hg38 package,
      load it and see the Note section in
      \code{?`\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}`}.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{available.genomes}},
  \link[GenomeInfoDb]{GenomeDescription-class},
  \link{BSgenome-utils},
  \link[Biostrings]{DNAString-class},
  \link[Biostrings]{DNAStringSet-class},
  \link[Biostrings]{MaskedDNAString-class},
  \code{\link{getSeq,BSgenome-method}},
  \code{\link{injectSNPs}},
  \link[XVector]{subseq,XVector-method},
  \code{\link[base]{rm}},
  \code{\link[base]{gc}}
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

## Lengths (i.e. number of nucleotides) of the single sequences:
seqlengths(Celegans)

## Load chromosome I from disk to memory (hence takes some time)
## and keep a reference to it:
chrI <- Celegans[["chrI"]]  # equivalent to Celegans$chrI

chrI

class(chrI)   # a DNAString instance
length(chrI)  # with 15080483 nucleotides

## Single sequence can be renamed:
seqnames(Celegans) <- sub("^chr", "", seqnames(Celegans))
seqlengths(Celegans)
Celegans$I
seqnames(Celegans) <- paste0("chr", seqnames(Celegans))

## Multiple sequences:
library(BSgenome.Rnorvegicus.UCSC.rn5)
rn5 <- BSgenome.Rnorvegicus.UCSC.rn5
rn5
seqnames(rn5)
rn5_chr1 <- rn5$chr1
mseqnames(rn5) 
rn5_random  <- Rnorvegicus$random
rn5_random
class(rn5_random)  # a DNAStringSet instance
## Character vector containing the description lines of the first
## 4 sequences in the original FASTA file:
names(rn5_random)[1:4]

## ---------------------------------------------------------------------
## PASS-BY-ADDRESS SEMANTIC, CACHING AND MEMORY USAGE
## ---------------------------------------------------------------------

## We want a message to be printed each time a sequence is removed
## from the cache:
options(verbose=TRUE)

gc()  # nothing seems to be removed from the cache
rm(rn5_chr1, rn5_random)
gc()  # rn5_chr1 and rn5_random are removed from the cache (they are
      # not in use anymore)

options(verbose=FALSE)

## Get the current amount of data in memory (in Mb):
mem0 <- gc()["Vcells", "(Mb)"]

system.time(rn5_chr2 <- rn5$chr2)  # read from disk
  
gc()["Vcells", "(Mb)"] - mem0  # 'rn5_chr2' occupies 20Mb in memory

system.time(tmp <- rn5$chr2)  # much faster! (sequence
                              # is in the cache)

gc()["Vcells", "(Mb)"] - mem0  # we're still using 20Mb (sequences
                               # have a pass-by-address semantic
                               # i.e. the sequence data are not
                               # duplicated)
  
## subseq() doesn't copy the sequence data either, hence it is very
## fast and memory efficient (but the returned object will hold a
## reference to 'rn5_chr2'):
y <- subseq(rn5_chr2, 10, 8000000) 
gc()["Vcells", "(Mb)"] - mem0

## We must remove all references to 'rn5_chr2' before it can be
## removed from the cache (so the 20Mb of memory used by this
## sequence are freed):
options(verbose=TRUE)
rm(rn5_chr2, tmp)
gc()

## Remember that 'y' holds a reference to 'rn5_chr2' too:
rm(y)
gc()

options(verbose=FALSE)
gc()["Vcells", "(Mb)"] - mem0
}

\keyword{methods}
\keyword{classes}
