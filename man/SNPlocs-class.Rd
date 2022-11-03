\name{SNPlocs-class}
\docType{class}

\alias{class:SNPlocs}
\alias{SNPlocs-class}
\alias{SNPlocs}

\alias{class:ODLT_SNPlocs}
\alias{ODLT_SNPlocs-class}
\alias{ODLT_SNPlocs}

\alias{class:OldFashionSNPlocs}
\alias{OldFashionSNPlocs-class}
\alias{OldFashionSNPlocs}

% accessors
\alias{provider,SNPlocs-method}
\alias{providerVersion,SNPlocs-method}
\alias{releaseDate,SNPlocs-method}
\alias{releaseName}
\alias{releaseName,SNPlocs-method}
\alias{referenceGenome}
\alias{referenceGenome,SNPlocs-method}
\alias{compatibleGenomes}
\alias{compatibleGenomes,SNPlocs-method}
\alias{organism,SNPlocs-method}
\alias{commonName,SNPlocs-method}
\alias{seqinfo,SNPlocs-method}
\alias{seqnames,SNPlocs-method}

% constructors
\alias{new_ODLT_SNPlocs}
\alias{newSNPlocs}

% displaying
\alias{show,SNPlocs-method}

% SNP extractors
\alias{snpcount}
\alias{snpcount,SNPlocs-method}
\alias{snpcount,ODLT_SNPlocs-method}
\alias{snpcount,OldFashionSNPlocs-method}
\alias{snplocs}
\alias{snplocs,SNPlocs-method}
\alias{snplocs,ODLT_SNPlocs-method}
\alias{snplocs,OldFashionSNPlocs-method}
\alias{snpsBySeqname}
\alias{snpsBySeqname,SNPlocs-method}
\alias{snpsBySeqname,ODLT_SNPlocs-method}
\alias{snpsBySeqname,OldFashionSNPlocs-method}
\alias{snpsByOverlaps}
\alias{snpsByOverlaps,SNPlocs-method}
\alias{snpsByOverlaps,ODLT_SNPlocs-method}
\alias{snpsByOverlaps,OldFashionSNPlocs-method}
\alias{snpsById}
\alias{snpsById,SNPlocs-method}
\alias{snpsById,ODLT_SNPlocs-method}
\alias{snpsById,OldFashionSNPlocs-method}

\alias{inferRefAndAltAlleles}

\title{SNPlocs objects}

\description{
  The SNPlocs class is a container for storing known SNP locations (of class
  \emph{snp}) for a given organism.

  SNPlocs objects are usually made in advance by a volunteer and made
  available to the Bioconductor community as \emph{SNPlocs data packages}.
  See \code{?\link{available.SNPs}} for how to get the list of
  \emph{SNPlocs and \link{XtraSNPlocs} data packages} curently available.

  The main focus of this man page is on how to extract SNPs from an
  SNPlocs object.
}

\usage{
snpcount(x)

snpsBySeqname(x, seqnames, ...)
\S4method{snpsBySeqname}{SNPlocs}(x, seqnames, drop.rs.prefix=FALSE, genome=NULL)

snpsByOverlaps(x, ranges, ...)
\S4method{snpsByOverlaps}{SNPlocs}(x, ranges, drop.rs.prefix=FALSE, ..., genome=NULL)

snpsById(x, ids, ...)
\S4method{snpsById}{SNPlocs}(x, ids, ifnotfound=c("error", "warning", "drop"), genome=NULL)

inferRefAndAltAlleles(gpos, genome)
}

\arguments{
  \item{x}{
    A SNPlocs object.
  }
  \item{seqnames}{
    The names of the sequences for which to get SNPs. Must be a subset of
    \code{seqlevels(x)}. NAs and duplicates are not allowed.
  }
  \item{...}{
    Additional arguments, for use in specific methods.

    Arguments passed to the \code{snpsByOverlaps} method for SNPlocs
    objects thru \code{...} are used internally in the call to
    \code{\link[IRanges]{subsetByOverlaps}()}.
    See \code{?IRanges::\link[IRanges]{subsetByOverlaps}} in the \pkg{IRanges}
    package and \code{?GenomicRanges::\link[GenomicRanges]{subsetByOverlaps}}
    in the \pkg{GenomicRanges} package for more information about the
    \code{subsetByOverlaps()} generic and its method for
    \link[GenomicRanges]{GenomicRanges} objects.
  }
  \item{drop.rs.prefix}{
    Should the \code{rs} prefix be dropped from the returned RefSNP ids?
    (RefSNP ids are stored in the \code{RefSNP_id} metadata column of the
    returned object.)
  }
  \item{genome}{
    For \code{snpsBySeqname}, \code{snpsByOverlaps}, and \code{snpsById}:

    \code{NULL} (the default), or a \link{BSgenome} object
    containing the sequences of the reference genome that corresponds
    to the SNP positions. See \code{inferRefAndAltAlleles} below for
    an alternative way to specify \code{genome}.

    If \code{genome} is supplied, then \code{inferRefAndAltAlleles} is
    called internally by \code{snpsBySeqname}, \code{snpsByOverlaps},
    or \code{snpsById} to \emph{infer} the reference allele (a.k.a.
    \emph{ref} allele) and alternate allele(s) (a.k.a. \emph{alt}
    allele(s)) for each SNP in the returned \link[GenomicRanges]{GPos}
    object. The inferred \emph{ref} allele and \emph{alt} allele(s) are
    returned in additional metadata columns \code{ref_allele} (character)
    and \code{alt_alleles} (CharacterList).

    For \code{inferRefAndAltAlleles}:

    A \link{BSgenome} object containing the sequences of the
    reference genome that corresponds to the SNP positions in
    \code{gpos}. Alternatively \code{genome} can be a single
    string containing the name of the reference genome, in which
    case it must be specified in a way that is accepted by the
    \code{\link{getBSgenome}} function (e.g. \code{"GRCh38"})
    and the corresponding BSgenome data package needs to be already
    installed (see \code{?\link{getBSgenome}} for the details).
  }
  \item{ranges}{
    One or more genomic regions of interest specified as a
    \link[GenomicRanges]{GRanges} or \link[GenomicRanges]{GPos} object.
    A single region of interest can be specified as a character string of
    the form \code{"ch14:5201-5300"}.
  }
  \item{ids}{
    The RefSNP ids to look up (a.k.a. rs ids). Can be integer or character
    vector, with or without the \code{"rs"} prefix. NAs are not allowed.
  }
  \item{ifnotfound}{
    What to do if SNP ids are not found.
  }
  \item{gpos}{
    A \link[GenomicRanges]{GPos} object containing SNPs. It must have a
    metadata column \code{alleles_as_ambig} like obtained when using
    any of the SNP extractor \code{snpsBySeqname}, \code{snpsByOverlaps},
    or \code{snpsById} on a SNPlocs object.
  }
}

\details{
  When the reference genome is specified via the \code{genome} argument,
  SNP extractors \code{snpsBySeqname}, \code{snpsByOverlaps}, and
  \code{snpsById} call \code{inferRefAndAltAlleles} internally to
  \emph{infer} the reference allele (a.k.a. \emph{ref} allele) and
  alternate allele(s) (a.k.a. \emph{alt} allele(s)) for each SNP.

  For each SNP the \emph{ref} allele is inferred from the actual
  nucleotide found in the reference genome at the SNP position.
  The \emph{alt} alleles are inferred from metadata column
  \code{alleles_as_ambig} and the \code{ref} allele. More precisely
  for each SNP the \emph{alt} alleles are considered to be the alleles
  in \code{alleles_as_ambig} minus the \emph{ref} allele.
}

\value{
  \code{snpcount} returns a named integer vector containing the number
  of SNPs for each sequence in the reference genome.

  \code{snpsBySeqname}, \code{snpsByOverlaps}, and \code{snpsById} return
  an \emph{unstranded} \link[GenomicRanges]{GPos} object with one element
  (genomic position) per SNP and the following metadata columns:
  \itemize{
    \item \code{RefSNP_id}: RefSNP ID (aka "rs id"). Character vector
          with no NAs and no duplicates.
    \item \code{alleles_as_ambig}: A character vector with no NAs
          containing the alleles for each SNP represented by an IUPAC
          nucleotide ambiguity code.
          See \code{?\link[Biostrings]{IUPAC_CODE_MAP}} in the
          \pkg{Biostrings} package for more information.
  }
  If the reference genome was specified (via the \code{genome} argument),
  the additional metadata columns are returned:
  \itemize{
    \item \code{genome_compat}: A logical vector indicating whether the
          alleles in \code{alleles_as_ambig} are consistent with the
          reference genome.
    \item \code{ref_allele}: A character vector containing the
          \emph{inferred} reference allele for each SNP.
    \item \code{alt_alleles}: A \link[IRanges]{CharacterList} object
          where each list element is a character vector containing the
          \emph{inferred} alternate allele(s) for the corresponding SNP.
  }

  Note that this \link[GenomicRanges]{GPos} object is \emph{unstranded}
  i.e. all the SNPs in it have their strand set to \code{"*"}.
  Alleles are always reported with respect to the \emph{positive} strand.

  If \code{ifnotfound="error"}, the object returned by \code{snpsById}
  is guaranteed to be \emph{parallel} to \code{ids}, that is, the i-th
  element in the \link[GenomicRanges]{GPos} object corresponds to the
  i-th element in \code{ids}.

  \code{inferRefAndAltAlleles} returns a \link[S4Vectors]{DataFrame} with
  one row per SNP in \code{gpos} and with columns \code{genome_compat}
  (logical), \code{ref_allele} (character), and \code{alt_alleles}
  (CharacterList).

}

\author{H. Pagès}

\seealso{
  \itemize{
    \item \code{\link{available.SNPs}}

    \item \link[GenomicRanges]{GPos} and \link[GenomicRanges]{GRanges}
          objects in the \pkg{GenomicRanges} package.

    \item \link{XtraSNPlocs} packages and objects for molecular variations
          of class other than \emph{snp} e.g. of class \emph{in-del},
          \emph{heterozygous}, \emph{microsatellite}, etc...

    \item \code{IRanges::\link[IRanges]{subsetByOverlaps}} in the
          \pkg{IRanges} package and
          \code{GenomicRanges::\link[GenomicRanges]{subsetByOverlaps}}
          in the \pkg{GenomicRanges} package for more information about the
          \code{subsetByOverlaps()} generic and its method for
          \link[GenomicRanges]{GenomicRanges} objects.

    \item \code{\link{injectSNPs}}

    \item \code{\link[Biostrings]{IUPAC_CODE_MAP}} in the \pkg{Biostrings}
          package.
  }
}

\examples{
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh38
snpcount(snps)

## ---------------------------------------------------------------------
## snpsBySeqname()
## ---------------------------------------------------------------------

## Get all SNPs located on chromosome 22 or MT:
snpsBySeqname(snps, c("22", "MT"))

## ---------------------------------------------------------------------
## snpsByOverlaps()
## ---------------------------------------------------------------------

## Get all SNPs overlapping some genomic region of interest:
snpsByOverlaps(snps, "X:3e6-33e6")

## With the regions of interest being all the known CDS for hg38
## located on chromosome 22 or MT (except for the chromosome naming
## convention, hg38 is the same as GRCh38):
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
my_cds <- cds(txdb)
seqlevels(my_cds, pruning.mode="coarse") <- c("chr22", "chrM")
seqlevelsStyle(my_cds)  # UCSC
seqlevelsStyle(snps)    # NCBI
seqlevelsStyle(my_cds) <- seqlevelsStyle(snps)
genome(my_cds) <- genome(snps)
my_snps <- snpsByOverlaps(snps, my_cds)
my_snps
table(my_snps \%within\% my_cds)

## ---------------------------------------------------------------------
## snpsById()
## ---------------------------------------------------------------------

## Lookup some RefSNP ids:
my_rsids <- c("rs10458597", "rs12565286", "rs7553394")
\dontrun{
  snpsById(snps, my_rsids)  # error, rs7553394 not found
}
## The following example uses more than 2GB of memory, which is more
## than what 32-bit Windows can handle:
is_32bit_windows <- .Platform$OS.type == "windows" &&
                    .Platform$r_arch == "i386"
if (!is_32bit_windows) {
    snpsById(snps, my_rsids, ifnotfound="drop")
}

## ---------------------------------------------------------------------
## Obtaining the ref allele and alt allele(s)
## ---------------------------------------------------------------------

## When the reference genome is specified (via the 'genome' argument),
## SNP extractors snpsBySeqname(), snpsByOverlaps(), and snpsById()
## call inferRefAndAltAlleles() internally to **infer** the ref allele
## and alt allele(s) for each SNP.
my_snps <- snpsByOverlaps(snps, "X:3e6-8e6", genome="GRCh38")
my_snps

## Most SNPs have only 1 alternate allele:
table(lengths(mcols(my_snps)$alt_alleles))

## SNPs with 2 alternate alleles:
my_snps[lengths(mcols(my_snps)$alt_alleles) == 2]

## SNPs with 3 alternate alleles:
my_snps[lengths(mcols(my_snps)$alt_alleles) == 3]

## Note that a small percentage of SNPs in dbSNP have alleles that
## are inconsistent with the reference genome (don't ask me why):
table(mcols(my_snps)$genome_compat)

## For the inconsistent SNPs, all the alleles reported by dbSNP
## are considered alternate alleles i.e. for each inconsistent SNP
## metadata columns "alleles_as_ambig" and "alt_alleles" represent
## the same set of nucleotides (the latter being just an expanded
## representation of the IUPAC ambiguity letter in the former):
my_snps[!mcols(my_snps)$genome_compat]
}

\keyword{methods}
\keyword{classes}
