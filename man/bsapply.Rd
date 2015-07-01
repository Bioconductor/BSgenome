\name{bsapply}

\alias{bsapply}

\title{bsapply}

\description{
  Apply a function to each chromosome in a genome.
}

\usage{
bsapply(BSParams, ...)
}

\arguments{
  \item{BSParams}{
    a BSParams object that holds the various parameters needed to
    configure the bsapply function
  }
  \item{...}{
    optional arguments to 'FUN'.
  }
}

\details{
  By default the exclude parameter is set to not exclude anything.  A
  popular option will probably be to set this to "rand" so that random
  bits of unassigned contigs are filtered out.
}

\value{
  If \code{BSParams} sets \code{simplify=FALSE}, an ordinary list is returned
  containing the results generated using the remaining BSParams specifications.
  If \code{BSParams} sets \code{simplify=TRUE}, an \code{sapply}-like
  simplification is performed on the results.
}

\author{Marc Carlson}

\seealso{
  \link{BSParams-class}, \link{BSgenome-class}, \link{BSgenome-utils}
}

\examples{
  ## Load the Worm genome:
  library("BSgenome.Celegans.UCSC.ce2")

  ## Count the alphabet frequencies for every chromosome but exclude
  ## mitochrondrial and scaffold ones:
  params <- new("BSParams", X = Celegans, FUN = alphabetFrequency,
  exclude = c("M", "_"))
  bsapply(params)

  ## Or we can do this same function with simplify = TRUE:
  params <- new("BSParams", X = Celegans, FUN = alphabetFrequency,
  exclude = c("M", "_"), simplify = TRUE)
  bsapply(params)


  ## Examples to show how we might look for a string (in this case an
  ## ebox motif) across the whole genome.  
  Ebox <- DNAStringSet("CACGTG")
  pdict0 <- PDict(Ebox)

  params <- new("BSParams", X = Celegans, FUN = countPDict, simplify = TRUE)
  bsapply(params, pdict = pdict0)

  params@FUN <- matchPDict
  bsapply(params, pdict = pdict0)

  ## And since its really overkill to use matchPDict to find a single pattern:
  params@FUN <- matchPattern
  bsapply(params, pattern = "CACGTG")


  ## Examples on how to use the masks
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
  genome <- BSgenome.Hsapiens.UCSC.hg38.masked
  ## I can make things verbose if I want to see the chromosomes getting processed.
  options(verbose=TRUE)
  ## For the 1st example, lets use default masks
  params <- new("BSParams", X = genome, FUN = alphabetFrequency,
  exclude = c(1:8,"M","X","_"), simplify = TRUE)
  bsapply(params)

  if (interactive()) {
    ## Set up the motifList to filter out all double T's and all double C's
    params@motifList <-c("TT","CC")
    bsapply(params)

    ## Get rid of the motifList
    params@motifList=as.character()
  }

  ##Enable all standard masks
  params@maskList <- c(RM=TRUE,TRF=TRUE)
  bsapply(params)

  ##Disable all standard masks
  params@maskList <- c(AGAPS=FALSE,AMB=FALSE)
  bsapply(params)
}


\keyword{manip}
