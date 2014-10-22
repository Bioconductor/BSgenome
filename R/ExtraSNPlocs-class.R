### =========================================================================
### ExtraSNPlocs objects
### -------------------------------------------------------------------------
###

setClass("ExtraSNPlocs",
    representation(
        ## Name of the ExtraSNPlocs data package where the ExtraSNPlocs
        ## object is defined.
        pkgname="character",

        ## OnDiskLongTable object containing the SNP table.
        snp_table="OnDiskLongTable",

        ## Provider of the SNPs (e.g. "dbSNP").
        provider="character",

        ## E.g. "dbSNP Human BUILD 141".
        provider_version="character",

        ## Official release date of the SNPs (e.g. "May 2014").
        release_date="character",

        ## Official release name of the SNPs (e.g. "dbSNP Human BUILD 141").
        release_name="character",

        ## URL to the place where the original SNP data was downloaded.
        download_url="character",

        ## Date the original SNP data was downloaded.
        download_date="character",

        ## Reference genome of the SNPs.
        reference_genome="GenomeDescription"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("provider", "ExtraSNPlocs", function(x) x@provider)

setMethod("providerVersion", "ExtraSNPlocs", function(x) x@provider_version)

setMethod("releaseDate", "ExtraSNPlocs", function(x) x@release_date)

setMethod("releaseName", "ExtraSNPlocs", function(x) x@release_name)

setMethod("referenceGenome", "ExtraSNPlocs", function(x) x@reference_genome)

setMethod("organism", "ExtraSNPlocs", function(x) organism(referenceGenome(x)))

setMethod("species", "ExtraSNPlocs", function(x) species(referenceGenome(x)))

setMethod("seqinfo", "ExtraSNPlocs", function(x) seqinfo(referenceGenome(x)))

setMethod("seqnames", "ExtraSNPlocs", function(x) seqnames(referenceGenome(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Not intended to be used directly.
newExtraSNPlocs <- function(pkgname, snp_table_dirpath,
                            provider, provider_version,
                            release_date, release_name,
                            download_url, download_date,
                            reference_genome)
{
    snp_table <- OnDiskLongTable(snp_table_dirpath)
    new("ExtraSNPlocs",
        pkgname=pkgname,
        snp_table=snp_table,
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        download_url=download_url,
        download_date=download_date,
        reference_genome=reference_genome)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

setMethod("show", "ExtraSNPlocs",
    function(object)
    {
        cat(class(object), " object for ", organism(object), " (",
            provider(object), ": ", releaseName(object), ")\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpcount() and snplocs()
###

setMethod("snpcount", "ExtraSNPlocs",
    function(x)
    {
        stop("NOT READY YET!")        
    }
)

### Returns a data frame (when 'as.GRanges=FALSE') or a GRanges object
### (when 'as.GRanges=TRUE').
setMethod("snplocs", "ExtraSNPlocs",
    function(x, seqname, as.GRanges=FALSE, caching=TRUE)
    {
        stop("NOT READY YET!")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpid2loc(), snpid2alleles(), and snpid2grange()
###

setMethod("snpid2loc", "ExtraSNPlocs",
    function(x, snpid, caching=TRUE)
    {
        stop("NOT READY YET!")
    }
)

setMethod("snpid2alleles", "ExtraSNPlocs",
    function(x, snpid, caching=TRUE)
    {
        stop("NOT READY YET!")
    }
)

setMethod("snpid2grange", "ExtraSNPlocs",
    function(x, snpid, caching=TRUE)
    {
        stop("NOT READY YET!")
    }
)

