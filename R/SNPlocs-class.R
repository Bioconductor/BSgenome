### =========================================================================
### SNPlocs objects
### -------------------------------------------------------------------------

setClass("SNPlocs",
    representation(
        ## provider of the SNPs (e.g. "dbSNP")
        provider="character",

        ## creation date (in compact format) of the flat files found at
        ## download_url (look at 2nd line of each file e.g.
        ## CREATED ON: 2012-06-08 14:53, and use the most recent date in
        ## case of mixed dates, e.g. "20120608")
        provider_version="character",

        ## official release date of the SNPs (e.g. "Nov 9, 2010")
        release_date="character",

        ## official release name of the SNPs (e.g. "Build 132")
        release_name="character",

        ## URL to the place where the original SNP data was downloaded
        download_url="character",

        ## date the original SNP data was downloaded
        download_date="character",

        ## reference genome of the SNPs
        reference_genome="GenomeDescription",

        ## named list of "sequence name translation tables" (one table per
        ## compatible genome and each table is represented by a named character
        ## vector)
        compatible_genomes="list",

        ## package name and absolute path to local directory where to find
        ## the serialized objects containing the SNPs
        data_pkgname="character",
        data_dirpath="character",

        .data_cache="environment"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Not intended to be used directly.
newSNPlocs <- function(provider, provider_version,
                       release_date, release_name,
                       download_url, download_date,
                       reference_genome, compatible_genomes,
                       data_pkgname, data_dirpath)
{
    new("SNPlocs",
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        download_url=download_url,
        download_date=download_date,
        reference_genome=reference_genome,
        compatible_genomes=compatible_genomes,
        data_pkgname=data_pkgname,
        data_dirpath=data_dirpath,
        .data_cache=new.env(parent=emptyenv()))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

