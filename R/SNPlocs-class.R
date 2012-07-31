### =========================================================================
### SNPlocs objects
### -------------------------------------------------------------------------

setClass("SNPlocs",
    representation(
        ## provider of the SNPs (e.g. "dbSNP")
        provider="character",

        ## official release date of the SNPs (e.g. "Nov 9, 2010")
        release_date="character",

        ## provider version of the SNPs (e.g. "Build 132")
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
newSNPlocs <- function(provider, release_date, release_name,
                       download_url, download_date,
                       reference_genome, compatible_genomes,
                       data_pkgname, data_dirpath)
{
    new("SNPlocs",
        provider=provider,
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

