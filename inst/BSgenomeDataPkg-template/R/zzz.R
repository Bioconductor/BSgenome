###
###
###

@BSGENOMEOBJNAME@ <- BSgenome(
    organism="@ORGANISM@",
    species="@SPECIES@",
    provider="@PROVIDER@",
    provider_version="@PROVIDERVERSION@",
    release_date="@RELEASEDATE@",
    release_name="@RELEASENAME@",
    source_url="@SOURCEURL@",
    seqnames=@SEQNAMES@,
    mseqnames=@MSEQNAMES@,
    package="@PKGNAME@",
    subdir="extdata"
)

