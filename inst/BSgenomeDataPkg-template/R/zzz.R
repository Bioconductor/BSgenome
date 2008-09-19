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
    seqs_pkg="@PKGNAME@",
    seqs_dir="extdata",
    nmask_per_seq=@NMASKPERSEQ@,
    masks_pkg="@PKGNAME@",
    masks_dir="data"
)

