###
###
###

.seqnames <- @SEQNAMES@
.mseqnames <- @MSEQNAMES@

@BSGENOMEOBJNAME@ <- BSgenome(
    organism="@ORGANISM@",
    species="@SPECIES@",
    provider="@PROVIDER@",
    provider_version="@PROVIDERVERSION@",
    release_date="@RELEASEDATE@",
    release_name="@RELEASENAME@",
    source_url="@SOURCEURL@",
    seqnames=.seqnames,
    mseqnames=.mseqnames,
    seqs_pkg="@PKGNAME@",
    seqlengths_dir="data",
    seqs_dir="extdata",
    nmask_per_seq=@NMASKPERSEQ@,
    masks_pkg="@PKGNAME@",
    masks_dir="data"
)

