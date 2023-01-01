### =========================================================================
### ODLT_SNPlocs objects
### -------------------------------------------------------------------------


setClass("ODLT_SNPlocs",
    contains="SNPlocs",
    representation(
        snp_table="OnDiskLongTable"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor
###

new_ODLT_SNPlocs <- function(provider, provider_version,
                             release_date, release_name,
                             source_data_url, download_date,
                             reference_genome, compatible_genomes,
                             data_pkgname, data_dirpath)
{
    snp_table <- OnDiskLongTable(data_dirpath)
    new("ODLT_SNPlocs",
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        source_data_url=source_data_url,
        download_date=download_date,
        reference_genome=reference_genome,
        compatible_genomes=compatible_genomes,
        data_pkgname=data_pkgname,
        snp_table=snp_table)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpcount()
###

setMethod("snpcount", "ODLT_SNPlocs",
    function(x)
    {
        spatial_index <- spatialIndex(x@snp_table)
        batch_seqnames <- seqnames(spatial_index)
        batches_per_seq <- runLength(batch_seqnames)
        batch_breakpoints <- breakpoints(x@snp_table)
        seq_breakpoints <- batch_breakpoints[cumsum(batches_per_seq)]
        setNames(S4Vectors:::diffWithInitialZero(seq_breakpoints),
                 runValue(batch_seqnames))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snplocs()
###
### Used internally for SNP injection. Not intended for the end user.
### Must return a 2-col data-frame-like object with columns "loc" (integer)
### and "alleles_as_ambig" (character).
###

setMethod("snplocs", "ODLT_SNPlocs",
    function(x, seqname)
    {
        df <- getBatchesBySeqnameFromOnDiskLongTable(x@snp_table, seqname)
        data.frame(loc=df[ , "pos"],
                   alleles_as_ambig=decode_bytes_as_letters(df[ , "alleles"]),
                   stringsAsFactors=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### inferRefAndAltAlleles()
###
### Infers the ref allele and alt allele(s) for each SNP in 'gpos'.
### 'gpos' must be a GPos derivative containing SNPs. It must have a metadata
### column "alleles_as_ambig" like obtained when using any of the SNP
### extractors snpsBySeqname(), snpsByOverlaps(), or snpsById() on a SNPlocs
### object.
### For each SNP the ref allele is inferred from the nucleotide
### found in reference genome 'genome' at the SNP position.
### The alt alleles are inferred from metadata column "alleles_as_ambig" and
### the ref allele. More precisely for each SNP the alt alleles are considered
### to be the alleles in "alleles_as_ambig" minus the ref allele.
### Return a DataFrame with one row per SNP in 'gpos' and with columns
### "genome_compat" (logical), "ref_allele" (character), and "alt_alleles"
### (CharacterList).
inferRefAndAltAlleles <- function(gpos, genome)
{
    ## Check 'gpos'.
    if(!is(gpos, "GPos"))
        stop(wmsg("'gpos' must be a GPos derivative"))

    ## Check metadata column "alleles_as_ambig".
    alleles_as_ambig <- mcols(gpos)$alleles_as_ambig
    if (!is.character(alleles_as_ambig))
        stop(wmsg("'gpos' must have metadata column \"alleles_as_ambig\" ",
                  "and it must be a character vector"))
    alleles <- IUPAC_CODE_MAP[alleles_as_ambig]
    if (anyNA(alleles))
        stop(wmsg("invalid metadata column \"alleles_as_ambig\""))
    alleles <- DNAStringSet(alleles)

    ## Check 'genome'.
    genome <- getBSgenome(genome)

    ## Compute 'ref_allele'.
    ## getSeq() will complain if 'seqinfo(gpos)' and 'seqinfo(genome)' use
    ## different genome names (e.g. GRCh38.p2 and GRCh38) so we get rid of
    ## the genome information in 'gpos'.
    genome(gpos) <- NA_character_
    ref_allele <- getSeq(genome, gpos)

    ## Compute 'genome_compat'.
    ## Note that a small percentage of SNPs in dbSNP have alleles that
    ## are inconsistent with the reference genome (don't ask me why).
    af <- alphabetFrequency(alleles, baseOnly=TRUE)[ , 1:4, drop=FALSE]
    ref_af <- alphabetFrequency(ref_allele, baseOnly=TRUE)[ , 1:4, drop=FALSE]
    genome_compat <- matrixStats::rowAlls(af >= ref_af)

    ## Compute 'alt_alleles'.
    gpos_len <- length(gpos)
    idx <- relist(as.logical(t((af - ref_af) >= 1L)),
                  PartitioningByEnd(seq_len(gpos_len) * 4L))
    alt_alleles <- rep.int(CharacterList(colnames(af)), gpos_len)[idx]

    DataFrame(genome_compat=genome_compat,
              ref_allele=as.character(ref_allele),
              alt_alleles=alt_alleles)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .as_naked_GPos() and .as_GPos()
###

.as_naked_GPos <- function(df, seqinfo)
{
    ans_seqnames <- df[ , "seqnames"]
    ans_pos <- IPos(df[ , "pos"])
    ans_strand <- Rle(strand("*"), nrow(df))
    GPos(ans_seqnames, ans_pos, ans_strand, seqinfo=seqinfo)
}

.as_GPos <- function(df, seqinfo, drop.rs.prefix=FALSE, genome=NULL)
{
    gpos <- .as_naked_GPos(df, seqinfo)
    alleles_as_ambig <- decode_bytes_as_letters(df[ , "alleles"])
    rowids <- df$rowids
    if (is.null(rowids)) {
        ans_mcols <- DataFrame(alleles_as_ambig=alleles_as_ambig)
    } else {
        if (!drop.rs.prefix && length(gpos) != 0L)
            rowids <- paste0("rs", rowids)
        ans_mcols <- DataFrame(RefSNP_id=rowids,
                               alleles_as_ambig=alleles_as_ambig)
    }
    mcols(gpos) <- ans_mcols
    if (!is.null(genome))
        mcols(gpos) <- cbind(ans_mcols, inferRefAndAltAlleles(gpos, genome))
    gpos
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SNP extractors: snpsBySeqname(), snpsByOverlaps(), snpsById()
###

.snpsBySeqname_ODLT_SNPlocs <- function(x, seqnames, drop.rs.prefix=FALSE,
                                                     genome=NULL)
{
    if (!isTRUEorFALSE(drop.rs.prefix))
        stop(wmsg("'drop.rs.prefix' must be TRUE or FALSE"))
    if (!is.null(genome)) {
        genome <- getBSgenome(genome)
        seqlevelsStyle(genome) <- "NCBI"
    }

    df <- getBatchesBySeqnameFromOnDiskLongTable(x@snp_table, seqnames,
                                                 with.rowids=TRUE)
    x_spatial_index <- spatialIndex(x@snp_table)
    x_seqinfo <- seqinfo(x_spatial_index)
    .as_GPos(df, x_seqinfo, drop.rs.prefix=drop.rs.prefix, genome=genome)
}

setMethod("snpsBySeqname", "ODLT_SNPlocs", .snpsBySeqname_ODLT_SNPlocs)

.snpsByOverlaps_ODLT_SNPlocs <- function(x, ranges, drop.rs.prefix=FALSE, ...,
                                                    genome=NULL)
{
    ranges <- normarg_ranges(ranges)
    if (!isTRUEorFALSE(drop.rs.prefix))
        stop(wmsg("'drop.rs.prefix' must be TRUE or FALSE"))
    if (!is.null(genome)) {
        genome <- getBSgenome(genome)
        seqlevelsStyle(genome) <- "NCBI"
    }

    dots <- list(...)
    if (isTRUE(dots$invert))
        stop(wmsg("snpsByOverlaps() does not support 'invert=TRUE'"))

    if (is.null(maxgap <- dots$maxgap))
        maxgap <- -1L
    if (is.null(minoverlap <- dots$minoverlap))
        minoverlap <- 0L
    df <- getBatchesByOverlapsFromOnDiskLongTable(x@snp_table, ranges,
                                                  maxgap=maxgap,
                                                  minoverlap=minoverlap,
                                                  with.rowids=TRUE)
    x_spatial_index <- spatialIndex(x@snp_table)
    x_seqinfo <- seqinfo(x_spatial_index)
    gpos0 <- .as_naked_GPos(df, x_seqinfo)
    idx <- which(overlapsAny(gpos0, ranges, ...))
    df <- df[idx, ]
    .as_GPos(df, x_seqinfo, drop.rs.prefix=drop.rs.prefix, genome=genome)
}

setMethod("snpsByOverlaps", "ODLT_SNPlocs", .snpsByOverlaps_ODLT_SNPlocs)

.snpsById_ODLT_SNPlocs <- function(x, ids,
                                   ifnotfound=c("error", "warning", "drop"),
                                   genome=NULL)
{
    user_rowids <- ids2rowids(ids)
    ifnotfound <- match.arg(ifnotfound)
    if (!is.null(genome)) {
        genome <- getBSgenome(genome)
        seqlevelsStyle(genome) <- "NCBI"
    }

    x_rowids_env <- get_rowids_env(x@snp_table)
    rowidx <- rowids2rowidx(user_rowids, ids, x_rowids_env, ifnotfound)

    df <- getRowsFromOnDiskLongTable(x@snp_table, rowidx[[1L]],
                                     with.rowids=FALSE)
    df$rowids <- rowidx[[2L]]
    x_spatial_index <- spatialIndex(x@snp_table)
    x_seqinfo <- seqinfo(x_spatial_index)
    .as_GPos(df, x_seqinfo, drop.rs.prefix=TRUE, genome=genome)
}

setMethod("snpsById", "ODLT_SNPlocs", .snpsById_ODLT_SNPlocs)

