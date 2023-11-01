### =========================================================================
### The BSgenomeForge functions
### -------------------------------------------------------------------------


.getMasksObjname <- function(seqnames)
{
    if (length(seqnames) == 0)
        return(character(0))
    paste0(seqnames, ".masks")
}

.saveObjectToRdsFile <- function(object, objname, destdir=".", verbose=TRUE)
{
    destfile <- file.path(destdir, paste0(objname, ".rds"))
    if (verbose)
        cat("Saving '", objname, "' object to compressed data file '",
            destfile, "' ... ", sep="")
    ## Using compress="xz" (instead of compress="gzip") would produce a .rds
    ## file that is about 20% smaller on disk but it would also take almost 3x
    ## longer to load it later on with load(). Tested on hg19 chr1 with R-2.14
    ## (2011-09-20 r57033). This is why we stick to compress="gzip".
    saveRDS(object, file=destfile, compress="gzip")
    if (verbose)
        cat("DONE\n")
}

.saveObjectToRdaFile <- function(object, objname, destdir=".", verbose=TRUE)
{
    assign(objname, object)
    destfile <- file.path(destdir, paste0(objname, ".rda"))
    if (verbose)
        cat("Saving '", objname, "' object to compressed data file '",
            destfile, "' ... ", sep="")
    ## Using compress="xz" (instead of compress="gzip") would produce a .rda
    ## file that is about 20% smaller on disk but it would also take almost 3x
    ## longer to load it later on with load(). Tested on hg19 chr1 with R-2.14
    ## (2011-09-20 r57033). This is why we stick to compress="gzip".
    save(list=objname, file=destfile, compress="gzip")
    if (verbose)
        cat("DONE\n")
    remove(list=objname)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "getSeqSrcpaths" function.
###

getSeqSrcpaths <- function(seqnames, prefix="", suffix=".fa", seqs_srcdir=".")
{
    if (!is.null(seqnames) && !is.character(seqnames))
        stop("'seqnames' must be a character vector (or NULL)")
    if (length(seqnames) == 0L) {
        warning("'seqnames' is empty")
        return(setNames(character(0), character(0)))
    }
    if (!isSingleString(prefix))
        stop("'prefix' must be a single string")
    if (!isSingleString(suffix))
        stop("'suffix' must be a single string")
    if (!isSingleString(seqs_srcdir))
        stop("'seqs_srcdir' must be a single string")
    srcfiles <- paste0(prefix, seqnames, suffix)
    ans <- file.path(seqs_srcdir, srcfiles)
    is_OK <- file.exists(ans)
    if (!all(is_OK)) {
        files_not_found <- paste(ans[!is_OK], collapse=", ")
        stop("file(s) not found: ", files_not_found)
    }
    names(ans) <- seqnames
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getSeqlengths(), forgeSeqlengthsRdsFile(), and forgeSeqlengthsRdaFile()
###

getSeqlengths <- function(seqnames, prefix="", suffix=".fa", seqs_srcdir=".",
                          genome=NA_character_)
{
    srcpaths <- getSeqSrcpaths(seqnames, prefix=prefix, suffix=suffix,
                               seqs_srcdir=seqs_srcdir)
    ans <- vapply(setNames(seq_along(seqnames), seqnames),
        function(i) {
            seqname <- seqnames[[i]]
            srcpath <- srcpaths[[seqname]]
            ans <- fasta.seqlengths(srcpath)
            if (i %% 200L == 0L) {
                ## Unfortunately, the files opened by fasta.seqlengths() only
                ## get closed at garbage collection time! So if we don't
                ## explicitly call gc() every once in a while in this loop,
                ## we run into the risk of reaching the maximum number of
                ## files that can be opened simultaneously on the system (this
                ## is OS-dependent). This is a huge drawback of using
                ## XVector::open_input_files() internally to open the files!
                gc()
            }
            if (length(ans) == 0L)
                stop("In file '", srcpath, "': no sequence found")
            if (length(ans) > 1L)
                warning("In file '", srcpath, "': ", length(ans),
                        " sequences found, using first sequence only")
            if (names(ans)[[1L]] != seqname)
                warning("In file '", srcpath, "': sequence description \"",
                        names(ans), "\" doesn't match user-specified ",
                        "sequence name \"", seqname, "\"")
            ans[[1L]]
        },
        integer(1),
        USE.NAMES=TRUE
    )
    if (!is.na(genome)) {
        si <- try(Seqinfo(genome=genome), silent=TRUE)
        if (inherits(si, "try-error")) {
            warning(wmsg("genome is unknown ('Seqinfo(genome=\"",
                         genome, "\")' failed) ==> unable to check ",
                         "the lengths of the sequences in the files"))
        } else {
            expected <- seqlengths(si)[names(ans)]
            if (!identical(expected, ans))
                stop(wmsg("the sequences in the files have lengths ",
                          "that don't match the lengths reported ",
                          "by 'Seqinfo(genome=\"", genome, "\")'"))
        }
    }
    ans
}

forgeSeqlengthsRdsFile <- function(seqnames, prefix="", suffix=".fa",
                                   seqs_srcdir=".", seqs_destdir=".",
                                   genome=NA_character_, verbose=TRUE)
{
    if (!isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    seqlengths <- getSeqlengths(seqnames, prefix=prefix, suffix=suffix,
                                seqs_srcdir=seqs_srcdir, genome=genome)
    .saveObjectToRdsFile(seqlengths, "seqlengths", destdir=seqs_destdir,
                         verbose=verbose)
}

forgeSeqlengthsRdaFile <- function(seqnames, prefix="", suffix=".fa",
                                   seqs_srcdir=".", seqs_destdir=".",
                                   genome=NA_character_, verbose=TRUE)
{
    if (!isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    seqlengths <- getSeqlengths(seqnames, prefix=prefix, suffix=suffix,
                                seqs_srcdir=seqs_srcdir, genome=genome)
    .saveObjectToRdaFile(seqlengths, "seqlengths", destdir=seqs_destdir,
                         verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### forgeSeqFiles()
###

.forgeRdsSeqFile <- function(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                             is.single.seq=TRUE, verbose=TRUE)
{
    if (!isSingleString(name))
        stop("'name' must be a single string")
    srcpath <- getSeqSrcpaths(name, prefix=prefix, suffix=suffix,
                              seqs_srcdir=seqs_srcdir)
    if (verbose)
        cat("Loading FASTA file '", srcpath, "' ... ", sep="")
    seq <- readDNAStringSet(srcpath, "fasta")
    if (verbose)
        cat("DONE\n")
    if (is.single.seq) {
        if (length(seq) == 0L)
            stop("file contains no DNA sequence")
        if (length(seq) > 1L)
            warning("file contains ", length(seq), " sequences, ",
                    "using the first sequence only")
        seq <- seq[[1L]] # now 'seq' is a DNAString object
    }
    .saveObjectToRdsFile(seq, name, destdir=seqs_destdir, verbose=verbose)
}

.forgeRdaSeqFile <- function(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                             is.single.seq=TRUE, verbose=TRUE)
{
    if (!isSingleString(name))
        stop("'name' must be a single string")
    srcpath <- getSeqSrcpaths(name, prefix=prefix, suffix=suffix,
                              seqs_srcdir=seqs_srcdir)
    if (verbose)
        cat("Loading FASTA file '", srcpath, "' in '", name,
            "' object ... ", sep="")
    seq <- readDNAStringSet(srcpath, "fasta")
    if (verbose)
        cat("DONE\n")
    if (is.single.seq) {
        if (length(seq) == 0L)
            stop("file contains no DNA sequence")
        if (length(seq) > 1L)
            warning("file contains ", length(seq), " sequences, ",
                    "using the first sequence only")
        seq <- seq[[1]] # now 'seq' is a DNAString object
    }
    .saveObjectToRdaFile(seq, name, destdir=seqs_destdir, verbose=verbose)
}

.forgeFastaRzFileFromFastaFiles <- function(seqnames, prefix, suffix,
                                            seqs_srcdir, seqs_destdir,
                                            ondisk_seq_format, verbose=TRUE)
{
    if (!is.character(seqnames))
        stop("'seqnames' must be a character vector")
    dest_filename <- "single_sequences.fa"
    dest_filepath <- file.path(seqs_destdir, dest_filename)
    for (seqname in seqnames) {
        srcpath <- getSeqSrcpaths(seqname, prefix=prefix, suffix=suffix,
                                  seqs_srcdir=seqs_srcdir)
        if (verbose)
            cat("Loading '", seqname, "' sequence from FASTA file '",
                srcpath, "' ... ", sep="")
        seq <- readDNAStringSet(srcpath, "fasta")
        if (verbose)
            cat("DONE\n")
        if (length(seq) == 0L)
            stop("file contains no DNA sequence")
        if (length(seq) > 1L) {
            warning("file contains ", length(seq), " sequences, ",
                    "using the first sequence only")
            seq <- seq[1L]
        }
        if (verbose)
            cat("Appending '", seqname, "' sequence to FASTA file '",
                dest_filepath, "' ... ", sep="")
        names(seq) <- seqname
        writeXStringSet(seq, dest_filepath, append=TRUE,
                        format="fasta", width=50L)
        if (verbose)
            cat("DONE\n")
    }

    ## Create the index.
    if (ondisk_seq_format == "fa") {
        ## "fa" format
        if (verbose)
            cat("Indexing FASTA file '", dest_filepath, "' ... ", sep="")
        indexFa(dest_filepath)
        if (verbose)
            cat("DONE\n")
    } else {
        ## "fa.rz" format
        if (verbose)
            cat("Compressing FASTA file '", dest_filepath, "' ... ", sep="")
        farz_filepath <- sprintf("%s.rz", dest_filepath)
        razip(dest_filepath, dest=farz_filepath)
        unlink(dest_filepath)
        if (verbose)
            cat("DONE\n")
        if (verbose)
            cat("Indexing compressed FASTA file '", farz_filepath,
                "' ... ", sep="")
        indexFa(farz_filepath)
        if (verbose)
            cat("DONE\n")
    }
}

.replace_non_ACGTN_with_N <- function(x)
{
    af <- alphabetFrequency(x)
    non_ACGTN <- setdiff(names(af), c(DNA_BASES, "N"))
    if (sum(af[non_ACGTN]) == 0L)
        return(x)
    warning(wmsg("DNA sequence contains letters not supported by UCSC 2bit ",
                 "format (the format\n  only supports A, C, G, T, and N). ",
                 "Replacing them with Ns."))
    old <- paste(non_ACGTN, collapse="")
    new <- paste(rep.int("N", length(non_ACGTN)), collapse="")
    chartr(old, new, x)
}

.forgeTwobitFileFromFastaFiles <- function(seqnames, prefix, suffix,
                                           seqs_srcdir, seqs_destdir,
                                           verbose=TRUE)
{
    if (!is.character(seqnames))
        stop("'seqnames' must be a character vector")
    dest_filename <- "single_sequences.2bit"
    dest_filepath <- file.path(seqs_destdir, dest_filename)
    seqs <- setNames(vector(mode="list", length=length(seqnames)), seqnames)
    for (seqname in seqnames) {
        srcpath <- getSeqSrcpaths(seqname, prefix=prefix, suffix=suffix,
                                  seqs_srcdir=seqs_srcdir)
        if (verbose)
            cat("Loading '", seqname, "' sequence from FASTA file '",
                srcpath, "' ... ", sep="")
        seq <- readDNAStringSet(srcpath, "fasta")
        if (verbose)
            cat("DONE\n")
        if (length(seq) == 0L)
            stop("file contains no DNA sequence")
        if (length(seq) > 1L) {
            warning("file contains ", length(seq), " sequences, ",
                    "using the first sequence only")
            seq <- seq[1L]
        }
        seq <- .replace_non_ACGTN_with_N(seq[[1L]])
        seqs[[seqname]] <- seq
    }
    seqs <- DNAStringSet(seqs)
    if (verbose)
        cat("Writing all sequences to '", dest_filepath, "' ... ", sep="")
    export(seqs, dest_filepath, format="2bit")
    if (verbose)
        cat("DONE\n")
}

.subsetTwobitFile <- function(seqfile_name, seqnames,
                              seqs_srcdir, seqs_destdir,
                              verbose=verbose)
{
    src_filepath <- file.path(seqs_srcdir, seqfile_name)
    dest_filename <- "single_sequences.2bit"
    dest_filepath <- file.path(seqs_destdir, dest_filename)
    if (verbose)
        cat("Loading '", src_filepath, "' ... ", sep="")
    seqs <- import(src_filepath)
    if (verbose)
        cat("DONE\n")
    idx <- match(seqnames, names(seqs))
    notfound_idx <- which(is.na(idx))
    if (length(notfound_idx) != 0L)
        stop("sequence(s) not found: ",
             paste(seqnames[notfound_idx], collapse=", "))
    seqs <- seqs[idx]
    if (verbose)
        cat("Writing sequences to '", dest_filepath, "' ... ", sep="")
    export(seqs, dest_filepath, format="2bit")
    if (verbose)
        cat("DONE\n")
}

.copyTwobitFile <- function(seqfile_name,
                            seqs_srcdir, seqs_destdir,
                            verbose=TRUE)
{
    src_filepath <- file.path(seqs_srcdir, seqfile_name)
    if (!file.exists(src_filepath))
        stop("File not found: ", src_filepath)
    dest_filename <- "single_sequences.2bit"
    dest_filepath <- file.path(seqs_destdir, dest_filename)
    if (verbose)
        cat("Copying '", src_filepath, "' to '", dest_filepath, "' ... ",
            sep="")
    if (!file.copy(src_filepath, dest_filepath,
                   copy.mode=FALSE, copy.date=FALSE)) {
        if (verbose)
            stop("FAILED")
        stop("Failed to copy '", src_filepath, "' to '", dest_filepath, "'")
    }
    if (verbose)
        cat("DONE\n")
}

.sortUCSCTwobitFile <- function(genome, seqs_dir, verbose=TRUE)
{
    if (verbose)
        cat("Getting chrom info from UCSC with 'getChromInfoFromUCSC(\"",
            genome, "\")' ... ", sep="")
    chrom_info <- getChromInfoFromUCSC(genome)
    if (verbose)
        cat("DONE\n")
    filename <- "single_sequences.2bit"
    filepath <- file.path(seqs_dir, filename)
    if (verbose)
        cat("Loading '", filepath, "' ... ", sep="")
    seqs <- import(filepath)
    if (verbose)
        cat("DONE\n")
    m <- match(chrom_info[ , "chrom"], names(seqs))
    if (nrow(chrom_info) != length(seqs) || anyNA(m))
        stop(wmsg("nb of sequences in 'chromInfo' table and 2bit file ",
                  "don't match for UCSC genome ", genome))
    if (verbose)
        cat("Sorting sequences as in 'getChromInfoFromUCSC(\"",
            genome, "\")' ... ", sep="")
    seqs <- seqs[m]
    if (verbose)
        cat("DONE\n")
    if (verbose)
        cat("Checking the sequence lengths ... ")
    if (!identical(chrom_info[ , "size"], width(seqs)))
        stop(wmsg("sequence lengths in 'chromInfo' table and 2bit file ",
                  "don't match for UCSC genome ", genome))
    if (verbose)
        cat("OK\n")
    if (verbose)
        cat("Writing sequences to '", filepath, "' ... ", sep="")
    export(seqs, filepath, format="2bit")
    if (verbose)
        cat("DONE\n")
}

forgeSeqFiles <- function(provider, genome,
                          seqnames, mseqnames=NULL,
                          seqfile_name=NA, prefix="", suffix=".fa",
                          seqs_srcdir=".", seqs_destdir=".",
                          ondisk_seq_format=c("2bit", "rds", "rda", "fa.rz", "fa"),
                          verbose=TRUE)
{
    if (length(seqnames) == 0L && is.na(seqfile_name))
        warning("'seqnames' is empty")
    if (!(length(mseqnames) == 0L || is.character(mseqnames)))
        stop("'mseqnames' must be a character vector (or NULL)")
    if (!isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    ondisk_seq_format <- match.arg(ondisk_seq_format)
    if (ondisk_seq_format == "rds") {          # "rds" format
        for (name in seqnames) {
            .forgeRdsSeqFile(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                             is.single.seq=TRUE, verbose=verbose)
        }
    } else if (ondisk_seq_format == "rda") {   # "rda" format
        for (name in seqnames) {
            .forgeRdaSeqFile(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                             is.single.seq=TRUE, verbose=verbose)
        }
    } else if (ondisk_seq_format == "2bit") {  # "2bit" format
        if (is.na(seqfile_name)) {
            .forgeTwobitFileFromFastaFiles(seqnames, prefix, suffix,
                                           seqs_srcdir, seqs_destdir,
                                           verbose=verbose)
        } else if (length(seqnames) != 0L) {
            .subsetTwobitFile(seqfile_name, seqnames,
                              seqs_srcdir, seqs_destdir,
                              verbose=verbose)
        } else {
            .copyTwobitFile(seqfile_name,
                            seqs_srcdir, seqs_destdir,
                            verbose=verbose)
            if (provider == "UCSC") {
                UCSC_genomes <- registered_UCSC_genomes()
                if (genome %in% UCSC_genomes[ , "genome"])
                    .sortUCSCTwobitFile(genome, seqs_destdir,
                                        verbose=verbose)
            }
        }
    } else {  # "fa" and "fa.rz" formats
        .forgeFastaRzFileFromFastaFiles(seqnames, prefix, suffix,
                                        seqs_srcdir, seqs_destdir,
                                        ondisk_seq_format, verbose=verbose)
    }
    for (name in mseqnames) {
        .forgeRdaSeqFile(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                         is.single.seq=FALSE, verbose=verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .copy_CITATION_file()
###

.copy_CITATION_file <- function(citation_file, pkgdir, verbose=TRUE)
{
    if (!file.exists(citation_file))
        stop("File not found: ", citation_file)
    dest_filename <- "CITATION"
    dest_filepath <- file.path(pkgdir, "inst", dest_filename)
    if (verbose)
        cat("Copying '", citation_file, "' to '", dest_filepath, "' ... ",
            sep="")
    if (!file.copy(citation_file, dest_filepath,
                   copy.mode=FALSE, copy.date=FALSE)) {
        if (verbose)
            stop("FAILED")
        stop("Failed to copy '", citation_file, "' to '", dest_filepath, "'")
    }
    if (verbose)
        cat("DONE\n")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### forgeMasksFiles()
###

## AGAPS is the mask of "assembly gaps" (inter-contig Ns).
## If 'filetype' is:
##        NA: then AGAPS is an empty mask;
##     "gap": then AGAPS is extracted from a UCSC "gap" file;
##     "agp": then AGAPS is extracted from an NCBI "agp" file.
## AGAPS is active by default.
.forge.AGAPSmask <- function(seqname, mask_width, masks_srcdir,
                             filetype, filename,
                             fileprefix, filesuffix)
{
    if (!isSingleStringOrNA(filetype))
        stop("'filetype' must be a single string or NA")
    if (!isSingleStringOrNA(filename))
        stop("'filename' must be a single string or NA")
    if (!isSingleStringOrNA(fileprefix))
        stop("'fileprefix' must be a single string or NA")
    if (!isSingleStringOrNA(filesuffix))
        stop("'filesuffix' must be a single string or NA")
    if (is.na(filetype)) {
        ans <- Mask(mask_width)
        desc(ans) <- "assembly gaps (empty)"
    } else {
        if (is.na(filename))
            filename <- paste0(fileprefix, seqname, filesuffix)
        filepath <- file.path(masks_srcdir, filename)
        if (filetype == "gap")
            ans <- read.gapMask(filepath, seqname=seqname, mask.width=mask_width)
        else
            ans <- read.agpMask(filepath, seqname=seqname, mask.width=mask_width)
    }
    active(ans) <- TRUE
    names(ans) <- "AGAPS"
    ans
}

## AMB is the mask of "intra-contig ambiguities".
## AMB is active by default.
.forge.AMBmask <- function(seq, AGAPSmask)
{
    active(AGAPSmask) <- TRUE
    masks(seq) <- AGAPSmask
    amb_letters <- names(IUPAC_CODE_MAP)[nchar(IUPAC_CODE_MAP) > 1]
    for (amb_letter in amb_letters)
        seq <- maskMotif(seq, amb_letter)
    ans <- collapse(masks(seq)[-1])
    desc(ans) <- "intra-contig ambiguities"
    if (isEmpty(ans))
        desc(ans) <- paste0(desc(ans), " (empty)")
    active(ans) <- TRUE
    names(ans) <- "AMB"
    ans
}

## RM is the "RepeatMasker" mask (from the RepeatMasker .out file).
## RM is NOT active by default.
.forge.RMmask <- function(seqname, mask_width, masks_srcdir,
                          filename, fileprefix, filesuffix)
{
    if (!isSingleStringOrNA(filename))
        stop("'filename' must be a single string or NA")
    if (!isSingleStringOrNA(fileprefix))
        stop("'fileprefix' must be a single string or NA")
    if (!isSingleStringOrNA(filesuffix))
        stop("'filesuffix' must be a single string or NA")
    if (is.na(filename))
        filename <- paste0(fileprefix, seqname, filesuffix)
    filepath <- file.path(masks_srcdir, filename)
    if (file.exists(filepath)) {
        ans <- read.rmMask(filepath, seqname=seqname, mask.width=mask_width)
        desc(ans) <- "RepeatMasker"
    } else {
        ans <- Mask(mask_width)
        desc(ans) <- "RepeatMasker (empty)"
    }
    active(ans) <- FALSE
    names(ans) <- "RM"
    ans
}

## TRF is the "Tandem Repeats Finder" mask (from the Tandem Repeats Finder
## .bed file).
## TRF is NOT active by default.
.forge.TRFmask <- function(seqname, mask_width, masks_srcdir,
                           filename, fileprefix, filesuffix)
{
    if (!isSingleStringOrNA(filename))
        stop("'filename' must be a single string or NA")
    if (!isSingleStringOrNA(fileprefix))
        stop("'fileprefix' must be a single string or NA")
    if (!isSingleStringOrNA(filesuffix))
        stop("'filesuffix' must be a single string or NA")
    if (is.na(filename))
        filename <- paste0(fileprefix, seqname, filesuffix)
    filepath <- file.path(masks_srcdir, filename)
    if (file.exists(filepath)) {
        ans <- read.trfMask(filepath, seqname=seqname, mask.width=mask_width)
        desc(ans) <- "Tandem Repeats Finder [period<=12]"
    } else {
        ans <- Mask(mask_width)
        desc(ans) <- "Tandem Repeats Finder [period<=12] (empty)"
    }
    active(ans) <- FALSE
    names(ans) <- "TRF"
    ans
}

.forgeMasksFile <- function(seqname, nmask_per_seq,
                            seqs_destdir=".",
                            ondisk_seq_format=c("2bit", "rda", "fa.rz", "fa"),
                            masks_srcdir=".", masks_destdir=".",
                            AGAPSfiles_type="gap", AGAPSfiles_name=NA,
                            AGAPSfiles_prefix="", AGAPSfiles_suffix="_gap.txt",
                            RMfiles_name=NA, RMfiles_prefix="", RMfiles_suffix=".fa.out",
                            TRFfiles_name=NA, TRFfiles_prefix="", TRFfiles_suffix=".bed",
                            verbose=TRUE)
{
    if (!isSingleString(seqname))
        stop("'seqname' must be a single string")
    if (!is.numeric(nmask_per_seq)
     || length(nmask_per_seq) != 1
     || !(nmask_per_seq %in% 0:4))
        stop("'nmask_per_seq' must be 0, 1, 2, 3 or 4")
    if (nmask_per_seq == 0)
        warning("forging an empty mask collection ('nmask_per_seq' is set to 0)")
    if (!isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    if (!isSingleString(masks_srcdir))
        stop("'masks_srcdir' must be a single string")
    if (!isSingleString(masks_destdir))
        stop("'masks_destdir' must be a single string")

    ## Load the sequence.
    ondisk_seq_format <- match.arg(ondisk_seq_format)
    if (ondisk_seq_format == "rda") {  # "rda" format
        seqfile <- file.path(seqs_destdir, paste0(seqname, ".rda"))
        load(seqfile)
        seq <- get(seqname)
        remove(list=seqname)
    } else if (ondisk_seq_format == "2bit") {  # "2bit" format
        twobit_filename <- "single_sequences.2bit"
        twobit_filepath <- file.path(seqs_destdir, twobit_filename)
        twobitfile <- TwoBitFile(twobit_filepath)
        which <- GRanges(seqname, IRanges(1L, seqlengths(twobitfile)[[seqname]]))
        seq <- import(twobitfile, which=which)[[1L]]
    } else {  # "fa" and "fa.rz" formats
        fa_filename <- "single_sequences.fa"
        if (ondisk_seq_format == "fa.rz")
            fa_filename <- paste0(fa_filename, ".rz")
        fa_filepath <- file.path(seqs_destdir, fa_filename)
        fafile <- FaFile(fa_filepath)
        ## TODO: Implement "[[" method for FaFile objects (in Rsamtools) and
        ## use it here (i.e. do 'seq <- fafile[[seqname]]' instead of the 2
        ## lines below).
        param <- GRanges(seqname, IRanges(1L, seqlengths(fafile)[[seqname]]))
        seq <- scanFa(fafile, param=param)[[1L]]
    }
    mask_width <- length(seq)

    ## Start with an empty mask collection (i.e. a MaskCollection of
    ## length 0).
    masks <- new("MaskCollection", width=mask_width)
    if (nmask_per_seq >= 1) {
        AGAPSmask <- .forge.AGAPSmask(seqname, mask_width, masks_srcdir,
                                      AGAPSfiles_type, AGAPSfiles_name,
                                      AGAPSfiles_prefix, AGAPSfiles_suffix)
        masks <- append(masks, AGAPSmask)
    }
    if (nmask_per_seq >= 2) {
        AMBmask <- .forge.AMBmask(seq, AGAPSmask)
        masks <- append(masks, AMBmask)
    }
    if (nmask_per_seq >= 3) {
        RMmask <- .forge.RMmask(seqname, mask_width, masks_srcdir,
                                RMfiles_name, RMfiles_prefix, RMfiles_suffix)
        masks <- append(masks, RMmask)
    }
    if (nmask_per_seq >= 4) {
        TRFmask <- .forge.TRFmask(seqname, mask_width, masks_srcdir,
                                  TRFfiles_name, TRFfiles_prefix, TRFfiles_suffix)
        masks <- append(masks, TRFmask)
    }
    objname <- .getMasksObjname(seqname)
    .saveObjectToRdaFile(masks, objname, destdir=masks_destdir, verbose=verbose)
}

forgeMasksFiles <- function(seqnames, nmask_per_seq,
                            seqs_destdir=".",
                            ondisk_seq_format=c("2bit", "rda", "fa.rz", "fa"),
                            masks_srcdir=".", masks_destdir=".",
                            AGAPSfiles_type="gap", AGAPSfiles_name=NA,
                            AGAPSfiles_prefix="", AGAPSfiles_suffix="_gap.txt",
                            RMfiles_name=NA, RMfiles_prefix="", RMfiles_suffix=".fa.out",
                            TRFfiles_name=NA, TRFfiles_prefix="", TRFfiles_suffix=".bed",
                            verbose=TRUE)
{
    if (length(seqnames) == 0L)
        warning("'seqnames' is empty")
    for (seqname in seqnames) {
        .forgeMasksFile(seqname, nmask_per_seq,
                        seqs_destdir=seqs_destdir,
                        ondisk_seq_format=ondisk_seq_format,
                        masks_srcdir=masks_srcdir, masks_destdir=masks_destdir,
                        AGAPSfiles_type=AGAPSfiles_type, AGAPSfiles_name=AGAPSfiles_name,
                        AGAPSfiles_prefix=AGAPSfiles_prefix, AGAPSfiles_suffix=AGAPSfiles_suffix,
                        RMfiles_name=RMfiles_name, RMfiles_prefix=RMfiles_prefix, RMfiles_suffix=RMfiles_suffix,
                        TRFfiles_name=TRFfiles_name, TRFfiles_prefix=TRFfiles_prefix, TRFfiles_suffix=TRFfiles_suffix,
                        verbose=verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The BSgenomeDataPkgSeed and MaskedBSgenomeDataPkgSeed classes and their
### low-level constructors.
###

setClass("BSgenomeDataPkgSeed",
    representation(
        Package="character",
        Title="character",
        Description="character",
        Version="character",
        Author="character",
        Maintainer="character",
        Suggests="character",
        License="character",
        organism="character",
        common_name="character",
        ## Should be accepted by 'Seqinfo(genome=genome)' e.g. "TAIR10.1",
        ## "hg38", "ce11" etc...
        genome="character",
        provider="character",
        provider_version="character", # deprecated in favor of 'genome'
        release_date="character",
        release_name="character",    # deprecated (not replaced by anything)
        source_url="character",
        organism_biocview="character",
        BSgenomeObjname="character",
        seqnames="character",        # a single string containing R source code
        circ_seqs="character",       # a single string containing R source code
        mseqnames="character",       # a single string containing R source code
        PkgDetails="character",
        SrcDataFiles="character",
        PkgExamples="character",
        seqfile_name="character",
        seqfiles_prefix="character",
        seqfiles_suffix="character",
        ondisk_seq_format="character",
        citation_file="character"    # a single string containing R source code
    ),
    prototype(
        Author="The Bioconductor Dev Team",
        Maintainer="Bioconductor Package Maintainer <maintainer@bioconductor.org>",
        Suggests="",
        License="Artistic-2.0",
        genome=NA_character_,
        provider_version=NA_character_,
        release_name=NA_character_,
        source_url="-- information not available --",
        seqnames=NA_character_,
        circ_seqs=NA_character_,
        mseqnames="NULL",            # equivalent to "character(0)"
        PkgDetails="",
        SrcDataFiles="-- information not available --",
        PkgExamples="",
        seqfile_name=NA_character_,
        seqfiles_prefix="",
        seqfiles_suffix=".fa",
        ondisk_seq_format="2bit",
        citation_file=NA_character_
    )
)

setClass("MaskedBSgenomeDataPkgSeed",
    representation(
        Package="character",
        Title="character",
        Description="character",
        Version="character",
        Author="character",
        Maintainer="character",
        RefPkgname="character",
        Suggests="character",
        License="character",
        source_url="character",
        organism_biocview="character",
        nmask_per_seq="integer",     # a single integer
        PkgDetails="character",
        SrcDataFiles="character",
        PkgExamples="character",
        AGAPSfiles_type="character",
        AGAPSfiles_name="character",
        AGAPSfiles_prefix="character",
        AGAPSfiles_suffix="character",
        RMfiles_name="character",
        RMfiles_prefix="character",
        RMfiles_suffix="character",
        TRFfiles_name="character",
        TRFfiles_prefix="character",
        TRFfiles_suffix="character"
    ),
    prototype(
        Author="The Bioconductor Dev Team",
        Maintainer="Bioconductor Package Maintainer <maintainer@bioconductor.org>",
        Suggests="",
        License="Artistic-2.0",
        source_url="-- information not available --",
        nmask_per_seq=0L,
        PkgDetails="",
        SrcDataFiles="-- information not available --",
        PkgExamples="",
        AGAPSfiles_type="gap",
        AGAPSfiles_name=as.character(NA),
        AGAPSfiles_prefix="",
        AGAPSfiles_suffix="_gap.txt",
        RMfiles_name=as.character(NA),
        RMfiles_prefix="",
        RMfiles_suffix=".fa.out",
        TRFfiles_name=as.character(NA),
        TRFfiles_prefix="",
        TRFfiles_suffix=".bed"
    )
)

### Generic transformation of a named list into an S4 object with automatic
### coercion of the list elements to the required types.
makeS4FromList <- function(Class, x)
{
    if (!is.list(x) || is.null(names(x)))
        stop("'x' must be a named list")
    explicit_slots <- getSlots(Class)[names(x)]
    if (any(is.na(explicit_slots))) {
        invalid_names <- setdiff(names(x), names(getSlots(Class)))
        stop("some names in 'x' are not valid ", Class, " slots (",
             paste(invalid_names, collapse=", "), ")")
    }
    y <- lapply(seq_len(length(x)),
                function(i) {
                    x_elt <- x[[i]]
                    if (is(x_elt, explicit_slots[i]))
                        return(x_elt)
                    as(x_elt, explicit_slots[i])
                })
    names(y) <- names(x)
    y$Class <- Class
    do.call(new, y)
}

BSgenomeDataPkgSeed <- function(x)
    makeS4FromList("BSgenomeDataPkgSeed", x)

MaskedBSgenomeDataPkgSeed <- function(x)
    makeS4FromList("MaskedBSgenomeDataPkgSeed", x)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeBSgenomeDataPkg" generic and methods.
###

setGeneric("forgeBSgenomeDataPkg", signature="x",
    function(x, seqs_srcdir=".", destdir=".", replace=FALSE, verbose=TRUE)
        standardGeneric("forgeBSgenomeDataPkg")
)

setMethod("forgeBSgenomeDataPkg", "BSgenomeDataPkgSeed",
    function(x, seqs_srcdir=".", destdir=".", replace=FALSE, verbose=TRUE)
    {
        if (!isTRUEorFALSE(replace))
            stop("'replace' must be TRUE or FALSE")

        ## The Biobase package is needed for createPackage().
        if (!requireNamespace("Biobase", quietly=TRUE))
            stop("Couldn't load the Biobase package. Please install ",
                 "the Biobase\n  package and try again.")
        template_path <- system.file("pkgtemplates", "BSgenome_datapkg",
                                     package="BSgenome")
        BSgenome_version <- installed.packages()['BSgenome','Version']
        if (is.na(x@genome)) {
            if (is.na(x@provider_version))
                stop("'genome' field is missing in seed file")
            warning("field 'provider_version' is deprecated ",
                    "in favor of 'genome'")
            x@genome <- x@provider_version
        }
        if (!is.na(x@release_name))
            warning("field 'release_name' is deprecated")
        seqnames <- x@seqnames
        if (!is.na(seqnames)) {
            .seqnames <- eval(parse(text=seqnames))
        } else {
            if (is.na(x@seqfile_name)) {
                .seqnames <- seqlevels(Seqinfo(genome=x@genome))
            } else {
                .seqnames <- NULL
            }
        }
        seqnames <- deparse1(.seqnames)
        circ_seqs <- x@circ_seqs
        if (!is.na(circ_seqs)) {
            .circ_seqs <- eval(parse(text=circ_seqs))
        } else {
            si <- Seqinfo(genome=x@genome)
            .circ_seqs <- seqlevels(si)[isCircular(si)]
        }
        circ_seqs <- deparse1(.circ_seqs)
        symvals <- list(
            PKGTITLE=x@Title,
            PKGDESCRIPTION=x@Description,
            PKGVERSION=x@Version,
            AUTHOR=x@Author,
            MAINTAINER=x@Maintainer,
            BSGENOMEVERSION=BSgenome_version,
            SUGGESTS=x@Suggests,
            LICENSE=x@License,
            ORGANISM=x@organism,
            COMMONNAME=x@common_name,
            GENOME=x@genome,
            PROVIDER=x@provider,
            RELEASEDATE=x@release_date,
            SOURCEURL=x@source_url,
            ORGANISMBIOCVIEW=x@organism_biocview,
            BSGENOMEOBJNAME=x@BSgenomeObjname,
            SEQNAMES=seqnames,
            CIRCSEQS=circ_seqs,
            MSEQNAMES=x@mseqnames,
            PKGDETAILS=x@PkgDetails,
            SRCDATAFILES=x@SrcDataFiles,
            PKGEXAMPLES=x@PkgExamples
        )
        ## Should never happen
        if (any(duplicated(names(symvals)))) {
            str(symvals)
            stop("'symvals' contains duplicated symbols")
        }
        ## All symvals should by single strings (non-NA)
        is_OK <- sapply(symvals, isSingleString)
        if (!all(is_OK)) {
            bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
            stop("values for symbols ", bad_syms, " are not single strings")
        }
        pkgdir <- file.path(destdir, x@Package)
        if (file.exists(pkgdir)) {
            if (replace) {
                unlink(pkgdir, recursive=TRUE)
            } else {
                stop("directory ", pkgdir, " exists. ",
                     "Use replace=TRUE to replace it.")
            }
        }
        Biobase::createPackage(x@Package, destdir, template_path, symvals)

        .mseqnames <- eval(parse(text=x@mseqnames))
        seqs_destdir <- file.path(pkgdir, "inst", "extdata")
        if (x@ondisk_seq_format == "rds") {
            ## Forge the "seqlengths.rds" file
            forgeSeqlengthsRdsFile(.seqnames,
                                   prefix=x@seqfiles_prefix,
                                   suffix=x@seqfiles_suffix,
                                   seqs_srcdir=seqs_srcdir,
                                   seqs_destdir=seqs_destdir,
                                   genome=x@genome,
                                   verbose=verbose)
        } else if (x@ondisk_seq_format == "rda") {
            ## Forge the "seqlengths.rda" file
            forgeSeqlengthsRdaFile(.seqnames,
                                   prefix=x@seqfiles_prefix,
                                   suffix=x@seqfiles_suffix,
                                   seqs_srcdir=seqs_srcdir,
                                   seqs_destdir=seqs_destdir,
                                   genome=x@genome,
                                   verbose=verbose)
        }
        ## Forge the sequence files (either "2bit", "rds", "rda", "fa.rz",
        ## or "fa").
        forgeSeqFiles(x@provider, x@genome,
                      .seqnames, mseqnames=.mseqnames,
                      seqfile_name=x@seqfile_name,
                      prefix=x@seqfiles_prefix,
                      suffix=x@seqfiles_suffix,
                      seqs_srcdir=seqs_srcdir,
                      seqs_destdir=seqs_destdir,
                      ondisk_seq_format=x@ondisk_seq_format,
                      verbose=verbose)
        if (!is.na(x@citation_file)) {
            citation_file <- eval(parse(text=x@citation_file))
            .copy_CITATION_file(citation_file, pkgdir, verbose=verbose)
        }
    }
)

setMethod("forgeBSgenomeDataPkg", "list",
    function(x, seqs_srcdir=".", destdir=".", replace=FALSE, verbose=TRUE)
    {
        y <- BSgenomeDataPkgSeed(x)
        forgeBSgenomeDataPkg(y, seqs_srcdir=seqs_srcdir, destdir=destdir,
                                replace=replace, verbose=verbose)
    }
)

.removeCommentLines <- function(infile=stdin(), outfile=stdout())
{
    if (is.character(infile)) {
        infile <- file(infile, "r")
        on.exit(close(infile))
    }
    if (is.character(outfile)) {
        outfile <- file(outfile, "w")
        on.exit({close(infile); close(outfile)})
    }
    while (TRUE) {
        lines <- readLines(infile, n=25000L)
        if (length(lines) == 0L)
            return()
        keep_it <- substr(lines, 1L, 1L) != "#"
        writeLines(lines[keep_it], outfile)
    }
}

### A "comment aware" version of read.dcf().
read.dcf2 <- function(file, ...)
{
    clean_file <- tempfile()
    .removeCommentLines(file, clean_file)
    on.exit(file.remove(clean_file))
    read.dcf(clean_file, ...)
}

### Return a named character vector.
.readSeedFile <- function(file, verbose=TRUE)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if (file.exists(file)) {
        ## Using 'x["isdir"][[1]]' is safer than using 'x$isdir' or
        ## 'x[["isdir"]]' because it will fail if "isdir" is not a defined
        ## column
        isdir <- file.info(file)["isdir"][[1]]
        if (isdir)
            stop("'", file, "' is a directory, not a seed file")
    } else {
        file0 <- file
        seed_dir <- system.file("extdata", "GentlemanLab", package="BSgenome")
        file <- file.path(seed_dir, file)
        if (!file.exists(file)) {
            file <- paste0(file, "-seed")
            if (!file.exists(file))
                stop("seed file '", file0, "' not found")
        }
        if (verbose)
            cat("Seed file '", file0, "' not found, using file '", file, "'\n",
                sep="")
    }
    ans <- read.dcf2(file)  # a character matrix
    if (nrow(ans) != 1)
        stop("seed file '", file, "' must have exactly 1 record")
    ans[1, , drop=TRUE]
}

setMethod("forgeBSgenomeDataPkg", "character",
    function(x, seqs_srcdir=".", destdir=".", replace=FALSE, verbose=TRUE)
    {
        y <- .readSeedFile(x, verbose=verbose)
        y <- as.list(y)
        if (missing(seqs_srcdir)) {
            seqs_srcdir <- y[["seqs_srcdir"]]
            if (is.null(seqs_srcdir))
                stop("'seqs_srcdir' argument is missing, and the ",
                     "'seqs_srcdir' field is missing in seed file")
        }
        y <- y[!(names(y) %in% "seqs_srcdir")]
        forgeBSgenomeDataPkg(y, seqs_srcdir=seqs_srcdir, destdir=destdir,
                                replace=replace, verbose=verbose)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeMaskedBSgenomeDataPkg" generic and methods.
###

setGeneric("forgeMaskedBSgenomeDataPkg", signature="x",
    function(x, masks_srcdir=".", destdir=".", verbose=TRUE)
        standardGeneric("forgeMaskedBSgenomeDataPkg")
)

setMethod("forgeMaskedBSgenomeDataPkg", "MaskedBSgenomeDataPkgSeed",
    function(x, masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        ## The Biobase package is needed for createPackage().
        if (!requireNamespace("Biobase", quietly=TRUE))
            stop("Couldn't load the Biobase package. Please install ",
                 "the Biobase\n  package and try again.")
        require(x@RefPkgname, character.only=TRUE) ||
            stop("the ", x@RefPkgname, " package is required")
        ref_envir <- as.environment(paste0("package:", x@RefPkgname))
        ref_bsgenome <- get(x@RefPkgname, envir=ref_envir)
        template_path <- system.file("pkgtemplates", "MaskedBSgenome_datapkg", package="BSgenome")
        BSgenome_version <- installed.packages()['BSgenome','Version']
        symvals <- list(
            PKGTITLE=x@Title,
            PKGDESCRIPTION=x@Description,
            PKGVERSION=x@Version,
            AUTHOR=x@Author,
            MAINTAINER=x@Maintainer,
            BSGENOMEVERSION=BSgenome_version,
            REFPKGNAME=x@RefPkgname,
            SUGGESTS=x@Suggests,
            LICENSE=x@License,
            ORGANISM=metadata(ref_bsgenome)$organism,
            COMMONNAME=metadata(ref_bsgenome)$common_name,
            GENOME=metadata(ref_bsgenome)$genome,
            PROVIDER=metadata(ref_bsgenome)$provider,
            RELEASEDATE=metadata(ref_bsgenome)$release_date,
            SOURCEURL=x@source_url,
            ORGANISMBIOCVIEW=x@organism_biocview,
            NMASKPERSEQ=as.character(x@nmask_per_seq),
            PKGDETAILS=x@PkgDetails,
            SRCDATAFILES=x@SrcDataFiles,
            PKGEXAMPLES=x@PkgExamples
        )
        ## Should never happen
        if (any(duplicated(names(symvals)))) {
            str(symvals)
            stop("'symvals' contains duplicated symbols")
        }
        ## All symvals should by single strings (non-NA)
        is_OK <- sapply(symvals, isSingleString)
        if (!all(is_OK)) {
            bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
            stop("values for symbols ", bad_syms, " are not single strings")
        }
        Biobase::createPackage(x@Package, destdir, template_path, symvals)
        ## Forge the "*.masks.rda" files
        seqs_destdir <- dirname(seqlengthsFilepath(ref_bsgenome@single_sequences))
        pkgdir <- file.path(destdir, x@Package)
        masks_destdir <- file.path(pkgdir, "inst", "extdata")
        forgeMasksFiles(seqnames(ref_bsgenome), x@nmask_per_seq,
                        seqs_destdir=seqs_destdir, ondisk_seq_format="2bit",
                        masks_srcdir=masks_srcdir, masks_destdir=masks_destdir,
                        AGAPSfiles_type=x@AGAPSfiles_type, AGAPSfiles_name=x@AGAPSfiles_name,
                        AGAPSfiles_prefix=x@AGAPSfiles_prefix, AGAPSfiles_suffix=x@AGAPSfiles_suffix,
                        RMfiles_name=x@RMfiles_name, RMfiles_prefix=x@RMfiles_prefix, RMfiles_suffix=x@RMfiles_suffix,
                        TRFfiles_name=x@TRFfiles_name, TRFfiles_prefix=x@TRFfiles_prefix, TRFfiles_suffix=x@TRFfiles_suffix,
                        verbose=verbose)
    }
)

setMethod("forgeMaskedBSgenomeDataPkg", "list",
    function(x, masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        y <- MaskedBSgenomeDataPkgSeed(x)
        forgeMaskedBSgenomeDataPkg(y,
            masks_srcdir=masks_srcdir, destdir=destdir,
            verbose=verbose)
    }
)

setMethod("forgeMaskedBSgenomeDataPkg", "character",
    function(x, masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        y <- .readSeedFile(x, verbose=verbose)
        y <- as.list(y)
        if (missing(masks_srcdir)) {
            masks_srcdir <- y[["masks_srcdir"]]
            if (is.null(masks_srcdir))
                stop("'masks_srcdir' argument is missing, and the ",
                     "'masks_srcdir' field is missing in seed file")
        }
        y <- y[!(names(y) %in% "masks_srcdir")]
        forgeMaskedBSgenomeDataPkg(y,
            masks_srcdir=masks_srcdir, destdir=destdir,
            verbose=verbose)
    }
)

