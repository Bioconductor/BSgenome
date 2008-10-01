### =========================================================================
### The BSgenomeForge functions
### -------------------------------------------------------------------------


.isSingleString <- function(x)
{
    is.character(x) && length(x) == 1 && !is.na(x)
}

.isSingleStringOrNA <- function(x)
{
    is.vector(x) && is.atomic(x) && length(x) == 1 && (is.character(x) || is.na(x))
}

.getMasksObjname <- function(seqnames)
{
    if (length(seqnames) == 0)
        return(character(0))
    paste(seqnames, ".masks", sep="")
}

.saveObject <- function(object, objname, destdir=".", verbose=TRUE)
{
    assign(objname, object)
    destfile <- file.path(destdir, paste(objname, ".rda", sep=""))
    if (verbose)
        cat("Saving '", objname, "' object to compressed data file '",
            destfile, "'... ", sep="")
    save(list=objname, file=destfile, compress=TRUE)
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
    if (length(seqnames) == 0) {
        warning("'seqnames' is empty")
        return(character(0))
    }
    if (!.isSingleString(prefix))
        stop("'prefix' must be a single string")
    if (!.isSingleString(suffix))
        stop("'suffix' must be a single string")
    if (!.isSingleString(seqs_srcdir))
        stop("'seqs_srcdir' must be a single string")
    srcfiles <- paste(prefix, seqnames, suffix, sep="")
    ans <- file.path(seqs_srcdir, srcfiles)
    is_OK <- file.exists(ans)
    if (!all(is_OK)) {
        missing_files <- paste(ans[!is_OK], collapse=", ")
        stop(missing_files, ": file(s) not found")
    }
    names(ans) <- seqnames
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "getSeqlengths" and "forgeSeqlengthsFile" functions.
###

getSeqlengths <- function(seqnames, prefix="", suffix=".fa", seqs_srcdir=".")
{
    if (length(seqnames) == 0) {
        warning("'seqnames' is empty")
        return(integer(0))
    }
    srcpaths <- getSeqSrcpaths(seqnames, prefix=prefix, suffix=suffix,
                               seqs_srcdir=seqs_srcdir)
    sapply(seqnames, function(seqname)
           {
               srcpath <- srcpaths[[seqname]]
               ans <- fasta.info(srcpath)
               if (length(ans) == 0)
                   stop("In file '", srcpath, "': no sequence found")
               if (length(ans) > 1)
                   warning("In file '", srcpath, "': ", length(ans),
                           " sequences found, using first sequence only")
               if (names(ans)[1] != seqname)
                   warning("In file '", srcpath, "': sequence description \"",
                           names(ans), "\" doesn't match user-specified ",
                           "sequence name \"", seqname, "\"")
               ans[[1]]
           },
           USE.NAMES=TRUE
    )
}

forgeSeqlengthsFile <- function(seqnames, prefix="", suffix=".fa",
                                seqs_srcdir=".", seqs_destdir=".",
                                verbose=TRUE)
{
    seqlengths <- getSeqlengths(seqnames, prefix=prefix, suffix=suffix,
                                seqs_srcdir=seqs_srcdir)
    if (!.isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    .saveObject(seqlengths, "seqlengths", destdir=seqs_destdir,
                verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeSeqFiles" function.
###

.forgeSeqFile <- function(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                          is.single.seq=TRUE, verbose=TRUE)
{
    if (!.isSingleString(name))
        stop("'name' must be a single string")
    srcpath <- getSeqSrcpaths(name, prefix=prefix, suffix=suffix,
                              seqs_srcdir=seqs_srcdir)
    if (verbose)
        cat("Loading FASTA file '", srcpath, "' in '", name, "' object... ", sep="")
    seq <- read.DNAStringSet(srcpath, "fasta")
    if (verbose)
        cat("DONE\n")
    if (is.single.seq) {
        if (length(seq) == 0)
            stop("file contains no DNA sequence")
        if (length(seq) > 1)
            warning("file contains ", length(seq), " sequences, ",
                    "using the first sequence only")
        seq <- seq[[1]] # now 'seq' is a DNAString object
    }
    .saveObject(seq, name, destdir=seqs_destdir, verbose=verbose)
}

forgeSeqFiles <- function(seqnames, mseqnames=NULL, prefix="", suffix=".fa",
                          seqs_srcdir=".", seqs_destdir=".", verbose=TRUE)
{
    if (length(seqnames) == 0) {
        warning("'seqnames' is empty")
    } else {
        ## just for the side effect of checking the arguments
        getSeqSrcpaths(seqnames, prefix=prefix, suffix=suffix,
                       seqs_srcdir=seqs_srcdir)
    }
    if (length(mseqnames) != 0) {
        if (!is.character(mseqnames))
            stop("'mseqnames' must be a character vector (or NULL)")
        ## just for the side effect of checking the arguments
        getSeqSrcpaths(mseqnames, prefix=prefix, suffix=suffix,
                       seqs_srcdir=seqs_srcdir)
    }
    if (!.isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    for (name in seqnames) {
        .forgeSeqFile(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                      is.single.seq=TRUE, verbose=verbose)
    }
    for (name in mseqnames) {
        .forgeSeqFile(name, prefix, suffix, seqs_srcdir, seqs_destdir,
                      is.single.seq=FALSE, verbose=verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeMasksFiles" function.
###

## AGAPS is the mask of "assembly gaps" (inter-contig Ns).
## If 'filetype' is:
##        NA: then AGAPS is an empty mask;
##     "gap": then AGAPS is extracted from a UCSC "gap" file;
##     "agp": then AGAPS is extracted from an NCBI "agp" file.
## AGAPS is active by default.
.forge.AGAPS.mask <- function(seqname, mask_width, masks_srcdir,
                              filetype, filename,
                              fileprefix, filesuffix)
{
    if (!.isSingleStringOrNA(filetype))
        stop("'filetype' must be a single string or NA")
    if (!.isSingleStringOrNA(filename))
        stop("'filename' must be a single string or NA")
    if (!.isSingleStringOrNA(fileprefix))
        stop("'fileprefix' must be a single string or NA")
    if (!.isSingleStringOrNA(filesuffix))
        stop("'filesuffix' must be a single string or NA")
    if (is.na(filetype)) {
        ans <- Mask(mask_width)
        desc(ans) <- "assembly gaps (empty)"
    } else {
        if (is.na(filename))
            filename <- paste(fileprefix, seqname, filesuffix, sep="")
        filepath <- file.path(masks_srcdir, filename)
        if (filetype == "gap")
            ans <- read.gapMask(filepath, mask_width, seqname=seqname)
        else
            ans <- read.agpMask(filepath, mask_width, seqname=seqname)
    }
    active(ans) <- TRUE
    names(ans) <- "AGAPS"
    ans
}

## AMB is the mask of "intra-contig ambiguities".
## AMB is active by default.
.forge.AMB.mask <- function(seq, AGAPS)
{
    active(AGAPS) <- TRUE
    masks(seq) <- AGAPS
    amb_letters <- names(IUPAC_CODE_MAP)[nchar(IUPAC_CODE_MAP) > 1]
    for (amb_letter in amb_letters)
        seq <- maskMotif(seq, amb_letter)
    ans <- reduce(masks(seq)[-1])
    desc(ans) <- "intra-contig ambiguities"
    if (isEmpty(ans))
        desc(ans) <- paste(desc(ans), "(empty)")
    active(ans) <- TRUE
    names(ans) <- "AMB"
    ans
}

## RM is the "RepeatMasker" mask (from the RepeatMasker .out file).
## RM is NOT active by default.
.forge.RM.mask <- function(seqname, mask_width, masks_srcdir)
{
    file <- file.path(masks_srcdir, paste(seqname, ".fa.out", sep=""))
    if (file.exists(file)) {
        ans <- read.rmMask(file, mask_width)
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
.forge.TRF.mask <- function(seqname, mask_width, masks_srcdir)
{
    file <- file.path(masks_srcdir, paste(seqname, ".bed", sep=""))
    if (file.exists(file)) {
        ans <- read.trfMask(file, mask_width)
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
                            seqs_destdir=".", masks_srcdir=".", masks_destdir=".",
                            AGAPSfiles.type="gap", AGAPSfiles.name=NA,
                            AGAPSfiles.prefix="", AGAPSfiles.suffix="_gap.txt",
                            verbose=TRUE)
{
    if (!.isSingleString(seqname))
        stop("'seqname' must be a single string")
    if (!is.numeric(nmask_per_seq)
     || length(nmask_per_seq) != 1
     || !(nmask_per_seq %in% 0:4))
        stop("'nmask_per_seq' must be 0, 1, 2, 3 or 4")
    if (nmask_per_seq == 0)
        warning("forging an empty mask collection ('nmask_per_seq' is set to 0)")
    if (!.isSingleString(seqs_destdir))
        stop("'seqs_destdir' must be a single string")
    if (!.isSingleString(masks_srcdir))
        stop("'masks_srcdir' must be a single string")
    if (!.isSingleString(masks_destdir))
        stop("'masks_destdir' must be a single string")

    ## Get the length of the sequence.
    seqfile <- file.path(seqs_destdir, paste(seqname, ".rda", sep=""))
    load(seqfile)
    seq <- get(seqname)
    mask_width <- length(seq)

    ## Start with an empty mask collection (i.e. a MaskCollection of
    ## length 0).
    masks <- new("MaskCollection", width=mask_width)
    if (nmask_per_seq >= 1) {
        AGAPS <- .forge.AGAPS.mask(seqname, mask_width, masks_srcdir,
                                   AGAPSfiles.type, AGAPSfiles.name,
                                   AGAPSfiles.prefix, AGAPSfiles.suffix)
        masks <- append(masks, AGAPS)
    }
    if (nmask_per_seq >= 2) {
        AMB <- .forge.AMB.mask(seq, AGAPS)
        masks <- append(masks, AMB)
    }
    remove(seq, list=seqname)
    if (nmask_per_seq >= 3) {
        RM <- .forge.RM.mask(seqname, mask_width, masks_srcdir)
        masks <- append(masks, RM)
    }
    if (nmask_per_seq >= 4) {
        TRF <- .forge.TRF.mask(seqname, mask_width, masks_srcdir)
        masks <- append(masks, TRF)
    }
    objname <- .getMasksObjname(seqname)
    .saveObject(masks, objname, destdir=masks_destdir, verbose=verbose)
}

forgeMasksFiles <- function(seqnames, nmask_per_seq,
                            seqs_destdir=".", masks_srcdir=".", masks_destdir=".",
                            AGAPSfiles.type="gap", AGAPSfiles.name=NA,
                            AGAPSfiles.prefix="", AGAPSfiles.suffix="_gap.txt",
                            verbose=TRUE)
{
    if (length(seqnames) == 0)
        warning("'seqnames' is empty")
    for (seqname in seqnames) {
        .forgeMasksFile(seqname, nmask_per_seq,
                        seqs_destdir=seqs_destdir,
                        masks_srcdir=masks_srcdir, masks_destdir=masks_destdir,
                        AGAPSfiles.type=AGAPSfiles.type, AGAPSfiles.name=AGAPSfiles.name,
                        AGAPSfiles.prefix=AGAPSfiles.prefix, AGAPSfiles.suffix=AGAPSfiles.suffix,
                        verbose=verbose)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "BSgenomeDataPkgSeed" class.
###

setClass(
    "BSgenomeDataPkgSeed",
    representation(
        Package="character",
        Title="character",
        Description="character",
        Version="character",
        Author="character",
        Maintainer="character",
        License="character",
        organism="character",
        species="character",
        provider="character",
        provider_version="character",
        release_date="character",
        release_name="character",
        source_url="character",
        organism_biocview="character",
        BSgenomeObjname="character",
        seqfiles.prefix="character",
        seqfiles.suffix="character",
        seqnames="character",         # a single string containing R source code
        mseqnames="character",        # a single string containing R source code
        nmask_per_seq="integer",      # a single integer
        AGAPSfiles.type="character",
        AGAPSfiles.name="character",
        AGAPSfiles.prefix="character",
        AGAPSfiles.suffix="character",
        PkgDetails="character",
        SrcDataFiles1="character",
        SrcDataFiles2="character",
        PkgExamples="character"
    ),
    prototype(
        Author="H. Pages",
        Maintainer="Biocore Team c/o BioC user list <bioconductor@stat.math.ethz.ch>",
        License="Artistic-2.0",
        source_url="-- information not available --",
        seqfiles.prefix="",
        seqfiles.suffix=".fa",
        seqnames="NULL",              # equivalent to "character(0)"
        mseqnames="NULL",
        nmask_per_seq=0L,
        AGAPSfiles.type="gap",
        AGAPSfiles.name=as.character(NA),
        AGAPSfiles.prefix="",
        AGAPSfiles.suffix="_gap.txt",
        PkgDetails="",
        SrcDataFiles1="-- information not available --",
        SrcDataFiles2="",
        PkgExamples=""
    )
)   


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "forgeBSgenomeDataPkg" generic and methods.
###

setGeneric("forgeBSgenomeDataPkg", signature="x",
    function(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)
        standardGeneric("forgeBSgenomeDataPkg")
)

setMethod("forgeBSgenomeDataPkg", "BSgenomeDataPkgSeed",
    function(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        template_path <- system.file("BSgenomeDataPkg-template", package="BSgenome")
        BSgenome_version <- installed.packages()['BSgenome','Version']
        .seqnames <- eval(parse(text=x@seqnames))
        symvals <- list(
            PKGTITLE=x@Title,
            PKGDESCRIPTION=x@Description,
            PKGVERSION=x@Version,
            AUTHOR=x@Author,
            MAINTAINER=x@Maintainer,
            BSGENOMEVERSION=BSgenome_version,
            LIC=x@License,
            ORGANISM=x@organism,
            SPECIES=x@species,
            PROVIDER=x@provider,
            PROVIDERVERSION=x@provider_version,
            RELEASEDATE=x@release_date,
            RELEASENAME=x@release_name,
            SOURCEURL=x@source_url,
            ORGANISMBIOCVIEW=x@organism_biocview,
            BSGENOMEOBJNAME=x@BSgenomeObjname,
            SEQNAMES=x@seqnames,
            MSEQNAMES=x@mseqnames,
            NMASKPERSEQ=as.character(x@nmask_per_seq),
            PKGDETAILS=x@PkgDetails,
            SRCDATAFILES1=x@SrcDataFiles1,
            SRCDATAFILES2=x@SrcDataFiles2,
            PKGEXAMPLES=x@PkgExamples
        )
        ## Should never happen
        if (any(duplicated(names(symvals)))) {
            str(symvals)
            stop("'symvals' contains duplicated symbols")
        }
        ## All symvals should by single strings (non-NA)
        is_OK <- sapply(symvals, function(val) {.isSingleString(val)})
        if (!all(is_OK)) {
            bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
            stop("values for symbols ", bad_syms, " are not single strings")
        }
        createPackage(x@Package, destdir, template_path, symvals)
        pkgdir <- file.path(destdir, x@Package)
        ## This sourcing defines the '.seqnames', '.mseqnames'
        ## and '.nmask_per_seq' variables
        source(file.path(pkgdir, "R", "zzz.R"), local=TRUE)
        ## Forge the "seqlengths.rda" file
        seqs_destdir <- file.path(pkgdir, "inst", "extdata")
        forgeSeqlengthsFile(.seqnames,
                            prefix=x@seqfiles.prefix, suffix=x@seqfiles.suffix,
                            seqs_srcdir=seqs_srcdir,
                            seqs_destdir=seqs_destdir,
                            verbose=verbose)
        ## Forge the sequence "*.rda" files
        forgeSeqFiles(.seqnames, mseqnames=.mseqnames,
                      prefix=x@seqfiles.prefix, suffix=x@seqfiles.suffix,
                      seqs_srcdir=seqs_srcdir,
                      seqs_destdir=seqs_destdir,
                      verbose=verbose)
        if (.nmask_per_seq > 0) {
            ## Forge the "*.masks.rda" files
            masks_destdir <- file.path(pkgdir, "inst", "extdata")
            forgeMasksFiles(.seqnames, .nmask_per_seq,
                            seqs_destdir=seqs_destdir,
                            masks_srcdir=masks_srcdir, masks_destdir=masks_destdir,
                            AGAPSfiles.type=x@AGAPSfiles.type, AGAPSfiles.name=x@AGAPSfiles.name,
                            AGAPSfiles.prefix=x@AGAPSfiles.prefix, AGAPSfiles.suffix=x@AGAPSfiles.suffix,
                            verbose=verbose)
        }
    }
)

setMethod("forgeBSgenomeDataPkg", "list",
    function(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        storage.mode(x$nmask_per_seq) <- "integer"
        x$Class <- "BSgenomeDataPkgSeed"
        y <- do.call("new", x)
        forgeBSgenomeDataPkg(y,
            seqs_srcdir=seqs_srcdir,
            masks_srcdir=masks_srcdir,
            destdir=destdir, verbose=verbose)
    }
)

.removeCommentsFromFile <- function(infile, outfile)
{
    if (!.isSingleString(infile))
        stop("'infile' must be a single string")
    if (!.isSingleString(outfile))
        stop("'outfile' must be a single string")
    if (file.exists(outfile))
        stop("file '", outfile, "' already exists")
    infile <- file(infile, "r")
    #on.exit(close(infile))
    outfile <- file(outfile, "w")
    #on.exit(close(outfile)) # doesn't seem to work
    while (TRUE) {
        text <- readLines(infile, n=1)
        if (length(text) == 0)
            break
        if (substr(text, 1, 1) != "#")
            writeLines(text, outfile)
    }
    close(outfile)
    close(infile)
}

### Return a named character vector.
.readSeedFile <- function(file, verbose=TRUE)
{
    if (!.isSingleString(file))
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
            file <- paste(file, "-seed", sep="")
            if (!file.exists(file))
                stop("seed file '", file0, "' not found")
        }
        if (verbose)
            cat("Seed file '", file0, "' not found, using file '", file, "'\n",
                sep="")
    }
    tmp_file <- file.path(tempdir(), "cleanseed9999.dcf")
    .removeCommentsFromFile(file, tmp_file)
    ans <- read.dcf(tmp_file)  # a character matrix
    file.remove(tmp_file)
    if (nrow(ans) != 1)
        stop("seed file '", file, "' must have exactly 1 record")
    ans[1, , drop=TRUE]
}

setMethod("forgeBSgenomeDataPkg", "character",
    function(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)
    {
        y <- as.list(.readSeedFile(x, verbose=verbose))
        if (missing(seqs_srcdir)) {
            seqs_srcdir <- y[["seqs_srcdir"]]
            if (is.null(seqs_srcdir))
                stop("'seqs_srcdir' argument is missing, and the ",
                     "'seqs_srcdir' field is missing in seed file")
        }
        if (missing(masks_srcdir) && !is.null(y[["masks_srcdir"]]))
            masks_srcdir <- y[["masks_srcdir"]]
        y <- y[!(names(y) %in% c("seqs_srcdir", "masks_srcdir"))]
        forgeBSgenomeDataPkg(y,
            seqs_srcdir=seqs_srcdir,
            masks_srcdir=masks_srcdir,
            destdir=destdir, verbose=verbose)
    }
)

