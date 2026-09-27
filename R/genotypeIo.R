#' @title Genotype I/O via GenotypeHandle
#' @description Read genotype data from various formats (VCF, plink1, plink2,
#'   GDS) and provide block-level genotype extraction without requiring format
#'   conversion.
#' @name pecotmr-genotype-io
#' @keywords internal
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom tools file_ext
#' @importFrom methods as
#' @include AllGenerics.R
NULL

# =============================================================================
# Main reader method -- returns a GenotypeHandle
# =============================================================================

# Abort on an unrecognized genotype format. A named helper so a switch() default
# branch stays a single call, keeping message construction out of the signaler.
# @noRd
.abortUnsupportedFormat <- function(format, context = NULL) {
    msg <- if (is.null(context)) {
        glue("Unsupported genotype format: {format}")
    } else {
        glue("Unsupported format in {context}: {format}")
    }
    abort(msg)
}

#' @rdname readGenotypes
#' @export
setMethod(
    "readGenotypes",
    signature(path = "character"),
    function(path, format = NULL, vcfArgs = list(), ...) {
        .genotypeExperiment(
            .readGenotypeHandle(path, format = format, vcfArgs = vcfArgs)
        )
    }
)

#' @rdname readGenotypes
#' @export
setMethod(
    "readGenotypes",
    signature(path = "missing"),
    function(path, format = NULL, vcfArgs = list(), ...) {
        .genotypeExperiment(
            GenotypeHandle(format = format, vcfArgs = vcfArgs, ...)
        )
    }
)

# The handle behind a panel. Internal: the handle is the seed layer, and the
# handle-construction machinery below calls this rather than readGenotypes()
# so it does not wrap and immediately unwrap a panel on every hop.
# @noRd
.readGenotypeHandle <- function(path, format = NULL, vcfArgs = list()) {
    if (is.null(format)) {
        format <- .h2DetectFormat(path)
    }
    switch(
        format,
        "gds" = .makeGdsHandle(path),
        "vcf" = .makeVcfHandle(path, vcfArgs = vcfArgs),
        "plink1" = .makePlink1Handle(path),
        "plink2" = .makePlink2Handle(path),
        .abortUnsupportedFormat(format)
    )
}

# =============================================================================
# Handle constructors -- read metadata, defer genotype loading
# =============================================================================

# Record each variant's original 1-based position in the genotype file. This
# lets @snpInfo be row-subset (e.g. to the range of a study's summary stats)
# while genotype reads that index BY FILE POSITION still resolve correctly:
# PLINK2's ReadList(variant_subset=) and the per-chromosome PLINK2 view inside
# sharded routing read fileIdx[snpIdx]; the by-id backends (plink1/gds/vcf) look
# up snpInfo$SNP[snpIdx] / $BP[snpIdx] and ignore fileIdx entirely. For a full
# (unsubset) handle fileIdx == seq_len(nrow), so reads are unchanged.
# @noRd
.withFileIdx <- function(snpInfo) {
    mutate(snpInfo, fileIdx = seq_len(nrow(snpInfo)))
}

# Restrict a GenotypeHandle's @snpInfo to `keep` (a logical mask or integer row
# indices into @snpInfo). Genotype reads stay correct because fileIdx carries
# each kept variant's original file position; everything else (path/format/
# pgenPtr/chromPaths/samples) is preserved. NULL-safe; a no-op when nothing is
# dropped. Handles built before the fileIdx column existed are NOT subset (the
# read path would be positional) -- return them unchanged.
# @noRd
# A handle trimmed to zero variants, every other property preserved.
#
# Separate from .subsetGenotypeHandle() because that declines to subset a
# legacy handle with no fileIdx column -- the read path would be positional,
# so it hands the full panel back rather than risk a wrong read. A
# zero-variant handle is never read for extraction, so that guard has nothing
# to protect here and the empty must win. NULL-safe and idempotent.
# @noRd
.emptyGenotypeHandle <- function(handle) {
    if (is.null(handle)) {
        return(NULL)
    }
    # Empty both axes. An emptied sketch references no LD, so its sample axis is
    # dead weight (~10k anonymous names, the bulk of a skipped-region file).
    # nSamples must go to 0 alongside sampleIds, or dm (variants x nSamples)
    # would disagree with the now-empty derived sample dimnames. Reached only
    # via .emptySketch (the empty / PIP-skip path), never on a surviving region.
    methods::initialize(
        handle,
        snpInfo = slice(getSnpInfo(handle), integer(0)),
        sampleIds = character(0),
        nSamples = 0L
    )
}

.subsetGenotypeHandle <- function(handle, keep) {
    if (is.null(handle)) {
        return(NULL)
    }
    si <- getSnpInfo(handle)
    keepIdx <- if (is.logical(keep)) which(keep) else as.integer(keep)
    if (length(keepIdx) >= nrow(si)) {
        return(handle)
    } # nothing dropped
    if (!is_in("fileIdx", names(si))) {
        return(handle)
    } # legacy handle: unsafe
    methods::initialize(handle, snpInfo = slice(si, keepIdx))
}

#' @importFrom checkmate assertFileExists
#' @keywords internal
.makeGdsHandle <- function(path) {
    if (!requireNamespace("SNPRelate", quietly = TRUE)) {
        abort("Package 'SNPRelate' is required for reading GDS files.")
    }
    if (!requireNamespace("gdsfmt", quietly = TRUE)) {
        abort("Package 'gdsfmt' is required for reading GDS files.")
    }
    assertFileExists(path, access = "r", .var.name = "GDS file")

    snpInfo <- .gdsSnpInfo(path)

    sampleIds <- .withGds(path, .gdsReadSampleIds)
    nSamples <- length(sampleIds)

    new(
        "GenotypeHandle",
        path = path,
        format = "gds",
        snpInfo = .withFileIdx(snpInfo),
        nSamples = as.integer(nSamples),
        sampleIds = sampleIds,
        pgenPtr = NULL
    )
}

#' @importFrom checkmate assertFileExists
#' @keywords internal
.makeVcfHandle <- function(path, vcfArgs = list()) {
    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
        abort("Package 'VariantAnnotation' is required for reading VCF files.")
    }
    assertFileExists(path, access = "r", .var.name = "VCF file")

    hdr <- VariantAnnotation::scanVcfHeader(path)
    sampleIds <- as.character(VariantAnnotation::samples(hdr))
    nSamples <- length(sampleIds)

    param <- VariantAnnotation::ScanVcfParam(
        fixed = c("ALT"),
        info = NA,
        geno = NA
    )
    vcf <- exec(
        VariantAnnotation::readVcf,
        path,
        param = param,
        !!!vcfArgs
    )
    rd <- rowRanges(vcf)

    # pecotmr convention: A1 = ALT (effect), A2 = REF
    snpInfo <- tibble(
        SNP = names(rd),
        CHR = as.character(seqnames(rd)),
        BP = as.integer(start(rd)),
        A1 = map_chr(rd$ALT, .gtFirstAllele),
        A2 = as.character(rd$REF)
    )

    new(
        "GenotypeHandle",
        path = normalizePath(path),
        format = "vcf",
        snpInfo = .withFileIdx(snpInfo),
        nSamples = as.integer(nSamples),
        sampleIds = sampleIds,
        pgenPtr = NULL
    )
}

#' @keywords internal
.makePlink1Handle <- function(path) {
    if (!requireNamespace("snpStats", quietly = TRUE)) {
        abort("Package 'snpStats' is required for reading plink1 files.")
    }
    stem <- .plink1RequireFiles(path)
    meta <- .plink1ReadMeta(stem)
    new(
        "GenotypeHandle",
        path = stem,
        format = "plink1",
        snpInfo = .withFileIdx(meta$snpInfo),
        nSamples = as.integer(meta$nSamples),
        sampleIds = meta$sampleIds,
        pgenPtr = NULL
    )
}

# Resolve the plink1 stem and assert its .bed/.bim/.fam all exist.
# @noRd
#' @importFrom checkmate assertFileExists
.plink1RequireFiles <- function(path) {
    stem <- .plinkStem(path)
    assertFileExists(
        str_c(stem, c(".bed", ".bim", ".fam")),
        access = "r",
        .var.name = "Plink file"
    )
    stem
}

# Read the .bim / .fam sidecars into snpInfo + sample metadata. Returns
# list(snpInfo, sampleIds, nSamples).
# @noRd
.plink1ReadMeta <- function(stem) {
    # col_types all-character so readr never coerces a column: a bim whose
    # allele column is uniformly "T"/"F" (all-A/T SNPs) would otherwise become
    # logical, silently corrupting A1/A2 (BP is cast to integer below).
    bim <- read_table(
        str_c(stem, ".bim"),
        col_names = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
        col_types = cols(.default = col_character())
    )
    fam <- read_table(
        str_c(stem, ".fam"),
        col_names = FALSE,
        col_types = cols(.default = col_character())
    )
    # plink1 bim: col5 = A1 (minor/effect), col6 = A2 (major/ref); matches the
    # pecotmr convention directly.
    snpInfo <- tibble(
        SNP = bim$SNP,
        CHR = as.character(bim$CHR),
        BP = as.integer(bim$BP),
        A1 = bim$A1,
        A2 = bim$A2
    )
    list(
        snpInfo = snpInfo,
        sampleIds = as.character(fam[[2]]),
        nSamples = nrow(fam)
    )
}

#' @keywords internal
.makePlink2Handle <- function(path) {
    if (!requireNamespace("pgenlibr", quietly = TRUE)) {
        abort("Package 'pgenlibr' is required for reading plink2 files.")
    }

    stem <- .plinkStem(path)

    # Use pecotmr's resolvePlink2Paths for robust path detection (.pvar.zst)
    paths <- resolvePlink2Paths(stem)

    # Use pecotmr's readPvar for robust .pvar/.pvar.zst handling via pgenlibr
    vi <- readPvar(paths$pvar)
    # readPvar returns: chrom, id, pos, A2 (REF), A1 (ALT) -- pecotmr convention
    snpInfo <- tibble(
        SNP = vi$id,
        CHR = as.character(vi$chrom),
        BP = as.integer(vi$pos),
        A1 = vi$A1,
        A2 = vi$A2
    )

    # Read sample IDs from .psam
    psam <- vroom(
        paths$psam,
        delim = "\t",
        show_col_types = FALSE,
        .name_repair = .psamHeaderNames
    )
    sampleIds <- as.character(psam$IID)

    pgen <- pgenlibr::NewPgen(paths$pgen)
    nSamples <- pgenlibr::GetRawSampleCt(pgen)

    new(
        "GenotypeHandle",
        path = stem,
        format = "plink2",
        snpInfo = .withFileIdx(snpInfo),
        nSamples = as.integer(nSamples),
        sampleIds = sampleIds,
        pgenPtr = pgen
    )
}

# =============================================================================
# Block genotype extraction -- dispatches by format
# =============================================================================

#' @title Extract Block Genotypes
#' @description Extract a genotype matrix for a subset of SNPs from a
#'   \code{GenotypeHandle}. Returns a \code{RangedSummarizedExperiment} with
#'   dosage assay in Bioconductor convention (variants x samples), variant
#'   metadata as \code{rowRanges} (GRanges), and sample IDs as \code{colData}.
#' @param handle A \code{GenotypeHandle} object.
#' @param snpIdx Integer vector of 1-based SNP indices into
#'   \code{handle@@snpInfo}.
#' @param meanImpute Logical, whether to mean-impute missing values. Default
#'   TRUE.
#' @return A \code{RangedSummarizedExperiment} with:
#'   \describe{
#'     \item{assay("dosage")}{Numeric matrix (variants x samples)}
#'     \item{rowRanges}{GRanges with A1, A2 metadata}
#'     \item{colData}{DataFrame with sampleId column}
#'   }
#' @importFrom purrr map2 map set_names
#' @keywords internal
extractBlockGenotypes <- function(handle, snpIdx, meanImpute = TRUE) {
    # One-file-per-chromosome handle: route by chromosome to the right file.
    # `.genotypeChromPaths` tolerates handles deserialized before the slot
    # existed (treated as single-file).
    if (length(.genotypeChromPaths(handle)) > 0L) {
        return(.extractBlockSharded(handle, snpIdx, meanImpute = meanImpute))
    }
    # An empty request is answered without reaching a reader. This is not
    # defensive tidiness: snpStats::read.plink(select.snps = character(0))
    # SEGFAULTS, taking the R session with it rather than raising a condition
    # tryCatch could see. The sharded branch above already returns early for
    # the same reason; this is the single-file path catching up.
    if (length(snpIdx) == 0L) {
        return(.emptyBlockSe(getSampleIds(handle)))
    }
    # Read ascending, then put the columns back in the requested order. See
    # .restoreRequestedOrder() for why this is not merely tidiness.
    ord <- order(.genotypeFilePos(handle, snpIdx))
    raw <- .extractBlockByFormat(handle, snpIdx[ord])
    if (is.null(raw)) {
        return(NULL)
    }
    ordered <- .restoreRequestedOrder(raw, ord)
    geno <- if (meanImpute) .meanImputeGeno(ordered) else ordered
    .blockGenotypesToSe(geno, handle, snpIdx)
}

# Position of each requested variant in the underlying file. snpInfo rows are
# not file positions once a handle has been row-subset; fileIdx carries the
# mapping. Ordering by FILE position is what makes "read ascending" line up
# with what the ID-selecting backends actually hand back.
# @noRd
.genotypeFilePos <- function(handle, snpIdx) {
    fileIdx <- getSnpInfo(handle)$fileIdx
    if (is.null(fileIdx)) snpIdx else fileIdx[snpIdx]
}

# Put a block's variant columns back into the order the caller asked for.
#
# snpStats (`select.snps=`) and SNPRelate (`snp.id=`) return variants in FILE
# order however they were asked for, while .blockGenotypesToSe() labels the
# block from snpInfo[snpIdx], i.e. the REQUESTED order. An unsorted request
# therefore produced a block whose names and dosages described different
# variants -- silently, since both are the right length. pgenlibr reads
# positionally and was unaffected.
#
# Reading ascending and permuting back is correct for every backend, including
# the ones that do honour the request: given a sorted request they return
# sorted output either way, so the inverse permutation restores the caller's
# order in both cases.
# @noRd
.restoreRequestedOrder <- function(geno, ord) {
    if (ncol(geno) != length(ord)) {
        abort(glue(
            "extractBlockGenotypes: backend returned {ncol(geno)} variant(s) ",
            "for a request of {length(ord)}; the block cannot be labelled."
        ))
    }
    geno[, order(ord), drop = FALSE]
}

# A genotype panel as a RangedSummarizedExperiment: variants x samples,
# dosages read lazily through the handle, optional per-sample covariates as
# colData. Nothing is read here -- `genotypeDelayedArray()` only describes
# the panel. Shared by QtlDataset's genotype experiment and the LD sketch,
# which are the same object with different provenance.
# @noRd
.genotypeExperiment <- function(genotypes, genotypeCovariates = NULL) {
    if (methods::is(genotypes, "RangedSummarizedExperiment")) {
        return(.genotypeExperimentCovariates(genotypes, genotypeCovariates))
    }
    dosage <- genotypeDelayedArray(genotypes)
    # An LD panel carries no per-sample covariates; a QTL dataset's genotype
    # experiment does.
    gCov <- if (is.null(genotypeCovariates)) {
        matrix(numeric(0), nrow = 0L, ncol = 0L)
    } else {
        as.matrix(genotypeCovariates)
    }
    cd <- .genotypeColData(gCov, colnames(dosage))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = dosage),
        rowRanges = .genotypeSnpRanges(genotypes, rownames(dosage)),
        colData = cd
    )
}

# Attach covariates to a panel the caller already built. Used when a public
# entry point is handed a panel rather than a handle: rebuilding the panel from
# the handle would discard any subsetting the caller had applied, since the
# DelayedArray's seed is the whole file either way.
# @noRd
.genotypeExperimentCovariates <- function(panel, genotypeCovariates) {
    gCov <- if (is.null(genotypeCovariates)) {
        matrix(numeric(0), nrow = 0L, ncol = 0L)
    } else {
        as.matrix(genotypeCovariates)
    }
    SummarizedExperiment::`colData<-`(
        panel,
        value = .genotypeColData(gCov, colnames(panel))
    )
}

# PLINK2 writes the .psam header line with a leading "#" on the first field.
# @noRd
.psamHeaderNames <- function(nms) str_remove(nms, "^#")

# Covariate rows carry the sample identity. An unnamed matrix is only
# interpretable when it already has one row per panel sample, in order.
# @noRd
.genotypeCovariateRownames <- function(gCov, sampleIds) {
    if (!is.null(rownames(gCov))) {
        return(gCov)
    }
    if (nrow(gCov) != length(sampleIds)) {
        abort(glue(
            "'genotypeCovariates' has {nrow(gCov)} rows but the panel ",
            "has {length(sampleIds)} samples; name its rows to align ",
            "them explicitly"
        ))
    }
    `rownames<-`(gCov, sampleIds)
}

# Per-sample covariates as a colData aligned to the panel's sample order.
# Samples the covariate matrix does not name get NA rather than being
# dropped: the assay still has a column for them.
# @noRd
.genotypeColData <- function(gCov, sampleIds) {
    empty <- S4Vectors::DataFrame(row.names = sampleIds)
    if (ncol(gCov) == 0L || nrow(gCov) == 0L) {
        return(empty)
    }
    named <- .genotypeCovariateRownames(gCov, sampleIds)
    S4Vectors::DataFrame(
        named[match(sampleIds, rownames(named)), , drop = FALSE],
        row.names = sampleIds
    )
}

# snpInfo as rowRanges for a genotype panel: width-1 ranges at each
# variant's position, alleles carried as mcols so plyranges predicates can
# reach them. Deliberately the same shape as the SummarizedExperiment
# `extractBlockGenotypes()` returns -- chr-prefixed seqnames, SNP / A1 / A2 --
# so a block read from this experiment and a block read through the handle
# describe variants the same way, and overlaps between the two actually hit.
# @noRd
.genotypeSnpRanges <- function(genotypes, variantIds) {
    si <- getSnpInfo(genotypes)
    `names<-`(
        GenomicRanges::GRanges(
            seqnames = withChrPrefix(as.character(si$CHR)),
            ranges = IRanges::IRanges(as.integer(si$BP), width = 1L),
            SNP = as.character(si$SNP),
            A1 = as.character(si$A1),
            A2 = as.character(si$A2)
        ),
        variantIds
    )
}

# Dispatch block extraction to the format-specific backend (samples x variants).
# @noRd
.extractBlockByFormat <- function(handle, snpIdx) {
    fmt <- getFormat(handle)
    switch(
        fmt,
        "gds" = .extractBlockGds(handle, snpIdx),
        "vcf" = .extractBlockVcf(handle, snpIdx),
        "plink1" = .extractBlockPlink1(handle, snpIdx),
        "plink2" = .extractBlockPlink2(handle, snpIdx),
        .abortUnsupportedFormat(fmt, "extractBlockGenotypes")
    )
}

# Wrap an extracted (samples x variants) dosage block into a variants x samples
# SummarizedExperiment with per-variant rowRanges + sample colData.
# @noRd
.blockGenotypesToSe <- function(geno, handle, snpIdx) {
    si <- slice(getSnpInfo(handle), snpIdx)
    chr <- str_c(
        "chr",
        str_remove(as.character(si$CHR), regex("^chr", ignore_case = TRUE))
    )
    rowRanges <- GRanges(
        seqnames = chr,
        ranges = IRanges(start = as.integer(si$BP), width = 1L),
        SNP = si$SNP,
        A1 = si$A1,
        A2 = si$A2
    )
    sampleIds <- getSampleIds(handle)
    # Transpose to Bioc convention: variants x samples.
    dosage <- `dimnames<-`(t(geno), list(si$SNP, sampleIds))
    SummarizedExperiment(
        assays = list(dosage = dosage),
        rowRanges = rowRanges,
        colData = DataFrame(sampleId = sampleIds, row.names = sampleIds)
    )
}

# Extract one chromosome's block from a sharded handle: reslice the handle to
# that chrom's file + SNP subset and delegate to extractBlockGenotypes.
# @noRd
.extractBlockForChrom <- function(
    chrom,
    posInReq,
    handle,
    snpIdx,
    unifiedChr,
    meanImpute
) {
    if (!is_in(chrom, names(getChromPaths(handle)))) {
        msg <- glue(
            "extractBlockGenotypes: no per-chromosome file for chromosome ",
            "'{chrom}' (have: ",
            "{str_flatten(names(handle@chromPaths), collapse = ', ')})."
        )
        abort(msg)
    }
    blockGlobal <- which(unifiedChr == chrom) # file-order global indices
    localIdx <- match(snpIdx[posInReq], blockGlobal)
    th <- methods::initialize(
        handle,
        path = handle@chromPaths[[chrom]],
        snpInfo = slice(handle@snpInfo, blockGlobal),
        pgenPtr = NULL,
        # Treated as single-file: the view addresses one chromosome's payload.
        chromPaths = character(0)
    )
    extractBlockGenotypes(th, localIdx, meanImpute = meanImpute)
}

# Extract a block from a one-file-per-chromosome (sharded) handle. The global
# snpIdx index the unified @snpInfo; we group them by chromosome, route each
# group to its per-chromosome payload via a transient single-file view (with
# the chromosome-local snpInfo so positional backends like PLINK2 stay valid),
# then row-bind across chromosomes (samples identical by construction) and
# restore the requested order. A single-chromosome request -- the common cis
# case -- returns its one SE directly.
#' @keywords internal
.extractBlockSharded <- function(handle, snpIdx, meanImpute = TRUE) {
    sampleIds <- getSampleIds(handle)
    if (length(snpIdx) == 0L) {
        return(.emptyBlockSe(sampleIds))
    }
    unifiedChr <- canonChrom(getSnpInfo(handle)$CHR)
    groups <- split(seq_along(snpIdx), unifiedChr[snpIdx])
    ses <- set_names(
        map2(
            names(groups),
            groups,
            .extractBlockForChrom,
            handle = handle,
            snpIdx = snpIdx,
            unifiedChr = unifiedChr,
            meanImpute = meanImpute
        ),
        names(groups)
    )
    if (length(ses) == 1L) {
        return(ses[[1L]])
    }
    .combineShardedSes(ses, groups)
}

# Empty-block SE (no variants) carrying the handle's sample dimension.
# @noRd
.emptyBlockSe <- function(sampleIds) {
    empty <- matrix(
        numeric(0),
        nrow = 0L,
        ncol = length(sampleIds),
        dimnames = list(character(0), sampleIds)
    )
    SummarizedExperiment(
        assays = list(dosage = empty),
        rowRanges = GRanges(),
        colData = DataFrame(sampleId = sampleIds, row.names = sampleIds)
    )
}

# Combine per-chromosome SEs at the assay/rowRanges level (rbind-ing the SEs
# trips on disjoint seqlevels) and restore the requested snpIdx order.
# @noRd
.combineShardedSes <- function(ses, groups) {
    ord <- order(unname(list_c(groups)))
    dosages <- map(ses, .seDosage)
    combinedDos <- exec(rbind, !!!dosages)[ord, , drop = FALSE]
    rowRangesList <- unname(map(ses, SummarizedExperiment::rowRanges))
    combinedGr <- suppressWarnings(exec(c, !!!rowRangesList))[ord]
    SummarizedExperiment(
        assays = list(dosage = combinedDos),
        rowRanges = combinedGr,
        colData = SummarizedExperiment::colData(ses[[1L]])
    )
}

#' @keywords internal
.extractBlockGds <- function(handle, snpIdx) {
    .withGds(
        .genotypeReadPath(handle),
        .gdsBlockGeno,
        handle = handle,
        snpIdx = snpIdx,
        allow.fork = TRUE
    )
}

#' @keywords internal
.extractBlockVcf <- function(handle, snpIdx) {
    si <- slice(getSnpInfo(handle), snpIdx)
    gr <- GRanges(
        seqnames = si$CHR,
        ranges = IRanges(start = si$BP, end = si$BP)
    )
    param <- VariantAnnotation::ScanVcfParam(
        which = gr,
        fixed = NA,
        info = NA,
        geno = "GT"
    )
    vcf <- VariantAnnotation::readVcf(
        .genotypeReadPath(handle),
        genome = "",
        param = param
    )
    gt <- VariantAnnotation::geno(vcf)$GT

    # Convert GT strings to ALT dosage (A1 dosage): one column per variant,
    # laid side by side rather than filled into a preallocated matrix.
    matrix(
        map_dbl(as.vector(t(gt)), .gtStringToDosage),
        nrow = ncol(gt),
        ncol = nrow(gt)
    )
}

#' @keywords internal
.extractBlockPlink1 <- function(handle, snpIdx) {
    snpIds <- getSnpInfo(handle)$SNP[snpIdx]
    pathStem <- .genotypeReadPath(handle)
    plinkData <- snpStats::read.plink(
        bed = str_c(pathStem, ".bed"),
        bim = str_c(pathStem, ".bim"),
        fam = str_c(pathStem, ".fam"),
        select.snps = snpIds
    )
    # snpStats as(x, "numeric") gives count of B allele (A2/bim col 6).
    # Flip to count A1 (bim col 5 / effect allele).
    `storage.mode<-`(2 - as(plinkData$genotypes, "numeric"), "double")
}

#' @keywords internal
#' @importFrom rlang try_fetch
.extractBlockPlink2 <- function(handle, snpIdx) {
    # pgenlibr::ReadList returns ALT dosage = A1 dosage in pecotmr convention.
    # The cached @pgenPtr does not survive saveRDS/readRDS (external pointers
    # become stale), so we re-open from getPath() on the fly if the cached
    # pointer errors out. Opening is cheap relative to dosage extraction.
    paths <- resolvePlink2Paths(.genotypeReadPath(handle))
    # `variant_subset` indexes the .pgen by FILE position. `snpIdx` is a
    # position into @snpInfo, which may have been row-subset; translate through
    # fileIdx to recover the true .pgen index. For a full handle fileIdx ==
    # seq_len(nrow), so this is a no-op. Older RDS handles predate the column ->
    # fall back to snpIdx.
    fileIdx <- getSnpInfo(handle)$fileIdx
    variantSubset <- if (is.null(fileIdx)) snpIdx else fileIdx[snpIdx]
    # A sharded handle routes through a transient view with pgenPtr = NULL (one
    # pgen per chromosome), and a deserialized pointer is stale; open a fresh
    # pgen up front in those cases rather than provoking a caught read error.
    ptr <- getPgenPtr(handle) %||% pgenlibr::NewPgen(paths$pgen)
    geno <- try_fetch(
        pgenlibr::ReadList(
            ptr,
            variant_subset = variantSubset,
            meanimpute = FALSE
        ),
        error = function(cnd) {
            reopened <- pgenlibr::NewPgen(paths$pgen)
            pgenlibr::ReadList(
                reopened,
                variant_subset = variantSubset,
                meanimpute = FALSE
            )
        }
    )
    `storage.mode<-`(geno, "double")
}

# =============================================================================
# LD correlation computation
# =============================================================================

# Samples x variants dosage matrix for a genotype block: extract via the
# GenotypeHandle pipeline and transpose out of the Bioc variants x samples
# layout. The shared form of the t(assay(extractBlockGenotypes(...))) idiom.
# @noRd
.dosageMatrix <- function(handle, snpIdx, meanImpute = TRUE) {
    block <- extractBlockGenotypes(handle, snpIdx, meanImpute = meanImpute)
    # NULL, not an error, when the reader declines the block: callers read
    # that as "no genotypes here" (computeLd answers with an identity).
    if (is.null(block)) {
        return(NULL)
    }
    t(SummarizedExperiment::assay(block, "dosage"))
}

# Open a GDS read-only, guarantee it is closed on exit, and return fn(gds).
# The shared form of the snpgdsOpen(...) + on.exit(snpgdsClose(gds)) idiom
# used across the GDS readers.
# @noRd
# Resource bracket: open the GDS at `path`, guarantee it is closed on exit, and
# run `fn(gds, ...)`. `fn` is a top-level function (not an inline closure); its
# per-call inputs are threaded through `...`.
.withGds <- function(path, fn, ..., readonly = TRUE, allow.fork = FALSE) {
    gds <- SNPRelate::snpgdsOpen(
        path,
        readonly = readonly,
        allow.fork = allow.fork
    )
    on.exit(SNPRelate::snpgdsClose(gds))
    fn(gds, ...)
}

#' @keywords internal
.computeBlockLdGds <- function(handle, snpIdx) {
    .withGds(
        .genotypeReadPath(handle),
        .gdsBlockLd,
        handle = handle,
        snpIdx = snpIdx,
        allow.fork = TRUE
    )
}

# =============================================================================
# Region filtering helper
# =============================================================================

#' @title Filter SNP Info by Region
#' @description Return 1-based indices into snpInfo for SNPs within a genomic
#'   region string.
#' @param snpInfo data.frame with CHR and BP columns.
#' @param region Character region string "chr:start-end".
#' @return Integer vector of matching SNP indices.
#' @keywords internal
.regionToSnpIdx <- function(snpInfo, region) {
    parsed <- parseRegion(region)
    chrMatch <- stripChrPrefix(as.character(snpInfo$CHR)) == parsed$chrom
    posMatch <- snpInfo$BP >= parsed$start & snpInfo$BP <= parsed$end
    which(chrMatch & posMatch)
}

#' @title Convert SNP Info to Variant Info
#' @description Convert GenotypeHandle snpInfo (uppercase columns) to pecotmr
#'   variant_info format (lowercase columns).
#' @param snpInfo data.frame with SNP, CHR, BP, A1, A2 columns.
#' @return data.frame with chrom, id, pos, A2, A1 columns.
#' @keywords internal
.snpInfoToVariantInfo <- function(snpInfo) {
    tibble(
        chrom = snpInfo$CHR,
        id = snpInfo$SNP,
        pos = snpInfo$BP,
        A2 = snpInfo$A2,
        A1 = snpInfo$A1
    )
}

# =============================================================================
# Helpers
# =============================================================================

#' @keywords internal
.meanImputeGeno <- function(geno) {
    naMask <- is.na(geno)
    if (!any(naMask)) {
        return(geno)
    }
    # One fill over the whole matrix rather than a copy per column: geno is
    # samples x variants, so rebuilding it column by column is the expensive
    # way to say this.
    means <- colMeans(geno, na.rm = TRUE)
    replace(geno, naMask, means[col(geno)[naMask]])
}

#' @keywords internal
.gdsSnpInfo <- function(gdsPath) {
    .withGds(gdsPath, .gdsReadSnpInfo)
}

# Map a lowercase file extension to a genotype format, or NULL if unrecognized.
# @noRd
.h2FormatFromExt <- function(ext) {
    switch(
        ext,
        "vcf" = "vcf",
        "bcf" = "vcf",
        "bed" = "plink1",
        "bim" = "plink1",
        "fam" = "plink1",
        "pgen" = "plink2",
        "pvar" = "plink2",
        "psam" = "plink2",
        "gds" = "gds",
        "rds" = "rds",
        "rdata" = "rds",
        "annot" = "ldsc_annot",
        "bw" = "bigwig",
        "bigwig" = "bigwig",
        NULL
    )
}

.h2DetectFormat <- function(path) {
    lpath <- str_to_lower(path)
    if (
        str_detect(lpath, "\\.vcf\\.gz$") ||
            str_detect(lpath, "\\.vcf\\.bgz$")
    ) {
        return("vcf")
    }
    if (str_detect(lpath, "\\.annot\\.gz$")) {
        return("ldsc_annot")
    }

    ext <- str_to_lower(file_ext(path))
    if (str_length(ext) > 0L) {
        detected <- .h2FormatFromExt(ext)
        if (!is.null(detected)) return(detected)
    }
    # Check for file stems, including dotted prefixes such as sample.EUR.chr21.
    if (
        file.exists(str_c(path, ".pgen")) || file.exists(str_c(path, ".pvar"))
    ) {
        return("plink2")
    }
    if (file.exists(str_c(path, ".bed")) || file.exists(str_c(path, ".bim"))) {
        return("plink1")
    }
    if (file.exists(str_c(path, ".gds"))) {
        return("gds")
    }
    if (str_length(ext) > 0L) {
        msg <- glue("Cannot detect format from extension: {ext}")
        abort(msg)
    }
    msg <- glue("Cannot detect genotype format for path: {path}")
    abort(msg)
}

#' @title Detect Plink File Stem
#' @description Given any plink file path, return the stem.
#' @param path Character, path to any plink file.
#' @return Character, file stem without extension.
#' @keywords internal
.plinkStem <- function(path) {
    # Only strip known plink extensions; leave other paths as-is (they may
    # already be the stem, e.g. "prefix.genotype" -> "prefix.genotype.bed")
    ext <- file_ext(path)
    plinkExts <- c("bed", "bim", "fam", "pgen", "pvar", "psam")
    if (is_in(str_to_lower(ext), plinkExts)) {
        file_path_sans_ext(path)
    } else {
        path
    }
}


# =============================================================================
# Format-specific file readers
# -----------------------------------------------------------------------------
# Low-level PLINK / VCF / GDS variant-metadata readers, the stochastic
# genotype sidecar helpers (.afreq / .stochastic_meta.tsv), and the
# top-level dispatchers loadGenotypeRegion + getRefVariantInfo that
# auto-detect the underlying format and route to the correct reader.
# =============================================================================

# read PLINK files

#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
#' @importFrom Rsamtools TabixFile seqnamesTabix scanTabix headerTabix
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom SummarizedExperiment assay
readBim <- function(bed) {
    bimf <- str_c(file_path_sans_ext(bed), ".bim")
    `colnames<-`(
        vroom(bimf, col_names = FALSE),
        c("chrom", "id", "gpos", "pos", "a1", "a0")
    )
}

#' @importFrom vroom vroom
#' @importFrom tools file_path_sans_ext
readFam <- function(bed) {
    famf <- str_c(file_path_sans_ext(bed), ".fam")
    return(vroom(famf, col_names = FALSE))
}

# open bed/bim/fam: A PLINK 1 .bed is a valid .pgen
openBed <- function(bed) {
    if (!requireNamespace("pgenlibr", quietly = TRUE)) {
        abort("Package 'pgenlibr' is required for this function.")
    }
    rawSCt <- nrow(readFam(bed))
    return(pgenlibr::NewPgen(bed, raw_sample_ct = rawSCt))
}

#' Read a PLINK2 allele frequency file (.afreq or .afreq.zst)
#'
#' @param prefix File prefix (without .afreq extension).
#' @return A data.frame with columns: chrom, id, A2 (REF), A1 (ALT), alt_freq,
#'   obs_ct. alt_freq is the frequency of the A1 (ALT/effect) allele.
#' @importFrom vroom vroom
#' @importFrom archive archive_read
#' @importFrom dplyr rename select
#' @examples
#' stem <- file.path(system.file("extdata", "ld_reference", "chr22",
#'   package = "pecotmr"), "protocol_example.LD.chr22")
#' readAfreq(stem)
#' @export
#' @importFrom checkmate assertString
readAfreq <- function(prefix) {
    assertString(prefix)
    afreqZst <- str_c(prefix, ".afreq.zst")
    afreqPlain <- str_c(prefix, ".afreq")
    if (file.exists(afreqZst)) {
        af <- vroom(
            archive_read(afreqZst, format = "raw", filter = "zstd"),
            delim = "\t",
            show_col_types = FALSE
        )
    } else if (file.exists(afreqPlain)) {
        af <- vroom(
            afreqPlain,
            delim = "\t",
            show_col_types = FALSE
        )
    } else {
        return(NULL)
    }
    # PLINK2 .afreq: REF = A2, ALT = A1, ALT_FREQS = A1 (effect allele)
    # frequency
    renamed <- rename(
        af,
        "chrom" = "#CHROM",
        "id" = "ID",
        "A2" = "REF",
        "A1" = "ALT",
        "alt_freq" = "ALT_FREQS",
        "obs_ct" = "OBS_CT"
    )
    # A stochastic-genotype .afreq also carries U_MIN / U_MAX, which
    # invertMinmaxScaling() needs for an exact min-max inversion.
    stochastic <- is_in("U_MIN", colnames(renamed))
    cols <- c(
        c("chrom", "id", "A2", "A1", "alt_freq", "obs_ct"),
        if (stochastic) c("u_min", "u_max")
    )
    out <- if (stochastic) {
        rename(renamed, "u_min" = "U_MIN", "u_max" = "U_MAX")
    } else {
        renamed
    }
    select(out, all_of(cols))
}

#' Read stochastic genotype sidecar metadata (U_MIN/U_MAX).
#'
#' Reads per-variant min/max values used to invert min-max [0,2] scaling
#' of stochastic genotype data. Supports two formats:
#' \itemize{
#'   \item \strong{afreq}: PLINK2 .afreq/.afreq.zst with U_MIN/U_MAX columns
#'     (read via \code{readAfreq}, which also returns allele frequencies).
#'   \item \strong{generic}: Tab-delimited file with columns id, u_min, u_max.
#' }
#'
#' @param path Path to the sidecar metadata file.
#' @param format One of \code{NULL} (auto-detect from extension),
#'   \code{"afreq"}, or \code{"generic"}. When \code{NULL}, files ending in
#'   \code{.afreq} or \code{.afreq.zst} are parsed as afreq; all others as
#'   generic.
#' @return A data.frame with columns \code{id}, \code{u_min}, \code{u_max}, or
#'   \code{NULL} if the file lacks U_MIN/U_MAX columns (afreq format) or doesn't
#'   exist.
#' @importFrom vroom vroom
#' @noRd
readStochasticMeta <- function(path, format = NULL) {
    if (!file.exists(path)) {
        return(NULL)
    }

    detected <- format %||%
        if (str_detect(path, "\\.afreq(\\.zst)?$")) {
            "afreq"
        } else {
            "generic"
        }
    resolved <- arg_match(detected, c("afreq", "generic"))

    if (resolved == "afreq") {
        # readAfreq expects a prefix, not a full path - strip the .afreq[.zst]
        # suffix
        prefix <- str_remove(path, "\\.afreq(\\.zst)?$")
        af <- readAfreq(prefix)
        if (is.null(af) || !all(is_in(c("u_min", "u_max"), colnames(af)))) {
            return(NULL)
        }
        return(select(af, all_of(c("id", "u_min", "u_max"))))
    }

    # Generic: expect tab-delimited with columns id, u_min, u_max
    meta <- vroom(path, delim = "\t", show_col_types = FALSE)
    required <- c("id", "u_min", "u_max")
    if (!all(is_in(required, colnames(meta)))) {
        msg <- glue(
            "Stochastic metadata file '{path}' must contain columns: ",
            "{str_flatten(required, collapse = ', ')}"
        )
        abort(msg)
    }
    select(meta, all_of(required))
}

#' Search for a stochastic genotype sidecar file alongside a genotype path.
#'
#' Looks for \code{.afreq}, \code{.afreq.zst}, and \code{.stochastic_meta.tsv}
#' files next to the given genotype path. For extension-based paths (VCF, GDS),
#' the extension is stripped first. For prefix-based paths (PLINK1/2), the
#' prefix is used directly.
#'
#' @param genotypePath Path to the genotype data (prefix or file path).
#' @return Path to the first sidecar file found, or \code{NULL}.
#' @noRd
findStochasticMeta <- function(genotypePath) {
    # Strip known genotype extensions to get the stem
    stem <- str_remove(
        genotypePath,
        "\\.(vcf|vcf\\.gz|bcf|gds|bed|bim|fam|pgen|pvar|psam)$"
    )
    candidates <- c(
        str_c(stem, ".afreq"),
        str_c(stem, ".afreq.zst"),
        str_c(stem, ".stochastic_meta.tsv")
    )
    found <- candidates[file.exists(candidates)]
    if (length(found) > 0) found[1] else NULL
}


#' Invert min-max [0,2] scaling to recover the original U matrix
#'
#' Stochastic genotype data is stored after min-max scaling: U_scaled = 2 * (U -
#' u_min) / (u_max - u_min). This function exactly inverts that transform using
#' the stored per-variant u_min and u_max values from a companion sidecar file
#' (.afreq or .stochastic_meta.tsv).
#'
#' The recovered U satisfies U'U/B ~ Wishart(B, R)/B, the correct distributional
#' property for LD-based fine-mapping with dynamic variance tracking.
#'
#' @param X Numeric matrix (B x p) of min-max scaled values in [0, 2].
#' @param uMin Numeric vector of per-variant minimum values before scaling.
#' @param uMax Numeric vector of per-variant maximum values before scaling.
#' @return Matrix of original U values with same dimensions.
#' @examples
#' X <- matrix(runif(12), 4, 3)
#' invertMinmaxScaling(X, uMin = rep(0, 3), uMax = rep(1, 3))
#' @export
invertMinmaxScaling <- function(X, uMin, uMax) {
    if (length(uMin) != ncol(X) || length(uMax) != ncol(X)) {
        msg <- glue(
            "Length of u_min/u_max ({length(uMin)}) must equal ",
            "ncol(X) ({ncol(X)})"
        )
        abort(msg)
    }
    span <- uMax - uMin
    # Monomorphic variant: the forward scaling was the identity.
    denom <- replace(span, span == 0, 1)
    # Invert: U_original = U_scaled * (u_max - u_min) / 2 + u_min
    sweep(sweep(X, 2, denom / 2, "*"), 2, uMin, "+")
}

# ---------- Internal helpers for PLINK2 format ----------

#' @importFrom checkmate assertFileExists
#' Resolve and validate PLINK2 file paths for a given prefix.
#' @return Named list with pgen, pvar, psam paths.
#' @noRd
resolvePlink2Paths <- function(prefix) {
    pgen <- str_c(prefix, ".pgen")
    assertFileExists(pgen, access = "r", .var.name = "PLINK2 .pgen file")
    # Prefer plain .pvar (fast, no extra deps); fall back to .pvar.zst
    pvar <- if (file.exists(str_c(prefix, ".pvar"))) {
        str_c(prefix, ".pvar")
    } else if (file.exists(str_c(prefix, ".pvar.zst"))) {
        str_c(prefix, ".pvar.zst")
    } else {
        msg <- glue("PLINK2 .pvar[.zst] file not found at prefix: {prefix}")
        abort(msg)
    }
    psam <- str_c(prefix, ".psam")
    assertFileExists(psam, access = "r", .var.name = "PLINK2 .psam file")
    list(pgen = pgen, pvar = pvar, psam = psam)
}

#' Read .pvar or .pvar.zst into a data.frame via pgenlibr.
#'
#' Uses pgenlibr::NewPvar() to parse the file (handles both plain .pvar and
#' zstd-compressed .pvar.zst natively, no external CLI required).
#'
#' @param pvarPath Path to .pvar or .pvar.zst file.
#' @return data.frame with columns: chrom, id, pos, A2 (REF), A1 (ALT).
#' @noRd
readPvar <- function(pvarPath) {
    if (!requireNamespace("pgenlibr", quietly = TRUE)) {
        abort("Package 'pgenlibr' is required.")
    }
    pvar <- pgenlibr::NewPvar(pvarPath)
    on.exit(pgenlibr::ClosePvar(pvar), add = TRUE)
    n <- pgenlibr::GetVariantCt(pvar)
    idx <- seq_len(n)
    tibble(
        chrom = map_chr(idx, .pvarChrom, pvar = pvar),
        id = map_chr(idx, .pvarId, pvar = pvar),
        pos = map_int(idx, .pvarPos, pvar = pvar),
        A2 = map_chr(idx, .pvarAlleleCode, pvar = pvar, k = 1L),
        A1 = map_chr(idx, .pvarAlleleCode, pvar = pvar, k = 2L)
    )
}

# PLINK1 .bim files come in a 6-column form and a 9-column extended form;
# the column meanings are positional, so the width selects the header.
# @noRd
.plinkVariantColnames <- function(n, snpFilePath) {
    if (n == 6) {
        return(c("chrom", "id", "gpos", "pos", "A1", "A2"))
    }
    if (n == 9) {
        return(c(
            "chrom",
            "id",
            "gpos",
            "pos",
            "A1",
            "A2",
            "variance",
            "allele_freq",
            "n_nomiss"
        ))
    }
    abort(glue(
        "Unexpected number of columns ({n}) in variant file: ",
        "{snpFilePath}"
    ))
}

#' Read variant metadata from either .bim or .pvar/.pvar.zst file.
#'
#' Auto-detects the format by extension and header, then returns a standardized
#' data.frame. For PLINK1 .bim files, assigns column names based on the number
#' of columns (6 or 9). For PLINK2 .pvar files, delegates to \code{readPvar()}.
#'
#' @param snpFilePath Path to .bim, .pvar, or .pvar.zst file.
#' @return data.frame with at minimum columns: chrom, id, pos, A2, A1. Extended
#'   .bim files (9 columns) also include: variance, allele_freq, n_nomiss.
#' @importFrom readr read_table cols col_character
#' @noRd
readVariantMetadata <- function(snpFilePath) {
    # Extension first, then the header line: a .pvar written without its
    # extension still announces itself with a #CHROM header.
    byExtension <- str_detect(snpFilePath, "\\.(pvar|pvar\\.zst)$")
    isPvar <- byExtension ||
        str_detect(read_lines(snpFilePath, n_max = 1), "^#CHROM")

    if (isPvar) {
        readPvar(snpFilePath)
    } else {
        df <- read_table(snpFilePath, col_names = FALSE, col_types = cols())
        `names<-`(df, .plinkVariantColnames(ncol(df), snpFilePath))
    }
}

#' Get variant information from any LD reference source
#'
#' Auto-detects the source type (PLINK2, PLINK1, VCF, GDS, or pre-computed LD
#' metadata) and returns variant metadata. For PLINK2, opens only the .pvar
#' file. For PLINK1, reads only the .bim file. For VCF and GDS, loads the full
#' file and extracts variant info.
#'
#' @param source Genotype file path/prefix or LD metadata file path.
#' @param region Region of interest: "chr:start-end" string or data.frame with
#'   chrom/start/end. If NULL, returns all variants.
#' @return A data.frame with columns: chrom, id, pos, A2, A1. May also include
#'   allele_freq, variance, n_nomiss depending on source.
#'
#' @importFrom vroom vroom
#' @examples
#' meta <- system.file("extdata", "ld_reference", "ld_meta_file.tsv",
#'   package = "pecotmr")
#' getRefVariantInfo(meta, region = "chr22:10000000-19000000")
#' @export
getRefVariantInfo <- function(source, region = NULL) {
    resolved <- resolveLdSource(source)
    dataPath <- .refDataPath(resolved, region)
    if (resolved$type == "plink2") {
        info <- .refInfoPlink2(dataPath)
    } else if (resolved$type == "plink1") {
        info <- .refInfoPlink1(dataPath)
    } else if (is_in(resolved$type, c("vcf", "gds"))) {
        return(.refInfoVcfGds(dataPath, region))
    } else {
        return(.refInfoPrecomputedLd(resolved, region))
    }
    # Region filter for plink2 / plink1 (the loader-based paths self-filter).
    .refRegionFilter(info, region)
}

# Resolve the per-chromosome data path when the source carries region metadata,
# else the resolved data path.
# @noRd
.refDataPath <- function(resolved, region) {
    usesMeta <- is_in(resolved$type, c("plink2", "plink1", "vcf", "gds")) &&
        !is.null(resolved$metaPath) &&
        !is.null(region)
    if (usesMeta) {
        return(resolveGenotypePathForRegion(resolved$metaPath, region))
    }
    resolved$dataPath
}

# plink2 variant info from the .pvar, with .afreq allele frequency merged in.
# @noRd
.refInfoPlink2 <- function(dataPath) {
    paths <- resolvePlink2Paths(dataPath)
    info <- readPvar(paths$pvar)
    afreq <- readAfreq(dataPath)
    mutate(
        info,
        !!!compact(list(
            allele_freq = if (!is.null(afreq)) {
                afreq$alt_freq[match(info$id, afreq$id)]
            }
        ))
    )
}

# plink1 variant info from the .bim (col5 = A1, col6 = A2).
# @noRd
.refInfoPlink1 <- function(dataPath) {
    bim <- readBim(str_c(dataPath, ".bed"))
    tibble(
        chrom = bim$chrom,
        id = bim$id,
        pos = bim$pos,
        A2 = bim$a0,
        A1 = bim$a1
    )
}

# VCF / GDS variant info via the genotype loader (already region-filtered), with
# allele frequency computed from the dosage matrix.
# @noRd
.refInfoVcfGds <- function(dataPath, region) {
    result <- loadGenotypeRegion(
        dataPath,
        region = region,
        returnVariantInfo = TRUE
    )
    mutate(
        result$variant_info,
        allele_freq = colMeans(result$X, na.rm = TRUE) / 2
    )
}

# Pre-computed LD variant info: read the per-intersection bim/pvar metadata
# (already region-filtered by getRegionalLdMeta).
# @noRd
.refInfoPrecomputedLd <- function(resolved, region) {
    bimPaths <- getRegionalLdMeta(
        resolved$metaPath,
        region
    )$intersections$bimFilePaths
    bind_rows(map(bimPaths, .refReadBimMeta)) |>
        mutate(id = normalizeVariantId(.data$id))
}

# Read one bim/pvar metadata file into a canonical variant-info data.frame,
# carrying the optional variance / allele_freq / n_nomiss columns when present.
# @noRd
.refReadBimMeta <- function(path) {
    df <- readVariantMetadata(path)
    out <- tibble(
        chrom = df$chrom,
        id = df$id,
        pos = df$pos,
        A2 = df$A2,
        A1 = df$A1
    )
    optional <- intersect(c("variance", "allele_freq", "n_nomiss"), names(df))
    mutate(
        out,
        !!!set_names(map(optional, .refMetaColumn, df = df), optional)
    )
}

# @noRd
.refMetaColumn <- function(col, df) {
    df[[col]]
}

# Filter plink2 / plink1 variant info to the requested region.
# @noRd
.refRegionFilter <- function(info, region) {
    if (is.null(region)) {
        return(info)
    }
    parsed <- parseRegion(region)
    infoChrom <- stripChrPrefix(info$chrom)
    inRegion <- if (is.data.frame(parsed) && nrow(parsed) > 1) {
        .refMultiRegionMask(infoChrom, info$pos, parsed)
    } else {
        infoChrom == as.character(parsed$chrom) &
            info$pos >= parsed$start &
            info$pos <= parsed$end
    }
    info[inRegion, , drop = FALSE]
}

# OR-mask across a multi-row (one row per chrom) parsed region.
# @noRd
.refMultiRegionMask <- function(infoChrom, pos, parsed) {
    reduce(
        map(
            seq_len(nrow(parsed)),
            .refRegionRowMask,
            infoChrom = infoChrom,
            pos = pos,
            parsed = parsed
        ),
        `|`,
        .init = rep(FALSE, length(infoChrom))
    )
}

# Variants falling inside parsed region row `r`.
# @noRd
.refRegionRowMask <- function(r, infoChrom, pos, parsed) {
    infoChrom == as.character(parsed$chrom[r]) &
        pos >= parsed$start[r] &
        pos <= parsed$end[r]
}

#' Match variant_info against a whitelist file, returning logical index. Uses
#' parse_variant_id() from misc.R to handle all variant ID formats.
#' @importFrom vroom vroom
#' @importFrom readr read_lines
#' @noRd
matchVariantsToKeep <- function(variantInfo, keepVariantsPath) {
    keepRaw <- try_fetch(
        as.data.frame(vroom(keepVariantsPath, show_col_types = FALSE)),
        error = function(cnd) NULL
    )
    if (
        !is.null(keepRaw) &&
            is_in("chrom", names(keepRaw)) &&
            is_in("pos", names(keepRaw))
    ) {
        keepVariants <- parseVariantId(keepRaw)
    } else {
        # Fall back to reading as single-column variant IDs
        ids <- read_lines(keepVariantsPath)
        keepVariants <- parseVariantId(ids)
    }
    viChrom <- canonChrom(variantInfo$chrom)
    hasAlleles <- is_in("A1", names(keepVariants)) &&
        is_in("A2", names(keepVariants)) &&
        !any(is.na(keepVariants$A1)) &&
        !any(is.na(keepVariants$A2))
    if (hasAlleles) {
        is_in(
            str_c(
                viChrom,
                variantInfo$pos,
                variantInfo$A2,
                variantInfo$A1,
                sep = ":"
            ),
            str_c(
                keepVariants$chrom,
                keepVariants$pos,
                keepVariants$A2,
                keepVariants$A1,
                sep = ":"
            )
        )
    } else {
        is_in(
            str_c(viChrom, variantInfo$pos, sep = ":"),
            str_c(keepVariants$chrom, keepVariants$pos, sep = ":")
        )
    }
}


#' Load genotype data for a specific region
#'
#' Auto-detects PLINK2 (.pgen/.pvar[.zst]/.psam), PLINK1 (.bed/.bim/.fam), VCF
#' (.vcf/.vcf.gz/.bcf), or GDS (.gds) format and loads genotype data via
#' \code{\link{readGenotypes}} and \code{\link{extractBlockGenotypes}}. If a
#' stochastic genotype sidecar file (.afreq or .stochastic_meta.tsv) is found
#' alongside the genotype file, non-integer dosages are automatically rescaled
#' using the stored U_MIN/U_MAX values.
#'
#' @param genotype Path to the genotype data file (without extension).
#' @param region The target region in the format "chr:start-end".
#' @param keepIndel Whether to keep indel SNPs.
#' @param keepVariantsPath Path to a file listing variants to keep.
#' @param returnVariantInfo If TRUE, return a list with X (dosage matrix) and
#'   variant_info (data.frame). If FALSE (default), return only the dosage
#'   matrix.
#' @param stochasticMetaPath Optional explicit path to a stochastic genotype
#'   sidecar file. If NULL (default), auto-detected via
#'   \code{findStochasticMeta}.
#' @param stochasticMetaFormat Optional format override for the sidecar file:
#'   \code{"afreq"} or \code{"generic"}. If NULL (default), auto-detected from
#'   file extension.
#' @return If return_variant_info is FALSE, a numeric dosage matrix
#'   (rows=samples, cols=variants). If TRUE, a list with elements X and
#'   variant_info.
#'
#' @examples
#' stem <- sub("[.]bed$", "",
#'   system.file("extdata", "toy_ref.bed", package = "pecotmr"))
#' loadGenotypeRegion(genotype = stem, region = "chr22:1-100000000")
#' @export
loadGenotypeRegion <- function(
    genotype,
    region = NULL,
    keepIndel = TRUE,
    keepVariantsPath = NULL,
    returnVariantInfo = FALSE,
    stochasticMetaPath = NULL,
    stochasticMetaFormat = NULL
) {
    handle <- .loadGenoHandle(genotype)
    handleSnpInfo <- getSnpInfo(handle)
    snpIdx <- .loadGenoSnpIdx(handleSnpInfo, region)
    # Samples x variants matrix (pecotmr convention); callers handle missing.
    extracted <- list(
        X = .dosageMatrix(handle, snpIdx, meanImpute = FALSE),
        variant_info = .loadGenoAttachAfreq(
            handle,
            .snpInfoToVariantInfo(slice(handleSnpInfo, snpIdx))
        )
    )
    filtered <- .loadGenoPostFilter(extracted, keepIndel, keepVariantsPath)
    result <- .loadGenoInvertStochastic(
        filtered,
        genotype,
        stochasticMetaPath,
        stochasticMetaFormat
    )
    if (returnVariantInfo) result else result$X
}

# Detect the genotype file format and open the matching GenotypeHandle.
# @noRd
.loadGenoHandle <- function(genotype) {
    if (str_detect(genotype, "\\.(vcf|vcf\\.gz|bcf)$")) {
        return(.readGenotypeHandle(genotype, format = "vcf"))
    }
    if (str_detect(genotype, "\\.gds$")) {
        return(.readGenotypeHandle(genotype, format = "gds"))
    }
    if (hasPlink2Files(genotype)) {
        return(.readGenotypeHandle(genotype, format = "plink2"))
    }
    if (hasPlink1Files(genotype)) {
        return(.readGenotypeHandle(genotype, format = "plink1"))
    }
    msg <- glue(
        "Genotype files not found at: {genotype}\n",
        "  Expected: .vcf/.vcf.gz/.bcf, .gds, or PLINK prefix ",
        "(.pgen/.pvar[.zst]/.psam or .bed/.bim/.fam)",
        .trim = FALSE
    )
    abort(msg)
}

# Resolve the SNP indices for a region (all variants when region is NULL).
# @noRd
.loadGenoSnpIdx <- function(handleSnpInfo, region) {
    if (is.null(region)) {
        return(seq_len(nrow(handleSnpInfo)))
    }
    snpIdx <- .regionToSnpIdx(handleSnpInfo, region)
    if (length(snpIdx) == 0) {
        msg <- glue("No SNPs found in the specified region {region}")
        abort(msg, class = "NoSnpsError")
    }
    snpIdx
}

# Attach allele frequency from the .afreq sidecar (plink2 only).
#' @importFrom dplyr left_join
#' @noRd
.loadGenoAttachAfreq <- function(handle, variantInfo) {
    if (getFormat(handle) != "plink2") {
        return(variantInfo)
    }
    afreq <- readAfreq(.genotypeReadPath(handle))
    if (is.null(afreq)) {
        return(variantInfo)
    }
    afreqCols <- intersect(c("id", "alt_freq", "obs_ct"), colnames(afreq))
    # left_join keeps every variant row + the .afreq order (merge sort = FALSE).
    left_join(
        variantInfo,
        select(afreq, all_of(afreqCols)),
        by = "id"
    )
}

# Apply the indel-drop and variant-whitelist post-filters to (X, variant_info).
# @noRd
.loadGenoPostFilter <- function(result, keepIndel, keepVariantsPath) {
    snpOnly <- if (keepIndel) {
        result
    } else {
        .loadGenoSubsetVariants(
            result,
            isSnpAlleles(result$variant_info$A1, result$variant_info$A2)
        )
    }
    if (is.null(keepVariantsPath)) {
        return(snpOnly)
    }
    .loadGenoSubsetVariants(
        snpOnly,
        matchVariantsToKeep(snpOnly$variant_info, keepVariantsPath)
    )
}

# Subset a loaded (X, variant_info) pair along the variant axis, keeping the
# dosage columns and the metadata rows in step.
# @noRd
.loadGenoSubsetVariants <- function(result, keep) {
    list_assign(
        result,
        X = result$X[, keep, drop = FALSE],
        variant_info = result$variant_info[keep, , drop = FALSE]
    )
}

# Detect stochastic genotype scaling and restore the original scale from a
# metadata sidecar; warn (only) when non-integer dosages have no sidecar.
# @noRd
.loadGenoInvertStochastic <- function(
    result,
    genotype,
    stochasticMetaPath,
    stochasticMetaFormat
) {
    metaPath <- stochasticMetaPath %||% findStochasticMeta(genotype)
    if (is.null(metaPath)) {
        .loadGenoWarnStochastic(result$X)
        return(result)
    }
    smeta <- readStochasticMeta(metaPath, format = stochasticMetaFormat)
    if (is.null(smeta)) {
        return(result)
    }
    idx <- match(colnames(result$X), smeta$id)
    matched <- !is.na(idx)
    if (!any(matched)) {
        return(result)
    }
    inverted <- `[<-`(
        result$X,
        ,
        matched,
        value = invertMinmaxScaling(
            result$X[, matched, drop = FALSE],
            smeta$u_min[idx[matched]],
            smeta$u_max[idx[matched]]
        )
    )
    restored <- list_assign(
        result,
        X = inverted,
        variant_info = mutate(
            result$variant_info,
            u_min = smeta$u_min[idx],
            u_max = smeta$u_max[idx]
        )
    )
    msg <- glue(
        "Stochastic genotype detected: restored original scale via ",
        "{basename(metaPath)}"
    )
    inform(msg)
    restored
}

# Warn when non-integer dosages are present but no stochastic sidecar was found.
# @noRd
.loadGenoWarnStochastic <- function(X) {
    if (all(X == round(X), na.rm = TRUE)) {
        return(invisible(NULL))
    }
    msg <- glue(
        "Non-integer genotype values detected but no stochastic metadata ",
        "sidecar found. Place a .afreq or .stochastic_meta.tsv file with ",
        "u_min/u_max columns alongside the genotype files to restore the ",
        "original scale."
    )
    warn(msg)
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# The first ALT allele of a variant (as character).
# @noRd
.gtFirstAllele <- function(x) {
    as.character(x)[1]
}

# The dosage assay matrix of one per-chromosome SummarizedExperiment.
# @noRd
.seDosage <- function(se) {
    SummarizedExperiment::assay(se, "dosage")
}

# ALT (A1) dosage from a VCF GT string; NA for missing ("./." or NA).
# @noRd
.gtStringToDosage <- function(x) {
    if (is.na(x) || x == "./.") {
        return(NA_real_)
    }
    alleles <- str_split(x, "[/|]")[[1L]]
    sum(alleles != "0")
}

# Variant `i`'s chromosome from an open pvar handle.
# @noRd
.pvarChrom <- function(i, pvar) {
    pgenlibr::GetVariantChrom(pvar, i)
}

# Variant `i`'s id from an open pvar handle.
# @noRd
.pvarId <- function(i, pvar) {
    pgenlibr::GetVariantId(pvar, i)
}

# Variant `i`'s base-pair position from an open pvar handle.
# @noRd
.pvarPos <- function(i, pvar) {
    pgenlibr::GetVariantPos(pvar, i)
}

# Allele code `k` (1 = A2/ref, 2 = A1/alt) of variant `i` from a pvar handle.
# @noRd
.pvarAlleleCode <- function(i, pvar, k) {
    pgenlibr::GetAlleleCode(pvar, i, k)
}

# ---- .withGds resource-body functions (run inside the open-GDS bracket) ---

# The sample ids stored in an open GDS.
# @noRd
.gdsReadSampleIds <- function(gds) {
    as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.id")))
}

# The dosage matrix for `handle`'s `snpIdx` variants (NULL when none selected).
# snpgdsGetGeno counts the first allele (labelled A1 in .gdsSnpInfo), so no flip
# is needed; it also handles non-contiguous SNP selection.
# @noRd
.gdsBlockGeno <- function(gds, handle, snpIdx) {
    snpIds <- getSnpInfo(handle)$SNP[snpIdx]
    geno <- SNPRelate::snpgdsGetGeno(
        gds,
        snp.id = snpIds,
        with.id = FALSE,
        verbose = FALSE
    )
    if (is.null(geno) || length(geno) == 0) {
        return(NULL)
    }
    `storage.mode<-`(geno, "double")
}

# The sample-LD matrix for `handle`'s `snpIdx` variants (NA correlations -> 0).
# @noRd
.gdsBlockLd <- function(gds, handle, snpIdx) {
    # snpgdsLDMat() returns variants in FILE order whatever order snp.id is
    # given in, and labels nothing. Ask in ascending order and permute back,
    # exactly as the block readers do, so the matrix agrees with its names.
    ord <- order(.genotypeFilePos(handle, snpIdx))
    snpIds <- getSnpInfo(handle)$SNP[snpIdx[ord]]
    ldMat <- SNPRelate::snpgdsLDMat(
        gds,
        snp.id = snpIds,
        method = "corr",
        slide = -1,
        verbose = FALSE
    )
    raw <- replace(ldMat$LD, is.na(ldMat$LD), 0)
    inv <- order(ord)
    ids <- getSnpInfo(handle)$SNP[snpIdx]
    `dimnames<-`(raw[inv, inv, drop = FALSE], list(ids, ids))
}

# The (SNP, CHR, BP, A1, A2) snpInfo frame from an open GDS. A1 = the first
# snp.allele (the allele snpgdsGetGeno counts), so dosage = count of A1.
# @noRd
.gdsReadSnpInfo <- function(gds) {
    snpId <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.id"))
    chr <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.chromosome"))
    pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.position"))
    allele <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.allele"))
    allelesSplit <- str_split(allele, "/")
    # `[` (single-bracket) so a missing second allele yields NA, not an error.
    a1 <- map_chr(allelesSplit, `[`, 1L)
    a2 <- map_chr(allelesSplit, `[`, 2L)
    tibble(
        SNP = snpId,
        CHR = as.character(chr),
        BP = as.integer(pos),
        A1 = a1,
        A2 = a2
    )
}
