#' @title Genotype I/O
#' @description Read genotype data in its own format (VCF, plink1, plink2,
#'   GDS, BGEN) as a lazily read genotype panel, without converting it:
#'   PLINK2 and BGEN files through a GenotypeHandle, VCF, GDS and PLINK1
#'   files through a GenotypeDosageSeed.
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
    function(path, format = NULL, vcfArgs = list()) {
        .readGenotypePanel(path, format = format, vcfArgs = vcfArgs)
    }
)

#' @rdname readGenotypes
#' @export
setMethod(
    "readGenotypes",
    signature(path = "missing"),
    function(
        path,
        plink1Prefix = NULL,
        plink2Prefix = NULL,
        bed = NULL,
        bim = NULL,
        fam = NULL,
        pgen = NULL,
        pvar = NULL,
        psam = NULL,
        ldMeta = NULL,
        block = NULL,
        genoMeta = NULL,
        chroms = NULL,
        format = NULL,
        vcfArgs = list()
    ) {
        # The keyword sources are declared here rather than absorbed by `...`
        # and forwarded, so a misspelled source is an error at the call
        # rather than a silently unset one.
        .genotypePanelFromSources(
            plink1Prefix = plink1Prefix,
            plink2Prefix = plink2Prefix,
            bed = bed,
            bim = bim,
            fam = fam,
            pgen = pgen,
            pvar = pvar,
            psam = psam,
            ldMeta = ldMeta,
            block = block,
            genoMeta = genoMeta,
            chroms = chroms,
            format = format,
            vcfArgs = vcfArgs
        )
    }
)

# The panel a keyword source names: exactly one source, checked as
# GenotypeHandle() checks it. A PLINK2 triplet is read through a handle, a
# PLINK1 one as a panel; a file an LD-meta row or a per-chromosome meta
# names, through whatever reads its format.
# @noRd
.genotypePanelFromSources <- function(
    plink1Prefix = NULL,
    plink2Prefix = NULL,
    bed = NULL,
    bim = NULL,
    fam = NULL,
    pgen = NULL,
    pvar = NULL,
    psam = NULL,
    ldMeta = NULL,
    block = NULL,
    genoMeta = NULL,
    chroms = NULL,
    format = NULL,
    vcfArgs = list()
) {
    flags <- .ghValidateArgs(
        bed = bed,
        bim = bim,
        fam = fam,
        pgen = pgen,
        pvar = pvar,
        psam = psam,
        ldMeta = ldMeta,
        block = block
    )
    sources <- .ghResolveSources(
        path = NULL,
        plink1Prefix = plink1Prefix,
        plink2Prefix = plink2Prefix,
        ldMeta = ldMeta,
        genoMeta = genoMeta,
        chroms = chroms,
        flags = flags
    )
    switch(
        names(sources)[which(sources)],
        plink1Prefix = .readGenotypePanel(plink1Prefix, format = "plink1"),
        plink2Prefix = .genotypeExperiment(.makePlink2Handle(plink2Prefix)),
        plink1Triplet = .readGenotypePanel(
            .plink1TripletStem(bed, bim, fam),
            format = "plink1"
        ),
        plink2Triplet = .genotypeExperiment(
            .genotypeHandleFromPlink2Triplet(pgen, pvar, psam)
        ),
        ldMeta = .ghLdPathToPanel(.ghResolveLdPath(ldMeta, block), vcfArgs),
        genoMeta = .chromMetaPanel(genoMeta, chroms, format)
    )
}

# The formats read as a panel directly, through a GenotypeDosageSeed; every
# other format is read through a GenotypeHandle.
# @noRd
.DOSAGE_PANEL_FORMATS <- c("vcf", "gds", "plink1")

# The panel of one genotype file or PLINK prefix. A VCF, GDS or PLINK1 file
# is read as a panel directly; a PLINK2 or BGEN file through a
# GenotypeHandle.
# @noRd
.readGenotypePanel <- function(path, format = NULL, vcfArgs = list()) {
    format <- format %||% .requireGenotypeFormat(path, "readGenotypes")
    if (!is_in(format, .DOSAGE_PANEL_FORMATS)) {
        return(.genotypeExperiment(.readGenotypeHandle(path, format = format)))
    }
    .withSourceProjection(.readDosagePanel(path, format, vcfArgs))
}

# The handle behind a panel of a format a GenotypeHandle reads: PLINK2 or
# BGEN. Internal: the handle is the seed layer, and the
# handle-construction machinery below calls this rather than readGenotypes()
# so it does not wrap and immediately unwrap a panel on every hop.
# @noRd
.readGenotypeHandle <- function(path, format = NULL) {
    format <- format %||% .requireGenotypeFormat(path, "readGenotypes")
    switch(
        format,
        "plink2" = .makePlink2Handle(path),
        "bgen" = .makeBgenHandle(path),
        "vcf" = ,
        "gds" = ,
        "plink1" = abort(c(
            glue("A {format} file is read as a genotype panel, not a handle."),
            i = "Read it with readGenotypes()."
        )),
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
# sharded routing read fileIdx[snpIdx]. For a full (unsubset) handle fileIdx
# == seq_len(nrow), so reads are unchanged.
# @noRd
.withFileIdx <- function(snpInfo) {
    mutate(snpInfo, fileIdx = seq_len(nrow(snpInfo)))
}

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
    # via .emptyLdRef (the empty / PIP-skip path), never on a surviving region.
    methods::initialize(
        handle,
        snpInfo = slice(snpInfo(handle), integer(0)),
        sampleIds = character(0),
        nSamples = 0L
    )
}

# Restrict a GenotypeHandle's @snpInfo to `keep` (a logical mask or integer row
# indices into @snpInfo). Genotype reads stay correct because fileIdx carries
# each kept variant's original file position; everything else (path/format/
# pgenPtr/chromPaths/samples) is preserved. NULL-safe; a no-op when nothing is
# dropped. Handles built before the fileIdx column existed are NOT subset (the
# read path would be positional) -- return them unchanged.
# @noRd
.subsetGenotypeHandle <- function(handle, keep) {
    if (is.null(handle)) {
        return(NULL)
    }
    si <- snpInfo(handle)
    keepIdx <- if (is.logical(keep)) which(keep) else as.integer(keep)
    if (length(keepIdx) >= nrow(si)) {
        return(handle)
    } # nothing dropped
    if (!is_in("fileIdx", names(si))) {
        return(handle)
    } # legacy handle: unsafe
    methods::initialize(handle, snpInfo = slice(si, keepIdx))
}

#' @keywords internal
.makePlink2Handle <- function(path) {
    if (!requireNamespace("pgenlibr", quietly = TRUE)) {
        abort("Package 'pgenlibr' is required for reading plink2 files.")
    }

    stem <- .genotypeStem(path)

    # Use pecotmr's resolvePlink2Paths for robust path detection (.pvar.zst)
    paths <- resolvePlink2Paths(stem)

    # Use pecotmr's readPvar for robust .pvar/.pvar.zst handling via pgenlibr
    vi <- readPvar(paths$pvar)
    # readPvar returns: chrom, id, pos, A2 (REF), A1 (ALT) -- pecotmr convention
    snpInfo <- tibble(
        SNP = vi$id,
        CHR = withChrPrefix(vi$chrom),
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
    # An empty request is answered without reaching a reader, which need not
    # accept one, as the sharded branch above already answers it.
    if (length(snpIdx) == 0L) {
        return(.emptyBlockSe(sampleIds(handle)))
    }
    # Read ascending, then put the columns back in the requested order. See
    # .restoreRequestedOrder() for why this is not merely tidiness.
    ord <- order(.genotypeFilePos(handle, snpIdx))
    raw <- .extractBlockByFormat(handle, snpIdx[ord])
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
    fileIdx <- snpInfo(handle)$fileIdx
    if (is.null(fileIdx)) snpIdx else fileIdx[snpIdx]
}

# Put a block's variant columns back into the order the caller asked for.
#
# .blockGenotypesToSe() labels the block from snpInfo[snpIdx], i.e. the
# REQUESTED order, so a backend returning variants in FILE order however they
# were asked for -- snpStats did, when it read PLINK1 -- would yield a block
# whose names and dosages describe different variants, silently, since both
# are the right length. Reading ascending and permuting back is correct for
# every backend, whether or not it honours the requested order: given a
# sorted request it returns sorted output either way, so the inverse
# permutation restores the caller's order.
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
    # A stochastic sketch's dosages are projections, so its real sample count
    # and allele frequencies are recorded from its sidecar here, once; and
    # where the genotypes come from is recorded for what needs the file
    # (.ldRefSource()).
    .withLdProjection(
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = .genotypeSnpRanges(genotypes, rownames(dosage)),
            colData = cd,
            metadata = list(genotypeSource = .genotypeSourceRecord(genotypes))
        ),
        .sketchProjectionRecord(genotypes)
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

# =============================================================================
# Panels saved by an earlier version
# -----------------------------------------------------------------------------
# A panel built today records where its genotypes come from
# (metadata$genotypeSource), and one over a VCF, GDS or PLINK1 file reads
# through that format's own backend. A panel saved before either reads
# through a GenotypeHandle's GhSeed, in a format a handle may no longer read,
# and a stochastic sketch saved then lacks the real sample count and allele
# frequencies it is read with today (.withLdProjection()). Such a panel is
# read again from its files, as readGenotypes() reads them, and cut to the
# variants and samples it held. Only the files' variant and sample lists and
# a sketch's sidecar are read, never genotypes, and nothing is written.
# =============================================================================

# `panel` as a panel built today: read again from its files when it was
# saved before the changes above, else `panel` itself -- as for one that
# reads through no handle (a VCF, GDS or PLINK1 panel built since, or one
# held in memory). A bare GenotypeHandle, which an LD reference was before
# it was a panel, is taken as the panel over it.
# @noRd
.updatePanel <- function(panel) {
    panel <- .asLdRef(panel)
    handle <- .ldRefHandle(panel)
    if (is.null(handle)) {
        return(panel)
    }
    source <- .genotypeSourceRecord(handle)
    recorded <- !is.null(S4Vectors::metadata(panel)$genotypeSource)
    if (recorded && !is_in(source$format, .DOSAGE_PANEL_FORMATS)) {
        return(panel)
    }
    .rereadPanel(panel, source)
}

# `panel` read again from `source`, the record of where its genotypes come
# from, and cut to the variants and samples it held, keeping its sample data
# (a QTL dataset's genotype covariates). Variants are matched on position
# and alleles, which identify them however their names were normalized when
# the panel was saved. A `pecotmrUnreadablePanel` error when the files
# cannot be read or no longer hold them all.
# @noRd
#' @importFrom rlang try_fetch
.rereadPanel <- function(panel, source) {
    fresh <- try_fetch(
        .readSourcePanel(source),
        error = function(cnd) {
            .abortUnreadablePanel(source, "they cannot be read", cnd)
        }
    )
    rows <- match(
        .grVariantIds(SummarizedExperiment::rowRanges(panel)),
        .grVariantIds(SummarizedExperiment::rowRanges(fresh))
    )
    cols <- match(colnames(panel), colnames(fresh))
    if (length(rows) != nrow(panel) || anyNA(rows) || anyNA(cols)) {
        .abortUnreadablePanel(
            source,
            "they no longer hold every variant and sample it did"
        )
    }
    SummarizedExperiment::`colData<-`(
        .panelCut(fresh, rows, cols),
        value = SummarizedExperiment::colData(panel)
    )
}

# The panel today's readers build over a genotype source: its per-chromosome
# files' panels stacked when it has them, else its one file's. A bundled
# example's portable path resolves under the installed package.
# @noRd
.readSourcePanel <- function(source) {
    if (length(source$chromPaths) > 0L) {
        return(.chromMetaPanel(
            map_chr(source$chromPaths, .resolveGenotypeResourcePath),
            format = source$format
        ))
    }
    .readGenotypePanel(
        .resolveGenotypeResourcePath(source$path),
        format = source$format
    )
}

# `panel` cut to rows `rows` and columns `cols`, each cut only when it does
# not keep everything in order.
# @noRd
.panelCut <- function(panel, rows, cols) {
    if (!identical(rows, seq_len(nrow(panel)))) {
        panel <- .ldRefSubset(panel, rows)
    }
    if (identical(cols, seq_len(ncol(panel)))) {
        return(panel)
    }
    panel[, cols]
}

# @noRd
.abortUnreadablePanel <- function(source, why, parent = NULL) {
    abort(
        c(
            glue(
                "A {source$format} genotype panel saved by an earlier pecotmr ",
                "has to be read again from its files at {source$path}, and ",
                "{why}."
            ),
            i = "Put the files back where they were, then update again."
        ),
        class = "pecotmrUnreadablePanel",
        parent = parent
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
    .snpInfoRanges(snpInfo(genotypes), variantIds)
}

# A panel's variants snpInfo-style, one row per rowRanges element: SNP (the
# file's own ids), CHR, BP, A1, A2. The inverse of .snpInfoRanges(), so
# genotypes from any backend answer alike.
# @noRd
.rangesSnpInfo <- function(gr) {
    mc <- S4Vectors::mcols(gr)
    tibble(
        SNP = as.character(mc$SNP),
        CHR = as.character(seqnames(gr)),
        BP = as.integer(start(gr)),
        A1 = as.character(mc$A1),
        A2 = as.character(mc$A2)
    )
}

# The rowRanges of variants described snpInfo-style (SNP, CHR, BP, A1, A2),
# named `variantIds`: a handle's, or a VCF, GDS or PLINK1 panel's rows. Every
# reader spells CHR in UCSC style as it reads the file.
# @noRd
.snpInfoRanges <- function(si, variantIds) {
    `names<-`(
        GenomicRanges::GRanges(
            seqnames = as.character(si$CHR),
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
    fmt <- genotypeFormat(handle)
    switch(
        fmt,
        "plink2" = .extractBlockPlink2(handle, snpIdx),
        "bgen" = .extractBlockBgen(handle, snpIdx),
        .abortUnsupportedFormat(fmt, "extractBlockGenotypes")
    )
}

# Wrap an extracted (samples x variants) dosage block into a variants x samples
# SummarizedExperiment with per-variant rowRanges + sample colData.
# @noRd
.blockGenotypesToSe <- function(geno, handle, snpIdx) {
    si <- slice(snpInfo(handle), snpIdx)
    rowRanges <- GRanges(
        seqnames = si$CHR,
        ranges = IRanges(start = as.integer(si$BP), width = 1L),
        SNP = si$SNP,
        A1 = si$A1,
        A2 = si$A2
    )
    sampleIds <- sampleIds(handle)
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
    if (!is_in(chrom, names(chromPaths(handle)))) {
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
    sampleIds <- sampleIds(handle)
    if (length(snpIdx) == 0L) {
        return(.emptyBlockSe(sampleIds))
    }
    unifiedChr <- snpInfo(handle)$CHR
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
#' @importFrom rlang try_fetch
.extractBlockPlink2 <- function(handle, snpIdx) {
    # pgenlibr::ReadList returns ALT dosage = A1 dosage in pecotmr convention.
    # The cached @pgenPtr does not survive saveRDS/readRDS (external pointers
    # become stale), so we re-open from path() on the fly if the cached
    # pointer errors out. Opening is cheap relative to dosage extraction.
    paths <- resolvePlink2Paths(.genotypeReadPath(handle))
    # `variant_subset` indexes the .pgen by FILE position. `snpIdx` is a
    # position into @snpInfo, which may have been row-subset; translate through
    # fileIdx to recover the true .pgen index. For a full handle fileIdx ==
    # seq_len(nrow), so this is a no-op. Older RDS handles predate the column ->
    # fall back to snpIdx.
    fileIdx <- snpInfo(handle)$fileIdx
    variantSubset <- if (is.null(fileIdx)) snpIdx else fileIdx[snpIdx]
    # A sharded handle routes through a transient view with pgenPtr = NULL (one
    # pgen per chromosome), and a deserialized pointer is stale; open a fresh
    # pgen up front in those cases rather than provoking a caught read error.
    ptr <- pgenPtr(handle) %||% pgenlibr::NewPgen(paths$pgen)
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
    t(SummarizedExperiment::assay(block, "dosage"))
}

# Resource bracket: open the SNPRelate GDS at `path` read-only, guarantee it
# is closed on exit, and run `fn(gds, ...)` -- the snpgdsOpen(...) +
# on.exit(snpgdsClose(gds)) idiom of the on-disk LD routine, which needs
# SNPRelate's own handle. `fn` is a top-level function (not an inline
# closure); its per-call inputs are threaded through `...`. A file a panel's
# GDSArray backend already holds open in this session is opened again rather
# than refused.
# @noRd
.withGds <- function(
    path,
    fn,
    fnArgs = list(),
    readonly = TRUE,
    allow.fork = FALSE
) {
    gds <- SNPRelate::snpgdsOpen(
        path,
        readonly = readonly,
        allow.duplicate = TRUE,
        allow.fork = allow.fork
    )
    on.exit(SNPRelate::snpgdsClose(gds))
    exec(fn, gds, !!!fnArgs)
}

#' @keywords internal
.computeBlockLdGds <- function(path, snpIds) {
    .withGds(
        path,
        .gdsBlockLd,
        fnArgs = list(snpIds = snpIds),
        allow.fork = TRUE
    )
}

# =============================================================================
# Region filtering helper
# =============================================================================

#' @title Filter SNP Info by Region
#' @description Return 1-based indices into snpInfo for SNPs within a genomic
#'   window.
#' @param snpInfo data.frame with CHR and BP columns.
#' @param region A \code{GRanges} window, of one range or several.
#' @return Integer vector of matching SNP indices.
#' @keywords internal
.regionToSnpIdx <- function(snpInfo, region) {
    which(.refMultiRegionMask(
        as.character(snpInfo$CHR),
        snpInfo$BP,
        .regionTable(region)
    ))
}

# Internal: a GRanges window argument as the chrom / start / end table the
# genotype readers match variants against, one row per range, its
# chromosome in UCSC style as the LD and genotype metadata are.
# @noRd
.regionTable <- function(region, what = "`region`") {
    region <- .assertRegionGRanges(region, what)
    tibble(
        chrom = as.character(GenomicRanges::seqnames(region)),
        start = GenomicRanges::start(region),
        end = GenomicRanges::end(region)
    )
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

# =============================================================================
# Genotype format detection
# -----------------------------------------------------------------------------
# The one place a genotype path's format is told. A file is told by its
# extension, in any case; a PLINK prefix by the PLINK file beside it, a .pgen
# before a .bed (PLINK2 can hold dosages, so it wins when both are there). A
# VCF, GDS or BGEN file is always named by its own path.
# =============================================================================

# The extensions naming each genotype format, as regular expressions.
# @noRd
.GENOTYPE_EXTENSIONS <- c(
    vcf = "vcf(\\.b?gz)?|bcf",
    gds = "gds",
    bgen = "bgen",
    plink1 = "bed|bim|fam",
    plink2 = "pgen|pvar(\\.zst)?|psam"
)

# The genotype format a file's extension names, or NULL when it names none.
# @noRd
.genotypeFileFormat <- function(path) {
    named <- str_detect(
        str_to_lower(path),
        str_c("\\.(", .GENOTYPE_EXTENSIONS, ")$")
    )
    if (any(named)) names(.GENOTYPE_EXTENSIONS)[which(named)[[1L]]] else NULL
}

# The genotype format of a file or PLINK prefix: by the file's extension, or
# for a prefix by the PLINK file beside it; NULL when neither tells.
# @noRd
.genotypePathFormat <- function(path) {
    byExtension <- .genotypeFileFormat(path)
    if (!is.null(byExtension)) {
        return(byExtension)
    }
    if (file.exists(str_c(path, ".pgen"))) {
        return("plink2")
    }
    if (file.exists(str_c(path, ".bed"))) {
        return("plink1")
    }
    NULL
}

# .genotypePathFormat(), or an error from `label` when the path names no
# genotype format.
# @noRd
.requireGenotypeFormat <- function(path, label) {
    format <- .genotypePathFormat(path)
    if (!is.null(format)) {
        return(format)
    }
    abort(c(
        glue(
            "{label}: '{path}' is neither a genotype file nor a PLINK prefix."
        ),
        i = str_c(
            "A genotype file is named by its extension (.vcf.gz, .vcf.bgz, ",
            ".gds, .bgen, .bed, .pgen); a PLINK prefix by the .pgen or .bed ",
            "beside it."
        )
    ))
}

# The stem of a genotype path: the path without the extension naming its
# format, which is what PLINK's readers and a file's sidecars (.afreq,
# .stochastic_meta.tsv, .sample) are named after. A PLINK prefix is its own
# stem.
# @noRd
.genotypeStem <- function(path) {
    str_remove(
        path,
        regex(
            str_c("\\.(", str_flatten(.GENOTYPE_EXTENSIONS, "|"), ")$"),
            ignore_case = TRUE
        )
    )
}

# =============================================================================
# Format-specific file readers
# -----------------------------------------------------------------------------
# Low-level PLINK / VCF / GDS variant-metadata readers, the stochastic
# genotype sidecar helpers (.afreq / .stochastic_meta.tsv), and the
# top-level dispatchers loadGenotypeRegion + refVariantInfo that
# auto-detect the underlying format and route to the correct reader.
# =============================================================================

# read PLINK files

#' @importFrom vroom vroom
#' @importFrom Rsamtools TabixFile seqnamesTabix scanTabix headerTabix
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom SummarizedExperiment assay
readBim <- function(bed) {
    bimf <- str_c(.genotypeStem(bed), ".bim")
    `colnames<-`(
        vroom(bimf, col_names = FALSE),
        c("chrom", "id", "gpos", "pos", "a1", "a0")
    )
}

#' Read a PLINK2 allele frequency file (.afreq or .afreq.zst)
#'
#' @param path File prefix (without .afreq extension).
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
readAfreq <- function(path) {
    assertString(path)
    afreqZst <- str_c(path, ".afreq.zst")
    afreqPlain <- str_c(path, ".afreq")
    if (file.exists(afreqZst)) {
        return(.readAfreqFile(afreqZst))
    }
    if (file.exists(afreqPlain)) {
        return(.readAfreqFile(afreqPlain))
    }
    NULL
}

# One PLINK2 .afreq table read from the file named, zstd-compressed when the
# name ends in .zst. readAfreq() finds that file from a prefix; a caller that
# holds the file itself (readStochasticMeta()'s format override) reads it
# here, so the file need not be named *.afreq.
# @noRd
.readAfreqFile <- function(file) {
    af <- if (str_detect(file, "\\.zst$")) {
        vroom(
            archive_read(file, format = "raw", filter = "zstd"),
            delim = "\t",
            show_col_types = FALSE
        )
    } else {
        vroom(file, delim = "\t", show_col_types = FALSE)
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
#'   \item \strong{afreq}: a PLINK2 .afreq table with U_MIN/U_MAX columns,
#'     whatever the file is named (zstd-compressed when it ends in .zst).
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
        # The file itself, under whatever name it has: going through
        # readAfreq()'s prefix would read <path>.afreq instead.
        af <- .readAfreqFile(path)
        if (!all(is_in(c("u_min", "u_max"), colnames(af)))) {
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
#' files next to the given genotype path. For extension-based paths (VCF, GDS,
#' BGEN, or a PLINK file), the extension is stripped first. For prefix-based
#' paths (PLINK1/2), the prefix is used directly.
#'
#' @param genotypePath Path to the genotype data (prefix or file path).
#' @return Path to the first sidecar file found, or \code{NULL}.
#' @noRd
findStochasticMeta <- function(genotypePath) {
    # Strip known genotype extensions to get the stem
    stem <- .genotypeStem(genotypePath)
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
    # Monomorphic variant: the forward step divided by 1 rather than 0, which
    # left the column all 0.
    denom <- replace(span, span == 0, 1)
    # Invert: U_original = U_scaled * (u_max - u_min) / 2 + u_min
    sweep(sweep(X, 2, denom / 2, "*"), 2, uMin, "+")
}

# ---------- Internal helpers for PLINK2 format ----------

#' Resolve and validate PLINK2 file paths for a given prefix.
#' @return Named list with pgen, pvar, psam paths.
#' @importFrom checkmate assertFileExists
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

#' Get variant information from any LD reference source
#'
#' Auto-detects the source type (PLINK2, PLINK1, VCF, GDS or BGEN, given
#' directly or through LD metadata whose rows point at genotype files) and
#' returns variant metadata. For PLINK2, opens only the .pvar file. For PLINK1,
#' reads only the .bim file. For VCF, GDS and BGEN, loads the full file and
#' extracts variant info.
#'
#' @param source Genotype file path/prefix or LD metadata file path.
#' @param which Optional \code{GRanges}: return only the variants that
#'   overlap it. If NULL, returns all variants.
#' @return A data.frame with columns: chrom, id, pos, A2, A1. May also include
#'   allele_freq, variance, n_nomiss depending on source.
#'
#' @importFrom vroom vroom
#' @examples
#' meta <- system.file("extdata", "ld_reference", "ld_meta_file.tsv",
#'   package = "pecotmr")
#' refVariantInfo(meta,
#'   which = GenomicRanges::GRanges("chr22:10000000-19000000"))
#' @export
refVariantInfo <- function(source, which = NULL) {
    resolved <- resolveLdSource(source)
    dataPath <- .refDataPath(resolved, which)
    if (resolved$type == "plink2") {
        info <- .refInfoPlink2(dataPath)
    } else if (resolved$type == "plink1") {
        info <- .refInfoPlink1(dataPath)
    } else {
        return(.refInfoVcfGds(dataPath, which))
    }
    # Region filter for plink2 / plink1 (the loader-based paths self-filter).
    .refRegionFilter(info, which)
}

# Resolve the per-chromosome data path when the source carries region metadata,
# else the resolved data path.
# @noRd
.refDataPath <- function(resolved, region) {
    usesMeta <- is_in(
        resolved$type,
        c("plink2", "plink1", "vcf", "gds", "bgen")
    ) &&
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
        chrom = withChrPrefix(.data$chrom),
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
        chrom = withChrPrefix(bim$chrom),
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
        which = region,
        returnVariantInfo = TRUE
    )
    mutate(
        result$variant_info,
        allele_freq = colMeans(result$X, na.rm = TRUE) / 2
    )
}

# Filter plink2 / plink1 variant info to the requested region.
# @noRd
.refRegionFilter <- function(info, region) {
    if (is.null(region)) {
        return(info)
    }
    inRegion <- .refMultiRegionMask(
        info$chrom,
        info$pos,
        .regionTable(region)
    )
    info[inRegion, , drop = FALSE]
}

# OR-mask across a region table's rows (.regionTable()), one per range.
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
    viChrom <- variantInfo$chrom
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
#' (.vcf.gz/.vcf.bgz, bgzipped with a tabix index), GDS (.gds) or BGEN (.bgen)
#' format, reads the file as a panel with \code{\link{readGenotypes}} and
#' loads the region's dosages. If a stochastic genotype sidecar file (.afreq
#' or .stochastic_meta.tsv) is found alongside the genotype file, non-integer
#' dosages are automatically rescaled using the stored U_MIN/U_MAX values.
#'
#' @param genotype Path to the genotype data file, or a PLINK prefix.
#' @param which Optional \code{GRanges}: load only the variants that overlap
#'   it. If NULL, loads every variant.
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
#' loadGenotypeRegion(genotype = stem,
#'   which = GenomicRanges::GRanges("chr22:1-100000000"))
#' @export
loadGenotypeRegion <- function(
    genotype,
    which = NULL,
    keepIndel = TRUE,
    keepVariantsPath = NULL,
    returnVariantInfo = FALSE,
    stochasticMetaPath = NULL,
    stochasticMetaFormat = NULL
) {
    panel <- .readGenotypePanel(
        genotype,
        format = .requireGenotypeFormat(genotype, "loadGenotypeRegion")
    )
    variants <- .rangesSnpInfo(SummarizedExperiment::rowRanges(panel))
    snpIdx <- .loadGenoSnpIdx(variants, which)
    # Samples x variants matrix (pecotmr convention), its columns named by
    # the file's own variant ids; callers handle missing.
    extracted <- list(
        X = `colnames<-`(
            .ldRefDosage(panel, snpIdx, meanImpute = FALSE),
            variants$SNP[snpIdx]
        ),
        variant_info = .loadGenoAttachAfreq(
            panel,
            .snpInfoToVariantInfo(slice(variants, snpIdx))
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

# Resolve the SNP indices for a region (all variants when region is NULL).
# @noRd
.loadGenoSnpIdx <- function(variants, region) {
    if (is.null(region)) {
        return(seq_len(nrow(variants)))
    }
    snpIdx <- .regionToSnpIdx(variants, region)
    if (length(snpIdx) == 0) {
        msg <- glue(
            "No SNPs found in the specified region ",
            "{str_flatten(as.character(region), ', ')}"
        )
        abort(msg, class = "NoSnpsError")
    }
    snpIdx
}

# Attach allele frequency from the .afreq sidecar (plink2 only).
#' @importFrom dplyr left_join
#' @noRd
.loadGenoAttachAfreq <- function(panel, variantInfo) {
    if (!identical(.ldRefFormat(panel), "plink2")) {
        return(variantInfo)
    }
    afreq <- readAfreq(.ldRefReadPath(panel))
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

# Restore a stochastic sketch's original scale from its metadata sidecar.
# Without one the dosages are taken as they are: a dosage need not be a whole
# number (imputed BGEN, VCF FORMAT/DS, PLINK2 dosages), so whether the values
# are whole says nothing about whether the file is a sketch.
# @noRd
.loadGenoInvertStochastic <- function(
    result,
    genotype,
    stochasticMetaPath,
    stochasticMetaFormat
) {
    metaPath <- stochasticMetaPath %||% findStochasticMeta(genotype)
    if (is.null(metaPath)) {
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

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# The dosage assay matrix of one per-chromosome SummarizedExperiment.
# @noRd
.seDosage <- function(se) {
    SummarizedExperiment::assay(se, "dosage")
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

# ---- GDS resource-body functions (run inside an open-GDS bracket) ---------

# The sample ids stored in an open GDS.
# @noRd
.gdsReadSampleIds <- function(gds) {
    as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "sample.id")))
}

# The sample-LD matrix of the SNPs `snpIds` (NA correlations -> 0).
# @noRd
.gdsBlockLd <- function(gds, snpIds) {
    # snpgdsLDMat() returns variants in FILE order whatever order snp.id is
    # given in, and labels nothing. Ask in file order and permute back, so the
    # matrix agrees with its names; the file's own ids are passed, since a
    # GDS may store them as integers.
    fileIds <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.id"))
    pos <- match(as.character(snpIds), as.character(fileIds))
    ord <- order(pos)
    ldMat <- SNPRelate::snpgdsLDMat(
        gds,
        snp.id = fileIds[pos[ord]],
        method = "corr",
        slide = -1,
        verbose = FALSE
    )
    raw <- replace(ldMat$LD, is.na(ldMat$LD), 0)
    inv <- order(ord)
    ids <- as.character(snpIds)
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
        CHR = withChrPrefix(chr),
        BP = as.integer(pos),
        A1 = a1,
        A2 = a2
    )
}

# =============================================================================
# BGEN genotype reading
# -----------------------------------------------------------------------------
# A "bgen" GenotypeHandle reads a BGEN file lazily, without converting it:
# its variants come from the .bgi index (bgenix's SQLite database), its
# samples from the file's sample block or a .sample file beside it, and a
# block's dosages are read record by record at the offsets the index gives.
#
# Every revision is read: layout 1 (BGEN 1.1: diploid, biallelic, three
# 16-bit probabilities per sample, zlib or no compression) and layout 2 (BGEN
# 1.2 and 1.3: phased or unphased, any ploidy, any bit depth, zlib, zstd or
# no compression). A record with more than two alleles becomes one variant
# per non-first allele. The dosage is the expected count of A1, the variant's
# second (or later) allele, with the first as A2, as bgen's own tools count
# it: from 0 to the sample's ploidy.
# =============================================================================

# The handle of a BGEN file. Only its header, its index and its sample ids
# are read here.
# @noRd
#' @importFrom checkmate assertFileExists
#' @importFrom purrr walk
.makeBgenHandle <- function(path) {
    walk(c("DBI", "RSQLite"), .bgenRequirePackage)
    assertFileExists(path, access = "r", .var.name = "BGEN file")
    path <- normalizePath(path)
    header <- .bgenReadHeader(path, withSamples = TRUE)
    sampleIds <- header$sampleIds %||% .bgenSampleFileIds(path)
    if (length(sampleIds) != header$nSamples) {
        abort(glue(
            "BGEN file {path} declares {header$nSamples} samples but ",
            "{length(sampleIds)} sample ids were found."
        ))
    }
    new(
        "GenotypeHandle",
        path = path,
        format = "bgen",
        snpInfo = .withFileIdx(.bgenSnpInfo(path)),
        nSamples = as.integer(header$nSamples),
        sampleIds = sampleIds,
        pgenPtr = NULL
    )
}

# @noRd
.bgenRequirePackage <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        abort(glue(
            "Package '{pkg}' is required for reading BGEN files, whose .bgi ",
            "index is an SQLite database."
        ))
    }
    invisible(NULL)
}

# A little-endian unsigned integer from raw bytes, as a double (a uint32 can
# exceed R's integer range).
# @noRd
.bgenUint <- function(bytes) {
    sum(as.numeric(bytes) * 256^(seq_along(bytes) - 1L))
}

# The header block: the sample count, the compression and layout flags and,
# when the file stores them and they are asked for, the sample ids. An
# unknown layout or compression is refused here, before any variant is
# read.
# @noRd
.bgenReadHeader <- function(path, withSamples = FALSE) {
    con <- file(path, "rb")
    on.exit(close(con))
    fixed <- readBin(con, "raw", 20L)
    if (length(fixed) < 20L) {
        abort(glue("{path} is too short to be a BGEN file."))
    }
    headerLength <- .bgenUint(fixed[5:8])
    magic <- fixed[17:20]
    if (!identical(magic, charToRaw("bgen")) && any(magic != as.raw(0L))) {
        abort(glue("{path} is not a BGEN file (no 'bgen' magic number)."))
    }
    readBin(con, "raw", headerLength - 20)
    flags <- as.integer(readBin(con, "raw", 4L))
    compression <- bitwAnd(flags[[1L]], 3L)
    layout <- bitwAnd(bitwShiftR(flags[[1L]], 2L), 15L)
    .bgenCheckFlags(path, compression, layout)
    nSamples <- .bgenUint(fixed[13:16])
    hasIds <- bitwAnd(flags[[4L]], 128L) > 0L
    list(
        nSamples = nSamples,
        compression = compression,
        layout = layout,
        sampleIds = if (withSamples && hasIds) {
            .bgenReadSampleBlock(con, nSamples)
        }
    )
}

# @noRd
.bgenCheckFlags <- function(path, compression, layout) {
    if (!is_in(layout, 1:2)) {
        abort(glue(
            "{path} has BGEN layout {layout}; the spec defines layouts 1 and 2."
        ))
    }
    if (compression > 2L) {
        abort(glue("{path} has an unknown BGEN compression ({compression})."))
    }
    if (layout == 1L && compression == 2L) {
        abort(glue(
            "{path} is a layout-1 (BGEN 1.1) file compressed with zstd, which ",
            "the spec does not allow."
        ))
    }
    invisible(NULL)
}

# The sample identifier block, which follows the header: its length, the
# sample count again, then each id as a 2-byte length and its bytes.
# @noRd
.bgenReadSampleBlock <- function(con, nSamples) {
    blockLength <- .bgenUint(readBin(con, "raw", 4L))
    block <- readBin(con, "raw", blockLength - 4)
    if (.bgenUint(block[1:4]) != nSamples) {
        abort("The BGEN sample block's sample count disagrees with its header.")
    }
    bgenSampleIdsCpp(block[-(1:4)], as.integer(nSamples))
}

# The sample ids of a BGEN file that does not store them, from the .sample
# file beside it (`<stem>.sample`): the first column, after the header line
# and the line of column types.
# @noRd
#' @importFrom readr read_table cols col_character
.bgenSampleFileIds <- function(path) {
    samplePath <- str_c(.genotypeStem(path), ".sample")
    if (!file.exists(samplePath)) {
        abort(c(
            glue("{path} stores no sample ids, and {samplePath} is missing."),
            i = "Put the Oxford-format .sample file beside the .bgen file."
        ))
    }
    tbl <- read_table(
        samplePath,
        col_types = cols(.default = col_character()),
        progress = FALSE
    )
    as.character(tbl[[1L]][-1L])
}

# The variant table: one row per (record, non-first allele), in file order,
# with each record's offset and length from the .bgi index and the allele the
# row reads. Built column-wise for the biallelic records, which are nearly
# all of them; only a multi-allelic record is read to learn its alleles.
# @noRd
#' @importFrom dplyr arrange bind_rows
#' @importFrom purrr pmap
.bgenSnpInfo <- function(path) {
    index <- .bgenReadIndex(str_c(path, ".bgi"))
    multi <- index$number_of_alleles > 2L
    rows <- .bgenBiallelicRows(index[!multi, ])
    if (!any(multi)) {
        return(rows)
    }
    warn(glue(
        "{sum(multi)} BGEN record(s) of {path} carry more than two alleles; ",
        "each is read as one biallelic variant per non-first allele."
    ))
    split <- pmap(
        list(
            chrom = index$chromosome[multi],
            pos = index$position[multi],
            alleles = .bgenRecordAlleles(
                path,
                index$file_start_position[multi],
                index$size_in_bytes[multi]
            ),
            offset = index$file_start_position[multi],
            size = index$size_in_bytes[multi]
        ),
        .bgenSplitRecordRows
    )
    arrange(
        bind_rows(rows, list_rbind(split)),
        .data$bgen_offset,
        .data$bgen_allele
    )
}

# The rows of the biallelic records: A2 the first allele, A1 the second. A
# record keeps its rsid as its id when it has one, else is named
# chr:pos:A2:A1.
# @noRd
.bgenBiallelicRows <- function(index) {
    ids <- formatVariantId(
        index$chromosome,
        index$position,
        index$allele1,
        index$allele2
    )
    named <- !is.na(index$rsid) & nzchar(index$rsid) & index$rsid != "."
    tibble(
        SNP = if_else(named, index$rsid, ids),
        CHR = withChrPrefix(index$chromosome),
        BP = as.integer(index$position),
        A1 = index$allele2,
        A2 = index$allele1,
        bgen_offset = as.numeric(index$file_start_position),
        bgen_size = as.numeric(index$size_in_bytes),
        bgen_allele = rep(2L, nrow(index))
    )
}

# A multi-allelic record's rows: A2 its first allele, A1 each other allele
# in turn, every row named chr:pos:A2:A1 so they stay distinct.
# @noRd
.bgenSplitRecordRows <- function(chrom, pos, alleles, offset, size) {
    others <- alleles[-1L]
    tibble(
        SNP = formatVariantId(chrom, pos, alleles[[1L]], others),
        CHR = withChrPrefix(chrom),
        BP = as.integer(pos),
        A1 = others,
        A2 = alleles[[1L]],
        bgen_offset = as.numeric(offset),
        bgen_size = as.numeric(size),
        bgen_allele = seq_along(others) + 1L
    )
}

# The .bgi index's Variant table, in file order.
# @noRd
.bgenReadIndex <- function(bgiPath) {
    if (!file.exists(bgiPath)) {
        abort(c(
            str_c(
                "A BGEN file is read through its .bgi index: ",
                bgiPath,
                " is missing."
            ),
            i = "Build it with bgenix (`bgenix -g <file>.bgen -index`)."
        ))
    }
    con <- DBI::dbConnect(
        RSQLite::SQLite(),
        bgiPath,
        flags = RSQLite::SQLITE_RO
    )
    on.exit(DBI::dbDisconnect(con))
    as_tibble(DBI::dbGetQuery(
        con,
        str_c(
            "SELECT chromosome, position, rsid, number_of_alleles, allele1, ",
            "allele2, file_start_position, size_in_bytes FROM Variant ",
            "ORDER BY file_start_position"
        )
    ))
}

# The alleles of the records at `offsets`, read from their identifying data:
# the .bgi index holds only the first two.
# @noRd
.bgenRecordAlleles <- function(path, offsets, sizes) {
    con <- file(path, "rb")
    on.exit(close(con))
    map2(offsets, sizes, .bgenReadRecordAlleles, con = con)
}

# Only a layout-2 record can carry more than two alleles, so only one is
# parsed here.
# @noRd
.bgenReadRecordAlleles <- function(offset, size, con) {
    record <- .bgenReadRecord(con, offset, size)
    .bgenParseIdentifying(record, layout = 2L)$alleles
}

# @noRd
.bgenReadRecord <- function(con, offset, size) {
    seek(con, where = offset)
    readBin(con, "raw", size)
}

# A record's identifying data: its alleles, and where its genotype block
# starts (1-based). A layout-1 record opens with its sample count and has two
# alleles without saying so; a layout-2 one states its allele count.
# @noRd
.bgenParseIdentifying <- function(record, layout) {
    at <- if (layout == 1L) 5L else 1L
    for (field in 1:3) {
        # variant id, rsid, chromosome: a 2-byte length and the bytes
        at <- at + 2L + .bgenUint(record[at:(at + 1L)])
    }
    at <- at + 4L
    nAlleles <- 2L
    if (layout == 2L) {
        nAlleles <- .bgenUint(record[at:(at + 1L)])
        at <- at + 2L
    }
    alleles <- character(nAlleles)
    for (k in seq_len(nAlleles)) {
        len <- .bgenUint(record[at:(at + 3L)])
        alleles[[k]] <- rawToChar(record[(at + 4L):(at + 3L + len)])
        at <- at + 4L + len
    }
    list(alleles = alleles, genotypeStart = at)
}

# The samples x variants dosages of the handle's variants at `snpIdx`: each
# record read once, its probabilities decompressed and decoded for the
# alleles its rows ask about.
# @noRd
.extractBlockBgen <- function(handle, snpIdx) {
    si <- slice(snpInfo(handle), snpIdx)
    path <- .genotypeReadPath(handle)
    header <- .bgenReadHeader(path)
    con <- file(path, "rb")
    on.exit(close(con))
    offsets <- unique(si$bgen_offset)
    blocks <- map(
        offsets,
        .bgenRecordDosage,
        si = si,
        con = con,
        header = header,
        nSamples = handle@nSamples
    )
    out <- exec(cbind, !!!blocks)
    # The blocks hold the rows grouped by record; put them back in `si`'s
    # order.
    rowOrder <- list_c(map(offsets, .bgenRowsOf, si = si))
    out[, order(rowOrder), drop = FALSE]
}

# @noRd
.bgenRowsOf <- function(offset, si) {
    which(si$bgen_offset == offset)
}

# One record's dosage columns, for the rows of `si` that read it.
# @noRd
.bgenRecordDosage <- function(offset, si, con, header, nSamples) {
    rows <- .bgenRowsOf(offset, si)
    record <- .bgenReadRecord(con, offset, si$bgen_size[[rows[[1L]]]])
    if (header$layout == 1L) {
        return(.bgenLayout1Dosage(record, header$compression, nSamples))
    }
    probabilities <- .bgenProbabilityData(record, header$compression)
    bgenDosageCpp(
        probabilities,
        as.integer(nSamples),
        as.integer(si$bgen_allele[rows])
    )
}

# A layout-1 record's dosage column: per sample three 16-bit probabilities,
# of the genotypes AA, AB and BB, each out of 32768, zlib-compressed or
# stored as they are. The second allele's expected count is P(AB) + 2 P(BB);
# a sample whose three are all zero is missing.
# @noRd
.bgenLayout1Dosage <- function(record, compression, nSamples) {
    at <- .bgenParseIdentifying(record, layout = 1L)$genotypeStart
    nBytes <- 6 * nSamples
    data <- if (compression == 0L) {
        record[at:(at + nBytes - 1)]
    } else {
        stored <- .bgenUint(record[at:(at + 3L)])
        memDecompress(record[(at + 4L):(at + 3L + stored)], type = "gzip")
    }
    if (length(data) != nBytes) {
        abort("A layout-1 BGEN record holds the wrong amount of data.")
    }
    p <- matrix(
        readBin(
            data,
            "integer",
            n = 3 * nSamples,
            size = 2L,
            signed = FALSE,
            endian = "little"
        ) /
            32768,
        nrow = 3L
    )
    dosage <- p[2L, ] + 2 * p[3L, ]
    matrix(replace(dosage, colSums(p) == 0, NA_real_), ncol = 1L)
}

# A layout-2 record's genotype probability data, uncompressed.
# @noRd
.bgenProbabilityData <- function(record, compression) {
    at <- .bgenParseIdentifying(record, layout = 2L)$genotypeStart
    stored <- .bgenUint(record[at:(at + 3L)])
    if (compression == 0L) {
        return(record[(at + 4L):(at + 3L + stored)])
    }
    uncompressed <- .bgenUint(record[(at + 4L):(at + 7L)])
    payload <- record[(at + 8L):(at + 3L + stored)]
    .bgenUncompress(payload, compression, uncompressed)
}

# zlib through base R; zstd through archive, since R need not be built with
# it.
# @noRd
#' @importFrom archive archive_read
.bgenUncompress <- function(payload, compression, uncompressed) {
    if (compression == 1L) {
        return(memDecompress(payload, type = "gzip"))
    }
    con <- archive_read(rawConnection(payload), format = "raw", filter = "zstd")
    open(con, "rb")
    on.exit(close(con))
    readBin(con, "raw", uncompressed)
}

# =============================================================================
# VCF, GDS and PLINK1 genotype panels
# -----------------------------------------------------------------------------
# A VCF, GDS or PLINK1 file is read as a genotype panel directly rather than
# through a GenotypeHandle: its variants and samples are read here, its
# dosages through a GenotypeDosageSeed, a block at a time, when something
# touches them. The
# panel records where it came from (metadata(panel)$genotypeSource), which is
# what anything needing the file reads -- a stochastic sketch's .afreq
# sidecar, the per-chromosome files of a manifest.
#
# A record with several ALT alleles becomes one variant per ALT allele, with a
# warning, its dosage the count of that allele, as for BGEN; a record with no
# ALT allele holds no variant. The dosage is FORMAT/DS where the file declares
# one, else the count from the GT calls. A SNPRelate GDS stores genotype codes
# counting the first of a SNP's two alleles, which is therefore A1, as
# SNPRelate reads it. A PLINK1 .bed's codes are read through PlinkMatrix as
# the count of A1, the .bim's fifth column.
# =============================================================================

# @noRd
.requireReaderPackage <- function(pkg, format) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        abort(glue("Package '{pkg}' is required for reading {format} files."))
    }
    invisible(NULL)
}

# The panel of a VCF, GDS or PLINK1 file, without the record of a stochastic
# sketch's sidecar (.withSourceProjection() adds it, once the panel is whole).
# @noRd
.readDosagePanel <- function(path, format, vcfArgs = list()) {
    switch(
        format,
        vcf = .vcfPanel(path, vcfArgs),
        gds = .gdsPanel(path),
        plink1 = .plinkBedPanel(path)
    )
}

# A panel over a GenotypeDosageSeed: its rows' variants as rowRanges, shaped
# as .genotypeExperiment() shapes a handle's, and where it came from.
# @noRd
.dosageSeedPanel <- function(seed, rows, format, path) {
    dosage <- DelayedArray::DelayedArray(seed)
    SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = dosage),
        rowRanges = .snpInfoRanges(rows, rownames(dosage)),
        colData = S4Vectors::DataFrame(row.names = colnames(dosage)),
        metadata = list(
            genotypeSource = list(
                format = format,
                path = path,
                chromPaths = character(0)
            )
        )
    )
}

# One row per (record, ALT allele), in file order: A2 the record's REF
# allele, A1 the ALT allele the row counts. A record keeps its own id when it
# has one and a single ALT allele; otherwise each of its rows is named
# chr:pos:A2:A1, so they stay distinct.
# @noRd
.recordVariantRows <- function(records, path, label) {
    nAlt <- lengths(records$alt)
    multi <- nAlt > 1L
    if (any(multi)) {
        warn(glue(
            "{sum(multi)} {label} record(s) of {path} carry more than one ",
            "ALT allele; each is read as one biallelic variant per ALT allele."
        ))
    }
    record <- rep(seq_len(nrow(records)), nAlt)
    a1 <- as.character(list_c(records$alt))
    a2 <- records$ref[record]
    named <- !is.na(records$id) & !multi
    tibble(
        SNP = if_else(
            named[record],
            records$id[record],
            formatVariantId(records$chrom[record], records$pos[record], a2, a1)
        ),
        CHR = withChrPrefix(records$chrom[record]),
        BP = as.integer(records$pos[record]),
        A1 = a1,
        A2 = a2,
        record = record,
        allele = sequence(nAlt)
    )
}

# ---- VCF -------------------------------------------------------------------

# A VCF's panel. Only its header and its records' fixed fields are read here.
# @noRd
.vcfPanel <- function(path, vcfArgs = list()) {
    .requireReaderPackage("VariantAnnotation", "VCF")
    file <- .vcfIndexedFile(path)
    header <- VariantAnnotation::scanVcfHeader(file)
    field <- .vcfDosageField(header, path(file))
    records <- .vcfRecords(file, vcfArgs)
    rows <- .recordVariantRows(records, path(file), "VCF")
    seed <- methods::new(
        "GenotypeDosageSeed",
        backend = file,
        field = if (field == "DS") "vcfDs" else "vcfGt",
        rowRecord = rows$record,
        rowAllele = rows$allele,
        records = GRanges(
            seqnames = records$chrom,
            ranges = IRanges(records$pos, records$end),
            posRank = .positionRank(records$chrom, records$pos)
        ),
        dn = list(
            normalizeVariantId(rows$SNP),
            as.character(VariantAnnotation::samples(header))
        )
    )
    .dosageSeedPanel(seed, rows, "vcf", path(file))
}

# The VCF at `path` with its tabix index, through which its records are read
# a block at a time. A VCF without one cannot be read lazily; a BCF cannot be
# read at all, Rsamtools having its BCF support disabled.
# @noRd
#' @importFrom checkmate assertFileExists
.vcfIndexedFile <- function(path) {
    assertFileExists(path, access = "r", .var.name = "VCF file")
    path <- normalizePath(path)
    if (str_ends(str_to_lower(path), "\\.bcf")) {
        abort(c(
            glue("{path} is a BCF file, which pecotmr cannot read."),
            i = str_c(
                "Rsamtools ships with its BCF support disabled; convert it ",
                "to a bgzipped, tabix-indexed VCF."
            )
        ))
    }
    index <- str_c(path, ".tbi")
    if (!file.exists(index)) {
        abort(c(
            glue(
                "A VCF is read a block at a time through its tabix index, ",
                "and {index} is missing."
            ),
            i = str_c(
                "Compress the VCF with bgzip and index it, e.g. with ",
                "Rsamtools::bgzip() and Rsamtools::indexTabix(format = 'vcf')."
            )
        ))
    }
    VariantAnnotation::VcfFile(path, index = index)
}

# The FORMAT field the dosages are read from: DS when the header declares
# it, with one value per ALT allele (Number=A) or one per record (Number=1),
# else GT.
# @noRd
.vcfDosageField <- function(header, path) {
    declared <- VariantAnnotation::geno(header)
    fields <- rownames(declared)
    if (is_in("DS", fields)) {
        number <- declared$Number[match("DS", fields)]
        if (!is_in(number, c("A", "1"))) {
            abort(glue(
                "{path} declares FORMAT/DS with Number={number}; a dosage ",
                "has one value per ALT allele (Number=A) or one per record ",
                "(Number=1)."
            ))
        }
        return("DS")
    }
    if (is_in("GT", fields)) {
        return("GT")
    }
    abort(glue(
        "{path} declares neither FORMAT/DS nor FORMAT/GT, so it holds no ",
        "genotypes to read."
    ))
}

# One row per record of the VCF, in file order: its position and the end of
# its REF allele, its own ID (NA where it has none), its REF allele and its
# ALT alleles.
# @noRd
.vcfRecords <- function(file, vcfArgs) {
    vcf <- exec(
        VariantAnnotation::readVcf,
        file,
        param = VariantAnnotation::ScanVcfParam(
            fixed = "ALT",
            info = NA,
            geno = NA
        ),
        !!!vcfArgs
    )
    rr <- SummarizedExperiment::rowRanges(vcf)
    chrom <- as.character(seqnames(rr))
    pos <- as.integer(start(rr))
    ref <- as.character(rr$REF)
    alt <- as.list(methods::as(rr$ALT, "CharacterList"))
    tibble(
        chrom = chrom,
        pos = pos,
        end = as.integer(end(rr)),
        id = .vcfRecordIds(names(rr), chrom, pos, ref, alt),
        ref = ref,
        alt = map(alt, .vcfAltAlleles)
    )
}

# Each record's own ID, or NA where it has none. VariantAnnotation names a
# record whose ID is "." from the line's own text, CHROM:POS_REF/ALT, and
# that name is no ID of the record's.
# @noRd
.vcfRecordIds <- function(ids, chrom, pos, ref, alt) {
    if (is.null(ids)) {
        return(rep(NA_character_, length(chrom)))
    }
    altText <- map_chr(alt, str_flatten, collapse = ",")
    synthetic <- str_c(
        chrom,
        ":",
        pos,
        "_",
        ref,
        "/",
        if_else(altText == "", ".", altText)
    )
    if_else(ids == synthetic | ids == ".", NA_character_, ids)
}

# A record's ALT alleles; a record whose ALT is "." has none.
# @noRd
.vcfAltAlleles <- function(alt) {
    alt[!is_in(alt, c("", "."))]
}

# Each record's rank among the records starting at its position, in file
# order: 1 for nearly all of them.
# @noRd
.positionRank <- function(chrom, pos) {
    tibble(chrom = chrom, pos = pos) |>
        mutate(rank = row_number(), .by = all_of(c("chrom", "pos"))) |>
        pull("rank")
}

# ---- GDS -------------------------------------------------------------------

# A GDS file's panel: a SNPRelate (SNP_ARRAY) or SeqArray (SEQ_ARRAY) file,
# told apart by the format its root node records.
# @noRd
#' @importFrom checkmate assertFileExists
.gdsPanel <- function(path) {
    walk(c("gdsfmt", "GDSArray"), .requireReaderPackage, format = "GDS")
    assertFileExists(path, access = "r", .var.name = "GDS file")
    path <- normalizePath(path)
    kind <- .withGdsFile(path, .gdsFileFormat)
    switch(
        kind,
        SNP_ARRAY = .snpGdsPanel(path),
        SEQ_ARRAY = .seqArrayPanel(path),
        abort(glue(
            "{path} is a GDS file of format '{kind}'; genotypes are read ",
            "from SNPRelate (SNP_ARRAY) and SeqArray (SEQ_ARRAY) files."
        ))
    )
}

# Resource bracket: open the GDS at `path` read-only with gdsfmt, close it on
# exit, and run `fn(gds, ...)`. A file GDSArray already holds open in this
# session is opened again rather than refused.
# @noRd
.withGdsFile <- function(path, fn, fnArgs = list()) {
    gds <- gdsfmt::openfn.gds(path, readonly = TRUE, allow.duplicate = TRUE)
    on.exit(gdsfmt::closefn.gds(gds))
    exec(fn, gds, !!!fnArgs)
}

# The kind of GDS file an open GDS is, as its root node records it.
# @noRd
.gdsFileFormat <- function(gds) {
    as.character(gdsfmt::get.attr.gdsn(gds$root)$FileFormat %||% "unknown")
}

# A node of an open GDS, read whole.
# @noRd
.gdsRead <- function(gds, name) {
    gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, name))
}

# A SNPRelate GDS's panel.
# @noRd
.snpGdsPanel <- function(path) {
    meta <- .withGdsFile(path, .snpGdsMeta)
    rows <- mutate(
        meta$snpInfo,
        SNP = as.character(.data$SNP),
        record = row_number(),
        allele = 1L
    )
    backend <- DelayedArray::seed(GDSArray::GDSArray(path, "genotype"))
    stored <- c(nrow(rows), length(meta$sampleIds))
    expected <- if (meta$snpFirst) stored else rev(stored)
    if (!identical(as.integer(dim(backend)), as.integer(expected))) {
        abort(glue(
            "The genotype matrix of {path} is ",
            "{str_flatten(dim(backend), ' x ')}, but the file lists ",
            "{nrow(rows)} SNPs and {length(meta$sampleIds)} samples."
        ))
    }
    seed <- methods::new(
        "GenotypeDosageSeed",
        backend = backend,
        field = if (meta$snpFirst) "snpGdsSnpOrder" else "snpGdsSampleOrder",
        rowRecord = rows$record,
        rowAllele = rows$allele,
        records = IRanges(rows$record, width = 1L),
        dn = list(normalizeVariantId(rows$SNP), meta$sampleIds)
    )
    .dosageSeedPanel(seed, rows, "gds", path)
}

# An open SNPRelate GDS's SNPs, its samples, and how its genotype matrix is
# stored: SNP by SNP unless its node says sample.order, as SNPRelate reads
# it.
# @noRd
.snpGdsMeta <- function(gds) {
    node <- gdsfmt::index.gdsn(gds, "genotype")
    list(
        snpInfo = .gdsReadSnpInfo(gds),
        sampleIds = .gdsReadSampleIds(gds),
        snpFirst = !is_in("sample.order", names(gdsfmt::get.attr.gdsn(node)))
    )
}

# A SeqArray GDS's panel.
# @noRd
.seqArrayPanel <- function(path) {
    meta <- .withGdsFile(path, .seqArrayMeta, fnArgs = list(path = path))
    rows <- .recordVariantRows(meta$records, path, "SeqArray GDS")
    counts <- as.integer(meta$counts)
    seed <- methods::new(
        "GenotypeDosageSeed",
        backend = DelayedArray::seed(GDSArray::GDSArray(path, meta$node)),
        field = meta$field,
        rowRecord = rows$record,
        rowAllele = rows$allele,
        records = IRanges(start = cumsum(counts) - counts + 1L, width = counts),
        dn = list(normalizeVariantId(rows$SNP), meta$sampleIds)
    )
    .dosageSeedPanel(seed, rows, "gds", path)
}

# An open SeqArray GDS's samples and records, and the node its dosages are
# read from, with how many of that node's columns (or 2-bit planes) each
# record takes.
# @noRd
.seqArrayMeta <- function(gds, path) {
    alleles <- str_split(as.character(.gdsRead(gds, "allele")), ",")
    c(
        list(
            sampleIds = as.character(.gdsRead(gds, "sample.id")),
            records = tibble(
                chrom = as.character(.gdsRead(gds, "chromosome")),
                pos = as.integer(.gdsRead(gds, "position")),
                id = .seqArrayIds(gds, length(alleles)),
                ref = map_chr(alleles, 1L),
                alt = map(alleles, .seqArrayAltAlleles)
            )
        ),
        .seqArrayDosageNode(gds, path)
    )
}

# Each record's own ID, or NA where it has none or the file stores none.
# @noRd
.seqArrayIds <- function(gds, nRecords) {
    node <- gdsfmt::index.gdsn(gds, "annotation/id", silent = TRUE)
    if (is.null(node)) {
        return(rep(NA_character_, nRecords))
    }
    ids <- as.character(gdsfmt::read.gdsn(node))
    if_else(is_in(ids, c("", ".")), NA_character_, ids)
}

# A record's ALT alleles: every allele after its first, the REF.
# @noRd
.seqArrayAltAlleles <- function(alleles) {
    .vcfAltAlleles(alleles[-1L])
}

# Where a SeqArray GDS's dosages are: its FORMAT/DS when it stores one, with
# one value per ALT allele or one per record, else its genotype calls.
# @noRd
.seqArrayDosageNode <- function(gds, path) {
    ds <- gdsfmt::index.gdsn(gds, "annotation/format/DS", silent = TRUE)
    if (!is.null(ds)) {
        number <- as.character(gdsfmt::get.attr.gdsn(ds)$Number %||% "")
        if (!is_in(number, c("A", "1"))) {
            abort(glue(
                "{path} stores FORMAT/DS with Number={number}; a dosage has ",
                "one value per ALT allele (Number=A) or one per record ",
                "(Number=1)."
            ))
        }
        return(list(
            node = "annotation/format/DS/data",
            field = "seqArrayDs",
            counts = .gdsRead(gds, "annotation/format/DS/@data")
        ))
    }
    .seqArrayCheckGenotypes(gds, path)
    list(
        node = "genotype/data",
        field = "seqArrayGt",
        counts = .gdsRead(gds, "genotype/@data")
    )
}

# The genotype calls are read from their 2-bit planes, so a file storing
# them otherwise, or storing calls beyond its ploidy apart from them, is
# refused rather than misread.
# @noRd
.seqArrayCheckGenotypes <- function(gds, path) {
    data <- gdsfmt::objdesp.gdsn(gdsfmt::index.gdsn(gds, "genotype/data"))
    if (!identical(data$trait, "Bit2")) {
        abort(glue(
            "{path} stores its genotype calls as {data$trait}; SeqArray ",
            "stores them in 2-bit planes, which is what is read."
        ))
    }
    extra <- gdsfmt::index.gdsn(gds, "genotype/extra.index", silent = TRUE)
    nExtra <- if (is.null(extra)) 0L else gdsfmt::objdesp.gdsn(extra)$dim[[2L]]
    if (nExtra > 0L) {
        abort(glue(
            "{path} stores {nExtra} genotype call(s) beyond its ploidy, ",
            "which are not read."
        ))
    }
    invisible(NULL)
}

# ---- PLINK1 ----------------------------------------------------------------

# A PLINK1 fileset's panel: its .bim and .fam are read here, its genotypes a
# block at a time through PlinkMatrix (.plinkBedSeedRead()). The stem is
# recorded absolute, as a VCF's or GDS file's path is, so the panel still
# reads from another working directory.
# @noRd
.plinkBedPanel <- function(path) {
    .requireReaderPackage("PlinkMatrix", "PLINK1")
    named <- .plink1RequireFiles(path)
    stem <- file.path(normalizePath(dirname(named)), basename(named))
    meta <- .plink1ReadMeta(stem)
    nVariants <- nrow(meta$snpInfo)
    seed <- methods::new(
        "GenotypeDosageSeed",
        backend = stem,
        field = "plinkBed",
        rowRecord = seq_len(nVariants),
        rowAllele = rep(1L, nVariants),
        records = IRanges(seq_len(nVariants), width = 1L),
        dn = list(normalizeVariantId(meta$snpInfo$SNP), meta$sampleIds)
    )
    .dosageSeedPanel(seed, meta$snpInfo, "plink1", stem)
}

# Resolve the plink1 stem and assert its .bed/.bim/.fam all exist.
# @noRd
#' @importFrom checkmate assertFileExists
.plink1RequireFiles <- function(path) {
    stem <- .genotypeStem(path)
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
#' @importFrom readr read_table cols col_character
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
        CHR = withChrPrefix(bim$CHR),
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

# ---- Per-chromosome VCF / GDS / PLINK1 sets --------------------------------

# The panel of a per-chromosome VCF, GDS or PLINK1 set: the shards' panels
# stacked -- DelayedArray binds their assays, so nothing is read -- recording
# each chromosome's file as a sharded GenotypeHandle does. Every shard must
# hold the same samples, in the same order.
# @noRd
.stackDosagePanels <- function(panels, format, metaPath, paths) {
    .stackCheckSamples(panels, paths)
    rowRangesList <- unname(map(panels, SummarizedExperiment::rowRanges))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(
            dosage = exec(
                BiocGenerics::rbind,
                !!!unname(map(panels, .seDosage))
            )
        ),
        rowRanges = suppressWarnings(exec(c, !!!rowRangesList)),
        colData = SummarizedExperiment::colData(panels[[1L]]),
        metadata = list(
            genotypeSource = list(
                format = format,
                path = metaPath,
                chromPaths = .panelChromPaths(panels)
            )
        )
    )
}

# @noRd
.stackCheckSamples <- function(panels, paths) {
    first <- colnames(panels[[1L]])
    same <- map_lgl(map(panels, colnames), identical, first)
    if (all(same)) {
        return(invisible(NULL))
    }
    abort(glue(
        "genoMeta: all per-chromosome files must have ",
        "identical sample IDs in the same order (mismatch at ",
        "'{paths[[which(!same)[[1L]]]]}')."
    ))
}

# Each chromosome's file, named by chromosome; a chromosome found in two
# files is refused, as for a sharded GenotypeHandle.
# @noRd
.panelChromPaths <- function(panels) {
    chromPaths <- list_c(map(panels, .panelShardPaths))
    counts <- table(names(chromPaths))
    walk2(names(counts), as.integer(counts), .checkChromNotSplit)
    chromPaths
}

# One shard panel's chromosomes, each mapped to its file.
# @noRd
.panelShardPaths <- function(panel) {
    chroms <- unique(.ldRefChrom(panel))
    set_names(rep(.ldRefSource(panel)$path, length(chroms)), chroms)
}
