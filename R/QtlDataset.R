# =============================================================================
# QtlDataset S4 class
# -----------------------------------------------------------------------------
# Single-study individual-level QTL container: a GenotypeHandle, a named
# list of per-context phenotype SummarizedExperiments, optional genotype
# covariates, and constructor-level QC knobs (mafCutoff, macCutoff,
# xvarCutoff, imissCutoff, keepSamples, keepVariants). Backed by lazy
# getGenotypes / getResidualizedGenotypes accessors that apply QC at
# extraction time. The entry point for individual-level fine-mapping
# (fineMappingPipeline), TWAS weight learning (twasWeightsPipeline), and
# multi-study composition (MultiStudyQtlDataset).
# =============================================================================

#' @include AllGenerics.R
NULL

#' @title QTL Dataset (individual-level data for one study)
#' @description S4 container for a single QTL study's regional data. Holds a
#'   genotype handle plus per-context \code{SummarizedExperiment} objects
#'   carrying molecular-trait measurements. Each context's SE has
#'   \code{rowRanges} describing per-trait genomic positions and \code{colData}
#'   carrying per-context phenotype covariates. A single matrix of
#'   genotype-derived covariates (e.g., ancestry PCs) applies across contexts.
#'
#' @slot study Character (length 1). Study identifier; used in collection
#'   classes to tag downstream \code{FineMappingResult} / \code{TwasWeights}
#'   entries.
#' @slot genotypes A \code{GenotypeHandle} for lazy access to genotype dosages.
#' @slot phenotypes Named list of \code{SummarizedExperiment} objects, one per
#'   QTL context. Each SE has rows = molecular traits with positions in
#'   \code{rowRanges(se)}, columns = samples, and per-context covariates in
#'   \code{colData(se)}. Different contexts may carry different subsets of
#'   traits (rows); traits shared across contexts must have identical
#'   \code{rowRanges} entries (enforced by validity).
#' @slot genotypeCovariates Numeric matrix (samples x covariates) of
#'   genotype-derived covariates applied uniformly across all contexts (e.g.,
#'   ancestry PCs).
#' @slot scaleResiduals Logical (length 1). Whether residualization accessors
#'   scale residuals to unit variance.
#' @slot mafCutoff Numeric (length 1). Minor allele frequency threshold;
#'   variants with \code{MAF < mafCutoff} are dropped at extraction time inside
#'   \code{getGenotypes()} / \code{getResidualizedGenotypes()}. Default 0 (no
#'   filter).
#' @slot macCutoff Numeric (length 1). Minor allele count threshold; converted
#'   to a MAF threshold using \code{max(mafCutoff, macCutoff / (2 * n))} where
#'   \code{n} is the post-narrowing sample count of the extracted block. Default
#'   0 (no filter).
#' @slot xvarCutoff Numeric (length 1). Per-variant genotype variance threshold;
#'   variants with column variance below this are dropped at extraction time.
#'   Default 0 (no filter).
#' @slot imissCutoff Numeric (length 1). Per-sample genotype-missingness
#'   threshold; samples with a missing-genotype rate above this are dropped at
#'   extraction time. Default 0 (no filter).
#' @slot keepSamples Character vector of sample identifiers to retain prior to
#'   per-block QC; intersected with the genotype handle's \code{sampleIds} and
#'   the \code{samples} argument of \code{getGenotypes()}. Length 0 means no
#'   restriction.
#' @slot keepVariants Character vector of variant identifiers to retain prior to
#'   per-block QC. Length 0 means no restriction.
#' @slot keepIndel Logical (length 1). When \code{FALSE}, indel variants
#'   (alleles that are not single nucleotides) are dropped at extraction time.
#'   Default \code{TRUE}.
#' @export
setClass(
    "QtlDataset",
    representation(
        study = "character",
        genotypes = "GenotypeHandle",
        phenotypes = "list",
        genotypeCovariates = "matrix",
        scaleResiduals = "logical",
        mafCutoff = "numeric",
        macCutoff = "numeric",
        xvarCutoff = "numeric",
        imissCutoff = "numeric",
        keepSamples = "character",
        keepVariants = "character",
        keepIndel = "logical"
    ),
    prototype = prototype(keepIndel = TRUE),
    validity = function(object) .validateQtlDataset(object)
)

# Validity for QtlDataset: scalar/cutoff slots, the phenotypes list, and
# cross-context trait-position consistency. Returns TRUE or an error vector.
# @noRd
.validateQtlDataset <- function(object) {
    errors <- c(
        .qtlValidateScalars(object),
        .qtlValidatePhenotypes(object),
        .qtlValidateTraitPositions(object)
    )
    if (length(errors) == 0) TRUE else errors
}

# study / scaleResiduals / keepIndel scalars + the four non-negative cutoffs.
# @noRd
.qtlValidateScalars <- function(object) {
    errors <- character()
    if (length(object@study) != 1L || str_length(object@study) == 0L) {
        errors <- c(
            errors,
            "'study' must be a single non-empty character string"
        )
    }
    if (length(object@scaleResiduals) != 1L) {
        errors <- c(errors, "'scaleResiduals' must be a single logical value")
    }
    if (length(object@keepIndel) != 1L || is.na(object@keepIndel)) {
        errors <- c(errors, "'keepIndel' must be a single logical value")
    }
    for (nm in c("mafCutoff", "macCutoff", "xvarCutoff", "imissCutoff")) {
        v <- methods::slot(object, nm)
        if (length(v) != 1L || is.na(v) || !is.finite(v) || v < 0) {
            errors <- c(
                errors,
                glue("'{nm}' must be a single finite non-negative numeric")
            )
        }
    }
    errors
}

# phenotypes must be a non-empty, uniquely + non-empty-named list of
# SummarizedExperiments.
# @noRd
.qtlValidatePhenotypes <- function(object) {
    errors <- character()
    if (length(object@phenotypes) == 0L) {
        errors <- c(errors, "'phenotypes' must not be empty")
    }
    contextNames <- names(object@phenotypes)
    if (
        is.null(contextNames) ||
            any(str_length(contextNames) == 0L, na.rm = TRUE) ||
            any(is.na(contextNames))
    ) {
        errors <- c(
            errors,
            "'phenotypes' must be a named list with non-empty names"
        )
    } else if (n_distinct(contextNames) < length(contextNames)) {
        errors <- c(errors, "context names in 'phenotypes' must be unique")
    }
    for (ctx in seq_along(object@phenotypes)) {
        se <- object@phenotypes[[ctx]]
        if (!methods::is(se, "SummarizedExperiment")) {
            errors <- c(
                errors,
                glue(
                    "phenotypes[[{ctx}]] must be a SummarizedExperiment ",
                    "(got {class(se)[[1L]]})"
                )
            )
        }
    }
    errors
}

# TRUE when two GRanges share canonical chrom + start + end.
# @noRd
.qtlSameRange <- function(prev, this) {
    isTRUE(all.equal(
        canonChrom(GenomicRanges::seqnames(prev)),
        canonChrom(GenomicRanges::seqnames(this))
    )) &&
        GenomicRanges::start(prev) == GenomicRanges::start(this) &&
        GenomicRanges::end(prev) == GenomicRanges::end(this)
}

# A trait shared across contexts must have consistent rowRanges everywhere.
# @noRd
.qtlValidateTraitPositions <- function(object) {
    allSe <- length(object@phenotypes) > 1L &&
        all(map_lgl(object@phenotypes, methods::is, "SummarizedExperiment"))
    if (!allSe) {
        return(character())
    }
    errors <- character()
    traitToRange <- list()
    for (ctx in seq_along(object@phenotypes)) {
        se <- object@phenotypes[[ctx]]
        rr <- SummarizedExperiment::rowRanges(se)
        ids <- rownames(se)
        if (length(rr) != length(ids)) {
            next
        }
        for (i in seq_along(ids)) {
            tid <- ids[[i]]
            prev <- traitToRange[[tid]]
            if (is.null(prev)) {
                traitToRange[[tid]] <- rr[i]
            } else if (!.qtlSameRange(prev, rr[i])) {
                errors <- c(
                    errors,
                    glue(
                        "trait '{tid}' has inconsistent rowRanges across ",
                        "contexts"
                    )
                )
            }
        }
    }
    errors
}

# =============================================================================
# QtlDataset constructor and accessors
# =============================================================================

#' @title Create a QtlDataset Object
#' @description Construct a \code{QtlDataset} S4 object containing one study's
#'   individual-level QTL data: a genotype handle and a named list of
#'   \code{SummarizedExperiment} objects (one per QTL context), plus
#'   genotype-derived covariates and a residual-scaling flag.
#' @param study Character (length 1). Study identifier.
#' @param genotypes A \code{GenotypeHandle}.
#' @param phenotypes Named list of \code{SummarizedExperiment} objects, keyed by
#'   context. Each SE must have \code{rowRanges} carrying trait positions and
#'   \code{colData} carrying per-context phenotype covariates.
#' @param genotypeCovariates Numeric matrix of genotype-derived covariates
#'   (e.g., ancestry PCs); rows are samples.
#' @param scaleResiduals Logical (length 1). Default \code{TRUE}.
#' @param mafCutoff Numeric (length 1). Minor allele frequency threshold;
#'   variants with \code{MAF < mafCutoff} are dropped at extraction time inside
#'   \code{getGenotypes()} / \code{getResidualizedGenotypes()}. Default 0 (no
#'   filter).
#' @param macCutoff Numeric (length 1). Minor allele count threshold; converted
#'   to a MAF threshold using \code{max(mafCutoff, macCutoff / (2 * n))} where
#'   \code{n} is the post-narrowing sample count of the extracted block. Default
#'   0 (no filter).
#' @param xvarCutoff Numeric (length 1). Per-variant genotype variance
#'   threshold; variants with column variance below this are dropped at
#'   extraction time. Default 0 (no filter).
#' @param imissCutoff Numeric (length 1). Per-sample genotype-missingness
#'   threshold; samples with a missing-genotype rate above this are dropped at
#'   extraction time. Default 0 (no filter).
#' @param keepSamples Character vector of sample identifiers to retain prior to
#'   per-block QC; intersected with the genotype handle's \code{sampleIds} and
#'   the \code{samples} argument of \code{getGenotypes()}. Length 0 means no
#'   restriction.
#' @param keepVariants Character vector of variant identifiers to retain prior
#'   to per-block QC. Length 0 means no restriction.
#' @param keepIndel Logical (length 1). When \code{FALSE}, variants whose
#'   alleles are not single nucleotides (indels) are dropped at extraction.
#'   Default \code{TRUE} (keep all variants).
#' @return A \code{QtlDataset} object.
#' @examples
#' gh <- new("GenotypeHandle", path = "toy.gds", format = "gds",
#'   snpInfo = data.frame(SNP = paste0("rs", 1:3), CHR = "1",
#'     BP = c(100L, 200L, 300L), A1 = "A", A2 = "G"),
#'   nSamples = 6L, sampleIds = paste0("s", 1:6), pgenPtr = NULL)
#' rng <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000L, width = 500L))
#' names(rng) <- "ENSG1"
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(expression = matrix(rnorm(6), 1, 6,
#'     dimnames = list("ENSG1", paste0("s", 1:6)))), rowRanges = rng)
#' QtlDataset(study = "s1", genotypes = gh, phenotypes = list(brain = se),
#'   genotypeCovariates = matrix(0, 6, 0))
#' @export
QtlDataset <- function(
    study,
    genotypes,
    phenotypes,
    genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0),
    scaleResiduals = TRUE,
    mafCutoff = 0,
    macCutoff = 0,
    xvarCutoff = 0,
    imissCutoff = 0,
    keepSamples = character(0),
    keepVariants = character(0),
    keepIndel = TRUE
) {
    obj <- new(
        "QtlDataset",
        study = as.character(study),
        genotypes = genotypes,
        phenotypes = phenotypes,
        genotypeCovariates = as.matrix(genotypeCovariates),
        scaleResiduals = isTRUE(scaleResiduals),
        mafCutoff = as.numeric(mafCutoff),
        macCutoff = as.numeric(macCutoff),
        xvarCutoff = as.numeric(xvarCutoff),
        imissCutoff = as.numeric(imissCutoff),
        keepSamples = as.character(keepSamples),
        keepVariants = as.character(keepVariants),
        keepIndel = isTRUE(keepIndel)
    )
    validObject(obj)
    obj
}

#' @rdname getStudy
#' @export
setMethod("getStudy", "QtlDataset", function(x) x@study)

#' @rdname getContexts
#' @export
setMethod("getContexts", "QtlDataset", function(x) names(x@phenotypes))

#' @rdname getGenotypeCovariates
#' @export
setMethod("getGenotypeCovariates", "QtlDataset", function(x) {
    x@genotypeCovariates
})

#' @rdname getScaleResiduals
#' @export
setMethod("getScaleResiduals", "QtlDataset", function(x) x@scaleResiduals)

# --- Internal: resolve the variant-selection region for the genotype handle.
# Returns a GRanges (one or more ranges). When `traitId` is supplied, expand
# each trait's rowRange by `cisWindow` bp and take the union span (per the
# multi-trait rule: `[min(start) - cisWindow, max(end) + cisWindow]`). When
# `region` is supplied it is taken literally and may contain multiple ranges
# (e.g. for joint multi-region extraction), optionally extended per-range by
# `cisWindow`. Exactly one of (traitId, region) may be supplied; if neither is,
# return NULL meaning "all variants in handle".
.qtlResolveVariantRegion <- function(
    x,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL
) {
    if (!is.null(traitId) && !is.null(region)) {
        abort("Specify either `traitId` or `region`, not both.")
    }
    if (is.null(traitId) && is.null(region)) {
        return(NULL)
    }
    if (!is.null(traitId)) {
        return(.qtlTraitRegion(x, traitId, cisWindow))
    }
    .qtlLiteralRegion(region, cisWindow)
}

# The union span (+/- cisWindow) of a trait's rowRanges across all contexts.
# Requires cisWindow and a single shared chromosome.
# @noRd
.qtlTraitRegion <- function(x, traitId, cisWindow) {
    if (is.null(cisWindow) || length(cisWindow) != 1L || cisWindow < 0) {
        msg <- glue(
            "`cisWindow` is required (and must be non-negative) when ",
            "`traitId` is specified."
        )
        abort(msg)
    }
    perTraitRanges <- list()
    for (ctxIdx in seq_along(x@phenotypes)) {
        se <- x@phenotypes[[ctxIdx]]
        hits <- match(traitId, rownames(se))
        hits <- hits[!is.na(hits)]
        if (length(hits) > 0) {
            rr <- SummarizedExperiment::rowRanges(se)
            perTraitRanges[[length(perTraitRanges) + 1L]] <- rr[hits]
        }
    }
    if (length(perTraitRanges) == 0L) {
        abort("None of the requested traitId values were found in any context.")
    }
    allRanges <- exec(c, !!!perTraitRanges)
    chrs <- unique(as.character(GenomicRanges::seqnames(allRanges)))
    if (length(chrs) != 1L) {
        msg <- glue(
            "Multi-trait variant extraction requires all selected traits ",
            "to share a chromosome ",
            "(got: {str_flatten(chrs, ', ')})."
        )
        abort(msg)
    }
    spanStart <- max(1L, min(GenomicRanges::start(allRanges)) - cisWindow)
    spanEnd <- max(GenomicRanges::end(allRanges)) + cisWindow
    GenomicRanges::GRanges(
        seqnames = chrs,
        ranges = IRanges::IRanges(start = spanStart, end = spanEnd)
    )
}

# A literal `region` GRanges, each range optionally extended by cisWindow.
# @noRd
.qtlLiteralRegion <- function(region, cisWindow) {
    if (!methods::is(region, "GRanges")) {
        abort("`region` must be a GRanges object.")
    }
    if (length(region) == 0L) {
        abort("`region` must contain at least one range.")
    }
    if (is.null(cisWindow)) {
        return(region)
    }
    if (length(cisWindow) != 1L || cisWindow < 0) {
        abort("`cisWindow` must be a single non-negative value.")
    }
    GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(region),
        ranges = IRanges::IRanges(
            start = pmax(1L, GenomicRanges::start(region) - cisWindow),
            end = GenomicRanges::end(region) + cisWindow
        )
    )
}

# Per-trait genomic position: each trait's OWN rowRanges span across contexts
# (union min-start to max-end), WITHOUT the cisWindow expansion. One range per
# traitId (a 0-width chrUn sentinel for a trait absent from every context), for
# threading trait-position provenance onto QtlFineMappingResult / TwasWeights.
# @noRd
.qtlTraitPos <- function(x, traitIds) {
    # Extract chrom/start/end as plain vectors per trait (union span across
    # contexts), then build ONE fresh GRanges at the end. Combining per-context
    # GRanges with do.call(c, .) can trip S4 seqinfo reconciliation in some
    # GenomeInfoDb builds, so we avoid it entirely.
    n <- length(traitIds)
    chrs <- rep("chrUn", n)
    starts <- rep(1L, n)
    ends <- rep(1L, n)
    for (i in seq_len(n)) {
        tid <- traitIds[[i]]
        st <- Inf
        en <- -Inf
        ch <- NA_character_
        for (se in x@phenotypes) {
            h <- match(tid, rownames(se))
            h <- h[!is.na(h)]
            if (length(h) == 0L) {
                next
            }
            rr <- SummarizedExperiment::rowRanges(se)[h]
            ch <- as.character(GenomicRanges::seqnames(rr))[1L]
            st <- min(st, GenomicRanges::start(rr))
            en <- max(en, GenomicRanges::end(rr))
        }
        if (!is.na(ch)) {
            chrs[i] <- ch
            starts[i] <- as.integer(st)
            ends[i] <- as.integer(en)
        }
    }
    gr <- GenomicRanges::GRanges(
        chrs,
        IRanges::IRanges(start = starts, end = pmax(ends, starts))
    )
    names(gr) <- traitIds
    gr
}

#' @rdname getTraitPosition
#' @export
setMethod("getTraitPosition", "QtlDataset", function(x, traitId = NULL, ...) {
    tids <- if (is.null(traitId)) {
        unique(unlist(map(x@phenotypes, rownames)))
    } else {
        as.character(traitId)
    }
    .qtlTraitPos(x, tids)
})

# Internal: map a GRanges region (one or more ranges) into 1-based snpIdx into
# handle@snpInfo. Indices are unioned across ranges in range order (first
# occurrence wins), so overlapping ranges contribute each variant once.
.qtlVariantIndices <- function(x, region = NULL) {
    handle <- x@genotypes
    if (is.null(region)) {
        return(seq_len(nrow(handle@snpInfo)))
    }
    snpInfo <- handle@snpInfo
    siChr <- canonChrom(snpInfo$CHR)
    bp <- as.integer(snpInfo$BP)
    rChr <- canonChrom(GenomicRanges::seqnames(region))
    rStart <- GenomicRanges::start(region)
    rEnd <- GenomicRanges::end(region)
    idx <- integer(0)
    for (i in seq_along(region)) {
        idx <- c(idx, which(siChr == rChr[i] & bp >= rStart[i] & bp <= rEnd[i]))
    }
    unique(idx)
}

# Internal: keepIndel slot read, tolerant of QtlDataset objects serialized
# before the slot existed (treat a missing slot as TRUE = keep indels).
.qtlKeepIndel <- function(x) {
    isTRUE(tryCatch(x@keepIndel, error = function(e) TRUE))
}

# Internal: return a copy of a QtlDataset with the supplied filter cutoffs /
# keep-lists REPLACING the stored slot values (NULL = leave the stored value
# untouched). This lets a pipeline accept per-call filter overrides as ordinary
# arguments instead of forcing callers to mutate @slots directly (which bypasses
# the class's validity checks). Applied against a validated copy.
.qtlApplyFilterOverrides <- function(
    data,
    mafCutoff = NULL,
    macCutoff = NULL,
    xvarCutoff = NULL,
    imissCutoff = NULL,
    keepIndel = NULL,
    keepSamples = NULL,
    keepVariants = NULL
) {
    if (!is.null(mafCutoff)) {
        data@mafCutoff <- as.numeric(mafCutoff)
    }
    if (!is.null(macCutoff)) {
        data@macCutoff <- as.numeric(macCutoff)
    }
    if (!is.null(xvarCutoff)) {
        data@xvarCutoff <- as.numeric(xvarCutoff)
    }
    if (!is.null(imissCutoff)) {
        data@imissCutoff <- as.numeric(imissCutoff)
    }
    if (!is.null(keepIndel)) {
        data@keepIndel <- as.logical(keepIndel)
    }
    if (!is.null(keepSamples)) {
        data@keepSamples <- as.character(keepSamples)
    }
    if (!is.null(keepVariants)) {
        data@keepVariants <- as.character(keepVariants)
    }
    methods::validObject(data)
    data
}

# Internal: extract the panel dosage block (samples x variants) for the
# requested region, narrow to the requested sample set, and apply lazy QC
# (per-sample imiss filter, then per-variant max(mafCutoff,
# macCutoff / (2 * n)) and xvarCutoff filters). Used by getGenotypes,
# getResidualizedGenotypes (via getGenotypes), and getMaf so all three
# share a single variant/sample selection result.
#
# Returns a list:
#   geno       : numeric matrix (kept samples x kept variants)
#   variantIds : character vector of kept variant IDs (= colnames(geno))
#   sampleIds  : character vector of kept sample IDs (= rownames(geno))
#   maf        : numeric vector of per-variant MAF for kept variants
#   af         : numeric vector of per-variant effect-allele (A1) frequency
#                for kept variants. Directional (NOT folded to the minor
#                allele): the frequency of the dosage-counted allele, which
#                is A1 by the same convention the marginal betas use. `maf`
#                is `pmin(af, 1 - af)`.
.qtlExtractBlock <- function(
    x,
    traitId = NULL,
    region = NULL,
    cisWindow = NULL,
    samples = NULL
) {
    gr <- .qtlResolveVariantRegion(
        x,
        traitId = traitId,
        region = region,
        cisWindow = cisWindow
    )
    snpIdx <- .qtlVariantIndices(x, gr)
    if (length(snpIdx) == 0L) {
        return(.qtlEmptyBlockAllSamples(x))
    }
    # Apply keepVariants + indel restrictions before materialization so we do
    # not extract dosage we would immediately drop.
    snpIdx <- .qtlNarrowSnpIdx(x, snpIdx)
    if (length(snpIdx) == 0L) {
        return(.qtlEmptyBlock())
    }
    dosage <- .dosageMatrix(x@genotypes, snpIdx, meanImpute = FALSE)
    keep <- .qtlResolveSamples(dosage, x, samples)
    if (length(keep) == 0L) {
        return(.qtlEmptyBlockNoSamples(dosage))
    }
    dosage <- dosage[keep, , drop = FALSE]
    # Per-sample missingness filter.
    if (x@imissCutoff > 0 && nrow(dosage) > 0L && ncol(dosage) > 0L) {
        dosage <- dosage[
            rowMeans(is.na(dosage)) <= x@imissCutoff,
            ,
            drop = FALSE
        ]
    }
    filtered <- .qtlVariantFilters(dosage, x)
    dosage <- .qtlMeanImpute(filtered$dosage)
    list(
        geno = dosage,
        variantIds = colnames(dosage),
        sampleIds = rownames(dosage),
        maf = filtered$maf,
        af = filtered$af
    )
}

# Empty block preserving the full panel sample set (no variants selected).
# @noRd
.qtlEmptyBlockAllSamples <- function(x) {
    list(
        geno = matrix(
            numeric(0),
            nrow = x@genotypes@nSamples,
            ncol = 0L,
            dimnames = list(x@genotypes@sampleIds, character(0))
        ),
        variantIds = character(0),
        sampleIds = x@genotypes@sampleIds,
        maf = numeric(0),
        af = numeric(0)
    )
}

# Fully empty block (no variants, no samples).
# @noRd
.qtlEmptyBlock <- function() {
    list(
        geno = matrix(
            numeric(0),
            nrow = 0L,
            ncol = 0L,
            dimnames = list(character(0), character(0))
        ),
        variantIds = character(0),
        sampleIds = character(0),
        maf = numeric(0),
        af = numeric(0)
    )
}

# Empty block preserving the selected variants (no samples survived).
# @noRd
.qtlEmptyBlockNoSamples <- function(dosage) {
    list(
        geno = dosage[integer(0), , drop = FALSE],
        variantIds = colnames(dosage),
        sampleIds = character(0),
        maf = rep(NA_real_, ncol(dosage)),
        af = rep(NA_real_, ncol(dosage))
    )
}

# Narrow the selected variant indices by keepVariants (matched by chrom/pos/
# allele) and, unless kept, by dropping indels. Done before materialization.
# @noRd
.qtlNarrowSnpIdx <- function(x, snpIdx) {
    if (length(x@keepVariants) > 0L) {
        snpAll <- as.character(x@genotypes@snpInfo$SNP[snpIdx])
        km <- matchVariants(snpAll, as.character(x@keepVariants))
        keepMask <- logical(length(snpAll))
        keepMask[km$idxA] <- TRUE
        snpIdx <- snpIdx[keepMask]
    }
    if (length(snpIdx) > 0L && !.qtlKeepIndel(x)) {
        si <- x@genotypes@snpInfo
        # which() (not the mask) so an NA mask drops the variant rather than
        # injecting an NA index.
        snpMask <- str_length(as.character(si$A1[snpIdx])) == 1L &
            str_length(as.character(si$A2[snpIdx])) == 1L
        snpIdx <- snpIdx[which(snpMask)]
    }
    snpIdx
}

# Resolve the sample set: panel samples intersected with keepSamples and the
# per-call `samples` arg.
# @noRd
.qtlResolveSamples <- function(dosage, x, samples) {
    keep <- rownames(dosage)
    if (length(x@keepSamples) > 0L) {
        keep <- intersect(keep, x@keepSamples)
    }
    if (!is.null(samples)) {
        keep <- intersect(keep, as.character(samples))
    }
    keep
}

# Per-variant MAF / MAC / X-variance filters against the post-narrowing sample
# count. MAF is computed from the un-imputed dosage (A1 = the effect allele
# gives directional `af`; `maf` folds to the minor allele). Returns
# list(dosage, maf, af).
# @noRd
.qtlVariantFilters <- function(dosage, x) {
    nSamp <- nrow(dosage)
    if (ncol(dosage) == 0L) {
        # nocov start (unreachable: zero-variant paths return early above)
        return(list(dosage = dosage, maf = numeric(0), af = numeric(0)))
        # nocov end
    }
    nObs <- colSums(!is.na(dosage))
    sumD <- colSums(dosage, na.rm = TRUE)
    p <- if_else(nObs > 0L, sumD / (2 * nObs), NA_real_)
    afVec <- p
    mafVec <- pmin(p, 1 - p)
    effectiveMaf <- max(
        x@mafCutoff,
        if (nSamp > 0L) x@macCutoff / (2 * nSamp) else 0
    )
    keepVarMask <- !is.na(mafVec) & mafVec >= effectiveMaf
    if (x@xvarCutoff > 0 && nSamp > 1L) {
        mu <- if_else(nObs > 0L, sumD / nObs, 0)
        centered <- sweep(dosage, 2L, mu, FUN = "-")
        centered[is.na(centered)] <- 0
        varVec <- colSums(centered * centered) / (nSamp - 1L)
        keepVarMask <- keepVarMask & varVec >= x@xvarCutoff
    }
    list(
        dosage = dosage[, keepVarMask, drop = FALSE],
        maf = mafVec[keepVarMask],
        af = afVec[keepVarMask]
    )
}

# Mean-impute remaining missing dosage cells (per column) so downstream linear
# algebra is well-defined; MAF was already computed pre-imputation.
# @noRd
.qtlMeanImpute <- function(dosage) {
    if (!anyNA(dosage)) {
        return(dosage)
    }
    for (j in seq_len(ncol(dosage))) {
        col <- dosage[, j]
        na <- is.na(col)
        if (any(na)) {
            col[na] <- mean(col[!na])
            dosage[, j] <- col
        }
    }
    dosage
}

#' @rdname getGenotypes
#' @export
setMethod(
    "getGenotypes",
    "QtlDataset",
    function(
        x,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        samples = NULL,
        ...
    ) {
        .qtlExtractBlock(
            x,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            samples = samples
        )$geno
    }
)

#' @rdname getMaf
#' @export
setMethod(
    "getMaf",
    "QtlDataset",
    function(x, region = NULL, cisWindow = NULL, samples = NULL, ...) {
        block <- .qtlExtractBlock(
            x,
            traitId = NULL,
            region = region,
            cisWindow = cisWindow,
            samples = samples
        )
        out <- block$maf
        names(out) <- block$variantIds
        out
    }
)

#' @rdname getAf
#' @export
setMethod(
    "getAf",
    "QtlDataset",
    function(
        x,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        samples = NULL,
        ...
    ) {
        block <- .qtlExtractBlock(
            x,
            traitId = traitId,
            region = region,
            cisWindow = cisWindow,
            samples = samples
        )
        out <- block$af
        names(out) <- block$variantIds
        out
    }
)

#' @rdname getPhenotypes
#' @export
setMethod(
    "getPhenotypes",
    "QtlDataset",
    function(
        x,
        contexts,
        traitId = NULL,
        region = NULL,
        naAction = c("keep", "drop", "impute"),
        outlierAction = c("keep", "drop"),
        outlierPvalThreshold = 1e-3,
        ...
    ) {
        naAction <- arg_match(naAction)
        outlierAction <- arg_match(outlierAction)
        .qtlValidateContexts(x, contexts)
        out <- x@phenotypes[contexts]
        out <- .qtlFilterPhenotypes(
            out,
            contexts,
            traitId,
            region,
            naAction,
            outlierAction,
            outlierPvalThreshold
        )
        if (length(contexts) == 1L) out[[1L]] else out
    }
)

# `contexts` is required and must all be known phenotype contexts.
# @noRd
.qtlValidateContexts <- function(x, contexts) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
        msg <- glue(
            "`contexts` is required for getPhenotypes(QtlDataset). Pass a ",
            "character vector of one or more context names; use ",
            "getContexts(x) to list the available contexts."
        )
        abort(msg)
    }
    available <- names(x@phenotypes)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
        msg <- glue(
            "Unknown context(s): {str_flatten(bad, ', ')}. ",
            "Available: {str_flatten(available, ', ')}"
        )
        abort(msg)
    }
}

# Restrict each context's SE to `traitId`, warning about traits absent in a
# context.
# @noRd
.qtlFilterTraits <- function(out, contexts, traitId) {
    filtered <- map(
        seq_along(out),
        .qtlFilterTraitSe,
        out = out,
        traitId = traitId
    )
    set_names(filtered, contexts)
}

# Apply the optional trait / region / NA / outlier filters to the per-context
# phenotype SEs.
# @noRd
.qtlFilterPhenotypes <- function(
    out,
    contexts,
    traitId,
    region,
    naAction,
    outlierAction,
    outlierPvalThreshold
) {
    if (!is.null(traitId)) {
        out <- .qtlFilterTraits(out, contexts, traitId)
    }
    if (!is.null(region)) {
        out <- set_names(
            map(out, .qtlSeInRegion, region = region),
            contexts
        )
    }
    if (naAction != "keep") {
        out <- set_names(
            map(out, .qtlApplyPhenoNaAction, naAction = naAction),
            contexts
        )
    }
    if (outlierAction != "keep") {
        out <- set_names(
            map(
                out,
                .qtlApplyPhenoOutliers,
                action = outlierAction,
                pvalThreshold = outlierPvalThreshold
            ),
            contexts
        )
    }
    out
}

# Internal: apply naAction to a SummarizedExperiment slice. SE assay rows
# are traits and columns are samples.
#   "drop"   -> drop samples (cols) where any selected trait is NA
#   "impute" -> mean-impute each trait (row) independently over its
#               non-NA sample values
# Operates jointly over the rows currently in `se` -- the caller is
# expected to have already subset the SE to the user's requested
# (traitId, region) subset.
.qtlApplyPhenoNaAction <- function(se, naAction) {
    assayName <- SummarizedExperiment::assayNames(se)[[1L]]
    Y <- SummarizedExperiment::assay(se, assayName)
    if (length(Y) == 0L) {
        return(se)
    }
    if (naAction == "drop") {
        keepSamp <- colSums(is.na(Y)) == 0L
        se <- se[, keepSamp, drop = FALSE]
    } else if (naAction == "impute") {
        if (anyNA(Y)) {
            for (j in seq_len(nrow(Y))) {
                row <- Y[j, ]
                na <- is.na(row)
                if (any(na)) {
                    obs <- row[!na]
                    row[na] <- if (length(obs) > 0L) mean(obs) else 0
                    Y[j, ] <- row
                }
            }
            SummarizedExperiment::assay(se, assayName) <- Y
        }
    }
    se
}

# Multivariate-outlier keep mask via Mahalanobis distance against a
# (preferably robust) centre / covariance estimate. Returns a logical
# vector of length nrow(Y); TRUE = keep, FALSE = drop.
#
# When the `robustbase` package is installed, the centre and covariance
# come from `robustbase::covMcd` (minimum-covariance-determinant) so the
# detector itself is resistant to the outliers it's trying to find.
# Without robustbase we fall back to `colMeans` / `cov` with a one-shot
# message; the test then still works but its estimates are pulled by
# the very outliers it should be flagging.
#
# Significance: per-sample chi-squared(p) p-value with Bonferroni
# correction over the sample count. A sample is flagged when its
# corrected p-value falls below `pvalThreshold`. With single-trait Y
# (ncol == 1) this reduces to the standard z-test on (y - center)/sd.
#
# Returns all-TRUE (no-op) when there are too few samples to support
# a covariance estimate (n < p + 2).
.qtlOutlierKeepMask <- function(Y, pvalThreshold) {
    Y <- as.matrix(Y)
    n <- nrow(Y)
    p <- ncol(Y)
    if (n == 0L || p == 0L) {
        return(rep(TRUE, n))
    }
    if (n < p + 2L) {
        msg <- glue(
            "outlier detection skipped: {n} samples < {p} traits + 2 ",
            "needed for a covariance estimate."
        )
        warn(msg)
        return(rep(TRUE, n))
    }
    if (requireNamespace("robustbase", quietly = TRUE)) {
        mcd <- tryCatch(robustbase::covMcd(Y), error = function(e) NULL)
        if (!is.null(mcd)) {
            ctr <- mcd$center
            covMat <- mcd$cov
        } else {
            ctr <- colMeans(Y)
            covMat <- stats::cov(Y)
        }
    } else {
        msg <- glue(
            "outlier detection: install 'robustbase' for an MCD-based ",
            "estimator; falling back to non-robust colMeans/cov."
        )
        inform(msg)
        ctr <- colMeans(Y)
        covMat <- stats::cov(Y)
    }
    invCov <- tryCatch(solve(covMat), error = function(e) MASS::ginv(covMat))
    Yc <- sweep(Y, 2L, ctr)
    d2 <- rowSums((Yc %*% invCov) * Yc)
    raw <- stats::pchisq(d2, df = p, lower.tail = FALSE)
    raw >= (pvalThreshold / n)
}

# Wrapper: apply the keep-mask to a SummarizedExperiment slice. SE
# columns are samples; transpose the assay (traits x samples) before
# calling .qtlOutlierKeepMask which expects samples x traits.
.qtlApplyPhenoOutliers <- function(se, action, pvalThreshold) {
    if (action == "keep") {
        return(se)
    }
    assayName <- SummarizedExperiment::assayNames(se)[[1L]]
    Y <- t(SummarizedExperiment::assay(se, assayName))
    keep <- .qtlOutlierKeepMask(Y, pvalThreshold)
    if (all(keep)) {
        return(se)
    }
    se[, keep, drop = FALSE]
}

#' @rdname getPhenotypeCovariates
#' @export
setMethod("getPhenotypeCovariates", "QtlDataset", function(x, contexts) {
    if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
        abort("`contexts` is required.")
    }
    available <- names(x@phenotypes)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
        msg <- glue("Unknown context(s): {str_flatten(bad, ', ')}")
        abort(msg)
    }
    out <- map(contexts, .qtlContextColData, x = x)
    names(out) <- contexts
    out
})

# Internal: residualize a numeric matrix Y (n x k) against a covariate
# matrix C (n x p) via pivoted QR decomposition. Adds an intercept column
# to C. When C is rank-deficient (e.g., union of all contexts' phenotype
# covariates includes collinear / duplicate columns), the pivoted QR drops
# the redundant columns automatically. Optionally rescales each residual
# column to unit standard deviation; constant-valued columns are left
# unchanged.
.qtlResidualizeQr <- function(Y, C, scaleResiduals = TRUE) {
    X <- if (is.null(C) || ncol(C) == 0L) {
        matrix(
            1,
            nrow = nrow(Y),
            ncol = 1L,
            dimnames = list(rownames(Y), "intercept")
        )
    } else {
        cbind(intercept = 1, C)
    }
    # `qr.resid` does not support LAPACK pivoted QR, so use `lm.fit`. It
    # handles rank-deficient designs gracefully via base-R's pivoted QR
    # internally -- same effect the LAPACK path was meant to deliver.
    res <- stats::lm.fit(x = X, y = Y)$residuals
    res <- as.matrix(res)
    rownames(res) <- rownames(Y)
    colnames(res) <- colnames(Y)
    if (isTRUE(scaleResiduals)) {
        sds <- apply(res, 2L, stats::sd, na.rm = TRUE)
        # `sds == 0` exact-zero test is unreliable for residuals coming out of
        # lm.fit on a constant Y: roundoff gives sd ~ 1e-16 instead of 0, and
        # dividing the (also-tiny) residuals by it amplifies floating-point
        # noise to unit-scale. Treat anything below sqrt(.Machine$double.eps)
        # as effectively zero (column is constant) and skip rescaling.
        nearZero <- !is.finite(sds) | sds < sqrt(.Machine$double.eps)
        sds[nearZero] <- 1
        res[, nearZero] <- 0
        res <- sweep(res, 2L, sds, FUN = "/")
    }
    res
}

# Resolve the phenotype-covariate selection for one context: NULL requested ->
# all available covariates; otherwise validate the requested names are present.
# @noRd
.qtlResolveOne <- function(ctx, requested, x) {
    se <- x@phenotypes[[ctx]]
    avail <- colnames(SummarizedExperiment::colData(se))
    if (is.null(requested)) {
        return(avail)
    }
    keep <- intersect(requested, avail)
    if (length(keep) != length(requested)) {
        missingNames <- setdiff(requested, avail)
        msg <- glue(
            "phenotypeCovariatesToResidualize: context '{ctx}' has no ",
            "covariate(s) named: {str_flatten(missingNames, ', ')}"
        )
        abort(msg)
    }
    keep
}

# Internal: validate and resolve the `*ToResidualize` argument against a
# set of contexts and the covariates actually present in those contexts'
# colData. Accepts either NULL (use all), a character vector (apply to all
# listed contexts), or a named list keyed by context. Returns a named list
# keyed by context giving the actual character vector of covariate names
# to use for that context (or character(0) if none). Errors when:
#   - a named-list key is not in `contexts`
#   - `contexts` contains entries missing from a supplied named-list
#     (per the rule: named-list keys must equal `contexts`)
#   - an explicitly requested name matches no actual covariate
.qtlResolvePhenoSelection <- function(x, contexts, toResidualize) {
    if (is.null(toResidualize)) {
        return(set_names(
            map(contexts, .qtlResolveOne, requested = NULL, x = x),
            contexts
        ))
    }
    if (is.list(toResidualize)) {
        return(.qtlPhenoSelectionList(x, contexts, toResidualize))
    }
    if (is.character(toResidualize)) {
        return(set_names(
            map(contexts, .qtlResolveOne, requested = toResidualize, x = x),
            contexts
        ))
    }
    msg <- glue(
        "phenotypeCovariatesToResidualize must be NULL, a character vector, ",
        "or a named list keyed by context."
    )
    abort(msg)
}

# List-form phenotypeCovariatesToResidualize: must be named with EXACTLY the
# `contexts` set. Resolves each context's selection.
# @noRd
.qtlPhenoSelectionList <- function(x, contexts, toResidualize) {
    toResNames <- names(toResidualize)
    if (
        is.null(toResNames) || any(str_length(toResNames) == 0L, na.rm = TRUE)
    ) {
        msg <- glue(
            "phenotypeCovariatesToResidualize: when supplied as a list, it ",
            "must be named with context names."
        )
        abort(msg)
    }
    badKeys <- setdiff(names(toResidualize), contexts)
    if (length(badKeys) > 0L) {
        msg <- glue(
            "phenotypeCovariatesToResidualize: list key(s) not in ",
            "`contexts`: {str_flatten(badKeys, ', ')}"
        )
        abort(msg)
    }
    missingKeys <- setdiff(contexts, names(toResidualize))
    if (length(missingKeys) > 0L) {
        msg <- glue(
            "phenotypeCovariatesToResidualize: list does not cover all ",
            "`contexts`. Per-context lists must have exactly the same context ",
            "set as `contexts`. Missing keys: ",
            "{str_flatten(missingKeys, ', ')}"
        )
        abort(msg)
    }
    set_names(
        map(contexts, .qtlResolveContext, toResidualize = toResidualize, x = x),
        contexts
    )
}

# Internal: validate the genotype-covariate selection vector. Returns
# character(0) when nothing selected, the resolved set otherwise.
.qtlResolveGenoSelection <- function(x, toResidualize) {
    avail <- colnames(x@genotypeCovariates)
    if (is.null(avail)) {
        avail <- character(0)
    }
    if (is.null(toResidualize)) {
        return(avail)
    }
    keep <- intersect(toResidualize, avail)
    if (length(keep) != length(toResidualize)) {
        missingNames <- setdiff(toResidualize, avail)
        msg <- glue(
            "genotypeCovariatesToResidualize: no covariate(s) named: ",
            "{str_flatten(missingNames, ', ')}"
        )
        abort(msg)
    }
    keep
}

# Internal: build the covariate matrix used for residualization, given a
# set of contexts, the resolved per-context phenotype selections, and the
# resolved genotype-covariate selection. Honors the inclusion flags. For
# `length(contexts) == 1` (per-context mode), the per-context phenotype
# covariates and the genotype covariates are taken with no cross-context
# alignment. For `length(contexts) >= 2` (joint mode), per-context
# phenotype covariates from all listed contexts are concatenated
# (prefixed with "{context}." to keep same-named columns distinct) and
# the sample set is intersected across all contributing matrices.
# Returns a single matrix (with rownames = sample IDs) or NULL.
.qtlBuildResidualizationDesign <- function(
    x,
    contexts,
    phenoSelection,
    genoSelection,
    includePheno,
    includeGeno
) {
    perContext <- if (includePheno) {
        .qtlPhenoCovBlocks(x, contexts, phenoSelection)
    } else {
        list()
    }
    gCov <- if (includeGeno && length(genoSelection) > 0L) {
        x@genotypeCovariates[, genoSelection, drop = FALSE]
    } else {
        matrix(numeric(0), nrow = 0, ncol = 0)
    }
    haveAny <- length(perContext) > 0L || (!is.null(gCov) && ncol(gCov) > 0L)
    if (!haveAny) {
        return(NULL)
    }
    .qtlAlignCovariates(perContext, gCov)
}

# Per-context phenotype covariate matrices (colData columns from
# phenoSelection), column names prefixed with the context.
# @noRd
.qtlPhenoCovBlocks <- function(x, contexts, phenoSelection) {
    perContext <- list()
    for (ctx in contexts) {
        keep <- phenoSelection[[ctx]]
        if (length(keep) == 0L) {
            next
        }
        se <- x@phenotypes[[ctx]]
        cd <- as.matrix(as.data.frame(SummarizedExperiment::colData(se)))
        cdMat <- cd[, keep, drop = FALSE]
        colnames(cdMat) <- str_c(ctx, ".", colnames(cdMat))
        if (is.null(rownames(cdMat))) {
            rownames(cdMat) <- as.character(
                rownames(SummarizedExperiment::colData(se))
            )
        }
        perContext[[ctx]] <- cdMat
    }
    perContext
}

# Intersect the covariate blocks to their common samples and column-bind them
# into one design matrix; NULL when no samples are shared.
# @noRd
.qtlAlignCovariates <- function(perContext, gCov) {
    sampleSets <- list()
    for (mat in perContext) {
        if (!is.null(rownames(mat))) {
            sampleSets[[length(sampleSets) + 1L]] <- rownames(mat)
        }
    }
    if (!is.null(gCov) && ncol(gCov) > 0L && !is.null(rownames(gCov))) {
        sampleSets[[length(sampleSets) + 1L]] <- rownames(gCov)
    }
    common <- if (length(sampleSets) == 0L) {
        character(0)
    } else {
        reduce(sampleSets, intersect)
    }
    if (length(common) == 0L) {
        return(NULL)
    }
    blocks <- list()
    for (mat in perContext) {
        blocks[[length(blocks) + 1L]] <- mat[common, , drop = FALSE]
    }
    if (!is.null(gCov) && ncol(gCov) > 0L) {
        blocks[[length(blocks) + 1L]] <- gCov[common, , drop = FALSE]
    }
    exec(cbind, !!!blocks)
}

# Internal: resolve missing values in the residualization covariate matrix
# `C` (samples x covariates) before it reaches `stats::lm.fit`, which does
# not tolerate NA in the design and would otherwise error. Two strategies:
#   - "impute" (default): replace each covariate column's NA cells with that
#     column's observed (non-NA) mean. A wholly-missing column has an
#     undefined mean and is filled with 0, so it contributes nothing to the
#     fit rather than poisoning every row.
#   - "drop": complete-case -- remove any sample (row) carrying an NA in any
#     covariate. The caller's downstream sample intersection then narrows the
#     response matrix (G or Y) to the retained samples.
# Returns the cleaned matrix (fewer rows possible under "drop"), or `C`
# unchanged when it is NULL or already NA-free.
.qtlHandleCovariateNa <- function(C, action = c("impute", "drop")) {
    action <- arg_match(action)
    if (is.null(C) || !anyNA(C)) {
        return(C)
    }
    if (action == "drop") {
        keep <- rowSums(is.na(C)) == 0L
        return(C[keep, , drop = FALSE])
    }
    for (j in seq_len(ncol(C))) {
        col <- C[, j]
        na <- is.na(col)
        if (any(na)) {
            mu <- mean(col[!na])
            col[na] <- if (is.finite(mu)) mu else 0
            C[, j] <- col
        }
    }
    C
}

# Internal: resolve a (convenience, precise) flag pair to a single boolean.
# `missing*` arguments are passed as the result of `missing()` evaluated in
# the calling method to detect whether the user explicitly set the value.
# Rules:
#   - both missing: returns TRUE (the documented default)
#   - only convenience set: returns convenience
#   - only precise set: returns precise
#   - both set: must agree, else error
.qtlResolveResidualizationFlag <- function(
    conveniencePassed,
    convenienceMissing,
    precisePassed,
    preciseMissing,
    convenienceName,
    preciseName
) {
    if (preciseMissing && convenienceMissing) {
        return(TRUE)
    }
    if (preciseMissing) {
        return(isTRUE(conveniencePassed))
    }
    if (convenienceMissing) {
        return(isTRUE(precisePassed))
    }
    if (isTRUE(conveniencePassed) != isTRUE(precisePassed)) {
        msg <- glue(
            "Conflicting values: `{convenienceName}` = {conveniencePassed} ",
            "and `{preciseName}` = {precisePassed}. ",
            "Set only one, or pass consistent values."
        )
        abort(msg)
    }
    isTRUE(precisePassed)
}

#' @rdname getResidualizedGenotypes
#' @export
setMethod(
    "getResidualizedGenotypes",
    "QtlDataset",
    function(
        x,
        contexts,
        traitId = NULL,
        region = NULL,
        cisWindow = NULL,
        samples = NULL,
        phenotypeCovariatesToResidualize = NULL,
        genotypeCovariatesToResidualize = NULL,
        residualizePhenotypeCovariates = TRUE,
        residualizeGenotypeCovariates = TRUE,
        residualizePhenotypeCovariatesFromGenotypes = NULL,
        residualizeGenotypeCovariatesFromGenotypes = NULL,
        covariateNaAction = c("impute", "drop"),
        ...
    ) {
        if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
            msg <- glue(
                "`contexts` is required for ",
                "getResidualizedGenotypes(QtlDataset). ",
                "Use getContexts(x) to list the available contexts. ",
                "Pass a single context for per-context mode or multiple ",
                "contexts for joint mode (sample intersection)."
            )
            abort(msg)
        }
        covariateNaAction <- arg_match(covariateNaAction)
        convPhenoMissing <- missing(residualizePhenotypeCovariates)
        convGenoMissing <- missing(residualizeGenotypeCovariates)
        precPhenoMissing <- missing(
            residualizePhenotypeCovariatesFromGenotypes
        ) ||
            is.null(residualizePhenotypeCovariatesFromGenotypes)
        precGenoMissing <- missing(
            residualizeGenotypeCovariatesFromGenotypes
        ) ||
            is.null(residualizeGenotypeCovariatesFromGenotypes)
        p <- as.list(environment())
        p$dots <- list(...)
        .qtlResidualizedGenotypesImpl(p)
    }
)

# Intersect a genotype matrix G and covariate design C to their common samples
# (no-op when C is NULL); errors if they share none.
# @noRd
.qtlAlignGC <- function(G, C, contexts) {
    if (is.null(C)) {
        return(list(G = G, C = C))
    }
    common <- intersect(rownames(G), rownames(C))
    if (length(common) == 0L) {
        msg <- glue(
            "No samples in common between the genotype matrix and the ",
            "covariate matrix for contexts: ",
            "{str_flatten(contexts, ', ')}"
        )
        abort(msg)
    }
    list(G = G[common, , drop = FALSE], C = C[common, , drop = FALSE])
}

# getResidualizedGenotypes worker: resolve the covariate inclusion flags +
# selections, extract genotypes, build + NA-handle the covariate design, and
# QR-residualize. `p` holds the setMethod args + the precomputed `missing()`
# booleans.
# @noRd
.qtlResidualizedGenotypesImpl <- function(p) {
    bad <- setdiff(p$contexts, names(p$x@phenotypes))
    if (length(bad) > 0L) {
        msg <- glue("Unknown context(s): {str_flatten(bad, ', ')}")
        abort(msg)
    }
    includePheno <- .qtlResolveResidualizationFlag(
        p$residualizePhenotypeCovariates,
        p$convPhenoMissing,
        p$residualizePhenotypeCovariatesFromGenotypes,
        p$precPhenoMissing,
        "residualizePhenotypeCovariates",
        "residualizePhenotypeCovariatesFromGenotypes"
    )
    includeGeno <- .qtlResolveResidualizationFlag(
        p$residualizeGenotypeCovariates,
        p$convGenoMissing,
        p$residualizeGenotypeCovariatesFromGenotypes,
        p$precGenoMissing,
        "residualizeGenotypeCovariates",
        "residualizeGenotypeCovariatesFromGenotypes"
    )
    phenoSel <- .qtlResolvePhenoSelection(
        p$x,
        p$contexts,
        p$phenotypeCovariatesToResidualize
    )
    genoSel <- .qtlResolveGenoSelection(p$x, p$genotypeCovariatesToResidualize)
    G <- getGenotypes(
        p$x,
        traitId = p$traitId,
        region = p$region,
        cisWindow = p$cisWindow,
        samples = p$samples
    )
    if (ncol(G) == 0L) {
        return(G)
    }
    C <- .qtlBuildResidualizationDesign(
        p$x,
        contexts = p$contexts,
        phenoSelection = phenoSel,
        genoSelection = genoSel,
        includePheno = includePheno,
        includeGeno = includeGeno
    )
    C <- .qtlHandleCovariateNa(C, p$covariateNaAction)
    aligned <- .qtlAlignGC(G, C, p$contexts)
    .qtlResidualizeQr(aligned$G, aligned$C, scaleResiduals = p$x@scaleResiduals)
}

#' @rdname getResidualizedPhenotypes
#' @export
setMethod(
    "getResidualizedPhenotypes",
    "QtlDataset",
    function(
        x,
        contexts,
        traitId = NULL,
        region = NULL,
        phenotypeCovariatesToResidualize = NULL,
        genotypeCovariatesToResidualize = NULL,
        residualizePhenotypeCovariates = TRUE,
        residualizeGenotypeCovariates = TRUE,
        residualizePhenotypeCovariatesFromPhenotypes = NULL,
        residualizeGenotypeCovariatesFromPhenotypes = NULL,
        naAction = c("keep", "drop", "impute"),
        covariateNaAction = c("impute", "drop"),
        outlierAction = c("keep", "drop"),
        outlierPvalThreshold = 1e-3,
        ...
    ) {
        if (missing(contexts) || is.null(contexts) || length(contexts) == 0L) {
            abort("`contexts` is required for getResidualizedPhenotypes().")
        }
        naAction <- arg_match(naAction)
        covariateNaAction <- arg_match(covariateNaAction)
        outlierAction <- arg_match(outlierAction)
        convPhenoMissing <- missing(residualizePhenotypeCovariates)
        convGenoMissing <- missing(residualizeGenotypeCovariates)
        precPhenoMissing <- missing(
            residualizePhenotypeCovariatesFromPhenotypes
        ) ||
            is.null(residualizePhenotypeCovariatesFromPhenotypes)
        precGenoMissing <- missing(
            residualizeGenotypeCovariatesFromPhenotypes
        ) ||
            is.null(residualizeGenotypeCovariatesFromPhenotypes)
        p <- as.list(environment())
        p$dots <- list(...)
        .qtlResidualizedPhenotypesImpl(p)
    }
)

# NA-handled raw phenotypes per context (re-wrapped to a list so single- and
# multi-context callers see the same shape).
# @noRd
.qtlResidPhenoY <- function(x, contexts, traitId, region, naAction) {
    Yraw <- getPhenotypes(
        x,
        contexts = contexts,
        traitId = traitId,
        region = region,
        naAction = naAction
    )
    if (length(contexts) == 1L) {
        Yraw <- set_names(list(Yraw), contexts)
    }
    Yraw
}

# Residualize one context's phenotypes against the covariate design (intersected
# to common samples) and drop residual-scale outliers.
# @noRd
.qtlResidualizeContextPheno <- function(
    se,
    C,
    ctx,
    outlierAction,
    outlierPvalThreshold,
    scaleResiduals
) {
    Y <- t(SummarizedExperiment::assay(se)) # samples x traits
    Cctx <- NULL
    if (!is.null(C)) {
        common <- intersect(rownames(Y), rownames(C))
        if (length(common) == 0L) {
            msg <- glue(
                "context '{ctx}': no samples shared between phenotype data ",
                "and the resolved covariate matrix."
            )
            abort(msg)
        }
        Y <- Y[common, , drop = FALSE]
        Cctx <- C[common, , drop = FALSE]
    }
    Yres <- .qtlResidualizeQr(Y, Cctx, scaleResiduals = scaleResiduals)
    # Outlier detection on the residualized scale.
    if (outlierAction != "keep") {
        keep <- .qtlOutlierKeepMask(Yres, outlierPvalThreshold)
        if (!all(keep)) {
            Yres <- Yres[keep, , drop = FALSE]
        }
    }
    Yres
}

# Resolve the phenotype/genotype covariate inclusion flags (convenience vs
# precise `*FromPhenotypes`) for getResidualizedPhenotypes.
# @noRd
.qtlResidPhenoFlags <- function(p) {
    list(
        includePheno = .qtlResolveResidualizationFlag(
            p$residualizePhenotypeCovariates,
            p$convPhenoMissing,
            p$residualizePhenotypeCovariatesFromPhenotypes,
            p$precPhenoMissing,
            "residualizePhenotypeCovariates",
            "residualizePhenotypeCovariatesFromPhenotypes"
        ),
        includeGeno = .qtlResolveResidualizationFlag(
            p$residualizeGenotypeCovariates,
            p$convGenoMissing,
            p$residualizeGenotypeCovariatesFromPhenotypes,
            p$precGenoMissing,
            "residualizeGenotypeCovariates",
            "residualizeGenotypeCovariatesFromPhenotypes"
        )
    )
}

# getResidualizedPhenotypes worker: resolve covariate inclusion + selections,
# NA-handle Y, build the covariate design, and per-context residualize +
# outlier-filter. `p` holds the setMethod args + precomputed missing() flags.
# @noRd
.qtlResidualizedPhenotypesImpl <- function(p) {
    bad <- setdiff(p$contexts, names(p$x@phenotypes))
    if (length(bad) > 0L) {
        msg <- glue("Unknown context(s): {str_flatten(bad, ', ')}")
        abort(msg)
    }
    flags <- .qtlResidPhenoFlags(p)
    includePheno <- flags$includePheno
    includeGeno <- flags$includeGeno
    phenoSel <- .qtlResolvePhenoSelection(
        p$x,
        p$contexts,
        p$phenotypeCovariatesToResidualize
    )
    genoSel <- .qtlResolveGenoSelection(p$x, p$genotypeCovariatesToResidualize)
    Yraw <- .qtlResidPhenoY(
        p$x,
        p$contexts,
        p$traitId,
        p$region,
        p$naAction
    )
    C <- .qtlBuildResidualizationDesign(
        p$x,
        contexts = p$contexts,
        phenoSelection = phenoSel,
        genoSelection = genoSel,
        includePheno = includePheno,
        includeGeno = includeGeno
    )
    C <- .qtlHandleCovariateNa(C, p$covariateNaAction)
    out <- set_names(
        map(p$contexts, .qtlResidualizeContext, Yraw = Yraw, C = C, p = p),
        p$contexts
    )
    if (length(p$contexts) == 1L) out[[1L]] else out
}


#' @rdname show-methods
#' @export
setMethod("show", "QtlDataset", function(object) {
    nCtx <- length(object@phenotypes)
    ctxNames <- names(object@phenotypes)
    totalTraits <- length(unique(unlist(
        map(object@phenotypes, rownames),
        use.names = FALSE
    )))
    cat(glue("QtlDataset for study '{object@study}'\n", .trim = FALSE))
    cat(glue(
        "  {nCtx} context(s): {str_flatten(ctxNames, ', ')}\n",
        .trim = FALSE
    ))
    cat(glue("  {totalTraits} unique traits across contexts\n", .trim = FALSE))
    cat(glue(
        "  Genotypes: {object@genotypes@format} @ {object@genotypes@path}\n",
        .trim = FALSE
    ))
    cat(glue(
        "  Genotype covariates: {ncol(object@genotypeCovariates)} cols\n",
        .trim = FALSE
    ))
    cat(glue("  Scale residuals: {object@scaleResiduals}\n", .trim = FALSE))
})

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# Restrict context `i`'s SE to the requested traits (warns about absent ones).
# @noRd
.qtlFilterTraitSe <- function(i, out, traitId) {
    se <- out[[i]]
    ctx <- names(out)[[i]]
    present <- intersect(traitId, rownames(se))
    missing <- setdiff(traitId, rownames(se))
    if (length(missing) > 0L) {
        msg <- glue(
            "context '{ctx}' is missing trait(s): ",
            "{str_flatten(missing, ', ')}"
        )
        warn(msg)
    }
    se[present, , drop = FALSE]
}

# One context's SE restricted to features overlapping `region`.
# @noRd
.qtlSeInRegion <- function(se, region) {
    rr <- SummarizedExperiment::rowRanges(se)
    se[IRanges::overlapsAny(rr, region), , drop = FALSE]
}

# One context's covariate matrix (colData of its phenotype SE).
# @noRd
.qtlContextColData <- function(ctx, x) {
    se <- x@phenotypes[[ctx]]
    cd <- SummarizedExperiment::colData(se)
    as.matrix(as.data.frame(cd))
}

# Resolve the per-context residualization spec for context `ctx`.
# @noRd
.qtlResolveContext <- function(ctx, toResidualize, x) {
    .qtlResolveOne(ctx, toResidualize[[ctx]], x)
}

# Residualize + filter one context's raw phenotype matrix against covariates C.
# @noRd
.qtlResidualizeContext <- function(ctx, Yraw, C, p) {
    .qtlResidualizeContextPheno(
        Yraw[[ctx]],
        C,
        ctx,
        p$outlierAction,
        p$outlierPvalThreshold,
        p$x@scaleResiduals
    )
}
