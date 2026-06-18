#' @title S4 Class Definitions
#' @description All S4 class definitions for pecotmr: genotype handles,
#'   summary statistics, LD references, annotations, fine-mapping results,
#'   TWAS weights, regional data, and heritability estimates.
#' @name pecotmr-classes
#' @keywords internal
#' @importFrom methods setClass setMethod new is validObject
#' @importFrom S4Vectors mcols mcols<-
NULL

# =============================================================================
# LD Block Definitions
# =============================================================================

#' @title LD Block Definitions
#' @description Defines approximately independent LD blocks for local
#'   estimation. Typically derived from Berisa & Pickrell (2016) or
#'   user-provided boundaries.
#' @slot blocks A \code{GRanges} object with one range per LD block.
#'   Metadata columns may include \code{blockId} (integer).
#' @slot genome Character string identifying the genome build (e.g., "hg19", "hg38").
#' @export
setClass("LdBlocks",
  representation(
    blocks = "GRanges",
    genome = "character"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@genome) != 1L)
      errors <- c(errors, "'genome' must be a single character string")
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Genotype Handle
# =============================================================================

#' @title Genotype Handle
#' @description Lightweight handle to genotype data in any supported format.
#'   Stores the file path, detected format, and cached SNP metadata. Used to
#'   defer reading genotypes until block-level extraction is needed.
#' @slot path Character, path to the genotype file (or stem for plink).
#' @slot format Character, one of "gds", "vcf", "plink1", "plink2".
#' @slot snpInfo A \code{data.frame} with columns \code{SNP}, \code{CHR},
#'   \code{BP}, \code{A1}, \code{A2}. Cached on first access.
#' @slot nSamples Integer, number of samples.
#' @slot sampleIds Character vector of sample identifiers.
#' @slot pgenPtr An external pointer for plink2 pgen handle, or NULL.
#' @export
setClass("GenotypeHandle",
  representation(
    path = "character",
    format = "character",
    snpInfo = "data.frame",
    nSamples = "integer",
    sampleIds = "character",
    pgenPtr = "ANY"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@path) != 1L)
      errors <- c(errors, "'path' must be a single character string")
    valid_formats <- c("gds", "vcf", "plink1", "plink2")
    if (!object@format %in% valid_formats)
      errors <- c(errors, paste("'format' must be one of:",
                                paste(valid_formats, collapse = ", ")))
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# GWAS Summary Statistics
# =============================================================================

# NOTE: the old single-study GwasSumStats class has been removed; see the
# DFrame-subclass replacement at the bottom of this file (search for
# `setClass("GwasSumStats", contains = "DFrame", ...)`).

# =============================================================================
# Annotation Data
# =============================================================================

#' @title Genomic Annotation Matrix
#' @description Container for SNP-level annotations used in stratified
#'   heritability analysis. Supports binary (0/1) and continuous annotations.
#'   Annotations are classified as baseline (always jointly fitted) or
#'   candidate (evaluated via score statistics).
#' @slot snpRanges A \code{GRanges} object with one range per SNP,
#'   defining genomic positions.
#' @slot annotations A numeric matrix (SNPs x annotations). Dense for
#'   small annotation counts, can be sparse (\code{dgCMatrix}) for large
#'   binary annotation sets.
#' @slot annotationMeta A \code{data.frame} with columns:
#'   \describe{
#'     \item{name}{Character, annotation name}
#'     \item{tier}{Character, one of "baseline" or "candidate"}
#'     \item{type}{Character, one of "binary" or "continuous"}
#'   }
#' @slot genome Character string for genome build.
#' @export
setClass("AnnotationMatrix",
  representation(
    snpRanges = "GRanges",
    annotations = "ANY",  # matrix or dgCMatrix
    annotationMeta = "data.frame",
    genome = "character"
  ),
  validity = function(object) {
    errors <- character()
    n_snp <- length(object@snpRanges)
    n_annot <- ncol(object@annotations)
    if (nrow(object@annotations) != n_snp)
      errors <- c(errors,
        "Number of rows in 'annotations' must match length of 'snpRanges'")
    required_meta_cols <- c("name", "tier", "type")
    if (!all(required_meta_cols %in% colnames(object@annotationMeta)))
      errors <- c(errors,
        "annotationMeta must have columns: name, tier, type")
    if (nrow(object@annotationMeta) != n_annot)
      errors <- c(errors,
        "Number of rows in 'annotationMeta' must match annotation count")
    valid_tiers <- c("baseline", "candidate")
    if (!all(object@annotationMeta$tier %in% valid_tiers))
      errors <- c(errors,
        "annotationMeta$tier must be 'baseline' or 'candidate'")
    valid_types <- c("binary", "continuous")
    if (!all(object@annotationMeta$type %in% valid_types))
      errors <- c(errors,
        "annotationMeta$type must be 'binary' or 'continuous'")
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# LD Statistic Classes (Virtual + Concrete Subclasses)
# =============================================================================

#' @title LD Statistic (Virtual Base Class)
#' @description Abstract container for pre-computed LD statistics. Subclasses
#'   provide method-specific representations: eigendecompositions (for
#'   LDER/HDL/sHDL) and LD score matrices (for S-LDSC/g-LDSC).
#' @slot ldBlocks An \code{LdBlocks} object defining the block structure.
#' @slot snpInfo A \code{data.frame} with columns \code{SNP}, \code{CHR},
#'   \code{BP}, \code{A1}, \code{A2}, and optionally \code{MAF}.
#' @slot nRef Integer, sample size of the LD reference panel.
#' @slot inSample Logical, whether the LD reference is from the same
#'   cohort as the GWAS (affects bias correction).
#' @slot genome Character string for genome build.
#' @export
setClass("LdStatistic",
  contains = "VIRTUAL",
  representation(
    ldBlocks = "LdBlocks",
    snpInfo = "data.frame",
    nRef = "integer",
    inSample = "logical",
    genome = "character"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@nRef) != 1L || object@nRef <= 0L)
      errors <- c(errors, "'nRef' must be a single positive integer")
    if (length(object@inSample) != 1L)
      errors <- c(errors, "'inSample' must be a single logical value")
    if (length(object@genome) != 1L || !nzchar(object@genome))
      errors <- c(errors, "'genome' must be a single non-empty character string")
    if (nrow(object@snpInfo) == 0L)
      errors <- c(errors, "'snpInfo' must have at least one row")
    if (length(errors) == 0) TRUE else errors
  }
)

#' @title Eigendecomposition-Based LD Statistic
#' @description Pre-computed per-block eigendecompositions of the LD
#'   correlation matrix. Used by LDER, HDL, and sHDL.
#' @slot eigenList A list of length \code{nBlocks}, each element a list
#'   with components:
#'   \describe{
#'     \item{values}{Numeric vector of eigenvalues}
#'     \item{vectors}{Numeric matrix of eigenvectors (SNPs x retained components)}
#'     \item{snpIdx}{Integer vector of SNP indices in \code{snpInfo}}
#'   }
#' @slot eigenvalueTruncation Numeric, proportion of variance retained
#'   (e.g., 0.9 for HDL's default). If 1.0, no truncation.
#' @export
setClass("LdEigen",
  contains = "LdStatistic",
  representation(
    eigenList = "list",
    eigenvalueTruncation = "numeric"
  ),
  validity = function(object) {
    parent_check <- getValidity(getClass("LdStatistic"))(object)
    errors <- if (isTRUE(parent_check)) character() else parent_check
    nBlocks <- length(object@ldBlocks@blocks)
    if (length(object@eigenList) != nBlocks)
      errors <- c(errors,
        "Length of 'eigenList' must match number of LD blocks")
    if (length(object@eigenvalueTruncation) != 1L ||
        object@eigenvalueTruncation <= 0 ||
        object@eigenvalueTruncation > 1)
      errors <- c(errors,
        "'eigenvalueTruncation' must be a single value in (0, 1]")
    if (length(errors) == 0) TRUE else errors
  }
)

#' @title LD Score-Based LD Statistic
#' @description Pre-computed LD scores for each SNP. Used by S-LDSC and
#'   g-LDSC. Supports both standard LD scores and annotation-stratified
#'   LD scores.
#' @slot ldScores A numeric matrix (SNPs x annotations+1). The first
#'   column is the base LD score (sum of r^2). Additional columns are
#'   annotation-stratified LD scores if annotations are provided.
#' @slot ldScoreWeights A numeric vector of regression weights for each SNP.
#' @slot ldMatrixList For g-LDSC: a list of per-block LD (R^2) matrices
#'   used to compute the FGLS residual covariance. NULL for S-LDSC.
#' @export
setClass("LdScore",
  contains = "LdStatistic",
  representation(
    ldScores = "matrix",
    ldScoreWeights = "numeric",
    ldMatrixList = "list"  # for g-LDSC; empty list for S-LDSC
  ),
  validity = function(object) {
    parent_check <- getValidity(getClass("LdStatistic"))(object)
    errors <- if (isTRUE(parent_check)) character() else parent_check
    if (nrow(object@ldScores) != nrow(object@snpInfo))
      errors <- c(errors,
        "Number of rows in 'ldScores' must match 'snpInfo'")
    if (length(object@ldScoreWeights) != nrow(object@snpInfo))
      errors <- c(errors,
        "Length of 'ldScoreWeights' must match 'snpInfo'")
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Result Classes
# =============================================================================

#' @title Heritability Estimate
#' @description Container for univariate heritability estimation results.
#'   Holds global, local, and annotation-stratified estimates.
#' @slot h2 Numeric, global SNP heritability estimate.
#' @slot h2Se Numeric, standard error of global h2.
#' @slot intercept Numeric, confounding intercept estimate (NA if method
#'   does not estimate one).
#' @slot interceptSe Numeric, SE of intercept.
#' @slot local A \code{data.frame} with per-block local heritability
#'   estimates (columns: \code{blockId}, \code{h2Local}, \code{h2LocalSe}).
#'   NULL if \code{local = FALSE}.
#' @slot enrichment A \code{data.frame} with baseline annotation enrichment
#'   estimates (columns: \code{annotation}, \code{tau}, \code{tauSe},
#'   \code{enrichment}, \code{enrichmentSe}, \code{enrichmentP},
#'   \code{propH2}, \code{propSnps}). NULL if unstratified.
#' @slot tauBlocks A numeric matrix (nBlocks x n_annotations) of per-block
#'   jackknife tau values. Required for Gazal tauStar standardization
#'   downstream. NULL if not available (e.g., unstratified analysis).
#' @slot scoreStats A list with score statistics for candidate annotations,
#'   suitable for input to \code{susieRss}. Contains:
#'   \describe{
#'     \item{z}{Numeric vector of z-scores for each candidate annotation}
#'     \item{R}{Correlation matrix of the score statistics}
#'     \item{annotationNames}{Character vector of candidate annotation names}
#'   }
#'   NULL if no candidate annotations provided.
#' @slot method Character string identifying the estimation method.
#' @slot nSnps Integer, number of SNPs used in estimation.
#' @slot traitName Character string for trait identifier.
#' @export
setClass("H2Estimate",
  representation(
    h2 = "numeric",
    h2Se = "numeric",
    intercept = "numeric",
    interceptSe = "numeric",
    local = "ANY",        # data.frame or NULL
    enrichment = "ANY",   # data.frame or NULL
    tauBlocks = "ANY",    # matrix or NULL
    scoreStats = "ANY",   # list or NULL
    method = "character",
    nSnps = "integer",
    traitName = "character"
  )
)

# =============================================================================
# LD Data (fine-mapping / colocalization input)
# =============================================================================

#' @title LD Data Container
#' @description S4 container for LD information. Stores either a pre-computed
#'   correlation matrix or a \code{GenotypeHandle} (or list of handles for
#'   mixture panels) for lazy genotype/correlation access.
#'
#' @slot correlation A correlation matrix, a list of per-block matrices
#'   (block-diagonal LD), or NULL if genotypes are available and R should
#'   be computed on demand.
#' @slot genotypeHandle A \code{GenotypeHandle}, a list of
#'   \code{GenotypeHandle}s (for mixture panels), or NULL when only
#'   pre-computed R is available.
#' @slot snpIdx Integer vector of 1-based SNP indices into the handle's
#'   \code{snpInfo}. NULL when correlation is pre-computed.
#' @slot variants A \code{GRanges} object with variant metadata (A1, A2,
#'   variant_id, and optionally allele_freq, variance, n_nomiss).
#' @slot blockMetadata An \code{LdBlocks} object or a \code{data.frame}
#'   with block boundary information.
#' @slot nRef Integer, reference panel sample size.
#' @slot mixtureWeights NULL when \code{genotypeHandle} is a single
#'   \code{GenotypeHandle}; a numeric vector of mixing proportions (one
#'   per panel, summing to 1) when \code{genotypeHandle} is a list of
#'   \code{GenotypeHandle}s. Used by \code{getCorrelation()} to compute
#'   a weighted-average mixture LD matrix; required whenever
#'   \code{genotypeHandle} is a list and \code{getCorrelation()} will be
#'   called.
#' @export
setClass("LdData",
  representation(
    correlation = "ANY",       # matrix, list of matrices, or NULL
    genotypeHandle = "ANY",    # GenotypeHandle, list of GenotypeHandles, or NULL
    snpIdx = "ANY",            # integer or NULL
    variants = "GRanges",
    blockMetadata = "ANY",     # LdBlocks or data.frame
    nRef = "integer",
    mixtureWeights = "ANY"     # numeric vector or NULL
  ),
  validity = function(object) {
    errors <- character()
    if (is.null(object@correlation) && is.null(object@genotypeHandle))
      errors <- c(errors,
        "At least one of 'correlation' or 'genotypeHandle' must be non-NULL")
    if (length(object@variants) == 0)
      errors <- c(errors, "'variants' must not be empty")
    if (!is.null(object@mixtureWeights)) {
      if (!is.list(object@genotypeHandle))
        errors <- c(errors,
          "'mixtureWeights' may only be set when 'genotypeHandle' is a list of GenotypeHandles")
      else {
        w <- object@mixtureWeights
        if (!is.numeric(w) || length(w) != length(object@genotypeHandle))
          errors <- c(errors,
            "'mixtureWeights' must be numeric of length equal to the genotypeHandle list")
        else if (any(w < 0) || abs(sum(w) - 1) > 1e-6)
          errors <- c(errors,
            "'mixtureWeights' must be non-negative and sum to 1")
      }
    }
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Fine-Mapping Result
# =============================================================================

#' @title Fine-Mapping Result Base Class
#' @description Virtual base class for fine-mapping result collections.
#'   Concrete subclasses (\code{QtlFineMappingResult},
#'   \code{GwasFineMappingResult}) carry a \code{DFrame} of per-fit rows
#'   and a shared \code{ldSketch} slot. Downstream pipelines should
#'   dispatch on \code{FineMappingResultBase} for behaviors that apply to
#'   either flavour, and on the concrete subclass when the tuple shape
#'   matters.
#' @slot ldSketch The LD reference \code{GenotypeHandle} the fits were
#'   computed against, or \code{NULL} when the fits were derived from
#'   individual-level data (no LD reference). Used downstream for
#'   cross-pipeline LD-sketch identity validation.
#' @export
setClass("FineMappingResultBase",
  contains = c("VIRTUAL", "DFrame"),
  representation(ldSketch = "ANY"))

#' @title QTL Fine-Mapping Result Collection
#' @description S4 collection of fine-mapping fits keyed by the identity
#'   tuple \code{(study, context, trait, method)}. Each entry is a
#'   \code{FineMappingEntry}. Implements the \code{DFrame}-subclass
#'   collection pattern: subsetting with \code{[} and concatenation via
#'   \code{c()} are inherited from \code{S4Vectors::DataFrame}.
#'
#'   Required columns: \code{study}, \code{context}, \code{trait},
#'   \code{method}, \code{entry}. The 4-tuple is unique. Each
#'   \code{entry} is a \code{FineMappingEntry}.
#' @export
setClass("QtlFineMappingResult",
  contains = "FineMappingResultBase",
  validity = function(object) {
    errors <- character()
    required <- c("study", "context", "trait", "method", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L)
      errors <- c(errors, paste("missing columns:",
                                paste(missingCols, collapse = ", ")))
    if (length(errors) == 0L) {
      if (length(object$entry) != nrow(object))
        errors <- c(errors,
          "length(entry) must equal nrow(.) for QtlFineMappingResult")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "FineMappingEntry"),
                          logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a FineMappingEntry")
      keyDf <- as.data.frame(object[, c("study", "context", "trait", "method")])
      if (anyDuplicated(keyDf))
        errors <- c(errors,
          "(study, context, trait, method) tuple uniqueness violated")
    }
    if (!is.null(object@ldSketch) &&
        !methods::is(object@ldSketch, "GenotypeHandle")) {
      errors <- c(errors,
        "'ldSketch' must be a GenotypeHandle or NULL")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

#' @title GWAS Fine-Mapping Result Collection
#' @description S4 collection of fine-mapping fits for one or more GWAS
#'   studies on a single LD block. Keyed by the identity tuple
#'   \code{(study, method)}; each entry is a \code{FineMappingEntry}.
#'
#'   Required columns: \code{study}, \code{method}, \code{entry}. The
#'   2-tuple is unique. The caller is expected to construct one
#'   \code{GwasFineMappingResult} per LD block (no in-class block
#'   indexing).
#' @export
setClass("GwasFineMappingResult",
  contains = "FineMappingResultBase",
  validity = function(object) {
    errors <- character()
    required <- c("study", "method", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L)
      errors <- c(errors, paste("missing columns:",
                                paste(missingCols, collapse = ", ")))
    if (length(errors) == 0L) {
      if (length(object$entry) != nrow(object))
        errors <- c(errors,
          "length(entry) must equal nrow(.) for GwasFineMappingResult")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "FineMappingEntry"),
                          logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a FineMappingEntry")
      keyDf <- as.data.frame(object[, c("study", "method")])
      if (anyDuplicated(keyDf))
        errors <- c(errors,
          "(study, method) tuple uniqueness violated")
    }
    if (!is.null(object@ldSketch) &&
        !methods::is(object@ldSketch, "GenotypeHandle")) {
      errors <- c(errors,
        "'ldSketch' must be a GenotypeHandle or NULL")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

# =============================================================================
# TWAS Weights
# =============================================================================

#' @title TWAS Weights Collection
#' @description S4 collection of TWAS weights keyed by the identity tuple
#'   \code{(study, context, trait, method)}. Each entry is a
#'   \code{TwasWeightsEntry} carrying one method's weights for one
#'   trait/context/study. Implements the \code{DFrame}-subclass
#'   collection pattern.
#'
#'   Required columns: \code{study}, \code{context}, \code{trait},
#'   \code{method}, \code{entry}. The 4-tuple is unique. Each
#'   \code{entry} is a \code{TwasWeightsEntry}.
#' @slot ldSketch The LD reference \code{GenotypeHandle} the weights were
#'   derived against, or \code{NULL} when the weights were learned from
#'   individual-level data. Used downstream for cross-pipeline
#'   LD-sketch identity validation.
#' @export
setClass("TwasWeights",
  contains = "DFrame",
  representation(ldSketch = "ANY"),
  validity = function(object) {
    errors <- character()
    required <- c("study", "context", "trait", "method", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L)
      errors <- c(errors, paste("missing columns:",
                                paste(missingCols, collapse = ", ")))
    if (length(errors) == 0L) {
      if (length(object$entry) != nrow(object))
        errors <- c(errors,
          "length(entry) must equal nrow(.) for TwasWeights")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "TwasWeightsEntry"),
                          logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a TwasWeightsEntry")
      keyDf <- as.data.frame(object[, c("study", "context", "trait", "method")])
      if (anyDuplicated(keyDf))
        errors <- c(errors,
          "(study, context, trait, method) tuple uniqueness violated")
    }
    if (!is.null(object@ldSketch) &&
        !methods::is(object@ldSketch, "GenotypeHandle")) {
      errors <- c(errors,
        "'ldSketch' must be a GenotypeHandle or NULL")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)


# =============================================================================
# QTL Dataset
# =============================================================================

#' @title QTL Dataset (individual-level data for one study)
#' @description S4 container for a single QTL study's regional data. Holds
#'   a genotype handle plus per-context \code{SummarizedExperiment} objects
#'   carrying molecular-trait measurements. Each context's SE has
#'   \code{rowRanges} describing per-trait genomic positions and
#'   \code{colData} carrying per-context phenotype covariates. A single
#'   matrix of genotype-derived covariates (e.g., ancestry PCs) applies
#'   across contexts.
#'
#' @slot study Character (length 1). Study identifier; used in collection
#'   classes to tag downstream \code{FineMappingResult} / \code{TwasWeights}
#'   entries.
#' @slot genotypes A \code{GenotypeHandle} for lazy access to genotype
#'   dosages.
#' @slot phenotypes Named list of \code{SummarizedExperiment} objects, one
#'   per QTL context. Each SE has rows = molecular traits with positions
#'   in \code{rowRanges(se)}, columns = samples, and per-context covariates
#'   in \code{colData(se)}. Different contexts may carry different subsets
#'   of traits (rows); traits shared across contexts must have identical
#'   \code{rowRanges} entries (enforced by validity).
#' @slot genotypeCovariates Numeric matrix (samples x covariates) of
#'   genotype-derived covariates applied uniformly across all contexts
#'   (e.g., ancestry PCs).
#' @slot scaleResiduals Logical (length 1). Whether residualization
#'   accessors scale residuals to unit variance.
#' @slot mafCutoff Numeric (length 1). Minor allele frequency threshold;
#'   variants with \code{MAF < mafCutoff} are dropped at extraction time
#'   inside \code{getGenotypes()} / \code{getResidualizedGenotypes()}.
#'   Default 0 (no filter).
#' @slot macCutoff Numeric (length 1). Minor allele count threshold;
#'   converted to a MAF threshold using
#'   \code{max(mafCutoff, macCutoff / (2 * n))} where \code{n} is the
#'   post-narrowing sample count of the extracted block. Default 0
#'   (no filter).
#' @slot xvarCutoff Numeric (length 1). Per-variant genotype variance
#'   threshold; variants with column variance below this are dropped at
#'   extraction time. Default 0 (no filter).
#' @slot imissCutoff Numeric (length 1). Per-sample genotype-missingness
#'   threshold; samples with a missing-genotype rate above this are
#'   dropped at extraction time. Default 0 (no filter).
#' @slot keepSamples Character vector of sample identifiers to retain
#'   prior to per-block QC; intersected with the genotype handle's
#'   \code{sampleIds} and the \code{samples} argument of
#'   \code{getGenotypes()}. Length 0 means no restriction.
#' @slot keepVariants Character vector of variant identifiers to retain
#'   prior to per-block QC. Length 0 means no restriction.
#' @export
setClass("QtlDataset",
  representation(
    study              = "character",
    genotypes          = "GenotypeHandle",
    phenotypes         = "list",
    genotypeCovariates = "matrix",
    scaleResiduals     = "logical",
    mafCutoff          = "numeric",
    macCutoff          = "numeric",
    xvarCutoff         = "numeric",
    imissCutoff        = "numeric",
    keepSamples        = "character",
    keepVariants       = "character"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@study) != 1L || !nzchar(object@study))
      errors <- c(errors, "'study' must be a single non-empty character string")
    if (length(object@scaleResiduals) != 1L)
      errors <- c(errors, "'scaleResiduals' must be a single logical value")
    for (nm in c("mafCutoff", "macCutoff", "xvarCutoff", "imissCutoff")) {
      v <- methods::slot(object, nm)
      if (length(v) != 1L || is.na(v) || !is.finite(v) || v < 0)
        errors <- c(errors, sprintf(
          "'%s' must be a single finite non-negative numeric", nm))
    }
    if (length(object@phenotypes) == 0L)
      errors <- c(errors, "'phenotypes' must not be empty")
    contextNames <- names(object@phenotypes)
    if (is.null(contextNames) || any(!nzchar(contextNames)) ||
        any(is.na(contextNames)))
      errors <- c(errors, "'phenotypes' must be a named list with non-empty names")
    else if (anyDuplicated(contextNames))
      errors <- c(errors, "context names in 'phenotypes' must be unique")
    for (ctx in seq_along(object@phenotypes)) {
      se <- object@phenotypes[[ctx]]
      if (!methods::is(se, "SummarizedExperiment")) {
        errors <- c(errors, sprintf(
          "phenotypes[[%d]] must be a SummarizedExperiment (got %s)",
          ctx, class(se)[[1L]]))
      }
    }
    # Trait-position consistency across shared traits in different contexts.
    if (length(object@phenotypes) > 1L &&
        all(vapply(object@phenotypes, methods::is, logical(1),
                   "SummarizedExperiment"))) {
      traitToRange <- list()
      for (ctx in seq_along(object@phenotypes)) {
        se <- object@phenotypes[[ctx]]
        rr <- SummarizedExperiment::rowRanges(se)
        ids <- rownames(se)
        if (length(rr) != length(ids)) next
        for (i in seq_along(ids)) {
          tid <- ids[[i]]
          prev <- traitToRange[[tid]]
          this <- rr[i]
          if (is.null(prev)) {
            traitToRange[[tid]] <- this
          } else {
            if (!isTRUE(all.equal(
                  as.character(GenomicRanges::seqnames(prev)),
                  as.character(GenomicRanges::seqnames(this))
                )) ||
                GenomicRanges::start(prev) != GenomicRanges::start(this) ||
                GenomicRanges::end(prev) != GenomicRanges::end(this)) {
              errors <- c(errors, sprintf(
                "trait '%s' has inconsistent rowRanges across contexts", tid))
            }
          }
        }
      }
    }
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Per-entry payload classes for FineMappingResult and TwasWeights
# =============================================================================

#' @title Fine-Mapping Entry (per-tuple payload)
#' @description S4 container for a single fine-mapping fit attached to a
#'   \code{FineMappingResult} row. One entry corresponds to one
#'   \code{(study, context, trait, method)} tuple.
#'
#'   For joint fits (e.g., multi-trait mvSuSiE or fSuSiE), multiple
#'   \code{FineMappingEntry} objects in the same \code{FineMappingResult}
#'   collection may carry references to the same underlying R fit object
#'   (R's copy-on-modify semantics keep this memory-efficient).
#' @slot variantIds Character vector of variant IDs in the fit.
#' @slot trimmedFit The method-specific fit object (SuSiE list, mvSuSiE
#'   object, fSuSiE object, etc.). May be shared by reference across
#'   joint-fit entries.
#' @slot topLoci A long-format \code{data.frame} with at minimum
#'   \code{variant_id} and \code{pip} columns; optional \code{cs},
#'   \code{coverage}, \code{betahat}, \code{sd}, \code{csLog10bf},
#'   \code{z}.
#' @slot sumstats A list of summary statistics used in the fit, or
#'   \code{NULL}.
#' @export
setClass("FineMappingEntry",
  representation(
    variantIds = "character",
    trimmedFit = "ANY",
    topLoci    = "data.frame",
    sumstats   = "ANY"
  ),
  validity = function(object) {
    errors <- character()
    if (nrow(object@topLoci) > 0L) {
      required <- c("variant_id", "pip")
      missingCols <- setdiff(required, colnames(object@topLoci))
      if (length(missingCols) > 0L)
        errors <- c(errors, paste("topLoci missing columns:",
                                  paste(missingCols, collapse = ", ")))
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

#' @title TWAS Weights Entry (per-tuple payload)
#' @description S4 container for one method's TWAS weights, attached to
#'   a \code{TwasWeights} row. One entry corresponds to one
#'   \code{(study, context, trait, method)} tuple.
#' @slot variantIds Character vector of variant IDs that have weights.
#' @slot weights Numeric vector (single-method, single-outcome) or
#'   matrix (multi-outcome).
#' @slot fits Optional method-specific fit object.
#' @slot cvPerformance Optional named list of CV metrics (\code{rsq},
#'   \code{pval}, etc.).
#' @slot standardized Logical (length 1). Whether the weights are on the
#'   standardized scale.
#' @slot dataType Data-type tag for downstream usage (e.g.,
#'   \code{"expression"}, \code{"splicing"}); may be \code{NULL}.
#' @export
setClass("TwasWeightsEntry",
  representation(
    variantIds    = "character",
    weights       = "ANY",
    fits          = "ANY",
    cvPerformance = "ANY",
    standardized  = "logical",
    dataType      = "ANY"
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@standardized) != 1L)
      errors <- c(errors, "'standardized' must be a single logical")
    if (is.matrix(object@weights) &&
        nrow(object@weights) != length(object@variantIds))
      errors <- c(errors,
        "nrow(weights) must equal length(variantIds)")
    if (length(errors) == 0L) TRUE else errors
  }
)

# =============================================================================
# MultiTaskQtlDataset
# =============================================================================

#' @title Multi-Task QTL Dataset
#' @description S4 container for a multi-task QTL analysis: a collection
#'   of individual-level \code{QtlDataset} studies, optionally combined
#'   with a \code{QtlSumStats} collection of summary-statistic-only
#'   studies. Used as the input to multi-study fine-mapping and
#'   colocboost-style analyses.
#'
#'   At least two studies must be present in total, counting
#'   \code{qtlDatasets} entries plus the studies in \code{sumStats}.
#'
#'   For traits that appear in more than one \code{qtlDatasets} entry,
#'   the per-trait genomic positions must agree across the entries
#'   (enforced by validity). No cross-checking is performed against
#'   \code{sumStats} variant positions, since summary statistics are
#'   already computed and cannot be re-aligned at construction time.
#' @slot qtlDatasets A named list of \code{QtlDataset} objects, keyed by
#'   study identifier.
#' @slot sumStats An optional \code{QtlSumStats} carrying additional
#'   summary-statistic-only studies. \code{NULL} when absent.
#' @export
setClass("MultiTaskQtlDataset",
  representation(
    qtlDatasets = "list",
    sumStats    = "ANY"
  ),
  validity = function(object) {
    errors <- character()
    # qtlDatasets elements
    if (!is.list(object@qtlDatasets) || length(object@qtlDatasets) == 0L) {
      errors <- c(errors, "'qtlDatasets' must be a non-empty named list")
    } else {
      nm <- names(object@qtlDatasets)
      if (is.null(nm) || any(!nzchar(nm)) || any(is.na(nm)))
        errors <- c(errors, "'qtlDatasets' must be a named list with non-empty names")
      else if (anyDuplicated(nm))
        errors <- c(errors, "names of 'qtlDatasets' must be unique")
      bad <- !vapply(object@qtlDatasets,
                    function(d) methods::is(d, "QtlDataset"), logical(1))
      if (any(bad))
        errors <- c(errors,
          "every element of 'qtlDatasets' must be a QtlDataset")
    }
    # sumStats: either NULL or a QtlSumStats
    if (!is.null(object@sumStats) && !methods::is(object@sumStats, "QtlSumStats")) {
      errors <- c(errors,
        "'sumStats' must be a QtlSumStats object or NULL")
    }
    # Total study count >= 2
    nQtl <- length(object@qtlDatasets)
    nSumstats <- if (is.null(object@sumStats)) 0L
                 else length(unique(as.character(object@sumStats$study)))
    if (length(errors) == 0L && (nQtl + nSumstats) < 2L) {
      errors <- c(errors, sprintf(
        "MultiTaskQtlDataset requires at least 2 studies in total (got ",
        "%d individual-level + %d summary-statistic = %d).",
        nQtl, nSumstats, nQtl + nSumstats))
    }
    # Pairwise trait-position consistency across qtlDatasets
    if (length(errors) == 0L && nQtl >= 2L) {
      traitRanges <- lapply(object@qtlDatasets, function(qd) {
        out <- list()
        for (ctx in seq_along(qd@phenotypes)) {
          se <- qd@phenotypes[[ctx]]
          rr <- SummarizedExperiment::rowRanges(se)
          ids <- rownames(se)
          for (i in seq_along(ids)) {
            tid <- ids[[i]]
            if (is.null(out[[tid]])) out[[tid]] <- rr[i]
          }
        }
        out
      })
      pairs <- utils::combn(seq_along(traitRanges), 2L)
      for (k in seq_len(ncol(pairs))) {
        i <- pairs[1L, k]; j <- pairs[2L, k]
        a <- traitRanges[[i]]; b <- traitRanges[[j]]
        shared <- intersect(names(a), names(b))
        for (tid in shared) {
          if (!isTRUE(all.equal(
                as.character(GenomicRanges::seqnames(a[[tid]])),
                as.character(GenomicRanges::seqnames(b[[tid]])))) ||
              GenomicRanges::start(a[[tid]]) != GenomicRanges::start(b[[tid]]) ||
              GenomicRanges::end(a[[tid]]) != GenomicRanges::end(b[[tid]])) {
            errors <- c(errors, sprintf(
              "trait '%s' has inconsistent rowRanges between studies '%s' and '%s'",
              tid, names(object@qtlDatasets)[i],
              names(object@qtlDatasets)[j]))
          }
        }
      }
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

# =============================================================================
# Summary Statistics collection classes (post-refactor; DFrame subclasses)
# =============================================================================

#' @title QTL Summary Statistics Collection
#' @description S4 collection of QTL summary statistics keyed by the
#'   identity tuple \code{(study, context, trait)}. Each entry holds a
#'   \code{GRanges} of summary statistics for that tuple. Class-level
#'   slots \code{ldSketch} (the LD reference \code{GenotypeHandle}) and
#'   \code{genome} (the genome build, a single character string) apply
#'   to every entry; the genome build must be uniform because all
#'   entries necessarily share the LD reference.
#'
#'   Required columns: \code{study}, \code{context}, \code{trait},
#'   \code{entry}. Optional columns include \code{varY} (numeric,
#'   per-tuple phenotype variance; \code{NA_real_} when unused). The
#'   3-tuple \code{(study, context, trait)} is unique. Each \code{entry}
#'   is a \code{GRanges} whose mcols carry the per-variant statistics
#'   (\code{SNP}, \code{A1}, \code{A2}, \code{Z}, \code{N}; plus
#'   optional \code{MAF}, \code{INFO}, \code{BETA}, \code{SE}, \code{P}).
#' @slot ldSketch A \code{GenotypeHandle} carrying the LD reference for
#'   downstream QC and RSS analysis.
#' @slot genome A single character string giving the genome build that
#'   the LD sketch and every entry are aligned to.
#' @title Summary Statistics Base Class
#' @description Virtual base class for summary-statistic collections.
#'   Concrete subclasses (\code{QtlSumStats}, \code{GwasSumStats}) carry
#'   a \code{DFrame} of per-entry rows plus shared slots \code{ldSketch},
#'   \code{genome}, and \code{qcInfo}. Downstream pipelines should
#'   dispatch on \code{SumStatsBase} for behaviors that apply to either
#'   flavour and on the concrete subclass when the tuple shape matters.
#' @slot ldSketch A \code{GenotypeHandle} carrying the LD reference for
#'   downstream QC and RSS analysis.
#' @slot genome A single character string giving the genome build that
#'   the LD sketch and every entry are aligned to.
#' @slot qcInfo A \code{list} recording which QC steps ran. Empty
#'   \code{list()} on construction; populated by \code{summaryStatsQc()}.
#'   Fine-mapping and TWAS-weights pipelines reject inputs where
#'   \code{length(getQcInfo(x)) == 0L} — the slot serves as both the
#'   gating flag and the audit trail.
#' @export
setClass("SumStatsBase",
  contains = c("VIRTUAL", "DFrame"),
  representation(
    ldSketch = "GenotypeHandle",
    genome   = "character",
    qcInfo   = "list"
  ))

#' @slot qcInfo A \code{list} recording which QC steps ran. Empty
#'   \code{list()} on construction; populated by \code{summaryStatsQc()}
#'   with a per-step audit record (filter names, drop counts, liftover
#'   target, RAISS settings, etc.). Fine-mapping and TWAS-weights
#'   pipelines reject inputs where \code{length(getQcInfo(x)) == 0L} — the
#'   slot serves as both the gating flag and the audit trail.
#' @export
setClass("QtlSumStats",
  contains = "SumStatsBase",
  validity = function(object) {
    errors <- character()
    required <- c("study", "context", "trait", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L)
      errors <- c(errors, paste("missing columns:",
                                paste(missingCols, collapse = ", ")))
    if (length(object@genome) != 1L || !nzchar(object@genome))
      errors <- c(errors,
        "'genome' slot must be a single non-empty character string")
    if (!is.list(object@qcInfo))
      errors <- c(errors, "'qcInfo' slot must be a list")
    if (length(errors) == 0L) {
      if (length(object$entry) != nrow(object))
        errors <- c(errors,
          "length(entry) must equal nrow(.) for QtlSumStats")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "GRanges"), logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a GRanges")
      keyDf <- as.data.frame(object[, c("study", "context", "trait")])
      if (anyDuplicated(keyDf))
        errors <- c(errors,
          "(study, context, trait) tuple uniqueness violated")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

#' @title GWAS Summary Statistics Collection
#' @description S4 collection of GWAS summary statistics keyed by
#'   \code{study}. Each entry holds a \code{GRanges} of summary statistics
#'   for that study. Class-level slots \code{ldSketch} (the LD reference
#'   \code{GenotypeHandle}) and \code{genome} (the genome build, a single
#'   character string) apply to every entry; the genome build must be
#'   uniform because all entries necessarily share the LD reference.
#'
#'   Required columns: \code{study}, \code{entry}. Optional columns
#'   include \code{varY} (numeric, phenotype variance for the
#'   sufficient-statistic interface; NA otherwise). \code{study} is
#'   unique. Each \code{entry} is a \code{GRanges} whose mcols carry the
#'   per-variant statistics (\code{SNP}, \code{A1}, \code{A2}, \code{Z},
#'   \code{N}; plus optional \code{MAF}, \code{INFO}, \code{BETA},
#'   \code{SE}, \code{P}).
#' @slot ldSketch A \code{GenotypeHandle} for the LD reference.
#' @slot genome A single character string giving the genome build that
#'   the LD sketch and every entry are aligned to.
#' @slot qcInfo A \code{list} recording which QC steps ran. Empty
#'   \code{list()} on construction; populated by \code{summaryStatsQc()}.
#'   Fine-mapping and TWAS-weights pipelines reject inputs where
#'   \code{length(getQcInfo(x)) == 0L} — the slot serves as both the
#'   gating flag and the audit trail.
#' @export
setClass("GwasSumStats",
  contains = "SumStatsBase",
  validity = function(object) {
    errors <- character()
    required <- c("study", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L)
      errors <- c(errors, paste("missing columns:",
                                paste(missingCols, collapse = ", ")))
    if (length(object@genome) != 1L || !nzchar(object@genome))
      errors <- c(errors,
        "'genome' slot must be a single non-empty character string")
    if (!is.list(object@qcInfo))
      errors <- c(errors, "'qcInfo' slot must be a list")
    if (length(errors) == 0L) {
      if (length(object$entry) != nrow(object))
        errors <- c(errors,
          "length(entry) must equal nrow(.) for GwasSumStats")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "GRanges"), logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a GRanges")
      if (anyDuplicated(as.character(object$study)))
        errors <- c(errors, "`study` must be unique")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)

# =============================================================================
# Show Methods
# =============================================================================

#' @export
setMethod("show", "GenotypeHandle", function(object) {
  cat(sprintf("GenotypeHandle [%s]\n", object@format))
  cat(sprintf("  Path: %s\n", object@path))
  cat(sprintf("  %d samples, %d SNPs\n",
              object@nSamples, nrow(object@snpInfo)))
})

# NOTE: the old single-study GwasSumStats show method has been removed;
# the new collection class's show method is defined alongside its class
# definition (see QtlSumStats / GwasSumStats section above).

#' @export
setMethod("show", "LdBlocks", function(object) {
  cat(sprintf("LdBlocks: %d blocks, genome build: %s\n",
              length(object@blocks), object@genome))
})

#' @export
setMethod("show", "AnnotationMatrix", function(object) {
  n_base <- sum(object@annotationMeta$tier == "baseline")
  n_cand <- sum(object@annotationMeta$tier == "candidate")
  n_bin <- sum(object@annotationMeta$type == "binary")
  n_cont <- sum(object@annotationMeta$type == "continuous")
  cat(sprintf("AnnotationMatrix: %d SNPs x %d annotations\n",
              nrow(object@annotations), ncol(object@annotations)))
  cat(sprintf("  Baseline: %d, Candidate: %d\n", n_base, n_cand))
  cat(sprintf("  Binary: %d, Continuous: %d\n", n_bin, n_cont))
  cat(sprintf("  Genome build: %s\n", object@genome))
})

#' @export
setMethod("show", "LdEigen", function(object) {
  cat(sprintf("LdEigen: %d SNPs across %d blocks\n",
              nrow(object@snpInfo), length(object@eigenList)))
  cat(sprintf("  Eigenvalue truncation: %.2f\n",
              object@eigenvalueTruncation))
  cat(sprintf("  Reference N: %d, In-sample: %s\n",
              object@nRef, object@inSample))
})

#' @export
setMethod("show", "LdScore", function(object) {
  n_scores <- ncol(object@ldScores)
  has_matrix <- length(object@ldMatrixList) > 0
  cat(sprintf("LdScore: %d SNPs, %d LD score columns\n",
              nrow(object@snpInfo), n_scores))
  cat(sprintf("  Full LD matrices: %s (needed for g-LDSC)\n", has_matrix))
  cat(sprintf("  Reference N: %d, In-sample: %s\n",
              object@nRef, object@inSample))
})

#' @export
setMethod("show", "H2Estimate", function(object) {
  cat(sprintf("H2Estimate for '%s' (method: %s)\n",
              object@traitName, object@method))
  cat(sprintf("  h2 = %.4f (SE = %.4f)\n", object@h2, object@h2Se))
  if (!is.na(object@intercept))
    cat(sprintf("  intercept = %.4f (SE = %.4f)\n",
                object@intercept, object@interceptSe))
  has_local <- !is.null(object@local)
  has_enrich <- !is.null(object@enrichment)
  has_tau_blocks <- !is.null(object@tauBlocks)
  cat(sprintf("  Local: %s, Enrichment: %s, tauBlocks: %s\n",
              has_local, has_enrich, has_tau_blocks))
  cat(sprintf("  N SNPs: %d\n", object@nSnps))
})

#' @export
setMethod("show", "LdData", function(object) {
  n_var <- length(object@variants)
  has_R <- !is.null(object@correlation)
  has_geno <- !is.null(object@genotypeHandle)
  r_type <- if (has_R && is.list(object@correlation)) "block-diagonal" else "single"
  cat(sprintf("LdData: %d variants\n", n_var))
  cat(sprintf("  Correlation: %s, Genotype handle: %s\n",
              if (has_R) r_type else "NULL",
              if (has_geno) "available" else "NULL"))
  cat(sprintf("  Reference N: %d\n", object@nRef))
})

#' @export
setMethod("show", "QtlFineMappingResult", function(object) {
  cat(sprintf("QtlFineMappingResult: %d entries\n", nrow(object)))
  if (nrow(object) > 0L) {
    cat(sprintf("  %d studies, %d contexts, %d traits, %d methods\n",
                length(unique(object$study)),
                length(unique(object$context)),
                length(unique(object$trait)),
                length(unique(object$method))))
  }
  ldSrc <- if (is.null(object@ldSketch)) "NULL (individual-level fit)"
           else sprintf("%s @ %s",
                         object@ldSketch@format,
                         object@ldSketch@path)
  cat(sprintf("  LD sketch: %s\n", ldSrc))
})

#' @export
setMethod("show", "GwasFineMappingResult", function(object) {
  cat(sprintf("GwasFineMappingResult: %d entries\n", nrow(object)))
  if (nrow(object) > 0L) {
    cat(sprintf("  %d studies, %d methods\n",
                length(unique(object$study)),
                length(unique(object$method))))
  }
  ldSrc <- if (is.null(object@ldSketch)) "NULL"
           else sprintf("%s @ %s",
                         object@ldSketch@format,
                         object@ldSketch@path)
  cat(sprintf("  LD sketch: %s\n", ldSrc))
})

#' @export
setMethod("show", "TwasWeights", function(object) {
  cat(sprintf("TwasWeights: %d entries\n", nrow(object)))
  if (nrow(object) > 0L) {
    cat(sprintf("  %d studies, %d contexts, %d traits, %d methods\n",
                length(unique(object$study)),
                length(unique(object$context)),
                length(unique(object$trait)),
                length(unique(object$method))))
  }
  ldSrc <- if (is.null(object@ldSketch)) "NULL (individual-level fit)"
           else sprintf("%s @ %s",
                         object@ldSketch@format,
                         object@ldSketch@path)
  cat(sprintf("  LD sketch: %s\n", ldSrc))
})

setMethod("show", "GwasSumStats", function(object) {
  cat(sprintf("GwasSumStats: %d studies, genome build %s\n",
              nrow(object), object@genome))
  cat(sprintf("  LD sketch: %s @ %s\n",
              object@ldSketch@format, object@ldSketch@path))
})

#' @export
setMethod("show", "FineMappingEntry", function(object) {
  nCs <- if (nrow(object@topLoci) > 0L && "cs" %in% names(object@topLoci))
           length(unique(object@topLoci$cs[object@topLoci$cs > 0]))
         else 0L
  cat(sprintf("FineMappingEntry: %d variants, %d credible sets\n",
              length(object@variantIds), nCs))
})

#' @export
setMethod("show", "TwasWeightsEntry", function(object) {
  cat(sprintf("TwasWeightsEntry: %d variants, standardized=%s\n",
              length(object@variantIds), object@standardized))
  hasCv <- !is.null(object@cvPerformance)
  cat(sprintf("  CV performance: %s\n", hasCv))
})

#' @export
setMethod("show", "MultiTaskQtlDataset", function(object) {
  nQtl <- length(object@qtlDatasets)
  ssEntries <- if (is.null(object@sumStats)) 0L
               else length(unique(as.character(object@sumStats$study)))
  cat(sprintf("MultiTaskQtlDataset: %d individual-level + %d sumstats studies\n",
              nQtl, ssEntries))
  if (nQtl > 0L) {
    cat(sprintf("  Individual-level studies: %s\n",
                paste(names(object@qtlDatasets), collapse = ", ")))
  }
  if (!is.null(object@sumStats)) {
    cat(sprintf("  Sumstats studies: %s\n",
                paste(unique(as.character(object@sumStats$study)),
                      collapse = ", ")))
  }
})

#' @export
setMethod("show", "QtlDataset", function(object) {
  nCtx <- length(object@phenotypes)
  ctxNames <- names(object@phenotypes)
  totalTraits <- length(unique(unlist(
    lapply(object@phenotypes, rownames), use.names = FALSE)))
  cat(sprintf("QtlDataset for study '%s'\n", object@study))
  cat(sprintf("  %d context(s): %s\n", nCtx,
              paste(ctxNames, collapse = ", ")))
  cat(sprintf("  %d unique traits across contexts\n", totalTraits))
  cat(sprintf("  Genotypes: %s\n",
              paste0(object@genotypes@format, " @ ", object@genotypes@path)))
  cat(sprintf("  Genotype covariates: %d cols\n",
              ncol(object@genotypeCovariates)))
  cat(sprintf("  Scale residuals: %s\n", object@scaleResiduals))
})

