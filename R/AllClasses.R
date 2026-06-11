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
#'   Metadata columns may include \code{block_id} (integer).
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

#' @title GWAS Summary Statistics
#' @description Standardized container for GWAS summary statistics. SNP
#'   positions are stored as a \code{GRanges} object, with effect sizes,
#'   standard errors, and sample sizes as metadata columns.
#' @slot sumstats A \code{GRanges} object with required metadata columns:
#'   \describe{
#'     \item{SNP}{Character, SNP identifier (rsID or chr:pos:ref:alt)}
#'     \item{A1}{Character, effect allele}
#'     \item{A2}{Character, non-effect allele}
#'     \item{Z}{Numeric, z-score (BETA/SE)}
#'     \item{N}{Numeric, sample size (per-SNP or constant)}
#'   }
#'   Optional metadata columns: \code{MAF}, \code{INFO}, \code{BETA},
#'   \code{SE}, \code{P}.
#' @slot genome Character string for genome build.
#' @slot traitName Character string for trait identifier.
#' @slot varY Numeric, phenotype variance. For observed-scale OLS on a
#'   centered 0/1 case-control trait, this is the \code{susieR}
#'   \code{sum(y^2) / (n - 1)} value after centering,
#'   \code{n / (n - 1) * phi * (1 - phi)}, where \code{phi = nCase / n}.
#'   Use it only with the full \code{bhat/shat/var_y} sufficient-statistic
#'   interface; z-score RSS analyses should leave it NULL.
#' @export
setClass("GwasSumStats",
  representation(
    sumstats = "GRanges",
    genome = "character",
    traitName = "character",
    varY = "ANY"  # numeric or NULL
  ),
  validity = function(object) {
    errors <- character()
    required_cols <- c("SNP", "A1", "A2", "Z", "N")
    mcols_names <- colnames(mcols(object@sumstats))
    missing <- setdiff(required_cols, mcols_names)
    if (length(missing) > 0)
      errors <- c(errors, paste("Missing required columns:",
                                paste(missing, collapse = ", ")))
    if (length(object@genome) != 1L)
      errors <- c(errors, "'genome' must be a single character string")
    if (length(errors) == 0) TRUE else errors
  }
)

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
#' @slot eigenList A list of length \code{n_blocks}, each element a list
#'   with components:
#'   \describe{
#'     \item{values}{Numeric vector of eigenvalues}
#'     \item{vectors}{Numeric matrix of eigenvectors (SNPs x retained components)}
#'     \item{snp_idx}{Integer vector of SNP indices in \code{snpInfo}}
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
    n_blocks <- length(object@ldBlocks@blocks)
    if (length(object@eigenList) != n_blocks)
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
#'   estimates (columns: \code{block_id}, \code{h2_local}, \code{h2_local_se}).
#'   NULL if \code{local = FALSE}.
#' @slot enrichment A \code{data.frame} with baseline annotation enrichment
#'   estimates (columns: \code{annotation}, \code{tau}, \code{tau_se},
#'   \code{enrichment}, \code{enrichment_se}, \code{enrichment_p},
#'   \code{prop_h2}, \code{prop_snps}). NULL if unstratified.
#' @slot tauBlocks A numeric matrix (n_blocks x n_annotations) of per-block
#'   jackknife tau values. Required for Gazal tau_star standardization
#'   downstream. NULL if not available (e.g., unstratified analysis).
#' @slot scoreStats A list with score statistics for candidate annotations,
#'   suitable for input to \code{susie_rss}. Contains:
#'   \describe{
#'     \item{z}{Numeric vector of z-scores for each candidate annotation}
#'     \item{R}{Correlation matrix of the score statistics}
#'     \item{annotation_names}{Character vector of candidate annotation names}
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
#' @export
setClass("LdData",
  representation(
    correlation = "ANY",       # matrix, list of matrices, or NULL
    genotypeHandle = "ANY",    # GenotypeHandle, list of GenotypeHandles, or NULL
    snpIdx = "ANY",            # integer or NULL
    variants = "GRanges",
    blockMetadata = "ANY",     # LdBlocks or data.frame
    nRef = "integer"
  ),
  validity = function(object) {
    errors <- character()
    if (is.null(object@correlation) && is.null(object@genotypeHandle))
      errors <- c(errors,
        "At least one of 'correlation' or 'genotypeHandle' must be non-NULL")
    if (length(object@variants) == 0)
      errors <- c(errors, "'variants' must not be empty")
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Fine-Mapping Result
# =============================================================================

#' @title Fine-Mapping Result
#' @description S4 container for fine-mapping output. Stores variant names,
#'   the trimmed model fit, and a long-format table of credible sets and PIPs.
#' @slot variantNames Character vector of variant IDs.
#' @slot trimmedFit List containing the method-specific trimmed fit.
#' @slot topLoci A \code{data.frame} in long format with columns:
#'   variant_id, method, coverage, cs, pip, and optionally betahat, sd,
#'   cs_log10bf, z.
#' @slot method Character string identifying the fine-mapping method.
#' @slot sumstats List of summary statistics used in the analysis, or NULL.
#' @export
setClass("FineMappingResult",
  representation(
    variantNames = "character",
    trimmedFit = "ANY",
    topLoci = "data.frame",
    method = "character",
    sumstats = "ANY"  # list or NULL
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@method) != 1L)
      errors <- c(errors, "'method' must be a single character string")
    if (nrow(object@topLoci) > 0) {
      required <- c("variant_id", "method")
      missing_cols <- setdiff(required, colnames(object@topLoci))
      if (length(missing_cols) > 0)
        errors <- c(errors, paste("topLoci missing columns:",
                                  paste(missing_cols, collapse = ", ")))
    }
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# TWAS Weights
# =============================================================================

#' @title TWAS Weights
#' @description S4 container for TWAS weight matrices.
#' @slot weights Named list of numeric matrices (variants x outcomes).
#' @slot variantIds Character vector of variant IDs (row names for all
#'   weight matrices).
#' @slot methods Character vector of method names (names of the weights
#'   list).
#' @slot fits Named list of model fit objects, or NULL.
#' @slot cvPerformance Named list of cross-validation performance
#'   metrics, or NULL.
#' @slot standardized Logical, whether weights are on standardized
#'   (correlation) scale. If TRUE, \code{harmonize_twas} skips the
#'   \code{sqrt(variance)} scaling step. Individual-level weights use
#'   FALSE (raw genotype scale); RSS weights use TRUE.
#' @export
setClass("TwasWeights",
  representation(
    weights = "list",
    variantIds = "character",
    methods = "character",
    fits = "ANY",           # list or NULL
    cvPerformance = "ANY",  # list or NULL
    standardized = "logical",
    molecularId = "character",  # gene/molecule name (length 0 or 1)
    dataType = "ANY"            # named list of data types per context, or NULL
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@standardized) != 1L)
      errors <- c(errors, "'standardized' must be a single logical value")
    if (length(object@methods) != length(object@weights))
      errors <- c(errors,
        "Length of 'methods' must match length of 'weights'")
    for (i in seq_along(object@weights)) {
      w <- object@weights[[i]]
      if (!is.null(w) && is.matrix(w) && nrow(w) != length(object@variantIds))
        errors <- c(errors, paste0(
          "Weight matrix '", object@methods[i],
          "' has ", nrow(w), " rows but variantIds has length ",
          length(object@variantIds)))
    }
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Allele QC Result
# =============================================================================

#' @title Allele QC Result
#' @description S4 container for the output of \code{matchRefPanel} /
#'   \code{alleleQc}. Carries the post-QC target variants alongside the full
#'   merge / flip / strand diagnostics needed by downstream callers that
#'   inspect what QC did.
#' @slot harmonizedData A \code{data.frame} of variants retained after
#'   allele harmonization, with reference-aligned A1/A2 and (when requested)
#'   sign-flipped effect columns.
#' @slot qcSummary A \code{data.frame} carrying per-variant QC diagnostics
#'   from the full merge: \code{variants_id_original}, \code{variants_id_qced},
#'   \code{exact_match}, \code{sign_flip}, \code{strand_flip}, \code{INDEL},
#'   \code{ID_match}, \code{keep}, etc.
#' @export
setClass("AlleleQcResult",
  representation(
    harmonizedData = "data.frame",
    qcSummary = "data.frame"
  )
)

# =============================================================================
# Summary-Statistics QC Result
# =============================================================================

#' @title Summary-Statistics QC Result
#' @description S4 container holding the output of \code{summaryStatsQc} and
#'   \code{.summary_stats_qc_single_study}. Carries the post-QC LD reference
#'   plus harmonized sumstats, a pre-imputation snapshot, and QC process
#'   metadata. Replaces the legacy list-of-named-fields return shape.
#' @slot ldData An \code{LdData} S4 object containing the post-QC LD
#'   reference (correlation and/or genotype), or NULL when QC produced no LD.
#' @slot rssInput List with \code{sumstats} (post-QC data.frame), \code{n},
#'   and \code{varY}.
#' @slot preprocess List with \code{sumstats} and \code{ldData} fields
#'   capturing the pre-imputation snapshot for downstream re-runs.
#' @slot outlierNumber Integer count of LD-mismatch outliers removed.
#' @slot skipped Single logical; TRUE when QC short-circuited.
#' @slot skipReason Character string explaining a skip; empty otherwise.
#' @export
setClass("QcResult",
  representation(
    ldData = "ANY",                   # LdData or NULL
    rssInput = "list",
    preprocess = "list",
    outlierNumber = "integer",
    skipped = "logical",
    skipReason = "character"
  ),
  validity = function(object) {
    errors <- character()
    if (!is.null(object@ldData) && !is(object@ldData, "LdData"))
      errors <- c(errors, "'ldData' must be an LdData object or NULL")
    if (length(object@skipped) != 1L)
      errors <- c(errors, "'skipped' must be a single logical value")
    if (length(object@outlierNumber) != 1L)
      errors <- c(errors, "'outlierNumber' must be a single integer")
    if (length(object@skipReason) > 1L)
      errors <- c(errors, "'skipReason' must be a single character string (or empty)")
    if (length(object@rssInput) > 0L) {
      required <- c("sumstats", "n", "varY")
      missing_keys <- setdiff(required, names(object@rssInput))
      if (length(missing_keys) > 0L)
        errors <- c(errors, paste0(
          "'rssInput' is missing key(s): ", paste(missing_keys, collapse = ", ")))
      if (!is.null(object@rssInput$sumstats) &&
          !is.data.frame(object@rssInput$sumstats))
        errors <- c(errors,
          "'rssInput$sumstats' must be a data.frame")
    }
    if (length(object@preprocess) > 0L) {
      pp_keys <- names(object@preprocess)
      if (!all(pp_keys %in% c("sumstats", "ldData")))
        errors <- c(errors,
          "'preprocess' may only contain 'sumstats' and 'ldData' keys")
      if (!is.null(object@preprocess$ldData) &&
          !is(object@preprocess$ldData, "LdData"))
        errors <- c(errors,
          "'preprocess$ldData' must be an LdData or NULL")
    }
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Regional Data (pipeline input)
# =============================================================================

#' @title Regional Association Data
#' @description S4 container for regional genotype/phenotype/covariate data.
#'   Residualized genotypes and phenotypes are computed lazily via accessors.
#' @slot genotypeMatrix Numeric matrix (samples x variants) of genotype
#'   dosages, with colnames as variant IDs and rownames as sample IDs.
#' @slot phenotypes Named list of phenotype matrices (per condition).
#' @slot covariates Named list of covariate matrices (per condition).
#' @slot scaleResiduals Logical, whether to scale residuals.
#' @slot maf Named list of MAF vectors (per condition).
#' @slot region A \code{GRanges} (single range) or NULL.
#' @slot droppedSamples Named list of dropped sample vectors.
#' @slot coordinates Phenotype coordinates, or NULL.
#' @export
setClass("RegionalData",
  representation(
    genotypeMatrix = "matrix",
    phenotypes = "list",
    covariates = "list",
    scaleResiduals = "logical",
    maf = "list",
    region = "ANY",           # GRanges or NULL
    droppedSamples = "list",
    coordinates = "ANY"       # data.frame or NULL
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@phenotypes) == 0)
      errors <- c(errors, "'phenotypes' must not be empty")
    if (length(object@covariates) != length(object@phenotypes))
      errors <- c(errors,
        "'covariates' and 'phenotypes' must have the same length")
    if (length(errors) == 0) TRUE else errors
  }
)

# =============================================================================
# Multivariate Regional Data
# =============================================================================

#' @title Multivariate Regional Association Data
#' @description S4 container for regional association data prepared for
#'   multivariate (joint-across-conditions) modeling. Unlike
#'   \code{RegionalData}, which carries a per-condition list of phenotype
#'   matrices, this class assumes all conditions are jointly observed in the
#'   same samples and packs the phenotypes into a single multivariate matrix
#'   (samples x conditions).
#' @slot genotypeMatrix Numeric matrix (samples x variants), rownames are
#'   sample IDs, colnames are variant IDs.
#' @slot Y Numeric matrix (samples x conditions) of residualized
#'   phenotypes after joining conditions and (optionally) filtering rows by
#'   minimum non-missing count.
#' @slot scaling Numeric vector of per-condition scaling factors
#'   (length = ncol(Y)).
#' @slot droppedSamples Character or list capturing sample IDs dropped
#'   during multivariate filtering.
#' @slot region A \code{GRanges} (single range) or NULL.
#' @slot coordinates A data.frame of phenotype coordinates, or NULL.
#' @export
setClass("MultivariateRegionalData",
  representation(
    genotypeMatrix = "matrix",
    Y = "matrix",
    scaling = "numeric",
    droppedSamples = "ANY",
    region = "ANY",
    coordinates = "ANY"
  ),
  validity = function(object) {
    errors <- character()
    if (nrow(object@genotypeMatrix) != nrow(object@Y))
      errors <- c(errors,
        "genotypeMatrix and Y must have the same number of rows")
    if (length(object@scaling) != ncol(object@Y))
      errors <- c(errors, "length(scaling) must equal ncol(Y)")
    if (length(errors) == 0) TRUE else errors
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

#' @export
setMethod("show", "GwasSumStats", function(object) {
  cat(sprintf("GwasSumStats for '%s'\n", object@traitName))
  cat(sprintf("  %d SNPs, genome build: %s\n",
              length(object@sumstats), object@genome))
  cat(sprintf("  Median N: %.0f\n",
              stats::median(S4Vectors::mcols(object@sumstats)$N)))
  has_maf <- "MAF" %in% colnames(S4Vectors::mcols(object@sumstats))
  cat(sprintf("  MAF available: %s\n", has_maf))
  if (!is.null(object@varY))
    cat(sprintf("  varY: %.4f\n", object@varY))
})

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
setMethod("show", "FineMappingResult", function(object) {
  n_cs <- if (nrow(object@topLoci) > 0 && "cs" %in% names(object@topLoci))
    length(unique(object@topLoci$cs[object@topLoci$cs > 0])) else 0L
  cat(sprintf("FineMappingResult [%s]: %d variants, %d credible sets\n",
              object@method, length(object@variantNames), n_cs))
})

#' @export
setMethod("show", "TwasWeights", function(object) {
  cat(sprintf("TwasWeights: %d methods, %d variants\n",
              length(object@methods), length(object@variantIds)))
  if (length(object@molecularId) > 0)
    cat(sprintf("  Molecular ID: %s\n", object@molecularId))
  cat(sprintf("  Methods: %s\n", paste(object@methods, collapse = ", ")))
  cat(sprintf("  Standardized: %s\n", object@standardized))
  has_cv <- !is.null(object@cvPerformance)
  cat(sprintf("  CV performance: %s\n", has_cv))
})

#' @export
setMethod("show", "RegionalData", function(object) {
  n_cond <- length(object@phenotypes)
  n_var <- ncol(object@genotypeMatrix)
  n_samp <- nrow(object@genotypeMatrix)
  cat(sprintf("RegionalData: %d conditions, %d variants, %d samples\n",
              n_cond, n_var, n_samp))
  cat(sprintf("  Scale residuals: %s\n", object@scaleResiduals))
})

#' @export
setMethod("show", "MultivariateRegionalData", function(object) {
  cat(sprintf("MultivariateRegionalData: %d conditions, %d variants, %d samples\n",
              ncol(object@Y), ncol(object@genotypeMatrix),
              nrow(object@genotypeMatrix)))
  if (!is.null(object@region))
    cat(sprintf("  Region: %s:%d-%d\n",
                as.character(GenomicRanges::seqnames(object@region))[1],
                GenomicRanges::start(object@region),
                GenomicRanges::end(object@region)))
})

#' @export
setMethod("show", "AlleleQcResult", function(object) {
  cat(sprintf("AlleleQcResult: %d harmonized variants (from %d scanned)\n",
              nrow(object@harmonizedData), nrow(object@qcSummary)))
})

#' @export
setMethod("show", "QcResult", function(object) {
  cat(sprintf("QcResult: %s\n",
              if (object@skipped) sprintf("skipped (%s)", object@skipReason) else "completed"))
  if (length(object@rssInput) > 0 && !is.null(object@rssInput$sumstats)) {
    cat(sprintf("  Sumstats: %d variants\n",
                nrow(object@rssInput$sumstats)))
  }
  if (!is.null(object@ldData)) {
    cat(sprintf("  LD: %d variants%s\n",
                length(getVariantIds(object@ldData)),
                if (hasGenotypes(object@ldData)) " (genotype-backed)" else " (correlation)"))
  }
  cat(sprintf("  Outliers removed: %d\n", object@outlierNumber))
})
