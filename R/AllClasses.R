#' @title S4 Class Definitions
#' @description All S4 class definitions for pecotmr: genotype handles,
#'   summary statistics, LD references, annotations, fine-mapping results,
#'   TWAS weights, regional data, and heritability estimates.
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
setClass("LDBlocks",
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
#' @slot snp_info A \code{data.frame} with columns \code{SNP}, \code{CHR},
#'   \code{BP}, \code{A1}, \code{A2}. Cached on first access.
#' @slot n_samples Integer, number of samples.
#' @slot sample_ids Character vector of sample identifiers.
#' @slot pgen_ptr An external pointer for plink2 pgen handle, or NULL.
#' @export
setClass("GenotypeHandle",
  representation(
    path = "character",
    format = "character",
    snp_info = "data.frame",
    n_samples = "integer",
    sample_ids = "character",
    pgen_ptr = "ANY"
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
#' @slot trait_name Character string for trait identifier.
#' @slot var_y Numeric, phenotype variance. For case-control studies this is
#'   \code{1 / (phi * (1 - phi))} where \code{phi = n_case / n}. For
#'   quantitative traits, typically 1 or NULL.
#' @export
setClass("GWASSumStats",
  representation(
    sumstats = "GRanges",
    genome = "character",
    trait_name = "character",
    var_y = "ANY"  # numeric or NULL
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
#' @slot snp_ranges A \code{GRanges} object with one range per SNP,
#'   defining genomic positions.
#' @slot annotations A numeric matrix (SNPs x annotations). Dense for
#'   small annotation counts, can be sparse (\code{dgCMatrix}) for large
#'   binary annotation sets.
#' @slot annotation_meta A \code{data.frame} with columns:
#'   \describe{
#'     \item{name}{Character, annotation name}
#'     \item{tier}{Character, one of "baseline" or "candidate"}
#'     \item{type}{Character, one of "binary" or "continuous"}
#'   }
#' @slot genome Character string for genome build.
#' @export
setClass("AnnotationMatrix",
  representation(
    snp_ranges = "GRanges",
    annotations = "ANY",  # matrix or dgCMatrix
    annotation_meta = "data.frame",
    genome = "character"
  ),
  validity = function(object) {
    errors <- character()
    n_snp <- length(object@snp_ranges)
    n_annot <- ncol(object@annotations)
    if (nrow(object@annotations) != n_snp)
      errors <- c(errors,
        "Number of rows in 'annotations' must match length of 'snp_ranges'")
    required_meta_cols <- c("name", "tier", "type")
    if (!all(required_meta_cols %in% colnames(object@annotation_meta)))
      errors <- c(errors,
        "annotation_meta must have columns: name, tier, type")
    if (nrow(object@annotation_meta) != n_annot)
      errors <- c(errors,
        "Number of rows in 'annotation_meta' must match annotation count")
    valid_tiers <- c("baseline", "candidate")
    if (!all(object@annotation_meta$tier %in% valid_tiers))
      errors <- c(errors,
        "annotation_meta$tier must be 'baseline' or 'candidate'")
    valid_types <- c("binary", "continuous")
    if (!all(object@annotation_meta$type %in% valid_types))
      errors <- c(errors,
        "annotation_meta$type must be 'binary' or 'continuous'")
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
#' @slot ld_blocks An \code{LDBlocks} object defining the block structure.
#' @slot snp_info A \code{data.frame} with columns \code{SNP}, \code{CHR},
#'   \code{BP}, \code{A1}, \code{A2}, and optionally \code{MAF}.
#' @slot n_ref Integer, sample size of the LD reference panel.
#' @slot in_sample Logical, whether the LD reference is from the same
#'   cohort as the GWAS (affects bias correction).
#' @slot genome Character string for genome build.
#' @export
setClass("LDStatistic",
  contains = "VIRTUAL",
  representation(
    ld_blocks = "LDBlocks",
    snp_info = "data.frame",
    n_ref = "integer",
    in_sample = "logical",
    genome = "character"
  )
)

#' @title Eigendecomposition-Based LD Statistic
#' @description Pre-computed per-block eigendecompositions of the LD
#'   correlation matrix. Used by LDER, HDL, and sHDL.
#' @slot eigen_list A list of length \code{n_blocks}, each element a list
#'   with components:
#'   \describe{
#'     \item{values}{Numeric vector of eigenvalues}
#'     \item{vectors}{Numeric matrix of eigenvectors (SNPs x retained components)}
#'     \item{snp_idx}{Integer vector of SNP indices in \code{snp_info}}
#'   }
#' @slot eigenvalue_truncation Numeric, proportion of variance retained
#'   (e.g., 0.9 for HDL's default). If 1.0, no truncation.
#' @export
setClass("LDEigen",
  contains = "LDStatistic",
  representation(
    eigen_list = "list",
    eigenvalue_truncation = "numeric"
  ),
  validity = function(object) {
    errors <- character()
    n_blocks <- length(object@ld_blocks@blocks)
    if (length(object@eigen_list) != n_blocks)
      errors <- c(errors,
        "Length of 'eigen_list' must match number of LD blocks")
    if (length(object@eigenvalue_truncation) != 1L ||
        object@eigenvalue_truncation <= 0 ||
        object@eigenvalue_truncation > 1)
      errors <- c(errors,
        "'eigenvalue_truncation' must be a single value in (0, 1]")
    if (length(errors) == 0) TRUE else errors
  }
)

#' @title LD Score-Based LD Statistic
#' @description Pre-computed LD scores for each SNP. Used by S-LDSC and
#'   g-LDSC. Supports both standard LD scores and annotation-stratified
#'   LD scores.
#' @slot ld_scores A numeric matrix (SNPs x annotations+1). The first
#'   column is the base LD score (sum of r^2). Additional columns are
#'   annotation-stratified LD scores if annotations are provided.
#' @slot ld_score_weights A numeric vector of regression weights for each SNP.
#' @slot ld_matrix_list For g-LDSC: a list of per-block LD (R^2) matrices
#'   used to compute the FGLS residual covariance. NULL for S-LDSC.
#' @export
setClass("LDScore",
  contains = "LDStatistic",
  representation(
    ld_scores = "matrix",
    ld_score_weights = "numeric",
    ld_matrix_list = "list"  # for g-LDSC; empty list for S-LDSC
  ),
  validity = function(object) {
    errors <- character()
    if (nrow(object@ld_scores) != nrow(object@snp_info))
      errors <- c(errors,
        "Number of rows in 'ld_scores' must match 'snp_info'")
    if (length(object@ld_score_weights) != nrow(object@snp_info))
      errors <- c(errors,
        "Length of 'ld_score_weights' must match 'snp_info'")
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
#' @slot h2_se Numeric, standard error of global h2.
#' @slot intercept Numeric, confounding intercept estimate (NA if method
#'   does not estimate one).
#' @slot intercept_se Numeric, SE of intercept.
#' @slot local A \code{data.frame} with per-block local heritability
#'   estimates (columns: \code{block_id}, \code{h2_local}, \code{h2_local_se}).
#'   NULL if \code{local = FALSE}.
#' @slot enrichment A \code{data.frame} with baseline annotation enrichment
#'   estimates (columns: \code{annotation}, \code{tau}, \code{tau_se},
#'   \code{enrichment}, \code{enrichment_se}, \code{enrichment_p},
#'   \code{prop_h2}, \code{prop_snps}). NULL if unstratified.
#' @slot tau_blocks A numeric matrix (n_blocks x n_annotations) of per-block
#'   jackknife tau values. Required for Gazal tau_star standardization
#'   downstream. NULL if not available (e.g., unstratified analysis).
#' @slot score_stats A list with score statistics for candidate annotations,
#'   suitable for input to \code{susie_rss}. Contains:
#'   \describe{
#'     \item{z}{Numeric vector of z-scores for each candidate annotation}
#'     \item{R}{Correlation matrix of the score statistics}
#'     \item{annotation_names}{Character vector of candidate annotation names}
#'   }
#'   NULL if no candidate annotations provided.
#' @slot method Character string identifying the estimation method.
#' @slot n_snps Integer, number of SNPs used in estimation.
#' @slot trait_name Character string for trait identifier.
#' @export
setClass("H2Estimate",
  representation(
    h2 = "numeric",
    h2_se = "numeric",
    intercept = "numeric",
    intercept_se = "numeric",
    local = "ANY",        # data.frame or NULL
    enrichment = "ANY",   # data.frame or NULL
    tau_blocks = "ANY",   # matrix or NULL
    score_stats = "ANY",  # list or NULL
    method = "character",
    n_snps = "integer",
    trait_name = "character"
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
#' @slot genotype_handle A \code{GenotypeHandle}, a list of
#'   \code{GenotypeHandle}s (for mixture panels), or NULL when only
#'   pre-computed R is available.
#' @slot snp_idx Integer vector of 1-based SNP indices into the handle's
#'   \code{snp_info}. NULL when correlation is pre-computed.
#' @slot variants A \code{GRanges} object with variant metadata (A1, A2,
#'   variant_id, and optionally allele_freq, variance, n_nomiss).
#' @slot block_metadata An \code{LDBlocks} object or a \code{data.frame}
#'   with block boundary information.
#' @slot n_ref Integer, reference panel sample size.
#' @export
setClass("LDData",
  representation(
    correlation = "ANY",       # matrix, list of matrices, or NULL
    genotype_handle = "ANY",   # GenotypeHandle, list of GenotypeHandles, or NULL
    snp_idx = "ANY",           # integer or NULL
    variants = "GRanges",
    block_metadata = "ANY",    # LDBlocks or data.frame
    n_ref = "integer"
  ),
  validity = function(object) {
    errors <- character()
    if (is.null(object@correlation) && is.null(object@genotype_handle))
      errors <- c(errors,
        "At least one of 'correlation' or 'genotype_handle' must be non-NULL")
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
#' @slot variant_names Character vector of variant IDs.
#' @slot trimmed_fit List containing the method-specific trimmed fit.
#' @slot top_loci A \code{data.frame} in long format with columns:
#'   variant_id, method, coverage, cs, pip, and optionally betahat, sd,
#'   cs_log10bf, z.
#' @slot method Character string identifying the fine-mapping method.
#' @slot sumstats List of summary statistics used in the analysis, or NULL.
#' @export
setClass("FineMappingResult",
  representation(
    variant_names = "character",
    trimmed_fit = "ANY",
    top_loci = "data.frame",
    method = "character",
    sumstats = "ANY"  # list or NULL
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@method) != 1L)
      errors <- c(errors, "'method' must be a single character string")
    if (nrow(object@top_loci) > 0) {
      required <- c("variant_id", "method")
      missing_cols <- setdiff(required, colnames(object@top_loci))
      if (length(missing_cols) > 0)
        errors <- c(errors, paste("top_loci missing columns:",
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
#' @slot variant_ids Character vector of variant IDs (row names for all
#'   weight matrices).
#' @slot methods Character vector of method names (names of the weights
#'   list).
#' @slot fits Named list of model fit objects, or NULL.
#' @slot cv_performance Named list of cross-validation performance
#'   metrics, or NULL.
#' @export
setClass("TWASWeights",
  representation(
    weights = "list",
    variant_ids = "character",
    methods = "character",
    fits = "ANY",           # list or NULL
    cv_performance = "ANY"  # list or NULL
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@methods) != length(object@weights))
      errors <- c(errors,
        "Length of 'methods' must match length of 'weights'")
    for (i in seq_along(object@weights)) {
      w <- object@weights[[i]]
      if (!is.null(w) && is.matrix(w) && nrow(w) != length(object@variant_ids))
        errors <- c(errors, paste0(
          "Weight matrix '", object@methods[i],
          "' has ", nrow(w), " rows but variant_ids has length ",
          length(object@variant_ids)))
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
#' @slot genotype_matrix Numeric matrix (samples x variants) of genotype
#'   dosages, with colnames as variant IDs and rownames as sample IDs.
#' @slot phenotypes Named list of phenotype matrices (per condition).
#' @slot covariates Named list of covariate matrices (per condition).
#' @slot scale_residuals Logical, whether to scale residuals.
#' @slot maf Named list of MAF vectors (per condition).
#' @slot region A \code{GRanges} (single range) or NULL.
#' @slot dropped_samples Named list of dropped sample vectors.
#' @slot Y_coordinates Phenotype coordinates, or NULL.
#' @export
setClass("RegionalData",
  representation(
    genotype_matrix = "matrix",
    phenotypes = "list",
    covariates = "list",
    scale_residuals = "logical",
    maf = "list",
    region = "ANY",           # GRanges or NULL
    dropped_samples = "list",
    Y_coordinates = "ANY"     # data.frame or NULL
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
# Show Methods
# =============================================================================

#' @export
setMethod("show", "GenotypeHandle", function(object) {
  cat(sprintf("GenotypeHandle [%s]\n", object@format))
  cat(sprintf("  Path: %s\n", object@path))
  cat(sprintf("  %d samples, %d SNPs\n",
              object@n_samples, nrow(object@snp_info)))
})

#' @export
setMethod("show", "GWASSumStats", function(object) {
  cat(sprintf("GWASSumStats for '%s'\n", object@trait_name))
  cat(sprintf("  %d SNPs, genome build: %s\n",
              length(object@sumstats), object@genome))
  cat(sprintf("  Median N: %.0f\n",
              stats::median(S4Vectors::mcols(object@sumstats)$N)))
  has_maf <- "MAF" %in% colnames(S4Vectors::mcols(object@sumstats))
  cat(sprintf("  MAF available: %s\n", has_maf))
  if (!is.null(object@var_y))
    cat(sprintf("  var_y: %.4f\n", object@var_y))
})

#' @export
setMethod("show", "LDBlocks", function(object) {
  cat(sprintf("LDBlocks: %d blocks, genome build: %s\n",
              length(object@blocks), object@genome))
})

#' @export
setMethod("show", "AnnotationMatrix", function(object) {
  n_base <- sum(object@annotation_meta$tier == "baseline")
  n_cand <- sum(object@annotation_meta$tier == "candidate")
  n_bin <- sum(object@annotation_meta$type == "binary")
  n_cont <- sum(object@annotation_meta$type == "continuous")
  cat(sprintf("AnnotationMatrix: %d SNPs x %d annotations\n",
              nrow(object@annotations), ncol(object@annotations)))
  cat(sprintf("  Baseline: %d, Candidate: %d\n", n_base, n_cand))
  cat(sprintf("  Binary: %d, Continuous: %d\n", n_bin, n_cont))
  cat(sprintf("  Genome build: %s\n", object@genome))
})

#' @export
setMethod("show", "LDEigen", function(object) {
  cat(sprintf("LDEigen: %d SNPs across %d blocks\n",
              nrow(object@snp_info), length(object@eigen_list)))
  cat(sprintf("  Eigenvalue truncation: %.2f\n",
              object@eigenvalue_truncation))
  cat(sprintf("  Reference N: %d, In-sample: %s\n",
              object@n_ref, object@in_sample))
})

#' @export
setMethod("show", "LDScore", function(object) {
  n_scores <- ncol(object@ld_scores)
  has_matrix <- length(object@ld_matrix_list) > 0
  cat(sprintf("LDScore: %d SNPs, %d LD score columns\n",
              nrow(object@snp_info), n_scores))
  cat(sprintf("  Full LD matrices: %s (needed for g-LDSC)\n", has_matrix))
  cat(sprintf("  Reference N: %d, In-sample: %s\n",
              object@n_ref, object@in_sample))
})

#' @export
setMethod("show", "H2Estimate", function(object) {
  cat(sprintf("H2Estimate for '%s' (method: %s)\n",
              object@trait_name, object@method))
  cat(sprintf("  h2 = %.4f (SE = %.4f)\n", object@h2, object@h2_se))
  if (!is.na(object@intercept))
    cat(sprintf("  intercept = %.4f (SE = %.4f)\n",
                object@intercept, object@intercept_se))
  has_local <- !is.null(object@local)
  has_enrich <- !is.null(object@enrichment)
  has_tau_blocks <- !is.null(object@tau_blocks)
  cat(sprintf("  Local: %s, Enrichment: %s, tau_blocks: %s\n",
              has_local, has_enrich, has_tau_blocks))
  cat(sprintf("  N SNPs: %d\n", object@n_snps))
})

#' @export
setMethod("show", "LDData", function(object) {
  n_var <- length(object@variants)
  has_R <- !is.null(object@correlation)
  has_geno <- !is.null(object@genotype_handle)
  r_type <- if (has_R && is.list(object@correlation)) "block-diagonal" else "single"
  cat(sprintf("LDData: %d variants\n", n_var))
  cat(sprintf("  Correlation: %s, Genotype handle: %s\n",
              if (has_R) r_type else "NULL",
              if (has_geno) "available" else "NULL"))
  cat(sprintf("  Reference N: %d\n", object@n_ref))
})

#' @export
setMethod("show", "FineMappingResult", function(object) {
  n_cs <- if (nrow(object@top_loci) > 0 && "cs" %in% names(object@top_loci))
    length(unique(object@top_loci$cs[object@top_loci$cs > 0])) else 0L
  cat(sprintf("FineMappingResult [%s]: %d variants, %d credible sets\n",
              object@method, length(object@variant_names), n_cs))
})

#' @export
setMethod("show", "TWASWeights", function(object) {
  cat(sprintf("TWASWeights: %d methods, %d variants\n",
              length(object@methods), length(object@variant_ids)))
  cat(sprintf("  Methods: %s\n", paste(object@methods, collapse = ", ")))
  has_cv <- !is.null(object@cv_performance)
  cat(sprintf("  CV performance: %s\n", has_cv))
})

#' @export
setMethod("show", "RegionalData", function(object) {
  n_cond <- length(object@phenotypes)
  n_var <- ncol(object@genotype_matrix)
  n_samp <- nrow(object@genotype_matrix)
  cat(sprintf("RegionalData: %d conditions, %d variants, %d samples\n",
              n_cond, n_var, n_samp))
  cat(sprintf("  Scale residuals: %s\n", object@scale_residuals))
})
