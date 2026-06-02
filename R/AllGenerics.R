#' @title S4 Generic Function Definitions
#' @description All S4 generic function definitions for pecotmr.
#' @include AllClasses.R
#' @importFrom methods setGeneric
NULL

# =============================================================================
# High-level estimation generic
# =============================================================================

#' @title Estimate SNP Heritability
#' @description Estimate SNP heritability from GWAS summary statistics using
#'   one of three methods: LDER, g-LDSC, or HDL/sHDL.
#' @param sumstats A \code{GWASSumStats} object.
#' @param ld_ref An \code{LDStatistic} object (method-appropriate subclass).
#' @param method Character, one of "lder", "gldsc", "hdl".
#' @param annotations An \code{AnnotationMatrix} object, or NULL for
#'   unstratified estimation.
#' @param local Logical, whether to compute per-block local estimates.
#' @param ... Additional method-specific arguments.
#' @return An \code{H2Estimate} object.
#' @export
setGeneric("estimateH2",
  function(sumstats, ld_ref, method = "lder", annotations = NULL,
           local = FALSE, ...)
    standardGeneric("estimateH2")
)

# =============================================================================
# LD score computation
# =============================================================================

#' @title Compute LD Scores
#' @description Compute LD scores from an LD reference, optionally
#'   stratified by annotations.
#' @param ld_ref An \code{LDStatistic} object.
#' @param annotations An \code{AnnotationMatrix} object, or NULL.
#' @param ... Additional arguments.
#' @return A numeric matrix of LD scores (SNPs x annotations+1).
#' @export
setGeneric("computeLdScores",
  function(ld_ref, annotations = NULL, ...)
    standardGeneric("computeLdScores")
)

# =============================================================================
# I/O generics
# =============================================================================

#' @title Read Genotype Data
#' @description Read genotype data from various formats (VCF, plink1,
#'   plink2, GDS) and return a \code{GenotypeHandle} for deferred
#'   genotype loading.
#' @param path Character, path to the genotype file.
#' @param format Character, one of "vcf", "plink1", "plink2", "gds".
#'   If NULL, inferred from file extension.
#' @param ... Additional arguments.
#' @return A \code{GenotypeHandle} object.
#' @export
setGeneric("readGenotypes",
  function(path, format = NULL, ...)
    standardGeneric("readGenotypes")
)

#' @title Read GWAS Summary Statistics
#' @description Read and standardize GWAS summary statistics from file.
#'   Optionally uses MungeSumStats for format detection and QC.
#' @param path Character, path to the summary statistics file.
#' @param trait_name Character, name for the trait.
#' @param genome Character, genome build (e.g., "hg19", "hg38").
#'   If NULL, inferred by MungeSumStats.
#' @param n Numeric, sample size (if not in the file).
#' @param use_mungesumstats Logical, whether to use MungeSumStats for
#'   standardization. Default TRUE if available.
#' @param ... Additional arguments.
#' @return A \code{GWASSumStats} object.
#' @export
setGeneric("readSumstats",
  function(path, trait_name = "trait", genome = NULL, n = NULL,
           use_mungesumstats = TRUE, ...)
    standardGeneric("readSumstats")
)

#' @title Read Annotations
#' @description Read genomic annotations from files (BED, BigWig,
#'   S-LDSC .annot format, or GRanges objects) and create an
#'   AnnotationMatrix.
#' @param paths Named character vector of file paths, or a named list
#'   of GRanges objects. Names become annotation names.
#' @param snp_ranges A \code{GRanges} object defining SNP positions.
#' @param annotation_meta A \code{data.frame} with annotation metadata
#'   (name, tier, type). If NULL, auto-detected from file format.
#' @param genome Character, genome build.
#' @param ... Additional arguments.
#' @return An \code{AnnotationMatrix} object.
#' @export
setGeneric("readAnnotations",
  function(paths, snp_ranges, annotation_meta = NULL,
           genome = "hg19", ...)
    standardGeneric("readAnnotations")
)

# =============================================================================
# Accessor generics
# =============================================================================

#' @title Get Local Estimates
#' @description Extract per-block local estimates from a result object.
#' @param object An \code{H2Estimate} object.
#' @return A \code{data.frame} of local estimates, or NULL.
#' @export
setGeneric("getLocal", function(object) standardGeneric("getLocal"))

#' @title Get Enrichment Estimates
#' @description Extract annotation enrichment estimates from a result object.
#' @param object An \code{H2Estimate} object.
#' @return A \code{data.frame} of enrichment estimates, or NULL.
#' @export
setGeneric("getEnrichment",
  function(object) standardGeneric("getEnrichment"))

#' @title Get Score Statistics
#' @description Extract score statistics for candidate annotations.
#' @param object An \code{H2Estimate} object.
#' @return A list with \code{z} and \code{R}, or NULL.
#' @export
setGeneric("getScoreStats",
  function(object) standardGeneric("getScoreStats"))

# =============================================================================
# GWASSumStats accessor generics
# =============================================================================

#' @title Get Z-scores
#' @description Extract z-score vector from a GWASSumStats object.
#' @param x A \code{GWASSumStats} object.
#' @return Numeric vector of z-scores.
#' @export
setGeneric("getZ", function(x) standardGeneric("getZ"))

#' @title Get Sample Sizes
#' @description Extract sample size vector from a GWASSumStats object.
#' @param x A \code{GWASSumStats} object.
#' @return Numeric vector of sample sizes.
#' @export
setGeneric("getN", function(x) standardGeneric("getN"))

#' @title Get Minor Allele Frequencies
#' @description Extract MAF vector from a GWASSumStats object.
#' @param x A \code{GWASSumStats} object.
#' @return Numeric vector of MAFs, or NULL if not available.
#' @export
setGeneric("getMaf", function(x) standardGeneric("getMaf"))

#' @title Get Number of SNPs
#' @description Get the number of SNPs in a GWASSumStats object.
#' @param x A \code{GWASSumStats} object.
#' @return Integer.
#' @export
setGeneric("nSnps", function(x) standardGeneric("nSnps"))

#' @title Subset by Chromosome
#' @description Extract a chromosome-specific subset of a GWASSumStats object.
#' @param x A \code{GWASSumStats} object.
#' @param chr Character, chromosome name (e.g., "1", "chr1").
#' @return A \code{GWASSumStats} object.
#' @export
setGeneric("subsetChr", function(x, chr) standardGeneric("subsetChr"))

#' @title Get Phenotype Variance
#' @description Extract phenotype variance from a GWASSumStats object.
#' @param x A \code{GWASSumStats} object.
#' @return Numeric phenotype variance, or NULL.
#' @export
setGeneric("getVarY", function(x) standardGeneric("getVarY"))

# =============================================================================
# LDData accessor generics
# =============================================================================

#' @title Get LD Correlation Matrix
#' @description Extract the LD correlation matrix from an \code{LDData}
#'   object. If only a genotype handle is available, recomputes R from
#'   genotypes on the fly.
#' @param x An \code{LDData} object.
#' @return A correlation matrix, or a list of per-block matrices.
#' @export
setGeneric("getCorrelation", function(x) standardGeneric("getCorrelation"))

#' @title Get Genotype Matrix
#' @description Extract genotype matrix from an \code{LDData} object via
#'   the genotype handle. Returns NULL if no handle is available.
#' @param x An \code{LDData} object.
#' @return A numeric matrix (samples x variants), a list of matrices
#'   for mixture panels, or NULL.
#' @export
setGeneric("getGenotypes", function(x) standardGeneric("getGenotypes"))

#' @title Check Genotype Availability
#' @description Check whether an \code{LDData} object has a genotype
#'   handle for extracting raw genotypes.
#' @param x An \code{LDData} object.
#' @return Logical.
#' @export
setGeneric("hasGenotypes", function(x) standardGeneric("hasGenotypes"))

#' @title Get Variant IDs
#' @description Extract variant ID vector from an \code{LDData} object.
#' @param x An \code{LDData} object.
#' @return Character vector of variant IDs.
#' @export
setGeneric("getVariantIds", function(x) standardGeneric("getVariantIds"))

#' @title Get Variant GRanges
#' @description Extract variant metadata as GRanges.
#' @param x An object with variant metadata.
#' @return A \code{GRanges} object.
#' @export
setGeneric("getVariantInfo", function(x) standardGeneric("getVariantInfo"))

#' @title Get Reference Panel
#' @description Extract reference panel metadata as a data.frame from
#'   an \code{LDData} object, including chrom and pos columns.
#' @param x An \code{LDData} object.
#' @return A data.frame with variant metadata including chrom, pos, A1, A2.
#' @export
setGeneric("getRefPanel", function(x) standardGeneric("getRefPanel"))

#' @title Get Block Metadata
#' @description Extract block metadata from an \code{LDData} object.
#' @param x An \code{LDData} object.
#' @return An \code{LDBlocks} or \code{data.frame}.
#' @export
setGeneric("getBlockMetadata",
  function(x) standardGeneric("getBlockMetadata"))

# =============================================================================
# RegionalData accessor generics
# =============================================================================

#' @title Get Residualized Genotypes
#' @description Compute residualized genotypes on demand for a given
#'   condition index.
#' @param x A \code{RegionalData} object.
#' @param condition Integer index of the condition.
#' @return A numeric matrix of residualized genotypes.
#' @export
setGeneric("getResidualX",
  function(x, condition = 1L) standardGeneric("getResidualX"))

#' @title Get Residualized Phenotypes
#' @description Compute residualized phenotypes on demand for a given
#'   condition index.
#' @param x A \code{RegionalData} object.
#' @param condition Integer index of the condition.
#' @return A numeric matrix of residualized phenotypes.
#' @export
setGeneric("getResidualY",
  function(x, condition = 1L) standardGeneric("getResidualY"))

#' @title Get Residual X Scalar
#' @description Compute per-variant SDs of residualized genotypes.
#' @param x A \code{RegionalData} object.
#' @param condition Integer index of the condition.
#' @return A numeric vector or 1 if not scaling.
#' @export
setGeneric("getResidualXScalar",
  function(x, condition = 1L) standardGeneric("getResidualXScalar"))

#' @title Get Residual Y Scalar
#' @description Compute per-column SDs of residualized phenotypes.
#' @param x A \code{RegionalData} object.
#' @param condition Integer index of the condition.
#' @return A numeric vector or 1 if not scaling.
#' @export
setGeneric("getResidualYScalar",
  function(x, condition = 1L) standardGeneric("getResidualYScalar"))

#' @title Get Per-Variant Variance
#' @description Per-variant variance of residualized genotypes for a
#'   condition.
#' @param x A \code{RegionalData} object.
#' @param condition Integer index of the condition.
#' @return A numeric vector (length = number of variants).
#' @export
setGeneric("getXVariance",
  function(x, condition = 1L) standardGeneric("getXVariance"))

#' @title Get Phenotype List
#' @description Extract the per-condition phenotype list from a
#'   \code{RegionalData}.
#' @param x A \code{RegionalData} object.
#' @return A named list of phenotype matrices.
#' @export
setGeneric("getPhenotypes", function(x) standardGeneric("getPhenotypes"))

#' @title Get Covariate List
#' @description Extract the per-condition covariate list from a
#'   \code{RegionalData}.
#' @param x A \code{RegionalData} object.
#' @return A named list of covariate matrices.
#' @export
setGeneric("getCovariates", function(x) standardGeneric("getCovariates"))

#' @title Get Genotype Matrix
#' @description Extract the raw genotype matrix from a
#'   \code{RegionalData} or \code{MultivariateRegionalData}.
#' @param x The object.
#' @return A numeric matrix (samples x variants).
#' @export
setGeneric("getGenotypeMatrix", function(x) standardGeneric("getGenotypeMatrix"))

#' @title Get Region Chromosome
#' @description Extract the chromosome name from a region-bearing S4 object.
#' @param x The object.
#' @return A single character string, or NULL.
#' @export
setGeneric("getChrom", function(x) standardGeneric("getChrom"))

#' @title Get Region Range
#' @description Extract the start/end positions from a region-bearing S4
#'   object as a character vector \code{c(start, end)}.
#' @param x The object.
#' @return A character vector of length 2, or NULL.
#' @export
setGeneric("getGrange", function(x) standardGeneric("getGrange"))

#' @title Get Multivariate Y Matrix
#' @description Extract the multivariate phenotype matrix from a
#'   \code{MultivariateRegionalData}.
#' @param x A \code{MultivariateRegionalData} object.
#' @return A numeric matrix (samples x conditions).
#' @export
setGeneric("getYMatrix", function(x) standardGeneric("getYMatrix"))

#' @title Get Y Scaling Factors
#' @description Per-condition scaling factors used for residualized
#'   multivariate phenotypes.
#' @param x A \code{MultivariateRegionalData} object.
#' @return A numeric vector (length = number of conditions).
#' @export
setGeneric("getYScalar", function(x) standardGeneric("getYScalar"))

# =============================================================================
# FineMappingResult accessor generics
# =============================================================================

#' @title Get PIP Values
#' @description Extract posterior inclusion probabilities.
#' @param x A \code{FineMappingResult} object.
#' @return A named numeric vector of PIPs.
#' @export
setGeneric("getPIP", function(x) standardGeneric("getPIP"))

#' @title Get Trimmed Fit
#' @description Extract the trimmed SuSiE fit from a FineMappingResult.
#' @param x A \code{FineMappingResult} object.
#' @return A list (trimmed SuSiE fit).
#' @export
setGeneric("getTrimmedFit", function(x) standardGeneric("getTrimmedFit"))

#' @title Get Variant Names
#' @description Extract variant names from a FineMappingResult.
#' @param x A \code{FineMappingResult} object.
#' @return Character vector of variant names.
#' @export
setGeneric("getVariantNames", function(x) standardGeneric("getVariantNames"))

#' @title Get Top Loci
#' @description Extract top loci data.frame from a FineMappingResult.
#' @param x A \code{FineMappingResult} object.
#' @return A data.frame of top loci.
#' @export
setGeneric("getTopLoci", function(x) standardGeneric("getTopLoci"))

#' @title Get Credible Sets
#' @description Extract credible set assignments.
#' @param x A \code{FineMappingResult} object.
#' @param coverage Numeric, coverage level to extract.
#' @return A data.frame of credible set information.
#' @export
setGeneric("getCS",
  function(x, coverage = 0.95) standardGeneric("getCS"))

#' @title Get Log Bayes Factors
#' @description Extract per-variant log Bayes factors from a fine-mapping result.
#'   Returns a data.frame with variant names and one column per effect (L1, L2, ...).
#' @param x A \code{FineMappingResult} object.
#' @return A data.frame with columns \code{variant_id} and one numeric column
#'   per effect containing log Bayes factors.
#' @export
setGeneric("getLBF", function(x) standardGeneric("getLBF"))

#' @title Get Per-Effect Fine-Mapping Summary
#' @description Extract per-effect information from a fine-mapping result:
#'   prior variance, credible set log BF, purity, coverage, and member variants.
#' @param x A \code{FineMappingResult} object.
#' @return A data.frame with one row per effect.
#' @export
setGeneric("getEffects", function(x) standardGeneric("getEffects"))

# =============================================================================
# TWASWeights accessor generics
# =============================================================================

#' @title Get TWAS Weight Matrices
#' @description Extract weight matrices from a TWASWeights object.
#' @param x A \code{TWASWeights} object.
#' @param method Character, specific method name. If NULL, returns all.
#' @return A matrix or named list of matrices.
#' @export
setGeneric("getWeights",
  function(x, method = NULL) standardGeneric("getWeights"))

#' @title Get Standardized Flag
#' @description Check whether weights are on the standardized (correlation) scale.
#' @param x A \code{TWASWeights} object.
#' @return Logical.
#' @export
setGeneric("getStandardized", function(x) standardGeneric("getStandardized"))

#' @title Get CV Performance
#' @description Extract cross-validation performance metrics.
#' @param x A \code{TWASWeights} object.
#' @param method Character, specific method name. If NULL, returns all.
#' @return A list or single element.
#' @export
setGeneric("getCVPerformance",
  function(x, method = NULL) standardGeneric("getCVPerformance"))

#' @title Get Model Fits
#' @description Extract fitted model objects from a TWASWeights object.
#' @param x A \code{TWASWeights} object.
#' @param method Character, specific method name. If NULL, returns all.
#' @return A list or single element.
#' @export
setGeneric("getFits",
  function(x, method = NULL) standardGeneric("getFits"))

#' @title Get Method Names
#' @description Extract method names from a TWASWeights object.
#' @param x A \code{TWASWeights} object.
#' @return Character vector.
#' @export
setGeneric("getMethodNames", function(x) standardGeneric("getMethodNames"))

#' @title Get Molecular ID
#' @description Extract molecular/gene identifier from a TWASWeights object.
#' @param x A \code{TWASWeights} object.
#' @return Character string (length 0 or 1).
#' @export
setGeneric("getMolecularId", function(x) standardGeneric("getMolecularId"))

#' @title Get Data Type
#' @description Extract data type metadata from a TWASWeights object.
#' @param x A \code{TWASWeights} object.
#' @return A named list of data types per context, or NULL.
#' @export
setGeneric("getDataType", function(x) standardGeneric("getDataType"))

# =============================================================================
# AlleleQCResult accessor generics
# =============================================================================

#' @title Get Harmonized Variant Data
#' @description Extract the post-QC, reference-harmonized variants from an
#'   \code{AlleleQCResult}.
#' @param x An \code{AlleleQCResult} object.
#' @return A \code{data.frame} of harmonized variants.
#' @export
setGeneric("getHarmonizedData", function(x) standardGeneric("getHarmonizedData"))

#' @title Get Allele QC Summary
#' @description Extract the full per-variant merge/flip/strand diagnostics
#'   produced by allele QC.
#' @param x An \code{AlleleQCResult} object.
#' @return A \code{data.frame} with the diagnostic columns.
#' @export
setGeneric("getQCSummary", function(x) standardGeneric("getQCSummary"))

# =============================================================================
# QCResult accessor generics
# =============================================================================

#' @title Get LD Data
#' @description Extract the post-QC LDData payload from a QCResult.
#' @param x A \code{QCResult} object.
#' @return An \code{LDData} object, or NULL when QC produced no LD reference.
#' @export
setGeneric("getLDData", function(x) standardGeneric("getLDData"))

#' @title Get RSS Input
#' @description Extract the post-QC summary-statistic record (sumstats, n, var_y).
#' @param x A \code{QCResult} object.
#' @return A list with \code{sumstats}, \code{n}, \code{var_y}.
#' @export
setGeneric("getRSSInput", function(x) standardGeneric("getRSSInput"))

#' @title Get Preprocess Snapshot
#' @description Extract the pre-imputation snapshot (\code{sumstats},
#'   \code{ld_data}) captured before any LD-mismatch QC or RAISS imputation.
#' @param x A \code{QCResult} object.
#' @return A list with \code{sumstats} and \code{ld_data}.
#' @export
setGeneric("getPreprocess", function(x) standardGeneric("getPreprocess"))

#' @title Get Outlier Number
#' @description Number of LD-mismatch outliers removed during QC.
#' @param x A \code{QCResult} object.
#' @return Integer count.
#' @export
setGeneric("getOutlierNumber", function(x) standardGeneric("getOutlierNumber"))

#' @title Is Skipped
#' @description Whether QC short-circuited (e.g. no signals, too few variants).
#' @param x A \code{QCResult} object.
#' @return Single logical.
#' @export
setGeneric("isSkipped", function(x) standardGeneric("isSkipped"))

#' @title Get Skip Reason
#' @description Why QC short-circuited; empty string if not skipped.
#' @param x A \code{QCResult} object.
#' @return Character scalar.
#' @export
setGeneric("getSkipReason", function(x) standardGeneric("getSkipReason"))

# =============================================================================
# VCF/BCF writer generic
# =============================================================================

#' Write summary statistics or fine-mapping results to VCF/BCF
#'
#' Creates a VCF object from GWAS summary statistics or fine-mapping results
#' and writes it to disk. Supports bgzipped VCF (.vcf.gz/.vcf.bgz) and
#' BCF (.bcf) output formats via VariantAnnotation and Rsamtools.
#'
#' @param x Input data: a \code{GWASSumStats} object, a
#'   \code{FineMappingResult} object, or a data.frame with columns
#'   \code{chrom}, \code{pos}, \code{ref}, \code{alt}.
#' @param output_path File path for output. Extension determines format:
#'   \code{.vcf.gz} or \code{.vcf.bgz} for bgzipped VCF,
#'   \code{.bcf} for BCF, \code{.vcf} for uncompressed VCF.
#' @param sample_name Name for the VCF sample column (default: trait name or
#'   method name from the S4 object).
#' @param ... Additional arguments passed to methods.
#' @return Invisible path to the written file.
#' @export
setGeneric("writeSumstatsVcf",
  function(x, output_path, sample_name = NULL, ...) standardGeneric("writeSumstatsVcf"))
