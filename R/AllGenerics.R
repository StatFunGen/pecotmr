#' @title S4 Generic Function Definitions
#' @description All S4 generic function definitions for pecotmr.
#' @name pecotmr-generics
#' @keywords internal
#' @importFrom methods setGeneric
NULL

# =============================================================================
# High-level estimation generic
# =============================================================================

#' @title Estimate SNP Heritability
#' @description Estimate SNP heritability from GWAS summary statistics using
#'   one of three methods: LDER, g-LDSC, or HDL/sHDL.
#' @param sumstats A \code{GwasSumStats} object.
#' @param ldRef An \code{LdStatistic} object (method-appropriate subclass).
#' @param method Character, one of "lder", "gldsc", "hdl".
#' @param annotations An \code{AnnotationMatrix} object, or NULL for
#'   unstratified estimation.
#' @param local Logical, whether to compute per-block local estimates.
#' @param ... Additional method-specific arguments.
#' @return An \code{H2Estimate} object.
#' @export
setGeneric("estimateH2",
  function(sumstats, ldRef, method = "lder", annotations = NULL,
           local = FALSE, ...)
    standardGeneric("estimateH2")
)

# =============================================================================
# LD score computation
# =============================================================================

#' @title Compute LD Scores
#' @description Compute LD scores from an LD reference, optionally
#'   stratified by annotations.
#' @param ldRef An \code{LdStatistic} object.
#' @param annotations An \code{AnnotationMatrix} object, or NULL.
#' @param ... Additional arguments.
#' @return A numeric matrix of LD scores (SNPs x annotations+1).
#' @export
setGeneric("computeLdScores",
  function(ldRef, annotations = NULL, ...)
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

#' @title Read Annotations
#' @description Read genomic annotations from files (BED, BigWig,
#'   S-LDSC .annot format, or GRanges objects) and create an
#'   AnnotationMatrix.
#' @param paths Named character vector of file paths, or a named list
#'   of GRanges objects. Names become annotation names.
#' @param snpRanges A \code{GRanges} object defining SNP positions.
#' @param annotationMeta A \code{data.frame} with annotation metadata
#'   (name, tier, type). If NULL, auto-detected from file format.
#' @param genome Character, genome build.
#' @param ... Additional arguments.
#' @return An \code{AnnotationMatrix} object.
#' @export
setGeneric("readAnnotations",
  function(paths, snpRanges, annotationMeta = NULL,
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
# GwasSumStats accessor generics
# =============================================================================

#' @title Get Z-scores
#' @description Extract z-score vector from a \code{GwasSumStats} or
#'   \code{QtlSumStats} entry, selected by its identity tuple.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments (e.g., \code{study} for
#'   \code{GwasSumStats}; \code{study}, \code{context}, \code{trait} for
#'   \code{QtlSumStats}).
#' @return Numeric vector of z-scores.
#' @export
setGeneric("getZ", function(x, ...) standardGeneric("getZ"))

#' @title Get Sample Sizes
#' @description Extract sample size vector from a \code{GwasSumStats} or
#'   \code{QtlSumStats} entry, selected by its identity tuple.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @return Numeric vector of sample sizes.
#' @export
setGeneric("getN", function(x, ...) standardGeneric("getN"))

#' @title Get Minor Allele Frequencies
#' @description Extract MAF vector from a GwasSumStats object.
#' @param x A \code{GwasSumStats} or \code{QtlDataset} object.
#' @param ... Class-specific selection arguments (e.g., \code{region},
#'   \code{cisWindow} for \code{QtlDataset}).
#' @return Numeric vector of MAFs, or NULL if not available.
#' @export
setGeneric("getMaf", function(x, ...) standardGeneric("getMaf"))

#' @title Get Number of SNPs
#' @description Number of SNPs in a \code{GwasSumStats} or
#'   \code{QtlSumStats} entry, selected by its identity tuple.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @return Integer.
#' @export
setGeneric("nSnps", function(x, ...) standardGeneric("nSnps"))

#' @title Subset by Chromosome
#' @description Extract a chromosome-specific subset of a GwasSumStats object.
#' @param x A \code{GwasSumStats} object.
#' @param chr Character, chromosome name (e.g., "1", "chr1").
#' @return A \code{GwasSumStats} object.
#' @export
setGeneric("subsetChr", function(x, chr) standardGeneric("subsetChr"))

#' @title Get Phenotype Variance
#' @description Extract phenotype variance from a \code{GwasSumStats} or
#'   \code{QtlSumStats} entry, selected by its identity tuple. Returns
#'   \code{NULL} when the entry has no \code{varY} recorded.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selection arguments.
#' @return Numeric phenotype variance, or NULL.
#' @export
setGeneric("getVarY", function(x, ...) standardGeneric("getVarY"))

#' @title Get a Single Summary-Statistic Entry or Embedded Collection
#' @description Behavior depends on the class of \code{x}:
#'   \describe{
#'     \item{For \code{GwasSumStats} / \code{QtlSumStats}}{Returns the
#'       per-variant \code{GRanges} of summary statistics for one entry,
#'       selected by its identity tuple (\code{study} for GWAS;
#'       \code{study}, \code{context}, \code{trait} for QTL).}
#'     \item{For \code{MultiStudyQtlDataset}}{Returns the embedded
#'       \code{QtlSumStats} collection (the summary-statistic-only
#'       studies), or \code{NULL} when absent. No selection arguments
#'       are accepted in this case.}
#'   }
#' @param x A \code{GwasSumStats}, \code{QtlSumStats}, or
#'   \code{MultiStudyQtlDataset} object.
#' @param ... Class-specific selection arguments (see above).
#' @return A \code{GRanges}, a \code{QtlSumStats}, or \code{NULL}.
#' @export
setGeneric("getSumStats", function(x, ...) standardGeneric("getSumStats"))

#' @title Get Standardized Sumstat Data Frame for One Tuple
#' @description Return a per-tuple summary-statistics \code{data.frame}
#'   in the standardized layout \code{variant_id, chrom, pos, A1, A2,
#'   z, beta, se, N, maf} (optional columns omitted when absent on the
#'   entry). Combines tuple-keyed row selection (\code{getSumStats})
#'   with mcols unpacking; replaces the pre-S4 idiom of pulling
#'   \code{S4Vectors::mcols(entry)$<col>} directly inside pipelines.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Class-specific selectors (\code{study} for
#'   \code{GwasSumStats}; \code{study}, \code{context}, \code{trait}
#'   for \code{QtlSumStats}) plus pass-throughs \code{require},
#'   \code{derive}, \code{keepChrPrefix} forwarded to the underlying
#'   unpacker.
#' @return A \code{data.frame}.
#' @export
setGeneric("getSumstatDf", function(x, ...) standardGeneric("getSumstatDf"))

#' @title Get the Embedded QtlDataset List
#' @description Return the named list of \code{QtlDataset} objects
#'   carried by a \code{MultiStudyQtlDataset}.
#' @param x A \code{MultiStudyQtlDataset} object.
#' @return A named list of \code{QtlDataset} objects.
#' @export
setGeneric("getQtlDatasets",
  function(x) standardGeneric("getQtlDatasets"))

#' @title Get the Genome Build
#' @description Return the genome build that the collection's LD sketch
#'   and every entry are aligned to. Because all entries in a
#'   \code{GwasSumStats} or \code{QtlSumStats} share the LD sketch, the
#'   genome build is a single value at the collection level.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Unused (present for method-signature compatibility).
#' @return Character (length 1).
#' @export
setGeneric("getGenome", function(x, ...) standardGeneric("getGenome"))

#' @title Get QC Audit Record
#' @description Return the audit record of QC steps applied to this
#'   collection. An empty \code{list()} (default on construction) means
#'   \code{\link{summaryStatsQc}} has not yet been run. Pipelines that
#'   require harmonized sumstats (\code{fineMappingPipeline},
#'   \code{twasWeightsPipeline}, and downstream consumers) reject inputs
#'   where \code{length(getQcInfo(x)) == 0L}.
#' @param x A \code{GwasSumStats} or \code{QtlSumStats} object.
#' @param ... Unused.
#' @return A \code{list} (possibly empty).
#' @export
setGeneric("getQcInfo", function(x, ...) standardGeneric("getQcInfo"))

#' @title Get LD Sketch
#' @description Return the \code{GenotypeHandle} carrying the LD
#'   reference for this collection. Defined on classes that embed an
#'   \code{ldSketch} slot: \code{GwasSumStats}, \code{QtlSumStats},
#'   \code{FineMappingResult}, \code{TwasWeights}. Returns \code{NULL}
#'   when the slot is unset (e.g. a \code{TwasWeights} fit from
#'   individual-level data via \code{QtlDataset}).
#' @param x An S4 object that carries an \code{ldSketch} slot.
#' @param ... Unused.
#' @return A \code{GenotypeHandle} or \code{NULL}.
#' @export
setGeneric("getLdSketch", function(x, ...) standardGeneric("getLdSketch"))

# =============================================================================
# LdData accessor generics
# =============================================================================

#' @title Get LD Correlation Matrix
#' @description Extract the LD correlation matrix from an \code{LdData}
#'   object. If only a genotype handle is available, recomputes R from
#'   genotypes on the fly.
#' @param x An \code{LdData} object.
#' @return A correlation matrix, or a list of per-block matrices.
#' @export
setGeneric("getCorrelation", function(x) standardGeneric("getCorrelation"))

#' @title Get Genotype Matrix
#' @description Extract a genotype matrix from an object that carries
#'   genotype data. For an \code{LdData}, returns the underlying genotype
#'   matrix via its handle (or \code{NULL} if no handle is available).
#'   For a \code{QtlDataset}, returns the genotype matrix for a selected
#'   set of traits or region (see method documentation for the
#'   per-class selection arguments).
#' @param x The object to extract from.
#' @param ... Class-specific selection arguments (e.g., \code{traitId},
#'   \code{region}, \code{cisWindow} for \code{QtlDataset}).
#' @return A numeric matrix, a list of matrices, or \code{NULL}.
#' @export
setGeneric("getGenotypes", function(x, ...) standardGeneric("getGenotypes"))

#' @title Check Genotype Availability
#' @description Check whether an \code{LdData} object has a genotype
#'   handle for extracting raw genotypes.
#' @param x An \code{LdData} object.
#' @return Logical.
#' @export
setGeneric("hasGenotypes", function(x) standardGeneric("hasGenotypes"))

#' @title Get Variant IDs
#' @description Extract variant ID vector from an object that carries one
#'   (e.g., \code{LdData}, \code{FineMappingEntry}, \code{TwasWeightsEntry})
#'   or from one entry of a collection class selected by its identity
#'   tuple.
#' @param x The object.
#' @param ... Class-specific selection arguments.
#' @return Character vector of variant IDs.
#' @export
setGeneric("getVariantIds", function(x, ...) standardGeneric("getVariantIds"))

#' @title Get Phenotype List
#' @description Extract phenotype data from an object that carries it.
#'   For a \code{QtlDataset}, the user can optionally select specific
#'   contexts, traits, or a region (see method documentation for the
#'   per-class selection arguments).
#' @param x The object to extract from.
#' @param ... Class-specific selection arguments (e.g., \code{contexts},
#'   \code{traitId}, \code{region}).
#' @return A named list of phenotype matrices or
#'   \code{SummarizedExperiment} objects.
#' @export
setGeneric("getPhenotypes", function(x, ...) standardGeneric("getPhenotypes"))
# =============================================================================
# FineMappingResult accessor generics
# =============================================================================

#' @title Get a Single Fine-Mapping Entry
#' @description Return the \code{FineMappingEntry} for one
#'   \code{(study, context, trait, method)} row of a
#'   \code{FineMappingResult} collection.
#' @param x A \code{FineMappingResult} object.
#' @param study,context,trait,method Single character identifiers. All
#'   required when the collection has more than one row; optional when
#'   the collection has a single row.
#' @return A \code{FineMappingEntry} object.
#' @export
setGeneric("getFineMappingResult",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL)
    standardGeneric("getFineMappingResult"))

#' @title Renormalize Fine-Mapping PIPs to a Variant Subset
#' @description Re-derive a \code{FineMappingEntry}'s PIPs (and the
#'   \code{topLoci} table) after restricting to a kept variant subset.
#'   For each effect the \code{lbf_variable} row is subset to the kept
#'   variants, renormalized via \code{lbfToAlpha()}, and the per-variant
#'   PIPs are recomputed as \code{1 - prod_l(1 - alpha[l, p])}.
#'
#'   The two scenarios this supports:
#'   \itemize{
#'     \item The user declined to impute missing variants in a GWAS
#'           \code{SumStats}, so a downstream fine-mapping result needs
#'           PIPs restricted to the GWAS-covered intersection.
#'     \item Colocalization between a GWAS \code{FineMappingResult} and a
#'           QTL \code{FineMappingResult} computed on different variant
#'           sets — the GWAS PIPs (or QTL PIPs) get renormalized to the
#'           common variant set.
#'   }
#'
#' @param x A \code{FineMappingEntry} or \code{FineMappingResultBase}.
#' @param keepVariants Character vector of variant IDs to keep. Intersected
#'   with the entry's own \code{variantIds}; an empty intersection raises
#'   an error.
#' @param ... Future expansion.
#' @return The same flavour of object with PIPs renormalized on the kept
#'   subset.
#' @export
setGeneric("adjustPips",
  function(x, keepVariants, ...) standardGeneric("adjustPips"))

#' @title Get PIP Values
#' @description Extract posterior inclusion probabilities from a single
#'   \code{FineMappingEntry} or from one entry of a
#'   \code{FineMappingResult} (selected by its identity tuple).
#' @param x A \code{FineMappingEntry} or \code{FineMappingResult}.
#' @param ... Class-specific selection arguments.
#' @return A named numeric vector of PIPs.
#' @export
setGeneric("getPip", function(x, ...) standardGeneric("getPip"))

#' @title Get SuSiE Fit
#' @description Extract the SuSiE fit object from a fine-mapping entry
#'   or result. The fit may be the trimmed view (when the pipeline ran
#'   with the default \code{trim = TRUE}) or the full untrimmed
#'   \code{susie()} return (when \code{trim = FALSE}).
#' @param x A \code{FineMappingEntry} or \code{FineMappingResult}.
#' @param ... Class-specific selection arguments.
#' @return A list (the SuSiE fit object).
#' @export
setGeneric("getSusieFit", function(x, ...) standardGeneric("getSusieFit"))

#' @title Get Marginal Effects
#' @description Extract per-variant marginal univariate effects from a
#'   fine-mapping entry or result. Returns a \code{data.frame} with
#'   identity columns (\code{variant_id, chrom, pos, A1, A2}), context
#'   (\code{N, MAF}), and the marginal effect columns
#'   (\code{beta, se, z, p}). Populated uniformly across the
#'   individual-level and RSS paths.
#' @param x A \code{FineMappingEntry} or \code{FineMappingResult}.
#' @param maxPval Optional numeric (length 1). When non-\code{NULL},
#'   filter rows where \code{p > maxPval}. Default \code{NULL}
#'   (no filter).
#' @param ... Class-specific selection arguments.
#' @return A \code{data.frame}.
#' @export
setGeneric("getMarginalEffects",
  function(x, maxPval = NULL, ...) standardGeneric("getMarginalEffects"))

#' @title Get Top Loci (posterior view)
#' @description Extract the per-variant posterior fine-mapping payload
#'   as either a \code{data.frame} (default) or a \code{GRanges}.
#'   Returns identity columns (\code{variant_id, chrom, pos, A1, A2}),
#'   context (\code{N, MAF}), the posterior effect columns
#'   (\code{beta = posterior_mean, se = posterior_sd}), \code{pip},
#'   and credible-set membership columns (\code{cs_95}, etc.).
#'   Rows are filtered by PIP by default — set \code{signalCutoff = 0}
#'   to return every variant.
#' @param x A \code{FineMappingEntry} or \code{FineMappingResult}.
#' @param type One of \code{"data.frame"} (default) or \code{"GRanges"}.
#' @param signalCutoff Numeric (length 1). Drop rows where
#'   \code{pip <= signalCutoff}. Default \code{0.025}. Use
#'   \code{signalCutoff = 0} to keep every variant.
#' @param ... Class-specific selection arguments.
#' @return A \code{data.frame} or a \code{GRanges}.
#' @export
setGeneric("getTopLoci",
  function(x, type = c("data.frame", "GRanges"),
           signalCutoff = 0.025, ...)
    standardGeneric("getTopLoci"))

#' @title Get Credible Sets
#' @description Extract credible set assignments at the requested coverage.
#' @param x A \code{FineMappingEntry} or \code{FineMappingResult}.
#' @param ... Class-specific selection arguments plus \code{coverage}.
#' @return A data.frame of credible set information.
#' @export
setGeneric("getCs", function(x, ...) standardGeneric("getCs"))

# =============================================================================
# TwasWeights accessor generics
# =============================================================================

#' @title Get TWAS Weights
#' @description Extract weights from a \code{TwasWeightsEntry} or from
#'   one entry of a \code{TwasWeights} collection.
#' @param x A \code{TwasWeightsEntry} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @return A numeric vector or matrix of weights.
#' @export
setGeneric("getWeights", function(x, ...) standardGeneric("getWeights"))

#' @title Get Standardized Flag
#' @description Check whether weights are on the standardized scale.
#' @param x A \code{TwasWeightsEntry} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @return Logical.
#' @export
setGeneric("getStandardized",
  function(x, ...) standardGeneric("getStandardized"))

#' @title Get CV Performance
#' @description Extract cross-validation performance metrics.
#' @param x A \code{TwasWeightsEntry} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @return Method-specific (typically a list).
#' @export
setGeneric("getCvPerformance",
  function(x, ...) standardGeneric("getCvPerformance"))

#' @title Get Model Fits
#' @description Extract fitted model objects.
#' @param x A \code{TwasWeightsEntry} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @return Method-specific (typically a list).
#' @export
setGeneric("getFits", function(x, ...) standardGeneric("getFits"))

#' @title Get Method Names
#' @description Extract method names from a collection class.
#' @param x A \code{FineMappingResult} or \code{TwasWeights} object.
#' @return Character vector.
#' @export
setGeneric("getMethodNames", function(x) standardGeneric("getMethodNames"))

#' @title Get Data Type
#' @description Extract the data-type tag.
#' @param x A \code{TwasWeightsEntry} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @return A character vector or NULL.
#' @export
setGeneric("getDataType", function(x, ...) standardGeneric("getDataType"))

# =============================================================================
# AlleleQcResult accessor generics
# QcResult accessor generics
# =============================================================================
# VCF/BCF writer generic
# =============================================================================

#' Write summary statistics or fine-mapping results to VCF/BCF
#'
#' Creates a VCF object from GWAS summary statistics or fine-mapping results
#' and writes it to disk. Supports bgzipped VCF (.vcf.gz/.vcf.bgz) and
#' BCF (.bcf) output formats via VariantAnnotation and Rsamtools.
#'
#' @param x Input data: a \code{GwasSumStats} object, a
#'   \code{FineMappingResult} object, or a data.frame with columns
#'   \code{chrom}, \code{pos}, \code{ref}, \code{alt}.
#' @param outputPath File path for output. Extension determines format:
#'   \code{.vcf.gz} or \code{.vcf.bgz} for bgzipped VCF,
#'   \code{.bcf} for BCF, \code{.vcf} for uncompressed VCF.
#' @param sampleName Name for the VCF sample column (default: trait name or
#'   method name from the S4 object).
#' @param ... Additional arguments passed to methods.
#' @return Invisible path to the written file.
#' @export
setGeneric("writeSumstatsVcf",
  function(x, outputPath, sampleName = NULL, ...) standardGeneric("writeSumstatsVcf"))

# =============================================================================
# QtlDataset accessor generics
# =============================================================================

#' @title Get Study Identifier
#' @description Return the study identifier carried by a \code{QtlDataset}.
#' @param x A \code{QtlDataset} object.
#' @return Character (length 1).
#' @export
setGeneric("getStudy", function(x) standardGeneric("getStudy"))

#' @title Get Context Names
#' @description Return the names of all contexts carried by an object
#'   (e.g., the keys of the \code{phenotypes} list on a \code{QtlDataset},
#'   or the unique \code{context} values of a \code{QtlSumStats}).
#' @param x The object.
#' @return Character vector of context names.
#' @export
setGeneric("getContexts", function(x) standardGeneric("getContexts"))

#' @title Get Unique Trait Names
#' @description Return the unique trait identifiers carried by a
#'   collection class (e.g., \code{QtlSumStats}).
#' @param x The object.
#' @return Character vector of unique trait names.
#' @export
setGeneric("getTraits", function(x) standardGeneric("getTraits"))

#' @title Get Residualized Genotypes
#' @description Residualize the genotype matrix against the per-context
#'   phenotype covariates and the genotype covariates, optionally
#'   subsetting variants to those falling within a trait's cis-window or
#'   an explicit region.
#' @param x A \code{QtlDataset} object.
#' @param ... Selection arguments: \code{traitId}, \code{region},
#'   \code{cisWindow}, \code{phenotypeCovariatesToRemove},
#'   \code{genotypeCovariatesToRemove}.
#' @return A numeric matrix (samples x variants).
#' @export
setGeneric("getResidualizedGenotypes",
  function(x, ...) standardGeneric("getResidualizedGenotypes"))

#' @title Get Residualized Phenotypes
#' @description Residualize the per-context phenotype matrices against
#'   the per-context phenotype covariates and the genotype covariates,
#'   for one or more requested contexts.
#' @param x A \code{QtlDataset} object.
#' @param ... Selection arguments: \code{contexts} (required),
#'   \code{traitId}, \code{region},
#'   \code{phenotypeCovariatesToRemove},
#'   \code{genotypeCovariatesToRemove}.
#' @return A named list of numeric matrices keyed by context.
#' @export
setGeneric("getResidualizedPhenotypes",
  function(x, ...) standardGeneric("getResidualizedPhenotypes"))

#' @title Get Per-Context Phenotype Covariates
#' @description Return per-context phenotype covariate matrices, taken
#'   from the \code{colData} of each context's \code{SummarizedExperiment}.
#' @param x A \code{QtlDataset} object.
#' @param contexts Character vector of context names (subset of
#'   \code{names(getPhenotypes(x))}).
#' @return A named list of matrices keyed by context.
#' @export
setGeneric("getPhenotypeCovariates",
  function(x, contexts) standardGeneric("getPhenotypeCovariates"))

#' @title Get Genotype Covariates
#' @description Return the single genotype-derived covariate matrix
#'   carried by a \code{QtlDataset} (e.g., ancestry PCs).
#' @param x A \code{QtlDataset} object.
#' @return Numeric matrix (samples x covariates).
#' @export
setGeneric("getGenotypeCovariates",
  function(x) standardGeneric("getGenotypeCovariates"))

#' @title Get scaleResiduals Flag
#' @description Whether residualization accessors scale residuals to unit
#'   variance.
#' @param x A \code{QtlDataset} object.
#' @return Logical (length 1).
#' @export
setGeneric("getScaleResiduals",
  function(x) standardGeneric("getScaleResiduals"))

# =============================================================================
# GenotypeHandle / LD-statistic / Annotation / LdData / H2Estimate accessors
# =============================================================================

#' @title Get SNP Info
#' @description Return the cached SNP metadata data.frame
#'   (columns: SNP, CHR, BP, A1, A2, optionally MAF).
#' @param x A \code{GenotypeHandle} or \code{LdStatistic}.
#' @return A data.frame.
#' @export
setGeneric("getSnpInfo", function(x) standardGeneric("getSnpInfo"))

#' @title Get Genotype Storage Format
#' @description Return the detected genotype storage format.
#' @param x A \code{GenotypeHandle}.
#' @return Character (length 1): one of "gds", "vcf", "plink1", "plink2".
#' @export
setGeneric("getFormat", function(x) standardGeneric("getFormat"))

#' @title Get File Path
#' @description Return the underlying genotype file path or stem.
#' @param x A \code{GenotypeHandle}.
#' @return Character (length 1).
#' @export
setGeneric("getPath", function(x) standardGeneric("getPath"))

#' @title Get Sample Identifiers
#' @description Return the sample-id vector.
#' @param x A \code{GenotypeHandle}.
#' @return Character vector.
#' @export
setGeneric("getSampleIds", function(x) standardGeneric("getSampleIds"))

#' @title Get plink2 pgen Pointer
#' @description Return the cached external pointer to the plink2 pgen
#'   handle (NULL when the handle is not pgen-backed).
#' @param x A \code{GenotypeHandle}.
#' @return An external pointer or NULL.
#' @export
setGeneric("getPgenPtr", function(x) standardGeneric("getPgenPtr"))

#' @title Get Sample Count
#' @description Return the number of samples carried by a
#'   \code{GenotypeHandle}.
#' @param x A \code{GenotypeHandle}.
#' @return Integer (length 1).
#' @export
setGeneric("getNSamples", function(x) standardGeneric("getNSamples"))

#' @title Get Per-Block Eigendecompositions
#' @description Return the per-block eigendecomposition list carried by
#'   an \code{LdEigen} object.
#' @param x An \code{LdEigen}.
#' @return List of per-block eigen decompositions.
#' @export
setGeneric("getEigenList", function(x) standardGeneric("getEigenList"))

#' @title Get LD Reference Panel Size
#' @description Return the reference-panel sample size used to compute
#'   an \code{LdStatistic} or carried by an \code{LdData}.
#' @param x An \code{LdStatistic} or \code{LdData}.
#' @return Integer (length 1).
#' @export
setGeneric("getNRef", function(x) standardGeneric("getNRef"))

#' @title Get In-Sample Flag
#' @description Whether the LD reference panel is from the same cohort
#'   as the GWAS (affects bias correction).
#' @param x An \code{LdStatistic}.
#' @return Logical (length 1).
#' @export
setGeneric("getInSample", function(x) standardGeneric("getInSample"))

#' @title Get LD Scores
#' @description Return the per-SNP LD score matrix carried by an
#'   \code{LdScore} object.
#' @param x An \code{LdScore}.
#' @return Numeric matrix (SNPs x annotations+1).
#' @export
setGeneric("getLdScores", function(x) standardGeneric("getLdScores"))

#' @title Get LD-Score Regression Weights
#' @description Return the per-SNP regression weights vector carried by
#'   an \code{LdScore} object.
#' @param x An \code{LdScore}.
#' @return Numeric vector.
#' @export
setGeneric("getLdScoreWeights",
  function(x) standardGeneric("getLdScoreWeights"))

#' @title Get Per-Block LD Matrix List
#' @description Return the list of per-block LD (R^2) matrices used for
#'   the FGLS residual covariance in g-LDSC.
#' @param x An \code{LdScore}.
#' @return List of matrices (empty list for S-LDSC).
#' @export
setGeneric("getLdMatrixList",
  function(x) standardGeneric("getLdMatrixList"))

#' @title Get LD Block Container
#' @description Return the \code{LdBlocks} object carried by an
#'   \code{LdStatistic}.
#' @param x An \code{LdStatistic}.
#' @return An \code{LdBlocks} object.
#' @export
setGeneric("getLdBlocks", function(x) standardGeneric("getLdBlocks"))

#' @title Get Annotation Matrix
#' @description Return the (SNPs x annotations) annotation matrix.
#' @param x An \code{AnnotationMatrix}.
#' @return Numeric matrix or dgCMatrix.
#' @export
setGeneric("getAnnotations",
  function(x) standardGeneric("getAnnotations"))

#' @title Get Annotation Metadata
#' @description Return the per-annotation metadata data.frame (columns
#'   \code{name}, \code{tier}, \code{type}).
#' @param x An \code{AnnotationMatrix}.
#' @return A data.frame.
#' @export
setGeneric("getAnnotationMeta",
  function(x) standardGeneric("getAnnotationMeta"))

#' @title Get SNP Ranges
#' @description Return the per-SNP \code{GRanges} carried by an
#'   \code{AnnotationMatrix}.
#' @param x An \code{AnnotationMatrix}.
#' @return A \code{GRanges} object.
#' @export
setGeneric("getSnpRanges", function(x) standardGeneric("getSnpRanges"))

#' @title Get LD Block Ranges
#' @description Return the per-block \code{GRanges} carried by an
#'   \code{LdBlocks} object.
#' @param x An \code{LdBlocks}.
#' @return A \code{GRanges} object.
#' @export
setGeneric("getBlocks", function(x) standardGeneric("getBlocks"))

#' @title Get GenotypeHandle from LdData
#' @description Return the \code{GenotypeHandle} (or list of handles for
#'   mixture panels) carried by an \code{LdData}.
#' @param x An \code{LdData}.
#' @return A \code{GenotypeHandle}, a list of them, or NULL.
#' @export
setGeneric("getGenotypeHandle",
  function(x) standardGeneric("getGenotypeHandle"))

#' @title Get Mixture Weights
#' @description Return the per-panel mixing proportions carried by an
#'   \code{LdData} when its \code{genotypeHandle} slot is a list of
#'   panels. NULL for single-panel objects.
#' @param x An \code{LdData}.
#' @return Numeric vector or NULL.
#' @export
setGeneric("getMixtureWeights",
  function(x) standardGeneric("getMixtureWeights"))

#' @title Get SNP Indices
#' @description Return the integer indices into the handle's snpInfo
#'   carried by an \code{LdData}.
#' @param x An \code{LdData}.
#' @return Integer vector or NULL.
#' @export
setGeneric("getSnpIdx", function(x) standardGeneric("getSnpIdx"))

#' @title Get Variant GRanges
#' @description Return the variant metadata \code{GRanges} of an
#'   \code{LdData}.
#' @param x An \code{LdData}.
#' @return A \code{GRanges}.
#' @export
setGeneric("getVariantInfo", function(x) standardGeneric("getVariantInfo"))

#' @title Get Block Metadata
#' @description Return the block metadata (\code{LdBlocks} or
#'   \code{data.frame}) carried by an \code{LdData}.
#' @param x An \code{LdData}.
#' @return An \code{LdBlocks} or \code{data.frame}.
#' @export
setGeneric("getBlockMetadata",
  function(x) standardGeneric("getBlockMetadata"))

#' @title Get Reference Panel (data.frame)
#' @description Flatten the variant \code{GRanges} of an \code{LdData}
#'   into a reference-panel data.frame.
#' @param x An \code{LdData}.
#' @return A data.frame.
#' @export
setGeneric("getRefPanel", function(x) standardGeneric("getRefPanel"))

#' @title Get Per-Block tau Matrix
#' @description Return the per-block jackknife tau matrix carried by an
#'   \code{H2Estimate}.
#' @param x An \code{H2Estimate}.
#' @return A numeric matrix or NULL.
#' @export
setGeneric("getTauBlocks", function(x) standardGeneric("getTauBlocks"))

#' @title Get Global SNP Heritability
#' @description Return the global SNP heritability estimate carried by
#'   an \code{H2Estimate}.
#' @param x An \code{H2Estimate}.
#' @return Numeric (length 1).
#' @export
setGeneric("getH2", function(x) standardGeneric("getH2"))
