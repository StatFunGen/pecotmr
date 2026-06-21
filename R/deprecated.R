# =============================================================================
# Centralized deprecation stubs.
#
# Every function in this file is removed from pecotmr's active surface but
# kept as an exported stub so callers see a helpful redirect via
# `.Deprecated()` rather than a silent "object not found" error.
#
# Stubs come in two flavors:
#   * Non-functional stubs (`invisible(NULL)` body): the legacy code path
#     no longer exists at all; the message tells the user what to use
#     instead.
#   * Functional stubs (delegate to a still-active internal helper): the
#     legacy behavior is preserved in the short term while the public name
#     is steered toward its successor.
#
# Group ordering: by deprecation target (which active API replaces them).
# =============================================================================


# -----------------------------------------------------------------------------
# Allele harmonization (matchRefPanel / alleleQc → summaryStatsQc)
# -----------------------------------------------------------------------------

#' (Deprecated) Match summary statistics to a reference panel
#'
#' \strong{Deprecated.} Allele harmonization now runs inside
#' \code{\link{summaryStatsQc}} as part of the SumStats QC pass.
#'
#' @param targetData Data frame of target variants.
#' @param refVariants Reference variant identifiers.
#' @param ... Forwarded to the internal helper.
#' @return Result of the internal helper.
#' @export
matchRefPanel <- function(targetData, refVariants, ...) {
  .Deprecated(new = "summaryStatsQc", package = "pecotmr",
    msg = paste(
      "matchRefPanel() is deprecated; allele harmonization now runs",
      "inside summaryStatsQc() as part of the SumStats QC pass."))
  .matchRefPanel(targetData = targetData, refVariants = refVariants, ...)
}

#' @rdname matchRefPanel
#' @export
alleleQc <- function(targetData, refVariants, ...) {
  .Deprecated(new = "summaryStatsQc", package = "pecotmr",
    msg = paste(
      "alleleQc() is deprecated; allele harmonization now runs inside",
      "summaryStatsQc() as part of the SumStats QC pass."))
  .matchRefPanel(targetData = targetData, refVariants = refVariants, ...)
}


# -----------------------------------------------------------------------------
# coloc / enloc / qtl-enrichment wrappers → colocPipeline / qtlEnrichmentPipeline
# -----------------------------------------------------------------------------

#' (Deprecated) xQTL GWAS Enrichment Analysis
#'
#' \strong{Deprecated.} Use \code{\link{qtlEnrichmentPipeline}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
xqtlEnrichmentWrapper <- function(...) {
  .Deprecated(new = "qtlEnrichmentPipeline", package = "pecotmr",
    msg = paste(
      "xqtlEnrichmentWrapper() has been removed. Use",
      "qtlEnrichmentPipeline() with FineMappingResult collections for",
      "the GWAS and the QTLs."))
  invisible(NULL)
}

#' (Deprecated) Low-level QTL Enrichment Estimator
#'
#' \strong{Deprecated.} Renamed to \code{\link{qtlEnrichment}} for
#' consistency with \code{\link{qtlEnrichmentPipeline}}. This wrapper
#' forwards every argument to \code{qtlEnrichment()} after emitting a
#' deprecation message.
#'
#' @param ... Forwarded to \code{\link{qtlEnrichment}}.
#' @return Result of \code{\link{qtlEnrichment}}.
#' @export
computeQtlEnrichment <- function(...) {
  .Deprecated(new = "qtlEnrichment", package = "pecotmr",
    msg = paste(
      "computeQtlEnrichment() has been renamed to qtlEnrichment().",
      "Update callers to use qtlEnrichment(); the call signature is",
      "unchanged."))
  qtlEnrichment(...)
}

#' (Deprecated) Colocalization Wrapper
#'
#' \strong{Deprecated.} Use \code{\link{colocPipeline}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
colocWrapper <- function(...) {
  .Deprecated(new = "colocPipeline", package = "pecotmr",
    msg = paste(
      "colocWrapper() has been removed. Use colocPipeline() with a",
      "QTL FineMappingResult plus a GwasSumStats or GWAS",
      "FineMappingResult."))
  invisible(NULL)
}

#' (Deprecated) Colocalization Post-Processor
#'
#' \strong{Deprecated.} Post-processing now happens inside
#' \code{\link{colocPipeline}}, which returns the cleaned, ranked
#' colocalization output directly.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
colocPostProcessor <- function(...) {
  .Deprecated(new = "colocPipeline", package = "pecotmr",
    msg = paste(
      "colocPostProcessor() has been removed. Post-processing is now",
      "internal to colocPipeline()."))
  invisible(NULL)
}

#' (Deprecated) Enrichment-Informed Colocalization Pipeline
#'
#' \strong{Deprecated.} \code{enlocPipeline} is now a thin alias for
#' \code{\link{colocPipeline}} with the \code{enrichment} and
#' \code{p12Max} arguments supplied. The enrichment-adjusted prior
#' \code{p12 * (1 + r)} (capped at \code{p12Max}) and the per-pair
#' \code{enrichment} / \code{p12Used} output columns are now produced by
#' \code{colocPipeline} directly when \code{enrichment} is non-NULL.
#'
#' @param qtlFineMappingResult,gwasInput,enrichment Forwarded to
#'   \code{\link{colocPipeline}}.
#' @param ... Forwarded to \code{\link{colocPipeline}}.
#' @return Data frame as returned by \code{\link{colocPipeline}} with
#'   \code{enrichment} and \code{p12Used} columns populated.
#' @export
enlocPipeline <- function(qtlFineMappingResult, gwasInput,
                          enrichment, ...) {
  .Deprecated(new = "colocPipeline", package = "pecotmr",
    msg = paste(
      "enlocPipeline() is deprecated. Call colocPipeline(",
      "..., enrichment = <data.frame>, p12Max = <max>) instead;",
      "the enrichment-informed prior adjustment is now integrated",
      "into colocPipeline()."))
  colocPipeline(qtlFineMappingResult = qtlFineMappingResult,
                gwasInput            = gwasInput,
                enrichment           = enrichment,
                ...)
}


# -----------------------------------------------------------------------------
# ctwas file-path loaders → readBim / ldLoader
# -----------------------------------------------------------------------------

#' (Deprecated) Load a PLINK bim file for cTWAS
#'
#' \strong{Deprecated.} Use \code{\link{readBim}}.
#'
#' @param bimFilePath Path to a \code{.bim} or \code{.bed} file.
#' @return A data.frame with legacy column names (chrom, id, GD, pos, A1, A2).
#' @export
ctwasBimfileLoader <- function(bimFilePath) {
  .Deprecated("readBim", package = "pecotmr",
              msg = "ctwasBimfileLoader() is deprecated. Use readBim() instead.")
  bedPath <- sub("\\.bim$", ".bed", bimFilePath)
  bim <- readBim(bedPath)
  data.frame(
    chrom = bim$chrom,
    id    = normalizeVariantId(bim$id),
    GD    = bim$gpos,
    pos   = bim$pos,
    A1    = bim$a1,
    A2    = bim$a0,
    stringsAsFactors = FALSE)
}

#' (Deprecated) Load cTWAS LD meta-data
#'
#' \strong{Deprecated.} Use \code{\link{ldLoader}} with its \code{ldInfo}
#' argument instead.
#'
#' @param ldMetaDataFile Path to the LD meta-data TSV file.
#' @param subsetRegionIds Optional character vector of region IDs
#'   (\code{"chrom_start_end"}) to subset to.
#' @return A list with \code{ldInfo} and \code{regionInfo} data.frames.
#' @importFrom vroom vroom
#' @export
getCtwasMetaData <- function(ldMetaDataFile, subsetRegionIds = NULL) {
  .Deprecated("ldLoader", package = "pecotmr",
              msg = "getCtwasMetaData() is deprecated. Use ldLoader() with ldInfo instead.")
  ldInfo <- as.data.frame(vroom(ldMetaDataFile))
  colnames(ldInfo)[1] <- "chrom"
  ldInfo$region_id <- paste(as.integer(stripChrPrefix(ldInfo$chrom)),
                            ldInfo$start, ldInfo$end, sep = "_")
  ldInfo$LD_file <- paste0(dirname(ldMetaDataFile), "/",
                           gsub(",.*$", "", ldInfo$path))
  ldInfo$SNP_file <- paste0(ldInfo$LD_file, ".bim")
  ldInfo <- ldInfo[, c("region_id", "LD_file", "SNP_file")]
  regionInfo <- ldInfo[, "region_id", drop = FALSE]
  regionInfo$chrom <- as.integer(gsub("\\_.*$", "", regionInfo$region_id))
  regionInfo$start <- as.integer(gsub("\\_.*$", "",
                                      sub("^.*?\\_", "", regionInfo$region_id)))
  regionInfo$stop <- as.integer(sub("^.*?\\_", "",
                                    sub("^.*?\\_", "", regionInfo$region_id)))
  regionInfo$region_id <- paste0(regionInfo$chrom, "_",
                                 regionInfo$start, "_",
                                 regionInfo$stop)
  regionInfo <- regionInfo[, c("chrom", "start", "stop", "region_id")]
  if (!is.null(subsetRegionIds)) {
    regionInfo <- regionInfo[regionInfo$region_id %in% subsetRegionIds, ]
  }
  list(ldInfo = ldInfo, regionInfo = regionInfo)
}


# -----------------------------------------------------------------------------
# Legacy regional-data loaders → QtlDataset / MultiStudyQtlDataset
# -----------------------------------------------------------------------------

#' (Deprecated) Load regional association data
#'
#' \strong{Deprecated.} Use \code{\link{QtlDataset}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalAssociationData <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalAssociationData() has been removed.",
      "Build a QtlDataset() directly from a GenotypeHandle and a named",
      "list of per-context SummarizedExperiment phenotypes; pass",
      "mafCutoff / macCutoff / xvarCutoff / imissCutoff /",
      "keepSamples / keepVariants to the constructor."))
  invisible(NULL)
}

#' (Deprecated) Load regional univariate data
#'
#' \strong{Deprecated.} Use \code{\link{QtlDataset}} with a single context.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalUnivariateData <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalUnivariateData() has been removed. Build a QtlDataset()",
      "with a single context entry in the phenotypes list."))
  invisible(NULL)
}

#' (Deprecated) Load regional data for regression modeling
#'
#' \strong{Deprecated.} Use \code{\link{QtlDataset}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalRegressionData <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalRegressionData() has been removed. Build a QtlDataset()",
      "directly; per-condition residualized genotype/phenotype views are",
      "available via getResidualizedGenotypes() and",
      "getResidualizedPhenotypes()."))
  invisible(NULL)
}

#' (Deprecated) Load and preprocess regional multivariate data
#'
#' \strong{Deprecated.} Use \code{\link{MultiStudyQtlDataset}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalMultivariateData <- function(...) {
  .Deprecated(new = "MultiStudyQtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalMultivariateData() has been removed. Build per-study",
      "QtlDataset() objects and combine them with MultiStudyQtlDataset();",
      "the multivariate-Y join is now a pipeline-side concern (mvSuSiE",
      "and mr.mash wrappers form the joint Y matrix from the QtlDataset",
      "list at use-time)."))
  invisible(NULL)
}

#' (Deprecated) Load regional functional association data
#'
#' \strong{Deprecated.} Use \code{\link{QtlDataset}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRegionalFunctionalData <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "loadRegionalFunctionalData() has been removed. Build a QtlDataset()",
      "directly; the previous `minMarkers` filter can be applied to the",
      "phenotype SummarizedExperiment list before constructor entry."))
  invisible(NULL)
}

#' (Deprecated) Load TWAS Weights from RDS Files
#'
#' \strong{Deprecated.} Construct a \code{\link{TwasWeights}} collection
#' directly with the \code{TwasWeights()} constructor.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadTwasWeights <- function(...) {
  .Deprecated(
    new = "TwasWeights",
    package = "pecotmr",
    msg = paste(
      "loadTwasWeights() has been removed. Construct TwasWeights",
      "collections directly with TwasWeights() and TwasWeightsEntry();",
      "file-path loaders are no longer part of pecotmr."))
  invisible(NULL)
}

#' (Deprecated) Load Summary Statistic Data
#'
#' \strong{Deprecated.} Build a \code{GRanges} of summary statistics in
#' your own code and pass it to \code{\link{GwasSumStats}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadRssData <- function(...) {
  .Deprecated(new = "GwasSumStats", package = "pecotmr",
    msg = paste(
      "loadRssData() has been removed. Build a GRanges of summary",
      "statistics in your own code (vignette('rss-qc') shows examples",
      "with MungeSumstats / Rsamtools / VariantAnnotation), then pass it",
      "to GwasSumStats()."))
  invisible(NULL)
}

#' (Deprecated) Load mixture regional data across multiple cohorts
#'
#' \strong{Deprecated.} Build per-study \code{\link{QtlDataset}} +
#' \code{\link{QtlSumStats}} and combine with
#' \code{\link{MultiStudyQtlDataset}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadMultitaskRegionalData <- function(...) {
  .Deprecated(new = "MultiStudyQtlDataset", package = "pecotmr",
    msg = paste(
      "loadMultitaskRegionalData() has been removed. Build per-study",
      "QtlDataset() objects for individual-level cohorts and a",
      "QtlSumStats() for summary-statistic cohorts, then combine with",
      "MultiStudyQtlDataset()."))
  invisible(NULL)
}

#' (Deprecated) Convert loaded regional data to individual-level inputs
#'
#' \strong{Deprecated.} Use \code{\link{QtlDataset}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
regionDataToIndInput <- function(...) {
  .Deprecated(new = "QtlDataset", package = "pecotmr",
    msg = paste(
      "regionDataToIndInput() has been removed alongside RegionalData.",
      "Build individual-level inputs directly via the QtlDataset()",
      "constructor."))
  invisible(NULL)
}

#' (Deprecated) Convert loaded regional data to RSS inputs
#'
#' \strong{Deprecated.} Use \code{\link{QtlSumStats}} or
#' \code{\link{GwasSumStats}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
regionDataToRssInput <- function(...) {
  .Deprecated(new = "QtlSumStats", package = "pecotmr",
    msg = paste(
      "regionDataToRssInput() has been removed alongside RegionalData.",
      "Build RSS inputs directly via QtlSumStats() / GwasSumStats() and",
      "run summaryStatsQc()."))
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Legacy TWAS / LD / RSS pipelines → S4 pipelines
# -----------------------------------------------------------------------------

#' (Deprecated) Harmonize TWAS panel data
#'
#' \strong{Deprecated.} TWAS-side panel harmonization is absorbed into
#' \code{\link{causalInferencePipeline}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
harmonizeTwas <- function(...) {
  .Deprecated(new = "causalInferencePipeline", package = "pecotmr",
    msg = paste(
      "harmonizeTwas() has been removed. Variant harmonization for TWAS",
      "is now an internal step of causalInferencePipeline()."))
  invisible(NULL)
}

#' (Deprecated) Harmonize GWAS sumstats against an LD panel
#'
#' \strong{Deprecated.} Use \code{\link{summaryStatsQc}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
harmonizeGwas <- function(...) {
  .Deprecated(new = "summaryStatsQc", package = "pecotmr",
    msg = paste(
      "harmonizeGwas() has been removed. GWAS / LD panel harmonization",
      "now runs inside summaryStatsQc()."))
  invisible(NULL)
}

#' (Deprecated) Load LD for a study, supporting single or mixture panels
#'
#' \strong{Deprecated.} Build a \code{GenotypeHandle} (or load via
#' \code{\link{loadLdMatrix}}) and pass it as the \code{ldSketch} slot of
#' your SumStats input. Mixture-LD support will be reintroduced in a
#' follow-up.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadStudyLd <- function(...) {
  .Deprecated(new = "GenotypeHandle", package = "pecotmr",
    msg = paste(
      "loadStudyLd() has been removed. For single-panel LD, build a",
      "GenotypeHandle (or use loadLdMatrix()) and pass it as the",
      "ldSketch slot of your SumStats. Mixture-LD support is deferred",
      "and will be reintroduced in a follow-up."))
  invisible(NULL)
}

#' (Deprecated) Univariate Analysis Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{fineMappingPipeline}} with
#' \code{methods = "susie"} (and the SuSiE-inf init chain via
#' \code{addSusieInf = TRUE}).
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
univariateAnalysisPipeline <- function(...) {
  .Deprecated(new = "fineMappingPipeline", package = "pecotmr",
    msg = paste(
      "univariateAnalysisPipeline() has been removed. Use",
      "fineMappingPipeline(qtlDataset, methods = 'susie') for",
      "fine-mapping and twasWeightsPipeline(qtlDataset, ...) for",
      "TWAS weights."))
  invisible(NULL)
}

#' (Deprecated) RSS Analysis Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{fineMappingPipeline}} on a
#' \code{\link{QtlSumStats}} or \code{\link{GwasSumStats}} after
#' \code{\link{summaryStatsQc}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
rssAnalysisPipeline <- function(...) {
  .Deprecated(new = "fineMappingPipeline", package = "pecotmr",
    msg = paste(
      "rssAnalysisPipeline() has been removed. Use",
      "fineMappingPipeline(sumStats, methods = 'susieRSS') after",
      "running summaryStatsQc() on the SumStats input."))
  invisible(NULL)
}

#' (Deprecated) Multivariate Analysis Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{fineMappingPipeline}} with
#' \code{methods = "mvsusie"}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
multivariateAnalysisPipeline <- function(...) {
  .Deprecated(new = "fineMappingPipeline", package = "pecotmr",
    msg = paste(
      "multivariateAnalysisPipeline() has been removed. Use",
      "fineMappingPipeline(data, methods = 'mvsusie') for",
      "individual-level joint analyses or",
      "fineMappingPipeline(data, methods = 'mvsusieRSS') for",
      "RSS-based joint analyses."))
  invisible(NULL)
}

#' (Deprecated) SuSiE-RSS Fine-mapping Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{fineMappingPipeline}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
susieRssPipeline <- function(...) {
  .Deprecated(new = "fineMappingPipeline", package = "pecotmr",
    msg = paste(
      "susieRssPipeline() has been removed. Use",
      "fineMappingPipeline(sumStats, methods = 'susieRSS') after",
      "running summaryStatsQc()."))
  invisible(NULL)
}

#' (Deprecated) TWAS Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{ctwasPipeline}} for the cTWAS
#' variant and \code{\link{causalInferencePipeline}} for per-tuple TWAS
#' Z / MR.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
twasPipeline <- function(...) {
  .Deprecated(new = "causalInferencePipeline", package = "pecotmr",
    msg = paste(
      "twasPipeline() has been removed. Use ctwasPipeline() for cTWAS",
      "and causalInferencePipeline() for per-tuple TWAS Z + MR."))
  invisible(NULL)
}

#' (Deprecated) Multivariate TWAS Weights Pipeline
#'
#' \strong{Deprecated.} Use \code{\link{twasWeightsPipeline}}, which now
#' handles both univariate and multivariate (mvSuSiE / mr.mash) weight
#' methods uniformly. Pass \code{methods = c("mrmash", "mvsusie")} to
#' fit them side-by-side; the per-(study, context, trait, method) output
#' shape is the same.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
twasMultivariateWeightsPipeline <- function(...) {
  .Deprecated(new = "twasWeightsPipeline", package = "pecotmr",
    msg = paste(
      "twasMultivariateWeightsPipeline() has been removed.",
      "Use twasWeightsPipeline(data, methods = c('mrmash', 'mvsusie'))",
      "to fit mr.mash and mvSuSiE side-by-side; multivariate dispatch",
      "and per-condition output partitioning are built into the unified",
      "pipeline."))
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Sumstat / phenotype I/O helpers → user-side reads
# -----------------------------------------------------------------------------
# Post-S4 refactor, GwasSumStats / QtlSumStats / QtlDataset constructors
# accept in-memory R objects directly. Users are expected to read their
# own files with MungeSumstats / vroom / Rsamtools / VariantAnnotation
# and pass the resulting in-memory objects in. The package no longer
# bundles file-format-specific TSV/tabix loaders.

#' (Deprecated) Standardize GWAS sumstats column names
#'
#' \strong{Deprecated.} Use \pkg{MungeSumstats}'s
#' \code{standardise_header()} directly (or your own column renaming) on
#' the data.frame before passing it to \code{\link{GwasSumStats}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
standardiseSumstatsColumns <- function(...) {
  .Deprecated(new = "MungeSumstats::standardise_header", package = "pecotmr",
    msg = paste(
      "standardiseSumstatsColumns() has been removed. Call",
      "MungeSumstats::standardise_header() directly on your data.frame",
      "and rename columns to pecotmr conventions in your own code",
      "before passing the data.frame to GwasSumStats()."))
  invisible(NULL)
}

#' (Deprecated) Load and filter tabular sumstats by region
#'
#' \strong{Deprecated.} Read your tabix-indexed or plain TSV files
#' yourself (Rsamtools::scanTabix, vroom, readr::read_delim), then build
#' a \code{GRanges} of sumstats and pass it to \code{\link{GwasSumStats}}.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadTsvRegion <- function(...) {
  .Deprecated(new = "GwasSumStats", package = "pecotmr",
    msg = paste(
      "loadTsvRegion() has been removed. Read your sumstats TSV (or",
      "tabix .gz) in your own code (Rsamtools::scanTabix() or vroom()",
      "for region subsetting), build a GRanges with the required mcols,",
      "and pass it to GwasSumStats()."))
  invisible(NULL)
}

#' (Deprecated) Split loaded TWAS weights into memory-bounded batches
#'
#' \strong{Deprecated.} The on-disk batching utility is no longer part of
#' pecotmr; construct \code{\link{TwasWeights}} collections in memory
#' and partition them in your own code if you need to bound memory
#' across many genes.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
batchLoadTwasWeights <- function(...) {
  .Deprecated(new = "TwasWeights", package = "pecotmr",
    msg = paste(
      "batchLoadTwasWeights() has been removed. Build TwasWeights",
      "collections in memory and partition them in your own code if",
      "needed; pecotmr no longer ships an on-disk batching utility."))
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Legacy mash file-reading helpers → QtlSumStats / GwasSumStats + mashPipeline
# -----------------------------------------------------------------------------
# Both helpers below pre-date the S4 sumstats refactor. They read tabix-
# indexed tensorQTL / R-format sumstats files directly from disk and
# assemble the per-region data.frame structure mash needs. The new
# workflow is: the user constructs QtlSumStats / GwasSumStats objects in
# their own code (loading files however they like) and passes a named
# list of those S4 objects to `mashPipeline()`.

#' (Deprecated) Load multi-trait tensorQTL summary statistics by region
#'
#' \strong{Deprecated.} Build \code{\link{QtlSumStats}} objects from your
#' tensorQTL output (read the files in your own code with
#' \pkg{Rsamtools} / \pkg{vroom} / \pkg{MungeSumstats}), then pass a
#' named list of them to \code{\link{mashPipeline}} (with the required
#' \code{strong} and \code{random} entries, and optional \code{null}).
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadMultitraitTensorqtlSumstat <- function(...) {
  .Deprecated(new = "mashPipeline", package = "pecotmr",
    msg = paste(
      "loadMultitraitTensorqtlSumstat() has been removed. Build",
      "QtlSumStats() objects from your tensorQTL files in your own code,",
      "then pass a named list (strong / random / optional null) to",
      "mashPipeline()."))
  invisible(NULL)
}

#' (Deprecated) Load multi-trait R-format summary statistics from a SuSiE fit
#'
#' \strong{Deprecated.} Use \code{\link{QtlSumStats}} (or
#' \code{\link{GwasSumStats}}) directly. The legacy file-path workflow
#' is no longer maintained; \code{summaryStatsQc()} handles allele
#' harmonization / LD-panel filtering on the in-memory S4 objects.
#'
#' @param ... Ignored.
#' @return \code{NULL} (invisibly).
#' @export
loadMultitraitRSumstat <- function(...) {
  .Deprecated(new = "mashPipeline", package = "pecotmr",
    msg = paste(
      "loadMultitraitRSumstat() has been removed. Build QtlSumStats() /",
      "GwasSumStats() objects from your inputs, run summaryStatsQc()",
      "for harmonization, then pass them to mashPipeline()."))
  invisible(NULL)
}
