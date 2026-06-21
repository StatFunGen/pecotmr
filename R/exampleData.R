#' @name gwas_sumstats_example
#'
#' @title Example GWAS Summary Statistics
#'
#' @docType data
#'
#' @description De-identified GWAS summary statistics for a single genomic
#' region. Sample names, variant positions, and identifiers have been
#' randomized; they do not correspond to any real locus or study.
#'
#' @format A data frame with 2,828 rows and 8 columns:
#'
#' \describe{
#'   \item{variant_id}{Character. Synthetic variant identifier
#'     (chrom:pos:A1:A2).}
#'   \item{chrom}{Character. Chromosome label.}
#'   \item{pos}{Integer. Genomic position (synthetic).}
#'   \item{A1}{Character. Effect allele.}
#'   \item{A2}{Character. Other allele.}
#'   \item{beta}{Numeric. GWAS effect size estimate.}
#'   \item{se}{Numeric. Standard error of the effect size.}
#'   \item{z}{Numeric. Z-score (beta / se).}
#' }
#'
#' @keywords data
#'
#' @examples
#' data(gwas_sumstats_example)
#' head(gwas_sumstats_example)
#'
NULL


#' @name eqtl_region_example
#'
#' @title Example eQTL Region Data (Individual-Level)
#'
#' @docType data
#'
#' @description De-identified individual-level eQTL data for a single genomic
#' region, containing a genotype matrix and residualized phenotype vector.
#' All sample names, variant positions, and identifiers are synthetic and do
#' not correspond to any real individuals or loci.
#'
#' @format A list with two elements:
#'
#' \describe{
#'   \item{X}{Numeric matrix (415 samples x 2,828 variants). Genotype dosage
#'     matrix with synthetic sample and variant names.}
#'   \item{y_res}{Named numeric vector (length 415). Residualized molecular
#'     phenotype values with synthetic sample names.}
#' }
#'
#' @keywords data
#'
#' @examples
#' data(eqtl_region_example)
#' dim(eqtl_region_example$X)
#' length(eqtl_region_example$y_res)
#'
NULL


#' @name gwas_finemapping_example
#'
#' @title Example GWAS Fine-Mapping Results (SuSiE)
#'
#' @docType data
#'
#' @description SuSiE RSS fine-mapping results for de-identified GWAS summary
#' statistics. Suitable for use with \code{\link{qtlEnrichment}}. All variant
#' identifiers are synthetic.
#'
#' @format A list of length 1, where the first element is a trimmed SuSiE
#' result list containing:
#'
#' \describe{
#'   \item{alpha}{Numeric matrix (L x p). Single-effect assignment
#'     probabilities.}
#'   \item{pip}{Named numeric vector (length 2,828). Posterior inclusion
#'     probabilities with synthetic variant names.}
#'   \item{V}{Numeric vector. Estimated prior variances for each single
#'     effect.}
#'   \item{sets}{List. Credible set information from SuSiE.}
#' }
#'
#' @keywords data
#'
#' @examples
#' data(gwas_finemapping_example)
#' length(gwas_finemapping_example[[1]]$pip)
#'
NULL


#' @name qtl_finemapping_example
#'
#' @title Example QTL Fine-Mapping Results (SuSiE)
#'
#' @docType data
#'
#' @description SuSiE fine-mapping results for de-identified eQTL
#' individual-level data. Stored as a nested list matching the format expected
#' by \code{\link{xqtl_enrichment_wrapper}}. All variant identifiers and
#' region/context names are synthetic.
#'
#' @format A nested list with structure
#' \code{[[region]][[context]]}, where each context contains:
#'
#' \describe{
#'   \item{susie_result_trimmed}{List. Trimmed SuSiE result with elements
#'     \code{alpha}, \code{pip}, \code{V}, and \code{sets}. (Legacy key; new
#'     code should use the \code{FineMappingResult} S4 object.)}
#'   \item{variantNames}{Character vector (length 2,828). Synthetic variant
#'     identifiers matching the variant names in the SuSiE result.
#'     (Legacy key; new code should use the \code{FineMappingResult} S4 object.)}
#' }
#'
#' @keywords data
#'
#' @examples
#' data(qtl_finemapping_example)
#' names(qtl_finemapping_example)
#' names(qtl_finemapping_example[["region_1"]][["context_1"]])
#'
NULL
#' @name multitraite_data
#'
#' @title Simulated Multi-condition Data for TWAS analysis
#'
#' @docType data
#'
#' @description Simulated data of a gene with multi-conditions (cell-type/tissues)
#' gene expression level matrix(Y) and genotype matrix(X) from 400 individuals,
#' plus mixure prior matrices, prior grid, as well as summary statistics from
#' univariate regression and GWAS summary statistics that is ready for use for
#' TWAS analysis. Genotype matrix is centered and scaled, expression matrix is
#' normalized.
#'
#' @format \code{multitraite_data} is a list with the following elements:
#'
#' \describe{
#'
#'   \item{X}{Centered and scaled n x p matrix of genotype, where n is the total
#'       number of individuals and p denotes the number of SNPs.}
#'
#'   \item{Y}{Normalized n x r matrix of residual for expression, where n is the
#'       total number of individuals and r is the total number of conditions
#'       (tissue/cell-types).}
#'
#'   \item{prior_matrices}{A list of data-driven covariance matrices.}
#'
#'   \item{prior_grid}{A vector of scaling factors to be used in fitting
#'         mr.mash model.}
#'
#'   \item{prior_matrices_cv}{A list of list containing data-driven covariance
#'         matrices for 5-fold cross validation.}
#'
#'   \item{prior_grid_cv}{A list of vectors of scaling factors for 5-fold
#'         cross validation via sample partition.}
#'
#'   \item{gwas_sumstats}{A data frame for GWAS summary statistics.}
#'
#'   \item{sumstat}{Summary statistics of Bhat and Sbhat from univariate
#'         regression for a gene.}
#'
#'    \item{sumstat_cv}{A list of 5 fold cross-validation summary statistics based
#'         on sample partition for a gene.}
#' }
#'
#' @keywords data
#'
#' @references
#' Morgante, F., Carbonetto, P., Wang, G., Zou, Y., Sarkar, A. & Stephens, M. (2023).
#'   A flexible empirical Bayes approach to multivariate multiple regression, and
#'   its improved accuracy in predicting multi-tissue gene expression from genotypes.
#'   PLoS Genetics 19(7): e1010539. https://doi.org/10.1371/journal.pgen.1010539
#'
#' @examples
#' data(multitraite_data)
#'
NULL


#' @name qtl_dataset_example
#'
#' @title Example QtlDataset (S4)
#'
#' @docType data
#'
#' @description A minimal but complete \code{\link{QtlDataset}} built from the
#' bundled \code{inst/extdata/toy_ref} PLINK1 panel (165 samples, 200 chr22
#' variants) plus a synthetic single-trait phenotype with two causal variants.
#' Intended as a self-contained input for the
#' \code{\link{fineMappingPipeline}}, \code{\link{twasWeightsPipeline}}, and
#' \code{\link{colocboostPipeline}} vignettes.
#'
#' @format A \code{QtlDataset} object: study = \code{"study1"}, single context
#'   \code{"brain"}, single trait \code{"ENSG_example"}. Genotype handle wraps
#'   the bundled toy PLINK1 reference; phenotypes are a single-row
#'   \code{SummarizedExperiment} with synthetic dosage-driven values.
#'
#' @keywords data
#'
#' @examples
#' data(qtl_dataset_example)
#' qtl_dataset_example
#'
NULL


#' @name qtl_sumstats_example
#'
#' @title Example QtlSumStats (S4)
#'
#' @docType data
#'
#' @description A pre-QC'd \code{\link{QtlSumStats}} collection covering the
#' same (study, context, trait) tuple as \code{\link{qtl_dataset_example}};
#' the per-variant Z / BETA / SE / N values were computed from the synthetic
#' phenotype by per-variant linear regression. The \code{ldSketch} slot
#' references the same toy PLINK1 panel.
#'
#' @format A \code{QtlSumStats} S4 collection with one row (study1, brain,
#'   ENSG_example) backed by a GRanges of 200 chr22 variants.
#'
#' @keywords data
#'
#' @examples
#' data(qtl_sumstats_example)
#' qtl_sumstats_example
#'
NULL


#' @name qtl_sumstats_multicontext_example
#'
#' @title Example multi-context QtlSumStats (S4) for mash demos
#'
#' @docType data
#'
#' @description A \code{\link{QtlSumStats}} collection covering one trait
#' (\code{"ENSG_example"}) across three synthetic contexts (\code{brain},
#' \code{blood}, \code{muscle}) on the same toy PLINK1 panel used by
#' \code{\link{qtl_sumstats_example}}. The signal pattern is wired for
#' mash pattern recovery: one variant is causal in all three contexts
#' (the shared eQTL), one is brain-only, one is blood-only, and muscle
#' carries only the shared signal. Per-variant \code{BETA}, \code{SE},
#' \code{Z}, \code{N}, and \code{MAF} are populated, so the bundle
#' works with both \code{inputScale = "beta"} (the default) and
#' \code{inputScale = "z"} paths through \code{\link{mashPipeline}}.
#'
#' @format A \code{QtlSumStats} S4 collection with three rows (one per
#'   context) backed by GRanges of 200 chr22 variants each.
#'
#' @keywords data
#'
#' @examples
#' data(qtl_sumstats_multicontext_example)
#' qtl_sumstats_multicontext_example
#' getContexts(qtl_sumstats_multicontext_example)
#'
NULL


#' @name gwas_sumstats_s4_example
#'
#' @title Example GwasSumStats (S4)
#'
#' @docType data
#'
#' @description A pre-QC'd \code{\link{GwasSumStats}} collection for one
#' synthetic trait (\code{"trait1"}, N = 50,000). The signal pattern shares
#' one causal variant with \code{\link{qtl_sumstats_example}} (for
#' colocalization demos) and adds a second GWAS-only causal variant.
#' The \code{ldSketch} slot references the same toy PLINK1 panel used
#' for the QTL example.
#'
#' @format A \code{GwasSumStats} S4 collection with one row (trait1)
#'   backed by a GRanges of 200 chr22 variants.
#'
#' @keywords data
#'
#' @examples
#' data(gwas_sumstats_s4_example)
#' gwas_sumstats_s4_example
#'
NULL


#' @title Resolve bundled-example GenotypeHandle paths
#'
#' @description The bundled S4 example objects (\code{qtl_dataset_example},
#' \code{qtl_sumstats_example}, \code{qtl_sumstats_multicontext_example},
#' \code{gwas_sumstats_s4_example}, \code{multi_study_qtl_dataset_example})
#' store a relative-style path to the bundled
#' \code{inst/extdata/toy_canonical} PLINK1 reference.
#' Call this helper once at the top of a vignette to re-point each
#' contained \code{GenotypeHandle} at the installed path resolved via
#' \code{system.file()}.
#'
#' @param x A bundled \code{QtlDataset}, \code{MultiStudyQtlDataset},
#'   \code{QtlSumStats}, or \code{GwasSumStats} example object.
#' @return The same object with \code{GenotypeHandle@@path} rewritten
#'   to the resolved install path.
#' @export
fixupExampleGenotypePaths <- function(x) {
  resolve <- function(handle) {
    bn <- basename(handle@path)
    # PLINK1 / PLINK2 are stems (no extension); resolve via the .bed
    # (or .pgen) sidecar and strip the extension back off.
    bedPath <- system.file("extdata", paste0(bn, ".bed"),
                            package = "pecotmr")
    if (nzchar(bedPath)) {
      handle@path <- sub("\\.bed$", "", bedPath)
      return(handle)
    }
    new <- system.file("extdata", bn, package = "pecotmr")
    if (!nzchar(new))
      stop("fixupExampleGenotypePaths: cannot resolve '", bn,
           "' under inst/extdata/")
    handle@path <- new
    handle
  }
  if (methods::is(x, "QtlDataset")) {
    x@genotypes <- resolve(x@genotypes)
  } else if (methods::is(x, "MultiStudyQtlDataset")) {
    for (nm in names(x@qtlDatasets))
      x@qtlDatasets[[nm]]@genotypes <-
        resolve(x@qtlDatasets[[nm]]@genotypes)
    if (!is.null(x@sumStats))
      x@sumStats@ldSketch <- resolve(x@sumStats@ldSketch)
  } else if (methods::is(x, "QtlSumStats") ||
             methods::is(x, "GwasSumStats")) {
    x@ldSketch <- resolve(x@ldSketch)
  } else {
    stop("fixupExampleGenotypePaths: unsupported class '",
         class(x)[[1L]], "'.")
  }
  x
}


#' @name multi_study_qtl_dataset_example
#'
#' @title Example MultiStudyQtlDataset (S4)
#'
#' @docType data
#'
#' @description A \code{\link{MultiStudyQtlDataset}} combining two
#' synthetic \code{QtlDataset}s (\code{study1} and \code{study2}) that
#' share the toy PLINK1 panel but have different causal variants. Use it
#' to demonstrate the multi-study dispatch of
#' \code{\link{fineMappingPipeline}}, \code{\link{twasWeightsPipeline}},
#' and \code{\link{colocboostPipeline}}.
#'
#' @format A \code{MultiStudyQtlDataset} with two embedded
#'   \code{QtlDataset}s (study1 and study2) and \code{sumStats = NULL}.
#'
#' @keywords data
#'
#' @examples
#' data(multi_study_qtl_dataset_example)
#' multi_study_qtl_dataset_example
#'
NULL
