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
#' statistics. Suitable for use with \code{\link{xqtl_enrichment_wrapper}} and
#' \code{\link{computeQtlEnrichment}}. All variant identifiers are synthetic.
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
