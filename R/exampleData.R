#' @name gwasSumStatsExample
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
#' data(gwasSumStatsExample)
#' head(gwasSumStatsExample)
#'
NULL


#' @name eqtlRegionExample
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
#'   \item{yRes}{Named numeric vector (length 415). Residualized molecular
#'     phenotype values with synthetic sample names.}
#' }
#'
#' @keywords data
#'
#' @examples
#' data(eqtlRegionExample)
#' dim(eqtlRegionExample$X)
#' length(eqtlRegionExample$yRes)
#'
NULL


#' @name gwasFineMappingExample
#'
#' @title Example GWAS Fine-Mapping Results (SuSiE)
#'
#' @docType data
#'
#' @description A \code{\link{GwasFineMappingResult}} S4 collection holding
#' SuSiE-RSS fine-mapping output for de-identified GWAS summary statistics (one
#' study), carrying the LD-reference \code{ldSketch}. Suitable as the GWAS input
#' to \code{\link{qtlEnrichmentPipeline}}. All variant identifiers are
#' synthetic.
#'
#' @format A \code{\link{GwasFineMappingResult}} object: a \code{DFrame}-backed
#' collection with one row carrying a \code{\link{FineMappingEntry}} payload
#' (the
#' variant ids, the SuSiE-RSS fit, and a per-variant \code{topLoci} table with
#' \code{variant_id}, \code{pip}, and \code{cs} columns for 2,828 variants) plus
#' a non-\code{NULL} \code{ldSketch} \code{\link{GenotypeHandle}}.
#'
#' @keywords data
#'
#' @examples
#' data(gwasFineMappingExample)
#' gwasFineMappingExample
#' head(getTopLoci(gwasFineMappingExample))
#'
NULL


#' @name qtlFineMappingExample
#'
#' @title Example QTL Fine-Mapping Results (SuSiE)
#'
#' @docType data
#'
#' @description A \code{\link{QtlFineMappingResult}} S4 collection holding SuSiE
#' fine-mapping output for de-identified individual-level eQTL data, for a
#' single
#' \code{(study, context, trait)} tuple. Suitable as the QTL input to
#' \code{\link{qtlEnrichmentPipeline}}. All variant identifiers and
#' region/context names are synthetic.
#'
#' @format A \code{\link{QtlFineMappingResult}} object: a \code{DFrame}-backed
#' collection with one row carrying a \code{\link{FineMappingEntry}} payload
#' (the
#' variant ids, the SuSiE fit, and a per-variant \code{topLoci} table with
#' \code{variant_id}, \code{pip}, and \code{cs} columns for 2,828 variants).
#'
#' @keywords data
#'
#' @examples
#' data(qtlFineMappingExample)
#' qtlFineMappingExample
#' head(getTopLoci(qtlFineMappingExample))
#'
NULL
#' @name multiTraitData
#'
#' @title Simulated Multi-condition Data for TWAS analysis
#'
#' @docType data
#'
#' @description Simulated data of a gene with multi-conditions
#' (cell-type/tissues)
#' gene expression level matrix(Y) and genotype matrix(X) from 400 individuals,
#' plus mixure prior matrices, prior grid, as well as summary statistics from
#' univariate regression and GWAS summary statistics that is ready for use for
#' TWAS analysis. Genotype matrix is centered and scaled, expression matrix is
#' normalized.
#'
#' @format \code{multiTraitData} is a list with the following elements:
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
#'   \item{priorMatrices}{A list of data-driven covariance matrices.}
#'
#'   \item{prior_grid}{A vector of scaling factors to be used in fitting
#'         mr.mash model.}
#'
#'   \item{priorMatricesCv}{A list of list containing data-driven covariance
#'         matrices for 5-fold cross validation.}
#'
#'   \item{prior_grid_cv}{A list of vectors of scaling factors for 5-fold
#'         cross validation via sample partition.}
#'
#'   \item{gwasSumStats}{A data frame for GWAS summary statistics.}
#'
#'   \item{sumstat}{Summary statistics of Bhat and Sbhat from univariate
#'         regression for a gene.}
#'
#'    \item{sumstat_cv}{A list of 5 fold cross-validation summary statistics
#'    based
#'         on sample partition for a gene.}
#' }
#'
#' @keywords data
#'
#' @references
#' Morgante, F., Carbonetto, P., Wang, G., Zou, Y., Sarkar, A. & Stephens, M.
#' (2023).
#'   A flexible empirical Bayes approach to multivariate multiple regression,
#'   and
#'   its improved accuracy in predicting multi-tissue gene expression from
#'   genotypes.
#'   PLoS Genetics 19(7): e1010539. https://doi.org/10.1371/journal.pgen.1010539
#'
#' @examples
#' data(multiTraitData)
#'
NULL


#' @name qtlDatasetExample
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
#' data(qtlDatasetExample)
#' qtlDatasetExample
#'
NULL


#' @name qtlSumStatsExample
#'
#' @title Example QtlSumStats (S4)
#'
#' @docType data
#'
#' @description A pre-QC'd \code{\link{QtlSumStats}} collection covering the
#' same (study, context, trait) tuple as \code{\link{qtlDatasetExample}};
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
#' data(qtlSumStatsExample)
#' qtlSumStatsExample
#'
NULL


#' @name qtlSumStatsMulticontextExample
#'
#' @title Example multi-context QtlSumStats (S4) for mash demos
#'
#' @docType data
#'
#' @description A \code{\link{QtlSumStats}} collection covering one trait
#' (\code{"ENSG_example"}) across three synthetic contexts (\code{brain},
#' \code{blood}, \code{muscle}) on the same toy PLINK1 panel used by
#' \code{\link{qtlSumStatsExample}}. The signal pattern is wired for
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
#' data(qtlSumStatsMulticontextExample)
#' qtlSumStatsMulticontextExample
#' getContexts(qtlSumStatsMulticontextExample)
#'
NULL


#' @name gwasSumStatsS4Example
#'
#' @title Example GwasSumStats (S4)
#'
#' @docType data
#'
#' @description A pre-QC'd \code{\link{GwasSumStats}} collection for one
#' synthetic trait (\code{"trait1"}, N = 50,000). The signal pattern shares
#' one causal variant with \code{\link{qtlSumStatsExample}} (for
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
#' data(gwasSumStatsS4Example)
#' gwasSumStatsS4Example
#'
NULL


#' @name multiStudyQtlDatasetExample
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
#' data(multiStudyQtlDatasetExample)
#' multiStudyQtlDatasetExample
#'
NULL


#' @name qtlFineMappingLbfExample
#'
#' @title Example QTL Fine-Mapping Results with LBF (SuSiE)
#'
#' @docType data
#'
#' @description A \code{\link{QtlFineMappingResult}} whose per-entry SuSiE fit
#' carries the full single-effect log-Bayes-factor matrix
#' (\code{lbf_variable}), posterior inclusion probabilities, and credible sets.
#' Unlike \code{\link{qtlFineMappingExample}} (a lightweight fit), this object
#' supports LBF-based downstream analysis such as \code{\link{colocPipeline}}.
#' Derived from the de-identified \code{protocol_example} toy data
#' (gene ENSG00000283047, chr22); all identifiers are synthetic.
#'
#' @format A \code{QtlFineMappingResult} with two rows (contexts
#'   \code{context1} and \code{context2}) for trait \code{ENSG00000283047},
#'   method \code{susie}. Each entry's \code{susieFit} carries
#'   \code{lbf_variable} (single effects x variants), \code{alpha}, \code{pip},
#'   and credible-set membership.
#'
#' @keywords data
#'
#' @examples
#' data(qtlFineMappingLbfExample)
#' qtlFineMappingLbfExample
#' head(getTopLoci(qtlFineMappingLbfExample))
#'
NULL


#' @name fsusieFineMappingExample
#'
#' @title Example Functional-SuSiE Fine-Mapping Result
#'
#' @docType data
#'
#' @description A \code{\link{QtlFineMappingResult}} produced by the functional
#' SuSiE (\code{fsusie}) method, for demonstrating the \code{fsusie} accessors
#' and weight extraction. Derived from the de-identified
#' \code{protocol_example} toy data; all identifiers are synthetic.
#'
#' @format A \code{QtlFineMappingResult} with one row, method \code{fsusie}.
#'
#' @keywords data
#'
#' @examples
#' data(fsusieFineMappingExample)
#' fsusieFineMappingExample
#'
NULL


#' @name mvsusieFineMappingExample
#'
#' @title Example Multivariate-SuSiE Fine-Mapping Result
#'
#' @docType data
#'
#' @description A \code{\link{QtlFineMappingResult}} produced by multivariate
#' SuSiE (\code{mvsusie}), for demonstrating the multivariate accessors and
#' weight extraction. Derived from the de-identified \code{protocol_example}
#' toy data; all identifiers are synthetic.
#'
#' @format A \code{QtlFineMappingResult} with two rows, method \code{mvsusie}.
#'
#' @keywords data
#'
#' @examples
#' data(mvsusieFineMappingExample)
#' mvsusieFineMappingExample
#'
NULL


#' @name twasWeightsExample
#'
#' @title Example TWAS Weights (Multiple Methods)
#'
#' @docType data
#'
#' @description A \code{\link{TwasWeights}} collection holding per-variant TWAS
#' weights learned by several regression / fine-mapping methods for one gene.
#' Derived from the de-identified \code{protocol_example} toy data
#' (gene ENSG00000130538, chr22); all identifiers are synthetic.
#'
#' @format A \code{TwasWeights} collection with 11 rows (one per method:
#'   \code{mrash}, \code{susie}, \code{susie_inf}, \code{enet}, \code{lasso},
#'   \code{mcp}, \code{scad}, \code{l0learn}, \code{bayes_r}, \code{bayes_c},
#'   \code{ensemble}).
#'
#' @keywords data
#'
#' @examples
#' data(twasWeightsExample)
#' twasWeightsExample
#'
NULL


#' @name ctwasWeightsExample
#'
#' @title Example cTWAS-formatted TWAS Weights
#'
#' @docType data
#'
#' @description A \code{\link{TwasWeights}} collection formatted for the cTWAS
#' workflow (weights, LD, and variant provenance aligned for LD-block
#' placement). Derived from the de-identified \code{protocol_example} toy data
#' (chr22); all identifiers are synthetic.
#'
#' @format A \code{TwasWeights} collection with one row, method \code{susie}.
#'
#' @keywords data
#'
#' @examples
#' data(ctwasWeightsExample)
#' ctwasWeightsExample
#'
NULL


#' @name mashInputExample
#'
#' @title Example mash Input Matrices (Multi-Condition)
#'
#' @docType data
#'
#' @description A named list of \code{mashr}-style input matrices (strong /
#' random / null subsets of effect sizes \code{.b}, standard errors \code{.s},
#' and z-scores \code{.z}) across eight brain cell-type conditions. Drives the
#' mash model workflow (\code{\link{mashModelFit}}, \code{\link{mashPosterior}})
#' via \code{\link{qtlSumStatsFromBetaMatrix}}. Derived from the de-identified
#' \code{protocol_example} toy data; all identifiers are synthetic.
#'
#' @format A list. The effect-size / standard-error / z-score elements
#'   (\code{strong.b}, \code{strong.s}, \code{strong.z}, and the \code{random.*}
#'   / \code{null.*} counterparts) are each a variants x 8-condition matrix
#'   (conditions \code{ALL, Ast, End, Exc, Inh, Mic, OPC, Oli}); the
#'   \code{.b} / \code{.s} / \code{.z} suffixes follow the \code{mashr}
#'   naming contract.
#'
#' @keywords data
#'
#' @examples
#' data(mashInputExample)
#' names(mashInputExample)
#' dim(mashInputExample$strong.b)
#'
NULL


#' @name mashPosteriorExample
#'
#' @title Example mash Posterior Bundle
#'
#' @docType data
#'
#' @description A small bundle of mash posterior outputs on a strong-variant
#' set: the posterior means / covariances, the original effect estimates, and a
#' toy credible-set assignment. Intended for the pairwise-contrast and
#' feature-score helpers. Derived from the de-identified \code{protocol_example}
#' toy data; all identifiers are synthetic.
#'
#' @format A list with elements \code{posterior} (a list of
#'   \code{PosteriorMean} and \code{PosteriorCov}), \code{orig} (a list of
#'   \code{bhat} and \code{sbhat}), and \code{fineMapping} (a data frame of
#'   variant ids, credible-set order, and pip).
#'
#' @keywords data
#'
#' @examples
#' data(mashPosteriorExample)
#' names(mashPosteriorExample)
#'
NULL


#' @name gwasFineMappingLbfExample
#'
#' @title Example GWAS Fine-Mapping Result with LBF (SuSiE)
#'
#' @docType data
#'
#' @description A \code{\link{GwasFineMappingResult}} whose SuSiE fit carries
#' the single-effect log-Bayes-factor matrix (\code{lbf_variable}), covering the
#' same synthetic region as \code{\link{qtlFineMappingLbfExample}}. The pair is
#' the intended input to \code{\link{colocPipeline}} (LBF-based
#' colocalization). Derived from the de-identified \code{protocol_example} toy
#' data; all identifiers are synthetic.
#'
#' @format A \code{GwasFineMappingResult} with one row (study \code{gwas1},
#'   method \code{susie}) whose \code{susieFit} carries \code{lbf_variable},
#'   \code{alpha}, and \code{pip}.
#'
#' @keywords data
#'
#' @examples
#' data(gwasFineMappingLbfExample)
#' gwasFineMappingLbfExample
#'
NULL


#' @name ctwasInputsExample
#'
#' @title Example Assembled cTWAS Inputs
#'
#' @docType data
#'
#' @description The assembled cTWAS input list produced by
#' \code{\link{assembleCtwasInputs}} over a two-LD-block chr22 grid with one
#' gene, ready for \code{\link{estCtwasParam}}. Self-contained (no external LD
#' files needed for parameter estimation). Derived from the de-identified
#' \code{protocol_example} toy data; all identifiers are synthetic.
#'
#' @format A list with the cTWAS input components \code{z_snp}, \code{z_gene},
#'   \code{weights}, \code{region_info}, \code{snp_map}, \code{LD_map}, and the
#'   LD / SNP-info loader closures.
#'
#' @keywords data
#'
#' @examples
#' data(ctwasInputsExample)
#' names(ctwasInputsExample)
#'
NULL


#' @name ctwasEstExample
#'
#' @title Example cTWAS Parameter-Estimation State
#'
#' @docType data
#'
#' @description The augmented cTWAS state returned by
#' \code{\link{estCtwasParam}} (the assembled inputs plus estimated
#' group priors and per-region data), ready for
#' \code{\link{screenCtwasRegions}}. Derived from the de-identified
#' \code{protocol_example} toy data; all identifiers are synthetic.
#'
#' @format A list carrying the \code{\link{ctwasInputsExample}} components plus
#'   the estimated \code{group_prior} / \code{group_prior_var} and per-region
#'   data.
#'
#' @keywords data
#'
#' @examples
#' data(ctwasEstExample)
#' names(ctwasEstExample)
#'
NULL


#' @name ctwasFinemapExample
#'
#' @title Example cTWAS Fine-Mapping Result
#'
#' @docType data
#'
#' @description The fine-mapped cTWAS result returned by
#' \code{\link{finemapCtwasRegions}}, ready for \code{\link{asCtwasResult}} or
#' \code{\link{mergeCtwasBoundaryRegions}}. Derived from the de-identified
#' \code{protocol_example} toy data; all identifiers are synthetic.
#'
#' @format A list of per-region cTWAS fine-mapping output (susie alpha / PIP
#'   tables and region provenance).
#'
#' @keywords data
#'
#' @examples
#' data(ctwasFinemapExample)
#' asCtwasResult(ctwasFinemapExample)
#'
NULL


#' @name ldEigenExample
#'
#' @title Example LD Eigen-Decomposition Reference
#'
#' @docType data
#'
#' @description A small synthetic \code{LdEigen} reference (per-block LD
#' eigen-decompositions over 20 SNPs in 2 blocks) for the heritability
#' estimators (\code{\link{estimateH2}} with \code{method = "lder"} or
#' \code{"hdl"}, and \code{\link{computeLdScores}}). All identifiers are
#' synthetic.
#'
#' @format An \code{LdEigen} object (extends \code{LdStatistic}) with an
#'   \code{eigenList} of per-block \code{values} / \code{vectors} /
#'   \code{snpIdx}
#'   over 20 SNPs in 2 LD blocks.
#'
#' @keywords data
#'
#' @examples
#' data(ldEigenExample)
#' ldEigenExample
#'
NULL


#' @name ldScoreExample
#'
#' @title Example LD-Score Reference
#'
#' @docType data
#'
#' @description A small synthetic \code{LdScore} reference (per-SNP LD scores,
#' weights, and per-block LD matrices over 20 SNPs in 2 blocks) for the
#' heritability estimators (\code{\link{estimateH2}} with
#' \code{method = "gldsc"}). All identifiers are synthetic.
#'
#' @format An \code{LdScore} object (extends \code{LdStatistic}) carrying
#'   \code{ldScores}, \code{ldScoreWeights}, and \code{ldMatrixList} over 20
#'   SNPs
#'   in 2 LD blocks.
#'
#' @keywords data
#'
#' @examples
#' data(ldScoreExample)
#' ldScoreExample
#'
NULL


#' @name h2EstimateExample
#'
#' @title Example Heritability Estimate
#'
#' @docType data
#'
#' @description A small synthetic \code{H2Estimate} (the output of
#' \code{\link{estimateH2}}) carrying the SNP heritability, intercept, and a
#' per-annotation enrichment table. Intended for the \code{H2Estimate}
#' accessors and \code{\link{h2EstimateToSldscTrait}}. All identifiers are
#' synthetic.
#'
#' @format An \code{H2Estimate} object with \code{h2}, \code{h2Se},
#'   \code{intercept}, an \code{enrichment} data frame, and \code{tauBlocks}.
#'
#' @keywords data
#'
#' @examples
#' data(h2EstimateExample)
#' h2EstimateExample
#'
NULL
