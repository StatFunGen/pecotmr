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
