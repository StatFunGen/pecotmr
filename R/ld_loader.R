#' Extract the LD or genotype matrix from an LdData S4 object.
#' @param ld An LdData object.
#' @param wantGenotype Logical; if TRUE, extract the genotype matrix
#'   (via \code{getGenotypes()}).
#' @return A matrix.
#' @noRd
extractLdMatrix <- function(ld, wantGenotype = FALSE) {
  if (!is(ld, "LdData")) stop("ld must be an LdData object")
  if (wantGenotype && hasGenotypes(ld)) {
    return(getGenotypes(ld))
  }
  getCorrelation(ld)
}

#' Create an LD loader for on-demand block-wise LD retrieval
#'
#' Constructs a loader function that retrieves per-block LD matrices on
#' demand. This avoids loading all blocks into memory simultaneously,
#' which is critical for genome-wide analyses with hundreds of blocks.
#'
#' Four modes are supported:
#'
#' \describe{
#'   \item{list mode (R)}{Pre-loaded list of LD correlation matrices.
#'     Simple but uses more memory. Set \code{R_list}.}
#'   \item{list mode (X)}{Pre-loaded list of genotype matrices (n x p_g).
#'     Set \code{X_list}.}
#'   \item{region mode}{Loads LD from a pecotmr metadata TSV file on the fly
#'     via \code{\link{load_LD_matrix}}. Memory-efficient for large datasets.
#'     Set \code{ld_meta_path} and \code{regions}.}
#'   \item{LD_info mode}{Loads pre-computed LD blocks from \code{.cor.xz}
#'     files listed in an \code{LD_info} data.frame (as returned by
#'     cTWAS meta-data utilities). Set \code{LD_info}.}
#' }
#'
#' @param rList List of G precomputed LD correlation matrices (p_g x p_g).
#' @param xList List of G genotype matrices (n x p_g).
#' @param ldMetaPath Path to a pecotmr LD metadata TSV file (as used by
#'   \code{\link{load_LD_matrix}}).
#' @param regions Character vector of G region strings (e.g.,
#'   \code{"chr22:17238266-19744294"}). Required when \code{ldMetaPath}
#'   is used.
#' @param ldInfo A data.frame with column \code{LD_file} (paths to
#'   genotype files or \code{.cor.xz} LD matrix files) and optionally
#'   \code{SNP_file} (paths to companion \code{.bim} files for pre-computed
#'   blocks; defaults to \code{paste0(LD_file, ".bim")} if absent).
#'   Genotype paths can be PLINK2 prefixes, PLINK1 prefixes, VCF files,
#'   or GDS files. As returned by cTWAS meta-data utilities.
#' @param returnGenotype Logical. When using region mode, return the
#'   genotype matrix X (\code{TRUE}) or LD correlation R (\code{FALSE},
#'   default).
#' @param maxVariants Integer or \code{NULL}. If set, randomly subsample
#'   blocks larger than this to control memory usage.
#'
#' @return A function \code{loader(g)} that, given a block index \code{g},
#'   returns the corresponding LD matrix or genotype matrix.
#'
#' @examples
#' # List mode with pre-computed LD
#' R1 <- diag(10)
#' R2 <- diag(15)
#' loader <- ldLoader(rList = list(R1, R2))
#' loader(1)  # returns R1
#' loader(2)  # returns R2
#'
#' @export
ldLoader <- function(rList = NULL, xList = NULL,
                     ldMetaPath = NULL, regions = NULL,
                     ldInfo = NULL,
                     returnGenotype = FALSE,
                     maxVariants = NULL) {
  # Validate: exactly one source
  nSources <- sum(!is.null(rList), !is.null(xList),
                  !is.null(ldMetaPath), !is.null(ldInfo))
  if (nSources != 1)
    stop("Provide exactly one of rList, xList, ldMetaPath, or ldInfo.")

  if (!is.null(rList)) {
    # List mode (R matrices)
    loader <- function(g) {
      R <- rList[[g]]
      if (!is.null(maxVariants) && ncol(R) > maxVariants) {
        keep <- sort(sample(ncol(R), maxVariants))
        R <- R[keep, keep]
      }
      R
    }
  } else if (!is.null(xList)) {
    # List mode (genotype matrices)
    loader <- function(g) {
      X <- xList[[g]]
      if (!is.null(maxVariants) && ncol(X) > maxVariants) {
        keep <- sort(sample(ncol(X), maxVariants))
        X <- X[, keep]
      }
      X
    }
  } else if (!is.null(ldMetaPath)) {
    # Region mode: load on the fly via loadLdMatrix()
    if (is.null(regions))
      stop("'regions' is required when using ldMetaPath.")

    loader <- function(g) {
      ld <- loadLdMatrix(ldMetaPath, region = regions[g],
                         returnGenotype = returnGenotype)
      mat <- extractLdMatrix(ld, wantGenotype = returnGenotype)
      if (!is.null(maxVariants) && ncol(mat) > maxVariants) {
        keep <- sort(sample(ncol(mat), maxVariants))
        if (returnGenotype || nrow(mat) > ncol(mat)) {
          mat <- mat[, keep]
        } else {
          mat <- mat[keep, keep]
        }
      }
      # Center and scale genotype matrices
      if (returnGenotype || nrow(mat) > ncol(mat)) {
        mat <- scale(mat)
        mat[is.na(mat)] <- 0
      }
      mat
    }
  } else {
    # ldInfo mode: load LD blocks by index from file paths
    # Supports all genotype formats (PLINK2, PLINK1, VCF, GDS) and
    # pre-computed .cor.xz + .bim/.pvar blocks
    if (!is.data.frame(ldInfo) || !"LD_file" %in% colnames(ldInfo))
      stop("ldInfo must be a data.frame with column 'LD_file'.")

    loader <- function(g) {
      ldPath <- ldInfo$LD_file[g]

      # Auto-detect format: genotype source or pre-computed block
      if (isGenotypeSource(ldPath)) {
        geno <- loadGenotypeRegion(ldPath)
        mat <- computeLd(geno)
      } else {
        # Pre-computed .cor.xz block
        snpFile <- if ("SNP_file" %in% colnames(ldInfo)) {
          ldInfo$SNP_file[g]
        } else {
          NULL  # let processLdMatrix auto-detect .bim/.pvar/.pvar.zst
        }
        ld <- processLdMatrix(ldPath, snpFile)
        mat <- extractLdMatrix(ld)
      }

      if (!is.null(maxVariants) && ncol(mat) > maxVariants) {
        keep <- sort(sample(ncol(mat), maxVariants))
        mat <- mat[keep, keep]
      }
      mat
    }
  }

  loader
}

