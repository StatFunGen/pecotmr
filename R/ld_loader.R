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
#' @param R_list List of G precomputed LD correlation matrices (p_g x p_g).
#' @param X_list List of G genotype matrices (n x p_g).
#' @param ld_meta_path Path to a pecotmr LD metadata TSV file (as used by
#'   \code{\link{load_LD_matrix}}).
#' @param regions Character vector of G region strings (e.g.,
#'   \code{"chr22:17238266-19744294"}). Required when \code{ld_meta_path}
#'   is used.
#' @param LD_info A data.frame with column \code{LD_file} (paths to
#'   \code{.cor.xz} LD matrix files) and optionally \code{SNP_file}
#'   (paths to companion \code{.bim} files; defaults to
#'   \code{paste0(LD_file, ".bim")} if absent). As returned by
#'   cTWAS meta-data utilities.
#' @param return_genotype Logical. When using region mode, return the
#'   genotype matrix X (\code{TRUE}) or LD correlation R (\code{FALSE},
#'   default).
#' @param max_variants Integer or \code{NULL}. If set, randomly subsample
#'   blocks larger than this to control memory usage.
#'
#' @return A function \code{loader(g)} that, given a block index \code{g},
#'   returns the corresponding LD matrix or genotype matrix.
#'
#' @examples
#' # List mode with pre-computed LD
#' R1 <- diag(10)
#' R2 <- diag(15)
#' loader <- ld_loader(R_list = list(R1, R2))
#' loader(1)  # returns R1
#' loader(2)  # returns R2
#'
#' @export
ld_loader <- function(R_list = NULL, X_list = NULL,
                      ld_meta_path = NULL, regions = NULL,
                      LD_info = NULL,
                      return_genotype = FALSE,
                      max_variants = NULL) {
  # Validate: exactly one source
  n_sources <- sum(!is.null(R_list), !is.null(X_list),
                   !is.null(ld_meta_path), !is.null(LD_info))
  if (n_sources != 1)
    stop("Provide exactly one of R_list, X_list, ld_meta_path, or LD_info.")

  if (!is.null(R_list)) {
    # List mode (R matrices)
    loader <- function(g) {
      R <- R_list[[g]]
      if (!is.null(max_variants) && ncol(R) > max_variants) {
        keep <- sort(sample(ncol(R), max_variants))
        R <- R[keep, keep]
      }
      R
    }
  } else if (!is.null(X_list)) {
    # List mode (genotype matrices)
    loader <- function(g) {
      X <- X_list[[g]]
      if (!is.null(max_variants) && ncol(X) > max_variants) {
        keep <- sort(sample(ncol(X), max_variants))
        X <- X[, keep]
      }
      X
    }
  } else if (!is.null(ld_meta_path)) {
    # Region mode: load on the fly via load_LD_matrix()
    if (is.null(regions))
      stop("'regions' is required when using ld_meta_path.")

    loader <- function(g) {
      ld <- load_LD_matrix(ld_meta_path, region = regions[g],
                           return_genotype = return_genotype)
      mat <- ld$LD_matrix
      if (!is.null(max_variants) && ncol(mat) > max_variants) {
        keep <- sort(sample(ncol(mat), max_variants))
        if (return_genotype || nrow(mat) > ncol(mat)) {
          mat <- mat[, keep]
        } else {
          mat <- mat[keep, keep]
        }
      }
      # Center and scale genotype matrices
      if (return_genotype || nrow(mat) > ncol(mat)) {
        mat <- scale(mat)
        mat[is.na(mat)] <- 0
      }
      mat
    }
  } else {
    # LD_info mode: load LD blocks by index from file paths
    # Supports all three formats:
    #   1. Pre-computed .cor.xz + .bim/.pvar (custom block format)
    #   2. PLINK1 prefix (.bed/.bim/.fam) — LD computed on the fly
    #   3. PLINK2 prefix (.pgen/.pvar/.psam) — LD computed on the fly
    if (!is.data.frame(LD_info) || !"LD_file" %in% colnames(LD_info))
      stop("LD_info must be a data.frame with column 'LD_file'.")

    loader <- function(g) {
      ld_path <- LD_info$LD_file[g]

      # Auto-detect format by checking what files exist
      if (has_plink2_files(ld_path)) {
        # PLINK2: load genotypes and compute LD
        geno <- load_genotype_region(ld_path)
        mat <- compute_LD(geno)
      } else if (has_plink1_files(ld_path)) {
        # PLINK1: load genotypes and compute LD
        geno <- load_genotype_region(ld_path)
        mat <- compute_LD(geno)
      } else {
        # Pre-computed .cor.xz block
        snp_file <- if ("SNP_file" %in% colnames(LD_info)) {
          LD_info$SNP_file[g]
        } else {
          NULL  # let process_LD_matrix auto-detect .bim/.pvar/.pvar.zst
        }
        ld <- process_LD_matrix(ld_path, snp_file)
        mat <- ld$LD_matrix
      }

      if (!is.null(max_variants) && ncol(mat) > max_variants) {
        keep <- sort(sample(ncol(mat), max_variants))
        mat <- mat[keep, keep]
      }
      mat
    }
  }

  loader
}
