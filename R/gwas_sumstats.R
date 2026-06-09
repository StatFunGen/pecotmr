#' @title Summary Statistics Handling
#' @description Functions for reading, validating, and constructing
#'   \code{GWASSumStats} objects from various input formats.
#' @name pecotmr-gwas-sumstats
#' @keywords internal
#' @importFrom GenomicRanges GRanges seqnames start
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom MungeSumstats format_sumstats
#' @include AllGenerics.R
NULL

# =============================================================================
# Constructor
# =============================================================================

#' @title Create a GWASSumStats Object
#' @description Construct a \code{GWASSumStats} from a data.frame with
#'   standardized column names.
#' @param data A data.frame with at minimum columns:
#'   SNP, CHR, BP, A1, A2, Z, N.
#' @param trait_name Character, name for the trait.
#' @param genome Character, genome build (e.g., "hg19", "hg38").
#' @param var_y Numeric, phenotype variance. For case-control studies this is
#'   \code{1 / (phi * (1 - phi))} where \code{phi = n_case / n}. NULL for
#'   quantitative traits.
#' @return A \code{GWASSumStats} object.
#' @export
GWASSumStats <- function(data, trait_name = "trait", genome = "hg19",
                          var_y = NULL) {
  data <- as.data.frame(data)

  required <- c("SNP", "CHR", "BP", "A1", "A2", "Z", "N")
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         ". Consider using readSumstats() with MungeSumStats for ",
         "automatic format detection.")
  }

  data$Z <- as.numeric(data$Z)
  data$N <- as.numeric(data$N)
  data$BP <- as.integer(data$BP)
  data$CHR <- as.character(data$CHR)

  complete <- complete.cases(data[, required])
  if (sum(!complete) > 0) {
    message(sprintf("Removed %d SNPs with missing values in required columns",
                    sum(!complete)))
    data <- data[complete, ]
  }

  chr <- data$CHR
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr <- paste0("chr", chr)

  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = data$BP, width = 1L)
  )

  mcols_data <- data[, c("SNP", "A1", "A2", "Z", "N")]
  optional_cols <- c("MAF", "INFO", "BETA", "SE", "P")
  for (col in optional_cols) {
    if (col %in% colnames(data)) {
      mcols_data[[col]] <- as.numeric(data[[col]])
    }
  }
  mcols(gr) <- DataFrame(mcols_data)

  new("GWASSumStats",
    sumstats = gr,
    genome = genome,
    trait_name = trait_name,
    var_y = var_y
  )
}

# =============================================================================
# Reader method
# =============================================================================

#' @rdname readSumstats
#' @export
setMethod("readSumstats",
  signature(path = "character"),
  function(path, trait_name = "trait", genome = NULL, n = NULL,
           use_mungesumstats = TRUE, ...) {

    if (use_mungesumstats && requireNamespace("MungeSumstats", quietly = TRUE)) {
      message("Standardizing summary statistics with MungeSumStats...")
      reformatted <- format_sumstats(
        path = path,
        ref_genome = genome,
        return_data = TRUE,
        log_folder_ind = FALSE,
        ...
      )
      dt <- as.data.frame(reformatted)

      col_map <- c(
        "SNP" = "SNP", "CHR" = "CHR", "BP" = "BP",
        "A1" = "A1", "A2" = "A2", "Z" = "Z", "N" = "N",
        "FRQ" = "MAF", "INFO" = "INFO", "BETA" = "BETA",
        "SE" = "SE", "P" = "P"
      )
      present <- intersect(names(col_map), colnames(dt))
      names(dt)[match(present, names(dt))] <- col_map[present]

      if (!"Z" %in% colnames(dt) && all(c("BETA", "SE") %in% colnames(dt))) {
        dt$Z <- dt$BETA / dt$SE
      }

      if (!is.null(n) && !"N" %in% colnames(dt)) {
        dt$N <- n
      }

      if (is.null(genome)) genome <- "hg19"

      return(GWASSumStats(dt, trait_name = trait_name, genome = genome))
    }

    message("Reading summary statistics directly (no MungeSumStats)...")
    dt <- as.data.frame(vroom(path, show_col_types = FALSE, ...))

    if (!is.null(n) && !"N" %in% colnames(dt)) {
      dt$N <- n
    }

    if (is.null(genome)) genome <- "hg19"

    GWASSumStats(dt, trait_name = trait_name, genome = genome)
  }
)

# =============================================================================
# Accessor methods
# =============================================================================

#' @rdname getZ
#' @export
setMethod("getZ", "GWASSumStats", function(x) {
  mcols(x@sumstats)$Z
})

#' @rdname getN
#' @export
setMethod("getN", "GWASSumStats", function(x) {
  mcols(x@sumstats)$N
})

#' @rdname getMaf
#' @export
setMethod("getMaf", "GWASSumStats", function(x) {
  mc <- mcols(x@sumstats)
  if ("MAF" %in% colnames(mc)) mc$MAF else NULL
})

#' @rdname nSnps
#' @export
setMethod("nSnps", "GWASSumStats", function(x) {
  length(x@sumstats)
})

#' @rdname subsetChr
#' @export
setMethod("subsetChr", "GWASSumStats", function(x, chr) {
  chr_name <- paste0("chr", sub("^chr", "", as.character(chr)))
  idx <- as.character(seqnames(x@sumstats)) == chr_name
  new("GWASSumStats",
    sumstats = x@sumstats[idx],
    genome = x@genome,
    trait_name = x@trait_name,
    var_y = x@var_y
  )
})

#' @rdname getVarY
#' @export
setMethod("getVarY", "GWASSumStats", function(x) {
  x@var_y
})

# =============================================================================
# Coercion
# =============================================================================

#' @title Convert GWASSumStats to data.frame
#' @description Extracts the genomic ranges and metadata columns into a plain
#'   data.frame with columns SNP, CHR, BP, A1, A2, Z, N (and any optional
#'   columns such as MAF, BETA, SE, P).
#' @param x A \code{GWASSumStats} object.
#' @param row.names Ignored (present for S3 generic compatibility).
#' @param optional Ignored.
#' @param ... Ignored.
#' @return A data.frame.
#' @method as.data.frame GWASSumStats
#' @export
as.data.frame.GWASSumStats <- function(x, row.names = NULL, optional = FALSE, ...) {
  gr <- x@sumstats
  mc <- as.data.frame(mcols(gr))
  mc$CHR <- as.character(seqnames(gr))
  mc$BP  <- start(gr)
  # Reorder: SNP, CHR, BP first, then remaining columns
  first_cols <- c("SNP", "CHR", "BP")
  rest_cols  <- setdiff(names(mc), first_cols)
  mc[, c(first_cols, rest_cols), drop = FALSE]
}

#' @title Convert load_rss_data Output to GWASSumStats
#' @description Converts the list returned by \code{load_rss_data}
#'   (with elements \code{sumstats}, \code{n}, \code{var_y}) into a
#'   \code{GWASSumStats} object.
#' @param rss_list A list with elements \code{sumstats} (data.frame),
#'   \code{n} (numeric or NULL), and \code{var_y} (numeric or NULL), as
#'   returned by \code{load_rss_data}.
#' @param trait_name Character, name for the trait.
#' @param genome Character, genome build (e.g., "hg19", "hg38").
#' @return A \code{GWASSumStats} object, or NULL if \code{rss_list$sumstats}
#'   has zero rows.
#' @export
rss_to_gwas_sumstats <- function(rss_list, trait_name = "trait",
                                  genome = "hg38") {
  ss <- rss_list$sumstats
  if (is.null(ss) || nrow(ss) == 0L) return(NULL)

  # Map pecotmr column names to GWASSumStats expected names
  col_map <- c(
    "chrom" = "CHR", "pos" = "BP", "variant_id" = "SNP",
    "z" = "Z", "beta" = "BETA", "se" = "SE"
  )
  for (old_name in names(col_map)) {
    new_name <- col_map[[old_name]]
    if (old_name %in% names(ss) && !(new_name %in% names(ss))) {
      names(ss)[names(ss) == old_name] <- new_name
    }
  }

  # Ensure SNP column exists
  if (!"SNP" %in% names(ss)) {
    if (all(c("CHR", "BP", "A2", "A1") %in% names(ss))) {
      ss$SNP <- paste(ss$CHR, ss$BP, ss$A2, ss$A1, sep = ":")
    } else {
      ss$SNP <- seq_len(nrow(ss))
    }
  }

  # Ensure N column exists
  if (!"N" %in% names(ss)) {
    ss$N <- rss_list$n
  }

  GWASSumStats(data = ss, trait_name = trait_name, genome = genome,
               var_y = rss_list$var_y)
}
