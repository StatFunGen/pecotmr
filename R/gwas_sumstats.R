#' @title Summary Statistics Handling
#' @description Functions for reading, validating, and constructing
#'   \code{GwasSumStats} objects from various input formats.
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

#' @title Create a GwasSumStats Object
#' @description Construct a \code{GwasSumStats} from a data.frame with
#'   standardized column names.
#' @param data A data.frame with at minimum columns:
#'   SNP, CHR, BP, A1, A2, Z, N.
#' @param traitName Character, name for the trait.
#' @param genome Character, genome build (e.g., "hg19", "hg38").
#' @param varY Numeric, phenotype variance. For observed-scale OLS on a
#'   centered 0/1 case-control trait, this is \code{n / (n - 1) * phi *
#'   (1 - phi)}, where \code{phi = nCase / n}. Use it only with the full
#'   \code{bhat/shat/var_y} sufficient-statistic interface; z-score RSS
#'   analyses should leave it NULL.
#' @return A \code{GwasSumStats} object.
#' @export
GwasSumStats <- function(data, traitName = "trait", genome = "hg19",
                         varY = NULL) {
  data <- as.data.frame(data)

  required <- c("SNP", "CHR", "BP", "A1", "A2", "Z", "N")
  missingCols <- setdiff(required, colnames(data))
  if (length(missingCols) > 0) {
    stop("Missing required columns: ", paste(missingCols, collapse = ", "),
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

  mcolsData <- data[, c("SNP", "A1", "A2", "Z", "N")]
  optionalCols <- c("MAF", "INFO", "BETA", "SE", "P")
  for (col in optionalCols) {
    if (col %in% colnames(data)) {
      mcolsData[[col]] <- as.numeric(data[[col]])
    }
  }
  mcols(gr) <- DataFrame(mcolsData)

  new("GwasSumStats",
    sumstats = gr,
    genome = genome,
    traitName = traitName,
    varY = varY
  )
}

# =============================================================================
# Reader method
# =============================================================================

#' @rdname readSumstats
#' @export
setMethod("readSumstats",
  signature(path = "character"),
  function(path, traitName = "trait", genome = NULL, n = NULL,
           useMungesumstats = TRUE, ...) {

    if (useMungesumstats && requireNamespace("MungeSumstats", quietly = TRUE)) {
      message("Standardizing summary statistics with MungeSumStats...")
      reformatted <- format_sumstats(
        path = path,
        ref_genome = genome,
        return_data = TRUE,
        log_folder_ind = FALSE,
        ...
      )
      dt <- as.data.frame(reformatted)

      colMap <- c(
        "SNP" = "SNP", "CHR" = "CHR", "BP" = "BP",
        "A1" = "A1", "A2" = "A2", "Z" = "Z", "N" = "N",
        "FRQ" = "MAF", "INFO" = "INFO", "BETA" = "BETA",
        "SE" = "SE", "P" = "P"
      )
      present <- intersect(names(colMap), colnames(dt))
      names(dt)[match(present, names(dt))] <- colMap[present]

      if (!"Z" %in% colnames(dt) && all(c("BETA", "SE") %in% colnames(dt))) {
        dt$Z <- dt$BETA / dt$SE
      }

      if (!is.null(n) && !"N" %in% colnames(dt)) {
        dt$N <- n
      }

      if (is.null(genome)) genome <- "hg19"

      return(GwasSumStats(dt, traitName = traitName, genome = genome))
    }

    message("Reading summary statistics directly (no MungeSumStats)...")
    dt <- as.data.frame(vroom(path, show_col_types = FALSE, ...))

    if (!is.null(n) && !"N" %in% colnames(dt)) {
      dt$N <- n
    }

    if (is.null(genome)) genome <- "hg19"

    GwasSumStats(dt, traitName = traitName, genome = genome)
  }
)

# =============================================================================
# Accessor methods
# =============================================================================

#' @rdname getZ
#' @export
setMethod("getZ", "GwasSumStats", function(x) {
  mcols(x@sumstats)$Z
})

#' @rdname getN
#' @export
setMethod("getN", "GwasSumStats", function(x) {
  mcols(x@sumstats)$N
})

#' @rdname getMaf
#' @export
setMethod("getMaf", "GwasSumStats", function(x) {
  mc <- mcols(x@sumstats)
  if ("MAF" %in% colnames(mc)) mc$MAF else NULL
})

#' @rdname nSnps
#' @export
setMethod("nSnps", "GwasSumStats", function(x) {
  length(x@sumstats)
})

#' @rdname subsetChr
#' @export
setMethod("subsetChr", "GwasSumStats", function(x, chr) {
  chrName <- paste0("chr", sub("^chr", "", as.character(chr)))
  idx <- as.character(seqnames(x@sumstats)) == chrName
  new("GwasSumStats",
    sumstats = x@sumstats[idx],
    genome = x@genome,
    traitName = x@traitName,
    varY = x@varY
  )
})

#' @rdname getVarY
#' @export
setMethod("getVarY", "GwasSumStats", function(x) {
  x@varY
})

# =============================================================================
# Coercion
# =============================================================================

#' @title Convert GwasSumStats to data.frame
#' @description Extracts the genomic ranges and metadata columns into a plain
#'   data.frame with columns SNP, CHR, BP, A1, A2, Z, N (and any optional
#'   columns such as MAF, BETA, SE, P).
#' @param x A \code{GwasSumStats} object.
#' @param row.names Ignored (present for S3 generic compatibility).
#' @param optional Ignored.
#' @param ... Ignored.
#' @return A data.frame.
#' @method as.data.frame GwasSumStats
#' @export
as.data.frame.GwasSumStats <- function(x, row.names = NULL, optional = FALSE, ...) {
  gr <- x@sumstats
  mc <- as.data.frame(mcols(gr))
  mc$CHR <- as.character(seqnames(gr))
  mc$BP  <- start(gr)
  # Reorder: SNP, CHR, BP first, then remaining columns
  firstCols <- c("SNP", "CHR", "BP")
  restCols  <- setdiff(names(mc), firstCols)
  mc[, c(firstCols, restCols), drop = FALSE]
}

#' @title Convert load_rss_data Output to GwasSumStats
#' @description Converts the list returned by \code{load_rss_data}
#'   (with elements \code{sumstats}, \code{n}, \code{var_y}) into a
#'   \code{GwasSumStats} object.
#' @param rssList A list with elements \code{sumstats} (data.frame),
#'   \code{n} (numeric or NULL), and \code{var_y} (numeric or NULL), as
#'   returned by \code{load_rss_data}.
#' @param traitName Character, name for the trait.
#' @param genome Character, genome build (e.g., "hg19", "hg38").
#' @return A \code{GwasSumStats} object, or NULL if \code{rssList$sumstats}
#'   has zero rows.
#' @export
rssToGwasSumstats <- function(rssList, traitName = "trait",
                              genome = "hg38") {
  ss <- rssList$sumstats
  if (is.null(ss) || nrow(ss) == 0L) return(NULL)

  # Map pecotmr column names to GwasSumStats expected names
  colMap <- c(
    "chrom" = "CHR", "pos" = "BP", "variant_id" = "SNP",
    "z" = "Z", "beta" = "BETA", "se" = "SE"
  )
  for (oldName in names(colMap)) {
    newName <- colMap[[oldName]]
    if (oldName %in% names(ss) && !(newName %in% names(ss))) {
      names(ss)[names(ss) == oldName] <- newName
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
    ss$N <- rssList$n
  }

  GwasSumStats(data = ss, traitName = traitName, genome = genome,
               varY = rssList$varY)
}

