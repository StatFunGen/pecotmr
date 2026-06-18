

#' Detect LD-Summary Statistic Mismatches
#'
#' Unified wrapper for detecting outlier variants due to LD-summary statistic
#' mismatches. Dispatches to either \code{\link{dentistSingleWindow}} or
#' \code{\link{slalom}} based on the \code{method} argument.
#'
#' @param zScore Numeric vector of z-scores.
#' @param R Square LD correlation matrix. Provide either \code{R} or \code{X}.
#' @param X Genotype matrix (samples x SNPs). If provided, LD is computed via
#'   \code{\link{computeLd}} and \code{nSample} defaults to \code{nrow(X)}.
#' @param nSample Number of samples in the LD reference panel. Required when
#'   \code{R} is provided and \code{method = "dentist"}; inferred from \code{X}
#'   when \code{X} is provided.
#' @param method Character string specifying the QC method: \code{"slalom"}
#'   (default) or \code{"dentist"}.
#' @param ldMethod Character string specifying the LD computation method when
#'   \code{X} is provided. One of \code{"sample"} (default), \code{"population"},
#'   or \code{"gcta"}. Ignored when \code{R} is provided directly.
#' @param ... Additional arguments passed to the underlying QC method
#'   (\code{\link{dentistSingleWindow}} or \code{\link{slalom}}).
#'
#' @return A data frame with at least a logical \code{outlier} column indicating
#'   which variants are identified as outliers. The remaining columns depend on
#'   the method used.
#'
#' @seealso \code{\link{dentistSingleWindow}}, \code{\link{slalom}},
#'   \code{\link{summaryStatsQc}}
#' @importFrom dplyr mutate row_number filter pull
#' @export
ldMismatchQc <- function(zScore, R = NULL, X = NULL, nSample = NULL,
                         method = c("slalom", "dentist"),
                         ldMethod = "sample", ...) {
  method <- match.arg(method)
  if (method == "dentist") {
    qcResults <- dentistSingleWindow(zScore, R = R, X = X, nSample = nSample,
                                     ldMethod = ldMethod, ...)
    return(qcResults)
  } else {
    qcResults <- slalom(zScore, R = R, X = X, ldMethod = ldMethod, ...)
    # Standardize output: slalom uses "outliers", rename to "outlier" for consistency
    result <- qcResults$data
    if ("outliers" %in% colnames(result) && !"outlier" %in% colnames(result)) {
      colnames(result)[colnames(result) == "outliers"] <- "outlier"
    }
    return(result)
  }
}

.resolveZMismatchQc <- function(zMismatchQc) {
  if (is.null(zMismatchQc)) return("none")
  match.arg(zMismatchQc, c("none", "slalom", "dentist"))
}

#' Kriging-style LD-consistency outlier QC
#'
#' Flags variants whose observed z-score is inconsistent with the value
#' predicted from its LD neighbours. For \code{z ~ N(0, R)} the leave-one-out
#' conditional distribution of \code{z_i} given the rest has mean
#' \code{-(1/Omega_ii) * Omega_{i,-i} z_{-i}} and variance \code{1/Omega_ii},
#' where \code{Omega = R^{-1}}. The standardized residual is ~\code{N(0,1)} when
#' the z-scores and LD are mutually consistent, so a large residual marks an
#' allele-flip / LD-mismatch outlier. RSS-only helper, opt-in via
#' \code{alleleFlipKriging}; never wired into \code{alleleQc()} /
#' \code{matchRefPanel()}.
#'
#' @param zScore Numeric vector of harmonized z-scores.
#' @param R Square LD correlation matrix aligned to \code{zScore}.
#' @param variantIds Optional variant IDs for the diagnostics table.
#' @param pThreshold Two-sided p-value cutoff for flagging an outlier
#'   (default \code{5e-8}).
#' @param ridge Small diagonal added to \code{R} before inversion for numerical
#'   stability (default \code{1e-3}).
#' @return A list with \code{outlier} (logical vector) and \code{diagnostics}
#'   (data frame of per-variant predicted z, residual, statistic, p-value, and
#'   outlier flag).
#' @importFrom stats pnorm
#' @export
krigingOutlierQc <- function(zScore, R, variantIds = NULL,
                             pThreshold = 5e-8, ridge = 1e-3) {
  zScore <- as.numeric(zScore)
  m <- length(zScore)
  if (is.null(R) || !is.matrix(R) || nrow(R) != m || ncol(R) != m) {
    stop("krigingOutlierQc requires a square LD matrix aligned to zScore.")
  }
  if (is.null(variantIds)) variantIds <- rownames(R)
  # Regularize so the precision matrix is well-defined for collinear panels.
  Omega <- solve(R + diag(ridge, m))
  d <- diag(Omega)
  omegaZ <- as.numeric(Omega %*% zScore)
  condMean <- -(omegaZ - d * zScore) / d
  condVar <- 1 / d
  residual <- zScore - condMean
  statistic <- residual / sqrt(condVar)
  pValue <- 2 * pnorm(-abs(statistic))
  outlier <- !is.na(pValue) & pValue < pThreshold
  list(
    outlier = outlier,
    diagnostics = data.frame(
      variant_id = if (is.null(variantIds)) seq_len(m) else variantIds,
      z = zScore, predicted = condMean, residual = residual,
      statistic = statistic, p_value = pValue, outlier = outlier,
      stringsAsFactors = FALSE
    )
  )
}

# =============================================================================
# summaryStatsQc — SumStats-input QC pipeline (replaces the previous
# data.frame/LdData/QcResult-based summaryStatsQc and rssBasicQc).
# =============================================================================

# Convert one entry's GRanges into a flat data.frame with the column shape that
# MungeSumstats and .matchRefPanel expect (lower-case chrom/pos plus the
# CapsCase mcols).
.entryGrangesToDf <- function(gr) {
  mc <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
  out <- data.frame(
    chrom = sub("^chr", "", as.character(GenomicRanges::seqnames(gr)),
                ignore.case = TRUE),
    pos   = GenomicRanges::start(gr),
    stringsAsFactors = FALSE)
  cbind(out, mc)
}

# Build a refVariants data.frame (chrom, pos, A1, A2, variant_id) from the
# ldSketch GenotypeHandle's snpInfo so .matchRefPanel can join by (chrom, pos).
.refVariantsFromSketch <- function(handle) {
  si <- getSnpInfo(handle)
  chr <- sub("^chr", "", as.character(si$CHR), ignore.case = TRUE)
  data.frame(
    chrom      = chr,
    pos        = as.integer(si$BP),
    A1         = as.character(si$A1),
    A2         = as.character(si$A2),
    variant_id = as.character(si$SNP),
    stringsAsFactors = FALSE)
}

# Reassemble a harmonized data.frame into a GRanges with the SumStats mcol
# shape (SNP, A1, A2, Z, N, ... optional MAF/INFO/BETA/SE/P kept if present).
.dfToEntryGranges <- function(df) {
  chr <- paste0("chr", sub("^chr", "", as.character(df$chrom),
                           ignore.case = TRUE))
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges   = IRanges::IRanges(start = as.integer(df$pos), width = 1L))
  if (!is.null(df$variant_id) && is.null(df$SNP)) df$SNP <- df$variant_id
  baseCols <- c("SNP", "A1", "A2", "Z", "N")
  optCols  <- c("MAF", "INFO", "BETA", "SE", "P")
  use <- intersect(c(baseCols, optCols), colnames(df))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(df[, use, drop = FALSE])
  gr
}

# Drop variants whose (chrom, pos) overlaps any user-supplied skipRegion.
# skipRegion may be a character vector of "chr:start-end" strings or a GRanges.
.applySkipRegion <- function(df, skipRegion) {
  if (is.null(skipRegion) || length(skipRegion) == 0L) return(df)
  if (is.character(skipRegion)) {
    parsed <- do.call(rbind, lapply(skipRegion, function(s) {
      m <- regmatches(s, regexec("^([^:]+):([0-9]+)-([0-9]+)$", s))[[1]]
      if (length(m) != 4L)
        stop("skipRegion entry must be 'chr:start-end'; got '", s, "'")
      data.frame(chrom = sub("^chr", "", m[2], ignore.case = TRUE),
                 start = as.integer(m[3]),
                 end   = as.integer(m[4]),
                 stringsAsFactors = FALSE)
    }))
  } else if (methods::is(skipRegion, "GRanges")) {
    parsed <- data.frame(
      chrom = sub("^chr", "", as.character(GenomicRanges::seqnames(skipRegion)),
                  ignore.case = TRUE),
      start = GenomicRanges::start(skipRegion),
      end   = GenomicRanges::end(skipRegion),
      stringsAsFactors = FALSE)
  } else {
    stop("skipRegion must be a character vector of 'chr:start-end' ",
         "strings or a GRanges.")
  }
  dropMask <- rep(FALSE, nrow(df))
  dfChr <- sub("^chr", "", as.character(df$chrom), ignore.case = TRUE)
  for (i in seq_len(nrow(parsed))) {
    dropMask <- dropMask |
      (dfChr == parsed$chrom[i] &
       df$pos >= parsed$start[i] &
       df$pos <= parsed$end[i])
  }
  df[!dropMask, , drop = FALSE]
}

# Run the curated MungeSumstats::format_sumstats() pass. Returns the cleaned
# data.frame in pecotmr CapsCase column convention plus a per-entry audit
# record describing what was applied. Caller is responsible for catching
# errors and reporting them with the entry identity.
.runMungeSumstatsFilter <- function(df, refGenome, useDbsnpRefCheck,
                                    removeIndels, removeStrandAmbiguous,
                                    mafCutoff, infoCutoff, nCutoff,
                                    convertRefGenome, mungeSumstatsArgs) {
  if (!requireNamespace("MungeSumstats", quietly = TRUE))
    stop("Package 'MungeSumstats' is required for summaryStatsQc(). ",
         "Install it from Bioconductor.")
  # MungeSumstats writes through a temp tsv when invoked with a path; pass the
  # in-memory data.frame via the file-write/read path it provides.
  tmpIn  <- tempfile(fileext = ".tsv.gz")
  on.exit(unlink(tmpIn), add = TRUE)
  data.table::fwrite(df, tmpIn, sep = "\t")

  baseArgs <- list(
    path                = tmpIn,
    ref_genome          = if (is.null(refGenome)) "GRCh38" else refGenome,
    convert_ref_genome  = convertRefGenome,
    on_ref_genome       = isTRUE(useDbsnpRefCheck),
    infer_eff_direction = isTRUE(useDbsnpRefCheck),
    allele_flip_check   = isTRUE(useDbsnpRefCheck),
    bi_allelic_filter   = isTRUE(useDbsnpRefCheck),
    strand_ambig_filter = isTRUE(removeStrandAmbiguous),
    drop_indels         = isTRUE(removeIndels),
    FRQ_filter          = mafCutoff,
    INFO_filter         = infoCutoff,
    N_std               = nCutoff,
    return_data         = TRUE,
    return_format       = "data.table",
    log_folder_ind      = FALSE,
    log_mungesumstats_msgs = FALSE)
  # Caller-supplied pass-through overrides anything we set above.
  for (nm in names(mungeSumstatsArgs)) baseArgs[[nm]] <- mungeSumstatsArgs[[nm]]

  before <- nrow(df)
  cleaned <- do.call(MungeSumstats::format_sumstats, baseArgs)
  cleaned <- as.data.frame(cleaned)
  # Restore lower-case chrom/pos for downstream pecotmr code.
  if ("CHR" %in% colnames(cleaned)) {
    cleaned$chrom <- sub("^chr", "", as.character(cleaned$CHR),
                          ignore.case = TRUE)
    cleaned$CHR <- NULL
  }
  if ("BP" %in% colnames(cleaned)) {
    cleaned$pos <- as.integer(cleaned$BP)
    cleaned$BP <- NULL
  }
  list(df = cleaned, droppedNVariants = before - nrow(cleaned))
}

# Apply the panel-vs-sumstats allele harmonization using the slim
# .matchRefPanel against the ldSketch's variant info.
.matchAgainstSketch <- function(df, ldSketch, matchMinProp) {
  refVariants <- .refVariantsFromSketch(ldSketch)
  flipCandidates <- c("Z", "BETA")
  colToFlip <- intersect(flipCandidates, colnames(df))
  if (length(colToFlip) == 0L)
    stop("summaryStatsQc: input entry must contain at least one of Z or BETA ",
         "before panel harmonization.")
  colToComplement <- intersect("MAF", colnames(df))
  if (!"A1" %in% colnames(df) || !"A2" %in% colnames(df))
    stop("summaryStatsQc: input entry must contain A1 and A2 columns.")
  res <- .matchRefPanel(
    targetData      = df,
    refVariants     = refVariants,
    colToFlip       = colToFlip,
    colToComplement = colToComplement,
    matchMinProp    = matchMinProp,
    removeUnmatched = TRUE)
  out <- res$harmonizedData
  if (!"chrom" %in% colnames(out) && "chr" %in% colnames(out))
    colnames(out)[colnames(out) == "chr"] <- "chrom"
  out
}

# Apply ldMismatchQc (SLALOM/DENTIST) against the LD sketch.
.applyLdMismatchQcToEntry <- function(df, ldSketch, method) {
  variantIds <- df$SNP
  if (is.null(variantIds) || any(is.na(variantIds)))
    stop("summaryStatsQc: ldMismatchQc requires SNP column on the entry.")
  # Extract the panel block for these variants.
  snpIdx <- match(variantIds, as.character(getSnpInfo(ldSketch)$SNP))
  if (anyNA(snpIdx))
    stop("summaryStatsQc: ", sum(is.na(snpIdx)), " variant(s) in entry are ",
         "absent from the ldSketch panel; harmonize / impute before ",
         "calling zMismatchQc.")
  block <- extractBlockGenotypes(ldSketch, snpIdx, meanImpute = TRUE)
  dosage <- t(SummarizedExperiment::assay(block, "dosage"))
  colnames(dosage) <- variantIds
  R <- computeLd(dosage, method = "sample")
  qc <- ldMismatchQc(zScore = df$Z, R = R, nSample = getNSamples(ldSketch),
                     method = method)
  list(df = df[!qc$outlier, , drop = FALSE], outliers = sum(qc$outlier))
}

# Per-entry SER-based pip-screen (skip if no signal above the cutoff).
.applyPipScreen <- function(df, n, cutoff) {
  if (cutoff <= 0) return(list(df = df, skipped = FALSE))
  effectiveCutoff <- if (cutoff < 0) 3 / nrow(df) else cutoff
  pip <- susieR::susie_ser(z = df$Z, n = n, coverage = NULL)$pip
  if (!any(pip > effectiveCutoff)) {
    return(list(df = df[FALSE, , drop = FALSE], skipped = TRUE,
                reason = sprintf("no signals above PIP threshold %g",
                                 effectiveCutoff)))
  }
  list(df = df, skipped = FALSE)
}

# Internal: per-entry pipeline. Returns the cleaned GRanges and an audit list.
.runEntrySummaryStatsQc <- function(gr, ldSketch, refGenome, opts) {
  entryAudit <- list()
  df <- .entryGrangesToDf(gr)
  entryAudit$variantsIn <- nrow(df)

  # 1. MungeSumstats variant-content pass.
  mungeResult <- .runMungeSumstatsFilter(
    df,
    refGenome             = refGenome,
    useDbsnpRefCheck      = opts$useDbsnpRefCheck,
    removeIndels          = opts$removeIndels,
    removeStrandAmbiguous = opts$removeStrandAmbiguous,
    mafCutoff             = opts$mafCutoff,
    infoCutoff            = opts$infoCutoff,
    nCutoff               = opts$nCutoff,
    convertRefGenome      = opts$convertRefGenome,
    mungeSumstatsArgs     = opts$mungeSumstatsArgs)
  df <- mungeResult$df
  entryAudit$mungeSumstatsDropped <- mungeResult$droppedNVariants

  # 2. keepVariants subset.
  if (length(opts$keepVariants) > 0L) {
    before <- nrow(df)
    df <- df[df$SNP %in% opts$keepVariants, , drop = FALSE]
    entryAudit$keepVariantsDropped <- before - nrow(df)
  }

  # 3. skipRegion drop.
  if (!is.null(opts$skipRegion) && length(opts$skipRegion) > 0L) {
    before <- nrow(df)
    df <- .applySkipRegion(df, opts$skipRegion)
    entryAudit$skipRegionDropped <- before - nrow(df)
  }

  # 4. Optional PIP screen.
  if (opts$pipCutoffToSkip != 0) {
    pip <- .applyPipScreen(df, n = opts$nForPip, cutoff = opts$pipCutoffToSkip)
    df <- pip$df
    entryAudit$pipScreenSkipped <- isTRUE(pip$skipped)
    if (isTRUE(pip$skipped)) entryAudit$pipScreenReason <- pip$reason
  }

  if (nrow(df) < 2L) {
    entryAudit$earlyExit <- "fewer than two variants after pre-harmonization QC"
    return(list(gr = .dfToEntryGranges(df), audit = entryAudit))
  }

  # 5. Panel-vs-sumstats allele harmonization.
  df <- .matchAgainstSketch(df, ldSketch, matchMinProp = opts$matchMinProp)
  entryAudit$matchedAgainstSketch <- nrow(df)

  # 6. Optional kriging prefilter.
  if (isTRUE(opts$alleleFlipKriging) && nrow(df) >= 2L) {
    snpIdx <- match(df$SNP, as.character(getSnpInfo(ldSketch)$SNP))
    block <- extractBlockGenotypes(ldSketch, snpIdx, meanImpute = TRUE)
    dosage <- t(SummarizedExperiment::assay(block, "dosage"))
    colnames(dosage) <- df$SNP
    R <- computeLd(dosage, method = "sample")
    kr <- krigingOutlierQc(df$Z, R, variantIds = df$SNP)
    nKr <- sum(kr$outlier)
    if (nKr > 0L) df <- df[!kr$outlier, , drop = FALSE]
    entryAudit$krigingOutliersDropped <- nKr
  }

  # 7. Optional LD-mismatch QC.
  if (!identical(opts$zMismatchQc, "none") && nrow(df) >= 2L) {
    ldQc <- .applyLdMismatchQcToEntry(df, ldSketch, opts$zMismatchQc)
    df <- ldQc$df
    entryAudit$ldMismatchOutliersDropped <- ldQc$outliers
    entryAudit$ldMismatchMethod          <- opts$zMismatchQc
  }

  # 8. Optional RAISS imputation against the ldSketch.
  if (isTRUE(opts$impute) && nrow(df) >= 1L) {
    refPanel <- .refVariantsFromSketch(ldSketch)
    refPanel <- refPanel[order(refPanel$pos), , drop = FALSE]

    knownVariantIds <- if (!is.null(df$SNP)) as.character(df$SNP)
                       else as.character(df$variant_id)
    knownZ <- data.frame(
      chrom      = as.character(df$chrom),
      pos        = as.integer(df$pos),
      variant_id = knownVariantIds,
      A1         = as.character(df$A1),
      A2         = as.character(df$A2),
      z          = as.numeric(df$Z),
      stringsAsFactors = FALSE)
    if ("N"    %in% colnames(df)) knownZ$n    <- as.numeric(df$N)
    if ("BETA" %in% colnames(df)) knownZ$beta <- as.numeric(df$BETA)
    if ("SE"   %in% colnames(df)) knownZ$se   <- as.numeric(df$SE)
    knownZ <- knownZ[order(knownZ$pos), , drop = FALSE]

    # Materialize the full panel dosage in panel-order matching refPanel.
    sketchSnpInfo <- getSnpInfo(ldSketch)
    block <- extractBlockGenotypes(
      ldSketch, seq_len(nrow(sketchSnpInfo)), meanImpute = TRUE)
    dosage <- t(SummarizedExperiment::assay(block, "dosage"))
    colnames(dosage) <- as.character(sketchSnpInfo$SNP)
    dosage <- dosage[, refPanel$variant_id, drop = FALSE]
    scaledDosage <- scale(dosage)
    scaledDosage[is.na(scaledDosage)] <- 0

    imputed <- raiss(
      refPanel       = refPanel,
      knownZscores   = knownZ,
      genotypeMatrix = scaledDosage,
      svdTol         = if (is.null(opts$imputeOpts$svdTol)) 1e-12
                       else opts$imputeOpts$svdTol,
      lamb           = if (is.null(opts$imputeOpts$lamb)) 0.01
                       else opts$imputeOpts$lamb,
      r2Threshold    = if (is.null(opts$imputeOpts$r2Threshold)) 0.6
                       else opts$imputeOpts$r2Threshold,
      minimumLd      = if (is.null(opts$imputeOpts$minimumLd)) 5
                       else opts$imputeOpts$minimumLd,
      verbose        = FALSE)
    if (!is.null(imputed) && !is.null(imputed$resultFilter)) {
      impDf <- imputed$resultFilter
      out <- data.frame(
        chrom = impDf$chrom,
        pos   = impDf$pos,
        SNP   = impDf$variant_id,
        A1    = impDf$A1,
        A2    = impDf$A2,
        Z     = impDf$z,
        stringsAsFactors = FALSE)
      if ("n"    %in% colnames(impDf)) out$N    <- impDf$n
      if ("beta" %in% colnames(impDf)) out$BETA <- impDf$beta
      if ("se"   %in% colnames(impDf)) out$SE   <- impDf$se
      if ("N" %in% colnames(out) && any(is.na(out$N)))
        out$N[is.na(out$N)] <- stats::median(out$N, na.rm = TRUE)
      entryAudit$raissTotalVariants    <- nrow(out)
      entryAudit$raissImputedVariants  <- nrow(out) - nrow(knownZ)
      df <- out
    } else {
      entryAudit$raissImputedVariants <- 0L
    }
  }

  entryAudit$variantsOut <- nrow(df)
  list(gr = .dfToEntryGranges(df), audit = entryAudit)
}

#' Run QC on a SumStats Collection
#'
#' Applies a single QC pass to a \code{QtlSumStats} or \code{GwasSumStats}
#' collection: delegates variant-content QC (column standardization,
#' indels, strand-ambiguous, MAF/INFO/N filters, p-value & effect-size
#' sanity checks, optional dbSNP / liftover) to
#' \code{MungeSumstats::format_sumstats()}, then runs pecotmr-specific
#' steps (\code{skipRegion}, optional PIP screen, panel harmonization
#' against the \code{ldSketch} via \code{.matchRefPanel}, optional
#' SLALOM/DENTIST LD-mismatch QC, optional RAISS imputation).
#'
#' The returned collection has its \code{qcInfo} slot populated with a
#' per-entry audit record (variant counts, drop counts at each step,
#' which filters fired, etc.). Fine-mapping and TWAS-weights pipelines
#' reject SumStats inputs where \code{length(getQcInfo(x)) == 0L}.
#'
#' Column-availability error contract: a non-zero \code{mafCutoff}
#' requires every entry to carry a \code{MAF} column; non-zero
#' \code{infoCutoff} requires \code{INFO}; non-zero \code{nCutoff}
#' requires \code{N}. Missing column with a non-zero cutoff is a hard
#' error.
#'
#' @param sumstats A \code{QtlSumStats} or \code{GwasSumStats}
#'   collection.
#' @param useDbsnpRefCheck One-shot opt-in to MungeSumstats's
#'   dbSNP-based reference-genome / allele-flip / biallelic-filter
#'   checks. When \code{TRUE}, sets \code{on_ref_genome},
#'   \code{infer_eff_direction}, \code{allele_flip_check}, and
#'   \code{bi_allelic_filter} all to \code{TRUE} simultaneously.
#'   Default \code{FALSE} (trust input alleles, lighter dependency
#'   footprint).
#' @param removeIndels Logical (length 1). When \code{TRUE}, drop
#'   indels. Default \code{FALSE} (match MungeSumstats default).
#' @param removeStrandAmbiguous Logical (length 1). When \code{TRUE},
#'   drop A/T and C/G strand-ambiguous variants. Default \code{TRUE}.
#' @param mafCutoff Numeric (length 1). MAF threshold (variants with
#'   \code{MAF < mafCutoff} are dropped). Default 0. Requires \code{MAF}
#'   column.
#' @param infoCutoff Numeric (length 1). INFO score threshold. Default
#'   0. Requires \code{INFO} column.
#' @param nCutoff Numeric (length 1). MungeSumstats \code{N_std} value;
#'   sample-size deviation threshold. Default 5.
#' @param keepVariants Optional character vector of variant IDs (SNP
#'   column) to retain prior to harmonization.
#' @param skipRegion Optional character vector of \code{"chr:start-end"}
#'   strings, or a \code{GRanges}, of regions to drop.
#' @param pipCutoffToSkip Numeric (length 1). When \code{!= 0}, run an
#'   LD-independent single-effect SER screen and skip the entry if no
#'   PIP exceeds the cutoff. \code{< 0} resolves to \code{3 / nVariants}.
#'   Default 0 (no screen).
#' @param zMismatchQc One of \code{"none"} (default), \code{"slalom"},
#'   \code{"dentist"}.
#' @param alleleFlipKriging Logical (length 1). Opt-in kriging
#'   LD-consistency prefilter run before SLALOM/DENTIST. Default
#'   \code{FALSE}.
#' @param impute Logical (length 1). Run RAISS imputation against the
#'   \code{ldSketch}. Default \code{FALSE}. (Note: RAISS against the
#'   sketch is not yet fully wired for the new path; the option is
#'   accepted but currently emits a warning and is skipped.)
#' @param imputeOpts Named list of RAISS parameters.
#' @param convertRefGenome Optional character (\code{"GRCh37"} or
#'   \code{"GRCh38"}) to liftover the sumstats to via MungeSumstats.
#'   \code{NULL} (default) skips liftover.
#' @param matchMinProp Minimum proportion of LD panel variants that must
#'   be matched by the sumstats; default 0.
#' @param mungeSumstatsArgs Optional named list of pass-through args to
#'   \code{MungeSumstats::format_sumstats()}. Any name supplied here
#'   overrides the value the curated knobs would have set.
#' @return A new \code{QtlSumStats} / \code{GwasSumStats} with cleaned
#'   entries and \code{qcInfo} populated.
#' @export
summaryStatsQc <- function(sumstats,
                           useDbsnpRefCheck       = FALSE,
                           removeIndels           = FALSE,
                           removeStrandAmbiguous  = TRUE,
                           mafCutoff              = 0,
                           infoCutoff             = 0,
                           nCutoff                = 5,
                           keepVariants           = NULL,
                           skipRegion             = NULL,
                           pipCutoffToSkip        = 0,
                           zMismatchQc            = c("none", "slalom",
                                                     "dentist"),
                           alleleFlipKriging      = FALSE,
                           impute                 = FALSE,
                           imputeOpts             = list(rcond = 0.01,
                                                        r2Threshold = 0.6,
                                                        minimumLd = 5,
                                                        lamb = 0.01),
                           convertRefGenome       = NULL,
                           matchMinProp           = 0,
                           mungeSumstatsArgs      = list()) {
  if (!methods::is(sumstats, "QtlSumStats") &&
      !methods::is(sumstats, "GwasSumStats")) {
    stop("summaryStatsQc requires a QtlSumStats or GwasSumStats input.")
  }
  zMismatchQc <- match.arg(zMismatchQc)

  # Column-availability checks across all entries.
  for (i in seq_len(nrow(sumstats))) {
    mc <- S4Vectors::mcols(sumstats$entry[[i]])
    cols <- colnames(mc)
    if (mafCutoff > 0 && !"MAF" %in% cols)
      stop("summaryStatsQc: mafCutoff > 0 requires every entry to carry a ",
           "MAF column; entry ", i, " does not.")
    if (infoCutoff > 0 && !"INFO" %in% cols)
      stop("summaryStatsQc: infoCutoff > 0 requires every entry to carry an ",
           "INFO column; entry ", i, " does not.")
  }

  opts <- list(
    useDbsnpRefCheck       = useDbsnpRefCheck,
    removeIndels           = removeIndels,
    removeStrandAmbiguous  = removeStrandAmbiguous,
    mafCutoff              = mafCutoff,
    infoCutoff             = infoCutoff,
    nCutoff                = nCutoff,
    keepVariants           = as.character(keepVariants),
    skipRegion             = skipRegion,
    pipCutoffToSkip        = pipCutoffToSkip,
    zMismatchQc            = zMismatchQc,
    alleleFlipKriging      = alleleFlipKriging,
    impute                 = impute,
    imputeOpts             = imputeOpts,
    convertRefGenome       = convertRefGenome,
    matchMinProp           = matchMinProp,
    mungeSumstatsArgs      = mungeSumstatsArgs,
    nForPip                = NULL)

  newEntries <- vector("list", nrow(sumstats))
  entryAudits <- vector("list", nrow(sumstats))
  for (i in seq_len(nrow(sumstats))) {
    opts$nForPip <- if ("N" %in% colnames(S4Vectors::mcols(sumstats$entry[[i]])))
      stats::median(S4Vectors::mcols(sumstats$entry[[i]])$N, na.rm = TRUE)
    else NULL
    result <- .runEntrySummaryStatsQc(
      gr        = sumstats$entry[[i]],
      ldSketch  = getLdSketch(sumstats),
      refGenome = getGenome(sumstats),
      opts      = opts)
    newEntries[[i]] <- result$gr
    entryAudits[[i]] <- result$audit
  }

  qcInfo <- list(
    timestamp        = NA_character_,
    options          = list(
      useDbsnpRefCheck      = useDbsnpRefCheck,
      removeIndels          = removeIndels,
      removeStrandAmbiguous = removeStrandAmbiguous,
      mafCutoff             = mafCutoff,
      infoCutoff            = infoCutoff,
      nCutoff               = nCutoff,
      zMismatchQc           = zMismatchQc,
      alleleFlipKriging     = alleleFlipKriging,
      impute                = impute,
      convertRefGenome      = convertRefGenome),
    entryAudit       = entryAudits,
    mungeSumstatsArgs = mungeSumstatsArgs)

  # Rebuild the SumStats with new entries and qcInfo.
  if (methods::is(sumstats, "GwasSumStats")) {
    GwasSumStats(
      study    = as.character(sumstats$study),
      entry    = newEntries,
      genome   = getGenome(sumstats),
      ldSketch = getLdSketch(sumstats),
      varY     = as.numeric(sumstats$varY),
      qcInfo   = qcInfo)
  } else {
    QtlSumStats(
      study    = as.character(sumstats$study),
      context  = as.character(sumstats$context),
      trait    = as.character(sumstats$trait),
      entry    = newEntries,
      genome   = getGenome(sumstats),
      ldSketch = getLdSketch(sumstats),
      varY     = as.numeric(sumstats$varY),
      qcInfo   = qcInfo)
  }
}
