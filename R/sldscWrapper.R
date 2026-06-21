# Stratified LD Score Regression (S-LDSC) post-processing wrappers around polyfun.
#
# This file provides the post-processing layer for the xqtl-protocol sLDSC pipeline:
# read polyfun outputs per trait, compute Gazal-style standardized tau* and the
# differential per-SNP heritability statistic (EnrichStat), and run DerSimonian-Laird
# random-effects meta-analysis across traits.
#
# Reference panel convention: all LD-derived quantities (baseline LD scores,
# target LD scores, regression weights, allele frequencies) must come from the
# same reference panel. Do not mix files from different panels (e.g. 1000G vs ADSP).
#
# MAF convention: by default we restrict to MAF > 5% per the sLDSC recommendation.
# Pass maf_cutoff = 0 to opt out (not recommended).
#
# Cross-type comparison: tau* (Gazal et al. 2017 standardization) is the
# cross-type comparable statistic. Use tau* to rank or meta-analyze annotations
# that mix binary and continuous types. E (proportion-based enrichment) is
# scale-dependent for continuous annotations and is only comparable within type.

# ---- internal helpers ----

.sldscStdCols <- c("CHR", "SNP", "BP", "CM", "A1", "A2", "MAF")

.sldscChromFromFilename <- function(f) {
  bn <- basename(f)
  m  <- regmatches(bn, regexec("\\.([0-9]+)\\.annot\\.gz$", bn))[[1]]
  if (length(m) >= 2) as.integer(m[2]) else NA_integer_
}

.sldscDetectAnnotCols <- function(filePath) {
  sample <- vroom(filePath, n_max = 5L, show_col_types = FALSE)
  setdiff(names(sample), .sldscStdCols)
}


#' @title Read S-LDSC outputs from polyfun for one trait/run
#'
#' @description Reads the regression outputs produced by `polyfun/ldsc.py` for a
#'   single polyfun run (one trait, one annotation set) and returns them as a
#'   tidy list ready for downstream standardization. Hides the underlying file
#'   formats; downstream code consumes only modeling quantities.
#'
#' @param prefix Character. Path prefix to the polyfun outputs for one trait/run.
#'   The function appends `.results`, `.log`, and `.part_delete` to this prefix.
#'   Example: `"/path/to/cwd/CAD_META.filtered.sumstats.gz"`.
#'
#' @return A named list. See `sldscPostprocessingPipeline` for components.
#'
#' @examples
#' \dontrun{
#' run <- readSldscTrait("/output/CAD_META.filtered.sumstats.gz")
#' run$tau["my_target_annotation"]
#' }
#'
#' @importFrom stats setNames var na.omit
#' @importFrom utils head
#' @importFrom vroom vroom
#' @export
readSldscTrait <- function(prefix) {
  resultsFile <- paste0(prefix, ".results")
  logFile     <- paste0(prefix, ".log")
  deleteFile  <- paste0(prefix, ".part_delete")

  for (f in c(resultsFile, logFile, deleteFile)) {
    if (!file.exists(f)) stop("readSldscTrait: missing file: ", f)
  }

  results <- vroom(resultsFile, show_col_types = FALSE)
  cats <- as.character(results$Category)

  logLines <- readLines(logFile, warn = FALSE)
  h2Line <- grep("Total Observed scale h2:", logLines, value = TRUE)
  if (length(h2Line) == 0L)
    stop("readSldscTrait: could not find 'Total Observed scale h2:' in ", logFile)
  h2g <- suppressWarnings(as.numeric(gsub(".*h2: (-?[0-9.eE+-]+).*", "\\1", h2Line[1])))
  if (is.na(h2g))
    stop("readSldscTrait: failed to parse h2g numeric from log line: ", h2Line[1])

  deleteValues <- as.matrix(vroom(deleteFile, show_col_types = FALSE))
  if (ncol(deleteValues) != length(cats)) {
    stop("readSldscTrait: .part_delete has ", ncol(deleteValues),
         " columns but .results has ", length(cats), " categories.")
  }
  colnames(deleteValues) <- cats

  list(
    categories     = cats,
    tau            = setNames(as.numeric(results$Coefficient),                 cats),
    tauSe         = setNames(as.numeric(results[["Coefficient_std_error"]]),  cats),
    enrichment     = setNames(as.numeric(results$Enrichment),                  cats),
    enrichmentSe  = setNames(as.numeric(results[["Enrichment_std_error"]]),   cats),
    enrichmentP   = setNames(as.numeric(results[["Enrichment_p"]]),           cats),
    propH2        = setNames(as.numeric(results[["Prop._h2"]]),               cats),
    propSnps      = setNames(as.numeric(results[["Prop._SNPs"]]),             cats),
    h2g            = h2g,
    tauBlocks     = deleteValues,
    nBlocks       = nrow(deleteValues)
  )
}


#' @title Compute per-annotation standard deviation, MAF-restricted
#'
#' @description Computes the standard deviation of each annotation column in the
#'   target annotation files, restricted to SNPs above a MAF cutoff via PLINK
#'   `.frq` files. Required for internal consistency with polyfun's regression,
#'   which operates on MAF > cutoff SNPs by default.
#'
#' @param targetAnnoDir Character. Directory containing target annotation files
#'   (one per chromosome) in polyfun's `.annot.gz` format.
#' @param frqfileDir Character or NULL. Directory containing PLINK `.frq` files
#'   for the reference panel. Required when `mafCutoff > 0`; the function
#'   errors if missing.
#' @param plinkName Character. Filename prefix of the `.frq` files
#'   (e.g. `"ADSP_chr"`). Files are expected at `{plinkName}{chr}.frq`.
#' @param mafCutoff Numeric, default `0.05`.
#' @param annotCols Character or integer vector, default NULL. Annotation columns
#'   to compute sd for. If NULL, all annotation columns are used.
#'
#' @return Named numeric vector of \eqn{sd_C} values, one per annotation.
#'
#' @importFrom stats setNames var
#' @export
computeSldscAnnotSd <- function(targetAnnoDir, frqfileDir = NULL,
                                plinkName = "ADSP_chr",
                                mafCutoff = 0.05, annotCols = NULL) {
  if (mafCutoff > 0 && (is.null(frqfileDir) || !dir.exists(frqfileDir))) {
    stop("computeSldscAnnotSd: mafCutoff = ", mafCutoff,
         " requires frqfileDir, but '", frqfileDir, "' is not a directory.")
  }
  if (!dir.exists(targetAnnoDir)) {
    stop("computeSldscAnnotSd: targetAnnoDir does not exist: ", targetAnnoDir)
  }

  annoFiles <- list.files(targetAnnoDir, pattern = "\\.annot\\.gz$", full.names = TRUE)
  if (length(annoFiles) == 0L)
    stop("computeSldscAnnotSd: no .annot.gz files in: ", targetAnnoDir)

  detected <- .sldscDetectAnnotCols(annoFiles[1])
  if (is.null(annotCols)) {
    colsUse <- detected
  } else if (is.numeric(annotCols)) {
    colsUse <- detected[annotCols]
  } else {
    colsUse <- annotCols
  }
  if (length(colsUse) == 0L)
    stop("computeSldscAnnotSd: no annotation columns to process.")

  num <- setNames(numeric(length(colsUse)), colsUse)
  den <- 0

  for (annoFile in annoFiles) {
    dat <- vroom(annoFile, show_col_types = FALSE)
    if (mafCutoff > 0) {
      chrom <- .sldscChromFromFilename(annoFile)
      if (is.na(chrom))
        stop("computeSldscAnnotSd: could not parse chromosome from: ", annoFile)
      frqFile <- file.path(frqfileDir, paste0(plinkName, chrom, ".frq"))
      if (!file.exists(frqFile))
        stop("computeSldscAnnotSd: .frq file not found: ", frqFile)
      frq <- vroom(frqFile, col_select = c("SNP", "MAF"), show_col_types = FALSE)
      dat <- merge(dat, frq, by = "SNP", all.x = FALSE, all.y = FALSE)
      dat <- dat[!is.na(dat$MAF) & dat$MAF > mafCutoff, ]
    }
    if (nrow(dat) <= 1L) next
    nMinus1 <- nrow(dat) - 1L
    for (col in colsUse) {
      vals <- as.numeric(dat[[col]])
      v <- var(vals, na.rm = TRUE)
      if (!is.na(v)) num[col] <- num[col] + nMinus1 * v
    }
    den <- den + nMinus1
  }

  if (den <= 0)
    stop("computeSldscAnnotSd: zero degrees of freedom after MAF filtering.")
  sqrt(num / den)
}


#' @title Reference-panel SNP count (the M_ref used to standardise tau*)
#'
#' @description `M_ref` is the number of SNPs in the REFERENCE PANEL over which
#'   heritability is partitioned in the sLDSC model
#'   (`h2(C) = sum_{j in M_ref} a_C(j) sum_{C'} tau_{C'} a_{C'}(j)`). It is
#'   panel-defined and is **not** the regression SNP set (HapMap3 ~1M) nor any
#'   HM3-subsetted target output:
#'   \itemize{
#'     \item `mafCutoff > 0` (Gazal/Finucane convention): count MAF > cutoff
#'       SNPs across all `.frq` files (the same set polyfun's `.l2.M_5_50` sums).
#'     \item `mafCutoff == 0` (all-M variant): count ALL SNPs across all
#'       `.frq` files (the same set polyfun's `.l2.M` sums).
#'   }
#'   `targetAnnoDir` is a fallback only, used when no `.frq` directory is
#'   given; that fallback counts `.l2.ldscore` rows and is WRONG when the target
#'   was HM3-subsetted (it then yields the regression SNP count, not M_ref).
#'
#' @param targetAnnoDir Character or NULL. Fallback only - directory of
#'   `.l2.ldscore` files. Used only when `frqfileDir` is unavailable.
#' @param frqfileDir Character or NULL. Directory of PLINK `.frq` files; the
#'   preferred (recommended) source of M_ref.
#' @param plinkName Character. Filename prefix of `.frq` files.
#' @param mafCutoff Numeric, default `0.05`.
#'
#' @return Scalar integer.
#'
#' @export
computeSldscMRef <- function(targetAnnoDir = NULL, frqfileDir = NULL,
                             plinkName = "ADSP_chr", mafCutoff = 0.05) {
  ## --- preferred path: count reference-panel SNPs from the .frq files ---
  if (!is.null(frqfileDir) && dir.exists(frqfileDir)) {
    pat <- paste0("^", gsub("([.])", "\\\\\\1", plinkName), "[0-9]+\\.frq$")
    frqFiles <- list.files(frqfileDir, pattern = pat, full.names = TRUE)
    if (length(frqFiles) == 0L)
      frqFiles <- list.files(frqfileDir, pattern = "\\.frq$", full.names = TRUE)
    if (length(frqFiles) > 0L) {
      total <- 0L
      for (f in frqFiles) {
        frq <- vroom(f, col_select = "MAF", show_col_types = FALSE)
        total <- total + if (mafCutoff > 0)
          sum(!is.na(frq$MAF) & frq$MAF > mafCutoff) else nrow(frq)
      }
      return(as.integer(total))
    }
  }

  ## --- fallback only (no .frq dir): count target .l2.ldscore rows ---
  if (mafCutoff > 0)
    stop("computeSldscMRef: mafCutoff = ", mafCutoff,
         " requires frqfileDir (to count MAF>cutoff reference-panel SNPs).")
  if (is.null(targetAnnoDir) || !dir.exists(targetAnnoDir))
    stop("computeSldscMRef: need frqfileDir, or targetAnnoDir as fallback.")
  files <- list.files(targetAnnoDir,
                      pattern = "\\.l2\\.ldscore\\.(gz|parquet)$",
                      full.names = TRUE)
  if (length(files) == 0L)
    stop("computeSldscMRef: no .frq files and no .l2.ldscore files found.")
  warning("computeSldscMRef: no .frq dir given; counting target .l2.ldscore ",
          "rows as M_ref. If the target was HM3-subsetted this UNDERCOUNTS the ",
          "reference panel and shrinks tau*. Pass frqfileDir instead.")
  total <- 0L
  for (f in files) {
    if (endsWith(f, ".parquet")) {
      if (!requireNamespace("arrow", quietly = TRUE))
        stop("computeSldscMRef: install 'arrow' to read .parquet files.")
      total <- total + nrow(arrow::read_parquet(f))
    } else {
      total <- total + nrow(vroom(f, show_col_types = FALSE))
    }
  }
  as.integer(total)
}


#' @title Detect whether each annotation is binary or continuous
#'
#' @description Inspects each annotation column and returns whether its values
#'   lie in \{0, 1\} (binary) or take other values (continuous).
#'
#' @param targetAnnoDir Character. Directory containing the target `.annot.gz`
#'   files (one per chromosome).
#' @param annotCols Character or integer vector, default NULL.
#'
#' @return Named logical vector: TRUE for binary, FALSE for continuous.
#'
#' @importFrom stats setNames
#' @export
isBinarySldscAnnot <- function(targetAnnoDir, annotCols = NULL) {
  annoFiles <- list.files(targetAnnoDir, pattern = "\\.annot\\.gz$", full.names = TRUE)
  if (length(annoFiles) == 0L)
    stop("isBinarySldscAnnot: no .annot.gz files in: ", targetAnnoDir)

  detected <- .sldscDetectAnnotCols(annoFiles[1])
  if (is.null(annotCols)) {
    colsUse <- detected
  } else if (is.numeric(annotCols)) {
    colsUse <- detected[annotCols]
  } else {
    colsUse <- annotCols
  }

  isBinary <- setNames(rep(TRUE, length(colsUse)), colsUse)

  for (f in annoFiles) {
    dat <- vroom(f, col_select = all_of(colsUse), show_col_types = FALSE)
    for (col in colsUse) {
      if (!isBinary[[col]]) next
      vals <- unique(na.omit(as.numeric(dat[[col]])))
      if (any(!(vals %in% c(0, 1)))) isBinary[[col]] <- FALSE
    }
    if (!any(isBinary)) break
  }

  isBinary
}


#' @title Standardize tau and compute EnrichStat for one polyfun run
#'
#' @description Applies the Gazal standardization
#'   \eqn{\tau^*_C = \tau_C \cdot sd_C \cdot M_{ref} / h^2_g} to the point and
#'   to each jackknife block. For `mode = "single"`, additionally computes
#'   EnrichStat and back-solves its standard error from polyfun's reported
#'   `Enrichment_p` using \eqn{|Z| = \Phi^{-1}(1 - p/2)}.
#'
#' @param traitData List from \code{\link{readSldscTrait}}.
#' @param sdAnnot Named numeric vector from \code{\link{computeSldscAnnotSd}}.
#' @param MRef Scalar from \code{\link{computeSldscMRef}}.
#' @param targetCategories Character vector or NULL. If NULL, intersects
#'   `traitData$categories` with `names(sdAnnot)`.
#' @param mode Character: `"single"` or `"joint"`.
#'
#' @return A list with `summary` (data frame), `tau_star_blocks` (matrix),
#'   `h2g`, `nBlocks`, `mode`.
#'
#' @importFrom stats qnorm var
#' @export
standardizeSldscTrait <- function(traitData, sdAnnot, MRef,
                                  targetCategories = NULL,
                                  mode = c("single", "joint")) {
  mode <- match.arg(mode)
  if (is.null(targetCategories))
    targetCategories <- intersect(traitData$categories, names(sdAnnot))
  if (length(targetCategories) == 0L)
    stop("standardizeSldscTrait: no target categories.")

  targetIdx <- match(targetCategories, traitData$categories)
  if (any(is.na(targetIdx)))
    stop("standardizeSldscTrait: missing categories: ",
         paste(targetCategories[is.na(targetIdx)], collapse = ", "))

  h2g <- traitData$h2g
  sdTarget <- as.numeric(sdAnnot[targetCategories])
  if (any(is.na(sdTarget) | sdTarget == 0))
    warning("standardizeSldscTrait: zero/NA sd for some targets; tau* will be NA/0.")

  tau        <- as.numeric(traitData$tau[targetCategories])
  tauSe      <- as.numeric(traitData$tauSe[targetCategories])
  blocksTarget   <- traitData$tauBlocks[, targetIdx, drop = FALSE]

  ts <- standardizeTauStar(tau, blocksTarget, sdTarget, MRef, h2g)
  tauStar    <- ts$tauStar
  tauStarSe  <- ts$tauStarSe

  tauStarBlocks <- sweep(blocksTarget, 2L, sdTarget * MRef / h2g, FUN = "*")

  summaryDf <- data.frame(
    target      = targetCategories,
    tau         = tau,
    tauSe      = tauSe,
    tauStar    = tauStar,
    tauStarSe = tauStarSe,
    stringsAsFactors = FALSE
  )

  if (mode == "single") {
    enrich    <- as.numeric(traitData$enrichment[targetCategories])
    enrichSe  <- as.numeric(traitData$enrichmentSe[targetCategories])
    enrichP   <- as.numeric(traitData$enrichmentP[targetCategories])
    pH2       <- as.numeric(traitData$propH2[targetCategories])
    pM        <- as.numeric(traitData$propSnps[targetCategories])

    diffRatio   <- (pH2 / pM) - (1 - pH2) / (1 - pM)
    enrichstat  <- (h2g / MRef) * diffRatio

    absZ <- qnorm(1 - enrichP / 2)
    enrichstatSe <- abs(enrichstat) / absZ
    enrichstatSe[!is.finite(absZ) | absZ <= 0] <- NA_real_

    summaryDf$enrichment    <- enrich
    summaryDf$enrichmentSe <- enrichSe
    summaryDf$enrichmentP  <- enrichP
    summaryDf$enrichstat    <- enrichstat
    summaryDf$enrichstatSe <- enrichstatSe
  }

  list(
    summary         = summaryDf,
    tau_star_blocks = tauStarBlocks,
    h2g             = h2g,
    nBlocks        = nrow(blocksTarget),
    mode            = mode
  )
}


#' @title Random-effects meta-analysis of S-LDSC quantities across traits
#'
#' @description DerSimonian-Laird random-effects meta-analysis of one S-LDSC
#'   quantity for one annotation across multiple traits.
#'
#' @details Per-trait \eqn{SE_i} sources:
#'   - `quantity = "tauStar"`: jackknife SE from per-block \eqn{\tau^*}.
#'   - `quantity = "enrichment"`: polyfun-reported `Enrichment_std_error`.
#'   - `quantity = "enrichstat"`: back-solved SE from polyfun's `Enrichment_p`.
#'
#' @param perTraitEstimates Named list of per-trait results (each with a
#'   `summary` data frame).
#' @param category Character. Annotation name to meta-analyze.
#' @param quantity Character: `"tauStar"`, `"enrichment"`, or `"enrichstat"`.
#'
#' @return List with `mean`, `se`, `p`, `nTraits`, `traitsUsed`, `tau2`.
#'
#' @importFrom stats pnorm
#' @export
metaSldscRandom <- function(perTraitEstimates, category,
                            quantity = c("tauStar", "enrichment", "enrichstat")) {
  quantity <- match.arg(quantity)
  colPairs <- list(
    tauStar   = c("tauStar",   "tauStarSe"),
    enrichment = c("enrichment", "enrichmentSe"),
    enrichstat = c("enrichstat", "enrichstatSe")
  )
  cols <- colPairs[[quantity]]
  traitNames <- names(perTraitEstimates)
  if (is.null(traitNames))
    traitNames <- as.character(seq_along(perTraitEstimates))

  means <- numeric(0); ses <- numeric(0); used <- character(0)
  for (i in seq_along(perTraitEstimates)) {
    pt <- perTraitEstimates[[i]]
    if (is.null(pt) || is.null(pt$summary)) next
    df <- pt$summary
    row <- df[df$target == category, , drop = FALSE]
    if (nrow(row) == 0L) next
    if (!all(cols %in% names(row))) next
    m <- as.numeric(row[[cols[1]]])[1]
    s <- as.numeric(row[[cols[2]]])[1]
    if (is.na(m) || is.na(s) || !is.finite(s) || s <= 0) next
    means <- c(means, m); ses <- c(ses, s); used <- c(used, traitNames[i])
  }

  if (length(means) < 2L) {
    return(list(mean = NA_real_, se = NA_real_, p = NA_real_,
                nTraits = length(means), traitsUsed = used,
                tau2 = NA_real_))
  }
  meta <- metaRandomEffects(means, ses)
  z    <- meta$mean / meta$se
  p    <- 2 * pnorm(-abs(z))
  list(
    mean        = meta$mean,
    se          = meta$se,
    p           = as.numeric(p),
    nTraits    = length(means),
    traitsUsed = used,
    tau2        = meta$tau2
  )
}


# Internal helper: assemble a wide per-trait summary frame with single + joint
# columns side by side.
.sldscAssembleTraitSummary <- function(singleDf, jointDf, targetCategories,
                                       isBinaryVec) {
  rows <- if (!is.null(singleDf)) singleDf$target else
          if (!is.null(jointDf))  jointDf$target  else targetCategories
  out <- data.frame(target = rows,
                    isBinary = unname(isBinaryVec[rows]),
                    stringsAsFactors = FALSE)

  addCols <- function(out, src, suffix) {
    colsToAdd <- c("tau", "tauSe", "tauStar", "tauStarSe",
                   "enrichment", "enrichmentSe", "enrichmentP",
                   "enrichstat", "enrichstatSe")
    suffixCap <- paste0(toupper(substring(suffix, 1, 1)),
                        substring(suffix, 2))
    for (c in colsToAdd) {
      newcol <- paste0(c, suffixCap)
      if (!is.null(src) && c %in% names(src)) {
        out[[newcol]] <- src[[c]][match(out$target, src$target)]
      } else {
        out[[newcol]] <- NA_real_
      }
    }
    out
  }
  out <- addCols(out, singleDf, "single")
  out <- addCols(out, jointDf,  "joint")
  out
}


# Internal helper: build a per-trait list view that metaSldscRandom can read.
# Each list element has a $summary frame with the requested mode's columns
# renamed to the canonical names (tauStar, tauStarSe, enrichment, ...).
.sldscViewForMeta <- function(perTrait, suffix) {
  lapply(perTrait, function(pt) {
    if (is.null(pt$summary)) return(NULL)
    df <- pt$summary
    colsHave <- c("tauStar", "tauStarSe", "enrichment", "enrichmentSe",
                  "enrichmentP", "enrichstat", "enrichstatSe")
    suffixCap <- paste0(toupper(substring(suffix, 1, 1)),
                        substring(suffix, 2))
    srcCols <- paste0(colsHave, suffixCap)
    avail    <- srcCols %in% names(df)
    if (!any(avail)) return(NULL)
    newDf <- data.frame(target = df$target, stringsAsFactors = FALSE)
    for (k in seq_along(colsHave)) {
      if (avail[k]) newDf[[colsHave[k]]] <- df[[srcCols[k]]]
    }
    list(summary = newDf)
  })
}


#' @title End-to-end S-LDSC post-processing across traits, single + joint in one pass
#'
#' @description Top-level orchestration. Reads polyfun outputs (one single-target
#'   run per target plus, when available, one joint run per trait), standardizes
#'   both modes, and runs the default random-effects meta across all traits.
#'
#' @param traitSinglePrefixes Named list. For each trait, a character vector
#'   of length \eqn{N} giving the polyfun output prefixes for the \eqn{N}
#'   single-target runs (order must match `targetCategories`).
#' @param traitJointPrefix Named character. For each trait, the polyfun output
#'   prefix for the joint run. Pass `NA` (or `""`) for a trait without a joint run.
#' @param targetAnnoDir Character. Directory of target `.annot.gz` files used
#'   for `sd_C` and binary detection (typically the joint-mode dir).
#' @param frqfileDir Character or NULL.
#' @param plinkName Character. Default `"ADSP_chr"`.
#' @param mafCutoff Numeric, default `0.05`.
#' @param targetCategories Character vector or NULL. Auto-detected from the
#'   first available run if NULL.
#' @param targetLabels Character vector or NULL. Optional user-friendly display
#'   names for the target annotations, same length and order as the resolved
#'   `targetCategories` (e.g. `c("quantile_eQTL", "eQTL")` to replace the
#'   polyfun `.results` names `c("ANNOT_1_0", "ANNOT_2_0")`). When given, every
#'   `target` column and `tau*`-block column name in the output is renamed;
#'   `params$target_categories` then holds the labels and
#'   `params$target_categories_orig` keeps the original polyfun names. When NULL
#'   (default), nothing is renamed - the original `.results` category names are
#'   used as before.
#'
#' @return List with `per_trait`, `meta` (three frames), `params`.
#'
#' @export
