#' @title ColocBoost multi-trait colocalization pipeline (S4)
#' @description Protocol-level multi-trait colocalization analysis using
#'   \pkg{colocboost}. Dispatches on the QTL input type:
#'   \itemize{
#'     \item \code{QtlDataset} — single-study, individual-level
#'           multi-context data. Per-context residualized X / Y are
#'           extracted from the dataset (filtering knobs on the
#'           constructor apply lazily inside the accessors).
#'     \item \code{QtlSumStats} — summary-statistic-only QTL data with a
#'           shared LD reference (\code{ldSketch}). Must already have
#'           been passed through \code{\link{summaryStatsQc}} (the
#'           pipeline rejects inputs whose \code{getQcInfo()} is empty).
#'     \item \code{MultiTaskQtlDataset} — a mixture of one or more
#'           individual-level \code{QtlDataset} studies and an optional
#'           \code{QtlSumStats} collection.
#'   }
#'   GWAS is optional and always passed separately as a
#'   \code{GwasSumStats} object (must also be QC'd).
#'
#'   \code{colocboostPipeline} does \strong{not} accept a
#'   \code{FineMappingResult} for either side; colocboost has its own
#'   variable-selection algorithm.
#'
#' @section QC contract:
#'   \itemize{
#'     \item Individual-level QC (MAF / MAC / X-variance / per-sample
#'           missingness, sample / variant restrictions) lives on the
#'           \code{QtlDataset} constructor and is applied lazily inside
#'           \code{getGenotypes()} / \code{getResidualizedGenotypes()}.
#'           The pipeline does \emph{not} run a separate
#'           individual-level QC pass.
#'     \item All summary-statistic QC (variant filters, harmonization
#'           against the \code{ldSketch}, LD-mismatch detection, RAISS
#'           imputation, etc.) lives in
#'           \code{\link{summaryStatsQc}}. The pipeline rejects any
#'           \code{QtlSumStats} or \code{GwasSumStats} where
#'           \code{length(getQcInfo(x)) == 0L}.
#'   }
#'
#' @section Analysis variants:
#'   \itemize{
#'     \item \code{xqtlColoc} (default \code{TRUE}): run a colocboost
#'           model over the QTL contexts only (individual-level inputs).
#'     \item \code{jointGwas} (default \code{FALSE}): run a non-focal
#'           colocboost model that combines all QTL contexts/studies
#'           with the supplied \code{gwasSumStats} studies.
#'     \item \code{separateGwas} (default \code{FALSE}): run one focal
#'           colocboost model per GWAS study, where the GWAS is the
#'           focal outcome.
#'   }
#'
#' @param qtlData One of \code{QtlDataset}, \code{QtlSumStats}, or
#'   \code{MultiTaskQtlDataset}.
#' @param gwasSumStats Optional \code{GwasSumStats} with the GWAS
#'   studies to colocalize against. \code{NULL} to skip GWAS
#'   colocalization.
#' @param contexts Optional character vector of context names to
#'   restrict the individual-level / QtlSumStats QTL analysis to. When
#'   \code{NULL} (default), every context present is used.
#' @param traitId Optional character vector of trait identifiers to
#'   restrict the analysis to. When supplied with an individual-level
#'   \code{QtlDataset} input, \code{cisWindow} is required (passed to
#'   \code{getResidualizedGenotypes} / \code{getPhenotypes} for the
#'   variant-window selection).
#' @param region Optional single-range \code{GRanges} describing the
#'   analysis window. Mutually exclusive with \code{traitId} (see the
#'   \code{QtlDataset} accessors).
#' @param cisWindow Optional cis window in basepairs; required with
#'   \code{traitId}, optional with \code{region}.
#' @param focalTrait Optional trait name; when supplied and present in
#'   the assembled outcome list, the colocboost xQTL-only run uses it
#'   as the focal outcome.
#' @param eventFilters Optional list of per-context event-name filter
#'   specifications (each a list with \code{type_pattern} and at least
#'   one of \code{valid_pattern} / \code{exclude_pattern}). Applied to
#'   the trait columns within each context's residualized phenotype
#'   block before the colocboost call.
#' @param xqtlColoc,jointGwas,separateGwas Logical flags selecting which
#'   colocboost variants to run.
#' @param ... Additional arguments forwarded to
#'   \code{\link[colocboost]{colocboost}} (e.g., \code{M}, \code{L},
#'   \code{output_level}).
#' @return A list with elements \code{xqtl_coloc}, \code{joint_gwas},
#'   \code{separate_gwas}, and \code{computing_time}.
#' @name colocboostPipeline
#' @importFrom methods is setGeneric setMethod
#' @importFrom S4Vectors mcols
#' @importFrom GenomicRanges seqnames start end GRanges
#' @importFrom IRanges IRanges
#' @export
NULL

# =============================================================================
# Generic
# =============================================================================

#' @rdname colocboostPipeline
#' @export
setGeneric("colocboostPipeline",
  function(qtlData, gwasSumStats = NULL, ...) {
    standardGeneric("colocboostPipeline")
  })

# =============================================================================
# Helpers (private)
# =============================================================================

# Run colocboost() with tryCatch + timing.
.cbRun <- function(label, args) {
  if (!requireNamespace("colocboost", quietly = TRUE)) {
    stop("The colocboost package is required for colocboostPipeline().")
  }
  t1 <- Sys.time()
  args <- Filter(Negate(is.null), args)
  res <- tryCatch(
    do.call(colocboost::colocboost, args),
    error = function(e) {
      message(label, " failed: ", conditionMessage(e))
      NULL
    }
  )
  list(result = res, time = Sys.time() - t1)
}

# Build the LD / X_ref slot of the colocboost call from a list of LD
# matrices. When any matrix is non-square it is treated as a samples x
# variants genotype reference and routed to X_ref; otherwise routed to LD.
.cbBuildLdArgs <- function(ldList) {
  ldList <- Filter(Negate(is.null), ldList)
  if (length(ldList) == 0L) return(list())
  isGeno <- any(vapply(ldList,
                       function(m) nrow(m) != ncol(m), logical(1)))
  if (isGeno) list(X_ref = ldList) else list(LD = ldList)
}

# Reject SumStats objects that have not been passed through
# summaryStatsQc(). Both QtlSumStats and GwasSumStats expose getQcInfo();
# an empty list (the constructor default) signals "no QC run".
.cbRequireSumStatsQc <- function(x, what) {
  if (is.null(x)) return(invisible(NULL))
  if (length(getQcInfo(x)) == 0L) {
    stop(
      sprintf(
        "%s must be passed through summaryStatsQc() before reaching ",
        what),
      "colocboostPipeline (getQcInfo() returned an empty list). ",
      "Call summaryStatsQc(x, ...) and pass the result.")
  }
  invisible(NULL)
}

# Apply user-supplied event filters to the columns of a single
# residualized phenotype matrix. Returns the filtered matrix, or NULL
# when no events pass.
.cbApplyEventFilters <- function(y, filters, context) {
  if (is.null(y) || length(filters) == 0L) return(y)
  events <- colnames(y)
  if (is.null(events) || length(events) == 0L) return(y)
  kept <- events
  for (filter in filters) {
    if (is.null(filter$type_pattern) ||
        (is.null(filter$valid_pattern) && is.null(filter$exclude_pattern))) {
      stop("Each event filter must specify type_pattern and at least ",
           "one of valid_pattern or exclude_pattern.")
    }
    typeEvents <- kept[grepl(filter$type_pattern, kept)]
    if (length(typeEvents) == 0L) next
    if (!is.null(filter$valid_pattern)) {
      validGroups <- unique(gsub(filter$type_pattern, "\\1",
                                 typeEvents[grepl(filter$valid_pattern,
                                                  typeEvents)]))
      if (length(validGroups) == 0L) {
        typeEvents <- character(0)
      } else {
        typeEvents <- events[grepl(paste(validGroups, collapse = "|"),
                                   events)]
      }
    }
    if (!is.null(filter$exclude_pattern)) {
      typeEvents <- typeEvents[!grepl(filter$exclude_pattern, typeEvents)]
    }
    if (length(typeEvents) == length(events)) {
      message("All events matching ", filter$type_pattern,
              " in ", context, " included in following analysis.")
    } else if (length(typeEvents) == 0L) {
      message("No events matching ", filter$type_pattern,
              " in ", context, " pass the filtering.")
      return(NULL)
    } else {
      excludeEvents <- paste(setdiff(events, typeEvents), collapse = ";")
      message("Some events, ", excludeEvents,
              " in ", context, " are removed.")
    }
    kept <- unique(c(kept[!grepl(filter$type_pattern, kept)], typeEvents))
  }
  keep <- intersect(colnames(y), kept)
  if (length(keep) == 0L) return(NULL)
  y[, keep, drop = FALSE]
}

# Materialise an individual-level QtlDataset into the colocboost
# (X, Y, dict_YX, outcome_names) bundle. Each context becomes one X /
# Y pair; the YA matrices are split into single-trait columns and
# dict_YX maps each split column back to its X. Returns NULL when no
# context survives selection / event filtering.
.cbIndividualBundle <- function(qd, contexts = NULL,
                                traitId = NULL,
                                region = NULL,
                                cisWindow = NULL,
                                samples = NULL,
                                eventFilters = NULL) {
  if (is.null(contexts) || length(contexts) == 0L) {
    contexts <- getContexts(qd)
  } else {
    available <- getContexts(qd)
    bad <- setdiff(contexts, available)
    if (length(bad) > 0L) {
      stop("Unknown context(s) for QtlDataset '", getStudy(qd),
           "': ", paste(bad, collapse = ", "),
           ". Available: ", paste(available, collapse = ", "))
    }
  }

  YperCtx <- list()
  XperCtx <- list()
  for (ctx in contexts) {
    yList <- tryCatch(
      getResidualizedPhenotypes(qd, contexts = ctx,
                                traitId = traitId, region = region),
      error = function(e) {
        message("colocboostPipeline: skipping context '", ctx,
                "' (residualized phenotypes unavailable: ",
                conditionMessage(e), ").")
        NULL
      })
    if (is.null(yList) || !ctx %in% names(yList)) next
    Y <- yList[[ctx]]
    Y <- .cbApplyEventFilters(Y, eventFilters, ctx)
    if (is.null(Y) || ncol(Y) == 0L) next
    X <- tryCatch(
      getResidualizedGenotypes(qd, contexts = ctx,
                               traitId = traitId, region = region,
                               cisWindow = cisWindow, samples = samples),
      error = function(e) {
        message("colocboostPipeline: skipping context '", ctx,
                "' (residualized genotypes unavailable: ",
                conditionMessage(e), ").")
        NULL
      })
    if (is.null(X) || ncol(X) == 0L) next
    common <- intersect(rownames(X), rownames(Y))
    if (length(common) == 0L) {
      message("colocboostPipeline: skipping context '", ctx,
              "' (no samples shared between residualized X and Y).")
      next
    }
    X <- X[common, , drop = FALSE]
    Y <- Y[common, , drop = FALSE]
    XperCtx[[ctx]] <- X
    YperCtx[[ctx]] <- Y
  }

  if (length(XperCtx) == 0L) return(NULL)

  # Deduplicate X matrices that are identical across contexts so dict_YX
  # can fan out to a smaller X set (mirrors the legacy behaviour).
  uniqueX <- list()
  xMatch <- integer(length(XperCtx))
  for (i in seq_along(XperCtx)) {
    matched <- names(uniqueX)[
      vapply(uniqueX, identical, logical(1), XperCtx[[i]])]
    if (length(matched) > 0L) {
      xMatch[[i]] <- match(matched[[1L]], names(uniqueX))
    } else {
      uniqueX[[names(XperCtx)[i]]] <- XperCtx[[i]]
      xMatch[[i]] <- length(uniqueX)
    }
  }

  # Split each Y matrix into single-trait columns. When the same trait
  # name appears in multiple contexts, qualify it with the context prefix
  # so colocboost sees distinct outcome names.
  allTraitNames <- unlist(lapply(YperCtx, colnames), use.names = FALSE)
  dupTraits <- unique(allTraitNames[
    duplicated(allTraitNames) | duplicated(allTraitNames, fromLast = TRUE)])
  YSplit <- list()
  dict <- matrix(integer(0), ncol = 2L)
  for (i in seq_along(YperCtx)) {
    Y <- YperCtx[[i]]
    ctx <- names(YperCtx)[i]
    for (j in seq_len(ncol(Y))) {
      tname <- colnames(Y)[j]
      if (is.null(tname) || is.na(tname) || tname == "") {
        tname <- paste0("outcome", length(YSplit) + 1L)
      }
      if (tname %in% dupTraits) {
        tname <- paste0(ctx, "_", tname)
      }
      if (tname %in% names(YSplit)) {
        tname <- make.unique(c(names(YSplit), tname))[length(YSplit) + 1L]
      }
      YSplit[[tname]] <- Y[, j, drop = FALSE]
      dict <- rbind(dict, c(length(YSplit), xMatch[[i]]))
    }
  }
  colnames(dict) <- c("Y", "X")

  list(
    X            = uniqueX,
    Y            = YSplit,
    dict_YX      = dict,
    outcomeNames = names(YSplit))
}

# Build a single (sumstat data.frame, LD correlation matrix) pair from a
# QtlSumStats / GwasSumStats entry. Returns NULL when the entry has no
# variants overlapping the ldSketch panel.
.cbSumstatPair <- function(gr, ldSketch, varY = NULL) {
  if (is.null(gr) || length(gr) == 0L) return(NULL)
  mc <- S4Vectors::mcols(gr)
  if (!"Z" %in% colnames(mc)) {
    stop("colocboostPipeline: sumstat entry is missing the Z mcol; ",
         "summaryStatsQc() should have written it.")
  }
  variantIds <- if ("SNP" %in% colnames(mc)) {
    as.character(mc$SNP)
  } else {
    formatVariantId(as.character(GenomicRanges::seqnames(gr)),
                    GenomicRanges::start(gr),
                    as.character(mc$A2),
                    as.character(mc$A1))
  }
  panelIds <- as.character(getSnpInfo(ldSketch)$SNP)
  panelMatch <- match(variantIds, panelIds)
  keep <- !is.na(panelMatch)
  if (!any(keep)) return(NULL)
  gr <- gr[keep]
  variantIds <- variantIds[keep]
  panelIdx <- panelMatch[keep]
  mc <- S4Vectors::mcols(gr)

  block <- extractBlockGenotypes(ldSketch, panelIdx, meanImpute = TRUE)
  X <- t(SummarizedExperiment::assay(block, "dosage"))
  colnames(X) <- variantIds
  R <- computeLd(X, method = "sample")

  N <- if ("N" %in% colnames(mc)) as.numeric(mc$N) else NA_real_
  ss <- data.frame(
    z       = as.numeric(mc$Z),
    n       = N,
    variant = variantIds,
    stringsAsFactors = FALSE)
  if (!is.null(varY) && !is.na(varY)) {
    ss$var_y <- varY
  }
  list(sumstat = ss, LD = R, variantIds = variantIds)
}

# Build the colocboost sumstat / LD bundle from a QtlSumStats and an
# optional contexts / traitId filter. Returns a list keyed by sumstat
# study label, where each entry has (sumstat, LD).
.cbQtlSumStatsBundle <- function(ss, contexts = NULL, traitId = NULL) {
  if (is.null(ss) || nrow(ss) == 0L) return(list())
  ldSketch <- .cbLdSketch(ss)
  keepRow <- rep(TRUE, nrow(ss))
  if (!is.null(contexts) && length(contexts) > 0L) {
    keepRow <- keepRow & as.character(ss$context) %in% contexts
  }
  if (!is.null(traitId) && length(traitId) > 0L) {
    keepRow <- keepRow & as.character(ss$trait) %in% traitId
  }
  if (!any(keepRow)) return(list())
  rows <- which(keepRow)

  bundle <- list()
  for (i in rows) {
    label <- paste(as.character(ss$study)[[i]],
                   as.character(ss$context)[[i]],
                   as.character(ss$trait)[[i]], sep = ":")
    pair <- .cbSumstatPair(
      gr       = ss$entry[[i]],
      ldSketch = ldSketch,
      varY     = if ("varY" %in% names(ss)) ss$varY[[i]] else NA_real_)
    if (!is.null(pair)) bundle[[label]] <- pair
  }
  bundle
}

# Same as .cbQtlSumStatsBundle for a GwasSumStats collection, keyed by
# study label.
.cbGwasSumStatsBundle <- function(gws) {
  if (is.null(gws) || nrow(gws) == 0L) return(list())
  ldSketch <- .cbLdSketch(gws)
  bundle <- list()
  for (i in seq_len(nrow(gws))) {
    label <- as.character(gws$study)[[i]]
    pair <- .cbSumstatPair(
      gr       = gws$entry[[i]],
      ldSketch = ldSketch,
      varY     = if ("varY" %in% names(gws)) gws$varY[[i]] else NA_real_)
    if (!is.null(pair)) bundle[[label]] <- pair
  }
  bundle
}

# Pipeline-side wrapper: delegate to the public getLdSketch generic.
.cbLdSketch <- function(x) getLdSketch(x)

# Compare two GenotypeHandles for the LD-sketch equality contract.
# Two handles match iff they have identical snpInfo (variant id, CHR,
# BP, A1, A2 in the same order) and identical sampleIds. Returns silently
# when matched; hard error on mismatch.
.cbRequireMatchingLdSketches <- function(qtlLd, gwasLd) {
  if (is.null(qtlLd) || is.null(gwasLd)) return(invisible(NULL))
  if (!methods::is(qtlLd, "GenotypeHandle") ||
      !methods::is(gwasLd, "GenotypeHandle")) {
    stop("colocboostPipeline: both ldSketch slots must be GenotypeHandle ",
         "objects for the cross-pipeline LD reference check.")
  }
  qSnps <- getSnpInfo(qtlLd)
  gSnps <- getSnpInfo(gwasLd)
  if (nrow(qSnps) != nrow(gSnps))
    stop("colocboostPipeline: QTL and GWAS ldSketch panels have different ",
         "variant counts (", nrow(qSnps), " vs ", nrow(gSnps),
         "). The two ldSketch GenotypeHandles must match exactly.")
  for (col in c("SNP", "CHR", "BP", "A1", "A2")) {
    if (!identical(as.character(qSnps[[col]]),
                   as.character(gSnps[[col]])))
      stop("colocboostPipeline: QTL and GWAS ldSketch panels differ in ",
           col, ". Use the same ldSketch on both sides.")
  }
  if (!identical(getSampleIds(qtlLd), getSampleIds(gwasLd)))
    stop("colocboostPipeline: QTL and GWAS ldSketch panels have different ",
         "sample sets. Use the same ldSketch on both sides.")
  invisible(NULL)
}

# Combine sumstat bundles into the (sumstat-list, LD-list, dict) shape
# colocboost expects. Deduplicates identical LD matrices so dict_sumstatLD
# can point multiple sumstats at one LD.
.cbMergeSumstatBundles <- function(bundles) {
  if (length(bundles) == 0L) {
    return(list(sumstat = list(), LD = list(),
                dict_sumstatLD = matrix(integer(0), ncol = 2L)))
  }
  ldUnique <- list()
  ldMatch <- integer(length(bundles))
  for (i in seq_along(bundles)) {
    ld <- bundles[[i]]$LD
    matched <- which(vapply(ldUnique, identical, logical(1), ld))
    if (length(matched) > 0L) {
      ldMatch[[i]] <- matched[[1L]]
    } else {
      ldUnique[[length(ldUnique) + 1L]] <- ld
      ldMatch[[i]] <- length(ldUnique)
    }
  }
  names(ldUnique) <- paste0("LD", seq_along(ldUnique))
  sumstat <- lapply(bundles, `[[`, "sumstat")
  names(sumstat) <- names(bundles)
  dict <- cbind(seq_along(bundles), ldMatch)
  colnames(dict) <- c("sumstat", "LD")
  list(sumstat = sumstat, LD = ldUnique, dict_sumstatLD = dict)
}

# Build an empty result skeleton consistent with what the per-method
# dispatch fills in.
.cbEmptyResult <- function() {
  list(
    xqtl_coloc     = NULL,
    joint_gwas     = NULL,
    separate_gwas  = NULL,
    computing_time = list(Analysis = list(
      xqtl_coloc     = NULL,
      joint_gwas     = NULL,
      separate_gwas  = NULL)))
}

# Shared dispatch: accepts a fully-prepared individual bundle (possibly
# NULL) plus a sumstat bundle (possibly empty) and runs the three
# colocboost variants the user requested.
.cbRunVariants <- function(individualBundle, sumstatBundle,
                           xqtlColoc, jointGwas, separateGwas,
                           focalTrait, dotArgs) {
  results <- .cbEmptyResult()
  hasInd <- !is.null(individualBundle)
  hasSs  <- length(sumstatBundle$sumstat) > 0L

  if (!hasInd && !hasSs) {
    message("colocboostPipeline: no QTL inputs remain after selection. ",
            "Nothing to run.")
    return(results)
  }

  # xQTL-only run
  if (isTRUE(xqtlColoc) && hasInd) {
    traits <- individualBundle$outcomeNames
    focalIdx <- if (!is.null(focalTrait) && focalTrait %in% traits)
      which(traits == focalTrait) else NULL
    message("====== Performing xQTL-only ColocBoost on ",
            length(individualBundle$Y), " contexts. =====")
    args <- c(
      list(X                  = individualBundle$X,
           Y                  = individualBundle$Y,
           dict_YX            = individualBundle$dict_YX,
           outcome_names      = traits,
           focal_outcome_idx  = focalIdx,
           output_level       = 2),
      dotArgs)
    run <- .cbRun("xQTL-only ColocBoost", args)
    results$xqtl_coloc <- run$result
    results$computing_time$Analysis$xqtl_coloc <- run$time
  }

  # Joint (non-focal) QTL + GWAS run
  if (isTRUE(jointGwas) && hasSs) {
    traits <- c(if (hasInd) individualBundle$outcomeNames else character(),
                names(sumstatBundle$sumstat))
    ldArgs <- .cbBuildLdArgs(sumstatBundle$LD)
    nContexts <- if (hasInd) length(individualBundle$Y) else 0L
    message("====== Performing non-focal GWAS-xQTL ColocBoost on ",
            nContexts, " contexts and ",
            length(sumstatBundle$sumstat), " GWAS. =====")
    args <- c(
      list(X                  = if (hasInd) individualBundle$X else NULL,
           Y                  = if (hasInd) individualBundle$Y else NULL,
           sumstat            = sumstatBundle$sumstat,
           dict_YX            = if (hasInd) individualBundle$dict_YX else NULL,
           dict_sumstatLD     = sumstatBundle$dict_sumstatLD,
           outcome_names      = traits,
           focal_outcome_idx  = NULL,
           output_level       = 2),
      ldArgs,
      dotArgs)
    run <- .cbRun("Joint GWAS ColocBoost", args)
    results$joint_gwas <- run$result
    results$computing_time$Analysis$joint_gwas <- run$time
  }

  # Separate (focal) per-GWAS runs
  if (isTRUE(separateGwas) && hasSs) {
    ssNames <- names(sumstatBundle$sumstat)
    separate <- vector("list", length(ssNames))
    names(separate) <- ssNames
    t1 <- Sys.time()
    for (i in seq_along(ssNames)) {
      study <- ssNames[[i]]
      ldIdx <- sumstatBundle$dict_sumstatLD[i, 2L]
      ldArgs <- .cbBuildLdArgs(sumstatBundle$LD[ldIdx])
      traits <- c(if (hasInd) individualBundle$outcomeNames else character(),
                  study)
      nContexts <- if (hasInd) length(individualBundle$Y) else 0L
      message("====== Performing focal GWAS-xQTL ColocBoost on ",
              nContexts, " contexts and ", study, " GWAS. =====")
      args <- c(
        list(X                  = if (hasInd) individualBundle$X else NULL,
             Y                  = if (hasInd) individualBundle$Y else NULL,
             sumstat            = sumstatBundle$sumstat[i],
             dict_YX            = if (hasInd) individualBundle$dict_YX else NULL,
             outcome_names      = traits,
             focal_outcome_idx  = length(traits),
             output_level       = 2),
        ldArgs,
        dotArgs)
      run <- .cbRun(paste("Separate GWAS ColocBoost for", study), args)
      separate[[study]] <- run$result
    }
    t2 <- Sys.time()
    results$separate_gwas <- separate
    results$computing_time$Analysis$separate_gwas <- list(
      total = t2 - t1,
      n_studies = length(ssNames),
      average = if (length(ssNames) > 0L) (t2 - t1) / length(ssNames) else NA)
  }

  results
}

# Top-level driver shared by all input methods. qtlPairs and gwasPairs
# are per-tuple lists of `list(sumstat, LD)` produced by the per-class
# bundle helpers; they are merged here so dict_sumstatLD can dedupe
# identical LD matrices across QTL and GWAS sides.
.cbDriver <- function(individualBundle, qtlPairs, gwasSumStats,
                      xqtlColoc, jointGwas, separateGwas,
                      focalTrait, dotArgs, qtlLdSketch = NULL) {
  if (!isTRUE(xqtlColoc) && !isTRUE(jointGwas) && !isTRUE(separateGwas)) {
    message("colocboostPipeline: no analysis flag is TRUE; nothing to do.")
    return(.cbEmptyResult())
  }

  combinedPairs <- qtlPairs
  if (!is.null(gwasSumStats)) {
    .cbRequireSumStatsQc(gwasSumStats, "gwasSumStats")
    if (!is.null(qtlLdSketch)) {
      .cbRequireMatchingLdSketches(qtlLdSketch, getLdSketch(gwasSumStats))
    }
    gwasPairs <- .cbGwasSumStatsBundle(gwasSumStats)
    for (label in names(gwasPairs)) {
      key <- label
      if (key %in% names(combinedPairs)) {
        key <- make.unique(
          c(names(combinedPairs), key))[length(combinedPairs) + 1L]
      }
      combinedPairs[[key]] <- gwasPairs[[label]]
    }
  }
  sumstatBundle <- .cbMergeSumstatBundles(combinedPairs)

  .cbRunVariants(individualBundle, sumstatBundle,
                 xqtlColoc, jointGwas, separateGwas,
                 focalTrait, dotArgs)
}

# =============================================================================
# Methods
# =============================================================================

#' @rdname colocboostPipeline
#' @export
setMethod("colocboostPipeline", "QtlDataset",
  function(qtlData, gwasSumStats = NULL,
           contexts = NULL,
           traitId = NULL, region = NULL, cisWindow = NULL,
           focalTrait = NULL,
           eventFilters = NULL,
           xqtlColoc = TRUE,
           jointGwas = FALSE,
           separateGwas = FALSE,
           samples = NULL,
           ...) {
    dotArgs <- list(...)
    indBundle <- .cbIndividualBundle(
      qtlData,
      contexts     = contexts,
      traitId      = traitId,
      region       = region,
      cisWindow    = cisWindow,
      samples      = samples,
      eventFilters = eventFilters)
    .cbDriver(indBundle, qtlPairs = list(), gwasSumStats,
              xqtlColoc, jointGwas, separateGwas,
              focalTrait, dotArgs)
  })

#' @rdname colocboostPipeline
#' @export
setMethod("colocboostPipeline", "QtlSumStats",
  function(qtlData, gwasSumStats = NULL,
           contexts = NULL,
           traitId = NULL, region = NULL, cisWindow = NULL,
           focalTrait = NULL,
           eventFilters = NULL,
           xqtlColoc = TRUE,
           jointGwas = FALSE,
           separateGwas = FALSE,
           ...) {
    if (!is.null(eventFilters)) {
      warning("eventFilters is ignored for QtlSumStats input (no per-trait ",
              "phenotype matrix to filter).")
    }
    .cbRequireSumStatsQc(qtlData, "qtlData")
    dotArgs <- list(...)
    qtlPairs <- .cbQtlSumStatsBundle(
      qtlData, contexts = contexts, traitId = traitId)
    .cbDriver(individualBundle = NULL,
              qtlPairs         = qtlPairs,
              gwasSumStats     = gwasSumStats,
              xqtlColoc        = xqtlColoc,
              jointGwas        = jointGwas,
              separateGwas     = separateGwas,
              focalTrait       = focalTrait,
              dotArgs          = dotArgs,
              qtlLdSketch      = getLdSketch(qtlData))
  })

#' @rdname colocboostPipeline
#' @export
setMethod("colocboostPipeline", "MultiTaskQtlDataset",
  function(qtlData, gwasSumStats = NULL,
           contexts = NULL,
           traitId = NULL, region = NULL, cisWindow = NULL,
           focalTrait = NULL,
           eventFilters = NULL,
           xqtlColoc = TRUE,
           jointGwas = FALSE,
           separateGwas = FALSE,
           samples = NULL,
           ...) {
    dotArgs <- list(...)

    # Aggregate the individual-level bundles across all QtlDataset
    # members. Per-study trait names are prefixed with "{study}:" so
    # colocboost sees distinct outcomes when two studies share a trait.
    qtlDatasets <- getQtlDatasets(qtlData)
    combinedX <- list()
    combinedY <- list()
    combinedDict <- matrix(integer(0), ncol = 2L)
    colnames(combinedDict) <- c("Y", "X")
    combinedOutcomes <- character()
    for (study in names(qtlDatasets)) {
      qd <- qtlDatasets[[study]]
      sub <- .cbIndividualBundle(
        qd,
        contexts     = contexts,
        traitId      = traitId,
        region       = region,
        cisWindow    = cisWindow,
        samples      = samples,
        eventFilters = eventFilters)
      if (is.null(sub)) next
      xOffset <- length(combinedX)
      yOffset <- length(combinedY)
      newXNames <- paste(study, names(sub$X), sep = ":")
      newYNames <- paste(study, sub$outcomeNames, sep = ":")
      names(sub$X) <- newXNames
      names(sub$Y) <- newYNames
      combinedX <- c(combinedX, sub$X)
      combinedY <- c(combinedY, sub$Y)
      shifted <- sub$dict_YX
      shifted[, "Y"] <- shifted[, "Y"] + yOffset
      shifted[, "X"] <- shifted[, "X"] + xOffset
      combinedDict <- rbind(combinedDict, shifted)
      combinedOutcomes <- c(combinedOutcomes, newYNames)
    }

    indBundle <- if (length(combinedX) > 0L) {
      list(X            = combinedX,
           Y            = combinedY,
           dict_YX      = combinedDict,
           outcomeNames = combinedOutcomes)
    } else NULL

    # Sumstat side: any QtlSumStats embedded in the MultiTaskQtlDataset.
    embeddedSs <- getSumStats(qtlData)
    qtlPairs <- list()
    qtlLdSketch <- NULL
    if (!is.null(embeddedSs)) {
      .cbRequireSumStatsQc(embeddedSs,
                           "MultiTaskQtlDataset@sumStats")
      qtlPairs <- .cbQtlSumStatsBundle(
        embeddedSs, contexts = contexts, traitId = traitId)
      qtlLdSketch <- getLdSketch(embeddedSs)
    }

    .cbDriver(indBundle, qtlPairs, gwasSumStats,
              xqtlColoc, jointGwas, separateGwas,
              focalTrait, dotArgs,
              qtlLdSketch = qtlLdSketch)
  })

#' @rdname colocboostPipeline
#' @export
setMethod("colocboostPipeline", "ANY",
  function(qtlData, gwasSumStats = NULL, ...) {
    stop("colocboostPipeline does not accept inputs of class '",
         class(qtlData)[[1L]], "'. Pass a QtlDataset, QtlSumStats, or ",
         "MultiTaskQtlDataset for QTL data.")
  })
