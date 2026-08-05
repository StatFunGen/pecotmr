# Filter rows of a z-score matrix by significance p-value cutoff.
# Returns integer indices of rows where any |z| exceeds the threshold.
# @noRd
filterBySignificance <- function(zMatrix, sigPCutoff) {
  zThreshold <- sqrt(qchisq(sigPCutoff, df = 1, lower.tail = FALSE))
  which(apply(zMatrix, 1, function(row) any(abs(row) >= zThreshold)))
}

# Coerce every column to numeric, then replace NaN/Inf/NA with `replaceWith`.
# @noRd
.mashReplaceValues <- function(df, replaceWith) {
  df <- df %>%
    mutate(across(everything(), as.numeric)) %>%
    mutate(across(everything(), ~ replace(., is.nan(.) | is.infinite(.) | is.na(.), replaceWith)))
}

# Coerce z-scores to a matrix (NaN/Inf/NA -> 0) and, when a missing-rate
# threshold is given, drop rows falling below it.
# @noRd
.mashProcessZ <- function(zData, filterByMissingRate) {
  zData <- as.matrix(.mashReplaceValues(zData, 0))

  if (!is.null(filterByMissingRate)) {
    proportionNonzero <- apply(zData, 1, function(row) mean(row != 0))
    zData <- zData[proportionNonzero >= filterByMissingRate, , drop = FALSE]
  }

  return(zData)
}

#' @importFrom vroom vroom
#' @export
filterInvalidSummaryStat <- function(datList, bhat = NULL, sbhat = NULL, z = NULL, btoz = FALSE, sigPCutoff = 1E-6, filterByMissingRate = 0.2) {
  # Function to process bhat, sbhat
  if (!is.null(bhat) && !is.null(sbhat) && all(c(bhat, sbhat) %in% names(datList))) {
    # If the element is a list with 'bhat' and 'sbhat'
    if (!is.null(datList[[bhat]]) && !is.null(datList[[sbhat]])) {
      datList[[bhat]] <- as.matrix(.mashReplaceValues(datList[[bhat]], 0))
      datList[[sbhat]] <- as.matrix(.mashReplaceValues(datList[[sbhat]], 1000))
      if (("null.b" %in% names(datList)) || ("random.b" %in% names(datList))) {
        if (!is.null(filterByMissingRate)) {
          proportionNonzero <- apply(datList[[bhat]], 1, function(row) {
            mean(row != 0)
          })
          datList[[bhat]] <- datList[[bhat]][proportionNonzero >= filterByMissingRate, ]
          datList[[sbhat]] <- datList[[sbhat]][proportionNonzero >= filterByMissingRate, ]
        }
      }
    }
  }
  # Function to filter strong signal using z score
  if (btoz) {
    if (any(grepl("\\.b$", bhat)) | any(grepl("\\.s$", sbhat))) {
      condition <- sub("\\.b$", "", bhat)
      if (!is.null(datList[[bhat]]) && !is.null(datList[[sbhat]])) {
        datList[[paste0(condition, ".z")]] <- as.matrix(datList[[bhat]] / datList[[sbhat]])
      } else {
        datList[paste0(condition, ".z")] <- list(NULL)
      }
    } else {
      if (!is.null(datList[[bhat]]) && !is.null(datList[[sbhat]])) {
        datList[["z"]] <- as.matrix(datList[[bhat]] / datList[[sbhat]])
      } else {
        datList["z"] <- list(NULL)
      }
    }
    if ("strong.z" %in% names(datList)) {
      if (!is.null(sigPCutoff)) {
        keepIndex <- filterBySignificance(datList$strong.z, sigPCutoff)
        datList[["strong.z"]] <- datList$strong.z[keepIndex, ]
        datList[["strong.b"]] <- datList$strong.b[keepIndex, ]
        datList[["strong.s"]] <- datList$strong.s[keepIndex, ]
      }
    }
  }
  # Function to process z-scores and filter directly
  if (!is.null(z)) {
    # Process each component if it exists
    for (comp in c("strong", "random", "null")) {
      if (!is.null(datList[[comp]]) && !is.null(datList[[comp]]$z)) {
        datList[[comp]]$z <- .mashProcessZ(datList[[comp]]$z, filterByMissingRate)
      }
    }

    # Apply significance cutoff to strong signals if applicable
    if (!is.null(datList$strong) && !is.null(datList$strong$z) && !is.null(sigPCutoff)) {
      keepIndex <- filterBySignificance(datList$strong$z, sigPCutoff)
      datList$strong$z <- datList$strong$z[keepIndex, , drop = FALSE]
    }
  }

  return(datList)
}

#' @importFrom purrr keep
#' @export
filterMixtureComponents <- function(conditionsToKeep, U, w = NULL, wCutoff = 1e-04) {
  # Identify conditions not to keep (to be removed)
  conditionsToFilter <- setdiff(colnames(U[[1]]), conditionsToKeep)
  sumW <- sum(w) # Original total sum of weights

  # Filter U by removing unwanted phenotypes (conditions)
  U <- lapply(U, function(mat, toKeep) {
    missingConditions <- setdiff(toKeep, colnames(mat))
    if (length(missingConditions) > 0) {
      stop(paste("Condition(s)", paste(missingConditions,
        collapse = ", "
      ), "not found in matrix"))
    }
    mat[toKeep, toKeep] # Keep only relevant conditions
  }, conditionsToKeep)

  # Remove matrices where all values are zero or weight is below cutoff
  keepNames <- names(keep(U, function(mat) !all(mat == 0)))
  if (!is.null(w)) {
    keepNames <- intersect(keepNames, names(w[w >= wCutoff]))
  }
  U <- U[keepNames]
  if (!is.null(w)) {
    w <- w[keepNames]
  }

  # Note: Matrices in U may contain very small values on the diagonal
  # even when contexts are not present, due to EM algorithm adjustments.
  # This makes the matrix not exactly zero, so it won't be removed even though it may not
  # have strong context relevance. This behavior arises because the algorithm attempts to
  # ensure matrices are full-rank, slightly changing initial values.

  # We cannot simply remove diagonal matrices as signals on the diagonal can be strong and relevant.
  # So we manually remove the U components that are driven by non-relevant contexts.
  U[conditionsToFilter] <- NULL
  w <- w[!names(w) %in% conditionsToFilter]

  # Recalculate the sum of remaining weights
  sumWnew <- sum(w)

  # Adjust weights to maintain the original sumW
  w <- (w / sumWnew) * sumW

  message(paste(length(U), "components of matrices remained after filtering."))

  return(list(U = U, w = w))
}


# Draw the random + null sub-samples used to estimate the null correlation.
# @noRd
.mashExtractOneData <- function(dat, nRandom, nNull) {
  if (is.null(dat)) {
    return(NULL)
  }

  if ("z" %in% names(dat)) {
    absZ <- abs(dat$z)
    zData <- dat$z
  } else {
    absZ <- abs(dat$bhat / dat$sbhat)
    zData <- NULL
  }

  sampleIdx <- 1:nrow(absZ)
  randomIdx <- sample(sampleIdx, min(nRandom, length(sampleIdx)), replace = FALSE)

  if (!is.null(zData)) {
    random <- list(z = zData[randomIdx, , drop = FALSE])
  } else {
    random <- list(
      bhat = dat$bhat[randomIdx, , drop = FALSE],
      sbhat = dat$sbhat[randomIdx, , drop = FALSE]
    )
  }

  null.id <- which(apply(absZ, 1, max) < 2)
  if (length(null.id) == 0) {
    warning(paste("no variants are included in the null dataset because absZ > 2 for all variants in", dat$region))
    null <- list()
  } else {
    if (length(null.id) < ncol(absZ)) {
      warning(paste("not enough null data to estimate null correlation in", dat$region))
      null <- list()
    } else {
      nullIdx <- sample(null.id, min(nNull, length(null.id)), replace = FALSE)
      if (!is.null(zData)) {
        null <- list(z = zData[nullIdx, , drop = FALSE])
      } else {
        null <- list(
          bhat = dat$bhat[nullIdx, , drop = FALSE],
          sbhat = dat$sbhat[nullIdx, , drop = FALSE]
        )
      }
    }
  }
  dat <- list(random = random, null = null)
  return(dat)
}

#' @export
mashRandNullSample <- function(dat, nRandom, nNull, excludeCondition, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (length(excludeCondition) > 0) {
    colsToCheck <- if ("z" %in% names(dat)) "z" else "bhat"
    if (!all(excludeCondition %in% colnames(dat[[colsToCheck]]))) {
      stop(paste("Error: excludeCondition are not present in", dat$region))
    }
    for (key in intersect(names(dat), c("z", "bhat", "sbhat"))) {
      keep <- setdiff(colnames(dat[[key]]), excludeCondition)
      dat[[key]] <- dat[[key]][, keep, drop = FALSE]
    }
  }

  result <- .mashExtractOneData(dat, nRandom, nNull)
  return(result)
}

#' @export
mergeMashData <- function(resData, oneData) {
  if (length(resData) == 0 || is.null(resData)) return(oneData)
  if (length(oneData) == 0 || is.null(oneData)) return(resData)

  combinedData <- lapply(names(oneData), function(d) {
    od <- oneData[[d]]
    rd <- resData[[d]]
    if (length(od) == 0 || is.null(od)) return(rd)
    if (is.null(rd) || length(rd) == 0) return(od)

    # bind_rows auto-aligns columns, filling missing with NA; replace with NaN
    rnRes <- rownames(as.data.frame(rd))
    rnOne <- rownames(as.data.frame(od))
    combined <- bind_rows(as.data.frame(rd), as.data.frame(od))
    combined[is.na(combined)] <- NaN
    rnAll <- make.names(c(rnRes, rnOne), unique = TRUE)
    rownames(combined) <- rnAll
    combined
  })
  names(combinedData) <- names(oneData)
  return(combinedData)
}

# Build variants x conditions (Bhat, Shat) matrices for ONE object plus the
# row indices of its "strong" variants. Class dispatch:
#   QtlSumStats / GwasSumStats -> .mashSumStatsToMatrices; strong = the single
#       most significant variant (max|z|) per condition, unioned.
#   FineMappingResultBase      -> pivot getMarginalEffects() into a
#       variants x context (beta, se) pair; strong = the lead (max PIP) variant
#       of each credible set in each condition (getCs()), unioned. Conditions
#       with no credible set contribute no strong variant.
# For z-scale QtlSumStats the returned Shat is 1, so downstream code that forms
# z = b / s recovers the z-scores uniformly across both scales.
# @noRd
.mashObjectMatrices <- function(obj, inputScale, coverage) {
  if (methods::is(obj, "QtlSumStats") || methods::is(obj, "GwasSumStats")) {
    mats <- .mashSumStatsToMatrices(obj, "mash input", inputScale = inputScale)
    z <- mats$b / mats$s
    # most significant variant (max|z|) per condition column, unioned
    strongRows <- sort(unique(apply(abs(z), 2L, which.max)))
    return(list(b = mats$b, s = mats$s, strongRows = strongRows))
  }
  if (methods::is(obj, "FineMappingResultBase")) {
    me <- getMarginalEffects(obj)
    cs <- getCs(obj, coverage = coverage)
    if (!all(c("variant_id", "context", "beta", "se") %in% names(me))) {
      stop("mashInput: getMarginalEffects() must return variant_id/context/",
           "beta/se columns; a FineMappingResult with >= 2 contexts is required.")
    }
    # A multi-method result would duplicate (variant, context) cells; pin the
    # first method so the pivot is unambiguous.
    if ("method" %in% names(me) && length(unique(me$method)) > 1L) {
      m1 <- me$method[[1L]]
      warning("mashInput: FineMappingResult carries multiple methods; ",
              "using '", m1, "'.")
      me <- me[me$method == m1, , drop = FALSE]
      if ("method" %in% names(cs)) cs <- cs[cs$method == m1, , drop = FALSE]
    }
    contexts <- unique(me$context)
    variants <- unique(me$variant_id)
    b <- matrix(NA_real_, length(variants), length(contexts),
                dimnames = list(variants, contexts))
    s <- b
    cell <- cbind(match(me$variant_id, variants), match(me$context, contexts))
    b[cell] <- me$beta
    s[cell] <- me$se
    # strong = lead (max PIP) variant of each credible set in each condition
    strongVar <- character(0)
    csCol <- grep("^cs_", names(cs), value = TRUE)
    if (nrow(cs) > 0L && length(csCol) > 0L && "pip" %in% names(cs)) {
      grp <- interaction(cs$context, cs[[csCol[[1L]]]], drop = TRUE)
      for (rows in split(seq_len(nrow(cs)), grp)) {
        strongVar <- c(strongVar, cs$variant_id[[rows[[which.max(cs$pip[rows])]]]])
      }
    }
    strongRows <- sort(match(unique(strongVar), variants))
    return(list(b = b, s = s, strongRows = strongRows[!is.na(strongRows)]))
  }
  stop("mashInput: each element of `objects` must be a QtlSumStats, ",
       "GwasSumStats, or FineMappingResult; got ",
       paste(class(obj), collapse = "/"), ".")
}

# Extract the strong / random / null partitions from ONE object as a flat
# list(strong.b, strong.s, random.b, random.s, null.b, null.s) of
# variants x conditions matrices. Random / null are drawn by the shared
# mashRandNullSample() over the object's (Bhat, Shat); strong is the
# deterministic class-specific selection from .mashObjectMatrices().
# @noRd
.mashObjectPartitions <- function(obj, nRandom, nNull, excludeCondition,
                                  coverage, inputScale, seed,
                                  independentVariants = NULL) {
  mats <- .mashObjectMatrices(obj, inputScale = inputScale, coverage = coverage)
  keepCols <- setdiff(colnames(mats$b), excludeCondition)
  if (length(keepCols) < 2L) {
    stop("mashInput: fewer than 2 conditions remain for an object (after ",
         "excludeCondition); mash operates across conditions and needs >= 2.")
  }

  # Random / null candidate pool: optionally restrict to LD-independent variants
  # so the background carries no LD-correlated SNPs (which would bias Vhat and
  # the mixture weights). Matching is delegated to matchVariants() -- proper
  # chrom/pos/allele matching (ref/alt flips OK, no strand flip), NOT a raw
  # string compare. Strong is ALWAYS drawn from the full set (every top signal
  # is kept regardless of LD pruning).
  poolB <- mats$b; poolS <- mats$s
  if (!is.null(independentVariants) && length(independentVariants) > 0L) {
    # QtlSumStats matrix rownames carry a "study::trait::" block prefix; strip
    # it to the bare variant id before matching.
    rawIds <- sub(".*::", "", rownames(mats$b))
    keepIdx <- matchVariants(rawIds, independentVariants, allowFlip = TRUE,
                             removeStrandAmbiguous = FALSE)$idxA
    if (length(keepIdx) == 0L) {
      warning("mashInput: no variants matched the independent-variant list; ",
              "the random/null background is empty for this object.")
    }
    poolB <- mats$b[keepIdx, , drop = FALSE]
    poolS <- mats$s[keepIdx, , drop = FALSE]
  }

  rn <- mashRandNullSample(list(bhat = poolB, sbhat = poolS),
                           nRandom = nRandom, nNull = nNull,
                           excludeCondition = excludeCondition, seed = seed)
  out <- list()
  if (length(mats$strongRows) > 0L) {
    out[["strong.b"]] <- mats$b[mats$strongRows, keepCols, drop = FALSE]
    out[["strong.s"]] <- mats$s[mats$strongRows, keepCols, drop = FALSE]
  }
  if (!is.null(rn$random) && length(rn$random) > 0L) {
    out[["random.b"]] <- rn$random$bhat
    out[["random.s"]] <- rn$random$sbhat
  }
  if (!is.null(rn$null) && length(rn$null) > 0L) {
    out[["null.b"]] <- rn$null$bhat
    out[["null.s"]] <- rn$null$sbhat
  }
  out
}

#' Assemble MASH strong / random / null input from S4 objects
#'
#' Unified, S4-native replacement for the legacy
#' \code{load_multitrait_*_sumstat} + \code{mash_ran_null_sample} assembly.
#' Consumes a list of already-constructed objects (one per region) and returns
#' the flat \code{variants x conditions} matrix list consumed by the MASH
#' mixture-prior / fit / posterior steps.
#'
#' For EACH object three partitions are extracted:
#' \describe{
#'   \item{strong}{The high-signal variants (deterministic, class-specific).
#'     \code{QtlSumStats}: the single most significant variant (\eqn{\max|z|})
#'     per condition, unioned. \code{FineMappingResult}: the lead variant
#'     (\eqn{\max} PIP) of each credible set in each condition, unioned
#'     (conditions with no credible set contribute nothing).}
#'   \item{random}{\code{nRandom} variants sampled uniformly at random --
#'     represents the genome-wide mixture of effects and drives the mixture
#'     weights.}
#'   \item{null}{\code{nNull} variants sampled from those with \eqn{\max|z|<2}
#'     -- the noise floor used to estimate the residual correlation (Vhat).}
#' }
#' Random and null are selected identically for both classes
#' (\code{\link{mashRandNullSample}} over the object's \code{Bhat}/\code{Shat}).
#' Partitions are merged across objects (rownames disambiguated by region name),
#' cleaned + z-derived by \code{\link{filterInvalidSummaryStat}} (\code{btoz}),
#' and the strong \code{XtX} cross-product appended.
#'
#' @param objects A named \code{list} of \code{\link{QtlSumStats}} and/or
#'   \code{\link{FineMappingResult}} objects, one per region. Names disambiguate
#'   rownames across regions (defaults to \code{region1}, \code{region2}, ...).
#'   For \code{QtlSumStats} inputs \code{\link{summaryStatsQc}} must have been
#'   run (the matrix builder rejects un-QC'd SumStats).
#' @param nRandom,nNull Per-object random / null sample sizes (default 10 each).
#' @param excludeCondition Character vector of condition (column) names to drop.
#' @param coverage Credible-set coverage for \code{FineMappingResult} strong
#'   selection (default 0.95).
#' @param zOnly When \code{TRUE} the returned partitions carry only \code{.z}
#'   (the \code{.b}/\code{.s} matrices are dropped after z is derived).
#' @param sigPCutoff Significance cutoff applied to the strong partition
#'   (default 1e-6).
#' @param inputScale Matrix scale for \code{QtlSumStats} inputs
#'   (\code{"auto"}/\code{"beta"}/\code{"z"}); ignored for \code{FineMappingResult}
#'   (always effect-size scale).
#' @param independentVariants Optional character vector of variant ids (e.g. an
#'   LD-pruned independent SNP list). When supplied, the \emph{random} and
#'   \emph{null} background of every object is restricted to variants that match
#'   this set, so the background carries no LD-correlated SNPs (which would bias
#'   the residual correlation and the mixture weights). Matching is delegated to
#'   \code{matchVariants()} (chrom/pos/allele aware, ref/alt flips tolerated),
#'   \emph{not} a raw string compare, so a chr-prefix / separator / allele-order
#'   difference still matches. The \emph{strong} partition is never filtered.
#' @param seed RNG seed for the random / null sampling (default 999).
#'
#' @return A flat \code{list}: \code{strong.b}, \code{strong.s}, \code{strong.z},
#'   \code{random.*}, \code{null.*} (each a \code{variants x conditions} matrix)
#'   and \code{XtX} (a \code{conditions x conditions} matrix). The \code{.b} /
#'   \code{.s} matrices are omitted when \code{zOnly = TRUE}.
#' @seealso \code{\link{mashRandNullSample}}, \code{\link{mergeMashData}},
#'   \code{\link{filterInvalidSummaryStat}}
#' @export
mashInput <- function(objects, nRandom = 10L, nNull = 10L,
                      excludeCondition = character(0), coverage = 0.95,
                      zOnly = FALSE, sigPCutoff = 1e-6,
                      inputScale = c("auto", "beta", "z"),
                      independentVariants = NULL, seed = 999L) {
  inputScale <- match.arg(inputScale)
  if (!is.null(independentVariants)) {
    independentVariants <- as.character(independentVariants)
  }
  if (!is.list(objects) || length(objects) == 0L) {
    stop("mashInput: `objects` must be a non-empty list of QtlSumStats ",
         "and/or FineMappingResult objects.")
  }
  if (is.null(names(objects)) || any(!nzchar(names(objects)))) {
    names(objects) <- paste0("region", seq_along(objects))
  }

  combined <- list()
  for (nm in names(objects)) {
    part <- .mashObjectPartitions(objects[[nm]], nRandom = nRandom,
                                  nNull = nNull,
                                  excludeCondition = excludeCondition,
                                  coverage = coverage, inputScale = inputScale,
                                  seed = seed,
                                  independentVariants = independentVariants)
    # Disambiguate rownames across regions before accumulating.
    part <- lapply(part, function(m) {
      if (!is.null(m) && nrow(m) > 0L) {
        rownames(m) <- paste(rownames(m), nm, sep = "_")
      }
      m
    })
    combined <- mergeMashData(combined, part)
  }

  # A single-object run leaves matrices; coerce to data.frame so the
  # dplyr-based filterInvalidSummaryStat accepts every partition uniformly.
  combined <- lapply(combined, function(m) {
    if (is.null(m)) return(NULL) # nocov  (partitions are always matrices here, never NULL)
    as.data.frame(m)
  })

  # Clean each partition and derive z = b / s. Order random/null before strong
  # so the strong significance filter (keyed on strong.z) runs exactly once.
  for (cond in c("random", "null", "strong")) {
    bKey <- paste0(cond, ".b"); sKey <- paste0(cond, ".s")
    if (!is.null(combined[[bKey]]) && !is.null(combined[[sKey]])) {
      combined <- filterInvalidSummaryStat(combined, bhat = bKey, sbhat = sKey,
                                           btoz = TRUE, sigPCutoff = sigPCutoff)
    }
  }

  # filterInvalidSummaryStat subsets the strong partition without drop = FALSE,
  # so a single surviving strong variant degrades to a vector; restore the
  # 1-row matrix shape (conditions preserved as column names).
  for (k in c("strong.b", "strong.s", "strong.z")) {
    v <- combined[[k]]
    if (!is.null(v) && is.null(dim(v))) {
      combined[[k]] <- matrix(v, nrow = 1L, dimnames = list(NULL, names(v)))
    }
  }

  # Strong XtX cross-product (conditions x conditions).
  if (!is.null(combined$strong.z) && nrow(as.matrix(combined$strong.z)) > 0L) {
    sz <- as.matrix(combined$strong.z)
    combined$XtX <- crossprod(sz) / nrow(sz)
  }

  if (zOnly) {
    combined[grep("\\.(b|s)$", names(combined), value = TRUE)] <- NULL
  }
  combined
}

#' @title Build a QtlSumStats from a Z-score matrix
#' @description Assemble a per-condition \code{\link{QtlSumStats}} from a
#'   \code{variants x conditions} Z-score matrix -- the input shape the mash
#'   pipeline uses when only Z is available. Each column is one condition, and
#'   conditions are distinguished by \code{context}, \code{trait}, or both: the
#'   columns may be different cell types / tissues (contexts), different
#'   molecular phenotypes (traits), or arbitrary context x trait pairs.
#'   Chromosome / position are decoded from the row (variant) identifiers via
#'   \code{\link{parseVariantId}} (with a synthetic-position fallback for ids
#'   that do not encode coordinates); \code{A1} / \code{A2} / \code{N} are
#'   placeholders because a Z-only input carries no alleles or sample sizes
#'   (mash reads only Z). A pass-through \code{qcInfo} record is set so the
#'   result clears the mash QC gate.
#' @param z Numeric matrix (variants x conditions). \code{rownames(z)} are
#'   variant ids (ideally \code{chr:pos:A2:A1}); \code{colnames(z)} label the
#'   conditions.
#' @param study Study identifier (recycled across conditions).
#' @param ldSketch A \code{\link{GenotypeHandle}} embedded in the collection, or
#'   \code{NULL} (default) -- mash operates across conditions per variant and
#'   needs no LD reference.
#' @param context Condition context label(s): a single value recycled across
#'   every column, or a length-\code{ncol(z)} vector (one per condition).
#'   Defaults to \code{colnames(z)} -- one context per column.
#' @param trait Condition trait label(s): a single value recycled across every
#'   column, or a length-\code{ncol(z)} vector. Default \code{"mash"}. Pass
#'   \code{colnames(z)} here (with a constant \code{context}) when the columns
#'   are traits rather than contexts.
#' @param genome Genome build. Default \code{"GRCh38"}.
#' @param n Placeholder per-variant sample size. Default \code{1000}.
#' @param a1,a2 Placeholder alleles. Defaults \code{"A"} / \code{"G"}.
#' @param role Tag stored in the \code{qcInfo} record. Default \code{"mash"}.
#' @return A \code{\link{QtlSumStats}} with one entry per condition (column).
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
qtlSumStatsFromZMatrix <- function(z, study, ldSketch = NULL,
                                   context = colnames(z), trait = "mash",
                                   genome = "GRCh38", n = 1000L,
                                   a1 = "A", a2 = "G", role = "mash") {
  if (!is.matrix(z) || !is.numeric(z)) {
    stop("qtlSumStatsFromZMatrix: `z` must be a numeric variants x conditions matrix.")
  }
  vids <- rownames(z)
  if (is.null(vids)) vids <- paste0("var", seq_len(nrow(z)))
  .qtlSumStatsFromMatrix(
    vids = vids, nCond = ncol(z), study = study, ldSketch = ldSketch,
    context = context, trait = trait, genome = genome, role = role,
    mcolFn = function(j) S4Vectors::DataFrame(
      SNP = vids, A1 = rep(a1, length(vids)), A2 = rep(a2, length(vids)),
      Z = as.numeric(z[, j]), N = rep(as.integer(n), length(vids))))
}

#' @title Build a QtlSumStats from Bhat / Shat (effect-size) Matrices
#' @description Assemble a per-condition \code{\link{QtlSumStats}} from an aligned
#'   pair of \code{variants x conditions} effect-size (\code{Bhat}) and
#'   standard-error (\code{Shat}) matrices -- the beta-scale (EE) counterpart of
#'   \code{\link{qtlSumStatsFromZMatrix}}. Each entry carries \code{BETA},
#'   \code{SE}, and (derived) \code{Z = BETA / SE} mcols, so the result feeds
#'   \code{\link{mashPipeline}} / \code{\link{mashModelFit}} on either scale
#'   (\code{inputScale = "beta"} or \code{"z"}). Chromosome / position are decoded
#'   from the row (variant) ids exactly as in \code{\link{qtlSumStatsFromZMatrix}}.
#' @param bhat Numeric matrix (variants x conditions) of effect sizes.
#'   \code{rownames(bhat)} are variant ids; \code{colnames(bhat)} label the
#'   conditions.
#' @param shat Numeric matrix of standard errors, aligned with \code{bhat}
#'   (identical dimensions and row/column order).
#' @param study Study identifier (recycled across conditions).
#' @param ldSketch A \code{\link{GenotypeHandle}} embedded in the collection, or
#'   \code{NULL} (default) -- mash operates across conditions per variant and
#'   needs no LD reference.
#' @param context,trait Condition labels; see \code{\link{qtlSumStatsFromZMatrix}}.
#'   Defaults \code{context = colnames(bhat)}, \code{trait = "mash"}.
#' @param genome Genome build. Default \code{"GRCh38"}.
#' @param n Placeholder per-variant sample size. Default \code{1000}.
#' @param a1,a2 Placeholder alleles. Defaults \code{"A"} / \code{"G"}.
#' @param role Tag stored in the \code{qcInfo} record. Default \code{"mash"}.
#' @return A \code{\link{QtlSumStats}} with one entry per condition (column).
#' @seealso \code{\link{qtlSumStatsFromZMatrix}}, \code{\link{mashModelFit}}
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
qtlSumStatsFromBetaMatrix <- function(bhat, shat, study, ldSketch = NULL,
                                      context = colnames(bhat), trait = "mash",
                                      genome = "GRCh38", n = 1000L,
                                      a1 = "A", a2 = "G", role = "mash") {
  if (!is.matrix(bhat) || !is.numeric(bhat)) {
    stop("qtlSumStatsFromBetaMatrix: `bhat` must be a numeric variants x conditions matrix.")
  }
  if (!is.matrix(shat) || !is.numeric(shat)) {
    stop("qtlSumStatsFromBetaMatrix: `shat` must be a numeric variants x conditions matrix.")
  }
  if (!identical(dim(bhat), dim(shat))) {
    stop(sprintf(paste0("qtlSumStatsFromBetaMatrix: `bhat` (%s) and `shat` (%s) ",
                        "must have identical dimensions."),
                 paste(dim(bhat), collapse = "x"), paste(dim(shat), collapse = "x")))
  }
  vids <- rownames(bhat)
  if (is.null(vids)) vids <- paste0("var", seq_len(nrow(bhat)))
  .qtlSumStatsFromMatrix(
    vids = vids, nCond = ncol(bhat), study = study, ldSketch = ldSketch,
    context = context, trait = trait, genome = genome, role = role,
    mcolFn = function(j) S4Vectors::DataFrame(
      SNP = vids, A1 = rep(a1, length(vids)), A2 = rep(a2, length(vids)),
      BETA = as.numeric(bhat[, j]), SE = as.numeric(shat[, j]),
      Z = as.numeric(bhat[, j] / shat[, j]),
      N = rep(as.integer(n), length(vids))))
}

# Internal: shared assembly for the z / beta matrix constructors. Recycles the
# context / trait labels, decodes chrom/pos from the variant ids (synthesising
# where they don't parse), builds one GRanges entry per condition with mcols
# from `mcolFn(j)`, and wraps the entries as a QtlSumStats.
# @noRd
.qtlSumStatsFromMatrix <- function(vids, nCond, study, ldSketch, context, trait,
                                   genome, role, mcolFn) {
  context <- .qszmRecycle(context, nCond, "context")
  trait   <- .qszmRecycle(trait,   nCond, "trait")
  # Decode chrom/pos from the variant ids; synthesise where they do not parse.
  parsed <- tryCatch(suppressWarnings(parseVariantId(vids)),
                     error = function(e) NULL)
  chrom <- if (!is.null(parsed)) as.character(parsed$chrom)
           else rep(NA_character_, length(vids))
  pos <- if (!is.null(parsed)) suppressWarnings(as.integer(parsed$pos))
         else rep(NA_integer_, length(vids))
  chrom[is.na(chrom) | !nzchar(chrom)] <- "chr1"
  pos[is.na(pos)] <- seq_along(pos)[is.na(pos)]
  entries <- lapply(seq_len(nCond), function(j) {
    gr <- GenomicRanges::GRanges(
      seqnames = chrom, ranges = IRanges::IRanges(start = pos, width = 1L))
    S4Vectors::mcols(gr) <- mcolFn(j)
    gr
  })
  QtlSumStats(
    study    = rep(as.character(study), nCond),
    context  = context,
    trait    = trait,
    entry    = entries,
    genome   = genome,
    ldSketch = ldSketch,
    qcInfo   = list(role = role, entryAudit = vector("list", nCond)))
}

# Internal: recycle a condition-label argument (context / trait) to one value
# per matrix column. Accepts a single value (recycled to every column) or a
# vector of length ncol; errors otherwise -- including NULL, which is how
# `context = colnames(x)` arrives when the matrix has no column names.
.qszmRecycle <- function(v, n, what) {
  if (is.null(v)) {
    stop(sprintf(paste0("qtlSumStats matrix constructor: `%s` is NULL; pass a ",
                        "length-1 or length-%d value (or give the matrix column ",
                        "names)."), what, n))
  }
  v <- as.character(v)
  if (length(v) == 1L) return(rep(v, n))
  if (length(v) == n) return(v)
  stop(sprintf(paste0("qtlSumStats matrix constructor: `%s` must be length 1 or ",
                      "ncol=%d, got %d."), what, n, length(v)))
}

# Internal: convert a single SumStats object (post-QC) into a (Bhat, Shat)
# pair of matrices keyed by context. Each row of the matrix corresponds to
# one (variantId × (study, trait)) cell from the SumStats entries; each
# column corresponds to a context (from QtlSumStats $context; from
# GwasSumStats $study, which is the per-study mash column).
#
# For QtlSumStats:
#   * Pivots entries on (study, trait) so each (study, trait) becomes a
#     block of rows and each context becomes a column. Missing
#     (study, trait, context) cells are filled with NA.
# For GwasSumStats:
#   * Each row of the collection is one study; we treat each study as a
#     mash "context" (single block of rows per study, columns = studies).
#     This is rarely used on its own but lets a flat GwasSumStats pass
#     through alongside (or instead of) a QtlSumStats without special
#     casing further upstream.
#
# Variant alignment within a (study, trait) block uses the entry's
# variant order; missing variants in any one context are filled with NA.
# NA in Bhat is mapped to 0 and NA in Shat is mapped to a large value
# (1000) inside mashr::mash_set_data via its `zero_Bhat_Shat_reset`
# pathway, matching the prior pipeline's handling of incomplete cells.
# @noRd
.mashSumStatsToMatrices <- function(x, role,
                                    inputScale = c("auto", "beta", "z")) {
  inputScale <- match.arg(inputScale)
  if (!methods::is(x, "QtlSumStats") && !methods::is(x, "GwasSumStats")) {
    stop(sprintf(
      "mashPipeline: '%s' input must be a QtlSumStats or GwasSumStats; got %s.",
      role, paste(class(x), collapse = "/")))
  }
  if (length(getQcInfo(x)) == 0L) {
    stop(sprintf(
      "mashPipeline: '%s' SumStats has no QC info (length(getQcInfo(x)) == 0L). ",
      role),
      "Run summaryStatsQc() on the SumStats before passing it to mashPipeline().")
  }
  if (nrow(x) == 0L) {
    stop(sprintf(
      "mashPipeline: '%s' SumStats has no entries (nrow == 0).", role))
  }

  isQtl <- methods::is(x, "QtlSumStats")
  if (isQtl) {
    studyCol <- as.character(x$study)
    traitCol <- as.character(x$trait)
    contextCol <- as.character(x$context)
    blockKeys <- paste(studyCol, traitCol, sep = "::")
    columnLabels <- unique(contextCol)
  } else {
    studyCol <- as.character(x$study)
    blockKeys <- studyCol
    contextCol <- studyCol
    columnLabels <- unique(studyCol)
  }

  # Per-entry GRanges (avoid `@`; use the public list-column accessor).
  entries <- x$entry

  # Resolve per-entry scale: which (Bhat, Shat) source to pull. mashr
  # expects one coherent convention per call:
  #   "beta" → Bhat = BETA, Shat = SE   (effect-size scale; standard)
  #   "z"    → Bhat = Z,    Shat = 1    (z-score scale)
  # "auto" picks "beta" when every entry has BETA + SE, else "z" if
  # every entry has Z. Mixed inputs (some entries missing BETA, others
  # missing Z) are a hard error.
  entryCaps <- lapply(seq_len(nrow(x)), function(i) {
    mc <- S4Vectors::mcols(entries[[i]])
    list(hasBetaSe = all(c("BETA", "SE") %in% colnames(mc)),
         hasZ      = "Z" %in% colnames(mc))
  })
  allHaveBetaSe <- all(vapply(entryCaps, `[[`, logical(1), "hasBetaSe"))
  allHaveZ      <- all(vapply(entryCaps, `[[`, logical(1), "hasZ"))
  resolvedScale <- switch(inputScale,
    beta = {
      if (!allHaveBetaSe)
        stop(sprintf(
          "mashPipeline: inputScale = 'beta' requires every '%s' entry to ",
          role),
          "carry both BETA and SE mcols.")
      "beta"
    },
    z = {
      if (!allHaveZ)
        stop(sprintf(
          "mashPipeline: inputScale = 'z' requires every '%s' entry to ",
          role), "carry a Z mcol.")
      "z"
    },
    auto = {
      if (allHaveBetaSe) "beta"
      else if (allHaveZ) "z"
      else stop(sprintf(
        "mashPipeline: '%s' SumStats has no usable scale — every entry ",
        role),
        "must carry (BETA, SE) or Z mcols.")
    })

  # Group rows of x by (study, trait) block; within each block, build a
  # variant × context matrix for Bhat and Shat.
  uniqBlocks <- unique(blockKeys)
  bhatBlocks <- vector("list", length(uniqBlocks))
  shatBlocks <- vector("list", length(uniqBlocks))
  variantRowNames <- vector("list", length(uniqBlocks))

  for (bi in seq_along(uniqBlocks)) {
    bkey <- uniqBlocks[[bi]]
    rowsInBlock <- which(blockKeys == bkey)

    # Variant universe for this block = union of variant IDs (SNP)
    # across the contexts in this block, preserving first-seen order.
    variantOrder <- character()
    perContextB  <- list()
    perContextSe <- list()
    requireCols <- if (resolvedScale == "beta") c("SNP", "BETA", "SE")
                   else                          c("SNP", "Z")
    for (rIdx in rowsInBlock) {
      df <- if (isQtl) {
        getSumstatDf(x,
                     study   = studyCol[[rIdx]],
                     context = contextCol[[rIdx]],
                     trait   = traitCol[[rIdx]],
                     require = requireCols)
      } else {
        getSumstatDf(x, study = studyCol[[rIdx]], require = requireCols)
      }
      snps <- df$variant_id
      newSnps <- setdiff(snps, variantOrder)
      variantOrder <- c(variantOrder, newSnps)
      ctx <- contextCol[[rIdx]]
      if (resolvedScale == "beta") {
        perContextB[[ctx]]  <- setNames(df$beta, snps)
        perContextSe[[ctx]] <- setNames(df$se,   snps)
      } else {
        perContextB[[ctx]]  <- setNames(df$z, snps)
        perContextSe[[ctx]] <- setNames(rep(1, length(snps)), snps)
      }
    }

    nVar <- length(variantOrder)
    bMat <- matrix(NA_real_, nrow = nVar, ncol = length(columnLabels),
                   dimnames = list(variantOrder, columnLabels))
    sMat <- matrix(NA_real_, nrow = nVar, ncol = length(columnLabels),
                   dimnames = list(variantOrder, columnLabels))
    for (ctx in names(perContextB)) {
      bMat[names(perContextB[[ctx]]), ctx] <- perContextB[[ctx]]
      sMat[names(perContextSe[[ctx]]), ctx] <- perContextSe[[ctx]]
    }
    # Disambiguate rownames across blocks to avoid silent dedup.
    rownames(bMat) <- paste(bkey, variantOrder, sep = "::")
    rownames(sMat) <- rownames(bMat)
    bhatBlocks[[bi]] <- bMat
    shatBlocks[[bi]] <- sMat
    variantRowNames[[bi]] <- rownames(bMat)
  }

  bhat <- do.call(rbind, bhatBlocks)
  shat <- do.call(rbind, shatBlocks)

  # Replace NAs: bhat NA -> 0, shat NA -> 1000 (the same convention as
  # filterInvalidSummaryStat() and mash_set_data()'s
  # zero_Bhat_Shat_reset; ensures missing-cell variants do not drive the
  # fit).
  bhat[is.na(bhat)] <- 0
  shat[is.na(shat) | shat <= 0] <- 1000

  list(b = bhat, s = shat)
}

#' Estimate mash covariance matrices and mixture weights from SumStats.
#'
#' Genome-wide pipeline that estimates the \pkg{mashr} canonical, PCA,
#' flash, and ED covariance matrices plus mixture weights between
#' contexts. The pipeline is memory-intensive and not gene-parallelizable
#' (see \code{dev/refactor-design.md}).
#'
#' @param sumStatsList A named \code{list} (or \code{S4Vectors::SimpleList})
#'   of \code{\link{QtlSumStats}} and/or \code{\link{GwasSumStats}}
#'   objects. The list MUST be named with at least \code{"strong"} and
#'   \code{"random"}; \code{"null"} is optional. Each element is one
#'   SumStats collection whose entries are pivoted internally into a
#'   variants \eqn{\times} contexts \eqn{Bhat} / \eqn{Shat} matrix pair.
#'   For \code{QtlSumStats} the columns of the resulting matrix are the
#'   \code{context} values; for \code{GwasSumStats} the columns are the
#'   \code{study} values. Every SumStats must have been processed by
#'   \code{\link{summaryStatsQc}} (the pipeline rejects inputs where
#'   \code{length(getQcInfo(x)) == 0L}).
#' @param alpha mash \code{alpha} parameter (passed to
#'   \code{mashr::mash_set_data()}).
#' @param residualCorrelation Optional residual correlation matrix. Used
#'   in place of the null-data-derived V matrix when no \code{"null"}
#'   entry is supplied (matches the prior pipeline).
#' @param nPcs Number of principal components for PCA-based covariance
#'   matrices. Defaults to \code{ncol(Bhat) - 1}.
#' @param setSeed Integer seed for reproducibility (default 999).
#'
#' @return A list with two elements: \code{U} (the combined list of mash
#'   covariance matrices: canonical + PCA + flash + ED) and \code{w}
#'   (the estimated mixture weights from \code{mashr::get_estimated_pi}).
#'
#' @section Behavioural notes (changes vs. the legacy contract):
#' The legacy \code{mashInput} list of six pre-built matrices
#' (\code{strong.b}/\code{strong.s}/\code{random.b}/\code{random.s}/
#' \code{null.b}/\code{null.s}) is replaced by a \code{list} of SumStats
#' objects. Construction of the per-partition matrices now happens
#' inside \code{mashPipeline} via \code{getSumStats()} /
#' \code{S4Vectors::mcols()} accessors — no \code{@@} slot access. The
#' mashr algorithm itself is unchanged: same \code{cov_canonical},
#' \code{cov_pca}, \code{cov_flash}, \code{cov_ed}, and
#' \code{mash(..., outputlevel = 1)} sequence as before.
#'
#' @export
