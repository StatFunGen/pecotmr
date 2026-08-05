# =============================================================================
# TwasWeights S4 class
# -----------------------------------------------------------------------------
# DFrame-subclass collection keyed by the identity tuple (study, context,
# trait, method). Each row holds a TwasWeightsEntry payload (variant ids
# + per-variant weight vector/matrix). Class-level slots:
#   * ldSketch   GenotypeHandle (NULL for individual-level fits, the
#                LD-sketch handle for RSS-derived weights).
# Constructor + accessors below. The twasWeights pipeline helpers
# (learnTwasWeights, CV, ensemble, etc.) follow at the bottom.
# =============================================================================

#' @include AllGenerics.R tupleSelectors.R
NULL

setClass("TwasWeights",
  contains = "DFrame",
  representation(ldSketch = "ANY"),
  validity = function(object) {
    errors <- character()
    required <- c("study", "context", "trait", "method", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L)
      errors <- c(errors, paste("missing columns:",
                                paste(missingCols, collapse = ", ")))
    if (length(errors) == 0L) {
      if (length(object$entry) != nrow(object))
        errors <- c(errors,
          "length(entry) must equal nrow(.) for TwasWeights")
      entryTypes <- vapply(object$entry,
                          function(e) methods::is(e, "TwasWeightsEntry"),
                          logical(1))
      if (!all(entryTypes))
        errors <- c(errors,
          "every element of the `entry` column must be a TwasWeightsEntry")
      # Optional `region` provenance (one genomic anchor per row; non-key).
      errors <- c(errors, .validateRegionColumn(object))
      errors <- c(errors, .validateTraitPosColumn(object))
      jointCols <- intersect(
        c("jointStudies", "jointContexts", "jointTraits"), names(object))
      for (jc in jointCols) {
        vals <- object[[jc]]
        if (!is.character(vals))
          errors <- c(errors, sprintf(
            "'%s' column must be character (got %s)", jc, class(vals)[[1L]]))
      }
      keyCols <- c("study", "context", "trait", "method", jointCols)
      # Extract key columns directly rather than via `object[, keyCols]`:
      # column-subsetting preserves the TwasWeights class while dropping the
      # required `entry` column, and older S4Vectors revalidates that
      # intermediate, spuriously failing with "missing columns: entry".
      keyDf <- as.data.frame(
        lapply(keyCols, function(cn) object[[cn]]),
        col.names = keyCols, stringsAsFactors = FALSE)
      if (anyDuplicated(keyDf))
        errors <- c(errors,
          "(study, context, trait, method[, joint*]) tuple uniqueness violated")
    }
    if (!is.null(object@ldSketch) &&
        !methods::is(object@ldSketch, "GenotypeHandle")) {
      errors <- c(errors,
        "'ldSketch' must be a GenotypeHandle or NULL")
    }
    if (length(errors) == 0L) TRUE else errors
  }
)


# =============================================================================
# QTL Dataset
# =============================================================================

#' @title QTL Dataset (individual-level data for one study)
#' @description S4 container for a single QTL study's regional data. Holds
#'   a genotype handle plus per-context \code{SummarizedExperiment} objects
#'   carrying molecular-trait measurements. Each context's SE has
#'   \code{rowRanges} describing per-trait genomic positions and
#'   \code{colData} carrying per-context phenotype covariates. A single
#'   matrix of genotype-derived covariates (e.g., ancestry PCs) applies
#'   across contexts.
#'
#' @slot study Character (length 1). Study identifier; used in collection
#'   classes to tag downstream \code{FineMappingResult} / \code{TwasWeights}
#'   entries.
#' @slot genotypes A \code{GenotypeHandle} for lazy access to genotype
#'   dosages.
#' @slot phenotypes Named list of \code{SummarizedExperiment} objects, one
#'   per QTL context. Each SE has rows = molecular traits with positions
#'   in \code{rowRanges(se)}, columns = samples, and per-context covariates
#'   in \code{colData(se)}. Different contexts may carry different subsets
#'   of traits (rows); traits shared across contexts must have identical
#'   \code{rowRanges} entries (enforced by validity).
#' @slot genotypeCovariates Numeric matrix (samples x covariates) of
#'   genotype-derived covariates applied uniformly across all contexts
#'   (e.g., ancestry PCs).
#' @slot scaleResiduals Logical (length 1). Whether residualization
#'   accessors scale residuals to unit variance.
#' @slot mafCutoff Numeric (length 1). Minor allele frequency threshold;
#'   variants with \code{MAF < mafCutoff} are dropped at extraction time
#'   inside \code{getGenotypes()} / \code{getResidualizedGenotypes()}.
#'   Default 0 (no filter).
#' @slot macCutoff Numeric (length 1). Minor allele count threshold;
#'   converted to a MAF threshold using
#'   \code{max(mafCutoff, macCutoff / (2 * n))} where \code{n} is the
#'   post-narrowing sample count of the extracted block. Default 0
#'   (no filter).
#' @slot xvarCutoff Numeric (length 1). Per-variant genotype variance
#'   threshold; variants with column variance below this are dropped at
#'   extraction time. Default 0 (no filter).
#' @slot imissCutoff Numeric (length 1). Per-sample genotype-missingness
#'   threshold; samples with a missing-genotype rate above this are
#'   dropped at extraction time. Default 0 (no filter).
#' @slot keepSamples Character vector of sample identifiers to retain
#'   prior to per-block QC; intersected with the genotype handle's
#'   \code{sampleIds} and the \code{samples} argument of
#'   \code{getGenotypes()}. Length 0 means no restriction.
#' @slot keepVariants Character vector of variant identifiers to retain
#'   prior to per-block QC. Length 0 means no restriction.
#' @export

# =============================================================================

#' @title Create a TwasWeights Collection Object
#' @description Construct a \code{TwasWeights} DFrame-subclass collection
#'   from per-tuple vectors and a list of \code{TwasWeightsEntry}
#'   payloads (one per tuple).
#' @param study Character vector of study identifiers. Use the sentinel
#'   \code{"joint"} for rows produced by a cross-study joint fit.
#' @param context Character vector of context labels. Use \code{"joint"}
#'   for rows produced by a cross-context joint fit.
#' @param trait Character vector of trait identifiers. Use \code{"joint"}
#'   for rows produced by a cross-trait joint fit.
#' @param method Character vector of TWAS weight method names.
#' @param entry List / \code{SimpleList} of \code{TwasWeightsEntry} objects.
#' @param jointStudies Optional character vector (length \code{length(study)})
#'   listing the semicolon-joined studies participating in each row's
#'   cross-study joint fit, or \code{NA_character_} for non-joint rows.
#'   When \code{NULL} (default) the column is omitted.
#' @param jointContexts Optional character vector for cross-context joints.
#'   Same shape as \code{jointStudies}.
#' @param jointTraits Optional character vector for cross-trait joints.
#'   Same shape as \code{jointStudies}.
#' @param region Optional \code{GRanges} (length \code{length(study)}) giving
#'   the genomic anchor of each row's trait -- the trait's own coordinates when
#'   built from a \code{QtlDataset}, or the summary-stat variant-span window
#'   when built from a \code{QtlSumStats}. Carried forward as provenance (e.g.
#'   for cTWAS LD-block placement); not part of the identity key. \code{NULL}
#'   (default) omits the column.
#' @param ldSketch An optional \code{GenotypeHandle}, or \code{NULL} for
#'   individual-level fits.
#' @return A \code{TwasWeights} object.
#' @export
TwasWeights <- function(study, context, trait, method, entry,
                        jointStudies = NULL,
                        jointContexts = NULL,
                        jointTraits = NULL,
                        region = NULL,
                        traitPos = NULL,
                        ldSketch = NULL) {
  n <- length(study)
  if (length(context) != n || length(trait) != n || length(method) != n ||
      length(entry) != n) {
    stop("`study`, `context`, `trait`, `method`, and `entry` must all ",
         "have the same length.")
  }
  cols <- list(
    study   = as.character(study),
    context = as.character(context),
    trait   = as.character(trait),
    method  = as.character(method),
    entry   = S4Vectors::SimpleList(entry)
  )
  for (pair in list(c("jointStudies", "jointStudies"),
                    c("jointContexts", "jointContexts"),
                    c("jointTraits", "jointTraits"))) {
    val <- get(pair[[1L]])
    if (is.null(val)) next
    if (length(val) != n)
      stop("`", pair[[1L]], "` must have the same length as `study`.")
    cols[[pair[[2L]]]] <- as.character(val)
  }
  cols <- .appendRegionCol(cols, region, n)
  cols <- .appendTraitPosCol(cols, traitPos, n)
  df <- do.call(S4Vectors::DataFrame,
                c(cols, list(check.names = FALSE)))
  obj <- new("TwasWeights", df, ldSketch = ldSketch)
  validObject(obj)
  obj
}

#' @rdname getRegion
#' @export
setMethod("getRegion", "TwasWeights", function(x, ...) .getRegionColumn(x))

#' @rdname getTraitPosition
#' @export
setMethod("getTraitPosition", "TwasWeights", function(x, ...) .getTraitPosColumn(x))

#' @title Get a Single TWAS Weights Entry
#' @description Return the \code{TwasWeightsEntry} for one
#'   \code{(study, context, trait, method)} row of a \code{TwasWeights}
#'   collection.
#' @param x A \code{TwasWeights} object.
#' @param study,context,trait,method Single character identifiers. All
#'   required when the collection has more than one row; optional when
#'   the collection has a single row.
#' @return A \code{TwasWeightsEntry} object.
#' @export
setGeneric("getTwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL)
    standardGeneric("getTwasWeights"))

#' @rdname getTwasWeights
#' @export
setMethod("getTwasWeights", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    idx <- .tupleSelectRow(x, study, context, trait, method,
                           cls = "TwasWeights")
    x$entry[[idx]]
  })

#' @rdname getWeights
#' @export
setMethod("getWeights", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getWeights(entry)
  })

#' @rdname getCvResult
#' @export
setMethod("getCvResult", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getCvResult(entry)
  })

#' @rdname getFits
#' @export
setMethod("getFits", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getFits(entry)
  })

#' @rdname getStandardized
#' @export
setMethod("getStandardized", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getStandardized(entry)
  })

#' @rdname getDataType
#' @export
setMethod("getDataType", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getDataType(entry)
  })

#' @rdname getVariantIds
#' @export
setMethod("getVariantIds", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getVariantIds(entry)
  })

#' @rdname getStudy
#' @export
setMethod("getStudy", "TwasWeights",
          function(x) unique(as.character(x$study)))

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "TwasWeights",
          function(x, ...) x@ldSketch)

#' @rdname getContexts
#' @export
setMethod("getContexts", "TwasWeights",
          function(x) unique(as.character(x$context)))

#' @rdname getTraits
#' @export
setMethod("getTraits", "TwasWeights",
          function(x) unique(as.character(x$trait)))

#' @rdname getMethodNames
#' @export
setMethod("getMethodNames", "TwasWeights",
          function(x) unique(as.character(x$method)))



#' @export
setMethod("show", "TwasWeights", function(object) {
  cat(sprintf("TwasWeights: %d entries\n", nrow(object)))
  if (nrow(object) > 0L) {
    cat(sprintf("  %d studies, %d contexts, %d traits, %d methods\n",
                length(unique(object$study)),
                length(unique(object$context)),
                length(unique(object$trait)),
                length(unique(object$method))))
  }
  ldSrc <- if (is.null(object@ldSketch)) "NULL (individual-level fit)"
           else sprintf("%s @ %s",
                         object@ldSketch@format,
                         object@ldSketch@path)
  cat(sprintf("  LD sketch: %s\n", ldSrc))
})



# =============================================================================
# TwasWeights pipeline helpers (learnTwasWeights + CV + ensemble + Mvsusie/Mrmash)
# =============================================================================

# Evaluate an expression while suppressing external package output.
# Catches both message() output (susieR, qgg) and Rprintf/cat stdout (mr.ash.alpha).
# @param expr An expression to evaluate.
# @return The result of evaluating expr.
# @noRd
.quietEval <- function(expr) {
  invisible(capture.output(
    result <- suppressMessages(expr),
    type = "output"
  ))
  result
}

# Rename a "_weights"/"Weights" suffix to the case-matching equivalent of `target`.
# Snake-case inputs get the underscored snake-case form (e.g. "lasso_weights" -> "lasso_predicted")
# and camelCase inputs get the CamelCase form (e.g. "lassoWeights" -> "lassoPredicted").
# Names without a recognized suffix are returned unchanged.
# @param x Character vector of names ending in "_weights" or "Weights".
# @param target A bare token such as "predicted" or "performance".
# @return Character vector with suffixes rewritten.
# @noRd
.renameSuffix <- function(x, target) {
  cap <- paste0(toupper(substr(target, 1, 1)), substr(target, 2, nchar(target)))
  x <- sub("_weights$", paste0("_", target), x)
  x <- sub("Weights$", cap, x)
  x
}

# Map short method names and presets to weightMethods lists.
# @param methods A character vector of short method names, or a preset string
#   ("default" or "fast_default").
# @return A named list suitable for the weightMethods parameter.
# @noRd
.twasMethodLookup <- function(methods) {
  # `fn` is the snake_case key used in weight method lists; `impl` is the
  # actual camelCase function name implemented by the package.
  methodMap <- list(
    susie = list(fn = "susie_weights", impl = "susieWeights", args = list(refine = FALSE, L = 20, L_greedy = 5)),
    susieAsh = list(fn = "susie_ash_weights", impl = "susieAshWeights", args = list()),
    susieInf = list(fn = "susie_inf_weights", impl = "susieInfWeights", args = list()),
    mrash = list(fn = "mrash_weights", impl = "mrashWeights", args = list(initPriorSd = TRUE, max.iter = 100)),
    enet = list(fn = "enet_weights", impl = "enetWeights", args = list()),
    lasso = list(fn = "lasso_weights", impl = "lassoWeights", args = list()),
    bayes_r = list(fn = "bayes_r_weights", impl = "bayesRWeights", args = list()),
    bayes_l = list(fn = "bayes_l_weights", impl = "bLassoWeights", args = list()),
    bayes_a = list(fn = "bayes_a_weights", impl = "bayesAWeights", args = list()),
    bayes_b = list(fn = "bayes_b_weights", impl = "bayesBWeights", args = list()),
    bayes_c = list(fn = "bayes_c_weights", impl = "bayesCWeights", args = list()),
    bayes_n = list(fn = "bayes_n_weights", impl = "bayesNWeights", args = list()),
    b_lasso = list(fn = "b_lasso_weights", impl = "bLassoWeights", args = list()),
    dpr_vb = list(fn = "dpr_vb_weights", impl = "dprVbWeights", args = list()),
    dpr_gibbs = list(fn = "dpr_gibbs_weights", impl = "dprGibbsWeights", args = list()),
    dpr_adaptive_gibbs = list(fn = "dpr_adaptive_gibbs_weights", impl = "dprAdaptiveGibbsWeights", args = list()),
    scad = list(fn = "scad_weights", impl = "scadWeights", args = list()),
    mcp = list(fn = "mcp_weights", impl = "mcpWeights", args = list()),
    l0learn = list(fn = "l0learn_weights", impl = "l0learnWeights", args = list()),
    mvsusie = list(fn = "mvsusie_weights", impl = "mvsusieWeights", args = list(L = 30, L_greedy = 5)),
    mrmash = list(fn = "mrmash_weights", impl = "mrmashWeights", args = list(canonicalPriorMatrices = TRUE)),
    fsusie = list(fn = "fsusie_weights", impl = "fsusieWeights", args = list())
  )

  # Handle presets
  fastDefault <- c("susie", "susieInf", "mrash", "enet", "lasso", "mcp", "scad", "l0learn")
  if (length(methods) == 1) {
    if (methods == "fast_default") {
      methods <- fastDefault
    } else if (methods == "default") {
      methods <- c(fastDefault, "bayes_r", "bayes_c")
    }
  }

  # Build reverse map: function name -> short name, so full names are accepted too
  fnToShort <- setNames(
    names(methodMap),
    vapply(methodMap, function(x) x$fn, character(1))
  )
  # Normalize any full function names to short names
  methods <- vapply(methods, function(m) {
    if (m %in% names(fnToShort)) fnToShort[[m]] else m
  }, character(1), USE.NAMES = FALSE)

  unknown <- setdiff(methods, names(methodMap))
  if (length(unknown) > 0) {
    stop(
      "Unknown TWAS method(s): ", paste(unknown, collapse = ", "),
      ". Available methods: ", paste(names(methodMap), collapse = ", ")
    )
  }

  result <- list()
  for (m in methods) {
    entry <- methodMap[[m]]
    args <- entry$args
    # Track the actual function implementation name so downstream dispatchers
    # can resolve snake_case keys to the camelCase implementation.
    attr(args, "impl") <- entry$impl
    result[[entry$fn]] <- args
  }
  result
}

# TRUE if `name` is a function visible in the search path or in namespace `ns`.
# @noRd
.functionExistsInNs <- function(name, ns) {
  exists(name, mode = "function") ||
    exists(name, mode = "function", envir = ns, inherits = FALSE)
}

# Resolve the actual function name for a method key. Honors an "impl" attribute
# on the per-method args list (set by .twasMethodLookup), and otherwise applies
# a snake_case -> camelCase transformation as a fallback for user-supplied
# weightMethods lists.
.resolveMethodFunction <- function(methodKey, methodArgs = NULL) {
  # Search pecotmr's namespace explicitly so this works equally well when the
  # function is called either from inside the package or from a user session.
  ns <- asNamespace("pecotmr")
  impl <- if (!is.null(methodArgs)) attr(methodArgs, "impl") else NULL
  if (!is.null(impl) && nzchar(impl) && .functionExistsInNs(impl, ns)) {
    return(impl)
  }
  # Direct match (e.g. caller already passed camelCase)
  if (.functionExistsInNs(methodKey, ns)) return(methodKey)
  # snake_case_weights -> camelCaseWeights
  parts <- strsplit(methodKey, "_", fixed = TRUE)[[1]]
  capRest <- paste0(toupper(substring(parts[-1], 1, 1)),
                    substring(parts[-1], 2))
  candidate <- paste0(parts[1], paste0(capRest, collapse = ""))
  if (.functionExistsInNs(candidate, ns)) return(candidate)
  methodKey
}

# Validate a Sample/Fold partition data frame: required columns, no sample in
# two folds, and (when `sampleNames` given) exact coverage of all samples.
# @noRd
.validateFoldPartition <- function(df, sampleNames) {
  if (!all(c("Sample", "Fold") %in% names(df)))
    stop("samplePartition must have columns `Sample` and `Fold`.")
  df$Sample <- as.character(df$Sample)
  dup <- unique(df$Sample[duplicated(df$Sample)])
  if (length(dup) > 0L)
    stop("Fold partition assigns sample(s) to more than one fold: ",
         paste(dup, collapse = ", "))
  if (!is.null(sampleNames)) {
    unknown <- setdiff(df$Sample, sampleNames)
    if (length(unknown) > 0L)
      stop("Fold partition references unknown sample(s): ",
           paste(unknown, collapse = ", "))
    uncovered <- setdiff(sampleNames, df$Sample)
    if (length(uncovered) > 0L)
      stop("Fold partition does not cover ", length(uncovered),
           " sample(s) (folds must partition all samples), e.g. ",
           paste(utils::head(uncovered, 5L), collapse = ", "))
  }
  df
}

# Normalize a cross-validation fold specification into the canonical
# samplePartition data.frame(Sample, Fold) used throughout the CV machinery, so
# callers can pass folds in any of three forms and downstream code has a single
# source of truth. Accepts:
#   * cvFolds an integer k: k-fold auto-partition (returned as NULL so the
#     partition is generated per (study, context, trait) downstream).
#   * cvFolds a list of vectors: each element defines one fold's SAMPLES, as
#     numeric column indices into `sampleNames` or character sample names.
#   * samplePartition a data.frame(Sample, Fold): used as-is.
# A list-form `cvFolds` and an explicit `samplePartition` are mutually exclusive.
# When `sampleNames` is supplied, the resolved partition is validated to
# reference only known samples, assign each sample to one fold, and cover every
# sample (a proper partition). Returns list(samplePartition, nFolds).
.normalizeCvFolds <- function(cvFolds = 0, samplePartition = NULL,
                              sampleNames = NULL) {
  isListFolds <- is.list(cvFolds) && !is.data.frame(cvFolds)
  if (isListFolds && !is.null(samplePartition))
    stop("Provide either a list-form `cvFolds` or an explicit ",
         "`samplePartition`, not both.")


  if (!is.null(samplePartition)) {
    df <- .validateFoldPartition(as.data.frame(samplePartition,
                                          stringsAsFactors = FALSE), sampleNames)
    return(list(samplePartition = df, nFolds = length(unique(df$Fold))))
  }

  if (isListFolds) {
    if (length(cvFolds) < 2L)
      stop("A list-form `cvFolds` must define at least 2 folds.")
    rows <- lapply(seq_along(cvFolds), function(k) {
      ids <- cvFolds[[k]]
      if (is.numeric(ids)) {
        if (is.null(sampleNames))
          stop("Numeric fold vectors require `sampleNames` to resolve ",
               "column indices.")
        if (any(ids < 1L | ids > length(sampleNames)))
          stop("Fold ", k, " has out-of-range sample column index/indices.")
        ids <- sampleNames[as.integer(ids)]
      } else {
        ids <- as.character(ids)
      }
      data.frame(Sample = ids, Fold = k, stringsAsFactors = FALSE)
    })
    df <- .validateFoldPartition(do.call(rbind, rows), sampleNames)
    return(list(samplePartition = df, nFolds = length(cvFolds)))
  }

  k <- suppressWarnings(as.integer(cvFolds))
  if (length(k) != 1L || is.na(k))
    stop("`cvFolds` must be a single integer, a list of fold vectors, or ",
         "paired with `samplePartition`.")
  list(samplePartition = NULL, nFolds = k)
}

# Identify non-zero-variance columns of X. Returns a logical vector.
#' @importFrom matrixStats colSds
#' @noRd
.nonzeroVarColumns <- function(X) {
  sds <- colSds(X, na.rm = TRUE)
  !is.na(sds) & sds != 0
}

# Embed a smaller weights matrix into a full-sized zero matrix matching X and Y dimensions.
# @param weightsMatrix The fitted weights (nrow = number of valid columns).
# @param validColumns Logical or character vector identifying which columns of X were used.
# @param XColnames Column names of the original X.
# @param YColnames Column names of Y.
# @noRd
.embedWeights <- function(weightsMatrix, validColumns, nColsX, nColsY,
                          XColnames = NULL, YColnames = NULL) {
  full <- matrix(0, nrow = nColsX, ncol = nColsY)
  if (!is.null(XColnames)) rownames(full) <- XColnames
  if (!is.null(YColnames)) colnames(full) <- YColnames
  full[validColumns, ] <- weightsMatrix
  full
}

.prepareSusieWeightMethods <- function(X, Y, weightMethods, fittedModels = NULL) {
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  if (is.null(fittedModels)) fittedModels <- list()
  hasSusie <- !is.null(weightMethods[["susie_weights"]])
  hasSusieInf <- !is.null(weightMethods[["susie_inf_weights"]])
  susieFit <- if (hasSusie) weightMethods[["susie_weights"]][["susieFit"]] else NULL
  susieInfFit <- if (hasSusieInf) weightMethods[["susie_inf_weights"]][["susieInfFit"]] else NULL
  if (is.null(susieFit)) susieFit <- fittedModels[["susie"]]
  if (is.null(susieInfFit)) susieInfFit <- fittedModels[["susieInf"]]

  if (!is.null(susieFit)) {
    susieFit <- .setFinemappingFitClass(susieFit, "susie")
  }
  if (!is.null(susieInfFit)) {
    susieInfFit <- .setFinemappingFitClass(susieInfFit, "susieInf")
  }

  if (hasSusie && hasSusieInf && ncol(Y) == 1 &&
      is.null(susieFit) && is.null(susieInfFit)) {
    fitArgNames <- c("susieFit", "susieInfFit", "retainFit")
    fits <- fitSusieInfThenSusie(
      X,
      Y[, 1],
      args = weightMethods[["susie_weights"]][setdiff(names(weightMethods[["susie_weights"]]), fitArgNames)],
      susieInfArgs = modifyList(
        list(convergence_method = "pip"),
        weightMethods[["susie_inf_weights"]][setdiff(names(weightMethods[["susie_inf_weights"]]), fitArgNames)]
      ),
      fittedModels = list(susie = susieFit, susieInf = susieInfFit)
    )
    susieFit <- fits[["susie"]]
    susieInfFit <- fits[["susieInf"]]
  }

  if (!is.null(susieInfFit) && hasSusieInf) {
    weightMethods[["susie_inf_weights"]][["susieInfFit"]] <- susieInfFit
  }
  if (!is.null(susieFit) && hasSusie) {
    weightMethods[["susie_weights"]][["susieFit"]] <- susieFit
  }
  if (hasSusie &&
      is.null(weightMethods[["susie_weights"]][["susieFit"]]) &&
      !is.null(susieInfFit)) {
    weightMethods[["susie_weights"]] <- prepareSusieFromInfArgs(weightMethods[["susie_weights"]], susieInfFit)
  }
  weightMethods
}

# Per-fold TWAS weight fit for the CV engine. `ctx` carries weightMethods,
# multivariateWeightMethods, cvArgs, retainFits, verbose. Weights are keyed by the
# canonical method key; captured fits keep the full method name.
# @noRd
.weightFitFold <- function(Xtr, Ytr, j, ctx) {
  weightMethods <- ctx$weightMethods
  multivariateWeightMethods <- ctx$multivariateWeightMethods
  cvArgs <- ctx$cvArgs; retainFits <- ctx$retainFits; verbose <- ctx$verbose
  foldWeightMethods <- .prepareSusieWeightMethods(Xtr, Ytr, weightMethods)
  weights <- list()
  fits <- list()
  for (method in names(foldWeightMethods)) {
    args <- foldWeightMethods[[method]]
    fnName <- .resolveMethodFunction(method, args)
    mk <- sub("_weights$|Weights$", "", method)
    capturedFit <- NULL
    if (method %in% multivariateWeightMethods) {
      # Per-fold priors bind to the fitter's camelCase args.
      if (!is.null(cvArgs$data_driven_prior_matrices_cv) &&
          method %in% c("mrmash_weights", "mrmashWeights")) {
        args$dataDrivenPriorMatrices <- cvArgs$data_driven_prior_matrices_cv[[j]]
      }
      if (!is.null(cvArgs$reweightedMixturePriorCv) &&
          method %in% c("mvsusie_weights", "mvsusieWeights")) {
        args$prior_variance <- cvArgs$reweightedMixturePriorCv[[j]]
      }
      if (isTRUE(retainFits) && "retainFit" %in% names(formals(fnName))) {
        args$retainFit <- TRUE
      }
      W <- if (verbose < 2) {
        .quietEval(do.call(fnName, c(list(X = Xtr, Y = Ytr), args)))
      } else {
        do.call(fnName, c(list(X = Xtr, Y = Ytr), args))
      }
      capturedFit <- attr(W, "fit")
      attr(W, "fit") <- NULL
      rownames(W) <- colnames(Xtr)
    } else {
      Wcols <- lapply(seq_len(ncol(Ytr)), function(k) {
        w <- if (verbose < 2) {
          .quietEval(do.call(fnName, c(list(X = Xtr, y = Ytr[, k]), args)))
        } else {
          do.call(fnName, c(list(X = Xtr, y = Ytr[, k]), args))
        }
        as.numeric(w)
      })
      W <- do.call(cbind, Wcols)
      rownames(W) <- colnames(Xtr)
    }
    weights[[mk]] <- W
    fits[[method]] <- capturedFit
  }
  list(weights = weights, fits = fits)
}

#' Cross-Validation for weights selection in Transcriptome-Wide Association Studies (TWAS)
#'
#' Performs cross-validation for TWAS, supporting both univariate and multivariate methods.
#' It can either create folds for cross-validation or use pre-defined sample partitions.
#' For multivariate methods, it applies the method to the entire Y matrix for each fold.
#'
#' @param X A matrix of samples by features, where each row represents a sample and each column a feature.
#' @param Y A matrix (or vector, which will be converted to a matrix) of samples by outcomes, where each row corresponds to a sample.
#' @param fold An optional integer specifying the number of folds for cross-validation.
#' If NULL, 'samplePartitions' must be provided.
#' @param samplePartitions An optional dataframe with predefined sample partitions,
#' containing columns 'Sample' (sample names) and 'Fold' (fold number). If NULL, 'fold' must be provided.
#' @param weightMethods A list of methods and their specific arguments, formatted as list(method1 = method1_args, method2 = method2_args), or alternatively a character vector of method names (eg, c("susie_weights", "enet_weights")) in which case default arguments will be used for all methods.
#' methods in the list can be either univariate (applied to each column of Y) or multivariate (applied to the entire Y matrix).
#' @param maxNumVariants An optional integer to set the randomly selected maximum number of variants to use for CV purpose, to save computing time.
#' @param variantsToKeep An optional integer to ensure that the listed variants are kept in the CV when there is a limit on the maxNumVariants to use.
#' @param numThreads The number of threads to use for parallel processing.
#'        If set to -1, the function uses all available cores.
#'        If set to 0 or 1, no parallel processing is performed.
#'        If set to 2 or more, parallel processing is enabled with that many threads.
#' @param verbose Integer controlling verbosity level: 0 = suppress all messages,
#'   1 = suppress external package messages (default),
#'   2 = show all messages including those from external packages.
#' @return A list with the following components:
#' \itemize{
#'   \item `samplePartition`: A dataframe showing the sample partitioning used in the cross-validation.
#'   \item `prediction`: A list of matrices with predicted Y values for each method and fold.
#'   \item `metrics`: A matrix with rows representing methods and columns for various metrics:
#'     \itemize{
#'       \item `corr`: Pearson's correlation between predicated and observed values.
#'       \item `adj_rsq`: Adjusted R-squared value (which indicates the proportion of variance explained by the model) that accounts for the number of predictors in the model.
#'       \item `pval`: P-value assessing the significance of the model's predictions.
#'       \item `RMSE`: Root Mean Squared Error, a measure of the model's prediction error.
#'       \item `MAE`: Mean Absolute Error, a measure of the average magnitude of errors in a set of predictions.
#'     }
#'   \item `timeElapsed`: The time taken to complete the cross-validation process.
#' }
#' @importFrom purrr map
#' @importFrom BiocParallel bplapply bpworkers MulticoreParam
#' @importFrom quadprog solve.QP
#' @export
twasWeightsCv <- function(X, Y, fold = NULL, samplePartitions = NULL, weightMethods = NULL, maxNumVariants = NULL, variantsToKeep = NULL, numThreads = 1, verbose = 1, retainFits = FALSE, ...) {
  if (is.character(weightMethods)) {
    weightMethods <- .twasMethodLookup(weightMethods)
  }
  if (!exists(".Random.seed") && verbose >= 1) {
    message("! No seed has been set. Please set seed for reproducable result. ")
  }
  cvArgs <- list(...)
  # Multivariate weight methods (snake + camel) are fit on the whole Y for a
  # fold; univariate methods are fit per Y column. fSuSiE is intentionally
  # absent -- it is functional and cannot be refit from a bare (X, y) fold
  # split, so its cross-validated predictions are supplied by fineMappingPipeline.
  multivariateWeightMethods <- c("mrmash_weights", "mvsusie_weights",
                                 "mrmashWeights", "mvsusieWeights")

  # Per-fold fit context passed to the shared engine's top-level fitter
  # (.weightFitFold). Weights are keyed by the canonical method key; captured
  # fits keep the full method name (foldFits back-compat).
  cvFitCtx <- list(weightMethods = weightMethods,
                   multivariateWeightMethods = multivariateWeightMethods,
                   cvArgs = cvArgs, retainFits = retainFits, verbose = verbose)

  # No weight methods: the caller only wants the fold partition.
  if (is.null(weightMethods)) {
    res <- .crossValidateWeights(
      X, Y, fold = fold, samplePartitions = samplePartitions,
      fitFold = .cvNoopFitFold,
      numThreads = numThreads, maxNumVariants = maxNumVariants,
      variantsToKeep = variantsToKeep, retainFits = retainFits, verbose = verbose)
    return(list(samplePartition = res$samplePartition))
  }

  .crossValidateWeights(
    X, Y, fold = fold, samplePartitions = samplePartitions,
    fitFold = .weightFitFold, fitFoldCtx = cvFitCtx, numThreads = numThreads,
    maxNumVariants = maxNumVariants, variantsToKeep = variantsToKeep,
    retainFits = retainFits, verbose = verbose)
}

# Fit one TWAS weight method by name against the filtered design matrix, embedding
# the fitted weights back into the full variant space. `ctx` carries the shared
# fit state (X, Y, Xfiltered, validColumns, retainFits, retainFitDetail, verbose).
# @noRd
.computeMethodWeights <- function(methodName, weightMethods, ctx) {
  X <- ctx$X; Y <- ctx$Y; Xfiltered <- ctx$Xfiltered
  validColumns <- ctx$validColumns; retainFits <- ctx$retainFits
  retainFitDetail <- ctx$retainFitDetail; verbose <- ctx$verbose
  shortName <- sub("_weights$", "", methodName)
  if (verbose >= 1) {
    message(sprintf("  Fitting %s ...", shortName))
    tic()
  }

  # Hardcoded vector of multivariate methods (accept both snake and camel).
  # fSuSiE is multivariate (variants x features weight matrix) but is never
  # refit here — fsusieWeights extracts from the supplied fsusieFit.
  multivariateWeightMethods <- c("mrmash_weights", "mvsusie_weights",
                                  "fsusie_weights",
                                  "mrmashWeights", "mvsusieWeights",
                                  "fsusieWeights")
  args <- weightMethods[[methodName]]
  fnName <- .resolveMethodFunction(methodName, args)

  # Only pass retainFit (or its legacy snake_case alias) to functions that accept it
  if (retainFits) {
    fnFormals <- names(formals(fnName))
    if ("retainFit" %in% fnFormals) {
      args$retainFit <- TRUE
    } else if ("retain_fit" %in% fnFormals) {
      args$retain_fit <- TRUE
    }
    # Propagate the slim/full payload choice to producers that support it
    # (mr.mash individual + RSS), unless the caller already set it per-method.
    if ("fitDetail" %in% fnFormals && is.null(args$fitDetail)) {
      args$fitDetail <- retainFitDetail
    }
  }

  methodFit <- NULL
  if (methodName %in% multivariateWeightMethods) {
    # Apply multivariate method
    weightsMatrix <- if (verbose < 2) {
      .quietEval(do.call(fnName, c(list(X = Xfiltered, Y = Y), args)))
    } else {
      do.call(fnName, c(list(X = Xfiltered, Y = Y), args))
    }
    if (retainFits) methodFit <- attr(weightsMatrix, "fit")
    if (nrow(weightsMatrix) != length(validColumns)) weightsMatrix <- weightsMatrix[names(validColumns), , drop = FALSE]
  } else {
    # Apply univariate method to each column of Y
    # Initialize it with zeros to avoid NA
    weightsMatrix <- matrix(0, nrow = ncol(Xfiltered), ncol = ncol(Y))

    for (k in 1:ncol(Y)) {
      weightsVector <- if (verbose < 2) {
        .quietEval(do.call(fnName, c(list(X = Xfiltered, y = Y[, k]), args)))
      } else {
        do.call(fnName, c(list(X = Xfiltered, y = Y[, k]), args))
      }
      if (retainFits && is.null(methodFit)) {
        methodFit <- attr(weightsVector, "fit")
      }
      if (is.matrix(weightsVector)) weightsVector <- weightsVector[, k]
      weightsMatrix[, k] <- weightsVector
    }
  }

  result <- .embedWeights(weightsMatrix, validColumns, ncol(X), ncol(Y), colnames(X), colnames(Y))
  if (!is.null(methodFit)) attr(result, "fit") <- methodFit
  if (verbose >= 1) {
    elapsed <- toc(quiet = TRUE)
    message(sprintf("  Fitting %s done in %.1fs", shortName, elapsed$toc - elapsed$tic))
  }
  return(result)
}

# Assemble the (study, context, trait, method, entry) row vectors for the
# TwasWeights collection from the fitted `weightsList`. `ctx` carries the shared
# identity + flags (study, context, trait, Y, retainFits, standardized, dataType).
# @noRd
.buildTwasWeightEntries <- function(weightsList, variantIds, ctx) {
  Y <- ctx$Y; study <- ctx$study; context <- ctx$context; trait <- ctx$trait
  retainFits <- ctx$retainFits; standardized <- ctx$standardized; dataType <- ctx$dataType
  studies   <- character(0)
  contexts  <- character(0)
  traits    <- character(0)
  methodsV  <- character(0)
  entries   <- list()
  for (m in names(weightsList)) {
    wMat <- weightsList[[m]]
    fitVal <- attr(wMat, "fit")
    attr(wMat, "fit") <- NULL
    shortMethod <- sub("(_weights|Weights)$", "", m)
    # When trait/context were supplied per-row (length == ncol(Y)), emit
    # one row per (method, outcome). Otherwise emit one row per method
    # and carry the (possibly multi-column) weights matrix as-is.
    perOutcome <- length(trait) == ncol(Y) &&
                  length(context) %in% c(1L, ncol(Y))
    if (perOutcome) {
      contextV <- if (length(context) == 1L) rep(context, ncol(Y)) else context
      studyV   <- if (length(study) == 1L) rep(study, ncol(Y)) else study
      for (k in seq_len(ncol(Y))) {
        studies  <- c(studies,  studyV[k])
        contexts <- c(contexts, contextV[k])
        traits   <- c(traits,   trait[k])
        methodsV <- c(methodsV, shortMethod)
        entries[[length(entries) + 1L]] <- TwasWeightsEntry(
          variantIds    = variantIds,
          weights       = wMat[, k],
          fits          = if (retainFits) fitVal else NULL,
          cvResult = NULL,
          standardized  = isTRUE(standardized),
          dataType      = dataType)
      }
    } else {
      studies  <- c(studies,  study[1L])
      contexts <- c(contexts, context[1L])
      traits   <- c(traits,   trait[1L])
      methodsV <- c(methodsV, shortMethod)
      wPayload <- if (ncol(wMat) == 1L) drop(wMat) else wMat
      entries[[length(entries) + 1L]] <- TwasWeightsEntry(
        variantIds    = variantIds,
        weights       = wPayload,
        fits          = if (retainFits) fitVal else NULL,
        cvResult = NULL,
        standardized  = isTRUE(standardized),
        dataType      = dataType)
    }
  }
  list(study = studies, context = contexts, trait = traits,
       method = methodsV, entry = entries)
}

#' Run multiple TWAS weight methods
#'
#' Applies specified weight methods to the datasets X and Y, returning weight matrices for each method.
#' Handles both univariate and multivariate methods, and filters out columns in X with zero standard error.
#' This function utilizes parallel processing to handle multiple methods.
#'
#' @param X A matrix of samples by features, where each row represents a sample and each column a feature.
#' @param Y A matrix (or vector, which will be converted to a matrix) of samples by outcomes, where each row corresponds to a sample.
#' @param weightMethods A list of methods and their specific arguments, formatted as list(method1 = method1_args, method2 = method2_args), or alternatively a character vector of method names (eg, c("susie_weights", "enet_weights")) in which case default arguments will be used for all methods.
#' methods in the list can be either univariate (applied to each column of Y) or multivariate (applied to the entire Y matrix).
#' @param numThreads The number of threads to use for parallel processing.
#'        If set to -1, the function uses all available cores.
#'        If set to 0 or 1, no parallel processing is performed.
#'        If set to 2 or more, parallel processing is enabled with that many threads.
#' @param fittedModels Optional named list of fitted SuSiE-family models.
#' @param retainFits If TRUE, retain fitted model objects as attributes on
#'   returned weight matrices when supported by the weight method.
#' @param verbose Integer controlling verbosity level: 0 = suppress all messages,
#'   1 = suppress external package messages (default),
#'   2 = show all messages including those from external packages.
#' @return A list where each element is named after a method and contains the weight matrix produced by that method.
#'
#' @export
#' @importFrom purrr map exec
#' @importFrom rlang !!!
#' @importFrom tictoc tic toc
learnTwasWeights <- function(X, Y, weightMethods,
                             study = "", context = "", trait = "",
                             numThreads = 1,
                             fittedModels = NULL,
                             retainFits = FALSE,
                             retainFitDetail = c("slim", "full"),
                             standardized = FALSE,
                             dataType = NULL,
                             ldSketch = NULL,
                             verbose = 1) {
  if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
    stop("X must be a matrix and Y must be a matrix or a vector.")
  }
  retainFitDetail <- match.arg(retainFitDetail)

  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }

  if (nrow(X) != nrow(Y)) {
    stop("The number of rows in X and Y must be the same.")
  }

  if (is.character(weightMethods)) {
    weightMethods <- .twasMethodLookup(weightMethods)
  }

  # Determine number of cores to use
  numCores <- ifelse(numThreads == -1,
    bpworkers(MulticoreParam()),
    numThreads)
  numCores <- min(numCores,
    bpworkers(MulticoreParam()))

  validColumns <- .nonzeroVarColumns(X)
  Xfiltered <- as.matrix(X[, validColumns, drop = FALSE])
  weightMethods <- .prepareSusieWeightMethods(
    Xfiltered, Y, weightMethods, fittedModels
  )

  # Shared fit context threaded to the top-level per-method worker + entry builder.
  ctx <- list(X = X, Y = Y, Xfiltered = Xfiltered, validColumns = validColumns,
              study = study, context = context, trait = trait,
              retainFits = retainFits, retainFitDetail = retainFitDetail,
              standardized = standardized, dataType = dataType, verbose = verbose)

  if (numCores >= 2) {
    bpParam <- MulticoreParam(workers = numCores,
                              RNGseed = 1L)
    weightsList <- bplapply(names(weightMethods),
      .computeMethodWeights, weightMethods, ctx, BPPARAM = bpParam)
  } else {
    weightsList <- names(weightMethods) %>% map(.computeMethodWeights, weightMethods, ctx)
  }
  names(weightsList) <- names(weightMethods)

  if (!is.null(colnames(X))) {
    weightsList <- lapply(weightsList, function(x) {
      fit <- attr(x, "fit")
      rownames(x) <- colnames(X)
      if (!is.null(fit)) attr(x, "fit") <- fit
      return(x)
    })
  }

  variantIds <- if (!is.null(colnames(X))) colnames(X) else paste0("variant_", seq_len(ncol(X)))
  traitLabels <- if (!is.null(colnames(Y))) colnames(Y) else paste0("outcome_", seq_len(ncol(Y)))

  # Build one TwasWeightsEntry per (method, trait/outcome) row. For multi-
  # outcome (multivariate) methods the per-method weights matrix has one
  # column per outcome, so the same row in the TwasWeights collection
  # carries the matrix for that method across outcomes via a single
  # `trait` value taken from the input `trait` arg (when length 1) or the
  # corresponding Y column name when `trait` matches `colnames(Y)`.

  rows <- .buildTwasWeightEntries(weightsList, variantIds, ctx)
  TwasWeights(
    study   = rows$study,
    context = rows$context,
    trait   = rows$trait,
    method  = rows$method,
    entry   = rows$entry,
    ldSketch = ldSketch)
}

#' Predict outcomes using TWAS weights
#'
#' This function takes a matrix of predictors (\code{X}) and a list of TWAS (transcriptome-wide
#' association studies) weights (\code{weightsList}), and calculates the predicted outcomes by
#' multiplying \code{X} by each set of weights in \code{weightsList}. The names of the elements
#' in the output list are derived from the names in \code{weightsList}, with "_weights" replaced
#' by "_predicted".
#'
#' @param X A matrix or data frame of predictors where each row is an observation and each
#' column is a variable.
#' @param weightsList A list of numeric vectors representing the weights for each predictor.
#' The names of the list elements should follow the pattern \code{[outcome]_weights}, where
#' \code{[outcome]} is the name of the outcome variable that the weights are associated with.
#'
#' @return A named list of numeric vectors, where each vector is the predicted outcome for the
#' corresponding set of weights in \code{weightsList}. The names of the list elements are
#' derived from the names in \code{weightsList} by replacing "_weights" with "_predicted".
#'
#' @export
#' @examples
#' # Assuming `X` is your matrix of predictors and `weightsList` is your list of weights:
#' predicted_outcomes <- twasPredict(X, weightsList)
#' print(predicted_outcomes)
twasPredict <- function(X, weightsList) {
  if (is(weightsList, "TwasWeights")) {
    # Per-row weights vector/matrix payloads. Use the method name as key
    # for compatibility with the legacy snake_case "<method>_predicted"
    # convention; ensembleWeights() rebinds the suffix.
    methodNames <- as.character(weightsList$method)
    wl <- setNames(
      lapply(seq_len(nrow(weightsList)), function(i) {
        getWeights(weightsList$entry[[i]])
      }),
      paste0(methodNames, "_weights"))
  } else {
    wl <- weightsList
  }
  setNames(lapply(wl, function(w) {
    if (!is.matrix(w)) w <- matrix(w, ncol = 1)
    X %*% w
  }), .renameSuffix(names(wl), "predicted"))
}

#' Estimate Sparsity from mr.ash Mixture Proportions
#'
#' Computes an empirical estimate of the proportion of non-zero effects
#' (sparsity) from the mr.ash fit. mr.ash fits a mixture model with a
#' point mass at zero (spike) plus continuous components (slab), and
#' learns the mixture proportions via variational EM. The sparsity
#' estimate \code{1 - pi[1]} is the empirical Bayes estimate of the
#' non-null proportion, which can be used as a data-driven prior for
#' the inclusion probability parameters (\code{pi} for bayesC,
#' \code{probIn} for BayesB) of spike-and-slab Bayesian methods.
#'
#' @param weightResults Named list of weight vectors or matrices as
#'   returned by \code{\link{learnTwasWeights}}. The mr.ash element should
#'   have a \code{"fit"} attribute containing the model fit object
#'   (set \code{retainFits = TRUE} in \code{learnTwasWeights} to obtain this).
#'
#' @return A scalar sparsity estimate (proportion of non-zero effects).
#' @export
estimateSparsity <- function(weightResults) {
  if (is(weightResults, "TwasWeights")) {
    # Method names on the new TwasWeights collection are bare tokens
    # ("mrash"), not the snake_case _weights suffix form.
    methods <- as.character(weightResults$method)
    idx <- which(methods == "mrash")
    if (length(idx) == 0L) {
      stop("mr.ash entry not found in TwasWeights. Run learnTwasWeights() ",
           "with retainFits = TRUE and ensure 'mrash' is in the method list.")
    }
    fit <- getFits(weightResults$entry[[idx[[1L]]]])
    if (is.null(fit) || is.null(fit$pi)) {
      stop("mr.ash fit object not found. Run learnTwasWeights() with retainFits = TRUE ",
           "and ensure mrash_weights is included.")
    }
  } else {
    w <- weightResults[["mrash_weights"]]
    if (is.null(w)) {
      stop("mr.ash weights ('mrash_weights') not found in weightResults.")
    }
    fit <- attr(w, "fit")
    if (is.null(fit) || is.null(fit$pi)) {
      stop("mr.ash fit object not found. Run learnTwasWeights() with retainFits = TRUE ",
           "and ensure mrash_weights is included.")
    }
  }

  # fit$pi[1] is the weight on the spike (sa2[1] = 0); 1 - pi[1] = non-null proportion
  return(1 - fit$pi[1])
}

