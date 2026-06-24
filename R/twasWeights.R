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
      jointCols <- intersect(
        c("jointStudies", "jointContexts", "jointTraits"), names(object))
      for (jc in jointCols) {
        vals <- object[[jc]]
        if (!is.character(vals))
          errors <- c(errors, sprintf(
            "'%s' column must be character (got %s)", jc, class(vals)[[1L]]))
      }
      keyCols <- c("study", "context", "trait", "method", jointCols)
      keyDf <- as.data.frame(object[, keyCols, drop = FALSE])
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
#' @param ldSketch An optional \code{GenotypeHandle}, or \code{NULL} for
#'   individual-level fits.
#' @return A \code{TwasWeights} object.
#' @export
TwasWeights <- function(study, context, trait, method, entry,
                        jointStudies = NULL,
                        jointContexts = NULL,
                        jointTraits = NULL,
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
  df <- do.call(S4Vectors::DataFrame,
                c(cols, list(check.names = FALSE)))
  obj <- new("TwasWeights", df, ldSketch = ldSketch)
  validObject(obj)
  obj
}

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

#' @rdname getCvPerformance
#' @export
setMethod("getCvPerformance", "TwasWeights",
  function(x, study = NULL, context = NULL, trait = NULL, method = NULL,
           ...) {
    entry <- getTwasWeights(x, study, context, trait, method)
    getCvPerformance(entry)
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
    mrmash = list(fn = "mrmash_weights", impl = "mrmashWeights", args = list()),
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

# Resolve the actual function name for a method key. Honors an "impl" attribute
# on the per-method args list (set by .twasMethodLookup), and otherwise applies
# a snake_case -> camelCase transformation as a fallback for user-supplied
# weightMethods lists.
.resolveMethodFunction <- function(methodKey, methodArgs = NULL) {
  # Search pecotmr's namespace explicitly so this works equally well when the
  # function is called either from inside the package or from a user session.
  ns <- asNamespace("pecotmr")
  fnExists <- function(name) {
    exists(name, mode = "function") ||
      exists(name, mode = "function", envir = ns, inherits = FALSE)
  }
  impl <- if (!is.null(methodArgs)) attr(methodArgs, "impl") else NULL
  if (!is.null(impl) && nzchar(impl) && fnExists(impl)) {
    return(impl)
  }
  # Direct match (e.g. caller already passed camelCase)
  if (fnExists(methodKey)) return(methodKey)
  # snake_case_weights -> camelCaseWeights
  parts <- strsplit(methodKey, "_", fixed = TRUE)[[1]]
  capRest <- paste0(toupper(substring(parts[-1], 1, 1)),
                    substring(parts[-1], 2))
  candidate <- paste0(parts[1], paste0(capRest, collapse = ""))
  if (fnExists(candidate)) return(candidate)
  methodKey
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

# Filter weight methods that produced all-zero weights from CV.
# Returns filtered weightMethods list and warns about removed methods.
# @noRd
.filterZeroWeightMethods <- function(weightMethods, twasWeightsRes) {
  if (is(twasWeightsRes, "TwasWeights")) {
    methodTokens <- as.character(twasWeightsRes$method)
    perMethodAllZero <- vapply(seq_len(nrow(twasWeightsRes)), function(i) {
      w <- getWeights(twasWeightsRes$entry[[i]])
      all(w == 0, na.rm = TRUE)
    }, logical(1))
    methodToZero <- tapply(perMethodAllZero, methodTokens, all)
    methodKeys <- names(weightMethods)
    methodBase <- sub("(_weights|Weights)$", "", methodKeys)
    isAllZero <- vapply(methodBase, function(mb) {
      if (mb %in% names(methodToZero)) isTRUE(methodToZero[[mb]]) else FALSE
    }, logical(1))
  } else {
    wl <- twasWeightsRes
    isAllZero <- vapply(wl, function(w) all(w == 0, na.rm = TRUE), logical(1))
  }
  removed <- names(weightMethods)[isAllZero]
  if (length(removed) > 0) {
    warning(sprintf(
      "Methods %s are removed from CV because all their weights are zeros.",
      paste(removed, collapse = ", ")
    ))
  }
  weightMethods[!isAllZero]
}

.susieWeightIntermediate <- function(fit, X) {
  keep <- intersect(c("mu", "lbf_variable", "X_column_scale_factors", "pip", "theta"), names(fit))
  intermediate <- fit[keep]
  if (!is.null(fit$sets$cs)) {
    intermediate$csVariants <- setNames(lapply(fit$sets$cs, function(L) colnames(X)[L]), names(fit$sets$cs))
    intermediate$csPurity <- .translateSusiePurity(fit$sets$purity)
  }
  intermediate
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
twasWeightsCv <- function(X, Y, fold = NULL, samplePartitions = NULL, weightMethods = NULL, maxNumVariants = NULL, variantsToKeep = NULL, numThreads = 1, verbose = 1, ...) {
  splitData <- function(X, Y, samplePartition, fold) {
    testIds <- samplePartition[which(samplePartition$Fold == fold), "Sample"]
    Xtrain <- X[!(rownames(X) %in% testIds), , drop = FALSE]
    Ytrain <- Y[!(rownames(Y) %in% testIds), , drop = FALSE]
    Xtest <- X[rownames(X) %in% testIds, , drop = FALSE]
    Ytest <- Y[rownames(Y) %in% testIds, , drop = FALSE]
    if (nrow(Xtrain) == 0 || nrow(Ytrain) == 0 || nrow(Xtest) == 0 || nrow(Ytest) == 0) {
      stop("Error: One of the datasets (train or test) has zero rows.")
    }
    return(list(Xtrain = Xtrain, Ytrain = Ytrain, Xtest = Xtest, Ytest = Ytest))
  }

  # Validation checks
  if (!is.null(fold) && (!is.numeric(fold) || fold <= 0)) {
    stop("Invalid value for 'fold'. It must be a positive integer.")
  }

  if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
    stop("X must be a matrix and Y must be a matrix or a vector.")
  }

  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
    if (verbose >= 1) message(paste("Y converted to matrix of", nrow(Y), "rows and", ncol(Y), "columns."))
  }

  if (nrow(X) != nrow(Y)) {
    stop("The number of rows in X and Y must be the same.")
  }
  if (!is.null(rownames(X)) && !is.null(rownames(Y))) {
    if (!identical(rownames(X), rownames(Y))) {
      rownames(X) <- rownames(Y)
    }
    sampleNames <- rownames(Y)
  } else if (!is.null(rownames(Y))) {
    sampleNames <- rownames(Y)
  } else if (!is.null(rownames(X))) {
    sampleNames <- rownames(X)
  } else {
    sampleNames <- paste0("sample_", 1:nrow(X))
  }
  if (is.null(rownames(X))) {
    rownames(X) <- sampleNames
  }
  if (is.null(rownames(Y))) {
    rownames(Y) <- sampleNames
  }

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("variable_", 1:ncol(X))
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("context_", 1:ncol(Y))
  }

  if (is.character(weightMethods)) {
    weightMethods <- .twasMethodLookup(weightMethods)
  }

  if (!exists(".Random.seed")) {
    if (verbose >= 1) message("! No seed has been set. Please set seed for reproducable result. ")
  }

  # Select variants if necessary
  if (!is.null(maxNumVariants) && ncol(X) > maxNumVariants) {
    if (!is.null(variantsToKeep) && length(variantsToKeep) > 0) {
      variantsToKeep <- intersect(variantsToKeep, colnames(X))
      remainingColumns <- setdiff(colnames(X), variantsToKeep)
      if (length(variantsToKeep) < maxNumVariants) {
        additionalColumns <- sample(remainingColumns, maxNumVariants - length(variantsToKeep), replace = FALSE)
        selectedColumns <- union(variantsToKeep, additionalColumns)
        if (verbose >= 1) message(sprintf(
          "Including %d specified variants and randomly selecting %d additional variants, for a total of %d variants out of %d for cross-validation purpose.",
          length(variantsToKeep), length(additionalColumns), length(selectedColumns), ncol(X)
        ))
      } else {
        selectedColumns <- sample(variantsToKeep, maxNumVariants, replace = FALSE)
        if (verbose >= 1) message(paste("Randomly selecting", length(selectedColumns), "out of", length(variantsToKeep), "input variants for cross validation purpose."))
      }
    } else {
      selectedColumns <- sort(sample(ncol(X), maxNumVariants, replace = FALSE))
      if (verbose >= 1) message(paste("Randomly selecting", length(selectedColumns), "out of", ncol(X), "variants for cross validation purpose."))
    }
    X <- X[, selectedColumns, drop = FALSE]
  }

  # Create or use provided folds
  if (!is.null(fold)) {
    if (!is.null(samplePartitions)) {
      if (fold != length(unique(samplePartitions$Fold))) {
        if (verbose >= 1) message(paste0(
          "fold number provided does not match with sample partition, performing ", length(unique(samplePartitions$Fold)),
          " fold cross validation based on provided sample partition. "
        ))
      }

      folds <- samplePartitions$Fold
      samplePartition <- samplePartitions
    } else {
      sampleIndices <- sample(nrow(X))
      folds <- cut(seq(1, nrow(X)), breaks = fold, labels = FALSE)
      samplePartition <- data.frame(Sample = sampleNames[sampleIndices], Fold = folds, stringsAsFactors = FALSE)
    }
  } else if (!is.null(samplePartitions)) {
    if (!all(samplePartitions$Sample %in% sampleNames)) {
      stop("Some samples in 'samplePartitions' do not match the samples in 'X' and 'Y'.")
    }
    samplePartition <- samplePartitions
    fold <- length(unique(samplePartition$Fold))
  } else {
    stop("Either 'fold' or 'samplePartitions' must be provided.")
  }

  st <- proc.time()
  if (is.null(weightMethods)) {
    return(list(samplePartition = samplePartition))
  } else {
    # Hardcoded vector of multivariate weightMethods (accept both snake and camel).
    # fSuSiE is excluded from the per-fold CV refit path: it is functional and
    # cannot be refit from a bare (X, y) fold split, so its cross-validated
    # predictions are supplied by fineMappingPipeline (FineMappingResult
    # cvResult) rather than recomputed here.
    multivariateWeightMethods <- c("mrmash_weights", "mvsusie_weights",
                                    "mrmashWeights", "mvsusieWeights")

    # Determine the number of cores to use
    numCores <- ifelse(numThreads == -1,
      bpworkers(MulticoreParam()),
      numThreads)
    numCores <- min(numCores,
      bpworkers(MulticoreParam()))

    cvArgs <- list(...)

    # Perform CV with parallel processing
    computeMethodPredictions <- function(j) {
      if (verbose >= 1) {
        message(sprintf("  CV fold %d/%d ...", j, fold))
        tic()
      }
      datSplit <- splitData(X, Y, samplePartition = samplePartition, fold = j)
      Xtrain <- datSplit$Xtrain
      Ytrain <- datSplit$Ytrain
      Xtest <- datSplit$Xtest
      Ytest <- datSplit$Ytest

      # Remove columns with zero variance in the training fold.
      # NOTE: Y was already NA-handled at the pipeline layer
      # (getResidualizedPhenotypes naAction = "drop"/"impute"), so there
      # is no need to mask X by Y NA rows here.
      validColumns <- .nonzeroVarColumns(Xtrain)
      Xtrain <- Xtrain[, validColumns, drop = FALSE]
      validColumns <- colnames(Xtrain)
      # Xtest <- Xtest[, validColumns, drop=FALSE]
      foldWeightMethods <- .prepareSusieWeightMethods(Xtrain, Ytrain, weightMethods)

      foldPreds <- setNames(lapply(names(foldWeightMethods), function(method) {
        args <- foldWeightMethods[[method]]
        fnName <- .resolveMethodFunction(method, args)

        if (method %in% multivariateWeightMethods) {
          # Apply multivariate method to entire Y for this fold
          if (!is.null(cvArgs$data_driven_prior_matrices_cv)) {
            if (method %in% c("mrmash_weights", "mrmashWeights")) {
              args$data_driven_prior_matrices <- cvArgs$data_driven_prior_matrices_cv[[j]]
            }
            if (method %in% c("mvsusie_weights", "mvsusieWeights")) {
              args$prior_variance <- cvArgs$reweightedMixturePriorCv[[j]]
            }
          }
          weightsMatrix <- if (verbose < 2) {
            .quietEval(do.call(fnName, c(list(X = Xtrain, Y = Ytrain), args)))
          } else {
            do.call(fnName, c(list(X = Xtrain, Y = Ytrain), args))
          }
          rownames(weightsMatrix) <- colnames(Xtrain)
          fullWeightsMatrix <- .embedWeights(weightsMatrix[validColumns, , drop = FALSE], validColumns, ncol(X), ncol(Y), colnames(X), colnames(Y))
          Ypred <- Xtest %*% fullWeightsMatrix
          rownames(Ypred) <- rownames(Xtest)
          return(Ypred)
        } else {
          Ypred <- sapply(1:ncol(Ytrain), function(k) {
            weights <- if (verbose < 2) {
              .quietEval(do.call(fnName, c(list(X = Xtrain, y = Ytrain[, k]), args)))
            } else {
              do.call(fnName, c(list(X = Xtrain, y = Ytrain[, k]), args))
            }
            fullWeights <- rep(0, ncol(X))
            names(fullWeights) <- colnames(X)
            fullWeights[validColumns] <- weights
            # Handle NAs in weights
            fullWeights[is.na(fullWeights)] <- 0
            Xtest %*% fullWeights
          })
          rownames(Ypred) <- rownames(Xtest)
          return(Ypred)
        }
      }), names(foldWeightMethods))
      if (verbose >= 1) {
        elapsed <- toc(quiet = TRUE)
        message(sprintf("  CV fold %d/%d done in %.1fs", j, fold, elapsed$toc - elapsed$tic))
      }
      foldPreds
    }

    if (numCores >= 2) {
      bpParam <- MulticoreParam(workers = numCores,
                                RNGseed = 1L)
      foldResults <- bplapply(1:fold,
        computeMethodPredictions, BPPARAM = bpParam)
    } else {
      foldResults <- map(1:fold, computeMethodPredictions)
    }

    # Reorganize into Ypred
    # After cross validation, each sample should have been in
    # test set at some point, and therefore has predicted value.
    # The prediction matrix is therefore exactly the same dimension as input Y
    Ypred <- setNames(lapply(weightMethods, function(x) `dimnames<-`(matrix(NA, nrow(Y), ncol(Y)), dimnames(Y))), names(weightMethods))
    for (j in seq_along(foldResults)) {
      for (method in names(weightMethods)) {
        Ypred[[method]][rownames(foldResults[[j]][[method]]), ] <- foldResults[[j]][[method]]
      }
    }

    names(Ypred) <- .renameSuffix(names(Ypred), "predicted")

    # Compute rsq, adj rsq, p-value, RMSE, and MAE for each method
    metricsTable <- list()

    for (m in names(weightMethods)) {
      metricsTable[[m]] <- matrix(NA, nrow = ncol(Y), ncol = 6)
      colnames(metricsTable[[m]]) <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
      rownames(metricsTable[[m]]) <- colnames(Y)

      for (r in 1:ncol(Y)) {
        methodPredictions <- Ypred[[.renameSuffix(m, "predicted")]][, r]
        actualValues <- Y[, r]
        # Remove missing values in the first place
        naIndx <- which(is.na(actualValues))
        if (length(naIndx) != 0) {
          methodPredictions <- methodPredictions[-naIndx]
          actualValues <- actualValues[-naIndx]
        }
        if (sd(methodPredictions) != 0) {
          lmFit <- lm(actualValues ~ methodPredictions)

          # Calculate raw correlation and and adjusted R-squared
          metricsTable[[m]][r, "corr"] <- cor(actualValues, methodPredictions)

          metricsTable[[m]][r, "rsq"] <- summary(lmFit)$r.squared
          metricsTable[[m]][r, "adj_rsq"] <- summary(lmFit)$adj.r.squared

          # Calculate p-value
          metricsTable[[m]][r, "pval"] <- summary(lmFit)$coefficients[2, 4]

          # Calculate RMSE
          residuals <- actualValues - methodPredictions
          metricsTable[[m]][r, "RMSE"] <- sqrt(mean(residuals^2))

          # Calculate MAE
          metricsTable[[m]][r, "MAE"] <- mean(abs(residuals))
        } else {
          metricsTable[[m]][r, ] <- NA
          if (verbose >= 1) message(paste0(
            "Predicted values for condition ", r, " using ", m,
            " have zero variance. Filling performance metric with NAs"
          ))
        }
      }
    }
    names(metricsTable) <- .renameSuffix(names(metricsTable), "performance")
    return(list(samplePartition = samplePartition, prediction = Ypred, performance = metricsTable, timeElapsed = proc.time() - st))
  }
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
                             standardized = FALSE,
                             dataType = NULL,
                             ldSketch = NULL,
                             verbose = 1) {
  if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
    stop("X must be a matrix and Y must be a matrix or a vector.")
  }

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

  computeMethodWeights <- function(methodName, weightMethods) {
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

  if (numCores >= 2) {
    bpParam <- MulticoreParam(workers = numCores,
                              RNGseed = 1L)
    weightsList <- bplapply(names(weightMethods),
      computeMethodWeights, weightMethods, BPPARAM = bpParam)
  } else {
    weightsList <- names(weightMethods) %>% map(computeMethodWeights, weightMethods)
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
  buildEntries <- function() {
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
            cvPerformance = NULL,
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
          cvPerformance = NULL,
          standardized  = isTRUE(standardized),
          dataType      = dataType)
      }
    }
    list(study = studies, context = contexts, trait = traits,
         method = methodsV, entry = entries)
  }

  rows <- buildEntries()
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

