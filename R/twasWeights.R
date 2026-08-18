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

setClass(
    "TwasWeights",
    contains = "DFrame",
    representation(ldSketch = "ANY"),
    validity = function(object) .validateTwasWeights(object)
)

# Validity for the TwasWeights collection: required columns, entry payloads,
# region/traitPos provenance, joint* column types, tuple uniqueness, and the
# optional ldSketch. Returns TRUE or a character vector of error messages.
# @noRd
.validateTwasWeights <- function(object) {
    errors <- .twasValidateRequiredCols(object)
    if (length(errors) == 0L) {
        errors <- .twasValidateColumns(object)
    }
    errors <- c(errors, .twasValidateLdSketch(object))
    if (length(errors) == 0L) TRUE else errors
}

# Required key/entry columns must all be present.
# @noRd
.twasValidateRequiredCols <- function(object) {
    required <- c("study", "context", "trait", "method", "entry")
    missingCols <- setdiff(required, names(object))
    if (length(missingCols) > 0L) {
        str_c("missing columns: ", str_flatten(missingCols, ", "))
    } else {
        character()
    }
}

# Column-level checks that run only once the required columns are present.
# @noRd
.twasValidateColumns <- function(object) {
    jointCols <- intersect(
        c("jointStudies", "jointContexts", "jointTraits"),
        names(object)
    )
    c(
        .twasValidateEntries(object),
        .validateRegionColumn(object),
        .validateTraitPosColumn(object),
        .twasValidateJointCols(object, jointCols),
        .twasValidateKeyUniqueness(object, jointCols)
    )
}

# `entry` column: one TwasWeightsEntry per row.
# @noRd
.twasValidateEntries <- function(object) {
    errors <- character()
    if (length(object$entry) != nrow(object)) {
        errors <- c(errors, "length(entry) must equal nrow(.) for TwasWeights")
    }
    entryOk <- map_lgl(object$entry, methods::is, "TwasWeightsEntry")
    if (!all(entryOk)) {
        errors <- c(
            errors,
            str_c(
                "every element of the `entry` column must be a ",
                "TwasWeightsEntry"
            )
        )
    }
    errors
}

# Any present joint* provenance columns must be character.
# @noRd
.twasValidateJointCols <- function(object, jointCols) {
    bad <- keep(jointCols, .twasColNotCharacter, object = object)
    map_chr(bad, .twasBadColMsg, object = object)
}

# (study, context, trait, method[, joint*]) tuple uniqueness.
# @noRd
.twasValidateKeyUniqueness <- function(object, jointCols) {
    keyCols <- c("study", "context", "trait", "method", jointCols)
    # Extract key columns directly rather than via `object[, keyCols]`:
    # column-subsetting preserves the TwasWeights class while dropping the
    # required `entry` column, and older S4Vectors revalidates that
    # intermediate, spuriously failing with "missing columns: entry".
    keyTbl <- as_tibble(set_names(
        map(keyCols, .twasColumn, object = object),
        keyCols
    ))
    if (nrow(distinct(keyTbl)) < nrow(keyTbl)) {
        str_c(
            "(study, context, trait, method[, joint*]) tuple uniqueness ",
            "violated"
        )
    } else {
        character()
    }
}

# Optional ldSketch slot must be a GenotypeHandle or NULL.
# @noRd
.twasValidateLdSketch <- function(object) {
    if (
        !is.null(object@ldSketch) &&
            !methods::is(object@ldSketch, "GenotypeHandle")
    ) {
        "'ldSketch' must be a GenotypeHandle or NULL"
    } else {
        character()
    }
}


# =============================================================================

#' @title Create a TwasWeights Collection Object
#' @description Construct a \code{TwasWeights} DFrame-subclass collection from
#'   per-tuple vectors and a list of \code{TwasWeightsEntry} payloads (one per
#'   tuple).
#' @param study Character vector of study identifiers. Use the sentinel
#'   \code{"joint"} for rows produced by a cross-study joint fit.
#' @param context Character vector of context labels. Use \code{"joint"} for
#'   rows produced by a cross-context joint fit.
#' @param trait Character vector of trait identifiers. Use \code{"joint"} for
#'   rows produced by a cross-trait joint fit.
#' @param method Character vector of TWAS weight method names.
#' @param entry List / \code{SimpleList} of \code{TwasWeightsEntry} objects.
#' @param jointStudies Optional character vector (length \code{length(study)})
#'   listing the semicolon-joined studies participating in each row's
#'   cross-study joint fit, or \code{NA_character_} for non-joint rows. When
#'   \code{NULL} (default) the column is omitted.
#' @param jointContexts Optional character vector for cross-context joints. Same
#'   shape as \code{jointStudies}.
#' @param jointTraits Optional character vector for cross-trait joints. Same
#'   shape as \code{jointStudies}.
#' @param region Optional \code{GRanges} (length \code{length(study)}) giving
#'   the genomic anchor of each row's trait -- the trait's own coordinates when
#'   built from a \code{QtlDataset}, or the summary-stat variant-span window
#'   when built from a \code{QtlSumStats}. Carried forward as provenance (e.g.
#'   for cTWAS LD-block placement); not part of the identity key. \code{NULL}
#'   (default) omits the column.
#' @param ldSketch An optional \code{GenotypeHandle}, or \code{NULL} for
#'   individual-level fits.
#' @param traitPos Optional per-row trait genomic anchor (a \code{GRanges} or
#'   \code{NULL}), carried forward as provenance; not part of the identity key.
#'   \code{NULL} (default) omits the column.
#' @return A \code{TwasWeights} object.
#' @examples
#' twe <- TwasWeightsEntry(variantIds = paste0("v", 1:4),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5), standardized = FALSE)
#' tw <- TwasWeights(study = "s1", context = "brain", trait = "gene1",
#'   method = "susie", entry = list(twe))
#' tw
#' @export
TwasWeights <- function(
    study,
    context,
    trait,
    method,
    entry,
    jointStudies = NULL,
    jointContexts = NULL,
    jointTraits = NULL,
    region = NULL,
    traitPos = NULL,
    ldSketch = NULL
) {
    n <- .twasCheckRowLengths(study, context, trait, method, entry)
    cols <- list(
        study = as.character(study),
        context = as.character(context),
        trait = as.character(trait),
        method = as.character(method),
        entry = S4Vectors::SimpleList(entry)
    )
    cols <- .twasAppendJointCols(
        cols,
        jointStudies,
        jointContexts,
        jointTraits,
        n
    )
    cols <- .appendRegionCol(cols, region, n)
    cols <- .appendTraitPosCol(cols, traitPos, n)
    dfArgs <- c(cols, list(check.names = FALSE))
    df <- exec(S4Vectors::DataFrame, !!!dfArgs)
    obj <- new("TwasWeights", df, ldSketch = ldSketch)
    validObject(obj)
    obj
}

# Require study/context/trait/method/entry to share one length; returns it.
# @noRd
.twasCheckRowLengths <- function(study, context, trait, method, entry) {
    n <- length(study)
    if (
        length(context) != n ||
            length(trait) != n ||
            length(method) != n ||
            length(entry) != n
    ) {
        msg <- glue(
            "`study`, `context`, `trait`, `method`, and `entry` must all ",
            "have the same length."
        )
        abort(msg)
    }
    n
}

# Append any supplied joint* provenance columns (each length n) as character.
# @noRd
.twasAppendJointCols <- function(
    cols,
    jointStudies,
    jointContexts,
    jointTraits,
    n
) {
    joint <- list(
        jointStudies = jointStudies,
        jointContexts = jointContexts,
        jointTraits = jointTraits
    )
    for (nm in names(joint)) {
        val <- joint[[nm]]
        if (is.null(val)) {
            next
        }
        if (length(val) != n) {
            msg <- glue("`{nm}` must have the same length as `study`.")
            abort(msg)
        }
        cols[[nm]] <- as.character(val)
    }
    cols
}

#' @rdname getRegion
#' @export
setMethod("getRegion", "TwasWeights", function(x, ...) .getRegionColumn(x))

#' @rdname getTraitPosition
#' @export
setMethod("getTraitPosition", "TwasWeights", function(x, ...) {
    .getTraitPosColumn(x)
})

#' @title Get a Single TWAS Weights Entry
#' @description Return the \code{TwasWeightsEntry} for one \code{(study,
#'   context, trait, method)} row of a \code{TwasWeights} collection.
#' @param x A \code{TwasWeights} object.
#' @param study,context,trait,method Single character identifiers. All required
#'   when the collection has more than one row; optional when the collection has
#'   a single row.
#' @return A \code{TwasWeightsEntry} object.
#' @examples
#' twe <- TwasWeightsEntry(variantIds = paste0("v", 1:4), weights = rep(0.1, 4),
#'   cvResult = list(rsq = 0.5), standardized = FALSE)
#' tw <- TwasWeights(study = "s1", context = "brain", trait = "g1",
#'   method = "susie", entry = list(twe))
#' getTwasWeights(tw, study = "s1", context = "brain", trait = "g1",
#'   method = "susie")
#' @export
setGeneric(
    "getTwasWeights",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        standardGeneric("getTwasWeights")
    }
)

#' @rdname getTwasWeights
#' @export
setMethod(
    "getTwasWeights",
    "TwasWeights",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        idx <- .tupleSelectRow(
            x,
            study,
            context,
            trait,
            method,
            cls = "TwasWeights"
        )
        x$entry[[idx]]
    }
)

#' @rdname getWeights
#' @export
setMethod(
    "getWeights",
    "TwasWeights",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        ...
    ) {
        entry <- getTwasWeights(x, study, context, trait, method)
        getWeights(entry)
    }
)

#' @rdname getCvResult
#' @export
setMethod(
    "getCvResult",
    "TwasWeights",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        ...
    ) {
        entry <- getTwasWeights(x, study, context, trait, method)
        getCvResult(entry)
    }
)

#' @rdname getFits
#' @export
setMethod(
    "getFits",
    "TwasWeights",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        ...
    ) {
        entry <- getTwasWeights(x, study, context, trait, method)
        getFits(entry)
    }
)

#' @rdname getStandardized
#' @export
setMethod(
    "getStandardized",
    "TwasWeights",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        entry <- getTwasWeights(x, study, context, trait, method)
        getStandardized(entry)
    }
)

#' @rdname getDataType
#' @export
setMethod(
    "getDataType",
    "TwasWeights",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        entry <- getTwasWeights(x, study, context, trait, method)
        getDataType(entry)
    }
)

#' @rdname getVariantIds
#' @export
setMethod(
    "getVariantIds",
    "TwasWeights",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        ...
    ) {
        entry <- getTwasWeights(x, study, context, trait, method)
        getVariantIds(entry)
    }
)

#' @rdname getStudy
#' @export
setMethod("getStudy", "TwasWeights", function(x) unique(as.character(x$study)))

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "TwasWeights", function(x, ...) x@ldSketch)

#' @rdname getContexts
#' @export
setMethod("getContexts", "TwasWeights", function(x) {
    unique(as.character(x$context))
})

#' @rdname getTraits
#' @export
setMethod("getTraits", "TwasWeights", function(x) unique(as.character(x$trait)))

#' @rdname getMethodNames
#' @export
setMethod("getMethodNames", "TwasWeights", function(x) {
    unique(as.character(x$method))
})


#' @rdname show-methods
#' @export
setMethod("show", "TwasWeights", function(object) {
    cat(glue("TwasWeights: {nrow(object)} entries\n", .trim = FALSE))
    if (nrow(object) > 0L) {
        cat(glue(
            "  {n_distinct(object$study)} studies, ",
            "{n_distinct(object$context)} contexts, ",
            "{n_distinct(object$trait)} traits, ",
            "{n_distinct(object$method)} methods\n",
            .trim = FALSE
        ))
    }
    ldSrc <- if (is.null(object@ldSketch)) {
        "NULL (individual-level fit)"
    } else {
        glue("{object@ldSketch@format} @ {object@ldSketch@path}")
    }
    cat(glue("  LD sketch: {ldSrc}\n", .trim = FALSE))
})


# =============================================================================
# TwasWeights pipeline helpers (learnTwasWeights + CV + ensemble +
# Mvsusie/Mrmash)
# =============================================================================

# Evaluate an expression while suppressing external package output.
# Catches both message() output (susieR, qgg) and Rprintf/cat stdout
# (mr.ash.alpha).
# @param expr An expression to evaluate.
# @return The result of evaluating expr.
# @noRd
.quietEval <- function(expr) {
    invisible(utils::capture.output(
        result <- suppressMessages(expr),
        type = "output"
    ))
    result
}

# Rename a "_weights"/"Weights" suffix to the case-matching equivalent of
# `target`. Snake-case inputs get the underscored snake-case form (e.g.
# "lasso_weights" -> "lasso_predicted") and camelCase inputs get the CamelCase
# form (e.g. "lassoWeights" -> "lassoPredicted"). Names without a recognized
# suffix are returned unchanged.
# @param x Character vector of names ending in "_weights" or "Weights".
# @param target A bare token such as "predicted" or "performance".
# @return Character vector with suffixes rewritten.
# @noRd
.renameSuffix <- function(x, target) {
    cap <- str_c(
        str_to_upper(str_sub(target, 1, 1)),
        str_sub(target, 2)
    )
    x <- str_replace(x, "_weights$", str_c("_", target))
    x <- str_replace(x, "Weights$", cap)
    x
}

# Method name/impl/args lookup for TWAS weight methods. `fn` is the
# snake_case key used in weight method lists; `impl` is the camelCase
# function implemented by the package.
# @noRd
.twasMethodMap <- list(
    susie = list(
        fn = "susie_weights",
        impl = "susieWeights",
        args = list(refine = FALSE, L = 20, L_greedy = 5)
    ),
    susieAsh = list(
        fn = "susie_ash_weights",
        impl = "susieAshWeights",
        args = list()
    ),
    susieInf = list(
        fn = "susie_inf_weights",
        impl = "susieInfWeights",
        args = list()
    ),
    mrash = list(
        fn = "mrash_weights",
        impl = "mrashWeights",
        args = list(initPriorSd = TRUE, max.iter = 100)
    ),
    enet = list(fn = "enet_weights", impl = "enetWeights", args = list()),
    lasso = list(
        fn = "lasso_weights",
        impl = "lassoWeights",
        args = list()
    ),
    bayes_r = list(
        fn = "bayes_r_weights",
        impl = "bayesRWeights",
        args = list()
    ),
    bayes_l = list(
        fn = "bayes_l_weights",
        impl = "bLassoWeights",
        args = list()
    ),
    bayes_a = list(
        fn = "bayes_a_weights",
        impl = "bayesAWeights",
        args = list()
    ),
    bayes_b = list(
        fn = "bayes_b_weights",
        impl = "bayesBWeights",
        args = list()
    ),
    bayes_c = list(
        fn = "bayes_c_weights",
        impl = "bayesCWeights",
        args = list()
    ),
    bayes_n = list(
        fn = "bayes_n_weights",
        impl = "bayesNWeights",
        args = list()
    ),
    b_lasso = list(
        fn = "b_lasso_weights",
        impl = "bLassoWeights",
        args = list()
    ),
    dpr_vb = list(
        fn = "dpr_vb_weights",
        impl = "dprVbWeights",
        args = list()
    ),
    dpr_gibbs = list(
        fn = "dpr_gibbs_weights",
        impl = "dprGibbsWeights",
        args = list()
    ),
    dpr_adaptive_gibbs = list(
        fn = "dpr_adaptive_gibbs_weights",
        impl = "dprAdaptiveGibbsWeights",
        args = list()
    ),
    scad = list(fn = "scad_weights", impl = "scadWeights", args = list()),
    mcp = list(fn = "mcp_weights", impl = "mcpWeights", args = list()),
    l0learn = list(
        fn = "l0learn_weights",
        impl = "l0learnWeights",
        args = list()
    ),
    mvsusie = list(
        fn = "mvsusie_weights",
        impl = "mvsusieWeights",
        args = list(L = 30, L_greedy = 5)
    ),
    mrmash = list(
        fn = "mrmash_weights",
        impl = "mrmashWeights",
        args = list(canonicalPriorMatrices = TRUE)
    ),
    fsusie = list(
        fn = "fsusie_weights",
        impl = "fsusieWeights",
        args = list()
    )
)

# Expand the `default` / `fast_default` preset strings to their method vectors.
# @noRd
.twasExpandPresets <- function(methods) {
    fastDefault <- c(
        "susie",
        "susieInf",
        "mrash",
        "enet",
        "lasso",
        "mcp",
        "scad",
        "l0learn"
    )
    if (length(methods) == 1) {
        if (methods == "fast_default") {
            return(fastDefault)
        }
        if (methods == "default") {
            return(c(fastDefault, "bayes_r", "bayes_c"))
        }
    }
    methods
}

# Map short method names and presets to weightMethods lists.
# @param methods A character vector of short method names, or a preset string
#   ("default" or "fast_default").
# @return A named list suitable for the weightMethods parameter.
# @importFrom purrr map_chr set_names
# @noRd
.twasMethodLookup <- function(methods) {
    methods <- .twasExpandPresets(methods)
    # Accept full function names too by mapping them back to short names.
    fnToShort <- set_names(names(.twasMethodMap), map_chr(.twasMethodMap, "fn"))
    methods <- map_chr(methods, .twasCanonicalShortName, fnToShort = fnToShort)
    unknown <- setdiff(methods, names(.twasMethodMap))
    if (length(unknown) > 0) {
        msg <- glue(
            "Unknown TWAS method(s): {str_flatten(unknown, ', ')}. ",
            "Available methods: {str_flatten(names(.twasMethodMap), ', ')}"
        )
        abort(msg)
    }
    # Track the impl name as an attr so dispatchers resolve snake_case -> impl.
    entries <- map(methods, .twasMethodArgsWithImpl)
    set_names(entries, map_chr(methods, .twasMethodFn))
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
    if (
        !is.null(impl) && str_length(impl) > 0L && .functionExistsInNs(impl, ns)
    ) {
        return(impl)
    }
    # Direct match (e.g. caller already passed camelCase)
    if (.functionExistsInNs(methodKey, ns)) {
        return(methodKey)
    }
    # snake_case_weights -> camelCaseWeights
    parts <- str_split(methodKey, "_")[[1]]
    capRest <- str_c(
        str_to_upper(str_sub(parts[-1], 1, 1)),
        str_sub(parts[-1], 2)
    )
    candidate <- str_c(parts[1], str_flatten(capRest))
    if (.functionExistsInNs(candidate, ns)) {
        return(candidate)
    }
    methodKey
}

# Validate a Sample/Fold partition data frame: required columns, no sample in
# two folds, and (when `sampleNames` given) exact coverage of all samples.
# @noRd
.validateFoldPartition <- function(df, sampleNames) {
    if (!all(is_in(c("Sample", "Fold"), names(df)))) {
        abort("samplePartition must have columns `Sample` and `Fold`.")
    }
    df$Sample <- as.character(df$Sample)
    dup <- unique(df$Sample[duplicated(df$Sample)])
    if (length(dup) > 0L) {
        msg <- glue(
            "Fold partition assigns sample(s) to more than one fold: ",
            "{str_flatten(dup, ', ')}"
        )
        abort(msg)
    }
    if (!is.null(sampleNames)) {
        unknown <- setdiff(df$Sample, sampleNames)
        if (length(unknown) > 0L) {
            msg <- glue(
                "Fold partition references unknown sample(s): ",
                "{str_flatten(unknown, ', ')}"
            )
            abort(msg)
        }
        uncovered <- setdiff(sampleNames, df$Sample)
        if (length(uncovered) > 0L) {
            msg <- glue(
                "Fold partition does not cover {length(uncovered)} ",
                "sample(s) (folds must partition all samples), e.g. ",
                "{str_flatten(utils::head(uncovered, 5L), ', ')}"
            )
            abort(msg)
        }
    }
    df
}

# Normalize a cross-validation fold specification into the canonical
# samplePartition data.frame(Sample, Fold) used throughout the CV machinery, so
# callers can pass folds in any of three forms and downstream code has a single
# source of truth. Accepts: * cvFolds an integer k: k-fold auto-partition
# (returned as NULL so the partition is generated per (study, context, trait)
# downstream). * cvFolds a list of vectors: each element defines one fold's
# SAMPLES, as numeric column indices into `sampleNames` or character sample
# names. * samplePartition a data.frame(Sample, Fold): used as-is. A list-form
# `cvFolds` and an explicit `samplePartition` are mutually exclusive. When
# `sampleNames` is supplied, the resolved partition is validated to reference
# only known samples, assign each sample to one fold, and cover every
# sample (a proper partition). Returns list(samplePartition, nFolds).
.normalizeCvFolds <- function(
    cvFolds = 0,
    samplePartition = NULL,
    sampleNames = NULL
) {
    isListFolds <- is.list(cvFolds) && !is.data.frame(cvFolds)
    if (isListFolds && !is.null(samplePartition)) {
        msg <- glue(
            "Provide either a list-form `cvFolds` or an explicit ",
            "`samplePartition`, not both."
        )
        abort(msg)
    }
    if (!is.null(samplePartition)) {
        df <- .validateFoldPartition(
            as_tibble(samplePartition),
            sampleNames
        )
        return(list(samplePartition = df, nFolds = n_distinct(df$Fold)))
    }
    if (isListFolds) {
        return(.twasListFoldsToPartition(cvFolds, sampleNames))
    }
    k <- suppressWarnings(as.integer(cvFolds))
    if (length(k) != 1L || is.na(k)) {
        msg <- glue(
            "`cvFolds` must be a single integer, a list of fold vectors, or ",
            "paired with `samplePartition`."
        )
        abort(msg)
    }
    list(samplePartition = NULL, nFolds = k)
}

# Convert a list-form `cvFolds` (>= 2 fold vectors) to a validated fold
# partition data.frame, resolving numeric column indices via `sampleNames`.
# @noRd
.twasListFoldsToPartition <- function(cvFolds, sampleNames) {
    if (length(cvFolds) < 2L) {
        abort("A list-form `cvFolds` must define at least 2 folds.")
    }
    rows <- map(
        seq_along(cvFolds),
        .twasFoldRowAt,
        cvFolds = cvFolds,
        sampleNames = sampleNames
    )
    df <- .validateFoldPartition(bind_rows(rows), sampleNames)
    list(samplePartition = df, nFolds = length(cvFolds))
}

# One fold's Sample/Fold rows. Numeric ids are resolved to sample names via
# `sampleNames` (required + range-checked); character ids are used as-is.
# @noRd
.twasFoldRow <- function(k, ids, sampleNames) {
    if (is.numeric(ids)) {
        if (is.null(sampleNames)) {
            msg <- glue(
                "Numeric fold vectors require `sampleNames` to resolve ",
                "column indices."
            )
            abort(msg)
        }
        if (any(ids < 1L | ids > length(sampleNames))) {
            msg <- glue(
                "Fold {k} has out-of-range sample column index/indices."
            )
            abort(msg)
        }
        ids <- sampleNames[as.integer(ids)]
    } else {
        ids <- as.character(ids)
    }
    tibble(Sample = ids, Fold = k)
}

# Identify non-zero-variance columns of X. Returns a logical vector.
#' @importFrom matrixStats colSds
#' @noRd
.nonzeroVarColumns <- function(X) {
    sds <- colSds(X, na.rm = TRUE)
    !is.na(sds) & sds != 0
}

# Embed a smaller weights matrix into a full-sized zero matrix matching X and Y
# dimensions.
# @param weightsMatrix The fitted weights (nrow = number of valid columns).
# @param validColumns Logical or character vector identifying which columns of
# X were used.
# @param XColnames Column names of the original X.
# @param YColnames Column names of Y.
# @noRd
.embedWeights <- function(
    weightsMatrix,
    validColumns,
    nColsX,
    nColsY,
    XColnames = NULL,
    YColnames = NULL
) {
    full <- matrix(0, nrow = nColsX, ncol = nColsY)
    if (!is.null(XColnames)) {
        rownames(full) <- XColnames
    }
    if (!is.null(YColnames)) {
        colnames(full) <- YColnames
    }
    full[validColumns, ] <- weightsMatrix
    full
}

# Resolve the incoming susie / susieInf fits for the weight run: prefer a fit
# carried on the method args, else one from `fittedModels`, tagging each with
# its fine-mapping class. Returns list(susieFit, susieInfFit, hasSusie,
# hasSusieInf).
# @noRd
.twasResolveSusieFits <- function(weightMethods, fittedModels) {
    hasSusie <- !is.null(weightMethods[["susie_weights"]])
    hasSusieInf <- !is.null(weightMethods[["susie_inf_weights"]])
    susieFit <- if (hasSusie) {
        weightMethods[["susie_weights"]][["susieFit"]]
    } else {
        NULL
    }
    susieInfFit <- if (hasSusieInf) {
        weightMethods[["susie_inf_weights"]][["susieInfFit"]]
    } else {
        NULL
    }
    if (is.null(susieFit)) {
        susieFit <- fittedModels[["susie"]]
    }
    if (is.null(susieInfFit)) {
        susieInfFit <- fittedModels[["susieInf"]]
    }
    if (!is.null(susieFit)) {
        susieFit <- .setFinemappingFitClass(susieFit, "susie")
    }
    if (!is.null(susieInfFit)) {
        susieInfFit <- .setFinemappingFitClass(susieInfFit, "susieInf")
    }
    list(
        susieFit = susieFit,
        susieInfFit = susieInfFit,
        hasSusie = hasSusie,
        hasSusieInf = hasSusieInf
    )
}

# Chained susieInf -> susie fit for the single-outcome case where both methods
# are requested and no fit was supplied. Returns list(susieFit, susieInfFit).
# @noRd
.twasChainedSusieFit <- function(X, Y, weightMethods) {
    fitArgNames <- c("susieFit", "susieInfFit", "retainFit")
    susieArgs <- weightMethods[["susie_weights"]][setdiff(
        names(weightMethods[["susie_weights"]]),
        fitArgNames
    )]
    # modifyList (not list_modify): a NULL user arg should UNSET the key here.
    susieInfArgs <- modifyList(
        list(convergence_method = "pip"),
        weightMethods[["susie_inf_weights"]][setdiff(
            names(weightMethods[["susie_inf_weights"]]),
            fitArgNames
        )]
    )
    fits <- fitSusieInfThenSusie(
        X,
        Y[, 1],
        args = susieArgs,
        susieInfArgs = susieInfArgs,
        fittedModels = list(susie = NULL, susieInf = NULL)
    )
    list(susieFit = fits[["susie"]], susieInfFit = fits[["susieInf"]])
}

# Write the resolved susie / susieInf fits back onto the method args, deriving
# the susie args from a susieInf fit when susie was requested without its own.
# @noRd
.twasWriteBackSusieFits <- function(weightMethods, r) {
    if (!is.null(r$susieInfFit) && r$hasSusieInf) {
        weightMethods[["susie_inf_weights"]][["susieInfFit"]] <- r$susieInfFit
    }
    if (!is.null(r$susieFit) && r$hasSusie) {
        weightMethods[["susie_weights"]][["susieFit"]] <- r$susieFit
    }
    if (
        r$hasSusie &&
            is.null(weightMethods[["susie_weights"]][["susieFit"]]) &&
            !is.null(r$susieInfFit)
    ) {
        weightMethods[["susie_weights"]] <- prepareSusieFromInfArgs(
            weightMethods[["susie_weights"]],
            r$susieInfFit
        )
    }
    weightMethods
}

.prepareSusieWeightMethods <- function(
    X,
    Y,
    weightMethods,
    fittedModels = NULL
) {
    if (is.vector(Y)) {
        Y <- matrix(Y, ncol = 1)
    }
    if (is.null(fittedModels)) {
        fittedModels <- list()
    }
    r <- .twasResolveSusieFits(weightMethods, fittedModels)
    if (
        r$hasSusie &&
            r$hasSusieInf &&
            ncol(Y) == 1 &&
            is.null(r$susieFit) &&
            is.null(r$susieInfFit)
    ) {
        fits <- .twasChainedSusieFit(X, Y, weightMethods)
        r$susieFit <- fits$susieFit
        r$susieInfFit <- fits$susieInfFit
    }
    .twasWriteBackSusieFits(weightMethods, r)
}

# Per-fold TWAS weight fit for the CV engine. `ctx` carries weightMethods,
# multivariateWeightMethods, cvArgs, retainFits, verbose. Weights are keyed by
# the canonical method key; captured fits keep the full method name.
# @noRd
.weightFitFold <- function(Xtr, Ytr, j, ctx) {
    foldWeightMethods <- .prepareSusieWeightMethods(Xtr, Ytr, ctx$weightMethods)
    weights <- list()
    fits <- list()
    for (method in names(foldWeightMethods)) {
        res <- .twasFoldMethodWeights(
            method,
            foldWeightMethods[[method]],
            Xtr,
            Ytr,
            j,
            ctx
        )
        weights[[res$mk]] <- res$W
        fits[[method]] <- res$fit
    }
    list(weights = weights, fits = fits)
}

# Per-fold priors bound to a multivariate fitter's camelCase args (mr.mash
# data-driven matrices / mvsusie reweighted mixture prior for fold `j`).
# @noRd
.twasFoldPriors <- function(args, method, j, cvArgs) {
    if (
        !is.null(cvArgs$data_driven_priorMatricesCv) &&
            is_in(method, c("mrmash_weights", "mrmashWeights"))
    ) {
        args$dataDrivenPriorMatrices <- cvArgs$data_driven_priorMatricesCv[[j]]
    }
    if (
        !is.null(cvArgs$reweightedMixturePriorCv) &&
            is_in(method, c("mvsusie_weights", "mvsusieWeights"))
    ) {
        args$prior_variance <- cvArgs$reweightedMixturePriorCv[[j]]
    }
    args
}

# One fold's multivariate weight fit; returns list(W, fit).
# @noRd
.twasFoldMultivariate <- function(method, fnName, args, Xtr, Ytr, j, ctx) {
    args <- .twasFoldPriors(args, method, j, ctx$cvArgs)
    if (isTRUE(ctx$retainFits) && is_in("retainFit", names(formals(fnName)))) {
        args$retainFit <- TRUE
    }
    callArgs <- c(list(X = Xtr, Y = Ytr), args)
    W <- if (ctx$verbose < 2) {
        .quietEval(exec(fnName, !!!callArgs))
    } else {
        exec(fnName, !!!callArgs)
    }
    capturedFit <- attr(W, "fit")
    attr(W, "fit") <- NULL
    rownames(W) <- colnames(Xtr)
    list(W = W, fit = capturedFit)
}

# One fold's univariate weight fit (per Y column, column-bound); no fit kept.
# @noRd
.twasFoldUnivariate <- function(fnName, args, Xtr, Ytr, ctx) {
    Wcols <- map(
        seq_len(ncol(Ytr)),
        .twasFitColWeight,
        ctx = ctx,
        fnName = fnName,
        Xtr = Xtr,
        Ytr = Ytr,
        args = args
    )
    W <- exec(cbind, !!!Wcols)
    rownames(W) <- colnames(Xtr)
    list(W = W, fit = NULL)
}

# Fit one method for one CV fold: dispatch to the multivariate or univariate
# path. Returns list(mk, W, fit) keyed by the canonical method key `mk`.
# @noRd
.twasFoldMethodWeights <- function(method, args, Xtr, Ytr, j, ctx) {
    fnName <- .resolveMethodFunction(method, args)
    mk <- str_remove(method, "_weights$|Weights$")
    fit <- if (is_in(method, ctx$multivariateWeightMethods)) {
        .twasFoldMultivariate(method, fnName, args, Xtr, Ytr, j, ctx)
    } else {
        .twasFoldUnivariate(fnName, args, Xtr, Ytr, ctx)
    }
    list(mk = mk, W = fit$W, fit = fit$fit)
}

#' Cross-Validation for weights selection in Transcriptome-Wide Association
#' Studies (TWAS)
#'
#' Performs cross-validation for TWAS, supporting both univariate and
#' multivariate methods. It can either create folds for cross-validation or use
#' pre-defined sample partitions. For multivariate methods, it applies the
#' method to the entire Y matrix for each fold.
#'
#' @param X A matrix of samples by features, where each row represents a sample
#'   and each column a feature.
#' @param Y A matrix (or vector, which will be converted to a matrix) of samples
#'   by outcomes, where each row corresponds to a sample.
#' @param fold An optional integer specifying the number of folds for
#'   cross-validation. If NULL, 'samplePartitions' must be provided.
#' @param samplePartitions An optional dataframe with predefined sample
#'   partitions, containing columns 'Sample' (sample names) and 'Fold' (fold
#'   number). If NULL, 'fold' must be provided.
#' @param weightMethods A list of methods and their specific arguments,
#'   formatted as list(method1 = method1_args, method2 = method2_args), or
#'   alternatively a character vector of method names (eg, c("susie_weights",
#'   "enet_weights")) in which case default arguments will be used for all
#'   methods. methods in the list can be either univariate (applied to each
#'   column of Y) or multivariate (applied to the entire Y matrix).
#' @param maxNumVariants An optional integer to set the randomly selected
#'   maximum number of variants to use for CV purpose, to save computing time.
#' @param variantsToKeep An optional integer to ensure that the listed variants
#'   are kept in the CV when there is a limit on the maxNumVariants to use.
#' @param numThreads The number of threads to use for parallel processing. If
#'   set to -1, the function uses all available cores. If set to 0 or 1, no
#'   parallel processing is performed. If set to 2 or more, parallel processing
#'   is enabled with that many threads.
#' @param verbose Integer controlling verbosity level: 0 = suppress all
#'   messages, 1 = suppress external package messages (default), 2 = show all
#'   messages including those from external packages.
#' @param retainFits Logical. Retain the per-fold / per-method fitted-model
#'   objects on the result. Default \code{FALSE}.
#' @param seed Integer or \code{NULL}. When supplied, seeds both the
#'   main-process RNG (fold partitioning, variant sub-sampling) via
#'   \code{set.seed} and the parallel fold-fitting RNG via the
#'   \code{BiocParallel} \code{RNGseed}, so results are reproducible even under
#'   multi-threading. \code{NULL} (default) leaves the session RNG untouched and
#'   uses the historical parallel default.
#' @param ... Additional arguments forwarded to the per-method weight learners.
#' @return A list with the following components:
#' \itemize{
#'   \item `samplePartition`: A dataframe showing the sample partitioning used
#'   in the cross-validation.
#'   \item `prediction`: A list of matrices with predicted Y values for each
#'   method and fold.
#'   \item `metrics`: A matrix with rows representing methods and columns for
#'   various metrics:
#'     \itemize{
#'       \item `corr`: Pearson's correlation between predicated and observed
#'       values.
#'       \item `adj_rsq`: Adjusted R-squared value (which indicates the
#'       proportion of variance explained by the model) that accounts for the
#'       number of predictors in the model.
#'       \item `pval`: P-value assessing the significance of the model's
#'       predictions.
#'       \item `RMSE`: Root Mean Squared Error, a measure of the model's
#'       prediction error.
#'       \item `MAE`: Mean Absolute Error, a measure of the average magnitude
#'       of errors in a set of predictions.
#'     }
#'   \item `timeElapsed`: The time taken to complete the cross-validation
#'   process.
#' }
#' @importFrom purrr map
#' @importFrom BiocParallel bplapply bpworkers MulticoreParam
#' @importFrom quadprog solve.QP
#' @examples
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:80]
#' Y <- multiTraitData$Y
#' twasWeightsCv(X, Y[, 1, drop = FALSE], fold = 3,
#'   weightMethods = list(susie_weights = list()))
#' @export
twasWeightsCv <- function(
    X,
    Y,
    fold = NULL,
    samplePartitions = NULL,
    weightMethods = NULL,
    maxNumVariants = NULL,
    variantsToKeep = NULL,
    numThreads = 1,
    verbose = 1,
    retainFits = FALSE,
    seed = NULL,
    ...
) {
    p <- as.list(environment())
    p$cvArgs <- list(...)
    .twasWeightsCvImpl(p)
}

# Multivariate weight methods (snake + camel) fit on the whole Y for a fold;
# univariate methods are fit per Y column. fSuSiE is intentionally absent -- it
# is functional and cannot be refit from a bare (X, y) fold split, so its
# cross-validated predictions are supplied by fineMappingPipeline.
# @noRd
.twasCvMultivariateMethods <- c(
    "mrmash_weights",
    "mvsusie_weights",
    "mrmashWeights",
    "mvsusieWeights"
)

# twasWeightsCv worker. `p` is the captured public arguments plus `cvArgs`
# (the `...`). With no weightMethods the caller only wants the fold partition.
# @noRd
.twasWeightsCvImpl <- function(p) {
    weightMethods <- if (is.character(p$weightMethods)) {
        .twasMethodLookup(p$weightMethods)
    } else {
        p$weightMethods
    }
    if (is.null(p$seed) && !exists(".Random.seed") && p$verbose >= 1) {
        inform(
            "! No seed set. Pass `seed=` or call set.seed() for reproducibility."
        )
    }
    if (is.null(weightMethods)) {
        res <- .crossValidateWeights(
            p$X,
            p$Y,
            fold = p$fold,
            samplePartitions = p$samplePartitions,
            fitFold = .cvNoopFitFold,
            numThreads = p$numThreads,
            maxNumVariants = p$maxNumVariants,
            variantsToKeep = p$variantsToKeep,
            retainFits = p$retainFits,
            verbose = p$verbose,
            seed = p$seed
        )
        return(list(samplePartition = res$samplePartition))
    }
    cvFitCtx <- list(
        weightMethods = weightMethods,
        multivariateWeightMethods = .twasCvMultivariateMethods,
        cvArgs = p$cvArgs,
        retainFits = p$retainFits,
        verbose = p$verbose
    )
    .crossValidateWeights(
        p$X,
        p$Y,
        fold = p$fold,
        samplePartitions = p$samplePartitions,
        fitFold = .weightFitFold,
        fitFoldCtx = cvFitCtx,
        numThreads = p$numThreads,
        maxNumVariants = p$maxNumVariants,
        variantsToKeep = p$variantsToKeep,
        retainFits = p$retainFits,
        verbose = p$verbose,
        seed = p$seed
    )
}

# Fit one TWAS weight method by name against the filtered design matrix,
# embedding the fitted weights back into the full variant space. `ctx` carries
# the shared fit state (X, Y, Xfiltered, validColumns, retainFits,
# retainFitDetail, verbose).
# @noRd
.computeMethodWeights <- function(methodName, weightMethods, ctx) {
    shortName <- str_remove(methodName, "_weights$")
    if (ctx$verbose >= 1) {
        msg <- glue("  Fitting {shortName} ...")
        inform(msg)
        tic()
    }
    args <- weightMethods[[methodName]]
    fnName <- .resolveMethodFunction(methodName, args)
    args <- .twasApplyRetainFit(
        args,
        fnName,
        ctx$retainFits,
        ctx$retainFitDetail
    )
    fit <- .twasFitWeightsMatrix(fnName, args, ctx, methodName)
    result <- .embedWeights(
        fit$weights,
        ctx$validColumns,
        ncol(ctx$X),
        ncol(ctx$Y),
        colnames(ctx$X),
        colnames(ctx$Y)
    )
    if (!is.null(fit$methodFit)) {
        attr(result, "fit") <- fit$methodFit
    }
    if (ctx$verbose >= 1) {
        elapsed <- toc(quiet = TRUE)
        secs <- sprintf("%.1f", elapsed$toc - elapsed$tic)
        msg <- glue("  Fitting {shortName} done in {secs}s")
        inform(msg)
    }
    result
}

# Multivariate weight methods (variants x features matrix), accepting both
# snake_case and camelCase keys. fSuSiE is multivariate but never refit here --
# fsusieWeights extracts from the supplied fsusieFit.
# @noRd
.twasMultivariateWeightMethods <- c(
    "mrmash_weights",
    "mvsusie_weights",
    "fsusie_weights",
    "mrmashWeights",
    "mvsusieWeights",
    "fsusieWeights"
)

# Add retainFit (or its legacy retain_fit alias) + fitDetail to `args` when the
# target weight function accepts them and fits are being retained.
# @noRd
.twasApplyRetainFit <- function(args, fnName, retainFits, retainFitDetail) {
    if (!retainFits) {
        return(args)
    }
    fnFormals <- names(formals(fnName))
    if (is_in("retainFit", fnFormals)) {
        args$retainFit <- TRUE
    } else if (is_in("retain_fit", fnFormals)) {
        args$retain_fit <- TRUE
    }
    if (is_in("fitDetail", fnFormals) && is.null(args$fitDetail)) {
        args$fitDetail <- retainFitDetail
    }
    args
}

# Dispatch weight fitting to the multivariate or per-column univariate path;
# returns list(weights, methodFit).
# @noRd
.twasFitWeightsMatrix <- function(fnName, args, ctx, methodName) {
    if (is_in(methodName, .twasMultivariateWeightMethods)) {
        .twasFitMultivariate(fnName, args, ctx)
    } else {
        .twasFitUnivariate(fnName, args, ctx)
    }
}

# Multivariate fit: one call producing the full variants x features matrix.
# @noRd
.twasFitMultivariate <- function(fnName, args, ctx) {
    call <- c(list(X = ctx$Xfiltered, Y = ctx$Y), args)
    weightsMatrix <- if (ctx$verbose < 2) {
        .quietEval(exec(fnName, !!!call))
    } else {
        exec(fnName, !!!call)
    }
    methodFit <- if (ctx$retainFits) attr(weightsMatrix, "fit") else NULL
    if (nrow(weightsMatrix) != length(ctx$validColumns)) {
        weightsMatrix <- weightsMatrix[names(ctx$validColumns), , drop = FALSE]
    }
    list(weights = weightsMatrix, methodFit = methodFit)
}

# Univariate fit: apply the method to each column of Y, filling a zero-init
# weights matrix.
# @noRd
.twasFitUnivariate <- function(fnName, args, ctx) {
    weightsMatrix <- matrix(0, nrow = ncol(ctx$Xfiltered), ncol = ncol(ctx$Y))
    methodFit <- NULL
    for (k in seq_len(ncol(ctx$Y))) {
        call <- c(list(X = ctx$Xfiltered, y = ctx$Y[, k]), args)
        weightsVector <- if (ctx$verbose < 2) {
            .quietEval(exec(fnName, !!!call))
        } else {
            exec(fnName, !!!call)
        }
        if (ctx$retainFits && is.null(methodFit)) {
            methodFit <- attr(weightsVector, "fit")
        }
        if (is.matrix(weightsVector)) {
            weightsVector <- weightsVector[, k]
        }
        weightsMatrix[, k] <- weightsVector
    }
    list(weights = weightsMatrix, methodFit = methodFit)
}

# Assemble the (study, context, trait, method, entry) row vectors for the
# TwasWeights collection from the fitted `weightsList`. `ctx` carries the shared
# identity + flags (study, context, trait, Y, retainFits, standardized,
# dataType).
# @noRd
.buildTwasWeightEntries <- function(weightsList, variantIds, ctx) {
    rows <- list_flatten(map(
        names(weightsList),
        .twasMethodRowsFor,
        weightsList = weightsList,
        variantIds = variantIds,
        ctx = ctx
    ))
    list(
        study = map_chr(rows, "study"),
        context = map_chr(rows, "context"),
        trait = map_chr(rows, "trait"),
        method = map_chr(rows, "method"),
        entry = map(rows, "entry")
    )
}

# One TwasWeightsEntry for a (variantIds, weights) pair with the shared flags.
# @noRd
.twasEntry <- function(variantIds, weights, fits, ctx) {
    TwasWeightsEntry(
        variantIds = variantIds,
        weights = weights,
        fits = fits,
        cvResult = NULL,
        standardized = isTRUE(ctx$standardized),
        dataType = ctx$dataType
    )
}

# Row-records for one fitted method. When trait/context were supplied per-row
# (length == ncol(Y)) emit one (method, outcome) row per Y column; otherwise a
# single row carrying the (possibly multi-column) weights matrix as-is.
# @noRd
.twasMethodRows <- function(m, wMat, variantIds, ctx) {
    fitVal <- attr(wMat, "fit")
    attr(wMat, "fit") <- NULL
    fits <- if (ctx$retainFits) fitVal else NULL
    shortMethod <- str_remove(m, "(_weights|Weights)$")
    nY <- ncol(ctx$Y)
    perOutcome <- length(ctx$trait) == nY &&
        is_in(length(ctx$context), c(1L, nY))
    if (!perOutcome) {
        wPayload <- if (ncol(wMat) == 1L) drop(wMat) else wMat
        return(list(list(
            study = ctx$study[1L],
            context = ctx$context[1L],
            trait = ctx$trait[1L],
            method = shortMethod,
            entry = .twasEntry(variantIds, wPayload, fits, ctx)
        )))
    }
    contextV <- if (length(ctx$context) == 1L) {
        rep(ctx$context, nY)
    } else {
        ctx$context
    }
    studyV <- if (length(ctx$study) == 1L) rep(ctx$study, nY) else ctx$study
    map(
        seq_len(nY),
        .twasMethodRowAt,
        studyV = studyV,
        contextV = contextV,
        ctx = ctx,
        shortMethod = shortMethod,
        variantIds = variantIds,
        wMat = wMat,
        fits = fits
    )
}

#' Run multiple TWAS weight methods
#'
#' Applies specified weight methods to the datasets X and Y, returning weight
#' matrices for each method. Handles both univariate and multivariate methods,
#' and filters out columns in X with zero standard error. This function utilizes
#' parallel processing to handle multiple methods.
#'
#' @param X A matrix of samples by features, where each row represents a sample
#'   and each column a feature.
#' @param Y A matrix (or vector, which will be converted to a matrix) of samples
#'   by outcomes, where each row corresponds to a sample.
#' @param weightMethods A list of methods and their specific arguments,
#'   formatted as list(method1 = method1_args, method2 = method2_args), or
#'   alternatively a character vector of method names (eg, c("susie_weights",
#'   "enet_weights")) in which case default arguments will be used for all
#'   methods. methods in the list can be either univariate (applied to each
#'   column of Y) or multivariate (applied to the entire Y matrix).
#' @param numThreads The number of threads to use for parallel processing. If
#'   set to -1, the function uses all available cores. If set to 0 or 1, no
#'   parallel processing is performed. If set to 2 or more, parallel processing
#'   is enabled with that many threads.
#' @param fittedModels Optional named list of fitted SuSiE-family models.
#' @param retainFits If TRUE, retain fitted model objects as attributes on
#'   returned weight matrices when supported by the weight method.
#' @param verbose Integer controlling verbosity level: 0 = suppress all
#'   messages, 1 = suppress external package messages (default), 2 = show all
#'   messages including those from external packages.
#' @param study Character. Study identity label recorded on the resulting
#'   weights.
#' @param context Character. Context identity label recorded on the resulting
#'   weights.
#' @param trait Character. Trait identity label recorded on the resulting
#'   weights.
#' @param retainFitDetail Character. Level of fit detail to retain:
#'   \code{"slim"} (default) or \code{"full"}.
#' @param standardized Logical. Whether the supplied \code{X} / \code{Y} are
#'   already standardized. Default \code{FALSE}.
#' @param dataType Character or \code{NULL}. Data-type label recorded on the
#'   weights (e.g. \code{"individual"}).
#' @param ldSketch A \code{GenotypeHandle} LD sketch to record on the weights,
#'   or \code{NULL}.
#' @param seed Integer or \code{NULL}. When supplied, seeds the main-process RNG
#'   via \code{set.seed} and the parallel method-fitting RNG via the
#'   \code{BiocParallel} \code{RNGseed}, for reproducibility under
#'   multi-threading. \code{NULL} (default) leaves the session RNG untouched.
#' @return A list where each element is named after a method and contains the
#'   weight matrix produced by that method.
#'
#' @examples
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:80]
#' Y <- multiTraitData$Y
#' learnTwasWeights(X, Y[, 1, drop = FALSE],
#'   weightMethods = list(susie_weights = list()))
#' @export
#' @importFrom purrr map exec
#' @importFrom rlang !!! abort warn inform arg_match cnd_signal .data
#' @importFrom glue glue
#' @importFrom tictoc tic toc
learnTwasWeights <- function(
    X,
    Y,
    weightMethods,
    study = "",
    context = "",
    trait = "",
    numThreads = 1,
    fittedModels = NULL,
    retainFits = FALSE,
    retainFitDetail = c("slim", "full"),
    standardized = FALSE,
    dataType = NULL,
    ldSketch = NULL,
    verbose = 1,
    seed = NULL
) {
    .learnTwasWeightsImpl(as.list(environment()))
}

# Validate X/Y shapes; coerce a vector Y to a one-column matrix. Returns Y.
# @noRd
.twasValidateXY <- function(X, Y) {
    if (!is.matrix(X) || (!is.matrix(Y) && !is.vector(Y))) {
        abort("X must be a matrix and Y must be a matrix or a vector.")
    }
    if (is.vector(Y)) {
        Y <- matrix(Y, ncol = 1)
    }
    if (nrow(X) != nrow(Y)) {
        abort("The number of rows in X and Y must be the same.")
    }
    Y
}

# Number of BiocParallel workers to use: -1 means all available, otherwise the
# requested count capped at what is available.
# @noRd
.twasResolveCores <- function(numThreads) {
    avail <- bpworkers(MulticoreParam())
    min(if (numThreads == -1) avail else numThreads, avail)
}

# Variant ids for the weight rows: colnames(X), or synthetic variant_i labels.
# @noRd
.twasVariantIds <- function(X) {
    if (!is.null(colnames(X))) {
        colnames(X)
    } else {
        str_c("variant_", seq_len(ncol(X)))
    }
}

# Fit every weight method (parallel when >= 2 cores, else serial map), keyed by
# method name.
# @noRd
.twasFitAllMethods <- function(weightMethods, ctx, numCores) {
    weightsList <- if (numCores >= 2) {
        bpParam <- .bpSeedParam(numCores, ctx$rngSeed)
        bplapply(
            names(weightMethods),
            .computeMethodWeights,
            weightMethods,
            ctx,
            BPPARAM = bpParam
        )
    } else {
        map(names(weightMethods), .computeMethodWeights, weightMethods, ctx)
    }
    names(weightsList) <- names(weightMethods)
    weightsList
}

# Set weight-matrix rownames to colnames(X), preserving any retained `fit` attr.
# @noRd
.twasApplyRownames <- function(weightsList, X) {
    if (is.null(colnames(X))) {
        return(weightsList)
    }
    map(weightsList, .twasSetRownames, X = X)
}

# learnTwasWeights worker: validate, resolve methods, fit each, and assemble the
# TwasWeights collection. `p` is the captured public arguments.
# @noRd
.learnTwasWeightsImpl <- function(p) {
    # Seed the main-process RNG (used by serial method fitting); the parallel
    # path is seeded via ctx$rngSeed in .twasFitAllMethods.
    if (!is.null(p$seed)) {
        set.seed(p$seed)
    }
    retainFitDetail <- p$retainFitDetail
    retainFitDetail <- arg_match(retainFitDetail, c("slim", "full"))
    Y <- .twasValidateXY(p$X, p$Y)
    weightMethods <- if (is.character(p$weightMethods)) {
        .twasMethodLookup(p$weightMethods)
    } else {
        p$weightMethods
    }
    validColumns <- .nonzeroVarColumns(p$X)
    Xfiltered <- as.matrix(p$X[, validColumns, drop = FALSE])
    weightMethods <- .prepareSusieWeightMethods(
        Xfiltered,
        Y,
        weightMethods,
        p$fittedModels
    )
    ctx <- list(
        X = p$X,
        Y = Y,
        Xfiltered = Xfiltered,
        validColumns = validColumns,
        study = p$study,
        context = p$context,
        trait = p$trait,
        retainFits = p$retainFits,
        retainFitDetail = retainFitDetail,
        standardized = p$standardized,
        dataType = p$dataType,
        verbose = p$verbose,
        rngSeed = p$seed
    )
    weightsList <- .twasFitAllMethods(
        weightMethods,
        ctx,
        .twasResolveCores(p$numThreads)
    )
    weightsList <- .twasApplyRownames(weightsList, p$X)
    rows <- .buildTwasWeightEntries(weightsList, .twasVariantIds(p$X), ctx)
    TwasWeights(
        study = rows$study,
        context = rows$context,
        trait = rows$trait,
        method = rows$method,
        entry = rows$entry,
        ldSketch = p$ldSketch
    )
}

#' Predict outcomes using TWAS weights
#'
#' This function takes a matrix of predictors (\code{X}) and a list of TWAS
#' (transcriptome-wide association studies) weights (\code{weightsList}), and
#' calculates the predicted outcomes by multiplying \code{X} by each set of
#' weights in \code{weightsList}. The names of the elements in the output list
#' are derived from the names in \code{weightsList}, with "_weights" replaced by
#' "_predicted".
#'
#' @param X A matrix or data frame of predictors where each row is an
#'   observation and each column is a variable.
#' @param weightsList A list of numeric vectors representing the weights for
#'   each predictor. The names of the list elements should follow the pattern
#'   \code{[outcome]_weights}, where \code{[outcome]} is the name of the outcome
#'   variable that the weights are associated with.
#'
#' @return A named list of numeric vectors, where each vector is the predicted
#'   outcome for the corresponding set of weights in \code{weightsList}. The
#'   names of the list elements are derived from the names in \code{weightsList}
#'   by replacing "_weights" with "_predicted".
#'
#' @export
#' @examples
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:4]
#' colnames(X) <- paste0("v", 1:4)
#' twe <- TwasWeightsEntry(variantIds = paste0("v", 1:4),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5), standardized = FALSE)
#' tw <- TwasWeights(study = "s1", context = "brain", trait = "g1",
#'   method = "susie", entry = list(twe))
#' twasPredict(X, tw)
twasPredict <- function(X, weightsList) {
    if (is(weightsList, "TwasWeights")) {
        # Per-row weights vector/matrix payloads. Use the method name as key
        # for compatibility with the legacy snake_case "<method>_predicted"
        # convention; ensembleWeights() rebinds the suffix.
        methodNames <- as.character(weightsList$method)
        wl <- set_names(
            map(weightsList$entry, getWeights),
            str_c(methodNames, "_weights")
        )
    } else {
        wl <- weightsList
    }
    set_names(
        map(wl, .twasPredictOne, X = X),
        .renameSuffix(names(wl), "predicted")
    )
}

#' Estimate Sparsity from mr.ash Mixture Proportions
#'
#' Computes an empirical estimate of the proportion of non-zero effects
#' (sparsity) from the mr.ash fit. mr.ash fits a mixture model with a point mass
#' at zero (spike) plus continuous components (slab), and learns the mixture
#' proportions via variational EM. The sparsity estimate \code{1 - pi[1]} is the
#' empirical Bayes estimate of the non-null proportion, which can be used as a
#' data-driven prior for the inclusion probability parameters (\code{pi} for
#' bayesC, \code{probIn} for BayesB) of spike-and-slab Bayesian methods.
#'
#' @param weightResults Named list of weight vectors or matrices as returned by
#'   \code{\link{learnTwasWeights}}. The mr.ash element should have a
#'   \code{"fit"} attribute containing the model fit object (set
#'   \code{retainFits = TRUE} in \code{learnTwasWeights} to obtain this).
#'
#' @return A scalar sparsity estimate (proportion of non-zero effects).
#' @examples
#' estimateSparsity(list(mrash_weights = structure(c(0.1, 0, 0.3),
#'   fit = list(pi = c(0.6, 0.2, 0.2)))))
#' @export
estimateSparsity <- function(weightResults) {
    if (is(weightResults, "TwasWeights")) {
        # Method names on the new TwasWeights collection are bare tokens
        # ("mrash"), not the snake_case _weights suffix form.
        methods <- as.character(weightResults$method)
        idx <- which(methods == "mrash")
        if (length(idx) == 0L) {
            msg <- glue(
                "mr.ash entry not found in TwasWeights. Run ",
                "learnTwasWeights() ",
                "with retainFits = TRUE and ensure 'mrash' is in the ",
                "method list."
            )
            abort(msg)
        }
        fit <- getFits(weightResults$entry[[idx[[1L]]]])
        if (is.null(fit) || is.null(fit$pi)) {
            msg <- glue(
                "mr.ash fit object not found. Run learnTwasWeights() with ",
                "retainFits = TRUE ",
                "and ensure mrash_weights is included."
            )
            abort(msg)
        }
    } else {
        w <- weightResults[["mrash_weights"]]
        if (is.null(w)) {
            abort(
                "mr.ash weights ('mrash_weights') not found in weightResults."
            )
        }
        fit <- attr(w, "fit")
        if (is.null(fit) || is.null(fit$pi)) {
            msg <- glue(
                "mr.ash fit object not found. Run learnTwasWeights() with ",
                "retainFits = TRUE ",
                "and ensure mrash_weights is included."
            )
            abort(msg)
        }
    }

    # fit$pi[1] is the weight on the spike (sa2[1] = 0); 1 - pi[1] = non-null
    # proportion
    return(1 - fit$pi[1])
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# TRUE when joint column `jc` of `object` is not a character vector.
# @noRd
.twasColNotCharacter <- function(jc, object) {
    !is.character(object[[jc]])
}

# Validation message for a non-character joint column `jc`.
# @noRd
.twasBadColMsg <- function(jc, object) {
    glue(
        "'{jc}' column must be character (got {class(object[[jc]])[[1L]]})"
    )
}

# Column `cn` of `object`, extracted via `[[` (preserves the required `entry`
# column that `object[, cn]` would drop).
# @noRd
.twasColumn <- function(cn, object) {
    object[[cn]]
}

# Canonical short method name: map a full function name back to its short name.
# @noRd
.twasCanonicalShortName <- function(m, fnToShort) {
    if (is_in(m, names(fnToShort))) fnToShort[[m]] else m
}

# The per-method args list for short name `m`, tagged with its `impl` attribute.
# @noRd
.twasMethodArgsWithImpl <- function(m) {
    args <- .twasMethodMap[[m]]$args
    attr(args, "impl") <- .twasMethodMap[[m]]$impl
    args
}

# The implementation function name for short method name `m`.
# @noRd
.twasMethodFn <- function(m) {
    .twasMethodMap[[m]]$fn
}

# The Sample/Fold rows for fold `k` of a list-form `cvFolds`.
# @noRd
.twasFoldRowAt <- function(k, cvFolds, sampleNames) {
    .twasFoldRow(k, cvFolds[[k]], sampleNames)
}

# One univariate fold's weight column for outcome `k` (quiet unless verbose).
# @noRd
.twasFitColWeight <- function(k, ctx, fnName, Xtr, Ytr, args) {
    callArgs <- c(list(X = Xtr, y = Ytr[, k]), args)
    w <- if (ctx$verbose < 2) {
        .quietEval(exec(fnName, !!!callArgs))
    } else {
        exec(fnName, !!!callArgs)
    }
    as.numeric(w)
}

# The row-records for method `m` in `weightsList` (keyed lookup + build).
# @noRd
.twasMethodRowsFor <- function(m, weightsList, variantIds, ctx) {
    .twasMethodRows(m, weightsList[[m]], variantIds, ctx)
}

# One (method, outcome) row-record for Y column `k`.
# @noRd
.twasMethodRowAt <- function(
    k,
    studyV,
    contextV,
    ctx,
    shortMethod,
    variantIds,
    wMat,
    fits
) {
    list(
        study = studyV[k],
        context = contextV[k],
        trait = ctx$trait[k],
        method = shortMethod,
        entry = .twasEntry(variantIds, wMat[, k], fits, ctx)
    )
}

# Set one weight matrix's rownames to colnames(X), preserving any `fit` attr.
# @noRd
.twasSetRownames <- function(x, X) {
    fit <- attr(x, "fit")
    rownames(x) <- colnames(X)
    if (!is.null(fit)) {
        attr(x, "fit") <- fit
    }
    x
}

# X %*% w for one weight vector/matrix (coerced to a 1-column matrix if needed).
# @noRd
.twasPredictOne <- function(w, X) {
    if (!is.matrix(w)) {
        w <- matrix(w, ncol = 1)
    }
    X %*% w
}
