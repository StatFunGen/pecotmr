# Check a token's arguments against the union of formals along the chain this
# run will actually take. Stricter than the constructor, which cannot know the
# input class: every summary-statistics chain ends in one of pecotmr's own
# solvers and enumerates its formals, so the RSS path is checkable even where
# the individual path is not (ncvreg, glmnet and RcppDPR all take `...`).
# @noRd
.twasCheckMethodArgsForInput <- function(methodArgs, inputKind) {
    for (token in names(methodArgs)) {
        accepted <- .twasChainAccepted(token, inputKind)
        if (is.null(accepted)) {
            next
        }
        .engineCheckExtra(
            methodArgs[[token]],
            accepted,
            glue("twasWeightsPipeline: method '{token}'"),
            str_flatten(.twasMethodChainFor(token, inputKind), " -> ")
        )
    }
    invisible(NULL)
}

# =============================================================================
# TwasWeights S4 class
# -----------------------------------------------------------------------------
# DFrame-subclass collection keyed by the identity tuple (study, context,
# trait, method). Each row holds a TwasWeightsRow payload (variant ids
# + per-variant weight vector/matrix). Class-level slots:
#   * ldSketch   GenotypeHandle (NULL for individual-level fits, the
#                LD-sketch handle for RSS-derived weights).
# Constructor + accessors below. The twasWeights pipeline helpers
# (learnTwasWeights, CV, ensemble, etc.) follow at the bottom.
# =============================================================================

#' @include AllGenerics.R tupleSelectors.R
NULL

#' @title TWAS Weights Collection
#' @description S4 collection of TWAS weights keyed by the identity tuple
#'   \code{(study, context, trait, method)}. Each entry is a
#'   \code{TwasWeightsRow} carrying one method's weights for one
#'   trait/context/study. Implements the \code{DFrame}-subclass collection
#'   pattern.
#'
#' Required columns: \code{study}, \code{context}, \code{trait}, \code{method},
#' \code{entry}. Each \code{entry} is a \code{TwasWeightsRow}.
#'
#' Optional columns \code{jointStudies}, \code{jointContexts},
#' \code{jointTraits} appear when the collection contains rows produced by a
#' \code{jointSpecification}-driven joint fit. For such a row, the corresponding
#' identity-tuple column carries the sentinel \code{"joint"} and the joint
#' column lists the semicolon-joined members of the joined axis. For non-joint
#' rows the joint columns are \code{NA_character_}. Tuple uniqueness is enforced
#' jointly across the identity-tuple columns and any present joint columns.
#' @slot ldSketch The LD reference genotype panel the weights were
#'   derived against, or \code{NULL} when the weights were learned from
#'   individual-level data. Used downstream for cross-pipeline LD-sketch
#'   identity validation.
#' @export
setClass(
    "TwasWeights",
    contains = "RangedTupleList",
    representation(ldSketch = "LdSketchOrNULL"),
    validity = function(object) .validateTwasWeights(object)
)

# Validity for the TwasWeights collection: required columns, entry payloads,
# region/traitPos provenance, joint* column types, tuple uniqueness, and the
# optional ldSketch. Returns TRUE or a character vector of error messages.
# @noRd
#' @importFrom checkmate makeAssertCollection assertNames
.validateTwasWeights <- function(object) {
    coll <- makeAssertCollection()
    assertNames(
        .tupleColumnNames(object),
        must.include = c("study", "context", "trait", "method"),
        what = "colnames",
        .var.name = "mcols",
        add = coll
    )
    # The checks below read those columns; running them on an object missing
    # them reports the consequence rather than the cause.
    if (!coll$isEmpty()) {
        return(coll$getMessages())
    }
    coll$push(.twasValidateColumns(object))
    coll$getMessages()
}

# Column-level checks that run only once the required columns are present.
# @noRd
.twasValidateColumns <- function(object) {
    jointCols <- intersect(
        c("jointStudies", "jointContexts", "jointTraits"),
        .tupleColumnNames(object)
    )
    c(
        .twasValidateEntries(object),
        .validateTraitPosColumn(object),
        .twasValidateJointCols(object, jointCols),
        .twasValidateKeyUniqueness(object, jointCols)
    )
}

# The per-entry payload columns. The variants and their weights are the
# elements now, so what is left to check is that the payload columns are
# present and parallel to them.
# @noRd
.twasValidateEntries <- function(object) {
    payload <- c("fits", "cvResult", "standardized", "dataType")
    missingCols <- setdiff(payload, .tupleColumnNames(object))
    if (length(missingCols) > 0L) {
        return(str_c(
            "missing entry payload columns: ",
            str_flatten(missingCols, ", ")
        ))
    }
    character()
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

# =============================================================================

#' @title Create a TwasWeights Collection Object
#' @description Construct a \code{TwasWeights} DFrame-subclass collection from
#'   per-tuple vectors and a list of \code{TwasWeightsRow} payloads (one per
#'   tuple).
#' @param study Character vector of study identifiers. Use the sentinel
#'   \code{"joint"} for rows produced by a cross-study joint fit.
#' @param context Character vector of context labels. Use \code{"joint"} for
#'   rows produced by a cross-context joint fit.
#' @param trait Character vector of trait identifiers. Use \code{"joint"} for
#'   rows produced by a cross-trait joint fit.
#' @param method Character vector of TWAS weight method names.
#' @param entry List / \code{SimpleList} of \code{TwasWeightsRow} objects.
#' @param jointStudies Optional character vector (length \code{length(study)})
#'   listing the semicolon-joined studies participating in each row's
#'   cross-study joint fit, or \code{NA_character_} for non-joint rows. When
#'   \code{NULL} (default) the column is omitted.
#' @param jointContexts Optional character vector for cross-context joints. Same
#'   shape as \code{jointStudies}.
#' @param jointTraits Optional character vector for cross-trait joints. Same
#'   shape as \code{jointStudies}.
#' @param ldSketch An optional genotype panel (see
#'   \code{\link{readGenotypes}}), or \code{NULL} for individual-level fits.
#' @param traitPos Optional per-row trait genomic anchor (a \code{GRanges} or
#'   \code{NULL}), carried forward as provenance; not part of the identity key.
#'   \code{NULL} (default) omits the column.
#' @return A \code{TwasWeights} object.
#' @examples
#' twe <- twasWeightsRow(variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5), standardized = FALSE)
#' tw <- TwasWeights(study = "s1", context = "brain", trait = "gene1",
#'   method = "susie", entry = list(twe))
#' tw
#' @importFrom checkmate assertCharacter assert checkList
#' @importFrom checkmate checkClass
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
    traitPos = NULL,
    ldSketch = NULL
) {
    assertCharacter(study, any.missing = FALSE)
    assertCharacter(context, any.missing = FALSE)
    assertCharacter(trait, any.missing = FALSE)
    assertCharacter(method, any.missing = FALSE)
    # `entry` is documented as "List / SimpleList"; SimpleList is S4 and
    # fails checkList, so this must be an or-combination.
    assert(
        checkList(entry),
        checkClass(entry, "SimpleList"),
        .var.name = "entry"
    )
    assertCharacter(jointStudies, null.ok = TRUE)
    assertCharacter(jointContexts, null.ok = TRUE)
    assertCharacter(jointTraits, null.ok = TRUE)
    n <- .twasCheckRowLengths(study, context, trait, method, entry)
    entry <- map(entry, .asTwRowPayload)
    .checkRowPayloads(entry, "TwasWeightsRow", "TWAS-weight")
    baseCols <- c(
        list(
            study = as.character(study),
            context = as.character(context),
            trait = as.character(trait),
            method = as.character(method)
        ),
        .twRowPayloadCols(entry)
    )
    withJoint <- .twasAppendJointCols(
        baseCols,
        jointStudies,
        jointContexts,
        jointTraits,
        n
    )
    cols <- .appendTraitPosCol(withJoint, traitPos, n)
    dfArgs <- c(cols, list(check.names = FALSE))
    # Each entry's variants become one ELEMENT; its weights become that
    # element's inner mcols and the rest of its payload becomes outer mcols.
    # A multi-seqname entry is split into one element per chromosome, with its
    # metadata row replicated alongside.
    split <- .rtlSplitBySeqname(map(entry, rowVariants))
    md <- exec(S4Vectors::DataFrame, !!!dfArgs)
    grl <- `mcols<-`(
        GenomicRanges::GRangesList(split$entry),
        value = md[split$fromIdx, , drop = FALSE]
    )
    obj <- new("TwasWeights", grl, ldSketch = .asLdSketch(ldSketch))
    validObject(obj)
    obj
}

# One entry -> its element: the variants as a GRanges, with the per-variant
# weights carried in the element's own mcols alongside the alleles.
# @noRd

# The non-variant half of each entry, as outer mcols columns. Objects of
# arbitrary shape (fits, cvResult, dataType) ride in SimpleLists.
# @noRd

# @noRd
.twPayloadStandardized <- function(p) {
    getStandardized(p)
}

# Every `entry` must be a TwasWeightsRow. Checked before the payload is
# unpacked, so a wrong type reports itself rather than surfacing as a missing
# method on whatever was passed instead.
# @noRd

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
    supplied <- compact(list(
        jointStudies = jointStudies,
        jointContexts = jointContexts,
        jointTraits = jointTraits
    ))
    walk2(names(supplied), supplied, .twasCheckJointColLength, n = n)
    c(cols, map(supplied, as.character))
}

# Each joint-provenance column must be one value per row.
# @noRd
.twasCheckJointColLength <- function(nm, val, n) {
    if (length(val) == n) {
        return(invisible(NULL))
    }
    msg <- glue("`{nm}` must have the same length as `study`.")
    abort(msg)
}

#' @rdname getRegion
#' @export
setMethod("getRegion", "TwasWeights", function(x) .getRegionColumn(x))

#' @rdname getTraitPosition
#' @export
setMethod("getTraitPosition", "TwasWeights", function(x) {
    .getTraitPosColumn(x)
})

#' @title Get a Single TWAS Weights Entry
#' @description Return the \code{TwasWeightsRow} for one \code{(study,
#'   context, trait, method)} row of a \code{TwasWeights} collection.
#' @param x A \code{TwasWeights} object.
#' @param study,context,trait,method Single character identifiers. All required
#'   when the collection has more than one row; optional when the collection has
#'   a single row.
#' @return A \code{TwasWeightsRow} object.
#' @examples
#' twe <- twasWeightsRow(
#'   variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)), weights = rep(0.1, 4),
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
        x[idx]
    }
)

# The weights of one row, read straight off that element's mcols.
# @noRd
.twasRowWeights <- function(i, x) {
    mcols(x[[i]])$weight
}

# Rebuild the TwasWeightsRow for a tuple from its stored element(s).
#
# The entry is a DERIVED VIEW, not stored state: the variants live in the
# element, the weights in that element's mcols, and the rest of the payload in
# the collection's mcols. A tuple split across chromosomes owns several
# elements, so they are stitched back into the single entry callers expect --
# the same contract getSumStats() keeps.
# @noRd

# The row payload a (study, context, trait, method) selector pins.
#
# Read directly rather than via getTwasWeights(): that accessor returns a
# single-row collection, so routing the per-field accessors through it would
# dispatch straight back into them.
# @noRd
.twrSelectRowParts <- function(x, study, context, trait, method) {
    idx <- .tupleSelectRow(
        x,
        study,
        context,
        trait,
        method,
        cls = "TwasWeights"
    )
    .twrRowParts(x, idx)
}

#' @rdname resolveWeights
#' @export
setMethod(
    "resolveWeights",
    "TwasWeights",
    function(
        x,
        study = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        .twrRowResolveWeights(
            .twrSelectRowParts(x, study, context, trait, method)
        )
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
        method = NULL
    ) {
        getWeights(.twrSelectRowParts(x, study, context, trait, method))
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
        method = NULL
    ) {
        getCvResult(.twrSelectRowParts(x, study, context, trait, method))
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
        method = NULL
    ) {
        getFits(.twrSelectRowParts(x, study, context, trait, method))
    }
)

#' @rdname getStandardized
#' @export
setMethod(
    "getStandardized",
    "TwasWeights",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        isTRUE(
            getStandardized(.twrSelectRowParts(
                x,
                study,
                context,
                trait,
                method
            ))
        )
    }
)

#' @rdname getDataType
#' @export
setMethod(
    "getDataType",
    "TwasWeights",
    function(x, study = NULL, context = NULL, trait = NULL, method = NULL) {
        getDataType(.twrSelectRowParts(x, study, context, trait, method))
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
        method = NULL
    ) {
        .twrPartsVariantIds(
            .twrSelectRowParts(x, study, context, trait, method)
        )
    }
)

#' @rdname getStudy
#' @export
setMethod("getStudy", "TwasWeights", function(x) unique(as.character(x$study)))

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "TwasWeights", function(x) x@ldSketch)

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
        .ldSketchLabel(object@ldSketch)
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
    str_replace(x, "_weights$", str_c("_", target)) |>
        str_replace("Weights$", cap)
}

# Method name/impl/args lookup for TWAS weight methods. `fn` is the
# snake_case key used in weight method lists; `impl` is the camelCase
# function implemented by the package.
# @noRd
.twasMethodMap <- list(
    susie = list(
        fn = "susie_weights",
        impl = "susieWeights",
        # No fitting defaults: susieWeights extracts from a supplied fit.
        args = list()
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
    bayesR = list(
        fn = "bayes_r_weights",
        impl = "bayesRWeights",
        args = list()
    ),
    bayesL = list(
        fn = "bayes_l_weights",
        impl = "bLassoWeights",
        args = list()
    ),
    bayesA = list(
        fn = "bayes_a_weights",
        impl = "bayesAWeights",
        args = list()
    ),
    bayesB = list(
        fn = "bayes_b_weights",
        impl = "bayesBWeights",
        args = list()
    ),
    bayesC = list(
        fn = "bayes_c_weights",
        impl = "bayesCWeights",
        args = list()
    ),
    bayesN = list(
        fn = "bayes_n_weights",
        impl = "bayesNWeights",
        args = list()
    ),
    bLasso = list(
        fn = "b_lasso_weights",
        impl = "bLassoWeights",
        args = list()
    ),
    dprVb = list(
        fn = "dpr_vb_weights",
        impl = "dprVbWeights",
        args = list()
    ),
    dprGibbs = list(
        fn = "dpr_gibbs_weights",
        impl = "dprGibbsWeights",
        args = list()
    ),
    dprAdaptiveGibbs = list(
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
        # No fitting defaults: mvsusieWeights extracts from a supplied fit.
        args = list()
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

# Expand the `default` / `fastDefault` preset strings to their method vectors.
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
        if (methods == "fastDefault") {
            return(fastDefault)
        }
        if (methods == "default") {
            return(c(fastDefault, "bayesR", "bayesC"))
        }
    }
    methods
}

# --- TWAS method-token registry ---------------------------------------------

# A token's kwargs are a FLAT bag that .splitMethodArgs routes in two
# directions: names matching the weight function's own formals go to it, and
# everything else is forwarded to the fitting engine. So the accepted set is
# the union of both, and a name in neither reaches nothing -- the same
# situation as ctwas, hence `filtered = TRUE`.
# @noRd
.twasTokenAccepted <- function(token) {
    caps <- .twasMethodCapabilities[[token]]
    impls <- unique(compact(list(caps$individualImpl, caps$sumstatImpl)))
    implNames <- unlist(map(impls, .twasImplFormals))
    engine <- .twasTokenEngineNames(token)
    if (is.null(engine)) {
        return(NULL)
    }
    sort(unique(c(implNames, engine)))
}

# The weight function's own formals, minus the `methodArgs` slot that carries
# everything destined for the engine.
# @noRd
.twasImplFormals <- function(impl) {
    fn <- tryCatch(
        get(impl, envir = asNamespace("pecotmr")),
        error = function(cnd) NULL
    )
    if (!is.function(fn)) {
        return(character())
    }
    setdiff(names(formals(fn)), c("methodArgs", "..."))
}

# The fitting engine a token reaches, per input class. The two differ for
# every method with a summary-statistics counterpart: the individual path
# goes to the external package, the RSS path to pecotmr's own solver. Written
# out because it cannot be read off the capability table, which records the
# pecotmr wrapper rather than the engine behind it.
# @noRd
.twasTokenEngines <- function(token) {
    individual <- list(
        mrash = "susieR::mr.ash",
        enet = "glmnet::cv.glmnet",
        lasso = "glmnet::cv.glmnet",
        bayesA = "qgg::gbayes",
        bayesC = "qgg::gbayes",
        bayesN = "qgg::gbayes",
        bayesR = "qgg::gbayes",
        bayesB = "BGLR::BGLR",
        bayesL = "BGLR::BGLR",
        bLasso = "BGLR::BGLR",
        dprVb = "RcppDPR::fit_model",
        dprGibbs = "RcppDPR::fit_model",
        dprAdaptiveGibbs = "RcppDPR::fit_model",
        scad = "ncvreg::cv.ncvreg",
        mcp = "ncvreg::cv.ncvreg",
        l0learn = "L0Learn::L0Learn.cvfit",
        mrmash = "mr.mashr::mr.mash"
    )
    sumstat <- list(
        mrash = "susieR::mr.ash.rss",
        lasso = "pecotmr::lassosumRss",
        scad = "pecotmr::penalizedRss",
        mcp = "pecotmr::penalizedRss",
        l0learn = "pecotmr::penalizedRss",
        dprGibbs = "pecotmr::sdpr",
        mrmash = "mr.mashr::mr.mash.rss",
        prsCs = "pecotmr::prsCs"
    )
    list(individual = individual[[token]], sumstat = sumstat[[token]])
}

# The engine for one input class, or NULL when the token has no path there.
# @noRd
.twasTokenEngineFor <- function(token, inputKind) {
    e <- .twasTokenEngines(token)
    if (identical(inputKind, "QtlDataset")) e$individual else e$sumstat
}

# Every engine a token can reach, across classes.
# @noRd
.twasTokenEngineNames <- function(token) {
    callees <- unlist(compact(unname(.twasTokenEngines(token))))
    if (length(callees) == 0L) {
        return(NULL)
    }
    # filtered = FALSE: .splitMethodArgs does not DROP unrecognised names, it
    # forwards them to the engine -- so an engine with `...` really would
    # accept them and nothing can be rejected across the union. The per-class
    # check below narrows to one engine, where this is no longer true.
    .engineAcceptedNames(callees, FALSE)
}

# The per-token kwargs carried by a normalized methodList, keyed by canonical
# token. The list is keyed by weight-function name, which is the spelling the
# fitters want but not the one the engine map uses.
# @noRd
.twasMethodListArgs <- function(methodList) {
    if (length(methodList) == 0L || is.null(names(methodList))) {
        return(list())
    }
    set_names(
        as.list(methodList),
        str_remove(names(methodList), "(_weights|Weights)$")
    )
}

# --- per-token argument chains ----------------------------------------------
#
# A token's kwargs are consumed incrementally: each function along the way
# takes the names matching its own formals and forwards the rest. So the set
# of legal names is the union of formals along the WHOLE chain, not the
# formals of any one function -- validating against the final engine alone
# rejects legitimate wrapper arguments such as mrmash's
# `canonicalPriorMatrices`, which belongs to a middle hop.
#
# These chains were established by EXECUTION -- tracing which functions a real
# fit enters -- not by reading the sources, which hides forwarding behind
# list_modify()/exec() and imported-without-:: calls. test_twasWeights.R
# re-traces them and fails if one drifts.
# @noRd
.twasMethodChains <- function() {
    list(
        individual = list(
            lasso = c("lassoWeights", "glmnetWeights", "glmnet::cv.glmnet"),
            enet = c("enetWeights", "glmnetWeights", "glmnet::cv.glmnet"),
            scad = c("scadWeights", "ncvregWeights", "ncvreg::cv.ncvreg"),
            mcp = c("mcpWeights", "ncvregWeights", "ncvreg::cv.ncvreg"),
            l0learn = c("l0learnWeights", "L0Learn::L0Learn.cvfit"),
            mrash = c("mrashWeights", "susieR::mr.ash"),
            mrmash = c(
                "mrmashWeights",
                "mrmashWrapper",
                "mrmashPriorConfig",
                "buildMrmashPriorMatrices",
                "mr.mashr::mr.mash"
            ),
            bayesA = c("bayesAWeights", "bayesAlphabetWeights", "qgg::gbayes"),
            bayesC = c("bayesCWeights", "bayesAlphabetWeights", "qgg::gbayes"),
            bayesN = c("bayesNWeights", "bayesAlphabetWeights", "qgg::gbayes"),
            bayesR = c("bayesRWeights", "bayesAlphabetWeights", "qgg::gbayes"),
            bayesB = c("bayesBWeights", "bglrWeights", "BGLR::BGLR"),
            bayesL = c("bLassoWeights", "bglrWeights", "BGLR::BGLR"),
            bLasso = c("bLassoWeights", "bglrWeights", "BGLR::BGLR"),
            dprVb = c("dprVbWeights", "dprWeights", "RcppDPR::fit_model"),
            dprGibbs = c(
                "dprGibbsWeights",
                "dprWeights",
                "RcppDPR::fit_model"
            ),
            dprAdaptiveGibbs = c(
                "dprAdaptiveGibbsWeights",
                "dprWeights",
                "RcppDPR::fit_model"
            )
        ),
        sumstat = list(
            scad = c(
                "scadRssWeights",
                ".penalizedRssWeights",
                ".rssShrinkGridWeights",
                "penalizedRss"
            ),
            mcp = c(
                "mcpRssWeights",
                ".penalizedRssWeights",
                ".rssShrinkGridWeights",
                "penalizedRss"
            ),
            l0learn = c(
                "l0learnRssWeights",
                ".rssShrinkGridWeights",
                "penalizedRss"
            ),
            lasso = c(
                "lassosumRssWeights",
                ".rssShrinkGridWeights",
                "lassosumRss"
            ),
            prsCs = c("prsCsWeights", "prsCs"),
            dprGibbs = c("sdprWeights", "sdpr"),
            mrash = c("mrashRssWeights", "susieR::mr.ash.rss"),
            mrmash = c(
                "mrmashRssWeights",
                "buildMrmashPriorMatrices",
                "mr.mashr::mr.mash.rss"
            )
        )
    )
}

# The chain for one token and input class, or NULL when the token has no path
# for that class.
# @noRd
.twasMethodChainFor <- function(token, inputKind) {
    chains <- .twasMethodChains()
    side <- if (identical(inputKind, "QtlDataset")) "individual" else "sumstat"
    chains[[side]][[token]]
}

# Argument names that are data or plumbing rather than settable options, so
# naming one is not a user choice the check should endorse.
# @noRd
.twasChainDataArgs <- function() {
    c(
        "X",
        "Y",
        "y",
        "stat",
        "LD",
        "bhat",
        "shat",
        "R",
        "n",
        "S0",
        "w0",
        "sumstats",
        "Bhat",
        "Shat",
        "methodArgs",
        "dotArgs",
        "...",
        "solverInput",
        "LDs",
        "config",
        "bInit",
        "vInit",
        "K"
    )
}

# One hop's formals. A hop is either an external "pkg::fn" or a pecotmr
# function; NA marks one that could not be read, which makes the whole chain
# uncheckable rather than silently narrower.
# @noRd
.twasChainHopFormals <- function(nm) {
    fn <- if (str_detect(nm, fixed("::"))) {
        .engineCallee(nm)
    } else {
        tryCatch(
            get(nm, envir = asNamespace("pecotmr")),
            error = function(cnd) NULL
        )
    }
    if (!is.function(fn)) {
        return(NA_character_)
    }
    names(formals(fn))
}

# Every name legal for a token on one input class: the union of formals along
# its chain. NULL when a hop takes `...`, since anything then reaches it.
# @noRd
.twasChainAccepted <- function(token, inputKind) {
    chain <- .twasMethodChainFor(token, inputKind)
    if (is.null(chain)) {
        return(NULL)
    }
    all <- map(chain, .twasChainHopFormals)
    if (any(map_lgl(all, .twasChainHopUnknown))) {
        return(NULL)
    }
    setdiff(unique(list_c(all)), .twasChainDataArgs())
}

# A hop contributes nothing checkable when it could not be read, or when it
# takes `...` and so accepts any name.
# @noRd
.twasChainHopUnknown <- function(fm) {
    (length(fm) == 1L && is.na(fm[[1L]])) || is_in("...", fm)
}


# The argument constructor for each weight-method token. Several tokens share
# an engine, so this maps token -> constructor rather than token -> package.
# Tokens absent from this list take no options at all; .twasTokenNoArgsReason
# says why, so the error can be specific.
# @noRd
.twasMethodCtors <- function() {
    list(
        mrash = mrashConfig,
        enet = glmnetConfig,
        lasso = glmnetConfig,
        bayesA = qggConfig,
        bayesC = qggConfig,
        bayesN = qggConfig,
        bayesR = qggConfig,
        bayesB = bglrConfig,
        bayesL = bglrConfig,
        bLasso = bglrConfig,
        dprVb = dprConfig,
        dprGibbs = dprConfig,
        dprAdaptiveGibbs = dprConfig,
        scad = ncvregConfig,
        mcp = ncvregConfig,
        l0learn = l0learnConfig,
        mrmash = mrmashConfig,
        prsCs = prsCsConfig
    )
}

# Entries may be keyed by the short token (`lasso`) or by the full weight
# function name (`lasso_weights`), the same two spellings .twasMethodLookup
# accepts. Canonicalize to the short form so the registry has one key per
# method and errors name the canonical spelling.
# @noRd
.twasCanonicalEntryNames <- function(entries) {
    if (length(entries) == 0L || is.null(names(entries))) {
        return(entries)
    }
    fnToShort <- set_names(
        names(.twasMethodMap),
        map_chr(.twasMethodMap, "fn")
    )
    set_names(
        entries,
        map_chr(names(entries), .twasCanonicalEntryName, fnToShort = fnToShort)
    )
}

# One entry name in canonical short form. .twasMethodMap covers the tokens
# with an individual-level fitter; a sumstat-only token such as prsCs is not
# in it, so a `<token>_weights` spelling is also recognised by stripping the
# suffix when what remains is a known method.
# @noRd
.twasCanonicalEntryName <- function(nm, fnToShort) {
    short <- .twasCanonicalShortName(nm, fnToShort)
    if (is_in(short, names(.twasMethodCtors()))) {
        return(short)
    }
    stripped <- str_remove(nm, "(_weights|Weights)$")
    if (is_in(stripped, names(.twasMethodCtors()))) {
        return(stripped)
    }
    short
}

# Why a token takes no options, so the error can say something useful rather
# than just "unknown method".
# @noRd
.twasTokenNoArgsReason <- function(token) {
    fromFit <- c("susie", "susieAsh", "susieInf", "mvsusie", "fsusie")
    if (is_in(token, fromFit)) {
        return(glue(
            "'{token}' extracts weights from a fit supplied via ",
            "`fineMappingResult`, so there is nothing to configure here. ",
            "Configure the fit in fineMappingPipeline() instead."
        ))
    }
    NULL
}

#' @title Per-Method Arguments For twasWeightsPipeline
#' @description Options for each TWAS weight method, keyed by method token.
#'   Naming a method here also selects it, so this is what
#'   \code{twasWeightsPipeline(methods = )} takes when you want to configure
#'   the fit; a plain character vector remains the shorthand for running
#'   methods with their defaults.
#'
#'   Each entry may be a plain list or the engine's own constructor --- a
#'   plain list is spliced into that constructor, so it gets the same
#'   checking either way. Several tokens share an engine and therefore share
#'   a constructor: \code{lasso} and \code{enet} both take
#'   \code{\link{glmnetConfig}}, \code{bayesA}/\code{bayesC}/\code{bayesN}/
#'   \code{bayesR} take \code{\link{qggConfig}}, and so on.
#'
#'   Methods that extract weights from a fit you supply --- \code{susie},
#'   \code{susieAsh}, \code{susieInf}, \code{mvsusie}, \code{fsusie} ---
#'   have nothing to configure and are rejected here; select them with the
#'   character form and configure the fit in
#'   \code{\link{fineMappingPipeline}}.
#' @param ... Named entries, one per method token.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' twasWeightsMethodsConfig(lasso = list(nfold = 10), mrmash = mrmashConfig())
#' @export
twasWeightsMethodsConfig <- function(...) {
    entries <- .twasCanonicalEntryNames(list(...))
    ctors <- .twasMethodCtors()
    optionless <- setdiff(names(entries), names(ctors))
    reasons <- compact(map(optionless, .twasTokenNoArgsReason))
    if (length(reasons) > 0L) {
        abort(glue(
            "twasWeightsMethodsConfig: ",
            "{str_flatten(map_chr(reasons, as.character), ' ')}"
        ))
    }
    .newNestedConfig(
        entries,
        set_names(
            map(names(ctors), .twasTokenCtor),
            names(ctors)
        ),
        "twasWeightsMethodsConfig",
        engines = set_names(
            map(names(ctors), .twasTokenEngineLabels),
            names(ctors)
        )
    )
}

# A token's summary-statistic constructor, where that path runs a different
# engine from the individual-level one. The four pecotmr solvers are real
# engines with their own formals, so each has its own constructor; a token
# whose two paths share an engine is absent here.
# @noRd
.twasMethodRssCtors <- function() {
    list(
        lasso = lassosumConfig,
        scad = penalizedRssConfig,
        mcp = penalizedRssConfig,
        l0learn = penalizedRssConfig,
        dprGibbs = sdprConfig
    )
}

# Every engine label a token's entry may declare. A method's two input paths
# can run different engines -- lasso fits with glmnet on individual data and
# with pecotmr's lassosum solver on summary statistics -- and either
# engine's constructor is a legitimate way to configure that method.
# @noRd
.twasTokenEngineLabels <- function(token) {
    ctors <- compact(list(
        .twasMethodCtors()[[token]],
        .twasMethodRssCtors()[[token]]
    ))
    unique(map_chr(ctors, .twasCtorEngineLabel, token = token))
}

# The engine label a constructor records, falling back to the token when the
# constructor cannot be built (its engine's package may be absent).
# @noRd
.twasCtorEngineLabel <- function(ctor, token) {
    rec <- tryCatch(ctor(), error = function(cnd) NULL)
    if (is.null(rec)) {
        return(token)
    }
    metadata(rec)$engine %||% token
}

# A token's constructor for the aggregator: validates the flat kwargs against
# the union of the weight function's formals and its engine's, since
# .splitMethodArgs routes them to both.
# @noRd
.twasTokenCtor <- function(token) {
    force(token)
    # Labelled with the ENGINE, not the token, so that passing the engine's
    # own constructor -- glmnetConfig() for either lasso or enet -- is
    # recognised as a match rather than reported as a mismatch.
    engine <- .twasTokenEngineLabel(token)
    function(...) {
        accepted <- .twasTokenAccepted(token)
        .newMethodConfig(
            NULL,
            defaults = list(),
            extra = list(...),
            label = glue("twasWeightsMethodsConfig: method '{token}'"),
            engine = engine,
            accepted = accepted,
            check = !is.null(accepted)
        )
    }
}

# The engine label a token's entries carry, taken from that token's engine
# constructor so the two agree.
# @noRd
.twasTokenEngineLabel <- function(token) {
    .twasCtorEngineLabel(.twasMethodCtors()[[token]], token)
}

# Map short method names and presets to weightMethods lists.
# @param methods A character vector of short method names, or a preset string
#   ("default" or "fastDefault").
# @return A named list suitable for the weightMethods parameter.
# @importFrom purrr map_chr set_names
# @noRd
.twasMethodLookup <- function(methods) {
    expanded <- .twasExpandPresets(methods)
    # Accept full function names too by mapping them back to short names.
    fnToShort <- set_names(names(.twasMethodMap), map_chr(.twasMethodMap, "fn"))
    shortNames <- map_chr(
        expanded,
        .twasCanonicalShortName,
        fnToShort = fnToShort
    )
    unknown <- setdiff(shortNames, names(.twasMethodMap))
    if (length(unknown) > 0) {
        msg <- glue(
            "unknown method token(s): {str_flatten(unknown, ', ')}. ",
            "Known tokens: {str_flatten(names(.twasMethodMap), ', ')}."
        )
        abort(msg)
    }
    # Track the impl name as an attr so dispatchers resolve snake_case -> impl.
    entries <- map(shortNames, .twasMethodArgsWithImpl)
    set_names(entries, map_chr(shortNames, .twasMethodFn))
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
    df <- mutate(df, Sample = as.character(.data$Sample))
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
    full <- matrix(
        0,
        nrow = nColsX,
        ncol = nColsY,
        dimnames = list(XColnames, YColnames)
    )
    full[validColumns, ] <- weightsMatrix
    full
}

# The canonical fine-mapping token behind a weight-method name, accepting the
# token itself (`susieInf`), its camelCase weight function (`susieInfWeights`)
# and its snake_case registry key (`susie_inf_weights`). Stripping the suffix
# is not enough on its own: that turns `susie_inf_weights` into `susie_inf`,
# which is not the registry's token.
# @noRd
.twasFmTokenFor <- function(method) {
    methodKeys <- map_chr(.twasFineMappingMethodAdapters, "methodKey")
    hit <- names(methodKeys)[methodKeys == method]
    if (length(hit) == 1L) {
        return(hit)
    }
    bare <- str_remove(method, "(_weights|Weights)$")
    if (is_in(bare, names(.twasFineMappingMethodAdapters))) {
        return(bare)
    }
    NA_character_
}

# The `weightMethods` name under which `token` was supplied, or NA when it was
# not requested. Every lookup keyed by method name goes through this, so one
# spelling cannot resolve where another silently misses.
# @noRd
.twasMethodNameFor <- function(weightMethods, token) {
    adapter <- .twasFineMappingMethodAdapters[[token]]
    spellings <- c(adapter$methodKey, adapter$weightFn, token)
    hit <- spellings[is_in(spellings, names(weightMethods))]
    if (length(hit) == 0L) {
        return(NA_character_)
    }
    hit[[1]]
}

# The arguments supplied for `token`, under whichever of its spellings the
# caller used as the `weightMethods` name.
# @noRd
.twasMethodArgsFor <- function(weightMethods, token) {
    nm <- .twasMethodNameFor(weightMethods, token)
    if (is.na(nm)) {
        return(NULL)
    }
    weightMethods[[nm]]
}

# The per-fold fits supplied for `token`, under any of its spellings. The gate
# and `.twasFoldFit()` must agree on this, or the gate can refuse fits the
# fold lookup would have found.
# @noRd
.twasFoldFitsFor <- function(fittedModelsCv, token) {
    if (is.null(fittedModelsCv)) {
        return(NULL)
    }
    adapter <- .twasFineMappingMethodAdapters[[token]]
    spellings <- c(token, adapter$methodKey, adapter$weightFn)
    hit <- spellings[is_in(spellings, names(fittedModelsCv))]
    if (length(hit) == 0L) {
        return(NULL)
    }
    fittedModelsCv[[hit[[1]]]]
}

# The SuSiE-family methods requested in `weightMethods`, as canonical tokens.
# @noRd
.twasSusieTokensRequested <- function(weightMethods) {
    if (is.null(weightMethods) || length(weightMethods) == 0L) {
        return(character(0))
    }
    keys <- if (is.character(weightMethods)) {
        weightMethods
    } else {
        names(weightMethods)
    }
    tokens <- map_chr(keys, .twasFmTokenFor)
    unique(tokens[!is.na(tokens)])
}

# Fail before any fitting work when a SuSiE-family method is requested without
# the fit it needs. These wrappers extract weights from an existing fit and
# never fine-map, so the run is already doomed; checking here reports the
# missing fit by name instead of surfacing it from inside the per-fold map.
# `available(token)` answers "is this token's fit present?" for the caller's
# own supply channel.
# @noRd
.twasRequireSusieFits <- function(weightMethods, available, fnLabel, how) {
    tokens <- .twasSusieTokensRequested(weightMethods)
    if (length(tokens) == 0L) {
        return(invisible(NULL))
    }
    missing <- tokens[!map_lgl(tokens, available)]
    if (length(missing) == 0L) {
        return(invisible(NULL))
    }
    missingStr <- str_flatten(missing, ", ")
    msg <- glue(
        "{fnLabel}: method(s) {missingStr} extract weights from an existing ",
        "fine-mapping fit and never run fine-mapping themselves, but no fit ",
        "was supplied for them. {how}"
    )
    abort(msg)
}

# Resolve the incoming susie / susieInf fits for the weight run: prefer a fit
# carried on the method args, else one from `fittedModels`, tagging each with
# its fine-mapping class. Returns list(susieFit, susieInfFit, hasSusie,
# hasSusieInf, susieName, susieInfName) -- the names being the spellings the
# caller actually used, so the write-back lands where the fitter will look.
# @noRd
.twasResolveSusieFits <- function(weightMethods, fittedModels) {
    susieName <- .twasMethodNameFor(weightMethods, "susie")
    susieInfName <- .twasMethodNameFor(weightMethods, "susieInf")
    hasSusie <- !is.na(susieName)
    hasSusieInf <- !is.na(susieInfName)
    argSusie <- if (hasSusie) {
        weightMethods[[susieName]][["susieFit"]]
    } else {
        NULL
    }
    argSusieInf <- if (hasSusieInf) {
        weightMethods[[susieInfName]][["susieInfFit"]]
    } else {
        NULL
    }
    # A fit the caller passed in the method arguments wins; otherwise the
    # pipeline's already-fitted model, if it produced one.
    susieRaw <- argSusie %||% fittedModels[["susie"]]
    susieInfRaw <- argSusieInf %||% fittedModels[["susieInf"]]
    susieFit <- if (is.null(susieRaw)) {
        NULL
    } else {
        .setFinemappingFitClass(susieRaw, "susie")
    }
    susieInfFit <- if (is.null(susieInfRaw)) {
        NULL
    } else {
        .setFinemappingFitClass(susieInfRaw, "susieInf")
    }
    list(
        susieFit = susieFit,
        susieInfFit = susieInfFit,
        hasSusie = hasSusie,
        hasSusieInf = hasSusieInf,
        susieName = susieName,
        susieInfName = susieInfName
    )
}

# Write the resolved susie / susieInf fits back onto the method args. Deriving
# susie's arguments from a susieInf fit belongs to fine-mapping, not here: it
# prepares a susie *fit* (model_init, unmappable_effects), which
# `susieWeights()` neither accepts nor runs.
# @noRd
.twasWriteBackSusieFits <- function(weightMethods, r) {
    withInf <- if (!is.null(r$susieInfFit) && r$hasSusieInf) {
        .twasSetMethodArg(
            weightMethods,
            r$susieInfName,
            "susieInfFit",
            r$susieInfFit
        )
    } else {
        weightMethods
    }
    if (is.null(r$susieFit) || !r$hasSusie) {
        return(withInf)
    }
    .twasSetMethodArg(withInf, r$susieName, "susieFit", r$susieFit)
}

# One method's argument list with `arg` set, leaving every other method and
# every other argument of that method untouched.
# @noRd
.twasSetMethodArg <- function(weightMethods, method, arg, value) {
    updated <- list_assign(
        weightMethods[[method]],
        !!!set_names(list(value), arg)
    )
    list_assign(weightMethods, !!!set_names(list(updated), method))
}

# Resolve the supplied susie / susieInf fits onto the method args. This never
# fine-maps: the SuSiE-family weight methods extract from a fit that
# fineMappingPipeline() produced, and a missing one is an error, not a cue to
# fit here.
# @noRd
.prepareSusieWeightMethods <- function(weightMethods, fittedModels = NULL) {
    if (is.null(fittedModels)) {
        fittedModels <- list()
    }
    r <- .twasResolveSusieFits(weightMethods, fittedModels)
    .twasWriteBackSusieFits(weightMethods, r)
}

# Per-fold TWAS weight fit for the CV engine. `ctx` carries weightMethods,
# multivariateWeightMethods, cvArgs, fitRetention, verbose. Weights are keyed by
# the canonical method key; captured fits keep the full method name.
# @noRd
.weightFitFold <- function(Xtr, Ytr, j, ctx) {
    foldWeightMethods <- .prepareSusieWeightMethods(ctx$weightMethods)
    methods <- names(foldWeightMethods)
    results <- map(
        methods,
        .twasFoldMethodWeightsAt,
        foldWeightMethods = foldWeightMethods,
        Xtr = Xtr,
        Ytr = Ytr,
        j = j,
        ctx = ctx
    )
    # Weights are keyed by the canonical method key, captured fits by the
    # full method name; a later method wins a shared key, as the keyed
    # assignment did.
    weights <- set_names(map(results, "W"), map_chr(results, "mk"))
    list(
        weights = weights[!duplicated(names(weights), fromLast = TRUE)],
        fits = set_names(map(results, "fit"), methods)
    )
}

# @noRd
.twasFoldMethodWeightsAt <- function(
    method,
    foldWeightMethods,
    Xtr,
    Ytr,
    j,
    ctx
) {
    .twasFoldMethodWeights(
        method,
        foldWeightMethods[[method]],
        Xtr,
        Ytr,
        j,
        ctx
    )
}

# Per-fold priors bound to a multivariate fitter's camelCase args (mr.mash
# data-driven matrices / mvsusie reweighted mixture prior for fold `j`).
# @noRd
.twasFoldPriors <- function(
    args,
    method,
    j,
    dataDrivenPriorMatricesCv,
    reweightedMixturePriorCv
) {
    list_assign(
        args,
        !!!compact(list(
            dataDrivenPriorMatrices = if (
                !is.null(dataDrivenPriorMatricesCv) &&
                    is_in(method, c("mrmash_weights", "mrmashWeights"))
            ) {
                dataDrivenPriorMatricesCv[[j]]
            },
            prior_variance = if (
                !is.null(reweightedMixturePriorCv) &&
                    is_in(method, c("mvsusie_weights", "mvsusieWeights"))
            ) {
                reweightedMixturePriorCv[[j]]
            }
        ))
    )
}

# Inject the fold's own fine-mapping fit for a SuSiE-family method. Those
# weight wrappers extract from a supplied fit and never fine-map, so a fold
# can only be scored if fineMappingPipeline's CV retained that fold's fit
# (`fittedModelsCv`, keyed method -> fold_<j>).
# @noRd
.twasFoldFit <- function(args, method, j, fittedModelsCv) {
    if (is.null(fittedModelsCv)) {
        return(args)
    }
    mk <- .twasFmTokenFor(method)
    if (is.na(mk)) {
        return(args)
    }
    adapter <- .twasFineMappingMethodAdapters[[mk]]
    perFold <- .twasFoldFitsFor(fittedModelsCv, mk)
    if (is.null(perFold)) {
        return(args)
    }
    key <- str_c("fold_", j)
    fit <- if (is_in(key, names(perFold))) {
        perFold[[key]]
    } else if (length(perFold) >= j) {
        perFold[[j]]
    } else {
        NULL
    }
    if (is.null(fit)) {
        return(args)
    }
    list_assign(args, !!!set_names(list(fit), adapter$fitArg))
}

# One fold's multivariate weight fit; returns list(W, fit).
# @noRd
.twasFoldMultivariate <- function(method, fnName, args, Xtr, Ytr, j, ctx) {
    withPriors <- .twasFoldPriors(
        args,
        method,
        j,
        ctx$dataDrivenPriorMatricesCv,
        ctx$reweightedMixturePriorCv
    )
    retaining <- !identical(ctx$fitRetention, "none") &&
        is_in("fitRetention", names(formals(fnName)))
    fitArgs <- list_assign(
        withPriors,
        !!!compact(list(
            fitRetention = if (retaining) ctx$fitRetention
        ))
    )
    callArgs <- .twasWeightCallArgs(fnName, list(X = Xtr, Y = Ytr), fitArgs)
    W <- if (ctx$verbose < 2) {
        .quietEval(exec(fnName, !!!callArgs))
    } else {
        exec(fnName, !!!callArgs)
    }
    capturedFit <- attr(W, "fit")
    bare <- `attr<-`(W, "fit", NULL)
    list(W = `rownames<-`(bare, colnames(Xtr)), fit = capturedFit)
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
    W <- `rownames<-`(exec(cbind, !!!Wcols), colnames(Xtr))
    list(W = W, fit = NULL)
}

# Fit one method for one CV fold: dispatch to the multivariate or univariate
# path. Returns list(mk, W, fit) keyed by the canonical method key `mk`.
# @noRd
.twasFoldMethodWeights <- function(method, args, Xtr, Ytr, j, ctx) {
    fnName <- .resolveMethodFunction(method, args)
    args <- .twasFoldFit(args, method, j, ctx$fittedModelsCv)
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
#' @param fitRetention How much of each per-method fit is kept on the result:
#'   \code{"none"} (default) keeps none, \code{"slim"} keeps the trimmed
#'   payload downstream needs, \code{"full"} keeps the whole fitted object.
#'   Only the mr.mash methods distinguish \code{"slim"} from \code{"full"};
#'   for every other method the two behave alike.
#' @param seed Integer or \code{NULL}. When supplied, seeds both the
#'   main-process RNG (fold partitioning, variant sub-sampling) and the
#'   parallel fold-fitting RNG via the \code{BiocParallel} \code{RNGseed}, so
#'   results are reproducible even under multi-threading. The main-process
#'   seed is scoped to the call, so the session RNG is left as it was found.
#'   \code{NULL} (default) does not seed at all and uses the historical
#'   parallel default.
#' @param fittedModelsCv Optional per-fold fine-mapping fits, as
#'   \code{method -> fold_<j> -> fit}, from a \code{fineMappingPipeline()}
#'   run with \code{cvFolds > 1}. SuSiE-family weight wrappers extract from a
#'   supplied fit and never fine-map, so this is what makes cross-validating
#'   them possible; it must have been produced on \code{samplePartitions},
#'   since a fit trained on different folds would leak held-out samples.
#' @param dataDrivenPriorMatricesCv Optional list, one element per fold, of
#'   data-driven prior matrices for the mr.mash learner.
#' @param reweightedMixturePriorCv Optional list, one element per fold, of
#'   reweighted mixture priors for the mvSuSiE learner.
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
#' @importFrom BiocParallel bplapply multicoreWorkers MulticoreParam
#' @importFrom quadprog solve.QP
#' @examples
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:80]
#' Y <- multiTraitData$Y
#' # A cross-validated method is refit on each fold's training rows. The
#' # SuSiE family never fine-maps, so it needs each fold's own fit, passed
#' # as `fittedModelsCv` from a fineMappingPipeline() run with cvFolds > 1.
#' twasWeightsCv(X, Y[, 1, drop = FALSE], fold = 3,
#'   weightMethods = list(lasso_weights = list()))
#' @importFrom checkmate assertDataFrame assertNumber assertInt
#' @importFrom checkmate assertFlag assertCount
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
    fitRetention = c("none", "slim", "full"),
    seed = NULL,
    dataDrivenPriorMatricesCv = NULL,
    reweightedMixturePriorCv = NULL,
    fittedModelsCv = NULL
) {
    # X / Y / fold are asserted downstream in .cvPrepareData; these are the
    # arguments nothing else checks.
    assertDataFrame(samplePartitions, null.ok = TRUE)
    # NOT assertCount: `Inf` is the "no cap" sentinel (jointEngine passes it
    # when cfg$crossValidation$maxVariants is unset), and
    # .cvSubsampleVariants relies on `ncol(X) <= maxNumVariants` being FALSE
    # for it.
    assertNumber(maxNumVariants, lower = 1, null.ok = TRUE)
    assertInt(numThreads)
    assertCount(verbose)
    fitRetention <- arg_match(fitRetention)
    assertInt(seed, null.ok = TRUE)
    .twasWeightsCvImpl(
        X = X,
        Y = Y,
        fold = fold,
        samplePartitions = samplePartitions,
        weightMethods = weightMethods,
        maxNumVariants = maxNumVariants,
        variantsToKeep = variantsToKeep,
        numThreads = numThreads,
        verbose = verbose,
        fitRetention = fitRetention,
        seed = seed,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        reweightedMixturePriorCv = reweightedMixturePriorCv,
        fittedModelsCv = fittedModelsCv
    )
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

# Refuse per-fold fits that were not produced on the folds being scored. A
# fit trained on a different split has seen some of this split's held-out
# samples, so its out-of-fold predictions are contaminated and the CV metrics
# come out optimistic -- silently. Fits carry the producer's partition
# fingerprint (see .cvPartitionKey); anything unstamped is refused too, since
# it cannot be shown to match.
# @noRd
.twasCheckFoldFitPartition <- function(fittedModelsCv, samplePartitions) {
    if (is.null(fittedModelsCv) || length(fittedModelsCv) == 0L) {
        return(invisible(NULL))
    }
    if (is.null(samplePartitions)) {
        msg <- glue(
            "twasWeightsCv: `fittedModelsCv` needs the fold partition those ",
            "fits were trained on. Pass the fine-mapping CV's ",
            "`samplePartition` as `samplePartitions`; a freshly drawn ",
            "partition would score each fold with a fit that saw its ",
            "held-out samples."
        )
        abort(msg)
    }
    want <- .cvPartitionKey(samplePartitions)
    for (m in names(fittedModelsCv)) {
        got <- attr(fittedModelsCv[[m]], "partitionKey")
        if (is.null(got)) {
            msg <- glue(
                "twasWeightsCv: the per-fold fits for '{m}' carry no ",
                "partition fingerprint, so they cannot be shown to match ",
                "`samplePartitions`. Take them from a fineMappingPipeline() ",
                "run with cvFolds > 1."
            )
            abort(msg)
        }
        if (!identical(got, want)) {
            msg <- glue(
                "twasWeightsCv: the per-fold fits for '{m}' were trained on ",
                "a different fold partition than the one being scored. Use ",
                "the fine-mapping CV's own `samplePartition`."
            )
            abort(msg)
        }
    }
    invisible(NULL)
}

# twasWeightsCv worker. With no weightMethods the caller only wants the fold
# partition.
# @noRd
.twasWeightsCvImpl <- function(
    X,
    Y,
    fold,
    samplePartitions,
    weightMethods,
    maxNumVariants,
    variantsToKeep,
    numThreads,
    verbose,
    fitRetention,
    seed,
    dataDrivenPriorMatricesCv,
    reweightedMixturePriorCv,
    fittedModelsCv
) {
    .twasCheckFoldFitPartition(fittedModelsCv, samplePartitions)
    .twasRequireSusieFits(
        weightMethods,
        available = function(tk) {
            !is.null(.twasFoldFitsFor(fittedModelsCv, tk))
        },
        fnLabel = "twasWeightsCv",
        how = str_c(
            "Cross-validation refits on each fold, so it needs that fold's ",
            "own fit: pass `fittedModelsCv = list(<token> = <fold fits>)` ",
            "from a fineMappingPipeline() run with cvFolds > 1."
        )
    )
    weightMethods <- if (is.character(weightMethods)) {
        .twasMethodLookup(weightMethods)
    } else {
        weightMethods
    }
    if (is.null(seed) && !exists(".Random.seed") && verbose >= 1) {
        inform(str_c(
            "! No seed set. Pass `seed=` or call ",
            "set.seed() for reproducibility."
        ))
    }
    if (is.null(weightMethods)) {
        res <- .crossValidateWeights(
            X,
            Y,
            fold = fold,
            samplePartitions = samplePartitions,
            fitFold = .cvNoopFitFold,
            numThreads = numThreads,
            maxNumVariants = maxNumVariants,
            variantsToKeep = variantsToKeep,
            fitRetention = fitRetention,
            verbose = verbose,
            seed = seed
        )
        return(list(samplePartition = res$samplePartition))
    }
    cvFitCtx <- list(
        weightMethods = weightMethods,
        multivariateWeightMethods = .twasCvMultivariateMethods,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        reweightedMixturePriorCv = reweightedMixturePriorCv,
        fittedModelsCv = fittedModelsCv,
        fitRetention = fitRetention,
        verbose = verbose
    )
    .crossValidateWeights(
        X,
        Y,
        fold = fold,
        samplePartitions = samplePartitions,
        fitFold = .weightFitFold,
        fitFoldCtx = cvFitCtx,
        numThreads = numThreads,
        maxNumVariants = maxNumVariants,
        variantsToKeep = variantsToKeep,
        fitRetention = fitRetention,
        verbose = verbose,
        seed = seed
    )
}

# Fit one TWAS weight method by name against the filtered design matrix,
# embedding the fitted weights back into the full variant space. `ctx` carries
# the shared fit state (X, Y, Xfiltered, validColumns, fitRetention,
# verbose).
# @noRd
.computeMethodWeights <- function(methodName, weightMethods, ctx) {
    shortName <- str_remove(methodName, "_weights$")
    if (ctx$verbose >= 1) {
        inform(glue("  Fitting {shortName} ..."))
        tic()
    }
    userArgs <- weightMethods[[methodName]]
    fnName <- .resolveMethodFunction(methodName, userArgs)
    fitArgs <- .twasApplyRetainFit(
        userArgs,
        fnName,
        ctx$fitRetention
    )
    fit <- .twasFitWeightsMatrix(fnName, fitArgs, ctx, methodName)
    embedded <- .embedWeights(
        fit$weights,
        ctx$validColumns,
        ncol(ctx$X),
        ncol(ctx$Y),
        colnames(ctx$X),
        colnames(ctx$Y)
    )
    result <- if (is.null(fit$methodFit)) {
        embedded
    } else {
        `attr<-`(embedded, "fit", fit$methodFit)
    }
    if (ctx$verbose >= 1) {
        elapsed <- toc(quiet = TRUE)
        secs <- sprintf("%.1f", elapsed$toc - elapsed$tic)
        inform(glue("  Fitting {shortName} done in {secs}s"))
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

# Pass the retention level down to the target weight function, when it takes
# one and the caller is keeping fits at all. A target that does not take
# `fitRetention` keeps nothing, which is what "no such argument" means.
# @noRd
.twasApplyRetainFit <- function(args, fnName, fitRetention) {
    if (identical(fitRetention, "none")) {
        return(args)
    }
    if (!is_in("fitRetention", names(formals(fnName)))) {
        return(args)
    }
    list_assign(args, fitRetention = args$fitRetention %||% fitRetention)
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
    call <- .twasWeightCallArgs(
        fnName,
        list(X = ctx$Xfiltered, Y = ctx$Y),
        args
    )
    fitted <- if (ctx$verbose < 2) {
        .quietEval(exec(fnName, !!!call))
    } else {
        exec(fnName, !!!call)
    }
    methodFit <- if (!identical(ctx$fitRetention, "none")) {
        attr(fitted, "fit")
    } else {
        NULL
    }
    weightsMatrix <- if (nrow(fitted) == length(ctx$validColumns)) {
        fitted
    } else {
        fitted[names(ctx$validColumns), , drop = FALSE]
    }
    list(weights = weightsMatrix, methodFit = methodFit)
}

# Univariate fit: apply the method to each column of Y, filling a zero-init
# weights matrix.
# @noRd
.twasFitUnivariate <- function(fnName, args, ctx) {
    columns <- map(
        seq_len(ncol(ctx$Y)),
        .twasUnivariateColumn,
        fnName = fnName,
        args = args,
        ctx = ctx
    )
    fits <- compact(map(columns, "fit"))
    list(
        weights = matrix(
            unname(list_c(map(columns, "weights"))),
            nrow = ncol(ctx$Xfiltered),
            ncol = ncol(ctx$Y)
        ),
        # The first outcome that carries one defines the method fit, which is
        # what "set it only while still NULL" produced.
        methodFit = if (length(fits) == 0L) NULL else fits[[1L]]
    )
}

# One outcome column's weights, plus the fit the method attached to them. A
# method that answers with a matrix is reporting every outcome at once, so
# this column's own slice is taken.
# @noRd
.twasUnivariateColumn <- function(k, fnName, args, ctx) {
    call <- .twasWeightCallArgs(
        fnName,
        list(X = ctx$Xfiltered, y = ctx$Y[, k]),
        args
    )
    w <- if (ctx$verbose < 2) {
        .quietEval(exec(fnName, !!!call))
    } else {
        exec(fnName, !!!call)
    }
    list(
        weights = if (is.matrix(w)) w[, k] else w,
        fit = if (!identical(ctx$fitRetention, "none")) {
            attr(w, "fit")
        } else {
            NULL
        }
    )
}

# Assemble the (study, context, trait, method, entry) row vectors for the
# TwasWeights collection from the fitted `weightsList`. `ctx` carries the shared
# identity + flags (study, context, trait, Y, fitRetention, standardized,
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

# One TwasWeightsRow for a (variantIds, weights) pair with the shared flags.
# @noRd
.twasEntry <- function(variantIds, weights, fits, ctx) {
    twasWeightsRow(
        variantIds = variantIds,
        weights = weights,
        fits = fits,
        cvResult = NULL,
        standardized = ctx$standardized,
        dataType = ctx$dataType
    )
}

# Row-records for one fitted method. When trait/context were supplied per-row
# (length == ncol(Y)) emit one (method, outcome) row per Y column; otherwise a
# single row carrying the (possibly multi-column) weights matrix as-is.
# @noRd
.twasMethodRows <- function(m, wMat, variantIds, ctx) {
    fitVal <- attr(wMat, "fit")
    wMat <- `attr<-`(wMat, "fit", NULL)
    fits <- if (!identical(ctx$fitRetention, "none")) fitVal else NULL
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
#' @param fittedModels Named list of fitted SuSiE-family models, keyed by
#'   token (\code{susie}, \code{susieInf}, \code{mvsusie}, \code{fsusie}).
#'   Required whenever a SuSiE-family weight method is requested: those
#'   methods extract weights from an existing fit and never fine-map, so a
#'   missing fit is an error. Run \code{\link{fineMappingPipeline}} to
#'   produce the fits.
#' @param fitRetention How much of each per-method fit is kept on the result:
#'   \code{"none"} (default) keeps none, \code{"slim"} keeps the trimmed
#'   payload downstream needs, \code{"full"} keeps the whole fitted object.
#'   Only the mr.mash methods distinguish \code{"slim"} from \code{"full"};
#'   for every other method the two behave alike.
#' @param verbose Integer controlling verbosity level: 0 = suppress all
#'   messages, 1 = suppress external package messages (default), 2 = show all
#'   messages including those from external packages.
#' @param study Character. Study identity label recorded on the resulting
#'   weights.
#' @param context Character. Context identity label recorded on the resulting
#'   weights.
#' @param trait Character. Trait identity label recorded on the resulting
#'   weights.
#' @param standardized Logical. Whether the supplied \code{X} / \code{Y} are
#'   already standardized. Default \code{FALSE}.
#' @param dataType Character or \code{NULL}. Data-type label recorded on the
#'   weights (e.g. \code{"individual"}).
#' @param ldSketch A genotype panel (see \code{\link{readGenotypes}}) to
#'   record on the weights as their LD sketch, or \code{NULL}.
#' @param seed Integer or \code{NULL}. When supplied, seeds the main-process
#'   RNG and the parallel method-fitting RNG via the \code{BiocParallel}
#'   \code{RNGseed}, for reproducibility under multi-threading. The
#'   main-process seed is scoped to the call, so the session RNG is left as it
#'   was found. \code{NULL} (default) does not seed at all.
#' @return A list where each element is named after a method and contains the
#'   weight matrix produced by that method.
#'
#' @examples
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:80]
#' Y <- multiTraitData$Y
#' # SuSiE-family methods extract weights from an existing fit and never
#' # fine-map themselves, so the fit is supplied via `fittedModels`.
#' fit <- susieR::susie(X, Y[, 1], L = 5)
#' learnTwasWeights(X, Y[, 1, drop = FALSE],
#'   weightMethods = list(susie_weights = list()),
#'   fittedModels = list(susie = fit))
#' @export
#' @importFrom purrr map exec
#' @importFrom rlang !!! abort warn inform arg_match cnd_signal .data
#' @importFrom glue glue
#' @importFrom tictoc tic toc
#' @importFrom checkmate assertString assertInt assertFlag assertCount
#' @importFrom checkmate assert checkList checkCharacter
learnTwasWeights <- function(
    X,
    Y,
    weightMethods,
    study = "",
    context = "",
    trait = "",
    numThreads = 1,
    fittedModels = NULL,
    fitRetention = c("none", "slim", "full"),
    standardized = FALSE,
    dataType = NULL,
    ldSketch = NULL,
    verbose = 1,
    seed = NULL
) {
    assertString(study)
    assertString(context)
    assertString(trait)
    assertInt(numThreads)
    fitRetention <- arg_match(fitRetention)
    assertFlag(standardized)
    assertString(dataType, null.ok = TRUE)
    assertCount(verbose)
    assertInt(seed, null.ok = TRUE)
    # weightMethods is documented as a named list OR a character vector.
    assert(
        checkList(weightMethods),
        checkCharacter(weightMethods),
        .var.name = "weightMethods"
    )
    .learnTwasWeightsImpl(
        X = X,
        Y = Y,
        weightMethods = weightMethods,
        study = study,
        context = context,
        trait = trait,
        numThreads = numThreads,
        fittedModels = fittedModels,
        fitRetention = fitRetention,
        standardized = standardized,
        dataType = dataType,
        ldSketch = ldSketch,
        verbose = verbose,
        seed = seed
    )
}

# Validate X/Y shapes; coerce a vector Y to a one-column matrix. Returns Y.
# @noRd
#' @importFrom checkmate assert assertMatrix checkAtomicVector checkMatrix
.twasValidateXY <- function(X, Y) {
    assertMatrix(X)
    assert(checkMatrix(Y), checkAtomicVector(Y), .var.name = "Y")
    if (is.vector(Y)) {
        Y <- matrix(Y, ncol = 1)
    }
    assertMatrix(Y, nrows = nrow(X))
    Y
}

# Number of BiocParallel workers to use: -1 means all available, otherwise the
# requested count capped at what is available.
# @noRd
.twasResolveCores <- function(numThreads) {
    avail <- multicoreWorkers()
    min(if (numThreads == -1) avail else numThreads, avail)
}

# Variant ids for the weight rows: colnames(X), or synthetic variant_i labels.
# @noRd
.twasVariantIds <- function(X) {
    # An unnamed genotype matrix carries no variant identity, and a synthetic
    # "variant_<i>" label does not create one -- it just defers the failure to
    # wherever the range is needed. A variant id renders (chrom, pos, ref, alt),
    # so the caller has to supply real ids as colnames.
    if (is.null(colnames(X))) {
        msg <- glue(
            "twasWeights: the genotype matrix has no colnames, so its ",
            "{ncol(X)} variants have no identity. Set colnames(X) to variant ",
            "ids of the form chrom:pos:ref:alt."
        )
        abort(msg)
    }
    colnames(X)
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
    set_names(weightsList, names(weightMethods))
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
# Whether a fine-mapping token's fit is already available -- on the method's
# own arguments, or in `fittedModels`.
# @noRd
.twasLearnFitAvailable <- function(tk, resolvedMethods, fittedModels) {
    adapter <- .twasFineMappingMethodAdapters[[tk]]
    args <- .twasMethodArgsFor(resolvedMethods, tk)
    !is.null(args[[adapter$fitArg]]) || !is.null(fittedModels[[tk]])
}

# Refuse a fine-mapping token that learnTwasWeights() has no fit for, naming
# both ways one can be supplied.
# @noRd
.twasRequireLearnFits <- function(resolvedMethods, fittedModels) {
    .twasRequireSusieFits(
        resolvedMethods,
        available = partial(
            .twasLearnFitAvailable,
            resolvedMethods = resolvedMethods,
            fittedModels = fittedModels
        ),
        fnLabel = "learnTwasWeights",
        how = str_c(
            "Pass it as `fittedModels = list(<token> = <fit>)`, or on the ",
            "method's own arguments; run fineMappingPipeline() first to ",
            "produce one."
        )
    )
}

.learnTwasWeightsImpl <- function(
    X,
    Y,
    weightMethods,
    study,
    context,
    trait,
    numThreads,
    fittedModels,
    fitRetention,
    standardized,
    dataType,
    ldSketch,
    verbose,
    seed
) {
    .applySeed(seed)
    # fitRetention arrives already matched: learnTwasWeights is the only
    # caller and validates it there, where the choices are declared.
    Y <- .twasValidateXY(X, Y)
    resolvedMethods <- if (is.character(weightMethods)) {
        .twasMethodLookup(weightMethods)
    } else {
        weightMethods
    }
    .twasRequireLearnFits(resolvedMethods, fittedModels)
    validColumns <- .nonzeroVarColumns(X)
    Xfiltered <- as.matrix(X[, validColumns, drop = FALSE])
    prepared <- .prepareSusieWeightMethods(resolvedMethods, fittedModels)
    ctx <- list(
        X = X,
        Y = Y,
        Xfiltered = Xfiltered,
        validColumns = validColumns,
        study = study,
        context = context,
        trait = trait,
        fitRetention = fitRetention,
        standardized = standardized,
        dataType = dataType,
        verbose = verbose,
        rngSeed = seed
    )
    weightsList <- .twasFitAllMethods(
        prepared,
        ctx,
        .twasResolveCores(numThreads)
    ) |>
        .twasApplyRownames(X)
    rows <- .buildTwasWeightEntries(weightsList, .twasVariantIds(X), ctx)
    TwasWeights(
        study = rows$study,
        context = rows$context,
        trait = rows$trait,
        method = rows$method,
        entry = rows$entry,
        ldSketch = ldSketch
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
#' colnames(X) <- sprintf("chr1:%d:A:G", 100L * (1:4))
#' twe <- twasWeightsRow(variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5), standardized = FALSE)
#' tw <- TwasWeights(study = "s1", context = "brain", trait = "g1",
#'   method = "susie", entry = list(twe))
#' twasPredict(X, tw)
#' @importFrom checkmate assert checkList checkClass
#' @importFrom checkmate assertList
twasPredict <- function(X, weightsList) {
    # The body branches on TwasWeights, so this is a list OR that S4 class.
    assert(
        checkList(weightsList),
        checkClass(weightsList, "TwasWeights"),
        .var.name = "weightsList"
    )
    if (is(weightsList, "TwasWeights")) {
        # Per-row weights vector/matrix payloads. Use the method name as key
        # for compatibility with the legacy snake_case "<method>_predicted"
        # convention; ensembleWeights() rebinds the suffix.
        methodNames <- as.character(weightsList$method)
        # The per-row weights live in each element's mcols now, so read them
        # off the elements rather than a stored `entry` column.
        wl <- set_names(
            map(seq_len(nrow(weightsList)), .twasRowWeights, x = weightsList),
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
#'   \code{fitRetention = "slim"} in \code{learnTwasWeights} to obtain
#'   this).
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
                "with fitRetention = \"slim\" and ensure 'mrash' is in ",
                "the method list."
            )
            abort(msg)
        }
        fit <- getFits(.twrRowParts(weightResults, idx[[1L]]))
        if (is.null(fit) || is.null(fit$pi)) {
            msg <- glue(
                "mr.ash fit object not found. Run learnTwasWeights() with ",
                "fitRetention = \"slim\" ",
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
                "fitRetention = \"slim\" ",
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
    !is.character(.tupleColumn(object, jc))
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
    .tupleColumn(object, cn)
}

# Canonical short method name: map a full function name back to its short name.
# @noRd
.twasCanonicalShortName <- function(m, fnToShort) {
    if (is_in(m, names(fnToShort))) fnToShort[[m]] else m
}

# The per-method args list for short name `m`, tagged with its `impl` attribute.
# @noRd
.twasMethodArgsWithImpl <- function(m) {
    `attr<-`(
        .twasMethodMap[[m]]$args,
        "impl",
        .twasMethodMap[[m]]$impl
    )
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

# Route a caller's per-method arguments into a weight function. A wrapper's
# own formals (a pre-fit, fitRetention, initPriorSd, ...) bind by name; any
# else is a tool option and goes in the wrapper's `methodArgs` list, so an
# unknown option errors inside the wrapper rather than vanishing. Wrappers
# with no `methodArgs` formal get everything by name, and R reports an unused
# argument -- which is the point.
# @noRd
.twasWeightCallArgs <- function(fnName, baseArgs, userArgs) {
    c(baseArgs, .splitMethodArgs(fnName, userArgs))
}

# One univariate fold's weight column for outcome `k` (quiet unless verbose).
# @noRd
.twasFitColWeight <- function(k, ctx, fnName, Xtr, Ytr, args) {
    callArgs <- .twasWeightCallArgs(fnName, list(X = Xtr, y = Ytr[, k]), args)
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
    named <- `rownames<-`(x, colnames(X))
    if (is.null(fit)) {
        return(named)
    }
    `attr<-`(named, "fit", fit)
}

# X %*% w for one weight vector/matrix (coerced to a 1-column matrix if needed).
# @noRd
.twasPredictOne <- function(w, X) {
    if (!is.matrix(w)) {
        w <- matrix(w, ncol = 1)
    }
    X %*% w
}
