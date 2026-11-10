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
# RangedTupleList collection keyed by (study, context, trait[, joint*]): one
# row per gene, holding every method fitted for it (see "Per-gene rows"
# below). Each element is the union of the row's methods' variants; the
# weights, fits, cross-validation and scale are columns of the collection's
# own mcols, kept per method. Class-level slots:
#   * ldReference   the LD references RSS-derived weights used, each stored
#                once; the ldReference mcols column names each row's (NA
#                for weights learned from individual-level data).
# Constructor + accessors below. The twasWeights pipeline helpers
# (learnTwasWeights, CV, ensemble, etc.) follow at the bottom.
# =============================================================================

#' @include AllGenerics.R tupleSelectors.R
NULL

#' @title TWAS Weights Collection
#' @description S4 collection of TWAS weights, one row per
#'   \code{(study, context, trait)} -- a gene -- holding the weights of every
#'   method fitted for it. Each element is the union of the row's methods'
#'   variants; a method fitted on fewer of them, or in another order, keeps
#'   its own. The weights form a variants x methods x weight-column tensor,
#'   stored dense for the methods most of whose weights are nonzero and
#'   sparse for the rest. One method of one gene is what a
#'   \code{\link{TwasWeightsRow}} holds: the constructor takes one per
#'   \code{(study, context, trait, method)}, and the accessors select one by
#'   \code{method}.
#'
#' Required \code{mcols()} columns: the identity \code{study},
#' \code{context} and \code{trait}, and each row's \code{methods} with
#' their weights (\code{weightTensor}), fits (\code{fits}, and the
#' per-fold fits of their cross-validation in \code{foldFits}),
#' \code{cvResult}, \code{standardized} and \code{dataType}.
#'
#' Optional columns \code{jointStudies}, \code{jointContexts},
#' \code{jointTraits} appear when the collection holds weights from a
#' \code{jointSpecification}-driven joint fit: they list the semicolon-joined
#' members of the joined axis, \code{NA_character_} for a row from no joint
#' fit. A row is unique on its identity together with any joint columns, so
#' one gene's weights from a joint fit and from fits of its own are two rows.
#' @slot ldReference A named list of the LD references the weights were
#'   derived against, each a genotype panel, a \code{GenotypeMixture} or an
#'   \code{LdMatrix}, stored once apiece; the \code{ldReference} column of
#'   \code{mcols} names the entry each row used (\code{NA} for weights
#'   learned from individual-level data). Used downstream to check that the
#'   rows a pipeline combines share an LD reference.
#' @slot cvPartition A named list of the cross-validation fold partitions the
#'   rows' results share, each with the order of the samples its out-of-fold
#'   predictions are held for, stored once apiece; the \code{cvPartition}
#'   column of \code{mcols} names the entry each row's result used
#'   (\code{NA} for a row without one). \code{cvResult()} returns a row's
#'   result with its partition and sample names in place.
#' @slot fits A named list of the methods' fits, and of the per-fold fits
#'   of their cross-validation, stored once apiece: the rows of a joint fit
#'   all name one entry. The \code{fits} and \code{foldFits} columns of
#'   \code{mcols} name each method's entries (\code{NA} for none);
#'   \code{methodFits()} and \code{cvResult()} return them in place.
#' @section Narrowed rows: A row's stored fit names the variants its weights
#'   were fitted on, where its engine keeps them or in a \code{variantIds}
#'   element. A row narrowed since keeps its fits as they were, so
#'   \code{cvResult()} and \code{methodFits()} warn, with condition class
#'   \code{pecotmrFitCoverage}, when a row holds fewer variants than they
#'   describe. So do the pipelines that read them:
#'   \code{causalInferencePipeline()} when it selects weights by their
#'   cross-validation, and \code{ctwasPipeline()} for rows whose SuSiE fit
#'   renormalizes their weights. The weights themselves are narrowed with
#'   the row and do not warn, and a row kept without a fit
#'   (\code{fitRetention = "none"}) has nothing to compare them with.
#' @export
setClass(
    "TwasWeights",
    contains = "RangedTupleList",
    representation(ldReference = "list", cvPartition = "list", fits = "list"),
    prototype(ldReference = list(), cvPartition = list(), fits = list()),
    validity = function(object) .validateTwasWeights(object)
)

# Validity for the TwasWeights collection: required columns, entry payloads,
# region/traitPos provenance, joint* column types, tuple uniqueness, and the
# optional ldReference. Returns TRUE or a character vector of error messages.
# @noRd
#' @importFrom checkmate makeAssertCollection assertNames
.validateTwasWeights <- function(object) {
    coll <- makeAssertCollection()
    assertNames(
        .tupleColumnNames(object),
        must.include = c("study", "context", "trait"),
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
        .validateRegionColumn(object),
        .twasValidateJointCols(object, jointCols),
        .twasValidateKeyUniqueness(object, jointCols),
        .twasValidateMethods(object),
        .ldRefDictProblems(object@ldReference, .ldRefKeys(object)),
        .cvPartitionProblems(object),
        .fitProblems(object, c("fits", "foldFits"))
    )
}

# A row holds each of its methods once.
# @noRd
.twasValidateMethods <- function(object) {
    methods <- map(as.list(.tupleColumn(object, "methods")), as.character)
    if (any(map_lgl(methods, .twasHasDuplicate))) {
        return("a row holds one set of weights per method")
    }
    character()
}

# @noRd
.twasHasDuplicate <- function(v) {
    anyDuplicated(v) > 0L
}

# The per-entry payload columns. The variants and their weights are the
# elements now, so what is left to check is that the payload columns are
# present and parallel to them.
# @noRd
.twasValidateEntries <- function(object) {
    payload <- c(
        "methods",
        "weightTensor",
        "fits",
        "foldFits",
        "cvResult",
        "standardized",
        "dataType"
    )
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

# (study, context, trait[, joint*]) uniqueness: one row per gene.
# @noRd
.twasValidateKeyUniqueness <- function(object, jointCols) {
    keyCols <- c("study", "context", "trait", jointCols)
    # Read column by column: `[` takes no column index on a collection.
    keyTbl <- as_tibble(set_names(
        map(keyCols, .twasColumn, object = object),
        keyCols
    ))
    if (nrow(distinct(keyTbl)) < nrow(keyTbl)) {
        str_c(
            "(study, context, trait[, joint*]) uniqueness violated: one row ",
            "per gene"
        )
    } else {
        character()
    }
}

# =============================================================================

#' @title Create a TwasWeights Collection Object
#' @description Construct a \code{TwasWeights} collection from per-entry
#'   vectors and a list of \code{TwasWeightsRow} payloads, one per (study,
#'   context, trait, method) entry. The entries of one gene -- one (study,
#'   context, trait) and joint group -- become one row, its methods an axis.
#' @param studyName Character vector of study identifiers. A joint fit's
#'   entries keep their own study, context and trait, and record the group in
#'   \code{jointStudies}, \code{jointContexts} or \code{jointTraits}.
#' @param context Character vector of context labels.
#' @param trait Character vector of trait identifiers.
#' @param method Character vector of TWAS weight method names. A snake_case
#'   name is stored as its camelCase token (\code{susie_inf} as
#'   \code{susieInf}); a gene holds one entry per method.
#' @param entry List / \code{SimpleList} of \code{TwasWeightsRow} objects.
#' @param jointStudies Optional character vector (one per entry) listing the
#'   semicolon-joined studies participating in each entry's cross-study joint
#'   fit, or \code{NA_character_} for non-joint entries. When \code{NULL}
#'   (default) the column is omitted.
#' @param jointContexts Optional character vector for cross-context joints. Same
#'   shape as \code{jointStudies}.
#' @param jointTraits Optional character vector for cross-trait joints. Same
#'   shape as \code{jointStudies}.
#' @param ldReference The LD reference the weights were derived against --
#'   a genotype panel (see \code{\link{readGenotypes}}), a
#'   \code{GenotypeMixture} or an \code{LdMatrix} -- used for every entry,
#'   or a list with one per entry; \code{NULL} for individual-level fits.
#'   Each distinct reference is stored once, and a gene's methods share one.
#' @param traitPos Optional per-entry trait genomic anchor (a \code{GRanges}
#'   or \code{NULL}), carried forward as provenance; not part of the
#'   identity key. \code{NULL} (default) omits the column.
#' @param region Optional per-entry analysis WINDOW (a \code{GRanges} or
#'   \code{NULL}) -- the interval each fit was run over, carried forward as
#'   provenance; not part of the identity key. A gene's methods share one.
#'   \code{NULL} (default) omits the column, and \code{\link{genomicRegion}}
#'   then derives the window from the span of the variants present.
#' @param metadata A list of free-form metadata about the collection,
#'   as \code{metadata()} holds on any Bioconductor vector. Default
#'   \code{list()}.
#' @return A \code{TwasWeights} object.
#' @examples
#' twe <- TwasWeightsRow(variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5),
#'   weightStandardized = FALSE)
#' tw <- TwasWeights(studyName = "s1", context = "brain", trait = "gene1",
#'   method = "susie", entry = list(twe))
#' tw
#' @importFrom checkmate assertCharacter assert checkList
#' @importFrom checkmate checkClass
#' @export
TwasWeights <- function(
    studyName,
    context,
    trait,
    method,
    entry,
    jointStudies = NULL,
    jointContexts = NULL,
    jointTraits = NULL,
    traitPos = NULL,
    region = NULL,
    ldReference = NULL,
    metadata = list()
) {
    .twasAssertArgs(
        studyName,
        context,
        trait,
        method,
        entry,
        jointStudies,
        jointContexts,
        jointTraits
    )
    n <- .twasCheckRowLengths(studyName, context, trait, method, entry)
    entry <- map(entry, .asTwRowPayload)
    .checkRowPayloads(entry, "TwasWeightsRow", "TWAS-weight")
    cells <- list(
        keys = .twasCellKeys(
            studyName,
            context,
            trait,
            list(
                jointStudies = jointStudies,
                jointContexts = jointContexts,
                jointTraits = jointTraits
            ),
            n
        ),
        method = .twasMethodToken(as.character(method)),
        entry = entry,
        traitPos = .twasCheckedRanges(traitPos, n, "traitPos"),
        region = .twasCheckedRanges(region, n, "region"),
        ldRefs = .ldRefPerTuple(ldReference, n)
    )
    .withMetadata(.twasNewValidated(cells), metadata)
}

# The identity and joint-axis arguments, checked. Mirrors
# .qfmrAssertArgs(); `entry` is documented as "List / SimpleList", and
# SimpleList is S4 and fails checkList, so that one is an or-combination.
# @noRd
.twasAssertArgs <- function(
    studyName,
    context,
    trait,
    method,
    entry,
    jointStudies,
    jointContexts,
    jointTraits
) {
    assertCharacter(studyName, any.missing = FALSE)
    assertCharacter(context, any.missing = FALSE)
    assertCharacter(trait, any.missing = FALSE)
    assertCharacter(method, any.missing = FALSE)
    assert(
        checkList(entry),
        checkClass(entry, "SimpleList"),
        .var.name = "entry"
    )
    assertCharacter(jointStudies, null.ok = TRUE)
    assertCharacter(jointContexts, null.ok = TRUE)
    assertCharacter(jointTraits, null.ok = TRUE)
}

# Every cell's identity: study, context and trait, and the joint axes
# supplied (each one value per cell).
# @noRd
.twasCellKeys <- function(studyName, context, trait, joints, n) {
    supplied <- compact(joints)
    walk2(names(supplied), supplied, .twasCheckJointColLength, n = n)
    as_tibble(c(
        list(
            study = as.character(studyName),
            context = as.character(context),
            trait = as.character(trait)
        ),
        map(supplied, as.character)
    ))
}

# A per-cell range argument (`traitPos`, `region`): NULL, or a GRanges with
# one range per cell.
# @noRd
.twasCheckedRanges <- function(gr, n, what) {
    if (is.null(gr)) {
        return(NULL)
    }
    if (!methods::is(gr, "GRanges")) {
        abort(glue("`{what}` must be a GRanges (one range per row) or NULL."))
    }
    if (length(gr) != n) {
        abort(glue("`{what}` must have the same length as `study`."))
    }
    gr
}

# Require study/context/trait/method/entry to share one length; returns it.
# @noRd
.twasCheckRowLengths <- function(studyName, context, trait, method, entry) {
    n <- length(studyName)
    if (
        length(context) != n ||
            length(trait) != n ||
            length(method) != n ||
            length(entry) != n
    ) {
        msg <- glue(
            "`studyName`, `context`, `trait`, `method`, and `entry` must all ",
            "have the same length."
        )
        abort(msg)
    }
    n
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

# =============================================================================
# Per-gene rows
# -----------------------------------------------------------------------------
# A TwasWeights row is one (study, context, trait[, joint*]) -- a gene -- and
# holds every method fitted for it. Its element is the union of its methods'
# variants. Their weights form a variants x methods x K tensor, K being a
# method's weight columns (one for a weight vector), held in two parts: a
# dense array for the methods most of whose weights are nonzero, an
# SVT_SparseArray for the rest, each method placed by its measured density. A
# method fitted on other variants than the union, or in another order,
# records its own as positions in the union. Fits, cross-validation, scale
# and data type are kept per method. One method of one row -- a cell -- is
# what a TwasWeightsRow holds, what the constructor takes and what the
# accessors read (.twrCellParts()).
# =============================================================================

# Weight methods are named by their camelCase token (susieInf, bayesR), the
# one spelling a row's method axis holds; a snake_case label such as
# learnTwasWeights() derives from a method function (susie_inf) maps onto it.
# @noRd
.twasMethodToken <- function(method) {
    str_replace_all(method, "_[a-z0-9]", .twasUpperTail)
}

# @noRd
.twasUpperTail <- function(x) {
    str_to_upper(str_sub(x, 2L))
}

# `expr`, without the warning SparseArray raises after its calls on builds
# linked to LLVM's libomp: it pauses the OpenMP runtime after every call, and
# libomp reports a pause of a runtime already paused as a failure. Nothing
# else is muffled.
# @noRd
.svtQuiet <- function(expr) {
    withCallingHandlers(expr, warning = .svtMuffleOmpPause)
}

# @noRd
.svtMuffleOmpPause <- function(w) {
    if (str_detect(conditionMessage(w), "omp_pause_resource\\(\\) failed")) {
        invokeRestart("muffleWarning")
    }
}

# ---- building rows from cells ------------------------------------------------

# The collection the cells `cells` build: one row per distinct identity, in
# order of first appearance. `cells` holds, one per cell, the identity
# (`keys`), the method, the TwasWeightsRow (`entry`), the trait position and
# region (GRanges, or NULL for none) and the LD reference (`ldRefs`).
# @noRd
.twasNewValidated <- function(cells) {
    byRow <- unname(split(
        seq_along(cells$method),
        as.integer(vctrs::vec_group_id(cells$keys))
    ))
    .twasRefuseRepeatedMethods(byRow, cells)
    first <- map_int(byRow, .twasFirstCell)
    stores <- map(
        byRow,
        .twasRowStore,
        entries = cells$entry,
        methods = cells$method
    )
    .twasRefuseMultiChromosome(stores)
    fitted <- .twasStoredFits(
        map(cells$entry, methodFits),
        map(cells$entry, cvResult)
    )
    cv <- .cvStoredColumn(fitted$cv)
    dict <- .ldRefDict(.twasRowLdRefs(cells$ldRefs, byRow, first))
    cols <- c(
        .twasRowCols(cells, byRow, first, stores, cv, fitted),
        list(ldReference = dict$keys, check.names = FALSE)
    )
    grl <- S4Vectors::`mcols<-`(
        .rtlBindElements(map(stores, "variants")),
        value = exec(S4Vectors::DataFrame, !!!cols)
    )
    new(
        "TwasWeights",
        .rtlNamed(grl),
        ldReference = dict$entries,
        cvPartition = cv$entries,
        fits = fitted$entries
    )
}

# The cells' fits `fits` and the per-fold fits of their CV `cvs` as the fits
# dictionary holds them, each distinct one once, with each cell's keys into
# it; and the cells' CV without the per-fold fits, which the dictionary now
# holds.
# @noRd
.twasStoredFits <- function(fits, cvs) {
    folds <- map(cvs, .cvFoldFits)
    dict <- .dictBuild(c(fits, folds), .FIT_DICT_PREFIX)
    n <- length(fits)
    list(
        fits = dict$keys[seq_len(n)],
        foldFits = dict$keys[n + seq_len(n)],
        entries = dict$entries,
        cv = map2(cvs, folds, .cvWithoutFoldFits)
    )
}

# A CV payload's per-fold fits (mr.mash keeps them), NULL for none.
# @noRd
.cvFoldFits <- function(cv) {
    if (is.list(cv) && !is.data.frame(cv)) cv[["foldFits"]]
}

# @noRd
.cvWithoutFoldFits <- function(cv, folds) {
    if (is.null(folds)) cv else list_assign(cv, foldFits = zap())
}

# @noRd
.twasFirstCell <- function(idx) {
    idx[[1L]]
}

# A gene holds one set of weights per method.
# @noRd
.twasRefuseRepeatedMethods <- function(byRow, cells) {
    repeated <- keep(byRow, .twasHasRepeatedMethod, methods = cells$method)
    if (length(repeated) == 0L) {
        return(invisible(NULL))
    }
    i <- repeated[[1L]][[1L]]
    key <- str_flatten(discard(map_chr(cells$keys, i), is.na), "|")
    abort(glue(
        "TwasWeights: '{key}' has more than one set of weights for one ",
        "method; a (study, context, trait, method[, joint*]) holds one."
    ))
}

# @noRd
.twasHasRepeatedMethod <- function(idx, methods) {
    anyDuplicated(methods[idx]) > 0L
}

# One gene's store: its methods' union of variants, their weights as a
# tensor, and their scale and data type by method.
# @noRd
.twasRowStore <- function(idx, entries, methods) {
    rowEntries <- set_names(entries[idx], methods[idx])
    union <- .twasUnionVariants(rowEntries)
    ids <- .grVariantIds(union)
    placed <- map(rowEntries, .twasPlacedWeights, ids = ids)
    list(
        variants = union,
        tensor = .twasTensor(placed, length(ids)),
        standardized = map_lgl(rowEntries, weightStandardized),
        dataType = map(rowEntries, dataType)
    )
}

# The variants of every method of a row, each once, in order of first
# appearance, carrying their identity columns only.
# @noRd
.twasUnionVariants <- function(rowEntries) {
    grs <- map(unname(rowEntries), .twasIdentityVariants)
    flat <- suppressWarnings(exec(c, !!!grs))
    flat[!duplicated(.grVariantIds(flat))]
}

# @noRd
.twasIdentityVariants <- function(e) {
    gr <- variants(e)
    md <- mcols(gr, use.names = FALSE)
    S4Vectors::`mcols<-`(
        gr,
        value = md[, colnames(md) != "weight", drop = FALSE]
    )
}

# One method's weights on the row's variants `ids`: its own variants' rows
# filled and zero elsewhere, with what reading them back needs. `positions`
# is NULL when the method holds the row's variants in their order.
# @noRd
.twasPlacedWeights <- function(e, ids) {
    own <- .grVariantIds(variants(e))
    hit <- match(own, ids)
    w <- .unnamedWeights(weights(e))
    list(
        w = if (is.null(w)) NULL else .twasFillRows(w, hit, length(ids)),
        positions = if (identical(hit, seq_along(ids))) NULL else hit,
        isMatrix = is.matrix(w),
        columns = colnames(w)
    )
}

# The row's methods' weights as a tensor in two parts, and where each method
# sits in it: `part` ("dense", "sparse", or "none" for a method without
# weights), `slice` along the method axis of that part, and `K`.
# @noRd
#' @importClassesFrom SparseArray SVT_SparseArray
.twasTensor <- function(placed, n) {
    K <- map_int(placed, .twasPlacedK)
    part <- rep("none", length(K))
    part[K > 0L] <- "sparse"
    part[K > 0L & map_lgl(placed, .twasIsDense)] <- "dense"
    slice <- integer(length(K))
    slice[part == "dense"] <- seq_len(sum(part == "dense"))
    slice[part == "sparse"] <- seq_len(sum(part == "sparse"))
    list(
        dense = .twasPartArray(placed[part == "dense"], n),
        sparse = .twasSparsePart(placed[part == "sparse"], n),
        part = unname(part),
        slice = unname(slice),
        K = unname(K),
        isMatrix = unname(map_lgl(placed, "isMatrix")),
        columns = discard(map(placed, "columns"), is.null),
        variants = discard(map(placed, "positions"), is.null)
    )
}

# `w` (a vector or a variants x K matrix) as rows `hit` of an n x K matrix of
# zeros.
# @noRd
.twasFillRows <- function(w, hit, n) {
    m <- if (is.matrix(w)) w else matrix(w, ncol = 1L)
    full <- matrix(0, n, ncol(m))
    full[hit, ] <- m
    full
}

# @noRd
.twasPlacedK <- function(p) {
    if (is.null(p$w)) 0L else ncol(p$w)
}

# Dense when at least two thirds of the weights are nonzero: a stored
# nonzero costs about 12 bytes in an SVT_SparseArray, a dense one 8.
# @noRd
.twasIsDense <- function(p) {
    nonzero <- sum(p$w != 0 | is.na(p$w))
    3 * nonzero >= 2 * length(p$w)
}

# The methods `placed` as one variants x methods x K array, K padded with
# zeros to the widest; NULL for none.
# @noRd
.twasPartArray <- function(placed, n) {
    if (length(placed) == 0L) {
        return(NULL)
    }
    K <- max(map_int(placed, .twasPlacedK))
    arr <- array(0, c(n, length(placed), K))
    for (j in seq_along(placed)) {
        arr[, j, seq_len(ncol(placed[[j]]$w))] <- placed[[j]]$w
    }
    arr
}

# @noRd
.twasSparsePart <- function(placed, n) {
    arr <- .twasPartArray(placed, n)
    if (is.null(arr)) {
        return(NULL)
    }
    .svtQuiet(methods::as(arr, "SVT_SparseArray"))
}

# Each row's LD reference: the one its cells share. A gene's methods were fit
# against one reference; cells on different ones are refused.
# @noRd
.twasRowLdRefs <- function(ldRefs, byRow, first) {
    ids <- .dictIdentities(ldRefs)
    ids[map_lgl(ldRefs, is.null)] <- NA_character_
    mixed <- keep(byRow, .twasRowDisagrees, values = ids)
    if (length(mixed) > 0L) {
        abort(glue(
            "TwasWeights: the methods of one (study, context, trait) were fit ",
            "against different LD references; a row holds one."
        ))
    }
    ldRefs[first]
}

# @noRd
.twasRowDisagrees <- function(idx, values) {
    length(unique(values[idx])) > 1L
}

# The outer columns of the rows `byRow` groups the cells into: the identity of
# each row's first cell, its methods and their payloads, and the trait
# position and region its cells share.
# @noRd
.twasRowCols <- function(cells, byRow, first, stores, cv, fitted) {
    n <- length(byRow)
    payload <- list(
        methods = IRanges::CharacterList(
            map(byRow, .twasPick, values = cells$method)
        ),
        weightTensor = S4Vectors::SimpleList(map(stores, "tensor")),
        fits = IRanges::CharacterList(
            map(byRow, .twasPick, values = fitted$fits)
        ),
        foldFits = IRanges::CharacterList(
            map(byRow, .twasPick, values = fitted$foldFits)
        ),
        cvResult = S4Vectors::SimpleList(
            map(byRow, .twasRowCv, cv = cv, methods = cells$method)
        ),
        cvPartition = IRanges::CharacterList(
            map(byRow, .twasPick, values = cv$keys)
        ),
        standardized = S4Vectors::SimpleList(map(stores, "standardized")),
        dataType = S4Vectors::SimpleList(map(stores, "dataType"))
    )
    # Assigned, not listed: a NULL in the list would still name a column.
    payload$jointIndex <- .twasRowJointIndex(cells$jointIndex, byRow)
    traitPos <- .twasRowRanges(cells$traitPos, byRow, "traitPos")
    region <- .twasRowRanges(cells$region, byRow, "region")
    c(as.list(cells$keys[first, ]), payload) |>
        .appendTraitPosCol(traitPos, n) |>
        .appendRegionCol(region, n)
}

# The rows' cells' joint-fit positions, or NULL (no column) when no cell is
# from a joint fit.
# @noRd
.twasRowJointIndex <- function(jointIndex, byRow) {
    if (is.null(jointIndex) || all(is.na(jointIndex))) {
        return(NULL)
    }
    IRanges::IntegerList(map(byRow, .twasPick, values = jointIndex))
}

# @noRd
.twasPick <- function(idx, values) {
    values[idx]
}

# A row's stored CV, by method.
# @noRd
.twasRowCv <- function(idx, cv, methods) {
    set_names(as.list(cv$cvResult)[idx], methods[idx])
}

# The range a row's cells share for `what` (trait position or region): the
# placeholder of a cell given none is unknown, and the cells that know one
# must agree. NULL when no cell has the column.
# @noRd
.twasRowRanges <- function(gr, byRow, what) {
    if (is.null(gr)) {
        return(NULL)
    }
    known <- !.twasIsPlaceholder(gr)
    labels <- replace(as.character(gr), !known, NA_character_)
    mixed <- keep(byRow, .twasRowDisagreesKnown, values = labels)
    if (length(mixed) > 0L) {
        abort(glue(
            "TwasWeights: the methods of one (study, context, trait) carry ",
            "different `{what}` values; a row holds one."
        ))
    }
    gr[map_int(byRow, .twasFirstKnown, known = known)]
}

# The placeholder range a combine pads a missing range column with.
# @noRd
.twasIsPlaceholder <- function(gr) {
    as.character(GenomicRanges::seqnames(gr)) == "chrUn" & width(gr) == 0L
}

# @noRd
.twasRowDisagreesKnown <- function(idx, values) {
    length(unique(discard(values[idx], is.na))) > 1L
}

# @noRd
.twasFirstKnown <- function(idx, known) {
    hit <- idx[known[idx]]
    if (length(hit) == 0L) idx[[1L]] else hit[[1L]]
}

# ---- reading cells -----------------------------------------------------------

# Every cell of `x`, one per (row, method), in row order: its row, its method
# and its row's identity columns.
# @noRd
.twasCells <- function(x) {
    md <- mcols(x, use.names = FALSE)
    methods <- map(as.list(md$methods), as.character)
    rows <- rep(seq_along(methods), lengths(methods))
    ident <- intersect(.TWAS_KEY_COLUMNS, colnames(md))
    identCols <- set_names(
        map(ident, .twasRowValues, md = md, rows = rows),
        ident
    )
    tibble(
        row = rows,
        method = list_c(methods, ptype = character()),
        !!!identCols
    )
}

# The identity columns a row is keyed by.
# @noRd
.TWAS_KEY_COLUMNS <- c(
    "study",
    "context",
    "trait",
    "jointStudies",
    "jointContexts",
    "jointTraits"
)

# @noRd
.twasRowValues <- function(cn, md, rows) {
    as.character(md[[cn]])[rows]
}

# One cell as a TwasWeightsRow: the method's variants and weights, and its
# fit, cross-validation, scale and data type.
# @noRd
.twrCellParts <- function(x, row, method) {
    md <- mcols(x, use.names = FALSE)
    tensor <- md$weightTensor[[row]]
    union <- .rtlGatherElements(x, row)
    positions <- tensor$variants[[method]] %||% seq_along(union)
    k <- match(method, as.character(md$methods[[row]]))
    w <- .twasCellWeights(tensor, k, method, positions)
    gr <- union[positions]
    new(
        "TwasWeightsRow",
        variants = S4Vectors::`mcols<-`(
            gr,
            value = `[[<-`(mcols(gr, use.names = FALSE), "weight", value = w)
        ),
        weights = w,
        methodFits = .fitEntry(md$fits[[row]][[k]], x),
        cvResult = .twasCellCv(x, md, row, k),
        weightStandardized = isTRUE(md$standardized[[row]][[method]]),
        weightsDataType = md$dataType[[row]][[method]]
    )
}

# Method `method`'s weights (the k-th of its row) on its own variants: a
# vector, or a variants x K matrix when that is what it was given. K = 0 is
# either no weights (NULL) or a matrix with no columns, such as a joint fit
# that kept no conditions.
# @noRd
.twasCellWeights <- function(tensor, k, method, positions) {
    K <- tensor$K[[k]]
    if (K == 0L) {
        if (!tensor$isMatrix[[k]]) {
            return(NULL)
        }
        return(matrix(numeric(0), length(positions), 0L))
    }
    j <- tensor$slice[[k]]
    block <- if (identical(tensor$part[[k]], "dense")) {
        tensor$dense[positions, j, seq_len(K), drop = FALSE]
    } else {
        .svtQuiet(methods::as(
            tensor$sparse[positions, j, seq_len(K), drop = FALSE],
            "array"
        ))
    }
    w <- matrix(as.numeric(block), length(positions), K)
    if (!tensor$isMatrix[[k]]) {
        return(as.numeric(w))
    }
    `colnames<-`(w, tensor$columns[[method]])
}

# The k-th method of row `row`'s cross-validation, its shared partition and
# per-fold fits put back.
# @noRd
.twasCellCv <- function(x, md, row, k) {
    stored <- md$cvResult[[row]][[k]]
    key <- md$cvPartition[[row]][[k]]
    cv <- .cvJoinPartition(stored, .dictEntryFor(key, .cvPartitionDict(x)))
    folds <- .fitEntry(md$foldFits[[row]][[k]], x)
    if (is.null(folds)) cv else list_assign(cv, foldFits = folds)
}

# ---- selecting cells ---------------------------------------------------------

# The cell (a row of .twasCells()) a (study, context, trait, method) selector
# pins. All four are needed unless `x` holds a single cell; the first match
# wins, and the joint axes do not select, as for any collection.
# @noRd
#' @importFrom checkmate assertVector
.twrSelectCell <- function(x, studyName, context, trait, method) {
    cells <- .twasCells(x)
    if (nrow(cells) == 0L) {
        abort("TwasWeights has no rows.")
    }
    selectors <- list(studyName, context, trait, method)
    if (any(map_lgl(selectors, is.null))) {
        if (nrow(cells) == 1L) {
            return(cells[1L, ])
        }
        abort(glue(
            "TwasWeights has {nrow(cells)} entries. Pass `studyName`, ",
            "`context`, `trait`, and `method` to select one."
        ))
    }
    walk(selectors, assertVector, len = 1L)
    hit <- .twasCellIndex(cells, studyName, context, trait, method)
    if (length(hit) == 0L) {
        abort(glue(
            "No entry for (study='{studyName}', context='{context}', ",
            "trait='{trait}', method='{method}')."
        ))
    }
    cells[hit[[1L]], ]
}

# Indices of the cells (rows of .twasCells()) a (study, context, trait,
# method) selector matches, the method in either spelling; empty for none.
# `joint`, a joint key (.noJointKey(), .jointKeyAt()), narrows them to the
# cells of that joint group; without one the joint axes do not select.
# @noRd
.twasCellIndex <- function(
    cells,
    studyName,
    context,
    trait,
    method,
    joint = NULL
) {
    which(
        cells$study == studyName &
            cells$context == context &
            cells$trait == trait &
            cells$method == .twasMethodToken(method) &
            .cellsInJointGroup(cells, joint)
    )
}

# Which cells belong to joint group `joint`; every cell when none is named.
# @noRd
.cellsInJointGroup <- function(cells, joint) {
    if (is.null(joint)) {
        return(rep(TRUE, nrow(cells)))
    }
    reduce(
        map(.JOINT_COLUMNS, .cellsJointColumnIn, cells = cells, joint = joint),
        `&`
    )
}

# @noRd
.cellsJointColumnIn <- function(col, cells, joint) {
    .jointMembersMatch(.cellsJointColumn(col, cells), joint[[col]])
}

# Joint column `col` of a cell table, NA where the cells come from no joint
# fit or the table has no such column.
# @noRd
.cellsJointColumn <- function(col, cells) {
    if (is_in(col, names(cells))) {
        as.character(cells[[col]])
    } else {
        rep(NA_character_, nrow(cells))
    }
}

# The constructor's inputs for cells `cells` of `x` (rows of .twasCells()).
# @noRd
.twasCellInputs <- function(x, cells) {
    md <- mcols(x, use.names = FALSE)
    rows <- cells$row
    list(
        keys = cells[intersect(.TWAS_KEY_COLUMNS, names(cells))],
        method = cells$method,
        entry = map2(rows, cells$method, .twrCellParts, x = x),
        traitPos = if (is_in("traitPos", colnames(md))) md$traitPos[rows],
        region = if (is_in("region", colnames(md))) md$region[rows],
        ldRefs = .rowLdRefs(x)[rows],
        jointIndex = .twasCellJointIndex(md, cells)
    )
}

# Each cell's position in its joint fit's condition axis (NA for a cell from
# no joint fit), or NULL when no row of the collection is from one.
# @noRd
.twasCellJointIndex <- function(md, cells) {
    if (!is_in("jointIndex", colnames(md))) {
        return(NULL)
    }
    map2_int(cells$row, cells$method, .twasCellJointIndexAt, md = md)
}

# @noRd
.twasCellJointIndexAt <- function(row, method, md) {
    k <- match(method, as.character(md$methods[[row]]))
    as.integer(md$jointIndex[[row]][[k]])
}

# `x` cut to its cells `cells`, rebuilt so each row keeps only those.
# @noRd
.twasKeepCells <- function(x, cells) {
    .withMetadata(
        .twasNewValidated(.twasCellInputs(x, cells)),
        S4Vectors::metadata(x)
    )
}

# ---- combining ---------------------------------------------------------------

# Collections `parts` as one: their cells regrouped, so the methods one gene
# has in different parts share its row. A non-NULL `ldReference` is every
# row's reference instead.
# @noRd
.twasCombine <- function(parts, ldReference) {
    inputs <- map(parts, .twasAllCellInputs)
    joined <- list(
        keys = bind_rows(map(inputs, "keys")),
        method = list_c(map(inputs, "method"), ptype = character()),
        entry = list_flatten(map(inputs, "entry")),
        traitPos = .twasJoinRanges(map(inputs, "traitPos"), inputs),
        region = .twasJoinRanges(map(inputs, "region"), inputs),
        ldRefs = list_flatten(map(inputs, "ldRefs")),
        jointIndex = .twasJoinJointIndex(inputs)
    )
    if (!is.null(ldReference)) {
        joined$ldRefs <- .ldRefPerTuple(ldReference, length(joined$method))
    }
    .withMetadata(.twasNewValidated(joined), S4Vectors::metadata(parts[[1L]]))
}

# The parts' cells' joint-fit positions joined, NA for a part's cells from no
# joint fit; NULL when no part has any.
# @noRd
.twasJoinJointIndex <- function(inputs) {
    if (all(map_lgl(map(inputs, "jointIndex"), is.null))) {
        return(NULL)
    }
    list_c(map(inputs, .twasInputJointIndex), ptype = integer())
}

# @noRd
.twasInputJointIndex <- function(input) {
    input$jointIndex %||% rep(NA_integer_, length(input$method))
}

# @noRd
.twasAllCellInputs <- function(x) {
    .twasCellInputs(x, .twasCells(x))
}

# The parts' range columns joined, a part without the column padded with
# placeholders as a combine pads it; NULL when no part has it.
# @noRd
.twasJoinRanges <- function(grs, inputs) {
    if (all(map_lgl(grs, is.null))) {
        return(NULL)
    }
    n <- map_int(inputs, .twasInputCells)
    padded <- map2(grs, n, .twasRangesOrPlaceholder)
    paddedRanges <- unname(padded)
    suppressWarnings(exec(c, !!!paddedRanges))
}

# @noRd
.twasInputCells <- function(input) {
    length(input$method)
}

# @noRd
.twasRangesOrPlaceholder <- function(gr, n) {
    gr %||% .placeholderRanges(n)
}

# ---- rebuilt elements --------------------------------------------------------

# `md` (the outer columns of rows `keep` of `x`) with each row's tensor
# realigned to its rebuilt element in `bare`: weights follow their variants by
# id, a variant the element no longer holds drops out of every method, and one
# it gained has no weight.
# @noRd
.twasRealignedMcols <- function(x, md, bare, keep) {
    old <- map(keep, .twasElementIds, grl = x)
    new <- map(seq_along(bare), .twasElementIds, grl = bare)
    md$weightTensor <- S4Vectors::SimpleList(pmap(
        list(
            as.list(md$weightTensor),
            map(as.list(md$methods), as.character),
            old,
            new
        ),
        .twasRealignTensor
    ))
    md
}

# @noRd
.twasElementIds <- function(i, grl) {
    .grVariantIds(grl[[i]])
}

# One row's tensor, its methods `methods`, on variants `new` instead of
# `old`.
# @noRd
.twasRealignTensor <- function(tensor, methods, old, new) {
    if (identical(old, new)) {
        return(tensor)
    }
    at <- match(new, old)
    dense <- .twasTakeRows(tensor$dense, at)
    sparse <- .twasTakeRows(.twasDensePart(tensor$sparse), at)
    list_assign(
        tensor,
        dense = dense,
        sparse = if (is.null(sparse)) {
            NULL
        } else {
            .svtQuiet(methods::as(sparse, "SVT_SparseArray"))
        },
        variants = discard(
            map(
                set_names(methods),
                .twasRemappedPositions,
                tensor = tensor,
                n = length(old),
                at = at
            ),
            is.null
        )
    )
}

# @noRd
.twasDensePart <- function(part) {
    if (is.null(part)) NULL else .svtQuiet(methods::as(part, "array"))
}

# Rows `at` of array `arr` (zero where `at` is NA).
# @noRd
.twasTakeRows <- function(arr, at) {
    if (is.null(arr)) {
        return(NULL)
    }
    d <- dim(arr)
    out <- array(0, c(length(at), d[[2L]], d[[3L]]))
    kept <- which(!is.na(at))
    out[kept, , ] <- arr[at[kept], , , drop = FALSE]
    out
}

# Method `method`'s positions among a row's new variants: its old ones (all
# `n`, in order, when it held the row's variants as they were) mapped through
# `at`, those dropped gone; NULL when it holds the new variants in order.
# @noRd
.twasRemappedPositions <- function(method, tensor, n, at) {
    old <- tensor$variants[[method]] %||% seq_len(n)
    moved <- match(old, at)
    moved <- moved[!is.na(moved)]
    if (identical(moved, seq_along(at))) NULL else moved
}

# ---- flat views --------------------------------------------------------------

# `x` as flattenTupleRanges() hands a TwasWeights out: one range per (variant,
# method, weight column), carrying that weight, the row's identity columns
# (dotted) and `.method`. `.weightColumn` names a weight-matrix column (an
# fSuSiE feature, an outcome) and is NA for a plain weight vector.
# @noRd
.twasFlatten <- function(x) {
    md <- mcols(x, use.names = FALSE)
    cells <- .twasCells(x)
    blocks <- list_flatten(map2(
        cells$row,
        cells$method,
        .twasFlatCell,
        x = x,
        md = md,
        keyCols = .rtlTupleKeyCols(x)
    ))
    if (length(blocks) == 0L) {
        return(.rtlGatherElements(x, seq_len(length(x))))
    }
    blockRanges <- unname(blocks)
    suppressWarnings(exec(c, !!!blockRanges))
}

# @noRd
.twasFlatCell <- function(row, method, x, md, keyCols) {
    gr <- variants(.twrCellParts(x, row, method))
    w <- mcols(gr, use.names = FALSE)$weight
    labels <- .twbWeightLabels(w)
    ident <- set_names(
        map(keyCols, .twasRowValue, md = md, row = row),
        map_chr(keyCols, .rtlDotName)
    )
    map(
        seq_along(labels),
        .twasFlatBlock,
        gr = gr,
        w = w,
        labels = labels,
        ident = ident,
        method = method
    )
}

# @noRd
.twasRowValue <- function(nm, md, row) {
    md[[nm]][[row]]
}

# One weight column of one entry as flat ranges.
# @noRd
.twasFlatBlock <- function(k, gr, w, labels, ident, method) {
    n <- length(gr)
    own <- mcols(gr, use.names = FALSE)
    added <- c(
        list(
            weight = if (is.null(w)) {
                rep(NA_real_, n)
            } else {
                .twbWeightColumn(k, w)
            }
        ),
        map(ident, rep, length.out = n),
        list(.method = rep(method, n), .weightColumn = rep(labels[[k]], n))
    )
    `mcols<-`(
        gr,
        value = cbind(
            own[, colnames(own) != "weight", drop = FALSE],
            DataFrame(added, check.names = FALSE)
        )
    )
}

# The columns a flat TwasWeights needs to be nested back, beyond the
# identity columns: the weights and where they belong.
# @noRd
.TWAS_FLAT_COLUMNS <- c("weight", ".method", ".weightColumn")

# `flat` (a flattened `template`, possibly filtered or mutated) nested back
# into `template`'s entries: each entry keeps the ranges carrying its row's
# identity and its method, with their weights. An entry left with none is
# kept, with no variants.
# @noRd
.twasNest <- function(flat, template) {
    md <- mcols(template, use.names = FALSE)
    keyCols <- .rtlTupleKeyCols(template)
    cells <- .twasCells(template)
    fmd <- mcols(flat, use.names = FALSE)
    wanted <- str_c(
        .rtlTupleKeys(md, keyCols, length(template))[cells$row],
        cells$method,
        sep = "\r"
    )
    have <- str_c(
        .rtlTupleKeys(fmd, map_chr(keyCols, .rtlDotName), length(flat)),
        as.character(fmd$.method),
        sep = "\r"
    )
    inputs <- .twasCellInputs(template, cells)
    inputs$entry <- map2(
        inputs$entry,
        wanted,
        .twasNestCell,
        flat = flat,
        have = have,
        drop = c(map_chr(keyCols, .rtlDotName), .TWAS_FLAT_COLUMNS)
    )
    .withMetadata(.twasNewValidated(inputs), S4Vectors::metadata(template))
}

# Entry `entry` rebuilt from its ranges in `flat`, its payload kept.
# @noRd
.twasNestCell <- function(entry, key, flat, have, drop) {
    picked <- flat[which(have == key)]
    pmd <- mcols(picked, use.names = FALSE)
    ids <- .grVariantIds(picked)
    first <- !duplicated(ids)
    gr <- `mcols<-`(
        picked[first],
        value = pmd[first, setdiff(colnames(pmd), drop), drop = FALSE]
    )
    w <- .twasNestWeights(weights(entry), pmd, ids, first)
    new(
        "TwasWeightsRow",
        variants = `mcols<-`(
            gr,
            value = `[[<-`(mcols(gr, use.names = FALSE), "weight", value = w)
        ),
        weights = w,
        methodFits = methodFits(entry),
        cvResult = cvResult(entry),
        weightStandardized = weightStandardized(entry),
        weightsDataType = dataType(entry)
    )
}

# An entry's weights from its flat ranges, in the shape `shape` (its weights
# before flattening) had. A matrix weight a predicate dropped from one column
# while the variant kept another is 0, as a method's weight is wherever it
# has none.
# @noRd
.twasNestWeights <- function(shape, pmd, ids, first) {
    if (is.null(shape)) {
        return(NULL)
    }
    if (is.null(dim(shape))) {
        return(as.numeric(pmd$weight[first]))
    }
    labels <- .twbWeightLabels(shape)
    kept <- ids[first]
    w <- matrix(0, length(kept), length(labels))
    w[cbind(match(ids, kept), match(pmd$.weightColumn, labels))] <- pmd$weight
    `colnames<-`(w, colnames(shape))
}

# ---- saved collections -------------------------------------------------------

# A collection saved with one row per (study, context, trait, method), its
# methods folded into gene rows; one saved per gene is left as it is.
# @noRd
.twasFoldMethodRows <- function(object) {
    cols <- .tupleColumnNames(object)
    if (!is_in("method", cols) || is_in("methods", cols)) {
        return(object)
    }
    .withMetadata(
        .twasNewValidated(.twasLegacyCellInputs(object)),
        S4Vectors::metadata(object)
    )
}

# The constructor's inputs for a saved per-method collection, one cell per
# saved row.
# @noRd
.twasLegacyCellInputs <- function(object) {
    md <- mcols(object, use.names = FALSE)
    keyCols <- intersect(.TWAS_KEY_COLUMNS, colnames(md))
    list(
        keys = as_tibble(set_names(
            map(keyCols, .twasMdCharacter, md = md),
            keyCols
        )),
        method = .twasMethodToken(as.character(md$method)),
        entry = map(seq_along(object), .twasLegacyRow, x = object, md = md),
        traitPos = if (is_in("traitPos", colnames(md))) md$traitPos,
        region = if (is_in("region", colnames(md))) md$region,
        ldRefs = .rowLdRefs(object)
    )
}

# @noRd
.twasMdCharacter <- function(cn, md) {
    as.character(md[[cn]])
}

# Saved row `i` as the TwasWeightsRow it was built from: its weights the
# element's `weight` column, its payload the row's columns.
# @noRd
.twasLegacyRow <- function(i, x, md) {
    gr <- .rtlGatherElements(x, i)
    mc <- mcols(gr, use.names = FALSE)
    w <- .unnamedWeights(mc$weight)
    new(
        "TwasWeightsRow",
        variants = S4Vectors::`mcols<-`(
            gr,
            value = `[[<-`(mc, "weight", value = w)
        ),
        weights = w,
        methodFits = md$fits[[i]],
        cvResult = .cvJoinPartition(md$cvResult[[i]], .mdCvPartition(x, md, i)),
        weightStandardized = isTRUE(md$standardized[[i]]),
        weightsDataType = md$dataType[[i]]
    )
}

# A gene's methods sit on one chromosome: its row is one element.
# @noRd
.twasRefuseMultiChromosome <- function(stores) {
    spans <- map_int(stores, .twasChromosomeCount)
    if (all(spans <= 1L)) {
        return(invisible(NULL))
    }
    abort(glue(
        "TwasWeights: row {which(spans > 1L)[[1L]]} holds variants on more ",
        "than one chromosome; a (study, context, trait) is one cis fit."
    ))
}

# @noRd
.twasChromosomeCount <- function(store) {
    n_distinct(as.character(GenomicRanges::seqnames(store$variants)))
}

#' @rdname traitPosition
#' @export
setMethod("traitPosition", "TwasWeights", function(x) {
    .getTraitPosColumn(x)
})

# Documented as `twasWeightsEntry`, not as `twasWeights`: man/ pages are
# files, and on a case-insensitive filesystem `twasWeights.Rd` and the
# constructor's `TwasWeights.Rd` are the same file, so whichever document()
# writes second silently destroys the other. Do not "simplify" this @rdname
# back to the generic's own name.
#' @title Get a Single TWAS Weights Entry
#' @description Return one \code{(study, context, trait, method)} entry of a
#'   \code{TwasWeights} collection: its gene's row, cut to that method.
#' @param x A \code{TwasWeights} object.
#' @param studyName,context,trait,method Single character identifiers. All
#'   required when the collection holds more than one entry; optional when
#'   it holds one. The method may be given in snake_case
#'   (\code{susie_inf}).
#' @return A \code{TwasWeights} holding that one entry; read it with
#'   \code{\link[=weights-methods]{weights}}, \code{\link{methodFits}} or
#'   \code{\link{cvResult}}.
#' @examples
#' twe <- TwasWeightsRow(
#'   variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)), weights = rep(0.1, 4),
#'   cvResult = list(rsq = 0.5), weightStandardized = FALSE)
#' tw <- TwasWeights(studyName = "s1", context = "brain", trait = "g1",
#'   method = "susie", entry = list(twe))
#' twasWeights(tw, studyName = "s1", context = "brain", trait = "g1",
#'   method = "susie")
#' @rdname twasWeightsEntry
#' @export
setGeneric(
    "twasWeights",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        standardGeneric("twasWeights")
    }
)

#' @rdname twasWeightsEntry
#' @export
setMethod(
    "twasWeights",
    "TwasWeights",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        .twasKeepCells(
            x,
            .twrSelectCell(x, studyName, context, trait, method)
        )
    }
)

# The cell payload that selector pins.
#
# Read directly rather than via twasWeights(): that accessor returns a
# single-cell collection, so routing the per-field accessors through it
# would dispatch straight back into them.
# @noRd
.twrSelectRowParts <- function(x, studyName, context, trait, method) {
    cell <- .twrSelectCell(x, studyName, context, trait, method)
    .twrCellParts(x, cell$row, cell$method)
}

#' @rdname resolveWeights
#' @export
setMethod(
    "resolveWeights",
    "TwasWeights",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        .twrRowResolveWeights(
            .twrSelectRowParts(x, studyName, context, trait, method)
        )
    }
)

# Adopts BiocGenerics' `weights` generic instead of a pecotmr-owned one. Its
# first formal is `object`, not `x`, and the inherited default returned NULL
# for these classes, so overriding it loses nothing. Documented as
# `weights-methods` so no \alias{weights} competes with stats::weights.
#' @title Get TWAS Weights
#' @description Extract weights from a \code{TwasWeightsRow} or from one entry
#'   of a \code{TwasWeights} collection.
#' @param object A \code{TwasWeightsRow} or \code{TwasWeights}.
#' @param ... Class-specific selection arguments.
#' @param studyName Character (length 1) or \code{NULL}. Restrict the selection
#'   to this study; \code{NULL} matches all studies.
#' @param context Character (length 1) or \code{NULL}. Restrict the selection to
#'   this context; \code{NULL} matches all contexts.
#' @param trait Character (length 1) or \code{NULL}. Restrict the selection to
#'   this trait; \code{NULL} matches all traits.
#' @param method Character (length 1) or \code{NULL}. Restrict the selection to
#'   this fine-mapping / weight method; \code{NULL} matches all methods.
#' @return A numeric vector or matrix of weights.
#' @examples
#' twe <- TwasWeightsRow(variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5),
#'   weightStandardized = FALSE)
#' weights(twe)
#' @importFrom BiocGenerics weights
#' @rdname weights-methods
#' @aliases weights-methods
#' @export
setMethod(
    "weights",
    "TwasWeights",
    function(
        object,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        ...
    ) {
        weights(
            .twrSelectRowParts(
                object,
                studyName,
                context,
                trait,
                method
            )
        )
    }
)

#' @rdname cvResult
#' @export
setMethod(
    "cvResult",
    "TwasWeights",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        parts <- .twrSelectRowParts(x, studyName, context, trait, method)
        cvResult(.fitCoverageChecked(parts, "cvResult()"))
    }
)

#' @rdname methodFits
#' @export
setMethod(
    "methodFits",
    "TwasWeights",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        parts <- .twrSelectRowParts(x, studyName, context, trait, method)
        methodFits(.fitCoverageChecked(parts, "methodFits()"))
    }
)

# A collection saved with one row per (study, context, trait, method) holds
# a gene's methods in rows of their own; they are folded into its row
# (.twasFoldMethodRows()), which also drops what saved rows carried that rows
# no longer do -- variant ids on the weights, a `fitVariants` column, each
# row's own CV partition. The inherited update runs either way, for the
# updates every collection shares.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "TwasWeights",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message("updateObject(object = 'TwasWeights')")
        }
        # A collection saved with one row per method is rebuilt per gene,
        # which stores everything as rows now hold it. Its LD references are
        # brought up to date first, as the rebuild reads them, and the
        # inherited update then runs on the rebuilt collection.
        folded <- .twasFoldMethodRows(.updateRowLdRefs(object))
        methods::callNextMethod(.twasShareFits(folded), ..., verbose = verbose)
    }
)

# A collection saved one row per gene before rows shared their fits holds
# each method's fit in its `fits` cell and mr.mash's per-fold fits in its
# cvResult. Both move to the fits dictionary, where the rows of a joint fit
# share one entry.
# @noRd
#' @importFrom methods slot<-
.twasShareFits <- function(object) {
    md <- mcols(object, use.names = FALSE)
    if (!is_in("methods", colnames(md)) || is_in("foldFits", colnames(md))) {
        return(object)
    }
    n <- lengths(md$methods)
    rowOf <- rep(seq_along(n), n)
    fitted <- .twasStoredFits(
        list_flatten(map(as.list(md$fits), unname)),
        list_flatten(map(as.list(md$cvResult), unname))
    )
    md$fits <- .twasKeysByRow(fitted$fits, rowOf)
    md$foldFits <- .twasKeysByRow(fitted$foldFits, rowOf)
    md$cvResult <- S4Vectors::SimpleList(map2(
        unname(split(fitted$cv, rowOf)),
        as.list(md$methods),
        set_names
    ))
    cols <- colnames(md)
    order <- append(setdiff(cols, "foldFits"), "foldFits", match("fits", cols))
    methods::slot(object, "fits", check = FALSE) <- fitted$entries
    S4Vectors::`mcols<-`(object, value = md[, order])
}

# @noRd
.twasKeysByRow <- function(keys, rowOf) {
    IRanges::CharacterList(unname(split(keys, rowOf)))
}

#' @rdname weightStandardized
#' @export
setMethod(
    "weightStandardized",
    "TwasWeights",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        isTRUE(
            weightStandardized(
                .twrSelectRowParts(x, studyName, context, trait, method)
            )
        )
    }
)

#' @rdname dataType
#' @export
setMethod(
    "dataType",
    "TwasWeights",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        dataType(.twrSelectRowParts(x, studyName, context, trait, method))
    }
)

# The assembled weights, each row taking the molecular type of the data it was
# fitted on (.rowDataTypes()). Applied to the assembled result, as a
# fine-mapping result's is (.fmRecordProvenance()), since joint rows are built
# on their own path.
# @noRd
.twasRecordDataTypes <- function(x, data) {
    if (is.null(x)) {
        return(x)
    }
    labels <- .rowDataTypes(x, data)
    if (all(is.na(labels))) {
        return(x)
    }
    `$<-`(
        x,
        "dataType",
        S4Vectors::SimpleList(map2(
            as.list(x$dataType),
            labels,
            .twasRowDataType
        ))
    )
}

# One row's per-method labels under one label, as they were if it has none.
# @noRd
.twasRowDataType <- function(row, label) {
    if (is.na(label)) {
        return(row)
    }
    set_names(rep(list(label), length(row)), names(row))
}

#' @rdname variantIds
#' @export
setMethod(
    "variantIds",
    "TwasWeights",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        .twrPartsVariantIds(
            .twrSelectRowParts(x, studyName, context, trait, method)
        )
    }
)

#' @rdname studyName
#' @export
setMethod("studyName", "TwasWeights", function(x) unique(as.character(x$study)))

#' @rdname ldReference
#' @export
setMethod("ldReference", "TwasWeights", function(x) {
    .ldRefsInUse(x)
})

#' @rdname contexts
#' @export
setMethod("contexts", "TwasWeights", function(x) {
    unique(as.character(x$context))
})

#' @rdname traitNames
#' @export
setMethod("traitNames", "TwasWeights", function(x) {
    unique(as.character(x$trait))
})

#' @rdname methodNames
#' @export
setMethod("methodNames", "TwasWeights", function(x) {
    unique(.twasCells(x)$method)
})

#' @rdname show-methods
#' @export
setMethod("show", "TwasWeights", function(object) {
    cells <- .twasCells(object)
    cat(glue(
        "TwasWeights: {nrow(object)} genes, {nrow(cells)} entries\n",
        .trim = FALSE
    ))
    if (nrow(object) > 0L) {
        cat(glue(
            "  {n_distinct(cells$study)} studies, ",
            "{n_distinct(cells$context)} contexts, ",
            "{n_distinct(cells$trait)} traits, ",
            "{n_distinct(cells$method)} methods\n",
            .trim = FALSE
        ))
    }
    ldSrc <- .ldRefsLabel(
        ldReference(object),
        "NULL (individual-level fit)"
    )
    cat(glue("  LD reference: {ldSrc}\n", .trim = FALSE))
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
        args = list(max.iter = 100)
    ),
    enet = list(fn = "enet_weights", impl = "enetWeights", args = list()),
    ridge = list(fn = "ridge_weights", impl = "ridgeWeights", args = list()),
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
        impl = "bayesLWeights",
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
    sparsevb = list(
        fn = "sparsevb_weights",
        impl = "sparsevbWeights",
        args = list()
    ),
    ser = list(
        fn = "ser_weights",
        impl = "serWeights",
        # No fitting defaults: serWeights extracts from a supplied fit.
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

# Which wrapper and engine each TWAS weight method actually reaches. Split
# by data path because the two share almost nothing: a method's individual
# and summary-statistic routes usually run different engines entirely, and
# only the token spelling is common.
# @noRd
.twasMethodChains <- function() {
    list(
        individual = .twasMethodChainsIndividual(),
        sumstat = .twasMethodChainsSumstat()
    )
}

# Individual-level data: X and y are in hand, so these are ordinary fitters.
# @noRd
.twasMethodChainsIndividual <- function() {
    list(
        lasso = c("lassoWeights", "glmnetWeights", "glmnet::cv.glmnet"),
        enet = c("enetWeights", "glmnetWeights", "glmnet::cv.glmnet"),
        ridge = c("ridgeWeights", "glmnetWeights", "glmnet::cv.glmnet"),
        scad = c("scadWeights", "ncvregWeights", "ncvreg::cv.ncvreg"),
        mcp = c("mcpWeights", "ncvregWeights", "ncvreg::cv.ncvreg"),
        l0learn = c("l0learnWeights", "L0Learn::L0Learn.cvfit"),
        sparsevb = c("sparsevbWeights", "sparsevb::svb.fit"),
        mrash = c("mrashWeights", "susieR::mr.ash"),
        mrmash = c(
            "mrmashWeights",
            "mrmashWrapper",
            "MrmashPriorParam",
            "buildMrmashPriorMatrices",
            "mr.mashr::mr.mash"
        ),
        bayesA = .alphabetChains("bayesA", "individual"),
        bayesB = .alphabetChains("bayesB", "individual"),
        bayesC = .alphabetChains("bayesC", "individual"),
        bayesL = .alphabetChains("bayesL", "individual"),
        bayesN = .alphabetChains("bayesN", "individual"),
        bayesR = .alphabetChains("bayesR", "individual"),
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
    )
}

# Summary statistics: each chain ends in pecotmr's own RSS solver or an
# upstream *.rss entry point, never the individual-level engine. pecotmr's
# solvers are reached through their LD-source paths (.penalizedRssPath(),
# .prsCsRun(), .sdprRun()), which take the options of penalizedRss() --
# lassosumRss()'s are a subset -- prsCs() and sdpr() under the same names, so
# a genotype reference can be used through its factor without forming the LD
# matrix. Every penalized method's path is .penalizedRssPath(): the penalty
# reaches it as one record its wrapper built, so only the path's own options
# (lambda, thr, maxiter) count as the method's.
# @noRd
.twasMethodChainsSumstat <- function() {
    list(
        scad = .penalizedRssChain("scadRssWeights"),
        mcp = .penalizedRssChain("mcpRssWeights"),
        l0learn = c(
            "l0learnRssWeights",
            ".rssShrinkGridWeights",
            ".penalizedRssPath"
        ),
        lasso = list(
            pecotmr = c(
                "lassosumRssWeights",
                ".rssShrinkGridWeights",
                ".penalizedRssPath"
            ),
            bigsnpr = c(
                "lassosumRssWeights",
                ".lassosum2Weights",
                "bigsnpr::snp_lassosum2"
            )
        ),
        enet = .penalizedRssChain("enetRssWeights"),
        ridge = .penalizedRssChain("ridgeRssWeights"),
        prsCs = c("prsCsWeights", ".prsCsRun"),
        # LDpred2's settings are all its wrappers' own arguments.
        ldpred2Auto = "ldpred2AutoWeights",
        ldpred2Grid = "ldpred2GridWeights",
        # A run across studies reads sdpr()'s names (.sdprCoupledSettings()),
        # its settings across studies, `sharing` and `rho`, among them.
        dprGibbs = c("sdprWeights", ".sdprRun", "sdpr"),
        bayesA = .alphabetChains("bayesA", "sumstat"),
        bayesB = .alphabetChains("bayesB", "sumstat"),
        bayesC = .alphabetChains("bayesC", "sumstat"),
        bayesL = .alphabetChains("bayesL", "sumstat"),
        bayesN = .alphabetChains("bayesN", "sumstat"),
        bayesR = .alphabetChains("bayesR", "sumstat"),
        mrash = c("mrashRssWeights", "susieR::mr.ash.rss"),
        mrmash = c(
            "mrmashRssWeights",
            "buildMrmashPriorMatrices",
            "mr.mashr::mr.mash.rss"
        )
    )
}

# A penalized RSS method's chain: its wrapper, the shared penalized dispatcher,
# the shrinkage grid and the path every penalty solves.
# @noRd
.penalizedRssChain <- function(wrapper) {
    c(
        wrapper,
        ".penalizedRssWeights",
        ".rssShrinkGridWeights",
        ".penalizedRssPath"
    )
}

# A Bayesian alphabet method's chains on one data path, one per backend:
# the method's wrapper, the shared dispatcher, the engine's adapter and the
# engine.
# @noRd
.alphabetChains <- function(method, path) {
    backends <- names(.alphabetBackends[[path]][[method]])
    set_names(
        map(backends, .alphabetChain, method = method, path = path),
        backends
    )
}

# @noRd
.alphabetChain <- function(backend, method, path) {
    c(
        .twasMethodCapabilities[[method]][[str_c(path, "Impl")]],
        .alphabetDispatchers[[path]],
        .alphabetFitters[[path]][[backend]],
        .alphabetEngineCallees[[path]][[backend]]
    )
}

# The function every alphabet method's wrapper hands its fit to, per path.
# @noRd
.alphabetDispatchers <- c(
    individual = "bayesAlphabetWeights",
    sumstat = "bayesAlphabetRssWeights"
)

# The engine each backend's adapter reaches, with any pecotmr hop between,
# per data path.
# @noRd
.alphabetEngineCallees <- list(
    individual = list(
        qgg = "qgg::gbayes",
        bglr = c(
            "bglrWeights",
            "BGLR::BGLR",
            "BGLR::Multitrait",
            ".bglrTermSettings"
        ),
        hibayes = "hibayes::ibrm"
    ),
    sumstat = list(
        qgg = "qgg::gbayes",
        hibayes = "hibayes::sbrm"
    )
)

# The chain for one token and input class, or NULL when the token has no path
# for that class. A method with several backends has a chain per backend,
# and this is its default backend's.
# @noRd
.twasMethodChainFor <- function(token, inputKind) {
    chains <- .twasMethodChainsOf(token, inputKind)
    if (!is.list(chains)) {
        return(chains)
    }
    chains[[1L]]
}

# Every chain a token has on one input class: a vector for a method with one
# implementation, a list keyed by backend for one with several.
# @noRd
.twasMethodChainsOf <- function(token, inputKind) {
    chains <- .twasMethodChains()
    side <- if (identical(inputKind, "QtlDataset")) "individual" else "sumstat"
    chains[[side]][[token]]
}

# Names twasWeightsPipeline() sets for each method from its own arguments.
# Refused as method options even where a `...` hop leaves the rest of a
# method's bundle unchecked (see .methodsEntryNameProblem).
# @noRd
.twasPipelineOwnedArgs <- function() {
    "backend"
}

# Argument names that are data or plumbing rather than settable options, so
# naming one is not a user choice the check should endorse.
# @noRd
.twasChainDataArgs <- function() {
    c(
        .twasPipelineOwnedArgs(),
        "X",
        "Y",
        "y",
        "stat",
        "ldReference",
        "LD",
        "ld",
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
        "spec",
        "label",
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
        try_fetch(
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
# its chain, or along every backend's chain for a method with several, whose
# own constructors check a record exactly. NULL when a hop takes `...`, since
# anything then reaches it.
# @noRd
.twasChainAccepted <- function(token, inputKind) {
    chain <- unique(unlist(.twasMethodChainsOf(token, inputKind)))
    if (is.null(chain)) {
        return(NULL)
    }
    perHop <- map(chain, .twasChainHopFormals)
    if (any(map_lgl(perHop, .twasChainHopUnknown))) {
        return(NULL)
    }
    setdiff(unique(list_c(perHop)), .twasChainDataArgs())
}

# A hop contributes nothing checkable when it could not be read, or when it
# takes `...` and so accepts any name.
# @noRd
.twasChainHopUnknown <- function(fm) {
    (length(fm) == 1L && is.na(fm[[1L]])) || is_in("...", fm)
}

# The implementation that runs one method on one input path.
# @noRd
.twasMethodImplFor <- function(token, inputKind) {
    info <- .twasMethodCapabilities[[token]]
    if (is.null(info)) {
        return(NULL)
    }
    if (identical(inputKind, "QtlDataset")) {
        info$individualImpl
    } else {
        info$sumstatImpl
    }
}

# The Options constructor for one implementation, read off that function's
# own `methodArgs` default rather than transcribed into a table here.
#
# Derived because a second table would be a second source of truth: the
# individual and summary-statistics paths of one method often reach
# DIFFERENT packages -- lasso is glmnet on individual data and lassosum on
# summary statistics, dprGibbs is RcppDPR then SDPR -- so a hand-written
# registry has to restate a pairing the implementation already declares.
# @noRd
#
# An implementation that runs on several engines (the Bayesian alphabet)
# defaults `methodArgs` to NULL and has a constructor per backend, from
# .alphabetBackends: this answers `backend`'s, or the default's when NULL.
# @noRd
.twasImplCtorName <- function(impl, backend = NULL) {
    nms <- .twasImplCtorNames(impl)
    if (length(nms) <= 1L) {
        return(nms)
    }
    unname(nms[[backend %||% names(nms)[[1L]]]])
}

# Every Options constructor an implementation takes: the one its `methodArgs`
# default names, or one per backend, keyed by backend. NULL for an
# implementation that takes none.
# @noRd
.twasImplCtorNames <- function(impl) {
    if (is.null(impl)) {
        return(NULL)
    }
    fn <- try_fetch(
        get(impl, envir = asNamespace("pecotmr")),
        error = function(cnd) NULL
    )
    if (is.null(fn)) {
        return(NULL)
    }
    default <- formals(fn)$methodArgs
    if (is.call(default)) {
        return(as.character(default[[1L]]))
    }
    .twasImplBackendCtors(impl)
}

# The constructors of an implementation with several backends, keyed by
# backend: found through the method and data path whose implementation it
# is. NULL for any other implementation.
# @noRd
.twasImplBackendCtors <- function(impl) {
    found <- compact(map(
        c("individual", "sumstat"),
        .twasPathBackendCtors,
        impl = impl
    ))
    if (length(found) == 0L) NULL else found[[1L]]
}

# The backend constructors of the method whose implementation on `path` is
# `impl`, or NULL.
# @noRd
.twasPathBackendCtors <- function(path, impl) {
    caps <- .twasMethodCapabilities
    hit <- names(caps)[map_lgl(
        caps,
        .twasCapUsesImpl,
        field = str_c(path, "Impl"),
        impl = impl
    )]
    if (length(hit) == 0L) {
        return(NULL)
    }
    ctors <- .methodBackends[[path]][[hit[[1L]]]]
    # A method whose engines take no Options constructor (LDpred2).
    if (all(is.na(ctors))) NULL else ctors
}

# @noRd
.twasCapUsesImpl <- function(cap, field, impl) {
    identical(cap[[field]], impl)
}

# The constructor itself. Deliberately separate from the name: the
# input-type index below needs the NAME and must not call the constructor,
# which would re-enter .newMethodOptions() and recurse forever.
# @noRd
.twasImplCtor <- function(impl) {
    nm <- .twasImplCtorName(impl)
    if (is.null(nm)) {
        return(NULL)
    }
    try_fetch(
        get(nm, envir = asNamespace("pecotmr")),
        error = function(cnd) NULL
    )
}

# One method's Options constructor on one input path -- its default
# backend's, for a method with several; NULL when the method does not run on
# that path at all.
# @noRd
.twasMethodCtorFor <- function(token, inputKind) {
    .twasImplCtor(.twasMethodImplFor(token, inputKind))
}

# Every Options constructor one method takes on one input path, by name.
# @noRd
.twasMethodCtorNamesFor <- function(token, inputKind) {
    .twasImplCtorNames(.twasMethodImplFor(token, inputKind))
}

# Every TWAS weight method that runs on one input path: the configurable
# ones with an implementation there, plus the fit-derived ones, whose
# availability is not ours to state -- a fit can only be extracted from
# where fineMappingPipeline() could produce it, so it comes from the
# fine-mapping registry.
# @noRd
.twasMethodsFor <- function(inputKind) {
    conf <- keep(
        .twasConfigurableMethods(),
        .twasMethodRunsOn,
        inputKind = inputKind
    )
    unique(c(
        conf,
        intersect(.twasFineMappingTokens(), .fmMethodsFor(inputKind))
    ))
}

# Whether configurable method `token` has an implementation on `inputKind`.
# @noRd
.twasMethodRunsOn <- function(token, inputKind) {
    !is.null(.twasMethodImplFor(token, inputKind))
}

# Every configurable method's name, independent of input path. The answer to
# "is this a method that takes options", which the token normalisers ask.
# @noRd
.twasConfigurableMethods <- function() {
    names(.twasMethodCapabilities)
}

# Why a token takes no options, so the error can say something useful rather
# than just "unknown method".
# @noRd
.twasTokenNoArgsReason <- function(token) {
    # Derived, not transcribed: this list was hand-written and had already
    # drifted from .twasFineMappingTokens(), which is what the character
    # path uses -- `ser` was in that one and missing here, so asking the
    # aggregator about it answered "unknown method" instead of saying there
    # is nothing to configure.
    if (is_in(token, .twasFineMappingTokens())) {
        return(glue(
            "'{token}' extracts weights from a fit supplied via ",
            "`fineMappingResult`, so there is nothing to configure here. ",
            "Configure the fit in fineMappingPipeline() instead."
        ))
    }
    NULL
}

# `backend` given as a method option. The pipeline owns it, keyed by method
# (.twasPipelineOwnedArgs), so the refusal says where it goes. The example
# uses the caller's own choice when it is one of the method's backends on
# this input.
# @noRd
.twasBackendArgHint <- function(token, entry, inputKind) {
    if (!is_in("backend", names(entry))) {
        return("")
    }
    side <- if (identical(inputKind, "QtlDataset")) "individual" else "sumstat"
    choices <- names(.methodBackends[[side]][[token]])
    if (!.methodHasBackends(token, side)) {
        return(glue(
            " `backend` is twasWeightsPipeline()'s own argument, and {token} ",
            "has one implementation on {inputKind}, so there is no backend ",
            "to choose."
        ))
    }
    given <- entry$backend
    example <- if (length(given) == 1L && is_in(given, choices)) {
        given
    } else {
        choices[[1L]]
    }
    glue(
        " Choose it with twasWeightsPipeline()'s own `backend` argument, ",
        "keyed by method: backend = c({token} = \"{example}\"). ",
        "{token} on {inputKind} offers {str_flatten(choices, ', ')}."
    )
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
# token itself (`susieInf`), its camelCase weight function (`susieInfWeights`),
# its snake_case registry key (`susie_inf_weights`) and that key's row label
# (`susie_inf`). Stripping the suffix is not enough on its own: that turns
# `susie_inf_weights` into `susie_inf`, which is not the registry's token.
# @noRd
.twasFmTokenFor <- function(method) {
    methodKeys <- map_chr(.twasFineMappingMethodAdapters, "methodKey")
    rowLabels <- str_remove(methodKeys, "_weights$")
    hit <- names(methodKeys)[which(
        methodKeys == method | rowLabels == method
    )]
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

# The per-fold weights supplied for `token`, under any of its spellings. The
# gate and .twasFoldWeightsW() must agree on this, or the gate can refuse
# weights the fold lookup would have found.
# @noRd
.twasFoldWeightsFor <- function(foldWeights, token) {
    if (is.null(foldWeights)) {
        return(NULL)
    }
    adapter <- .twasFineMappingMethodAdapters[[token]]
    spellings <- c(token, adapter$methodKey, adapter$weightFn)
    hit <- spellings[is_in(spellings, names(foldWeights))]
    if (length(hit) == 0L) {
        return(NULL)
    }
    foldWeights[[hit[[1]]]]
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
# `available(token, ...)` answers "is this token's fit present?" for the
# caller's own supply channel; `...` carries what that channel needs.
# @noRd
.twasRequireSusieFits <- function(weightMethods, available, fnLabel, how, ...) {
    tokens <- .twasSusieTokensRequested(weightMethods)
    if (length(tokens) == 0L) {
        return(invisible(NULL))
    }
    missing <- tokens[!map_lgl(tokens, available, ...)]
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

# One method's argument list with `arg` set, leaving every other method and
# every other argument of that method untouched.
# @noRd
.twasSetMethodArg <- function(weightMethods, method, arg, value) {
    argValue <- set_names(list(value), arg)
    updated <- list_assign(weightMethods[[method]], !!!argValue)
    methodEntry <- set_names(list(updated), method)
    list_assign(weightMethods, !!!methodEntry)
}

# Resolve the supplied fine-mapping fits onto the method args, for every
# token with an adapter. This never fine-maps: the SuSiE-family weight methods
# extract from a fit that fineMappingPipeline() produced, and a missing one is
# an error, not a cue to fit here.
# @noRd
.prepareSusieWeightMethods <- function(weightMethods, fittedModels = NULL) {
    reduce(
        names(.twasFineMappingMethodAdapters),
        .twasResolveTokenFit,
        fittedModels = fittedModels %||% list(),
        .init = weightMethods
    )
}

# One token's fit on its method's args, tagged with its fine-mapping class: a
# fit the caller put there wins, else the one in `fittedModels`. It is written
# under the spelling the caller used, where the extractor will look.
# @noRd
.twasResolveTokenFit <- function(weightMethods, token, fittedModels) {
    nm <- .twasMethodNameFor(weightMethods, token)
    if (is.na(nm)) {
        return(weightMethods)
    }
    fitArg <- .twasFineMappingMethodAdapters[[token]]$fitArg
    fit <- weightMethods[[nm]][[fitArg]] %||% fittedModels[[token]]
    if (is.null(fit)) {
        return(weightMethods)
    }
    .twasSetMethodArg(
        weightMethods,
        nm,
        fitArg,
        .setFinemappingFitClass(fit, token)
    )
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
.twasFoldPriors <- function(args, method, j, dataDrivenPriorMatricesCv) {
    priorUpdates <- compact(list(
        dataDrivenPriorMatrices = if (
            !is.null(dataDrivenPriorMatricesCv) &&
                is_in(method, c("mrmash_weights", "mrmashWeights"))
        ) {
            dataDrivenPriorMatricesCv[[j]]
        }
    ))
    list_assign(args, !!!priorUpdates)
}

# A SuSiE-family method's weights for fold `j`, from fineMappingPipeline's CV
# (`foldWeights`: method -> variants x folds, rows named by variant, or
# variants x folds x outcomes): those methods extract from a fit and never
# fine-map, so a fold can only be scored with the weights its own fit gave.
# Aligned to the fold's training variants (zero where the fit had none), one
# column per outcome, or NULL when the method is not a fine-mapping one or
# has no weights for the fold.
# @noRd
.twasFoldWeightsW <- function(method, j, Xtr, Ytr, foldWeights) {
    mk <- .twasFmTokenFor(method)
    fw <- if (is.na(mk)) NULL else .twasFoldWeightsFor(foldWeights, mk)
    key <- str_c("fold_", j)
    if (is.null(fw) || !is_in(key, colnames(fw))) {
        return(NULL)
    }
    if (is.null(rownames(fw))) {
        abort(glue(
            "twasWeightsCv: the per-fold weights for '{mk}' do not name their ",
            "variants, so they cannot be aligned to X. Take them through ",
            "cvResult(), which names them."
        ))
    }
    perOutcome <- .twasFoldOutcomeWeights(fw, key, Ytr, mk)
    hit <- match(rownames(fw), colnames(Xtr))
    W <- matrix(0, ncol(Xtr), ncol(Ytr), dimnames = list(colnames(Xtr), NULL))
    W[hit[!is.na(hit)], ] <- perOutcome[!is.na(hit), , drop = FALSE]
    W
}

# Fold `key`'s weights as variants x outcomes. A variants x folds matrix
# holds one weight per variant, which every outcome shares, as it does for a
# univariate fit; a variants x folds x outcomes array, as mvSuSiE gives, holds
# each outcome's own, matched to the columns of Y by name.
# @noRd
.twasFoldOutcomeWeights <- function(fw, key, Ytr, mk) {
    if (length(dim(fw)) != 3L) {
        return(matrix(fw[, key], nrow(fw), ncol(Ytr)))
    }
    outcomes <- colnames(Ytr)
    missing <- setdiff(outcomes, dimnames(fw)[[3L]])
    if (length(missing) > 0L) {
        missingStr <- str_flatten(missing, ", ")
        abort(glue(
            "twasWeightsCv: the per-fold weights for '{mk}' have no outcome ",
            "{missingStr}. Name the third axis by the columns of Y."
        ))
    }
    matrix(fw[, key, outcomes], nrow(fw), length(outcomes))
}

# One fold's multivariate weight fit; returns list(W, fit).
# @noRd
.twasFoldMultivariate <- function(method, fnName, args, Xtr, Ytr, j, ctx) {
    withPriors <- .twasFoldPriors(
        args,
        method,
        j,
        ctx$dataDrivenPriorMatricesCv
    )
    retaining <- !identical(ctx$fitRetention, "none")
    retentionArg <- compact(list(
        fitRetention = if (retaining) ctx$fitRetention
    ))
    fitArgs <- list_assign(withPriors, !!!retentionArg)
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
    mk <- str_remove(method, "_weights$|Weights$")
    given <- .twasFoldWeightsW(method, j, Xtr, Ytr, ctx$foldWeights)
    if (!is.null(given)) {
        return(list(mk = mk, W = given, fit = NULL))
    }
    fnName <- .resolveMethodFunction(method, args)
    wholeY <- is_in(method, ctx$multivariateWeightMethods) ||
        .twasFitsWholeY(fnName)
    fit <- if (wholeY) {
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
#' @param BPPARAM A \code{BiocParallelParam} the folds run on;
#'   \code{SerialParam()} (default) runs them one at a time. Each draws from
#'   its own random-number stream, seeded by \code{seed}, else by the param's
#'   own \code{RNGseed}, else by one draw from the session RNG, so a run
#'   gives the same numbers however many workers it has.
#' @param verbose Integer controlling verbosity level: 0 = suppress all
#'   messages, 1 = suppress external package messages (default), 2 = show all
#'   messages including those from external packages.
#' @param fitRetention How much of each per-method fit is kept on the result:
#'   \code{"none"} (default) keeps none, \code{"slim"} keeps the trimmed
#'   payload downstream needs, \code{"full"} keeps the whole fitted object.
#'   Every method keeps a fit, which names the variants it was fitted on.
#'   \code{"slim"} and \code{"full"} differ for mr.mash, glmnet
#'   (\code{lasso}, \code{enet}, \code{ridge}), ncvreg (\code{scad},
#'   \code{mcp}), \code{l0learn}, the Bayesian alphabet and \code{mrash};
#'   for the others the two behave alike.
#' @param seed Integer or \code{NULL}. When supplied, seeds both the
#'   main-process RNG (fold partitioning, variant sub-sampling) and the
#'   streams the folds draw from (see \code{BPPARAM}). The main-process seed
#'   is scoped to the call, so the session RNG is left as it was found.
#'   \code{NULL} (default) seeds neither.
#' @param foldWeights Optional per-fold weights of fine-mapping methods, as
#'   \code{method -> variants x folds} matrix (rows named by variant, columns
#'   \code{fold_<j>}): the \code{foldWeights} of the \code{cvResult} a
#'   \code{fineMappingPipeline()} run with \code{cvFolds > 1} records. Every
#'   column of \code{Y} shares those weights; a method with weights of its own
#'   for each outcome, such as mvSuSiE, takes a variants x folds x outcomes
#'   array instead, its third axis named by the columns of \code{Y}.
#'   SuSiE-family weight wrappers extract from a fit and never fine-map, so a
#'   fold's own weights are what make cross-validating them possible; they
#'   must have been produced on \code{samplePartitions}, since weights
#'   trained on different folds would leak held-out samples.
#' @param dataDrivenPriorMatricesCv Optional list, one element per fold, of
#'   data-driven prior matrices for the mr.mash learner.
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
#' @importFrom BiocParallel SerialParam
#' @importFrom quadprog solve.QP
#' @examples
#' data(multiTraitData)
#' X <- multiTraitData$X[, 1:80]
#' Y <- multiTraitData$Y
#' # A cross-validated method is refit on each fold's training rows. The
#' # SuSiE family never fine-maps, so it needs each fold's own weights,
#' # passed as `foldWeights` from a fineMappingPipeline() run with
#' # cvFolds > 1.
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
    BPPARAM = SerialParam(),
    verbose = 1,
    fitRetention = c("none", "slim", "full"),
    seed = NULL,
    dataDrivenPriorMatricesCv = NULL,
    foldWeights = NULL
) {
    # X / Y / fold are asserted downstream in .cvPrepareData; these are the
    # arguments nothing else checks.
    assertDataFrame(samplePartitions, null.ok = TRUE)
    # NOT assertCount: `Inf` is the "no cap" sentinel (jointEngine passes it
    # when cfg$crossValidation$maxVariants is unset), and
    # .cvSubsampleVariants relies on `ncol(X) <= maxNumVariants` being FALSE
    # for it.
    assertNumber(maxNumVariants, lower = 1, null.ok = TRUE)
    .assertBpparam(BPPARAM)
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
        BPPARAM = BPPARAM,
        verbose = verbose,
        fitRetention = fitRetention,
        seed = seed,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        foldWeights = foldWeights
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
# fingerprint (see .cvPartitionKey); a fit carrying none is refused too, since
# it cannot be shown to match.
# @noRd
.twasCheckFoldWeightsPartition <- function(foldWeights, samplePartitions) {
    if (is.null(foldWeights) || length(foldWeights) == 0L) {
        return(invisible(NULL))
    }
    if (is.null(samplePartitions)) {
        msg <- glue(
            "twasWeightsCv: `foldWeights` needs the fold partition those ",
            "weights were trained on. Pass the fine-mapping CV's ",
            "`samplePartition` as `samplePartitions`; a freshly drawn ",
            "partition would score each fold with weights that saw its ",
            "held-out samples."
        )
        abort(msg)
    }
    want <- .cvPartitionKey(samplePartitions)
    for (m in names(foldWeights)) {
        got <- attr(foldWeights[[m]], "partitionKey")
        if (is.null(got)) {
            msg <- glue(
                "twasWeightsCv: the per-fold weights for '{m}' carry no ",
                "partition fingerprint, so they cannot be shown to match ",
                "`samplePartitions`. Take them from a fineMappingPipeline() ",
                "run with cvFolds > 1."
            )
            abort(msg)
        }
        if (!identical(got, want)) {
            msg <- glue(
                "twasWeightsCv: the per-fold weights for '{m}' were trained ",
                "on a different fold partition than the one being scored. ",
                "Use the fine-mapping CV's own `samplePartition`."
            )
            abort(msg)
        }
    }
    invisible(NULL)
}

# What .weightFitFold() needs to refit one fold, or NULL on the
# partition-only path where .cvNoopFitFold() reads nothing.
# @noRd
.twasCvFoldCtx <- function(
    weightMethods,
    dataDrivenPriorMatricesCv,
    foldWeights,
    fitRetention,
    verbose
) {
    if (is.null(weightMethods)) {
        return(NULL)
    }
    list(
        weightMethods = weightMethods,
        multivariateWeightMethods = .twasCvMultivariateMethods,
        dataDrivenPriorMatricesCv = dataDrivenPriorMatricesCv,
        foldWeights = foldWeights,
        fitRetention = fitRetention,
        verbose = verbose
    )
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
    BPPARAM,
    verbose,
    fitRetention,
    seed,
    dataDrivenPriorMatricesCv,
    foldWeights
) {
    weightMethods <- .twasCvResolveMethods(
        weightMethods,
        foldWeights,
        samplePartitions,
        seed,
        verbose
    )
    # No methods means "partition only": the no-op fold fit walks the same
    # engine so the returned partition is the one a real run would use.
    noMethods <- is.null(weightMethods)
    res <- .crossValidateWeights(
        X,
        Y,
        fold = fold,
        samplePartitions = samplePartitions,
        fitFold = if (noMethods) .cvNoopFitFold else .weightFitFold,
        fitFoldCtx = .twasCvFoldCtx(
            weightMethods,
            dataDrivenPriorMatricesCv,
            foldWeights,
            fitRetention,
            verbose
        ),
        BPPARAM = BPPARAM,
        maxNumVariants = maxNumVariants,
        variantsToKeep = variantsToKeep,
        fitRetention = fitRetention,
        verbose = verbose,
        seed = seed
    )
    if (noMethods) list(samplePartition = res$samplePartition) else res
}

# Check the fold-fit inputs and turn a character `weightMethods` into the
# method record the fold fitter takes. Also the one place that notices an
# unseeded call, since this is where the partition is about to be drawn.
# @noRd
.twasCvResolveMethods <- function(
    weightMethods,
    foldWeights,
    samplePartitions,
    seed,
    verbose
) {
    .twasCheckFoldWeightsPartition(foldWeights, samplePartitions)
    .twasRequireSusieFits(
        weightMethods,
        available = .twasFoldWeightsAvailable,
        fnLabel = "twasWeightsCv",
        how = str_c(
            "Cross-validation refits on each fold, so it needs the weights ",
            "that fold's own fit gave: pass `foldWeights = list(<token> = ",
            "<fold weights>)` from a fineMappingPipeline() run with ",
            "cvFolds > 1."
        ),
        foldWeights = foldWeights
    )
    if (is.null(seed) && !exists(".Random.seed") && verbose >= 1) {
        inform(str_c(
            "! No seed set. Pass `seed=` or call ",
            "set.seed() for reproducibility."
        ))
    }
    .twasAsMethodRecord(weightMethods)
}

# Whether `foldWeights` carries per-fold weights for `tk`.
# @noRd
.twasFoldWeightsAvailable <- function(tk, foldWeights) {
    !is.null(.twasFoldWeightsFor(foldWeights, tk))
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
        ctx$fitRetention,
        ctx$retainX,
        fnName
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

# Pass the retention level down to the target weight function when the caller
# is keeping fits at all; every weight method takes `fitRetention`. `retainX`
# goes only where it can apply: a method's own level may be lower than the
# run's, and most methods' fits hold no copy of X.
# @noRd
.twasApplyRetainFit <- function(args, fitRetention, retainX, fnName) {
    if (identical(fitRetention, "none")) {
        return(args)
    }
    level <- args$fitRetention %||% fitRetention
    keepingX <- retainX && identical(level, "full") && .twasKeepsX(fnName)
    list_assign(
        args,
        fitRetention = level,
        !!!compact(list(retainX = if (keepingX) TRUE))
    )
}

# Whether a weight function's fit can hold a copy of X, as its own formal says.
# @noRd
.twasKeepsX <- function(fnName) {
    is_in("retainX", names(formals(get0(fnName, mode = "function"))))
}

# `retainX = TRUE` with no method whose fit has a copy of X would keep nothing.
# @noRd
.twasAssertKeepsX <- function(weightMethods) {
    fns <- map2_chr(names(weightMethods), weightMethods, .resolveMethodFunction)
    if (!any(map_lgl(fns, .twasKeepsX))) {
        abort(glue(
            "learnTwasWeights: `retainX = TRUE` needs a method whose fit has ",
            "a copy of X; the fits of {str_flatten(fns, ', ')} hold none."
        ))
    }
    invisible(NULL)
}

# Dispatch weight fitting to the multivariate or per-column univariate path;
# returns list(weights, methodFit).
# @noRd
.twasFitWeightsMatrix <- function(fnName, args, ctx, methodName) {
    if (
        is_in(methodName, .twasMultivariateWeightMethods) ||
            .twasFitsWholeY(fnName)
    ) {
        .twasFitMultivariate(fnName, args, ctx)
    } else {
        .twasFitUnivariate(fnName, args, ctx)
    }
}

# Whether weight function `fnName` fits every column of Y at once, as it does
# when it takes Y: a multivariate method, or a univariate one fitted across a
# group's studies (.lassoCoupledWeights()), which a method's arguments name
# as their implementation (.resolveMethodFunction()).
# @noRd
.twasFitsWholeY <- function(fnName) {
    is_in("Y", names(formals(match.fun(fnName))))
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
# this column's own slice is taken. An outcome missing on some samples, as a
# cross-study group's condition is on the other studies' samples, is fitted
# on the samples observing it (.twasObservedColumn()).
# @noRd
.twasUnivariateColumn <- function(k, fnName, args, ctx) {
    if (anyNA(ctx$Y[, k])) {
        return(.twasObservedColumn(
            ctx$Xfiltered,
            ctx$Y[, k],
            k,
            fnName,
            args,
            ctx
        ))
    }
    w <- .twasColumnFit(fnName, ctx$Xfiltered, ctx$Y[, k], args, ctx)
    list(
        weights = if (is.matrix(w)) w[, k] else w,
        fit = .twasColumnFitRecord(w, ctx)
    )
}

# Outcome `y` (column `k`) fitted on the samples observing it alone, over the
# variants that vary among them; the rest -- a variant its study lacks,
# constant there -- take weight 0. The weights come back over X's columns.
# @noRd
.twasObservedColumn <- function(X, y, k, fnName, args, ctx) {
    observed <- !is.na(y)
    held <- X[observed, , drop = FALSE]
    varying <- .nonzeroVarColumns(held)
    w <- .twasColumnFit(
        fnName,
        held[, varying, drop = FALSE],
        y[observed],
        args,
        ctx
    )
    own <- if (is.matrix(w)) w[, min(k, ncol(w))] else w
    list(
        weights = replace(numeric(ncol(X)), varying, unname(own)),
        fit = .twasColumnFitRecord(w, ctx)
    )
}

# One weight function's fit of `y` on `X`, quiet below verbose level 2.
# @noRd
.twasColumnFit <- function(fnName, X, y, args, ctx) {
    call <- .twasWeightCallArgs(fnName, list(X = X, y = y), args)
    if (ctx$verbose < 2) {
        return(.quietEval(exec(fnName, !!!call)))
    }
    exec(fnName, !!!call)
}

# The fit a column's weights carry, kept unless the run keeps none.
# @noRd
.twasColumnFitRecord <- function(w, ctx) {
    if (identical(ctx$fitRetention, "none")) {
        return(NULL)
    }
    attr(w, "fit")
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
    TwasWeightsRow(
        variantIds = variantIds,
        weights = weights,
        methodFits = fits,
        cvResult = NULL,
        weightStandardized = ctx$standardized,
        weightsDataType = ctx$dataType
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
#' @param BPPARAM A \code{BiocParallelParam} the methods run on;
#'   \code{SerialParam()} (default) runs them one at a time. Each draws from
#'   its own random-number stream, seeded by \code{seed}, else by the param's
#'   own \code{RNGseed}, else by one draw from the session RNG, so a run
#'   gives the same numbers however many workers it has.
#' @param fittedModels Named list of fitted SuSiE-family models, keyed by
#'   token (\code{susie}, \code{susieInf}, \code{mvsusie}, \code{fsusie}).
#'   Required whenever a SuSiE-family weight method is requested: those
#'   methods extract weights from an existing fit and never fine-map, so a
#'   missing fit is an error. Run \code{\link{fineMappingPipeline}} to
#'   produce the fits.
#' @param fitRetention How much of each per-method fit is kept on the result:
#'   \code{"none"} (default) keeps none, \code{"slim"} keeps the trimmed
#'   payload downstream needs, \code{"full"} keeps the whole fitted object but
#'   for any copy of \code{X}.
#'   Every method keeps a fit, which names the variants it was fitted on.
#'   \code{"slim"} and \code{"full"} differ for mr.mash, glmnet
#'   (\code{lasso}, \code{enet}, \code{ridge}), ncvreg (\code{scad},
#'   \code{mcp}), \code{l0learn}, the Bayesian alphabet and \code{mrash};
#'   for the others the two behave alike.
#' @param retainX Logical. Keep a fit's copy of \code{X}, for the methods whose
#'   fit has one (\code{lasso}, \code{enet}, \code{ridge}, \code{scad},
#'   \code{mcp}, \code{mrash} and the DPR methods), at
#'   \code{fitRetention = "full"} (default \code{FALSE}). An error at any
#'   other \code{fitRetention}, and when no method in \code{weightMethods} has
#'   such a copy.
#' @param verbose Integer controlling verbosity level: 0 = suppress all
#'   messages, 1 = suppress external package messages (default), 2 = show all
#'   messages including those from external packages.
#' @param studyName Character. Study identity label recorded on the resulting
#'   weights.
#' @param context Character. Context identity label recorded on the resulting
#'   weights.
#' @param trait Character. Trait identity label recorded on the resulting
#'   weights.
#' @param weightStandardized Logical. Whether the supplied \code{X} / \code{Y}
#'   are already standardized. Default \code{FALSE}.
#' @param weightsDataType Character or \code{NULL}. The molecular type recorded
#'   on the weights (e.g. \code{"expression"}).
#' @param ldReference An LD reference (a genotype panel, see
#'   \code{\link{readGenotypes}}, a \code{GenotypeMixture} or an
#'   \code{LdMatrix}) to record on the weights, or \code{NULL}.
#' @param seed Integer or \code{NULL}. When supplied, seeds the main-process
#'   RNG and the streams the methods draw from (see \code{BPPARAM}). The
#'   main-process seed is scoped to the call, so the session RNG is left as it
#'   was found. \code{NULL} (default) seeds neither.
#' @return A \code{TwasWeights} with one row, for the (\code{studyName},
#'   \code{context}, \code{trait}) gene, holding each method's weights, fit
#'   and cross-validation. Its methods are named by camelCase token: a
#'   method given as \code{susie_inf_weights} is \code{"susieInf"}.
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
    studyName = "",
    context = "",
    trait = "",
    BPPARAM = SerialParam(),
    fittedModels = NULL,
    fitRetention = c("none", "slim", "full"),
    retainX = FALSE,
    weightStandardized = FALSE,
    weightsDataType = NULL,
    ldReference = NULL,
    verbose = 1,
    seed = NULL
) {
    .twasAssertRowLabels(studyName, context, trait)
    .assertBpparam(BPPARAM)
    fitRetention <- arg_match(fitRetention)
    .assertRetainX(retainX, fitRetention)
    assertFlag(weightStandardized)
    assertString(weightsDataType, null.ok = TRUE)
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
        studyName = studyName,
        context = context,
        trait = trait,
        BPPARAM = BPPARAM,
        fittedModels = fittedModels,
        fitRetention = fitRetention,
        retainX = retainX,
        weightStandardized = weightStandardized,
        weightsDataType = weightsDataType,
        ldReference = ldReference,
        verbose = verbose,
        seed = seed
    )
}

# The labels learnTwasWeights() gives its rows.
# @noRd
.twasAssertRowLabels <- function(studyName, context, trait) {
    assertString(studyName)
    assertString(context)
    assertString(trait)
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

# Fit every weight method on BPPARAM's workers, keyed by method name.
# @noRd
.twasFitAllMethods <- function(weightMethods, ctx, BPPARAM) {
    weightsList <- .bpMap(
        names(weightMethods),
        .computeMethodWeights,
        weightMethods,
        ctx,
        BPPARAM = BPPARAM,
        seed = ctx$rngSeed
    )
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
        available = .twasLearnFitAvailable,
        fnLabel = "learnTwasWeights",
        how = str_c(
            "Pass it as `fittedModels = list(<token> = <fit>)`, or on the ",
            "method's own arguments; run fineMappingPipeline() first to ",
            "produce one."
        ),
        resolvedMethods = resolvedMethods,
        fittedModels = fittedModels
    )
}

.learnTwasWeightsImpl <- function(
    X,
    Y,
    weightMethods,
    studyName,
    context,
    trait,
    BPPARAM,
    fittedModels,
    fitRetention,
    retainX,
    weightStandardized,
    weightsDataType,
    ldReference,
    verbose,
    seed
) {
    .applySeed(seed)
    # fitRetention arrives already matched: learnTwasWeights is the only
    # caller and validates it there, where the choices are declared.
    Y <- .twasValidateXY(X, Y)
    resolvedMethods <- .twasAsMethodRecord(weightMethods)
    .twasRequireLearnFits(resolvedMethods, fittedModels)
    if (retainX) {
        .twasAssertKeepsX(resolvedMethods)
    }
    validColumns <- .nonzeroVarColumns(X)
    Xfiltered <- as.matrix(X[, validColumns, drop = FALSE])
    prepared <- .prepareSusieWeightMethods(resolvedMethods, fittedModels)
    ctx <- list(
        X = X,
        Y = Y,
        Xfiltered = Xfiltered,
        validColumns = validColumns,
        study = studyName,
        context = context,
        trait = trait,
        fitRetention = fitRetention,
        retainX = retainX,
        standardized = weightStandardized,
        dataType = weightsDataType,
        verbose = verbose,
        rngSeed = seed
    )
    weightsList <- .twasFitAllMethods(prepared, ctx, BPPARAM) |>
        .twasApplyRownames(X)
    rows <- .buildTwasWeightEntries(weightsList, .twasVariantIds(X), ctx)
    .twasLearnResult(rows, ldReference)
}

# The learn path's column-wise row record as a TwasWeights collection. NOT
# the pipeline's same-shaped helper, which takes a LIST of per-row records
# and map_chr()s across them -- a different input entirely, and sharing a
# name with it would let collate order decide which one a call reaches.
# @noRd
.twasLearnResult <- function(rows, ldReference) {
    TwasWeights(
        studyName = rows$study,
        context = rows$context,
        trait = rows$trait,
        method = rows$method,
        entry = rows$entry,
        ldReference = ldReference
    )
}

# `weightMethods` takes the method record or a character shorthand for it.
# @noRd
.twasAsMethodRecord <- function(weightMethods) {
    if (is.character(weightMethods)) {
        return(.twasMethodLookup(weightMethods))
    }
    weightMethods
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
#' twe <- TwasWeightsRow(variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
#'   weights = rep(0.1, 4), cvResult = list(rsq = 0.5),
#'   weightStandardized = FALSE)
#' tw <- TwasWeights(studyName = "s1", context = "brain", trait = "g1",
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
        # Every (gene, method) cell's weights, keyed by the method name for
        # the legacy "<method>_predicted" convention; ensembleWeights()
        # rebinds the suffix.
        cells <- .twasCells(weightsList)
        wl <- set_names(
            map2(cells$row, cells$method, .twasCellWeightsAt, x = weightsList),
            str_c(cells$method, "_weights")
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
#' data-driven prior for the inclusion probability, \code{pi}, of the
#' spike-and-slab Bayesian methods (bayesB and bayesC).
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
        cells <- .twasCells(weightResults)
        idx <- which(cells$method == "mrash")
        if (length(idx) == 0L) {
            msg <- glue(
                "mr.ash entry not found in TwasWeights. Run ",
                "learnTwasWeights() ",
                "with fitRetention = \"slim\" and ensure 'mrash' is in ",
                "the method list."
            )
            abort(msg)
        }
        fit <- methodFits(
            .twrCellParts(weightResults, cells$row[[idx[[1L]]]], "mrash")
        )
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

# Column `cn` of `object`'s mcols.
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
# own formals (a pre-fit, fitRetention, ...) bind by name; any
# else is a tool option and goes in the wrapper's `methodArgs` list, so an
# unknown option errors inside the wrapper rather than vanishing. Wrappers
# with no `methodArgs` formal get everything by name, and R reports an unused
# argument -- which is the point.
# @noRd
.twasWeightCallArgs <- function(fnName, baseArgs, userArgs) {
    c(baseArgs, .splitMethodArgs(fnName, userArgs))
}

# One univariate fold's weight column for outcome `k` (quiet unless verbose),
# fitted on its observed training samples alone when it is missing on some
# (.twasObservedColumn()).
# @noRd
.twasFitColWeight <- function(k, ctx, fnName, Xtr, Ytr, args) {
    if (anyNA(Ytr[, k])) {
        return(.twasObservedColumn(Xtr, Ytr[, k], k, fnName, args, ctx)$weights)
    }
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

# The weights of method `method` of row `row` of `x`.
# @noRd
.twasCellWeightsAt <- function(row, method, x) {
    weights(.twrCellParts(x, row, method))
}

# X %*% w for one weight vector/matrix (coerced to a 1-column matrix if needed).
# @noRd
.twasPredictOne <- function(w, X) {
    if (!is.matrix(w)) {
        w <- matrix(w, ncol = 1)
    }
    X %*% w
}
