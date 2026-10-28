# =============================================================================
# QtlFineMappingResult S4 class
# -----------------------------------------------------------------------------
# RangedTupleList collection keyed by the identity tuple (study, context,
# trait, method). Each element is one fit's variants, with their topLoci
# values as its mcols; the fit and its cross-validation are columns of the
# collection's own mcols. Class-level slots:
#   * ldReference   the LD references RSS-derived fits used, each stored
#                once; the ldReference mcols column names each row's (NA
#                for an individual-level fit).
# Optional columns jointStudies / jointContexts / jointTraits carry
# semicolon-joined member identities for cross-axis joint fits emitted by
# the joint dispatchers.
# =============================================================================

#' @include AllClasses.R tupleSelectors.R
NULL

#' @title QTL Fine-Mapping Result Collection
#' @description S4 collection of fine-mapping fits for one or more QTL
#'   studies, one row per \code{(study, context, trait, method)} fit. Each
#'   element is a fit's variants: a \code{GRanges} on one chromosome whose
#'   \code{mcols} hold the per-variant values \code{\link{topLoci}} reads.
#'
#' Required \code{mcols()} columns: the identity tuple \code{study},
#' \code{context}, \code{trait} and \code{method}, which is unique, and
#' each row's stored fit \code{susieFit} and cross-validation result
#' \code{cvResult}. The constructor also records \code{fitVariants}, the
#' variants each fit was made with. There is no LD-block axis: a trait's fit
#' is one row wherever its window falls.
#'
#' Optional columns \code{jointStudies}, \code{jointContexts},
#' \code{jointTraits} appear when the collection contains rows produced by a
#' \code{jointSpecification}-driven joint fit, carrying the semicolon-joined
#' member identities of the joined axis.
#' @export
setClass(
    "QtlFineMappingResult",
    contains = "FineMappingResultBase",
    validity = function(object) .validateQtlFineMappingResult(object)
)

# ---- QtlFineMappingResult validity helpers ---------------------------------

# Collect all contract violations (empty vector = valid). The per-entry checks
# run only once the required columns are present; FineMappingResultBase's own
# validity checks the LD references.
# @noRd
#' @importFrom checkmate makeAssertCollection assertNames checkNames
.validateQtlFineMappingResult <- function(object) {
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
    coll$push(.qfmrCheckEntries(object))
    coll$getMessages()
}

# Payload-column + region/traitPos + joint-column + tuple-uniqueness
# contract.
# @noRd
.qfmrCheckEntries <- function(object) {
    jointCols <- intersect(
        c("jointStudies", "jointContexts", "jointTraits"),
        .tupleColumnNames(object)
    )
    c(
        .qfmrCheckEntryLength(object),
        .validateTraitPosColumn(object),
        .validateRegionColumn(object),
        .qfmrCheckJointCols(object, jointCols),
        .qfmrCheckTupleUniqueness(object, jointCols)
    )
}

# The variants and their topLoci ARE the elements now, so what is left
# to check is that the fit payload columns are present and parallel.
# @noRd
.qfmrCheckEntryLength <- function(object) {
    res <- checkNames(
        .tupleColumnNames(object),
        must.include = c("susieFit", "cvResult"),
        what = "colnames"
    )
    if (isTRUE(res)) {
        return(NULL)
    }
    str_c("missing entry payload columns: ", res)
}

# Each present joint* column must be character.
# @noRd
.qfmrCheckJointCols <- function(object, jointCols) {
    list_c(compact(map(jointCols, .qfmrJointColError, object = object)))
}

# @noRd
.qfmrJointColError <- function(jc, object) {
    vals <- .tupleColumn(object, jc)
    if (!is.character(vals)) {
        return(glue(
            "'{jc}' column must be character (got {class(vals)[[1L]]})"
        ))
    }
    NULL
}

# (study, context, trait, method[, joint*]) tuple uniqueness.
# @noRd
.qfmrCheckTupleUniqueness <- function(object, jointCols) {
    keyCols <- c("study", "context", "trait", "method", jointCols)
    # Read column by column: `[` takes no column index on a collection.
    keyTbl <- as_tibble(set_names(
        map(keyCols, .qfmrColOf, object = object),
        keyCols
    ))
    if (nrow(distinct(keyTbl)) < nrow(keyTbl)) {
        return(str_c(
            "(study, context, trait, method[, joint*]) tuple uniqueness ",
            "violated"
        ))
    }
    NULL
}

# @noRd
.qfmrColOf <- function(cn, object) {
    .tupleColumn(object, cn)
}

# The identity-tuple vectors and the payload list are parallel: one entry per
# (study, context, trait, method). A length mismatch would otherwise surface as
# a recycling artefact deep inside the DataFrame assembly.
# @noRd
.qfmrCheckTupleLengths <- function(studyName, context, trait, method, entry) {
    n <- length(studyName)
    ok <- length(context) == n &&
        length(trait) == n &&
        length(method) == n &&
        length(entry) == n
    if (ok) {
        return(invisible(NULL))
    }
    abort(glue(
        "`studyName`, `context`, `trait`, `method`, and `entry` must all ",
        "have the same length."
    ))
}


#' @title Create a QtlFineMappingResult Collection
#' @description Construct a \code{QtlFineMappingResult} collection from
#'   per-tuple vectors and a list of \code{FineMappingRow} payloads (one per
#'   tuple). The optional \code{ldReference} records the LD
#'   reference each RSS-derived fit used; leave it \code{NULL} (the default)
#'   for individual-level fits.
#' @param studyName Character vector of study identifiers (per tuple). Use the
#'   sentinel \code{"joint"} for rows produced by a cross-study joint fit.
#' @param context Character vector of context labels (per tuple). Use
#'   \code{"joint"} for rows produced by a cross-context joint fit.
#' @param trait Character vector of trait identifiers (per tuple). Use
#'   \code{"joint"} for rows produced by a cross-trait joint fit.
#' @param method Character vector of fine-mapping method names (per tuple).
#' @param entry List / \code{SimpleList} of \code{FineMappingRow} objects.
#' @param jointStudies Optional character vector (length \code{length(study)})
#'   listing the semicolon-joined studies participating in each row's
#'   cross-study joint fit, or \code{NA_character_} for non-joint rows. When
#'   \code{NULL} (default) the column is omitted.
#' @param jointContexts Optional character vector for cross-context joints. Same
#'   shape as \code{jointStudies}.
#' @param jointTraits Optional character vector for cross-trait joints. Same
#'   shape as \code{jointStudies}.
#' @param dataType Optional assay / modality label ("bulk_rnaseq",
#'   "methylation"), recorded as a per-row provenance column. Length 1 or
#'   one per row; \code{NULL} (default) omits the column.
#' @param dataSource Optional per-row label for the kind of input the fit saw,
#'   \code{"individual"} or \code{"sumstat"}. Where it is absent, the result
#'   writers read a row's kind from its LD reference: a row with none saw
#'   individual-level data.
#' @param weightStandardized Optional per-row flag: whether the posterior
#'   effect is on the STANDARDIZED genotype scale. The posterior effect IS the
#'   prediction weight, so it carries the same scale question as a TWAS
#'   weight. \code{\link{fineMappingPipeline}} always passes \code{TRUE}:
#'   it standardizes genotypes before fitting, and the posterior mean is
#'   reported on the scale of the \code{X} the fit was given. The argument
#'   exists for a collection assembled by hand, where only the caller knows.
#' @param ldReference The LD reference RSS-derived fits were computed
#'   against -- a genotype panel (see \code{\link{readGenotypes}}), a
#'   \code{GenotypeMixture} or an \code{LdMatrix} -- used for every row, or
#'   a list with one per row; \code{NULL} for individual-level fits. Each
#'   distinct reference is stored once.
#' @param traitPos Optional per-row trait genomic anchor (a \code{GRanges} or
#'   \code{NULL}), carried forward as provenance; not part of the identity key.
#'   \code{NULL} (default) omits the column.
#' @param region Optional per-row analysis WINDOW (a \code{GRanges} or
#'   \code{NULL}) -- the interval each row's fit was run over, carried
#'   forward as provenance; not part of the identity key. \code{NULL}
#'   (default) omits the column, and \code{\link{genomicRegion}} then
#'   derives the window from the span of the variants present.
#' @param metadata A list of free-form metadata about the collection,
#'   as \code{metadata()} holds on any Bioconductor vector. Default
#'   \code{list()}.
#' @return A \code{QtlFineMappingResult} object.
#' @examples
#' tl <- data.frame(variant_id = paste0("chr1:", 100 * 1:3, ":A:G"),
#'   pip = c(0.9, 0.5, 0.1), cs = c(1L, 1L, NA))
#' fe <- FineMappingRow(
#'   variantIds = tl$variant_id, susieFit = list(), topLoci = tl)
#' QtlFineMappingResult(studyName = "s1", context = "brain", trait = "g1",
#'   method = "susie", entry = list(fe))
#' @importFrom checkmate assertCharacter assert checkList
#' @importFrom checkmate checkClass
#' @export
QtlFineMappingResult <- function(
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
    dataType = NULL,
    dataSource = NULL,
    weightStandardized = NULL,
    ldReference = NULL,
    metadata = list()
) {
    .qfmrAssertArgs(
        studyName,
        context,
        trait,
        method,
        entry,
        jointStudies,
        jointContexts,
        jointTraits
    )
    payloads <- .qfmrPayloads(entry)
    cols <- .qfmrMetadataCols(
        studyName,
        context,
        trait,
        method,
        payloads,
        jointStudies,
        jointContexts,
        jointTraits,
        traitPos,
        region,
        dataType,
        dataSource,
        weightStandardized
    )
    .withMetadata(
        .rtlNewValidated("QtlFineMappingResult", payloads, cols, ldReference),
        metadata
    )
}

# The entries as fine-mapping row payloads, each on one chromosome.
# @noRd
.qfmrPayloads <- function(entry) {
    payloads <- map(entry, .asFmRowPayload)
    .checkRowPayloads(payloads, "FineMappingRow", "fine-mapping")
    .qfmrRefuseMultiChromosome(payloads)
    payloads
}

# A QTL fit is one row per (study, context, trait, method), on one chromosome,
# so an entry whose variants span several -- a trans fit -- has no row to go
# in: split per chromosome, it would repeat its tuple and fail as a
# uniqueness violation. Said here, before the split, instead.
# @noRd
.qfmrRefuseMultiChromosome <- function(payloads) {
    bad <- .multiChromosomeEntries(payloads)
    if (length(bad) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "QtlFineMappingResult: `entry` element(s) ",
        "{str_flatten(utils::head(bad, 5L), ', ')} span more than one ",
        "chromosome. A row is one (study, context, trait, method) fit on one ",
        "chromosome, so a fit across chromosomes has no row to go in."
    ))
}

# Every argument's type and length contract, checked before anything is
# built.
# @noRd
.qfmrAssertArgs <- function(
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
    .qfmrCheckTupleLengths(studyName, context, trait, method, entry)
}

# The collection's outer mcols: the identity tuple, each entry's fit / cv
# payload, and any joint-provenance or trait-position columns.
# @noRd
.qfmrMetadataCols <- function(
    studyName,
    context,
    trait,
    method,
    payloads,
    jointStudies,
    jointContexts,
    jointTraits,
    traitPos,
    region,
    dataType,
    dataSource,
    weightStandardized
) {
    n <- length(studyName)
    baseCols <- list(
        study = as.character(studyName),
        context = as.character(context),
        trait = as.character(trait),
        method = as.character(method),
        susieFit = S4Vectors::SimpleList(map(payloads, susieFit)),
        cvResult = S4Vectors::SimpleList(map(payloads, cvResult)),
        fitVariants = IRanges::CharacterList(map(payloads, .fmRowFitVariants))
    )
    withJoint <- .qfmrAppendJointCols(
        baseCols,
        jointStudies,
        jointContexts,
        jointTraits,
        n
    )
    .appendFmProvenanceCols(
        .appendRegionCol(
            .appendTraitPosCol(withJoint, traitPos, n),
            region,
            n
        ),
        dataType,
        dataSource,
        weightStandardized,
        n
    )
}

# Each joint-provenance column must be one value per row.
# @noRd
.qfmrCheckJointColLength <- function(nm, val, n) {
    if (length(val) == n) {
        return(invisible(NULL))
    }
    msg <- glue("`{nm}` must have the same length as `study`.")
    abort(msg)
}

# Append any supplied joint* provenance columns (each must match
# length(study)).
# @noRd
.qfmrAppendJointCols <- function(
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
    walk2(names(supplied), supplied, .qfmrCheckJointColLength, n = n)
    c(cols, map(supplied, as.character))
}

# The single row a (study, context, trait, method) selector pins.
# @noRd
.qfmrSelectRowIndex <- function(x, studyName, context, trait, method) {
    .tupleSelectRow(
        x,
        studyName,
        context,
        trait,
        method,
        cls = "QtlFineMappingResult"
    )
}

#' @rdname fineMappingResult
setMethod(
    "fineMappingResult",
    "QtlFineMappingResult",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        blockId = NULL
    ) {
        .qfmrRefuseBlockId(blockId, "fineMappingResult(QtlFineMappingResult)")
        x[.qfmrSelectRowIndex(x, studyName, context, trait, method)]
    }
)

# A QTL fit is one row, named by (study, context, trait, method), so there is
# no block for `blockId` to pick; one is refused rather than ignored. NA
# counts as unset, as it does for the GWAS-side refusal of context / trait.
# @noRd
.qfmrRefuseBlockId <- function(blockId, caller) {
    if (is.null(blockId) || (length(blockId) == 1L && is.na(blockId))) {
        return(invisible(NULL))
    }
    abort(glue(
        "{caller}: `blockId` names an LD block of a GWAS fit. A QTL fit is ",
        "one row; select it with `studyName`, `context`, `trait` and `method`."
    ))
}

# Derived collection-level accessors (delegate to entry-level methods).

#' @rdname pip
#' @export
setMethod(
    "pip",
    "QtlFineMappingResult",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        idx <- .qfmrSelectRowIndex(x, studyName, context, trait, method)
        .fmrRowPip(.fmrRowParts(x, idx))
    }
)

# Row selector for the base delegating accessors (credibleSets / topLoci /
# marginalEffects / susieFit / variantIds now live on
# FineMappingResultBase, AllClasses.R).
setMethod(
    ".fmrSelectEntry",
    "QtlFineMappingResult",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL,
        blockId = NULL
    ) {
        # The shared FineMappingResultBase accessors pass the whole selector
        # set, and only GwasFineMappingResult is block-indexed. Refusing a
        # non-NULL `blockId` is the point of naming it: a selector absorbed
        # by `...` would silently return the study/context/trait/method match.
        if (!is.null(blockId)) {
            abort(glue(
                "QtlFineMappingResult has no LD blocks: select with ",
                "study / context / trait / method. `blockId` ",
                "applies to GwasFineMappingResult only."
            ))
        }
        .fmrRowParts(
            x,
            .qfmrSelectRowIndex(x, studyName, context, trait, method)
        )
    }
)

#' @rdname cvResult
#' @export
setMethod(
    "cvResult",
    "QtlFineMappingResult",
    function(
        x,
        studyName = NULL,
        context = NULL,
        trait = NULL,
        method = NULL
    ) {
        parts <- .fmrRowParts(
            x,
            .qfmrSelectRowIndex(x, studyName, context, trait, method)
        )
        cvResult(.fitCoverageChecked(parts, "cvResult()"))
    }
)

# topLoci / marginalEffects / susieFit / variantIds are defined once
# on FineMappingResultBase (AllClasses.R), routing through .fmrSelectEntry.

#' @rdname contexts
#' @export
setMethod("contexts", "QtlFineMappingResult", function(x) {
    unique(as.character(x$context))
})

#' @rdname traitNames
#' @export
setMethod("traitNames", "QtlFineMappingResult", function(x) {
    unique(as.character(x$trait))
})
#' @rdname show-methods
#' @export
setMethod("show", "QtlFineMappingResult", function(object) {
    cat(glue("QtlFineMappingResult: {nrow(object)} entries\n", .trim = FALSE))
    if (nrow(object) > 0L) {
        cat(glue(
            "  {n_distinct(object$study)} studies, ",
            "{n_distinct(object$context)} contexts, ",
            "{n_distinct(object$trait)} traits, ",
            "{n_distinct(object$method)} methods\n",
            .trim = FALSE
        ))
    }
    ldSrc <- .ldRefsLabel(
        ldReference(object),
        "NULL (individual-level fit)"
    )
    cat(glue("  LD reference: {ldSrc}\n", .trim = FALSE))
})
