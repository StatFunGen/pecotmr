# =============================================================================
# H2Estimate S4 class
# -----------------------------------------------------------------------------
# RangedTupleList collection of heritability estimates, one row each, keyed by
# the identity tuple (study, context, trait, method): a genome-wide estimate
# from GWAS summary statistics, one per study, with no context or trait, or a
# trait's cis estimate from QTL summary statistics. A cis estimate's element
# is the window it covers; a genome-wide estimate covers its whole reference,
# more chromosomes than an element may span, so its element holds no ranges.
# Produced by estimateH2() and consumed by h2EstimateToSldscTrait().
# =============================================================================

#' @include AllGenerics.R RangedTupleList.R
NULL

#' @title Heritability Estimates
#' @description Heritability estimates from \code{\link{estimateH2}}, one row
#'   each, keyed by the identity tuple \code{(study, context, trait,
#'   method)}: a genome-wide estimate from GWAS summary statistics, one per
#'   study, whose \code{context} and \code{trait} are \code{NA}, or a trait's
#'   cis estimate from QTL summary statistics. A cis estimate's element is the
#'   window it covers. A genome-wide estimate covers its whole LD reference,
#'   more chromosomes than an element may span, so its element holds no
#'   ranges.
#' @details Besides the identity, \code{mcols()} hold each estimate's values:
#'   \describe{
#'     \item{\code{h2}, \code{h2_se}}{The heritability and its standard
#'       error.}
#'     \item{\code{intercept}, \code{intercept_se}}{The fitted intercept and
#'       its standard error.}
#'     \item{\code{p_value}}{The test of no heritability.}
#'     \item{\code{n_snps}}{The variants the estimate used.}
#'     \item{\code{n_components}}{The LD components or eigenpairs a cis
#'       estimate used.}
#'     \item{\code{scale}}{FastHer's fitted phenotypic variance relative to
#'       its standardized value.}
#'     \item{\code{enrichment}}{The per-annotation enrichment table of a
#'       stratified estimate (\code{\link{stratifiedHeritabilityEnrichment}}).}
#'     \item{\code{jackknife}}{The per-block jackknife matrix of the
#'       annotation coefficients (\code{\link{annotationJackknifeCoefs}}).}
#'     \item{\code{score_stats}}{The candidate annotations' score statistics
#'       (\code{\link{scoreStats}}).}
#'   }
#'   A value a method does not report is \code{NA}, and an object it does not
#'   produce is \code{NULL}. The accessors return one value per estimate, in
#'   row order.
#' @seealso \code{\link{estimateH2}}, which builds one.
#' @export
setClass(
    "H2Estimate",
    contains = "RangedTupleList",
    validity = function(object) .validateH2Estimate(object)
)

# The identity of an estimate.
# @noRd
.H2_ID_COLS <- c("study", "context", "trait", "method")

# The values an estimate carries, and the type each is held as.
# @noRd
.H2_VALUE_COLS <- c(
    h2 = "double",
    h2_se = "double",
    intercept = "double",
    intercept_se = "double",
    p_value = "double",
    n_snps = "integer",
    n_components = "integer",
    scale = "double"
)

# The objects an estimate carries, each NULL when it has none.
# @noRd
.H2_OBJECT_COLS <- c("enrichment", "jackknife", "score_stats")

# ---- H2Estimate validity helpers -------------------------------------------

# Collect all contract violations (TRUE when there are none). The type and
# uniqueness checks read the columns, so they run once every one is present.
# @noRd
#' @importFrom checkmate makeAssertCollection assertNames
.validateH2Estimate <- function(object) {
    md <- mcols(object, use.names = FALSE)
    coll <- makeAssertCollection()
    assertNames(
        colnames(md),
        must.include = c(.H2_ID_COLS, names(.H2_VALUE_COLS), .H2_OBJECT_COLS),
        what = "colnames",
        .var.name = "columns",
        add = coll
    )
    if (!coll$isEmpty()) {
        return(coll$getMessages())
    }
    problems <- c(
        list_c(
            map2(
                names(.H2_VALUE_COLS),
                unname(.H2_VALUE_COLS),
                .h2TypeProblem,
                md = md
            ),
            ptype = character()
        ),
        list_c(
            map(.H2_OBJECT_COLS, .h2ObjectColumnProblem, md = md),
            ptype = character()
        ),
        .h2UniquenessProblem(md)
    )
    if (length(problems) == 0L) TRUE else problems
}

# @noRd
.h2TypeProblem <- function(cn, type, md) {
    got <- typeof(md[[cn]])
    if (identical(got, type)) {
        return(character(0))
    }
    as.character(glue("the `{cn}` column must be {type} (got {got})."))
}

# @noRd
.h2ObjectColumnProblem <- function(cn, md) {
    if (is.list(md[[cn]])) {
        return(character(0))
    }
    as.character(glue(
        "the `{cn}` column must be a list, holding each row's object or NULL."
    ))
}

# One row per estimate: the identity names it.
# @noRd
.h2UniquenessProblem <- function(md) {
    if (anyDuplicated(.rtlTupleKeys(md, .H2_ID_COLS, nrow(md))) == 0L) {
        return(character(0))
    }
    str_c(
        "(study, context, trait, method) tuple uniqueness violated: each ",
        "estimate takes one row."
    )
}

# ---- Columns ---------------------------------------------------------------

# Column `cn`, one entry per estimate.
# @noRd
.h2Column <- function(x, cn) {
    mcols(x, use.names = FALSE)[[cn]]
}

# An H2Estimate from one element per estimate, with the estimates' `values`
# (one vector per value column) and `objects` (one list per object column).
# The result takes its genome from `seqinfoOf`.
# @noRd
.h2EstimateNew <- function(elements, values, objects, seqinfoOf) {
    dataFrameArgs <- c(values, map(objects, I), list(check.names = FALSE))
    md <- exec(S4Vectors::DataFrame, !!!dataFrameArgs)
    grl <- .rtlRestoreSeqinfo(.rtlBindElements(elements), seqinfoOf)
    new("H2Estimate", .rtlNamed(S4Vectors::`mcols<-`(grl, value = md)))
}

# The object columns of `n` estimates that carry none.
# @noRd
.h2NoObjects <- function(n) {
    set_names(
        rep(list(rep(list(NULL), n)), length(.H2_OBJECT_COLS)),
        .H2_OBJECT_COLS
    )
}


# =============================================================================
# Accessors
# =============================================================================

#' @rdname heritabilityEstimate
#' @export
setMethod("heritabilityEstimate", "H2Estimate", function(x) {
    .h2Column(x, "h2")
})

#' @rdname heritabilityEstimateSe
#' @export
setMethod("heritabilityEstimateSe", "H2Estimate", function(x) {
    .h2Column(x, "h2_se")
})

#' @rdname heritabilityIntercept
#' @export
setMethod("heritabilityIntercept", "H2Estimate", function(x) {
    .h2Column(x, "intercept")
})

#' @rdname heritabilityInterceptSe
#' @export
setMethod("heritabilityInterceptSe", "H2Estimate", function(x) {
    .h2Column(x, "intercept_se")
})

#' @rdname nSnps
#' @export
setMethod("nSnps", "H2Estimate", function(x) .h2Column(x, "n_snps"))

#' @rdname studyName
#' @export
setMethod("studyName", "H2Estimate", function(x) {
    unique(as.character(x$study))
})

#' @rdname traitNames
#' @export
setMethod("traitNames", "H2Estimate", function(x) {
    traits <- as.character(x$trait)
    unique(traits[!is.na(traits)])
})

#' @rdname methodNames
#' @export
setMethod("methodNames", "H2Estimate", function(x) {
    unique(as.character(x$method))
})

#' @rdname annotationJackknifeCoefs
#' @export
setMethod("annotationJackknifeCoefs", "H2Estimate", function(x) {
    .h2Column(x, "jackknife")
})

#' @rdname scoreStats
#' @export
setMethod("scoreStats", "H2Estimate", function(object) {
    .h2Column(object, "score_stats")
})

#' @rdname stratifiedHeritabilityEnrichment
#' @export
setMethod("stratifiedHeritabilityEnrichment", "H2Estimate", function(object) {
    md <- mcols(object, use.names = FALSE)
    ids <- as_tibble(as.data.frame(md[.H2_ID_COLS]))
    list_rbind(c(
        list(ids[0L, ]),
        map(
            seq_along(object),
            .h2EnrichmentWithId,
            tables = md$enrichment,
            ids = ids
        )
    ))
})

# Estimate `k`'s enrichment table, each row led by its identity; NULL when it
# has none.
# @noRd
.h2EnrichmentWithId <- function(k, tables, ids) {
    table <- tables[[k]]
    if (is.null(table)) {
        return(NULL)
    }
    bind_cols(ids[rep(k, nrow(table)), ], as_tibble(table))
}


# =============================================================================
# updateObject
# =============================================================================

# Each slot an H2Estimate was saved under before its slots were renamed, and
# its name then.
# @noRd
.H2_RENAMED_SLOTS <- c(
    h2 = "estimate",
    h2Se = "estimateSe",
    local = "localBlocks",
    tauBlocks = "annotationJackknifeCoefs"
)

# An H2Estimate saved when it held one estimate in slots is rebuilt as a
# collection of one.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "H2Estimate",
    function(object, ..., verbose = FALSE) {
        if (isTRUE(verbose)) {
            message("updateObject(object = 'H2Estimate')")
        }
        if (.h2SavedInSlots(object)) {
            object <- .h2FromSlots(.renameSavedSlots(object, .H2_RENAMED_SLOTS))
        }
        methods::callNextMethod(object, ..., verbose = verbose)
    }
)

# @noRd
.h2SavedInSlots <- function(object) {
    !is.null(attr(object, "estimate", exact = TRUE)) ||
        !is.null(attr(object, "h2", exact = TRUE))
}

# The estimate an H2Estimate held in slots, as one row. The per-block local
# estimates some held are dropped, with a warning: the genome-wide
# estimators make one estimate per study.
# @noRd
.h2FromSlots <- function(object) {
    if (!is.null(methods::slot(object, "localBlocks"))) {
        warn(str_c(
            "updateObject: this H2Estimate's per-block local estimates are ",
            "dropped; the genome-wide estimators make one estimate per ",
            "study."
        ))
    }
    .h2EstimateNew(
        list(GenomicRanges::GRanges()),
        list(
            study = methods::slot(object, "traitName"),
            context = NA_character_,
            trait = NA_character_,
            method = methods::slot(object, "method"),
            h2 = methods::slot(object, "estimate"),
            h2_se = methods::slot(object, "estimateSe"),
            intercept = methods::slot(object, "intercept"),
            intercept_se = methods::slot(object, "interceptSe"),
            p_value = NA_real_,
            n_snps = methods::slot(object, "nSnps"),
            n_components = NA_integer_,
            scale = NA_real_
        ),
        list(
            enrichment = list(methods::slot(object, "enrichment")),
            jackknife = list(methods::slot(object, "annotationJackknifeCoefs")),
            score_stats = list(methods::slot(object, "scoreStats"))
        ),
        GenomicRanges::GRanges()
    )
}


# =============================================================================
# Show
# =============================================================================

#' @rdname show-methods
#' @export
setMethod("show", "H2Estimate", function(object) {
    cat(glue("H2Estimate: {length(object)} estimate(s)\n", .trim = FALSE))
    if (length(object) == 0L) {
        return(invisible(object))
    }
    cols <- c(.H2_ID_COLS, "h2", "h2_se", "intercept", "p_value", "n_snps")
    md <- mcols(object, use.names = FALSE)[cols]
    print(as_tibble(as.data.frame(md)), n = 10L)
    invisible(object)
})
