# =============================================================================
# H2Estimate S4 class
# -----------------------------------------------------------------------------
# RangedTupleList collection of heritability estimates, one row each, keyed by
# the identity tuple (study, context, trait): a genome-wide estimate from GWAS
# summary statistics, one per study, with no context or trait, or a trait's
# cis estimate from QTL data. A cis estimate's element is the window it
# covers; a genome-wide estimate covers its whole reference, more chromosomes
# than an element may span, so its element holds no ranges. The settings the
# estimates were fitted under are one H2EstimateParam for the whole
# collection, in the `param` slot. Produced by estimateH2(),
# estimateStratifiedH2() and estimateLocalH2(), and consumed by
# h2EstimateToSldscTrait().
# =============================================================================

#' @include AllGenerics.R RangedTupleList.R H2EstimateParam.R
NULL

#' @title Heritability Estimates
#' @description Heritability estimates from \code{\link{estimateH2}},
#'   \code{\link{estimateStratifiedH2}} and
#'   \code{\link{estimateLocalH2}}, one row each, keyed by the identity tuple
#'   \code{(study, context, trait)}: a genome-wide estimate from GWAS summary
#'   statistics, one per study, whose \code{context} and \code{trait} are
#'   \code{NA}, or a trait's cis estimate from QTL data. A cis estimate's
#'   element is the window it covers. A genome-wide estimate covers its whole
#'   LD reference, more chromosomes than an element may span, so its element
#'   holds no ranges. Every estimate in one object was fitted under the same
#'   settings, an \code{\link{H2EstimateParam}} that
#'   \code{\link{h2EstimateParam}} returns, so \code{c()} refuses to combine
#'   objects made under different settings.
#' @details Besides the identity, \code{mcols()} hold each estimate's values:
#'   \describe{
#'     \item{\code{h2}, \code{h2_se}}{The heritability and its standard
#'       error.}
#'     \item{\code{h2_lower}, \code{h2_upper}}{Its 95\% interval, where the
#'       settings give one: under the likelihood with its likelihood-ratio
#'       test (HDL-L's settings, the defaults), from inverting that test, the
#'       intercept refitted at each heritability; an end at the bound of the
#'       heritability's range where the likelihood never falls far enough.}
#'     \item{\code{h2_liability}, \code{h2_liability_se},
#'       \code{h2_liability_lower}, \code{h2_liability_upper}}{The same on
#'       the liability scale (Lee et al. 2011), for a genome-wide estimate of
#'       a case/control study that records its population prevalence K
#'       (\code{prevalence} in \code{\link{GwasSumStats}}): each times
#'       \eqn{K^2 (1 - K)^2 / (P (1 - P) \phi(t)^2)}, with t the liability
#'       threshold \eqn{\Phi^{-1}(1 - K)} and P the share of cases among the
#'       samples the fit's N counts. That is one half where N is the
#'       effective sample size \eqn{4 / (1 / n_{case} + 1 / n_{control})}, as
#'       \code{\link{summaryStatsQc}} sets it by default (Grotzinger et al.
#'       2023), and \code{nCase / (nCase + nControl)} where N counts every
#'       sample. \code{NA} for any other estimate, a quantitative trait's
#'       among them.}
#'     \item{\code{intercept}, \code{intercept_se}}{The fitted intercept and
#'       its standard error.}
#'     \item{\code{p_value}}{The test of no heritability.}
#'     \item{\code{n_snps}}{The variants the estimate used.}
#'     \item{\code{n_components}}{The LD components or eigenpairs a cis
#'       estimate used.}
#'     \item{\code{scale}}{The fitted phenotypic variance relative to its
#'       standardized value, where the settings profile it out.}
#'     \item{\code{reference_overlap}}{The share of the study's samples taken
#'       to be in the LD reference (\code{referenceOverlap} in
#'       \code{\link{estimateLocalH2}}), where the settings read it.}
#'     \item{\code{enrichment}}{The per-annotation enrichment table of a
#'       stratified estimate (\code{\link{stratifiedHeritabilityEnrichment}}).}
#'     \item{\code{jackknife}}{The per-block jackknife matrix of the
#'       annotation coefficients (\code{\link{annotationJackknifeCoefs}}).}
#'     \item{\code{score_stats}}{The candidate annotations' score statistics
#'       (\code{\link{scoreStats}}).}
#'   }
#'   A value the settings do not report is \code{NA}, and an object they do
#'   not produce is \code{NULL}. The accessors return one value per estimate,
#'   in row order.
#' @seealso \code{\link{estimateH2}}, \code{\link{estimateStratifiedH2}}
#'   and \code{\link{estimateLocalH2}}, which build one.
#' @references
#'   Lee SH, Wray NR, Goddard ME, Visscher PM (2011). Estimating missing
#'   heritability for disease from genome-wide association studies. Am J Hum
#'   Genet, 88:294-305.
#'
#'   Grotzinger AD, de la Fuente J, \enc{Privé}{Prive} F, Nivard MG,
#'   Tucker-Drob EM (2023). Pervasive downward bias in estimates of
#'   liability-scale heritability in genome-wide association study
#'   meta-analysis: a simple solution. Biol Psychiatry, 93:29-36.
#' @export
setClass(
    "H2Estimate",
    contains = "RangedTupleList",
    representation(param = "H2EstimateParam"),
    prototype(param = H2EstimateParam()),
    validity = function(object) .validateH2Estimate(object)
)

# The identity of an estimate.
# @noRd
.H2_ID_COLS <- c("study", "context", "trait")

# The values an estimate carries, and the type each is held as.
# @noRd
.H2_VALUE_COLS <- c(
    h2 = "double",
    h2_se = "double",
    h2_lower = "double",
    h2_upper = "double",
    h2_liability = "double",
    h2_liability_se = "double",
    h2_liability_lower = "double",
    h2_liability_upper = "double",
    intercept = "double",
    intercept_se = "double",
    p_value = "double",
    n_snps = "integer",
    n_components = "integer",
    scale = "double",
    reference_overlap = "double"
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
        "(study, context, trait) tuple uniqueness violated: each estimate ",
        "takes one row."
    )
}

# ---- Columns ---------------------------------------------------------------

# Column `cn`, one entry per estimate.
# @noRd
.h2Column <- function(x, cn) {
    mcols(x, use.names = FALSE)[[cn]]
}

# An H2Estimate from one element per estimate, with the estimates' `values`
# (one vector per value column), `objects` (one list per object column) and
# the settings `param` they were fitted under. The result takes its genome
# from `seqinfoOf`.
# @noRd
.h2EstimateNew <- function(elements, values, objects, seqinfoOf, param) {
    dataFrameArgs <- c(values, map(objects, I), list(check.names = FALSE))
    md <- exec(S4Vectors::DataFrame, !!!dataFrameArgs)
    grl <- .rtlRestoreSeqinfo(.rtlBindElements(elements), seqinfoOf)
    new(
        "H2Estimate",
        .rtlNamed(S4Vectors::`mcols<-`(grl, value = md)),
        param = param
    )
}

# The object columns of `n` estimates that carry none.
# @noRd
.h2NoObjects <- function(n) {
    set_names(
        rep(list(rep(list(NULL), n)), length(.H2_OBJECT_COLS)),
        .H2_OBJECT_COLS
    )
}

# The settings of estimates combined by c(), an H2Estimate's or an
# RgEstimate's: every part's must be the same, as one object records one set.
# @noRd
.h2CombineParam <- function(parts, fn) {
    params <- map(parts, h2EstimateParam)
    same <- map_lgl(params, identical, params[[1L]])
    if (all(same)) {
        return(params[[1L]])
    }
    abort(c(
        glue(
            "{fn}: the {class(parts[[1L]])[[1L]]} objects were fitted under ",
            "different settings, and one object records one set."
        ),
        i = "Compare them as separate objects (see h2EstimateParam())."
    ))
}


# =============================================================================
# Accessors
# =============================================================================

#' @rdname heritabilityEstimate
#' @export
setMethod("heritabilityEstimate", "H2Estimate", function(x) {
    .h2Column(x, "h2")
})

#' @rdname H2Estimate-class
#' @export
setMethod("h2EstimateParam", "H2Estimate", function(x) {
    x@param
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
# collection of one; one saved before estimates recorded the reference
# overlap, a heritability interval or the liability scale has none to record;
# one saved with a `method` column has it replaced by the settings that
# method ran under.
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
        object <- .h2WithParamSlot(.h2WithLiabilityColumns(
            .h2WithIntervalColumns(.h2WithOverlapColumn(object))
        ))
        methods::callNextMethod(object, ..., verbose = verbose)
    }
)

# @noRd
.h2WithIntervalColumns <- function(object) {
    md <- mcols(object, use.names = FALSE)
    if (is_in("h2_lower", colnames(md))) {
        return(object)
    }
    none <- rep(NA_real_, nrow(md))
    S4Vectors::`mcols<-`(
        object,
        value = cbind(
            md,
            S4Vectors::DataFrame(h2_lower = none, h2_upper = none)
        )
    )
}

# @noRd
.h2WithLiabilityColumns <- function(object) {
    md <- mcols(object, use.names = FALSE)
    if (is_in("h2_liability", colnames(md))) {
        return(object)
    }
    none <- rep(NA_real_, nrow(md))
    S4Vectors::`mcols<-`(
        object,
        value = cbind(
            md,
            S4Vectors::DataFrame(
                h2_liability = none,
                h2_liability_se = none,
                h2_liability_lower = none,
                h2_liability_upper = none
            )
        )
    )
}

# @noRd
.h2WithOverlapColumn <- function(object) {
    md <- mcols(object, use.names = FALSE)
    if (is_in("reference_overlap", colnames(md))) {
        return(object)
    }
    S4Vectors::`mcols<-`(
        object,
        value = cbind(
            md,
            S4Vectors::DataFrame(reference_overlap = rep(NA_real_, nrow(md)))
        )
    )
}

# An H2Estimate saved with a `method` column, before the settings were
# recorded: the column goes, and the settings its one method ran under fill
# the `param` slot. A saved LDER estimate does not say whether its reference
# was in-sample, and is given the out-of-sample eigenvalue taper.
# @noRd
.h2WithParamSlot <- function(object) {
    md <- mcols(object, use.names = FALSE)
    if (!is_in("method", colnames(md))) {
        return(object)
    }
    method <- unique(as.character(md$method))
    stratified <- any(!map_lgl(md$enrichment, is.null))
    # A slot is an attribute; `attr<-` called as a function copies the
    # object, where `slot<-` would change the caller's in place.
    withParam <- `attr<-`(
        object,
        "param",
        .h2SavedMethodSettings(method, stratified)
    )
    S4Vectors::`mcols<-`(
        withParam,
        value = md[setdiff(colnames(md), "method")]
    )
}

# The settings of the one method a saved H2Estimate names.
# @noRd
.h2SavedMethodSettings <- function(method, stratified) {
    local <- c("lava", "hdlL", "fastHer", "greml")
    genomeWide <- c("sldsc", "gldsc", "lder", "hdl")
    if (length(method) != 1L || !is_in(method, c(local, genomeWide))) {
        abort(c(
            glue(
                "updateObject: this H2Estimate holds estimates by ",
                "{str_flatten(method, ', ')}, and one object now records one ",
                "set of settings."
            ),
            i = str_c(
                "Rebuild the estimates with estimateH2(), ",
                "estimateStratifiedH2() or estimateLocalH2()."
            )
        ))
    }
    if (is_in(method, local)) {
        return(.h2LocalSettings(method))
    }
    .h2GenomeWideSettings(method, stratified, inSample = FALSE)
}

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
    enrichment <- methods::slot(object, "enrichment")
    .h2EstimateNew(
        list(GenomicRanges::GRanges()),
        list(
            study = methods::slot(object, "traitName"),
            context = NA_character_,
            trait = NA_character_,
            h2 = methods::slot(object, "estimate"),
            h2_se = methods::slot(object, "estimateSe"),
            h2_lower = NA_real_,
            h2_upper = NA_real_,
            h2_liability = NA_real_,
            h2_liability_se = NA_real_,
            h2_liability_lower = NA_real_,
            h2_liability_upper = NA_real_,
            intercept = methods::slot(object, "intercept"),
            intercept_se = methods::slot(object, "interceptSe"),
            p_value = NA_real_,
            n_snps = methods::slot(object, "nSnps"),
            n_components = NA_integer_,
            scale = NA_real_,
            reference_overlap = NA_real_
        ),
        list(
            enrichment = list(enrichment),
            jackknife = list(methods::slot(object, "annotationJackknifeCoefs")),
            score_stats = list(methods::slot(object, "scoreStats"))
        ),
        GenomicRanges::GRanges(),
        .h2SavedMethodSettings(
            methods::slot(object, "method"),
            !is.null(enrichment)
        )
    )
}


# =============================================================================
# Show
# =============================================================================

#' @rdname show-methods
#' @export
setMethod("show", "H2Estimate", function(object) {
    cat(glue("H2Estimate: {length(object)} estimate(s)\n", .trim = FALSE))
    cat(glue("settings: {.h2SettingsText(object@param)}\n", .trim = FALSE))
    if (length(object) == 0L) {
        return(invisible(object))
    }
    cols <- c(.H2_ID_COLS, "h2", "h2_se", "intercept", "p_value", "n_snps")
    md <- mcols(object, use.names = FALSE)[cols]
    print(as_tibble(as.data.frame(md)), n = 10L)
    invisible(object)
})

# The core settings on one line, in the order H2EstimateParam() takes them.
# @noRd
.h2SettingsText <- function(param) {
    str_flatten(
        c(
            param@statistic,
            param@criterion,
            str_c(param@noise, " noise"),
            str_c("reference ", param@referenceNoise),
            str_c(param@sampleSize, " n"),
            param@inference
        ),
        " / "
    )
}
