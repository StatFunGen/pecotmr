# =============================================================================
# RgEstimate S4 class
# -----------------------------------------------------------------------------
# RangedTupleList collection of local genetic correlations, one row per pair
# of conditions, keyed by the two conditions' identity tuples (study_1,
# context_1, trait_1, study_2, context_2, trait_2). A pair's element is the
# window its fit covers. The settings the pairs were fitted under are one
# H2EstimateParam for the whole collection, in the `param` slot, and whether
# each pair is reported given its group's other conditions, in the
# `conditional` slot. Produced by estimateLocalRg().
# =============================================================================

#' @include AllGenerics.R RangedTupleList.R H2EstimateParam.R H2Estimate.R
NULL

#' @title Local Genetic Correlations
#' @description Local genetic correlations from
#'   \code{\link{estimateLocalRg}}, one row per pair of conditions, keyed by
#'   the two conditions' identity tuples: \code{study_1}, \code{context_1},
#'   \code{trait_1} and \code{study_2}, \code{context_2}, \code{trait_2}. A
#'   GWAS side has no context or trait, so those are \code{NA}. A pair's
#'   element is the window its fit covers. Every pair in one object was
#'   fitted under the same settings, an \code{\link{H2EstimateParam}} that
#'   \code{\link{h2EstimateParam}} returns, and in one reading: each pair on
#'   its own, or given the other conditions of its group
#'   (\code{conditional = TRUE}); \code{c()} refuses to combine objects that
#'   differ in either.
#' @details Besides the identity, \code{mcols()} hold each pair's values,
#'   read as columns (\code{x$rg}):
#'   \describe{
#'     \item{\code{rg}, \code{rg_se}}{The genetic correlation and its
#'       standard error.}
#'     \item{\code{rg_lower}, \code{rg_upper}}{Its 95\% confidence
#'       interval.}
#'     \item{\code{covariance}, \code{covariance_se}}{The genetic
#'       covariance, on the heritability's scale, and its standard error.}
#'     \item{\code{h2_1}, \code{h2_2}}{Each condition's heritability in the
#'       window, from the fit the covariance rests on.}
#'     \item{\code{covariance_liability}, \code{covariance_liability_se},
#'       \code{h2_1_liability}, \code{h2_2_liability}}{The same with a
#'       case/control GWAS side on the liability scale, where its
#'       \code{GwasSumStats} records its population \code{prevalence}: its
#'       heritability times the factor \code{\linkS4class{H2Estimate}}'s
#'       \code{h2_liability} takes, and the covariance and its standard error
#'       times that factor's root, as ldsc's genetic covariance takes it (a
#'       quantitative side stays on its own scale). \code{NA} for a side, or
#'       a pair, with nothing to convert. rg does not change.}
#'     \item{\code{noise_correlation}}{The correlation of the two conditions'
#'       statistics under no effect: the noise they share through shared
#'       samples.}
#'     \item{\code{noise_source}}{Where that came from: \code{"window"},
#'       fitted in the window; \code{"supplied"}; \code{"fallback"}, supplied
#'       because the window fit failed; \code{"tied"}, the two phenotypes'
#'       correlation on their shared donors less the genetic covariance (a
#'       \code{QtlDataset}); or \code{"zero"}.}
#'     \item{\code{p_value}}{The test of no genetic covariance.}
#'     \item{\code{r2_1}, \code{r2_2}}{With \code{conditional = TRUE}, the
#'       share of each condition's genetic variance that the other conditions
#'       of its unit explain together, with their 95\% intervals
#'       (\code{r2_1_lower}, \code{r2_1_upper} and so on); \code{NA}
#'       otherwise.}
#'     \item{\code{n_snps}}{The variants the pair's fit used.}
#'     \item{\code{n_components}}{The LD components or eigenpairs it used.}
#'     \item{\code{rg_loglik}}{A list column: under the likelihood, the
#'       pair's profile log-likelihood at 101 values of rg evenly spaced on
#'       [-1, 1], every other parameter refitted, relative to its largest (what
#'       \code{\link{summarizeLocalRg}}'s shared summary adds up across
#'       genes); \code{NULL} under the other settings.}
#'   }
#'   A value the settings do not report is \code{NA}.
#' @param x An \code{RgEstimate}.
#' @return \code{h2EstimateParam()}: the \code{\link{H2EstimateParam}} the
#'   pairs were fitted under. \code{isConditional()}: whether each pair is
#'   given the other conditions of its group (\code{TRUE}) or on its own.
#' @seealso \code{\link{estimateLocalRg}}, which builds one.
#' @examples
#' if (requireNamespace("PlinkMatrix", quietly = TRUE)) {
#'     data(qtlSumStatsMulticontextExample)
#'     x <- estimateLocalRg(qtlSumStatsMulticontextExample)
#'     x$rg
#'     h2EstimateParam(x)
#' }
#' @export
setClass(
    "RgEstimate",
    contains = "RangedTupleList",
    representation(param = "H2EstimateParam", conditional = "logical"),
    prototype(param = H2EstimateParam(), conditional = FALSE),
    validity = function(object) .validateRgEstimate(object)
)

# The identity of a pair: the two conditions' (study, context, trait).
# @noRd
.RG_ID_COLS <- c(
    "study_1",
    "context_1",
    "trait_1",
    "study_2",
    "context_2",
    "trait_2"
)

# The values a pair carries, and the type each is held as.
# @noRd
.RG_VALUE_COLS <- c(
    rg = "double",
    rg_se = "double",
    rg_lower = "double",
    rg_upper = "double",
    covariance = "double",
    covariance_se = "double",
    h2_1 = "double",
    h2_2 = "double",
    covariance_liability = "double",
    covariance_liability_se = "double",
    h2_1_liability = "double",
    h2_2_liability = "double",
    noise_correlation = "double",
    noise_source = "character",
    p_value = "double",
    r2_1 = "double",
    r2_1_lower = "double",
    r2_1_upper = "double",
    r2_2 = "double",
    r2_2_lower = "double",
    r2_2_upper = "double",
    n_snps = "integer",
    n_components = "integer"
)

# Where a pair's shared noise came from.
# @noRd
.RG_NOISE_SOURCES <- c("window", "supplied", "fallback", "tied", "zero")

# What a pair holds as an object: its log-likelihood on the rg grid, NULL
# where the settings have no likelihood.
# @noRd
.RG_OBJECT_COLS <- "rg_loglik"

# The rg grid a pair's log-likelihood is kept on.
# @noRd
.RG_CURVE_GRID <- seq(-1, 1, length.out = 101L)

# ---- RgEstimate validity helpers -------------------------------------------

# Collect all contract violations (TRUE when there are none). The type,
# level and uniqueness checks read the columns, so they run once every one is
# present.
# @noRd
#' @importFrom checkmate makeAssertCollection assertNames assertFlag
.validateRgEstimate <- function(object) {
    md <- mcols(object, use.names = FALSE)
    coll <- makeAssertCollection()
    assertFlag(object@conditional, .var.name = "conditional", add = coll)
    assertNames(
        colnames(md),
        must.include = c(.RG_ID_COLS, names(.RG_VALUE_COLS), .RG_OBJECT_COLS),
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
                names(.RG_VALUE_COLS),
                unname(.RG_VALUE_COLS),
                .h2TypeProblem,
                md = md
            ),
            ptype = character()
        ),
        .rgNoiseSourceProblem(md),
        .h2ObjectColumnProblem(.RG_OBJECT_COLS, md),
        .rgCurveProblem(md),
        .rgUniquenessProblem(md)
    )
    if (length(problems) == 0L) TRUE else problems
}

# @noRd
.rgNoiseSourceProblem <- function(md) {
    sources <- md$noise_source
    if (!is.character(sources)) {
        return(character(0))
    }
    unknown <- setdiff(sources[!is.na(sources)], .RG_NOISE_SOURCES)
    if (length(unknown) == 0L) {
        return(character(0))
    }
    as.character(glue(
        "the `noise_source` column holds {str_flatten(unknown, ', ')}; it ",
        "takes {str_flatten(.RG_NOISE_SOURCES, ', ')} or NA."
    ))
}

# Each pair's log-likelihood is NULL or one value per point of the rg grid.
# @noRd
.rgCurveProblem <- function(md) {
    curves <- md$rg_loglik
    if (!is.list(curves) || all(map_lgl(curves, .rgIsCurve))) {
        return(character(0))
    }
    glue(
        "the `rg_loglik` column holds, per pair, NULL or the log-likelihood ",
        "at each of the {length(.RG_CURVE_GRID)} points of the rg grid."
    )
}

# @noRd
.rgIsCurve <- function(curve) {
    is.null(curve) ||
        (is.numeric(curve) && length(curve) == length(.RG_CURVE_GRID))
}

# One row per pair: the two identities name it.
# @noRd
.rgUniquenessProblem <- function(md) {
    if (anyDuplicated(.rtlTupleKeys(md, .RG_ID_COLS, nrow(md))) == 0L) {
        return(character(0))
    }
    str_c(
        "pair uniqueness violated: each pair of (study, context, trait) ",
        "tuples takes one row."
    )
}

# An RgEstimate from one element per pair, the pairs' `values` (one vector per
# column, identity first), the settings `param` they were fitted under and
# whether they are given their groups' other conditions. The result takes its
# genome from `seqinfoOf`.
# @noRd
.rgEstimateNew <- function(
    elements,
    values,
    seqinfoOf,
    param,
    conditional = FALSE
) {
    objects <- map(values[intersect(names(values), .RG_OBJECT_COLS)], I)
    dataFrameArgs <- c(
        values[setdiff(names(values), .RG_OBJECT_COLS)],
        objects,
        list(check.names = FALSE)
    )
    md <- exec(S4Vectors::DataFrame, !!!dataFrameArgs)
    grl <- .rtlRestoreSeqinfo(.rtlBindElements(elements), seqinfoOf)
    new(
        "RgEstimate",
        .rtlNamed(S4Vectors::`mcols<-`(grl, value = md)),
        param = param,
        conditional = conditional
    )
}

# The reading of RgEstimates combined by c(): every part's must be the same,
# as one object holds one.
# @noRd
.rgCombineConditional <- function(parts, fn) {
    flags <- map_lgl(parts, isConditional)
    if (all(flags == flags[[1L]])) {
        return(flags[[1L]])
    }
    abort(c(
        glue(
            "{fn}: some of the RgEstimate objects report each pair given its ",
            "group's other conditions (conditional = TRUE) and some do not, ",
            "and one object holds one reading."
        ),
        i = "Compare them as separate objects."
    ))
}

#' @rdname RgEstimate-class
#' @export
setMethod("h2EstimateParam", "RgEstimate", function(x) {
    x@param
})

#' @rdname RgEstimate-class
#' @export
setMethod("isConditional", "RgEstimate", function(x) {
    x@conditional
})

#' @rdname show-methods
#' @export
setMethod("show", "RgEstimate", function(object) {
    cat(glue("RgEstimate: {length(object)} pair(s)\n", .trim = FALSE))
    cat(glue("settings: {.rgSettingsText(object@param)}\n", .trim = FALSE))
    if (object@conditional) {
        cat("each pair given its group's other conditions\n")
    }
    if (length(object) == 0L) {
        return(invisible(object))
    }
    cols <- c(.RG_ID_COLS, "rg", "rg_lower", "rg_upper", "p_value")
    md <- mcols(object, use.names = FALSE)[cols]
    print(as_tibble(as.data.frame(md)), n = 10L)
    invisible(object)
})

# The core settings on one line, with the shared noise the pairs fit.
# @noRd
.rgSettingsText <- function(param) {
    str_c(.h2SettingsText(param), " / cross noise ", param@crossNoise)
}
