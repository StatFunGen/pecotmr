# =============================================================================
# JointGroup S4 hierarchy + dispatch scaffolding
# -----------------------------------------------------------------------------
# The intermediate contract for the unified joint-analysis engine (see
# dev/jointSpecification-s4-refactor.md). Every enumerator emits a list of
# `JointGroup`s; every fitter consumes one. The grammar/parsing half of
# jointSpecification.R and the auto-detection paths both funnel through this.
#
#   JointGroup (VIRTUAL)        the conditions fitted jointly: a data.frame with
#                               one row per fitted condition (= per Y/Z column),
#                               carrying its (study, context, trait) identity.
#     IndividualJointGroup      design = individual-level (X, Y)
#     SumStatsJointGroup        design = summary-statistic (Z, R, N)
#
# The OUTPUT row identity is DERIVED from `conditions`: an axis that takes one
# value across all conditions is fixed (that value); an axis that varies is
# collapsed to "joint" with the distinct members recorded in jointStudies /
# jointContexts / jointTraits. So cross-context / cross-trait / cross-study are
# the single-varying-axis case and composed is the >1-varying-axis case --
# uniformly, with the actual fitted tuples preserved (composed loses nothing).
#
#   JointDispatchCell           one row of the wiring table: (pattern, dataForm)
#                               -> enumerator + minGroup
#   JointPipeline (VIRTUAL)     pipeline marker carrying per-pipeline config
#     FmJointPipeline           fine-mapping  -> QtlFineMappingResult
#     TwasJointPipeline         twas weights  -> TwasWeights
#
# Construction is validated (new() runs validity), so an enumerator cannot emit
# a malformed group and a mistyped dispatch cell fails at package load.
# =============================================================================

#' @include AllGenerics.R
NULL

# ---- JointGroup virtual base ------------------------------------------------
setClass(
    "JointGroup",
    contains = "VIRTUAL",
    # one row per condition (Y/Z column)
    representation(conditions = "data.frame"),
    validity = function(object) {
        errors <- character()
        if (
            !all(is_in(
                c("study", "context", "trait"),
                names(object@conditions)
            ))
        ) {
            errors <- c(
                errors,
                "'conditions' must have columns 'study', 'context', 'trait'"
            )
        } else if (nrow(object@conditions) < 1L) {
            errors <- c(errors, "a group needs >= 1 condition (Y/Z column)")
        }
        if (length(errors) == 0L) TRUE else errors
    }
)

# ---- IndividualJointGroup ---------------------------------------------------
# `pos` is the per-condition functional position (one per Y column), set only by
# the cross-trait enumerator for fsusie (functional SuSiE over the trait
# domain);
# empty for every other pattern/method.
setClass(
    "IndividualJointGroup",
    contains = "JointGroup",
    representation(X = "matrix", Y = "matrix", traitPos = "numeric"),
    validity = function(object) {
        errors <- character()
        if (nrow(object@X) != nrow(object@Y)) {
            errors <- c(errors, "X and Y must share the sample (row) dimension")
        }
        if (ncol(object@Y) != nrow(object@conditions)) {
            errors <- c(errors, "ncol(Y) must equal nrow(conditions)")
        }
        if (
            length(object@traitPos) > 0L &&
                length(object@traitPos) != ncol(object@Y)
        ) {
            errors <- c(
                errors,
                str_c(
                    "when set, 'traitPos' must have one entry per TRAIT ",
                    "(Y column), not per variant"
                )
            )
        }
        if (length(errors) == 0L) TRUE else errors
    }
)

# ---- SumStatsJointGroup -----------------------------------------------------
setClass(
    "SumStatsJointGroup",
    contains = "JointGroup",
    # The LD REFERENCE, not the LD matrix: R is variants x variants, and
    # every group in an enumeration would hold its own dense copy long before
    # any of them is fitted. The sketch is a lazy panel, so the matrix is
    # derived once per group at fit time -- which is what the univariate RSS
    # path in fineMappingPipeline already does.
    representation(Z = "matrix", ldSketch = "LdSketchOrNULL", N = "numeric"),
    validity = function(object) {
        errors <- character()
        if (ncol(object@Z) != nrow(object@conditions)) {
            errors <- c(errors, "ncol(Z) must equal nrow(conditions)")
        }
        if (is.null(rownames(object@Z))) {
            errors <- c(
                errors,
                str_c(
                    "'Z' must carry the variant ids as rownames: they are ",
                    "what the LD matrix is derived over"
                )
            )
        }
        if (length(errors) == 0L) TRUE else errors
    }
)

# ---- JointDispatchCell ------------------------------------------------------
setClass(
    "JointDispatchCell",
    representation(
        pattern = "character", # context / trait / study / composed (a label)
        dataForm = "character", # individual / sumstats
        enumerate = "function", # (data, scope, args) -> list<JointGroup>
        minGroup = "integer"
    ), # smallest fittable condition count (joint cells
    # use >= 2; the univariate cell uses 1)
    validity = function(object) {
        errors <- character()
        if (
            length(object@dataForm) != 1L ||
                !is_in(object@dataForm, c("individual", "sumstats"))
        ) {
            errors <- c(errors, "'dataForm' must be 'individual' or 'sumstats'")
        }
        if (length(object@minGroup) != 1L || object@minGroup < 1L) {
            errors <- c(errors, "'minGroup' must be a single integer >= 1")
        }
        if (length(errors) == 0L) TRUE else errors
    }
)

# ---- Pipeline markers -------------------------------------------------------
# Not empty: the `config` list carries the per-pipeline parameter tail
# (coverage/cvFolds/samplePartition/fitFullData/retainFit/... for fm;
# retainFit/retainFitDetail/cvFolds/... for twas), and dispatch on the concrete
# class selects the result type via `construct()`.
setClass("JointPipeline", contains = "VIRTUAL", representation(config = "list"))

setClass("FmJointPipeline", contains = "JointPipeline")
setClass("TwasJointPipeline", contains = "JointPipeline")

# ---- Internal accessors -----------------------------------------------------
# These classes are engine internals -- none of them appear in NAMESPACE, and
# none is useful to a user on its own -- so they get INTERNAL accessors rather
# than exported `get*` generics. The point of the no-`@` rule is encapsulation:
# keeping every slot read in the file that defines the layout is what buys
# that, and exporting a `getMinGroup()` for a dispatch cell nobody can reach
# would only enlarge the public API.

# @noRd
.jgConditions <- function(g) g@conditions

# @noRd
.jgX <- function(g) g@X

# @noRd
.jgY <- function(g) g@Y

# The genomic midpoint of each TRAIT, one per Y column -- fsusie's functional
# domain coordinate, not a variant position. Named apart from `pos`, which
# everywhere else in the package means the position of a variant.
# @noRd
.jgTraitPos <- function(g) g@traitPos

# @noRd
.jgZ <- function(g) g@Z

# @noRd
.jgLdSketch <- function(g) g@ldSketch

# The variants a summary-statistics group covers, in Z's row order -- which is
# the order the derived LD matrix comes back in.
# @noRd
.jgVariantIds <- function(g) rownames(.jgZ(g))

# The group's LD matrix, derived from its sketch.
#
# Callers must derive ONCE per group and pass the result down: the per-
# condition entry builder runs once per Z column, so calling this from there
# would recompute the same matrix for every condition.
# @noRd
.jgLdMatrix <- function(g) {
    .ldFromSketch(
        .jgLdSketch(g),
        .jgVariantIds(g),
        label = "jointEngine"
    )
}

# @noRd
.jgN <- function(g) g@N

# @noRd
.jcPattern <- function(cell) cell@pattern

# @noRd
.jcDataForm <- function(cell) cell@dataForm

# @noRd
.jcEnumerate <- function(cell) cell@enumerate

# @noRd
.jcMinGroup <- function(cell) cell@minGroup

# @noRd
.jpConfig <- function(pipeline) pipeline@config
