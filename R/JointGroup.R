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
setClass("JointGroup",
  contains = "VIRTUAL",
  representation(conditions = "data.frame"),  # one row per condition (Y/Z column)
  validity = function(object) {
    errors <- character()
    if (!all(c("study", "context", "trait") %in% names(object@conditions))) {
      errors <- c(errors,
        "'conditions' must have columns 'study', 'context', 'trait'")
    } else if (nrow(object@conditions) < 1L) {
      errors <- c(errors, "a group needs >= 1 condition (Y/Z column)")
    }
    if (length(errors) == 0L) TRUE else errors
  })

# ---- IndividualJointGroup ---------------------------------------------------
# `pos` is the per-condition functional position (one per Y column), set only by
# the cross-trait enumerator for fsusie (functional SuSiE over the trait domain);
# empty for every other pattern/method.
setClass("IndividualJointGroup",
  contains = "JointGroup",
  representation(X = "matrix", Y = "matrix", pos = "numeric"),
  validity = function(object) {
    errors <- character()
    if (nrow(object@X) != nrow(object@Y))
      errors <- c(errors, "X and Y must share the sample (row) dimension")
    if (ncol(object@Y) != nrow(object@conditions))
      errors <- c(errors, "ncol(Y) must equal nrow(conditions)")
    if (length(object@pos) > 0L && length(object@pos) != ncol(object@Y))
      errors <- c(errors, "when set, 'pos' must have one entry per Y column")
    if (length(errors) == 0L) TRUE else errors
  })

# ---- SumStatsJointGroup -----------------------------------------------------
setClass("SumStatsJointGroup",
  contains = "JointGroup",
  representation(Z = "matrix", R = "matrix", N = "numeric"),
  validity = function(object) {
    errors <- character()
    if (nrow(object@R) != ncol(object@R))
      errors <- c(errors, "'R' (LD) must be square")
    if (nrow(object@Z) != nrow(object@R))
      errors <- c(errors, "'Z' rows (variants) must match the 'R' dimension")
    if (ncol(object@Z) != nrow(object@conditions))
      errors <- c(errors, "ncol(Z) must equal nrow(conditions)")
    if (length(errors) == 0L) TRUE else errors
  })

# ---- JointDispatchCell ------------------------------------------------------
setClass("JointDispatchCell",
  representation(
    pattern   = "character",   # context / trait / study / composed (a label)
    dataForm  = "character",   # individual / sumstats
    enumerate = "function",    # (data, scope, args) -> list<JointGroup>
    minGroup  = "integer"),    # smallest fittable condition count (joint cells
                               # use >= 2; the univariate cell uses 1)
  validity = function(object) {
    errors <- character()
    if (length(object@dataForm) != 1L ||
        !object@dataForm %in% c("individual", "sumstats"))
      errors <- c(errors, "'dataForm' must be 'individual' or 'sumstats'")
    if (length(object@minGroup) != 1L || object@minGroup < 1L)
      errors <- c(errors, "'minGroup' must be a single integer >= 1")
    if (length(errors) == 0L) TRUE else errors
  })

# ---- Pipeline markers -------------------------------------------------------
# Not empty: the `config` list carries the per-pipeline parameter tail
# (coverage/cvFolds/samplePartition/fitFullData/retainFit/... for fm;
# retainFit/retainFitDetail/cvFolds/... for twas), and dispatch on the concrete
# class selects the result type via `construct()`.
setClass("JointPipeline",
  contains = "VIRTUAL",
  representation(config = "list"))

setClass("FmJointPipeline",   contains = "JointPipeline")
setClass("TwasJointPipeline", contains = "JointPipeline")
