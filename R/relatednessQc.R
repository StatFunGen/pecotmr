#' @title Column Names In A Relatedness Table
#' @description Which columns of the \code{relatedness} data frame hold the
#'   pair identifiers and the relatedness measure. Built by this constructor
#'   rather than passed as a bare list so a misspelled field fails at the call
#'   site instead of silently leaving the default in place.
#' @param iid1,iid2 Column names for the first and second individual ID.
#'   Defaults \code{"IID1"} and \code{"IID2"}.
#' @param fid1,fid2 Column names for the first and second family ID, or
#'   \code{NULL} (default) when the table has none.
#' @param value Column name for the relatedness measure. Default
#'   \code{"PI_HAT"}.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' relatednessColumns(value = "KINSHIP")
#' @export
relatednessColumns <- function(
    iid1 = "IID1",
    iid2 = "IID2",
    fid1 = NULL,
    fid2 = NULL,
    value = "PI_HAT"
) {
    # No `...`: every field is pecotmr's own, so R's argument matching is the
    # check and an unknown name is an "unused argument" error.
    .newMethodConfig(
        NULL,
        defaults = list(
            iid1 = iid1,
            iid2 = iid2,
            fid1 = fid1,
            fid2 = fid2,
            value = value
        ),
        extra = list(),
        label = "relatednessColumns"
    )
}

#' @title Graph Pre-Pruning Controls For filterRelatedness
#' @description Bounds on the graph-based pre-pruning that runs before
#'   plinkQC, and on the iterative cleanup that runs after it.
#' @param maxComponentSize Largest connected component left untouched by
#'   pre-pruning. Default \code{20}.
#' @param reduceFraction Fraction of highest-degree nodes removed per
#'   pre-pruning iteration. Default \code{0.05}.
#' @param maxIterations Maximum plinkQC iterations used to resolve any
#'   remaining related pairs. Default \code{20}.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' relatednessPruning(maxComponentSize = 50L)
#' @export
relatednessPruning <- function(
    maxComponentSize = 20L,
    reduceFraction = 0.05,
    maxIterations = 20L
) {
    .newMethodConfig(
        NULL,
        defaults = list(
            maxComponentSize = maxComponentSize,
            reduceFraction = reduceFraction,
            maxIterations = maxIterations
        ),
        extra = list(),
        label = "relatednessPruning"
    )
}

#' @title Additional Arguments For plinkQC::relatednessFilter
#' @description The secondary-criterion settings \code{filterRelatedness}
#'   forwards to \code{plinkQC::relatednessFilter}, plus anything else that
#'   function accepts. Names are plinkQC's own and \strong{are checked}
#'   against \code{plinkQC::relatednessFilter}'s live formals, so the check
#'   cannot drift from the upstream signature. See
#'   \code{\link{MethodConfig}} for when checking is not possible.
#'
#'   The pair identifiers, the relatedness column and the threshold are
#'   \emph{not} settable here: \code{filterRelatedness} derives those from
#'   its own \code{columns} and \code{relatednessThreshold} arguments, which
#'   its graph pre-pruning uses as well, and injects them.
#' @param otherCriterion Optional data frame of additional filtering criteria.
#' @param otherCriterionTh Threshold for the additional criterion.
#' @param otherCriterionThDirection Direction of the threshold comparison.
#'   Default \code{"ge"}.
#' @param otherCriterionIID Column name for the individual ID in
#'   \code{otherCriterion}. Default \code{"IID"}.
#' @param otherCriterionMeasure Column name for the criterion measure.
#' @param ... Any other \code{plinkQC::relatednessFilter} argument.
#' @return A \code{\link{MethodConfig}} object.
#' @examples
#' plinkQcConfig(otherCriterionThDirection = "le")
#' @export
plinkQcConfig <- function(
    otherCriterion = NULL,
    otherCriterionTh = NULL,
    otherCriterionThDirection = "ge",
    otherCriterionIID = "IID",
    otherCriterionMeasure = NULL,
    ...
) {
    .newMethodConfig(
        "plinkQC::relatednessFilter",
        defaults = list(
            otherCriterion = otherCriterion,
            otherCriterionTh = otherCriterionTh,
            otherCriterionThDirection = otherCriterionThDirection,
            otherCriterionIID = otherCriterionIID,
            otherCriterionMeasure = otherCriterionMeasure
        ),
        extra = list(...),
        label = "plinkQcConfig",
        engine = "plinkQC::relatednessFilter"
    )
}

#' Filter related individuals from a study
#'
#' Iterative greedy algorithm that removes related individuals exceeding a
#' kinship threshold. First reduces large connected components via graph-based
#' pruning (removing highest-degree nodes), then applies
#' \code{plinkQC::relatednessFilter} iteratively until no related pairs remain.
#'
#' @param relatedness A data.frame of pairwise relatedness estimates (e.g. KING
#'   .kin0 output). Must contain columns for IID1, IID2, and relatedness value.
#' @param relatednessThreshold Kinship threshold above which individuals are
#'   considered related (default 0.0625, i.e. 2nd degree).
#' @param analysisType One of \code{"maximizeUnrelated"} (default) or
#'   \code{"maximizeCases"}. The latter preserves cases in case-control
#'   studies.
#' @param columns Which columns of \code{relatedness} hold the pair
#'   identifiers and the relatedness measure. Build with
#'   \code{\link{relatednessColumns}}.
#' @param phenoData A data.frame with columns \code{IID} and the column named by
#'   \code{phenoCol}. Required when \code{analysisType = "maximizeCases"}.
#' @param phenoCol Column name for the phenotype (default "pheno"). Expected to
#'   be binary (1 = case, 0 = control).
#' @param pruning Bounds on the graph pre-pruning and the iterative cleanup.
#'   Build with \code{\link{relatednessPruning}}.
#' @param plinkQcArgs Additional arguments forwarded to
#'   \code{plinkQC::relatednessFilter}. Build with
#'   \code{\link{plinkQcConfig}}.
#' @param verbose Logical, print progress messages (default FALSE).
#' @return A character vector of individual IDs to exclude.
#' @examples
#' rel <- data.frame(IID1 = c("s1", "s2"), IID2 = c("s2", "s3"),
#'   value = c(0.5, 0.1))
#' filterRelatedness(rel, relatednessThreshold = 0.2,
#'   columns = relatednessColumns(value = "value"))
#' @export
filterRelatedness <- function(
    relatedness,
    relatednessThreshold = 0.0625,
    analysisType = c("maximizeUnrelated", "maximizeCases"),
    columns = relatednessColumns(),
    phenoData = NULL,
    phenoCol = "pheno",
    pruning = relatednessPruning(),
    plinkQcArgs = plinkQcConfig(),
    verbose = FALSE
) {
    .relatednessRequirePackages()
    .assertMethodConfig(columns, "relatednessColumns", "columns")
    .assertMethodConfig(pruning, "relatednessPruning", "pruning")
    .assertMethodConfig(plinkQcArgs, "plinkQcConfig", "plinkQc")
    analysisType <- arg_match(analysisType)
    relatedness <- as_tibble(relatedness)
    if (analysisType == "maximizeCases" && is.null(phenoData)) {
        abort("Must provide phenoData when analysisType is 'maximizeCases'")
    }
    # Phase 1: graph-based pre-pruning of large components.
    highRelatedIndiv <- .relatednessPrune(
        relatedness,
        columns,
        relatednessThreshold,
        pruning,
        verbose
    )
    kin <- .relatednessRemovePruned(relatedness, highRelatedIndiv, columns)
    # Phase 2: plinkQC-based filtering (analysis-type dependent). The pair
    # columns and the threshold are pecotmr's own -- the graph pruning above
    # reads the same values -- so they are injected rather than settable in
    # `plinkQc`.
    pqArgs <- .relatednessBuildPlinkQcArgs(
        plinkQcArgs,
        relatednessThreshold,
        columns,
        verbose
    )
    filtered <- .relatednessPhase2(
        kin,
        pqArgs,
        analysisType,
        phenoData,
        phenoCol,
        columns
    )
    # Phase 3: iterative cleanup + combine with the graph-pruned individuals.
    cleaned <- .relatednessIterativeCleanup(
        filtered$kin,
        filtered$allExclude,
        pqArgs,
        pruning,
        verbose,
        columns,
        relatednessThreshold
    )
    allExclude <- unique(c(cleaned, highRelatedIndiv))
    .relatednessReport(allExclude, verbose, relatednessThreshold)
    allExclude
}

# Phase-2 dispatch: maximizeUnrelated runs plinkQC directly; maximizeCases
# preserves cases. Returns list(allExclude, kin).
# @noRd
.relatednessPhase2 <- function(
    kin,
    pqArgs,
    analysisType,
    phenoData,
    phenoCol,
    columns
) {
    if (analysisType == "maximizeUnrelated") {
        return(list(
            allExclude = .relatednessRunPlinkQc(kin, pqArgs)$IID,
            kin = kin
        ))
    }
    .relatednessMaximizeCases(kin, pqArgs, phenoData, phenoCol, columns)
}

# @noRd
.relatednessRequirePackages <- function() {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        abort("Package 'igraph' is required for filterRelatedness")
    }
    if (!requireNamespace("plinkQC", quietly = TRUE)) {
        abort("Package 'plinkQC' is required for filterRelatedness")
    }
}

# Size of the largest component, or 0 when the graph has none. Guards
# max(integer(0)), which warns and returns -Inf -- with no related pairs the
# loop below must simply not run.
# @noRd
.relatednessLargestComponent <- function(workingComp) {
    if (length(workingComp$csize) == 0L) 0L else max(workingComp$csize)
}

# Graph pre-pruning: iteratively remove the highest-degree nodes of any
# component larger than maxComponentSize. Returns the pruned individuals.
# @noRd
.relatednessPrune <- function(
    relatedness,
    columns,
    relatednessThreshold,
    pruning,
    verbose
) {
    relatedPairs <- filter(
        relatedness,
        .data[[columns$value]] >= relatednessThreshold
    )
    edges <- select(
        relatedPairs,
        all_of(c(columns$iid1, columns$iid2))
    )
    # igraph requires a base data.frame (it sets row names on the input).
    workingGraph <- igraph::graph_from_data_frame(
        as.data.frame(edges),
        directed = FALSE
    )
    .relatednessPruneStep(
        workingGraph,
        character(0),
        pruning,
        verbose
    )
}

# One pruning round: stop once the largest component fits, otherwise drop the
# chosen nodes and recurse on the smaller graph. Each round removes a
# fraction of the largest component, so the recursion is shallow.
# @noRd
.relatednessPruneStep <- function(graph, removed, pruning, verbose) {
    comp <- igraph::components(graph)
    if (.relatednessLargestComponent(comp) <= pruning$maxComponentSize) {
        return(removed)
    }
    .relatednessPruneMessage(comp, verbose, pruning)
    nodesToRemove <- .relatednessNodesToRemove(graph, comp, pruning)
    .relatednessPruneStep(
        igraph::delete_vertices(graph, nodesToRemove),
        c(removed, nodesToRemove),
        pruning,
        verbose
    )
}

# @noRd
.relatednessPruneMessage <- function(workingComp, verbose, pruning) {
    if (verbose) {
        pct <- round(pruning$reduceFraction * 100)
        msg <- glue(
            "Largest component has {max(workingComp$csize)} individuals. ",
            "Removing top {pct}% highest-degree nodes."
        )
        inform(msg)
    }
    invisible(NULL)
}

# The highest-degree nodes to remove across all over-sized components.
# @noRd
.relatednessNodesToRemove <- function(workingGraph, workingComp, pruning) {
    largeCompIds <- which(workingComp$csize > pruning$maxComponentSize)
    list_c(map(
        largeCompIds,
        .relatednessCompNodesToRemove,
        workingGraph = workingGraph,
        membership = workingComp$membership,
        pruning = pruning
    ))
}

# @noRd
.relatednessCompNodesToRemove <- function(
    compId,
    workingGraph,
    membership,
    pruning
) {
    compNodes <- igraph::V(workingGraph)[membership == compId]
    compDegrees <- igraph::degree(workingGraph, v = compNodes)
    numToRemove <- ceiling(length(compNodes) * pruning$reduceFraction)
    names(sort(compDegrees, decreasing = TRUE))[seq_len(numToRemove)]
}

# Drop the pre-pruned individuals from the relatedness data.
# @noRd
.relatednessRemovePruned <- function(relatedness, highRelatedIndiv, columns) {
    filter(
        relatedness,
        !is_in(.data[[columns$iid1]], highRelatedIndiv) &
            !is_in(.data[[columns$iid2]], highRelatedIndiv)
    )
}

# @noRd
.relatednessBuildPlinkQcArgs <- function(
    plinkQcArgs,
    relatednessThreshold,
    columns,
    verbose
) {
    # pecotmr's own settings translated into plinkQC's spelling, then the
    # user's own plinkQC arguments on top. The user cannot reach the derived
    # names: plinkQcConfig() would have rejected a duplicate anyway, and these
    # must agree with what the graph pruning used.
    derived <- list(
        relatednessTh = relatednessThreshold,
        relatednessIID1 = columns$iid1,
        relatednessIID2 = columns$iid2,
        relatednessFID1 = columns$fid1,
        relatednessFID2 = columns$fid2,
        relatednessRelatedness = columns$value,
        verbose = verbose
    )
    user <- as.list(plinkQcArgs)
    clash <- intersect(names(user), names(derived))
    if (length(clash) > 0L) {
        # Without this the duplicate surfaces as R's bare "matched by multiple
        # actual arguments" from inside exec(), which says nothing about where
        # the other value came from.
        abort(glue(
            "filterRelatedness: {str_flatten(clash, ', ')} ",
            "{if (length(clash) == 1L) 'is' else 'are'} derived from ",
            "`relatednessThreshold` and `columns`, which the graph pruning ",
            "reads as well, so {if (length(clash) == 1L) 'it' else 'they'} ",
            "cannot also be set in `plinkQc`."
        ))
    }
    c(derived, user)
}

# maximizeCases: preserve cases, preferentially remove controls. Returns
# list(allExclude, kin) (kin is restricted to phenotyped individuals).
# @noRd
.relatednessMaximizeCases <- function(
    kin,
    pqArgs,
    phenoData,
    phenoCol,
    columns
) {
    iid1Col <- columns$iid1
    iid2Col <- columns$iid2
    relatedIndividuals <- unique(c(kin[[iid1Col]], kin[[iid2Col]]))
    related <- as_tibble(phenoData) |>
        filter(!is.na(.data[[phenoCol]])) |>
        filter(is_in(.data$IID, relatedIndividuals))
    relatedCases <- related |>
        filter(.data[[phenoCol]] == 1) |>
        pull("IID")
    relatedControls <- related |>
        filter(.data[[phenoCol]] == 0) |>
        pull("IID")
    kin <- filter(
        kin,
        is_in(.data[[iid1Col]], related$IID) &
            is_in(.data[[iid2Col]], related$IID)
    )
    # Step 1: filter among cases.
    caseKin <- filter(
        kin,
        is_in(.data[[iid1Col]], relatedCases) &
            is_in(.data[[iid2Col]], relatedCases)
    )
    relCases <- .relatednessRunPlinkQc(caseKin, pqArgs)
    casesKeep <- setdiff(relatedCases, relCases$IID)
    # Step 2: remove controls related to retained cases.
    controlsExclude <- .relatednessControlsToExclude(
        kin,
        casesKeep,
        relatedControls,
        columns
    )
    # Step 3: filter among the remaining controls.
    controlsKeep <- setdiff(relatedControls, controlsExclude)
    controlKin <- filter(
        kin,
        is_in(.data[[iid1Col]], controlsKeep) &
            is_in(.data[[iid2Col]], controlsKeep)
    )
    relControls <- .relatednessRunPlinkQc(controlKin, pqArgs)
    list(
        allExclude = c(relCases$IID, controlsExclude, relControls$IID),
        kin = kin
    )
}

# Controls related to a retained case (row order preserved; a case--control
# edge excludes the control, mirroring the original per-row if / else-if).
# @noRd
.relatednessControlsToExclude <- function(
    kin,
    casesKeep,
    relatedControls,
    columns
) {
    iid1 <- kin[[columns$iid1]]
    iid2 <- kin[[columns$iid2]]
    mask1 <- is_in(iid1, casesKeep) & is_in(iid2, relatedControls)
    mask2 <- is_in(iid2, casesKeep) & is_in(iid1, relatedControls)
    contrib <- case_when(
        mask1 ~ iid2,
        mask2 ~ iid1,
        .default = NA_character_
    )
    contrib[!is.na(contrib)]
}

# Iteratively re-run plinkQC on the still-related pairs until none remain or
# maxIterations is hit. Returns the accumulated exclusion set.
# @noRd
.relatednessIterativeCleanup <- function(
    kin,
    allExclude,
    pqArgs,
    pruning,
    verbose,
    columns,
    relatednessThreshold
) {
    final <- .relatednessCleanupStep(
        kin = kin,
        allExclude = allExclude,
        iter = 0L,
        pruning = pruning,
        pqArgs = pqArgs,
        columns = columns,
        relatednessThreshold = relatednessThreshold,
        verbose = verbose
    )
    if (nrow(final$remaining) > 0) {
        iters <- pruning$maxIterations
        msg <- glue(
            "After {iters} iterations, {nrow(final$remaining)} ",
            "related pairs remain."
        )
        warn(msg)
    }
    final$allExclude
}

# One cleanup round: re-run plinkQC on whatever is still related, add its
# exclusions, and recurse until nothing is related or the iteration cap is
# reached. Returns the accumulated exclusions and what is still related.
# @noRd
.relatednessCleanupStep <- function(
    kin,
    allExclude,
    iter,
    pruning,
    pqArgs,
    columns,
    relatednessThreshold,
    verbose
) {
    remaining <- .relatednessRemaining(
        kin,
        allExclude,
        columns,
        relatednessThreshold
    )
    if (nrow(remaining) == 0 || iter >= pruning$maxIterations) {
        return(list(allExclude = allExclude, remaining = remaining))
    }
    if (verbose) {
        msg <- glue(
            "Iteration {iter + 1L}: {nrow(remaining)} related pairs ",
            "remaining."
        )
        inform(msg)
    }
    additional <- .relatednessRunPlinkQc(remaining, pqArgs)
    .relatednessCleanupStep(
        kin = kin,
        allExclude = c(allExclude, additional$IID),
        iter = iter + 1L,
        pruning = pruning,
        pqArgs = pqArgs,
        columns = columns,
        relatednessThreshold = relatednessThreshold,
        verbose = verbose
    )
}

# The still-related pairs above threshold after excluding `allExclude`.
# @noRd
.relatednessRemaining <- function(
    kin,
    allExclude,
    columns,
    relatednessThreshold
) {
    remaining <- filter(
        kin,
        !is_in(.data[[columns$iid1]], allExclude) &
            !is_in(.data[[columns$iid2]], allExclude)
    )
    filter(remaining, .data[[columns$value]] > relatednessThreshold)
}

# @noRd
.relatednessReport <- function(allExclude, verbose, relatednessThreshold) {
    if (verbose) {
        msg <- glue(
            "{length(allExclude)} individuals excluded at kinship ",
            "threshold {relatednessThreshold}"
        )
        inform(msg)
    }
    invisible(NULL)
}

# Run plinkQC::relatednessFilter with the pre-bound column names + thresholds
# (`args`), returning its $failIDs.
# @noRd
.relatednessRunPlinkQc <- function(relDf, args) {
    # plinkQC requires a base data.frame (it sets row names on the input).
    rfArgs <- c(list(relatedness = as.data.frame(relDf)), as.list(args))
    exec(plinkQC::relatednessFilter, !!!rfArgs)$failIDs
}
