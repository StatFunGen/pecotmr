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
#' @param relatednessIid1 Column name for first individual ID (default "IID1").
#' @param relatednessIid2 Column name for second individual ID (default "IID2").
#' @param relatednessFid1 Column name for first family ID (default NULL).
#' @param relatednessFid2 Column name for second family ID (default NULL).
#' @param relatednessValue Column name for the relatedness measure (default
#'   "PI_HAT").
#' @param phenoData A data.frame with columns \code{IID} and the column named by
#'   \code{phenoCol}. Required when \code{analysisType = "maximizeCases"}.
#' @param phenoCol Column name for the phenotype (default "pheno"). Expected to
#'   be binary (1 = case, 0 = control).
#' @param otherCriterion Optional data.frame with additional filtering criteria
#'   (passed to \code{plinkQC::relatednessFilter}).
#' @param otherCriterionThreshold Threshold for additional criterion.
#' @param otherCriterionDirection Direction for threshold comparison (default
#'   "ge").
#' @param otherCriterionIid Column name for individual ID in criterion data
#'   (default "IID").
#' @param otherCriterionMeasure Column name for the criterion measure.
#' @param maxComponentSize Maximum component size before graph-based pre-pruning
#'   (default 20).
#' @param reduceFraction Fraction of highest-degree nodes to remove per
#'   iteration during pre-pruning (default 0.05).
#' @param maxIterations Maximum plinkQC iterations for resolving remaining
#'   related pairs (default 20).
#' @param verbose Logical, print progress messages (default FALSE).
#' @return A character vector of individual IDs to exclude.
#' @examples
#' rel <- data.frame(IID1 = c("s1", "s2"), IID2 = c("s2", "s3"),
#'   value = c(0.5, 0.1))
#' filterRelatedness(rel, relatednessIid1 = "IID1", relatednessIid2 = "IID2",
#'   relatednessValue = "value", relatednessThreshold = 0.2)
#' @export
filterRelatedness <- function(
    relatedness,
    relatednessThreshold = 0.0625,
    analysisType = c("maximizeUnrelated", "maximizeCases"),
    relatednessIid1 = "IID1",
    relatednessIid2 = "IID2",
    relatednessFid1 = NULL,
    relatednessFid2 = NULL,
    relatednessValue = "PI_HAT",
    phenoData = NULL,
    phenoCol = "pheno",
    otherCriterion = NULL,
    otherCriterionThreshold = NULL,
    otherCriterionDirection = "ge",
    otherCriterionIid = "IID",
    otherCriterionMeasure = NULL,
    maxComponentSize = 20L,
    reduceFraction = 0.05,
    maxIterations = 20L,
    verbose = FALSE
) {
    .relatednessRequirePackages()
    analysisType <- arg_match(analysisType)
    relatedness <- as_tibble(relatedness)
    if (analysisType == "maximizeCases" && is.null(phenoData)) {
        abort("Must provide phenoData when analysisType is 'maximizeCases'")
    }
    # Phase 1: graph-based pre-pruning of large components.
    highRelatedIndiv <- .relatednessPrune(
        relatedness,
        relatednessValue,
        relatednessThreshold,
        relatednessIid1,
        relatednessIid2,
        maxComponentSize,
        reduceFraction,
        verbose
    )
    kin <- .relatednessRemovePruned(
        relatedness,
        highRelatedIndiv,
        relatednessIid1,
        relatednessIid2
    )
    # Phase 2: plinkQC-based filtering (analysis-type dependent).
    plinkqcArgs <- .relatednessBuildPlinkqcArgs(
        otherCriterion,
        relatednessThreshold,
        relatednessIid1,
        relatednessIid2,
        otherCriterionThreshold,
        otherCriterionDirection,
        relatednessFid1,
        relatednessFid2,
        relatednessValue,
        otherCriterionIid,
        otherCriterionMeasure,
        verbose
    )
    filtered <- .relatednessPhase2(
        kin,
        plinkqcArgs,
        analysisType,
        phenoData,
        phenoCol,
        relatednessIid1,
        relatednessIid2
    )
    # Phase 3: iterative cleanup + combine with the graph-pruned individuals.
    cleaned <- .relatednessIterativeCleanup(
        filtered$kin,
        filtered$allExclude,
        plinkqcArgs,
        maxIterations,
        verbose,
        relatednessIid1,
        relatednessIid2,
        relatednessValue,
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
    plinkqcArgs,
    analysisType,
    phenoData,
    phenoCol,
    relatednessIid1,
    relatednessIid2
) {
    if (analysisType == "maximizeUnrelated") {
        return(list(
            allExclude = .relatednessRunPlinkqc(kin, plinkqcArgs)$IID,
            kin = kin
        ))
    }
    .relatednessMaximizeCases(
        kin,
        plinkqcArgs,
        phenoData,
        phenoCol,
        relatednessIid1,
        relatednessIid2
    )
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
    relatednessValue,
    relatednessThreshold,
    relatednessIid1,
    relatednessIid2,
    maxComponentSize,
    reduceFraction,
    verbose
) {
    relatedPairs <- filter(
        relatedness,
        .data[[relatednessValue]] >= relatednessThreshold
    )
    edges <- select(
        relatedPairs,
        all_of(c(relatednessIid1, relatednessIid2))
    )
    # igraph requires a base data.frame (it sets row names on the input).
    workingGraph <- igraph::graph_from_data_frame(
        as.data.frame(edges),
        directed = FALSE
    )
    .relatednessPruneStep(
        workingGraph,
        character(0),
        maxComponentSize,
        reduceFraction,
        verbose
    )
}

# One pruning round: stop once the largest component fits, otherwise drop the
# chosen nodes and recurse on the smaller graph. Each round removes a
# fraction of the largest component, so the recursion is shallow.
# @noRd
.relatednessPruneStep <- function(
    graph,
    removed,
    maxComponentSize,
    reduceFraction,
    verbose
) {
    comp <- igraph::components(graph)
    if (.relatednessLargestComponent(comp) <= maxComponentSize) {
        return(removed)
    }
    .relatednessPruneMessage(comp, verbose, reduceFraction)
    nodesToRemove <- .relatednessNodesToRemove(
        graph,
        comp,
        maxComponentSize,
        reduceFraction
    )
    .relatednessPruneStep(
        igraph::delete_vertices(graph, nodesToRemove),
        c(removed, nodesToRemove),
        maxComponentSize,
        reduceFraction,
        verbose
    )
}

# @noRd
.relatednessPruneMessage <- function(workingComp, verbose, reduceFraction) {
    if (verbose) {
        msg <- glue(
            "Largest component has {max(workingComp$csize)} individuals. ",
            "Removing top {round(reduceFraction * 100)}% ",
            "highest-degree nodes."
        )
        inform(msg)
    }
    invisible(NULL)
}

# The highest-degree nodes to remove across all over-sized components.
# @noRd
.relatednessNodesToRemove <- function(
    workingGraph,
    workingComp,
    maxComponentSize,
    reduceFraction
) {
    largeCompIds <- which(workingComp$csize > maxComponentSize)
    list_c(map(
        largeCompIds,
        .relatednessCompNodesToRemove,
        workingGraph = workingGraph,
        membership = workingComp$membership,
        reduceFraction = reduceFraction
    ))
}

# @noRd
.relatednessCompNodesToRemove <- function(
    compId,
    workingGraph,
    membership,
    reduceFraction
) {
    compNodes <- igraph::V(workingGraph)[membership == compId]
    compDegrees <- igraph::degree(workingGraph, v = compNodes)
    numToRemove <- ceiling(length(compNodes) * reduceFraction)
    names(sort(compDegrees, decreasing = TRUE))[seq_len(numToRemove)]
}

# Drop the pre-pruned individuals from the relatedness data.
# @noRd
.relatednessRemovePruned <- function(
    relatedness,
    highRelatedIndiv,
    relatednessIid1,
    relatednessIid2
) {
    filter(
        relatedness,
        !is_in(.data[[relatednessIid1]], highRelatedIndiv) &
            !is_in(.data[[relatednessIid2]], highRelatedIndiv)
    )
}

# @noRd
.relatednessBuildPlinkqcArgs <- function(
    otherCriterion,
    relatednessThreshold,
    relatednessIid1,
    relatednessIid2,
    otherCriterionThreshold,
    otherCriterionDirection,
    relatednessFid1,
    relatednessFid2,
    relatednessValue,
    otherCriterionIid,
    otherCriterionMeasure,
    verbose
) {
    list(
        otherCriterion = otherCriterion,
        relatednessTh = relatednessThreshold,
        relatednessIID1 = relatednessIid1,
        relatednessIID2 = relatednessIid2,
        otherCriterionTh = otherCriterionThreshold,
        otherCriterionThDirection = otherCriterionDirection,
        relatednessFID1 = relatednessFid1,
        relatednessFID2 = relatednessFid2,
        relatednessRelatedness = relatednessValue,
        otherCriterionIID = otherCriterionIid,
        otherCriterionMeasure = otherCriterionMeasure,
        verbose = verbose
    )
}

# maximizeCases: preserve cases, preferentially remove controls. Returns
# list(allExclude, kin) (kin is restricted to phenotyped individuals).
# @noRd
.relatednessMaximizeCases <- function(
    kin,
    plinkqcArgs,
    phenoData,
    phenoCol,
    relatednessIid1,
    relatednessIid2
) {
    relatedIndividuals <- unique(c(
        kin[[relatednessIid1]],
        kin[[relatednessIid2]]
    ))
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
        is_in(.data[[relatednessIid1]], related$IID) &
            is_in(.data[[relatednessIid2]], related$IID)
    )
    # Step 1: filter among cases.
    caseKin <- filter(
        kin,
        is_in(.data[[relatednessIid1]], relatedCases) &
            is_in(.data[[relatednessIid2]], relatedCases)
    )
    relCases <- .relatednessRunPlinkqc(caseKin, plinkqcArgs)
    casesKeep <- setdiff(relatedCases, relCases$IID)
    # Step 2: remove controls related to retained cases.
    controlsExclude <- .relatednessControlsToExclude(
        kin,
        casesKeep,
        relatedControls,
        relatednessIid1,
        relatednessIid2
    )
    # Step 3: filter among the remaining controls.
    controlsKeep <- setdiff(relatedControls, controlsExclude)
    controlKin <- filter(
        kin,
        is_in(.data[[relatednessIid1]], controlsKeep) &
            is_in(.data[[relatednessIid2]], controlsKeep)
    )
    relControls <- .relatednessRunPlinkqc(controlKin, plinkqcArgs)
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
    relatednessIid1,
    relatednessIid2
) {
    iid1 <- kin[[relatednessIid1]]
    iid2 <- kin[[relatednessIid2]]
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
    plinkqcArgs,
    maxIterations,
    verbose,
    relatednessIid1,
    relatednessIid2,
    relatednessValue,
    relatednessThreshold
) {
    remainingArgs <- list(
        relatednessIid1 = relatednessIid1,
        relatednessIid2 = relatednessIid2,
        relatednessValue = relatednessValue,
        relatednessThreshold = relatednessThreshold
    )
    final <- .relatednessCleanupStep(
        kin = kin,
        allExclude = allExclude,
        iter = 0L,
        maxIterations = maxIterations,
        plinkqcArgs = plinkqcArgs,
        remainingArgs = remainingArgs,
        verbose = verbose
    )
    if (nrow(final$remaining) > 0) {
        msg <- glue(
            "After {maxIterations} iterations, {nrow(final$remaining)} ",
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
    maxIterations,
    plinkqcArgs,
    remainingArgs,
    verbose
) {
    remaining <- exec(.relatednessRemaining, kin, allExclude, !!!remainingArgs)
    if (nrow(remaining) == 0 || iter >= maxIterations) {
        return(list(allExclude = allExclude, remaining = remaining))
    }
    if (verbose) {
        msg <- glue(
            "Iteration {iter + 1L}: {nrow(remaining)} related pairs ",
            "remaining."
        )
        inform(msg)
    }
    additional <- .relatednessRunPlinkqc(remaining, plinkqcArgs)
    .relatednessCleanupStep(
        kin = kin,
        allExclude = c(allExclude, additional$IID),
        iter = iter + 1L,
        maxIterations = maxIterations,
        plinkqcArgs = plinkqcArgs,
        remainingArgs = remainingArgs,
        verbose = verbose
    )
}

# The still-related pairs above threshold after excluding `allExclude`.
# @noRd
.relatednessRemaining <- function(
    kin,
    allExclude,
    relatednessIid1,
    relatednessIid2,
    relatednessValue,
    relatednessThreshold
) {
    remaining <- filter(
        kin,
        !is_in(.data[[relatednessIid1]], allExclude) &
            !is_in(.data[[relatednessIid2]], allExclude)
    )
    filter(remaining, .data[[relatednessValue]] > relatednessThreshold)
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
.relatednessRunPlinkqc <- function(relDf, args) {
    # plinkQC requires a base data.frame (it sets row names on the input).
    rfArgs <- c(list(relatedness = as.data.frame(relDf)), args)
    exec(plinkQC::relatednessFilter, !!!rfArgs)$failIDs
}
