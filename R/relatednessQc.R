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
#' @param analysisType One of \code{"maximize_unrelated"} (default) or
#'   \code{"maximize_cases"}. The latter preserves cases in case-control
#'   studies.
#' @param relatednessIid1 Column name for first individual ID (default "IID1").
#' @param relatednessIid2 Column name for second individual ID (default "IID2").
#' @param relatednessFid1 Column name for first family ID (default NULL).
#' @param relatednessFid2 Column name for second family ID (default NULL).
#' @param relatednessValue Column name for the relatedness measure (default
#'   "PI_HAT").
#' @param phenoData A data.frame with columns \code{IID} and the column named by
#'   \code{phenoCol}. Required when \code{analysisType = "maximize_cases"}.
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
    analysisType = c("maximize_unrelated", "maximize_cases"),
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
    analysisType <- match.arg(analysisType)
    p <- as.list(environment())
    p$relatedness <- as.data.frame(relatedness)
    if (analysisType == "maximize_cases" && is.null(phenoData)) {
        stop("Must provide phenoData when analysisType is 'maximize_cases'")
    }
    # Phase 1: graph-based pre-pruning of large components.
    highRelatedIndiv <- .relatednessPrune(p)
    kin <- .relatednessRemovePruned(p$relatedness, highRelatedIndiv, p)
    # Phase 2: plinkQC-based filtering (analysis-type dependent).
    plinkqcArgs <- .relatednessBuildPlinkqcArgs(p)
    filtered <- .relatednessPhase2(kin, plinkqcArgs, analysisType, p)
    # Phase 3: iterative cleanup + combine with the graph-pruned individuals.
    allExclude <- .relatednessIterativeCleanup(
        filtered$kin,
        filtered$allExclude,
        plinkqcArgs,
        p
    )
    allExclude <- unique(c(allExclude, highRelatedIndiv))
    .relatednessReport(allExclude, p)
    allExclude
}

# Phase-2 dispatch: maximize_unrelated runs plinkQC directly; maximize_cases
# preserves cases. Returns list(allExclude, kin).
# @noRd
.relatednessPhase2 <- function(kin, plinkqcArgs, analysisType, p) {
    if (analysisType == "maximize_unrelated") {
        return(list(
            allExclude = .relatednessRunPlinkqc(kin, plinkqcArgs)$IID,
            kin = kin
        ))
    }
    .relatednessMaximizeCases(kin, plinkqcArgs, p)
}

# @noRd
.relatednessRequirePackages <- function() {
    # nocov start
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required for filterRelatedness")
    }
    if (!requireNamespace("plinkQC", quietly = TRUE)) {
        stop("Package 'plinkQC' is required for filterRelatedness")
    }
    # nocov end
}

# Graph pre-pruning: iteratively remove the highest-degree nodes of any
# component larger than maxComponentSize. Returns the pruned individuals.
# @noRd
.relatednessPrune <- function(p) {
    relatedPairs <- p$relatedness[
        p$relatedness[[p$relatednessValue]] >= p$relatednessThreshold,
    ]
    edges <- relatedPairs[, c(p$relatednessIid1, p$relatednessIid2)]
    workingGraph <- igraph::graph_from_data_frame(edges, directed = FALSE)
    workingComp <- igraph::components(workingGraph)
    highRelatedIndiv <- character(0)
    while (max(workingComp$csize) > p$maxComponentSize) {
        .relatednessPruneMessage(workingComp, p)
        nodesToRemove <- .relatednessNodesToRemove(workingGraph, workingComp, p)
        highRelatedIndiv <- c(highRelatedIndiv, nodesToRemove)
        workingGraph <- igraph::delete_vertices(workingGraph, nodesToRemove)
        workingComp <- igraph::components(workingGraph)
    }
    highRelatedIndiv
}

# @noRd
.relatednessPruneMessage <- function(workingComp, p) {
    if (p$verbose) {
        message(
            "Largest component has ",
            max(workingComp$csize),
            " individuals. Removing top ",
            round(p$reduceFraction * 100),
            "% highest-degree nodes."
        )
    }
    invisible(NULL)
}

# The highest-degree nodes to remove across all over-sized components.
# @noRd
.relatednessNodesToRemove <- function(workingGraph, workingComp, p) {
    largeCompIds <- which(workingComp$csize > p$maxComponentSize)
    unlist(map(
        largeCompIds,
        .relatednessCompNodesToRemove,
        workingGraph = workingGraph,
        membership = workingComp$membership,
        reduceFraction = p$reduceFraction
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
.relatednessRemovePruned <- function(relatedness, highRelatedIndiv, p) {
    relatedness[
        !(relatedness[[p$relatednessIid1]] %in% highRelatedIndiv) &
            !(relatedness[[p$relatednessIid2]] %in% highRelatedIndiv),
    ]
}

# @noRd
.relatednessBuildPlinkqcArgs <- function(p) {
    list(
        otherCriterion = p$otherCriterion,
        relatednessTh = p$relatednessThreshold,
        relatednessIID1 = p$relatednessIid1,
        relatednessIID2 = p$relatednessIid2,
        otherCriterionTh = p$otherCriterionThreshold,
        otherCriterionThDirection = p$otherCriterionDirection,
        relatednessFID1 = p$relatednessFid1,
        relatednessFID2 = p$relatednessFid2,
        relatednessRelatedness = p$relatednessValue,
        otherCriterionIID = p$otherCriterionIid,
        otherCriterionMeasure = p$otherCriterionMeasure,
        verbose = p$verbose
    )
}

# maximize_cases: preserve cases, preferentially remove controls. Returns
# list(allExclude, kin) (kin is restricted to phenotyped individuals).
# @noRd
.relatednessMaximizeCases <- function(kin, plinkqcArgs, p) {
    phenoData <- as.data.frame(p$phenoData)
    phenoData <- phenoData[!is.na(phenoData[[p$phenoCol]]), ]
    relatedIndividuals <- unique(c(
        kin[[p$relatednessIid1]],
        kin[[p$relatednessIid2]]
    ))
    phenoData <- phenoData[phenoData$IID %in% relatedIndividuals, ]
    relatedCases <- phenoData$IID[phenoData[[p$phenoCol]] == 1]
    relatedControls <- phenoData$IID[phenoData[[p$phenoCol]] == 0]
    kin <- kin[
        kin[[p$relatednessIid1]] %in%
            phenoData$IID &
            kin[[p$relatednessIid2]] %in% phenoData$IID,
    ]
    # Step 1: filter among cases.
    caseKin <- kin[
        kin[[p$relatednessIid1]] %in%
            relatedCases &
            kin[[p$relatednessIid2]] %in% relatedCases,
    ]
    relCases <- .relatednessRunPlinkqc(caseKin, plinkqcArgs)
    casesKeep <- setdiff(relatedCases, relCases$IID)
    # Step 2: remove controls related to retained cases.
    controlsExclude <- .relatednessControlsToExclude(
        kin,
        casesKeep,
        relatedControls,
        p
    )
    # Step 3: filter among the remaining controls.
    controlsKeep <- setdiff(relatedControls, controlsExclude)
    controlKin <- kin[
        kin[[p$relatednessIid1]] %in%
            controlsKeep &
            kin[[p$relatednessIid2]] %in% controlsKeep,
    ]
    relControls <- .relatednessRunPlinkqc(controlKin, plinkqcArgs)
    list(
        allExclude = c(relCases$IID, controlsExclude, relControls$IID),
        kin = kin
    )
}

# Controls related to a retained case (row order preserved; a case--control
# edge excludes the control, mirroring the original per-row if / else-if).
# @noRd
.relatednessControlsToExclude <- function(kin, casesKeep, relatedControls, p) {
    iid1 <- kin[[p$relatednessIid1]]
    iid2 <- kin[[p$relatednessIid2]]
    mask1 <- iid1 %in% casesKeep & iid2 %in% relatedControls
    mask2 <- iid2 %in% casesKeep & iid1 %in% relatedControls
    contrib <- ifelse(mask1, iid2, ifelse(mask2, iid1, NA_character_))
    contrib[!is.na(contrib)]
}

# Iteratively re-run plinkQC on the still-related pairs until none remain or
# maxIterations is hit. Returns the accumulated exclusion set.
# @noRd
.relatednessIterativeCleanup <- function(kin, allExclude, plinkqcArgs, p) {
    remaining <- .relatednessRemaining(kin, allExclude, p)
    iter <- 0L
    while (nrow(remaining) > 0 && iter < p$maxIterations) {
        if (p$verbose) {
            message(
                "Iteration ",
                iter + 1L,
                ": ",
                nrow(remaining),
                " related pairs remaining."
            )
        }
        additional <- .relatednessRunPlinkqc(remaining, plinkqcArgs)
        allExclude <- c(allExclude, additional$IID)
        remaining <- .relatednessRemaining(kin, allExclude, p)
        iter <- iter + 1L
    }
    if (nrow(remaining) > 0) {
        warning(
            "After ",
            p$maxIterations,
            " iterations, ",
            nrow(remaining),
            " related pairs remain."
        )
    }
    allExclude
}

# The still-related pairs above threshold after excluding `allExclude`.
# @noRd
.relatednessRemaining <- function(kin, allExclude, p) {
    remaining <- kin[
        !(kin[[p$relatednessIid1]] %in% allExclude) &
            !(kin[[p$relatednessIid2]] %in% allExclude),
    ]
    remaining[remaining[[p$relatednessValue]] > p$relatednessThreshold, ]
}

# @noRd
.relatednessReport <- function(allExclude, p) {
    if (p$verbose) {
        message(
            length(allExclude),
            " individuals excluded at kinship threshold ",
            p$relatednessThreshold
        )
    }
    invisible(NULL)
}

# Run plinkQC::relatednessFilter with the pre-bound column names + thresholds
# (`args`), returning its $failIDs.
# @noRd
.relatednessRunPlinkqc <- function(relDf, args) {
    do.call(
        plinkQC::relatednessFilter,
        c(list(relatedness = relDf), args)
    )$failIDs
}
