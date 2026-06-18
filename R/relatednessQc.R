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
#'   \code{"maximize_cases"}. The latter preserves cases in case-control studies.
#' @param relatednessIid1 Column name for first individual ID (default "IID1").
#' @param relatednessIid2 Column name for second individual ID (default "IID2").
#' @param relatednessFid1 Column name for first family ID (default NULL).
#' @param relatednessFid2 Column name for second family ID (default NULL).
#' @param relatednessValue Column name for the relatedness measure
#'   (default "PI_HAT").
#' @param phenoData A data.frame with columns \code{IID} and the column named
#'   by \code{phenoCol}. Required when \code{analysisType = "maximize_cases"}.
#' @param phenoCol Column name for the phenotype (default "pheno"). Expected
#'   to be binary (1 = case, 0 = control).
#' @param otherCriterion Optional data.frame with additional filtering criteria
#'   (passed to \code{plinkQC::relatednessFilter}).
#' @param otherCriterionThreshold Threshold for additional criterion.
#' @param otherCriterionDirection Direction for threshold comparison
#'   (default "ge").
#' @param otherCriterionIid Column name for individual ID in criterion data
#'   (default "IID").
#' @param otherCriterionMeasure Column name for the criterion measure.
#' @param maxComponentSize Maximum component size before graph-based
#'   pre-pruning (default 20).
#' @param reduceFraction Fraction of highest-degree nodes to remove per
#'   iteration during pre-pruning (default 0.05).
#' @param maxIterations Maximum plinkQC iterations for resolving remaining
#'   related pairs (default 20).
#' @param verbose Logical, print progress messages (default FALSE).
#' @return A character vector of individual IDs to exclude.
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
    verbose = FALSE) {

  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required for filterRelatedness")
  if (!requireNamespace("plinkQC", quietly = TRUE))
    stop("Package 'plinkQC' is required for filterRelatedness")

  analysisType <- match.arg(analysisType)
  relatedness <- as.data.frame(relatedness)

  if (analysisType == "maximize_cases" && is.null(phenoData))
    stop("Must provide phenoData when analysisType is 'maximize_cases'")

  # --- Phase 1: Graph-based pre-pruning of large components ----
  relatedPairs <- relatedness[relatedness[[relatednessValue]] >= relatednessThreshold, ]
  edges <- relatedPairs[, c(relatednessIid1, relatednessIid2)]
  workingGraph <- igraph::graph_from_data_frame(edges, directed = FALSE)
  workingComp <- igraph::components(workingGraph)

  highRelatedIndiv <- character(0)

  while (max(workingComp$csize) > maxComponentSize) {
    if (verbose) {
      message("Largest component has ", max(workingComp$csize),
              " individuals. Removing top ", round(reduceFraction * 100),
              "% highest-degree nodes.")
    }
    largeCompIds <- which(workingComp$csize > maxComponentSize)
    nodesToRemove <- character(0)

    for (compId in largeCompIds) {
      compNodes <- igraph::V(workingGraph)[workingComp$membership == compId]
      compDegrees <- igraph::degree(workingGraph, v = compNodes)
      numToRemove <- ceiling(length(compNodes) * reduceFraction)
      highDegreeNodes <- names(sort(compDegrees, decreasing = TRUE))[seq_len(numToRemove)]
      nodesToRemove <- c(nodesToRemove, highDegreeNodes)
    }

    highRelatedIndiv <- c(highRelatedIndiv, nodesToRemove)
    workingGraph <- igraph::delete_vertices(workingGraph, nodesToRemove)
    workingComp <- igraph::components(workingGraph)
  }

  # Remove pre-pruned individuals from the relatedness data
  kin <- relatedness[
    !(relatedness[[relatednessIid1]] %in% highRelatedIndiv) &
    !(relatedness[[relatednessIid2]] %in% highRelatedIndiv), ]

  # --- Phase 2: plinkQC-based filtering ----
  runPlinkqc <- function(relDf) {
    plinkQC::relatednessFilter(
      relatedness = relDf,
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
    )$failIDs
  }

  if (analysisType == "maximize_unrelated") {
    rel <- runPlinkqc(kin)
    allExclude <- rel$IID

  } else {
    # maximize_cases: preserve cases, preferentially remove controls
    phenoData <- as.data.frame(phenoData)
    phenoData <- phenoData[!is.na(phenoData[[phenoCol]]), ]

    relatedIndividuals <- unique(c(kin[[relatednessIid1]], kin[[relatednessIid2]]))
    phenoData <- phenoData[phenoData$IID %in% relatedIndividuals, ]

    relatedCases <- phenoData$IID[phenoData[[phenoCol]] == 1]
    relatedControls <- phenoData$IID[phenoData[[phenoCol]] == 0]

    kin <- kin[
      kin[[relatednessIid1]] %in% phenoData$IID &
      kin[[relatednessIid2]] %in% phenoData$IID, ]

    # Step 1: Filter among cases
    caseKin <- kin[
      kin[[relatednessIid1]] %in% relatedCases &
      kin[[relatednessIid2]] %in% relatedCases, ]
    relCases <- runPlinkqc(caseKin)
    casesKeep <- setdiff(relatedCases, relCases$IID)

    # Step 2: Remove controls related to retained cases
    controlsExclude <- character(0)
    for (i in seq_len(nrow(kin))) {
      iid1 <- kin[[relatednessIid1]][i]
      iid2 <- kin[[relatednessIid2]][i]
      if (iid1 %in% casesKeep && iid2 %in% relatedControls) {
        controlsExclude <- c(controlsExclude, iid2)
      } else if (iid2 %in% casesKeep && iid1 %in% relatedControls) {
        controlsExclude <- c(controlsExclude, iid1)
      }
    }

    # Step 3: Filter among remaining controls
    controlsKeep <- setdiff(relatedControls, controlsExclude)
    controlKin <- kin[
      kin[[relatednessIid1]] %in% controlsKeep &
      kin[[relatednessIid2]] %in% controlsKeep, ]
    relControls <- runPlinkqc(controlKin)

    allExclude <- c(relCases$IID, controlsExclude, relControls$IID)
  }

  # --- Phase 3: Iterative cleanup ----
  remaining <- kin[
    !(kin[[relatednessIid1]] %in% allExclude) &
    !(kin[[relatednessIid2]] %in% allExclude), ]
  remaining <- remaining[remaining[[relatednessValue]] > relatednessThreshold, ]

  iter <- 0L
  while (nrow(remaining) > 0 && iter < maxIterations) {
    if (verbose)
      message("Iteration ", iter + 1L, ": ", nrow(remaining), " related pairs remaining.")
    additional <- runPlinkqc(remaining)
    allExclude <- c(allExclude, additional$IID)
    remaining <- kin[
      !(kin[[relatednessIid1]] %in% allExclude) &
      !(kin[[relatednessIid2]] %in% allExclude), ]
    remaining <- remaining[remaining[[relatednessValue]] > relatednessThreshold, ]
    iter <- iter + 1L
  }

  if (nrow(remaining) > 0)
    warning("After ", maxIterations, " iterations, ",
            nrow(remaining), " related pairs remain.")

  # Combine with graph-pruned individuals
  allExclude <- unique(c(allExclude, highRelatedIndiv))

  if (verbose)
    message(length(allExclude), " individuals excluded at kinship threshold ",
            relatednessThreshold)

  allExclude
}
