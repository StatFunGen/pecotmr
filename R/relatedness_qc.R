#' Filter related individuals from a study
#'
#' Iterative greedy algorithm that removes related individuals exceeding a
#' kinship threshold. First reduces large connected components via graph-based
#' pruning (removing highest-degree nodes), then applies
#' \code{plinkQC::relatednessFilter} iteratively until no related pairs remain.
#'
#' @param relatedness A data.frame of pairwise relatedness estimates (e.g. KING
#'   .kin0 output). Must contain columns for IID1, IID2, and relatedness value.
#' @param relatedness_threshold Kinship threshold above which individuals are
#'   considered related (default 0.0625, i.e. 2nd degree).
#' @param analysis_type One of \code{"maximize_unrelated"} (default) or
#'   \code{"maximize_cases"}. The latter preserves cases in case-control studies.
#' @param relatedness_iid1 Column name for first individual ID (default "IID1").
#' @param relatedness_iid2 Column name for second individual ID (default "IID2").
#' @param relatedness_fid1 Column name for first family ID (default NULL).
#' @param relatedness_fid2 Column name for second family ID (default NULL).
#' @param relatedness_value Column name for the relatedness measure
#'   (default "PI_HAT").
#' @param pheno_data A data.frame with columns \code{IID} and the column named
#'   by \code{pheno_col}. Required when \code{analysis_type = "maximize_cases"}.
#' @param pheno_col Column name for the phenotype (default "pheno"). Expected
#'   to be binary (1 = case, 0 = control).
#' @param other_criterion Optional data.frame with additional filtering criteria
#'   (passed to \code{plinkQC::relatednessFilter}).
#' @param other_criterion_threshold Threshold for additional criterion.
#' @param other_criterion_direction Direction for threshold comparison
#'   (default "ge").
#' @param other_criterion_iid Column name for individual ID in criterion data
#'   (default "IID").
#' @param other_criterion_measure Column name for the criterion measure.
#' @param max_component_size Maximum component size before graph-based
#'   pre-pruning (default 20).
#' @param reduce_fraction Fraction of highest-degree nodes to remove per
#'   iteration during pre-pruning (default 0.05).
#' @param max_iterations Maximum plinkQC iterations for resolving remaining
#'   related pairs (default 20).
#' @param verbose Logical, print progress messages (default FALSE).
#' @return A character vector of individual IDs to exclude.
#' @export
filter_relatedness <- function(
    relatedness,
    relatedness_threshold = 0.0625,
    analysis_type = c("maximize_unrelated", "maximize_cases"),
    relatedness_iid1 = "IID1",
    relatedness_iid2 = "IID2",
    relatedness_fid1 = NULL,
    relatedness_fid2 = NULL,
    relatedness_value = "PI_HAT",
    pheno_data = NULL,
    pheno_col = "pheno",
    other_criterion = NULL,
    other_criterion_threshold = NULL,
    other_criterion_direction = "ge",
    other_criterion_iid = "IID",
    other_criterion_measure = NULL,
    max_component_size = 20L,
    reduce_fraction = 0.05,
    max_iterations = 20L,
    verbose = FALSE) {

  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required for filter_relatedness")
  if (!requireNamespace("plinkQC", quietly = TRUE))
    stop("Package 'plinkQC' is required for filter_relatedness")

  analysis_type <- match.arg(analysis_type)
  relatedness <- as.data.frame(relatedness)

  if (analysis_type == "maximize_cases" && is.null(pheno_data))
    stop("Must provide pheno_data when analysis_type is 'maximize_cases'")

  # --- Phase 1: Graph-based pre-pruning of large components ----
  related_pairs <- relatedness[relatedness[[relatedness_value]] >= relatedness_threshold, ]
  edges <- related_pairs[, c(relatedness_iid1, relatedness_iid2)]
  working_graph <- igraph::graph_from_data_frame(edges, directed = FALSE)
  working_comp <- igraph::components(working_graph)

  high_related_indiv <- character(0)

  while (max(working_comp$csize) > max_component_size) {
    if (verbose) {
      message("Largest component has ", max(working_comp$csize),
              " individuals. Removing top ", round(reduce_fraction * 100),
              "% highest-degree nodes.")
    }
    large_comp_ids <- which(working_comp$csize > max_component_size)
    nodes_to_remove <- character(0)

    for (comp_id in large_comp_ids) {
      comp_nodes <- igraph::V(working_graph)[working_comp$membership == comp_id]
      comp_degrees <- igraph::degree(working_graph, v = comp_nodes)
      num_to_remove <- ceiling(length(comp_nodes) * reduce_fraction)
      high_degree_nodes <- names(sort(comp_degrees, decreasing = TRUE))[seq_len(num_to_remove)]
      nodes_to_remove <- c(nodes_to_remove, high_degree_nodes)
    }

    high_related_indiv <- c(high_related_indiv, nodes_to_remove)
    working_graph <- igraph::delete_vertices(working_graph, nodes_to_remove)
    working_comp <- igraph::components(working_graph)
  }

  # Remove pre-pruned individuals from the relatedness data
  kin <- relatedness[
    !(relatedness[[relatedness_iid1]] %in% high_related_indiv) &
    !(relatedness[[relatedness_iid2]] %in% high_related_indiv), ]

  # --- Phase 2: plinkQC-based filtering ----
  run_plinkqc <- function(rel_df) {
    plinkQC::relatednessFilter(
      relatedness = rel_df,
      otherCriterion = other_criterion,
      relatednessTh = relatedness_threshold,
      relatednessIID1 = relatedness_iid1,
      relatednessIID2 = relatedness_iid2,
      otherCriterionTh = other_criterion_threshold,
      otherCriterionThDirection = other_criterion_direction,
      relatednessFID1 = relatedness_fid1,
      relatednessFID2 = relatedness_fid2,
      relatednessRelatedness = relatedness_value,
      otherCriterionIID = other_criterion_iid,
      otherCriterionMeasure = other_criterion_measure,
      verbose = verbose
    )$failIDs
  }

  if (analysis_type == "maximize_unrelated") {
    rel <- run_plinkqc(kin)
    all_exclude <- rel$IID

  } else {
    # maximize_cases: preserve cases, preferentially remove controls
    pheno_data <- as.data.frame(pheno_data)
    pheno_data <- pheno_data[!is.na(pheno_data[[pheno_col]]), ]

    related_individuals <- unique(c(kin[[relatedness_iid1]], kin[[relatedness_iid2]]))
    pheno_data <- pheno_data[pheno_data$IID %in% related_individuals, ]

    related_cases <- pheno_data$IID[pheno_data[[pheno_col]] == 1]
    related_controls <- pheno_data$IID[pheno_data[[pheno_col]] == 0]

    kin <- kin[
      kin[[relatedness_iid1]] %in% pheno_data$IID &
      kin[[relatedness_iid2]] %in% pheno_data$IID, ]

    # Step 1: Filter among cases
    case_kin <- kin[
      kin[[relatedness_iid1]] %in% related_cases &
      kin[[relatedness_iid2]] %in% related_cases, ]
    rel_cases <- run_plinkqc(case_kin)
    cases_keep <- setdiff(related_cases, rel_cases$IID)

    # Step 2: Remove controls related to retained cases
    controls_exclude <- character(0)
    for (i in seq_len(nrow(kin))) {
      iid1 <- kin[[relatedness_iid1]][i]
      iid2 <- kin[[relatedness_iid2]][i]
      if (iid1 %in% cases_keep && iid2 %in% related_controls) {
        controls_exclude <- c(controls_exclude, iid2)
      } else if (iid2 %in% cases_keep && iid1 %in% related_controls) {
        controls_exclude <- c(controls_exclude, iid1)
      }
    }

    # Step 3: Filter among remaining controls
    controls_keep <- setdiff(related_controls, controls_exclude)
    control_kin <- kin[
      kin[[relatedness_iid1]] %in% controls_keep &
      kin[[relatedness_iid2]] %in% controls_keep, ]
    rel_controls <- run_plinkqc(control_kin)

    all_exclude <- c(rel_cases$IID, controls_exclude, rel_controls$IID)
  }

  # --- Phase 3: Iterative cleanup ----
  remaining <- kin[
    !(kin[[relatedness_iid1]] %in% all_exclude) &
    !(kin[[relatedness_iid2]] %in% all_exclude), ]
  remaining <- remaining[remaining[[relatedness_value]] > relatedness_threshold, ]

  iter <- 0L
  while (nrow(remaining) > 0 && iter < max_iterations) {
    if (verbose)
      message("Iteration ", iter + 1L, ": ", nrow(remaining), " related pairs remaining.")
    additional <- run_plinkqc(remaining)
    all_exclude <- c(all_exclude, additional$IID)
    remaining <- kin[
      !(kin[[relatedness_iid1]] %in% all_exclude) &
      !(kin[[relatedness_iid2]] %in% all_exclude), ]
    remaining <- remaining[remaining[[relatedness_value]] > relatedness_threshold, ]
    iter <- iter + 1L
  }

  if (nrow(remaining) > 0)
    warning("After ", max_iterations, " iterations, ",
            nrow(remaining), " related pairs remain.")

  # Combine with graph-pruned individuals
  all_exclude <- unique(c(all_exclude, high_related_indiv))

  if (verbose)
    message(length(all_exclude), " individuals excluded at kinship threshold ",
            relatedness_threshold)

  all_exclude
}
