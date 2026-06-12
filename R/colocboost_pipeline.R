#' Build LD/X_ref arguments for colocboost based on data type.
#'
#' When LD matrices are genotype X (non-square, rows=samples, cols=variants),
#' passes them as X_ref to colocboost. Otherwise passes as LD (correlation).
#'
#' @param ldList A list of matrices (correlation R or genotype X).
#' @param subset Optional index vector to subset ldList (e.g., from dict_sumstatLD).
#' @return A named list with either `LD = ...` or `X_ref = ...`.
#' @noRd
buildLdArgs <- function(ldList, subset = NULL) {
  if (!is.null(subset)) ldList <- ldList[subset]
  # Detect: if any matrix is non-square, it's genotype X (samples x variants)
  isGeno <- any(sapply(ldList, function(m) nrow(m) != ncol(m)))
  if (isGeno) list(X_ref = ldList) else list(LD = ldList)
}

#' Run colocboost with tryCatch and timing.
#' @importFrom colocboost colocboost
#' @noRd
.runColocboost <- function(label, ...) {
  t1 <- Sys.time()
  res <- tryCatch(
    colocboostAnalysis(...),
    error = function(e) {
      message(label, " failed: ", conditionMessage(e))
      NULL
    }
  )
  list(result = res, time = Sys.time() - t1)
}

.cbCallColocboost <- function(args, dots) {
  if (!requireNamespace("colocboost", quietly = TRUE)) {
    stop("The colocboost package is required for colocboostAnalysis().")
  }
  do.call(colocboost, c(args, dots))
}

#' Convert loaded regional data to ColocBoost inputs
#'
#' @param regionData A list returned by \code{load_multitask_regional_data()}.
#' @return A structured list containing \code{colocboost_input},
#'   \code{qc_input}, and \code{source_info}.
#' @export
regionDataToColocboostInput <- function(regionData) {
  indRecordsFromInput <- function(input) {
    X <- .cbAsNamedList(input$X, "individual")
    Y <- .cbAsNamedList(input$Y, "individual")
    contexts <- intersect(names(X), names(Y))
    records <- list()
    for (context in contexts) {
      if (is.null(X[[context]]) || .cbYNcol(Y[[context]]) == 0) next
      records[[context]] <- list(
        X = X[[context]],
        Y = Y[[context]],
        maf = .cbListValue(input$maf, context),
        X_variance = .cbListValue(input$X_variance, context)
      )
    }
    records
  }

  indInput <- regionDataToIndInput(regionData)
  rssInput <- regionDataToRssInput(regionData)

  indRecords <- indRecordsFromInput(indInput)
  indArgs <- .cbFormatIndividual(indRecords)

  # Wrap each (rss_input, LD_data) pair as a QcResult (with no QC applied)
  # so .cbFormatSumstat consumes a uniform shape regardless of whether the
  # records came from summary_stats_qc or directly from regionData.
  sumstatRecords <- lapply(names(rssInput$rss_input), function(study) {
    QcResult(
      ldData = rssInput$LD_data[[study]],
      rssInput = rssInput$rss_input[[study]],
      preprocess = list(),
      outlierNumber = 0L,
      skipped = FALSE,
      skipReason = ""
    )
  })
  names(sumstatRecords) <- names(rssInput$rss_input)
  sumstatArgs <- .cbFormatSumstat(sumstatRecords)

  outcomeNames <- c(indArgs$outcome_names, names(sumstatArgs$sumstat))
  indArgs$outcome_names <- NULL
  colocboostInput <- .cbMergeArgs(indArgs, sumstatArgs)
  if (length(outcomeNames) > 0) colocboostInput$outcome_names <- outcomeNames

  list(
    colocboost_input = Filter(Negate(is.null), colocboostInput),
    qc_input = list(
      individual = indInput[c("X", "Y", "maf", "X_variance")],
      sumstat = rssInput[c("rss_input", "LD_data")]
    ),
    source_info = list(individual = indInput$source_info,
                       sumstat = rssInput$source_info)
  )
}

#' ColocBoost analysis with optional pipeline QC
#'
#' This wrapper keeps the direct \code{colocboost()} argument surface. All
#' ColocBoost inputs and model parameters are supplied through \code{...}. When
#' no QC options are requested, the call is passed directly to
#' \code{colocboost()}. When QC options are requested, the wrapper
#' inspects named \code{X}/\code{Y} and/or \code{sumstat}/\code{LD}/\code{X_ref}
#' arguments in \code{...}, runs the relevant reusable QC step, and then calls
#' ColocBoost on the cleaned inputs. If the required named inputs are not
#' available, QC is skipped with a warning and the original ColocBoost call is
#' used.
#'
#' @details
#' Use \code{colocboostAnalysis()} the same way you would use
#' \code{colocboost()}: pass the native ColocBoost arguments by
#' name, for example \code{X}, \code{Y}, \code{sumstat}, \code{LD},
#' \code{X_ref}, \code{dict_YX}, \code{dict_sumstatLD},
#' \code{outcome_names}, \code{focal_outcome_idx}, \code{effect_est},
#' \code{effect_se}, \code{effect_n}, \code{M}, and other ColocBoost model or
#' post-processing options. These arguments are forwarded unchanged unless one
#' or more QC controls are requested.
#'
#' Individual-level QC is only attempted when at least one individual QC control
#' is non-\code{NULL} and named \code{X} and \code{Y} inputs are available in
#' \code{...}. Summary-statistic QC is only attempted when \code{zMismatchQc},
#' \code{pip_cutoff_to_skip_sumstat}, \code{impute = TRUE}, or
#' \code{LD_reference_info} is supplied and named \code{sumstat} plus either
#' \code{LD}, \code{X_ref}, or \code{LD_reference_info} are available.
#' \code{zMismatchQc = "none"} means run basic allele/variant harmonization
#' only; it does not run SLALOM/DENTIST
#' LD-mismatch QC. RAISS imputation is controlled separately by
#' \code{impute = TRUE}.
#'
#' If no QC controls are supplied, this function is a thin direct call to
#' \code{colocboost(...)}.
#' When QC removes outcomes, \code{outcome_names} and \code{focal_outcome_idx}
#' are updated to match the post-QC outcome order. If the requested focal outcome
#' is removed by QC, \code{focal_outcome_idx} is set to \code{NULL} with a
#' warning.
#'
#' @param ... Arguments passed to \code{colocboost()}, including
#'   data inputs such as \code{X}, \code{Y}, \code{sumstat}, \code{LD},
#'   \code{X_ref}, \code{dict_YX}, \code{dict_sumstatLD},
#'   \code{outcome_names}, and all ColocBoost model/post-processing options.
#'   QC can only inspect inputs that are supplied by name.
#' @param missingRateThresh,mafCutoff,xvarCutoff,ldReferenceMetaFile,pipCutoffToSkipInd
#'   Individual-level QC controls. If all are \code{NULL}, individual-level QC
#'   is not run.
#' @param keepIndel,pipCutoffToSkipSumstat,zMismatchQc,impute,imputeOpts
#'   Summary-statistic QC controls. \code{zMismatchQc = "none"} runs
#'   basic allele harmonization without
#'   LD-mismatch outlier detection. Imputation is only run when
#'   \code{impute = TRUE}.
#' @param ldReferenceInfo Optional LD reference information for
#'   summary-statistic QC. This is only needed when the native \code{LD} matrix
#'   row/column names or \code{X_ref} column names are missing or are not
#'   parseable genomic variant IDs. It can be a .bim/.pvar/.pvar.zst file path,
#'   a data.frame with variant metadata, or a \code{load_LD_matrix()} result.
#'   This is a QC-only argument and is not passed to
#'   \code{colocboost()}.
#' @param variantConvention Allele order used by native ColocBoost-style
#'   \code{sumstat$variant} and LD/X_ref names when deriving QC inputs:
#'   \code{"A2_A1"} for pecotmr canonical \code{chr:pos:A2:A1}, or
#'   \code{"A1_A2"} for \code{chr:pos:A1:A2}.
#' @return The object returned by \code{colocboost()}.
#' @examples
#' \dontrun{
#' # Direct ColocBoost call without QC.
#' fit <- colocboostAnalysis(X = X, Y = Y, M = 500)
#'
#' # Summary-statistic input with basic allele/variant harmonization only.
#' fit <- colocboostAnalysis(sumstat = sumstat, LD = LD,
#'                           zMismatchQc = "none", M = 500)
#'
#' # Summary-statistic input with LD-mismatch QC and RAISS imputation.
#' fit <- colocboostAnalysis(sumstat = sumstat, LD = LD,
#'                           zMismatchQc = "slalom", impute = TRUE)
#'
#' # Use richer LD metadata from load_LD_matrix() for QC, while still passing
#' # ColocBoost's native LD input.
#' ldData <- load_LD_matrix(ldMetaFile, region)
#' fit <- colocboostAnalysis(sumstat = sumstat, LD = getCorrelation(ldData),
#'                           ldReferenceInfo = ldData, zMismatchQc = "none")
#'
#' # Individual-level input with explicit genotype QC thresholds.
#' fit <- colocboostAnalysis(X = X, Y = Y,
#'                           missingRateThresh = 0.1,
#'                           mafCutoff = 0.0005)
#' }
#' @export
colocboostAnalysis <- function(...,
                               # individual QC
                               missingRateThresh = NULL,
                               mafCutoff = NULL,
                               xvarCutoff = NULL,
                               ldReferenceMetaFile = NULL,
                               pipCutoffToSkipInd = NULL,
                               # sumstat QC
                               keepIndel = TRUE,
                               pipCutoffToSkipSumstat = NULL,
                               zMismatchQc = NULL,
                               impute = FALSE,
                               imputeOpts = list(rcond = 0.01, R2_threshold = 0.6,
                                                 minimum_ld = 5, lamb = 0.01),
                               ldReferenceInfo = NULL,
                               variantConvention = c("A2_A1", "A1_A2")) {
  variantConvention <- match.arg(variantConvention)
  directArgs <- list(...)
  preQcDataOutcomes <- .cbColocboostOutcomeNames(directArgs, preferSupplied = FALSE)
  preQcDisplayOutcomes <- .cbColocboostOutcomeNames(directArgs, preferSupplied = TRUE)
  if (!is.null(zMismatchQc)) zMismatchQc <- .resolveZMismatchQc(zMismatchQc)

  individualQcRequested <- !is.null(missingRateThresh) ||
    !is.null(mafCutoff) || !is.null(xvarCutoff) ||
    !is.null(ldReferenceMetaFile) || !is.null(pipCutoffToSkipInd)
  sumstatQcRequested <- !is.null(zMismatchQc) || isTRUE(impute) ||
    !is.null(pipCutoffToSkipSumstat) || !is.null(ldReferenceInfo)
  qcRequested <- individualQcRequested || sumstatQcRequested
  if (!qcRequested) {
    return(.cbCallColocboost(directArgs, list()))
  }

  X <- directArgs$X
  Y <- directArgs$Y
  sumstat <- directArgs$sumstat
  LD <- directArgs$LD
  X_ref <- directArgs$X_ref
  dict_YX <- directArgs$dict_YX
  dict_sumstatLD <- directArgs$dict_sumstatLD

  qcSkipMessages <- character()
  individualQcInput <- NULL
  if (individualQcRequested) {
    if (!is.null(X) && !is.null(Y)) {
      individualQcInput <- .cbIndividualQcInputFromColocboost(X, Y, dict_YX)
    } else {
      qcSkipMessages <- c(qcSkipMessages,
                          "Individual-level QC requested but named X and Y were not both supplied.")
    }
  }

  sumstatQcInput <- NULL
  if (sumstatQcRequested) {
    if (!is.null(sumstat) && (!is.null(LD) || !is.null(X_ref) || !is.null(ldReferenceInfo))) {
      sumstatQcInput <- tryCatch(
        .cbSumstatQcInputFromColocboost(
          sumstat, LD, X_ref, dict_sumstatLD,
          LD_reference_info = ldReferenceInfo,
          variant_convention = variantConvention
        ),
        error = function(e) {
          qcSkipMessages <<- c(
            qcSkipMessages,
            paste("Summary-statistic QC input could not be prepared:", conditionMessage(e))
          )
          NULL
        }
      )
      if (!is.null(sumstatQcInput)) {
        if (length(sumstatQcInput$skip_reasons) > 0) {
          qcSkipMessages <- c(qcSkipMessages, sumstatQcInput$skip_reasons)
        }
        if (length(sumstatQcInput$rss_input) == 0) {
          sumstatQcInput <- NULL
        }
      }
    } else {
      qcSkipMessages <- c(qcSkipMessages,
                          "Summary-statistic QC requested but named sumstat plus LD, X_ref, or ldReferenceInfo were not supplied.")
    }
  }
  if (is.null(individualQcInput) && is.null(sumstatQcInput)) {
    warning("QC requested but required QC inputs are unavailable. Calling colocboost() directly. ",
            paste(qcSkipMessages, collapse = " "))
    return(.cbCallColocboost(directArgs, list()))
  }
  if (length(qcSkipMessages) > 0) {
    warning(paste(qcSkipMessages, collapse = " "), " Skipping unavailable QC branch.")
  }

  qcArgs <- tryCatch({
    args <- list()
    if (!is.null(individualQcInput)) {
      message("QC track: processing individual-level inputs before ColocBoost.")
      ind <- qcIndividualData(
        X = individualQcInput$X,
        Y = individualQcInput$Y,
        maf = individualQcInput$maf,
        XVariance = individualQcInput$X_variance,
        missingRateThresh = missingRateThresh,
        mafCutoff = mafCutoff,
        xvarCutoff = .cbDefault(xvarCutoff, 0),
        ldReferenceMetaFile = ldReferenceMetaFile,
        keepIndel = keepIndel,
        pipCutoffToSkip = .cbDefault(pipCutoffToSkipInd, 0)
      )
      args <- .cbMergeArgs(args, .cbFormatIndividual(ind))
    }
    if (!is.null(sumstatQcInput) && length(sumstatQcInput$rss_input) > 0) {
      message("QC track: processing summary-statistic inputs before ColocBoost.")
      sumstatQc <- summaryStatsQc(
        rssInput = sumstatQcInput$rss_input,
        ldData = sumstatQcInput$LD_data,
        keepIndel = keepIndel,
        pipCutoffToSkip = .cbDefault(pipCutoffToSkipSumstat, 0),
        zMismatchQc = if (is.null(zMismatchQc)) "none" else zMismatchQc,
        impute = impute,
        imputeOpts = imputeOpts
      )
      args <- .cbMergeArgs(args, .cbFormatSumstat(sumstatQc))
    }
    args
  }, error = function(e) {
    warning("QC requested but skipped: ", conditionMessage(e),
            ". Calling colocboost() directly.")
    NULL
  })

  if (is.null(qcArgs) || length(qcArgs) == 0) {
    return(.cbCallColocboost(directArgs, list()))
  }
  mergedArgs <- .cbMergeArgs(directArgs, qcArgs)
  if (!is.null(qcArgs$LD)) mergedArgs$X_ref <- NULL
  if (!is.null(qcArgs$X_ref)) mergedArgs$LD <- NULL
  postQcDataOutcomes <- .cbColocboostOutcomeNames(mergedArgs, preferSupplied = FALSE)
  if (length(postQcDataOutcomes) > 0) {
    mergedArgs$outcome_names <- .cbResolveQcOutcomeNames(
      preQcDataOutcomes,
      preQcDisplayOutcomes,
      postQcDataOutcomes
    )
    mergedArgs$focal_outcome_idx <- .cbRemapFocalOutcomeIdx(
      focalOutcomeIdx = directArgs$focal_outcome_idx,
      preQcDataOutcomes = preQcDataOutcomes,
      preQcDisplayOutcomes = preQcDisplayOutcomes,
      postQcDataOutcomes = postQcDataOutcomes,
      postQcDisplayOutcomes = mergedArgs$outcome_names
    )
  }
  .cbCallColocboost(Filter(Negate(is.null), mergedArgs), list())
}

#' Multi-trait colocalization analysis protocol pipeline
#'
#' This function performs protocol-level multi-trait colocalization using
#' ColocBoost. It accepts loaded regional data, performs QC once, then runs the
#' requested xQTL-only, joint GWAS, and separate GWAS analyses.
#'
#' @param regionData A region data loaded from \code{load_regional_data}.
#' @param focalTrait Name of trait if perform focaled ColocBoost
#' @param eventFilters A list of pattern for filtering events based on context names. Example: for sQTL, list(type_pattern = ".*clu_(\\d+_[+-?]).*",valid_pattern = "clu_(\\d+_[+-?]):PR:",exclude_pattern = "clu_(\\d+_[+-?]):IN:")
#' @param mafCutoff A scalar to remove variants with maf < mafCutoff, dafault is 0.005.
#' @param pipCutoffToSkipInd A vector of cutoff values for skipping analysis based on PIP values for each context. Default is 0.
#' @param pipCutoffToSkipSumstat A vector of cutoff values for skipping analysis based on PIP values for each sumstat Default is 0.
#' @param zMismatchQc Quality control method to use. Options are "none",
#'   "slalom", or "dentist". \code{NULL} is treated as \code{"none"} for
#'   basic-only summary-stat preprocessing.
#' @param impute Logical; if TRUE, performs imputation for outliers identified in the analysis (default: TRUE).
#' @param imputeOpts A list of imputation options including rcond, R2_threshold, and minimum_ld (default: list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5)).
#'
#'
#' @return A list containing the individual_data and sumstat_data after QC:
#' individual_data contains the following components if exist
#' \itemize{
#'   \item Y: A list of residualized phenotype values for all tasks.
#'   \item X: A list of residualized genotype matrices all tasks.
#' }
#' sumstat_data contains the following components if exist
#' \itemize{
#'   \item sumstats: A list of summary statistics f or the matched LD_info, each sublist contains sumstats, n, var_y from \code{load_rss_data}.
#'   \item LD_info: A list of LD information, each sublist contains LD_variants, LD_matrix, ref_panel  \code{load_LD_matrix}.
#' }
#'
#' @importFrom susieR susie_rss
#' @importFrom purrr imap map_int
#' @export
colocboostPipeline <- function(
  regionData,
  focalTrait = NULL,
  eventFilters = NULL,
  # - analysis
  xqtlColoc = TRUE,
  jointGwas = FALSE,
  separateGwas = FALSE,
  # - individual QC
  mafCutoff = 0.0005,
  pipCutoffToSkipInd = 0,
  # - sumstat QC
  keepIndel = TRUE,
  pipCutoffToSkipSumstat = 0,
  zMismatchQc = NULL,
  impute = TRUE,
  imputeOpts = list(
    rcond = 0.01, R2_threshold = 0.6,
    minimum_ld = 5, lamb = 0.01
  ),
  ...
) {
  # - internal function by filtering events based on eventFilters
  filterEvents <- function(events, filters, condition) {
    # filters is a list of filter specifications
    # Each filter spec must have:
    #   type_pattern: pattern to identify event type
    #   And at least ONE of:
    #   valid_pattern: pattern that must exist in group
    #   exclude_pattern: pattern to exclude

    filteredEvents <- events
    for (filter in filters) {
      if (is.null(filter$type_pattern) ||
        (is.null(filter$valid_pattern) && is.null(filter$exclude_pattern))) {
        stop("Each filter must specify type_pattern and at least one of valid_pattern or exclude_pattern")
      }
      # Get events of this type
      typeEvents <- filteredEvents[grepl(filter$type_pattern, filteredEvents)]

      if (length(typeEvents) == 0) next
      # Apply valid pattern if specified
      if (!is.null(filter$valid_pattern)) {
        validGroups <- unique(gsub(
          filter$type_pattern, "\\1",
          typeEvents[grepl(filter$valid_pattern, typeEvents)]
        ))
        if (length(validGroups) > 0) {
          typeEvents <- events[grepl(paste(validGroups, collapse = "|"), typeEvents)]
        } else {
          typeEvents <- character(0)
        }
      }
      # Apply exclusions if specified
      if (!is.null(filter$exclude_pattern)) {
        typeEvents <- typeEvents[!grepl(filter$exclude_pattern, typeEvents)]
      }
      if (length(typeEvents) == length(events)) {
        message(paste("All events matching", filter$type_pattern, "in", condition, "included in following analysis."))
      } else if (length(typeEvents) == 0) {
        message(paste("No events matching", filter$type_pattern, "in", condition, "pass the filtering."))
        return(NULL)
      } else {
        excludeEvents <- paste0(setdiff(events, typeEvents), collapse = ";")
        message(paste("Some events,", excludeEvents, "in", condition, "are removed."))
      }
      # Update events list
      filteredEvents <- unique(c(
        filteredEvents[!grepl(filter$type_pattern, filteredEvents)],
        typeEvents
      ))
    }

    return(filteredEvents)
  }

  # - extract contexts and studies from region data, handling both pre- and post-QC
  extractContextsStudies <- function(regionData, phenotypesInit = NULL) {
    individualData <- regionData$individual_data
    sumstatData <- regionData$sumstat_data
    phenotypes <- list("individual_contexts" = NULL, "sumstat_studies" = NULL)

    # Extract individual contexts
    if (!is.null(individualData)) {
      if (is.null(phenotypesInit)) {
        # Pre-QC: individualData is a RegionalData (S4)
        phenotypes$individual_contexts <- names(getPhenotypes(individualData))
      } else {
        nullY <- which(sapply(individualData$Y, is.null))
        if (length(nullY) == 0) {
          message("All individual data pass QC steps.")
          phenotypes$individual_contexts <- names(individualData$Y)
        } else if (length(nullY) < length(individualData$Y)) {
          message(paste(
            "Skipping follow-up analysis for individual traits",
            paste(names(individualData$Y)[nullY], collapse = ";"), "after QC."
          ))
          phenotypes$individual_contexts <- names(individualData$Y)[-nullY]
        } else {
          message("No individual data pass QC.")
        }
      }
    } else {
      message(if (is.null(phenotypesInit)) "No individual data in this region!" else "No individual data pass QC.")
    }

    # Extract sumstat studies
    if (!is.null(sumstatData)) {
      if (is.null(phenotypesInit)) {
        phenotypes$sumstat_studies <- unlist(sapply(sumstatData$sumstats, names))
      } else {
        phenotypes$sumstat_studies <- names(sumstatData$sumstats)
        if (length(phenotypesInit$sumstat_studies) == length(phenotypes$sumstat_studies)) {
          message("All sumstat studies pass QC steps.")
        } else {
          message(paste(
            "Skipping follow-up analysis for sumstat studies",
            paste(setdiff(phenotypesInit$sumstat_studies, phenotypes$sumstat_studies), collapse = ";"), "after QC."
          ))
        }
      }
    } else {
      if (is.null(phenotypesInit)) {
        message("No sumstat data in this region!")
      } else if (!is.null(phenotypesInit$sumstat_studies)) {
        message(paste(
          "Skipping follow-up analysis for sumstat studies",
          paste(phenotypesInit$sumstat_studies, collapse = ";"), "after QC."
        ))
      } else {
        message("No sumstat data pass QC.")
      }
    }

    return(phenotypes)
  }

  ####### ========= resolve defaults ======== #######
  zMismatchQc <- .resolveZMismatchQc(zMismatchQc)

  ####### ========= initial output results before QC ======== #######
  analysisResults <- list("xqtl_coloc" = NULL, "joint_gwas" = NULL, "separate_gwas" = NULL)
  analysisResults$computing_time <- list("QC" = NULL, "Analysis" = list("xqtl_coloc" = NULL, "joint_gwas" = NULL, "separate_gwas" = NULL))
  if (!xqtlColoc & !jointGwas & !separateGwas) {
    message("No colocalization has been performed!")
    return(analysisResults)
  }
  phenotypesInit <- extractContextsStudies(regionData)
  if (is.null(phenotypesInit$individual_contexts) & is.null(phenotypesInit$sumstat_studies)) {
    return(analysisResults)
  }
  if (!is.null(phenotypesInit$individual_contexts)) {
    analysisResults$xqtl_coloc <- list(NULL)
  }
  if (!is.null(phenotypesInit$sumstat_studies)) {
    analysisResults$joint_gwas <- list(NULL)
    if (length(phenotypesInit$sumstat_studies) > 1) {
      analysisResults$separate_gwas <- vector("list", length(phenotypesInit$sumstat_studies)) %>% setNames(phenotypesInit$sumstat_studies)
    } else {
      analysisResults$separate_gwas[[1]] <- list(NULL)
      names(analysisResults$separate_gwas) <- phenotypesInit$sumstat_studies
    }
  }

  ####### ========= Filtering events before QC =========== #########
  if (!is.null(eventFilters) & !is.null(regionData$individual_data)) {
    indData <- regionData$individual_data
    YList <- getPhenotypes(indData)
    YNames <- names(YList)
    YFiltered <- lapply(seq_along(YList), function(i) {
      y <- YList[[i]]
      events <- colnames(y)
      condition <- YNames[i]
      filteredEvents <- filterEvents(events, eventFilters, condition)
      if (is.null(filteredEvents)) {
        return(NULL)
      }
      y[, filteredEvents, drop = FALSE]
    }) %>% setNames(YNames)
    # Drop conditions whose events were entirely filtered out so the
    # RegionalData validity is preserved; ones to drop are remembered for
    # downstream QC messaging via a synthetic NULL-Y list.
    keepCond <- !vapply(YFiltered, is.null, logical(1))
    if (!any(keepCond)) {
      regionData$individual_data <- NULL
    } else {
      YClean <- YFiltered[keepCond]
      # Attach a record of dropped conditions for extractContextsStudies()
      # to surface the post-QC "Skipping follow-up analysis" message.
      droppedNames <- names(YFiltered)[!keepCond]
      mafList <- indData@maf
      YCoords <- indData@coordinates
      regionData$individual_data <- RegionalData(
        genotypeMatrix = getGenotypeMatrix(indData),
        phenotypes = YClean,
        covariates = getCovariates(indData)[keepCond],
        scaleResiduals = indData@scaleResiduals,
        maf = if (length(mafList) == length(keepCond)) mafList[keepCond] else mafList,
        region = indData@region,
        droppedSamples = indData@droppedSamples,
        coordinates = if (!is.null(YCoords)) YCoords[keepCond] else NULL
      )
      if (length(droppedNames) > 0) {
        attr(regionData$individual_data, "filtered_out_contexts") <- droppedNames
      }
    }
  }

  ####### ========= QC for the regionData ======== ########
  t01 <- Sys.time()
  regionData <- qcRegionalData(regionData,
    mafCutoff = mafCutoff,
    pipCutoffToSkipInd = pipCutoffToSkipInd,
    keepIndel = keepIndel,
    pipCutoffToSkipSumstat = pipCutoffToSkipSumstat,
    zMismatchQc = zMismatchQc,
    impute = impute,
    imputeOpts = imputeOpts
  )
  phenotypesQc <- extractContextsStudies(regionData, phenotypesInit = phenotypesInit)
  if (!is.null(phenotypesInit$sumstat_studies) &&
      is.null(phenotypesQc$sumstat_studies)) {
    message("No valid summary statistic studies remain after validation.")
  }
  t02 <- Sys.time()
  analysisResults$computing_time$QC <- t02 - t01

  ####### ========= convert QC'd regional data to ColocBoost input ======== ########
  colocboostInput <- regionDataToColocboostInput(regionData)$colocboost_input
  X <- colocboostInput$X
  Y <- colocboostInput$Y
  dict_YX <- colocboostInput$dict_YX
  sumstats <- colocboostInput$sumstat
  dict_sumstatLD <- colocboostInput$dict_sumstatLD
  ldMat <- colocboostInput$LD
  if (is.null(ldMat)) ldMat <- colocboostInput$X_ref


  ####### ========= streamline three types of analyses ======== ########
  if (is.null(X) & is.null(sumstats)) {
    message("No data pass QC and will not perform analyses.")
    return(analysisResults)
  }
  # - run xQTL-only version of ColocBoost
  if (xqtlColoc & !is.null(X)) {
    message(paste("====== Performing xQTL-only ColocBoost on", length(Y), "contexts. ====="))
    traits <- names(Y)
    focalOutcomeIdx <- if (!is.null(focalTrait) && focalTrait %in% traits) which(traits == focalTrait) else NULL
    cbRes <- .runColocboost("xQTL-only ColocBoost",
      X = X, Y = Y, dict_YX = dict_YX,
      outcome_names = traits, focal_outcome_idx = focalOutcomeIdx,
      output_level = 2, ...
    )
    analysisResults["xqtl_coloc"] <- list(cbRes$result)
    analysisResults$computing_time$Analysis$xqtl_coloc <- cbRes$time
  }
  # - run joint GWAS no focaled version of ColocBoost
  if (jointGwas & !is.null(sumstats)) {
    message(paste("====== Performing non-focaled version GWAS-xQTL ColocBoost on", length(Y), "contexts and", length(sumstats), "GWAS. ====="))
    traits <- c(names(Y), names(sumstats))
    ldArgs <- buildLdArgs(ldMat)
    cbRes <- do.call(.runColocboost, c(list("Joint GWAS ColocBoost",
      X = X, Y = Y, sumstat = sumstats,
      dict_YX = dict_YX, dict_sumstatLD = dict_sumstatLD,
      outcome_names = traits, focal_outcome_idx = NULL,
      output_level = 2), ldArgs, list(...)))
    analysisResults["joint_gwas"] <- list(cbRes$result)
    analysisResults$computing_time$Analysis$joint_gwas <- cbRes$time
  }
  # - run focaled version of ColocBoost for each GWAS
  if (separateGwas & !is.null(sumstats)) {
    t31 <- Sys.time()
    resGwasSeparate <- analysisResults$separate_gwas
    for (iGwas in 1:nrow(dict_sumstatLD)) {
      dict <- dict_sumstatLD[iGwas, ]
      currentStudy <- names(sumstats)[dict[1]]
      message(paste("====== Performing focaled version GWAS-xQTL ColocBoost on", length(Y), "contexts and ", currentStudy, "GWAS. ====="))
      traits <- c(names(Y), currentStudy)
      ldArgsSep <- buildLdArgs(ldMat, subset = dict[2])
      cbRes <- do.call(.runColocboost, c(
        list(paste("Separate GWAS ColocBoost for", currentStudy),
          X = X, Y = Y, sumstat = sumstats[dict[1]],
          dict_YX = dict_YX,
          outcome_names = traits, focal_outcome_idx = length(traits),
          output_level = 2), ldArgsSep, list(...)))
      resGwasSeparate[currentStudy] <- list(cbRes$result)
    }
    t32 <- Sys.time()
    analysisResults$separate_gwas <- resGwasSeparate
    analysisResults$computing_time$Analysis$separate_gwas <- list("total" = t32 - t31, "n_studies" = nrow(dict_sumstatLD), "average" = (t32 - t31) / nrow(dict_sumstatLD))
  }

  return(analysisResults)
}

#' Initial QC for the region data loaded from \code{load_regional_data}
#'
#' This compatibility wrapper converts loaded regional data to reusable individual
#' and RSS inputs, runs the shared QC helpers once, and returns the historical
#' post-QC structure consumed by \code{colocboostPipeline()}.
#'
#' @param regionData A region data loaded from \code{load_regional_data}.
#' @param mafCutoff A scalar to remove variants with maf < mafCutoff.
#' @param pipCutoffToSkipInd A vector of cutoff values for skipping individual contexts.
#' @param keepIndel Logical; if \code{FALSE}, remove indel variants during
#'   summary-statistic allele harmonization.
#' @param pipCutoffToSkipSumstat A vector of cutoff values for skipping summary-stat studies.
#' @param zMismatchQc Quality control method to use. Options are "none",
#'   "slalom", or "dentist". \code{NULL} is treated as \code{"none"} for
#'   basic-only summary-stat preprocessing.
#' @param impute Logical; if TRUE, performs imputation when required metadata are available.
#' @param imputeOpts A list of imputation options.
#' @return A list containing post-QC \code{individual_data} and \code{sumstat_data}.
#' @export
qcRegionalData <- function(regionData,
                           # - individual
                           mafCutoff = 0.0005,
                           pipCutoffToSkipInd = 0,
                           # - sumstat
                           keepIndel = TRUE,
                           pipCutoffToSkipSumstat = 0,
                           zMismatchQc = NULL,
                           impute = FALSE,
                           imputeOpts = list(rcond = 0.01, R2_threshold = 0.6, minimum_ld = 5, lamb = 0.01)) {
  zMismatchQc <- .resolveZMismatchQc(zMismatchQc)
  qcedIndividualToRegionData <- function(indQc) {
    if (is.null(indQc) || length(indQc) == 0) return(NULL)
    list(
      X = lapply(indQc, `[[`, "X"),
      Y = lapply(indQc, `[[`, "Y")
    )
  }
  qcedSumstatToRegionData <- function(sumstatQc) {
    if (is.null(sumstatQc) || length(sumstatQc) == 0) return(NULL)
    if (is(sumstatQc, "QcResult")) {
      sumstatQc <- list(study1 = sumstatQc)
    }
    sumstats <- lapply(sumstatQc, getRssInput)
    LD_data <- list()
    LD_match <- character()
    ldVariantIndex <- list()
    for (study in names(sumstatQc)) {
      ldObj <- getLdData(sumstatQc[[study]])
      variantKey <- paste(if (is.null(ldObj)) "" else getVariantIds(ldObj), collapse = ",")
      if (variantKey %in% names(ldVariantIndex)) {
        LD_match <- c(LD_match, ldVariantIndex[[variantKey]])
      } else {
        LD_data[[study]] <- ldObj
        ldVariantIndex[[variantKey]] <- study
        LD_match <- c(LD_match, study)
      }
    }
    list(sumstats = sumstats, LD_data = LD_data, LD_match = LD_match)
  }

  individualData <- NULL
  indInput <- regionDataToIndInput(regionData)
  if (isTRUE(indInput$source_info$has_individual)) {
    indQc <- qcIndividualData(
      X = indInput$X,
      Y = indInput$Y,
      maf = indInput$maf,
      XVariance = indInput$X_variance,
      mafCutoff = mafCutoff,
      pipCutoffToSkip = pipCutoffToSkipInd
    )
    individualData <- qcedIndividualToRegionData(indQc)
    # If eventFilters dropped any pre-QC contexts entirely, surface them as
    # NULL-Y entries so downstream extractContextsStudies() emits the
    # "Skipping follow-up analysis for individual traits ..." message.
    droppedCtx <- attr(regionData$individual_data, "filtered_out_contexts")
    if (!is.null(droppedCtx) && length(droppedCtx) > 0 && !is.null(individualData)) {
      for (ctx in droppedCtx) {
        individualData$X[[ctx]] <- NULL
        individualData$Y[[ctx]] <- list(NULL)[[1]]
      }
      # NULL inserts via [[ removed entries; re-insert as explicit NULL.
      for (ctx in droppedCtx) {
        if (!ctx %in% names(individualData$Y)) {
          individualData$Y <- c(individualData$Y, stats::setNames(list(NULL), ctx))
        }
        if (!ctx %in% names(individualData$X)) {
          individualData$X <- c(individualData$X, stats::setNames(list(NULL), ctx))
        }
      }
    }
  }

  sumstatData <- NULL
  rssInput <- regionDataToRssInput(regionData)
  if (isTRUE(rssInput$source_info$has_sumstat)) {
    sumstatQc <- summaryStatsQc(
      rssInput = rssInput$rss_input,
      ldData = rssInput$LD_data,
      keepIndel = keepIndel,
      pipCutoffToSkip = pipCutoffToSkipSumstat,
      zMismatchQc = zMismatchQc,
      impute = impute,
      imputeOpts = imputeOpts
    )
    sumstatData <- qcedSumstatToRegionData(sumstatQc)
  }

  list(individual_data = individualData, sumstat_data = sumstatData)
}

#' Run reusable individual-level QC
#'
#' @param X Genotype matrix or named list of genotype matrices.
#' @param Y Phenotype vector/matrix or named list of phenotype matrices.
#' @param maf Optional MAF vector or named list.
#' @param XVariance Optional variant variance vector or named list.
#' @param missingRateThresh Maximum missing genotype rate.
#' @param mafCutoff Minimum MAF cutoff.
#' @param xvarCutoff Minimum genotype variance cutoff.
#' @param ldReferenceMetaFile Optional LD reference metadata file.
#' @param keepIndel Whether indel variants are kept during LD-reference
#'   filtering.
#' @param pipCutoffToSkip Optional single-effect PIP cutoff.
#' @return A named list of cleaned context-level \code{X}/\code{Y} records, or
#'   one cleaned record for matrix inputs.
#' @noRd
qcIndividualData <- function(X, Y, maf = NULL, XVariance = NULL,
                             missingRateThresh = NULL,
                             mafCutoff = 0.0005,
                             xvarCutoff = 0,
                             ldReferenceMetaFile = NULL,
                             keepIndel = TRUE,
                             pipCutoffToSkip = 0,
                             context = NULL) {
  qcOne <- function(X, Y, maf = NULL, XVariance = NULL, context = NULL,
                    pipCutoffToSkip = 0) {
    if (is.null(X) || is.null(Y)) return(NULL)
    if (is.null(colnames(X))) stop("X must have variant colnames for individual-level QC.")
    if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
    if (is.null(colnames(Y))) colnames(Y) <- .cbDefault(context, paste0("outcome", seq_len(ncol(Y))))
    if (!is.null(context)) colnames(Y) <- paste0(context, "_", colnames(Y))

    message("QC track: starting individual-level QC for ", .cbDefault(context, "individual data"), ".")
    originalVariants <- colnames(X)
    if (!is.null(names(maf)) || length(maf) == length(originalVariants)) {
      if (is.null(names(maf))) names(maf) <- originalVariants
    }
    if (!is.null(names(XVariance)) || length(XVariance) == length(originalVariants)) {
      if (is.null(names(XVariance))) names(XVariance) <- originalVariants
    }
    if (!is.null(ldReferenceMetaFile)) {
      referenceFilter <- filterVariantsByLdReference(originalVariants, ldReferenceMetaFile,
                                                     keepIndel = keepIndel)
      X <- X[, referenceFilter$data, drop = FALSE]
      if (!is.null(names(maf))) maf <- maf[colnames(X)]
      if (!is.null(names(XVariance))) XVariance <- XVariance[colnames(X)]
    }
    X <- filterX(X, missingRateThresh = missingRateThresh,
                 mafThresh = mafCutoff, varThresh = xvarCutoff,
                 maf = maf, xVariance = XVariance)
    if (!is.null(names(maf))) maf <- maf[colnames(X)]
    if (!is.null(names(XVariance))) XVariance <- XVariance[colnames(X)]

    if (!is.null(pipCutoffToSkip) && pipCutoffToSkip != 0) {
      cutoff <- pipCutoffToSkip
      if (cutoff < 0) cutoff <- 3 / ncol(X)
      keepY <- logical(ncol(Y))
      for (j in seq_len(ncol(Y))) {
        observed <- !is.na(Y[, j])
        if (sum(observed) < 2) next
        pip <- susieR::susie(X[observed, , drop = FALSE], Y[observed, j],
                             L = 1, max_iter = 100)$pip
        keepY[j] <- any(pip > cutoff)
      }
      if (!any(keepY)) {
        message("QC track: skipping individual context ", context,
                ". No outcomes passed PIP threshold ", cutoff, ".")
        return(NULL)
      }
      Y <- Y[, keepY, drop = FALSE]
    }
    message("QC track: retained ", ncol(X), " variants and ", ncol(Y),
            " outcome(s) for individual context ", .cbDefault(context, "input"), ".")
    list(X = X, Y = Y, maf = maf, X_variance = XVariance)
  }

  if (is.list(X) && !is.matrix(X) && !is.data.frame(X)) {
    X <- .cbAsNamedList(X, "individual")
    Y <- .cbAsNamedList(Y, "individual")
    contexts <- intersect(names(X), names(Y))
    if (length(contexts) == 0) stop("No matched X/Y contexts for individual-level QC.")
    cutoffs <- .cbNamedCutoff(pipCutoffToSkip, contexts, "pipCutoffToSkipInd")
    out <- list()
    for (context in contexts) {
      res <- qcOne(
        X[[context]], Y[[context]],
        maf = .cbListValue(maf, context),
        XVariance = .cbListValue(XVariance, context),
        context = context,
        pipCutoffToSkip = cutoffs[[context]]
      )
      if (!is.null(res)) out[[context]] <- res
    }
    return(out)
  }
  qcOne(
    X = X, Y = Y, maf = maf, XVariance = XVariance,
    pipCutoffToSkip = pipCutoffToSkip,
    context = context
  )
}


##### Generic ColocBoost helper functions #####

.cbDefault <- function(x, y) if (is.null(x)) y else x

.cbMergeArgs <- function(x, y) {
  for (nm in names(y)) {
    if (!is.null(y[[nm]])) x[[nm]] <- y[[nm]]
  }
  x
}

.cbAsNamedList <- function(x, defaultName) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.matrix(x) && !is.data.frame(x)) return(x)
  stats::setNames(list(x), defaultName)
}

.cbListValue <- function(x, name, default = NULL) {
  if (is.null(x)) return(default)
  if (is.list(x) && !is.matrix(x) && !is.data.frame(x)) {
    if (name %in% names(x)) return(x[[name]])
    return(default)
  }
  x
}

.cbNamedCutoff <- function(x, namesToFill, argName) {
  if (length(x) == 1 && is.null(names(x))) {
    return(stats::setNames(rep(x, length(namesToFill)), namesToFill))
  }
  if (!is.null(names(x))) {
    missing <- setdiff(namesToFill, names(x))
    if (length(missing) > 0) x <- c(x, stats::setNames(rep(0, length(missing)), missing))
    return(x[namesToFill])
  }
  if (length(x) == length(namesToFill)) {
    return(stats::setNames(x, namesToFill))
  }
  stop(argName, " must be a scalar, named vector, or match the number of inputs.")
}

.cbYNcol <- function(y) {
  if (is.null(y)) return(0L)
  if (is.null(dim(y))) return(as.integer(length(y) > 0))
  ncol(y)
}

##### Outcome helper functions #####

.cbColocboostOutcomeNames <- function(args, preferSupplied = TRUE) {
  if (isTRUE(preferSupplied) && !is.null(args$outcome_names)) {
    return(as.character(args$outcome_names))
  }
  Y <- args$Y
  yOutcomes <- character()
  if (!is.null(Y)) {
    if (is.data.frame(Y)) Y <- as.matrix(Y)
    if (is.matrix(Y)) {
      yOutcomes <- colnames(Y)
      if (is.null(yOutcomes)) yOutcomes <- paste0("Y", seq_len(ncol(Y)))
    } else if (is.atomic(Y) && !is.list(Y)) {
      yOutcomes <- "Y1"
    } else {
      yNames <- names(Y)
      yCols <- unlist(Map(function(y, nm) {
        if (is.null(y)) return(character())
        if (is.data.frame(y)) y <- as.matrix(y)
        if (!is.null(dim(y))) {
          cn <- colnames(y)
          if (is.null(cn) || any(is.na(cn) | cn == "")) {
            if (!is.null(nm) && !is.na(nm) && nm != "" && ncol(y) == 1) return(nm)
            return(paste0(.cbDefault(nm, "Y"), "_", seq_len(ncol(y))))
          }
          return(as.character(cn))
        }
        if (!is.null(nm) && !is.na(nm) && nm != "") return(nm)
        "Y"
      }, Y, .cbDefault(yNames, rep("", length(Y)))), use.names = FALSE)
      if (length(yCols) != length(Y) || anyDuplicated(yCols) == 0) {
        yOutcomes <- yCols
      } else {
        if (is.null(yNames) || any(is.na(yNames) | yNames == "")) {
          yNames <- paste0("Y", seq_along(Y))
        }
        yOutcomes <- yNames
      }
    }
  }

  sumstat <- args$sumstat
  effectEst <- args$effect_est
  sumstatOutcomes <- character()
  if (!is.null(sumstat)) {
    if (is.data.frame(sumstat)) {
      sumstatOutcomes <- "sumstat1"
    } else {
      ssNames <- names(sumstat)
      if (is.null(ssNames) || any(is.na(ssNames) | ssNames == "")) {
        ssNames <- paste0("sumstat", seq_along(sumstat))
      }
      sumstatOutcomes <- as.character(ssNames)
    }
  } else if (!is.null(effectEst)) {
    effectEst <- as.matrix(effectEst)
    effectNames <- colnames(effectEst)
    if (is.null(effectNames)) effectNames <- paste0("sumstat", seq_len(ncol(effectEst)))
    sumstatOutcomes <- as.character(effectNames)
  }
  c(as.character(yOutcomes), sumstatOutcomes)
}

.cbResolveQcOutcomeNames <- function(preQcDataOutcomes,
                                     preQcDisplayOutcomes,
                                     postQcDataOutcomes) {
  if (length(postQcDataOutcomes) == 0) return(character())
  if (length(preQcDataOutcomes) == length(preQcDisplayOutcomes) &&
      length(preQcDataOutcomes) > 0 &&
      !anyDuplicated(preQcDataOutcomes)) {
    idx <- match(postQcDataOutcomes, preQcDataOutcomes)
    if (all(!is.na(idx))) {
      return(preQcDisplayOutcomes[idx])
    }
  }
  postQcDataOutcomes
}

.cbRemapFocalOutcomeIdx <- function(focalOutcomeIdx,
                                    preQcDataOutcomes,
                                    preQcDisplayOutcomes,
                                    postQcDataOutcomes,
                                    postQcDisplayOutcomes) {
  if (is.null(focalOutcomeIdx)) return(NULL)
  if (length(focalOutcomeIdx) != 1 || is.na(focalOutcomeIdx)) {
    warning("focal_outcome_idx must be a single non-missing index. Passing it through unchanged.")
    return(focalOutcomeIdx)
  }
  focalOutcomeIdx <- as.integer(focalOutcomeIdx)
  if (length(postQcDataOutcomes) == 0) return(NULL)
  if (focalOutcomeIdx < 1 ||
      focalOutcomeIdx > max(length(preQcDisplayOutcomes), length(preQcDataOutcomes))) {
    warning("focal_outcome_idx is outside the pre-QC outcome range. Passing it through unchanged.")
    return(focalOutcomeIdx)
  }

  focalDisplay <- if (focalOutcomeIdx <= length(preQcDisplayOutcomes)) {
    preQcDisplayOutcomes[[focalOutcomeIdx]]
  } else {
    NULL
  }
  focalData <- if (focalOutcomeIdx <= length(preQcDataOutcomes)) {
    preQcDataOutcomes[[focalOutcomeIdx]]
  } else {
    NULL
  }

  candidates <- integer()
  for (needle in unique(c(focalDisplay, focalData))) {
    if (is.null(needle) || is.na(needle) || needle == "") next
    candidates <- c(
      candidates,
      which(postQcDisplayOutcomes == needle),
      which(postQcDataOutcomes == needle)
    )
    suffixPattern <- paste0("_", needle)
    candidates <- c(
      candidates,
      which(endsWith(postQcDisplayOutcomes, suffixPattern)),
      which(endsWith(postQcDataOutcomes, suffixPattern))
    )
  }
  candidates <- unique(candidates)
  if (length(candidates) > 0) {
    return(candidates[[1]])
  }
  if (length(preQcDataOutcomes) == length(postQcDataOutcomes) ||
      length(preQcDisplayOutcomes) == length(postQcDisplayOutcomes)) {
    return(focalOutcomeIdx)
  }

  focalLabel <- .cbDefault(focalDisplay, .cbDefault(focalData, focalOutcomeIdx))
  warning("focal_outcome_idx refers to outcome ", focalLabel,
          ", which is not present after QC. Setting focal_outcome_idx to NULL.")
  NULL
}

##### Individual-level ColocBoost helper functions #####

.cbIndividualQcInputFromColocboost <- function(X, Y, dict_YX = NULL) {
  bindYForQc <- function(YList, yNames = NULL) {
    if (is.null(yNames)) yNames <- names(YList)
    if (is.null(yNames)) yNames <- rep("", length(YList))
    mats <- Map(function(y, nm) {
      if (is.null(dim(y))) y <- matrix(y, ncol = 1)
      if (is.null(colnames(y))) {
        colnames(y) <- if (!is.null(nm) && !is.na(nm) && nm != "") nm else paste0("Y", seq_len(ncol(y)))
      }
      y
    }, YList, yNames)
    do.call(cbind, mats)
  }

  XList <- .cbAsNamedList(X, "X1")
  YList <- .cbAsNamedList(Y, "Y1")
  matched <- intersect(names(XList), names(YList))
  if (length(matched) > 0 && is.null(dict_YX)) {
    return(list(X = XList[matched], Y = YList[matched]))
  }

  if (!is.null(dict_YX)) {
    dict <- as.matrix(dict_YX)
    if (ncol(dict) < 2) stop("dict_YX must have at least two columns.")
    XQc <- list()
    YQc <- list()
    for (xIdx in unique(dict[, 2])) {
      if (is.na(xIdx) || xIdx < 1 || xIdx > length(XList)) next
      yIdx <- dict[dict[, 2] == xIdx, 1]
      yIdx <- yIdx[!is.na(yIdx) & yIdx >= 1 & yIdx <= length(YList)]
      if (length(yIdx) == 0) next
      context <- names(XList)[xIdx]
      if (is.null(context) || is.na(context) || context == "") context <- paste0("X", xIdx)
      XQc[[context]] <- XList[[xIdx]]
      YQc[[context]] <- bindYForQc(YList[yIdx], names(YList)[yIdx])
    }
    if (length(XQc) > 0) return(list(X = XQc, Y = YQc))
  }

  if (length(XList) == 1 && length(YList) > 0) {
    context <- names(XList)[1]
    if (is.null(context) || is.na(context) || context == "") context <- "X1"
    return(list(
      X = stats::setNames(list(XList[[1]]), context),
      Y = stats::setNames(list(bindYForQc(YList, names(YList))), context)
    ))
  }

  if (length(XList) == length(YList)) {
    if (is.null(names(XList))) names(XList) <- paste0("X", seq_along(XList))
    names(YList) <- names(XList)
    return(list(X = XList, Y = YList))
  }

  list(X = XList, Y = YList)
}

.cbFormatIndividual <- function(ind) {
  if (length(ind) == 0) return(list())
  if (!is.null(ind$X) && !is.null(ind$Y)) {
    ind <- list(individual = ind)
  }
  ind <- Filter(function(x) {
    !is.null(x$X) && .cbYNcol(x$Y) > 0
  }, ind)
  if (length(ind) == 0) return(list())
  X <- lapply(ind, `[[`, "X")
  Y <- lapply(ind, `[[`, "Y")
  uniqueX <- list()
  XMatch <- integer(length(X))
  for (i in seq_along(X)) {
    matched <- names(uniqueX)[vapply(uniqueX, identical, logical(1), X[[i]])]
    if (length(matched) > 0) {
      XMatch[[i]] <- match(matched[[1]], names(uniqueX))
    } else {
      uniqueX[[names(ind)[i]]] <- X[[i]]
      XMatch[[i]] <- length(uniqueX)
    }
  }
  YSplit <- list()
  dict <- matrix(nrow = 0, ncol = 2)
  yIndex <- 0L
  allOutcomeNames <- unlist(lapply(Y, colnames), use.names = FALSE)
  duplicatedOutcomeNames <- unique(allOutcomeNames[
    duplicated(allOutcomeNames) | duplicated(allOutcomeNames, fromLast = TRUE)
  ])
  for (i in seq_along(Y)) {
    y <- Y[[i]]
    if (is.null(dim(y))) y <- matrix(y, ncol = 1)
    for (j in seq_len(ncol(y))) {
      yIndex <- yIndex + 1L
      outcome <- colnames(y)[j]
      if (is.null(outcome) || is.na(outcome) || outcome == "") {
        outcome <- paste0("outcome", yIndex)
      }
      context <- names(ind)[i]
      if (outcome %in% duplicatedOutcomeNames &&
          !is.null(context) && !is.na(context) && context != "") {
        outcome <- paste0(context, "_", outcome)
      }
      if (outcome %in% names(YSplit)) {
        outcome <- make.unique(c(names(YSplit), outcome))[length(YSplit) + 1]
      }
      YSplit[[outcome]] <- y[, j, drop = FALSE]
      dict <- rbind(dict, c(yIndex, XMatch[[i]]))
    }
  }
  colnames(dict) <- c("Y", "X")
  list(X = uniqueX, Y = YSplit, dict_YX = dict, outcome_names = names(YSplit))
}

##### Summary-statistic ColocBoost helper functions #####

.cbSumstatQcInputFromColocboost <- function(sumstat, LD, X_ref, dict_sumstatLD,
                                            LD_reference_info = NULL,
                                            variant_convention = c("A2_A1", "A1_A2")) {
  isLdData <- function(x) {
    is(x, "LdData")
  }
  asReferenceInfoList <- function(x) {
    if (is.null(x)) return(NULL)
    if (isLdData(x) || is.data.frame(x) || is.character(x)) {
      return(list(LD = x))
    }
    if (is.list(x) && length(x) > 0) return(x)
    stop("LD_reference_info must be a .bim/.pvar path, data.frame, load_LD_matrix() result, or a list of these.")
  }
  referenceInfoForIndex <- function(referenceInfo, index) {
    if (is.null(referenceInfo)) return(NULL)
    if (length(referenceInfo) == 1) return(referenceInfo[[1]])
    referenceInfo[[min(index, length(referenceInfo))]]
  }
  validateSumstatForQc <- function(sumstat, study) {
    if (is.null(sumstat) || !is.data.frame(sumstat)) {
      return(paste0("Summary-statistic QC for study ", study,
                    " requires each sumstat input to be a data.frame."))
    }
    missingCols <- setdiff(c("variant", "z"), colnames(sumstat))
    if (length(missingCols) > 0) {
      return(paste0("Summary-statistic QC for study ", study,
                    " requires sumstat columns: ", paste(missingCols, collapse = ", "), "."))
    }
    NULL
  }
  ldVariantNames <- function(ld) {
    if (!is.matrix(ld)) return(NULL)
    if (nrow(ld) == ncol(ld)) rownames(ld) else colnames(ld)
  }

  variant_convention <- match.arg(variant_convention)
  sumstat <- .cbAsNamedList(sumstat, "sumstat")
  usingXRef <- is.null(LD) && !is.null(X_ref)
  ldSource <- .cbAsNamedList(if (!is.null(LD)) LD else X_ref, "LD")
  referenceInfo <- asReferenceInfoList(LD_reference_info)
  if (is.null(dict_sumstatLD)) {
    dict_sumstatLD <- cbind(seq_along(sumstat), pmin(seq_along(sumstat), length(ldSource)))
  }
  rssInput <- list()
  LD_data <- list()
  skipReasons <- character()
  for (i in seq_along(sumstat)) {
    study <- names(sumstat)[i]
    ldIndex <- dict_sumstatLD[i, 2]
    validationReason <- validateSumstatForQc(sumstat[[i]], study)
    if (!is.null(validationReason)) {
      skipReasons <- c(skipReasons, validationReason)
      next
    }
    ldMat <- ldSource[[ldIndex]]
    refInfo <- referenceInfoForIndex(referenceInfo, ldIndex)
    if (is.null(refInfo)) {
      message("QC track: checking LD/X_ref variant names for summary-stat study ", study, ".")
      refIds <- ldVariantNames(ldMat)
      if (is.null(refIds)) {
        skipReasons <- c(
          skipReasons,
          paste0("Summary-statistic QC for study ", study,
                 " requires LD row/column names or X_ref column names parseable as genomic variant IDs; ",
                 "provide LD_reference_info for QC. Skipping summary-statistic QC for this study.")
        )
        next
      }
      ldData <- .cbMakeLdData(
        ldMat,
        isGenotype = usingXRef,
        variantConvention = variant_convention
      )
      if (is.null(ldData)) {
        skipReasons <- c(
          skipReasons,
          paste0("Summary-statistic QC for study ", study,
                 " could not parse LD/X_ref names as genomic variant IDs; ",
                 "provide LD_reference_info for QC. Skipping summary-statistic QC for this study.")
        )
        next
      }
      message("QC track: LD/X_ref names are parseable for summary-stat study ", study, ".")
    } else if (isLdData(refInfo)) {
      message("QC track: using supplied LD_reference_info LD data for summary-stat study ", study, ".")
      ldData <- refInfo
    } else {
      message("QC track: using supplied LD_reference_info variant metadata for summary-stat study ", study, ".")
      ldData <- .cbMakeLdData(
        ldMat,
        isGenotype = usingXRef,
        referenceInfo = refInfo,
        variantConvention = variant_convention
      )
    }
    parsed <- tryCatch(
      .cbParseVariantsForQc(sumstat[[i]]$variant, variant_convention),
      error = function(e) {
        skipReasons <<- c(
          skipReasons,
          paste0("Summary-statistic QC for study ", study,
                 " could not parse sumstat$variant as genomic variant IDs: ",
                 conditionMessage(e), " Skipping summary-statistic QC for this study.")
        )
        NULL
      }
    )
    if (is.null(parsed)) next
    variantId <- formatVariantId(parsed$chrom, parsed$pos, parsed$A2, parsed$A1)
    ss <- data.frame(parsed, z = sumstat[[i]]$z,
                     variant_id = variantId,
                     stringsAsFactors = FALSE)
    n <- if ("n" %in% colnames(sumstat[[i]])) unique(sumstat[[i]]$n)[1] else NULL
    varY <- if ("var_y" %in% colnames(sumstat[[i]])) unique(sumstat[[i]]$var_y)[1] else 1
    rssInput[[study]] <- list(sumstats = ss, n = n, var_y = varY)
    LD_data[[study]] <- ldData
  }
  list(rss_input = rssInput, LD_data = LD_data, skip_reasons = skipReasons)
}

.cbFormatSumstat <- function(sumstatQc) {
  validSumstatEntry <- function(ssDf, minVariants = 2) {
    if (is.null(ssDf) || !is.data.frame(ssDf)) return(FALSE)
    if (nrow(ssDf) < minVariants) return(FALSE)
    if (all(is.na(ssDf$z))) return(FALSE)
    if (all(ssDf$n <= 0 | is.na(ssDf$n))) return(FALSE)
    TRUE
  }
  filterValidSumstats <- function(sumstats, ldMat, minVariants = 2) {
    dedupeLd <- function(ldMat, studies) {
      uniqueLd <- list()
      LD_match <- character()
      for (study in studies) {
        ld <- ldMat[[study]]
        matched <- names(uniqueLd)[vapply(uniqueLd, identical, logical(1), ld)]
        if (length(matched) > 0) {
          LD_match <- c(LD_match, matched[[1]])
        } else {
          uniqueLd[[study]] <- ld
          LD_match <- c(LD_match, study)
        }
      }
      list(LD_mat = uniqueLd, LD_match = LD_match)
    }

    validIdx <- vapply(sumstats, validSumstatEntry, logical(1),
                       minVariants = minVariants)
    if (!any(validIdx)) {
      message("No valid summary statistic studies remain after validation.")
      return(NULL)
    }
    removed <- names(sumstats)[!validIdx]
    if (length(removed) > 0) {
      message("Removed invalid sumstat studies: ", paste(removed, collapse = ", "))
    }
    sumstats <- sumstats[validIdx]
    ldMat <- ldMat[validIdx]
    deduped <- dedupeLd(ldMat, names(sumstats))
    ldMat <- deduped$LD_mat
    LD_match <- deduped$LD_match
    dict_sumstatLD <- cbind(seq_along(sumstats), match(LD_match, names(ldMat)))
    list(sumstats = sumstats, LD_mat = ldMat, LD_match = LD_match,
         dict_sumstatLD = dict_sumstatLD)
  }
  if (length(sumstatQc) == 0) return(list())
  if (is(sumstatQc, "QcResult")) {
    sumstatQc <- list(sumstat = sumstatQc)
  }
  sumstat <- lapply(sumstatQc, function(x) {
    ss <- getRssInput(x)$sumstats
    variantId <- if ("variant_id" %in% colnames(ss)) {
      ss$variant_id
    } else {
      formatVariantId(ss$chrom, ss$pos, ss$A2, ss$A1)
    }
    data.frame(z = ss$z, n = getRssInput(x)$n,
               variant = normalizeVariantId(variantId),
               stringsAsFactors = FALSE)
  })
  ldMat <- lapply(sumstatQc, function(x) {
    ld <- getLdData(x)
    if (is.null(ld)) return(NULL)
    if (hasGenotypes(ld)) getGenotypes(ld) else getCorrelation(ld)
  })
  filtered <- filterValidSumstats(sumstat, ldMat)
  if (is.null(filtered)) return(list())
  c(
    list(
      sumstat = filtered$sumstats,
      dict_sumstatLD = filtered$dict_sumstatLD
    ),
    buildLdArgs(filtered$LD_mat)
  )
}


##### LD/reference QC helper functions #####

.cbMakeLdData <- function(ld, isGenotype = NULL,
                          referenceInfo = NULL,
                          variantConvention = c("A2_A1", "A1_A2")) {
  referenceInfoToRefPanel <- function(referenceInfo) {
    if (is.character(referenceInfo) && length(referenceInfo) == 1) {
      if (!file.exists(referenceInfo)) {
        stop("LD_reference_info file does not exist: ", referenceInfo)
      }
      referenceInfo <- readVariantMetadata(referenceInfo)
    }
    if (!is.data.frame(referenceInfo)) {
      stop("LD_reference_info must be a .bim/.pvar path or data.frame when it is not a load_LD_matrix() result.")
    }
    info <- as.data.frame(referenceInfo, stringsAsFactors = FALSE)
    names(info) <- sub("^#", "", names(info))
    names(info) <- sub("^chr$", "chrom", names(info), ignore.case = TRUE)
    names(info) <- sub("^CHROM$", "chrom", names(info))
    names(info) <- sub("^POS$", "pos", names(info))
    names(info) <- sub("^ID$", "id", names(info))
    names(info) <- sub("^REF$", "A2", names(info))
    names(info) <- sub("^ALT$", "A1", names(info))
    names(info) <- sub("^a0$", "A2", names(info))
    names(info) <- sub("^a1$", "A1", names(info))
    names(info) <- sub("^rsid$", "id", names(info), ignore.case = TRUE)

    if ("variants" %in% names(info) && !"variant_id" %in% names(info)) {
      info$variant_id <- normalizeVariantId(as.character(info$variants))
    }
    if (!"variant_id" %in% names(info)) {
      missingCols <- setdiff(c("chrom", "pos", "A2", "A1"), names(info))
      if (length(missingCols) > 0) {
        stop("LD_reference_info must contain variant_id or columns chrom, pos, A2, A1. Missing: ",
             paste(missingCols, collapse = ", "), ".")
      }
      info$variant_id <- normalizeVariantId(formatVariantId(info$chrom, info$pos, info$A2, info$A1))
    } else {
      info$variant_id <- normalizeVariantId(as.character(info$variant_id))
    }
    if (!all(c("chrom", "pos", "A2", "A1") %in% names(info))) {
      parsed <- parseVariantId(info$variant_id)
      info$chrom <- parsed$chrom
      info$pos <- parsed$pos
      info$A2 <- parsed$A2
      info$A1 <- parsed$A1
    }
    keepCols <- intersect(c("chrom", "id", "pos", "A2", "A1", "variant_id",
                            "allele_freq", "variance", "n_nomiss"), names(info))
    info[, keepCols, drop = FALSE]
  }

  alignRefPanelToLd <- function(refPanel, ldNames, nVariants) {
    if (!is.null(ldNames)) {
      matchIdx <- rep(NA_integer_, length(ldNames))
      if ("id" %in% names(refPanel)) {
        matchIdx <- match(ldNames, refPanel$id)
      }
      if (any(is.na(matchIdx))) {
        variantMatch <- match(normalizeVariantId(ldNames), refPanel$variant_id)
        matchIdx[is.na(matchIdx)] <- variantMatch[is.na(matchIdx)]
      }
      if (all(!is.na(matchIdx))) {
        refPanel <- refPanel[matchIdx, , drop = FALSE]
        rownames(refPanel) <- NULL
        return(refPanel)
      }
      if (length(ldNames) == nrow(refPanel)) {
        message("QC track: LD_reference_info could not be matched by LD names; using LD_reference_info row order.")
        return(refPanel[seq_along(ldNames), , drop = FALSE])
      }
      stop("LD_reference_info could not be matched to LD/X_ref names. ",
           "Provide an id/variant_id column matching LD names, or provide rows in LD matrix order.")
    }
    if (nrow(refPanel) != nVariants) {
      stop("LD_reference_info has ", nrow(refPanel), " variants but LD/X_ref has ",
           nVariants, " columns. Provide LD_reference_info in LD matrix order or with matching LD names.")
    }
    message("QC track: LD/X_ref names are missing; using LD_reference_info row order.")
    refPanel
  }

  variantConvention <- match.arg(variantConvention)
  if (is.null(isGenotype)) isGenotype <- is.matrix(ld) && nrow(ld) != ncol(ld)
  refIds <- if (is.matrix(ld) && nrow(ld) == ncol(ld)) rownames(ld) else colnames(ld)

  if (!is.null(referenceInfo)) {
    refPanel <- referenceInfoToRefPanel(referenceInfo)
    refPanel <- alignRefPanelToLd(refPanel, refIds, ncol(ld))
    variantIds <- refPanel$variant_id
    if (is.matrix(ld) && nrow(ld) == ncol(ld)) {
      rownames(ld) <- colnames(ld) <- variantIds
    } else if (is.matrix(ld)) {
      colnames(ld) <- variantIds
    }
    refPanel$chrom <- as.character(refPanel$chrom)
    variantsGr <- .refPanelToGranges(refPanel)
    corr <- if (isTRUE(isGenotype)) cor(ld) else ld
    bm <- .inferSingleLdBlockMetadata(refPanel)
    return(LdData(
      correlation = corr,
      variants = variantsGr,
      blockMetadata = bm
    ))
  }

  if (is.null(refIds)) {
    return(NULL)
  }

  parsed <- if (!is.null(refIds)) {
    tryCatch(
      .cbParseVariantsForQc(refIds, variantConvention),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  if (is.null(parsed)) return(NULL)
  variantIds <- if (!is.null(parsed)) formatVariantId(parsed$chrom, parsed$pos, parsed$A2, parsed$A1) else NULL
  if (!is.null(variantIds)) {
    if (is.matrix(ld) && nrow(ld) == ncol(ld)) {
      rownames(ld) <- colnames(ld) <- variantIds
    } else if (is.matrix(ld)) {
      colnames(ld) <- variantIds
    }
    parsed$variant_id <- variantIds
  }
  parsed$chrom <- as.character(parsed$chrom)
  variantsGr <- .refPanelToGranges(parsed)
  corr <- if (isTRUE(isGenotype)) cor(ld) else ld
  bm <- if (!is.null(parsed)) .inferSingleLdBlockMetadata(parsed) else data.frame()
  LdData(
    correlation = corr,
    variants = variantsGr,
    blockMetadata = bm
  )
}

.cbParseVariantsForQc <- function(ids, variantConvention = c("A2_A1", "A1_A2")) {
  variantConvention <- match.arg(variantConvention)
  parsed <- parseVariantId(ids)
  if (any(is.na(parsed$chrom)) || any(is.na(parsed$pos)) ||
      any(is.na(parsed$A1)) || any(is.na(parsed$A2))) {
    stop("QC requires variant IDs parseable as chr:pos:A2:A1 by default. ",
         "If the input uses chr:pos:A1:A2, set variant_convention = 'A1_A2'.")
  }
  if (identical(variantConvention, "A1_A2")) {
    parsed <- data.frame(chrom = parsed$chrom, pos = parsed$pos,
                         A2 = parsed$A1, A1 = parsed$A2,
                         stringsAsFactors = FALSE)
  }
  parsed
}
