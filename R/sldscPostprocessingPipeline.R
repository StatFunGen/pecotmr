#' @title sLDSC Postprocessing Pipeline
#' @description Postprocess polyfun's per-trait sLDSC outputs (one
#'   single-target run per target annotation, plus an optional joint
#'   run) into a single results object with per-trait tau*, EnrichStat
#'   with back-solved jackknife SE, and a DerSimonian-Laird random-
#'   effects meta-analysis across traits.
#' @param traitSinglePrefixes Named list of file prefixes for the
#'   single-target polyfun runs (one entry per trait; each value is a
#'   length-N character vector of `<dir>/<trait>` prefixes, one per
#'   target annotation).
#' @param traitJointPrefix Named list of file prefixes for the joint
#'   polyfun runs (one entry per trait; each value a `<dir>/<trait>`
#'   prefix into the joint LD-score dir). Pass an empty list to skip
#'   the joint branch.
#' @param targetAnnoDir Directory containing the target `.annot.gz`
#'   files used for sd_C and binary detection (typically the joint dir).
#' @param frqfileDir Optional directory of `.frq` files for the MAF
#'   cutoff. Pass \code{NULL} to skip MAF filtering.
#' @param plinkName File-name prefix of the PLINK reference panel
#'   (default \code{"ADSP_chr"}; combined per-chromosome as
#'   \code{paste0(plinkName, chrom)}).
#' @param mafCutoff Numeric MAF cutoff applied via the `.frq` files.
#'   Default \code{0.05}. Set to \code{0} to opt out.
#' @param targetCategories Optional character vector of target
#'   annotation names to retain. Auto-detected from the joint run when
#'   \code{NULL}.
#' @param targetLabels Optional display names, same length / order as
#'   \code{targetCategories}, applied to every output column / tau*
#'   block colname.
#' @return A list with \code{per_trait} (per-trait standardised tables),
#'   meta tables (\code{tau_star_meta}, \code{E_meta},
#'   \code{enrich_stat_meta}), and a \code{params} record of the call
#'   options.
#' @export
sldscPostprocessingPipeline <- function(traitSinglePrefixes,
                                        traitJointPrefix,
                                        targetAnnoDir,
                                        frqfileDir = NULL,
                                        plinkName = "ADSP_chr",
                                        mafCutoff = 0.05,
                                        targetCategories = NULL,
                                        targetLabels = NULL) {
  traitNames <- names(traitSinglePrefixes)
  if (is.null(traitNames))
    stop("sldscPostprocessingPipeline: traitSinglePrefixes must be a named list.")

  message("[sldsc] Computing M_ref...")
  MRef <- computeSldscMRef(targetAnnoDir = targetAnnoDir,
                           frqfileDir = frqfileDir,
                           plinkName = plinkName,
                           mafCutoff = mafCutoff)
  message(sprintf("[sldsc]   M_ref = %d (MAF cutoff %g)", MRef, mafCutoff))

  message("[sldsc] Computing per-annotation sd...")
  sdAnnotFull <- computeSldscAnnotSd(targetAnnoDir = targetAnnoDir,
                                     frqfileDir = frqfileDir,
                                     plinkName = plinkName,
                                     mafCutoff = mafCutoff)
  message(sprintf("[sldsc]   sd computed for %d annotation columns",
                  length(sdAnnotFull)))

  message("[sldsc] Detecting binary vs continuous annotations...")
  isBinaryFull <- isBinarySldscAnnot(targetAnnoDir = targetAnnoDir)

  # Polyfun renames target columns to `<col>_0` (file_idx=0 in --ref-ld-chr);
  # mirror that suffix so intersect() with pivotRun$categories matches.
  names(sdAnnotFull)  <- paste0(names(sdAnnotFull),  "_0")
  names(isBinaryFull) <- paste0(names(isBinaryFull), "_0")

  # Auto-detect target categories from a representative run.
  if (is.null(targetCategories)) {
    pivotRun <- NULL
    if (!is.null(traitJointPrefix) && length(traitJointPrefix) > 0) {
      jp <- traitJointPrefix[[1]]
      if (is.character(jp) && length(jp) == 1L && !is.na(jp) && nzchar(jp)) {
        pivotRun <- tryCatch(readSldscTrait(jp), error = function(e) NULL)
      }
    }
    if (is.null(pivotRun) &&
        length(traitSinglePrefixes) > 0L &&
        length(traitSinglePrefixes[[1]]) > 0L) {
      pivotRun <- tryCatch(readSldscTrait(traitSinglePrefixes[[1]][1]),
                           error = function(e) NULL)
    }
    if (is.null(pivotRun))
      stop("sldscPostprocessingPipeline: cannot auto-detect targetCategories.")
    targetCategories <- intersect(pivotRun$categories, names(sdAnnotFull))
    # Fallback when .annot.gz names + "_0" do not match polyfun's .results
    # Category. Happens when the pipeline ran with --snp-list, since
    # ldsc.py --l2 hardcodes the LD score column to "L2" (single annot) or
    # "<annot>L2" (joint), instead of preserving the .annot.gz names.
    # Trust polyfun's invariant that target categories occupy the first
    # length(sdAnnotFull) rows of .results, then rename positionally.
    if (length(targetCategories) == 0L) {
      nTarget    <- length(sdAnnotFull)
      nBaseline  <- length(pivotRun$categories) - nTarget
      oldNames   <- names(sdAnnotFull)
      targetCategories <- pivotRun$categories[seq_len(nTarget)]
      names(sdAnnotFull)  <- targetCategories
      names(isBinaryFull) <- targetCategories
      baselinePreview <- if (nBaseline > 0L)
        paste(head(pivotRun$categories[-seq_len(nTarget)], 3), collapse = ", ")
      else "(none)"
      message(sprintf(paste0(
        "[sldsc] sdAnnot/isBinary names did not match polyfun .results categories;\n",
        "        falling back to positional rename (target = first %d rows of .results)\n",
        "        target  (%d): %s -> %s\n",
        "        baseline (%d): %s%s"),
        nTarget,
        nTarget, paste(oldNames, collapse = ", "),
                 paste(targetCategories, collapse = ", "),
        nBaseline, baselinePreview,
        if (nBaseline > 3L) ", ..." else ""))
    }
    message(sprintf("[sldsc] Auto-detected %d target categories", length(targetCategories)))
  }

  baselineCategories <- character(0)
  if (!is.null(traitJointPrefix) && length(traitJointPrefix) > 0L) {
    jp <- traitJointPrefix[[1]]
    if (is.character(jp) && length(jp) == 1L && !is.na(jp) && nzchar(jp)) {
      pivot <- tryCatch(readSldscTrait(jp), error = function(e) NULL)
      if (!is.null(pivot))
        baselineCategories <- setdiff(pivot$categories, targetCategories)
    }
  }
  if (length(baselineCategories) > 0L) {
    msgHead <- paste(head(baselineCategories, 5), collapse = ", ")
    msgTail <- if (length(baselineCategories) > 5) ", ..." else ""
    message(sprintf("[sldsc] Detected %d baseline annotations: %s%s",
                    length(baselineCategories), msgHead, msgTail))
  } else {
    message("[sldsc] No baseline annotations detected (joint-run prefix missing or unreadable).")
  }

  sdAnnot <- sdAnnotFull[targetCategories]
  isBinary <- if (length(isBinaryFull) > 0L) isBinaryFull[targetCategories] else
              setNames(rep(FALSE, length(targetCategories)), targetCategories)

  message(sprintf("[sldsc] Standardizing %d traits...", length(traitNames)))
  perTrait <- list()

  for (trait in traitNames) {
    # ---- single-mode ----
    singleSummaries <- list()
    singleBlocks    <- list()
    singleH2gs      <- numeric(0)
    singPrefs <- traitSinglePrefixes[[trait]]
    for (i in seq_along(targetCategories)) {
      catName <- targetCategories[i]
      if (i > length(singPrefs)) break
      pref <- singPrefs[i]
      run  <- tryCatch(readSldscTrait(pref), error = function(e) {
        warning(sprintf("[sldsc] Failed to read single %s for %s: %s",
                        catName, trait, e$message)); NULL
      })
      if (is.null(run)) next
      std <- tryCatch(
        standardizeSldscTrait(run, sdAnnot[catName], MRef,
                              targetCategories = catName, mode = "single"),
        error = function(e) {
          warning(sprintf("[sldsc] Failed to standardize single %s for %s: %s",
                          catName, trait, e$message)); NULL
        })
      if (is.null(std)) next
      singleSummaries[[catName]] <- std$summary
      singleBlocks[[catName]]    <- std$tau_star_blocks
      singleH2gs                 <- c(singleH2gs, std$h2g)
    }
    singleDf <- if (length(singleSummaries) > 0L)
                  do.call(rbind, singleSummaries) else NULL
    if (!is.null(singleDf)) rownames(singleDf) <- NULL
    blocksSingle <- if (length(singleBlocks) > 0L) do.call(cbind, singleBlocks) else NULL

    # ---- joint-mode ----
    jointDf       <- NULL
    blocksJoint   <- NULL
    jointH2g      <- NA_real_
    nBlocksTrait  <- NA_integer_
    if (!is.null(traitJointPrefix) && trait %in% names(traitJointPrefix)) {
      jp <- traitJointPrefix[[trait]]
      if (is.character(jp) && length(jp) == 1L && !is.na(jp) && nzchar(jp)) {
        run <- tryCatch(readSldscTrait(jp), error = function(e) {
          warning(sprintf("[sldsc] Failed to read joint for %s: %s",
                          trait, e$message)); NULL
        })
        if (!is.null(run)) {
          std <- tryCatch(
            standardizeSldscTrait(run, sdAnnot, MRef,
                                  targetCategories = targetCategories,
                                  mode = "joint"),
            error = function(e) {
              warning(sprintf("[sldsc] Failed to standardize joint for %s: %s",
                              trait, e$message)); NULL
            })
          if (!is.null(std)) {
            jointDf       <- std$summary
            blocksJoint   <- std$tau_star_blocks
            jointH2g      <- std$h2g
            nBlocksTrait  <- std$nBlocks
          }
        }
      }
    }

    summaryWide <- .sldscAssembleTraitSummary(singleDf, jointDf,
                                              targetCategories, isBinary)
    perTrait[[trait]] <- list(
      summary                = summaryWide,
      tau_star_blocks_single = blocksSingle,
      tau_star_blocks_joint  = blocksJoint,
      h2g                    = if (!is.na(jointH2g)) jointH2g
                               else if (length(singleH2gs) > 0L) median(singleH2gs)
                               else NA_real_,
      nBlocks               = nBlocksTrait
    )
  }

  message("[sldsc] Running random-effects meta across traits...")
  ptViewSingle <- .sldscViewForMeta(perTrait, "single")
  ptViewJoint  <- .sldscViewForMeta(perTrait, "joint")

  buildTable <- function(quantity, view, label) {
    rows <- list()
    for (cat in targetCategories) {
      m <- metaSldscRandom(view, cat, quantity)
      rows[[cat]] <- data.frame(
        target    = cat,
        isBinary = unname(isBinary[cat]),
        mean      = m$mean,
        se        = m$se,
        p         = m$p,
        nTraits  = m$nTraits,
        stringsAsFactors = FALSE
      )
    }
    df <- do.call(rbind, rows)
    rownames(df) <- NULL
    nmOld <- c("mean", "se", "p")
    nmNew <- paste0(label, toupper(substring(nmOld, 1, 1)),
                    substring(nmOld, 2))
    names(df)[names(df) %in% nmOld] <- nmNew
    df
  }

  metaTauStarSingle <- buildTable("tauStar",   ptViewSingle, "single")
  metaTauStarJoint  <- buildTable("tauStar",   ptViewJoint,  "joint")
  metaESingle       <- buildTable("enrichment", ptViewSingle, "single")
  metaEsSingle      <- buildTable("enrichstat", ptViewSingle, "single")

  # Combine tauStar single + joint into one wide frame.
  metaTauStar <- metaTauStarSingle
  ord <- match(metaTauStar$target, metaTauStarJoint$target)
  metaTauStar$jointMean <- metaTauStarJoint$jointMean[ord]
  metaTauStar$jointSe   <- metaTauStarJoint$jointSe[ord]
  metaTauStar$jointP    <- metaTauStarJoint$jointP[ord]

  # Two-channel enrichment meta: effect/SE from E meta, p from EnrichStat meta.
  metaEnrichment <- metaESingle
  metaEnrichment$singleP <- metaEsSingle$singleP[match(metaEnrichment$target,
                                                          metaEsSingle$target)]

  # Pure EnrichStat meta (separate frame).
  metaEnrichstat <- metaEsSingle

  res <- list(
    per_trait = perTrait,
    meta = list(
      tauStar   = metaTauStar,
      enrichment = metaEnrichment,
      enrichstat = metaEnrichstat
    ),
    params = list(
      maf_cutoff          = mafCutoff,
      M_ref               = MRef,
      target_categories   = targetCategories,
      n_baseline          = length(baselineCategories),
      baseline_categories = baselineCategories,
      trait_names         = traitNames
    )
  )

  # Optional: relabel target categories to user-friendly display names.
  # If targetLabels is NULL we keep the original polyfun .results names.
  if (!is.null(targetLabels)) {
    targetLabels <- as.character(targetLabels)
    if (length(targetLabels) != length(targetCategories))
      stop(sprintf(paste0("sldscPostprocessingPipeline: targetLabels has length %d but ",
                          "there are %d target categories (%s)."),
                   length(targetLabels), length(targetCategories),
                   paste(targetCategories, collapse = ", ")))
    relab    <- setNames(targetLabels, targetCategories)
    relabVec <- function(x) { y <- unname(relab[x]); y[is.na(y)] <- x[is.na(y)]; y }

    for (t in names(res$per_trait)) {
      pt <- res$per_trait[[t]]
      if (!is.null(pt$summary) && "target" %in% names(pt$summary))
        res$per_trait[[t]]$summary$target <- relabVec(pt$summary$target)
      for (bn in c("tau_star_blocks_single", "tau_star_blocks_joint")) {
        b <- pt[[bn]]
        if (!is.null(b) && !is.null(colnames(b)))
          colnames(res$per_trait[[t]][[bn]]) <- relabVec(colnames(b))
      }
    }
    for (mn in names(res$meta)) {
      if (!is.null(res$meta[[mn]]) && "target" %in% names(res$meta[[mn]]))
        res$meta[[mn]]$target <- relabVec(res$meta[[mn]]$target)
    }
    res$params$target_categories_orig <- res$params$target_categories
    res$params$target_categories      <- unname(relab[targetCategories])
    message(sprintf("[sldsc] Relabeled target categories: %s",
                    paste(sprintf("%s -> %s", targetCategories,
                                  unname(relab[targetCategories])), collapse = ", ")))
  }

  res
}
