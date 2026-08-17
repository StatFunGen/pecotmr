#' @title sLDSC Postprocessing Pipeline
#' @description Postprocess polyfun's per-trait sLDSC outputs (already loaded
#'   into an \code{\link{SldscData}} object) into a single results object with
#'   per-trait tau*, EnrichStat with back-solved jackknife SE, and a
#'   DerSimonian-Laird random-effects meta-analysis across traits. All file I/O
#'   is done up front by the reader functions (\code{\link{readSldscAnnot}},
#'   \code{\link{readSldscFrq}}, \code{\link{readSldscTrait}}); this pipeline is
#'   pure computation over the in-memory \code{SldscData}.
#' @param sldscData An \code{\link{SldscData}} object bundling the annotation
#'   table, the reference-panel allele frequencies, and the per-trait
#'   single/joint polyfun runs.
#' @param mafCutoff Numeric MAF cutoff applied via the object's frq table.
#'   Default \code{0.05}. Set to \code{0} to opt out (requires frq data when
#'   \code{> 0}).
#' @param targetCategories Optional character vector of target annotation names
#'   to retain. Auto-detected from the joint run (or first single run) when
#'   \code{NULL}.
#' @param targetLabels Optional display names, same length / order as
#'   \code{targetCategories}, applied to every output column / tau* block
#'   colname.
#' @return A list with \code{per_trait} (per-trait standardised tables), meta
#'   tables (\code{tauStar}, \code{enrichment}, \code{enrichstat}), and a
#'   \code{params} record of the call options.
#' @seealso \code{\link{SldscData}}, \code{\link{readSldscAnnot}},
#'   \code{\link{readSldscFrq}}, \code{\link{readSldscTrait}}
#' @importFrom stats median
#' @importFrom methods is
#' @include SldscData.R
#' @examples
#' mkRun <- function(cats) {
#'   n <- length(cats)
#'   list(categories = cats, tau = setNames(rep(1e-7, n), cats),
#'     tauSe = setNames(rep(3e-8, n), cats),
#'     enrichment = setNames(rep(2, n), cats),
#'     enrichmentSe = setNames(rep(0.4, n), cats),
#'     enrichmentP = setNames(rep(0.01, n), cats),
#'     propH2 = setNames(rep(0.2, n), cats),
#'     propSnps = setNames(rep(0.1, n), cats), h2g = 0.3,
#'     tauBlocks = matrix(1e-7, 10, n, dimnames = list(NULL, cats)),
#'     nBlocks = 10L)
#' annot <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   annot_A = c(1, 0, 1, 0, 1, 0), annot_B = c(2.1, 1.8, 2.5, 1.9, 2.3, 2))
#' frq <- data.frame(CHR = c(1, 1, 1, 2, 2, 2), SNP = paste0("rs", 1:6),
#'   MAF = rep(0.2, 6))
#' mkTrait <- function() {
#'   list(single = list(mkRun(c("annot_A_0", "baselineLD_0")),
#'     mkRun(c("annot_B_0", "baselineLD_0"))),
#'     joint = mkRun(c("annot_A_0", "annot_B_0", "baselineLD_0")))
#' }
#' traits <- setNames(list(mkTrait(), mkTrait()), c("traitX", "traitY"))
#' sd <- SldscData(annot = annot, frq = frq, traits = traits)
#' sldscPostprocessingPipeline(sldscData = sd)
#' }
#' @export
sldscPostprocessingPipeline <- function(
    sldscData,
    mafCutoff = 0.05,
    targetCategories = NULL,
    targetLabels = NULL
) {
    traitNames <- .sldscValidate(sldscData)
    ref <- .sldscComputeRefStats(sldscData, mafCutoff)
    tgt <- .sldscResolveTargets(sldscData, traitNames, targetCategories, ref)
    baselineCategories <- .sldscBaselineCategories(
        sldscData,
        traitNames,
        tgt$targetCategories
    )
    sub <- .sldscTargetStats(
        tgt$sdAnnotFull,
        tgt$isBinaryFull,
        tgt$targetCategories
    )
    nTraits <- length(traitNames)
    msg <- glue("[sldsc] Standardizing {nTraits} traits...")
    inform(msg)
    ctx <- list(
        sldscData = sldscData,
        MRef = ref$MRef,
        sdAnnot = sub$sdAnnot,
        isBinary = sub$isBinary,
        targetCategories = tgt$targetCategories
    )
    perTrait <- set_names(
        map(traitNames, .sldscStandardizeOneTrait, ctx = ctx),
        traitNames
    )
    meta <- .sldscMetaTables(perTrait, sub$isBinary, tgt$targetCategories)
    res <- .sldscAssembleResult(
        perTrait,
        meta,
        tgt$targetCategories,
        mafCutoff,
        ref$MRef,
        baselineCategories,
        traitNames
    )
    .sldscRelabel(res, targetLabels, tgt$targetCategories)
}

# Require an SldscData with at least one trait; returns the trait names.
# @noRd
.sldscValidate <- function(sldscData) {
    if (!is(sldscData, "SldscData")) {
        msg <- glue(
            "sldscPostprocessingPipeline: `sldscData` must be an ",
            "SldscData object."
        )
        abort(msg)
    }
    traitNames <- getTraitNames(sldscData)
    if (length(traitNames) == 0L) {
        abort("sldscPostprocessingPipeline: SldscData has no traits.")
    }
    traitNames
}

# M_ref, per-annotation sd, and binary/continuous flags. Target column names get
# polyfun's `_0` (file_idx=0) suffix so intersect() with a run's categories
# matches. Returns list(MRef, sdAnnotFull, isBinaryFull).
# @noRd
.sldscComputeRefStats <- function(sldscData, mafCutoff) {
    inform("[sldsc] Computing M_ref...")
    MRef <- computeSldscMRef(sldscData, mafCutoff = mafCutoff)
    msg <- glue("[sldsc]   M_ref = {MRef} (MAF cutoff {mafCutoff})")
    inform(msg)
    inform("[sldsc] Computing per-annotation sd...")
    sdAnnotFull <- computeSldscAnnotSd(sldscData, mafCutoff = mafCutoff)
    nSd <- length(sdAnnotFull)
    msg <- glue("[sldsc]   sd computed for {nSd} annotation columns")
    inform(msg)
    inform("[sldsc] Detecting binary vs continuous annotations...")
    isBinaryFull <- isBinarySldscAnnot(sldscData)
    names(sdAnnotFull) <- str_c(names(sdAnnotFull), "_0")
    names(isBinaryFull) <- str_c(names(isBinaryFull), "_0")
    list(MRef = MRef, sdAnnotFull = sdAnnotFull, isBinaryFull = isBinaryFull)
}

# Resolve the target categories: keep the user's set, else auto-detect from a
# pivot run (falling back to a positional rename when the `_0` names don't match
# polyfun's .results). Returns list(targetCategories, sdAnnotFull,
# isBinaryFull).
# @noRd
.sldscResolveTargets <- function(sldscData, traitNames, targetCategories, ref) {
    if (!is.null(targetCategories)) {
        return(list(
            targetCategories = targetCategories,
            sdAnnotFull = ref$sdAnnotFull,
            isBinaryFull = ref$isBinaryFull
        ))
    }
    pivotRun <- .sldscPivotRun(sldscData, traitNames[1])
    if (is.null(pivotRun)) {
        abort(
            "sldscPostprocessingPipeline: cannot auto-detect targetCategories."
        )
    }
    detected <- intersect(pivotRun$categories, names(ref$sdAnnotFull))
    out <- if (length(detected) == 0L) {
        .sldscFallbackRename(pivotRun, ref)
    } else {
        list(
            targetCategories = detected,
            sdAnnotFull = ref$sdAnnotFull,
            isBinaryFull = ref$isBinaryFull
        )
    }
    nDetected <- length(out$targetCategories)
    msg <- glue("[sldsc] Auto-detected {nDetected} target categories")
    inform(msg)
    out
}

# A representative run for auto-detection: the first trait's joint run, else its
# first single run.
# @noRd
.sldscPivotRun <- function(sldscData, trait1) {
    pivotRun <- getTraitRun(sldscData, trait1, "joint")
    if (is.null(pivotRun)) {
        pivotRun <- getTraitRun(sldscData, trait1, "single", 1L)
    }
    pivotRun
}

# Positional-rename fallback: trust polyfun's invariant that target categories
# occupy the first length(sdAnnotFull) rows of .results. Returns the renamed
# list(targetCategories, sdAnnotFull, isBinaryFull).
# @noRd
.sldscFallbackRename <- function(pivotRun, ref) {
    sdAnnotFull <- ref$sdAnnotFull
    isBinaryFull <- ref$isBinaryFull
    nTarget <- length(sdAnnotFull)
    nBaseline <- length(pivotRun$categories) - nTarget
    oldNames <- names(sdAnnotFull)
    targetCategories <- pivotRun$categories[seq_len(nTarget)]
    names(sdAnnotFull) <- targetCategories
    names(isBinaryFull) <- targetCategories
    .sldscFallbackMessage(
        pivotRun,
        oldNames,
        targetCategories,
        nTarget,
        nBaseline
    )
    list(
        targetCategories = targetCategories,
        sdAnnotFull = sdAnnotFull,
        isBinaryFull = isBinaryFull
    )
}

# @noRd
.sldscFallbackMessage <- function(
    pivotRun,
    oldNames,
    targetCategories,
    nTarget,
    nBaseline
) {
    baselinePreview <- if (nBaseline > 0L) {
        str_flatten(head(pivotRun$categories[-seq_len(nTarget)], 3), ", ")
    } else {
        "(none)"
    }
    oldStr <- str_flatten(oldNames, ", ")
    targetStr <- str_flatten(targetCategories, ", ")
    baselineMore <- if (nBaseline > 3L) ", ..." else ""
    msg <- glue(
        "[sldsc] sdAnnot/isBinary names did not match polyfun .results ",
        "categories;\n",
        "        falling back to positional rename (target = first ",
        "{nTarget} rows of .results)\n",
        "        target  ({nTarget}): {oldStr} -> {targetStr}\n",
        "        baseline ({nBaseline}): {baselinePreview}{baselineMore}",
        .trim = FALSE
    )
    inform(msg)
}

# Baseline categories = the first trait's joint categories minus the targets.
# @noRd
.sldscBaselineCategories <- function(sldscData, traitNames, targetCategories) {
    jointPivot <- getTraitRun(sldscData, traitNames[1], "joint")
    baselineCategories <- if (!is.null(jointPivot)) {
        setdiff(jointPivot$categories, targetCategories)
    } else {
        character(0)
    }
    .sldscBaselineMessage(baselineCategories)
    baselineCategories
}

# @noRd
.sldscBaselineMessage <- function(baselineCategories) {
    if (length(baselineCategories) == 0L) {
        msg <- glue(
            "[sldsc] No baseline annotations detected (no joint run on ",
            "the first trait)."
        )
        inform(msg)
        return(invisible(NULL))
    }
    msgTail <- if (length(baselineCategories) > 5) ", ..." else ""
    nBase <- length(baselineCategories)
    baseStr <- str_flatten(head(baselineCategories, 5), ", ")
    msg <- glue(
        "[sldsc] Detected {nBase} baseline annotations: {baseStr}{msgTail}"
    )
    inform(msg)
}

# Subset the per-annotation sd + binary flags to the target categories.
# @noRd
.sldscTargetStats <- function(sdAnnotFull, isBinaryFull, targetCategories) {
    isBinary <- if (length(isBinaryFull) > 0L) {
        isBinaryFull[targetCategories]
    } else {
        set_names(rep(FALSE, length(targetCategories)), targetCategories)
    }
    list(sdAnnot = sdAnnotFull[targetCategories], isBinary = isBinary)
}

# Standardize one trait (single + joint modes) into a per-trait record.
# @noRd
.sldscStandardizeOneTrait <- function(trait, ctx) {
    single <- .sldscTraitSingle(trait, ctx)
    joint <- .sldscTraitJoint(trait, ctx)
    summaryWide <- .sldscAssembleTraitSummary(
        single$singleDf,
        joint$jointDf,
        ctx$targetCategories,
        ctx$isBinary
    )
    list(
        summary = summaryWide,
        tau_star_blocks_single = single$blocksSingle,
        tau_star_blocks_joint = joint$blocksJoint,
        h2g = .sldscTraitH2g(joint$jointH2g, single$singleH2gs),
        nBlocks = joint$nBlocks
    )
}

# Trait h2g: the joint estimate when present, else the median of the single
# estimates.
# @noRd
.sldscTraitH2g <- function(jointH2g, singleH2gs) {
    if (!is.na(jointH2g)) {
        jointH2g
    } else if (length(singleH2gs) > 0L) {
        median(singleH2gs)
    } else {
        NA_real_
    }
}

# Single-mode standardization over the target categories (one run each). Returns
# list(singleDf, blocksSingle, singleH2gs).
# @noRd
.sldscTraitSingle <- function(trait, ctx) {
    singleRuns <- getTraitRun(ctx$sldscData, trait, "single")
    if (is.null(singleRuns)) {
        singleRuns <- list()
    }
    nRun <- min(length(ctx$targetCategories), length(singleRuns))
    stds <- compact(map(
        seq_len(nRun),
        .sldscStandardizeSingle,
        trait = trait,
        ctx = ctx
    ))
    .sldscCollectSingle(stds)
}

# Standardize the i-th single run (NULL + warning on failure).
# @noRd
.sldscStandardizeSingle <- function(i, trait, ctx) {
    catName <- ctx$targetCategories[i]
    std <- tryCatch(
        standardizeSldscTrait(
            ctx$sldscData,
            trait,
            mode = "single",
            idx = i,
            sdAnnot = ctx$sdAnnot[catName],
            MRef = ctx$MRef,
            targetCategories = catName
        ),
        error = function(e) {
            eMsg <- e$message
            msg <- glue(
                "[sldsc] Failed to standardize single {catName} for ",
                "{trait}: {eMsg}"
            )
            warn(msg)
            NULL
        }
    )
    if (is.null(std)) {
        return(NULL)
    }
    list(
        catName = catName,
        summary = std$summary,
        blocks = std$tau_star_blocks,
        h2g = std$h2g
    )
}

# Collect the per-category single results into (singleDf, blocksSingle,
# singleH2gs).
# @noRd
.sldscCollectSingle <- function(stds) {
    if (length(stds) == 0L) {
        return(list(
            singleDf = NULL,
            blocksSingle = NULL,
            singleH2gs = numeric(0)
        ))
    }
    catNames <- map_chr(stds, "catName")
    singleDf <- bind_rows(map(stds, "summary"))
    blocksList <- set_names(map(stds, "blocks"), catNames)
    list(
        singleDf = singleDf,
        blocksSingle = exec(cbind, !!!blocksList),
        singleH2gs = map_dbl(stds, "h2g")
    )
}

# Joint-mode standardization (when the trait has a joint run). Returns
# list(jointDf, blocksJoint, jointH2g, nBlocks).
# @noRd
.sldscTraitJoint <- function(trait, ctx) {
    empty <- list(
        jointDf = NULL,
        blocksJoint = NULL,
        jointH2g = NA_real_,
        nBlocks = NA_integer_
    )
    if (is.null(getTraitRun(ctx$sldscData, trait, "joint"))) {
        return(empty)
    }
    std <- tryCatch(
        standardizeSldscTrait(
            ctx$sldscData,
            trait,
            mode = "joint",
            sdAnnot = ctx$sdAnnot,
            MRef = ctx$MRef,
            targetCategories = ctx$targetCategories
        ),
        error = function(e) {
            eMsg <- e$message
            msg <- glue(
                "[sldsc] Failed to standardize joint for {trait}: {eMsg}"
            )
            warn(msg)
            NULL
        }
    )
    if (is.null(std)) {
        return(empty)
    }
    list(
        jointDf = std$summary,
        blocksJoint = std$tau_star_blocks,
        jointH2g = std$h2g,
        nBlocks = std$nBlocks
    )
}

# Random-effects meta across traits -> list(tauStar, enrichment, enrichstat).
# @noRd
.sldscMetaTables <- function(perTrait, isBinary, targetCategories) {
    inform("[sldsc] Running random-effects meta across traits...")
    ptViewSingle <- .sldscViewForMeta(perTrait, "single")
    ptViewJoint <- .sldscViewForMeta(perTrait, "joint")
    metaTauStarSingle <- .sldscBuildTable(
        "tauStar",
        ptViewSingle,
        "single",
        isBinary,
        targetCategories
    )
    metaTauStarJoint <- .sldscBuildTable(
        "tauStar",
        ptViewJoint,
        "joint",
        isBinary,
        targetCategories
    )
    metaESingle <- .sldscBuildTable(
        "enrichment",
        ptViewSingle,
        "single",
        isBinary,
        targetCategories
    )
    metaEsSingle <- .sldscBuildTable(
        "enrichstat",
        ptViewSingle,
        "single",
        isBinary,
        targetCategories
    )
    list(
        tauStar = .sldscCombineTauStar(metaTauStarSingle, metaTauStarJoint),
        enrichment = .sldscCombineEnrichment(metaESingle, metaEsSingle),
        enrichstat = metaEsSingle
    )
}

# Combine tauStar single + joint into one wide frame (joint aligned by target).
# @noRd
.sldscCombineTauStar <- function(metaTauStarSingle, metaTauStarJoint) {
    metaTauStar <- metaTauStarSingle
    ord <- match(metaTauStar$target, metaTauStarJoint$target)
    metaTauStar$jointMean <- metaTauStarJoint$jointMean[ord]
    metaTauStar$jointSe <- metaTauStarJoint$jointSe[ord]
    metaTauStar$jointP <- metaTauStarJoint$jointP[ord]
    metaTauStar
}

# Two-channel enrichment meta: effect/SE from E, p-value from EnrichStat.
# @noRd
.sldscCombineEnrichment <- function(metaESingle, metaEsSingle) {
    metaEnrichment <- metaESingle
    metaEnrichment$singleP <- metaEsSingle$singleP[match(
        metaEnrichment$target,
        metaEsSingle$target
    )]
    metaEnrichment
}

# Assemble the pipeline result (per_trait + meta + params).
# @noRd
.sldscAssembleResult <- function(
    perTrait,
    meta,
    targetCategories,
    mafCutoff,
    MRef,
    baselineCategories,
    traitNames
) {
    list(
        per_trait = perTrait,
        meta = meta,
        params = list(
            maf_cutoff = mafCutoff,
            M_ref = MRef,
            target_categories = targetCategories,
            n_baseline = length(baselineCategories),
            baseline_categories = baselineCategories,
            trait_names = traitNames
        )
    )
}

# Optionally relabel target categories to user-friendly display names (keeps the
# polyfun names when targetLabels is NULL).
# @noRd
.sldscRelabel <- function(res, targetLabels, targetCategories) {
    if (is.null(targetLabels)) {
        return(res)
    }
    targetLabels <- as.character(targetLabels)
    if (length(targetLabels) != length(targetCategories)) {
        nLabels <- length(targetLabels)
        nTargets <- length(targetCategories)
        targetStr <- str_flatten(targetCategories, ", ")
        msg <- glue(
            "sldscPostprocessingPipeline: targetLabels has length ",
            "{nLabels} but there are {nTargets} target categories ",
            "({targetStr})."
        )
        abort(msg)
    }
    relab <- set_names(targetLabels, targetCategories)
    res$per_trait <- .sldscRelabelPerTrait(res$per_trait, relab)
    res$meta <- .sldscRelabelMeta(res$meta, relab)
    res$params$target_categories_orig <- res$params$target_categories
    res$params$target_categories <- unname(relab[targetCategories])
    .sldscRelabelMessage(targetCategories, relab)
    res
}

# @noRd
.sldscRelabelPerTrait <- function(perTrait, relab) {
    for (t in names(perTrait)) {
        perTrait[[t]] <- .sldscRelabelOneTrait(perTrait[[t]], relab)
    }
    perTrait
}

# @noRd
.sldscRelabelOneTrait <- function(pt, relab) {
    if (!is.null(pt$summary) && is_in("target", names(pt$summary))) {
        pt$summary$target <- .sldscRelabVec(pt$summary$target, relab)
    }
    for (bn in c("tau_star_blocks_single", "tau_star_blocks_joint")) {
        b <- pt[[bn]]
        if (!is.null(b) && !is.null(colnames(b))) {
            colnames(pt[[bn]]) <- .sldscRelabVec(colnames(b), relab)
        }
    }
    pt
}

# @noRd
.sldscRelabelMeta <- function(meta, relab) {
    for (mn in names(meta)) {
        if (!is.null(meta[[mn]]) && is_in("target", names(meta[[mn]]))) {
            meta[[mn]]$target <- .sldscRelabVec(meta[[mn]]$target, relab)
        }
    }
    meta
}

# @noRd
.sldscRelabelMessage <- function(targetCategories, relab) {
    labels <- unname(relab[targetCategories])
    pairs <- glue("{targetCategories} -> {labels}")
    pairsStr <- str_flatten(pairs, ", ")
    msg <- glue("[sldsc] Relabeled target categories: {pairsStr}")
    inform(msg)
}

# Build a per-category meta table for one quantity/view, with `label`-prefixed
# mean/se/p columns and an isBinary flag.
# @noRd
.sldscBuildTable <- function(
    quantity,
    view,
    label,
    isBinary,
    targetCategories
) {
    rows <- list()
    for (category in targetCategories) {
        m <- metaSldscRandom(view, category, quantity)
        rows[[category]] <- tibble(
            target = category,
            isBinary = unname(isBinary[category]),
            mean = m$mean,
            se = m$se,
            p = m$p,
            nTraits = m$nTraits
        )
    }
    df <- bind_rows(rows)
    nmOld <- c("mean", "se", "p")
    nmNew <- str_c(label, str_to_upper(str_sub(nmOld, 1, 1)), str_sub(nmOld, 2))
    names(df)[is_in(names(df), nmOld)] <- nmNew
    df
}

# Relabel a target-category vector via the `relab` map, leaving unmapped values
# unchanged.
# @noRd
.sldscRelabVec <- function(x, relab) {
    y <- unname(relab[x])
    y[is.na(y)] <- x[is.na(y)]
    y
}
