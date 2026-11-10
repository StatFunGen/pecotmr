#' @include MethodParam.R
NULL

#' @title Causal Inference Pipeline (TWAS-Z + Mendelian Randomization)
#' @description Per-region pipeline that pairs QTL-derived weight vectors
#'   (\code{\link{TwasWeights}} and/or a QTL \code{\link{QtlFineMappingResult}})
#'   with one or more GWAS studies (\code{\link{GwasSumStats}}) to produce
#'   per-tuple TWAS Z-scores and, when fine-mapping is supplied, Wald-ratio
#'   Mendelian Randomization estimates over the QTL credible sets.
#'
#'   Input combinations:
#'   \itemize{
#'     \item \code{twasWeights} alone (no \code{fineMappingResult}):
#'           TWAS Z only, no MR.
#'     \item \code{fineMappingResult} alone (no \code{twasWeights}):
#'           TWAS Z derived from the SuSiE-style coefficients in each
#'           row's top loci (\code{\link{topLoci}}); plus MR.
#'     \item both: TWAS Z computed from \code{twasWeights};
#'           MR computed from \code{fineMappingResult}.
#'   }
#'
#'   MR is computed only from \code{fineMappingResult}: for every gene and
#'   fine-mapping method it holds, against every GWAS study it overlaps,
#'   whichever TWAS methods were chosen and whether or not
#'   \code{twasWeights} holds the gene at all. A fine-mapping method among
#'   the TWAS weights (such as \code{susie}) is a TWAS method like any other,
#'   and MR never reads it.
#'
#' @section LD reference compatibility check: Per pair of rows scored
#'   together: a QTL input row (TwasWeights or QtlFineMappingResult) that
#'   carries an \code{ldReference} must be the same LD reference as the
#'   \code{gwasSumStats} rows its region overlaps, and the rows of one GWAS
#'   study that a QTL tuple overlaps must share one reference, since its TWAS
#'   z reads one LD matrix. Rows that are never scored together may sit on
#'   different references.
#'   The two must be the same kind of reference. Two genotype panels must
#'   share their samples and allele orientation; two \code{LdMatrix} objects
#'   must agree on the LD they share; two \code{GenotypeMixture} objects must
#'   average the same panels with the same weights. The two need NOT carry
#'   the same variants: QC-ing the two sides separately trims each reference
#'   to its own surviving variants, and a partial overlap only warns. A
#'   mismatch of any of those, or no shared variant at all, is a hard error. A
#'   QTL input row with no reference (the fit was learned from
#'   individual-level data) skips the validation.
#'
#' @section Output shape:
#' A long-format \code{GRanges} of two kinds of row, told apart by
#' \code{analysis}: a TWAS row per chosen \code{(qtlStudy, context, trait,
#' method)} weight set and GWAS study, positioned at the span of its
#' variants the GWAS holds, and an MR row per fine-mapping row and GWAS
#' study, positioned at the span of the fine-mapping row's variants. mcols:
#' \describe{
#'   \item{\code{qtlStudy}, \code{context}, \code{trait},
#'         \code{method}, \code{gwasStudy}}{Identity columns. On an MR row,
#'         \code{method} is the fine-mapping method.}
#'   \item{\code{jointStudies}, \code{jointContexts},
#'         \code{jointTraits}}{Present when a scored weight set came from a
#'         joint fit: the studies, contexts or traits the fit joined
#'         (\code{";"}-joined), as on its \code{twasWeights} or
#'         \code{fineMappingResult} row, \code{NA} for a weight set from no
#'         joint fit and on the \code{combined.*} rows. They tell apart the
#'         weight sets of one \code{(qtlStudy, context, trait, method)}
#'         fitted in different joint groups.}
#'   \item{\code{analysis}}{\code{"twas"} or \code{"mr"}.}
#'   \item{\code{twasZ}, \code{twasPval}}{TWAS Z and p-value; \code{NA} on
#'         MR rows.}
#'   \item{\code{waldRatio}, \code{waldRatioSe}, \code{mrPval}}{The
#'         IVW-aggregated Wald-ratio MR estimate, standard error, and
#'         p-value; \code{NA} on TWAS rows, and on an MR row whose estimate
#'         \code{mrParam}'s \code{pvalCutoff} withheld or that had no
#'         instrument.}
#'   \item{\code{nIV}}{Number of instrumental variables used in the
#'         MR aggregation.}
#'   \item{\code{Q}, \code{I2}, \code{nCs}}{The credible-set heterogeneity
#'         and the number of credible sets, from \code{method = "csAware"}
#'         MR; \code{NA} otherwise.}
#' }
#'
#' @param gwasSumStats A \code{\link{GwasSumStats}} object. Must be QC'd
#'   (\code{length(qcOptions(x)) > 0L}). A study may span many rows (one per
#'   chromosome or LD block); each QTL tuple is scored once per study, over the
#'   statistics of the rows its window overlaps. A tuple's variants outside
#'   every loaded block of a study it overlaps are an error where the blocks'
#'   ranges were given (\code{regionGiven}, see
#'   \code{\linkS4class{GwasSumStats}}), since they lie in a block that was not
#'   loaded, and a warning for more than one in a row past a range computed
#'   from the GWAS variants.
#' @param twasWeights Optional \code{\link{TwasWeights}} carrying per-(study,
#'   context, trait, method) weights. When supplied, drives the TWAS-Z
#'   computation.
#' @param fineMappingResult Optional \code{\link{QtlFineMappingResult}}. When
#'   supplied, drives the MR computation, for every row it holds, and (when
#'   \code{twasWeights = NULL}) the TWAS-Z weights via the SuSiE-style
#'   coefficients on each entry's \code{topLoci}.
#' @param weightSelectionParam Which TWAS methods' weights are scored for each
#'   gene, built with \code{\link{WeightSelectionParam}}: by default each
#'   gene's best method by cross-validated r-squared, or \code{"all"}, or the
#'   methods named; and the \code{cutoff} / \code{pvalCutoff} gates every
#'   choice is bounded by. The choice reads the weights' \code{cvResult}
#'   (\code{twasWeights}, else \code{fineMappingResult}), never a GWAS study.
#'   Ports the legacy \code{twas_pipeline} \code{pick_best_model}. The
#'   metrics describe each fit on all of its variants, so a choice that reads
#'   them warns about a row narrowed since its fit (condition class
#'   \code{pecotmrFitCoverage}).
#' @param mrParam How Mendelian randomization is run, built with
#'   \code{\link{MrParam}}: \code{method} (\code{"ivwPerVariant"}, the
#'   default, or \code{"csAware"}), the PIP threshold that method reads
#'   (\code{pipCutoff} or \code{cpipCutoff}), and \code{pvalCutoff}: a gene's
#'   MR against a GWAS study is estimated only when one of the gene's TWAS
#'   p-values for that study is below it (default \code{1}, no gate). A gene
#'   with no TWAS row for the study is not gated.
#' @param combineMethods Optional character vector forwarded to
#'   \code{\link{combinePValues}} for cross-method combination of the TWAS
#'   rows of each \code{(qtlStudy, context, trait, gwasStudy)} group.
#'   \code{NULL} (default) skips combination.
#' @param alleleFlip Logical, default \code{TRUE}. When TRUE, match QTL variants
#'   to the GWAS by (chrom, pos) with ref/alt swaps recognized and the exposure
#'   effect / weight sign-flipped accordingly; when FALSE, match on exact
#'   alleles only, so a ref/alt swap is treated as a distinct variant.
#' @return A \code{GRanges} as described above, its \code{seqinfo()} naming
#'   the genome build the inputs name (inputs naming different builds are
#'   refused).
#' @examples
#' data(qtlDatasetExample)
#' data(gwasSumStatsS4Example)
#' tw <- twasWeightsPipeline(qtlDatasetExample,
#'   methods = "lasso", cisWindow = 1e6)
#' causalInferencePipeline(
#'   gwasSumStats = gwasSumStatsS4Example, twasWeights = tw)
#' @export
causalInferencePipeline <- function(
    gwasSumStats,
    twasWeights = NULL,
    fineMappingResult = NULL,
    weightSelectionParam = WeightSelectionParam(),
    mrParam = MrParam(),
    combineMethods = NULL,
    alleleFlip = TRUE
) {
    .assertMethodParam(
        weightSelectionParam,
        "WeightSelectionParam"
    )
    .assertMethodParam(mrParam, "MrParam")
    build <- .inputsBuild(
        list(gwasSumStats, twasWeights, fineMappingResult),
        "causalInferencePipeline"
    )
    res <- .cipRun(
        gwasSumStats = gwasSumStats,
        twasWeights = twasWeights,
        fineMappingResult = fineMappingResult,
        combineMethods = combineMethods,
        weightSelectionParam = weightSelectionParam,
        alleleFlip = alleleFlip,
        mrParam = mrParam
    )
    .withGenome(res, build)
}

# Orchestrate the causal-inference pipeline: validate; choose each gene's TWAS
# weights (twasWeights, else the fine-mapping result); score them against
# every GWAS study they overlap; run MR for every fine-mapping row; finalize.
# @noRd
.cipRun <- function(
    gwasSumStats,
    twasWeights,
    fineMappingResult,
    combineMethods,
    weightSelectionParam,
    alleleFlip,
    mrParam
) {
    .cipValidateInputs(gwasSumStats, twasWeights, fineMappingResult)
    .cipCheckLdRefs(gwasSumStats, twasWeights, fineMappingResult)
    weightSource <- twasWeights %||% fineMappingResult
    cells <- .cipChosenCells(weightSource, weightSelectionParam)
    twasTuples <- compact(map(
        seq_len(nrow(cells)),
        .cipTwasTuple,
        cells = cells,
        weightSource = weightSource
    ))
    mrTuples <- .cipMrTuples(fineMappingResult)
    .cipCheckBlockReach(c(twasTuples, mrTuples), gwasSumStats)
    gwasIndex <- .rowSpans(gwasSumStats)
    twasRows <- list_flatten(map(
        twasTuples,
        .cipScoreTwasTuple,
        gwasSumStats = gwasSumStats,
        gwasIndex = gwasIndex,
        alleleFlip = alleleFlip
    ))
    mrRows <- list_flatten(map(
        mrTuples,
        .cipScoreMrTuple,
        gwasSumStats = gwasSumStats,
        gwasIndex = gwasIndex,
        twasPvals = .cipTwasPvals(twasRows),
        alleleFlip = alleleFlip,
        mrParam = mrParam
    ))
    if (length(twasRows) + length(mrRows) == 0L) {
        abort(
            "causalInferencePipeline: no (qtl, gwas) tuples produced a result."
        )
    }
    .cipFinalize(c(twasRows, mrRows), combineMethods)
}

# Validate the object classes + QC state of the pipeline inputs.
# @noRd
#' @importFrom checkmate assertClass checkClass
.cipValidateInputs <- function(gwasSumStats, twasWeights, fineMappingResult) {
    assertClass(gwasSumStats, "GwasSumStats")
    if (length(qcOptions(gwasSumStats)) == 0L) {
        msg <- glue(
            "causalInferencePipeline: gwasSumStats has no QC record ",
            "(qcOptions() is empty). Call summaryStatsQc() first."
        )
        abort(msg)
    }
    if (is.null(twasWeights) && is.null(fineMappingResult)) {
        msg <- glue(
            "causalInferencePipeline: at least one of `twasWeights` or ",
            "`fineMappingResult` must be supplied."
        )
        abort(msg)
    }
    assertClass(twasWeights, "TwasWeights", null.ok = TRUE)
    # Composed rather than asserted: rejecting the GWAS-side class here is
    # deliberate, and checkClass alone cannot say so.
    res <- checkClass(
        fineMappingResult,
        "QtlFineMappingResult",
        null.ok = TRUE
    )
    if (!isTRUE(res)) {
        msg <- glue(
            "`fineMappingResult` {res} (causalInferencePipeline does not ",
            "accept GWAS-side fine mapping for the QTL slot)."
        )
        abort(msg)
    }
    invisible(NULL)
}

# Require the twasWeights / fineMappingResult rows to share an LD reference
# with each GWAS row they are scored against -- the rows whose regions
# overlap; rows that never meet may sit on different references.
# @noRd
.cipCheckLdRefs <- function(gwasSumStats, twasWeights, fineMappingResult) {
    if (!is.null(twasWeights)) {
        .cipRequireMatchingLdRefs(twasWeights, gwasSumStats, "twasWeights")
    }
    if (!is.null(fineMappingResult)) {
        .cipRequireMatchingLdRefs(
            fineMappingResult,
            gwasSumStats,
            "fineMappingResult"
        )
    }
    invisible(NULL)
}

# The LD reference a (QTL tuple, GWAS study) pair's TWAS z reads: the one the
# study's `rows` -- those the tuple overlaps -- share. The z folds their
# variants into a single quadratic form, so they cannot sit on different
# references; that is refused rather than read from one of them. The study's
# other rows are not read and may.
# @noRd
.cipRowsLdRef <- function(gwasSumStats, rows, gStudy) {
    .rowsSharedLdRef(
        gwasSumStats,
        rows,
        glue(
            "causalInferencePipeline: the rows of GWAS study '{gStudy}' ",
            "a QTL tuple overlaps"
        )
    )
}

# The (gene, method) weight sets of `weightSource` the run scores: every one it
# holds, narrowed by `weightSelectionParam` (.selectWeightCells()). An error
# when the source holds none, or the choice leaves none.
# @noRd
.cipChosenCells <- function(weightSource, weightSelectionParam) {
    cells <- .weightCells(weightSource, weightSelectionParam)
    if (nrow(cells) == 0L) {
        abort(glue(
            "causalInferencePipeline: no QTL tuples to score (the supplied ",
            "twasWeights / fineMappingResult collections are empty)."
        ))
    }
    .selectionCoverageCheck(
        list(weightSource),
        weightSelectionParam,
        "causalInferencePipeline()"
    )
    chosen <- .selectWeightCells(
        cells,
        weightSelectionParam,
        "causalInferencePipeline"
    )
    if (nrow(chosen) == 0L) {
        abort(glue(
            "causalInferencePipeline: every QTL tuple was filtered out by ",
            "WeightSelectionParam(cutoff = {weightSelectionParam$cutoff}, ",
            "pvalCutoff = {weightSelectionParam$pvalCutoff}) ",
            "(no method cleared the CV cutoffs)."
        ))
    }
    chosen
}

# Chosen cell `i`'s TWAS inputs: its identity, its weights, and the variants
# the GWAS is read at. NULL when it has no usable weights.
# @noRd
.cipTwasTuple <- function(i, cells, weightSource) {
    tuple <- .cipCellTuple(cells, i)
    weightsInfo <- .cipExtractWeights(weightSource, tuple)
    if (is.null(weightsInfo)) {
        return(NULL)
    }
    list(
        tuple = tuple,
        weightsInfo = weightsInfo,
        variants = .variantIdsToGRanges(
            unique(weightsInfo$variantIds),
            what = "causalInferencePipeline variant id"
        )
    )
}

# Row `i` of a weight-set table (.weightCellIds()) as a tuple identity: its
# (study, context, trait, method), its joint key, and for a TwasWeights cell
# the row it sits in.
# @noRd
.cipCellTuple <- function(cells, i) {
    list(
        qStudy = cells$study[[i]],
        qContext = cells$context[[i]],
        qTrait = cells$trait[[i]],
        qMethod = cells$method[[i]],
        qJoint = map(.cellsJointColumns(cells), .cipValueAt, i = i),
        qRow = if (is_in("row", names(cells))) cells$row[[i]]
    )
}

# @noRd
.cipValueAt <- function(values, i) {
    values[[i]]
}

# One MR input per fine-mapping row -- each (study, context, trait, method) it
# holds, a multi-block one as one -- with the variants it reads the GWAS at.
# None without a fine-mapping result.
# @noRd
.cipMrTuples <- function(fineMappingResult) {
    if (is.null(fineMappingResult)) {
        return(list())
    }
    keys <- .weightCellIds(fineMappingResult)
    map(
        seq_len(nrow(keys)),
        .cipMrTuple,
        keys = keys,
        fineMappingResult = fineMappingResult
    )
}

# @noRd
.cipMrTuple <- function(i, keys, fineMappingResult) {
    tuple <- .cipCellTuple(keys, i)
    fmrEntry <- .cipFmrRow(fineMappingResult, tuple)
    list(
        tuple = tuple,
        fmrEntry = fmrEntry,
        variants = .variantIdsToGRanges(
            unique(.grVariantIds(unlist(fmrEntry, use.names = FALSE))),
            what = "causalInferencePipeline variant id"
        )
    )
}

# Each tuple's variants against the LD blocks of the GWAS studies they
# overlap (.checkGwasBlockReach()).
# @noRd
.cipCheckBlockReach <- function(tuples, gwasSumStats) {
    .checkGwasBlockReach(
        set_names(map(tuples, "variants"), map_chr(tuples, .cipTupleLabel)),
        gwasSumStats,
        "causalInferencePipeline",
        "QTL tuples"
    )
}

# @noRd
.cipTupleLabel <- function(input) {
    t <- input$tuple
    .fmrTupleLabel(
        "QTL",
        c(
            list(
                study = t$qStudy,
                context = t$qContext,
                trait = t$qTrait,
                method = t$qMethod
            ),
            t$qJoint
        )
    )
}

# A weight set's joint key for a message, ", jointStudies = a;b"; empty for a
# set from no joint fit.
# @noRd
.cipJointLabel <- function(tuple) {
    set <- discard(tuple$qJoint, is.na)
    if (length(set) == 0L) "" else str_c(", ", .twasJointKeyLabel(set))
}

# The rows of each GWAS study that `window` overlaps, split by study. A study
# split by chromosome or LD block has many rows, and a tuple is scored once
# per study, over just the rows its window -- the span of its variants --
# overlaps. `gwasIndex` is .rowSpans(gwasSumStats), built once for the run.
# @noRd
.cipStudyRows <- function(gwasIndex, gwasSumStats, window) {
    rows <- unique(gwasIndex$rows[.overlapsOnShared(gwasIndex$spans, window)])
    studies <- as.character(gwasSumStats$study)[rows]
    split(rows, factor(studies, levels = unique(studies)))
}

# A GWAS study's statistics inside `window`, over its rows `rows`.
# @noRd
.cipGwasDf <- function(gwasSumStats, rows, window, gStudy) {
    .entryToSumstatDf(
        .ssStitchElements(gwasSumStats, rows, window),
        require = c("SNP", "Z"),
        derive = "none",
        keepChrPrefix = TRUE,
        label = glue("GwasSumStats[{gStudy}]")
    )
}

# Score one TWAS tuple against every GWAS study it overlaps -> a list of TWAS
# records (empty when it overlaps no GWAS row).
# @noRd
.cipScoreTwasTuple <- function(input, gwasSumStats, gwasIndex, alleleFlip) {
    window <- range(input$variants)
    compact(imap(
        .cipStudyRows(gwasIndex, gwasSumStats, window),
        .cipScoreTwasPair,
        input = input,
        gwasSumStats = gwasSumStats,
        window = window,
        alleleFlip = alleleFlip
    ))
}

# Score one (TWAS tuple, GWAS study) pair -> a TWAS record, or NULL when the
# TWAS Z cannot be computed (too little overlap). Only the study's statistics
# inside the window are read, so a gene whose weights cross an LD-block
# boundary is scored once, over every block it reaches.
# @noRd
.cipScoreTwasPair <- function(
    rows,
    gStudy,
    input,
    gwasSumStats,
    window,
    alleleFlip
) {
    gwasLd <- .cipRowsLdRef(gwasSumStats, rows, gStudy)
    twasOut <- .cipComputeTwasZ(
        weights = input$weightsInfo$weights,
        variantIds = input$weightsInfo$variantIds,
        gwasDf = .cipGwasDf(gwasSumStats, rows, window, gStudy),
        gwasLd = gwasLd,
        alleleFlip = alleleFlip,
        label = .cipPairLabel(input$tuple, gStudy)
    )
    if (is.null(twasOut)) {
        return(NULL)
    }
    .cipRecord(input$tuple, gStudy, "twas", twasOut, .cipEmptyMr(), twasOut)
}

# Run MR for one fine-mapping row against every GWAS study it overlaps -> one
# MR record per study, whatever its estimate came to. `twasPvals` holds the
# run's TWAS p-values by gene and study (.cipTwasPvals()), which mrParam's
# pvalCutoff reads.
# @noRd
.cipScoreMrTuple <- function(
    input,
    gwasSumStats,
    gwasIndex,
    twasPvals,
    alleleFlip,
    mrParam
) {
    window <- range(input$variants)
    imap(
        .cipStudyRows(gwasIndex, gwasSumStats, window),
        .cipScoreMrPair,
        input = input,
        gwasSumStats = gwasSumStats,
        window = window,
        twasPvals = twasPvals,
        alleleFlip = alleleFlip,
        mrParam = mrParam
    )
}

# @noRd
.cipScoreMrPair <- function(
    rows,
    gStudy,
    input,
    gwasSumStats,
    window,
    twasPvals,
    alleleFlip,
    mrParam
) {
    mrOut <- .cipRunMr(
        input$fmrEntry,
        .cipGwasDf(gwasSumStats, rows, window, gStudy),
        twasPvals[[.cipGeneKey(input$tuple, gStudy)]],
        alleleFlip = alleleFlip,
        mrParam = mrParam
    )
    .cipRecord(
        input$tuple,
        gStudy,
        "mr",
        .cipEmptyTwas(),
        mrOut,
        .cipWindowSpan(window)
    )
}

# A gene's key against GWAS study `gStudy`, as .cipTwasPvals() keys them.
# @noRd
.cipGeneKey <- function(tuple, gStudy) {
    str_c(tuple$qStudy, tuple$qContext, tuple$qTrait, gStudy, sep = "\r")
}

# The TWAS p-values of `records`, by gene and GWAS study (.cipGeneKey()).
# @noRd
.cipTwasPvals <- function(records) {
    if (length(records) == 0L) {
        return(list())
    }
    keys <- str_c(
        map_chr(records, "qtlStudy"),
        map_chr(records, "context"),
        map_chr(records, "trait"),
        map_chr(records, "gwasStudy"),
        sep = "\r"
    )
    split(map_dbl(records, "twasPval"), keys)
}

# Identify a (QTL tuple, GWAS study) pair in diagnostics.
# @noRd
.cipPairLabel <- function(tuple, gStudy) {
    glue(
        "QTL (study='{tuple$qStudy}', context='{tuple$qContext}', ",
        "trait='{tuple$qTrait}', method='{tuple$qMethod}') x GWAS ",
        "(study='{gStudy}')"
    )
}

# Run MR for a fine-mapping row against a GWAS study -> its estimate, or the
# all-NA one when mrParam's pvalCutoff withholds it: when none of the gene's
# TWAS p-values for the study (`genePvals`) is below the cutoff. A gene with
# no TWAS row for the study is not gated; a cutoff of 1 or more gates none.
# @noRd
.cipRunMr <- function(fmrEntry, gdf, genePvals, alleleFlip, mrParam) {
    # The bundle's terminal; the MR engines below take plain scalars.
    mrMethod <- mrParam$method
    mrPipCutoff <- mrParam$pipCutoff
    mrCpipCutoff <- mrParam$cpipCutoff
    mrPvalCutoff <- mrParam$pvalCutoff
    gateOpen <- mrPvalCutoff >= 1 ||
        length(genePvals) == 0L ||
        any(!is.na(genePvals) & genePvals < mrPvalCutoff)
    if (!gateOpen) {
        return(.cipEmptyMr())
    }
    if (mrMethod == "csAware") {
        .cipComputeMrCsAware(
            fmrEntry = fmrEntry,
            gwasDf = gdf,
            cpipCutoff = mrCpipCutoff,
            alleleFlip = alleleFlip
        )
    } else {
        .cipComputeMr(
            fmrEntry = fmrEntry,
            gwasDf = gdf,
            pipCutoff = mrPipCutoff,
            alleleFlip = alleleFlip
        )
    }
}

# The all-NA MR result: a TWAS row's, and an MR row's when its estimate is
# withheld.
# @noRd
.cipEmptyMr <- function() {
    list(
        waldRatio = NA_real_,
        waldRatioSe = NA_real_,
        mrPval = NA_real_,
        nIV = NA_integer_,
        Q = NA_real_,
        I2 = NA_real_,
        nCs = NA_integer_
    )
}

# The all-NA TWAS result an MR row carries.
# @noRd
.cipEmptyTwas <- function() {
    list(Z = NA_real_, pval = NA_real_)
}

# Where an MR row sits: the span of its fine-mapping row's variants.
# @noRd
.cipWindowSpan <- function(window) {
    list(
        chrom = as.character(seqnames(window))[[1L]],
        startPos = start(window)[[1L]],
        endPos = end(window)[[1L]]
    )
}

# One result record: a tuple's identity against GWAS study `gStudy`, which
# `analysis` it is, its TWAS and MR outputs (each all-NA on the other kind of
# row), and the span `at` (chrom, startPos, endPos) it is placed at.
# @noRd
.cipRecord <- function(tuple, gStudy, analysis, twasOut, mrOut, at) {
    list(
        qtlStudy = tuple$qStudy,
        context = tuple$qContext,
        trait = tuple$qTrait,
        method = tuple$qMethod,
        jointStudies = tuple$qJoint$jointStudies,
        jointContexts = tuple$qJoint$jointContexts,
        jointTraits = tuple$qJoint$jointTraits,
        gwasStudy = gStudy,
        analysis = analysis,
        twasZ = twasOut$Z,
        twasPval = twasOut$pval,
        waldRatio = mrOut$waldRatio,
        waldRatioSe = mrOut$waldRatioSe,
        mrPval = mrOut$mrPval,
        nIV = mrOut$nIV,
        Q = mrOut$Q %||% NA_real_,
        I2 = mrOut$I2 %||% NA_real_,
        nCs = mrOut$nCs %||% NA_integer_,
        chrom = at$chrom,
        startPos = at$startPos,
        endPos = at$endPos
    )
}

# The records as the result GRanges, with the TWAS rows' p-values combined
# across methods when asked.
# @noRd
.cipFinalize <- function(records, combineMethods) {
    perMethod <- .cipDfToGranges(.cipRowsToDf(records))
    if (is.null(combineMethods)) {
        return(perMethod)
    }
    .cipCombineAcrossMethods(perMethod, methods = combineMethods)
}

# =============================================================================
# Internal helpers
# =============================================================================

# Check the QTL-side rows against the GWAS rows they overlap for
# reference-panel compatibility. Thin wrapper over the shared
# `.requireOverlapLdRefs` helper (R/ld.R).
.cipRequireMatchingLdRefs <- function(qtl, gwasSumStats, label) {
    .requireOverlapLdRefs(
        qtl,
        gwasSumStats,
        pipelineName = "causalInferencePipeline",
        label = label
    )
}

# =============================================================================
# Method choice, shared with ctwasPipeline (WeightSelectionParam)
# =============================================================================

# One row per (gene, method, joint group) weight set of `source`: a
# TwasWeights' cells (with the `row` each sits in), or a fine-mapping result's
# distinct (study, context, trait, method, joint key) rows, a multi-block one
# once. A weight set's joint key (NA for none) tells the fits of one gene and
# method in different joint groups apart.
# @noRd
.weightCellIds <- function(source) {
    if (methods::is(source, "TwasWeights")) {
        cells <- .twasCells(source)
        return(tibble(
            study = cells$study,
            context = cells$context,
            trait = cells$trait,
            method = cells$method,
            !!!.cellsJointColumns(cells),
            row = cells$row
        ))
    }
    distinct(.weightRowIds(source))
}

# Each row of fine-mapping result `x` as its (study, context, trait, method)
# and joint key.
# @noRd
.weightRowIds <- function(x) {
    tibble(
        study = as.character(x$study),
        context = as.character(x$context),
        trait = as.character(x$trait),
        method = as.character(x$method),
        !!!set_names(
            map(.JOINT_COLUMNS, .fmrIdentityColumn, x = x),
            .JOINT_COLUMNS
        )
    )
}

# A cell table's three joint columns, NA-filled (.cellsJointColumn()).
# @noRd
.cellsJointColumns <- function(cells) {
    set_names(
        map(.JOINT_COLUMNS, .cellsJointColumn, cells = cells),
        .JOINT_COLUMNS
    )
}

# One key per weight set of a weight-set table (.weightCellIds()), to match
# weight sets with each other only.
# @noRd
.cellKey <- function(cells) {
    exec(
        str_c,
        cells$study,
        cells$context,
        cells$trait,
        cells$method,
        !!!map(.cellsJointColumns(cells), coalesce, ""),
        sep = "\r"
    )
}

# `source`'s weight sets (.weightCellIds()) with what WeightSelectionParam
# reads: each one's CV `rsq` (its `metric`) and `pval` (its `pvalMetric`), and
# whether its weights are all `zero`. They are read only when the choice needs
# them -- "best", or a cutoff (.selectionReadsCv()).
# @noRd
.weightCells <- function(source, weightSelectionParam) {
    cells <- .weightCellIds(source)
    if (nrow(cells) == 0L || !.selectionReadsCv(weightSelectionParam)) {
        return(mutate(cells, rsq = NA_real_, pval = NA_real_, zero = FALSE))
    }
    parts <- map(
        seq_len(nrow(cells)),
        .weightCellPart,
        cells = cells,
        source = source
    )
    mutate(
        cells,
        rsq = map_dbl(
            parts,
            .weightCvValue,
            which = weightSelectionParam$metric
        ),
        pval = map_dbl(
            parts,
            .weightCvValue,
            which = weightSelectionParam$pvalMetric
        ),
        zero = map_lgl(parts, .weightsAllZero)
    )
}

# Whether `weightSelectionParam`'s choice reads the CV metrics: "best" ranks
# by them, and a cutoff gates on them.
# @noRd
.selectionReadsCv <- function(weightSelectionParam) {
    identical(weightSelectionParam$methods, "best") ||
        weightSelectionParam$cutoff > 0 ||
        is.finite(weightSelectionParam$pvalCutoff)
}

# The CV metrics describe each fit on all of its variants, so a choice that
# reads them warns once about the rows of `sources` narrowed since their fits.
# @noRd
.selectionCoverageCheck <- function(sources, weightSelectionParam, caller) {
    if (.selectionReadsCv(weightSelectionParam)) {
        .fitCoverageCheckedAll(sources, caller)
    }
    invisible(NULL)
}

# Cell `i`'s payload: a TwasWeightsRow, or a FineMappingRow.
# @noRd
.weightCellPart <- function(i, cells, source) {
    if (methods::is(source, "TwasWeights")) {
        return(.twrCellParts(source, cells$row[[i]], cells$method[[i]]))
    }
    .fmrRowParts(source, .cipFmrRowIndex(source, .cipCellTuple(cells, i)))
}

# The first of the CV metrics `which` that payload `part` records, NA when it
# records none. A cvResult holds its metrics as a named vector or a one-row
# matrix; one without a `metrics` entry is read as the metrics themselves.
# @noRd
#' @importFrom rlang try_fetch
.weightCvValue <- function(part, which) {
    cv <- try_fetch(cvResult(part), error = function(cnd) NULL)
    metrics <- .cvMetricsVector(
        if (is.list(cv) && !is.null(cv[["metrics"]])) cv[["metrics"]] else cv
    )
    hit <- intersect(which, names(metrics))
    if (length(hit) == 0L) {
        return(NA_real_)
    }
    value <- suppressWarnings(as.numeric(metrics[[hit[[1L]]]]))
    if (length(value) == 0L) NA_real_ else value[[1L]]
}

# A metrics record as a named vector: a matrix's first row.
# @noRd
.cvMetricsVector <- function(metrics) {
    if (is.matrix(metrics) && nrow(metrics) > 0L) {
        return(metrics[1L, , drop = TRUE])
    }
    metrics
}

# Whether payload `part`'s weights are all zero, or absent.
# @noRd
.weightsAllZero <- function(part) {
    w <- if (methods::is(part, "TwasWeightsRow")) {
        .twrRowResolveWeights(part)$weights
    } else {
        .fmrRowResolveWeights(part)$weights
    }
    length(w) == 0L || all(w == 0, na.rm = TRUE)
}

# The cells of `cells` (.weightCells()) `weightSelectionParam` keeps. Its
# `jointAxes` narrows them to the kinds of joint fit it names, and its cutoffs
# bound every choice; then "best" keeps each gene's best weight set
# (.bestCellIn()), "all" every set, and named methods those. A name the cells
# do not hold is an error, as a misspelling would otherwise score nothing.
# @noRd
.selectWeightCells <- function(cells, weightSelectionParam, caller) {
    methods <- weightSelectionParam$methods
    kept <- .cellsOfKinds(cells, weightSelectionParam$jointAxes, caller)
    eligible <- .cellsByGene(
        kept,
        .cellsPassingIn,
        cutoff = weightSelectionParam$cutoff,
        pvalCutoff = weightSelectionParam$pvalCutoff
    )
    if (identical(methods, "all")) {
        return(eligible)
    }
    if (identical(methods, "best")) {
        return(.cellsByGene(eligible, .bestCellIn))
    }
    .assertCellMethods(methods, kept$method, caller)
    eligible[is_in(eligible$method, methods), ]
}

# The weight sets of the kinds of joint fit `jointAxes` names (.jointKinds()),
# every set when it is NULL. A choice keeping none is an error, as a kind the
# weights do not hold would otherwise score nothing.
# @noRd
.cellsOfKinds <- function(cells, jointAxes, caller) {
    if (is.null(jointAxes) || nrow(cells) == 0L) {
        return(cells)
    }
    kinds <- .jointKinds(.cellsJointColumns(cells))
    kept <- cells[is_in(kinds, jointAxes), ]
    if (nrow(kept) > 0L) {
        return(kept)
    }
    abort(glue(
        "{caller}: WeightSelectionParam(jointAxes = ) names ",
        "{str_flatten(jointAxes, ', ')}, which the weights do not hold ",
        "(they hold {str_flatten(unique(kinds), ', ')})."
    ))
}

# `cells` narrowed gene by gene: `keep(idx, cells, ...)` gives the rows of a
# gene's rows `idx` to keep. Row order is kept.
# @noRd
.cellsByGene <- function(cells, keep, ...) {
    if (nrow(cells) == 0L) {
        return(cells)
    }
    gene <- str_c(cells$study, cells$context, cells$trait, sep = "\r")
    byGene <- split(seq_len(nrow(cells)), factor(gene, levels = unique(gene)))
    kept <- list_c(map(byGene, keep, cells = cells, ...), ptype = integer())
    cells[sort(kept), ]
}

# The rows `idx` of one gene that clear the r-squared and p-value gates. A
# gate is set only by a cutoff above 0, or a finite p-value cutoff; and a gene
# none of whose methods has an r-squared keeps every row, there being nothing
# to gate it by (summary-statistic weights).
# @noRd
.cellsPassingIn <- function(idx, cells, cutoff, pvalCutoff) {
    rsq <- cells$rsq[idx]
    if (all(is.na(rsq))) {
        return(idx)
    }
    pval <- cells$pval[idx]
    byRsq <- cutoff <= 0 | (!is.na(rsq) & rsq >= cutoff)
    byPval <- !is.finite(pvalCutoff) | (!is.na(pval) & pval < pvalCutoff)
    idx[byRsq & byPval]
}

# Gene `idx`'s best row: the highest-r-squared method whose weights are not
# all zero, else the highest-r-squared one. Ties go to the first. A gene with
# no r-squared keeps every row, there being nothing to rank them by.
# @noRd
.bestCellIn <- function(idx, cells) {
    rsq <- cells$rsq[idx]
    if (all(is.na(rsq))) {
        return(idx)
    }
    ranked <- idx[order(rsq, decreasing = TRUE)]
    withWeights <- ranked[!cells$zero[ranked]]
    if (length(withWeights) > 0L) withWeights[[1L]] else ranked[[1L]]
}

# Named methods must be ones the weights hold.
# @noRd
.assertCellMethods <- function(methods, available, caller) {
    unknown <- setdiff(methods, available)
    if (length(unknown) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "{caller}: WeightSelectionParam(methods =) names ",
        "{str_flatten(unknown, ', ')}, which the weights do not hold ",
        "(they hold {str_flatten(unique(available), ', ')})."
    ))
}

# The first row of fine-mapping result `fmr` with weight set `tuple`'s full
# identity (.cipCellTuple()): its (study, context, trait, method) and joint
# key. The weight sets are read off `fmr`, so there is one.
# @noRd
.cipFmrRowIndex <- function(fmr, tuple) {
    .matchTupleRows(
        fmr,
        c(
            list(
                study = tuple$qStudy,
                context = tuple$qContext,
                trait = tuple$qTrait,
                method = tuple$qMethod
            ),
            tuple$qJoint
        )
    )[[1L]]
}

# That row, as a one-row fine-mapping result.
# @noRd
.cipFmrRow <- function(fmr, tuple) {
    fmr[.cipFmrRowIndex(fmr, tuple)]
}

# Fetch the per-tuple weight-source entry (the cell's TwasWeightsRow, or its
# row of a fine-mapping result) and resolve it to an aligned (variantIds,
# weights) pair as `resolveWeights` resolves one -- the SuSiE-style posterior
# effect for a fine-mapping result. Returns NULL when it has no usable
# weights.
.cipExtractWeights <- function(weightSource, tuple) {
    wr <- if (methods::is(weightSource, "TwasWeights")) {
        .twrRowResolveWeights(
            .twrCellParts(weightSource, tuple$qRow, tuple$qMethod)
        )
    } else {
        .fmrRowResolveWeights(
            .fmrRowParts(weightSource, .cipFmrRowIndex(weightSource, tuple))
        )
    }
    single <- .singleWeightColumn(wr$weights)
    which <- glue(
        "(study={tuple$qStudy}, context={tuple$qContext}, ",
        "trait={tuple$qTrait}, method={tuple$qMethod}{.cipJointLabel(tuple)})"
    )
    if (is.null(single)) {
        warn(glue(
            "causalInferencePipeline: the weights of {which} are a matrix ",
            "with {ncol(wr$weights)} columns, and a tuple is scored with one ",
            "weight vector; skipping."
        ))
        return(NULL)
    }
    if (length(wr$variantIds) == 0L) {
        warn(glue(
            "causalInferencePipeline: no weights resolved for {which}; ",
            "skipping. The row carries variants but no usable weight, or the ",
            "two are not the same length."
        ))
        return(NULL)
    }
    list_assign(wr, weights = single)
}

# Compute the per-tuple TWAS Z from a single GwasSumStats tuple's
# unpacked data.frame (produced by as.data.frame upstream). Returns
# NULL when the overlap is too small.
.cipComputeTwasZ <- function(
    weights,
    variantIds,
    gwasDf,
    gwasLd,
    alleleFlip = TRUE,
    label = "this (QTL, GWAS) pair"
) {
    # Match QTL weights to the GWAS sumstats by (chrom, pos, allele) rather than
    # by raw id string, so chr-prefix / separator / allele-swap differences are
    # reconciled instead of silently dropped. sign = -1 marks an allele swap;
    # weights[idxA] * sign brings the weight into the GWAS allele coding, which
    # the GWAS z and the GWAS-panel LD already share.
    m <- matchVariants(variantIds, gwasDf$variant_id, allowFlip = alleleFlip)
    # Reconciliation here is one-sided and unbiased -- subsetting both w and R
    # keeps Z = w'z / sqrt(w'Rw) honest, it only costs power -- but the caller
    # cannot see how much power was lost unless the drop is reported.
    .cipWarnVariantDrop(length(variantIds), length(m$idxA), label)
    if (length(m$idxA) < 2L) {
        return(NULL)
    }

    wSub <- weights[m$idxA] * m$sign
    zSub <- gwasDf$z[m$idxB]
    gwasIds <- gwasDf$variant_id[m$idxB]
    res <- .cipTwasZ(wSub, zSub, gwasLd, gwasIds)
    zMat <- res$Z
    zVal <- as.numeric(zMat[1L, "Z"])
    pVal <- as.numeric(zMat[1L, "pval"])
    # Position the row at the variant span.
    chrom <- gwasDf$chrom[[m$idxB[[1L]]]]
    startPos <- min(gwasDf$pos[m$idxB])
    endPos <- max(gwasDf$pos[m$idxB])
    list(
        Z = zVal,
        pval = pVal,
        chrom = chrom,
        startPos = startPos,
        endPos = endPos
    )
}

# twasZ() of weights `w` against GWAS z-scores `z` over the GWAS variants
# `gwasIds` (`w` and `z` in their order), needing w'Rw only: a genotype
# reference gives its LD factor and an LdMatrix its packed triangle, so the
# p x p LD is never formed.
# @noRd
.cipTwasZ <- function(w, z, gwasLd, gwasIds) {
    label <- "causalInferencePipeline"
    if (.ldIsGenotypeRef(gwasLd)) {
        return(twasZ(
            weights = w,
            z = z,
            X = .ldFactorFromRef(gwasLd, gwasIds, label)
        ))
    }
    weights <- .twasZPrepWeights(unname(w), z)
    ld <- .ldPackedFromRef(gwasLd, gwasIds, label)
    .twasZStatistics(
        weights,
        z,
        crossprod(weights, .ldSourceTimes(ld, weights))
    )
}

# Report how many QTL weight variants failed to match the GWAS sumstats.
# Dropping below two matches yields no TWAS Z at all, which used to be a silent
# NULL, so that case is called out separately.
# @noRd
.cipWarnVariantDrop <- function(nWeights, nMatched, label) {
    if (nMatched >= nWeights) {
        return(invisible(NULL))
    }
    dropped <- nWeights - nMatched
    if (nMatched < 2L) {
        msg <- glue(
            "causalInferencePipeline: {label} matched only {nMatched} of ",
            "{nWeights} weight variant(s) to the GWAS sumstats; ",
            "at least 2 are needed for a TWAS Z, so no result is reported."
        )
    } else {
        msg <- glue(
            "causalInferencePipeline: {label} dropped {dropped} of ",
            "{nWeights} weight variant(s) absent from the GWAS sumstats; ",
            "the TWAS Z uses the remaining {nMatched}."
        )
    }
    warn(msg)
}

# Compute the Wald-ratio IVW MR estimate for a single tuple. Uses the
# FineMappingRow's topLoci as the instrumental variable source: each
# variant with PIP > pipCutoff contributes one ratio = beta_y / beta_x.
# Returns list(waldRatio, waldRatioSe, mrPval, nIV) with NA fields when
# no IVs survive.
.cipComputeMr <- function(fmrEntry, gwasDf, pipCutoff, alleleFlip = TRUE) {
    iv <- .cipMrInstruments(fmrEntry, pipCutoff)
    if (is.null(iv)) {
        return(.cipEmptyMrRatio())
    }
    m <- matchVariants(iv$ivVars, gwasDf$variant_id, allowFlip = alleleFlip)
    if (length(m$idxA) == 0L) {
        return(.cipEmptyMrRatio())
    }
    gIdx <- m$idxB
    gZ <- gwasDf$z[gIdx]
    gN <- .cipGwasCol(gwasDf$N, gIdx)
    gMaf <- .cipGwasCol(gwasDf$maf, gIdx)
    # Align the QTL exposure effect to the GWAS allele coding (sign = -1 on an
    # allele swap) so the Wald ratio betaY / betaX has the correct sign, and
    # put it on the same scale the outcome lands on.
    exposure <- .cipMrExposureScale(
        iv$betaX[m$idxA] * m$sign,
        iv$seX[m$idxA],
        .cipGwasHasScale(gMaf, gN)
    )
    betaX <- exposure$beta
    seX <- exposure$se
    betaY <- .cipZToBeta(gZ, gMaf, gN)
    seY <- .cipZToSe(gZ, gMaf, gN)
    ratio <- betaY / betaX
    # Standard Wald-ratio SE via delta method.
    rSe <- sqrt((seY / betaX)^2 + (betaY * seX / betaX^2)^2)
    # IVW pooling of the per-instrument Wald ratios.
    pooled <- .ivwPool(ratio, rSe)
    if (pooled$n == 0L) {
        return(.cipEmptyMrRatio())
    }
    list(
        waldRatio = pooled$effect,
        waldRatioSe = pooled$se,
        mrPval = pooled$pval,
        nIV = pooled$n
    )
}

# The all-NA per-variant Wald-ratio MR result (4 fields).
# @noRd
.cipEmptyMrRatio <- function() {
    list(
        waldRatio = NA_real_,
        waldRatioSe = NA_real_,
        mrPval = NA_real_,
        nIV = 0L
    )
}

# The per-variant top-loci table of one entry. topLoci() returns ranges and
# mcols() is exactly the pre-conversion table -- as.data.frame() on the
# ranges would add the standard seqnames/start/end/width/strand columns the
# MR column matching does not expect. Both instrument builders read through
# here, so the equivalence is stated once and the tests have one seam to
# substitute a controlled table at.
# @noRd
.cipTopLociTable <- function(fmrEntry) {
    as.data.frame(S4Vectors::mcols(topLoci(fmrEntry)))
}

# Instrumental variables from a FineMappingRow's topLoci: variants with PIP >
# pipCutoff. Returns list(ivVars, betaX, seX) or NULL when none qualify.
# @noRd
.cipMrInstruments <- function(fmrEntry, pipCutoff) {
    tl <- .cipTopLociTable(fmrEntry)
    if (is.null(tl) || nrow(tl) == 0L) {
        return(NULL)
    }
    cols <- .cipTlCols(tl)
    .cipRequireMrCols(cols, c("pip", "beta", "se"), "the Wald-ratio path")
    keep <- !is.na(tl[[cols$pip[[1L]]]]) & tl[[cols$pip[[1L]]]] > pipCutoff
    ivVars <- as.character(tl$variant_id)[keep]
    if (length(ivVars) == 0L) {
        return(NULL)
    }
    list(
        ivVars = ivVars,
        betaX = as.numeric(tl[[cols$beta[[1L]]]])[keep],
        seX = as.numeric(tl[[cols$se[[1L]]]])[keep]
    )
}

# Index a possibly-absent GWAS column, filling NA when the column is missing.
# @noRd
.cipGwasCol <- function(col, gIdx) {
    if (!is.null(col)) col[gIdx] else rep(NA_real_, length(gIdx))
}

# Cochran's Q-based I-squared heterogeneity statistic. Ported from
# the legacy mr.R::calcI2. Q = 0 or near-zero -> I2 = 0; clipped to
# [0, 1] to stay in the usual heterogeneity convention.
# @noRd
.cipCalcI2 <- function(Q, nGroups) {
    if (!is.finite(Q) || Q <= 1e-3 || nGroups <= 1L) {
        return(0)
    }
    i2 <- (Q - (nGroups - 1)) / Q
    max(0, min(1, i2))
}

# CS-aware Wald-ratio MR estimate for a single tuple. Adapted from the
# legacy mr.R::mrAnalysis but operating on a `FineMappingRow`'s topLoci
# data.frame + a GWAS `GRanges` instead of the old data.frame-only API.
#
# Logic:
#   1. Group topLoci variants by credible set (column `cs` in topLoci, or
#      first column matching ^cs -- cs_95 on a stored result, whose
#      "<method>_<idx>" labels are parsed to the set index).
#   2. Per CS, compute cumulative PIP; drop CSs with cpip < cpipCutoff.
#   3. Per surviving CS, compute a PIP-weighted composite Wald ratio
#      `composite_bhat = sum((bhatY/bhatX * pip) / cpip)` with the
#      delta-method composite SE.
#   4. IVW-pool composite Walds across CSs with weights wv = 1/se^2.
#   5. Report (metaEff, metaSe, metaPval, Q, I2, nCs, nIv) where the
#      result list keys are renamed (waldRatio/waldRatioSe/mrPval/nIV) so
#      the calling pipeline emits consistent column names.
#
# Returns list(waldRatio, waldRatioSe, mrPval, nIV, Q, I2, nCs) with NA
# fields when no usable CS survives.
# @noRd
.cipComputeMrCsAware <- function(
    fmrEntry,
    gwasDf,
    cpipCutoff,
    alleleFlip = TRUE
) {
    inst <- .cipCsAwareInstruments(fmrEntry, gwasDf, alleleFlip)
    if (is.null(inst)) {
        return(.cipEmptyMrCs())
    }
    comp <- .cipCsComposites(inst, cpipCutoff)
    if (length(comp$bhat) == 0L) {
        return(.cipEmptyMrCs())
    }
    .cipCsAwareResult(comp, inst)
}

# The all-NA CS-aware MR result (7 fields).
# @noRd
.cipEmptyMrCs <- function() {
    list(
        waldRatio = NA_real_,
        waldRatioSe = NA_real_,
        mrPval = NA_real_,
        nIV = 0L,
        Q = NA_real_,
        I2 = NA_real_,
        nCs = 0L
    )
}

# Extract + GWAS-align the CS-aware exposure instruments. Returns the aligned,
# unit-SE-standardized (cs, pip, bhatX, sbhatX, bhatY, sbhatY) or NULL.
# @noRd
.cipCsAwareInstruments <- function(fmrEntry, gwasDf, alleleFlip) {
    tl <- .cipTopLociTable(fmrEntry)
    if (is.null(tl) || nrow(tl) == 0L) {
        return(NULL)
    }
    cols <- .cipTlCols(tl)
    .cipRequireMrCols(
        cols,
        c("pip", "beta", "se", "cs"),
        "the credible-set-aware path"
    )
    raw <- .cipCsAwareRaw(tl, cols)
    if (is.null(raw)) {
        return(NULL)
    }
    .cipCsAwareAlign(raw, gwasDf, alleleFlip)
}

# Pull the CS / PIP / beta / se / variant columns from topLoci, dropping rows
# with NA cs / non-positive cs / NA pip / NA beta. Returns NULL when none valid.
# @noRd
.cipCsAwareRaw <- function(tl, cols) {
    csCol <- tl[[cols$cs[[1L]]]]
    # A stored topLoci table labels each variant's set "<method>_<idx>"
    # (cs_95); a caller-built table may carry the bare integer instead.
    cs <- if (is.character(csCol)) .csLabelIndex(csCol) else as.integer(csCol)
    pip <- as.numeric(tl[[cols$pip[[1L]]]])
    bhatX <- as.numeric(tl[[cols$beta[[1L]]]])
    sbhatX <- as.numeric(tl[[cols$se[[1L]]]])
    vids <- as.character(tl$variant_id)
    ok <- !is.na(cs) & cs > 0L & !is.na(pip) & !is.na(bhatX) & !is.na(sbhatX)
    if (!any(ok)) {
        return(NULL)
    }
    list(
        cs = cs[ok],
        pip = pip[ok],
        bhatX = bhatX[ok],
        sbhatX = sbhatX[ok],
        vids = vids[ok]
    )
}

# Match the exposure variants to the GWAS sumstats (allele-aware), derive the
# outcome beta/se, and put the exposure on the same scale. Returns NULL when
# nothing matches.
#
# The exposure used to be standardized to unit SE unconditionally, which left
# csAware dividing a real-scale outcome by a z-scale exposure whenever the
# GWAS did carry maf and N -- the mirror of the ivwPerVariant mismatch. It is
# conditional now, on the same test the outcome conversion uses.
# @noRd
.cipCsAwareAlign <- function(raw, gwasDf, alleleFlip) {
    m <- matchVariants(raw$vids, gwasDf$variant_id, allowFlip = alleleFlip)
    if (length(m$idxA) == 0L) {
        return(NULL)
    }
    gIdx <- m$idxB
    gZ <- gwasDf$z[gIdx]
    gN <- .cipGwasCol(gwasDf$N, gIdx)
    gMaf <- .cipGwasCol(gwasDf$maf, gIdx)
    exposure <- .cipMrExposureScale(
        raw$bhatX[m$idxA] * m$sign,
        raw$sbhatX[m$idxA],
        .cipGwasHasScale(gMaf, gN)
    )
    list(
        cs = raw$cs[m$idxA],
        pip = raw$pip[m$idxA],
        bhatX = exposure$beta,
        sbhatX = exposure$se,
        bhatY = .cipZToBeta(gZ, gMaf, gN),
        sbhatY = .cipZToSe(gZ, gMaf, gN)
    )
}

# Per-CS PIP-weighted composite Wald ratios (dropping CSs below cpipCutoff or
# with a degenerate SE). Returns list(bhat, sbhat) vectors.
# @noRd
.cipCsComposites <- function(inst, cpipCutoff) {
    csIds <- sort(unique(inst$cs))
    comps <- compact(map(
        csIds,
        .cipCsComposite,
        inst = inst,
        cpipCutoff = cpipCutoff
    ))
    list(bhat = map_dbl(comps, "bhat"), sbhat = map_dbl(comps, "sbhat"))
}

# One CS's composite Wald ratio + delta-method SE, or NULL when the CS is
# dropped (cpip below cutoff, or non-finite / non-positive SE).
# @noRd
.cipCsComposite <- function(cid, inst, cpipCutoff) {
    inCs <- which(inst$cs == cid)
    cpip <- sum(inst$pip[inCs])
    if (!is.finite(cpip) || cpip < cpipCutoff) {
        return(NULL)
    }
    pNorm <- inst$pip[inCs] / cpip
    bx <- inst$bhatX[inCs]
    sx <- inst$sbhatX[inCs]
    by <- inst$bhatY[inCs]
    sy <- inst$sbhatY[inCs]
    cBhat <- sum((by / bx) * pNorm)
    second <- sum(
        ((by / bx)^2 + (sy^2 / bx^2) + ((by^2 * sx^2) / bx^4)) * pNorm
    )
    cSe <- sqrt(max(0, second - cBhat^2))
    if (!is.finite(cBhat) || !is.finite(cSe) || cSe <= 0) {
        return(NULL)
    }
    list(bhat = cBhat, sbhat = cSe)
}

# IVW-meta the per-CS composite Walds into the final CS-aware MR result.
# @noRd
.cipCsAwareResult <- function(comp, inst) {
    pooled <- .ivwPool(comp$bhat, comp$sbhat)
    list(
        waldRatio = pooled$effect,
        waldRatioSe = pooled$se,
        mrPval = pooled$pval,
        nIV = length(inst$cs),
        Q = pooled$Q,
        I2 = .cipCalcI2(pooled$Q, pooled$n),
        nCs = length(comp$bhat)
    )
}


# Whether the GWAS can be put on an effect-size scale at all, which needs
# both maf and N. All-or-nothing per call, matching how the conversions below
# behave.
#
# Stated once because three callers have to agree on it: the outcome beta, the
# outcome se, and the EXPOSURE scaling. While the exposure had no such test,
# a z-scale outcome was divided by a raw-scale QTL effect and the Wald ratio
# came out inflated by 1/seX -- answers like -1094 +/- 1105 for a tight
# instrument.
# @noRd
.cipGwasHasScale <- function(maf, n) {
    !any(is.na(maf)) && !any(is.na(n))
}

# Derive beta / se from z using maf + n via the shared zToBetaSe() (model-exact
# se = 1/sqrt(2*p*q*(N + z^2)), beta = z*se). Fall back to z as a beta surrogate
# / se = 1 when maf or n is unavailable.
.cipZToBeta <- function(z, maf, n) {
    if (!.cipGwasHasScale(maf, n)) {
        return(z)
    }
    .zToBetaSe(z, maf, n)$beta
}
.cipZToSe <- function(z, maf, n) {
    if (!.cipGwasHasScale(maf, n)) {
        return(rep(1, length(z)))
    }
    .zToBetaSe(z, maf, n)$se
}

# The exposure on whatever scale the outcome ended up on. With maf and N both
# sides are real effect sizes, so the Wald ratio is an effect-size ratio and
# the raw QTL estimate is the one worth having. Without them the outcome is a
# z with unit se, so the exposure standardizes to z as well.
#
# Both MR methods route through this, so neither can end up on a scale the
# other is not on -- which is how ivwPerVariant paired a z outcome with a raw
# exposure, and csAware a real-scale outcome with a standardized one.
# @noRd
.cipMrExposureScale <- function(beta, se, gwasHasScale) {
    if (gwasHasScale) {
        return(list(beta = beta, se = se))
    }
    list(beta = beta / se, se = rep(1, length(beta)))
}

# Fixed-effect inverse-variance-weighted (IVW) pooling of per-instrument
# effects. Weights 1/se^2, drops non-finite / non-positive weights, and
# returns the pooled effect, its SE (1/sqrt(sum w)), two-tailed p-value,
# Cochran's Q, and the number pooled (n = 0 when nothing is poolable).
.ivwPool <- function(effect, se) {
    invVar <- 1 / se^2
    ok <- is.finite(invVar) & invVar > 0
    if (!any(ok)) {
        return(list(
            effect = NA_real_,
            se = NA_real_,
            pval = NA_real_,
            Q = NA_real_,
            n = 0L
        ))
    }
    effect <- effect[ok]
    w <- invVar[ok]
    metaEff <- sum(w * effect) / sum(w)
    metaSe <- 1 / sqrt(sum(w))
    list(
        effect = metaEff,
        se = metaSe,
        pval = .zToPvalue(metaEff / metaSe),
        Q = sum(w * (effect - metaEff)^2),
        n = length(effect)
    )
}

# Resolve the topLoci column aliases used across the CIP MR helpers. Returns the
# first matching column name for each role (character(0) if absent); `cs` also
# falls back to any column whose name starts with "cs" (e.g. cs_0.95).
# Abort naming the columns Mendelian randomization needs but the exposure does
# not have.
#
# These used to return NULL, so a result lacking `se` -- routine for a table
# converted from another tool, or read back from a published file whose source
# recorded only z -- produced "no (qtl, gwas) tuples produced a result" with
# nothing pointing at the cause. A column that is ABSENT is a structural
# property of the input and identical for every tuple, so it is worth stopping
# over; a column that is present but empty for one tuple stays a skip.
# @noRd
.cipRequireMrCols <- function(cols, need, label) {
    absent <- need[map_int(cols[need], length) == 0L]
    if (length(absent) == 0L) {
        return(invisible(NULL))
    }
    abort(glue(
        "causalInferencePipeline: {label} needs ",
        "{str_flatten(need, ', ')} on the exposure, and this result has no ",
        "{str_flatten(absent, ', ')}. A result read from a published table ",
        "carries beta/se only if the source recorded them, or if af and n ",
        "let .zToBetaSe() derive them."
    ))
}

.cipTlCols <- function(tl) {
    cn <- colnames(tl)
    exact <- intersect("cs", cn)
    # No bare `cs` column: fall back to the first cs_<coverage> variant.
    csCand <- cn[str_detect(cn, "^cs")]
    cs <- if (length(exact) > 0L) {
        exact
    } else if (length(csCand) > 0L) {
        csCand[[1L]]
    } else {
        exact
    }
    list(
        pip = intersect(c("pip", "PIP"), cn),
        beta = intersect(
            c("betahat", "beta", "bhat_x", "posterior_effect"),
            cn
        ),
        se = intersect(
            c("sebetahat", "se", "sbhat_x", "posterior_se"),
            cn
        ),
        cs = cs
    )
}

# Convert the accumulated list of row records to a flat data.frame.
.cipRowsToDf <- function(rows) {
    bind_rows(rows)
}

# Convert the assembled result data.frame to a GRanges with mcols.
.cipDfToGranges <- function(df) {
    gr <- GenomicRanges::GRanges(
        seqnames = df$chrom,
        ranges = IRanges::IRanges(
            start = as.integer(df$startPos),
            end = as.integer(df$endPos)
        )
    )
    mcols <- select(
        df,
        all_of(c(
            "qtlStudy",
            "context",
            "trait",
            "method",
            keep(.JOINT_COLUMNS, .cipColumnHasValues, df = df),
            "gwasStudy",
            "analysis",
            "twasZ",
            "twasPval",
            "waldRatio",
            "waldRatioSe",
            "mrPval",
            "nIV",
            "Q",
            "I2",
            "nCs"
        ))
    )
    S4Vectors::`mcols<-`(gr, value = S4Vectors::DataFrame(mcols))
}

# Whether a joint column of the result records holds any value: a joint
# column is reported only when some scored weight set came from a joint fit.
# @noRd
.cipColumnHasValues <- function(col, df) {
    is_in(col, names(df)) && any(!is.na(df[[col]]))
}

# Combine TWAS p-values across method for each (qtlStudy, context,
# trait, gwasStudy) group of TWAS rows. Appends one row per group with a
# combined p-value and methodName = "combined.<methodToken>". Uses
# combinePValues() with the cross-method correlation set to the identity
# (we have no cross-method covariance available downstream).
# @noRd
.cipConcat <- function(pieces) {
    if (length(pieces) == 0L) {
        return(list())
    }
    list_c(pieces)
}

# One combined row: the group's first TWAS row, carrying the combined p-value
# for method `m` in place of its own, and placed where that row is (`.at`).
# Its own Z does not carry over, so it is blanked rather than inherited; a
# TWAS row has no MR to carry, and pooling rows of several fits, it has no
# joint group of its own.
# @noRd
.cipCombinedRow <- function(m, mc, rows, cp) {
    mutate(
        slice(mc, rows[[1L]]),
        method = str_c("combined.", m),
        twasZ = NA_real_,
        twasPval = as.numeric(cp$results[[m]]$pval),
        across(any_of(.JOINT_COLUMNS), .cipNoJointGroup),
        .at = rows[[1L]]
    )
}

# @noRd
.cipNoJointGroup <- function(values) {
    rep(NA_character_, length(values))
}

# One group's combined rows, or none when it holds a single method (nothing
# to combine across).
# @noRd
.cipCombinedRowsFor <- function(rows, mc, methods) {
    if (length(rows) < 2L) {
        return(list())
    }
    cp <- combinePValues(
        pvals = as.numeric(mc$twasPval[rows]),
        zScores = as.numeric(mc$twasZ[rows]),
        methods = methods
    )
    map(methods, .cipCombinedRow, mc = mc, rows = rows, cp = cp)
}

.cipCombineAcrossMethods <- function(gr, methods) {
    mc <- as_tibble(as.data.frame(S4Vectors::mcols(gr)))
    twas <- which(mc$analysis == "twas")
    key <- str_c(mc$qtlStudy, mc$context, mc$trait, mc$gwasStudy, sep = "||")
    groups <- split(twas, key[twas])
    extras <- .cipConcat(map(
        groups,
        .cipCombinedRowsFor,
        mc = mc,
        methods = methods
    ))
    if (length(extras) == 0L) {
        return(gr)
    }
    newMcs <- bind_rows(extras)
    newGr <- S4Vectors::`mcols<-`(
        gr[newMcs$.at],
        value = S4Vectors::DataFrame(select(newMcs, -".at"))
    )
    c(gr, newGr)
}


# =============================================================================
# Unified TWAS Z-statistic
# =============================================================================

# Internal: build the K x K covariance W^T R W. Uses the SVD path when the
# triplet (V, D, nSketch) is supplied, otherwise the R / X path. Aligns LD
# rows/cols to the rownames of W when both are named; falls back to
# positional alignment otherwise.
.twasZCovY <- function(
    weights,
    R = NULL,
    X = NULL,
    V = NULL,
    D = NULL,
    nSketch = NULL
) {
    rn <- rownames(weights)
    if (!is.null(V) && !is.null(D) && !is.null(nSketch)) {
        return(.twasZCovYSvd(weights, V, D, nSketch, rn))
    }
    .twasZCovYR(weights, R, X, rn)
}

# SVD path: covY = (V^TW * sqrt(Lambda))^T (V^TW * sqrt(Lambda)) with
# Lambda_i = D_i^2 / (nSketch - 1).
# @noRd
.twasZCovYSvd <- function(weights, V, D, nSketch, rn) {
    vSub <- .twasZAlignV(V, rn, nrow(weights))
    Lambda <- D^2 / (nSketch - 1)
    VtW <- crossprod(vSub, weights)
    list(covY = crossprod(VtW * sqrt(Lambda)))
}

# Row-align V to the weights (by name when available, else positionally).
# @noRd
#' @importFrom checkmate assertMatrix
.twasZAlignV <- function(V, rn, nW) {
    if (!is.null(rownames(V)) && !is.null(rn)) {
        idx <- match(rn, rownames(V))
        if (anyNA(idx)) {
            msg <- glue(
                "twasZ: V is missing rows for {sum(is.na(idx))} ",
                "variant(s) named in weights."
            )
            abort(msg)
        }
        return(V[idx, , drop = FALSE])
    }
    assertMatrix(V, nrows = nW, .var.name = "V (positional alignment)")
    V
}

# R / X path: covY = W^T R W. From a genotype matrix X the LD is never formed:
# W^T R W = crossprod(F W) / (n - 1) for F, X's LD factor
# (.ldFactorFromDosage()), which is computeLd()'s LD exactly.
# @noRd
.twasZCovYR <- function(weights, R, X, rn) {
    if (!is.null(R)) {
        rSub <- .twasZAlignR(R, rn, nrow(weights))
        return(list(covY = crossprod(weights, rSub) %*% weights))
    }
    if (is.null(X)) {
        abort("twasZ: provide R, X, or the (V, D, nSketch) SVD triplet.")
    }
    factor <- .ldFactorFromDosage(.twasZAlignX(X, rn, nrow(weights)))
    list(covY = crossprod(factor %*% weights) / (nrow(factor) - 1))
}

# Column-align X to the weights (by name when available, else positionally).
# @noRd
#' @importFrom checkmate assertMatrix
.twasZAlignX <- function(X, rn, nW) {
    if (!is.null(colnames(X)) && !is.null(rn)) {
        idx <- match(rn, colnames(X))
        if (anyNA(idx)) {
            msg <- glue(
                "twasZ: X is missing columns for {sum(is.na(idx))} ",
                "variant(s) named in weights."
            )
            abort(msg)
        }
        return(X[, idx, drop = FALSE])
    }
    assertMatrix(X, ncols = nW, .var.name = "X (positional alignment)")
    X
}

# Symmetric-align R to the weights (by name when available, else positionally).
# @noRd
#' @importFrom checkmate assertMatrix
.twasZAlignR <- function(R, rn, nW) {
    if (!is.null(rownames(R)) && !is.null(rn)) {
        idx <- match(rn, rownames(R))
        if (anyNA(idx)) {
            msg <- glue(
                "twasZ: R is missing rows for {sum(is.na(idx))} ",
                "variant(s) named in weights."
            )
            abort(msg)
        }
        return(R[idx, idx, drop = FALSE])
    }
    assertMatrix(R, nrows = nW, .var.name = "R (positional alignment)")
    R
}

#' Calculate TWAS Z-Statistics for One or More Methods / Contexts
#'
#' Unified TWAS Z-statistic: accepts a weight vector (single method/context) or
#' a (variants x K) weight matrix, computes the per-tuple TWAS Z-score and
#' two-sided p-value, and optionally delegates cross-tuple p-value combination
#' to \code{\link{combinePValues}}.
#'
#' For each column k of \code{weights}:
#' \itemize{
#'   \item \code{stat_k = w_k^T z}
#'   \item \code{denom_k = w_k^T R w_k}
#'   \item \code{Z_k = stat_k / sqrt(denom_k)}, \code{p_k = 2 * (1 -
#'   Phi(|Z_k|))}
#' }
#' When \code{combineMethods} is non-NULL and K >= 2, the cross-tuple
#' correlation matrix \code{rho_{i,j} = covY_{i,j} / sqrt(covY_{i,i} *
#' covY_{j,j})} is constructed once and forwarded to
#' \code{combinePValues} as the \code{R} argument. When K == 1, the
#' combined p-value trivially equals the per-tuple p-value.
#'
#' The SVD path (\code{V}, \code{D}, \code{nSketch}) lets the caller avoid
#' materializing the full LD matrix: \code{covY = (V^TW * sqrt(Lambda))^T (V^TW
#' * sqrt(Lambda))} with \code{Lambda_i = D_i^2 / (nSketch - 1)}. Use the
#' \code{R} path when an LD correlation matrix is already available; the
#' \code{X} path works from a genotype matrix directly, as
#' \code{covY = (F W)^T (F W) / (n - 1)} for its standardized dosages
#' \code{F}, never forming \code{R}.
#'
#' @param weights Numeric vector of weights (single tuple) or a numeric matrix
#'   with one column per tuple (method / context). When a vector, the column
#'   name defaults to \code{"method1"}.
#' @param z Numeric vector of GWAS Z-scores aligned to the rows of
#'   \code{weights}.
#' @param R Optional LD correlation matrix.
#' @param X Optional genotype matrix (samples x variants) standing in for
#'   \code{R} when \code{R} is missing; its LD is taken as
#'   \code{\link{computeLd}}'s, without forming it.
#' @param V,D,nSketch SVD components of the LD sketch (right-singular vectors,
#'   singular values, panel sample size). Supplying all three selects the SVD
#'   path.
#' @param combineMethods Optional character vector of method names to forward to
#'   \code{\link{combinePValues}} for cross-tuple combination. \code{NULL}
#'   (default) skips combination.
#' @return A list with:
#' \describe{
#'   \item{Z}{A \code{K x 2} numeric matrix with columns
#'     \code{c("Z", "pval")}; rownames are the column names of
#'     \code{weights}.}
#'   \item{combined}{Output of \code{combinePValues} (or its trivial
#'     K=1 equivalent) when \code{combineMethods} is non-NULL,
#'     otherwise \code{NULL}.}
#' }
#' @seealso \code{\link{combinePValues}} for the combination method menu.
#' @importFrom stats pnorm
#' @examples
#' data(eqtlRegionExample)
#' X <- eqtlRegionExample$X[, 1:20]
#' w <- setNames(rnorm(20) * 0.1, colnames(X))
#' twasZ(weights = w, z = rnorm(20), R = cor(X))
#' @export
#' @importFrom checkmate assertCharacter assertCount
twasZ <- function(
    weights,
    z,
    R = NULL,
    X = NULL,
    V = NULL,
    D = NULL,
    nSketch = NULL,
    combineMethods = NULL
) {
    assertCount(nSketch, positive = TRUE, null.ok = TRUE)
    assertCharacter(combineMethods, null.ok = TRUE)
    weights <- .twasZPrepWeights(weights, z)
    covY <- .twasZCovY(
        weights = weights,
        R = R,
        X = X,
        V = V,
        D = D,
        nSketch = nSketch
    )$covY
    .twasZStatistics(weights, z, covY, combineMethods)
}

# twasZ()'s statistics from the weights' covariance under the LD,
# covY = W' R W, however that was computed.
# @noRd
.twasZStatistics <- function(weights, z, covY, combineMethods = NULL) {
    ySd <- sqrt(diag(covY))
    stats <- as.numeric(crossprod(weights, as.numeric(z)))
    zVec <- stats / ySd
    pVec <- .zToPvalue(zVec)
    zMatrix <- `rownames<-`(
        cbind(Z = zVec, pval = pVec),
        colnames(weights)
    )
    combined <- .twasZCombine(
        combineMethods,
        ncol(weights),
        pVec,
        zVec,
        covY,
        ySd,
        zMatrix
    )
    list(Z = zMatrix, combined = combined)
}

# Coerce a weight vector to a one-column matrix, validate the class / dims, and
# default the column names.
# @noRd
#' @importFrom checkmate assertMatrix
.twasZPrepWeights <- function(weights, z) {
    if (is.numeric(weights) && is.null(dim(weights))) {
        nm <- if (!is.null(names(weights))) names(weights) else NULL
        weights <- matrix(weights, ncol = 1L, dimnames = list(nm, "method1"))
    }
    if (!is.matrix(weights)) {
        abort("`weights` must be a numeric vector or a matrix.")
    }
    named <- if (is.null(colnames(weights))) {
        `colnames<-`(weights, str_c("method", seq_len(ncol(weights))))
    } else {
        weights
    }
    assertMatrix(named, nrows = length(z), .var.name = "weights")
    named
}

# Optional cross-tuple p-value combination via combinePValues (K == 1 uses the
# trivial identity form). Returns NULL when combineMethods is NULL.
# @noRd
.twasZCombine <- function(combineMethods, K, pVec, zVec, covY, ySd, zMatrix) {
    if (is.null(combineMethods)) {
        return(NULL)
    }
    combineMethods <- as.character(combineMethods)
    if (K == 1L) {
        return(.twasZCombineSingle(combineMethods, pVec, zMatrix))
    }
    methodNames <- rownames(zMatrix)
    sig <- `dimnames<-`(
        covY / tcrossprod(ySd, ySd),
        list(methodNames, methodNames)
    )
    combinePValues(
        pvals = set_names(pVec, methodNames),
        zScores = set_names(zVec, methodNames),
        methods = combineMethods,
        R = sig
    )
}

# K == 1 combination: every requested method trivially takes the single tuple's
# p-value, with a 1x1 aligned correlation matrix.
# @noRd
.twasZCombineSingle <- function(combineMethods, pVec, zMatrix) {
    perMethod <- set_names(
        map(combineMethods, .cipSingleMethodPval, pVec = pVec),
        combineMethods
    )
    list(
        input = list(
            nPvalsIn = 1L,
            nZScoresIn = 1L,
            nValid = 1L,
            Raligned = matrix(
                1.0,
                1L,
                1L,
                dimnames = list(rownames(zMatrix), rownames(zMatrix))
            )
        ),
        results = perMethod
    )
}

# ---- map/apply helpers (lambda-free callbacks) ---------------------------

# One method's single-tuple record (K == 1: it just takes the lone p-value).
# @noRd
.cipSingleMethodPval <- function(m, pVec) {
    list(method = m, pval = as.numeric(pVec[[1L]]))
}

#' @rdname WeightSelectionParam
#' @aliases WeightSelectionParam-class
#' @exportClass WeightSelectionParam
setClass(
    "WeightSelectionParam",
    contains = "MethodParam",
    slots = c(
        methods = "character",
        jointAxes = "character_OR_NULL",
        cutoff = "numeric",
        pvalCutoff = "numeric",
        metric = "character",
        pvalMetric = "character"
    )
)

# "best" and "all" each stand alone; anything else names methods.
setValidity("WeightSelectionParam", function(object) {
    coll <- makeAssertCollection()
    assertCharacter(
        object@methods,
        min.len = 1L,
        min.chars = 1L,
        any.missing = FALSE,
        unique = TRUE,
        .var.name = "methods",
        add = coll
    )
    if (
        length(object@methods) > 1L &&
            any(is_in(object@methods, c("best", "all")))
    ) {
        coll$push(str_c(
            "methods: \"best\" and \"all\" stand alone; to report several ",
            "methods, name them"
        ))
    }
    if (!is.null(object@jointAxes)) {
        assertSubset(
            object@jointAxes,
            .JOINT_KINDS,
            empty.ok = FALSE,
            .var.name = "jointAxes",
            add = coll
        )
    }
    problems <- coll$getMessages()
    if (length(problems) == 0L) TRUE else problems
})

#' @title TWAS Weight-Selection Settings
#' @description Which TWAS weight sets
#'   \code{\link{causalInferencePipeline}} and \code{\link{ctwasPipeline}}
#'   use for each gene (each
#'   \code{(study, context, trait)}), from the cross-validation metrics on
#'   the weights' \code{cvResult}. A gene can hold several weight sets of one
#'   method, one per kind of joint fit -- for instance mvsusie fitted across
#'   its contexts and across studies -- and each is a candidate of its own.
#' @param methods \code{"best"} (default), \code{"all"}, or the names of the
#'   methods to use. \code{"best"} takes each gene's weight set -- its method
#'   and kind of joint fit -- with the highest \code{metric} among those
#'   clearing the cutoffs whose weights are not all zero. A gene whose sets
#'   carry no \code{metric} keeps them all, as under \code{"all"}, there being
#'   nothing to rank them by; so does every gene of a source with only one
#'   set per gene, which then simply uses it. The choice depends on the
#'   weights alone, never on a GWAS study, so the two pipelines pick the same
#'   weight set for a gene.
#' @param jointAxes Optional character vector: the kinds of joint fit whose
#'   weight sets are candidates, narrowing them before \code{methods} applies.
#'   A kind is \code{"none"}, a fit of the tuple alone, or the axes the fit
#'   joined, \code{"+"}-joined in study, context, trait order:
#'   \code{"study"}, \code{"context"}, \code{"trait"}, \code{"context+trait"}
#'   and so on, as the weights' \code{jointStudies} / \code{jointContexts} /
#'   \code{jointTraits} record. For instance
#'   \code{WeightSelectionParam(methods = "mvsusie", jointAxes = "study")}
#'   keeps only mvsusie fitted across studies. \code{NULL} (default) keeps
#'   every kind.
#' @param cutoff Numeric. Minimum cross-validated r-squared a method must
#'   reach to be used, whatever \code{methods} is. Default \code{0} (no
#'   r-squared gate).
#' @param pvalCutoff Numeric. A method is used only when its
#'   \code{cvResult} p-value metric is below this, whatever \code{methods}
#'   is. Default \code{Inf} (no gate). Where a gate is set, a method without
#'   the metric fails it, unless no method of its gene has the r-squared.
#' @param metric Character. Which \code{cvResult} metric is the r-squared
#'   used for the cutoff and for ranking; typically \code{"rsq"} or
#'   \code{"adj_rsq"}. Default \code{"rsq"}.
#' @param pvalMetric Character vector of candidate \code{cvResult} metric
#'   names for the p-value gate; the first one present in a tuple's metrics
#'   is used. Default \code{c("adj_rsq_pval", "pval")}.
#' @return A \code{WeightSelectionParam} object, a \code{\link{MethodParam}}.
#' @examples
#' WeightSelectionParam(cutoff = 0.01, metric = "adj_rsq")
#' # Every method, or just the two named.
#' WeightSelectionParam(methods = "all")
#' WeightSelectionParam(methods = c("lasso", "mrash"))
#' # Only weights fitted across studies.
#' WeightSelectionParam(jointAxes = "study")
#' @importFrom checkmate assertCharacter makeAssertCollection
#' @export
WeightSelectionParam <- function(
    methods = "best",
    jointAxes = NULL,
    cutoff = 0,
    pvalCutoff = Inf,
    metric = "rsq",
    pvalMetric = c("adj_rsq_pval", "pval")
) {
    new(
        "WeightSelectionParam",
        methods = methods,
        jointAxes = jointAxes,
        cutoff = cutoff,
        pvalCutoff = pvalCutoff,
        metric = metric,
        pvalMetric = pvalMetric
    )
}

#' @rdname MrParam
#' @aliases MrParam-class
#' @exportClass MrParam
setClass(
    "MrParam",
    contains = "MethodParam",
    slots = c(
        method = "character",
        pipCutoff = "numeric",
        cpipCutoff = "numeric",
        pvalCutoff = "numeric"
    )
)

#' @title Mendelian-Randomization Settings
#' @description How \code{\link{causalInferencePipeline}} runs MR from a
#'   tuple's \code{topLoci} table.
#' @param method \code{"ivwPerVariant"} (default) filters \code{topLoci}
#'   variants by \code{pip > pipCutoff} and IVW-pools their Wald ratios.
#'   \code{"csAware"} groups variants by credible set and computes a
#'   PIP-weighted composite Wald ratio per set, using \code{cpipCutoff}.
#' @param pipCutoff Numeric. PIP threshold for a variant to be used as an
#'   instrument. Default \code{0.5}. Used only by \code{"ivwPerVariant"}.
#' @param cpipCutoff Numeric. Cumulative-PIP cutoff for retaining a credible
#'   set. Default \code{0.5}. Used only by \code{"csAware"}.
#' @param pvalCutoff Numeric. A gene's MR against a GWAS study is estimated
#'   only when one of the gene's TWAS p-values for that study is below this;
#'   otherwise its MR row's estimate is \code{NA}. A gene with no TWAS row for
#'   the study is not gated. Default \code{1} (no gate).
#' @return A \code{MrParam} object, a \code{\link{MethodParam}}.
#' @examples
#' MrParam(method = "csAware", cpipCutoff = 0.8)
#' @export
MrParam <- function(
    method = c("ivwPerVariant", "csAware"),
    pipCutoff = 0.5,
    cpipCutoff = 0.5,
    pvalCutoff = 1
) {
    method <- arg_match(method)
    .mrAssertApplicable(method, pipCutoff, cpipCutoff)
    new(
        "MrParam",
        method = method,
        pipCutoff = pipCutoff,
        cpipCutoff = cpipCutoff,
        pvalCutoff = pvalCutoff
    )
}

# Each method reads exactly one of the two PIP cutoffs. Setting the other to
# something other than its default is a request that would be ignored, so it
# fails here rather than silently doing nothing.
# @noRd
.mrAssertApplicable <- function(method, pipCutoff, cpipCutoff) {
    unused <- if (method == "csAware") {
        if (!isTRUE(all.equal(pipCutoff, 0.5))) "pipCutoff" else NULL
    } else {
        if (!isTRUE(all.equal(cpipCutoff, 0.5))) "cpipCutoff" else NULL
    }
    if (is.null(unused)) {
        return(invisible(NULL))
    }
    other <- if (method == "csAware") "cpipCutoff" else "pipCutoff"
    abort(glue(
        "MrParam: `{unused}` is not read by method = \"{method}\"; that ",
        "method uses `{other}`. Set the one its method reads, or change ",
        "`method`."
    ))
}
