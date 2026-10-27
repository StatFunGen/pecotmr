# =============================================================================
# ColocBoostResult S4 class
# -----------------------------------------------------------------------------
# The output of colocboostPipeline. One ELEMENT per confidence set (CoS),
# holding that set's member variants; the per-set statistics live in outer
# mcols and the per-variant colocalization probabilities in inner mcols.
#
# Sibling of ColocResult, NOT a reuse of it. Both share the RangedTupleList
# shape, but the statistics are different in kind and forcing colocboost into
# the coloc schema would store claims that were never computed:
#
#   * coloc is PAIRWISE -- one QTL credible set against one GWAS credible set.
#     A colocboost CoS colocalizes an arbitrary-sized SET of outcomes (verified:
#     a four-trait run yields cos_id "cos1:y1_y2_y3" over three outcomes), so
#     there is no (qtlCs, gwasCs) pair to key on. Expanding one CoS into all
#     pairs of its members would fabricate pairwise quantities.
#   * ColocResult REQUIRES PP.H0.abf..PP.H4.abf. colocboost computes no
#     five-hypothesis decomposition at all; its quantities are cos_npc,
#     min_npc_outcome, per-outcome npc and purity.
#   * The variant layers differ in meaning. SNP.PP.H4 is conditional on a pair
#     ("given these two colocalize, this is the shared causal variant"); vcp is
#     a marginal per-variant colocalization probability over the region.
#
# What the two DO share is the projection surface, so both extend
# ColocResultBase and both answer colocPairs() and colocVariants(). Each keeps
# its own LD references. A colocboost run matches every summary-statistic
# outcome against its own row's reference, so this class keeps a named list of
# them -- each reference once -- and outcomeInfo records the one every outcome
# used.
#
#   element      the member variants of one CoS
#   inner mcols  per-variant: vcp
#   outer mcols  per-CoS: cosId, analysis, gwasStudy, outcomes (a
#                CharacterList -- native variable arity), purity, cosNpc,
#                minNpcOutcome, topVariable, topVariableVcp, focalOutcome,
#                isColocalized
#   slots        ldReference (a named list), outcomeInfo, regionVcp,
#                computingTime
#
# Uncolocalized (outcome-specific) sets are elements too, flagged
# isColocalized = FALSE, rather than living in a second container: they have
# the same shape and the same views apply.
# =============================================================================

#' @include AllClasses.R RangedTupleList.R ColocResult.R
NULL

#' @title ColocBoost Result
#' @description A collection of ColocBoost results, one element per confidence
#'   set (CoS). Project it with \code{\link{colocPairs}},
#'   \code{\link{colocVariants}} or \code{\link{colocBoostOutcomes}}.
#' @slot ldReference A named list of the LD references the summary statistics
#'   were matched against, each a genotype panel, a \code{GenotypeMixture} or
#'   an \code{LdMatrix}, and each stored once however many outcomes used it.
#'   Empty when no summary statistics were analysed.
#' @slot outcomeInfo A data frame mapping each outcome \code{name} to its
#'   \code{study}, \code{context}, \code{trait}, \code{dataForm} and
#'   \code{ldReference}, the name of the \code{ldReference} entry its
#'   summary statistics were matched against (\code{NA} for an
#'   individual-level outcome, whose LD comes from its own genotypes).
#'   ColocBoost is given outcome names as bare labels, and those labels are
#'   lossy -- a trait is context-qualified only when ambiguous -- so the
#'   mapping is recorded when the names are minted rather than parsed back
#'   out later.
#' @slot regionVcp A \code{GRanges} of every variant in the analysed region
#'   carrying the marginal \code{vcp}. Collection-level: it spans the region,
#'   not the elements, so it cannot live in \code{mcols}.
#' @slot computingTime A list of per-analysis timings.
#' @seealso \code{\link{colocboostPipeline}},
#'   \code{\linkS4class{ColocResult}}
#' @export
setClass(
    "ColocBoostResult",
    contains = "ColocResultBase",
    representation(
        ldReference = "list",
        outcomeInfo = "data.frame",
        regionVcp = "ANY",
        computingTime = "list"
    ),
    prototype(
        ldReference = list(),
        regionVcp = NULL,
        computingTime = list()
    )
)

methods::setValidity("ColocBoostResult", function(object) {
    .validateColocBoostResult(object)
})

# @noRd
.cbrRequiredCols <- function() {
    c(
        "cosId",
        "analysis",
        "outcomes",
        "purity",
        "cosNpc",
        "nVariables",
        "isColocalized"
    )
}

# @noRd
#' @importFrom checkmate makeAssertCollection assertNames checkNames
#' @importFrom checkmate checkSubset
.validateColocBoostResult <- function(object) {
    coll <- makeAssertCollection()
    assertNames(
        colnames(mcols(object, use.names = FALSE)) %||% character(0),
        must.include = .cbrRequiredCols(),
        what = "colnames",
        .var.name = "mcols",
        add = coll
    )
    # The checks below read those columns; running them on an object missing
    # them reports the consequence rather than the cause.
    if (!coll$isEmpty()) {
        return(coll$getMessages())
    }
    coll$push(.cbrCheckVcpColumn(object))
    coll$push(.cbrCheckOutcomeInfo(object))
    # The outcomes, not the sets, are what match against a reference, so
    # their keys live in outcomeInfo rather than in mcols.
    coll$push(.ldRefDictProblems(
        object@ldReference,
        object@outcomeInfo$ldReference
    ))
    coll$getMessages()
}

# The per-variant layer is the point of the class, exactly as SNP.PP.H4 is for
# ColocResult: an element without it reproduces the gap this closes. The
# elements share one variant metadata table, so they all have the column or
# none does.
# @noRd
.cbrCheckVcpColumn <- function(object) {
    if (length(object) == 0L || .rtlHasVariantColumn(object, "vcp")) {
        return(NULL)
    }
    str_c(
        "element(s) ",
        str_flatten(seq_along(object), ", "),
        " have no vcp column in their variant metadata"
    )
}

# Every outcome named on an element must be resolvable in outcomeInfo, or the
# identity join silently drops rows -- which looks like "this CoS colocalized
# fewer traits" rather than like a bookkeeping bug.
# @noRd
.cbrCheckOutcomeInfo <- function(object) {
    info <- object@outcomeInfo
    if (nrow(info) == 0L) {
        return(NULL)
    }
    cols <- checkNames(
        colnames(info),
        must.include = c("name", "study", "context", "trait", "dataForm"),
        what = "colnames"
    )
    if (!isTRUE(cols)) {
        return(str_c("outcomeInfo is missing columns: ", cols))
    }
    if (length(object) == 0L) {
        return(NULL)
    }
    named <- unique(unlist(mcols(object, use.names = FALSE)$outcomes))
    resolved <- checkSubset(named, as.character(info$name))
    if (isTRUE(resolved)) {
        return(NULL)
    }
    str_c("outcome(s) not in outcomeInfo: ", resolved)
}

# ---- accessors --------------------------------------------------------------

#' @rdname show-methods
#' @export
setMethod("show", "ColocBoostResult", function(object) {
    cat(class(object), "with", length(object), "confidence set(s)\n")
    if (length(object) == 0L) {
        return(invisible(NULL))
    }
    md <- mcols(object, use.names = FALSE)
    coloc <- sum(md$isColocalized)
    cat("  colocalized :", coloc, "\n")
    cat("  outcome-only:", length(object) - coloc, "\n")
    cat("  analyses    :", str_flatten(unique(md$analysis), ", "), "\n")
    cat("  variants    :", sum(lengths(object)), "across all sets\n")
    if (coloc > 0L) {
        cat(
            "  max cos_npc :",
            format(max(md$cosNpc[md$isColocalized], na.rm = TRUE), digits = 4),
            "\n"
        )
    }
    invisible(NULL)
})

#' @rdname ldReference
#' @export
setMethod("ldReference", "ColocBoostResult", function(x) x@ldReference)

# A ColocBoostResult keys its outcomes, not its sets, to its references, so
# its update is its own: the one reference an earlier version recorded was
# the one every summary-statistic outcome was matched against.
#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "ColocBoostResult",
    function(object, ..., verbose = FALSE) {
        methods::callNextMethod(
            .cbrUpdateLdRefs(object),
            ...,
            verbose = verbose
        )
    }
)

# @noRd
.cbrUpdateLdRefs <- function(object) {
    old <- methods::slot(object, "ldReference")
    info <- .cbrWithLdKeys(object@outcomeInfo)
    if (is.list(old) && !is.object(old)) {
        return(`slot<-`(object, "outcomeInfo", value = info))
    }
    dict <- .ldRefDict(list(.asLdRef(old)))
    key <- if (length(dict$refs) == 0L) NA_character_ else dict$keys[[1L]]
    info <- mutate(
        info,
        ldReference = if_else(.data$dataForm == "sumstats", key, NA_character_)
    )
    updated <- `slot<-`(object, "ldReference", check = FALSE, value = dict$refs)
    `slot<-`(updated, "outcomeInfo", value = info)
}

# ---- construction -----------------------------------------------------------

#' @title Build a ColocBoostResult
#' @description Assemble a \code{\linkS4class{ColocBoostResult}} from the raw
#'   objects \code{colocboost::colocboost()} returns. Callers normally get one
#'   from \code{\link{colocboostPipeline}} rather than building it directly.
#' @param results A named list of \code{colocboost} objects, one per analysis
#'   (\code{xqtl_coloc}, \code{joint_gwas}, or one entry per GWAS study for
#'   \code{separate_gwas}). \code{NULL} entries -- a run that failed or was not
#'   requested -- are skipped.
#' @param analysis Character vector, parallel to \code{results}, naming the
#'   analysis each element came from.
#' @param gwasStudy Optional character vector, parallel to \code{results},
#'   naming the GWAS study for the per-study \code{separate_gwas} runs.
#' @param outcomeInfo Data frame mapping outcome \code{name} to \code{study},
#'   \code{context}, \code{trait} and \code{dataForm}, and optionally to
#'   \code{ldReference}, the name of the \code{ldReference} entry the
#'   outcome's summary statistics were matched against. Without that column
#'   every outcome gets \code{NA}.
#' @param ldReference Optional named list of the LD references the summary
#'   statistics were matched against: genotype panels (see
#'   \code{\link{readGenotypes}}), \code{GenotypeMixture} or \code{LdMatrix}
#'   objects.
#' @param computingTime Optional list of per-analysis timings.
#' @param metadata A list of free-form metadata about the collection,
#'   as \code{metadata()} holds on any Bioconductor vector. Default
#'   \code{list()}.
#' @return A \code{ColocBoostResult}. The outcome-specific (uncolocalized)
#'   sets are kept as elements flagged \code{isColocalized = FALSE}, so "no
#'   colocalization here" stays distinct from "no signal here";
#'   \code{x[x$isColocalized]} keeps the colocalized ones.
#' @examples
#' # Built from `colocboost` objects, which colocboostPipeline() supplies;
#' # colocboostResultExample is one such result, from a real run.
#' data(colocboostResultExample)
#' colocboostResultExample
#' # Its unit is a SET of colocalized outcomes rather than a pair, so one row
#' # can carry any number of them.
#' colocboostResultExample$outcomes
#' colocPairs(colocboostResultExample)$nOutcomes
#' # A run that failed or was not requested arrives as NULL and is skipped,
#' # giving an empty result rather than an error.
#' nrow(ColocBoostResult(list(NULL), "xqtl_coloc"))
#' @seealso \code{\linkS4class{ColocResult}},
#'   \code{\link{colocboostPipeline}}
#' @export
ColocBoostResult <- function(
    results,
    analysis,
    gwasStudy = NA_character_,
    outcomeInfo = NULL,
    ldReference = NULL,
    computingTime = list(),
    metadata = list()
) {
    results <- as.list(results)
    n <- length(results)
    analysis <- rep_len(as.character(analysis), n)
    gwasStudy <- rep_len(as.character(gwasStudy), n)
    if (is.null(outcomeInfo)) {
        outcomeInfo <- .cbEmptyOutcomeInfo()
    }
    rows <- map(
        seq_len(n),
        .cbrRowsFor,
        results = results,
        analysis = analysis,
        gwasStudy = gwasStudy
    ) |>
        list_flatten() |>
        unname()
    assembled <- .cbrAssemble(
        rows,
        outcomeInfo = .cbrWithLdKeys(as.data.frame(outcomeInfo)),
        regionVcp = .cbrRegionVcp(results),
        ldReference = .cbrAsLdRefList(ldReference),
        computingTime = computingTime
    )
    .withMetadata(assembled, metadata)
}

# The `ldReference` argument as the slot holds it: an empty list for NULL, or
# the named list with each entry normalised like any LD-reference input (a
# bare GenotypeHandle becomes a panel). A lone reference is refused rather
# than given a made-up name, since the name is what outcomeInfo keys on.
# @noRd
.cbrAsLdRefList <- function(ldReference) {
    if (is.null(ldReference)) {
        return(list())
    }
    if (!is.list(ldReference) || is.object(ldReference)) {
        abort(str_c(
            "`ldReference` must be NULL or a named list of LD references, ",
            "e.g. list(qtl = panel)."
        ))
    }
    map(ldReference, .asLdRef)
}

# outcomeInfo with its `ldReference` column, NA throughout when the caller's
# frame has none: no recorded reference is not the same as a shared one.
# @noRd
.cbrWithLdKeys <- function(info) {
    if (is_in("ldReference", colnames(info))) {
        return(info)
    }
    mutate(info, ldReference = NA_character_)
}

# Every set (colocalized, then uncolocalized) of one colocboost object, as a
# list of one-row records.
# @noRd
.cbrRowsFor <- function(i, results, analysis, gwasStudy) {
    res <- results[[i]]
    if (is.null(res)) {
        return(list())
    }
    c(
        .cbrColocalizedRows(res, analysis[[i]], gwasStudy[[i]]),
        .cbrUncolocalizedRows(res, analysis[[i]], gwasStudy[[i]])
    )
}

# @noRd
.cbrColocalizedRows <- function(res, analysis, gwasStudy) {
    d <- res$cos_details
    ids <- names(d$cos$cos_index)
    if (is.null(ids) || length(ids) == 0L) {
        return(list())
    }
    map(
        ids,
        .cbrOneCosRow,
        res = res,
        analysis = analysis,
        gwasStudy = gwasStudy
    )
}

# @noRd
.cbrOneCosRow <- function(id, res, analysis, gwasStudy) {
    d <- res$cos_details
    idx <- d$cos$cos_index[[id]]
    vars <- as.character(d$cos$cos_variables[[id]])
    # cos_vcp spans every variant in the region; the element is the SET, so it
    # is indexed down to the members rather than carrying the region.
    vcpAll <- d$cos_vcp[[id]]
    vcp <- if (is.null(vcpAll)) rep(NA_real_, length(idx)) else vcpAll[idx]
    summaryRow <- .cbrSummaryFor(res, id)
    list(
        variants = tibble(variant_id = vars, vcp = as.numeric(vcp)),
        meta = tibble(
            cosId = id,
            analysis = analysis,
            gwasStudy = gwasStudy,
            purity = .cbrPurity(d$cos_purity, id),
            cosNpc = .cbrNum(d$cos_npc[[id]]),
            minNpcOutcome = .cbrNum(d$cos_min_npc_outcome[[id]]),
            nVariables = length(idx),
            topVariable = .cbrTopVariable(d$cos_top_variables, id),
            topVariableVcp = .cbrNum(summaryRow$top_variable_vcp),
            # focal_outcome is `character` (the trait) when a focal outcome was
            # set and `logical FALSE` when it was not; left alone it would
            # coerce the whole column on rbind, so it is normalized to
            # character with NA meaning "no focal outcome".
            focalOutcome = .cbrFocalOutcome(summaryRow$focal_outcome),
            isColocalized = TRUE
        ),
        outcomes = as.character(d$cos_outcomes$outcome_name[[id]])
    )
}

# Outcome-specific (uncolocalized) sets, kept as elements so "no
# colocalization" stays distinguishable from "no signal".
# @noRd
.cbrUncolocalizedRows <- function(res, analysis, gwasStudy) {
    d <- res$ucos_details
    if (is.null(d)) {
        return(list())
    }
    ids <- names(d$ucos$ucos_index)
    if (is.null(ids) || length(ids) == 0L) {
        return(list())
    }
    map(
        ids,
        .cbrOneUcosRow,
        d = d,
        analysis = analysis,
        gwasStudy = gwasStudy
    )
}

# @noRd
.cbrOneUcosRow <- function(id, d, analysis, gwasStudy) {
    idx <- d$ucos$ucos_index[[id]]
    vars <- as.character(d$ucos$ucos_variables[[id]])
    list(
        variants = tibble(
            variant_id = vars,
            vcp = rep(NA_real_, length(vars))
        ),
        meta = tibble(
            cosId = id,
            analysis = analysis,
            gwasStudy = gwasStudy,
            purity = .cbrPurity(d$ucos_purity, id),
            cosNpc = NA_real_,
            minNpcOutcome = NA_real_,
            nVariables = length(idx),
            topVariable = .cbrTopVariable(d$ucos_top_variables, id),
            topVariableVcp = NA_real_,
            focalOutcome = NA_character_,
            isColocalized = FALSE
        ),
        outcomes = as.character(d$ucos_outcomes$outcome_name[[id]])
    )
}

# The cos_summary row for one set, or an empty tibble when absent.
# @noRd
.cbrSummaryFor <- function(res, id) {
    raw <- res$cos_summary
    if (is.null(raw) || nrow(raw) == 0L || !is_in("cos_id", colnames(raw))) {
        return(tibble(
            top_variable = NA_character_,
            top_variable_vcp = NA_real_,
            focal_outcome = NA
        ))
    }
    s <- as.data.frame(raw, stringsAsFactors = FALSE)
    hit <- which(as.character(s$cos_id) == id)
    if (length(hit) == 0L) {
        return(tibble(
            top_variable = NA_character_,
            top_variable_vcp = NA_real_,
            focal_outcome = NA
        ))
    }
    s[hit[[1L]], , drop = FALSE]
}

# Purity of one set. colocboost reports purity as a list of three SQUARE
# matrices over sets (min / max / median absolute correlation); a set's own
# purity is the diagonal entry, and the off-diagonal is between-set
# correlation. Indexing the list by id -- which is what it looks like it wants
# -- silently yields NA.
#
# min_abs_cor is used because that is what `minAbsCorr` means everywhere else
# in the package (fineMappingPipeline, colocCredibleSets).
# @noRd
.cbrPurity <- function(purity, id) {
    m <- purity$min_abs_cor
    if (is.null(m) || is.null(rownames(m)) || !is_in(id, rownames(m))) {
        return(NA_real_)
    }
    as.numeric(m[id, id])
}

# Top variable of one set. Both cos_top_variables and ucos_top_variables are
# data frames whose ROWNAMES are the set ids, not per-id lists.
# @noRd
.cbrTopVariable <- function(tv, id) {
    if (is.null(tv) || is.null(rownames(tv)) || !is_in(id, rownames(tv))) {
        return(NA_character_)
    }
    as.character(tv[id, "top_variables"])
}

# @noRd
.cbrNum <- function(v) {
    if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v)[[1L]]
}

# `FALSE` means "no focal outcome was set"; anything else is the trait name.
# @noRd
.cbrFocalOutcome <- function(v) {
    if (is.null(v) || length(v) == 0L || is.logical(v)) {
        return(NA_character_)
    }
    as.character(v)[[1L]]
}

# The region-wide marginal vcp, as a GRanges. Taken from the first run that
# carries one: every analysis in a call scores the same region, so they agree.
# @noRd
.cbrRegionVcp <- function(results) {
    for (res in results) {
        if (is.null(res) || is.null(res$vcp)) {
            next
        }
        ids <- names(res$vcp)
        if (is.null(ids)) {
            next
        }
        bare <- .variantIdsToGRanges(ids, what = "colocboost vcp names")
        gr <- S4Vectors::`mcols<-`(
            bare,
            value = cbind(
                mcols(bare, use.names = FALSE),
                S4Vectors::DataFrame(vcp = as.numeric(res$vcp))
            )
        )
        return(gr)
    }
    NULL
}

# @noRd
.cbrAssemble <- function(
    rows,
    outcomeInfo,
    regionVcp,
    ldReference,
    computingTime
) {
    elements <- map(rows, .cbrElementFor)
    grl <- S4Vectors::`mcols<-`(
        .rtlBindElements(elements),
        value = .cbrMcolsFor(rows)
    )
    new(
        "ColocBoostResult",
        .rtlNamed(grl),
        ldReference = ldReference,
        outcomeInfo = outcomeInfo,
        regionVcp = regionVcp,
        computingTime = computingTime
    )
}

# @noRd
.cbrElementFor <- function(row) {
    v <- row$variants
    bare <- .variantIdsToGRanges(
        as.character(v$variant_id),
        what = "colocboost variant name"
    )
    # Appended, not assigned: the A1 / A2 columns .variantIdsToGRanges()
    # attaches are the variant identity, and overwriting mcols would leave the
    # element unable to name its own variants.
    gr <- S4Vectors::`mcols<-`(
        bare,
        value = cbind(
            mcols(bare, use.names = FALSE),
            S4Vectors::DataFrame(vcp = as.numeric(v$vcp))
        )
    )
    gr
}

# @noRd
.cbrMcolsFor <- function(rows) {
    # A CharacterList column, so a set of any size lives in one row instead of
    # being flattened to a delimited string a caller has to re-split.
    if (length(rows) == 0L) {
        return(cbind(
            S4Vectors::DataFrame(.cbrEmptyMeta()),
            S4Vectors::DataFrame(outcomes = IRanges::CharacterList())
        ))
    }
    cbind(
        S4Vectors::DataFrame(
            bind_rows(map(rows, "meta")),
            check.names = FALSE
        ),
        S4Vectors::DataFrame(
            outcomes = IRanges::CharacterList(map(rows, "outcomes"))
        )
    )
}

# @noRd
.cbrEmptyMeta <- function() {
    tibble(
        cosId = character(0),
        analysis = character(0),
        gwasStudy = character(0),
        purity = numeric(0),
        cosNpc = numeric(0),
        minNpcOutcome = numeric(0),
        nVariables = integer(0),
        topVariable = character(0),
        topVariableVcp = numeric(0),
        focalOutcome = character(0),
        isColocalized = logical(0)
    )
}

# ---- views ------------------------------------------------------------------

#' @rdname colocViews
#' @export
setMethod(
    "colocPairs",
    "ColocBoostResult",
    function(x) {
        .asViewDataFrame(.cbrPairsTable(x = x))
    }
)

# colocPairs() on a ColocBoostResult as a tibble.
# @noRd
.cbrPairsTable <- function(x) {
    # Named "pairs" for symmetry with ColocResult, but a ColocBoost row is a
    # SET, not a pair: `outcomes` holds however many outcomes colocalized.
    md <- mcols(x, use.names = FALSE)
    if (is.null(md) || length(x) == 0L) {
        return(tibble())
    }
    # `drop = FALSE` as in the sibling views: a single kept column would
    # arrive as a vector and be named after the subscript text.
    flat <- as.data.frame(
        md[, setdiff(colnames(md), "outcomes"), drop = FALSE]
    )
    mutate(
        as_tibble(flat, .name_repair = "minimal"),
        outcomes = map_chr(as.list(md$outcomes), str_flatten, collapse = "; "),
        nOutcomes = lengths(md$outcomes)
    )
}

#' @rdname colocViews
#' @export
setMethod(
    "colocVariants",
    "ColocBoostResult",
    function(x) {
        .withGenome(.variantTableGRanges(.cbrLongVariants(x)), .oneBuild(x))
    }
)

#' @rdname colocViews
#' @export
setMethod(
    "colocPooledVariants",
    "ColocBoostResult",
    function(x) {
        long <- .cbrLongVariants(x)
        if (nrow(long) == 0L) {
            return(.withGenome(.variantTableGRanges(long), .oneBuild(x)))
        }
        # Pooling across sets is noisy-OR only: unlike coloc, ColocBoost's sets
        # are not per-block slices of one signal, so there is no
        # mutually-exclusive sum stage to apply first.
        pooledTbl <- long |>
            group_by(across(all_of(c("analysis", "variant_id")))) |>
            summarise(vcp = .crNoisyOr(.data$vcp), .groups = "drop")
        .withGenome(.variantTableGRanges(pooledTbl), .oneBuild(x))
    }
)

# One row per (set, variant).
# @noRd
.cbrLongVariants <- function(x) {
    if (length(x) == 0L) {
        return(tibble())
    }
    md <- as.data.frame(
        mcols(x, use.names = FALSE)[,
            setdiff(colnames(mcols(x, use.names = FALSE)), "outcomes"),
            drop = FALSE
        ]
    )
    n <- lengths(x)
    setIdx <- rep(seq_len(length(x)), n)
    flat <- unlist(x, use.names = FALSE)
    mutate(
        as_tibble(md[setIdx, , drop = FALSE], .name_repair = "minimal"),
        variant_id = .grVariantIds(flat),
        vcp = as.numeric(mcols(flat, use.names = FALSE)$vcp)
    )
}

#' @title ColocBoost Outcome View
#' @description One row per (confidence set, colocalized outcome), joined to
#'   the identity the outcome was built from. This is what the
#'   \code{outcomeInfo} slot exists for: ColocBoost is handed outcome names as
#'   bare labels, so without the recorded mapping a caller can only get the
#'   label back, never the study / context / trait behind it.
#' @param x A \code{\linkS4class{ColocBoostResult}}.
#' @param ... Ignored.
#' @return A \code{DataFrame} with one row per (set, outcome), the
#'   set-level columns repeated, plus \code{outcome}, \code{study},
#'   \code{context},
#'   \code{trait}, \code{dataForm} and \code{ldReference} (the name of the
#'   \code{\link{ldReference}} entry the outcome's summary statistics were
#'   matched against).
#' @examples
#' # See colocboostPipeline(); a ColocBoostResult carries the mapping it needs.
#' methods::existsMethod("colocBoostOutcomes", "ColocBoostResult")
#' @export
setGeneric("colocBoostOutcomes", function(x, ...) {
    standardGeneric("colocBoostOutcomes")
})

#' @rdname colocBoostOutcomes
#' @export
setMethod(
    "colocBoostOutcomes",
    "ColocBoostResult",
    function(x) {
        .asViewDataFrame(.cbrOutcomesTable(x = x))
    }
)

# colocBoostOutcomes() as a tibble.
# @noRd
.cbrOutcomesTable <- function(x) {
    if (length(x) == 0L) {
        return(tibble())
    }
    md <- mcols(x, use.names = FALSE)
    outcomes <- as.list(md$outcomes)
    reps <- lengths(outcomes)
    flat <- as.data.frame(md[, setdiff(colnames(md), "outcomes"), drop = FALSE])
    out <- mutate(
        as_tibble(
            flat[rep(seq_len(nrow(flat)), reps), , drop = FALSE],
            .name_repair = "minimal"
        ),
        outcome = unlist(outcomes, use.names = FALSE)
    )
    info <- x@outcomeInfo
    if (nrow(info) == 0L) {
        return(out)
    }
    left_join(
        out,
        as_tibble(info, .name_repair = "minimal"),
        by = c(outcome = "name")
    )
}

#' @title Per-Analysis Timings
#' @description The wall-clock timings colocboostPipeline recorded for each
#'   analysis it ran.
#' @param x A \code{\linkS4class{ColocBoostResult}}.
#' @param ... Ignored.
#' @return A \code{DataFrame}, one row per analysis that ran:
#'   \code{analysis} (\code{xqtl_coloc}, \code{joint_gwas} or
#'   \code{separate_gwas}), \code{nRuns} (1, or one per GWAS study for
#'   \code{separate_gwas}) and \code{seconds}, their total wall-clock time.
#'   No rows when none were recorded.
#' @examples
#' methods::existsMethod("computingTime", "ColocBoostResult")
#' @export
setGeneric("computingTime", function(x, ...) {
    standardGeneric("computingTime")
})

#' @rdname computingTime
#' @export
setMethod("computingTime", "ColocBoostResult", function(x) {
    .cbrTimingTable(x@computingTime)
})

# The recorded timings, one row per analysis that ran.
# @noRd
.cbrTimingTable <- function(timing) {
    runs <- compact(timing$Analysis %||% list())
    S4Vectors::DataFrame(
        analysis = names(runs) %||% character(0),
        nRuns = map_int(unname(runs), .cbrTimingRuns),
        seconds = map_dbl(unname(runs), .cbrTimingSeconds)
    )
}

# A per-study sweep records its run count; every other analysis is one run.
# @noRd
.cbrTimingRuns <- function(t) {
    if (is.list(t)) as.integer(t$n_studies %||% NA_integer_) else 1L
}

# @noRd
.cbrTimingSeconds <- function(t) {
    total <- if (is.list(t)) t$total else t
    if (is.null(total)) {
        return(NA_real_)
    }
    if (inherits(total, "difftime")) {
        return(as.numeric(total, units = "secs"))
    }
    as.numeric(total)
}

#' @title Region-Wide Variant Colocalization Probabilities
#' @description The marginal \code{vcp} for every variant in the analysed
#'   region, not just the confidence-set members. Distinct from
#'   \code{\link{colocVariants}}, which reports the per-set \code{vcp} of
#'   members only.
#' @param x A \code{\linkS4class{ColocBoostResult}}.
#' @param ... Ignored.
#' @return A \code{GRanges} carrying a \code{vcp} column, or \code{NULL}.
#' @examples
#' methods::existsMethod("regionVcp", "ColocBoostResult")
#' @export
setGeneric("regionVcp", function(x, ...) standardGeneric("regionVcp"))

#' @rdname regionVcp
#' @export
setMethod("regionVcp", "ColocBoostResult", function(x) x@regionVcp)
