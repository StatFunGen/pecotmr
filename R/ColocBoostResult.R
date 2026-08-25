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
# What the two DO share is the LD panel they were computed against and the
# projection surface, so both extend ColocResultBase -- which carries the
# ldSketch slot and its accessor -- and both answer getColocPairs() and
# getColocVariants().
#
#   element      the member variants of one CoS
#   inner mcols  per-variant: vcp
#   outer mcols  per-CoS: cosId, analysis, gwasStudy, outcomes (a
#                CharacterList -- native variable arity), purity, cosNpc,
#                minNpcOutcome, topVariable, topVariableVcp, focalOutcome,
#                isColocalized
#   slots        ldSketch, outcomeInfo, regionVcp, computingTime
#
# Uncolocalized (outcome-specific) sets are elements too, flagged
# isColocalized = FALSE, rather than living in a second container: they have
# the same shape and the same views apply.
# =============================================================================

#' @include AllClasses.R RangedTupleList.R ColocResult.R
NULL

#' @title ColocBoost Result
#' @description A collection of ColocBoost results, one element per confidence
#'   set (CoS). Project it with \code{\link{getColocPairs}},
#'   \code{\link{getColocVariants}} or \code{\link{getColocBoostOutcomes}}.
#' @slot outcomeInfo A data frame mapping each outcome \code{name} to its
#'   \code{study}, \code{context}, \code{trait} and \code{dataForm}. ColocBoost
#'   is given outcome names as bare labels, and those labels are lossy -- a
#'   trait is context-qualified only when ambiguous -- so the mapping is
#'   recorded when the names are minted rather than parsed back out later.
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
        outcomeInfo = "data.frame",
        regionVcp = "ANY",
        computingTime = "list"
    ),
    prototype(
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
.validateColocBoostResult <- function(object) {
    errors <- .cbrCheckRequiredCols(object)
    if (length(errors) == 0L) {
        errors <- c(
            .cbrCheckVcpColumn(object),
            .cbrCheckOutcomeInfo(object)
        )
    }
    if (length(errors) == 0L) TRUE else errors
}

# @noRd
.cbrCheckRequiredCols <- function(object) {
    md <- mcols(object, use.names = FALSE)
    have <- if (is.null(md)) character(0) else colnames(md)
    missingCols <- setdiff(.cbrRequiredCols(), have)
    if (length(missingCols) > 0L) {
        return(str_c("missing columns: ", str_flatten(missingCols, ", ")))
    }
    NULL
}

# The per-variant layer is the point of the class, exactly as SNP.PP.H4 is for
# ColocResult: an element without it reproduces the gap this closes.
# @noRd
.cbrCheckVcpColumn <- function(object) {
    if (length(object) == 0L) {
        return(NULL)
    }
    bad <- which(!map_lgl(as.list(object), .cbrHasVcp))
    if (length(bad) == 0L) {
        return(NULL)
    }
    str_c(
        "element(s) ",
        str_flatten(bad, ", "),
        " have no vcp column in their variant metadata"
    )
}

# @noRd
.cbrHasVcp <- function(g) {
    md <- mcols(g, use.names = FALSE)
    !is.null(md) && is_in("vcp", colnames(md))
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
    missingCols <- setdiff(
        c("name", "study", "context", "trait", "dataForm"),
        colnames(info)
    )
    if (length(missingCols) > 0L) {
        return(str_c(
            "outcomeInfo is missing columns: ",
            str_flatten(missingCols, ", ")
        ))
    }
    if (length(object) == 0L) {
        return(NULL)
    }
    named <- unique(unlist(mcols(object, use.names = FALSE)$outcomes))
    unknown <- setdiff(named, as.character(info$name))
    if (length(unknown) == 0L) {
        return(NULL)
    }
    str_c(
        "outcome(s) not present in outcomeInfo: ",
        str_flatten(utils::head(unknown, 5L), ", ")
    )
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
#'   \code{context}, \code{trait} and \code{dataForm}.
#' @param ldSketch Optional genotype panel (see
#'   \code{\link{readGenotypes}}).
#' @param computingTime Optional list of per-analysis timings.
#' @param includeUncolocalized Keep the outcome-specific (uncolocalized) sets
#'   as elements flagged \code{isColocalized = FALSE}. Default \code{TRUE}:
#'   dropping them at write time would make "no colocalization here" and "no
#'   signal here" indistinguishable downstream.
#' @return A \code{ColocBoostResult}.
#' @examples
#' # Built from `colocboost` objects, which colocboostPipeline() supplies;
#' # colocboostResultExample is one such result, from a real run.
#' data(colocboostResultExample)
#' colocboostResultExample
#' # Its unit is a SET of colocalized outcomes rather than a pair, so one row
#' # can carry any number of them.
#' colocboostResultExample$outcomes
#' getColocPairs(colocboostResultExample)$nOutcomes
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
    ldSketch = NULL,
    computingTime = list(),
    includeUncolocalized = TRUE
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
        gwasStudy = gwasStudy,
        includeUncolocalized = includeUncolocalized
    )
    rows <- unlist(rows, recursive = FALSE, use.names = FALSE)
    .cbrAssemble(
        rows,
        outcomeInfo = as.data.frame(outcomeInfo),
        regionVcp = .cbrRegionVcp(results),
        ldSketch = .asLdSketch(ldSketch),
        computingTime = computingTime
    )
}

# Every set (colocalized, then uncolocalized) of one colocboost object, as a
# list of one-row records.
# @noRd
.cbrRowsFor <- function(i, results, analysis, gwasStudy, includeUncolocalized) {
    res <- results[[i]]
    if (is.null(res)) {
        return(list())
    }
    out <- .cbrColocalizedRows(res, analysis[[i]], gwasStudy[[i]])
    if (isTRUE(includeUncolocalized)) {
        out <- c(out, .cbrUncolocalizedRows(res, analysis[[i]], gwasStudy[[i]]))
    }
    out
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
    s <- res$cos_summary
    if (is.null(s) || nrow(s) == 0L || !is_in("cos_id", colnames(s))) {
        return(tibble(
            top_variable = NA_character_,
            top_variable_vcp = NA_real_,
            focal_outcome = NA
        ))
    }
    s <- as.data.frame(s, stringsAsFactors = FALSE)
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
# in the package (fineMappingPipeline, getColocCredibleSets).
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

# @noRd

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
        gr <- .variantIdsToGRanges(ids, what = "colocboost vcp names")
        mcols(gr) <- cbind(
            mcols(gr, use.names = FALSE),
            S4Vectors::DataFrame(vcp = as.numeric(res$vcp))
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
    ldSketch,
    computingTime
) {
    elements <- map(rows, .cbrElementFor)
    grl <- GenomicRanges::GRangesList(elements)
    mcols(grl) <- .cbrMcolsFor(rows)
    obj <- new(
        "ColocBoostResult",
        grl,
        ldSketch = .asLdSketch(ldSketch),
        outcomeInfo = outcomeInfo,
        regionVcp = regionVcp,
        computingTime = computingTime
    )
    validObject(obj)
    obj
}

# @noRd
.cbrElementFor <- function(row) {
    v <- row$variants
    gr <- .variantIdsToGRanges(
        as.character(v$variant_id),
        what = "colocboost variant name"
    )
    # Appended, not assigned: the A1 / A2 columns .variantIdsToGRanges()
    # attaches are the variant identity, and overwriting mcols would leave the
    # element unable to name its own variants.
    mcols(gr) <- cbind(
        mcols(gr, use.names = FALSE),
        S4Vectors::DataFrame(vcp = as.numeric(v$vcp))
    )
    gr
}

# @noRd
.cbrMcolsFor <- function(rows) {
    if (length(rows) == 0L) {
        md <- S4Vectors::DataFrame(.cbrEmptyMeta())
        md$outcomes <- IRanges::CharacterList()
        return(md)
    }
    meta <- bind_rows(map(rows, "meta"))
    md <- S4Vectors::DataFrame(meta, check.names = FALSE)
    # A CharacterList column, so a set of any size lives in one row instead of
    # being flattened to a delimited string a caller has to re-split.
    md$outcomes <- IRanges::CharacterList(map(rows, "outcomes"))
    md
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
setMethod("getColocPairs", "ColocBoostResult", function(x, ...) {
    # Named "pairs" for symmetry with ColocResult, but a ColocBoost row is a
    # SET, not a pair: `outcomes` holds however many outcomes colocalized.
    md <- mcols(x, use.names = FALSE)
    if (is.null(md) || length(x) == 0L) {
        return(tibble())
    }
    flat <- as.data.frame(md[, setdiff(colnames(md), "outcomes")])
    out <- as_tibble(flat, .name_repair = "minimal")
    out$outcomes <- map_chr(as.list(md$outcomes), str_flatten, collapse = "; ")
    out$nOutcomes <- lengths(md$outcomes)
    out
})

#' @rdname colocViews
#' @export
setMethod(
    "getColocVariants",
    "ColocBoostResult",
    function(x, pooled = FALSE, ...) {
        long <- .cbrLongVariants(x)
        if (!isTRUE(pooled) || nrow(long) == 0L) {
            return(long)
        }
        # Pooling across sets is noisy-OR only: unlike coloc, ColocBoost's sets
        # are not per-block slices of one signal, so there is no
        # mutually-exclusive sum stage to apply first.
        long |>
            group_by(across(all_of(c("analysis", "variant_id")))) |>
            summarise(vcp = .crNoisyOr(.data$vcp), .groups = "drop")
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
    out <- as_tibble(md[setIdx, , drop = FALSE], .name_repair = "minimal")
    out$variant_id <- .grVariantIds(flat)
    out$vcp <- as.numeric(mcols(flat, use.names = FALSE)$vcp)
    out
}

#' @title ColocBoost Outcome View
#' @description One row per (confidence set, colocalized outcome), joined to
#'   the identity the outcome was built from. This is what the
#'   \code{outcomeInfo} slot exists for: ColocBoost is handed outcome names as
#'   bare labels, so without the recorded mapping a caller can only get the
#'   label back, never the study / context / trait behind it.
#' @param x A \code{\linkS4class{ColocBoostResult}}.
#' @param ... Ignored.
#' @return A tibble with one row per (set, outcome) and the set-level columns
#'   repeated, plus \code{outcome}, \code{study}, \code{context},
#'   \code{trait} and \code{dataForm}.
#' @examples
#' # See colocboostPipeline(); a ColocBoostResult carries the mapping it needs.
#' methods::existsMethod("getColocBoostOutcomes", "ColocBoostResult")
#' @export
setGeneric("getColocBoostOutcomes", function(x, ...) {
    standardGeneric("getColocBoostOutcomes")
})

#' @rdname getColocBoostOutcomes
#' @export
setMethod("getColocBoostOutcomes", "ColocBoostResult", function(x, ...) {
    if (length(x) == 0L) {
        return(tibble())
    }
    md <- mcols(x, use.names = FALSE)
    outcomes <- as.list(md$outcomes)
    reps <- lengths(outcomes)
    flat <- as.data.frame(md[, setdiff(colnames(md), "outcomes"), drop = FALSE])
    out <- as_tibble(
        flat[rep(seq_len(nrow(flat)), reps), , drop = FALSE],
        .name_repair = "minimal"
    )
    out$outcome <- unlist(outcomes, use.names = FALSE)
    info <- x@outcomeInfo
    if (nrow(info) == 0L) {
        return(out)
    }
    left_join(
        out,
        as_tibble(info, .name_repair = "minimal"),
        by = c(outcome = "name")
    )
})

# ---- coercion ---------------------------------------------------------------

#' @title Coerce a ColocBoostResult to a data frame
#' @description Returns the set-level view. New code should call the view it
#'   wants (\code{\link{colocViews}},
#'   \code{\link{getColocBoostOutcomes}}).
#' @param x A \code{ColocBoostResult}.
#' @param row.names,optional Ignored; present for generic compatibility.
#' @param ... Ignored.
#' @return A data frame with one row per confidence set.
#' @examples
#' methods::existsMethod("as.data.frame", "ColocBoostResult")
#' @export
setMethod(
    "as.data.frame",
    "ColocBoostResult",
    function(x, row.names = NULL, optional = FALSE, ...) {
        as.data.frame(getColocPairs(x))
    }
)

#' @title Per-Analysis Timings
#' @description The wall-clock timings colocboostPipeline recorded for each
#'   analysis it ran.
#' @param x A \code{\linkS4class{ColocBoostResult}}.
#' @param ... Ignored.
#' @return A list of timings, empty when none were recorded.
#' @examples
#' methods::existsMethod("getComputingTime", "ColocBoostResult")
#' @export
setGeneric("getComputingTime", function(x, ...) {
    standardGeneric("getComputingTime")
})

#' @rdname getComputingTime
#' @export
setMethod("getComputingTime", "ColocBoostResult", function(x, ...) {
    x@computingTime
})

#' @title Region-Wide Variant Colocalization Probabilities
#' @description The marginal \code{vcp} for every variant in the analysed
#'   region, not just the confidence-set members. Distinct from
#'   \code{\link{getColocVariants}}, which reports the per-set \code{vcp} of
#'   members only.
#' @param x A \code{\linkS4class{ColocBoostResult}}.
#' @param ... Ignored.
#' @return A \code{GRanges} carrying a \code{vcp} column, or \code{NULL}.
#' @examples
#' methods::existsMethod("getRegionVcp", "ColocBoostResult")
#' @export
setGeneric("getRegionVcp", function(x, ...) standardGeneric("getRegionVcp"))

#' @rdname getRegionVcp
#' @export
setMethod("getRegionVcp", "ColocBoostResult", function(x, ...) x@regionVcp)
