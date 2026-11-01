# =============================================================================
# ColocResult S4 class
# -----------------------------------------------------------------------------
# The output of colocPipeline. One ELEMENT per scored group -- a QTL row
# against a GWAS row in one block, the credible-set pairs one coloc.bf_bf call
# scores -- holding the group's aligned variants once. Every pair of a group
# is scored over the same aligned variants, so they are stored once rather
# than once per pair.
#
# Why a class and not the flat tibble it replaces: the pair / variant /
# credible-set / gene "tiers" are not four shapes of table, they are four VIEWS
# over one representation. Storing any single granularity forces the coarser
# statistics to be repeated down the finer rows, and makes the finer ones
# unrecoverable -- which is exactly how the per-variant layer got lost when
# process_coloc_results() was dropped.
#
#   element      the aligned variants of one scored group (A1, A2)
#   outer mcols  per group: identity tuple of both sides, blockId, the LD
#                reference key, and `pairs` -- one table per group, one row
#                per credible-set pair: qtlCs, gwasCs, nSnps, PP.H0..PP.H4,
#                retained alpha mass (both sides), enrichment / p12Used, and
#                the pair's SNP.PP.H4, aligned to the group's variants
#   slots        ldReference, for the credible-set view's purity calculation
#
# Nothing is filtered at write time (spec 5.3): every testable pair and every
# variant is kept, and all narrowing is a view decision. Purity in particular
# depends on the coverage the caller asks for, so it cannot be precomputed
# without freezing coverage at construction.
#
# Column naming follows the package line: pecotmr-owned metadata is camelCase
# (gwasStudy, blockId, nSnps), upstream names are preserved verbatim, so
# coloc's own PP.H0.abf..PP.H4.abf and SNP.PP.H4 stay exactly as coloc writes
# them -- downstream code greps for those.
# =============================================================================

#' @include AllClasses.R RangedTupleList.R
NULL

#' @title Colocalization Result Base
#' @description Virtual base shared by \code{\linkS4class{ColocResult}} (from
#'   \code{coloc}) and \code{\linkS4class{ColocBoostResult}} (from
#'   \code{colocboost}).
#'
#'   The two store genuinely different statistics -- coloc is pairwise and
#'   reports \code{PP.H0}--\code{PP.H4}, while ColocBoost reports a
#'   normalized posterior over an arbitrary-sized set of outcomes -- so neither
#'   can be expressed in the other's schema. What they share is the shape (one
#'   element per colocalization unit, its variants inside) and therefore the
#'   projection surface: both answer \code{\link{colocPairs}} and
#'   \code{\link{colocVariants}}, so code that only needs "which variants
#'   colocalized, and how strongly" works against either.
#'
#'   Each subclass carries its own LD reference, and both answer
#'   \code{\link{ldReference}}: a \code{ColocResult} holds the one reference
#'   its pairs were computed against, a \code{ColocBoostResult} one per
#'   reference its summary statistics were matched against.
#' @export
setClass("ColocResultBase", contains = c("VIRTUAL", "RangedTupleList"))

#' @title Colocalization Result
#' @description A collection of colocalization results, one element per scored
#'   group: a QTL row against a GWAS row in one block, whose credible-set
#'   pairs coloc scores together over the same aligned variants. Build the
#'   four views with \code{\link{colocPairs}} (one row per pair),
#'   \code{\link{colocVariants}}, \code{\link{colocCredibleSets}} and
#'   \code{\link{colocGenes}}.
#' @details Each element holds its group's variants once. Its identity and
#'   block are columns of \code{mcols()}, and its pairs the \code{pairs}
#'   column: one table per group, one row per credible-set pair, with the
#'   pair's \code{SNP.PP.H4} aligned to the group's variants. So
#'   \code{length()} counts groups; \code{colocPairs()} has a row per pair.
#'   The \code{ldReference} is what \code{colocCredibleSets()} uses to
#'   compute credible-set purity, each group's from its own reference.
#' @slot ldReference A named list of the LD references the pairs were computed
#'   against, each a genotype panel, a \code{GenotypeMixture} or an
#'   \code{LdMatrix}, stored once apiece; the \code{ldReference} column of
#'   \code{mcols} names the entry each group used (\code{NA} for none).
#' @seealso \code{\link{colocPipeline}}
#' @export
setClass(
    "ColocResult",
    contains = "ColocResultBase",
    representation(ldReference = "list"),
    prototype(ldReference = list())
)

#' @rdname ldReference
#' @export
setMethod("ldReference", "ColocResult", function(x) .ldRefsInUse(x))

methods::setValidity("ColocResult", function(object) {
    .validateColocResult(object)
})

# @noRd
#' @importFrom checkmate makeAssertCollection assertNames checkNames
.validateColocResult <- function(object) {
    coll <- makeAssertCollection()
    assertNames(
        colnames(mcols(object, use.names = FALSE)) %||% character(0),
        must.include = .crRequiredCols(),
        what = "colnames",
        .var.name = "mcols",
        add = coll
    )
    # The checks below read those columns; running them on an object missing
    # them reports the consequence rather than the cause.
    if (!coll$isEmpty()) {
        return(coll$getMessages())
    }
    coll$push(.ldRefDictProblems(object@ldReference, .ldRefKeys(object)))
    coll$push(.crCheckPairColumns(object))
    coll$getMessages()
}

# The columns every pair of a group shares: both sides' identity and the
# block.
# @noRd
.crGroupCols <- function() {
    c(.crGeneCols(), "blockId")
}

# @noRd
.crRequiredCols <- function() {
    c(.crGroupCols(), "pairs")
}

# @noRd
.crPpCols <- function() {
    str_c("PP.H", 0:4, ".abf")
}

# The columns each pair's row must carry. The posteriors are a distinct
# contract from the identity columns, so a missing one keeps its own wording.
# The per-variant layer is the whole point of the class, so a pair without
# SNP.PP.H4 -- or with one not aligned to its group's variants -- is rejected
# rather than tolerated.
# @noRd
.crCheckPairColumns <- function(object) {
    pairs <- mcols(object, use.names = FALSE)$pairs
    if (!methods::is(pairs, "SplitDataFrameList")) {
        return("`pairs` must be a SplitDataFrameList, one table per element")
    }
    flat <- unlist(pairs, use.names = FALSE)
    res <- checkNames(
        colnames(flat),
        must.include = c(.crPpCols(), "qtlCs", "gwasCs", "SNP.PP.H4"),
        what = "colnames"
    )
    if (!isTRUE(res)) {
        return(str_c("missing posterior columns: ", res))
    }
    .crCheckPairVariants(object, flat$SNP.PP.H4, elementNROWS(pairs))
}

# @noRd
.crCheckPairVariants <- function(object, snpPp, pairsPerGroup) {
    bad <- which(lengths(snpPp) != rep(lengths(object), pairsPerGroup))
    if (length(bad) == 0L) {
        return(NULL)
    }
    str_c(
        "pair(s) ",
        str_flatten(bad, ", "),
        " have SNP.PP.H4 not aligned to their element's variants"
    )
}

# ---- accessors --------------------------------------------------------------

#' @rdname show-methods
#' @export
setMethod("show", "ColocResult", function(object) {
    pairs <- .crPairsTable(object)
    cat(
        class(object),
        "with",
        nrow(pairs),
        "colocalized pair(s) in",
        length(object),
        "scored group(s)\n"
    )
    if (length(object) == 0L) {
        return(invisible(NULL))
    }
    cat("  studies     :", str_flatten(unique(pairs$study), ", "), "\n")
    cat("  paired with :", str_flatten(unique(pairs$gwasStudy), ", "), "\n")
    cat("  variants    :", sum(lengths(object)), "across all groups\n")
    cat(
        "  max PP.H4   :",
        format(max(pairs$PP.H4.abf, na.rm = TRUE), digits = 4),
        "\n"
    )
    invisible(NULL)
})

# ---- construction -----------------------------------------------------------

#' @title Build a ColocResult
#' @description Assemble a \code{\linkS4class{ColocResult}} from a pair-level
#'   table and the per-pair variant tables that go with it. Callers normally
#'   get one from \code{\link{colocPipeline}} rather than building it directly.
#' @details Pairs that share their identity columns, block, LD reference and
#'   variants are one scored group, an element of the result; the group's
#'   variants are stored once.
#' @param pairs A data frame with one row per tested pair, carrying at least
#'   the identity columns of both sides (\code{study}, \code{context},
#'   \code{trait}, \code{method} and \code{gwasStudy}, \code{gwasContext},
#'   \code{gwasTrait}, \code{gwasMethod}), \code{blockId}, \code{qtlCs},
#'   \code{gwasCs}, \code{nSnps} and \code{PP.H0.abf} through
#'   \code{PP.H4.abf}. A GWAS side has no context or trait axis, so those two
#'   columns are \code{NA} for it.
#' @param variants A list, parallel to \code{pairs}' rows, of per-pair data
#'   frames with a \code{variant_id} column and a \code{SNP.PP.H4} column,
#'   all with the same columns.
#' @param ldReference Optional LD reference -- a genotype panel (see
#'   \code{\link{readGenotypes}}), a \code{GenotypeMixture} or an
#'   \code{LdMatrix} -- for every pair, or a list with one per pair; used by
#'   \code{\link{colocCredibleSets}} to compute purity. Each distinct
#'   reference is stored once.
#' @param metadata A list of free-form metadata about the collection,
#'   as \code{metadata()} holds on any Bioconductor vector. Default
#'   \code{list()}.
#' @return A \code{ColocResult}.
#' @examples
#' pairs <- data.frame(
#'     study = "s1", context = "c1", trait = "g1", method = "susie",
#'     gwasStudy = "G1", gwasContext = NA_character_,
#'     gwasTrait = NA_character_, gwasMethod = "susie",
#'     blockId = "chr1_1_1000",
#'     qtlCs = 1L, gwasCs = 1L, nSnps = 2L,
#'     PP.H0.abf = 0.1, PP.H1.abf = 0.1, PP.H2.abf = 0.1,
#'     PP.H3.abf = 0.1, PP.H4.abf = 0.6
#' )
#' variants <- list(data.frame(
#'     variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
#'     SNP.PP.H4 = c(0.7, 0.3)
#' ))
#' ColocResult(pairs, variants)
#' @export
ColocResult <- function(
    pairs,
    variants,
    ldReference = NULL,
    metadata = list()
) {
    pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
    if (!is.list(variants) || length(variants) != nrow(pairs)) {
        msg <- glue(
            "`variants` must be a list running parallel to `pairs`' rows ",
            "(got {length(variants)} vs {nrow(pairs)})."
        )
        abort(msg)
    }
    tables <- .crVariantTables(variants)
    dict <- .ldRefDict(.ldRefPerTuple(ldReference, nrow(pairs)))
    pairs$ldReference <- dict$keys
    groups <- .crGroups(pairs, tables)
    grl <- S4Vectors::`mcols<-`(
        .rtlBindElements(map(tables[groups$first], .crGroupVariants)),
        value = .crGroupMcols(pairs, tables, groups)
    )
    obj <- new("ColocResult", .rtlNamed(grl), ldReference = dict$entries)
    .withMetadata(obj, metadata)
}

# Each pair's variant table as a data frame with a variant_id and a
# SNP.PP.H4 column, all with the same columns (the per-variant columns a
# group's pairs carry side by side).
# @noRd
.crVariantTables <- function(variants) {
    tables <- map(variants, as.data.frame, stringsAsFactors = FALSE)
    required <- colnames(.crEmptyVariants())
    missingCols <- setdiff(
        required,
        reduce(map(tables, colnames), intersect, .init = required)
    )
    if (length(missingCols) > 0L) {
        msg <- glue(
            "each `variants` element needs columns ",
            "{str_flatten(missingCols, ', ')}."
        )
        abort(msg)
    }
    if (n_distinct(map_chr(tables, .crColumnKey)) > 1L) {
        abort("every `variants` element must carry the same columns.")
    }
    tables
}

# @noRd
.crColumnKey <- function(v) {
    str_flatten(colnames(v), "\r")
}

# Which scored group each pair is in: pairs that share their identity
# columns, block, LD reference and variants, in order of first appearance.
# `first` is each group's first pair, `order` the pairs group by group.
# @noRd
.crGroups <- function(pairs, tables) {
    keys <- str_c(
        .rtlTupleKeys(pairs, c(.crGroupCols(), "ldReference"), nrow(pairs)),
        map_chr(tables, .crVariantSetKey),
        sep = "\r"
    )
    groupKeys <- unique(keys)
    of <- match(keys, groupKeys)
    list(of = of, first = match(groupKeys, keys), order = order(of))
}

# @noRd
.crVariantSetKey <- function(v) {
    rlang::hash(as.character(v$variant_id))
}

# One group's variants as ranges with their alleles, parsed from the ids:
# a variant id is a RENDERING of (range, alleles) and never the stored
# identity, which is what lets subsetByOverlaps() reach coloc results at all.
# @noRd
.crGroupVariants <- function(v) {
    .variantIdsToGRanges(as.character(v$variant_id), what = "variant_id")
}

# The groups' outer columns: the identity, block and LD reference key their
# pairs share, then `pairs`.
# @noRd
.crGroupMcols <- function(pairs, tables, groups) {
    shared <- intersect(c(.crGroupCols(), "ldReference"), colnames(pairs))
    md <- exec(
        S4Vectors::DataFrame,
        !!!c(
            as.list(pairs[groups$first, shared, drop = FALSE]),
            list(check.names = FALSE)
        )
    )
    # Assigned, not passed to DataFrame(), which would spread the tables.
    md$pairs <- .crPairsColumn(pairs, tables, groups)
    md
}

# Each group's pairs as one table: every other column of `pairs`, and each
# pair's per-variant columns -- SNP.PP.H4 and any other -- as list columns
# aligned to the group's variants.
# @noRd
.crPairsColumn <- function(pairs, tables, groups) {
    ord <- groups$order
    own <- setdiff(colnames(pairs), c(.crGroupCols(), "ldReference"))
    perVariant <- setdiff(
        colnames(if (length(tables) > 0L) tables[[1L]] else .crEmptyVariants()),
        "variant_id"
    )
    flat <- exec(
        S4Vectors::DataFrame,
        !!!c(
            as.list(pairs[ord, own, drop = FALSE]),
            set_names(
                map(perVariant, .crPerVariantList, tables = tables[ord]),
                perVariant
            ),
            list(check.names = FALSE)
        )
    )
    S4Vectors::splitAsList(
        flat,
        factor(groups$of[ord], levels = seq_along(groups$first))
    )
}

# Column `cn` of each table, as one list column entry per pair.
# @noRd
.crPerVariantList <- function(cn, tables) {
    values <- unname(map(tables, cn))
    flat <- if (length(values) == 0L) numeric(0) else list_c(values)
    IRanges::relist(
        flat,
        IRanges::PartitioningByEnd(cumsum(lengths(values)))
    )
}

# Pivot coloc.bf_bf's WIDE $results (one SNP.PP.H4.rowK column per $summary
# row) into one variant table per pair.
#
# Verified against coloc: the K-th results column corresponds positionally to
# the K-th summary row. This pivot happens once, here, because this is the only
# place that knows the row-to-pair correspondence -- downstream views must
# never have to re-derive it.
# @noRd
.crPivotColocResults <- function(results, nPairs) {
    if (is.null(results) || nrow(results) == 0L) {
        return(rep(list(.crEmptyVariants()), nPairs))
    }
    results <- as.data.frame(results, stringsAsFactors = FALSE)
    map(seq_len(nPairs), .crPivotOne, results = results)
}

# @noRd
.crPivotOne <- function(k, results) {
    col <- .crPairPpColumn(k, colnames(results))
    if (is.null(col)) {
        return(.crEmptyVariants())
    }
    tibble(
        variant_id = as.character(results$snp),
        SNP.PP.H4 = as.numeric(results[[col]])
    )
}

# The column of coloc's per-variant results holding pair `k`'s SNP.PP.H4, or
# NULL when it holds none. coloc names it SNP.PP.H4.row<k> when it scored
# several pairs and SNP.PP.H4.abf when it scored one (coloc.bf_bf).
# @noRd
.crPairPpColumn <- function(k, cols) {
    found <- intersect(
        c(str_c("SNP.PP.H4.row", k), if (k == 1L) "SNP.PP.H4.abf"),
        cols
    )
    if (length(found) == 0L) NULL else found[[1L]]
}

# @noRd
.crEmptyVariants <- function() {
    tibble(variant_id = character(0), SNP.PP.H4 = numeric(0))
}

# ---- element operations and saved objects ----------------------------------

# `rowsMd`, the outer columns of `x`'s rows `keep`, with each pair's
# per-variant columns following its group's variants into `bare`, the
# rebuilt elements, by variant id: an element operation -- endoapply(), sort(),
# a filter() -- may drop or reorder a group's variants.
# @noRd
.crRealignedMcols <- function(x, rowsMd, bare, keep) {
    old <- unname(.crElementIds(x)[keep])
    new <- unname(.crElementIds(bare))
    if (identical(old, new)) {
        return(rowsMd)
    }
    rowsMd$pairs <- .crRealignPairs(rowsMd$pairs, map2(new, old, match))
    rowsMd
}

# `pairs` with each pair's per-variant columns taken at its group's positions
# `pos` of the rebuilt variants among the old.
# @noRd
.crRealignPairs <- function(pairs, pos) {
    flat <- unlist(pairs, use.names = FALSE)
    at <- pos[rep(seq_along(pos), elementNROWS(pairs))]
    listCols <- names(keep(as.list(flat), methods::is, "List"))
    IRanges::relist(
        reduce(listCols, .crRealignColumn, .init = flat, at = at),
        pairs
    )
}

# @noRd
.crRealignColumn <- function(flat, cn, at) {
    values <- map2(as.list(flat[[cn]]), at, `[`)
    flat[[cn]] <- IRanges::relist(
        list_c(values, ptype = numeric()),
        IRanges::PartitioningByEnd(cumsum(lengths(values)))
    )
    flat
}

# A ColocResult saved with one element per credible-set pair -- each holding
# the pair's variants and SNP.PP.H4 -- folded into one element per scored
# group (ColocResult()). A current one is returned as it is.
# @noRd
.crFoldSavedPairs <- function(object) {
    object <- .renameSavedSlots(object, c(ldSketch = "ldReference"))
    md <- attr(object, "elementMetadata")
    if (is_in("pairs", colnames(md))) {
        return(object)
    }
    elements <- as.list(as(object, "CompressedGRangesList"))
    pairs <- as.data.frame(md, optional = TRUE)
    pairs$ldReference <- NULL
    folded <- ColocResult(
        pairs,
        map(elements, .crSavedVariants),
        ldReference = .crSavedPairRefs(object, md),
        metadata = S4Vectors::metadata(object)
    )
    `attr<-`(folded, "gwasFineMapping", attr(object, "gwasFineMapping"))
}

# A saved pair's variant table: its ids, rendered from ranges and alleles,
# and its per-variant columns.
# @noRd
.crSavedVariants <- function(el) {
    mc <- as.data.frame(mcols(el, use.names = FALSE), optional = TRUE)
    as_tibble(
        c(
            list(variant_id = .grVariantIds(el)),
            as.list(mc[setdiff(colnames(mc), c("A1", "A2"))])
        ),
        .name_repair = "minimal"
    )
}

# The LD reference each saved pair used: its key's entry, or the one
# reference an object saved before references were per row held.
# @noRd
.crSavedPairRefs <- function(object, md) {
    old <- methods::slot(object, "ldReference")
    if (!is.list(old) || is.object(old)) {
        return(rep(list(old), nrow(md)))
    }
    keys <- md$ldReference %||% rep(NA_character_, nrow(md))
    map(as.character(keys), .dictEntryFor, entries = old)
}

#' @rdname updateObject-methods
#' @export
setMethod(
    "updateObject",
    "ColocResult",
    function(object, ..., verbose = FALSE) {
        methods::callNextMethod(
            .crFoldSavedPairs(object),
            ...,
            verbose = verbose
        )
    }
)

# ---- views ------------------------------------------------------------------

# The columns that identify a gene-level unit: everything that is fixed within
# one first-side trait tested against one second-side trait. Both sides carry
# the full (study, context, trait, method) tuple, since either may be a QTL
# collection whose rows differ only on context or trait -- grouping on the
# study alone would pool two distinct molecular phenotypes into one unit.
# @noRd
.crGeneCols <- function() {
    c(
        "study",
        "context",
        "trait",
        "method",
        "gwasStudy",
        "gwasContext",
        "gwasTrait",
        "gwasMethod"
    )
}

#' @rdname colocViews
#' @export
setMethod(
    "colocPairs",
    "ColocResult",
    function(x) {
        .asViewDataFrame(.crPairsTable(x = x))
    }
)

# colocPairs() as a tibble, for the views built on it: one row per pair, its
# group's identity and block, its own columns, then its group's LD reference
# key -- the per-variant columns left in `pairs`.
# @noRd
.crPairsTable <- function(x) {
    md <- mcols(x, use.names = FALSE)
    if (is.null(md)) {
        return(tibble())
    }
    group <- .crPairGroups(x)
    flat <- unlist(md$pairs, use.names = FALSE)
    own <- flat[, !map_lgl(as.list(flat), methods::is, "List"), drop = FALSE]
    shared <- intersect(.crGroupCols(), colnames(md))
    bind_cols(
        as_tibble(
            as.data.frame(md[group, shared, drop = FALSE], optional = TRUE),
            .name_repair = "minimal"
        ),
        as_tibble(
            as.data.frame(own, optional = TRUE),
            .name_repair = "minimal"
        ),
        tibble(ldReference = as.character(md$ldReference)[group]),
        .name_repair = "minimal"
    )
}

# The group each pair is in, pair by pair.
# @noRd
#' @importFrom S4Vectors elementNROWS
.crPairGroups <- function(x) {
    rep(seq_along(x), elementNROWS(mcols(x, use.names = FALSE)$pairs))
}

# Each pair's SNP.PP.H4, pair by pair, aligned to its group's variants.
# @noRd
.crPairPp <- function(x) {
    unlist(mcols(x, use.names = FALSE)$pairs, use.names = FALSE)$SNP.PP.H4
}

# Each element's variant ids, element by element.
# @noRd
.crElementIds <- function(grl) {
    split(
        .grVariantIds(unlist(grl, use.names = FALSE)),
        factor(rep(seq_along(grl), lengths(grl)), levels = seq_along(grl))
    )
}

#' @rdname colocViews
#' @export
setMethod("colocVariants", "ColocResult", function(x) {
    # Ranged, placed by each row's variant id, so the view narrows with
    # subsetByOverlaps().
    .withGenome(.variantTableGRanges(.crLongVariants(x)), .oneBuild(x))
})

#' @rdname colocViews
#' @export
setMethod("colocPooledVariants", "ColocResult", function(x) {
    long <- .crLongVariants(x)
    pooled <- if (nrow(long) == 0L) long else .crPoolVariants(long)
    .withGenome(.variantTableGRanges(pooled), .oneBuild(x))
})

# One row per (pair, variant): the pair's columns repeated across its group's
# variants, plus colocPp -- the posterior that THIS variant is the shared
# causal one, PP.H4 for the pair times the variant's share within it.
# @noRd
.crLongVariants <- function(x) {
    if (length(x) == 0L) {
        return(tibble())
    }
    pairs <- .crPairsTable(x)
    group <- .crPairGroups(x)
    snpPp <- as.numeric(unlist(.crPairPp(x), use.names = FALSE))
    pairs[rep(seq_len(nrow(pairs)), lengths(x)[group]), , drop = FALSE] |>
        mutate(
            variant_id = list_c(
                .crElementIds(x)[group],
                ptype = character()
            ),
            SNP.PP.H4 = snpPp,
            colocPp = .data$PP.H4.abf * snpPp
        )
}

# Pool per-variant posteriors to the gene tier, by the section 3.5 rule: sum
# within a QTL credible set (its per-block / per-GWAS-CS results are mutually
# exclusive), then noisy-OR across QTL credible sets (independent signals).
# @noRd
.crPoolVariants <- function(long) {
    byCs <- long |>
        group_by(across(all_of(c(.crGeneCols(), "qtlCs", "variant_id")))) |>
        summarise(csPp = sum(.data$colocPp), .groups = "drop")
    .crWarnOverflow(byCs$csPp, "variant")
    byCs |>
        mutate(csPp = pmin(1, .data$csPp)) |>
        group_by(across(all_of(c(.crGeneCols(), "variant_id")))) |>
        summarise(colocPp = .crNoisyOr(.data$csPp), .groups = "drop")
}

# 1 - prod(1 - p): the probability at least one independent signal colocalizes.
# @noRd
.crNoisyOr <- function(p) {
    1 - prod(1 - p)
}

# A sum above 1 is not clipped silently: it means several GWAS credible sets
# are competing for one QTL signal, which is diagnostic rather than noise.
# @noRd
.crWarnOverflow <- function(sums, tier) {
    over <- sum(sums > 1 + 1e-8, na.rm = TRUE)
    if (over == 0L) {
        return(invisible(NULL))
    }
    msg <- glue(
        "getColoc*(): {over} {tier}-level sum(s) exceeded 1 before clipping; ",
        "several GWAS credible sets are competing for one QTL signal. The ",
        "pair-level view shows which."
    )
    warn(msg)
    invisible(NULL)
}

#' @rdname colocViews
#' @export
setMethod(
    "colocGenes",
    "ColocResult",
    function(x) {
        .asViewDataFrame(.crGenesTable(x = x))
    }
)

# colocGenes() as a tibble.
# @noRd
.crGenesTable <- function(x) {
    pairs <- .crPairsTable(x)
    if (nrow(pairs) == 0L) {
        return(tibble())
    }
    byCs <- pairs |>
        group_by(across(all_of(c(.crGeneCols(), "qtlCs")))) |>
        summarise(
            csPp = sum(.data$PP.H4.abf),
            nPairs = dplyr::n(),
            .groups = "drop"
        )
    .crWarnOverflow(byCs$csPp, "gene")
    byCs |>
        mutate(csPp = pmin(1, .data$csPp)) |>
        group_by(across(all_of(.crGeneCols()))) |>
        summarise(
            PP.H4 = .crNoisyOr(.data$csPp),
            nQtlCs = dplyr::n(),
            nPairs = sum(.data$nPairs),
            .groups = "drop"
        )
}

# ---- credible-set view ------------------------------------------------------

#' @rdname colocViews
#' @export
setMethod(
    "colocCredibleSets",
    "ColocResult",
    function(
        x,
        coverage = 0.95,
        minPp4 = NULL,
        minAbsCorr = 0.8,
        requireMaxH4 = FALSE
    ) {
        sets <- .crCredibleSets(
            x = x,
            coverage = coverage,
            minPp4 = minPp4,
            minAbsCorr = minAbsCorr,
            requireMaxH4 = requireMaxH4
        )
        .withGenome(.crSetsAsList(sets), .oneBuild(x))
    }
)

# The coloc credible sets: one record per set, its members' ranges and its
# set-level row, filtered on the pair and then on purity.
# @noRd
.crCredibleSets <- function(
    x,
    coverage = 0.95,
    minPp4 = NULL,
    minAbsCorr = 0.8,
    requireMaxH4 = FALSE
) {
    view <- list(
        x = x,
        pairs = .crPairsTable(x),
        group = .crPairGroups(x),
        pp = as.list(.crPairPp(x))
    )
    kept <- .crPairFilter(view$pairs, minPp4, requireMaxH4)
    # Filtering happens BEFORE any LD work, so a stricter threshold costs
    # strictly less -- which is the reason purity lives on the accessor
    # rather than being precomputed at construction.
    sets <- compact(map(kept, .crCsForPair, view = view, coverage = coverage))
    if (length(sets) == 0L) {
        return(list())
    }
    purity <- map_dbl(sets, .crSetPurity, refs = ldReference(x))
    sets <- map2(sets, purity, .crWithPurity)
    if (is.null(minAbsCorr)) {
        return(sets)
    }
    # An unmeasurable purity does not drop a set: with no LD reference there is
    # no evidence the set is impure, and silently discarding it would make the
    # view depend on whether a sketch happened to be attached.
    keep(sets, .crPureEnough, minAbsCorr = minAbsCorr)
}

# The sets as one GRangesList: each set's members its element, its set-level
# columns that element's mcols.
# @noRd
.crSetsAsList <- function(sets) {
    if (length(sets) == 0L) {
        return(GenomicRanges::GRangesList())
    }
    S4Vectors::`mcols<-`(
        .rtlBindElements(map(sets, "members")),
        value = S4Vectors::DataFrame(
            bind_rows(map(sets, "row")),
            check.names = FALSE
        )
    )
}

# Purity of one set: the minimum absolute correlation among its members, one
# LD extraction per set (section 3.7), against its own pair's reference. NA
# when the pair has no LD reference to extract from -- reported honestly
# rather than defaulted to a passing value.
# @noRd
.crSetPurity <- function(set, refs) {
    .crPurityOne(.grVariantIds(set$members), set$row$ldReference, refs)
}

# @noRd
.crWithPurity <- function(set, purity) {
    list_assign(set, row = mutate(set$row, purity = purity))
}

# @noRd
.crPureEnough <- function(set, minAbsCorr) {
    p <- set$row$purity
    is.na(p) || p >= minAbsCorr
}

# Which pairs -- rows of the colocPairs() table `md` -- survive the
# pair-level filters. Both conditions read columns already stored on the
# pair, so neither costs anything to evaluate.
# @noRd
.crPairFilter <- function(md, minPp4, requireMaxH4) {
    if (nrow(md) == 0L) {
        return(integer(0))
    }
    byPp4 <- if (is.null(minPp4)) {
        rep(TRUE, nrow(md))
    } else {
        md$PP.H4.abf >= minPp4
    }
    byMaxH4 <- if (isTRUE(requireMaxH4)) {
        pp <- as.matrix(md[, .crPpCols(), drop = FALSE])
        max.col(pp, ties.method = "first") == 5L
    } else {
        TRUE
    }
    ok <- byPp4 & byMaxH4
    which(ok & !is.na(ok))
}

# The coloc credible set of one pair: the smallest prefix of variants, ordered
# by SNP.PP.H4, whose cumulative mass reaches `coverage`.
#
# This is a THIRD variant set -- section 3.7 -- not guaranteed to be a subset of
# either input credible set, which is why nothing about it can be inherited
# from the QTL or GWAS fine-mapping.
#
# `view` holds the result `x`, its colocPairs() table `pairs`, each pair's
# group and each pair's SNP.PP.H4 (`pp`).
# @noRd
.crCsForPair <- function(i, view, coverage) {
    g <- view$x[[view$group[[i]]]]
    if (length(g) == 0L) {
        return(NULL)
    }
    pp <- as.numeric(view$pp[[i]])
    ord <- order(pp, decreasing = TRUE)
    cum <- cumsum(pp[ord])
    nKeep <- .crFirstAtLeast(cum, coverage)
    members <- ord[seq_len(nKeep)]
    ids <- .grVariantIds(g[members])
    list(
        members = S4Vectors::`mcols<-`(
            g[members],
            value = cbind(
                mcols(g[members], use.names = FALSE),
                S4Vectors::DataFrame(SNP.PP.H4 = pp[members])
            )
        ),
        row = mutate(
            view$pairs[i, , drop = FALSE],
            csSize = length(members),
            csCoverage = cum[[nKeep]],
            leadVariant = ids[[1L]],
            leadPp = pp[ord][[1L]]
        )
    )
}

# Index of the first cumulative value reaching `target`; the whole vector when
# it never does (the set cannot cover more than everything).
#
# The tolerance absorbs accumulated floating-point error, not statistical
# slack: 0.6 + 0.3 is 0.8999... in binary, so a bare `>=` would take a third
# variant into a 0.9-coverage set whose first two already carry exactly 0.9.
# @noRd
.crFirstAtLeast <- function(cum, target, tol = 1e-10) {
    hit <- which(cum >= target - tol)
    if (length(hit) == 0L) length(cum) else hit[[1L]]
}

# @noRd
#' @importFrom rlang try_fetch
.crPurityOne <- function(ids, key, refs) {
    ldReference <- .dictEntryFor(key, refs)
    if (is.null(ldReference)) {
        return(NA_real_)
    }
    # A singleton set has no pair to correlate, and susie treats it as pure.
    if (length(ids) < 2L) {
        return(1)
    }
    ld <- try_fetch(
        .ldFromRef(
            ldReference,
            ids,
            label = "colocCredibleSets",
            onMissing = "drop"
        ),
        error = function(cnd) NULL
    )
    if (is.null(ld) || nrow(ld) < 2L) {
        return(NA_real_)
    }
    min(abs(ld[upper.tri(ld)]))
}
