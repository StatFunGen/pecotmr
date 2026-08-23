# =============================================================================
# ColocResult S4 class
# -----------------------------------------------------------------------------
# The output of colocPipeline / colocboostPipeline. One ELEMENT per tested
# (QTL credible set x GWAS credible set x block) pair, holding that pair's
# aligned variants; the pair-level statistics live in outer mcols and the
# per-variant ones in inner mcols.
#
# Why a class and not the flat tibble it replaces: the pair / variant /
# credible-set / gene "tiers" are not four shapes of table, they are four VIEWS
# over one representation. Storing any single granularity forces the coarser
# statistics to be repeated down the finer rows, and makes the finer ones
# unrecoverable -- which is exactly how the per-variant layer got lost when
# process_coloc_results() was dropped.
#
#   element      the aligned variants of one pair
#   inner mcols  per-variant: SNP.PP.H4
#   outer mcols  per-pair: identity tuple, blockId, PP.H0..PP.H4, nSnps,
#                retained alpha mass (both sides), enrichment / p12Used
#   slots        ldSketch, for the credible-set view's purity calculation
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
#'   projection surface: both answer \code{\link{getColocPairs}} and
#'   \code{\link{getColocVariants}}, so code that only needs "which variants
#'   colocalized, and how strongly" works against either.
#' @export
setClass("ColocResultBase", contains = c("VIRTUAL", "RangedTupleList"))

#' @title Colocalization Result
#' @description A collection of colocalization results, one element per tested
#'   (QTL credible set, GWAS credible set, block) pair. Build the four views
#'   with \code{\link{getColocPairs}}, \code{\link{getColocVariants}},
#'   \code{\link{getColocCredibleSets}} and \code{\link{getColocGenes}}.
#' @slot ldSketch The LD reference \code{GenotypeHandle} the underlying fits
#'   were computed against, or \code{NULL}. Used by
#'   \code{getColocCredibleSets()} to compute credible-set purity.
#' @seealso \code{\link{colocPipeline}}
#' @export
setClass(
    "ColocResult",
    contains = "ColocResultBase",
    representation(ldSketch = "ANY"),
    prototype(ldSketch = NULL)
)

methods::setValidity("ColocResult", function(object) {
    .validateColocResult(object)
})

# @noRd
.validateColocResult <- function(object) {
    errors <- .crCheckRequiredCols(object)
    if (length(errors) == 0L) {
        errors <- c(
            .crCheckPpColumns(object),
            .crCheckVariantColumn(object)
        )
    }
    if (length(errors) == 0L) TRUE else errors
}

# @noRd
.crRequiredCols <- function() {
    c(
        "study",
        "context",
        "trait",
        "method",
        "gwasStudy",
        "gwasMethod",
        "blockId",
        "qtlCs",
        "gwasCs",
        "nSnps"
    )
}

# @noRd
.crPpCols <- function() {
    str_c("PP.H", 0:4, ".abf")
}

# @noRd
.crCheckRequiredCols <- function(object) {
    md <- mcols(object, use.names = FALSE)
    have <- if (is.null(md)) character(0) else colnames(md)
    missingCols <- setdiff(.crRequiredCols(), have)
    if (length(missingCols) > 0L) {
        return(str_c("missing columns: ", str_flatten(missingCols, ", ")))
    }
    NULL
}

# @noRd
.crCheckPpColumns <- function(object) {
    md <- mcols(object, use.names = FALSE)
    missingCols <- setdiff(.crPpCols(), colnames(md))
    if (length(missingCols) > 0L) {
        return(str_c(
            "missing posterior columns: ",
            str_flatten(missingCols, ", ")
        ))
    }
    NULL
}

# The per-variant layer is the whole point of the class, so an element without
# it is rejected rather than tolerated -- a silently absent SNP.PP.H4 would
# reproduce the exact gap this class exists to close.
# @noRd
.crCheckVariantColumn <- function(object) {
    if (length(object) == 0L) {
        return(NULL)
    }
    bad <- which(!map_lgl(as.list(object), .crHasSnpPp))
    if (length(bad) == 0L) {
        return(NULL)
    }
    str_c(
        "element(s) ",
        str_flatten(bad, ", "),
        " have no SNP.PP.H4 column in their variant metadata"
    )
}

# @noRd
.crHasSnpPp <- function(g) {
    md <- mcols(g, use.names = FALSE)
    !is.null(md) && is_in("SNP.PP.H4", colnames(md))
}

# ---- accessors --------------------------------------------------------------

#' @rdname getLdSketch
#' @export
setMethod("getLdSketch", "ColocResult", function(x, ...) x@ldSketch)

#' @rdname show-methods
#' @export
setMethod("show", "ColocResult", function(object) {
    cat(class(object), "with", length(object), "colocalized pair(s)\n")
    if (length(object) == 0L) {
        return(invisible(NULL))
    }
    md <- mcols(object, use.names = FALSE)
    cat("  QTL studies :", str_flatten(unique(md$study), ", "), "\n")
    cat("  GWAS studies:", str_flatten(unique(md$gwasStudy), ", "), "\n")
    cat("  variants    :", sum(lengths(object)), "across all pairs\n")
    cat(
        "  max PP.H4   :",
        format(max(md$PP.H4.abf, na.rm = TRUE), digits = 4),
        "\n"
    )
    invisible(NULL)
})

# ---- construction -----------------------------------------------------------

#' @title Build a ColocResult
#' @description Assemble a \code{\linkS4class{ColocResult}} from a pair-level
#'   table and the per-pair variant tables that go with it. Callers normally
#'   get one from \code{\link{colocPipeline}} rather than building it directly.
#' @param pairs A data frame with one row per tested pair, carrying at least
#'   the identity columns (\code{study}, \code{context}, \code{trait},
#'   \code{method}, \code{gwasStudy}, \code{gwasMethod}), \code{blockId},
#'   \code{qtlCs}, \code{gwasCs}, \code{nSnps} and \code{PP.H0.abf} through
#'   \code{PP.H4.abf}.
#' @param variants A list, parallel to \code{pairs}' rows, of per-pair data
#'   frames with a \code{variant_id} column and a \code{SNP.PP.H4} column.
#' @param ldSketch Optional \code{GenotypeHandle} for the LD reference, used by
#'   \code{\link{getColocCredibleSets}} to compute purity.
#' @return A \code{ColocResult}.
#' @examples
#' pairs <- data.frame(
#'     study = "s1", context = "c1", trait = "g1", method = "susie",
#'     gwasStudy = "G1", gwasMethod = "susie", blockId = "chr1_1_1000",
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
ColocResult <- function(pairs, variants, ldSketch = NULL) {
    pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
    if (!is.list(variants) || length(variants) != nrow(pairs)) {
        msg <- glue(
            "`variants` must be a list running parallel to `pairs`' rows ",
            "(got {length(variants)} vs {nrow(pairs)})."
        )
        abort(msg)
    }
    elements <- map(variants, .crVariantsToGRanges)
    grl <- GenomicRanges::GRangesList(elements)
    # Set unconditionally, including at zero rows: an empty result still has to
    # carry the column schema, or it fails validity and a caller reading
    # `result$PP.H4.abf` gets NULL exactly when there is nothing to report.
    mcols(grl) <- exec(
        S4Vectors::DataFrame,
        !!!c(as.list(pairs), list(check.names = FALSE))
    )
    obj <- new("ColocResult", grl, ldSketch = ldSketch)
    validObject(obj)
    obj
}

# One pair's variant table -> a GRanges carrying SNP.PP.H4 in its mcols.
#
# The ids are parsed back into ranges rather than stored as strings, per the
# rule that a variant id is a RENDERING of (range, alleles) and never the
# stored identity -- that is what lets subsetByOverlaps() reach coloc results
# at all.
# @noRd
.crVariantsToGRanges <- function(v) {
    v <- as.data.frame(v, stringsAsFactors = FALSE)
    missingCols <- setdiff(c("variant_id", "SNP.PP.H4"), colnames(v))
    if (length(missingCols) > 0L) {
        msg <- glue(
            "each `variants` element needs columns ",
            "{str_flatten(missingCols, ', ')}."
        )
        abort(msg)
    }
    gr <- .variantIdsToGRanges(as.character(v$variant_id), what = "variant_id")
    # The A1 / A2 columns .variantIdsToGRanges() attaches ARE the variant
    # identity: range plus alleles is what an id renders, so overwriting mcols
    # wholesale here would leave the element unable to name its own variants.
    # Append instead.
    keep <- setdiff(colnames(v), "variant_id")
    extra <- exec(
        S4Vectors::DataFrame,
        !!!c(as.list(v[keep]), list(check.names = FALSE))
    )
    mcols(gr) <- cbind(mcols(gr, use.names = FALSE), extra)
    gr
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
    col <- str_c("SNP.PP.H4.row", k)
    if (!is_in(col, colnames(results))) {
        return(.crEmptyVariants())
    }
    tibble(
        variant_id = as.character(results$snp),
        SNP.PP.H4 = as.numeric(results[[col]])
    )
}

# @noRd
.crEmptyVariants <- function() {
    tibble(variant_id = character(0), SNP.PP.H4 = numeric(0))
}

# ---- views ------------------------------------------------------------------

# The columns that identify a gene-level unit: everything that is fixed within
# one QTL molecular trait tested against one GWAS study.
# @noRd
.crGeneCols <- function() {
    c("study", "context", "trait", "method", "gwasStudy", "gwasMethod")
}

#' @rdname colocViews
#' @export
setMethod("getColocPairs", "ColocResult", function(x, ...) {
    md <- mcols(x, use.names = FALSE)
    if (is.null(md) || length(x) == 0L) {
        return(tibble())
    }
    as_tibble(as.data.frame(md), .name_repair = "minimal")
})

#' @rdname colocViews
#' @export
setMethod("getColocVariants", "ColocResult", function(x, pooled = FALSE, ...) {
    long <- .crLongVariants(x)
    if (!isTRUE(pooled) || nrow(long) == 0L) {
        return(long)
    }
    .crPoolVariants(long)
})

# One row per (pair, variant): the pair's metadata repeated across its
# variants, plus colocPp -- the posterior that THIS variant is the shared
# causal one, PP.H4 for the pair times the variant's share within it.
# @noRd
.crLongVariants <- function(x) {
    if (length(x) == 0L) {
        return(tibble())
    }
    md <- as.data.frame(mcols(x, use.names = FALSE))
    n <- lengths(x)
    pairIdx <- rep(seq_len(length(x)), n)
    flat <- unlist(x, use.names = FALSE)
    out <- as_tibble(md[pairIdx, , drop = FALSE], .name_repair = "minimal")
    out$variant_id <- .grVariantIds(flat)
    out$SNP.PP.H4 <- as.numeric(mcols(flat, use.names = FALSE)$SNP.PP.H4)
    out$colocPp <- out$PP.H4.abf * out$SNP.PP.H4
    out
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
setMethod("getColocGenes", "ColocResult", function(x, ...) {
    pairs <- getColocPairs(x)
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
})

# ---- credible-set view ------------------------------------------------------

#' @rdname colocViews
#' @export
setMethod(
    "getColocCredibleSets",
    "ColocResult",
    function(
        x,
        coverage = 0.95,
        minPp4 = NULL,
        minAbsCorr = 0.8,
        requireMaxH4 = FALSE,
        ...
    ) {
        keep <- .crPairFilter(x, minPp4, requireMaxH4)
        if (length(keep) == 0L) {
            return(tibble())
        }
        # Filtering happens BEFORE any LD work, so a stricter threshold costs
        # strictly less -- which is the reason purity lives on the accessor
        # rather than being precomputed at construction.
        sets <- map(keep, .crCsForPair, x = x, coverage = coverage)
        sets <- compact(sets)
        if (length(sets) == 0L) {
            return(tibble())
        }
        out <- bind_rows(sets)
        out$purity <- .crPuritiesFor(out, getLdSketch(x))
        .crApplyPurityFilter(out, minAbsCorr)
    }
)

# Which pairs survive the pair-level filters. Both conditions read columns
# already stored on the pair, so neither costs anything to evaluate.
# @noRd
.crPairFilter <- function(x, minPp4, requireMaxH4) {
    if (length(x) == 0L) {
        return(integer(0))
    }
    md <- as.data.frame(mcols(x, use.names = FALSE))
    ok <- rep(TRUE, nrow(md))
    if (!is.null(minPp4)) {
        ok <- ok & md$PP.H4.abf >= minPp4
    }
    if (isTRUE(requireMaxH4)) {
        pp <- as.matrix(md[, .crPpCols(), drop = FALSE])
        ok <- ok & (max.col(pp, ties.method = "first") == 5L)
    }
    which(ok & !is.na(ok))
}

# The coloc credible set of one pair: the smallest prefix of variants, ordered
# by SNP.PP.H4, whose cumulative mass reaches `coverage`.
#
# This is a THIRD variant set -- section 3.7 -- not guaranteed to be a subset of
# either input credible set, which is why nothing about it can be inherited
# from the QTL or GWAS fine-mapping.
# @noRd
.crCsForPair <- function(i, x, coverage) {
    g <- x[[i]]
    if (length(g) == 0L) {
        return(NULL)
    }
    pp <- as.numeric(mcols(g, use.names = FALSE)$SNP.PP.H4)
    ord <- order(pp, decreasing = TRUE)
    cum <- cumsum(pp[ord])
    nKeep <- .crFirstAtLeast(cum, coverage)
    members <- ord[seq_len(nKeep)]
    ids <- .grVariantIds(g[members])
    md <- as.data.frame(mcols(x, use.names = FALSE))[i, , drop = FALSE]
    row <- as_tibble(md, .name_repair = "minimal")
    row$csSize <- length(members)
    row$csCoverage <- cum[[nKeep]]
    row$leadVariant <- ids[[1L]]
    row$leadPp <- pp[ord][[1L]]
    row$csVariants <- list(ids)
    row
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

# Purity per credible set: the minimum absolute correlation among its members,
# one LD extraction per set (section 3.7). NA when there is no LD reference to
# extract from -- reported honestly rather than defaulted to a passing value.
# @noRd
.crPuritiesFor <- function(out, ldSketch) {
    if (is.null(ldSketch)) {
        return(rep(NA_real_, nrow(out)))
    }
    map_dbl(out$csVariants, .crPurityOne, ldSketch = ldSketch)
}

# @noRd
.crPurityOne <- function(ids, ldSketch) {
    # A singleton set has no pair to correlate, and susie treats it as pure.
    if (length(ids) < 2L) {
        return(1)
    }
    ld <- tryCatch(
        .ldFromSketch(
            ldSketch,
            ids,
            label = "getColocCredibleSets",
            onMissing = "drop"
        ),
        error = function(e) NULL
    )
    if (is.null(ld) || nrow(ld) < 2L) {
        return(NA_real_)
    }
    min(abs(ld[upper.tri(ld)]))
}

# @noRd
.crApplyPurityFilter <- function(out, minAbsCorr) {
    if (is.null(minAbsCorr)) {
        return(out)
    }
    # An unmeasurable purity does not drop a set: with no LD reference there is
    # no evidence the set is impure, and silently discarding it would make the
    # view depend on whether a sketch happened to be attached.
    keep <- is.na(out$purity) | out$purity >= minAbsCorr
    out[keep, , drop = FALSE]
}

# ---- coercion ---------------------------------------------------------------

#' @title Coerce a ColocResult to a data frame
#' @description Returns the pair-level view, which is the flat table
#'   \code{colocPipeline()} used to return directly. Provided so existing
#'   consumers keep working; new code should call the view it actually wants
#'   (see \code{\link{colocViews}}).
#' @param x A \code{ColocResult}.
#' @param row.names,optional Ignored; present for generic compatibility.
#' @param ... Ignored.
#' @return A data frame with one row per tested pair.
#' @examples
#' pairs <- data.frame(
#'     study = "s1", context = "c1", trait = "g1", method = "susie",
#'     gwasStudy = "G1", gwasMethod = "susie", blockId = "chr1_1_1000",
#'     qtlCs = 1L, gwasCs = 1L, nSnps = 2L,
#'     PP.H0.abf = 0.1, PP.H1.abf = 0.1, PP.H2.abf = 0.1,
#'     PP.H3.abf = 0.1, PP.H4.abf = 0.6
#' )
#' variants <- list(data.frame(
#'     variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
#'     SNP.PP.H4 = c(0.7, 0.3)
#' ))
#' as.data.frame(ColocResult(pairs, variants))
#' @export
setMethod(
    "as.data.frame",
    "ColocResult",
    function(x, row.names = NULL, optional = FALSE, ...) {
        as.data.frame(getColocPairs(x))
    }
)
