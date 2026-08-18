# Multiply a signed-effect column by the per-row sign flip (lambda-free across
# helper for .overlapRelabelGwas).
# @noRd
.overlapApplySign <- function(x, sign) x * sign

# Prefix a top-loci frame's NON-key columns with `pfx` (key columns unchanged).
# @noRd
.overlapPrefixNonKey <- function(df, pfx, keyCols) {
    set_names(
        df,
        if_else(
            is_in(names(df), keyCols),
            names(df),
            str_c(pfx, names(df))
        )
    )
}

# The zero-row merged frame returned when either side has no signal: an
# identity-only join on variant_id (GWAS contributes only its non-coord cols).
# @noRd
.overlapEmptyMerge <- function(qtlTl, gwasTl, coordCols, keyCols) {
    qp <- .overlapPrefixNonKey(slice(qtlTl, 0), "qtl_", keyCols)
    gp <- .overlapPrefixNonKey(
        select(slice(gwasTl, 0), all_of(setdiff(names(gwasTl), coordCols))),
        "gwas_",
        keyCols
    )
    # A collection with no signal above the cutoff yields an identity-only
    # getTopLoci frame (study/context/trait/method) that carries no
    # variant_id; restore the join key so the empty join still resolves
    # instead of erroring in inner_join()'s `by` check.
    if (!is_in("variant_id", names(qp))) {
        qp$variant_id <- character(0)
    }
    if (!is_in("variant_id", names(gp))) {
        gp$variant_id <- character(0)
    }
    inner_join(qp, gp, by = "variant_id")
}

# Convert the merged frame to the requested output type.
# @noRd
.overlapFinish <- function(df, type) {
    if (type == "GRanges") .overlapToGRanges(df) else df
}

#' Overlap QTL and GWAS top loci by allele-aware variant matching
#'
#' Intersect the top-loci tables of a \code{QtlFineMappingResult} and a
#' \code{GwasFineMappingResult} on shared variants, matched with pecotmr's
#' allele-aware \code{matchVariants} (handling strand flips / ref-alt swaps
#' rather than naive id equality). The GWAS side is harmonized to the QTL
#' orientation: its signed effect columns (\code{beta}, \code{z},
#' \code{conditional_effect}) are sign-flipped and its effect-allele frequency
#' (\code{af}) is complemented wherever a swap occurred. The result keeps the
#' variant key columns once (from the QTL, the reference orientation) and
#' prefixes every other column \code{qtl_} / \code{gwas_}. A variant shared
#' across several QTL contexts and/or GWAS studies yields one row per (QTL entry
#' x GWAS entry) pair (a wide cross-product per variant).
#'
#' @param qtl A \code{QtlFineMappingResult}.
#' @param gwas A \code{GwasFineMappingResult}.
#' @param signalCutoff PIP cutoff forwarded to \code{\link{getTopLoci}} for both
#'   inputs. Default 0.025.
#' @param type \code{"data.frame"} (default) or \code{"GRanges"}.
#' @param ... Ignored.
#' @return A \code{tibble} (or \code{GRanges}) keyed on the QTL variant
#'   (\code{variant_id, chrom, pos, A1, A2}) with all other columns prefixed
#'   \code{qtl_} / \code{gwas_}. Zero rows when there is no allele-aware
#'   overlap.
#' @seealso \code{\link{getTopLoci}}, \code{matchVariants}
#' @examples
#' data(qtlFineMappingLbfExample)
#' data(gwasFineMappingLbfExample)
#' overlapTopLoci(qtlFineMappingLbfExample, gwasFineMappingLbfExample)
#' @include AllGenerics.R AllClasses.R QtlFineMappingResult.R
#'   GwasFineMappingResult.R
#' @export
setGeneric("overlapTopLoci", function(qtl, gwas, ...) {
    standardGeneric("overlapTopLoci")
})

#' @rdname overlapTopLoci
#' @export
setMethod(
    "overlapTopLoci",
    signature("QtlFineMappingResult", "GwasFineMappingResult"),
    function(
        qtl,
        gwas,
        signalCutoff = 0.025,
        type = c("data.frame", "GRanges"),
        ...
    ) {
        type <- arg_match(type)
        keyCols <- c("variant_id", "chrom", "pos", "A1", "A2")
        coordCols <- c("chrom", "pos", "A1", "A2")
        qtlTl <- as_tibble(getTopLoci(qtl, signalCutoff = signalCutoff))
        gwasTl <- as_tibble(getTopLoci(gwas, signalCutoff = signalCutoff))
        if (nrow(qtlTl) == 0L || nrow(gwasTl) == 0L) {
            return(.overlapEmptyReturn(qtlTl, gwasTl, coordCols, keyCols, type))
        }
        # Allele-aware correspondence between the unique QTL + GWAS variant
        # sets. target = GWAS, ref = QTL, so `sign` is the flip applied to the
        # GWAS side.
        vmap <- .overlapVariantMap(qtlTl, gwasTl)
        if (is.null(vmap)) {
            return(.overlapEmptyReturn(qtlTl, gwasTl, coordCols, keyCols, type))
        }
        g <- .overlapRelabelGwas(gwasTl, vmap, coordCols)
        if (nrow(g) == 0L) {
            # nocov start
            return(.overlapEmptyReturn(qtlTl, gwasTl, coordCols, keyCols, type))
            # nocov end
        }
        # Wide cross-product per shared variant, variant key kept once (QTL).
        merged <- inner_join(
            .overlapPrefixNonKey(qtlTl, "qtl_", keyCols),
            .overlapPrefixNonKey(g, "gwas_", keyCols),
            by = "variant_id"
        )
        # Restore key-column order (put the variant key columns first).
        ordered <- c(
            intersect(keyCols, names(merged)),
            setdiff(names(merged), keyCols)
        )
        .overlapFinish(select(merged, all_of(ordered)), type)
    }
)

# Finish an empty-overlap result (an empty prefixed merge) in the chosen type.
# @noRd
.overlapEmptyReturn <- function(qtlTl, gwasTl, coordCols, keyCols, type) {
    .overlapFinish(
        .overlapEmptyMerge(qtlTl, gwasTl, coordCols, keyCols),
        type
    )
}

# Allele-aware unique-variant correspondence (GWAS -> canonical QTL id + sign),
# or NULL when nothing matches.
# @noRd
.overlapVariantMap <- function(qtlTl, gwasTl) {
    uq <- unique(qtlTl$variant_id)
    ug <- unique(gwasTl$variant_id)
    m <- matchVariants(ug, uq, allowFlip = TRUE)
    if (length(m$idxA) == 0L) {
        return(NULL)
    }
    tibble(
        gwas_vid = ug[m$idxA],
        canon_vid = uq[m$idxB],
        .sign = as.numeric(m$sign)
    )
}

# Relabel every GWAS row to the canonical (QTL-orientation) variant, applying
# the swap to signed effects + effect-allele frequency, and drop the coord /
# helper columns.
# @noRd
.overlapRelabelGwas <- function(gwasTl, vmap, coordCols) {
    g <- inner_join(gwasTl, vmap, by = c("variant_id" = "gwas_vid"))
    signedCols <- c("beta", "z", "conditional_effect")
    g <- mutate(g, across(any_of(signedCols), .overlapApplySign, g$.sign))
    if (is_in("af", names(g))) {
        g$af <- if_else(g$.sign < 0 & !is.na(g$af), 1 - g$af, g$af)
    }
    g$variant_id <- g$canon_vid
    select(g, all_of(setdiff(names(g), c(coordCols, "canon_vid", ".sign"))))
}

# Build a GRanges from an overlap table: variants as width-1 ranges, all other
# columns (qtl_* / gwas_*) as mcols.
# @noRd
.overlapToGRanges <- function(df) {
    if (is.null(df) || nrow(df) == 0L) {
        return(GenomicRanges::GRanges())
    }
    p <- parseVariantId(df$variant_id)
    gr <- GenomicRanges::GRanges(
        seqnames = str_c("chr", p$chrom),
        ranges = IRanges::IRanges(start = p$pos, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(df, check.names = FALSE)
    gr
}
