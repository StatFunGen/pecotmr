#' Overlap QTL and GWAS top loci by allele-aware variant matching
#'
#' Intersect the top-loci tables of a \code{QtlFineMappingResult} and a
#' \code{GwasFineMappingResult} on shared variants, matched with pecotmr's
#' allele-aware \code{\link{matchVariants}} (handling strand flips / ref-alt
#' swaps rather than naive id equality). The GWAS side is harmonized to the QTL
#' orientation: its signed effect columns (\code{beta}, \code{z},
#' \code{conditional_effect}) are sign-flipped and its effect-allele frequency
#' (\code{af}) is complemented wherever a swap occurred. The result keeps the
#' variant key columns once (from the QTL, the reference orientation) and
#' prefixes every other column \code{qtl_} / \code{gwas_}. A variant shared
#' across several QTL contexts and/or GWAS studies yields one row per
#' (QTL entry x GWAS entry) pair (a wide cross-product per variant).
#'
#' @param qtl A \code{QtlFineMappingResult}.
#' @param gwas A \code{GwasFineMappingResult}.
#' @param signalCutoff PIP cutoff forwarded to \code{\link{getTopLoci}} for both
#'   inputs. Default 0.025.
#' @param type \code{"data.frame"} (default) or \code{"GRanges"}.
#' @param ... Ignored.
#' @return A \code{data.frame} (or \code{GRanges}) keyed on the QTL variant
#'   (\code{variant_id, chrom, pos, A1, A2}) with all other columns prefixed
#'   \code{qtl_} / \code{gwas_}. Zero rows when there is no allele-aware overlap.
#' @seealso \code{\link{getTopLoci}}, \code{\link{matchVariants}}
#' @include AllGenerics.R AllClasses.R QtlFineMappingResult.R GwasFineMappingResult.R
#' @export
setGeneric("overlapTopLoci",
  function(qtl, gwas, ...) standardGeneric("overlapTopLoci"))

#' @rdname overlapTopLoci
#' @export
setMethod("overlapTopLoci",
  signature("QtlFineMappingResult", "GwasFineMappingResult"),
  function(qtl, gwas, signalCutoff = 0.025,
           type = c("data.frame", "GRanges"), ...) {
    type      <- match.arg(type)
    keyCols   <- c("variant_id", "chrom", "pos", "A1", "A2")
    coordCols <- c("chrom", "pos", "A1", "A2")

    qtlTl  <- as.data.frame(getTopLoci(qtl,  signalCutoff = signalCutoff))
    gwasTl <- as.data.frame(getTopLoci(gwas, signalCutoff = signalCutoff))

    prefixNonKey <- function(df, pfx)
      stats::setNames(df, ifelse(names(df) %in% keyCols, names(df),
                                 paste0(pfx, names(df))))
    # The GWAS side contributes only variant_id (the join key) + its non-key
    # columns; the shared coordinate columns come from the QTL side.
    emptyMerge <- function() {
      qp <- prefixNonKey(qtlTl[0, , drop = FALSE], "qtl_")
      gp <- prefixNonKey(
        gwasTl[0, setdiff(names(gwasTl), coordCols), drop = FALSE], "gwas_")
      # A collection with no signal above the cutoff yields an identity-only
      # getTopLoci frame (study/context/trait/method) that carries no
      # variant_id; restore the join key so the empty merge still resolves
      # instead of erroring in merge()'s `by` check.
      if (!"variant_id" %in% names(qp)) qp$variant_id <- character(0)
      if (!"variant_id" %in% names(gp)) gp$variant_id <- character(0)
      merge(qp, gp, by = "variant_id")
    }
    finish <- function(df) if (type == "GRanges") .overlapToGRanges(df) else df

    if (nrow(qtlTl) == 0L || nrow(gwasTl) == 0L) return(finish(emptyMerge()))

    # Allele-aware correspondence between the unique QTL and GWAS variant sets.
    # target = GWAS, ref = QTL, so `sign` is the flip applied to the GWAS side.
    uq <- unique(qtlTl$variant_id)
    ug <- unique(gwasTl$variant_id)
    m  <- matchVariants(ug, uq, allowFlip = TRUE)
    if (length(m$idxA) == 0L) return(finish(emptyMerge()))
    map <- data.frame(gwas_vid = ug[m$idxA], canon_vid = uq[m$idxB],
                      .sign = as.numeric(m$sign), stringsAsFactors = FALSE)

    # Relabel every GWAS row to the canonical (QTL-orientation) variant, and
    # apply the swap to its signed effects + effect-allele frequency.
    g <- merge(gwasTl, map, by.x = "variant_id", by.y = "gwas_vid")
    # Unreachable defensive guard: map$gwas_vid is drawn from gwasTl$variant_id,
    # so this inner merge can never drop every row (length(m$idxA) > 0 here).
    if (nrow(g) == 0L) return(finish(emptyMerge())) # nocov
    for (cc in intersect(c("beta", "z", "conditional_effect"), names(g)))
      g[[cc]] <- g[[cc]] * g$.sign
    if ("af" %in% names(g))
      g$af <- ifelse(g$.sign < 0 & !is.na(g$af), 1 - g$af, g$af)
    g$variant_id <- g$canon_vid
    g <- g[, setdiff(names(g), c(coordCols, "canon_vid", ".sign")), drop = FALSE]

    # Wide cross-product per shared variant, variant key kept once (QTL side).
    merged <- merge(prefixNonKey(qtlTl, "qtl_"), prefixNonKey(g, "gwas_"),
                    by = "variant_id")
    # Restore key-column order (merge puts the `by` column first).
    ordered <- c(intersect(keyCols, names(merged)),
                 setdiff(names(merged), keyCols))
    finish(merged[, ordered, drop = FALSE])
  })

# Build a GRanges from an overlap table: variants as width-1 ranges, all other
# columns (qtl_* / gwas_*) as mcols.
# @noRd
.overlapToGRanges <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(GenomicRanges::GRanges())
  p <- parseVariantId(df$variant_id)
  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", p$chrom),
    ranges = IRanges::IRanges(start = p$pos, width = 1L))
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(df, check.names = FALSE)
  gr
}
