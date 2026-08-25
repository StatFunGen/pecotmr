# =============================================================================
# LD block coercion
# -----------------------------------------------------------------------------
# LD block boundaries are a GRanges. There was once an LdBlocks class wrapping
# one alongside a `genome` string, but GRanges already carries the build in
# seqinfo(), and the only consumer -- LdStatistic -- has its own genome slot,
# so the build was stored twice. What remains is the coercion the splitters
# need, which already accepted a bare GRanges among its input forms.
# =============================================================================

#' @include AllGenerics.R
NULL

# Coerce an LD-block specification to the GRanges the splitters consume.
#
# Accepts a `GRanges`, a data.frame with chrom/start/end (plus an optional
# `blockId`), or a path to a delimited file holding those columns.
# @noRd
.asLdBlockRanges <- function(x) {
    if (methods::is(x, "GRanges")) {
        return(x)
    }
    if (is.character(x) && length(x) == 1L) {
        x <- .readManifest(x)
    }
    if (!is.data.frame(x)) {
        msg <- glue(
            "`ldBlocks` must be a GRanges, a data.frame with ",
            "chrom/start/end, or a path to one (got {class(x)[[1L]]})."
        )
        abort(msg)
    }
    .ldBlockTableToRanges(x)
}

# One block per row. `blockId` is carried into mcols when present so the key
# survives; without it the splitter falls back to the block's coordinates.
# @noRd
.ldBlockTableToRanges <- function(df) {
    missingCols <- setdiff(c("chrom", "start", "end"), colnames(df))
    if (length(missingCols) > 0L) {
        msg <- glue(
            "`ldBlocks` table is missing required column(s): ",
            "{str_flatten(missingCols, ', ')}."
        )
        abort(msg)
    }
    gr <- GenomicRanges::GRanges(
        seqnames = as.character(df$chrom),
        ranges = IRanges::IRanges(
            start = as.integer(df$start),
            end = as.integer(df$end)
        )
    )
    if (is_in("blockId", colnames(df))) {
        mcols(gr)$blockId <- as.character(df$blockId)
    }
    gr
}
