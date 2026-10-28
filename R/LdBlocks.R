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

# Coerce an LD-block specification to the GRanges the splitters consume: a
# GRanges, or the path of a BED file, read as BED. The rules are those of
# every window argument (.asRegionRanges()).
# @noRd
.asLdBlockRanges <- function(x) {
    .asRegionRanges(x, "`ldBlocks`")
}
