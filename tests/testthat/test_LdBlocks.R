# Tests for LD-block coercion.
#
# There was once an LdBlocks class wrapping a GRanges beside a `genome`
# string; it was retired because GRanges already carries the build in
# seqinfo() and its only consumer, LdStatistic, has its own genome slot. What
# remains is .asLdBlockRanges(), the one place an LD-block specification
# becomes the GRanges the splitters consume.

# @noRd
.ldb_blocks <- function() {
    b <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 800L), c(500L, 1200L))
    )
    names(b) <- c("b1", "b2")
    b
}

test_that(".asLdBlockRanges passes a GRanges through unchanged", {
    expect_identical(.asLdBlockRanges(.ldb_blocks()), .ldb_blocks())
})

test_that(".asLdBlockRanges builds ranges from a block table", {
    df <- data.frame(
        chrom = "chr1",
        start = c(1L, 800L),
        end = c(500L, 1200L),
        stringsAsFactors = FALSE
    )
    gr <- .asLdBlockRanges(df)
    expect_s4_class(gr, "GRanges")
    expect_equal(length(gr), 2L)
    expect_equal(start(gr), c(1L, 800L))
})

test_that(".asLdBlockRanges carries blockId from the table into mcols", {
    # Without this the splitter would fall back to coordinate keys, and the
    # caller's own block names would be silently discarded.
    df <- data.frame(
        chrom = "chr1",
        start = c(1L, 800L),
        end = c(500L, 1200L),
        blockId = c("a", "b"),
        stringsAsFactors = FALSE
    )
    expect_equal(mcols(.asLdBlockRanges(df))$blockId, c("a", "b"))
})

test_that(".asLdBlockRanges reads a block table from a path", {
    df <- data.frame(
        chrom = "chr1",
        start = c(1L, 800L),
        end = c(500L, 1200L),
        stringsAsFactors = FALSE
    )
    tf <- withr::local_tempfile(fileext = ".tsv")
    readr::write_tsv(df, tf, progress = FALSE)
    expect_equal(length(.asLdBlockRanges(tf)), 2L)
})

test_that(".asLdBlockRanges names the columns a block table is missing", {
    expect_error(
        .asLdBlockRanges(data.frame(chrom = "chr1", start = 1L)),
        "end"
    )
})

test_that(".asLdBlockRanges rejects a type it cannot interpret", {
    expect_error(.asLdBlockRanges(1L), "must be a GRanges")
})
