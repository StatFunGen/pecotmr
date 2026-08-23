# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (LdBlocks) ===

test_that("LdBlocks constructs and validates correctly", {
    obj <- make_test_ldblocks()
    expect_s4_class(obj, "LdBlocks")
    expect_equal(length(obj@blocks), 2)
    expect_equal(obj@genome, "hg19")
    expect_true(methods::validObject(obj))
})


test_that("LdBlocks rejects genome of length != 1", {
    blocks_gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 1, end = 5000)
    )
    expect_error(
        methods::validObject(
            new("LdBlocks", blocks = blocks_gr, genome = c("hg19", "hg38"))
        ),
        "genome.*single"
    )
})

# show() smoke test, moved here from test_showMethods.R so the test
# tree mirrors R/.
test_that("show(LdBlocks) does not error", {
    expect_output(show(make_test_ldblocks()), "LdBlocks")
})


# .asLdBlockRanges: the one place an LD-block specification becomes the GRanges
# the splitters consume. `LdBlocks` is the canonical type the rest of the
# package passes around; the looser forms exist so a caller holding a plain
# block table does not have to construct one first.

# @noRd
.ldb_blocks <- function() {
    b <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 800L), c(500L, 1200L))
    )
    names(b) <- c("b1", "b2")
    b
}

test_that(".asLdBlockRanges unwraps an LdBlocks", {
    lb <- new("LdBlocks", blocks = .ldb_blocks(), genome = "hg38")
    expect_identical(.asLdBlockRanges(lb), .ldb_blocks())
})

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
    expect_error(.asLdBlockRanges(1L), "must be an LdBlocks")
})
