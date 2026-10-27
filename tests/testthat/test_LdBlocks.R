# Tests for LD-block coercion.
#
# There was once an LdBlocks class wrapping a GRanges beside a `genome`
# string; it was retired because GRanges already carries the build in
# seqinfo() and its only consumer, LdStatistic, has its own genome slot. What
# remains is .asLdBlockRanges(), the one place an LD-block specification
# becomes the GRanges the splitters consume: a GRanges, or a BED file.

# @noRd
.ldb_blocks <- function() {
    b <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 800L), c(500L, 1200L))
    )
    names(b) <- c("b1", "b2")
    b
}

# @noRd
.ldb_bedFile <- function(lines) {
    path <- withr::local_tempfile(
        fileext = ".bed",
        .local_envir = parent.frame()
    )
    writeLines(lines, path)
    path
}

test_that(".asLdBlockRanges passes a GRanges through unchanged", {
    expect_identical(.asLdBlockRanges(.ldb_blocks()), .ldb_blocks())
})

test_that(".asLdBlockRanges reads a BED file with BED's 0-based start", {
    gr <- .asLdBlockRanges(.ldb_bedFile(c("chr1\t0\t500", "chr1\t799\t1200")))
    expect_s4_class(gr, "GRanges")
    expect_equal(start(gr), c(1L, 800L))
    expect_equal(end(gr), c(500L, 1200L))
    expect_null(mcols(gr)$blockId)
})

test_that(".asLdBlockRanges takes a BED file's fourth column as the names", {
    # Without this the splitter would fall back to coordinate keys, and the
    # file's own block names would be silently discarded.
    path <- .ldb_bedFile(c("chr1\t0\t500\ta", "chr1\t799\t1200\tb"))
    expect_equal(mcols(.asLdBlockRanges(path))$blockId, c("a", "b"))
})

test_that(".asLdBlockRanges passes over BED comments and a header line", {
    # An LD-block file often names its columns, as LDetect's `chr start stop`
    # does; BED's columns are positional, so the header line is skipped.
    path <- .ldb_bedFile(c(
        "track name=blocks",
        "# LD blocks",
        "chr\tstart\tstop",
        "chr1\t0\t500",
        ""
    ))
    expect_identical(as.character(.asLdBlockRanges(path)), "chr1:1-500")
})

test_that(".asLdBlockRanges refuses a file that is not BED", {
    expect_error(
        .asLdBlockRanges(.ldb_bedFile(c("chr1\t0\t500", "chr1 799 1200"))),
        "is not a BED file"
    )
    expect_error(.asLdBlockRanges(.ldb_bedFile(character(0))), "not a BED")
})

test_that(".asLdBlockRanges refuses a table, which R calls give 1-based", {
    # Given in R a window is a GRanges, so a table -- whose start could be
    # either base -- is refused rather than read.
    df <- data.frame(chrom = "chr1", start = c(1L, 800L), end = c(500L, 1200L))
    expect_error(
        .asLdBlockRanges(df),
        "`ldBlocks` must be a GRanges or the path of a BED file"
    )
})

test_that(".asLdBlockRanges names a path that does not exist", {
    expect_error(
        .asLdBlockRanges("chr1:1-500"),
        "no file exists at 'chr1:1-500'"
    )
})

test_that(".asLdBlockRanges rejects a type it cannot interpret", {
    expect_error(.asLdBlockRanges(1L), "must be a GRanges")
})
