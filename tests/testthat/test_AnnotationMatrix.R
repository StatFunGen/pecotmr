context("AnnotationMatrix")

# Migrated from test_h2Annotations.R: the AnnotationMatrix S4 class + validity,
# its constructor, the getGenome accessor, and the getBaseline /
# getCandidates tier accessors — all now defined in R/AnnotationMatrix.R.
# The class is a RangedSummarizedExperiment, so the annotation matrix,
# its per-column metadata and its ranges are read with assay(),
# colData() and rowRanges() rather than through accessors of our own.
# Annotation *readers* stay in test_h2Annotations.R. Fixtures come from
# helper-h2Classes.R.

test_that("AnnotationMatrix validates dimensions and meta", {
    gr <- make_test_granges(10)
    meta <- make_test_annotation_meta()
    mat <- matrix(0, nrow = 10, ncol = 3)

    obj <- AnnotationMatrix(mat, gr, meta)
    expect_s4_class(obj, "AnnotationMatrix")
    expect_true(methods::validObject(obj))
})

test_that("AnnotationMatrix rejects row mismatch", {
    gr <- make_test_granges(10)
    meta <- make_test_annotation_meta()
    mat <- matrix(0, nrow = 5, ncol = 3) # wrong number of rows

    expect_error(AnnotationMatrix(mat, gr, meta), "rows.*must match")
})

test_that("AnnotationMatrix rejects column mismatch with meta", {
    gr <- make_test_granges(10)
    meta <- make_test_annotation_meta() # 3 annotations
    mat <- matrix(0, nrow = 10, ncol = 2) # only 2 columns

    # Constructor errors when colnames assignment fails (dimnames mismatch)
    expect_error(AnnotationMatrix(mat, gr, meta))
})

test_that("AnnotationMatrix rejects invalid tier values", {
    gr <- make_test_granges(10)
    meta <- data.frame(
        name = "x",
        tier = "invalid_tier",
        type = "binary",
        stringsAsFactors = FALSE
    )
    mat <- matrix(0, nrow = 10, ncol = 1)

    expect_error(AnnotationMatrix(mat, gr, meta), "tier.*baseline.*candidate")
})

test_that("AnnotationMatrix rejects invalid type values", {
    gr <- make_test_granges(10)
    meta <- data.frame(
        name = "x",
        tier = "baseline",
        type = "ordinal",
        stringsAsFactors = FALSE
    )
    mat <- matrix(0, nrow = 10, ncol = 1)

    expect_error(AnnotationMatrix(mat, gr, meta), "type.*binary.*continuous")
})

test_that("AnnotationMatrix() constructor creates object from matrix", {
    n <- 10
    gr <- make_test_granges(n)
    meta <- make_test_annotation_meta()
    mat <- matrix(runif(n * 3), nrow = n, ncol = 3)

    obj <- AnnotationMatrix(mat, gr, meta, genome = "hg38")
    expect_s4_class(obj, "AnnotationMatrix")
    expect_equal(nrow(assay(obj, "annotations")), n)
    expect_equal(ncol(assay(obj, "annotations")), 3)
    expect_equal(getGenome(obj), "hg38")
})

test_that("AnnotationMatrix() sets column names from annotation_meta", {
    n <- 10
    gr <- make_test_granges(n)
    meta <- make_test_annotation_meta()
    mat <- matrix(0, nrow = n, ncol = 3) # No colnames set on mat

    obj <- AnnotationMatrix(mat, gr, meta)
    expect_equal(
        colnames(assay(obj, "annotations")),
        c("base", "enhancer", "promoter")
    )
})

test_that("AnnotationMatrix rejects a non-data.frame annotationMeta", {
    gr <- make_test_granges(10)
    expect_error(
        AnnotationMatrix(matrix(0, 10, 3), gr, annotationMeta = list(a = 1)),
        "annotationMeta must be a data.frame"
    )
})

test_that("AnnotationMatrix rejects annotationMeta missing required columns", {
    gr <- make_test_granges(10)
    bad_meta <- data.frame(foo = c("a", "b", "c"), stringsAsFactors = FALSE)
    expect_error(
        AnnotationMatrix(matrix(0, 10, 3), gr, annotationMeta = bad_meta),
        "must have columns: name, tier, type"
    )
})

test_that("AnnotationMatrix accessors round-trip the stored slots", {
    n <- 10
    gr <- make_test_granges(n)
    meta <- make_test_annotation_meta()
    mat <- matrix(runif(n * 3), nrow = n, ncol = 3)
    obj <- AnnotationMatrix(mat, gr, meta, genome = "hg38")

    expect_equal(dim(assay(obj, "annotations")), c(n, 3L))
    # colData is a DFrame keyed by the annotation names, so the plain
    # fixture compares as data after dropping those rownames.
    cd <- SummarizedExperiment::colData(obj)
    expect_equal(rownames(cd), meta$name)
    bare <- as.data.frame(cd)
    rownames(bare) <- NULL
    expect_equal(bare, meta)
    # The build now lives in seqinfo() rather than a slot, so the stored
    # ranges carry it while the input ranges did not.
    want <- gr
    GenomeInfoDb::genome(want) <- "hg38"
    expect_equal(rowRanges(obj), want)
    expect_equal(getGenome(obj), "hg38")
})

test_that("getBaseline() subsets to baseline-tier only", {
    n <- 10
    gr <- make_test_granges(n)
    meta <- make_test_annotation_meta() # 1 baseline, 2 candidate
    mat <- matrix(0, nrow = n, ncol = 3)

    obj <- AnnotationMatrix(mat, gr, meta)
    baseline <- getBaseline(obj)

    expect_s4_class(baseline, "AnnotationMatrix")
    expect_equal(ncol(assay(baseline, "annotations")), 1)
    expect_true(all(SummarizedExperiment::colData(baseline)$tier == "baseline"))
})

test_that("getCandidates() subsets to candidate-tier only", {
    n <- 10
    gr <- make_test_granges(n)
    meta <- make_test_annotation_meta() # 1 baseline, 2 candidate
    mat <- matrix(0, nrow = n, ncol = 3)

    obj <- AnnotationMatrix(mat, gr, meta)
    cand <- getCandidates(obj)

    expect_s4_class(cand, "AnnotationMatrix")
    expect_equal(ncol(assay(cand, "annotations")), 2)
    expect_true(all(SummarizedExperiment::colData(cand)$tier == "candidate"))
})

# show() smoke test, moved here from test_showMethods.R so the test
# tree mirrors R/.
test_that("show(AnnotationMatrix) does not error", {
    am <- AnnotationMatrix(
        matrix(0, nrow = 10, ncol = 1),
        make_test_granges(10),
        data.frame(
            name = "base",
            tier = "baseline",
            type = "binary",
            stringsAsFactors = FALSE
        )
    )
    expect_output(show(am), "AnnotationMatrix")
})

test_that("the constructor rejects a column/metadata-row mismatch", {
    # Rows-vs-ranges is checked above; this is the other axis -- one metadata
    # row per annotation COLUMN.
    n <- 10
    gr <- make_test_granges(n)
    mat <- matrix(runif(n * 4), nrow = n, ncol = 4)
    expect_error(
        AnnotationMatrix(mat, gr, make_test_annotation_meta()),
        "column\\(s\\) for 3 metadata row\\(s\\)"
    )
})

test_that("validity names the metadata columns it requires", {
    n <- 10
    obj <- AnnotationMatrix(
        matrix(runif(n * 3), nrow = n, ncol = 3),
        make_test_granges(n),
        make_test_annotation_meta()
    )
    bad <- obj
    SummarizedExperiment::colData(bad)$tier <- NULL
    expect_error(
        methods::validObject(bad),
        "must have columns: name, tier, type"
    )
})
