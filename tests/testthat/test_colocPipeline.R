context("colocPipeline")

# ===========================================================================
# Strategy: mock coloc::coloc.bf_bf to return a tiny fake summary, then
# drive the QTL/GWAS pairing loop on a small fixture. Mock
# fineMappingPipeline so the GwasSumStats input path also runs without the
# heavy susie fits.
# ===========================================================================

.cp_makeHandle <- function(snp_n = 6L, n_samples = 30L, sample_prefix = "s") {
    new(
        "GenotypeHandle",
        path = "/tmp/sketch.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(snp_n)),
            CHR = rep("1", snp_n),
            BP = seq(100L, by = 100L, length.out = snp_n),
            A1 = rep("A", snp_n),
            A2 = rep("G", snp_n),
            stringsAsFactors = FALSE
        ),
        nSamples = n_samples,
        sampleIds = paste0(sample_prefix, seq_len(n_samples)),
        pgenPtr = NULL
    )
}

.cp_makeFmEntry <- function(
    variant_ids = paste0("chr1:", 100 * (1:5), ":A:G"),
    withLbf = TRUE,
    n_eff = 2L
) {
    pip <- seq(0.9, by = -0.15, length.out = length(variant_ids))
    n <- length(variant_ids)
    tl <- data.frame(
        variant_id = variant_ids,
        chrom = rep("1", n),
        pos = as.integer(100 * seq_len(n)),
        A1 = rep("G", n),
        A2 = rep("A", n),
        N = rep(1000, n),
        MAF = rep(0.1, n),
        marginal_beta = rep(0.1, n),
        marginal_se = rep(0.05, n),
        marginal_z = rep(2.0, n),
        marginal_p = rep(0.05, n),
        pip = pip,
        posterior_mean = rep(0.05, n),
        posterior_sd = rep(0.02, n),
        stringsAsFactors = FALSE
    )
    fit <- list(
        alpha = matrix(
            1 / length(variant_ids),
            nrow = n_eff,
            ncol = length(variant_ids),
            dimnames = list(NULL, variant_ids)
        ),
        pip = setNames(pip, variant_ids),
        V = rep(0.05, n_eff)
    )
    if (withLbf) {
        fit$lbf_variable <- matrix(
            rnorm(n_eff * length(variant_ids)),
            nrow = n_eff,
            ncol = length(variant_ids),
            dimnames = list(NULL, variant_ids)
        )
    }
    fineMappingRow(variantIds = variant_ids, susieFit = fit, topLoci = tl)
}

.cp_makeQtlFmr <- function(
    tuples = list(c("Q1", "c1", "t1", "susie")),
    entries = NULL,
    with_sketch = TRUE
) {
    if (is.null(entries)) {
        entries <- replicate(
            length(tuples),
            .cp_makeFmEntry(),
            simplify = FALSE
        )
    }
    QtlFineMappingResult(
        study = map_chr(tuples, 1L),
        context = map_chr(tuples, 2L),
        trait = map_chr(tuples, 3L),
        method = map_chr(tuples, 4L),
        entry = entries,
        ldSketch = if (with_sketch) .cp_makeHandle() else NULL
    )
}

.cp_makeGwasFmr <- function(
    tuples = list(c("G1", "susie")),
    entries = NULL,
    with_sketch = TRUE
) {
    if (is.null(entries)) {
        entries <- replicate(
            length(tuples),
            .cp_makeFmEntry(),
            simplify = FALSE
        )
    }
    GwasFineMappingResult(
        study = map_chr(tuples, 1L),
        method = map_chr(tuples, 2L),
        entry = entries,
        ldSketch = if (with_sketch) .cp_makeHandle() else NULL
    )
}

.cp_makeGwasSumstats <- function(study = "G1", qc = TRUE) {
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = 5L),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:5)),
        A1 = rep("A", 5),
        A2 = rep("G", 5),
        Z = rnorm(5),
        N = rep(1000L, 5)
    )
    GwasSumStats(
        study = study,
        entry = list(gr),
        genome = "hg19",
        ldSketch = .cp_makeHandle(),
        qcInfo = if (qc) list(step1 = "ok") else list()
    )
}

.cp_makeQtlSumstats <- function(study = "Q1", qc = TRUE) {
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = 5L),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:5)),
        A1 = rep("A", 5),
        A2 = rep("G", 5),
        Z = rnorm(5),
        N = rep(1000L, 5)
    )
    QtlSumStats(
        study = study,
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .cp_makeHandle(),
        qcInfo = if (qc) list(step1 = "ok") else list()
    )
}

# The enrichment lookup takes each side's identity list, as the scoring loop
# builds it.
.cp_side <- function(study, context = NA_character_, trait = NA_character_) {
    list(study = study, context = context, trait = trait)
}

.cp_mockColocBfBf <- function() {
    function(qLbf, gLbf, p1, p2, p12, ...) {
        list(
            summary = data.frame(
                idx1 = 1L,
                idx2 = 1L,
                nSnps = ncol(qLbf),
                PP.H0.abf = 0.1,
                PP.H1.abf = 0.2,
                PP.H2.abf = 0.2,
                PP.H3.abf = 0.2,
                PP.H4.abf = 0.3,
                stringsAsFactors = FALSE
            )
        )
    }
}

# ===========================================================================
# Input-type validation
# ===========================================================================

test_that("colocPipeline: rejects a non-fine-mapping first side", {
    expect_error(
        colocPipeline(
            qtlFineMappingResult = "no",
            gwasInput = .cp_makeGwasFmr()
        ),
        "must be a QtlFineMappingResult or a GwasFineMappingResult"
    )
})

test_that("colocPipeline: rejects an unusable second side", {
    expect_error(
        colocPipeline(qtlFineMappingResult = .cp_makeQtlFmr(), gwasInput = 42L),
        "must be a fine-mapping result"
    )
})

test_that("colocPipeline: rejects un-QCd GwasSumStats input", {
    qfmr <- .cp_makeQtlFmr()
    gss <- .cp_makeGwasSumstats(qc = FALSE)
    expect_error(
        colocPipeline(qtlFineMappingResult = qfmr, gwasInput = gss),
        "has no QC record"
    )
})

# ===========================================================================
# .colocRequireMatchingLdSketches
# ===========================================================================

test_that(".colocRequireMatchingLdSketches: NULL qtl-side ldSketch is allowed", {
    qfmr <- .cp_makeQtlFmr(with_sketch = FALSE)
    gfmr <- .cp_makeGwasFmr()
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gfmr
    ))
    expect_s4_class(out, "ColocResult")
})

test_that(".colocRequireMatchingLdSketches: non-NULL qtl + NULL gwas errors", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr(with_sketch = FALSE)
    expect_error(
        colocPipeline(qtlFineMappingResult = qfmr, gwasInput = gfmr),
        "ldSketch is NULL"
    )
})

test_that(".colocRequireMatchingLdSketches: sample set mismatch errors", {
    qfmr <- .cp_makeQtlFmr()
    otherPanel <- .cp_makeHandle(sample_prefix = "other")
    gfmr <- GwasFineMappingResult(
        study = "G1",
        method = "susie",
        entry = list(.cp_makeFmEntry()),
        ldSketch = otherPanel
    )
    expect_error(
        colocPipeline(qtlFineMappingResult = qfmr, gwasInput = gfmr),
        "different sample sets"
    )
})

test_that(".colocRequireMatchingLdSketches: differently trimmed panels pass", {
    # The two sides QC'd separately keep different subsets of one LD
    # reference; the panels overlap, so the check warns rather than aborting.
    expect_warning(
        expect_null(pecotmr:::.colocRequireMatchingLdSketches(
            .cp_makeHandle(snp_n = 6L),
            .cp_makeHandle(snp_n = 7L)
        )),
        "share 6 variant"
    )
})

# ===========================================================================
# End-to-end with mocked coloc.bf_bf
# ===========================================================================

test_that("colocPipeline: returns one row per (QTL tuple, GWAS tuple) pair", {
    qfmr <- .cp_makeQtlFmr(
        tuples = list(
            c("Q1", "c1", "t1", "susie"),
            c("Q1", "c2", "t1", "susie")
        )
    )
    gfmr <- .cp_makeGwasFmr(tuples = list(c("G1", "susie"), c("G2", "susie")))
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gfmr
    ))
    expect_equal(nrow(out), 4L) # 2 QTL tuples * 2 GWAS tuples
    expect_setequal(out$study, "Q1")
    expect_setequal(out$context, c("c1", "c2"))
    expect_setequal(out$gwasStudy, c("G1", "G2"))
})

# ===========================================================================
# Either side may be a QTL or a GWAS fine-mapping result
# ===========================================================================

test_that("colocPipeline: pairs two QTL fine-mapping results", {
    qfmr <- .cp_makeQtlFmr()
    # The two second-side rows share (study, method, block) and differ only on
    # trait -- keyed on the GWAS 2-tuple, one would silently replace the other.
    other <- .cp_makeQtlFmr(
        tuples = list(
            c("Q2", "c1", "t1", "susie"),
            c("Q2", "c1", "t2", "susie")
        )
    )
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = other
    ))
    expect_equal(nrow(out), 2L)
    expect_setequal(out$gwasStudy, "Q2")
    expect_setequal(out$gwasContext, "c1")
    expect_setequal(out$gwasTrait, c("t1", "t2"))
})

test_that("colocPipeline: pairs two GWAS fine-mapping results", {
    first <- .cp_makeGwasFmr(tuples = list(c("G1", "susie")))
    second <- .cp_makeGwasFmr(tuples = list(c("G2", "susie")))
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = first,
        gwasInput = second
    ))
    expect_equal(nrow(out), 1L)
    expect_equal(as.character(out$study), "G1")
    expect_equal(as.character(out$gwasStudy), "G2")
    # Neither side has a context or trait axis, so all four are reported as NA
    # rather than invented.
    expect_true(all(is.na(c(
        out$context,
        out$trait,
        out$gwasContext,
        out$gwasTrait
    ))))
})

test_that("colocPipeline: a GWAS first side names itself in warnings", {
    first <- .cp_makeGwasFmr()
    second <- .cp_makeGwasFmr(tuples = list(c("G2", "susie")))
    local_mocked_bindings(
        coloc.bf_bf = function(...) stop("synthetic test failure"),
        .package = "coloc"
    )
    expect_warning(
        colocPipeline(qtlFineMappingResult = first, gwasInput = second),
        "GWAS \\(study='G1', method='susie'\\)"
    )
})

test_that("colocPipeline: resolves QtlSumStats by fine-mapping it", {
    first <- .cp_makeGwasFmr()
    qss <- .cp_makeQtlSumstats()
    resolved <- .cp_makeQtlFmr()
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    local_mocked_bindings(
        fineMappingPipeline = function(data, methods, ...) resolved,
        .package = "pecotmr"
    )
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = first,
        gwasInput = qss,
        returnGwasFineMapping = TRUE,
        adjustPips = FALSE
    ))
    expect_equal(as.character(out$gwasTrait), "t1")
    expect_identical(attr(out, "gwasFineMapping"), resolved)
})

test_that("colocPipeline: rejects an un-QCd QtlSumStats second side", {
    expect_error(
        colocPipeline(
            qtlFineMappingResult = .cp_makeQtlFmr(),
            gwasInput = .cp_makeQtlSumstats(qc = FALSE)
        ),
        "has no QC record"
    )
})

test_that("colocPipeline: resolves GwasSumStats via fineMappingPipeline (mocked)", {
    qfmr <- .cp_makeQtlFmr()
    gss <- .cp_makeGwasSumstats()
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    local_mocked_bindings(
        fineMappingPipeline = function(data, methods, ...) .cp_makeGwasFmr(),
        .package = "pecotmr"
    )
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gss
    ))
    expect_s4_class(out, "ColocResult")
    expect_gte(nrow(out), 1L)
})

test_that("colocPipeline: returnGwasFineMapping attaches the resolved FMR", {
    qfmr <- .cp_makeQtlFmr()
    gss <- .cp_makeGwasSumstats()
    resolvedGfmr <- .cp_makeGwasFmr()
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    local_mocked_bindings(
        fineMappingPipeline = function(data, methods, ...) resolvedGfmr,
        .package = "pecotmr"
    )
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gss,
        returnGwasFineMapping = TRUE,
        adjustPips = FALSE
    ))
    expect_identical(attr(out, "gwasFineMapping"), resolvedGfmr)
})

test_that("colocPipeline: coloc.bf_bf failure surfaces as warning and skip", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr()
    local_mocked_bindings(
        coloc.bf_bf = function(...) stop("synthetic test failure"),
        .package = "coloc"
    )
    expect_warning(
        out <- colocPipeline(qtlFineMappingResult = qfmr, gwasInput = gfmr),
        "coloc.bf_bf failed"
    )
    expect_equal(nrow(out), 0L)
})

test_that("colocPipeline: disjoint QTL and GWAS variant sets error", {
    # Returning an empty table here is indistinguishable from "ran fine, no
    # colocalization" -- a materially different scientific conclusion -- so
    # fully disjoint inputs are refused rather than silently emptied.
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr(
        entries = list(.cp_makeFmEntry(
            variant_ids = paste0("chr9:", 100 * (1:5), ":A:G")
        ))
    )
    expect_error(
        colocPipeline(qtlFineMappingResult = qfmr, gwasInput = gfmr),
        "share no variant|disjoint"
    )
})

test_that("colocPipeline: empty result has the documented schema", {
    qfmr <- .cp_makeQtlFmr()
    # Build a GWAS FMR whose entry has no usable LBF -> pre-extract returns
    # empty. The variant OVERLAPS the QTL fixture: a disjoint pair is now an
    # error (see the disjoint test below), so it cannot be used to reach the
    # empty-schema path.
    shared <- "chr1:100:A:G"
    emptyFit <- list(
        alpha = matrix(1, 1, 1),
        pip = 1,
        V = 0,
        lbf_variable = matrix(NA_real_, 1, 1)
    )
    e <- fineMappingRow(
        variantIds = shared,
        susieFit = emptyFit,
        topLoci = data.frame(
            variant_id = shared,
            pip = 1,
            stringsAsFactors = FALSE
        )
    )
    gfmr <- GwasFineMappingResult(
        study = "G1",
        method = "susie",
        entry = list(e),
        ldSketch = .cp_makeHandle()
    )
    out <- suppressWarnings(
        colocPipeline(qtlFineMappingResult = qfmr, gwasInput = gfmr)
    )
    expect_equal(nrow(out), 0L)
    expect_setequal(
        colnames(out),
        c(
            "study",
            "context",
            "trait",
            "method",
            "gwasStudy",
            "gwasContext",
            "gwasTrait",
            "gwasMethod",
            "blockId",
            "qtlCs",
            "gwasCs",
            "idx1",
            "idx2",
            "nSnps",
            "PP.H0.abf",
            "PP.H1.abf",
            "PP.H2.abf",
            "PP.H3.abf",
            "PP.H4.abf",
            "qtlRetainedMass",
            "gwasRetainedMass"
        )
    )
})

# ===========================================================================
# Internal helpers
# ===========================================================================

test_that(".colocExtractLbfFromEntry: entry without trimmedFit returns NULL with warning", {
    e <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = NULL,
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.1,
            stringsAsFactors = FALSE
        )
    )
    expect_warning(
        out <- pecotmr:::.colocExtractLbfFromEntry(
            e,
            ColocLbfFilterParam()
        ),
        "has no trimmedFit"
    )
    expect_null(out)
})

test_that(".colocExtractLbfFromEntry: filterLbfCs subsets by cs_index", {
    fit <- list(
        alpha = matrix(
            0,
            3,
            4,
            dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * (1:4)))
        ),
        pip = setNames(
            c(0.9, 0.1, 0.5, 0.2),
            sprintf("chr1:%d:A:G", 100L * (1:4))
        ),
        V = c(0.1, 0.1, 0.1),
        lbf_variable = matrix(
            1:12,
            3,
            4,
            dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * (1:4)))
        ),
        sets = list(cs_index = c(1L, 3L))
    ) # keep effects 1 and 3
    e <- fineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        susieFit = fit,
        topLoci = data.frame(
            variant_id = sprintf("chr1:%d:A:G", 100L * (1:4)),
            pip = c(0.9, 0.1, 0.5, 0.2),
            stringsAsFactors = FALSE
        )
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        lbfFilterArgs = ColocLbfFilterParam(
            filterLbfCs = TRUE,
            secondary = NULL
        )
    )
    expect_equal(nrow(out$lbf), 2L)
})

test_that(".colocAlignLbf: aligned matrices share the common variant set", {
    q <- matrix(
        0,
        2,
        4,
        dimnames = list(
            NULL,
            c("chr1:10:A:G", "chr1:20:A:G", "chr1:30:A:G", "chr1:40:A:G")
        )
    )
    g <- matrix(
        0,
        2,
        3,
        dimnames = list(NULL, c("chr1:20:A:G", "chr1:30:A:G", "chr1:50:A:G"))
    )
    aligned <- pecotmr:::.colocAlignLbf(q, g)
    expect_setequal(colnames(aligned$qtl), c("chr1:20:A:G", "chr1:30:A:G"))
    expect_setequal(colnames(aligned$gwas), c("chr1:20:A:G", "chr1:30:A:G"))
})

test_that(".colocAlignLbf aligns non-autosomal columns across a chr-prefix difference", {
    # chrX ids: alignVariantNames' numeric-chrom regex could not parse these and
    # fell back to a raw intersect (no overlap); tuple matching resolves them.
    q <- matrix(
        0,
        2,
        2,
        dimnames = list(NULL, c("chrX:100:A:G", "chrX:200:C:T"))
    )
    g <- matrix(0, 2, 2, dimnames = list(NULL, c("X:100:A:G", "X:200:C:T")))
    aligned <- pecotmr:::.colocAlignLbf(q, g)
    expect_false(is.null(aligned))
    expect_equal(ncol(aligned$qtl), 2L)
    expect_identical(colnames(aligned$qtl), colnames(aligned$gwas))
})

test_that(".colocStandardiseRow: fills in missing PP columns with NA", {
    row <- data.frame(idx1 = 1L, stringsAsFactors = FALSE)
    out <- pecotmr:::.colocStandardiseRow(row)
    expect_true(all(
        c("idx2", "nSnps", paste0("PP.H", 0:4, ".abf")) %in% colnames(out)
    ))
})


context("encoloc")

# The file-path colocalization / enrichment wrappers (colocWrapper,
# xqtlEnrichmentWrapper, colocPostProcessor) and their helpers
# (filterAndOrderColocResults, calculateCumsum, calculate_purity,
# processColocResults, extract_ld_for_variants) have been removed in
# favor of the S4 colocPipeline / qtlEnrichmentPipeline / enlocPipeline
# entry points. The previous tests against the file-path wrappers no
# longer apply; they relied on mocking rssAnalysisPipeline (also
# removed). New tests for colocPipeline / qtlEnrichmentPipeline /
# enlocPipeline live alongside their pipeline implementations.

# ===========================================================================
# .colocLookupEnrichment (formerly .enlocLookupEnrichment, now shared)
# ===========================================================================

test_that(".colocLookupEnrichment: returns the value for a (gwasStudy, qtlStudy, qtlContext) hit", {
    enr <- data.frame(
        gwasStudy = c("G1", "G2"),
        qtlStudy = c("Q1", "Q1"),
        qtlContext = c("c1", "c1"),
        enrichment = c(2.0, 3.5),
        stringsAsFactors = FALSE
    )
    # Returns the matching ROW, so the caller can read the enloc priors off
    # it rather than only the enrichment factor.
    expect_equal(
        pecotmr:::.colocLookupEnrichment(
            enr,
            .cp_side("G2"),
            .cp_side("Q1", "c1")
        ),
        2L
    )
})

test_that(".colocLookupEnrichment: returns NA when no row matches", {
    enr <- data.frame(
        gwasStudy = "G1",
        qtlStudy = "Q1",
        qtlContext = "c1",
        enrichment = 2.0,
        stringsAsFactors = FALSE
    )
    expect_true(is.na(pecotmr:::.colocLookupEnrichment(
        enr,
        .cp_side("ghost"),
        .cp_side("Q1", "c1")
    )))
    # qtlStudy mismatch also a miss.
    expect_true(is.na(pecotmr:::.colocLookupEnrichment(
        enr,
        .cp_side("G1"),
        .cp_side("Qghost", "c1")
    )))
})

# ===========================================================================
# .colocEmptyResult(enriched = TRUE) — formerly .enlocEmptyResult
# ===========================================================================

test_that(".colocEmptyResult(enriched=TRUE): includes enrichment + p12Used schema", {
    out <- pecotmr:::.colocEmptyResult(enriched = TRUE)
    expect_s4_class(out, "ColocResult")
    expect_equal(nrow(out), 0L)
    expect_true(all(c("enrichment", "p12Used") %in% colnames(out)))
})


# ===========================================================================
# Additional coverage (appended)
# ===========================================================================

# Minimal canonical topLoci (variant_id + pip) for skeletal entries.
.cp_tl <- function(vids, pip = rep(0.1, length(vids))) {
    data.frame(
        variant_id = as.character(vids),
        pip = pip,
        stringsAsFactors = FALSE
    )
}

# --- enrichment= argument validation ---------------------------------------

test_that("colocPipeline: rejects a non-data.frame enrichment", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr()
    expect_error(
        colocPipeline(
            qtlFineMappingResult = qfmr,
            gwasInput = gfmr,
            enrichment = "not a data frame"
        ),
        "Must be of type 'data.frame'"
    )
})

test_that("colocPipeline: rejects enrichment missing required columns", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr()
    bad <- data.frame(
        gwasStudy = "G1",
        qtlStudy = "Q1",
        stringsAsFactors = FALSE
    ) # missing qtlContext + enrichment
    expect_error(
        colocPipeline(
            qtlFineMappingResult = qfmr,
            gwasInput = gfmr,
            enrichment = bad
        ),
        "Colnames must include the elements"
    )
})

# --- enrichment end-to-end (lookup + p12 scaling + output columns) ----------

test_that("colocPipeline: enrichment mode takes its priors from the table", {
    qfmr <- .cp_makeQtlFmr() # Q1 / c1 / t1 / susie
    gfmr <- .cp_makeGwasFmr() # G1 / susie
    enr <- data.frame(
        gwasStudy = "G1",
        qtlStudy = "Q1",
        qtlContext = "c1",
        enrichment = 2.0,
        colocP1 = 3e-4,
        colocP2 = 7e-4,
        colocP12 = 2.5e-5,
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    out <- suppressWarnings(
        colocPipeline(
            qtlFineMappingResult = qfmr,
            gwasInput = gfmr,
            enrichment = enr
        )
    )
    # Enrichment mode scores every (QTL effect, GWAS effect) pair itself
    # rather than calling coloc.bf_bf, so the fixture's 2 x 2 effects give
    # four rows; the mock used elsewhere collapses them to one.
    expect_equal(nrow(out), 4L)
    expect_true(all(c("enrichment", "p12Used") %in% colnames(out)))
    expect_equal(unique(out$enrichment), 2.0)
    # qtlEnrichmentPipeline derives p12 = P_eqtl * expit(a0 + a1) the way
    # fastenloc does; the pipeline uses that, rather than scaling a default.
    expect_equal(unique(out$p12Used), 2.5e-5)
})

test_that("enrichment mode reports RCP and LCP alongside the hypotheses", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr()
    enr <- data.frame(
        gwasStudy = "G1",
        qtlStudy = "Q1",
        qtlContext = "c1",
        enrichment = 2.0,
        colocP1 = 3e-4,
        colocP2 = 7e-4,
        colocP12 = 2.5e-5,
        stringsAsFactors = FALSE
    )
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gfmr,
        enrichment = enr
    ))
    pairs <- getColocPairs(out)
    expect_true(all(c("RCP", "LCP") %in% colnames(pairs)))
    # The two names for the same quantities: fastenloc's RCP is the posterior
    # of one shared causal variant, and LCP adds the distinct-variant case.
    expect_equal(pairs$RCP, pairs$PP.H4.abf)
    expect_equal(pairs$LCP, pairs$PP.H3.abf + pairs$PP.H4.abf)
    probs <- pairs$PP.H0.abf +
        pairs$PP.H1.abf +
        pairs$PP.H2.abf +
        pairs$PP.H3.abf +
        pairs$PP.H4.abf
    expect_equal(probs, rep(1, nrow(pairs)))
})

test_that("colocPipeline: an enrichment table without the priors is refused", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr()
    enr <- data.frame(
        gwasStudy = "G1",
        qtlStudy = "Q1",
        qtlContext = "c1",
        enrichment = 2.0,
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    expect_error(
        suppressWarnings(colocPipeline(
            qtlFineMappingResult = qfmr,
            gwasInput = gfmr,
            enrichment = enr
        )),
        "missing colocP1, colocP2, colocP12"
    )
})

test_that("colocPipeline: enrichment joins on the second side's trait", {
    # A QTL second side puts two traits under one study. Joining on the study
    # alone would give both the first row's factor.
    qfmr <- .cp_makeQtlFmr()
    other <- .cp_makeQtlFmr(
        tuples = list(
            c("Q2", "c1", "t1", "susie"),
            c("Q2", "c1", "t2", "susie")
        )
    )
    enr <- data.frame(
        gwasStudy = c("Q2", "Q2"),
        gwasContext = c("c1", "c1"),
        gwasTrait = c("t1", "t2"),
        qtlStudy = c("Q1", "Q1"),
        qtlContext = c("c1", "c1"),
        enrichment = c(1.0, 3.0),
        colocP1 = c(3e-4, 3e-4),
        colocP2 = c(7e-4, 7e-4),
        colocP12 = c(1e-5, 3e-5),
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = other,
        enrichment = enr
    ))
    pairs <- getColocPairs(out)
    # Four effect pairs per trait now that enrichment mode scores each one,
    # so the check is that each trait carries its own factor throughout.
    expect_equal(
        unique(pairs$enrichment[pairs$gwasTrait == "t1"]),
        1.0
    )
    expect_equal(
        unique(pairs$enrichment[pairs$gwasTrait == "t2"]),
        3.0
    )
})

test_that("colocPipeline: rejects an enrichment table with repeated keys", {
    enr <- data.frame(
        gwasStudy = c("G1", "G1"),
        qtlStudy = c("Q1", "Q1"),
        qtlContext = c("c1", "c1"),
        enrichment = c(2.0, 3.0),
        stringsAsFactors = FALSE
    )
    expect_error(
        colocPipeline(
            qtlFineMappingResult = .cp_makeQtlFmr(),
            gwasInput = .cp_makeGwasFmr(),
            enrichment = enr
        ),
        "repeated"
    )
})

test_that("colocPipeline: an enrichment miss warns and drops to plain priors", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr()
    enr <- data.frame(
        gwasStudy = "OTHER", # no row matches gwasStudy 'G1'
        qtlStudy = "Q1",
        qtlContext = "c1",
        enrichment = 2.0,
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    w <- testthat::capture_warnings(
        out <- colocPipeline(
            qtlFineMappingResult = qfmr,
            gwasInput = gfmr,
            enrichment = enr
        )
    )
    expect_match(w, "no enrichment entry", all = FALSE)
    expect_equal(nrow(out), 4L)
    # A miss means there is no enrichment estimate for this pair, so it is
    # scored with the unconditional priors rather than a fabricated factor
    # of zero applied to a default.
    expect_true(is.na(unique(out$enrichment)))
    expect_equal(unique(out$p12Used), 5e-6)
})

# --- returnGwasFineMapping attached on the empty-pair early return (~205) ----

test_that("colocPipeline: attaches gwasFineMapping when no pairs survive", {
    qfmr <- .cp_makeQtlFmr()
    gss <- .cp_makeGwasSumstats()
    # Resolved GWAS FMR whose only entry has no usable LBF (V filtered to 0
    # rows) -> .colocPreextractGwasLbf returns an empty list -> early return.
    emptyFit <- list(
        alpha = matrix(0, 1, 1),
        pip = c(v1 = 0),
        V = 0,
        lbf_variable = matrix(NA_real_, 1, 1)
    )
    e <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = emptyFit,
        topLoci = .cp_tl("chr1:100:A:G", pip = 0)
    )
    resolved <- GwasFineMappingResult(
        study = "G1",
        method = "susie",
        entry = list(e),
        ldSketch = .cp_makeHandle()
    )
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    local_mocked_bindings(
        fineMappingPipeline = function(data, methods, ...) resolved,
        .package = "pecotmr"
    )
    out <- suppressWarnings(
        colocPipeline(
            qtlFineMappingResult = qfmr,
            gwasInput = gss,
            returnGwasFineMapping = TRUE,
            adjustPips = FALSE
        )
    )
    expect_equal(nrow(out), 0L)
    expect_identical(attr(out, "gwasFineMapping"), resolved)
})

# --- per-pair loop skips: no usable QTL LBF (~223) and no overlap (~233) -----

test_that("colocPipeline: skips a QTL entry with no usable LBF", {
    badFit <- list(
        alpha = matrix(0, 1, 1),
        pip = c(v1 = 0),
        V = 0,
        lbf_variable = matrix(NA_real_, 1, 1)
    )
    badEntry <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = badFit,
        topLoci = .cp_tl("chr1:100:A:G", pip = 0)
    )
    goodEntry <- .cp_makeFmEntry()
    qfmr <- .cp_makeQtlFmr(
        tuples = list(
            c("Q1", "c1", "t1", "susie"),
            c("Q1", "c2", "t1", "susie")
        ),
        entries = list(badEntry, goodEntry)
    )
    gfmr <- .cp_makeGwasFmr()
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    out <- suppressWarnings(
        colocPipeline(
            qtlFineMappingResult = qfmr,
            gwasInput = gfmr,
            adjustPips = FALSE
        )
    )
    # bad entry (c1) skipped; only the good entry (c2) yields a row.
    expect_equal(nrow(out), 1L)
    expect_equal(unique(out$context), "c2")
})

test_that("colocPipeline: skips a pair with no shared variants", {
    qEntry <- .cp_makeFmEntry(
        variant_ids = paste0("chr1:", 100 * (1:5), ":A:G")
    )
    gEntry <- .cp_makeFmEntry(
        variant_ids = paste0("chr2:", 100 * (1:5), ":A:G")
    )
    qfmr <- .cp_makeQtlFmr(entries = list(qEntry))
    gfmr <- .cp_makeGwasFmr(entries = list(gEntry))
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    out <- suppressWarnings(
        colocPipeline(
            qtlFineMappingResult = qfmr,
            gwasInput = gfmr,
            adjustPips = FALSE
        )
    )
    expect_equal(nrow(out), 0L)
})

# --- .colocFilterCsByConcentration (~324-331) -------------------------------

test_that(".colocFilterCsByConcentration: keeps narrow CS, drops diffuse ones", {
    fit <- list(alpha = matrix(0.1, nrow = 3, ncol = 10)) # 10 variants
    # maxSize = ncol * coverage * concentration = 10 * 0.5 * 0.5 = 2.5
    local_mocked_bindings(
        susie_get_cs = function(s, coverage = 0.5, dedup = TRUE, ...) {
            list(
                cs = list(
                    L1 = c(1L, 2L), # size 2  -> keep
                    L2 = 1:8, # size 8  -> drop
                    L3 = 3L
                )
            )
        }, # size 1  -> keep
        .package = "pecotmr"
    )
    keep <- pecotmr:::.colocFilterCsByConcentration(
        fit,
        coverage = 0.5,
        concentration = 0.5
    )
    expect_setequal(keep, c(1, 3))
})

# --- .colocExtractLbfFromEntry branches -------------------------------------

test_that(".colocExtractLbfFromEntry: stacks fSuSiE lBF list into a matrix", {
    m1 <- matrix(
        rnorm(8),
        2,
        4,
        dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    m2 <- matrix(
        rnorm(8),
        2,
        4,
        dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    fit <- list(fsusie_result = list(lBF = list(m1, m2)))
    e <- fineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        susieFit = fit,
        topLoci = .cp_tl(sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        lbfFilterArgs = ColocLbfFilterParam(
            filterLbfCs = FALSE,
            secondary = NULL,
            concentration = 0.5
        )
    )
    expect_equal(nrow(out$lbf), 4L) # 2 + 2 stacked
    expect_setequal(colnames(out$lbf), sprintf("chr1:%d:A:G", 100L * (1:4)))
})

test_that(".colocExtractLbfFromEntry: stacks nested fSuSiE lBF (fit[[1]] path)", {
    m1 <- matrix(
        rnorm(8),
        2,
        4,
        dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    fit <- list(list(fsusie_result = list(lBF = list(m1))))
    e <- fineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        susieFit = fit,
        topLoci = .cp_tl(sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        lbfFilterArgs = ColocLbfFilterParam(
            filterLbfCs = FALSE,
            secondary = NULL,
            concentration = 0.5
        )
    )
    expect_equal(nrow(out$lbf), 2L)
})

test_that(".colocExtractLbfFromEntry: warns + NULL when fit carries no LBF slot", {
    fit <- list(notLbf = list(a = 1)) # fit[[1]] is a list -> $ access is safe
    e <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = fit,
        topLoci = .cp_tl("chr1:100:A:G")
    )
    expect_warning(
        out <- pecotmr:::.colocExtractLbfFromEntry(
            e,
            lbfFilterArgs = ColocLbfFilterParam(
                filterLbfCs = FALSE,
                secondary = NULL,
                concentration = 0.5
            )
        ),
        "no lbf_variable"
    )
    expect_null(out)
})

test_that(".colocExtractLbfFromEntry: warns + NULL on an empty LBF matrix", {
    fit <- list(
        lbf_variable = matrix(
            numeric(0),
            nrow = 0,
            ncol = 3,
            dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * (1:3)))
        )
    )
    e <- fineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
        susieFit = fit,
        topLoci = .cp_tl(sprintf("chr1:%d:A:G", 100L * (1:3)))
    )
    expect_warning(
        out <- pecotmr:::.colocExtractLbfFromEntry(
            e,
            lbfFilterArgs = ColocLbfFilterParam(
                filterLbfCs = FALSE,
                secondary = NULL,
                concentration = 0.5
            )
        ),
        "LBF matrix is empty"
    )
    expect_null(out)
})

test_that(".colocExtractLbfFromEntry: secondary CS filter subsets rows", {
    fit <- list(
        lbf_variable = matrix(
            1:12,
            3,
            4,
            dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * (1:4)))
        ),
        alpha = matrix(0.25, 3, 4)
    )
    e <- fineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        susieFit = fit,
        topLoci = .cp_tl(sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    local_mocked_bindings(
        .colocFilterCsByConcentration = function(fit, coverage, concentration) {
            2L
        },
        .package = "pecotmr"
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        lbfFilterArgs = ColocLbfFilterParam(
            filterLbfCs = FALSE,
            secondary = 0.95,
            concentration = 0.5
        )
    )
    expect_equal(nrow(out$lbf), 1L)
    expect_equal(as.numeric(out$lbf), c(2, 5, 8, 11)) # row 2 of matrix(1:12, 3, 4)
})

test_that(".colocExtractLbfFromEntry: assigns colnames from variantIds when fit lacks them", {
    fit <- list(lbf_variable = matrix(1:8, 2, 4)) # no dimnames
    e <- fineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        susieFit = fit,
        topLoci = .cp_tl(sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        lbfFilterArgs = ColocLbfFilterParam(
            filterLbfCs = FALSE,
            secondary = NULL,
            concentration = 0.5
        )
    )
    expect_equal(colnames(out$lbf), sprintf("chr1:%d:A:G", 100L * (1:4)))
})

test_that(".colocExtractLbfFromEntry: NULL when every variant column is NA-named", {
    fit <- list(
        lbf_variable = matrix(
            1:8,
            2,
            4,
            dimnames = list(NULL, rep(NA_character_, 4))
        )
    )
    # variantIds length (1) != ncol (4) so colnames are NOT back-filled, and the
    # NA-named columns are dropped, leaving a 0-column matrix.
    e <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = fit,
        topLoci = .cp_tl("chr1:100:A:G")
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        lbfFilterArgs = ColocLbfFilterParam(
            filterLbfCs = FALSE,
            secondary = NULL,
            concentration = 0.5
        )
    )
    expect_null(out)
})

# --- .colocAlignLbf no-overlap branch (~479) --------------------------------

test_that(".colocAlignLbf: returns NULL when there are no shared variants", {
    q <- matrix(
        0,
        2,
        3,
        dimnames = list(NULL, paste0("chr1:", 100 * (1:3), ":A:G"))
    )
    g <- matrix(
        0,
        2,
        3,
        dimnames = list(NULL, paste0("chr2:", 100 * (1:3), ":A:G"))
    )
    expect_null(pecotmr:::.colocAlignLbf(q, g))
})

test_that(".colocEffectRetainedMass reports NA when no reconciliation ran", {
    fit <- list(alpha = matrix(0.25, nrow = 2L, ncol = 4L))
    expect_equal(.colocEffectRetainedMass(fit, 2L), c(NA_real_, NA_real_))
    fit$retained_mass <- c(0.8, 0.3)
    expect_equal(.colocEffectRetainedMass(fit, 2L), c(0.8, 0.3))
    # A length mismatch is a broken parallel, not something to recycle.
    expect_equal(.colocEffectRetainedMass(fit, 3L), rep(NA_real_, 3L))
})

test_that(".colocSelectLbfRows returns indices into the unfiltered rows", {
    lbf <- matrix(0, nrow = 4L, ncol = 3L)
    fit <- list(
        sets = list(cs_index = c(2L, 4L)),
        V = c(0.5, 0, 0.2, 0)
    )
    expect_equal(
        .colocSelectLbfRows(
            lbf,
            fit,
            ColocLbfFilterParam(filterLbfCs = TRUE)
        ),
        c(2L, 4L)
    )
    expect_equal(
        .colocSelectLbfRows(
            lbf,
            fit,
            ColocLbfFilterParam(filterLbfCs = FALSE)
        ),
        c(1L, 3L)
    )
})

test_that(".colocSelectLbfRows does not prefix-match sets_secondary", {
    # `$` would resolve fit$sets to sets_secondary here; `[[` must not.
    lbf <- matrix(0, nrow = 3L, ncol = 2L)
    fit <- list(sets_secondary = list(cs_index = 1L), V = c(1, 1, 0))
    expect_equal(
        .colocSelectLbfRows(
            lbf,
            fit,
            ColocLbfFilterParam(filterLbfCs = TRUE)
        ),
        seq_len(3L)
    )
})

test_that(".colocPickAt indexes per effect and tolerates a missing index", {
    mass <- c(0.9, 0.4, 0.1)
    expect_equal(.colocPickAt(mass, c(3L, 1L), 2L), c(0.1, 0.9))
    expect_equal(.colocPickAt(mass, NULL, 2L), rep(NA_real_, 2L))
    expect_equal(.colocPickAt(mass, c(1L, 99L), 2L), c(0.9, NA_real_))
    expect_equal(.colocPickAt(NULL, 1L, 1L), NA_real_)
})

test_that("colocPipeline reports the retained mass of the scored effects", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr()
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gfmr,
        adjustPips = TRUE
    ))
    skip_if(nrow(out) == 0L, "fixture produced no coloc pairs")
    expect_true(is_in("qtlRetainedMass", colnames(out)))
    expect_true(all(
        is.na(out$qtlRetainedMass) |
            (out$qtlRetainedMass >= 0 & out$qtlRetainedMass <= 1)
    ))
    expect_true(all(
        is.na(out$gwasRetainedMass) |
            (out$gwasRetainedMass >= 0 & out$gwasRetainedMass <= 1)
    ))
})

test_that("colocPipeline publishes coloc's variant count as nSnps", {
    # coloc.bf_bf returns `nsnps`; without the rename .colocStandardiseRow()
    # invents an all-NA `nSnps` beside it and the real count never surfaces.
    sm <- data.frame(nsnps = 12L, PP.H4.abf = 0.5)
    expect_equal(.colocRenameNsnps(sm)$nSnps, 12L)
    expect_false(is_in("nsnps", colnames(.colocRenameNsnps(sm))))
    already <- data.frame(nSnps = 7L)
    expect_equal(.colocRenameNsnps(already)$nSnps, 7L)
})


test_that("gwasFineMapping is refused when the GWAS is already fine-mapped", {
    # The bundle configures the fine-mapping run colocPipeline does on a raw
    # GwasSumStats. Handed an already-fine-mapped GWAS there is no such run,
    # so a setting given here would be dropped rather than honoured.
    data(qtlFineMappingExample, gwasFineMappingExample)
    expect_error(
        colocPipeline(
            qtlFineMappingExample,
            gwasFineMappingExample,
            gwasFineMappingArgs = GwasFineMappingParam(methods = "susieInf")
        ),
        "there is no run to configure"
    )
    # The default bundle is not a request, so it passes silently.
    expect_no_error(pecotmr:::.colocResolveGwasFmr(
        gwasFineMappingExample,
        GwasFineMappingParam()
    ))
    expect_no_error(pecotmr:::.colocResolveGwasFmr(
        gwasFineMappingExample,
        list()
    ))
})

test_that("an explicit blockId is preferred over the derived range key", {
    # blockId keys the external block manifest, so it carries the true block
    # BOUNDARIES rather than the span of whichever variants survived QC.
    data(gwasFineMappingExample)
    expect_equal(
        pecotmr:::.colocGwasBlockIds(gwasFineMappingExample),
        "region_1"
    )
    # Without the column the identity falls back to the variant span.
    g <- gwasFineMappingExample
    S4Vectors::mcols(g)$blockId <- NULL
    expect_equal(
        pecotmr:::.colocGwasBlockIds(g),
        pecotmr:::.rtlRangeKeys(g)
    )
    expect_match(pecotmr:::.colocGwasBlockIds(g), "^chr22_")
})

test_that("an enrichment axis neither side has matches on NA", {
    # `==` evaluates to NA against an absent axis, which would drop every row;
    # the NA-wanted case has to match NA values explicitly instead.
    f <- pecotmr:::.colocEnrichmentColumnMatches
    enr <- data.frame(a = c("x", NA))
    expect_equal(f("a", enr, list(a = NA_character_)), c(FALSE, TRUE))
    expect_equal(f("a", enr, list(a = "x")), c(TRUE, FALSE))
})


test_that("PIP adjustment is skipped when either side has no rows", {
    # Intersecting variants across an empty side would empty the other, so
    # the inputs are passed through untouched instead.
    data(qtlFineMappingExample, gwasFineMappingExample)
    # It now returns just the two (possibly adjusted) inputs, not a bundle.
    emptyQtl <- list(
        qtlFineMappingResult = qtlFineMappingExample[0],
        gwasFmr = gwasFineMappingExample
    )
    expect_identical(
        pecotmr:::.colocMaybeAdjustPips(
            adjustPips = TRUE,
            qtlFineMappingResult = emptyQtl$qtlFineMappingResult,
            gwasFmr = emptyQtl$gwasFmr
        ),
        emptyQtl
    )
    # ...and it is skipped outright when not requested.
    notAsked <- list(
        qtlFineMappingResult = qtlFineMappingExample,
        gwasFmr = gwasFineMappingExample
    )
    expect_identical(
        pecotmr:::.colocMaybeAdjustPips(
            adjustPips = FALSE,
            qtlFineMappingResult = notAsked$qtlFineMappingResult,
            gwasFmr = notAsked$gwasFmr
        ),
        notAsked
    )
})

test_that("pre-extracting LBF from an empty GWAS result yields no blocks", {
    data(gwasFineMappingExample)
    expect_equal(
        pecotmr:::.colocPreextractGwasLbf(
            gwasFineMappingExample[0],
            ColocLbfFilterParam()
        ),
        list()
    )
})

test_that("ColocPriorParam carries the priors and rejects a typo", {
    pr <- ColocPriorParam(p12 = 1e-5)
    expect_s4_class(pr, "ColocPriorParam")
    expect_equal(pr$p12, 1e-5)
    expect_equal(pr$p1, 1e-4)
    expect_error(ColocPriorParam(p13 = 1e-5), "unused argument")
})

test_that("ColocLbfFilterParam carries the filter settings", {
    lf <- ColocLbfFilterParam(filterLbfCs = TRUE, concentration = 0.25)
    expect_true(lf$filterLbfCs)
    expect_equal(lf$concentration, 0.25)
    # NULL secondary means "primary sets only" and is not forwarded.
    expect_false(is_in("secondary", names(ColocLbfFilterParam())))
    expect_error(ColocLbfFilterParam(concentrations = 0.25), "unused argument")
})

test_that("ColocOptions is checked against coloc.bf_bf", {
    skip_if_not_installed("coloc")
    expect_s4_class(ColocOptions(), "MethodOptions")
    expect_output(show(ColocOptions()), "checked against coloc::coloc.bf_bf")
    expect_error(ColocOptions(overlapMin = 0.5), "unknown argument")
})

test_that("priors cannot be set twice", {
    skip_if_not_installed("coloc")
    # p1/p2/p12 are read by the enrichment adjustment as well as forwarded,
    # so allowing them in methodArgs would let the two disagree.
    expect_error(
        pecotmr:::.colocEngineArgs(ColocOptions(p12 = 1e-5), ColocPriorParam()),
        "set through `priors`"
    )
    expect_silent(
        pecotmr:::.colocEngineArgs(ColocOptions(), ColocPriorParam())
    )
})

test_that("colocPipeline refuses bare lists where a constructor is due", {
    expect_error(
        pecotmr:::.colocAssertGroups(
            list(p1 = 1e-4),
            ColocLbfFilterParam(),
            ColocOptions()
        ),
        "must be built with ColocPriorParam\\(\\)"
    )
    expect_error(
        pecotmr:::.colocAssertGroups(
            ColocPriorParam(),
            list(filterLbfCs = TRUE),
            ColocOptions()
        ),
        "must be built with ColocLbfFilterParam\\(\\)"
    )
})

test_that("priors and enrichment cannot both be given", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr()
    enr <- data.frame(
        gwasStudy = "G1",
        qtlStudy = "Q1",
        qtlContext = "c1",
        enrichment = 2.0,
        colocP1 = 3e-4,
        colocP2 = 7e-4,
        colocP12 = 2.5e-5,
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    # fastenloc skips enrichment silently in this situation; pecotmr refuses,
    # because scoring with one set of priors while enriching with another is
    # never what the caller meant.
    expect_error(
        colocPipeline(
            qtlFineMappingResult = qfmr,
            gwasInput = gfmr,
            priors = ColocPriorParam(p12 = 1e-5),
            enrichment = enr
        ),
        "cannot both be given"
    )
    # Either alone is fine.
    expect_no_error(suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gfmr,
        priors = ColocPriorParam(p12 = 1e-5)
    )))
    expect_no_error(suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gfmr,
        enrichment = enr
    )))
})

test_that("GwasFineMappingParam carries only what the inline fit can use", {
    g <- GwasFineMappingParam()
    expect_s4_class(g, "GwasFineMappingParam")
    expect_setequal(
        names(g),
        c(
            "methods",
            "credibleSetArgs",
            "rssArgs",
            "panelFilterArgs",
            "addSusieInf",
            "fitRetention"
        )
    )
    # Nested bundles survive: the inline fit gets a real credibleSetArgs.
    expect_s4_class(
        GwasFineMappingParam(
            credibleSetArgs = CredibleSetParam(coverage = 0.9)
        )$credibleSetArgs,
        "CredibleSetParam"
    )
    expect_equal(
        GwasFineMappingParam(
            credibleSetArgs = CredibleSetParam(coverage = 0.9)
        )$credibleSetArgs$coverage,
        0.9
    )
    # The fineMappingPipeline settings that cannot apply are absent, not
    # carried-and-ignored: it IS the fine-mapping call, so there is no
    # resume cache; CV is refused on sumstats; residualization has nothing
    # to regress out.
    for (nm in c("fineMappingResult", "crossValidation", "residualization")) {
        expect_false(nm %in% names(formals(GwasFineMappingParam)), label = nm)
    }
    # Every field it does carry is a real fineMappingPipeline(GwasSumStats)
    # argument.
    accepted <- names(formals(
        getMethod(
            "fineMappingPipeline",
            "GwasSumStats"
        )@.Data
    ))
    expect_true(all(is_in(
        setdiff(names(formals(GwasFineMappingParam)), "methods"),
        c(
            accepted,
            "credibleSetArgs",
            "rssArgs",
            "panelFilterArgs",
            "addSusieInf",
            "fitRetention"
        )
    )))
    expect_error(GwasFineMappingParam(fitRetention = "none"), "must be one of")
})

test_that("GwasFineMappingParam accessors round-trip a nested bundle", {
    g <- GwasFineMappingParam()
    expect_equal(getFineMappingMethods(g), "susie")
    expect_s4_class(getCredibleSetArgs(g), "CredibleSetParam")

    # Why this class has accessors at all: once CredibleSetParam is
    # settable, the container must let the edited record back in, or the
    # nested settings would be readable but not writable.
    g2 <- setCredibleSetArgs(g, setCoverage(getCredibleSetArgs(g), 0.8))
    expect_equal(getCoverage(getCredibleSetArgs(g2)), 0.8)
    expect_equal(getCoverage(getCredibleSetArgs(g)), 0.95)

    # The nested slots are typed, so the wrong Param is refused rather than
    # surfacing as a missing field somewhere downstream.
    expect_error(
        setRssArgs(g, CredibleSetParam()),
        "not valid for @.+rssArgs"
    )
    expect_s4_class(
        setRssArgs(g, SusieRssParam(serFallback = TRUE)),
        "GwasFineMappingParam"
    )
})

test_that("GwasFineMappingParam accessors read and replace every field", {
    # Dispatching each accessor is what exercises the generics declared in
    # AllGenerics.R; the constructor alone never reaches them.
    p <- GwasFineMappingParam()
    expect_equal(getFineMappingMethods(p), "susie")
    expect_s4_class(getRssArgs(p), "SusieRssParam")
    expect_s4_class(getPanelFilterArgs(p), "PanelFilterParam")
    expect_s4_class(getCredibleSetArgs(p), "CredibleSetParam")
    expect_true(getAddSusieInf(p))
    expect_equal(getFitRetention(p), "slim")

    # Each setter returns a new record and leaves the original alone.
    expect_equal(
        getFineMappingMethods(setFineMappingMethods(p, "susieInf")),
        "susieInf"
    )
    expect_false(getAddSusieInf(setAddSusieInf(p, FALSE)))
    expect_equal(getFitRetention(setFitRetention(p, "full")), "full")
    expect_s4_class(
        getPanelFilterArgs(setPanelFilterArgs(p, PanelFilterParam())),
        "PanelFilterParam"
    )
    expect_equal(getFineMappingMethods(p), "susie")
    expect_true(getAddSusieInf(p))
})
