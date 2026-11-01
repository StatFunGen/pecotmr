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
        path = "/tmp/sketch",
        format = "plink2",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * seq_len(snp_n)),
            CHR = rep("chr1", snp_n),
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
    FineMappingRow(variantIds = variant_ids, susieFit = fit, topLoci = tl)
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
        studyName = map_chr(tuples, 1L),
        context = map_chr(tuples, 2L),
        trait = map_chr(tuples, 3L),
        method = map_chr(tuples, 4L),
        entry = entries,
        ldReference = if (with_sketch) .cp_makeHandle() else NULL
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
        studyName = map_chr(tuples, 1L),
        method = map_chr(tuples, 2L),
        entry = entries,
        ldReference = if (with_sketch) .cp_makeHandle() else NULL
    )
}

.cp_makeGwasSumstats <- function(studyName = "G1", qc = TRUE) {
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
        studyName = studyName,
        entry = list(gr),
        genome = "hg19",
        ldReference = .cp_makeHandle(),
        qcOptions = if (qc) list(step1 = "ok") else list()
    )
}

.cp_makeQtlSumstats <- function(studyName = "Q1", qc = TRUE) {
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
        studyName = studyName,
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .cp_makeHandle(),
        qcOptions = if (qc) list(step1 = "ok") else list()
    )
}

# The enrichment lookup takes each side's identity list, as the scoring loop
# builds it.
.cp_side <- function(
    studyName,
    context = NA_character_,
    trait = NA_character_
) {
    list(study = studyName, context = context, trait = trait)
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
# .colocRequireMatchingLdRefs
# ===========================================================================

test_that(".colocRequireMatchingLdRefs: NULL qtl-side ldReference is allowed", {
    qfmr <- .cp_makeQtlFmr(with_sketch = FALSE)
    gfmr <- .cp_makeGwasFmr()
    local_mocked_bindings(coloc.bf_bf = .cp_mockColocBfBf(), .package = "coloc")
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gfmr
    ))
    expect_s4_class(out, "ColocResult")
})

test_that(".colocRequireMatchingLdRefs: non-NULL qtl + NULL gwas errors", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- .cp_makeGwasFmr(with_sketch = FALSE)
    expect_error(
        colocPipeline(qtlFineMappingResult = qfmr, gwasInput = gfmr),
        "ldReference is NULL"
    )
})

test_that(".colocRequireMatchingLdRefs: sample set mismatch errors", {
    qfmr <- .cp_makeQtlFmr()
    otherPanel <- .cp_makeHandle(sample_prefix = "other")
    gfmr <- GwasFineMappingResult(
        studyName = "G1",
        method = "susie",
        entry = list(.cp_makeFmEntry()),
        ldReference = otherPanel
    )
    expect_error(
        colocPipeline(qtlFineMappingResult = qfmr, gwasInput = gfmr),
        "different sample sets"
    )
})

test_that(".colocRequireMatchingLdRefs: differently trimmed panels pass", {
    # The two sides QC'd separately keep different subsets of one LD
    # reference; the panels overlap, so the check warns rather than aborting.
    qfmr <- pecotmr:::.withRowLdRefs(
        .cp_makeQtlFmr(),
        .cp_makeHandle(snp_n = 6L)
    )
    gfmr <- pecotmr:::.withRowLdRefs(
        .cp_makeGwasFmr(),
        .cp_makeHandle(snp_n = 7L)
    )
    expect_warning(
        expect_null(pecotmr:::.colocRequireMatchingLdRefs(qfmr, gfmr)),
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

test_that("colocPipeline carries the sides' genome build and refuses two", {
    qfmr <- .cp_makeQtlFmr()
    gfmr <- GenomeInfoDb::`genome<-`(.cp_makeGwasFmr(), value = "hg38")
    # Per-variant results too, so the pair holds variants for seqinfo to
    # describe: a collection of empty elements has no sequence to name a
    # build for.
    summaryOnly <- .cp_mockColocBfBf()
    local_mocked_bindings(
        coloc.bf_bf = function(qLbf, gLbf, ...) {
            res <- summaryOnly(qLbf, gLbf, ...)
            # The LBF matrices carry a `null` column beside the variants.
            ids <- colnames(qLbf)[str_detect(colnames(qLbf), ":")]
            res$results <- data.frame(
                snp = ids,
                SNP.PP.H4.row1 = 1 / length(ids)
            )
            res
        },
        .package = "coloc"
    )
    out <- suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gfmr
    ))
    expect_gt(sum(lengths(out)), 0L)
    expect_identical(unique(unname(GenomeInfoDb::genome(out))), "hg38")
    expect_identical(
        unique(unname(GenomeInfoDb::genome(colocVariants(out)))),
        "hg38"
    )
    hg19 <- GenomeInfoDb::`genome<-`(qfmr, value = "hg19")
    expect_error(
        colocPipeline(qtlFineMappingResult = hg19, gwasInput = gfmr),
        "colocPipeline: the inputs name different genome builds"
    )
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
    # Like the real pipeline, the stand-in keeps each row's block name, which
    # is what the pairing finds a region's fit by.
    local_mocked_bindings(
        fineMappingPipeline = function(data, methods, ...) {
            fmr <- .cp_makeGwasFmr()
            fmr$blockId <- as.character(data$blockId)
            fmr
        },
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
    e <- FineMappingRow(
        variantIds = shared,
        susieFit = emptyFit,
        topLoci = data.frame(
            variant_id = shared,
            pip = 1,
            stringsAsFactors = FALSE
        )
    )
    gfmr <- GwasFineMappingResult(
        studyName = "G1",
        method = "susie",
        entry = list(e),
        ldReference = .cp_makeHandle()
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
            "gwasRetainedMass",
            "ldReference"
        )
    )
})

# ===========================================================================
# Internal helpers
# ===========================================================================

test_that(".colocExtractLbfFromEntry: refuses an entry with no retained fit", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = NULL,
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.1,
            stringsAsFactors = FALSE
        )
    )
    # An absent fit means an absent LBF, and coloc.bf_bf has nothing to score.
    # That is a refusal, not a skip: warning per tuple produced an empty
    # result that read as "no colocalization" rather than "this input cannot
    # be colocalized".
    expect_error(
        pecotmr:::.colocExtractLbfFromEntry(e, ColocLbfFilterParam()),
        "no retained fit"
    )
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
    e <- FineMappingRow(
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
        colocLbfFilterParam = ColocLbfFilterParam(
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

# ===========================================================================
# .colocLookupEnrichment
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
# .colocEmptyResult(enriched = TRUE)
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

test_that("colocPipeline: enrichment takes table priors and reports RCP/LCP", {
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
    # RCP and LCP are reported alongside the hypotheses.
    pairs <- colocPairs(out)
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
    pairs <- colocPairs(out)
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

# --- returnGwasFineMapping attached on the empty-pair early return ----------

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
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = emptyFit,
        topLoci = .cp_tl("chr1:100:A:G", pip = 0)
    )
    resolved <- GwasFineMappingResult(
        studyName = "G1",
        method = "susie",
        entry = list(e),
        ldReference = .cp_makeHandle()
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

# --- per-pair loop skips: no usable QTL LBF and no overlap -------------------

test_that("colocPipeline: skips a QTL entry with no usable LBF", {
    badFit <- list(
        alpha = matrix(0, 1, 1),
        pip = c(v1 = 0),
        V = 0,
        lbf_variable = matrix(NA_real_, 1, 1)
    )
    badEntry <- FineMappingRow(
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

# --- .colocFilterCsByConcentration -----------------------------------------

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
    e <- FineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        susieFit = fit,
        topLoci = .cp_tl(sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        colocLbfFilterParam = ColocLbfFilterParam(
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
    e <- FineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        susieFit = fit,
        topLoci = .cp_tl(sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        colocLbfFilterParam = ColocLbfFilterParam(
            filterLbfCs = FALSE,
            secondary = NULL,
            concentration = 0.5
        )
    )
    expect_equal(nrow(out$lbf), 2L)
})

test_that(".colocExtractLbfFromEntry: refuses a fit carrying no LBF slot", {
    fit <- list(notLbf = list(a = 1))
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = fit,
        topLoci = .cp_tl("chr1:100:A:G")
    )
    expect_error(
        pecotmr:::.colocExtractLbfFromEntry(
            e,
            colocLbfFilterParam = ColocLbfFilterParam(
                filterLbfCs = FALSE,
                secondary = NULL,
                concentration = 0.5
            )
        ),
        "no per-effect LBF matrix"
    )
})

test_that(".colocExtractLbfFromEntry: reads a bare fSuSiE lBF list", {
    # The shape an fSuSiE fit actually carries. coloc used to look only for
    # lbf_variable / fsusie_result$lBF and skipped these silently, while
    # lbf() read them fine -- the two extractors had drifted apart.
    data(fsusieFineMappingExample)
    fit <- susieFit(fsusieFineMappingExample)
    expect_type(fit$lBF, "list")
    m <- pecotmr:::.colocLbfMatrix(fit, "probe")
    expect_true(is.matrix(m))
    expect_identical(nrow(m), length(fit$lBF))
    expect_identical(m, pecotmr:::.asLbfMatrix(fit))
})

test_that(".colocExtractLbfFromEntry: refuses an empty LBF matrix", {
    fit <- list(
        lbf_variable = matrix(
            numeric(0),
            nrow = 0,
            ncol = 3,
            dimnames = list(NULL, sprintf("chr1:%d:A:G", 100L * (1:3)))
        )
    )
    e <- FineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
        susieFit = fit,
        topLoci = .cp_tl(sprintf("chr1:%d:A:G", 100L * (1:3)))
    )
    expect_error(
        pecotmr:::.colocExtractLbfFromEntry(
            e,
            colocLbfFilterParam = ColocLbfFilterParam(
                filterLbfCs = FALSE,
                secondary = NULL,
                concentration = 0.5
            )
        ),
        "no per-effect LBF matrix"
    )
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
    e <- FineMappingRow(
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
        colocLbfFilterParam = ColocLbfFilterParam(
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
    e <- FineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:4)),
        susieFit = fit,
        topLoci = .cp_tl(sprintf("chr1:%d:A:G", 100L * (1:4)))
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        colocLbfFilterParam = ColocLbfFilterParam(
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
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = fit,
        topLoci = .cp_tl("chr1:100:A:G")
    )
    out <- pecotmr:::.colocExtractLbfFromEntry(
        e,
        colocLbfFilterParam = ColocLbfFilterParam(
            filterLbfCs = FALSE,
            secondary = NULL,
            concentration = 0.5
        )
    )
    expect_null(out)
})

# --- .colocAlignLbf no-overlap branch --------------------------------------

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
            gwasFineMappingParam = GwasFineMappingParam(methods = "susieInf")
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
    # Refused at construction now -- ColocOptions() owns the check, so the
    # deeper .colocEngineArgs() guard is never reached for this input.
    expect_error(ColocOptions(p12 = 1e-5), "supplied by pecotmr")
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
    # Either alone is fine; the same table alone is the enrichment-mode run
    # tested above.
    expect_no_error(suppressWarnings(colocPipeline(
        qtlFineMappingResult = qfmr,
        gwasInput = gfmr,
        priors = ColocPriorParam(p12 = 1e-5)
    )))
})

test_that("GwasFineMappingParam carries only what the inline fit can use", {
    g <- GwasFineMappingParam()
    expect_s4_class(g, "GwasFineMappingParam")
    expect_setequal(
        names(g),
        c(
            "methods",
            "credibleSetParam",
            "susieRssParam",
            "panelFilterParam",
            "initializeWithSusieInf",
            "fitRetention"
        )
    )
    # Nested bundles survive: the inline fit gets a real Param back.
    expect_s4_class(
        GwasFineMappingParam(
            credibleSetParam = CredibleSetParam(coverage = 0.9)
        )$credibleSetParam,
        "CredibleSetParam"
    )
    expect_equal(
        GwasFineMappingParam(
            credibleSetParam = CredibleSetParam(coverage = 0.9)
        )$credibleSetParam$coverage,
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
            # Listed explicitly because the S4 method's visible formals
            # are only `(data, ...)` -- the per-input-kind arguments are
            # documented on the generic -- so there is nothing in `accepted`
            # to match them against. The names themselves now agree with the
            # pipeline's own; the earlier split that spelled these
            # `credibleSetArgs` / `rssArgs` / `panelFilterArgs` /
            # `addSusieInf` on pipelines was retired once the
            # accessor-shadowing conflict was measured not to arise.
            "credibleSetParam",
            "susieRssParam",
            "panelFilterParam",
            "initializeWithSusieInf",
            "fitRetention"
        )
    )))
    expect_error(GwasFineMappingParam(fitRetention = "none"), "must be one of")
})

test_that("GwasFineMappingParam accessors round-trip a nested bundle", {
    g <- GwasFineMappingParam()
    expect_equal(fineMappingMethods(g), "susie")
    expect_s4_class(credibleSetParam(g), "CredibleSetParam")

    # The container must let an edited nested record back in, or the nested
    # settings would be readable but not writable.
    g2 <- g
    cs <- g2$credibleSetParam
    cs$coverage <- 0.8
    g2$credibleSetParam <- cs
    expect_equal(csCoverage(credibleSetParam(g2)), 0.8)
    expect_equal(csCoverage(credibleSetParam(g)), 0.95)

    # The nested slots are typed, so the wrong Param is refused rather than
    # surfacing as a missing field somewhere downstream.
    expect_error(
        g$susieRssParam <- CredibleSetParam(),
        "susieRssParam"
    )
    g3 <- g
    g3$susieRssParam <- SusieRssParam(serFallback = TRUE)
    expect_s4_class(g3, "GwasFineMappingParam")
})

test_that("GwasFineMappingParam fields read and replace, each one", {
    # Dispatching each accessor is what exercises the generics declared in
    # AllGenerics.R; the constructor alone never reaches them.
    p <- GwasFineMappingParam()
    expect_equal(fineMappingMethods(p), "susie")
    expect_s4_class(susieRssParam(p), "SusieRssParam")
    expect_s4_class(panelFilterParam(p), "PanelFilterParam")
    expect_s4_class(credibleSetParam(p), "CredibleSetParam")
    expect_true(initializeWithSusieInf(p))
    expect_equal(fitRetention(p), "slim")

    # Replacement edits a copy and leaves the original alone.
    edited <- p
    edited$methods <- "susieInf"
    edited$initializeWithSusieInf <- FALSE
    edited$fitRetention <- "full"
    edited$panelFilterParam <- PanelFilterParam()
    expect_equal(fineMappingMethods(edited), "susieInf")
    expect_false(initializeWithSusieInf(edited))
    expect_equal(fitRetention(edited), "full")
    expect_s4_class(panelFilterParam(edited), "PanelFilterParam")
    expect_equal(fineMappingMethods(p), "susie")
    expect_true(initializeWithSusieInf(p))
})

# ===========================================================================
# Per-pair LD references
# ===========================================================================

# An LdMatrix over five variants at 100 bp spacing on `chrom`.
.cp_ldm <- function(rho, chrom = "chr1") {
    ids <- sprintf("%s:%d:A:G", chrom, 100L * (1:5))
    LdMatrix(`dimnames<-`(rho^abs(outer(1:5, 1:5, "-")), list(ids, ids)))
}

test_that("coloc refuses an overlapping pair on two different references", {
    qfmr <- pecotmr:::.withRowLdRefs(.cp_makeQtlFmr(), .cp_ldm(0.2))
    gfmr <- pecotmr:::.withRowLdRefs(.cp_makeGwasFmr(), .cp_ldm(0.8))
    expect_error(
        pecotmr:::.colocRequireMatchingLdRefs(qfmr, gfmr),
        "disagree on the LD"
    )
})

test_that("coloc lets rows that never pair sit on different references", {
    # The QTL row is on chr1 and the GWAS row on chr2, so coloc never scores
    # them against each other and their references need not match.
    qfmr <- pecotmr:::.withRowLdRefs(.cp_makeQtlFmr(), .cp_ldm(0.2))
    chr2 <- .cp_makeFmEntry(variant_ids = paste0("chr2:", 100 * (1:5), ":A:G"))
    gfmr <- pecotmr:::.withRowLdRefs(
        .cp_makeGwasFmr(entries = list(chr2)),
        .cp_ldm(0.8, chrom = "chr2")
    )
    expect_null(pecotmr:::.colocRequireMatchingLdRefs(qfmr, gfmr))
})

# ===========================================================================
# A QTL region over several GWAS LD blocks is one colocalization
# ===========================================================================

# One GWAS study fine-mapped in two LD blocks: chr1:100-500 in block
# chr1_1_550 and chr1:600-1000 in block chr1_551_2000, two effects each.
.cp_blockIds <- c("chr1_1_550", "chr1_551_2000")
.cp_blockWindows <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(c(1L, 551L), c(550L, 2000L))
)

.cp_makeBlockedGwasFmr <- function(
    blockId = .cp_blockIds,
    method = "susie",
    region = .cp_blockWindows
) {
    GwasFineMappingResult(
        studyName = rep("G1", 2L),
        method = rep(method, 2L),
        entry = list(
            .cp_makeFmEntry(paste0("chr1:", 100 * (1:5), ":A:G")),
            .cp_makeFmEntry(paste0("chr1:", 100 * (6:10), ":A:G"))
        ),
        blockId = blockId,
        region = region,
        ldReference = .cp_makeHandle(snp_n = 10L)
    )
}

# A QTL fit over `positions` (multiples of 100 on chr1).
.cp_makeQtlOver <- function(positions, trait = "t1") {
    QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = trait,
        method = "susie",
        entry = list(.cp_makeFmEntry(paste0("chr1:", positions, ":A:G"))),
        ldReference = .cp_makeHandle(snp_n = 10L)
    )
}

test_that(".colocCollapseBlocks stacks blocks, each confined to its own", {
    recs <- pecotmr:::.colocPreextractGwasLbf(
        .cp_makeBlockedGwasFmr(),
        ColocLbfFilterParam()
    )
    expect_length(recs, 2L)
    one <- pecotmr:::.colocCollapseBlocks(recs[1L])
    expect_identical(one, recs[[1L]])
    merged <- pecotmr:::.colocCollapseBlocks(recs)
    blockA <- paste0("chr1:", 100 * (1:5), ":A:G")
    blockB <- paste0("chr1:", 100 * (6:10), ":A:G")
    expect_identical(colnames(merged$lbf), c(blockA, blockB))
    # Block A's effects carry no evidence about block B's variants, and the
    # reverse.
    expect_true(all(merged$lbf[1:2, blockB] == -Inf))
    expect_true(all(merged$lbf[3:4, blockA] == -Inf))
    expect_equal(merged$lbf[1:2, blockA], recs[[1L]]$lbf)
    # Block B's effects are numbered on from block A's two.
    expect_identical(merged$effect, 1:4)
    expect_identical(merged$nEffects, 4L)
    # The blockIds spell coordinates, so the merged block is their union.
    expect_identical(merged$blockId, "chr1_1_2000")
    expect_match(merged$label, "chr1_1_2000")
})

test_that(".colocMergedBlockId names the stretch a fit's blocks cover", {
    blocks <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 551L), c(550L, 2000L))
    )
    spans <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(100L, 600L), c(500L, 1000L))
    )
    # The stretch comes from the blocks' windows, whatever the ids say.
    expect_identical(
        pecotmr:::.colocMergedBlockId(.cp_blockIds, blocks),
        "chr1_1_2000"
    )
    expect_identical(
        pecotmr:::.colocMergedBlockId(c("blockA", "blockB"), spans),
        "chr1_100_1000"
    )
    # One block keeps its own id, whatever it is.
    expect_identical(
        pecotmr:::.colocMergedBlockId("blockA", blocks[1L]),
        "blockA"
    )
    expect_identical(pecotmr:::.colocMergedBlockId("chr1", spans[1L]), "chr1")
})

test_that(".colocMergedLdRef keeps a shared reference and drops a split one", {
    ref <- .cp_makeHandle()
    same <- list(
        list(ldKey = "ld1", ldReference = ref),
        list(ldKey = "ld1", ldReference = ref)
    )
    expect_identical(pecotmr:::.colocMergedLdRef(same), ref)
    split <- list(
        list(ldKey = "ld1", ldReference = ref),
        list(ldKey = "ld2", ldReference = ref)
    )
    expect_null(pecotmr:::.colocMergedLdRef(split))
})

test_that("colocPipeline scores a region over two blocks as one fit", {
    skip_if_not_installed("coloc")
    set.seed(3)
    res <- suppressWarnings(colocPipeline(
        .cp_makeQtlOver(100 * (1:10)),
        gwasInput = .cp_makeBlockedGwasFmr(),
        adjustPips = FALSE
    ))
    pairs <- colocPairs(res)
    expect_identical(unique(pairs$blockId), "chr1_1_2000")
    # Every QTL effect is scored over its whole region, against all four GWAS
    # effects of the one collapsed fit.
    expect_true(all(pairs$nSnps == 10L))
    expect_setequal(pairs$gwasCs, 1:4)
    pp <- as.matrix(pairs[, str_c("PP.H", 0:4, ".abf")])
    expect_false(anyNA(pp))
    expect_equal(unname(rowSums(pp)), rep(1, nrow(pairs)))
    # A region inside one block reads that block alone, under its own id.
    inside <- suppressWarnings(colocPipeline(
        .cp_makeQtlOver(100 * (1:5)),
        gwasInput = .cp_makeBlockedGwasFmr(),
        adjustPips = FALSE
    ))
    expect_identical(unique(colocPairs(inside)$blockId), "chr1_1_550")
    expect_true(all(colocPairs(inside)$nSnps == 5L))
})

test_that("colocPipeline keeps two methods' block fits apart", {
    skip_if_not_installed("coloc")
    set.seed(4)
    gfmr <- combineFineMappingResults(
        .cp_makeBlockedGwasFmr(method = "susie"),
        .cp_makeBlockedGwasFmr(method = "ser")
    )
    res <- suppressWarnings(colocPipeline(
        .cp_makeQtlOver(100 * (1:10)),
        gwasInput = gfmr,
        adjustPips = FALSE
    ))
    pairs <- colocPairs(res)
    byMethod <- split(pairs$gwasCs, pairs$gwasMethod)
    expect_setequal(names(byMethod), c("susie", "ser"))
    expect_setequal(byMethod$susie, 1:4)
    expect_setequal(byMethod$ser, 1:4)
})

# A GwasSumStats split into the two blocks above, and a fineMappingPipeline
# stand-in that records what it was given and fits each row it was handed.
.cp_makeBlockedGwasSumstats <- function() {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(100L * (1:10), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:10)),
        A1 = "A",
        A2 = "G",
        Z = seq_len(10L) / 10,
        N = 1000L
    )
    GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .cp_makeHandle(snp_n = 10L),
        qcOptions = list(step1 = "ok"),
        region = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(c(1L, 551L), c(550L, 2000L))
        )
    )
}

.cp_recordingFineMapper <- function(seen) {
    function(data, methods, ...) {
        seen$data <- data
        GwasFineMappingResult(
            studyName = as.character(data$study),
            method = rep("susie", nrow(data)),
            entry = map(
                as.list(data),
                function(gr) .cp_makeFmEntry(S4Vectors::mcols(gr)$SNP)
            ),
            blockId = as.character(data$blockId),
            ldReference = .cp_makeHandle(snp_n = 10L)
        )
    }
}

test_that("colocPipeline fine-maps each set of blocks a region reaches once", {
    seen <- new.env()
    scored <- new.env()
    scored$gwas <- character(0)
    local_mocked_bindings(
        fineMappingPipeline = .cp_recordingFineMapper(seen),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        coloc.bf_bf = function(qLbf, gLbf, p1, p2, p12, ...) {
            scored$gwas <- c(scored$gwas, str_flatten(colnames(gLbf), ","))
            .cp_mockColocBfBf()(qLbf, gLbf, p1, p2, p12)
        },
        .package = "coloc"
    )
    # The first region reaches no GWAS row, so the PIP reconciliation drops
    # it, and the regions after it move up a row; each must still find its
    # own fit.
    elsewhere <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "elsewhere",
        method = "susie",
        entry = list(.cp_makeFmEntry(paste0("chr2:", 100 * (1:5), ":A:G")))
    )
    qtl <- combineFineMappingResults(
        elsewhere,
        .cp_makeQtlOver(100 * (1:5), trait = "inside"),
        .cp_makeQtlOver(100 * (4:8), trait = "across")
    )
    res <- suppressMessages(suppressWarnings(colocPipeline(
        qtl,
        gwasInput = .cp_makeBlockedGwasSumstats(),
        returnGwasFineMapping = TRUE
    )))
    # Two sets: block A alone, and both blocks stitched together.
    expect_identical(
        as.character(seen$data$blockId),
        c("chr1_1_550", "chr1_1_2000")
    )
    expect_identical(unname(lengths(seen$data)), c(5L, 10L))
    # Each set's range is the stretch its blocks' given ranges cover.
    expect_identical(
        as.character(seen$data$region),
        c("chr1:1-550", "chr1:1-2000")
    )
    expect_identical(seen$data$regionGiven, c(TRUE, TRUE))
    # Each region is scored against its own set's fit and no other.
    pairs <- colocPairs(res)
    expect_identical(
        pairs$blockId[match(c("inside", "across"), pairs$trait)],
        c("chr1_1_550", "chr1_1_2000")
    )
    expect_length(scored$gwas, 2L)
    expect_identical(
        S4Vectors::mcols(attr(res, "gwasFineMapping"))$blockId,
        c("chr1_1_550", "chr1_1_2000")
    )
})

test_that("colocPipeline needs one LD reference across stitched blocks", {
    gws <- pecotmr:::.withRowLdRefs(
        .cp_makeBlockedGwasSumstats(),
        list(.cp_ldm(0.2), .cp_ldm(0.5))
    )
    expect_error(
        pecotmr:::.colocResolveGwasFmr(
            gws,
            GwasFineMappingParam(),
            .cp_makeQtlOver(100 * (4:8))
        ),
        "the rows of GWAS study 'G1' a QTL region reaches use 2"
    )
})

test_that("colocPipeline fine-maps nothing when no region reaches the GWAS", {
    resolved <- pecotmr:::.colocResolveGwasFmr(
        .cp_makeBlockedGwasSumstats(),
        GwasFineMappingParam(),
        QtlFineMappingResult(
            studyName = "Q1",
            context = "c1",
            trait = "t1",
            method = "susie",
            entry = list(.cp_makeFmEntry(paste0("chr2:", 100 * (1:5), ":A:G")))
        )
    )
    expect_equal(nrow(resolved$fmr), 0L)
    expect_equal(nrow(resolved$pairing), 0L)
})

test_that(".colocWithNull adds coloc's null column at log BF 0", {
    m <- matrix(c(1, -Inf), 1L, dimnames = list(NULL, c("a", "b")))
    expect_identical(
        pecotmr:::.colocWithNull(m),
        cbind(m, null = 0)
    )
})

test_that("enloc reads SuSiE's natural-log Bayes factors as such", {
    # With a flat GWAS effect the colocalization mass follows the QTL PIPs,
    # which must be SuSiE's softmax over natural-log Bayes factors.
    lbf <- c(2, 1, 0)
    aligned <- list(
        qtl = matrix(lbf, 1L, dimnames = list(NULL, c("a", "b", "c"))),
        gwas = matrix(0, 1L, 3L, dimnames = list(NULL, c("a", "b", "c")))
    )
    out <- pecotmr:::.colocEnlocOneEffectPair(
        1L,
        1L,
        aligned,
        list(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
    )
    expect_equal(
        unname(out$scp / sum(out$scp)),
        exp(lbf) / sum(exp(lbf))
    )
})

test_that("the enrichment lookup matches on method when the table has it", {
    enr <- data.frame(
        gwasStudy = "G1",
        gwasMethod = c("susie", "ser"),
        qtlStudy = "Q1",
        qtlContext = "c1",
        qtlMethod = "susie",
        enrichment = c(1, 2)
    )
    expect_null(pecotmr:::.colocValidateEnrichmentKeys(enr))
    q <- list(study = "Q1", context = "c1", trait = "t1", method = "susie")
    gwas <- function(m) {
        list(study = "G1", context = NA, trait = NA, method = m)
    }
    expect_identical(
        pecotmr:::.colocLookupEnrichment(enr, gwas("ser"), q),
        2L
    )
    expect_identical(
        pecotmr:::.colocLookupEnrichment(enr, gwas("susie"), q),
        1L
    )
    expect_true(is.na(
        pecotmr:::.colocLookupEnrichment(enr, gwas("susieInf"), q)
    ))
})

test_that("a stitched set of blocks keeps its case and control counts", {
    gws <- .cp_makeBlockedGwasSumstats()
    withCounts <- GwasSumStats(
        studyName = "G1",
        entry = list(sumStats(gws)),
        genome = "hg19",
        ldReference = .cp_makeHandle(snp_n = 10L),
        nCase = 100,
        nControl = 200,
        qcOptions = list(step1 = "ok"),
        region = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(c(1L, 551L), c(550L, 2000L))
        )
    )
    sets <- pecotmr:::.colocQtlRowSets(
        withCounts,
        .cp_makeQtlOver(100 * (4:8))
    )
    u <- pecotmr:::.colocUnionSumStats(withCounts, sets)
    expect_equal(nrow(u), 1L)
    expect_equal(u$nCase, 100)
    expect_equal(u$nControl, 200)
    expect_false("nSample" %in% colnames(S4Vectors::mcols(u)))
})

test_that("a stitched set is given its range only when every block was", {
    gws <- .cp_makeBlockedGwasSumstats()
    sets <- pecotmr:::.colocQtlRowSets(gws, .cp_makeQtlOver(100 * (4:8)))
    given <- pecotmr:::.colocUnionSumStats(gws, sets)
    expect_identical(as.character(given$region), "chr1:1-2000")
    expect_true(given$regionGiven)
    # The same blocks with their ranges computed: the set's range is the
    # stretch their spans cover, and it is computed too.
    computed <- GwasSumStats(
        studyName = c("G1", "G1"),
        entry = as.list(gws),
        genome = "hg19",
        ldReference = .cp_makeHandle(snp_n = 10L),
        qcOptions = list(step1 = "ok"),
        blockId = c("bA", "bB")
    )
    u <- pecotmr:::.colocUnionSumStats(
        computed,
        pecotmr:::.colocQtlRowSets(computed, .cp_makeQtlOver(100 * (4:8)))
    )
    expect_identical(u$blockId, "chr1_100_1000")
    expect_identical(as.character(u$region), "chr1:100-1000")
    expect_false(u$regionGiven)
})

# ===========================================================================
# A region reaching past the GWAS blocks that were loaded
# ===========================================================================

test_that("colocPipeline refuses a region reaching past given blocks", {
    skip_if_not_installed("coloc")
    # The blocks' ids spell their ranges, so chr1:2100 lies in a block that
    # was not loaded.
    expect_error(
        colocPipeline(
            .cp_makeQtlOver(c(100 * (1:10), 2100)),
            gwasInput = .cp_makeBlockedGwasFmr(),
            adjustPips = FALSE
        ),
        "1 variant(s) past the end of block chr1:551-2000 of GWAS study",
        fixed = TRUE
    )
    # Statistics are checked before any block is fine-mapped.
    expect_error(
        colocPipeline(
            .cp_makeQtlOver(c(100 * (4:8), 2100)),
            gwasInput = .cp_makeBlockedGwasSumstats()
        ),
        "lie in blocks that were not loaded"
    )
})

test_that("colocPipeline warns of a region running past computed blocks", {
    # Labels rather than coordinates: the blocks' ranges are the spans of
    # their variants, chr1:100-500 and chr1:600-1000.
    gfmr <- .cp_makeBlockedGwasFmr(blockId = c("bA", "bB"), region = NULL)
    expect_warning(
        pecotmr:::.colocCheckBlockReach(.cp_makeQtlOver(100 * (1:12)), gfmr),
        "2 variant\\(s\\) past the end of block bB \\(chr1:600-1000\\)"
    )
    # One variant past a computed range may only be one the GWAS lacks.
    expect_no_warning(
        pecotmr:::.colocCheckBlockReach(.cp_makeQtlOver(100 * (1:11)), gfmr)
    )
    # A region meeting no block of the study is not checked against it.
    expect_no_warning(
        pecotmr:::.colocCheckBlockReach(.cp_makeQtlOver(100 * (30:35)), gfmr)
    )
    # Nor is a second side that has no blocks.
    expect_null(
        pecotmr:::.colocCheckBlockReach(
            .cp_makeQtlOver(100 * (1:12)),
            .cp_makeQtlOver(100 * (1:5))
        )
    )
})

# ===========================================================================
# LBF column names come from the fit's own variants
# ===========================================================================

test_that("an unnamed LBF matrix is named by the variants of its fit", {
    vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    tl <- data.frame(variant_id = vids, pip = c(0.5, 0.3, 0.2))
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "g1",
        method = "susie",
        entry = list(FineMappingRow(
            vids,
            list(pip = set_names(tl$pip, vids)),
            tl
        ))
    )
    # The row holds two of the three variants the fit's columns are, found
    # by the names the fit gives them.
    parts <- pecotmr:::.fmrRowParts(endoapply(res, utils::head, 2L), 1L)
    named <- pecotmr:::.colocLbfNamesFrom(matrix(c(1, 2, 3), 1L), parts)
    expect_equal(colnames(named), vids)
    # A matrix that names its own columns keeps them; one matching no
    # variant list stays unnamed.
    own <- matrix(1:3, 1L, dimnames = list(NULL, c("a", "b", "c")))
    expect_identical(pecotmr:::.colocLbfNamesFrom(own, parts), own)
    expect_null(
        colnames(pecotmr:::.colocLbfNamesFrom(matrix(1:2, 1L), parts))
    )
})

# ===========================================================================
# Narrowed rows: a fit scored whole warns
# ===========================================================================

test_that("fits scored as supplied warn about narrowed rows on both sides", {
    data(qtlFineMappingExample, gwasFineMappingExample)
    qtl <- endoapply(qtlFineMappingExample, utils::head, 2000L)
    gwas <- endoapply(gwasFineMappingExample, utils::head, 2000L)
    cnd <- tryCatch(
        pecotmr:::.colocMaybeAdjustPips(
            adjustPips = FALSE,
            qtlFineMappingResult = qtl,
            gwasFmr = gwas
        ),
        warning = identity
    )
    expect_s3_class(cnd, "pecotmrFitCoverage")
    expect_match(conditionMessage(cnd), "colocPipeline\\(\\): 2 row\\(s\\)")
    # An empty side scores nothing, so nothing is read.
    expect_no_warning(pecotmr:::.colocMaybeAdjustPips(
        adjustPips = FALSE,
        qtlFineMappingResult = qtl[0],
        gwasFmr = gwas
    ))
})

test_that("colocPipeline warns about narrowed fits only when not adjusting", {
    skip_if_not_installed("coloc")
    data(qtlFineMappingLbfExample, gwasFineMappingLbfExample)
    narrowed <- endoapply(qtlFineMappingLbfExample, utils::head, 150L)
    asSupplied <- .fitCoverageMessages(colocPipeline(
        narrowed,
        gwasInput = gwasFineMappingLbfExample,
        adjustPips = FALSE
    ))
    expect_length(asSupplied, 1L)
    expect_match(asSupplied, "colocPipeline\\(\\): 2 row\\(s\\)")
    # The adjustment restricts each fit to the variants its row keeps.
    expect_length(
        .fitCoverageMessages(colocPipeline(
            narrowed,
            gwasInput = gwasFineMappingLbfExample
        )),
        0L
    )
})
