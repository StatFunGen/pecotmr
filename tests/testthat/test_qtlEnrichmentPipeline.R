context("qtlEnrichmentPipeline")

# ===========================================================================
# Strategy: mock qtlEnrichment so the pipeline runs end-to-end on a
# small fixture, but the heavy mixture-of-enrichment estimator never fires.
# ===========================================================================

# The internal builders take an identity list keyed by the collection's own
# columns; the axes a GWAS collection lacks are NA.
.qep_gwasIdent <- function(
    studyName,
    context = NA_character_,
    trait = NA_character_
) {
    list(study = studyName, context = context, trait = trait)
}

.qep_qtlIdent <- function(studyName, context) {
    list(study = studyName, context = context)
}

.qep_makeHandle <- function(
    snp_n = 6L,
    n_samples = 30L,
    path = "/tmp/sketch.gds"
) {
    new(
        "GenotypeHandle",
        path = path,
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
        sampleIds = paste0("s", seq_len(n_samples)),
        pgenPtr = NULL
    )
}

.qep_makeFmEntry <- function(
    variant_ids = sprintf("chr1:%d:A:G", 100L * (1:5)),
    pip = seq(0.9, by = -0.15, length.out = 5L),
    alpha = NULL
) {
    if (is.null(alpha)) {
        alpha <- matrix(
            1 / length(variant_ids),
            nrow = 1,
            ncol = length(variant_ids)
        )
    }
    tl <- data.frame(
        variant_id = variant_ids,
        pip = pip,
        stringsAsFactors = FALSE
    )
    fit <- list(alpha = alpha, pip = setNames(pip, variant_ids), V = 0.1)
    FineMappingRow(variantIds = variant_ids, susieFit = fit, topLoci = tl)
}

.qep_makeGwasFmr <- function(
    studies = "G1",
    n_blocks = 1L,
    with_sketch = TRUE
) {
    entries <- vector("list", n_blocks)
    studyVec <- character(0)
    methodVec <- character(0)
    for (k in seq_len(n_blocks)) {
        # Different variants per block to avoid duplication.
        ids <- sprintf("chr1:%d:A:G", 100L * ((k - 1L) * 3L + (1:3)))
        entries[[k]] <- .qep_makeFmEntry(
            variant_ids = ids,
            pip = c(0.5, 0.2, 0.1)
        )
        studyVec <- c(studyVec, studies)
        methodVec <- c(methodVec, "susie")
    }
    GwasFineMappingResult(
        studyName = studyVec,
        method = methodVec,
        entry = entries,
        ldReference = if (with_sketch) .qep_makeHandle() else NULL
    )
}

.qep_makeQtlFmr <- function(
    contexts = "c1",
    traits = "t1",
    with_sketch = TRUE
) {
    n <- length(contexts) * length(traits)
    studies <- rep("Q1", n)
    ctx <- rep(contexts, length.out = n)
    trs <- rep(traits, each = length(contexts))[seq_len(n)]
    methods <- rep("susie", n)
    entries <- replicate(
        n,
        .qep_makeFmEntry(variant_ids = sprintf("chr1:%d:A:G", 100L * (1:5))),
        simplify = FALSE
    )
    QtlFineMappingResult(
        studyName = studies,
        context = ctx,
        trait = trs,
        method = methods,
        entry = entries,
        ldReference = if (with_sketch) .qep_makeHandle() else NULL
    )
}

# Mock that returns a plausible enrichment list.
# The shape qtlEnrichment() really returns: the C++ estimator's own field
# names (src/qtl_enrichment.h), not the pipeline's output columns. A mock that
# invented the output names is what let a total field-name mismatch --
# every estimate NA on real data -- sit undetected.
#
# `value` is the enrichment the pipeline should REPORT, so it is encoded here
# as the log-odds the estimator would have produced for it.
.qep_mockEnrichment <- function(value = 1.5) {
    logOdds <- log1p(value)
    function(gwasPip, susieQtlRegions, ...) {
        list(
            "Intercept" = -6.5,
            "sd (intercept)" = 0.5,
            "Enrichment (no shrinkage)" = logOdds,
            "Enrichment (w/ shrinkage)" = logOdds,
            "sd (no shrinkage)" = 0.2,
            "sd (w/ shrinkage)" = 0.1,
            "Alternative (coloc) p1" = 1e-4,
            "Alternative (coloc) p2" = 1e-4,
            "Alternative (coloc) p12" = 5e-6,
            "Effective MI rounds" = 25,
            unused_xqtl_variants = list()
        )
    }
}

# ===========================================================================
# End-to-end against the real estimator (no mock)
# ===========================================================================

test_that("qtlEnrichmentPipeline: the real estimator fills the value columns", {
    # Unmocked on purpose. The C++ kernel's field names are the contract, and
    # a mismatch between them and this pipeline's output columns is invisible
    # to a mocked run -- which is how every estimate came back NA on real data
    # while the mocked tests stayed green.
    ids <- sprintf("chr1:%d:A:G", 100L * seq_len(50L))
    gwasPip <- rep(0.01, 50L)
    gwasPip[c(5L, 20L, 35L)] <- c(0.8, 0.6, 0.9)
    alpha <- matrix(0.001, nrow = 2L, ncol = 50L)
    alpha[1L, 5L] <- 0.95
    alpha[2L, 20L] <- 0.95
    alpha <- alpha / rowSums(alpha)
    sketch <- .qep_makeHandle()
    gfmr <- GwasFineMappingResult(
        studyName = "G1",
        method = "susie",
        entry = list(FineMappingRow(
            variantIds = ids,
            susieFit = list(
                alpha = matrix(1 / 50, nrow = 1L, ncol = 50L),
                pip = setNames(gwasPip, ids),
                V = 0.5
            ),
            topLoci = data.frame(variant_id = ids, pip = gwasPip)
        )),
        ldReference = sketch
    )
    qfmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(
            variantIds = ids,
            susieFit = list(
                alpha = alpha,
                pip = setNames(colSums(alpha), ids),
                V = c(0.5, 0.3)
            ),
            topLoci = data.frame(variant_id = ids, pip = colSums(alpha))
        )),
        ldReference = sketch
    )
    suppressWarnings(capture.output(
        out <- qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = qfmr,
            qtlEnrichmentParam = QtlEnrichmentParam(impN = 5),
            seed = 1L
        )
    ))
    expect_equal(nrow(out), 1L)
    valueCols <- names(pecotmr:::.enrNaEnrichment())
    expect_true(all(valueCols %in% colnames(out)))
    expect_true(all(is.finite(unlist(out[, valueCols]))))
    # The multiplicative factor and the log-odds must agree.
    expect_equal(out$enrichment, expm1(out$enrichmentLogOdds))
})

# ===========================================================================
# Either side may be a QTL or a GWAS fine-mapping result
# ===========================================================================

test_that("qtlEnrichmentPipeline: keys a QTL outcome side by its trait", {
    # Two molecular traits over the same variants with DIFFERENT pips. Keyed by
    # study alone they collide -- the pipeline would abort on conflicting PIPs
    # -- and one estimate would stand in for both traits.
    outcome <- QtlFineMappingResult(
        studyName = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t2"),
        method = c("susie", "susie"),
        entry = list(
            .qep_makeFmEntry(pip = c(0.9, 0.05, 0.02, 0.02, 0.01)),
            .qep_makeFmEntry(pip = c(0.1, 0.7, 0.1, 0.05, 0.05))
        ),
        ldReference = .qep_makeHandle()
    )
    local_mocked_bindings(
        .qtlEnrichmentCoded = .qep_mockEnrichment(2.0),
        .package = "pecotmr"
    )
    out <- qtlEnrichmentPipeline(
        gwasFineMappingResult = outcome,
        qtlFineMappingResult = .qep_makeQtlFmr()
    )
    expect_equal(nrow(out), 2L)
    expect_setequal(out$gwasStudy, "Q1")
    expect_setequal(out$gwasContext, "c1")
    expect_setequal(out$gwasTrait, c("t1", "t2"))
})

test_that("qtlEnrichmentPipeline: pairs two GWAS collections", {
    local_mocked_bindings(
        .qtlEnrichmentCoded = .qep_mockEnrichment(2.0),
        .package = "pecotmr"
    )
    out <- qtlEnrichmentPipeline(
        gwasFineMappingResult = .qep_makeGwasFmr(studies = "G1"),
        qtlFineMappingResult = .qep_makeGwasFmr(studies = "G2")
    )
    expect_equal(nrow(out), 1L)
    expect_equal(out$gwasStudy, "G1")
    expect_equal(out$qtlStudy, "G2")
    # Neither side has a context or trait axis, so all three read NA.
    expect_true(all(is.na(c(
        out$gwasContext,
        out$gwasTrait,
        out$qtlContext
    ))))
})

test_that("qtlEnrichmentPipeline: a QTL outcome may carry no ldReference", {
    # The RSS-derived requirement is a GWAS-side contract; an individual-level
    # QTL outcome has no panel to check.
    local_mocked_bindings(
        .qtlEnrichmentCoded = .qep_mockEnrichment(2.0),
        .package = "pecotmr"
    )
    out <- qtlEnrichmentPipeline(
        gwasFineMappingResult = .qep_makeQtlFmr(with_sketch = FALSE),
        qtlFineMappingResult = .qep_makeQtlFmr()
    )
    expect_equal(nrow(out), 1L)
})

# ===========================================================================
# Input-type validation
# ===========================================================================

test_that("qtlEnrichmentPipeline: rejects non-GwasFineMappingResult gwasFmr", {
    qfmr <- .qep_makeQtlFmr()
    expect_error(
        qtlEnrichmentPipeline(
            gwasFineMappingResult = "no",
            qtlFineMappingResult = qfmr
        ),
        "Must inherit from class 'GwasFineMappingResult'"
    )
})

test_that("qtlEnrichmentPipeline: rejects non-QtlFineMappingResult qtlFmr", {
    gfmr <- .qep_makeGwasFmr()
    expect_error(
        qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = "no"
        ),
        "Must inherit from class 'GwasFineMappingResult'"
    )
})

test_that("qtlEnrichmentPipeline: NULL ldReference on the GWAS side errors", {
    gfmr <- .qep_makeGwasFmr(with_sketch = FALSE)
    qfmr <- .qep_makeQtlFmr()
    expect_error(
        qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = qfmr
        ),
        "row must have an ldReference"
    )
})

test_that("qtlEnrichmentPipeline: ldReference mismatch errors", {
    # Build the QTL with a sketch carrying a different sample set.
    gfmr <- .qep_makeGwasFmr()
    qSketch <- .qep_makeHandle()
    qSketch@sampleIds <- paste0("z", seq_len(nSamples(qSketch)))
    qfmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(.qep_makeFmEntry()),
        ldReference = qSketch
    )
    expect_error(
        qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = qfmr
        ),
        "different sample sets"
    )
})

# ===========================================================================
# Per-study / per-context iteration via mocked qtlEnrichment
# ===========================================================================

test_that("qtlEnrichmentPipeline: returns one row per (gwasStudy, qtlStudy, qtlContext) triple", {
    gfmr <- .qep_makeGwasFmr()
    qfmr <- .qep_makeQtlFmr(contexts = c("c1", "c2"))
    local_mocked_bindings(
        .qtlEnrichmentCoded = .qep_mockEnrichment(2.0),
        .package = "pecotmr"
    )
    out <- qtlEnrichmentPipeline(
        gwasFineMappingResult = gfmr,
        qtlFineMappingResult = qfmr
    )
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2L) # 1 GWAS study * 1 QTL study * 2 contexts
    expect_setequal(out$gwasStudy, "G1")
    expect_setequal(out$qtlStudy, "Q1")
    expect_setequal(out$qtlContext, c("c1", "c2"))
    expect_equal(out$enrichment, c(2.0, 2.0))
})

test_that("qtlEnrichmentPipeline: distinguishes two QTL studies that share a context label", {
    # Build a QtlFineMappingResult with two studies (Q1, Q2) both
    # tagging the same context "shared_ctx". Per-study filtering must
    # keep them separate (a context-only filter would merge them and
    # produce a single enrichment row instead of two).
    e1 <- .qep_makeFmEntry(variant_ids = sprintf("chr1:%d:A:G", 100L * (1:5)))
    e2 <- .qep_makeFmEntry(variant_ids = sprintf("chr1:%d:A:G", 100L * (1:5)))
    qfmr <- QtlFineMappingResult(
        studyName = c("Q1", "Q2"),
        context = c("shared_ctx", "shared_ctx"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2),
        ldReference = .qep_makeHandle()
    )
    gfmr <- .qep_makeGwasFmr()
    capturedRegions <- list()
    local_mocked_bindings(
        .qtlEnrichmentCoded = function(gwasPip, susieQtlRegions, ...) {
            capturedRegions[[length(capturedRegions) + 1L]] <<- susieQtlRegions
            .qep_mockEnrichment(2.0)()
        },
        .package = "pecotmr"
    )
    out <- qtlEnrichmentPipeline(
        gwasFineMappingResult = gfmr,
        qtlFineMappingResult = qfmr
    )
    expect_equal(nrow(out), 2L)
    expect_setequal(out$qtlStudy, c("Q1", "Q2"))
    # Each call to qtlEnrichment sees exactly one region (the per-study one),
    # not both regions merged together.
    expect_true(all(lengths(capturedRegions) == 1L))
})

test_that("qtlEnrichmentPipeline: qtlEnrichment failure produces a warning + skip", {
    gfmr <- .qep_makeGwasFmr()
    qfmr <- .qep_makeQtlFmr()
    local_mocked_bindings(
        .qtlEnrichmentCoded = function(...) stop("synthetic failure"),
        .package = "pecotmr"
    )
    expect_warning(
        out <- qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = qfmr
        ),
        "qtlEnrichment failed"
    )
    expect_equal(nrow(out), 0L)
})

test_that("qtlEnrichmentPipeline: empty input collections yield the empty schema", {
    # Build a GwasFineMappingResult whose entries have empty fits so the PIP
    # vector is empty.
    emptyEntry <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(), # no pip -> .enrBuildGwasPipVector returns numeric(0)
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.1,
            stringsAsFactors = FALSE
        )
    )
    gfmr <- GwasFineMappingResult(
        studyName = "G1",
        method = "susie",
        entry = list(emptyEntry),
        ldReference = .qep_makeHandle()
    )
    qfmr <- .qep_makeQtlFmr()
    expect_warning(
        out <- qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = qfmr
        ),
        "no usable PIPs"
    )
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 0L)
    expect_setequal(
        colnames(out),
        c(
            "gwasStudy",
            "gwasContext",
            "gwasTrait",
            "gwasMethod",
            "qtlStudy",
            "qtlContext",
            "qtlMethod",
            names(pecotmr:::.enrNaEnrichment())
        )
    )
})

# ===========================================================================
# Internal helpers: .enrBuildGwasPipVector + .enrBuildQtlRegionsList
# ===========================================================================

test_that(".enrBuildGwasPipVector: extracts pip per study", {
    gfmr <- .qep_makeGwasFmr()
    out <- pecotmr:::.enrBuildGwasPipVector(gfmr, .qep_gwasIdent("G1"))
    expect_equal(length(out), 3L)
    expect_setequal(names(out), sprintf("chr1:%d:A:G", 100L * (1:3)))
})

test_that(".enrBuildGwasPipVector: deduplicates identical PIPs across blocks", {
    # Two rows under the same study, same method, different LD blocks
    # (distinct blockIds auto-supplied by the constructor) — the
    # genome-wide multi-block shape. Both rows share v1 with the same
    # PIP, so dedup must collapse it.
    e1 <- .qep_makeFmEntry(
        variant_ids = c("chr1:100:A:G", "chr1:200:A:G"),
        pip = c(0.5, 0.2)
    )
    e2 <- .qep_makeFmEntry(
        variant_ids = c("chr1:100:A:G", "chr1:300:A:G"),
        pip = c(0.5, 0.4)
    )
    g <- GwasFineMappingResult(
        studyName = c("G1", "G1"),
        method = c("susie", "susie"),
        entry = list(e1, e2),
        ldReference = .qep_makeHandle()
    )
    out <- pecotmr:::.enrBuildGwasPipVector(g, .qep_gwasIdent("G1"))
    expect_setequal(
        names(out),
        c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G")
    )
})

test_that(".enrBuildGwasPipVector: conflicting PIPs across blocks errors", {
    # The shared variant chr1:100 appears in two OVERLAPPING blocks with
    # different PIPs — .enrBuildGwasPipVector must refuse to merge them
    # silently. The two rows must span different ranges, because row identity
    # is (study, method, range): two rows over the same span are the same
    # block by definition, not two blocks.
    e1 <- .qep_makeFmEntry(
        variant_ids = c("chr1:100:A:G", "chr1:200:A:G"),
        pip = c(0.5, 0.1),
        alpha = matrix(c(0.5, 0.1), nrow = 1)
    )
    e2 <- .qep_makeFmEntry(
        variant_ids = c("chr1:100:A:G", "chr1:300:A:G"),
        pip = c(0.8, 0.1),
        alpha = matrix(c(0.8, 0.1), nrow = 1)
    )
    g <- GwasFineMappingResult(
        studyName = c("G1", "G1"),
        method = c("susie", "susie"),
        entry = list(e1, e2),
        ldReference = .qep_makeHandle()
    )
    expect_error(
        pecotmr:::.enrBuildGwasPipVector(g, .qep_gwasIdent("G1"), "error"),
        "conflicting PIPs"
    )
})

# ===========================================================================
# qtlEnrichmentPipeline(): a variant repeated across GWAS blocks
# ===========================================================================

# Two GWAS blocks sharing their boundary variant chr1:200, to which each
# block's fit gives its own PIP. `region` gives the blocks' windows.
.qep_boundaryGwasFmr <- function(region = NULL, pips = c(0.3, 0.6)) {
    e1 <- .qep_makeFmEntry(
        variant_ids = c("chr1:100:A:G", "chr1:200:A:G"),
        pip = c(0.5, pips[[1L]])
    )
    e2 <- .qep_makeFmEntry(
        variant_ids = c("chr1:200:A:G", "chr1:300:A:G"),
        pip = c(pips[[2L]], 0.1)
    )
    GwasFineMappingResult(
        studyName = c("G1", "G1"),
        method = c("susie", "susie"),
        entry = list(e1, e2),
        region = region,
        ldReference = .qep_makeHandle()
    )
}

.qep_windows <- function(starts, ends) {
    GenomicRanges::GRanges("chr1", IRanges::IRanges(starts, ends))
}

# The GWAS PIP vector qtlEnrichmentPipeline() hands the estimator.
.qep_capturedGwasPip <- function(gwasFmr, ...) {
    captured <- NULL
    mock <- .qep_mockEnrichment()
    local_mocked_bindings(
        .qtlEnrichmentCoded = function(gwasPip, susieQtlRegions, ...) {
            # The pipeline hands the PIPs with their variants; named here as
            # the variants' ids.
            variants <- list(...)$gwasVariants
            captured <<- set_names(
                gwasPip,
                pecotmr:::.variantTableIds(variants)
            )
            mock(gwasPip, susieQtlRegions, ...)
        },
        .package = "pecotmr"
    )
    qtlEnrichmentPipeline(gwasFmr, .qep_makeQtlFmr(), verbose = FALSE, ...)
    captured
}

test_that("qtlEnrichmentPipeline: keeps the largest of conflicting PIPs", {
    pip <- .qep_capturedGwasPip(.qep_boundaryGwasFmr())
    expect_equal(
        pip,
        c(`chr1:100:A:G` = 0.5, `chr1:200:A:G` = 0.6, `chr1:300:A:G` = 0.1)
    )
})

test_that("qtlEnrichmentPipeline: duplicatePips = 'error' refuses conflicts", {
    expect_error(
        .qep_capturedGwasPip(.qep_boundaryGwasFmr(), duplicatePips = "error"),
        "'chr1:200:A:G' appears with conflicting PIPs"
    )
    expect_error(
        .qep_capturedGwasPip(.qep_boundaryGwasFmr(), duplicatePips = "first"),
        "must be one of"
    )
})

test_that("qtlEnrichmentPipeline: 'blockEnd' keeps the block ending at the variant", {
    # The upstream block's window ends at chr1:200; the downstream one starts
    # there.
    pip <- .qep_capturedGwasPip(
        .qep_boundaryGwasFmr(.qep_windows(c(1, 200), c(200, 400))),
        duplicatePips = "blockEnd"
    )
    expect_equal(pip[["chr1:200:A:G"]], 0.3)
    expect_length(pip, 3L)
})

test_that("qtlEnrichmentPipeline: 'blockEnd' refuses a copy it cannot choose", {
    # No window given, or none ending at the variant: nothing to decide by.
    for (region in list(NULL, .qep_windows(c(1, 150), c(250, 400)))) {
        expect_error(
            .qep_capturedGwasPip(
                .qep_boundaryGwasFmr(region),
                duplicatePips = "blockEnd"
            ),
            "cannot choose between them: 0 of its blocks"
        )
    }
})

test_that("qtlEnrichmentPipeline: agreeing copies merge under every policy", {
    for (policy in c("max", "error", "blockEnd")) {
        pip <- .qep_capturedGwasPip(
            .qep_boundaryGwasFmr(pips = c(0.3, 0.3)),
            duplicatePips = policy
        )
        expect_equal(
            pip,
            c(`chr1:100:A:G` = 0.5, `chr1:200:A:G` = 0.3, `chr1:300:A:G` = 0.1),
            info = policy
        )
    }
})

test_that(".enrBuildQtlRegionsList: returns per-entry fit shapes for a (study, context) hit", {
    qfmr <- .qep_makeQtlFmr(contexts = c("c1", "c2"))
    out <- pecotmr:::.enrBuildQtlRegionsList(qfmr, .qep_qtlIdent("Q1", "c1"))
    expect_equal(length(out), 1L)
    expect_true(!is.null(out[[1L]]$alpha))
    expect_true(!is.null(out[[1L]]$pip))
})

test_that(".enrBuildQtlRegionsList: returns empty list when the (study, context) tuple is absent", {
    qfmr <- .qep_makeQtlFmr(contexts = c("c1", "c2"))
    # Correct context but wrong study -> no hit, even though context exists.
    expect_equal(
        length(pecotmr:::.enrBuildQtlRegionsList(
            qfmr,
            .qep_qtlIdent("Q_ghost", "c1")
        )),
        0L
    )
    # Correct study but wrong context.
    expect_equal(
        length(pecotmr:::.enrBuildQtlRegionsList(
            qfmr,
            .qep_qtlIdent("Q1", "c_ghost")
        )),
        0L
    )
})

# ===========================================================================
# qtlEnrichment() — kernel wrapper + real C++ integration
# These tests deliberately do NOT mock qtlEnrichmentCpp so the C++
# kernel in src/qtl_enrichment.cpp gets coverage. The wrapper itself
# (R/qtlEnrichmentPipeline.R::qtlEnrichment) is exercised here directly.
# ===========================================================================

# Build a small (gwasPip, susieQtlRegions) fixture with a sparse causal
# signal at known indices so the C++ enrichment routine has something
# meaningful to compute.
.qep_makeRealKernelInputs <- function(
    seed = 42,
    nSnps = 50,
    causalIdx = c(5, 20, 35),
    # Trimmed to the number of causal variants actually asked for; callers
    # that pass a shorter causalIdx would otherwise warn on the assignment.
    causalPips = c(0.8, 0.6, 0.9)[seq_along(causalIdx)],
    maxNumSingleEffects = 2L
) {
    set.seed(seed)
    variantNames <- paste0("1:", seq_len(nSnps), ":A:G")
    gwasPip <- rep(0.01, nSnps)
    gwasPip[causalIdx] <- causalPips
    names(gwasPip) <- variantNames

    alpha <- matrix(1 / nSnps, nrow = maxNumSingleEffects, ncol = nSnps)
    alpha[1, ] <- 0.001
    alpha[1, causalIdx[1]] <- 0.95
    alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])
    alpha[2, ] <- 0.001
    alpha[2, causalIdx[2]] <- 0.95
    alpha[2, ] <- alpha[2, ] / sum(alpha[2, ])
    pip <- colSums(alpha)
    names(pip) <- variantNames
    susieFits <- list(
        fit1 = list(pip = pip, alpha = alpha, prior_variance = c(0.5, 0.3))
    )
    list(
        gwasPip = gwasPip,
        susieQtlRegions = susieFits,
        variantNames = variantNames
    )
}

test_that("qtlEnrichment: real C++ kernel returns the expected keys (numGwas + piQtl supplied)", {
    fx <- .qep_makeRealKernelInputs()
    res <- qtlEnrichment(
        gwasPip = fx$gwasPip,
        susieQtlRegions = fx$susieQtlRegions,
        numGwas = 5000,
        piQtl = 0.5,
        lambda = 1,
        impN = 5,
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE
    )
    expect_type(res, "list")
    # Flat, and addressable by name: the estimates used to sit inside an
    # unnamed first element, where every by-name read of them returned NULL.
    en <- res
    expectedKeys <- c(
        "Intercept",
        "Enrichment (no shrinkage)",
        "Enrichment (w/ shrinkage)",
        "sd (no shrinkage)",
        "sd (w/ shrinkage)",
        "Alternative (coloc) p1",
        "Alternative (coloc) p2",
        "Alternative (coloc) p12"
    )
    expect_setequal(intersect(expectedKeys, names(en)), expectedKeys)
    expect_true(all(is.finite(unlist(en[expectedKeys]))))
})

test_that("qtlEnrichment: a single MI round is estimable, not NaN", {
    # With one round the Bessel-corrected dispersion is 0/0, so no
    # |x - mean| <= 3 * NaN comparison held, the outlier filter dropped the
    # only round, and dividing by (m - 1) = 0 made every estimate NaN.
    fx <- .qep_makeRealKernelInputs()
    res <- qtlEnrichment(
        gwasPip = fx$gwasPip,
        susieQtlRegions = fx$susieQtlRegions,
        numGwas = 5000,
        piQtl = 0.5,
        impN = 1,
        BPPARAM = BiocParallel::SerialParam(),
        verbose = FALSE,
        seed = 1L
    )
    expect_equal(res[["Effective MI rounds"]], 1)
    expect_true(all(is.finite(c(
        res[["Enrichment (w/ shrinkage)"]],
        res[["sd (w/ shrinkage)"]],
        res[["Intercept"]],
        res[["Alternative (coloc) p12"]]
    ))))
})

test_that("qtlEnrichment: numGwas omitted -> estimates piGwas from data + warns", {
    fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
    expect_warning(
        res <- qtlEnrichment(
            gwasPip = fx$gwasPip,
            susieQtlRegions = fx$susieQtlRegions,
            piQtl = 0.5,
            impN = 5,
            BPPARAM = BiocParallel::SerialParam(),
            verbose = FALSE
        ),
        "numGwas is not provided"
    )
    expect_type(res, "list")
})

test_that("qtlEnrichment: piQtl omitted -> estimates from susieQtlRegions + warns", {
    fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
    expect_warning(
        res <- qtlEnrichment(
            gwasPip = fx$gwasPip,
            susieQtlRegions = fx$susieQtlRegions,
            numGwas = 3000,
            impN = 5,
            BPPARAM = BiocParallel::SerialParam(),
            verbose = FALSE
        ),
        "piQtl is not provided"
    )
    expect_type(res, "list")
})

test_that("qtlEnrichment: errors when piGwas resolves to zero", {
    fx <- .qep_makeRealKernelInputs()
    zeroGwas <- rep(0, length(fx$gwasPip))
    names(zeroGwas) <- names(fx$gwasPip)
    expect_error(
        qtlEnrichment(
            gwasPip = zeroGwas,
            susieQtlRegions = fx$susieQtlRegions,
            numGwas = 5000,
            piQtl = 0.5,
            BPPARAM = BiocParallel::SerialParam(),
            verbose = FALSE
        ),
        "No association signal found"
    )
})

test_that("qtlEnrichment: errors when piQtl resolves to zero", {
    fx <- .qep_makeRealKernelInputs()
    expect_error(
        qtlEnrichment(
            gwasPip = fx$gwasPip,
            susieQtlRegions = fx$susieQtlRegions,
            numGwas = 5000,
            piQtl = 0,
            BPPARAM = BiocParallel::SerialParam(),
            verbose = FALSE
        ),
        "No QTL associated"
    )
})

test_that("qtlEnrichment: errors when gwasPip has no names", {
    fx <- .qep_makeRealKernelInputs()
    unnamed <- unname(fx$gwasPip)
    expect_error(
        qtlEnrichment(
            gwasPip = unnamed,
            susieQtlRegions = fx$susieQtlRegions,
            numGwas = 5000,
            piQtl = 0.5,
            BPPARAM = BiocParallel::SerialParam(),
            verbose = FALSE
        ),
        "Variant names are missing in gwasPip"
    )
})

test_that("qtlEnrichment: errors when susieQtlRegions$pip lacks names", {
    fx <- .qep_makeRealKernelInputs()
    fx$susieQtlRegions$fit1$pip <- unname(fx$susieQtlRegions$fit1$pip)
    expect_error(
        qtlEnrichment(
            gwasPip = fx$gwasPip,
            susieQtlRegions = fx$susieQtlRegions,
            numGwas = 5000,
            piQtl = 0.5,
            BPPARAM = BiocParallel::SerialParam(),
            verbose = FALSE
        ),
        "Variant names are missing in susieQtlRegions"
    )
})

test_that("qtlEnrichment: tracks unmatched QTL variants, either alignNames", {
    fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
    # Inject a couple of variant IDs into the QTL fit that don't exist
    # in the GWAS PIP vector.
    newNames <- names(fx$susieQtlRegions$fit1$pip)
    newNames[1:2] <- c("1:9999:A:G", "1:9998:A:G")
    names(fx$susieQtlRegions$fit1$pip) <- newNames
    colnames(fx$susieQtlRegions$fit1$alpha) <- newNames
    # alignNames = FALSE is the cheap set-membership branch
    # qtlEnrichmentPipeline uses once it has aligned the QTL names against
    # the union GWAS panel; both branches report the two injected ids.
    for (alignNames in c(TRUE, FALSE)) {
        res <- qtlEnrichment(
            gwasPip = fx$gwasPip,
            susieQtlRegions = fx$susieQtlRegions,
            numGwas = 3000,
            piQtl = 0.5,
            impN = 5,
            BPPARAM = BiocParallel::SerialParam(),
            verbose = FALSE,
            alignNames = alignNames
        )
        expect_true(
            all(
                c("1:9999:A:G", "1:9998:A:G") %in%
                    unlist(res$unused_xqtl_variants)
            ),
            info = str_c("alignNames = ", alignNames)
        )
    }
})


# One-hot QTL fixtures for the enrichment kernel. Every alpha row puts all of
# its mass on one variant, so each imputation round samples the same variants
# whatever the random stream: the estimates are exactly reproducible, and do
# not depend on how a standard library implements std::discrete_distribution.

.qeh_ids <- function(n = 200L) {
    sprintf("chr1:%d:A:G", 100L * seq_len(n))
}

# A GWAS PIP vector over .qeh_ids(): background 0.001, with high PIPs on the
# last ten variants.
.qeh_gwasPip <- function() {
    stats::setNames(
        c(
            rep(0.001, 190),
            0.9,
            0.8,
            0.5,
            0.3,
            0.2,
            0.1,
            0.05,
            0.02,
            0.01,
            0.6
        ),
        .qeh_ids()
    )
}

# A region over variants `idx` of the GWAS panel whose effects each sit on the
# variant at position `hits[l]` within the region. `ids` relabels the region's
# variants, e.g. to put some outside the GWAS panel.
.qeh_oneHotRegion <- function(
    idx,
    hits,
    priorVariance = c(0.5, 0.2)[seq_along(hits)],
    ids = .qeh_ids()[idx]
) {
    alpha <- matrix(0, length(hits), length(idx))
    alpha[cbind(seq_along(hits), hits)] <- 1
    list(
        alpha = alpha,
        pip = stats::setNames(as.numeric(colSums(alpha) > 0), ids),
        prior_variance = priorVariance
    )
}

# Three regions: r1's effects sit on variants 191 (GWAS PIP 0.9) and 152, r2's
# on 192 (0.8) and 195 (0.2), r3's on two background variants.
.qeh_oneHotRegions <- function() {
    list(
        r1 = .qeh_oneHotRegion(150:200, c(42, 3)),
        r2 = .qeh_oneHotRegion(170:200, c(23, 26)),
        r3 = .qeh_oneHotRegion(1:60, c(5, 40))
    )
}

# ===========================================================================
# qtlEnrichment(): pinned estimates, null fits, fit validation, threads
# Built on the one-hot fixtures above, whose sampling is
# deterministic.
# ===========================================================================

.qep_oneHotRun <- function(
    regions = .qeh_oneHotRegions(),
    numGwas = 5000,
    piQtl = 0.01,
    ...
) {
    qtlEnrichment(
        .qeh_gwasPip(),
        regions,
        numGwas = numGwas,
        piQtl = piQtl,
        impN = 3,
        seed = 1L,
        verbose = FALSE,
        ...
    )
}

test_that("qtlEnrichment: one-hot fixture reproduces the pinned estimates", {
    # Recorded from the kernel before its fits were integer-coded and null
    # fits allowed. Without shrinkage the shrunk and unshrunk a1 coincide.
    pinned <- c(
        `Intercept` = -21.896577386392181,
        `sd (intercept)` = 11.186832846318783,
        `Enrichment (no shrinkage)` = 19.284745626108926,
        `Enrichment (w/ shrinkage)` = 19.284745626108926,
        `sd (no shrinkage)` = 13.735002410315799,
        `sd (w/ shrinkage)` = 13.735002410315799,
        `Alternative (coloc) p1` = 3.0624743716206618e-10,
        `Alternative (coloc) p2` = 0.0093161918064975369,
        `Alternative (coloc) p12` = 0.00068380819350246172,
        `Effective MI rounds` = 3
    )
    res <- .qep_oneHotRun(lambda = 0)
    expect_equal(unlist(res[names(pinned)]), pinned, tolerance = 1e-10)
})

test_that("qtlEnrichment: one-hot fixture reproduces the pinned shrunk a1", {
    keys <- c(
        "Enrichment (no shrinkage)",
        "Enrichment (w/ shrinkage)",
        "sd (no shrinkage)",
        "sd (w/ shrinkage)"
    )
    single <- .qep_oneHotRun(lambda = 1)
    double <- .qep_oneHotRun(
        lambda = 1,
        doubleShrinkage = TRUE,
        besselCorrection = FALSE
    )
    expect_equal(
        unname(unlist(single[keys])),
        c(
            19.284745626108926,
            0.10168582132370414,
            13.735002410315799,
            0.99736008393071596
        ),
        tolerance = 1e-10
    )
    expect_equal(
        unname(unlist(double[keys])),
        c(
            19.284745626108926,
            0.050977308844208112,
            13.735002410315799,
            0.70617157949200415
        ),
        tolerance = 1e-10
    )
})

test_that("qtlEnrichment: a null QTL fit samples nothing instead of aborting", {
    keys <- c(
        "Intercept",
        "Enrichment (no shrinkage)",
        "sd (no shrinkage)",
        "Alternative (coloc) p12"
    )
    base <- .qep_oneHotRun(lambda = 0)
    # Every effect of prior variance 0: the kernel used to stop on this.
    nullFit <- .qeh_oneHotRegion(61:80, c(1, 2), priorVariance = c(0, 0))
    res <- .qep_oneHotRun(
        c(.qeh_oneHotRegions(), list(n = nullFit)),
        lambda = 0
    )
    expect_equal(unlist(res[keys]), unlist(base[keys]))
    expect_equal(res[["Null QTL fits"]], 1L)
    # A fit already trimmed to its non-null effects has a 0-row alpha.
    trimmed <- list(
        alpha = matrix(0, 0L, 3L),
        pip = stats::setNames(rep(0, 3), .qeh_ids()[61:63]),
        prior_variance = numeric(0)
    )
    res <- .qep_oneHotRun(
        c(.qeh_oneHotRegions(), list(n = trimmed)),
        lambda = 0
    )
    expect_equal(unlist(res[keys]), unlist(base[keys]))
})

test_that("qtlEnrichment: refuses a tuple made only of null fits", {
    nullFit <- .qeh_oneHotRegion(1:10, c(1, 2), priorVariance = c(0, 0))
    expect_error(.qep_oneHotRun(list(a = nullFit)), "nothing to sample")
})

test_that("qtlEnrichment: refuses malformed fits, naming them", {
    ok <- .qeh_oneHotRegion(1:10, c(1, 2))
    run <- function(bad) .qep_oneHotRun(list(good = ok, bad = bad))
    # A null_weight fit carries one more alpha column than it has variants.
    extraColumn <- ok
    extraColumn$alpha <- cbind(ok$alpha, 0)
    expect_error(run(extraColumn), "bad: alpha has 11 columns for 10 pip")
    shortV <- ok
    shortV$prior_variance <- 0.5
    expect_error(run(shortV), "bad: alpha has 2 rows for 1 prior variances")
    infV <- ok
    infV$prior_variance <- c(Inf, 0.2)
    expect_error(run(infV), "bad: prior variances must be finite")
    nanAlpha <- ok
    nanAlpha$alpha[1, 2] <- NaN
    expect_error(run(nanAlpha), "bad: alpha must be finite and non-negative")
    negativeAlpha <- ok
    negativeAlpha$alpha[1, 1:2] <- c(1.1, -0.1)
    expect_error(
        run(negativeAlpha),
        "bad: alpha must be finite and non-negative"
    )
    offSum <- ok
    offSum$alpha[1, ] <- offSum$alpha[1, ] / 2
    expect_error(run(offSum), "bad: an alpha row of positive prior variance")
})

test_that("qtlEnrichment: estimates do not depend on the thread count", {
    skip_on_os("windows")
    fx <- .qep_makeRealKernelInputs()
    run <- function(bp) {
        qtlEnrichment(
            fx$gwasPip,
            fx$susieQtlRegions,
            numGwas = 5000,
            piQtl = 0.5,
            impN = 6,
            BPPARAM = bp,
            seed = 7L,
            verbose = FALSE
        )
    }
    expect_identical(
        run(BiocParallel::SerialParam()),
        run(BiocParallel::MulticoreParam(2))
    )
})

test_that("qtlEnrichment: verbose = FALSE is silent", {
    expect_silent(.qep_oneHotRun())
})

test_that("qtlEnrichment: reports per-round diagnostics", {
    # r4's second effect sits on a variant the GWAS vector does not carry:
    # one of the eight variants sampled each round.
    r4 <- .qeh_oneHotRegion(
        1:10,
        c(3, 7),
        ids = c(.qeh_ids()[1:5], sprintf("chr2:%d:C:T", 1:5))
    )
    regions <- c(.qeh_oneHotRegions(), list(r4 = r4))
    res <- .qep_oneHotRun(regions, lambda = 0)
    expect_equal(res[["Missing proportion per MI round"]], rep(1 / 8, 3))
    expect_identical(res[["EM converged per MI round"]], rep(TRUE, 3))
    expect_identical(res[["Null QTL fits"]], 0L)
    expect_identical(
        res$unused_xqtl_variants$r4,
        sprintf("chr2:%d:C:T", 1:5)
    )
    expect_null(res$unused_xqtl_variants$r1)
    msgs <- testthat::capture_messages(qtlEnrichment(
        .qeh_gwasPip(),
        regions,
        numGwas = 5000,
        piQtl = 0.01,
        impN = 3,
        seed = 1L,
        verbose = TRUE
    ))
    expect_match(
        paste(msgs, collapse = ""),
        "missing from the GWAS PIP vector: 0.125"
    )
})

# An R transcription of fastenloc's enrichment estimator
# (github.com/xqwen/fastenloc, src/controller.cc: run_EM and the
# multiple-imputation post-processing of enrich_est), to check the C++ kernel
# against the original method. It takes the variants each imputation round
# samples directly, so a test can feed it what a one-hot fixture (above)
# makes every round of the kernel sample.

# One round's EM. `gwasPip` and `eqtlSample` run over every variant fastenloc
# reads, GWAS and QTL alike; a QTL variant the GWAS file lacks has GWAS PIP 0.
.flr_runEm <- function(gwasPip, eqtlSample, pGwas, pEqtl, totalSnp) {
    a0 <- log(pGwas / (1 - pGwas))
    a1 <- 0
    r0 <- exp(a0)
    r1 <- exp(a0 + a1)
    rm <- pGwas / (1 - pGwas)
    val <- gwasPip
    val[val == 1] <- 1 - 1e-8
    odds <- val / (1 - val)
    sampled <- eqtlSample == 1
    repeat {
        e0g0 <- (1 - pGwas) * (1 - pEqtl)
        e0g1 <- (1 - pEqtl) * pGwas
        e1g0 <- (1 - pGwas) * pEqtl
        e1g1 <- pGwas * pEqtl
        post0 <- r0 * (odds[!sampled] / rm)
        post0 <- post0 / (1 + post0)
        post1 <- r1 * (odds[sampled] / rm)
        post1 <- post1 / (1 + post1)
        e0g1 <- e0g1 + sum(post0)
        e0g0 <- e0g0 + sum(1 - post0)
        e1g1 <- e1g1 + sum(post1)
        e1g0 <- e1g0 + sum(1 - post1)
        e0g0 <- e0g0 + totalSnp - (e0g0 + e0g1 + e1g0 + e1g1)
        a1New <- log(e1g1 * e0g0 / (e1g0 * e0g1))
        if (abs(a1New - a1) < 0.01) {
            # fastenloc keeps the previous iteration's a0 here.
            return(c(
                a0 = a0,
                a1 = a1New,
                v0 = 1 / e0g1 + 1 / e0g0,
                v1 = 1 / e0g0 + 1 / e1g0 + 1 / e1g1 + 1 / e0g1
            ))
        }
        a1 <- a1New
        a0 <- log(e0g1 / e0g0)
        r0 <- exp(a0)
        r1 <- exp(a0 + a1)
    }
}

# The post-processing of the rounds' estimates (a matrix with rows a0, a1, v0,
# v1): outlier removal on the shrunk a1, which fastenloc's command line always
# applies; Rubin's rules; the intercept from the unshrunk a1; double
# shrinkage; and the coloc priors from the intercept and the shrunk a1.
# `lambda` is fastenloc's -s, where 0 means no shrinkage (a prior variance of
# 1e16). fastenloc takes the outliers' population sd in one pass; two passes
# give the same value without its cancellation.
.flr_combine <- function(rounds, pGwas, pEqtl, lambda) {
    pv <- if (lambda == 0) 1e16 else 1 / lambda
    shrunk <- rounds["a1", ] * pv / (pv + rounds["v1", ])
    mu <- mean(shrunk)
    sdev <- sqrt(mean((shrunk - mu)^2))
    rounds <- rounds[, !(abs(shrunk - mu) > 3 * sdev), drop = FALSE]
    m <- ncol(rounds)
    a0 <- rounds["a0", ]
    a1 <- rounds["a1", ]
    v1 <- rounds["v1", ]
    postA1 <- a1 * pv / (pv + v1)
    postVar <- 1 / (1 / pv + 1 / v1)
    sd0 <- sqrt(mean(rounds["v0", ]) + stats::var(a0) * (m + 1) / m)
    sd1 <- sqrt(mean(v1) + stats::var(a1) * (m + 1) / m)
    sd1Shrink <- sqrt(mean(postVar) + stats::var(postA1) * (m + 1) / m)
    intercept <- log(pGwas / (1 + pEqtl * exp(mean(a1)) - pEqtl - pGwas))
    a1Shrink <- pv * mean(postA1) / (pv + sd1Shrink^2)
    sd1Shrink <- sqrt(1 / (1 / pv + 1 / sd1Shrink^2))
    c(
        `Intercept` = intercept,
        `sd (intercept)` = sd0,
        `Enrichment (no shrinkage)` = mean(a1),
        `Enrichment (w/ shrinkage)` = a1Shrink,
        `sd (no shrinkage)` = sd1,
        `sd (w/ shrinkage)` = sd1Shrink,
        `Alternative (coloc) p1` = (1 - pEqtl) *
            exp(intercept) /
            (1 + exp(intercept)),
        `Alternative (coloc) p2` = pEqtl / (1 + exp(intercept + a1Shrink)),
        `Alternative (coloc) p12` = pEqtl *
            exp(intercept + a1Shrink) /
            (1 + exp(intercept + a1Shrink)),
        `Effective MI rounds` = m
    )
}

# fastenloc on qtlEnrichment()'s inputs, every round sampling the variant ids
# `sampled`. Its variants are the GWAS ones plus the QTL ones they lack, its
# total is -total_variants raised to their number, and its priors are the
# GWAS and QTL PIP mass over that total.
.flr_reference <- function(
    gwasPip,
    regions,
    sampled,
    numGwas = NULL,
    lambda = 1,
    impN = 3
) {
    qtlIds <- unique(unlist(
        lapply(regions, function(x) names(x$pip)),
        use.names = FALSE
    ))
    variants <- c(names(gwasPip), setdiff(qtlIds, names(gwasPip)))
    pip <- c(unname(gwasPip), rep(0, length(variants) - length(gwasPip)))
    totalSnp <- max(if (is.null(numGwas)) 0 else numGwas, length(variants))
    pGwas <- sum(gwasPip) / totalSnp
    pEqtl <- sum(unlist(lapply(regions, `[[`, "pip"))) / totalSnp
    round <- .flr_runEm(
        pip,
        as.integer(variants %in% sampled),
        pGwas,
        pEqtl,
        totalSnp
    )
    rounds <- matrix(
        round,
        nrow = length(round),
        ncol = impN,
        dimnames = list(names(round), NULL)
    )
    .flr_combine(rounds, pGwas, pEqtl, lambda)
}

# The variants every round samples from one-hot regions: the hit of each
# effect with positive prior variance.
.flr_oneHotSamples <- function(regions) {
    unlist(
        lapply(regions, function(x) {
            alpha <- x$alpha[x$prior_variance > 0, , drop = FALSE]
            names(x$pip)[max.col(alpha, ties.method = "first")]
        }),
        use.names = FALSE
    )
}

# ===========================================================================
# qtlEnrichment(): fastenloc semantics
# Checked against the R transcription of fastenloc above. fastenloc's EM keeps the previous iteration's
# intercept at convergence, so only sd (intercept) is not compared.
# ===========================================================================

# The one-hot regions plus r4, whose effects sit on a GWAS variant and on
# chr2:2, a variant the GWAS vector lacks; chr2:1 and chr2:3-5 are QTL-only
# variants no round samples.
.qep_qtlOnlyRegions <- function() {
    r4 <- .qeh_oneHotRegion(
        1:10,
        c(3, 7),
        ids = c(.qeh_ids()[1:5], sprintf("chr2:%d:C:T", 1:5))
    )
    c(.qeh_oneHotRegions(), list(r4 = r4))
}

.qep_parityKeys <- function() {
    c(
        "Intercept",
        "Enrichment (no shrinkage)",
        "Enrichment (w/ shrinkage)",
        "sd (no shrinkage)",
        "sd (w/ shrinkage)",
        "Alternative (coloc) p1",
        "Alternative (coloc) p2",
        "Alternative (coloc) p12",
        "Effective MI rounds"
    )
}

test_that("qtlEnrichment: matches fastenloc, QTL-only variants included", {
    keys <- .qep_parityKeys()
    for (regions in list(.qeh_oneHotRegions(), .qep_qtlOnlyRegions())) {
        for (numGwas in list(5000, NULL)) {
            for (lambda in c(1, 0.5, 0)) {
                info <- paste0(
                    length(regions),
                    " regions, numGwas = ",
                    deparse(numGwas),
                    ", lambda = ",
                    lambda
                )
                reference <- .flr_reference(
                    .qeh_gwasPip(),
                    regions,
                    .flr_oneHotSamples(regions),
                    numGwas = numGwas,
                    lambda = lambda
                )
                res <- suppressWarnings(qtlEnrichment(
                    .qeh_gwasPip(),
                    regions,
                    numGwas = numGwas,
                    lambda = lambda,
                    impN = 3,
                    doubleShrinkage = TRUE,
                    besselCorrection = FALSE,
                    seed = 1L,
                    verbose = FALSE
                ))
                expect_equal(
                    unlist(res[keys]),
                    reference[keys],
                    tolerance = 1e-8,
                    info = info
                )
            }
        }
    }
})

test_that("qtlEnrichment: the intercept uses the unshrunk enrichment", {
    res <- .qep_oneHotRun(lambda = 1)
    piGwas <- sum(.qeh_gwasPip()) / 5000
    a1 <- res[["Enrichment (no shrinkage)"]]
    expect_equal(
        res[["Intercept"]],
        log(piGwas / (1 + 0.01 * exp(a1) - 0.01 - piGwas))
    )
    a01 <- res[["Intercept"]] + res[["Enrichment (w/ shrinkage)"]]
    expect_equal(res[["Alternative (coloc) p12"]], 0.01 * plogis(a01))
})

test_that("qtlEnrichment: a QTL-only variant counts as a GWAS PIP of 0", {
    # Padding gwasPip with zeros for the QTL-only variants, at a fixed
    # variant total, changes nothing but the diagnostics.
    regions <- .qep_qtlOnlyRegions()
    padded <- c(
        .qeh_gwasPip(),
        stats::setNames(rep(0, 5), sprintf("chr2:%d:C:T", 1:5))
    )
    keys <- c(.qep_parityKeys(), "sd (intercept)")
    asGiven <- .qep_oneHotRun(regions, lambda = 1)
    res <- qtlEnrichment(
        padded,
        regions,
        numGwas = 5000,
        piQtl = 0.01,
        impN = 3,
        seed = 1L,
        verbose = FALSE
    )
    expect_equal(unlist(res[keys]), unlist(asGiven[keys]), tolerance = 1e-10)
    expect_equal(res[["Missing proportion per MI round"]], rep(0, 3))
    expect_equal(asGiven[["Missing proportion per MI round"]], rep(1 / 8, 3))
})

test_that("qtlEnrichment: estimates do not depend on variant or region order", {
    keys <- c(.qep_parityKeys(), "sd (intercept)")
    regions <- .qep_qtlOnlyRegions()
    base <- .qep_oneHotRun(regions, lambda = 1)
    gwasPip <- .qeh_gwasPip()
    res <- qtlEnrichment(
        rev(gwasPip),
        rev(regions),
        numGwas = 5000,
        piQtl = 0.01,
        impN = 3,
        seed = 1L,
        verbose = FALSE
    )
    expect_equal(unlist(res[keys]), unlist(base[keys]), tolerance = 1e-10)
})

test_that("qtlEnrichment: the variant total is the GWAS and QTL-only variants", {
    regions <- .qep_qtlOnlyRegions()
    # 200 GWAS variants and 5 QTL-only ones.
    expect_warning(
        res <- qtlEnrichment(
            .qeh_gwasPip(),
            regions,
            piQtl = 0.01,
            impN = 3,
            seed = 1L,
            verbose = FALSE
        ),
        "numGwas is not provided. Using the 205 distinct variants"
    )
    at205 <- .qep_oneHotRun(regions, numGwas = 205)
    keys <- .qep_parityKeys()
    expect_equal(unlist(res[keys]), unlist(at205[keys]))
    expect_warning(
        below <- .qep_oneHotRun(regions, numGwas = 100),
        "numGwas \\(100\\) is smaller than the 205 distinct variants"
    )
    expect_equal(unlist(below[keys]), unlist(at205[keys]))
})

test_that("qtlEnrichment: piQtl defaults to the summed QTL PIPs over the total", {
    regions <- .qep_qtlOnlyRegions()
    expect_warning(
        estimated <- qtlEnrichment(
            .qeh_gwasPip(),
            regions,
            numGwas = 5000,
            impN = 3,
            seed = 1L,
            verbose = FALSE
        ),
        "piQtl is not provided"
    )
    # Four regions of two one-hot effects: 8 units of PIP.
    supplied <- .qep_oneHotRun(regions, piQtl = 8 / 5000)
    keys <- c(.qep_parityKeys(), "sd (intercept)")
    expect_equal(unlist(estimated[keys]), unlist(supplied[keys]))
})

test_that("qtlEnrichment: refuses priors that do not sum to below 1", {
    expect_error(
        .qep_oneHotRun(piQtl = 0.9999),
        "must sum to less than 1. Supply numGwas"
    )
})

test_that("qtlEnrichment: refuses a variant total too small for the data", {
    # Every GWAS variant is sampled, so no variant is left that is neither
    # QTL nor causal: the EM's (0, 0) cell is negative.
    ids <- .qeh_ids(2L)
    fit <- list(
        alpha = diag(2),
        pip = stats::setNames(c(1, 1), ids),
        prior_variance = c(0.5, 0.2)
    )
    expect_error(
        suppressWarnings(qtlEnrichment(
            stats::setNames(c(0.5, 0.5), ids),
            list(a = fit),
            piQtl = 0.1,
            impN = 3,
            seed = 1L,
            verbose = FALSE
        )),
        "variant total is too small for the data in 3 of 3 MI round"
    )
})

# ===========================================================================
# qtlEnrichment(): aligning QTL variant names to the GWAS vector
# ===========================================================================

# `regions` with the pip names at positions `at` of region `name` replaced.
.qep_renameIds <- function(regions, name, at, ids) {
    nm <- names(regions[[name]]$pip)
    nm[at] <- ids
    names(regions[[name]]$pip) <- nm
    regions
}

.qep_estimates <- function(res) {
    unlist(res[c(.qep_parityKeys(), "sd (intercept)")])
}

test_that("qtlEnrichment: aligns swapped and reformatted QTL ids to GWAS ids", {
    regions <- .qeh_oneHotRegions()
    # r1's effects sit on chr1:19100 and chr1:15200 (positions 42 and 3),
    # r2's first on chr1:19200 (position 23).
    renamed <- .qep_renameIds(
        regions,
        "r1",
        c(42, 3),
        c("chr1:19100:G:A", "1_15200_A_G")
    )
    renamed <- .qep_renameIds(renamed, "r2", 23, "chr1:19200:G:A")
    base <- .qep_oneHotRun(regions)
    aligned <- .qep_oneHotRun(renamed)
    expect_equal(.qep_estimates(aligned), .qep_estimates(base))
    expect_true(all(vapply(aligned$unused_xqtl_variants, is.null, TRUE)))
    # Unaligned, the three names match nothing and count as QTL-only.
    unaligned <- .qep_oneHotRun(renamed, alignNames = FALSE)
    expect_setequal(
        unlist(unaligned$unused_xqtl_variants),
        c("chr1:19100:G:A", "1_15200_A_G", "chr1:19200:G:A")
    )
    expect_false(isTRUE(all.equal(
        .qep_estimates(unaligned),
        .qep_estimates(base)
    )))
})

test_that("qtlEnrichment: an unparseable id leaves the others allele-aware", {
    # An rsID among a tuple's names used to turn the whole match into exact
    # string identity, leaving the swapped ids unmatched.
    withRsid <- .qep_renameIds(.qeh_oneHotRegions(), "r1", 10, "rs999")
    swapped <- .qep_renameIds(withRsid, "r1", 42, "chr1:19100:G:A")
    res <- .qep_oneHotRun(swapped)
    expect_identical(res$unused_xqtl_variants$r1, "rs999")
    expect_equal(
        .qep_estimates(res),
        .qep_estimates(.qep_oneHotRun(withRsid))
    )
})

test_that("qtlEnrichment: an exactly spelled GWAS id beats a swapped twin", {
    # The GWAS vector carries chr1:19100 in both orientations, with
    # different PIPs; each QTL spelling must take its own.
    gwasPip <- c(.qeh_gwasPip(), `chr1:19100:G:A` = 0.05)
    run <- function(regions, alignNames) {
        qtlEnrichment(
            gwasPip,
            regions,
            numGwas = 5000,
            piQtl = 0.01,
            impN = 3,
            seed = 1L,
            verbose = FALSE,
            alignNames = alignNames
        )
    }
    regions <- .qeh_oneHotRegions()
    swapped <- .qep_renameIds(regions, "r1", 42, "chr1:19100:G:A")
    for (x in list(regions, swapped)) {
        expect_equal(
            .qep_estimates(run(x, TRUE)),
            .qep_estimates(run(x, FALSE))
        )
    }
    expect_false(isTRUE(all.equal(
        .qep_estimates(run(regions, TRUE)),
        .qep_estimates(run(swapped, TRUE))
    )))
})

test_that("qtlEnrichment: matches a tuple's ids allele-aware in one call", {
    calls <- 0L
    local_mocked_bindings(
        .matchVariantTables = function(a, b, ...) {
            calls <<- calls + 1L
            list(idxA = integer(0), idxB = integer(0), sign = numeric(0))
        },
        .package = "pecotmr"
    )
    # Every name spelled as a GWAS id: nothing is left to match.
    .qep_oneHotRun()
    expect_identical(calls, 0L)
    # Swapped ids in all three regions: one call for the tuple.
    renamed <- .qep_renameIds(.qeh_oneHotRegions(), "r1", 42, "chr1:19100:G:A")
    renamed <- .qep_renameIds(renamed, "r2", 26, "chr1:19500:G:A")
    renamed <- .qep_renameIds(renamed, "r3", 5, "chr1:500:G:A")
    .qep_oneHotRun(renamed)
    expect_identical(calls, 1L)
})

test_that("qtlEnrichment: refuses gwasPip naming a variant twice", {
    gwasPip <- .qeh_gwasPip()
    expect_error(
        qtlEnrichment(
            c(gwasPip, gwasPip[1]),
            .qeh_oneHotRegions(),
            numGwas = 5000,
            piQtl = 0.01,
            verbose = FALSE
        ),
        "names variant 'chr1:100:A:G' more than once"
    )
})

# ===========================================================================
# qtlEnrichment(): verbose messages
# ===========================================================================

test_that("qtlEnrichment: verbose=TRUE emits 'Estimated piGwas'/'Estimated piQtl' messages", {
    # numGwas = NULL and piQtl = NULL both trigger data-estimation paths;
    # with verbose = TRUE each emits a message.
    fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
    msgs <- testthat::capture_messages(
        suppressWarnings(qtlEnrichment(
            gwasPip = fx$gwasPip,
            susieQtlRegions = fx$susieQtlRegions,
            impN = 5,
            BPPARAM = BiocParallel::SerialParam(),
            verbose = TRUE
        ))
    )
    allMsgs <- paste(msgs, collapse = "\n")
    expect_match(allMsgs, "Estimated piGwas")
    expect_match(allMsgs, "Estimated piQtl")
})

# ===========================================================================
# .enrFlattenEnrichment(): shape coercion + NA fallbacks
# ===========================================================================

test_that(".enrFlattenEnrichment: reads the estimator's own field names", {
    out <- pecotmr:::.enrFlattenEnrichment(.qep_mockEnrichment(2.0)())
    # enrichment is expm1 of the log-odds, so colocPipeline's
    # p12 * (1 + enrichment) is the enloc prior p12 * exp(a1).
    expect_equal(out$enrichment, 2.0)
    expect_equal(out$enrichmentLogOdds, log(3))
    expect_equal(out$enrichmentSe, 0.1)
    expect_equal(out$enrichmentSeNoShrinkage, 0.2)
    expect_equal(out$intercept, -6.5)
    expect_equal(out$colocP12, 5e-6)
    expect_equal(out$effectiveMiRounds, 25)
})

test_that(".enrFlattenEnrichment: a zero log-odds leaves p12 unscaled", {
    mock <- .qep_mockEnrichment(2.0)()
    mock[["Enrichment (w/ shrinkage)"]] <- 0
    expect_equal(pecotmr:::.enrFlattenEnrichment(mock)$enrichment, 0)
})

test_that(".enrFlattenEnrichment: an unrecognized shape warns, not NAs", {
    # Silent all-NA here is exactly the failure that hid a field-name
    # mismatch between this pipeline and the C++ estimator.
    expect_warning(
        out <- pecotmr:::.enrFlattenEnrichment(1.5),
        "returned no 'Enrichment \\(w/ shrinkage\\)' field"
    )
    expect_equal(out, pecotmr:::.enrNaEnrichment())
    expect_warning(
        out <- pecotmr:::.enrFlattenEnrichment(list(enrichment = 2.0)),
        "returned no 'Enrichment"
    )
    expect_true(is.na(out$enrichment))
})


# ===========================================================================
# .enrBuildGwasPipVector(): empty-index + length-mismatch skip branches
# ===========================================================================

test_that(".enrBuildGwasPipVector: unknown study yields numeric(0)", {
    gfmr <- .qep_makeGwasFmr()
    expect_identical(
        pecotmr:::.enrBuildGwasPipVector(gfmr, .qep_gwasIdent("GHOST")),
        numeric(0)
    )
})

test_that(".enrBuildGwasPipVector: skips a block whose pip length disagrees with variantIds", {
    # Unnamed pip of length 2 against 3 variant ids -> ids come from
    # variantIds and length(ids) != length(pip) -> the block is skipped,
    # leaving no pieces -> numeric(0).
    badEntry <- FineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        susieFit = list(pip = c(0.1, 0.2)),
        topLoci = data.frame(
            variant_id = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
            pip = c(0.5, 0.3, 0.2),
            stringsAsFactors = FALSE
        )
    )
    g <- GwasFineMappingResult(
        studyName = "G1",
        method = "susie",
        entry = list(badEntry),
        ldReference = .qep_makeHandle()
    )
    expect_identical(
        pecotmr:::.enrBuildGwasPipVector(g, .qep_gwasIdent("G1")),
        numeric(0)
    )
})


# ===========================================================================
# .enrBuildQtlRegionsList(): incomplete-fit / no-prior / unnamed-pip branches
# ===========================================================================

test_that(".enrBuildQtlRegionsList: skips an entry whose fit lacks alpha/pip", {
    badEntry <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.1,
            stringsAsFactors = FALSE
        )
    )
    qfmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(badEntry),
        ldReference = NULL
    )
    expect_equal(
        length(pecotmr:::.enrBuildQtlRegionsList(
            qfmr,
            .qep_qtlIdent("Q1", "c1")
        )),
        0L
    )
})

test_that(".enrBuildQtlRegionsList: names an unnamed pip from the entry's variant ids", {
    unnamedPipFit <- list(
        alpha = matrix(c(0.6, 0.4), nrow = 1),
        pip = c(0.6, 0.4), # unnamed -> names assigned from variantIds
        V = 0.1
    )
    entry <- FineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
        susieFit = unnamedPipFit,
        topLoci = data.frame(
            variant_id = c("chr1:100:A:G", "chr1:200:A:G"),
            pip = c(0.6, 0.4),
            stringsAsFactors = FALSE
        )
    )
    qfmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(entry),
        ldReference = NULL
    )
    out <- pecotmr:::.enrBuildQtlRegionsList(qfmr, .qep_qtlIdent("Q1", "c1"))
    expect_equal(length(out), 1L)
    expect_equal(names(out[[1L]]$pip), c("chr1:100:A:G", "chr1:200:A:G"))
    expect_equal(out[[1L]]$prior_variance, 0.1)
})


# ===========================================================================
# qtlEnrichmentPipeline(): no-triples error, alignTuple cache reuse,
# and the per-tuple "no usable QTL regions" skip
# ===========================================================================

test_that("qtlEnrichmentPipeline: empty QTL collection errors with the no-triples message", {
    gfmr <- .qep_makeGwasFmr()
    qfmrEmpty <- QtlFineMappingResult(
        studyName = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        entry = list(),
        ldReference = NULL
    )
    expect_error(
        qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = qfmrEmpty
        ),
        "pairs to compute"
    )
})

test_that("qtlEnrichmentPipeline: alignTuple cache is reused across GWAS studies", {
    # Two GWAS studies -> each QTL tuple's alignment is computed for the
    # first study (cache miss) and served from cache for the second
    # (alignTuple early-return branch).
    e1 <- .qep_makeFmEntry(
        variant_ids = sprintf("chr1:%d:A:G", 100L * (1:3)),
        pip = c(0.5, 0.2, 0.1)
    )
    e2 <- .qep_makeFmEntry(
        variant_ids = sprintf("chr1:%d:A:G", 100L * (1:3)),
        pip = c(0.4, 0.3, 0.2)
    )
    gfmr <- GwasFineMappingResult(
        studyName = c("G1", "G2"),
        method = c("susie", "susie"),
        entry = list(e1, e2),
        ldReference = .qep_makeHandle()
    )
    qfmr <- .qep_makeQtlFmr()
    local_mocked_bindings(
        .qtlEnrichmentCoded = .qep_mockEnrichment(2.0),
        .package = "pecotmr"
    )
    out <- qtlEnrichmentPipeline(
        gwasFineMappingResult = gfmr,
        qtlFineMappingResult = qfmr
    )
    expect_equal(nrow(out), 2L) # 2 GWAS studies * 1 QTL tuple
    expect_setequal(out$gwasStudy, c("G1", "G2"))
})

test_that("qtlEnrichmentPipeline: a tuple with no usable QTL regions warns and is skipped", {
    gfmr <- .qep_makeGwasFmr()
    emptyQtlEntry <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.1,
            stringsAsFactors = FALSE
        )
    )
    qfmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(emptyQtlEntry),
        ldReference = .qep_makeHandle()
    )
    expect_warning(
        out <- qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = qfmr
        ),
        "no usable regions"
    )
    expect_equal(nrow(out), 0L)
})

test_that(".enrRunEnrichment turns a stored alignment failure into a warning", {
    cnd <- simpleCondition("boom")
    class(cnd) <- c("simpleError", "error", "condition")
    # The alignment failure was captured earlier and stored in place of the
    # aligned data; re-signalling it here routes it through the same handler
    # a live failure would take.
    expect_warning(
        out <- pecotmr:::.enrRunEnrichment(
            1L,
            outcome = NULL,
            outcomeAt = NULL,
            k = 1L,
            prep = list(aligned = list(cnd), annotations = list(NULL)),
            qtlEnrichmentParam = QtlEnrichmentParam(),
            BPPARAM = BiocParallel::SerialParam(),
            seed = NULL,
            verbose = FALSE,
            gwasFineMappingResult = NULL,
            gwasTuples = data.frame(study = "g"),
            qtlFineMappingResult = NULL,
            qtlTuples = data.frame(study = "q")
        ),
        "qtlEnrichment failed"
    )
    expect_null(out)
})

test_that(".enrBuildQtlRegionsList reads prior variance under either name", {
    pip <- c(v1 = 0.5, v2 = 0.5)
    run <- function(fit) {
        local_mocked_bindings(
            .enrMatchRows = function(qtlFmr, ident) 1L,
            .fmrFitsAt = function(x, idx) list(fit),
            .package = "pecotmr"
        )
        pecotmr:::.enrBuildQtlRegionsList("fmr", "ident")
    }
    base <- list(alpha = matrix(0.5, 1L, 2L), pip = pip)
    # `V` is susie's name for it; `prior_variance` is the older one. Either
    # is accepted, and a fit carrying neither is skipped rather than kept
    # with a NULL prior.
    expect_equal(run(c(base, list(V = 0.7)))[[1L]]$prior_variance, 0.7)
    expect_equal(
        run(c(base, list(prior_variance = 0.3)))[[1L]]$prior_variance,
        0.3
    )
    expect_length(run(base), 0L)
})

# =============================================================================
# QtlEnrichmentParam()
# =============================================================================

test_that("QtlEnrichmentParam carries the estimator's tunable settings", {
    a <- QtlEnrichmentParam()
    expect_s4_class(a, "QtlEnrichmentParam")
    # NULL defaults mean "estimate from the data", which qtlEnrichment says by
    # the argument being absent rather than explicitly NULL.
    expect_equal(
        names(a),
        c("lambda", "impN", "doubleShrinkage", "besselCorrection")
    )
    expect_equal(a$lambda, 1.0)
    expect_equal(a$impN, 25)
    expect_false(a$doubleShrinkage)
    expect_true(a$besselCorrection)
})

test_that("QtlEnrichmentParam forwards numGwas and piQtl once set", {
    a <- QtlEnrichmentParam(numGwas = 5000, piQtl = 0.01)
    expect_equal(a$numGwas, 5000)
    expect_equal(a$piQtl, 0.01)
})

test_that("QtlEnrichmentParam rejects a name the estimator has no formal for", {
    # No `...`: R's own argument matching is the check for a bundle pecotmr
    # owns outright.
    expect_error(QtlEnrichmentParam(shrinkage = 2), "unused argument")
    # alignNames is the pipeline's invariant, not a user setting.
    expect_error(QtlEnrichmentParam(alignNames = TRUE), "unused argument")
})

test_that("every QtlEnrichmentParam field is a qtlEnrichment formal", {
    # The bundle is spliced straight into qtlEnrichment, so a field it has no
    # formal for would be an error at a call site far from the constructor.
    expect_true(all(
        names(formals(QtlEnrichmentParam)) %in% names(formals(qtlEnrichment))
    ))
})

test_that("qtlEnrichmentPipeline refuses a bare list for qtlEnrichmentParam", {
    expect_error(
        qtlEnrichmentPipeline(
            gwasFineMappingResult = NULL,
            qtlFineMappingResult = NULL,
            qtlEnrichmentParam = list(impN = 5)
        ),
        "must be built with QtlEnrichmentParam"
    )
})

test_that("qtlEnrichmentPipeline refuses overlapping rows on two references", {
    ids <- sprintf("chr1:%d:A:G", 100L * (1:9))
    ldm <- function(rho) {
        LdMatrix(`dimnames<-`(rho^abs(outer(1:9, 1:9, "-")), list(ids, ids)))
    }
    gfmr <- pecotmr:::.withRowLdRefs(.qep_makeGwasFmr(), ldm(0.2))
    qfmr <- pecotmr:::.withRowLdRefs(.qep_makeQtlFmr(), ldm(0.8))
    expect_error(
        pecotmr:::.enrValidateInputs(gfmr, qfmr),
        "disagree on the LD"
    )
})

test_that("qtlEnrichmentPipeline keeps each fine-mapping method's fits apart", {
    # The same GWAS blocks fine-mapped by two methods are two outcome PIP
    # vectors, and the same QTL regions by two methods two annotations.
    gwas <- .qep_makeGwasFmr(n_blocks = 2L)
    gwas2 <- gwas
    S4Vectors::mcols(gwas2)$method <- "ser"
    gfmr <- combineFineMappingResults(gwas, gwas2)
    qtl <- .qep_makeQtlFmr()
    qtl2 <- qtl
    S4Vectors::mcols(qtl2)$method <- "ser"
    qfmr <- combineFineMappingResults(qtl, qtl2)
    seen <- new.env()
    seen$regions <- integer(0)
    local_mocked_bindings(
        .qtlEnrichmentCoded = function(gwasPip, susieQtlRegions, ...) {
            seen$regions <- c(seen$regions, length(susieQtlRegions))
            .qep_mockEnrichment(2.0)(gwasPip, susieQtlRegions, ...)
        },
        .package = "pecotmr"
    )
    out <- qtlEnrichmentPipeline(
        gwasFineMappingResult = gfmr,
        qtlFineMappingResult = qfmr,
        verbose = FALSE
    )
    expect_equal(nrow(out), 4L)
    expect_setequal(
        str_c(out$gwasMethod, out$qtlMethod, sep = "/"),
        c("susie/susie", "susie/ser", "ser/susie", "ser/ser")
    )
    # Each annotation is one method's one region, counted once.
    expect_true(all(seen$regions == 1L))
    # The GWAS outcome is both blocks of its fit.
    ident <- list(
        study = "G1",
        context = NA_character_,
        trait = NA_character_,
        method = "ser"
    )
    expect_length(pecotmr:::.enrBuildGwasPipVector(gfmr, ident), 6L)
})

# ===========================================================================
# A narrowed row's PIPs are named by the variants of its fit
# ===========================================================================

test_that("enrichment names a narrowed row's PIPs by its fit's variants", {
    vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    pip <- c(0.5, 0.3, 0.2)
    fit <- list(pip = set_names(pip, vids), alpha = matrix(pip, 1L), V = 1)
    tl <- data.frame(variant_id = vids, pip = pip)
    row <- FineMappingRow(vids, fit, tl)
    gwas <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(row)
    )
    expect_equal(
        pecotmr:::.enrRowPipVector(1L, endoapply(gwas, utils::head, 2L)),
        set_names(pip, vids)
    )
    qtl <- QtlFineMappingResult(
        studyName = "q1",
        context = "c1",
        trait = "g1",
        method = "susie",
        entry = list(row)
    )
    region <- pecotmr:::.enrRowRegion(1L, endoapply(qtl, utils::head, 2L))
    expect_equal(region$pip, set_names(pip, vids))
    expect_equal(region$prior_variance, 1)
    # A narrowed row whose fit names none has nothing to name them by.
    unnamed <- QtlFineMappingResult(
        studyName = "q1",
        context = "c1",
        trait = "g1",
        method = "susie",
        entry = list(FineMappingRow(
            vids,
            list(pip = pip, alpha = fit$alpha, V = 1),
            tl
        ))
    )
    expect_error(
        pecotmr:::.enrRowRegion(1L, endoapply(unnamed, utils::head, 2L)),
        "does not name its variables"
    )
})

test_that("qtlEnrichmentPipeline warns about narrowed fits on both sides", {
    data(gwasFineMappingExample, qtlFineMappingExample)
    cnd <- .fitCoverageCondition(qtlEnrichmentPipeline(
        gwasFineMappingResult = endoapply(
            gwasFineMappingExample,
            utils::head,
            2000L
        ),
        qtlFineMappingResult = endoapply(
            qtlFineMappingExample,
            utils::head,
            2000L
        ),
        verbose = FALSE
    ))
    expect_s3_class(cnd, "pecotmrFitCoverage")
    expect_match(
        conditionMessage(cnd),
        "qtlEnrichmentPipeline\\(\\): 2 row\\(s\\)"
    )
})
