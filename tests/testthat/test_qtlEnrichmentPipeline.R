context("qtlEnrichmentPipeline")

# ===========================================================================
# Strategy: mock qtlEnrichment so the pipeline runs end-to-end on a
# small fixture, but the heavy mixture-of-enrichment estimator never fires.
# ===========================================================================

# The internal builders take an identity list keyed by the collection's own
# columns; the axes a GWAS collection lacks are NA.
.qep_gwasIdent <- function(
    study,
    context = NA_character_,
    trait = NA_character_
) {
    list(study = study, context = context, trait = trait)
}

.qep_qtlIdent <- function(study, context) {
    list(study = study, context = context)
}

.qep_makeHandle <- function(
    snp_n = 6L,
    n_samples = 30L,
    path = "/tmp/sketch.gds"
) {
    new(
        "GenotypeHandle",
        path = path,
        format = "gds",
        snpInfo = data.frame(
            SNP = sprintf("chr1:%d:A:G", 100L * (seq_len(snp_n))),
            CHR = rep("1", snp_n),
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
    fineMappingRow(variantIds = variant_ids, susieFit = fit, topLoci = tl)
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
        study = studyVec,
        method = methodVec,
        entry = entries,
        ldSketch = if (with_sketch) .qep_makeHandle() else NULL
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
        study = studies,
        context = ctx,
        trait = trs,
        method = methods,
        entry = entries,
        ldSketch = if (with_sketch) .qep_makeHandle() else NULL
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
        study = "G1",
        method = "susie",
        entry = list(fineMappingRow(
            variantIds = ids,
            susieFit = list(
                alpha = matrix(1 / 50, nrow = 1L, ncol = 50L),
                pip = setNames(gwasPip, ids),
                V = 0.5
            ),
            topLoci = data.frame(variant_id = ids, pip = gwasPip)
        )),
        ldSketch = sketch
    )
    qfmr <- QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fineMappingRow(
            variantIds = ids,
            susieFit = list(
                alpha = alpha,
                pip = setNames(colSums(alpha), ids),
                V = c(0.5, 0.3)
            ),
            topLoci = data.frame(variant_id = ids, pip = colSums(alpha))
        )),
        ldSketch = sketch
    )
    suppressWarnings(capture.output(
        out <- qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = qfmr,
            impN = 5,
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
        study = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t2"),
        method = c("susie", "susie"),
        entry = list(
            .qep_makeFmEntry(pip = c(0.9, 0.05, 0.02, 0.02, 0.01)),
            .qep_makeFmEntry(pip = c(0.1, 0.7, 0.1, 0.05, 0.05))
        ),
        ldSketch = .qep_makeHandle()
    )
    local_mocked_bindings(
        qtlEnrichment = .qep_mockEnrichment(2.0),
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
        qtlEnrichment = .qep_mockEnrichment(2.0),
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

test_that("qtlEnrichmentPipeline: a QTL outcome side may carry no ldSketch", {
    # The RSS-derived requirement is a GWAS-side contract; an individual-level
    # QTL outcome has no panel to check.
    local_mocked_bindings(
        qtlEnrichment = .qep_mockEnrichment(2.0),
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
        "must be a GwasFineMappingResult"
    )
})

test_that("qtlEnrichmentPipeline: rejects non-QtlFineMappingResult qtlFmr", {
    gfmr <- .qep_makeGwasFmr()
    expect_error(
        qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = "no"
        ),
        "must be a QtlFineMappingResult"
    )
})

test_that("qtlEnrichmentPipeline: NULL ldSketch on the GWAS side errors", {
    gfmr <- .qep_makeGwasFmr(with_sketch = FALSE)
    qfmr <- .qep_makeQtlFmr()
    expect_error(
        qtlEnrichmentPipeline(
            gwasFineMappingResult = gfmr,
            qtlFineMappingResult = qfmr
        ),
        "must have a non-NULL ldSketch"
    )
})

test_that("qtlEnrichmentPipeline: ldSketch mismatch errors", {
    # Build the QTL with a sketch carrying a different sample set.
    gfmr <- .qep_makeGwasFmr()
    qSketch <- .qep_makeHandle()
    qSketch@sampleIds <- paste0("z", seq_len(getNSamples(qSketch)))
    qfmr <- QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(.qep_makeFmEntry()),
        ldSketch = qSketch
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
        qtlEnrichment = .qep_mockEnrichment(2.0),
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
        study = c("Q1", "Q2"),
        context = c("shared_ctx", "shared_ctx"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2),
        ldSketch = .qep_makeHandle()
    )
    gfmr <- .qep_makeGwasFmr()
    capturedRegions <- list()
    local_mocked_bindings(
        qtlEnrichment = function(gwasPip, susieQtlRegions, ...) {
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
        qtlEnrichment = function(...) stop("synthetic failure"),
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
    emptyEntry <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(), # no pip -> .enrBuildGwasPipVector returns numeric(0)
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.1,
            stringsAsFactors = FALSE
        )
    )
    gfmr <- GwasFineMappingResult(
        study = "G1",
        method = "susie",
        entry = list(emptyEntry),
        ldSketch = .qep_makeHandle()
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
            "qtlStudy",
            "qtlContext",
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
        study = c("G1", "G1"),
        method = c("susie", "susie"),
        entry = list(e1, e2),
        ldSketch = .qep_makeHandle()
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
        study = c("G1", "G1"),
        method = c("susie", "susie"),
        entry = list(e1, e2),
        ldSketch = .qep_makeHandle()
    )
    expect_error(
        pecotmr:::.enrBuildGwasPipVector(g, .qep_gwasIdent("G1")),
        "conflicting PIPs"
    )
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
# These tests deliberately do NOT mock qtlEnrichmentRcpp so the C++
# kernel in src/qtl_enrichment.cpp gets coverage. The wrapper itself
# (R/qtlEnrichmentPipeline.R::qtlEnrichment) is exercised here directly
# rather than via the deprecated `computeQtlEnrichment` shim (which has
# skip_on_covr()).
# ===========================================================================

# Build a small (gwasPip, susieQtlRegions) fixture with a sparse causal
# signal at known indices so the C++ enrichment routine has something
# meaningful to compute.
.qep_makeRealKernelInputs <- function(
    seed = 42,
    nSnps = 50,
    causalIdx = c(5, 20, 35),
    causalPips = c(0.8, 0.6, 0.9),
    L = 2L
) {
    set.seed(seed)
    variantNames <- paste0("1:", seq_len(nSnps), ":A:G")
    gwasPip <- rep(0.01, nSnps)
    gwasPip[causalIdx] <- causalPips
    names(gwasPip) <- variantNames

    alpha <- matrix(1 / nSnps, nrow = L, ncol = nSnps)
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
        numThreads = 1,
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
        numThreads = 1,
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
            numThreads = 1,
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
            numThreads = 1,
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
            piQtl = 0.5,
            numThreads = 1,
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
            numThreads = 1,
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
            numThreads = 1,
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
            numThreads = 1,
            verbose = FALSE
        ),
        "Variant names are missing in susieQtlRegions"
    )
})

test_that("qtlEnrichment: tracks unmatched QTL variants in the output", {
    fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
    # Inject a couple of variant IDs into the QTL fit that don't exist
    # in the GWAS PIP vector.
    newNames <- names(fx$susieQtlRegions$fit1$pip)
    newNames[1:2] <- c("1:9999:A:G", "1:9998:A:G")
    names(fx$susieQtlRegions$fit1$pip) <- newNames
    colnames(fx$susieQtlRegions$fit1$alpha) <- newNames
    res <- qtlEnrichment(
        gwasPip = fx$gwasPip,
        susieQtlRegions = fx$susieQtlRegions,
        numGwas = 3000,
        piQtl = 0.5,
        impN = 5,
        numThreads = 1,
        verbose = FALSE
    )
    expect_true("unused_xqtl_variants" %in% names(res))
})


# ===========================================================================
# qtlEnrichment(): verbose messages + alignNames = FALSE branch
# ===========================================================================

test_that("qtlEnrichment: verbose=TRUE emits 'Estimated piGwas'/'Estimated piQtl' messages", {
    # numGwas = NULL and piQtl = NULL both trigger data-estimation paths;
    # with verbose = TRUE each emits a message (lines 391 and 407).
    fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
    msgs <- testthat::capture_messages(
        suppressWarnings(qtlEnrichment(
            gwasPip = fx$gwasPip,
            susieQtlRegions = fx$susieQtlRegions,
            impN = 5,
            numThreads = 1,
            verbose = TRUE
        ))
    )
    allMsgs <- paste(msgs, collapse = "\n")
    expect_match(allMsgs, "Estimated piGwas")
    expect_match(allMsgs, "Estimated piQtl")
})

test_that("qtlEnrichment: alignNames=FALSE recomputes only the unmatched set", {
    # Exercises the cheap set-membership branch (lines 438-444) used by
    # qtlEnrichmentPipeline after it has already aligned QTL names against
    # the union GWAS panel. Two injected names are absent from gwasPip, so
    # the unmatched-variant accumulation (lines 441-442) runs.
    fx <- .qep_makeRealKernelInputs(nSnps = 30, causalIdx = c(5, 15))
    newNames <- names(fx$susieQtlRegions$fit1$pip)
    newNames[1:2] <- c("1:9999:A:G", "1:9998:A:G")
    names(fx$susieQtlRegions$fit1$pip) <- newNames
    colnames(fx$susieQtlRegions$fit1$alpha) <- newNames
    res <- qtlEnrichment(
        gwasPip = fx$gwasPip,
        susieQtlRegions = fx$susieQtlRegions,
        numGwas = 3000,
        piQtl = 0.5,
        impN = 5,
        numThreads = 1,
        verbose = FALSE,
        alignNames = FALSE
    )
    expect_true("unused_xqtl_variants" %in% names(res))
    expect_true(any(
        c("1:9999:A:G", "1:9998:A:G") %in%
            unlist(res$unused_xqtl_variants)
    ))
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
    # getVariantIds and length(ids) != length(pip) -> the block is skipped,
    # leaving no pieces -> numeric(0).
    badEntry <- fineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        susieFit = list(pip = c(0.1, 0.2)),
        topLoci = data.frame(
            variant_id = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
            pip = c(0.5, 0.3, 0.2),
            stringsAsFactors = FALSE
        )
    )
    g <- GwasFineMappingResult(
        study = "G1",
        method = "susie",
        entry = list(badEntry),
        ldSketch = .qep_makeHandle()
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
    badEntry <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.1,
            stringsAsFactors = FALSE
        )
    )
    qfmr <- QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(badEntry),
        ldSketch = NULL
    )
    expect_equal(
        length(pecotmr:::.enrBuildQtlRegionsList(
            qfmr,
            .qep_qtlIdent("Q1", "c1")
        )),
        0L
    )
})

test_that(".enrBuildQtlRegionsList: skips a fit with no V and no prior_variance", {
    noVfit <- list(
        alpha = matrix(c(0.6, 0.4), nrow = 1),
        pip = setNames(c(0.6, 0.4), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    entry <- fineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
        susieFit = noVfit,
        topLoci = data.frame(
            variant_id = c("chr1:100:A:G", "chr1:200:A:G"),
            pip = c(0.6, 0.4),
            stringsAsFactors = FALSE
        )
    )
    qfmr <- QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(entry),
        ldSketch = NULL
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
        pip = c(0.6, 0.4), # unnamed -> names assigned from getVariantIds
        V = 0.1
    )
    entry <- fineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
        susieFit = unnamedPipFit,
        topLoci = data.frame(
            variant_id = c("chr1:100:A:G", "chr1:200:A:G"),
            pip = c(0.6, 0.4),
            stringsAsFactors = FALSE
        )
    )
    qfmr <- QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(entry),
        ldSketch = NULL
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
        study = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        entry = list(),
        ldSketch = NULL
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
    # (alignTuple early-return branch, line 119).
    e1 <- .qep_makeFmEntry(
        variant_ids = sprintf("chr1:%d:A:G", 100L * (1:3)),
        pip = c(0.5, 0.2, 0.1)
    )
    e2 <- .qep_makeFmEntry(
        variant_ids = sprintf("chr1:%d:A:G", 100L * (1:3)),
        pip = c(0.4, 0.3, 0.2)
    )
    gfmr <- GwasFineMappingResult(
        study = c("G1", "G2"),
        method = c("susie", "susie"),
        entry = list(e1, e2),
        ldSketch = .qep_makeHandle()
    )
    qfmr <- .qep_makeQtlFmr()
    local_mocked_bindings(
        qtlEnrichment = .qep_mockEnrichment(2.0),
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
    emptyQtlEntry <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.1,
            stringsAsFactors = FALSE
        )
    )
    qfmr <- QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(emptyQtlEntry),
        ldSketch = .qep_makeHandle()
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
