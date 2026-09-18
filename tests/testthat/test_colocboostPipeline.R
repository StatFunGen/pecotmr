context("colocboostPipeline")

# ============================================================================
# Legacy colocboostPipeline tests removed for the post-S4-refactor API.
# ============================================================================
#
# This file previously contained ~100 test_that blocks targeting the legacy
# colocboost pipeline API (rssData/ldData list-shape inputs, RegionalData
# fixtures, qcRegionalData driver, colocboostAnalysis adapters, etc.). Every
# block exercised functions and classes that no longer exist in pecotmr:
#
#   Functions removed:
#     - colocboostAnalysis()         (direct colocboost adapter)
#     - qcRegionalData()             (legacy QC driver)
#     - regionDataToIndInput()       (now a .Deprecated() no-op)
#     - regionDataToRssInput()       (now a .Deprecated() no-op)
#     - regionDataToColocboostInput()
#     - rssAnalysisPipeline()        (replaced by fineMappingPipeline)
#     - rssBasicQc()                 (folded into summaryStatsQc)
#     - loadRssData()                (replaced by SumStats constructors)
#     - getrssinput(), getlddata(), getoutliernumber()
#     - colocWrapper(), xqtlEnrichmentWrapper(), colocPostProcessor()
#       (now .Deprecated() no-ops returning NULL)
#     - .runColocboost()             (replaced by internal .cbRun(label, args))
#     - buildLdArgs()                (replaced by internal .cbBuildLdArgs())
#
#   Classes removed:
#     - RegionalData, MultivariateRegionalData
#     - QcResult, AlleleQcResult
#
# The replacement API has a fundamentally different contract:
#
#   colocboostPipeline() is now an S4 generic dispatching on the QTL input
#   class. Signatures live in R/colocboostPipeline.R:
#
#     setMethod("colocboostPipeline", "QtlDataset",         ...)
#     setMethod("colocboostPipeline", "QtlSumStats",        ...)
#     setMethod("colocboostPipeline", "MultiStudyQtlDataset", ...)
#
#   - QTL data is supplied as a QtlDataset / QtlSumStats / MultiStudyQtlDataset
#     (DFrame-based S4 objects with an ldSketch slot for sumstats inputs and
#     getResidualizedGenotypes() / getResidualizedPhenotypes() accessors for
#     individual-level inputs).
#   - GWAS is supplied separately via the gwasSumStats = GwasSumStats(...)
#     argument.
#   - Individual-level QC (MAF / X-variance / sample missingness / event
#     selection) lives on the QtlDataset constructor and is applied lazily
#     by its accessors. There is no separate qcRegionalData() pass.
#   - All summary-statistic QC lives in summaryStatsQc(). The pipeline
#     rejects QtlSumStats / GwasSumStats whose getQcInfo() is empty.
#
# Rewriting the legacy tests in place would require fabricating new
# QtlDataset / QtlSumStats / GwasSumStats / MultiStudyQtlDataset fixtures
# and asserting against a different result shape -- i.e. inventing new
# coverage rather than porting existing coverage. That is out of scope
# for this legacy-cleanup pass.
#
# New tests for the S4 colocboostPipeline() methods, summaryStatsQc(), and
# the QtlDataset / QtlSumStats / GwasSumStats / MultiStudyQtlDataset
# constructors should be added in dedicated files alongside this one
# (e.g. test_colocboost_pipeline_qtl_dataset.R,
# test_colocboost_pipeline_qtl_sumstats.R,
# test_colocboost_pipeline_multitask.R). See:
#
#   - R/colocboostPipeline.R    (the new S4 generic + methods)
#   - R/allClasses.R             (QtlDataset, QtlSumStats, GwasSumStats,
#                                 MultiStudyQtlDataset, FineMappingRow)
#   - R/allMethods.R             (constructors)
#   - R/sumstatsQc.R            (summaryStatsQc, which is now the only
#                                 summary-statistic QC entry point)
#   - tests/testthat/test_sumstatsQc.R  (existing QC coverage)
# ============================================================================

# Sentinel test so the testthat context is non-empty.
test_that("colocboostPipeline is exported as an S4 generic", {
    expect_true("colocboostPipeline" %in% getNamespaceExports("pecotmr"))
    expect_true(methods::isGeneric("colocboostPipeline"))
})


context("colocboostPipeline (S4 dispatch)")

# ===========================================================================
# Strategy
# ----------------------------------------------------------------------------
# colocboost::colocboost is the heavy compute; we mock it to return a stub
# result so the pipeline orchestration runs end-to-end. The helpers
# (.cbBuildLdArgs, .cbMergeSumstatBundles,
# .cbRequireMatchingLdSketches, .cbEmptyResult) are exercised directly.
# ===========================================================================

.cbp_makeHandle <- function(snp_n = 6L, n_samples = 30L, sample_prefix = "s") {
    new(
        "GenotypeHandle",
        path = "/tmp/cb.gds",
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
        sampleIds = paste0(sample_prefix, seq_len(n_samples)),
        pgenPtr = NULL
    )
}

.cbp_mockExtractor <- function(seed = 11, n_samples = 30L) {
    function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(seed)
        panel <- matrix(
            rbinom(n_samples * nrow(getSnpInfo(handle)), 2, 0.3),
            nrow = n_samples,
            ncol = nrow(getSnpInfo(handle)),
            dimnames = list(getSampleIds(handle), getSnpInfo(handle)$SNP)
        )
        sub <- panel[, snpIdx, drop = FALSE]
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", getSnpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = getSnpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = getSnpInfo(handle)$SNP[snpIdx],
            A1 = getSnpInfo(handle)$A1[snpIdx],
            A2 = getSnpInfo(handle)$A2[snpIdx]
        )
        cd <- S4Vectors::DataFrame(
            sampleId = getSampleIds(handle),
            row.names = getSampleIds(handle)
        )
        dosage <- t(sub)
        rownames(dosage) <- getSnpInfo(handle)$SNP[snpIdx]
        colnames(dosage) <- getSampleIds(handle)
        SummarizedExperiment::SummarizedExperiment(
            assays = list(dosage = dosage),
            rowRanges = rr,
            colData = cd
        )
    }
}

.cbp_makeSe <- function(traits = c("ENSG_A", "ENSG_B"), n_samples = 30L) {
    rng <- GenomicRanges::GRanges(
        seqnames = rep("chr1", length(traits)),
        ranges = IRanges::IRanges(
            start = seq(1000L, by = 1000L, length.out = length(traits)),
            width = 500L
        )
    )
    names(rng) <- traits
    set.seed(0)
    expr <- matrix(
        rnorm(length(traits) * n_samples),
        nrow = length(traits),
        ncol = n_samples,
        dimnames = list(traits, paste0("s", seq_len(n_samples)))
    )
    cd <- S4Vectors::DataFrame(
        sex = rep(c(0, 1), length.out = n_samples),
        age = seq_len(n_samples),
        row.names = paste0("s", seq_len(n_samples))
    )
    SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = expr),
        rowRanges = rng,
        colData = cd
    )
}

.cbp_makeQtlDataset <- function(
    contexts = "brain",
    traits = c("ENSG_A", "ENSG_B")
) {
    gh <- .cbp_makeHandle()
    phen <- setNames(
        lapply(contexts, function(.) .cbp_makeSe(traits = traits)),
        contexts
    )
    QtlDataset(
        study = "study1",
        genotypes = gh,
        phenotypes = phen,
        genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0)
    )
}

.cbp_makeQtlSumStats <- function(qc = TRUE) {
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
        study = "Q1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .cbp_makeHandle(),
        qcInfo = if (qc) list(step1 = "ok") else list()
    )
}

.cbp_makeGwasSumStats <- function(qc = TRUE) {
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
        study = "G1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .cbp_makeHandle(),
        qcInfo = if (qc) list(step1 = "ok") else list()
    )
}

# ===========================================================================
# Internal helpers (run-anywhere)
# ===========================================================================

test_that(".cbBuildLdArgs: square matrices route to LD", {
    R1 <- diag(4)
    R2 <- diag(4)
    res <- pecotmr:::.cbBuildLdArgs(list(R1, R2))
    expect_true("LD" %in% names(res))
    expect_false("X_ref" %in% names(res))
})

test_that(".cbBuildLdArgs: non-square matrices route to X_ref", {
    X1 <- matrix(0, 10, 4)
    res <- pecotmr:::.cbBuildLdArgs(list(X1))
    expect_true("X_ref" %in% names(res))
})

test_that(".cbBuildLdArgs: empty list returns empty list", {
    expect_equal(pecotmr:::.cbBuildLdArgs(list()), list())
})

test_that(".cbScreenSpec enforces one metric and threads the right spec", {
    expect_error(pecotmr:::.cbScreenSpec(0.5, 5, 0, 0), "one signal screen")
    expect_equal(
        pecotmr:::.cbScreenSpec(0, 5, 0, 0),
        list(metric = "absZ", cutoff = 5)
    )
    expect_equal(
        pecotmr:::.cbScreenSpec(0, 0, 100, 0),
        list(metric = "bf", cutoff = 100)
    )
    expect_equal(pecotmr:::.cbScreenSpec(0.3, 0, 0, 0), 0.3) # legacy pip scalar
    expect_equal(
        pecotmr:::.cbScreenSpec(c(brain = 0.3), 0, 0, 0),
        c(brain = 0.3)
    ) # context-named pass-through
})

test_that(".cbResolveCutoff passes a screen object through uniformly", {
    sc <- list(metric = "absZ", cutoff = 5)
    expect_identical(pecotmr:::.cbResolveCutoff(sc, "brain"), sc)
    expect_identical(pecotmr:::.cbResolveCutoff(sc, "blood"), sc)
    expect_equal(pecotmr:::.cbResolveCutoff(c(brain = 0.3), "brain"), 0.3)
    expect_equal(pecotmr:::.cbResolveCutoff(c(brain = 0.3), "blood"), 0) # unlisted ctx
    expect_equal(pecotmr:::.cbResolveCutoff(0.5, "any"), 0.5)
})

test_that("colocboostPipeline(QtlDataset): enabling two screen metrics errors", {
    qd <- .cbp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
    # .cbScreenSpec runs before the bundle/engine, so this fails fast.
    expect_error(
        colocboostPipeline(qd, pipCutoffToSkip = 0.5, bfCutoffToSkip = 100),
        "one signal screen"
    )
})

test_that(".cbRequireSumStatsQc: un-QCd input errors", {
    ss <- .cbp_makeQtlSumStats(qc = FALSE)
    expect_error(
        pecotmr:::.cbRequireSumStatsQc(ss, "qtlData"),
        "summaryStatsQc"
    )
})

test_that(".cbRequireSumStatsQc: NULL input is a no-op", {
    expect_silent(pecotmr:::.cbRequireSumStatsQc(NULL, "x"))
})

test_that(".cbRequireSumStatsQc: QCd input passes", {
    ss <- .cbp_makeQtlSumStats(qc = TRUE)
    expect_silent(pecotmr:::.cbRequireSumStatsQc(ss, "x"))
})


test_that(".cbRequireMatchingLdSketches: NULL sides are allowed", {
    expect_silent(pecotmr:::.cbRequireMatchingLdSketches(
        NULL,
        .cbp_makeHandle()
    ))
    expect_silent(pecotmr:::.cbRequireMatchingLdSketches(
        .cbp_makeHandle(),
        NULL
    ))
})

test_that(".cbRequireMatchingLdSketches: accepts differently trimmed panels", {
    # Separate QC of the two sides trims one shared LD reference to two
    # overlapping-but-unequal variant sets; the check reports, not refuses.
    expect_warning(
        expect_null(pecotmr:::.cbRequireMatchingLdSketches(
            .cbp_makeHandle(snp_n = 4L),
            .cbp_makeHandle(snp_n = 5L)
        )),
        "share 4 variant"
    )
})

test_that(".cbRequireMatchingLdSketches: sample-set mismatch errors", {
    expect_error(
        pecotmr:::.cbRequireMatchingLdSketches(
            .cbp_makeHandle(sample_prefix = "a"),
            .cbp_makeHandle(sample_prefix = "b")
        ),
        "different sample sets"
    )
})

test_that(".cbMergeSumstatBundles: empty input gives empty dict", {
    res <- pecotmr:::.cbMergeSumstatBundles(list())
    expect_equal(length(res$sumstat), 0L)
    expect_equal(length(res$LD), 0L)
    expect_equal(nrow(res$dict_sumstatLD), 0L)
})

test_that(".cbMergeSumstatBundles: identical LD matrices are deduplicated", {
    R <- diag(3)
    bundles <- list(
        a = list(sumstat = data.frame(z = 1:3), LD = R),
        b = list(sumstat = data.frame(z = 4:6), LD = R)
    )
    res <- pecotmr:::.cbMergeSumstatBundles(bundles)
    expect_equal(length(res$LD), 1L)
    expect_equal(unique(res$dict_sumstatLD[, 2L]), 1L)
})

test_that(".cbEmptyResult: is the internal per-analysis accumulator", {
    # NOT the pipeline's return schema any more -- colocboostPipeline returns a
    # ColocBoostResult. This list is the accumulator .cbRunVariants fills in and
    # .cbToResultObject then flattens, so the four slots still have to be there
    # for the three analyses and their timings to have somewhere to land.
    res <- pecotmr:::.cbEmptyResult()
    expect_true(all(
        c("xqtl_coloc", "joint_gwas", "separate_gwas", "computing_time") %in%
            names(res)
    ))
})

# ===========================================================================
# colocboostPipeline(QtlDataset)
# ===========================================================================

test_that("colocboostPipeline(QtlDataset): runs xqtl-only ColocBoost with mocked engine", {
    qd <- .cbp_makeQtlDataset()
    capturedArgs <- NULL
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            capturedArgs <<- list(...)
            list(stub = TRUE)
        },
        .package = "colocboost"
    )
    out <- suppressMessages(
        colocboostPipeline(
            qd,
            xqtlColoc = TRUE,
            jointGwas = FALSE,
            separateGwas = FALSE
        )
    )
    # The pipeline now returns a ColocBoostResult, so "the xQTL-only analysis
    # ran" is read off the recorded timing rather than off a retained raw
    # object -- the timing is written whether or not the run produced sets.
    expect_s4_class(out, "ColocBoostResult")
    expect_true(is_in("xqtl_coloc", names(getComputingTime(out)$Analysis)))
    expect_true("X" %in% names(capturedArgs))
})

test_that("colocboostPipeline(QtlDataset): no QTL bundle and no GWAS returns empty result", {
    qd <- .cbp_makeQtlDataset()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    # An impossible traitId narrows the phenotype matrix to zero columns;
    # the bundle ends up NULL and the driver short-circuits.
    out <- suppressMessages(suppressWarnings(
        colocboostPipeline(
            qd,
            contexts = "brain",
            traitId = "ENSG_DOES_NOT_EXIST",
            xqtlColoc = TRUE
        )
    ))
    expect_null(out$xqtl_coloc)
})

# ===========================================================================
# colocboostPipeline(QtlSumStats)
# ===========================================================================

test_that("colocboostPipeline(QtlSumStats): un-QCd input rejected", {
    ss <- .cbp_makeQtlSumStats(qc = FALSE)
    expect_error(
        colocboostPipeline(ss, xqtlColoc = TRUE),
        "summaryStatsQc"
    )
})


# ===========================================================================
# colocboostPipeline(ANY)
# ===========================================================================

test_that("colocboostPipeline(ANY): unsupported input class errors", {
    expect_error(
        colocboostPipeline(matrix(0, 3, 3)),
        "does not accept inputs of class"
    )
})

# ===========================================================================
# Driver: jointGwas / separateGwas paths via QtlSumStats + GwasSumStats
# ===========================================================================

test_that("colocboostPipeline: jointGwas merges qtl + gwas sumstats and runs once", {
    ss <- .cbp_makeQtlSumStats()
    gs <- .cbp_makeGwasSumStats()
    capturedArgs <- NULL
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            capturedArgs <<- list(...)
            list(jointly_run = TRUE)
        },
        .package = "colocboost"
    )
    out <- suppressMessages(
        colocboostPipeline(
            ss,
            gwasSumStats = gs,
            xqtlColoc = FALSE,
            jointGwas = TRUE,
            separateGwas = FALSE
        )
    )
    expect_true(is_in("joint_gwas", names(getComputingTime(out)$Analysis)))
})

test_that("colocboostPipeline: separateGwas runs once per merged sumstat study", {
    ss <- .cbp_makeQtlSumStats()
    gs <- .cbp_makeGwasSumStats()
    callCount <- 0
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            callCount <<- callCount + 1L
            list(round = callCount)
        },
        .package = "colocboost"
    )
    out <- suppressMessages(
        colocboostPipeline(
            ss,
            gwasSumStats = gs,
            xqtlColoc = FALSE,
            jointGwas = FALSE,
            separateGwas = TRUE
        )
    )
    # Driver merges QTL + GWAS sumstats into a single bundle and the
    # separate-loop iterates over every merged study label (Q1:c1:t1 + G1).
    expect_equal(callCount, 2L)
    expect_true(
        is_in("separate_gwas", names(getComputingTime(out)$Analysis))
    )
})

test_that("colocboostPipeline: no analysis flag set emits a message and returns empty", {
    ss <- .cbp_makeQtlSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- suppressMessages(
        colocboostPipeline(
            ss,
            xqtlColoc = FALSE,
            jointGwas = FALSE,
            separateGwas = FALSE
        )
    )
    expect_null(out$xqtl_coloc)
    expect_null(out$joint_gwas)
    expect_null(out$separate_gwas)
})

test_that("colocboostPipeline: GWAS ldSketch mismatch errors during the driver", {
    ss <- .cbp_makeQtlSumStats()
    # Build a GwasSumStats whose ldSketch has a different sample set.
    gh_diff <- .cbp_makeHandle(sample_prefix = "z")
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
    gs <- GwasSumStats(
        study = "G1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = gh_diff,
        qcInfo = list(step1 = "ok")
    )
    # .cbQtlSumStatsBundle reads the qtl sketch via extractBlockGenotypes
    # before the LD-sketch mismatch check fires further down the driver, so
    # mock the extractor here too.
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    expect_error(
        suppressMessages(
            colocboostPipeline(
                ss,
                gwasSumStats = gs,
                xqtlColoc = FALSE,
                jointGwas = TRUE
            )
        ),
        "different sample sets"
    )
})

# ===========================================================================
# GWAS case/control: optional nCase/nControl columns + effective-N wiring
# ===========================================================================

test_that("GwasSumStats: nCase/nControl are optional columns (absent by default)", {
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = 5L),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:5)),
        A1 = "A",
        A2 = "G",
        Z = rnorm(5),
        N = rep(1000L, 5)
    )
    base <- list(
        study = "G1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .cbp_makeHandle(),
        qcInfo = list(ok = 1)
    )
    g0 <- do.call(GwasSumStats, base)
    expect_false(any(
        c("nCase", "nControl") %in% colnames(S4Vectors::mcols(g0))
    ))
    g1 <- do.call(GwasSumStats, c(base, list(nCase = 500, nControl = 1500)))
    expect_true(all(c("nCase", "nControl") %in% colnames(S4Vectors::mcols(g1))))
    expect_equal(g1$nCase, 500)
    expect_equal(g1$nControl, 1500)
})

test_that("colocboost GWAS bundle: effective N for case/control, per-variant N otherwise", {
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = 5L),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:5)),
        A1 = "A",
        A2 = "G",
        Z = rnorm(5),
        N = rep(1000L, 5)
    )
    base <- list(
        study = "G1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .cbp_makeHandle(),
        qcInfo = list(ok = 1)
    )
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    # case/control -> effective N = 4 / (1/500 + 1/1500) = 1500
    gcc <- do.call(GwasSumStats, c(base, list(nCase = 500, nControl = 1500)))
    bcc <- pecotmr:::.cbGwasSumStatsBundle(gcc)
    expect_true(all(bcc[["G1"]]$sumstat$n == 4 / (1 / 500 + 1 / 1500)))
    # quantitative (no nCase/nControl) -> per-variant N (1000)
    bq <- pecotmr:::.cbGwasSumStatsBundle(do.call(GwasSumStats, base))
    expect_true(all(bq[["G1"]]$sumstat$n == 1000L))
})

# ===========================================================================
# pipCutoffToSkip: per-context single-trait (L=1 SuSiE) outcome skip
# ===========================================================================

test_that(".cbPipSkipOutcomes: keeps signal outcomes, drops noise, honours cutoff", {
    skip_if_not_installed("susieR")
    set.seed(1)
    n <- 200L
    p <- 20L
    X <- matrix(
        rbinom(n * p, 2, 0.3),
        n,
        p,
        dimnames = list(paste0("s", 1:n), sprintf("chr1:%d:A:G", 100L * (1:p)))
    )
    Y <- cbind(
        sig = X[, 1] * 1.5 + rnorm(n, sd = 0.3), # strong signal at v1
        noise = rnorm(n)
    ) # null
    # cutoff 0 -> no-op
    expect_identical(pecotmr:::.cbPipSkipOutcomes(X, Y, 0), Y)
    # cutoff 0.5 -> keep the signal outcome, drop the noise outcome
    kept <- pecotmr:::.cbPipSkipOutcomes(X, Y, 0.5)
    expect_equal(colnames(kept), "sig")
    # all-noise -> NULL (whole context would be skipped)
    Yn <- cbind(n1 = rnorm(n), n2 = rnorm(n))
    expect_null(pecotmr:::.cbPipSkipOutcomes(X, Yn, 0.5))
})

test_that(".cbResolveCutoff: scalar applies to all; named vector is per-context", {
    expect_equal(pecotmr:::.cbResolveCutoff(0.5, "brain"), 0.5)
    expect_equal(
        pecotmr:::.cbResolveCutoff(c(brain = 0.3, blood = 0.7), "blood"),
        0.7
    )
    expect_equal(pecotmr:::.cbResolveCutoff(c(brain = 0.3), "missing"), 0)
    expect_equal(pecotmr:::.cbResolveCutoff(NULL, "brain"), 0)
})

# ===========================================================================
# Additional coverage: MultiStudy method, engine-failure path, multi-context
# bundle building, and sumstat-bundle helper early returns.
# ===========================================================================

.cbp_makeMultiStudy <- function() {
    MultiStudyQtlDataset(
        qtlDatasets = list(
            study1 = .cbp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
        ),
        sumStats = .cbp_makeQtlSumStats()
    )
}

test_that("colocboostPipeline(MultiStudyQtlDataset): combines per-study bundles + embedded sumstats", {
    mt <- .cbp_makeMultiStudy()
    capturedArgs <- NULL
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            capturedArgs <<- list(...)
            list(stub = TRUE)
        },
        .package = "colocboost"
    )
    out <- suppressMessages(colocboostPipeline(
        mt,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE
    ))
    expect_s4_class(out, "ColocBoostResult")
    expect_true(is_in("xqtl_coloc", names(getComputingTime(out)$Analysis)))
    # The individual study's outcome is prefixed "study1:" in the combined bundle.
    expect_true(any(grepl("study1:", names(capturedArgs$Y))))
})

test_that(".cbRun: an engine failure is caught -> message + NULL", {
    qd <- .cbp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) stop("engine boom"),
        .package = "colocboost"
    )
    out <- suppressMessages(colocboostPipeline(qd, xqtlColoc = TRUE))
    expect_null(out$xqtl_coloc) # .cbRun caught (123-124)
})

test_that(".cbIndividualBundle: multi-context bundle names + prefixes outcomes", {
    qd <- .cbp_makeQtlDataset(contexts = c("brain", "liver"), traits = "ENSG_A")
    capturedArgs <- NULL
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            capturedArgs <<- list(...)
            list(stub = TRUE)
        },
        .package = "colocboost"
    )
    out <- suppressMessages(colocboostPipeline(qd, xqtlColoc = TRUE))
    expect_true(is_in("xqtl_coloc", names(getComputingTime(out)$Analysis)))
    # Two contexts -> two context-prefixed outcomes (covers the xMatch + naming).
    expect_gte(length(capturedArgs$Y), 2L)
})

test_that("colocboostPipeline(QtlDataset): unknown context errors", {
    qd <- .cbp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    expect_error(
        colocboostPipeline(qd, contexts = "ghost", xqtlColoc = TRUE),
        "Unknown context"
    ) # 209
})

test_that("colocboostPipeline(QtlDataset): pipCutoffToSkip dropping every outcome -> empty", {
    qd <- .cbp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- suppressMessages(suppressWarnings(
        colocboostPipeline(qd, xqtlColoc = TRUE, pipCutoffToSkip = 0.9999)
    ))
    expect_null(out$xqtl_coloc) # 249-253 skip -> empty
})

test_that(".cbPipSkipOutcomes: an outcome with < 2 observations is skipped", {
    set.seed(4)
    X <- matrix(
        rnorm(60),
        30,
        2,
        dimnames = list(paste0("s", 1:30), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    Y <- cbind(a = c(1, rep(NA, 29)), b = rnorm(30)) # col a: 1 obs (< 2)
    res <- pecotmr:::.cbPipSkipOutcomes(X, Y, 0.5)
    expect_true(is.null(res) || is.matrix(res)) # col a -> next (180)
})

test_that(".cbSumstatPair: NULL / empty df -> NULL", {
    expect_null(pecotmr:::.cbSumstatPair(NULL, .cbp_makeHandle())) # 317
    expect_null(pecotmr:::.cbSumstatPair(data.frame(), .cbp_makeHandle()))
})

test_that(".cbQtlSumStatsBundle: NULL / context / trait filters and empty result", {
    ss <- .cbp_makeQtlSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    expect_equal(pecotmr:::.cbQtlSumStatsBundle(NULL), list()) # 359
    expect_length(pecotmr:::.cbQtlSumStatsBundle(ss, contexts = "c1"), 1L) # 363
    expect_length(pecotmr:::.cbQtlSumStatsBundle(ss, traitId = "t1"), 1L) # 366
    expect_equal(pecotmr:::.cbQtlSumStatsBundle(ss, contexts = "ghost"), list()) # 368
})

test_that(".cbGwasSumStatsBundle: NULL -> empty list", {
    expect_equal(pecotmr:::.cbGwasSumStatsBundle(NULL), list()) # 390
})

test_that(".cbSumstatPair: varY attaches var_y; NA variant ids fall back to chr:pos", {
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    h <- .cbp_makeHandle() # panel SNPs v1..v6
    df <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        z = c(1, -1, 0.5),
        N = rep(1000, 3),
        stringsAsFactors = FALSE
    )
    pair <- pecotmr:::.cbSumstatPair(df, h, varY = 0.7)
    expect_true("var_y" %in% names(pair$sumstat)) # 350
    expect_equal(unique(pair$sumstat$var_y), 0.7)
    # NA variant_id -> formatVariantId fallback (322-323); no panel overlap -> NULL (330)
    dfNA <- data.frame(
        variant_id = NA_character_,
        chrom = "chr1",
        pos = 999999L,
        A2 = "G",
        A1 = "A",
        z = 1,
        N = 1000,
        stringsAsFactors = FALSE
    )
    expect_null(pecotmr:::.cbSumstatPair(dfNA, h))
})

test_that(".cbSumstatPair canonicalizes variant ids to chr-prefixed for name alignment", {
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    h <- .cbp_makeHandle(snp_n = 3L)
    h@snpInfo$SNP <- c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G") # chr-prefixed panel
    # sumstat carries the same variants without the "chr" prefix; the pipeline
    # should canonicalize them so the sumstat / LD ids align by name with other
    # sources (previously the returned ids kept the caller's convention).
    df <- data.frame(
        variant_id = c("1:100:A:G", "1:200:A:G", "1:300:A:G"),
        z = c(1, -1, 0.5),
        N = rep(1000, 3),
        stringsAsFactors = FALSE
    )
    pair <- pecotmr:::.cbSumstatPair(df, h)
    expect_equal(
        pair$sumstat$variant,
        c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G")
    )
    expect_identical(colnames(pair$LD), pair$sumstat$variant)
})

test_that(".cbFlipPairToCanonical flips z + LD for swapped variants, relabels to canonical", {
    ss <- data.frame(
        z = c(2, -1),
        n = c(1000, 1000),
        variant = c("chr1:100:A:G", "chr1:200:C:T"),
        stringsAsFactors = FALSE
    )
    LD <- matrix(c(1, 0.5, 0.5, 1), 2, dimnames = list(ss$variant, ss$variant))
    p <- list(sumstat = ss, LD = LD, variantIds = ss$variant)
    canonical <- c("chr1:100:G:A", "chr1:200:C:T") # variant 1 swapped, 2 identical
    out <- pecotmr:::.cbFlipPairToCanonical(p, canonical)
    expect_equal(out$sumstat$variant, c("chr1:100:G:A", "chr1:200:C:T"))
    expect_equal(out$sumstat$z, c(-2, -1)) # v1 z flipped; v2 unchanged
    expect_equal(out$LD["chr1:100:G:A", "chr1:200:C:T"], -0.5) # one endpoint flipped
    expect_equal(unname(diag(out$LD)), c(1, 1)) # diagonal preserved
})

test_that(".cbFlipMatrixToCanonical negates residualized dosage for swapped columns", {
    X <- matrix(
        c(1, -1, 2, 0.5, -0.5, 1),
        nrow = 3,
        dimnames = list(paste0("s", 1:3), c("chr1:100:A:G", "chr1:200:C:T"))
    )
    canonical <- c("chr1:100:G:A", "chr1:200:C:T") # col 1 swapped, col 2 identical
    out <- pecotmr:::.cbFlipMatrixToCanonical(X, canonical)
    expect_equal(colnames(out), c("chr1:100:G:A", "chr1:200:C:T"))
    expect_equal(unname(out[, "chr1:100:G:A"]), -c(1, -1, 2)) # negated
    expect_equal(unname(out[, "chr1:200:C:T"]), c(0.5, -0.5, 1))
})

test_that(".cbHarmonizeAlleles aligns opposite-coded sumstats to one canonical (invariance)", {
    mkPair <- function(v, z) {
        ss <- data.frame(z = z, n = 1000, variant = v, stringsAsFactors = FALSE)
        list(
            sumstat = ss,
            LD = matrix(1, 1, 1, dimnames = list(v, v)),
            variantIds = v
        )
    }
    # Same locus, opposite ref/alt coding across two sumstats, same underlying z.
    pairs <- list(A = mkPair("chr1:100:A:G", 3), B = mkPair("chr1:100:G:A", 3))
    h <- pecotmr:::.cbHarmonizeAlleles(NULL, pairs)
    expect_equal(h$pairs$A$sumstat$variant, "chr1:100:A:G") # first-seen = canonical
    expect_equal(h$pairs$B$sumstat$variant, "chr1:100:A:G") # B relabeled to it
    expect_equal(h$pairs$A$sumstat$z, 3) # A already canonical
    expect_equal(h$pairs$B$sumstat$z, -3) # B flipped to match
})

test_that("colocboostPipeline(MultiStudyQtlDataset): a study with no usable bundle is skipped", {
    mt <- .cbp_makeMultiStudy() # qd (study1, ENSG_A) + ss (Q1, t1)
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) list(stub = TRUE),
        .package = "colocboost"
    )
    # traitId="t1" matches the embedded sumstats but not the QtlDataset -> its
    # per-study bundle is NULL and skipped (684); the sumstat side still runs.
    out <- suppressMessages(suppressWarnings(
        colocboostPipeline(mt, traitId = "t1", xqtlColoc = TRUE)
    ))
    expect_s4_class(out, "ColocBoostResult")
})


# ---------------------------------------------------------------------------
# RSS panel filters.
#
# colocboostPipeline never runs summary-statistic QC itself -- that lives in
# summaryStatsQc() -- so these cutoffs act at analysis time, against the LD
# reference panel, on both the QTL and the GWAS sumstat sides.
# ---------------------------------------------------------------------------

# @noRd
.cbf_qcd <- function() {
    data(qtlSumStatsExample, envir = environment())
    suppressMessages(summaryStatsQc(qtlSumStatsExample))
}

# @noRd
.cbf_n <- function(ss, ...) {
    b <- suppressMessages(.cbQtlSumStatsBundle(
        ss,
        cutoffs = .panelCutoffs(list(...))
    ))
    if (length(b) == 0L) 0L else length(b[[1L]]$variantIds)
}

test_that("colocboost sumstat bundle filters nothing by default", {
    ss <- .cbf_qcd()
    expect_equal(.cbf_n(ss), sum(lengths(ss)))
})

test_that("colocboost sumstat bundle drops panel-rare variants", {
    ss <- .cbf_qcd()
    full <- .cbf_n(ss)
    loose <- .cbf_n(ss, mafCutoff = 0.05)
    tight <- .cbf_n(ss, mafCutoff = 0.2)
    expect_lt(loose, full)
    expect_lt(tight, loose)
})

test_that("colocboost RSS cutoffs match .panelVariantFilter", {
    ss <- .cbf_qcd()
    ids <- normalizeVariantId(
        getSumStatsDf(
            ss,
            study = ss$study[[1L]],
            context = ss$context[[1L]],
            trait = ss$trait[[1L]],
            require = "Z"
        )$variant_id
    )
    for (cut in c(0.05, 0.2)) {
        expect_equal(
            .cbf_n(ss, mafCutoff = cut),
            length(.panelVariantFilter(
                getLdSketch(ss),
                ids,
                mafCutoff = cut
            )),
            label = str_c("mafCutoff ", cut)
        )
    }
})

test_that("colocboost RSS treats MAC as a MAF equivalent", {
    ss <- .cbf_qcd()
    nSamp <- ncol(getLdSketch(ss))
    expect_equal(
        .cbf_n(ss, macCutoff = 0.1 * 2 * nSamp),
        .cbf_n(ss, mafCutoff = 0.1)
    )
})

test_that("colocboost RSS honours a missingness cutoff", {
    ss <- .cbf_qcd()
    expect_lt(.cbf_n(ss, imissCutoff = 0), sum(lengths(ss)))
    expect_equal(.cbf_n(ss, imissCutoff = 1), sum(lengths(ss)))
})

test_that(".cbSumstatPair keeps sumstat rows aligned to the LD matrix", {
    # The row filter indexes df by the FULL-length variant id vector, so the
    # filtered set has to be passed to the LD build separately rather than by
    # shrinking that vector.
    ss <- .cbf_qcd()
    df <- getSumStatsDf(
        ss,
        study = ss$study[[1L]],
        context = ss$context[[1L]],
        trait = ss$trait[[1L]],
        require = "Z"
    )
    sketch <- getLdSketch(ss)
    pair <- suppressMessages(.cbSumstatPair(
        df = df,
        ldSketch = sketch,
        cutoffs = list(mafCutoff = 0.2, macCutoff = 0, imissCutoff = 1)
    ))
    expect_lt(length(pair$variantIds), nrow(df))
    expect_equal(nrow(pair$sumstat), length(pair$variantIds))
    expect_equal(pair$sumstat$variant, pair$variantIds)
    expect_equal(nrow(pair$LD), length(pair$variantIds))
    expect_equal(rownames(pair$LD), pair$variantIds)
    # The surviving z values are the originals for those variants.
    orig <- df$z[is_in(normalizeVariantId(df$variant_id), pair$variantIds)]
    expect_equal(pair$sumstat$z, orig)
})

test_that(".cbSumstatPair returns NULL when a cutoff removes everything", {
    ss <- .cbf_qcd()
    df <- getSumStatsDf(
        ss,
        study = ss$study[[1L]],
        context = ss$context[[1L]],
        trait = ss$trait[[1L]],
        require = "Z"
    )
    expect_null(suppressMessages(.cbSumstatPair(
        df = df,
        ldSketch = getLdSketch(ss),
        cutoffs = list(mafCutoff = 0.99, macCutoff = 0, imissCutoff = 1)
    )))
})


# ===========================================================================
# Outcome naming and canonical allele flipping
# ===========================================================================

test_that("outcome names are defaulted, context-qualified and de-duplicated", {
    # Three distinct problems: an unnamed column, the same trait appearing in
    # two contexts, and a name that clashes with one already assigned.
    f <- pecotmr:::.cbTraitName
    expect_equal(f(NULL, "brain", character(0), c("a", "b")), "outcome3")
    expect_equal(f("", "brain", character(0), "a"), "outcome2")
    # Same trait in more than one context gets the context prefix.
    expect_equal(f("G1", "brain", "G1", character(0)), "brain_G1")
    # A clash with an already-assigned name is made unique.
    expect_equal(f("G1", "brain", character(0), "G1"), "G1.1")
    # Nothing to fix.
    expect_equal(f("G1", "brain", character(0), character(0)), "G1")
})

test_that("flipping to canonical drops a matrix with no shared variants", {
    # Zero columns but the sample rows preserved, so the result still binds
    # against the other contexts rather than collapsing.
    m <- matrix(
        1:4,
        2L,
        2L,
        dimnames = list(NULL, c("chr9:1:A:G", "chr9:2:C:T"))
    )
    out <- pecotmr:::.cbFlipMatrixToCanonical(m, "chr1:100:A:G")
    expect_equal(dim(out), c(2L, 0L))
})

test_that("flipping a sumstat/LD pair with no shared variants yields NULL", {
    # NULL, not an empty pair: there is nothing to colocalize, and an empty
    # sumstat would look like a fitted-but-null result.
    pair <- list(
        sumstat = tibble(variant = "chr9:1:A:G", z = 1),
        LD = diag(1)
    )
    expect_null(pecotmr:::.cbFlipPairToCanonical(pair, "chr1:100:A:G"))
})


test_that("outcome info is the empty frame when nothing contributed", {
    # Same shape as a populated one, so the downstream bind and join still
    # work on a run that produced no outcomes.
    out <- pecotmr:::.cbOutcomeInfo(list(), NULL, NULL)
    expect_equal(nrow(out), 0L)
    expect_equal(
        colnames(out),
        colnames(pecotmr:::.cbEmptyOutcomeInfo())
    )
})

test_that("a MultiStudyQtlDataset with no embedded sumstats yields no pairs", {
    # The sumstats arm is optional; its absence is an empty bundle rather
    # than an error, so a purely individual-level collection still runs.
    data(multiStudyQtlDatasetExample)
    out <- pecotmr:::.cbMultiStudySumstats(
        multiStudyQtlDatasetExample,
        contexts = NULL,
        traitId = NULL
    )
    expect_equal(out$qtlPairs, list())
    expect_null(out$qtlLdSketch)
})

test_that("the signal screen passes Y through when no cutoff applies", {
    # A NULL / absent screen spec means "do not screen", which is different
    # from a screen that everything failed.
    Y <- matrix(1:4, 2L, 2L, dimnames = list(c("s1", "s2"), c("t1", "t2")))
    expect_identical(
        pecotmr:::.cbApplyScreen(NULL, Y, "brain", NULL),
        Y
    )
})

# ---------------------------------------------------------------------------
# Per-context skips: a context that cannot supply usable X/Y drops out with an
# explanation, rather than reaching colocboost as an empty or misaligned pair.
# ---------------------------------------------------------------------------

.cbp_mat <- function(rowNames, colNames) {
    matrix(
        0,
        nrow = length(rowNames),
        ncol = length(colNames),
        dimnames = list(rowNames, colNames)
    )
}

test_that(".cbBuildContextXY skips a context with no genotypes", {
    local_mocked_bindings(
        .cbResidualizedY = function(...) .cbp_mat(c("s1", "s2"), "f1"),
        .cbResidualizedX = function(...) NULL,
        .package = "pecotmr"
    )
    expect_null(pecotmr:::.cbBuildContextXY("c1", list()))
})

test_that(".cbBuildContextXY skips a context with no shared samples", {
    local_mocked_bindings(
        .cbResidualizedY = function(...) .cbp_mat(c("s1", "s2"), "f1"),
        .cbResidualizedX = function(...) {
            .cbp_mat(c("z9", "z8"), c("chr1:1:A:G", "chr1:2:C:T"))
        },
        .package = "pecotmr"
    )
    # X and Y are each non-empty, but they describe disjoint sample sets.
    expect_message(
        res <- pecotmr:::.cbBuildContextXY("c1", list()),
        "no samples shared between residualized X and Y"
    )
    expect_null(res)
})

test_that(".cbResidualizedX reports why genotypes were unavailable", {
    local_mocked_bindings(
        getResidualizedGenotypes = function(...) stop("kaboom"),
        .package = "pecotmr"
    )
    # The underlying message is carried through so the skip is diagnosable.
    expect_message(
        res <- pecotmr:::.cbResidualizedX(
            NULL,
            "c1",
            NULL,
            NULL,
            NULL,
            NULL
        ),
        "residualized genotypes unavailable: kaboom"
    )
    expect_null(res)
})

test_that(".cbApplyScreen keeps the outcomes that clear the screen", {
    Y <- .cbp_mat(c("s1", "s2"), c("f1", "f2"))
    local_mocked_bindings(
        .cbPipSkipOutcomes = function(X, Y, cutoff) Y[, 1, drop = FALSE],
        .package = "pecotmr"
    )
    out <- pecotmr:::.cbApplyScreen(NULL, Y, "c1", 0.5)
    expect_equal(ncol(out), 1L)
    expect_equal(colnames(out), "f1")
})

test_that(".cbToResultObject skips a separate_gwas study that produced none", {
    raw <- list(
        xqtl_coloc = NULL,
        joint_gwas = NULL,
        separate_gwas = list(G1 = NULL, G2 = .cbr_fake()),
        computing_time = list()
    )
    x <- pecotmr:::.cbToResultObject(raw, .cbr_info())
    # G1 contributes no row at all; its key must not survive as an empty one.
    expect_equal(nrow(x), 1L)
    expect_equal(as.character(x$gwasStudy), "G2")
    expect_equal(as.character(x$analysis), "separate_gwas")
})

test_that(".cbRunXqtlOnly passes a focal outcome through as an index", {
    local_mocked_bindings(
        .cbRun = function(label, args) {
            list(result = list(focal = args$focal_outcome_idx), time = 0)
        },
        .package = "pecotmr"
    )
    bundle <- list(
        outcomeNames = c("tA", "tB"),
        Y = list(1, 2),
        X = list(),
        dict_YX = NULL
    )
    empty <- pecotmr:::.cbMergeSumstatBundles(list())
    run <- suppressMessages(
        pecotmr:::.cbRunXqtlOnly(bundle, empty, TRUE, "tB", list())
    )
    # colocboost wants a position, not a name.
    expect_equal(run$result$focal, 2L)
    absent <- suppressMessages(
        pecotmr:::.cbRunXqtlOnly(bundle, empty, TRUE, "nope", list())
    )
    expect_null(absent$result$focal)
    none <- suppressMessages(
        pecotmr:::.cbRunXqtlOnly(bundle, empty, TRUE, NULL, list())
    )
    expect_null(none$result$focal)
})

# A summary-level QTL side used to make xqtlColoc a silent no-op: the run was
# gated on having an individual bundle, so the default flags returned an empty
# ColocBoostResult with all-NULL timings and no message.
test_that(".cbRunXqtlOnly runs on a summary-level QTL side alone", {
    captured <- NULL
    local_mocked_bindings(
        .cbRun = function(label, args) {
            captured <<- args
            list(result = list(ran = TRUE), time = 0)
        },
        .package = "pecotmr"
    )
    ssBundle <- pecotmr:::.cbMergeSumstatBundles(list(
        qtlA = list(sumstat = list(z = 1), LD = diag(2)),
        qtlB = list(sumstat = list(z = 2), LD = diag(2))
    ))
    run <- suppressMessages(
        pecotmr:::.cbRunXqtlOnly(NULL, ssBundle, FALSE, "qtlB", list())
    )
    expect_equal(run$result$ran, TRUE)
    expect_null(captured$X)
    expect_null(captured$Y)
    expect_null(captured$dict_YX)
    expect_equal(names(captured$sumstat), c("qtlA", "qtlB"))
    expect_equal(captured$outcome_names, c("qtlA", "qtlB"))
    # focalTrait is honored on the summary-level side too.
    expect_equal(captured$focal_outcome_idx, 2L)
    # Identical LD matrices dedupe to one, so the dict points both at it.
    expect_equal(unname(captured$dict_sumstatLD[, "LD"]), c(1L, 1L))
})

test_that(".cbRunXqtlOnly combines individual and summary-level QTL sides", {
    captured <- NULL
    local_mocked_bindings(
        .cbRun = function(label, args) {
            captured <<- args
            list(result = NULL, time = 0)
        },
        .package = "pecotmr"
    )
    bundle <- list(
        outcomeNames = c("tA", "tB"),
        Y = list(1, 2),
        X = list(),
        dict_YX = NULL
    )
    ssBundle <- pecotmr:::.cbMergeSumstatBundles(list(
        qtlC = list(sumstat = list(z = 3), LD = diag(2))
    ))
    suppressMessages(
        pecotmr:::.cbRunXqtlOnly(bundle, ssBundle, TRUE, "qtlC", list())
    )
    expect_equal(captured$outcome_names, c("tA", "tB", "qtlC"))
    expect_equal(captured$focal_outcome_idx, 3L)
})

test_that(".cbRunVariants: xqtlColoc runs on a QTL-only sumstat bundle", {
    called <- character(0)
    local_mocked_bindings(
        .cbRunXqtlOnly = function(...) {
            called <<- c(called, "xqtl")
            list(result = list(stub = TRUE), time = 1)
        },
        .cbOutcomeInfo = function(...) pecotmr:::.cbEmptyOutcomeInfo(),
        .package = "pecotmr"
    )
    merged <- pecotmr:::.cbMergeSumstatBundles(list(
        qtlA = list(sumstat = list(z = 1), LD = diag(2)),
        gwasG = list(sumstat = list(z = 2), LD = diag(2))
    ))
    qtlOnly <- pecotmr:::.cbMergeSumstatBundles(list(
        qtlA = list(sumstat = list(z = 1), LD = diag(2))
    ))
    out <- suppressMessages(pecotmr:::.cbRunVariants(
        NULL,
        merged,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE,
        focalTrait = NULL,
        dotArgs = list(),
        qtlSumstatBundle = qtlOnly
    ))
    expect_equal(called, "xqtl")
    expect_false(is.null(getComputingTime(out)$Analysis$xqtl_coloc))
})

test_that(".cbRunVariants warns instead of silently skipping an analysis", {
    local_mocked_bindings(
        .cbOutcomeInfo = function(...) pecotmr:::.cbEmptyOutcomeInfo(),
        .package = "pecotmr"
    )
    # Individual-level QTL side, no sumstats anywhere: the two GWAS variants
    # cannot run, and the caller is told so rather than getting a bare empty.
    ind <- list(
        outcomeNames = "tA",
        Y = list(1),
        X = list(),
        dict_YX = NULL
    )
    warnings <- capture_warnings(
        suppressMessages(pecotmr:::.cbRunVariants(
            ind,
            pecotmr:::.cbMergeSumstatBundles(list()),
            xqtlColoc = FALSE,
            jointGwas = TRUE,
            separateGwas = TRUE,
            focalTrait = NULL,
            dotArgs = list()
        ))
    )
    expect_length(warnings, 2L)
    expect_match(warnings[[1L]], "jointGwas = TRUE was requested")
    expect_match(warnings[[2L]], "separateGwas = TRUE was requested")
})

test_that(".cbRunVariants warns when xqtlColoc has only GWAS sumstats", {
    local_mocked_bindings(
        .cbOutcomeInfo = function(...) pecotmr:::.cbEmptyOutcomeInfo(),
        .package = "pecotmr"
    )
    gwasOnly <- pecotmr:::.cbMergeSumstatBundles(list(
        gwasG = list(sumstat = list(z = 1), LD = diag(2))
    ))
    out <- expect_warning(
        suppressMessages(pecotmr:::.cbRunVariants(
            NULL,
            gwasOnly,
            xqtlColoc = TRUE,
            jointGwas = FALSE,
            separateGwas = FALSE,
            focalTrait = NULL,
            dotArgs = list(),
            qtlSumstatBundle = pecotmr:::.cbMergeSumstatBundles(list())
        )),
        "xqtlColoc = TRUE was requested"
    )
    expect_null(getComputingTime(out)$Analysis$xqtl_coloc)
})

test_that(".cbAppendGwasPairs disambiguates a colliding study key", {
    local_mocked_bindings(
        .cbRequireSumStatsQc = function(...) invisible(NULL),
        .cbGwasSumStatsBundle = function(gwasSumStats, cutoffs) {
            list(dup = "GWAS")
        },
        .package = "pecotmr"
    )
    out <- pecotmr:::.cbAppendGwasPairs(
        list(dup = "QTL", other = "X"),
        "notNull",
        qtlLdSketch = NULL
    )
    # The QTL pair keeps its key; the GWAS pair is suffixed rather than
    # overwriting it.
    expect_equal(names(out), c("dup", "other", "dup.1"))
    expect_equal(out[["dup"]], "QTL")
    expect_equal(out[["dup.1"]], "GWAS")
})
