context("colocboostPipeline")

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
# (.cbBuildLdArgs, .cbMergeSumstatBundles, .cbResultLdRefs, .cbEmptyResult)
# are exercised directly.
# ===========================================================================

.cbp_mockExtractor <- function(seed = 11, n_samples = 30L) {
    function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(seed)
        panel <- matrix(
            rbinom(n_samples * nrow(snpInfo(handle)), 2, 0.3),
            nrow = n_samples,
            ncol = nrow(snpInfo(handle)),
            dimnames = list(sampleIds(handle), snpInfo(handle)$SNP)
        )
        sub <- panel[, snpIdx, drop = FALSE]
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", snpInfo(handle)$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = snpInfo(handle)$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = snpInfo(handle)$SNP[snpIdx],
            A1 = snpInfo(handle)$A1[snpIdx],
            A2 = snpInfo(handle)$A2[snpIdx]
        )
        cd <- S4Vectors::DataFrame(
            sampleId = sampleIds(handle),
            row.names = sampleIds(handle)
        )
        dosage <- t(sub)
        rownames(dosage) <- snpInfo(handle)$SNP[snpIdx]
        colnames(dosage) <- sampleIds(handle)
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
        map(contexts, function(.) .cbp_makeSe(traits = traits)),
        contexts
    )
    QtlDataset(
        studyName = "study1",
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
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .cbp_makeHandle(),
        qcOptions = if (qc) list(step1 = "ok") else list()
    )
}

.cbp_makeGwasSumStats <- function(qc = TRUE, ldReference = .cbp_makeHandle()) {
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
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = ldReference,
        qcOptions = if (qc) list(step1 = "ok") else list()
    )
}

# An LdMatrix over the fixtures' five variants, so a side can be put on a
# reference other than the shared genotype panel.
.cbp_makeLdMatrix <- function(rho = 0.2) {
    ids <- sprintf("chr1:%d:A:G", 100L * (1:5))
    r <- rho^abs(outer(1:5, 1:5, "-"))
    LdMatrix(`dimnames<-`(r, list(ids, ids)))
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

test_that(".cbBuildLdArgs: factors route to X_ref, whatever their shape", {
    X1 <- matrix(0, 10, 4)
    res <- pecotmr:::.cbBuildLdArgs(list(X1), "X_ref")
    expect_identical(names(res), "X_ref")
    # A factor with as many rows as variants is still a factor.
    square <- pecotmr:::.cbBuildLdArgs(list(diag(4)), "X_ref")
    expect_identical(names(square), "X_ref")
})

test_that(".cbBuildLdArgs: empty list returns empty list", {
    expect_equal(pecotmr:::.cbBuildLdArgs(list()), list())
})

test_that(".cbResidualizedY skips a context whose residualization fails", {
    # A failing context is skipped with a message naming it, so one bad
    # context does not end the whole run.
    local_mocked_bindings(
        .qtlResidualizedPhenotypes = function(...) {
            stop("no samples shared")
        },
        .package = "pecotmr"
    )
    expect_message(
        out <- pecotmr:::.cbResidualizedY(NULL, "brain", NULL, NULL),
        "skipping context 'brain' (residualized phenotypes unavailable)",
        fixed = TRUE
    )
    expect_null(out)
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
    expect_identical(res$ldArg, "LD")
})

test_that(".cbMergeSumstatBundles records factors and turns a mix into LD", {
    ids <- c("v1", "v2")
    f <- matrix(c(1, -1, 0.5, -0.5), 2, dimnames = list(NULL, ids))
    factors <- list(
        a = list(sumstat = data.frame(z = 1:2), LD = f, isFactor = TRUE),
        b = list(sumstat = data.frame(z = 3:4), LD = f, isFactor = TRUE)
    )
    res <- pecotmr:::.cbMergeSumstatBundles(factors)
    expect_identical(res$ldArg, "X_ref")
    expect_equal(length(res$LD), 1L)
    # colocboost takes LD or X_ref, never both, so a run whose two sides use
    # different kinds of reference gets the factor as its LD matrix.
    ld <- matrix(c(1, 0.3, 0.3, 1), 2, dimnames = list(ids, ids))
    mixed <- list_assign(
        factors,
        b = list(sumstat = data.frame(z = 3:4), LD = ld, isFactor = FALSE)
    )
    res <- pecotmr:::.cbMergeSumstatBundles(mixed)
    expect_identical(res$ldArg, "LD")
    expect_equal(res$LD[[1L]], crossprod(f) / (nrow(f) - 1L))
    expect_identical(res$LD[[2L]], ld)
    expect_equal(unname(res$dict_sumstatLD[, "LD"]), 1:2)
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
    expect_true(is_in("xqtl_coloc", computingTime(out)$analysis))
    expect_true("X" %in% names(capturedArgs))
})

test_that("colocboostPipeline(QtlDataset): no QTL bundle and no GWAS returns empty result", {
    qd <- .cbp_makeQtlDataset()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    # A signal screen no outcome clears skips every context; the bundle
    # ends up NULL and the driver short-circuits.
    out <- suppressMessages(suppressWarnings(
        colocboostPipeline(
            qd,
            contexts = "brain",
            signalScreenParam = SignalScreenParam(absZ = 1e6),
            xqtlColoc = TRUE
        )
    ))
    # No confidence set, and no analysis ran to be timed.
    expect_s4_class(out, "ColocBoostResult")
    expect_length(out, 0L)
    expect_equal(nrow(computingTime(out)), 0L)
})

test_that("colocboostPipeline refuses a trait selection that matches nothing", {
    # As fineMappingPipeline and twasWeightsPipeline do, rather than
    # returning an empty result.
    expect_error(
        colocboostPipeline(.cbp_makeQtlDataset(), traitId = "nope"),
        "colocboostPipeline(QtlDataset): no traits selected.",
        fixed = TRUE
    )
    expect_error(
        colocboostPipeline(.cbp_makeQtlSumStats(), traitId = "nope"),
        "colocboostPipeline(QtlSumStats): no entries matched the ",
        fixed = TRUE
    )
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
    expect_true(is_in("joint_gwas", computingTime(out)$analysis))
})

test_that("colocboostPipeline: separateGwas runs once per GWAS study", {
    ss <- .cbp_makeQtlSumStats()
    gs <- .cbp_makeGwasSumStats()
    seen <- list()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            seen[[length(seen) + 1L]] <<- list(...)
            list(round = length(seen))
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
    # ONE run per GWAS study, not per merged sumstat label. The QTL study is
    # a non-focal outcome of that run, not a focal "GWAS" of its own: the
    # merged bundle holds both sides, so iterating it ran the QTL as its own
    # focal GWAS and reported "0 contexts".
    expect_equal(length(seen), 1L)
    call <- seen[[1L]]
    expect_equal(length(call$sumstat), 2L)
    expect_equal(call$outcome_names, c("Q1:c1:t1", "G1"))
    # The focal outcome is the GWAS, and it is last.
    expect_equal(call$focal_outcome_idx, 2L)
    expect_equal(call$outcome_names[[call$focal_outcome_idx]], "G1")
    # The QTL side's LD travels with it, so the dict covers both sumstats.
    expect_equal(nrow(call$dict_sumstatLD), 2L)
    expect_true(
        is_in("separate_gwas", computingTime(out)$analysis)
    )
})

test_that("separateGwas refuses a QTL-only sumstat input", {
    ss <- .cbp_makeQtlSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    # Guarded on GWAS presence, not on "any summary statistics": the old
    # `hasSs` guard let a QTL-only input through and ran one focal pass per
    # QTL study.
    warnings <- character()
    withCallingHandlers(
        suppressMessages(colocboostPipeline(
            ss,
            xqtlColoc = FALSE,
            jointGwas = FALSE,
            separateGwas = TRUE
        )),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_true(any(str_detect(warnings, "GWAS summary statistics")))
})

test_that("colocboostPipeline: no analysis flag set is an error", {
    ss <- .cbp_makeQtlSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    # All three FALSE is the one combination that is not a request. It used
    # to run silently and hand back an empty result.
    expect_error(
        colocboostPipeline(
            ss,
            xqtlColoc = FALSE,
            jointGwas = FALSE,
            separateGwas = FALSE
        ),
        "no analysis selected"
    )
    # Any single flag is still a valid request -- they are independent, not
    # a mutually exclusive set.
    for (flag in c("xqtlColoc", "jointGwas", "separateGwas")) {
        args <- set_names(
            list(FALSE, FALSE, FALSE),
            c("xqtlColoc", "jointGwas", "separateGwas")
        )
        args[[flag]] <- TRUE
        expect_silent(exec(
            pecotmr:::.cbAssertAnyRun,
            !!!set_names(unname(args), NULL)
        ))
    }
})

test_that("a GWAS on a different LD reference runs, and both are recorded", {
    # The two sides need not share a reference: a GWAS panel with another
    # sample set is matched against on its own rather than refused.
    ss <- .cbp_makeQtlSumStats()
    gs <- .cbp_makeGwasSumStats(
        ldReference = .cbp_makeHandle(sample_prefix = "z")
    )
    seen <- list()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            seen[[length(seen) + 1L]] <<- list(...)
            list(jointly_run = TRUE)
        },
        .package = "colocboost"
    )
    out <- suppressMessages(
        colocboostPipeline(
            ss,
            gwasSumStats = gs,
            xqtlColoc = FALSE,
            jointGwas = TRUE
        )
    )
    expect_equal(length(seen), 1L)
    expect_equal(nrow(seen[[1L]]$dict_sumstatLD), 2L)
    # Each outcome keys to the reference it was matched against.
    refs <- ldReference(out)
    expect_length(refs, 2L)
    info <- out@outcomeInfo
    keys <- info$ldReference[match(c("Q1:c1:t1", "G1"), info$name)]
    expect_identical(refs[[keys[[1L]]]], ldReference(ss)[[1L]])
    expect_identical(refs[[keys[[2L]]]], ldReference(gs)[[1L]])
})

test_that("a panel and an LdMatrix on opposite sides reach colocboost as LD", {
    ss <- .cbp_makeQtlSumStats()
    ldm <- .cbp_makeLdMatrix()
    gs <- .cbp_makeGwasSumStats(ldReference = ldm)
    seen <- list()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            seen[[length(seen) + 1L]] <<- list(...)
            list(jointly_run = TRUE)
        },
        .package = "colocboost"
    )
    out <- suppressMessages(
        colocboostPipeline(
            ss,
            gwasSumStats = gs,
            xqtlColoc = TRUE,
            jointGwas = TRUE
        )
    )
    nSumstats <- map_int(seen, function(call) length(call$sumstat))
    xqtl <- seen[[which(nSumstats == 1L)]]
    joint <- seen[[which(nSumstats == 2L)]]
    # The QTL-only run has one kind of reference and keeps the panel's factor.
    expect_null(xqtl$LD)
    factor <- xqtl$X_ref[[1L]]
    # colocboost takes LD or X_ref but not both, so the joint run gets the
    # panel's factor as its LD matrix alongside the LdMatrix's own.
    expect_null(joint$X_ref)
    dict <- joint$dict_sumstatLD
    expect_equal(
        joint$LD[[dict[1L, "LD"]]],
        crossprod(factor) / (nrow(factor) - 1L)
    )
    expect_equal(
        unname(joint$LD[[dict[2L, "LD"]]]),
        0.2^abs(outer(1:5, 1:5, "-"))
    )
    expect_identical(
        unname(ldReference(out)),
        list(ldReference(ss)[[1L]], ldm)
    )
})

test_that(".cbResultLdRefs keys each outcome to its pair's reference", {
    panel <- ldReference(.cbp_makeQtlSumStats())[[1L]]
    ldm <- .cbp_makeLdMatrix()
    pairs <- list(
        `Q1:c1:t1` = list(ldReference = panel),
        `Q1:c2:t1` = list(ldReference = panel),
        G1 = list(ldReference = ldm)
    )
    out <- pecotmr:::.cbResultLdRefs(pairs)
    # Each distinct reference once, however many outcomes used it.
    expect_identical(out$refs, list(ld1 = panel, ld2 = ldm))
    expect_equal(
        out$keys,
        c(`Q1:c1:t1` = "ld1", `Q1:c2:t1` = "ld1", G1 = "ld2")
    )
    # No pairs, no references.
    none <- pecotmr:::.cbResultLdRefs(list())
    expect_identical(none$refs, list())
    expect_length(none$keys, 0L)
})

# ===========================================================================
# GWAS case/control: optional nCase/nControl columns + effective-N wiring
# ===========================================================================

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
        ldReference = .cbp_makeHandle(),
        qcOptions = list(ok = 1)
    )
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    # case/control -> effective N = 4 / (1/500 + 1/1500) = 1500
    gcc <- exec(GwasSumStats, !!!c(base, list(nCase = 500, nControl = 1500)))
    bcc <- pecotmr:::.cbGwasSumStatsBundle(gcc)
    expect_true(all(bcc[["G1"]]$sumstat$n == 4 / (1 / 500 + 1 / 1500)))
    # quantitative (no nCase/nControl) -> per-variant N (1000)
    bq <- pecotmr:::.cbGwasSumStatsBundle(exec(GwasSumStats, !!!base))
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
        dimnames = list(
            paste0("s", seq_len(n)),
            sprintf("chr1:%d:A:G", 100L * seq_len(p))
        )
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

test_that(".cbResolveCutoff: screen/scalar uniform, named vector per context", {
    # A resolved screen object applies to every context alike.
    sc <- list(metric = "absZ", cutoff = 5)
    expect_identical(pecotmr:::.cbResolveCutoff(sc, "brain"), sc)
    expect_identical(pecotmr:::.cbResolveCutoff(sc, "blood"), sc)
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

test_that("colocboostPipeline(MultiStudy) refuses a trait no study carries", {
    expect_error(
        colocboostPipeline(.cbp_makeMultiStudy(), traitId = "nope"),
        "colocboostPipeline(MultiStudyQtlDataset): no traits selected.",
        fixed = TRUE
    )
    # A trait only the embedded sumstats carry still selects something.
    expect_silent(pecotmr:::.cbMsAssertTraitsSelected(
        .cbp_makeMultiStudy(),
        NULL,
        "t1"
    ))
})

test_that("colocboostPipeline(MultiStudy): context selection skips a study", {
    # As in fineMappingPipeline: a context only some studies have is skipped
    # for the others, with a warning, and one no study has is refused.
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) list(stub = TRUE),
        .package = "colocboost"
    )
    expect_warning(
        out <- suppressMessages(colocboostPipeline(
            .cbp_makeMultiStudy(),
            contexts = "c1",
            xqtlColoc = TRUE
        )),
        "study 'study1' is missing requested context\\(s\\): c1"
    )
    # Only the summary statistics' outcome reached colocboost.
    expect_identical(out@outcomeInfo$name, "Q1:c1:t1")
    expect_error(
        colocboostPipeline(.cbp_makeMultiStudy(), contexts = "ghost"),
        "no study has: ghost"
    )
})

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
    expect_true(is_in("xqtl_coloc", computingTime(out)$analysis))
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
    expect_message(
        out <- colocboostPipeline(qd, xqtlColoc = TRUE),
        "xQTL-only ColocBoost failed"
    )
    # The run was attempted and timed, and contributes no confidence set.
    expect_true(is_in("xqtl_coloc", computingTime(out)$analysis))
    expect_length(out, 0L)
})

test_that(".cbIndividualBundle: multi-context outcomes get context prefixes", {
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
    expect_true(is_in("xqtl_coloc", computingTime(out)$analysis))
    # ENSG_A is in both contexts, so each outcome is prefixed with its context
    # to tell the two apart, and colocboost is handed those names.
    prefixed <- c("brain_ENSG_A", "liver_ENSG_A")
    expect_identical(names(capturedArgs$Y), prefixed)
    expect_identical(capturedArgs$outcome_names, prefixed)
    # The contexts share genotypes and samples, so their X deduplicates to the
    # first context's, and dict_YX maps both outcomes to it (the xMatch).
    expect_identical(names(capturedArgs$X), "brain")
    expect_identical(capturedArgs$dict_YX, cbind(Y = 1:2, X = c(1L, 1L)))
    # Each prefixed name still resolves to the context and trait behind it.
    info <- out@outcomeInfo
    expect_identical(info$name, prefixed)
    expect_identical(info$context, c("brain", "liver"))
    expect_identical(info$trait, c("ENSG_A", "ENSG_A"))
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
    )
})

test_that("colocboostPipeline(QtlDataset): pipCutoffToSkip dropping every outcome -> empty", {
    qd <- .cbp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- suppressMessages(suppressWarnings(
        colocboostPipeline(
            qd,
            xqtlColoc = TRUE,
            signalScreenParam = SignalScreenParam(pip = 0.9999)
        )
    ))
    # The screen leaves no outcome, so nothing runs and no set is reported.
    expect_length(out, 0L)
    expect_equal(nrow(computingTime(out)), 0L)
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
    expect_false(is_in("a", colnames(res)))
})

test_that(".cbSumstatPair: NULL / empty df -> NULL", {
    expect_null(pecotmr:::.cbSumstatPair(NULL, .cbp_makeHandle()))
    expect_null(pecotmr:::.cbSumstatPair(data.frame(), .cbp_makeHandle()))
})

test_that(".cbQtlSumStatsBundle: NULL / context / trait filters and empty result", {
    ss <- .cbp_makeQtlSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    expect_equal(pecotmr:::.cbQtlSumStatsBundle(NULL), list())
    expect_length(pecotmr:::.cbQtlSumStatsBundle(ss, contexts = "c1"), 1L)
    expect_length(pecotmr:::.cbQtlSumStatsBundle(ss, traitId = "t1"), 1L)
    expect_equal(pecotmr:::.cbQtlSumStatsBundle(ss, contexts = "ghost"), list())
})

test_that(".cbGwasSumStatsBundle: NULL -> empty list", {
    expect_equal(pecotmr:::.cbGwasSumStatsBundle(NULL), list())
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
    expect_true("var_y" %in% names(pair$sumstat))
    expect_equal(unique(pair$sumstat$var_y), 0.7)
    # NA variant_id -> formatVariantId fallback; no panel overlap -> NULL
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

test_that(".cbFlipPairToCanonical negates a factor's swapped columns", {
    ss <- data.frame(
        z = c(2, -1),
        n = c(1000, 1000),
        variant = c("chr1:100:A:G", "chr1:200:C:T"),
        stringsAsFactors = FALSE
    )
    f <- matrix(c(1, -1, 0, 0.5, 0.5, -1), nrow = 3)
    p <- list(sumstat = ss, LD = f, variantIds = ss$variant, isFactor = TRUE)
    canonical <- c("chr1:100:G:A", "chr1:200:C:T")
    out <- pecotmr:::.cbFlipPairToCanonical(p, canonical)
    expect_true(out$isFactor)
    expect_identical(colnames(out$LD), canonical)
    expect_equal(unname(out$LD), cbind(-f[, 1], f[, 2]))
    # The factor's LD is the flipped LD.
    expect_equal(
        unname(crossprod(out$LD)),
        unname(crossprod(f) * outer(c(-1, 1), c(-1, 1)))
    )
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

test_that("colocboostPipeline(MultiStudy): a study with no bundle is skipped", {
    mt <- .cbp_makeMultiStudy() # qd (study1, ENSG_A) + ss (Q1, t1)
    seen <- list()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            seen[[length(seen) + 1L]] <<- list(...)
            list(stub = TRUE)
        },
        .package = "colocboost"
    )
    # traitId="t1" matches the embedded sumstats but not the QtlDataset -> its
    # per-study bundle is NULL and skipped; the sumstat side still runs.
    out <- suppressMessages(suppressWarnings(
        colocboostPipeline(mt, traitId = "t1", xqtlColoc = TRUE)
    ))
    # One xQTL-only run over the summary statistics alone: no individual-level
    # data reached colocboost, and no "study1:" outcome was named.
    expect_length(seen, 1L)
    call <- seen[[1L]]
    expect_false(any(is_in(c("X", "Y", "dict_YX"), names(call))))
    expect_identical(call$outcome_names, "Q1:c1:t1")
    expect_identical(names(call$sumstat), "Q1:c1:t1")
    # The result times that run and records only the summary-level outcome,
    # keyed to the embedded sumstats' LD reference.
    expect_identical(computingTime(out)$analysis, "xqtl_coloc")
    expect_length(out, 0L)
    info <- out@outcomeInfo
    expect_identical(info$name, "Q1:c1:t1")
    expect_identical(info$study, "Q1")
    expect_identical(info$dataForm, "sumstats")
    expect_length(ldReference(out), 1L)
    expect_identical(
        ldReference(out)[[info$ldReference]],
        ldReference(sumStats(mt))[[1L]]
    )
})


# ---------------------------------------------------------------------------
# RSS panel filters.
#
# colocboostPipeline never runs summary-statistic QC itself -- that lives in
# summaryStatsQc() -- so these cutoffs act at analysis time, against the LD
# reference panel, on both the QTL and the GWAS sumstat sides.
# ---------------------------------------------------------------------------

# summaryStatsQc() of the bundled example, run once for the tests below.
.cbf_ss <- local({
    data(qtlSumStatsExample, envir = environment())
    suppressMessages(summaryStatsQc(qtlSumStatsExample))
})

# @noRd
.cbf_n <- function(ss, ...) {
    b <- suppressMessages(.cbQtlSumStatsBundle(
        ss,
        cutoffs = .panelCutoffs(PanelFilterParam(...))
    ))
    if (length(b) == 0L) 0L else length(b[[1L]]$variantIds)
}

test_that("colocboost sumstat bundle filters nothing by default", {
    ss <- .cbf_ss
    expect_equal(.cbf_n(ss), sum(lengths(ss)))
})

test_that("colocboost sumstat bundle drops panel-rare variants", {
    ss <- .cbf_ss
    full <- .cbf_n(ss)
    loose <- .cbf_n(ss, mafCutoff = 0.05)
    tight <- .cbf_n(ss, mafCutoff = 0.2)
    expect_lt(loose, full)
    expect_lt(tight, loose)
})

test_that("colocboost RSS cutoffs match .panelVariantFilter", {
    ss <- .cbf_ss
    ids <- normalizeVariantId(
        .sumStatsDf(
            ss,
            studyName = ss$study[[1L]],
            context = ss$context[[1L]],
            trait = ss$trait[[1L]],
            require = "Z"
        )$variant_id
    )
    for (cut in c(0.05, 0.2)) {
        expect_equal(
            .cbf_n(ss, mafCutoff = cut),
            length(.panelVariantFilter(
                ldReference(ss)[[1L]],
                ids,
                PanelFilterParam(mafCutoff = cut)
            )),
            label = str_c("mafCutoff ", cut)
        )
    }
})

test_that("colocboost RSS treats MAC as a MAF equivalent", {
    ss <- .cbf_ss
    nSamp <- ncol(ldReference(ss)[[1L]])
    expect_equal(
        .cbf_n(ss, macCutoff = 0.1 * 2 * nSamp),
        .cbf_n(ss, mafCutoff = 0.1)
    )
})

test_that("colocboost RSS honours a missingness cutoff", {
    ss <- .cbf_ss
    expect_lt(.cbf_n(ss, imissCutoff = 0), sum(lengths(ss)))
    expect_equal(.cbf_n(ss, imissCutoff = 1), sum(lengths(ss)))
})

test_that(".cbSumstatPair keeps sumstat rows aligned to the LD matrix", {
    # The row filter indexes df by the FULL-length variant id vector, so the
    # filtered set has to be passed to the LD build separately rather than by
    # shrinking that vector.
    ss <- .cbf_ss
    df <- .sumStatsDf(
        ss,
        studyName = ss$study[[1L]],
        context = ss$context[[1L]],
        trait = ss$trait[[1L]],
        require = "Z"
    )
    sketch <- ldReference(ss)[[1L]]
    pair <- suppressMessages(.cbSumstatPair(
        df = df,
        ldReference = sketch,
        cutoffs = list(mafCutoff = 0.2, macCutoff = 0, imissCutoff = 1)
    ))
    expect_lt(length(pair$variantIds), nrow(df))
    expect_equal(nrow(pair$sumstat), length(pair$variantIds))
    expect_equal(pair$sumstat$variant, pair$variantIds)
    # A genotype reference gives a factor of its LD: variants on the columns.
    expect_true(pair$isFactor)
    expect_equal(ncol(pair$LD), length(pair$variantIds))
    expect_equal(colnames(pair$LD), pair$variantIds)
    # The surviving z values are the originals for those variants.
    orig <- df$z[is_in(normalizeVariantId(df$variant_id), pair$variantIds)]
    expect_equal(pair$sumstat$z, orig)
})

test_that("colocboost reads a panel's factor and its LdMatrix alike", {
    # A genotype reference reaches colocboost as X_ref, a factor of its LD;
    # an LdMatrix as LD. The same LD gives the same colocalization.
    skip_if_not_installed("colocboost")
    data(qtlSumStatsMulticontextExample)
    ss <- suppressMessages(summaryStatsQc(qtlSumStatsMulticontextExample))
    panel <- ldReference(ss)[[1L]]
    withLd <- QtlSumStats(
        studyName = ss$study,
        context = ss$context,
        trait = ss$trait,
        entry = map(seq_len(length(ss)), .collectionEntry, x = ss),
        genome = "hg19",
        ldReference = LdMatrix(
            unname(computeLd(panel)),
            variants = SummarizedExperiment::rowRanges(panel)
        ),
        qcOptions = qcOptions(ss)
    )
    set.seed(1)
    viaPanel <- suppressWarnings(suppressMessages(colocboostPipeline(ss)))
    set.seed(1)
    viaMatrix <- suppressWarnings(suppressMessages(colocboostPipeline(withLd)))
    expect_gt(nrow(as.data.frame(viaPanel)), 0L)
    expect_equal(as.data.frame(viaPanel), as.data.frame(viaMatrix))
})

test_that(".cbSumstatPair returns NULL when a cutoff removes everything", {
    ss <- .cbf_ss
    df <- .sumStatsDf(
        ss,
        studyName = ss$study[[1L]],
        context = ss$context[[1L]],
        trait = ss$trait[[1L]],
        require = "Z"
    )
    expect_null(suppressMessages(.cbSumstatPair(
        df = df,
        ldReference = ldReference(ss)[[1L]],
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
        traitId = NULL,
        cisWindow = NULL
    )
    expect_equal(out$qtlPairs, list())
    expect_null(out$qtlLdReference)
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
        .qtlResidualizedGenotypes = function(...) stop("kaboom"),
        .package = "pecotmr"
    )
    # The cause is chained via `parent`, not flattened into the text, so the
    # rendered message still shows it AND it stays reachable as a condition.
    # catch_cnd() rather than expect_message(), which returns NULL.
    cnd <- rlang::catch_cnd(
        pecotmr:::.cbResidualizedX(
            NULL,
            "c1",
            NULL,
            NULL,
            NULL,
            NULL,
            genotypeFilter = pecotmr:::.qtlResolveFilter()
        ),
        classes = "message"
    )
    expect_match(conditionMessage(cnd), "residualized genotypes unavailable")
    expect_match(conditionMessage(cnd), "kaboom")
    expect_match(conditionMessage(cnd$parent), "kaboom")
    expect_null(suppressMessages(
        pecotmr:::.cbResidualizedX(
            NULL,
            "c1",
            NULL,
            NULL,
            NULL,
            NULL,
            genotypeFilter = pecotmr:::.qtlResolveFilter()
        )
    ))
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
    # A focal outcome is honored on the summary-level side too.
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
        focal = NULL,
        methodArgs = ColocboostOptions(),
        qtlSumstatBundle = qtlOnly
    ))
    expect_equal(called, "xqtl")
    expect_true(is_in("xqtl_coloc", computingTime(out)$analysis))
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
            focal = NULL,
            methodArgs = ColocboostOptions()
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
            focal = NULL,
            methodArgs = ColocboostOptions(),
            qtlSumstatBundle = pecotmr:::.cbMergeSumstatBundles(list())
        )),
        "xqtlColoc = TRUE was requested"
    )
    expect_false(is_in("xqtl_coloc", computingTime(out)$analysis))
})

test_that("ColocboostOptions() options reach the engine by name", {
    # Regression: the bundle was merged with c(), which appends a MethodOptions
    # (a SimpleList) as ONE opaque element instead of splicing its entries.
    # exec() then handed that object to colocboost as a positional argument,
    # so every option set through ColocboostOptions() was silently lost.
    seen <- NULL
    local_mocked_bindings(
        colocboost = function(...) {
            seen <<- list(...)
            list(stub = TRUE)
        },
        .package = "colocboost"
    )
    local_mocked_bindings(
        .cbOutcomeInfo = function(...) pecotmr:::.cbEmptyOutcomeInfo(),
        .package = "pecotmr"
    )
    qtlOnly <- pecotmr:::.cbMergeSumstatBundles(list(
        qtlA = list(sumstat = list(z = 1), LD = diag(2))
    ))
    suppressMessages(pecotmr:::.cbRunVariants(
        NULL,
        qtlOnly,
        xqtlColoc = TRUE,
        jointGwas = FALSE,
        separateGwas = FALSE,
        focal = NULL,
        methodArgs = ColocboostOptions(M = 3L, stop_thresh = 0.02),
        qtlSumstatBundle = qtlOnly
    ))
    expect_equal(seen$M, 3L)
    expect_equal(seen$stop_thresh, 0.02)
    # ...and nothing rode along unnamed, which is how the record used to travel
    expect_false(any(names(seen) == ""))
})

test_that(".cbAppendGwasPairs disambiguates a colliding study key", {
    local_mocked_bindings(
        .cbRequireSumStatsQc = function(...) invisible(NULL),
        .cbGwasSumStatsBundle = function(gwasSumStats, window, cutoffs) {
            list(dup = "GWAS")
        },
        .package = "pecotmr"
    )
    out <- pecotmr:::.cbAppendGwasPairs(
        list(dup = "QTL", other = "X"),
        "notNull"
    )
    # The QTL pair keeps its key; the GWAS pair is suffixed rather than
    # overwriting it.
    expect_equal(names(out), c("dup", "other", "dup.1"))
    expect_equal(out[["dup"]], "QTL")
    expect_equal(out[["dup.1"]], "GWAS")
})

test_that(".cbConcat answers an empty list for no pieces", {
    expect_identical(pecotmr:::.cbConcat(list()), list())
})

test_that("colocboostPipeline forwards residualizationParam to both sides", {
    # It accepted no bundle at all, so it ran on the dataset's own settings
    # with no per-call say -- unlike fineMappingPipeline() and
    # twasWeightsPipeline(). Both SIDES matter here: colocboost consumes X and
    # Y together, so a bundle that reached only the genotypes would
    # residualize the two against different covariate sets.
    data(qtlDatasetExample, envir = environment())
    origG <- pecotmr:::.qtlResidualizedGenotypes
    origP <- pecotmr:::.qtlResidualizedPhenotypes
    seen <- new.env(parent = emptyenv())
    note <- function(which, ra) {
        assign(
            which,
            c(
                get0(which, seen, ifnotfound = character(0)),
                as.character(ra$scaleResiduals %||% "unset")
            ),
            envir = seen
        )
    }
    local_mocked_bindings(
        .qtlResidualizedGenotypes = function(
            x,
            ...,
            residualizationParam = ResidualizationParam()
        ) {
            note("geno", residualizationParam)
            origG(x, ..., residualizationParam = residualizationParam)
        },
        .qtlResidualizedPhenotypes = function(
            x,
            ...,
            residualizationParam = ResidualizationParam()
        ) {
            note("pheno", residualizationParam)
            origP(x, ..., residualizationParam = residualizationParam)
        },
        .package = "pecotmr"
    )
    suppressWarnings(suppressMessages(colocboostPipeline(
        qtlDatasetExample,
        cisWindow = 1e6,
        residualizationParam = ResidualizationParam(scaleResiduals = FALSE)
    )))
    expect_equal(unique(get("geno", envir = seen)), "FALSE")
    expect_equal(unique(get("pheno", envir = seen)), "FALSE")
})

# ===========================================================================
# One outcome per GWAS study; the focal outcome
# ===========================================================================

# One GWAS study over three rows: two LD blocks on chr1 (100-300, 400-500)
# and one on chr2.
.cbp_makeBlockedGwas <- function() {
    gr <- GenomicRanges::GRanges(
        c(rep("chr1", 5L), "chr2"),
        IRanges::IRanges(c(100L * (1:5), 100L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = c(sprintf("chr1:%d:A:G", 100L * (1:5)), "chr2:100:A:G"),
        A1 = rep("A", 6L),
        A2 = rep("G", 6L),
        Z = seq_len(6L) / 10,
        N = rep(1000L, 6L)
    )
    GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .cbp_makeHandle(),
        qcOptions = list(ok = 1),
        region = GenomicRanges::GRanges(
            c("chr1", "chr1", "chr2"),
            IRanges::IRanges(c(1L, 301L, 1L), c(300L, 1000L, 1000L))
        )
    )
}

# Two GWAS studies over the fixtures' five chr1 variants.
.cbp_makeTwoGwas <- function() {
    gr <- .cbp_makeGwasSumStats()[[1L]]
    GwasSumStats(
        studyName = c("G1", "G2"),
        entry = list(gr, gr),
        genome = "hg19",
        ldReference = .cbp_makeHandle(),
        qcOptions = list(ok = 1)
    )
}

test_that("a GWAS study is one outcome over the rows the QTL data reach", {
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    gws <- .cbp_makeBlockedGwas()
    expect_equal(nrow(gws), 3L)
    chr1 <- sprintf("chr1:%d:A:G", 100L * (1:5))
    # A window reaching both chr1 blocks gets their union as one outcome; the
    # chr2 row is not read.
    both <- pecotmr:::.cbGwasSumStatsBundle(
        gws,
        window = GenomicRanges::GRanges("chr1", IRanges::IRanges(150L, 450L))
    )
    expect_named(both, "G1")
    expect_setequal(both$G1$sumstat$variant, chr1)
    expect_equal(both$G1$study, "G1")
    # A window inside one block reads that block alone.
    one <- pecotmr:::.cbGwasSumStatsBundle(
        gws,
        window = GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, 200L))
    )
    expect_setequal(one$G1$sumstat$variant, chr1[1:3])
    # A study none of whose rows the QTL data reach is left out, and said so.
    expect_message(
        none <- pecotmr:::.cbGwasSumStatsBundle(
            gws,
            window = GenomicRanges::GRanges("chr3", IRanges::IRanges(1L, 9L))
        ),
        "leaving out GWAS G1"
    )
    expect_length(none, 0L)
})

test_that("a GWAS study's overlapped rows must share one LD reference", {
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    gws <- pecotmr:::.withRowLdRefs(
        .cbp_makeBlockedGwas(),
        list(.cbp_makeHandle(), .cbp_makeLdMatrix(), .cbp_makeHandle())
    )
    expect_error(
        pecotmr:::.cbGwasSumStatsBundle(
            gws,
            window = GenomicRanges::GRanges(
                "chr1",
                IRanges::IRanges(150L, 450L)
            )
        ),
        "GWAS study 'G1' the QTL data overlap use 2 different LD references"
    )
})

test_that(".cbAssertFocalOutcome takes a study name or a named QTL vector", {
    ok <- list(
        NULL,
        NA,
        "G1",
        c(trait = "ENSG1"),
        c(study = "S1", context = "brain", trait = "ENSG1")
    )
    for (x in ok) {
        expect_invisible(pecotmr:::.cbAssertFocalOutcome(x))
    }
    bad <- list(
        1,
        c("G1", "G2"),
        c(gene = "ENSG1"),
        c("brain", trait = "ENSG1"),
        c(trait = NA_character_),
        character(0),
        list(trait = "ENSG1")
    )
    for (x in bad) {
        expect_error(pecotmr:::.cbAssertFocalOutcome(x), "focalOutcome")
    }
})

test_that(".cbResolveFocal: the lone GWAS by default, a named one, or none", {
    one <- .cbp_makeGwasSumStats()
    two <- .cbp_makeTwoGwas()
    qtlInfo <- tibble::tibble(
        name = c("brain_ENSG1", "liver_ENSG1", "Q1:c1:t1"),
        study = c("S1", "S1", "Q1"),
        context = c("brain", "liver", "c1"),
        trait = c("ENSG1", "ENSG1", "t1")
    )
    resolve <- function(focal, gws, keys) {
        pecotmr:::.cbResolveFocal(focal, gws, qtlInfo, keys)
    }
    keys <- c(G1 = "G1", G2 = "G2")
    # One study: the GWAS is the joint run's focal outcome, unless NA.
    expect_equal(resolve(NULL, one, keys[1]), list(xqtl = NULL, joint = "G1"))
    expect_equal(resolve(NA, one, keys[1]), list(xqtl = NULL, joint = NULL))
    # ...and only when it reached the runs.
    expect_null(resolve(NULL, one, character(0))$joint)
    # Several studies: none, unless one is named.
    expect_equal(resolve(NULL, two, keys), list(xqtl = NULL, joint = NULL))
    expect_equal(resolve("G2", two, keys)$joint, "G2")
    expect_error(resolve("G3", two, keys), "names no GWAS study")
    expect_error(resolve("G1", NULL, character(0)), "no gwasSumStats")
    expect_error(resolve("G2", two, keys[1]), "none of its rows overlap")
    # A QTL outcome, by as many parts as it takes, is focal in both runs.
    expect_equal(
        resolve(c(context = "brain", trait = "ENSG1"), two, keys),
        list(xqtl = "brain_ENSG1", joint = "brain_ENSG1")
    )
    expect_equal(resolve(c(trait = "t1"), two, keys)$xqtl, "Q1:c1:t1")
    expect_error(resolve(c(trait = "ENSG1"), two, keys), "matches 2 QTL")
    expect_error(resolve(c(trait = "nope"), two, keys), "matches no QTL")
})

test_that("the joint run is focal on a lone GWAS study by default", {
    seen <- list()
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        colocboost = function(...) {
            seen[[length(seen) + 1L]] <<- list(...)
            list(stub = TRUE)
        },
        .package = "colocboost"
    )
    # Built once: the runs below differ only in their focal outcome.
    qtl <- .cbp_makeQtlSumStats()
    oneGwas <- .cbp_makeGwasSumStats()
    twoGwas <- .cbp_makeTwoGwas()
    run <- function(gws, ...) {
        seen <<- list()
        suppressMessages(colocboostPipeline(
            qtl,
            gwasSumStats = gws,
            xqtlColoc = FALSE,
            jointGwas = TRUE,
            ...
        ))
        seen[[1L]]
    }
    call <- run(oneGwas)
    expect_equal(call$outcome_names, c("Q1:c1:t1", "G1"))
    expect_equal(call$focal_outcome_idx, 2L)
    off <- run(oneGwas, focalOutcome = NA)
    expect_null(off$focal_outcome_idx)
    toQtl <- run(oneGwas, focalOutcome = c(trait = "t1"))
    expect_equal(toQtl$focal_outcome_idx, 1L)
    # Two studies: non-focal unless one is named.
    expect_null(run(twoGwas)$focal_outcome_idx)
    named <- run(twoGwas, focalOutcome = "G2")
    expect_equal(named$outcome_names[[named$focal_outcome_idx]], "G2")
})

test_that("a QTL sumstat outcome records the tuple it came from", {
    pairs <- list(
        `Q1:c1:t1` = list(study = "Q1", context = "c1", trait = "t1"),
        G1 = list(study = "G1")
    )
    bundle <- list(sumstat = list(`Q1:c1:t1` = 1, G1 = 2))
    info <- pecotmr:::.cbOutcomeInfo(
        NULL,
        bundle,
        FALSE,
        c(`Q1:c1:t1` = "ld1", G1 = "ld1"),
        pecotmr:::.cbPairsInfo(pairs)
    )
    expect_equal(info$study, c("Q1", "G1"))
    expect_equal(info$context, c("c1", NA))
    expect_equal(info$trait, c("t1", NA))
})

test_that("a focalOutcome vector with no QTL outcomes to match says so", {
    noQtl <- tibble::tibble(
        name = character(0),
        study = character(0),
        context = character(0),
        trait = character(0)
    )
    expect_error(
        pecotmr:::.cbResolveFocal(
            c(trait = "t1"),
            .cbp_makeGwasSumStats(),
            noQtl,
            c(G1 = "G1")
        ),
        "The QTL outcomes are: none"
    )
})

test_that("an outcome whose variants its LD reference lacks is left out", {
    # The fixtures' LdMatrix holds chr1:100-500; these sit at 1100-1500.
    ids <- sprintf("chr1:%d:A:G", 1100L + 100L * (0:4))
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(1100L + 100L * (0:4), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = ids,
        A1 = rep("A", 5L),
        A2 = rep("G", 5L),
        Z = seq_len(5L) / 10,
        N = rep(1000L, 5L)
    )
    qss <- QtlSumStats(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .cbp_makeLdMatrix(),
        qcOptions = list(ok = 1)
    )
    expect_null(pecotmr:::.cbQtlEntryPair(1L, qss, NULL))
    gws <- GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .cbp_makeLdMatrix(),
        qcOptions = list(ok = 1)
    )
    expect_length(pecotmr:::.cbGwasSumStatsBundle(gws), 0L)
})

test_that(".cbQtlWindow places the QTL variants it can parse", {
    asX <- function(ids) {
        list(X = list(matrix(0, 2L, length(ids), dimnames = list(NULL, ids))))
    }
    # Nothing placed -- no QTL side, or only rsIDs -- leaves every GWAS row in.
    expect_null(pecotmr:::.cbQtlWindow(NULL, list()))
    expect_null(pecotmr:::.cbQtlWindow(asX(c("rs1", "rs2")), list()))
    window <- pecotmr:::.cbQtlWindow(
        asX(c("1:100:A:G", "rs2", "chr1:400:C:T")),
        list()
    )
    expect_equal(as.character(GenomicRanges::seqnames(window)), "chr1")
    expect_equal(GenomicRanges::start(window), 100L)
    expect_equal(GenomicRanges::end(window), 400L)
})

test_that(".cbQtlVariants places the variants every QTL outcome reads", {
    asX <- function(ids) {
        list(X = list(matrix(0, 2L, length(ids), dimnames = list(NULL, ids))))
    }
    expect_null(pecotmr:::.cbQtlVariants(NULL, list()))
    expect_null(pecotmr:::.cbQtlVariants(asX("rs1"), list()))
    placed <- pecotmr:::.cbQtlVariants(
        asX(c("1:100:A:G", "rs2")),
        list(q1 = list(variantIds = c("chr1:100:A:G", "chr1:700:C:T")))
    )
    expect_identical(as.character(placed), c("chr1:100", "chr1:700"))
})

test_that("colocboost refuses QTL data reaching past a given block", {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(100L * (1:5), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:5)),
        A1 = "A",
        A2 = "G",
        Z = 1,
        N = 100L
    )
    gws <- GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        qcOptions = list(ok = 1),
        region = GenomicRanges::GRanges("chr1:1-550")
    )
    qtl <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(300L, 600L), width = 1L)
    )
    expect_error(
        pecotmr:::.cbAppendGwasPairs(list(), gws, qtlVariants = qtl),
        "the QTL data: 1 variant(s) past the end of block chr1:1-550",
        fixed = TRUE
    )
    # Without QTL variants to place there is nothing to check.
    expect_null(pecotmr:::.cbCheckBlockReach(NULL, gws))
})

test_that("colocboostPipeline(QtlSumStats) lays a cisWindow around traitPos", {
    expect_error(
        colocboostPipeline(.cbp_makeQtlSumStats(), cisWindow = 1000),
        "colocboostPipeline: `cisWindow` is laid around each tuple's trait",
        fixed = TRUE
    )
    # Naming a context reaches the same check: the window is laid on the
    # tuples colocboost's own selection keeps.
    expect_error(
        colocboostPipeline(
            .cbp_makeQtlSumStats(),
            contexts = "c1",
            cisWindow = 1000
        ),
        "colocboostPipeline: `cisWindow` is laid around each tuple's trait",
        fixed = TRUE
    )
})

test_that("colocboostPipeline(MSQD) lays cisWindow on its summary statistics", {
    local_mocked_bindings(
        extractBlockGenotypes = .cbp_mockExtractor(),
        .package = "pecotmr"
    )
    ss <- .cbp_makeQtlSumStats()
    ss$traitPos <- GenomicRanges::GRanges("chr1:300")
    mt <- MultiStudyQtlDataset(
        qtlDatasets = list(
            study1 = .cbp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
        ),
        sumStats = ss
    )
    pairs <- pecotmr:::.cbMultiStudySumstats(mt, NULL, NULL, 100)$qtlPairs
    expect_identical(as.character(pairs[[1L]]$span), "chr1:200-400")
    # A summary-statistics study without a trait position is an error, as it
    # is on its own.
    expect_error(
        colocboostPipeline(.cbp_makeMultiStudy(), cisWindow = 100),
        "colocboostPipeline: `cisWindow` is laid around each tuple's trait",
        fixed = TRUE
    )
})

# ===========================================================================
# .cbToResultObject
# ===========================================================================

test_that("separate_gwas runs are keyed by gwasStudy", {
    # The three analyses are not parallel: separate_gwas is a named list one
    # level deeper, so without gwasStudy its sets would be indistinguishable.
    raw <- list(
        xqtl_coloc = NULL,
        joint_gwas = NULL,
        separate_gwas = list(G1 = .cbr_fake(), G2 = .cbr_fake()),
        computing_time = list()
    )
    x <- .cbToResultObject(raw, .cbr_info())
    expect_equal(nrow(x), 2L)
    expect_equal(x$analysis, rep("separate_gwas", 2L))
    expect_setequal(x$gwasStudy, c("G1", "G2"))
})
