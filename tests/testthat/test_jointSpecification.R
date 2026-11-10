# Tests for the joint-specification grammar and input-argument parsers
# (R/jointSpecification.R).

# The (X, Y) a joint fit is on, as the engine makes them: its phenotypes first,
# then, if it can be fit at all, its genotypes on the samples they leave.
.js_residGeno <- function(
    data,
    contexts,
    traits,
    cisWindow,
    region,
    residualizationParam,
    genotypeFilter,
    samples
) {
    if (is.null(region)) {
        return(pecotmr:::.fmResidGeno(
            data,
            contexts = contexts,
            traitId = traits,
            cisWindow = cisWindow,
            samples = samples,
            residualizationParam = residualizationParam,
            genotypeFilter = genotypeFilter
        ))
    }
    pecotmr:::.fmResidGeno(
        data,
        contexts = contexts,
        region = region,
        samples = samples,
        residualizationParam = residualizationParam,
        genotypeFilter = genotypeFilter
    )
}

.js_crossContextXy <- function(
    data,
    tid,
    scopedContexts,
    cisWindow,
    verbose,
    label,
    region = NULL,
    residualizationParam = ResidualizationParam(),
    genotypeFilter = pecotmr:::.qtlResolveFilter(),
    phenotypeFilter = NULL
) {
    planned <- pecotmr:::.crossContextPlan(
        data,
        tid,
        scopedContexts,
        verbose,
        label,
        residualizationParam,
        phenotypeFilter
    )
    if (is.null(planned)) {
        return(NULL)
    }
    pecotmr:::.crossContextAssemble(
        .js_residGeno(
            data,
            planned$perTraitContexts,
            tid,
            cisWindow,
            region,
            residualizationParam,
            genotypeFilter,
            planned$fitSamples
        ),
        planned$Yres,
        planned$perTraitContexts,
        planned$fitSamples,
        planned$report
    )
}

.js_crossTraitXy <- function(
    data,
    cx,
    scopedTraits,
    cisWindow,
    verbose,
    label,
    studyName,
    region = NULL,
    residualizationParam = ResidualizationParam(),
    genotypeFilter = pecotmr:::.qtlResolveFilter(),
    phenotypeFilter = NULL
) {
    planned <- pecotmr:::.crossTraitPlan(
        data,
        cx,
        scopedTraits,
        studyName,
        verbose,
        label,
        residualizationParam,
        phenotypeFilter
    )
    if (is.null(planned)) {
        return(NULL)
    }
    pecotmr:::.crossTraitAssemble(
        .js_residGeno(
            data,
            cx,
            planned$traitsHere,
            cisWindow,
            region,
            residualizationParam,
            genotypeFilter,
            planned$fitSamples
        ),
        planned$allY,
        planned$traitsHere,
        planned$se
    )
}

.js_composedXy <- function(
    data,
    scope,
    studyName,
    cisWindow,
    verbose,
    label,
    region = NULL,
    residualizationParam = ResidualizationParam(),
    genotypeFilter = pecotmr:::.qtlResolveFilter(),
    phenotypeFilter = NULL
) {
    planned <- pecotmr:::.composedPlan(
        data,
        scope,
        studyName,
        verbose,
        label,
        residualizationParam,
        phenotypeFilter
    )
    if (is.null(planned)) {
        return(NULL)
    }
    pecotmr:::.composedAlignXy(
        list(
            X = .js_residGeno(
                data,
                planned$allContexts,
                planned$allTraits,
                cisWindow,
                region,
                residualizationParam,
                genotypeFilter,
                planned$fitSamples
            ),
            YresList = planned$YresList
        ),
        planned$tuples
    )
}

# -----------------------------------------------------------------------------
# Fixture builders
# -----------------------------------------------------------------------------

.js_makeGenotypeHandle <- function(snp_n = 5L) {
    new(
        "GenotypeHandle",
        path = "/tmp/test",
        format = "plink2",
        snpInfo = data.frame(
            SNP = paste0("rs", seq_len(snp_n)),
            CHR = rep("chr1", snp_n),
            BP = seq(100L, by = 100L, length.out = snp_n),
            A1 = rep("A", snp_n),
            A2 = rep("G", snp_n),
            stringsAsFactors = FALSE
        ),
        nSamples = 5L,
        sampleIds = paste0("s", seq_len(5)),
        pgenPtr = NULL
    )
}

.js_makeSe <- function(
    traits = c("ENSG1", "ENSG2"),
    samples = paste0("s", seq_len(5))
) {
    rng <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq.int(100L, by = 100L, length.out = length(traits)),
            width = 50L
        )
    )
    names(rng) <- traits
    expr <- matrix(
        rnorm(length(traits) * length(samples)),
        nrow = length(traits),
        ncol = length(samples),
        dimnames = list(traits, samples)
    )
    cd <- S4Vectors::DataFrame(
        sex = rep(c("M", "F"), length.out = length(samples)),
        row.names = samples
    )
    SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = expr),
        rowRanges = rng,
        colData = cd
    )
}

.js_makeQtlDataset <- function(
    studyName = "s1",
    contexts = c("brain", "liver"),
    traits = c("ENSG1", "ENSG2")
) {
    phenos <- setNames(
        map(contexts, function(cx) .js_makeSe(traits = traits)),
        contexts
    )
    QtlDataset(
        studyName = studyName,
        genotypes = .js_makeGenotypeHandle(),
        phenotypes = phenos,
        genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0)
    )
}

.js_makeQtlSumStats <- function(
    studies = c("ssA", "ssB"),
    contexts = c("DLPFC"),
    traits = c("ENSG3")
) {
    rows <- expand.grid(
        study = studies,
        context = contexts,
        trait = traits,
        stringsAsFactors = FALSE
    )
    entries <- map(seq_len(nrow(rows)), function(i) {
        gr <- GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(start = c(100L, 200L), width = 1L)
        )
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
            SNP = c("rs1", "rs2"),
            A1 = c("A", "G"),
            A2 = c("G", "A"),
            Z = c(1.5, -1.0),
            N = c(1000L, 1000L)
        )
        gr
    })
    QtlSumStats(
        studyName = rows$study,
        context = rows$context,
        trait = rows$trait,
        entry = entries,
        genome = "hg19",
        ldReference = .js_makeGenotypeHandle()
    )
}

# The default dataset, built on first use and shared by the parser tests
# that only read its names.
delayedAssign(".jsQd", .js_makeQtlDataset())

# -----------------------------------------------------------------------------
# Scope helpers
# -----------------------------------------------------------------------------

test_that("a merged row's window spans its known windows, on one chromosome", {
    unknown <- GenomicRanges::GRanges(
        c("chrUn", "chrUn"),
        IRanges::IRanges(c(1L, 1L), width = 0L)
    )
    expect_identical(pecotmr:::.twasWindowSpan(1:2, unknown), unknown[1L])
    twoChrom <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(c(1L, 1L), width = 10L)
    )
    expect_error(
        pecotmr:::.twasWindowSpan(1:2, twoChrom),
        "windows on 2 chromosomes"
    )
})

test_that(".spListStudies / .spListContexts / .spListTraits: QtlDataset", {
    qd <- .js_makeQtlDataset(
        studyName = "S1",
        contexts = c("brain", "liver"),
        traits = c("ENSG_A", "ENSG_B")
    )
    expect_equal(pecotmr:::.spListStudies(qd), "S1")
    expect_setequal(pecotmr:::.spListContexts(qd), c("brain", "liver"))
    expect_setequal(
        pecotmr:::.spListTraits(qd, studyName = "S1", context = "brain"),
        c("ENSG_A", "ENSG_B")
    )
    expect_equal(pecotmr:::.spStudyDataForm(qd, "S1"), "individual")
})

test_that(".sp* helpers: MultiStudyQtlDataset combines individual + sumstats", {
    qd1 <- .js_makeQtlDataset(
        studyName = "indA",
        contexts = "brain",
        traits = "ENSG_A"
    )
    qd2 <- .js_makeQtlDataset(
        studyName = "indB",
        contexts = "liver",
        traits = "ENSG_A"
    )
    ss <- .js_makeQtlSumStats(
        studies = "ssC",
        contexts = "DLPFC",
        traits = "ENSG_A"
    )
    mt <- MultiStudyQtlDataset(
        qtlDatasets = list(indA = qd1, indB = qd2),
        sumStats = ss
    )
    expect_setequal(pecotmr:::.spListStudies(mt), c("indA", "indB", "ssC"))
    expect_setequal(pecotmr:::.spListContexts(mt), c("brain", "liver", "DLPFC"))
    expect_equal(pecotmr:::.spStudyDataForm(mt, "indA"), "individual")
    expect_equal(pecotmr:::.spStudyDataForm(mt, "ssC"), "sumstats")
    expect_error(pecotmr:::.spStudyDataForm(mt, "missing"), "not in")
})

# -----------------------------------------------------------------------------
# parseJointSpecification
# -----------------------------------------------------------------------------

test_that("parseJointSpecification: NULL returns empty list", {
    qd <- .jsQd
    expect_equal(pecotmr:::parseJointSpecification(NULL, qd), list())
    expect_silent(pecotmr:::parseJointSpecification(NULL, qd))
})

test_that("parseJointSpecification: auto-wraps a single char vector", {
    qd <- .jsQd
    out <- pecotmr:::parseJointSpecification("context", qd)
    expect_equal(length(out), 1L)
    expect_equal(out[[1L]]$axes, "context")
    expect_null(out[[1L]]$scope)
})

test_that("parseJointSpecification: accepts a list of specs", {
    qd <- .jsQd
    out <- pecotmr:::parseJointSpecification(
        list("context", c("context", "trait")),
        qd
    )
    expect_equal(length(out), 2L)
    expect_equal(out[[2L]]$axes, c("context", "trait"))
})

test_that("parseJointSpecification: accepts scope-restricted spec", {
    qd1 <- .js_makeQtlDataset(studyName = "A")
    qd2 <- .js_makeQtlDataset(studyName = "B")
    mt <- MultiStudyQtlDataset(qtlDatasets = list(A = qd1, B = qd2))
    out <- pecotmr:::parseJointSpecification(
        list(list(axes = "context", scope = list(study = c("A")))),
        mt
    )
    expect_equal(out[[1L]]$scope$study, "A")
})

test_that("parseJointSpecification: rejects unknown axes", {
    qd <- .jsQd
    expect_error(
        pecotmr:::parseJointSpecification(c("context", "bogus"), qd),
        "unknown axes"
    )
})

test_that("parseJointSpecification: rejects duplicate axes", {
    qd <- .jsQd
    expect_error(
        pecotmr:::parseJointSpecification(c("context", "context"), qd),
        "duplicate axes"
    )
})

test_that("parseJointSpecification: rejects unknown scope keys", {
    qd <- .jsQd
    expect_error(
        pecotmr:::parseJointSpecification(
            list(list(axes = "context", scope = list(bogus = "x"))),
            qd
        ),
        "unknown scope key"
    )
})

test_that("parseJointSpecification: rejects scope values absent from data", {
    qd <- .js_makeQtlDataset(studyName = "S1")
    expect_error(
        pecotmr:::parseJointSpecification(
            list(list(axes = "context", scope = list(study = "NotPresent"))),
            qd
        ),
        "scope\\$study contains values not in data"
    )
})

test_that("parseJointSpecification: rejects unknown spec elements", {
    qd <- .jsQd
    expect_error(
        pecotmr:::parseJointSpecification(
            list(list(axes = "context", scope = NULL, bogus = "x")),
            qd
        ),
        "unknown element"
    )
})

# -----------------------------------------------------------------------------
# parseContexts
# -----------------------------------------------------------------------------

test_that("parseContexts: NULL passes through", {
    qd <- .jsQd
    expect_null(pecotmr:::parseContexts(NULL, qd))
})

test_that("parseContexts: vector intersects with each study's contexts", {
    qd1 <- .js_makeQtlDataset(studyName = "A", contexts = c("brain", "liver"))
    qd2 <- .js_makeQtlDataset(studyName = "B", contexts = c("brain"))
    mt <- MultiStudyQtlDataset(qtlDatasets = list(A = qd1, B = qd2))
    expect_warning(
        out <- pecotmr:::parseContexts(c("brain", "liver"), mt),
        "B.*liver"
    )
    expect_setequal(out$A, c("brain", "liver"))
    expect_equal(out$B, "brain")
})

test_that("parseContexts: named-list form requires valid studies", {
    qd <- .js_makeQtlDataset(studyName = "A")
    expect_error(
        pecotmr:::parseContexts(list(B = "brain"), qd),
        "unknown studies"
    )
})

test_that("parseContexts: named-list rejects unknown contexts", {
    qd <- .js_makeQtlDataset(studyName = "A", contexts = c("brain", "liver"))
    expect_error(
        pecotmr:::parseContexts(list(A = "bogus"), qd),
        "unknown contexts"
    )
})

test_that("parseContexts: list form fills unmentioned studies with defaults", {
    qd1 <- .js_makeQtlDataset(studyName = "A", contexts = c("brain", "liver"))
    qd2 <- .js_makeQtlDataset(studyName = "B", contexts = c("brain", "liver"))
    mt <- MultiStudyQtlDataset(qtlDatasets = list(A = qd1, B = qd2))
    out <- pecotmr:::parseContexts(list(A = "brain"), mt)
    expect_equal(out$A, "brain")
    expect_setequal(out$B, c("brain", "liver"))
})

# -----------------------------------------------------------------------------
# parseTraitIds
# -----------------------------------------------------------------------------

test_that("parseTraitIds: NULL passes through", {
    qd <- .jsQd
    expect_null(pecotmr:::parseTraitIds(NULL, qd))
})

test_that("parseTraitIds: vector form returns the vector", {
    qd <- .js_makeQtlDataset(traits = c("X", "Y"))
    expect_equal(pecotmr:::parseTraitIds(c("X", "Y"), qd), c("X", "Y"))
})

test_that("parseTraitIds: study-keyed list validates per study", {
    qd1 <- .js_makeQtlDataset(studyName = "A", traits = c("ENSG_A"))
    qd2 <- .js_makeQtlDataset(studyName = "B", traits = c("ENSG_B"))
    mt <- MultiStudyQtlDataset(qtlDatasets = list(A = qd1, B = qd2))
    expect_error(
        pecotmr:::parseTraitIds(list(A = "ENSG_B"), mt),
        "unknown traits"
    )
    out <- pecotmr:::parseTraitIds(list(A = "ENSG_A", B = "ENSG_B"), mt)
    expect_equal(out$A, "ENSG_A")
    expect_equal(out$B, "ENSG_B")
})

test_that("parseTraitIds: doubly-nested study->context validates per context", {
    qd <- .js_makeQtlDataset(
        studyName = "A",
        contexts = c("brain", "liver"),
        traits = c("ENSG_A", "ENSG_B")
    )
    out <- pecotmr:::parseTraitIds(
        list(A = list(brain = "ENSG_A", liver = "ENSG_B")),
        qd
    )
    expect_equal(out$A$brain, "ENSG_A")
    expect_equal(out$A$liver, "ENSG_B")
    expect_error(
        pecotmr:::parseTraitIds(list(A = list(bogus = "ENSG_A")), qd),
        "unknown contexts"
    )
})

# -----------------------------------------------------------------------------
# Pipeline-level wiring: jointSpecification accepted on all three methods
# -----------------------------------------------------------------------------

test_that("fineMappingPipeline(QtlDataset): composed joint dispatcher wired", {
    # Composed dispatcher accepts axes = c("context", "trait") for individual-
    # level input. Reaching the composed builder proves the dispatcher was
    # invoked (rather than the previous stub error).
    qd <- .js_makeQtlDataset(
        studyName = "A",
        contexts = c("brain", "liver"),
        traits = c("ENSG_A", "ENSG_B")
    )
    local_mocked_bindings(
        .composedPlan = function(...) abort("composed builder reached"),
        .package = "pecotmr"
    )
    expect_error(
        fineMappingPipeline(
            qd,
            methods = "mvsusie",
            cisWindow = 1000L,
            jointSpecification = list(c("context", "trait"))
        ),
        "composed builder reached"
    )
})

test_that("fineMappingPipeline(QtlDataset): composed joint refuses 'study'", {
    qd <- .js_makeQtlDataset(
        studyName = "A",
        contexts = c("brain", "liver"),
        traits = "ENSG_A"
    )
    expect_error(
        fineMappingPipeline(
            qd,
            methods = "mvsusie",
            cisWindow = 1000L,
            jointSpecification = list(c("study", "context"))
        ),
        "axes including 'study' require sumstats"
    )
})

test_that("fineMappingPipeline(QtlDataset): study axis on individual errors", {
    qd <- .js_makeQtlDataset(
        studyName = "A",
        contexts = c("brain", "liver"),
        traits = "ENSG_A"
    )
    expect_error(
        fineMappingPipeline(
            qd,
            methods = "mvsusie",
            cisWindow = 1000L,
            jointSpecification = "study"
        ),
        "needs more than one study"
    )
})

test_that("fineMappingPipeline(QtlDataset): bad jointSpec errors before fit", {
    qd <- .js_makeQtlDataset(studyName = "A")
    expect_error(
        fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            jointSpecification = "bogus_axis"
        ),
        "unknown axes"
    )
})

test_that("twasWeightsPipeline(QtlDataset): study axis on individual errors", {
    qd <- .js_makeQtlDataset(
        studyName = "A",
        contexts = c("brain", "liver"),
        traits = "ENSG_A"
    )
    expect_error(
        twasWeightsPipeline(
            qd,
            methods = "mrmash",
            cisWindow = 1000L,
            jointSpecification = "study"
        ),
        "needs more than one study"
    )
})

# ============================================================================
# Joint dispatchers (merged from former tests/testthat/test_jointDispatchers.R)
# ============================================================================

context("joint dispatchers (fineMappingPipeline / twasWeightsPipeline)")

# ============================================================================
# Strategy: each joint-dispatcher function is exercised by driving
# fineMappingPipeline / twasWeightsPipeline through the user-facing
# `jointSpecification` argument and mocking the underlying fitters
# (mvsusieRss, mrmashWeights, mrmashRssWeights, ...). The mocks return tiny
# stub objects so postprocessing builds plausible result rows.
# ============================================================================

# -----------------------------------------------------------------------------
# Fixture builders
# -----------------------------------------------------------------------------

.jd_makeHandle <- function(snp_n = 5L, n_samples = 30L) {
    new(
        "GenotypeHandle",
        path = "/tmp/jd",
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

.jd_mockExtractor <- function(seed = 11, n_samples = 30L) {
    function(handle, snpIdx, sampleIdx = NULL) {
        set.seed(seed)
        panel <- matrix(
            rbinom(n_samples * nrow(snpInfo(handle)), 2, 0.3),
            nrow = n_samples,
            ncol = nrow(snpInfo(handle))
        )
        panel[sampleIdx %||% seq_len(nrow(panel)), snpIdx, drop = FALSE]
    }
}

# Multi-(study, context, trait) QtlSumStats. Every row carries the same SNP
# order (5 variants) so jointCrossContext / jointCrossTrait / jointCrossStudy
# can stack Z columns without alignment problems.
.jd_makeQtlSumStats <- function(
    studies = "Q1",
    contexts = c("c1", "c2"),
    traits = "t1"
) {
    rows <- expand.grid(
        study = studies,
        context = contexts,
        trait = traits,
        stringsAsFactors = FALSE
    )
    makeGr <- function() {
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
        gr
    }
    QtlSumStats(
        studyName = rows$study,
        context = rows$context,
        trait = rows$trait,
        entry = map(seq_len(nrow(rows)), function(.) makeGr()),
        genome = "hg19",
        ldReference = .jd_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
}

# -----------------------------------------------------------------------------
# Mocks for SuSiE / mvsusie / mr.mash families
# -----------------------------------------------------------------------------

.jd_mockMvsusie <- function() {
    function(X, Y, prior_variance, coverage, standardize = TRUE) {
        list(token = "mvsusie", n_X_cols = ncol(X), n_Y_cols = ncol(Y))
    }
}

.jd_mockMvsusieRss <- function() {
    # The real mvsusie_rss passes its `...` on to mvsusie, which is how L and
    # L_greedy (the pipeline's CredibleSetParam) arrive; recorded for tests.
    function(Z, R, N, prior_variance, coverage, L = NULL, L_greedy = NULL) {
        list(
            token = "mvsusieRss",
            nVariants = nrow(Z),
            nOutcomes = ncol(Z),
            L = L,
            L_greedy = L_greedy
        )
    }
}

.jd_mockMixturePrior <- function() {
    function(R, ...) list(R = R)
}

# A stub postprocessor that returns a tiny FineMappingRow. Mirrors the
# `.fmp_mockPostprocess` shape from test_fineMappingPipeline.R.
.jd_mockPostprocess <- function() {
    function(
        fit,
        method,
        dataX,
        dataY,
        coverage,
        secondaryCoverage,
        signalCutoff,
        minAbsCorr,
        csInput = NULL,
        af = NULL,
        region = NULL,
        conditionIdx = NULL,
        ...
    ) {
        if (is.matrix(dataX)) {
            vids <- colnames(dataX)
        } else if (is.list(dataY) && !is.null(dataY$z)) {
            vids <- names(dataY$z)
        } else {
            vids <- "chr1:1200:A:G"
        }
        if (is.null(vids)) {
            vids <- "chr1:1200:A:G"
        }
        FineMappingRow(
            variantIds = vids,
            susieFit = list(method = method, payload = fit),
            topLoci = data.frame(
                variant_id = vids,
                pip = seq(0.9, by = -0.1, length.out = length(vids)),
                stringsAsFactors = FALSE
            )
        )
    }
}

.jd_mockMrmashWeights <- function() {
    function(X, Y, ...) {
        w <- matrix(
            0,
            nrow = ncol(X),
            ncol = ncol(Y),
            dimnames = list(colnames(X), colnames(Y))
        )
        w
    }
}

.jd_mockMrmashRssWeights <- function() {
    function(stat, LD, ...) {
        nCols <- if (is.matrix(stat$z)) ncol(stat$z) else 1L
        nVars <- if (is.matrix(stat$z)) nrow(stat$z) else length(stat$z)
        w <- matrix(0, nrow = nVars, ncol = nCols)
        rownames(w) <- if (is.matrix(stat$z)) {
            rownames(stat$z)
        } else {
            stat$variantNames
        }
        if (is.matrix(stat$z) && !is.null(colnames(stat$z))) {
            colnames(w) <- colnames(stat$z)
        }
        w
    }
}

# =============================================================================
# fineMappingPipeline joint dispatch: QtlSumStats
# =============================================================================

test_that("fineMappingPipeline(QtlSumStats): context joint per study x trait", {
    ss <- .jd_makeQtlSumStats(
        studies = "Q1",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        .fmPostprocessOne = .jd_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie_rss = .jd_mockMvsusieRss(),
        create_mixture_prior = .jd_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            ss,
            methods = "mvsusie",
            jointSpecification = "context",
            credibleSetParam = CredibleSetParam(maxNumSingleEffects = 3L)
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(nrow(res), 2L) # per-context rows
    expect_setequal(as.character(res$context), c("c1", "c2")) # REAL contexts
    expect_true(all(grepl("c1;c2|c2;c1", as.character(res$jointContexts))))
    # maxNumSingleEffects reaches the joint RSS fit, as it does susie_rss.
    expect_equal(susieFit(res[1, ])$payload$L, 3L)
})

test_that("fineMappingPipeline(QtlSumStats): one context skips context joint", {
    ss <- .jd_makeQtlSumStats(studies = "Q1", contexts = "c1", traits = "t1")
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        .fmPostprocessOne = .jd_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie_rss = .jd_mockMvsusieRss(),
        create_mixture_prior = .jd_mockMixturePrior(),
        .package = "mvsusieR"
    )
    expect_error(
        suppressMessages(
            fineMappingPipeline(
                ss,
                methods = "mvsusie",
                jointSpecification = "context"
            )
        ),
        "no joint fits produced"
    )
})

test_that("fineMappingPipeline(QtlSumStats): trait joint per study x context", {
    ss <- .jd_makeQtlSumStats(
        studies = "Q1",
        contexts = "c1",
        traits = c("t1", "t2")
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        .fmPostprocessOne = .jd_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie_rss = .jd_mockMvsusieRss(),
        create_mixture_prior = .jd_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            ss,
            methods = "mvsusie",
            jointSpecification = "trait"
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(nrow(res), 2L) # per-trait rows
    expect_setequal(as.character(res$trait), c("t1", "t2"))
})

test_that("fineMappingPipeline(QtlSumStats): trait joint refuses fsusie", {
    ss <- .jd_makeQtlSumStats(
        studies = "Q1",
        contexts = "c1",
        traits = c("t1", "t2")
    )
    expect_error(
        fineMappingPipeline(
            ss,
            methods = "fsusie",
            jointSpecification = "trait"
        ),
        "fsusie"
    )
})

test_that("fineMappingPipeline(QtlSumStats): study joint per context x trait", {
    ss <- .jd_makeQtlSumStats(
        studies = c("Q1", "Q2"),
        contexts = "c1",
        traits = "t1"
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        .fmPostprocessOne = .jd_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie_rss = .jd_mockMvsusieRss(),
        create_mixture_prior = .jd_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            ss,
            methods = "mvsusie",
            jointSpecification = "study"
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(nrow(res), 2L) # per-study rows
    expect_setequal(as.character(res$study), c("Q1", "Q2"))
})

test_that("fineMappingPipeline(QtlSumStats): study+context joint fits", {
    ss <- .jd_makeQtlSumStats(
        studies = c("Q1", "Q2"),
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        .fmPostprocessOne = .jd_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie_rss = .jd_mockMvsusieRss(),
        create_mixture_prior = .jd_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            ss,
            methods = "mvsusie",
            jointSpecification = list(c("study", "context"))
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_true("jointStudies" %in% pecotmr:::.tupleColumnNames(res))
    expect_true("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
    expect_equal(nrow(res), 4L) # study x context
    expect_setequal(as.character(res$study), c("Q1", "Q2")) # both vary -> real
    expect_setequal(as.character(res$context), c("c1", "c2"))
})

test_that("fineMappingPipeline(QtlSumStats): composed joint rejects fsusie", {
    ss <- .jd_makeQtlSumStats(
        studies = c("Q1", "Q2"),
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    expect_error(
        fineMappingPipeline(
            ss,
            methods = "fsusie",
            jointSpecification = list(c("study", "context"))
        ),
        "fsusie"
    )
})

# =============================================================================
# twasWeightsPipeline joint dispatch: QtlDataset
# =============================================================================

.jd_makeSe <- function(traits = c("t1", "t2"), n_samples = 30L, starts = NULL) {
    if (is.null(starts)) {
        starts <- seq(1000L, by = 1000L, length.out = length(traits))
    }
    rng <- GenomicRanges::GRanges(
        seqnames = rep("chr1", length(traits)),
        ranges = IRanges::IRanges(start = starts, width = 500L)
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

.jd_makeQtlDataset <- function(
    studyName = "Q1",
    contexts = c("c1", "c2"),
    traits = c("t1", "t2")
) {
    phen <- setNames(
        map(contexts, function(.) .jd_makeSe(traits = traits)),
        contexts
    )
    QtlDataset(
        studyName = studyName,
        genotypes = .jd_makeHandle(),
        phenotypes = phen,
        genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0)
    )
}

test_that("twasWeightsPipeline(QtlDataset): context joint mr.mash per trait", {
    qd <- .jd_makeQtlDataset(
        studyName = "Q1",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        mrmashWeights = .jd_mockMrmashWeights(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        twasWeightsPipeline(
            qd,
            methods = "mrmash",
            cisWindow = 1000L,
            jointSpecification = "context"
        )
    )
    expect_s4_class(res, "TwasWeights")
    expect_equal(nrow(res), 2L) # per-context rows
    expect_setequal(as.character(res$context), c("c1", "c2"))
    expect_true("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
})

test_that("twasWeightsPipeline(QtlDataset): one context skips context joint", {
    qd <- .jd_makeQtlDataset(
        studyName = "Q1",
        contexts = "c1",
        traits = c("t1", "t2")
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        mrmashWeights = .jd_mockMrmashWeights(),
        .package = "pecotmr"
    )
    expect_error(
        suppressMessages(
            twasWeightsPipeline(
                qd,
                methods = "mrmash",
                cisWindow = 1000L,
                jointSpecification = "context"
            )
        ),
        "no joint fits produced|context"
    )
})

test_that("twasWeightsPipeline(QtlDataset): trait joint mr.mash per context", {
    qd <- .jd_makeQtlDataset(
        studyName = "Q1",
        contexts = "c1",
        traits = c("t1", "t2")
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        mrmashWeights = .jd_mockMrmashWeights(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        twasWeightsPipeline(
            qd,
            methods = "mrmash",
            cisWindow = 1000L,
            jointSpecification = "trait"
        )
    )
    expect_s4_class(res, "TwasWeights")
    expect_equal(nrow(res), 2L) # per-trait rows
    expect_setequal(as.character(res$trait), c("t1", "t2"))
    expect_true("jointTraits" %in% pecotmr:::.tupleColumnNames(res))
})

test_that("twasWeightsPipeline(QtlDataset): context+trait joint fits", {
    qd <- .jd_makeQtlDataset(
        studyName = "Q1",
        contexts = c("c1", "c2"),
        traits = c("t1", "t2")
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        mrmashWeights = .jd_mockMrmashWeights(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        twasWeightsPipeline(
            qd,
            methods = "mrmash",
            cisWindow = 1000L,
            jointSpecification = list(c("context", "trait"))
        )
    )
    expect_s4_class(res, "TwasWeights")
    expect_equal(nrow(res), 4L) # context x trait
    expect_setequal(as.character(res$context), c("c1", "c2"))
    expect_setequal(as.character(res$trait), c("t1", "t2"))
})

test_that("twasWeightsPipeline(QtlDataset): composed joint with study errors", {
    qd <- .jd_makeQtlDataset(
        studyName = "Q1",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    expect_error(
        twasWeightsPipeline(
            qd,
            methods = "mrmash",
            cisWindow = 1000L,
            jointSpecification = list(c("study", "context"))
        ),
        "require sumstats|requires sumstats"
    )
})

# =============================================================================
# twasWeightsPipeline joint dispatch: QtlSumStats
# =============================================================================

test_that("twasWeightsPipeline(QtlSumStats): context joint per study x trait", {
    ss <- .jd_makeQtlSumStats(
        studies = "Q1",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        mrmashRssWeights = .jd_mockMrmashRssWeights(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        twasWeightsPipeline(
            ss,
            methods = "mrmash",
            jointSpecification = "context"
        )
    )
    expect_s4_class(res, "TwasWeights")
    expect_equal(nrow(res), 2L) # per-context rows
    expect_setequal(as.character(res$context), c("c1", "c2"))
    expect_true("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
})

test_that("twasWeightsPipeline(QtlSumStats): trait joint per study x context", {
    ss <- .jd_makeQtlSumStats(
        studies = "Q1",
        contexts = "c1",
        traits = c("t1", "t2")
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        mrmashRssWeights = .jd_mockMrmashRssWeights(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        twasWeightsPipeline(
            ss,
            methods = "mrmash",
            jointSpecification = "trait"
        )
    )
    expect_s4_class(res, "TwasWeights")
    expect_equal(nrow(res), 2L) # per-trait rows
    expect_setequal(as.character(res$trait), c("t1", "t2"))
})

test_that("twasWeightsPipeline(QtlSumStats): study joint per context x trait", {
    ss <- .jd_makeQtlSumStats(
        studies = c("Q1", "Q2"),
        contexts = "c1",
        traits = "t1"
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        mrmashRssWeights = .jd_mockMrmashRssWeights(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        twasWeightsPipeline(
            ss,
            methods = "mrmash",
            jointSpecification = "study"
        )
    )
    expect_s4_class(res, "TwasWeights")
    expect_equal(nrow(res), 2L) # per-study rows
    expect_setequal(as.character(res$study), c("Q1", "Q2"))
})

test_that("twasWeightsPipeline(QtlSumStats): study+context joint fits", {
    ss <- .jd_makeQtlSumStats(
        studies = c("Q1", "Q2"),
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    local_mocked_bindings(
        .handleDosage = .jd_mockExtractor(),
        mrmashRssWeights = .jd_mockMrmashRssWeights(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        twasWeightsPipeline(
            ss,
            methods = "mrmash",
            jointSpecification = list(c("study", "context"))
        )
    )
    expect_s4_class(res, "TwasWeights")
    expect_true("jointStudies" %in% pecotmr:::.tupleColumnNames(res))
    expect_true("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
    expect_equal(nrow(res), 4L) # study x context
    expect_setequal(as.character(res$study), c("Q1", "Q2"))
    expect_setequal(as.character(res$context), c("c1", "c2"))
})

# ============================================================================
# Additional coverage: scope helpers, parser error branches, scope resolution,
# X/Y/Z builders' skip paths, multi-region merges, and MultiStudy dispatchers.
# ============================================================================

# -----------------------------------------------------------------------------
# .sp* helpers: QtlSumStats / MultiStudy branches + unsupported-class errors
# -----------------------------------------------------------------------------

test_that(".sp* helpers: QtlSumStats study/context/trait listing + dataForm", {
    ss <- .js_makeQtlSumStats(
        studies = c("A", "B"),
        contexts = c("c1", "c2"),
        traits = c("t1", "t2")
    )
    expect_setequal(pecotmr:::.spListStudies(ss), c("A", "B"))
    expect_setequal(pecotmr:::.spListContexts(ss), c("c1", "c2")) # study=NULL
    expect_setequal(pecotmr:::.spListContexts(ss, "A"), c("c1", "c2"))
    expect_setequal(
        pecotmr:::.spListTraits(ss, studyName = "A", context = "c1"),
        c("t1", "t2")
    )
    expect_equal(pecotmr:::.spStudyDataForm(ss, "A"), "sumstats")
    expect_error(
        pecotmr:::.spStudyDataForm(ss, "missing"),
        "not in QtlSumStats"
    )
})

test_that(".spStudyDataForm: QtlDataset wrong study errors", {
    qd <- .js_makeQtlDataset(studyName = "S1")
    expect_error(pecotmr:::.spStudyDataForm(qd, "wrong"), "not in QtlDataset")
})

test_that(".sp* helpers: MultiStudy per-study context/trait routing", {
    qd1 <- .js_makeQtlDataset(
        studyName = "indA",
        contexts = "brain",
        traits = c("g1", "g2")
    )
    ss <- .js_makeQtlSumStats(
        studies = "ssC",
        contexts = "DLPFC",
        traits = "g3"
    )
    mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qd1), sumStats = ss)
    expect_setequal(pecotmr:::.spListContexts(mt, "ssC"), "DLPFC") # ss branch
    expect_setequal(pecotmr:::.spListContexts(mt, "indA"), "brain")
    expect_setequal(pecotmr:::.spListContexts(mt), c("brain", "DLPFC")) # all
    expect_equal(pecotmr:::.spListContexts(mt, "nope"), character(0))
    # The sumstats branch.
    expect_setequal(pecotmr:::.spListTraits(mt, studyName = "ssC"), "g3")
    expect_setequal(
        pecotmr:::.spListTraits(mt, studyName = "indA"),
        c("g1", "g2")
    )
    expect_equal(pecotmr:::.spListTraits(mt, studyName = "nope"), character(0))
    # study = NULL aggregates ALL traits across individual + sumstats components
    # (Regression guard: a present sumStats slot must not shadow the
    # QtlDatasets.)
    expect_setequal(pecotmr:::.spListTraits(mt), c("g1", "g2", "g3"))
})

test_that(".spListTraits: no study on a sumstats-free MultiStudy aggregates", {
    qdA <- .js_makeQtlDataset(
        studyName = "indA",
        contexts = "brain",
        traits = c("g1", "g2")
    )
    qdB <- .js_makeQtlDataset(
        studyName = "indB",
        contexts = "liver",
        traits = "g3"
    )
    mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qdA, indB = qdB))
    # Aggregated over the studies.
    expect_setequal(pecotmr:::.spListTraits(mt), c("g1", "g2", "g3"))
})

test_that(".sp* helpers: unsupported class raises a labelled error", {
    expect_error(pecotmr:::.spListStudies(42), "unsupported class")
    expect_error(pecotmr:::.spStudyDataForm(42, "x"), "unsupported class")
    expect_error(pecotmr:::.spListContexts(42), "unsupported class")
    expect_error(pecotmr:::.spListTraits(42), "unsupported class")
})

# -----------------------------------------------------------------------------
# parseJointSpecification / parseContexts / parseTraitIds error branches
# -----------------------------------------------------------------------------

test_that("parseJointSpecification: malformed inputs error", {
    qd <- .jsQd
    expect_error(pecotmr:::parseJointSpecification(42, qd), "must be NULL")
    expect_error(
        pecotmr:::parseJointSpecification(list(list(scope = NULL)), qd),
        "missing `axes`"
    )
    expect_error(
        pecotmr:::parseJointSpecification(list(42), qd),
        "character vector or a named list"
    )
    expect_error(
        pecotmr:::parseJointSpecification(list(list(axes = character(0))), qd),
        "non-empty character"
    )
    expect_error(
        pecotmr:::parseJointSpecification(
            list(list(axes = "context", scope = c("a", "b"))),
            qd
        ),
        "named list"
    )
    expect_error(
        pecotmr:::parseJointSpecification(
            list(list(axes = "context", scope = list(study = integer(0)))),
            qd
        ),
        "non-empty character vector"
    )
})

test_that("parseContexts: malformed inputs error", {
    qd <- .js_makeQtlDataset(studyName = "S1", contexts = c("brain", "liver"))
    expect_error(pecotmr:::parseContexts(character(0), qd), "non-empty")
    expect_error(
        pecotmr:::parseContexts(list(brain = "x"), qd),
        "unknown studies"
    )
    expect_error(
        pecotmr:::parseContexts(list(S1 = character(0)), qd),
        "non-empty character vector"
    )
    # Unnamed.
    expect_error(pecotmr:::parseContexts(list("brain"), qd), "named list")
    expect_error(pecotmr:::parseContexts(42, qd), "`contexts` must be NULL")
})

test_that("parseContexts: vector form refuses a context no study has", {
    # With one study a missing context can only be a mistake, and across
    # studies one that none of them has is the same: an error, not a warning.
    qd <- .js_makeQtlDataset(studyName = "S1", contexts = c("brain", "liver"))
    expect_error(
        pecotmr:::parseContexts(c("brain", "absent"), qd),
        "no study has: absent"
    )
})

test_that("parseTraitIds: malformed inputs error", {
    qd <- .js_makeQtlDataset(
        studyName = "S1",
        contexts = "brain",
        traits = c("ENSG1", "ENSG2")
    )
    expect_error(pecotmr:::parseTraitIds(character(0), qd), "non-empty")
    expect_error(pecotmr:::parseTraitIds(42, qd), "`traitId` must be NULL")
    expect_error(pecotmr:::parseTraitIds(list("ENSG1"), qd), "named by study")
    expect_error(
        pecotmr:::parseTraitIds(list(nope = "ENSG1"), qd),
        "unknown studies"
    )
    expect_error(
        pecotmr:::parseTraitIds(list(S1 = character(0)), qd),
        "non-empty character vector"
    )
    expect_error(
        pecotmr:::parseTraitIds(list(S1 = list("brain")), qd),
        "named by context"
    )
    expect_error(
        pecotmr:::parseTraitIds(list(S1 = list(brain = character(0))), qd),
        "non-empty character vector"
    )
    expect_error(
        pecotmr:::parseTraitIds(list(S1 = list(brain = "nope")), qd),
        "unknown traits"
    )
    expect_error(
        pecotmr:::parseTraitIds(list(S1 = 42), qd),
        "character vector or a named list"
    )
})
# (The doubly-nested study->context success path is covered by the existing
#  "parseTraitIds: doubly-nested study->context validates per context" test.)

# -----------------------------------------------------------------------------
# .fmResolveSpecScope: scope$study / scope$context / scope$trait + filters
# -----------------------------------------------------------------------------

test_that(".fmResolveSpecScope: scope, contexts and traitIds intersect", {
    qd <- .js_makeQtlDataset(
        studyName = "S1",
        contexts = c("brain", "liver"),
        traits = c("ENSG1", "ENSG2")
    )
    spec <- list(scope = list(study = "S1", context = "brain", trait = "ENSG1"))
    out <- pecotmr:::.fmResolveSpecScope(spec, qd)
    expect_equal(out$studies, "S1")
    expect_equal(out$contexts$S1, "brain")
    expect_equal(out$traits$S1, "ENSG1")
    # Named-list contexts + study-keyed list traitIds filters.
    out2 <- pecotmr:::.fmResolveSpecScope(
        list(scope = NULL),
        qd,
        contexts = list(S1 = "liver"),
        traitIds = list(S1 = "ENSG2")
    )
    expect_equal(out2$contexts$S1, "liver")
    expect_equal(out2$traits$S1, "ENSG2")
})

test_that(".fmResolveSpecScope: a per-context traitId gives the union", {
    qd <- .js_makeQtlDataset(
        studyName = "S1",
        contexts = c("brain", "liver"),
        traits = c("ENSG1", "ENSG2")
    )
    nested <- list(S1 = list(brain = "ENSG1", liver = "ENSG2"))
    out <- pecotmr:::.fmResolveSpecScope(
        list(scope = NULL),
        qd,
        traitIds = nested
    )
    # A joint fit spans the study's contexts: it takes every trait selected
    # in any of them, and the cell runner drops the pairs left out.
    expect_setequal(out$traits$S1, c("ENSG1", "ENSG2"))
    expect_identical(out$traitSelection, nested)
    plain <- pecotmr:::.fmResolveSpecScope(
        list(scope = NULL),
        qd,
        traitIds = list(S1 = "ENSG1")
    )
    expect_null(plain$traitSelection)
})

test_that("a parsed selection is read per study and context", {
    nested <- list(A = list(brain = "g1"), B = "g2")
    expect_equal(pecotmr:::.selectedTraits(nested, "A", "brain"), "g1")
    # Unnamed contexts and studies keep everything.
    expect_null(pecotmr:::.selectedTraits(nested, "A", "liver"))
    expect_null(pecotmr:::.selectedTraits(nested, "C", "brain"))
    expect_equal(pecotmr:::.selectedTraits(nested, "B", "liver"), "g2")
    expect_equal(pecotmr:::.selectedTraits("g3", "A", "brain"), "g3")
    expect_true(pecotmr:::.traitSelectedAt("g9", NULL, "A", "brain"))
    expect_false(pecotmr:::.traitSelectedAt("g2", nested, "A", "brain"))
    # A study's traits across its contexts: the union, unless a context in
    # play is not named, which keeps everything.
    both <- list(A = list(brain = "g1", liver = c("g2", "g1")))
    expect_setequal(
        pecotmr:::.selectedStudyTraits(both, "A", c("brain", "liver")),
        c("g1", "g2")
    )
    expect_null(pecotmr:::.selectedStudyTraits(
        nested,
        "A",
        c("brain", "liver")
    ))
    expect_equal(pecotmr:::.selectedStudyTraits(nested, "A", "brain"), "g1")
    expect_true(pecotmr:::.hasNestedTraitSelection(nested))
    expect_false(pecotmr:::.hasNestedTraitSelection(list(A = "g1")))
    expect_false(pecotmr:::.hasNestedTraitSelection("g1"))
})

test_that("each study's share of traitId parses back to the same selection", {
    nested <- list(A = list(brain = "g1"), B = "g2")
    expect_equal(
        pecotmr:::.studyTraitIdArg(nested, "A"),
        list(A = list(brain = "g1"))
    )
    expect_equal(pecotmr:::.studyTraitIdArg(nested, "B"), "g2")
    expect_null(pecotmr:::.studyTraitIdArg(nested, "C"))
    expect_equal(pecotmr:::.studyTraitIdArg("g3", "A"), "g3")
})

test_that(".sumStatsRowsSelected applies study and per-context selections", {
    rows <- data.frame(
        study = c("A", "A", "B", "B"),
        context = c("brain", "liver", "brain", "brain"),
        trait = c("g1", "g1", "g1", "g2")
    )
    expect_equal(
        pecotmr:::.sumStatsRowsSelected(
            rows,
            list(A = "liver", B = "brain"),
            NULL
        ),
        c(FALSE, TRUE, TRUE, TRUE)
    )
    expect_equal(
        pecotmr:::.sumStatsRowsSelected(rows, NULL, list(B = "g2")),
        c(TRUE, TRUE, FALSE, TRUE)
    )
    expect_equal(
        pecotmr:::.sumStatsRowsSelected(
            rows,
            NULL,
            list(A = list(liver = "g2"))
        ),
        c(TRUE, FALSE, TRUE, TRUE)
    )
})

test_that("a study's selected traits follow traitId, else region", {
    # .js_makeSe places ENSG1 at chr1:100-149 and ENSG2 at chr1:200-249.
    qd <- .js_makeQtlDataset(studyName = "A")
    sel <- function(cx, traitId = NULL, region = NULL) {
        pecotmr:::.contextSelectedTraits(cx, qd, traitId, region)
    }
    expect_equal(sel("brain"), c("ENSG1", "ENSG2"))
    expect_equal(sel("brain", traitId = c("ENSG2", "ENSG9")), "ENSG2")
    # A selection that leaves the context out keeps all of its traits.
    expect_equal(
        sel("liver", traitId = list(A = list(brain = "ENSG1"))),
        c("ENSG1", "ENSG2")
    )
    # region picks traits only when no trait was named at all.
    nearEnsg1 <- GenomicRanges::GRanges("chr1:120-130")
    expect_equal(sel("brain", region = nearEnsg1), "ENSG1")
    expect_equal(sel("brain", traitId = "ENSG2", region = nearEnsg1), "ENSG2")
    has <- function(contexts = NULL, traitId = NULL, region = NULL) {
        pecotmr:::.studyHasSelectedTraits(qd, contexts, traitId, region)
    }
    expect_true(has(traitId = "ENSG1"))
    expect_false(has(traitId = "ENSG9"))
    expect_false(has(region = GenomicRanges::GRanges("chr1:5000-6000")))
    # None of the study's contexts selected leaves it no trait either.
    expect_false(has(contexts = list(A = character(0))))
    expect_true(has(contexts = list(A = "liver"), traitId = "ENSG2"))
})

test_that("a multi-study run is refused only when no component has a trait", {
    msqd <- MultiStudyQtlDataset(
        qtlDatasets = list(A = .js_makeQtlDataset(studyName = "A")),
        sumStats = .js_makeQtlSumStats(studies = "ssA")
    )
    check <- function(traitId) {
        pecotmr:::.msAssertTraitsSelected(
            msqd,
            list(contexts = NULL, traitId = traitId, region = NULL),
            "aPipeline"
        )
    }
    # A trait the individual-level study carries, or one only the summary
    # statistics carry: either is enough.
    expect_silent(check("ENSG1"))
    expect_silent(check("ENSG3"))
    expect_error(
        check("ENSG9"),
        "aPipeline(MultiStudyQtlDataset): no traits selected.",
        fixed = TRUE
    )
})

# -----------------------------------------------------------------------------
# .buildJointSumstatZMatrix: SNP-order mismatch error
# -----------------------------------------------------------------------------

test_that(".buildJointSumstatZMatrix: entries in different SNP orders error", {
    df <- data.frame(
        study = "S",
        context = c("c1", "c2"),
        trait = "t1",
        stringsAsFactors = FALSE
    )
    calls <- 0L
    local_mocked_bindings(
        .sumStatsDf = function(x, studyName, context, trait, require, ...) {
            calls <<- calls + 1L
            vid <- if (calls == 1L) {
                c("chr1:100:A:G", "chr1:200:A:G")
            } else {
                c("chr1:200:A:G", "chr1:100:A:G")
            } # reordered
            data.frame(
                variant_id = vid,
                z = c(1, 2),
                N = c(100, 100),
                stringsAsFactors = FALSE
            )
        },
        .package = "pecotmr"
    )
    expect_error(
        pecotmr:::.buildJointSumstatZMatrix(
            df,
            c(1L, 2L),
            c("c1", "c2"),
            errorLabel = "TESTLABEL"
        ),
        "identical SNP order"
    )
})

# -----------------------------------------------------------------------------
# Individual X/Y builders: the skip paths (return NULL / message branches)
# -----------------------------------------------------------------------------

test_that("a cross-context (X, Y): skips a trait in < 2 contexts", {
    se1 <- .js_makeSe(traits = "g1") # g1 present
    se0 <- .js_makeSe(traits = "other") # g1 absent
    local_mocked_bindings(
        .qtlContextSe = function(data, context) {
            if (identical(contexts, "c1")) se1 else se0
        },
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_crossContextXy(
        NULL,
        "g1",
        c("c1", "c2"),
        cisWindow = 1000L,
        verbose = 1,
        label = "X"
    )))
})

test_that("a cross-context (X, Y): region path + complete-case skip", {
    se <- .js_makeSe(traits = "g1", samples = paste0("s", 1:6))
    samp <- paste0("s", 1:6)
    region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10000))
    local_mocked_bindings(
        .qtlContextSe = function(data, context) se,
        .fmResidGeno = function(
            x,
            contexts,
            traitId = NULL,
            cisWindow = NULL,
            region = NULL,
            ...
        ) {
            matrix(
                0,
                6,
                2,
                dimnames = list(samp, c("chr1:100:A:G", "chr1:200:A:G"))
            )
        },
        .fmResidPheno = function(x, contexts, traitId = NULL, ...) {
            ym <- function(v) matrix(v, 6, 1, dimnames = list(samp, "g1"))
            list(
                c1 = ym(c(NA, NA, NA, NA, NA, 1)), # mostly NA -> < 2 complete
                c2 = ym(rnorm(6))
            )
        },
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_crossContextXy(
        NULL,
        "g1",
        c("c1", "c2"),
        cisWindow = NULL,
        verbose = 1,
        label = "X",
        region = region
    )))
})

test_that("a cross-context (X, Y): too few shared samples skips", {
    se <- .js_makeSe(traits = "g1", samples = paste0("s", 1:6))
    local_mocked_bindings(
        .qtlContextSe = function(data, context) se,
        .fmResidGeno = function(
            x,
            contexts,
            traitId = NULL,
            cisWindow = NULL,
            region = NULL,
            ...
        ) {
            matrix(
                0,
                1,
                2,
                dimnames = list("zz", c("chr1:100:A:G", "chr1:200:A:G"))
            )
        }, # disjoint sample
        .fmResidPheno = function(x, contexts, traitId = NULL, ...) {
            list(
                c1 = matrix(0, 6, 1, dimnames = list(paste0("s", 1:6), "g1")),
                c2 = matrix(0, 6, 1, dimnames = list(paste0("s", 1:6), "g1"))
            )
        },
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_crossContextXy(
        NULL,
        "g1",
        c("c1", "c2"),
        cisWindow = 1000L,
        verbose = 1,
        label = "X"
    )))
})

test_that(".crossContextFitSamples skips when too few samples are shared", {
    Yres <- list(
        c1 = matrix(0, 2, 1, dimnames = list(c("s1", "s2"), "g1")),
        c2 = matrix(0, 2, 1, dimnames = list(c("s2", "s3"), "g1"))
    )
    expect_message(
        res <- pecotmr:::.crossContextFitSamples(
            Yres,
            c("c1", "c2"),
            list(verbose = 1, label = "lbl", tid = "g1")
        ),
        "trait 'g1' has too few shared samples across contexts"
    )
    expect_null(res)
})

test_that(".msStudySpecs narrows a joint spec's scope to one component", {
    skip_if_not_installed("pgenlibr")
    qd <- .fitSamplesDataset()
    spec <- function(scope) list(axes = "context", scope = scope)
    # No scope, or one the component holds all of, stands as it is.
    expect_identical(
        pecotmr:::.msStudySpecs(list(spec(NULL)), qd),
        list(spec(NULL))
    )
    # What it holds is kept, as re-parsing refuses what the data lacks.
    expect_identical(
        pecotmr:::.msStudySpecs(
            list(spec(list(context = c("brain", "elsewhere")))),
            qd
        ),
        list(spec(list(context = "brain")))
    )
    # A scope naming only what other components hold leaves this one no spec.
    expect_length(
        pecotmr:::.msStudySpecs(list(spec(list(study = "elsewhere"))), qd),
        0L
    )
    expect_length(pecotmr:::.msStudySpecs(NULL, qd), 0L)
})

test_that(".fmTraitsInRegion keeps overlapping traits; NULL keeps all", {
    se <- .js_makeSe(traits = c("g1", "g2")) # g1@~100, g2@~200
    region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(90, 160))
    expect_equal(pecotmr:::.fmTraitsInRegion(se, c("g1", "g2"), region), "g1")
    both <- GenomicRanges::GRanges("chr1", IRanges::IRanges(90, 260))
    expect_setequal(
        pecotmr:::.fmTraitsInRegion(se, c("g1", "g2"), both),
        c("g1", "g2")
    )
    expect_equal(
        pecotmr:::.fmTraitsInRegion(se, c("g1", "g2"), NULL),
        c("g1", "g2")
    ) # NULL region -> unchanged
})

test_that("a cross-trait (X, Y): < 2 traits, region, complete skip", {
    se2 <- .js_makeSe(traits = c("g1", "g2"), samples = paste0("s", 1:6))
    samp <- paste0("s", 1:6)
    # < 2 scoped traits in the context -> NULL.
    local_mocked_bindings(
        .qtlContextSe = function(data, context) se2,
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_crossTraitXy(
        NULL,
        "cx",
        "g1",
        cisWindow = 1000L,
        verbose = 1,
        label = "X",
        studyName = "S"
    )))
    # region path + < 2 complete cases.
    local_mocked_bindings(
        .qtlContextSe = function(data, context) se2,
        .fmResidGeno = function(
            x,
            contexts,
            traitId = NULL,
            cisWindow = NULL,
            region = NULL,
            ...
        ) {
            matrix(
                0,
                6,
                2,
                dimnames = list(samp, c("chr1:100:A:G", "chr1:200:A:G"))
            )
        },
        .fmResidPheno = function(x, contexts, traitId = NULL, ...) {
            cbind(g1 = c(NA, NA, NA, NA, NA, 1), g2 = rnorm(6)) |>
                `rownames<-`(samp)
        },
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_crossTraitXy(
        NULL,
        "cx",
        c("g1", "g2"),
        cisWindow = NULL,
        verbose = 1,
        label = "X",
        studyName = "S",
        region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 9999))
    )))
})

test_that("a composed (X, Y): skip branches and single-context wrap", {
    se <- .js_makeSe(traits = c("g1", "g2"), samples = paste0("s", 1:6))
    samp <- paste0("s", 1:6)
    scope <- list(contexts = list(S = "c1"), traits = list(S = c("g1", "g2")))
    # Single context: YresList wrap branch, then a valid 2-tuple build.
    local_mocked_bindings(
        .qtlContextSe = function(data, context) se,
        .fmResidGeno = function(
            x,
            contexts,
            traitId = NULL,
            cisWindow = NULL,
            region = NULL,
            ...
        ) {
            matrix(
                0,
                6,
                2,
                dimnames = list(samp, c("chr1:100:A:G", "chr1:200:A:G"))
            )
        },
        .fmResidPheno = function(x, contexts, traitId = NULL, ...) {
            matrix(rnorm(12), 6, 2, dimnames = list(samp, c("g1", "g2")))
        },
        .package = "pecotmr"
    )
    out <- suppressMessages(.js_composedXy(
        NULL,
        scope,
        "S",
        cisWindow = 1000L,
        verbose = 1,
        label = "X"
    ))
    expect_equal(ncol(out$Y), 2L)
    expect_setequal(colnames(out$Y), c("c1:g1", "c1:g2"))
    # < 2 tuples -> NULL.
    scope1 <- list(contexts = list(S = "c1"), traits = list(S = "g1"))
    local_mocked_bindings(
        .qtlContextSe = function(data, context) {
            .js_makeSe(traits = "g1")
        },
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_composedXy(
        NULL,
        scope1,
        "S",
        cisWindow = 1000L,
        verbose = 1,
        label = "X"
    )))
})

# -----------------------------------------------------------------------------
# .enumerateComposedSumstatGroups: empty scope + all-axes (no complement)
# -----------------------------------------------------------------------------

test_that(".enumerateComposedSumstatGroups: empty scope and no complement", {
    df <- data.frame(
        study = "S",
        context = c("c1", "c2"),
        trait = "t1",
        stringsAsFactors = FALSE
    )
    empty <- pecotmr:::.enumerateComposedSumstatGroups(
        list(axes = c("context", "trait")),
        df,
        list(studies = character(0), contexts = list(), traits = list())
    )
    expect_null(empty)
    # axes = all three -> complement empty -> a single "__all__" block.
    gi <- pecotmr:::.enumerateComposedSumstatGroups(
        list(axes = c("study", "context", "trait")),
        df,
        list(
            studies = "S",
            contexts = list(S = c("c1", "c2")),
            traits = list(S = "t1")
        )
    )
    expect_equal(length(gi$groups), 1L)
    expect_equal(names(gi$groups), "__all__")
})

# -----------------------------------------------------------------------------
# .defaultJointSpec
# -----------------------------------------------------------------------------

test_that(".defaultJointSpec: trait beats context; single axes -> empty", {
    expect_equal(pecotmr:::.defaultJointSpec(3L, 2L)[[1L]]$axes, "trait")
    expect_equal(pecotmr:::.defaultJointSpec(2L, 1L)[[1L]]$axes, "context")
    expect_equal(pecotmr:::.defaultJointSpec(1L, 1L), list())
})

# -----------------------------------------------------------------------------
# Multi-region merges: .fmMergeResultsByKey / .twasMergeResultsByKey
# -----------------------------------------------------------------------------

.js_fmEntry <- function(vid = "chr1:100:A:G") {
    FineMappingRow(
        variantIds = vid,
        susieFit = list(),
        topLoci = data.frame(
            variant_id = vid,
            pip = 0.9,
            stringsAsFactors = FALSE
        )
    )
}

test_that(".fmMergeResultsByKey: merges per-region entries by (s,c,t,method)", {
    mk <- function() {
        QtlFineMappingResult(
            studyName = "S",
            context = "c1",
            trait = "t1",
            method = "mvsusie",
            entry = list(.js_fmEntry())
        )
    }
    local_mocked_bindings(
        .fmMergeEntries = function(entries) entries[[1L]],
        .package = "pecotmr"
    )
    out <- pecotmr:::.fmMergeResultsByKey(list(mk(), mk()))
    expect_s4_class(out, "QtlFineMappingResult")
    expect_equal(nrow(out), 1L)
    # n == 0 short-circuit returns the (empty) base unchanged.
    empty <- QtlFineMappingResult(
        studyName = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        entry = list()
    )
    expect_equal(nrow(pecotmr:::.fmMergeResultsByKey(list(empty, empty))), 0L)
})

test_that(".fmMergeResultsByKey keeps joint groups apart, and their places", {
    # One tuple in two joint fits per region; the regions agree on the first
    # fit's position for it and disagree on the second's.
    region <- function(jointIndex) {
        fm <- QtlFineMappingResult(
            studyName = c("S", "S"),
            context = c("c1", "c1"),
            trait = c("t1", "t1"),
            method = c("mvsusie", "mvsusie"),
            entry = list(.js_fmEntry(), .js_fmEntry()),
            jointStudies = c(NA, "S;T"),
            jointContexts = c("c1;c2", NA)
        )
        fm$jointIndex <- jointIndex
        fm
    }
    seen <- list()
    local_mocked_bindings(
        .fmMergeEntries = function(entries) {
            seen[[length(seen) + 1L]] <<- length(entries)
            entries[[1L]]
        },
        .package = "pecotmr"
    )
    out <- pecotmr:::.fmMergeResultsByKey(list(region(1:2), region(c(1L, 1L))))
    expect_equal(nrow(out), 2L)
    expect_equal(as.character(out$jointStudies), c(NA, "S;T"))
    # Each row merged its own group's two regions.
    expect_equal(unlist(seen), c(2L, 2L))
    expect_equal(out$jointIndex, c(1L, NA))
    # No region records one: no column.
    plain <- function() {
        QtlFineMappingResult(
            studyName = "S",
            context = "c1",
            trait = "t1",
            method = "mvsusie",
            entry = list(.js_fmEntry()),
            jointContexts = "c1;c2"
        )
    }
    merged <- pecotmr:::.fmMergeResultsByKey(list(plain(), plain()))
    expect_false("jointIndex" %in% colnames(S4Vectors::mcols(merged)))
})

test_that(".twasMergeResultsByKey: merges per-region TwasWeights entries", {
    mk <- function() {
        TwasWeights(
            studyName = "S",
            context = "c1",
            trait = "t1",
            method = "lasso",
            entry = list(TwasWeightsRow(
                variantIds = "chr1:100:A:G",
                weights = 0.5
            ))
        )
    }
    out <- pecotmr:::.twasMergeResultsByKey(list(mk(), mk()), c("r1", "r2"))
    expect_s4_class(out, "TwasWeights")
    expect_equal(nrow(out), 1L)
    empty <- TwasWeights(
        studyName = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        entry = list()
    )
    expect_equal(
        length(
            pecotmr:::.twasMergeResultsByKey(
                list(empty),
                "r1"
            )$method
        ),
        0L
    )
})

# -----------------------------------------------------------------------------
# Multi-region QtlDataset dispatch (xRegions length 2 -> merge path)
# -----------------------------------------------------------------------------

test_that(".twasDispatchJointSpecsQtlDataset: two blocks merge by key", {
    r1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(50, 250))
    r2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(300, 500))
    parsed <- list(list(axes = "context", scope = NULL))
    mkRegionRes <- function() {
        TwasWeights(
            studyName = c("Q1", "Q1"),
            context = c("c1", "c2"),
            trait = c("t1", "t1"),
            method = c("mrmash", "mrmash"),
            entry = list(
                TwasWeightsRow(variantIds = "chr1:100:A:G", weights = 0.5),
                TwasWeightsRow(variantIds = "chr1:100:A:G", weights = 0.5)
            )
        )
    }
    local_mocked_bindings(
        .twasJointSpecRegionWork = function(...) {
            pecotmr:::.qtlWorkValue(mkRegionRes())
        },
        .package = "pecotmr"
    )
    res <- pecotmr:::.twasDispatchJointSpecsQtlDataset(
        parsed,
        "mrmash",
        marker = NULL,
        ctx = list(xRegions = list(r1, r2)),
        contexts = NULL,
        traitIds = NULL
    )
    expect_s4_class(res, "TwasWeights")
    expect_equal(nrow(res), 2L) # per-context, merged regions
})

# -----------------------------------------------------------------------------
# MultiStudy dispatchers: .fmDispatchJointSpecsMultiStudy /
# .twasDispatchJointSpecsMultiStudy (leaf dispatchers mocked)
# -----------------------------------------------------------------------------

test_that(".fmDispatchJointSpecsMultiStudy: study specs go to sumstats", {
    qd <- .js_makeQtlDataset(
        studyName = "indA",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    ss <- .js_makeQtlSumStats(studies = "ssC", contexts = "c1", traits = "t1")
    mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qd), sumStats = ss)
    parsed <- list(
        list(axes = "context", scope = NULL),
        list(axes = "study", scope = NULL)
    )
    qdRes <- QtlFineMappingResult(
        studyName = "indA",
        context = "c1",
        trait = "t1",
        method = "mvsusie",
        entry = list(.js_fmEntry())
    )
    ssRes <- QtlFineMappingResult(
        studyName = "ssC",
        context = "c1",
        trait = "t1",
        method = "mvsusie",
        entry = list(.js_fmEntry()),
        ldReference = .js_makeGenotypeHandle()
    )
    local_mocked_bindings(
        .fmDispatchJointSpecsQtlDataset = function(...) qdRes,
        .fmDispatchJointSpecsQtlSumStats = function(...) ssRes,
        .package = "pecotmr"
    )
    out <- suppressMessages(pecotmr:::.fmDispatchJointSpecsMultiStudy(
        parsed,
        "mvsusie",
        list(),
        list(
            data = mt,
            verbose = 1,
            xRegions = list(NULL),
            genotypeFilterParam = GenotypeFilterParam(),
            credibleSetParam = CredibleSetParam(
                coverage = 0.95,
                secondaryCoverage = 0.5,
                signalCutoff = 0.1,
                minAbsCorr = 0.5
            )
        )
    ))
    expect_s4_class(out, "QtlFineMappingResult")
    expect_equal(nrow(out), 2L) # indA row + ssC row
})

test_that(".fmDispatchJointSpecsMultiStudy: a study spec joins the studies", {
    qdA <- .js_makeQtlDataset(
        studyName = "indA",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    qdB <- .js_makeQtlDataset(
        studyName = "indB",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    # No sumstats: the individual-level studies are the cross-study fit.
    mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qdA, indB = qdB))
    acrossStudies <- QtlFineMappingResult(
        studyName = c("indA", "indB"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("mvsusie", "mvsusie"),
        entry = list(.js_fmEntry(), .js_fmEntry())
    )
    seen <- new.env(parent = emptyenv())
    local_mocked_bindings(
        .fmDispatchStudySpecsIndividual = function(studySpecs, ...) {
            seen$axes <- map(studySpecs, "axes")
            acrossStudies
        },
        .package = "pecotmr"
    )
    parsed <- list(
        list(axes = "study", scope = NULL),
        list(axes = c("study", "context"), scope = NULL)
    )
    expect_message(
        out <- pecotmr:::.fmDispatchJointSpecsMultiStudy(
            parsed,
            "mvsusie",
            list(),
            list(data = mt, verbose = 1, xRegions = list(NULL))
        ),
        "composed specs including the study axis produce no result"
    )
    # Only the study-axis-only spec reaches the individual-level studies.
    expect_equal(seen$axes, list("study"))
    expect_identical(out, acrossStudies)
})

test_that(".fmDispatchStudySpecsIndividual: one study, or no mvsusie: NULL", {
    qdA <- .js_makeQtlDataset(studyName = "indA", traits = "t1")
    ss <- .js_makeQtlSumStats(studies = "ssC", contexts = "c1", traits = "t1")
    oneStudy <- MultiStudyQtlDataset(
        qtlDatasets = list(indA = qdA),
        sumStats = ss
    )
    spec <- list(list(axes = "study", scope = NULL))
    expect_null(pecotmr:::.fmDispatchStudySpecsIndividual(
        spec,
        "mvsusie",
        list(),
        list(data = oneStudy, verbose = 0)
    ))
    qdB <- .js_makeQtlDataset(studyName = "indB", traits = "t1")
    twoStudies <- MultiStudyQtlDataset(
        qtlDatasets = list(indA = qdA, indB = qdB)
    )
    expect_message(
        out <- pecotmr:::.fmDispatchStudySpecsIndividual(
            spec,
            "fsusie",
            list(),
            list(data = twoStudies, verbose = 1)
        ),
        "fsusie has no cross-study fit"
    )
    expect_null(out)
})

test_that(".twasDispatchJointSpecsMultiStudy: routes and rbinds by axis", {
    qd <- .jd_makeQtlDataset(
        studyName = "indA",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    ss <- .jd_makeQtlSumStats(studies = "ssC", contexts = "c1", traits = "t1")
    mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qd), sumStats = ss)
    parsed <- list(
        list(axes = "context", scope = NULL),
        list(axes = "study", scope = NULL)
    )
    qdRes <- TwasWeights(
        studyName = "indA",
        context = "c1",
        trait = "t1",
        method = "mrmash",
        entry = list(TwasWeightsRow(
            variantIds = "chr1:100:A:G",
            weights = 0.5
        ))
    )
    ssRes <- TwasWeights(
        studyName = "ssC",
        context = "c1",
        trait = "t1",
        method = "mrmash",
        entry = list(TwasWeightsRow(
            variantIds = "chr1:100:A:G",
            weights = 0.5
        )),
        ldReference = .jd_makeHandle()
    )
    local_mocked_bindings(
        .twasStudyJointResult = function(...) qdRes,
        .twasDispatchJointSpecsQtlSumStats = function(...) ssRes,
        .package = "pecotmr"
    )
    out <- suppressMessages(pecotmr:::.twasDispatchJointSpecsMultiStudy(
        parsed,
        "mrmash",
        "mrmash",
        list(data = mt, verbose = 1)
    ))
    expect_s4_class(out, "TwasWeights")
    expect_equal(nrow(out), 2L)
})

test_that(".twasDispatchJointSpecsMultiStudy: study spec, no sumStats: NULL", {
    qdA <- .jd_makeQtlDataset(
        studyName = "indA",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    qdB <- .jd_makeQtlDataset(
        studyName = "indB",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    mt <- MultiStudyQtlDataset(qtlDatasets = list(indA = qdA, indB = qdB))
    parsed <- list(list(axes = "study", scope = NULL))
    expect_message(
        out <- pecotmr:::.twasDispatchJointSpecsMultiStudy(
            parsed,
            "mrmash",
            "mrmash",
            list(data = mt, verbose = 1)
        ),
        "no sumStats slot"
    )
    expect_null(out)
})

# ============================================================================
# Mop-up: remaining .sp* empties, X/Y builder skip paths, region-missing keys
# in merges, FM multi-region merge.
# ============================================================================

test_that(".sp* QtlDataset: mismatched study / absent context return empty", {
    qd <- .js_makeQtlDataset(
        studyName = "S1",
        contexts = "brain",
        traits = c("g1", "g2")
    )
    expect_equal(pecotmr:::.spListContexts(qd, "wrong"), character(0))
    expect_equal(pecotmr:::.spListTraits(qd, studyName = "wrong"), character(0))
    expect_equal(pecotmr:::.spListTraits(qd, context = "nope"), character(0))
})

test_that("a cross-trait (X, Y): disjoint X/Y samples skip", {
    se <- .js_makeSe(traits = c("g1", "g2"), samples = paste0("s", 1:6))
    local_mocked_bindings(
        .qtlContextSe = function(data, context) se,
        .fmResidGeno = function(
            x,
            contexts,
            traitId = NULL,
            cisWindow = NULL,
            region = NULL,
            ...
        ) {
            matrix(
                0,
                1,
                2,
                dimnames = list("zz", c("chr1:100:A:G", "chr1:200:A:G"))
            )
        }, # disjoint sample
        .fmResidPheno = function(x, contexts, traitId = NULL, ...) {
            matrix(0, 6, 2, dimnames = list(paste0("s", 1:6), c("g1", "g2")))
        },
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_crossTraitXy(
        NULL,
        "cx",
        c("g1", "g2"),
        cisWindow = 1000L,
        verbose = 1,
        label = "X",
        studyName = "S"
    )))
})

test_that("a composed (X, Y): disjoint, missing or NA rows skip", {
    samp <- paste0("s", 1:6)
    se <- .js_makeSe(traits = c("g1", "g2"), samples = samp)
    scope <- list(
        contexts = list(S = c("c1", "c2")),
        traits = list(S = c("g1", "g2"))
    )
    # (a) disjoint X samples -> < 2 common -> NULL.
    local_mocked_bindings(
        .qtlContextSe = function(data, context) se,
        .fmResidGeno = function(
            x,
            contexts,
            traitId = NULL,
            cisWindow = NULL,
            region = NULL,
            ...
        ) {
            matrix(
                0,
                1,
                2,
                dimnames = list("zz", c("chr1:100:A:G", "chr1:200:A:G"))
            )
        },
        .fmResidPheno = function(x, contexts, traitId = NULL, ...) {
            setNames(
                map(c("c1", "c2"), function(.) {
                    matrix(
                        rnorm(12),
                        6,
                        2,
                        dimnames = list(samp, c("g1", "g2"))
                    )
                }),
                c("c1", "c2")
            )
        },
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_composedXy(
        NULL,
        scope,
        "S",
        cisWindow = 1000L,
        verbose = 1,
        label = "X"
    )))
    # (b) one context's Y lacks the trait column -> tuple skipped -> < 2 yCols.
    local_mocked_bindings(
        .qtlContextSe = function(data, context) se,
        .fmResidGeno = function(
            x,
            contexts,
            traitId = NULL,
            cisWindow = NULL,
            region = NULL,
            ...
        ) {
            matrix(
                0,
                6,
                2,
                dimnames = list(samp, c("chr1:100:A:G", "chr1:200:A:G"))
            )
        },
        .fmResidPheno = function(x, contexts, traitId = NULL, ...) {
            list(
                c1 = matrix(0, 6, 1, dimnames = list(samp, "g1")),
                c2 = matrix(0, 6, 1, dimnames = list(samp, "zzz"))
            )
        }, # no g1/g2
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_composedXy(
        NULL,
        list(contexts = list(S = c("c1", "c2")), traits = list(S = "g1")),
        "S",
        cisWindow = 1000L,
        verbose = 1,
        label = "X"
    )))
    # (c) two valid columns but < 2 complete rows (NA) -> NULL.
    local_mocked_bindings(
        .qtlContextSe = function(data, context) se,
        .fmResidGeno = function(
            x,
            contexts,
            traitId = NULL,
            cisWindow = NULL,
            region = NULL,
            ...
        ) {
            matrix(
                0,
                6,
                2,
                dimnames = list(samp, c("chr1:100:A:G", "chr1:200:A:G"))
            )
        },
        .fmResidPheno = function(x, contexts, traitId = NULL, ...) {
            list(
                c1 = matrix(
                    c(NA, NA, NA, NA, NA, 1),
                    6,
                    1,
                    dimnames = list(samp, "g1")
                ),
                c2 = matrix(
                    c(NA, NA, NA, NA, NA, 1),
                    6,
                    1,
                    dimnames = list(samp, "g1")
                )
            )
        },
        .package = "pecotmr"
    )
    expect_null(suppressMessages(.js_composedXy(
        NULL,
        list(contexts = list(S = c("c1", "c2")), traits = list(S = "g1")),
        "S",
        cisWindow = 1000L,
        verbose = 1,
        label = "X"
    )))
})

test_that(".fmMergeResultsByKey: a key a later region lacks adds nothing", {
    twoRow <- QtlFineMappingResult(
        studyName = c("S", "S"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("mvsusie", "mvsusie"),
        entry = list(.js_fmEntry("chr1:100:A:G"), .js_fmEntry("chr1:200:A:G"))
    )
    oneRow <- QtlFineMappingResult(
        studyName = "S",
        context = "c1",
        trait = "t1",
        method = "mvsusie",
        entry = list(.js_fmEntry("chr1:100:A:G"))
    ) # missing the c2 key
    seen <- 0L
    local_mocked_bindings(
        .fmMergeEntries = function(entries) {
            seen <<- seen + length(entries)
            entries[[1L]]
        },
        .package = "pecotmr"
    )
    out <- pecotmr:::.fmMergeResultsByKey(list(twoRow, oneRow))
    expect_equal(nrow(out), 2L)
    expect_equal(seen, 3L) # 2 for c1 row + 1 for c2 row
})

test_that(".fmDispatchJointSpecsQtlDataset: two region blocks are merged", {
    qd <- .jd_makeQtlDataset(
        studyName = "Q1",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    r1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(50, 250))
    r2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(300, 500))
    mkRes <- function() {
        QtlFineMappingResult(
            studyName = c("Q1", "Q1"),
            context = c("c1", "c2"),
            trait = c("t1", "t1"),
            method = c("mvsusie", "mvsusie"),
            entry = list(
                .js_fmEntry("chr1:100:A:G"),
                .js_fmEntry("chr1:100:A:G")
            )
        )
    }
    local_mocked_bindings(
        .fmJointSpecRegionWork = function(...) {
            pecotmr:::.qtlWorkValue(mkRes())
        },
        .fmMergeEntries = function(entries) entries[[1L]],
        .package = "pecotmr"
    )
    res <- pecotmr:::.fmDispatchJointSpecsQtlDataset(
        list(list(axes = "context", scope = NULL)),
        "mvsusie",
        list(),
        list(data = qd, verbose = 0, xRegions = list(r1, r2))
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(nrow(res), 2L) # the two region blocks merged by key
})

test_that(".fmDispatchJointSpecsQtlDataset: one region direct; all-NULL NULL", {
    qd <- .jd_makeQtlDataset(
        studyName = "Q1",
        contexts = c("c1", "c2"),
        traits = "t1"
    )
    res1 <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "mvsusie",
        entry = list(.js_fmEntry("chr1:100:A:G"))
    )
    local_mocked_bindings(
        .fmJointSpecRegionWork = function(...) pecotmr:::.qtlWorkValue(res1),
        .package = "pecotmr"
    )
    out <- pecotmr:::.fmDispatchJointSpecsQtlDataset(
        list(list(axes = "context", scope = NULL)),
        "mvsusie",
        list(),
        list(data = qd, verbose = 0, cisWindow = 1000L, xRegions = list(NULL))
    )
    expect_identical(out, res1) # length-1 short-circuit
    local_mocked_bindings(
        .fmJointSpecRegionWork = function(...) pecotmr:::.qtlWorkValue(NULL),
        .package = "pecotmr"
    )
    expect_null(pecotmr:::.fmDispatchJointSpecsQtlDataset(
        list(list(axes = "context", scope = NULL)),
        "mvsusie",
        list(),
        list(data = qd, verbose = 0, cisWindow = 1000L, xRegions = list(NULL))
    )) # all-NULL
})

test_that(".twasMergeResultsByKey keeps an entry whichever regions hold it", {
    twoRow <- TwasWeights(
        studyName = c("S", "S"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("lasso", "lasso"),
        entry = list(
            TwasWeightsRow(variantIds = "chr1:100:A:G", weights = 0.1),
            TwasWeightsRow(variantIds = "chr1:200:A:G", weights = 0.2)
        )
    )
    oneRow <- TwasWeights(
        studyName = "S",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(TwasWeightsRow(
            variantIds = "chr1:300:A:G",
            weights = 0.3
        ))
    )
    out <- pecotmr:::.twasMergeResultsByKey(list(twoRow, oneRow), c("rA", "rB"))
    expect_equal(nrow(out), 2L) # c2 kept from the first region
    # Keys came from the first region alone, so one only a later region
    # held was lost.
    out <- pecotmr:::.twasMergeResultsByKey(list(oneRow, twoRow), c("rA", "rB"))
    expect_setequal(as.character(out$context), c("c1", "c2"))
})

test_that(".twasMergeResultsByKey concatenates one entry across regions", {
    one <- function(pos, w) {
        TwasWeights(
            studyName = "S",
            context = "c1",
            trait = "t1",
            method = "lasso",
            entry = list(TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", pos),
                weights = w
            )),
            region = GenomicRanges::GRanges(
                "chr1",
                IRanges::IRanges(min(pos), max(pos))
            )
        )
    }
    out <- pecotmr:::.twasMergeResultsByKey(
        list(one(c(100L, 200L), c(0.1, 0.2)), one(c(300L, 400L), c(0.3, 0.4))),
        c("rA", "rB")
    )
    expect_equal(nrow(out), 1L)
    m <- pecotmr:::.twrCellParts(out, 1L, "lasso")
    expect_equal(names(methodFits(m)), c("rA", "rB"))
    expect_equal(unname(weights(m)), c(0.1, 0.2, 0.3, 0.4))
    # The merged window spans both regions' windows; it used to be dropped.
    expect_equal(start(genomicRegion(out)), 100L)
    expect_equal(end(genomicRegion(out)), 400L)
})

test_that(".twasMergeResultsByKey keeps two joint groups of one gene apart", {
    group <- function(jointContexts, w) {
        TwasWeights(
            studyName = "S",
            context = "c1",
            trait = "t1",
            method = "mrmash",
            entry = list(TwasWeightsRow(
                variantIds = "chr1:100:A:G",
                weights = w
            )),
            jointContexts = jointContexts
        )
    }
    perRegion <- function(w) {
        combineTwasWeights(group("c1,c2", w), group("c1,c3", w + 1))
    }
    out <- pecotmr:::.twasMergeResultsByKey(
        list(perRegion(0.1), perRegion(0.2)),
        c("rA", "rB")
    )
    expect_equal(nrow(out), 2L)
})

test_that(".mergeEntryForKey gives NULL for a region lacking the tuple", {
    # NULL, not an empty row: a region that never fitted the tuple must not
    # dilute the merged entry. Only the fine-mapping merge keys rows this way.
    data(qtlFineMappingExample)
    fm <- qtlFineMappingExample
    key <- pecotmr:::.mergeResultKeyOf(fm)[[1L]]
    expect_false(is.null(pecotmr:::.mergeEntryForKey(fm, key)))
    expect_null(pecotmr:::.mergeEntryForKey(fm, "nope\rx\ry\rz"))
})


# ---------------------------------------------------------------------------
# .buildJointSumstatZMatrix: LD-panel cutoffs on the joint RSS path.
#
# Without these, mafCutoff would filter univariate RSS fits but silently not
# joint ones -- an inconsistency that is hard to spot in results, since both
# return a well-formed fit either way.
# ---------------------------------------------------------------------------

# @noRd
.bjz_qcd <- function() {
    data(qtlSumStatsMulticontextExample, envir = environment())
    suppressWarnings(suppressMessages(
        summaryStatsQc(qtlSumStatsMulticontextExample)
    ))
}
# Built once on first use: summaryStatsQc() dominates these tests.
delayedAssign(".bjzSs", .bjz_qcd())

# @noRd
.bjz_build <- function(ss, ...) {
    suppressMessages(.buildJointSumstatZMatrix(
        ss,
        seq_len(nrow(ss)),
        as.character(ss$context),
        "probe",
        ldReference = ldReference(ss)[[1L]],
        cutoffs = .panelCutoffs(PanelFilterParam(...))
    ))
}

test_that(".buildJointSumstatZMatrix filters nothing by default", {
    ss <- .bjzSs
    jz <- .bjz_build(ss)
    expect_equal(length(jz$variantIds), sum(lengths(ss)) / nrow(ss))
    expect_equal(ncol(jz$Z), nrow(ss))
})

test_that(".buildJointSumstatZMatrix drops panel-rare variants", {
    ss <- .bjzSs
    full <- .bjz_build(ss)
    tight <- .bjz_build(ss, mafCutoff = 0.2)
    expect_lt(length(tight$variantIds), length(full$variantIds))
    # Variants are rows; every context survives as a column.
    expect_equal(nrow(tight$Z), length(tight$variantIds))
    expect_equal(ncol(tight$Z), ncol(full$Z))
    expect_equal(colnames(tight$Z), colnames(full$Z))
})

test_that(".buildJointSumstatZMatrix agrees with .panelVariantFilter", {
    ss <- .bjzSs
    ids <- .bjz_build(ss)$variantIds
    for (cut in c(0.05, 0.2)) {
        expect_equal(
            length(.bjz_build(ss, mafCutoff = cut)$variantIds),
            length(.panelVariantFilter(
                ldReference(ss)[[1L]],
                ids,
                PanelFilterParam(mafCutoff = cut)
            )),
            label = str_c("mafCutoff ", cut)
        )
    }
})

test_that(".buildJointSumstatZMatrix keeps Z aligned to its ids", {
    # Each context's own frame is narrowed to the filtered set before the
    # shared-SNP-order check; otherwise that check fires on a length mismatch
    # the filter itself created.
    ss <- .bjzSs
    full <- .bjz_build(ss)
    tight <- .bjz_build(ss, mafCutoff = 0.2)
    expect_equal(rownames(tight$Z), tight$variantIds)
    keep <- is_in(full$variantIds, tight$variantIds)
    expect_equal(tight$Z[, 1L], full$Z[keep, 1L])
    expect_equal(tight$nVec, full$nVec)
})

test_that(".buildJointSumstatZMatrix honours a missingness cutoff", {
    ss <- .bjzSs
    expect_lt(
        length(.bjz_build(ss, imissCutoff = 0)$variantIds),
        length(.bjz_build(ss)$variantIds)
    )
})

test_that(".fmScopeUserLimit answers NULL for an absent limit", {
    expect_null(pecotmr:::.fmScopeUserLimit(NULL, "s1"))
})

test_that(".msFoldStudyResults answers NULL when every fold is empty", {
    expect_null(pecotmr:::.msFoldStudyResults(list(), NULL))
    expect_null(pecotmr:::.msFoldStudyResults(list(NULL, NULL), NULL))
})

test_that("the joint builders filter genotypes on the samples fit", {
    skip_if_not_installed("pgenlibr")
    qd <- .fitSamplesDataset()
    rp <- ResidualizationParam()
    filter <- pecotmr:::.qtlResolveFilter(GenotypeFilterParam(mafCutoff = 0.05))
    pf <- pecotmr:::.jointPhenotypeFilter(list())
    window <- .fitSamplesWindow()
    # g2 misses five samples, which a fit on it cannot use.
    contexts <- .js_crossContextXy(
        qd,
        "g2",
        c("brain", "liver"),
        NULL,
        0,
        "lbl",
        region = window,
        residualizationParam = rp,
        genotypeFilter = filter,
        phenotypeFilter = pf
    )
    traits <- .js_crossTraitXy(
        qd,
        "brain",
        c("g1", "g2"),
        NULL,
        0,
        "lbl",
        "s",
        region = window,
        residualizationParam = rp,
        genotypeFilter = filter,
        phenotypeFilter = pf
    )
    composed <- .js_composedXy(
        qd,
        list(
            contexts = list(s = c("brain", "liver")),
            traits = list(s = c("g1", "g2"))
        ),
        "s",
        NULL,
        0,
        "lbl",
        region = window,
        residualizationParam = rp,
        genotypeFilter = filter,
        phenotypeFilter = pf
    )
    fitSamples <- rownames(contexts$Y)
    expect_length(fitSamples, 25L)
    for (xy in list(contexts, traits, composed)) {
        expect_setequal(rownames(xy$Y), fitSamples)
        .expectFitOn(
            xy$X,
            fitSamples,
            pecotmr:::.fmResidGeno(
                qd,
                rp,
                filter,
                contexts = "brain",
                region = window,
                samples = fitSamples
            )
        )
    }
})

test_that("a joint fit with too few samples to fit reads no genotypes", {
    skip_if_not_installed("pgenlibr")
    qd <- .fitSamplesDataset()
    local_mocked_bindings(
        .handleDosage = function(...) abort("read genotypes"),
        .fmResidPheno = function(data, contexts, traitId, ...) {
            y <- matrix(
                c(1, NA),
                ncol = length(traitId),
                nrow = 2L,
                dimnames = list(c("s1", "s2"), traitId)
            )
            if (length(contexts) == 1L) {
                return(y)
            }
            set_names(rep(list(y), length(contexts)), contexts)
        },
        .package = "pecotmr"
    )
    filter <- pecotmr:::.qtlResolveFilter()
    expect_null(.js_crossContextXy(
        qd,
        "g1",
        c("brain", "liver"),
        NULL,
        1,
        "lbl",
        region = .fitSamplesWindow(),
        residualizationParam = ResidualizationParam(),
        genotypeFilter = filter
    )) |>
        suppressMessages()
    expect_null(.js_crossTraitXy(
        qd,
        "brain",
        c("g1", "g2"),
        NULL,
        0,
        "lbl",
        "s",
        region = .fitSamplesWindow(),
        residualizationParam = ResidualizationParam(),
        genotypeFilter = filter
    ))
    expect_null(.js_composedXy(
        qd,
        list(
            contexts = list(s = "brain"),
            traits = list(s = c("g1", "g2"))
        ),
        "s",
        NULL,
        0,
        "lbl",
        region = .fitSamplesWindow(),
        residualizationParam = ResidualizationParam(),
        genotypeFilter = filter
    ))
})

test_that("the joint builders read phenotypes under the run's phenotype QC", {
    # The cross-context, cross-trait and composed assemblies each read their
    # phenotypes through .fmResidPheno(); the run's resolved filter has to
    # reach it, or a joint fit would ignore the QC every other fit applies.
    seen <- list()
    samp <- c("s1", "s2", "s3")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(e = matrix(0, 2, 3, dimnames = list(c("t1", "t2"), samp)))
    )
    local_mocked_bindings(
        .crossContextPerTrait = function(...) c("c1", "c2"),
        .crossTraitTooFew = function(...) FALSE,
        .qtlContextSe = function(data, context) se,
        .fmResidGeno = function(...) {
            matrix(0, 3, 1, dimnames = list(samp, "v1"))
        },
        .fmResidPheno = function(
            data,
            contexts,
            traitId,
            phenotypeFilter,
            ...
        ) {
            seen[[length(seen) + 1L]] <<- phenotypeFilter
            matrix(0, 3, 1, dimnames = list(samp, "t1"))
        },
        .crossContextFitSamples = function(...) samp,
        .composedFitSamples = function(...) samp,
        .composedTuples = function(...) {
            list(
                list(context = "c1", trait = "t1"),
                list(context = "c2", trait = "t2")
            )
        },
        .crossContextAssemble = function(...) "cc",
        .crossTraitAssemble = function(...) "ct",
        .package = "pecotmr"
    )
    pf <- list(naAction = "impute", outlierAction = "drop")
    common <- list(
        residualizationParam = ResidualizationParam(),
        genotypeFilter = pecotmr:::.qtlResolveFilter(),
        phenotypeFilter = pf
    )
    .js_crossContextXy(
        NULL,
        "t1",
        c("c1", "c2"),
        1e5,
        0,
        "lbl",
        residualizationParam = common$residualizationParam,
        genotypeFilter = common$genotypeFilter,
        phenotypeFilter = pf
    )
    .js_crossTraitXy(
        NULL,
        "c1",
        c("t1", "t2"),
        1e5,
        0,
        "lbl",
        "S",
        residualizationParam = common$residualizationParam,
        genotypeFilter = common$genotypeFilter,
        phenotypeFilter = pf
    )
    pecotmr:::.composedPlan(
        NULL,
        NULL,
        "S",
        0,
        "lbl",
        common$residualizationParam,
        pf
    )
    expect_length(seen, 3L)
    for (got in seen) {
        expect_identical(got, pf)
    }
})

test_that(".twasMergeResultsByKey keeps a position its regions agree on", {
    region <- function(pos, idx) {
        pecotmr:::.withJointIndex(
            TwasWeights(
                studyName = "S",
                context = "c1",
                trait = "t1",
                method = "mrmash",
                entry = list(
                    TwasWeightsRow(sprintf("chr1:%d:A:G", pos), 0.1)
                ),
                jointContexts = "c1;c2"
            ),
            idx
        )
    }
    out <- pecotmr:::.twasMergeResultsByKey(
        list(region(100L, 2L), region(200L, 2L)),
        c("rA", "rB")
    )
    expect_identical(as.list(S4Vectors::mcols(out)$jointIndex), list(2L))
    # Regions whose screens kept different conditions leave it unknown.
    out <- pecotmr:::.twasMergeResultsByKey(
        list(region(100L, 1L), region(200L, 2L)),
        c("rA", "rB")
    )
    expect_null(S4Vectors::mcols(out)$jointIndex)
})

test_that("a merged row's joint position is NA where no result holds it", {
    local_mocked_bindings(.mergeResultKeyOf = function(r) "a")
    expect_true(is.na(pecotmr:::.fmJointIndexForKey(NULL, "b")))
})

test_that("a composed spec with the study axis says which studies sit out", {
    composed <- list(list(axes = c("study", "context")))
    expect_message(
        pecotmr:::.twasMultiStudyWarnExcluded(
            composed,
            list(A = NULL, B = NULL),
            1L
        ),
        "excluding individual-level studies \\(A, B\\)"
    )
    # The study axis alone, or a quiet run, says nothing.
    expect_silent(pecotmr:::.twasMultiStudyWarnExcluded(
        list(list(axes = "study")),
        list(A = NULL),
        1L
    ))
    expect_silent(pecotmr:::.twasMultiStudyWarnExcluded(
        composed,
        list(A = NULL),
        0L
    ))
})
