context("fineMappingPipeline")

# ===========================================================================
# Strategy
# ----------------------------------------------------------------------------
# Mock the SuSiE fitters (.fmFitSusieIndiv / .fmFitSusieRss) and the post-
# processor (.fmPostprocessOne) so the pipeline orchestration runs end-to-
# end without firing real susieR / susie_rss / postprocess_finemapping_fits
# calls. The fixture uses a small in-memory QtlDataset (mocked
# extractBlockGenotypes) and small QtlSumStats / GwasSumStats collections
# (with mocked extractBlockGenotypes for the LD sketch).
# ===========================================================================

.fmp_makeHandle <- function(snp_n = 6L, n_samples = 40L) {
    new(
        "GenotypeHandle",
        path = "/tmp/fmsketch",
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

.fmp_mockExtractor <- function(seed = 3, n_samples = 40L) {
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

.fmp_makeSe <- function(
    traits = c("ENSG_A", "ENSG_B"),
    n_samples = 40L,
    starts = NULL
) {
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

.fmp_makeQtlDataset <- function(
    contexts = "brain",
    traits = c("ENSG_A", "ENSG_B")
) {
    gh <- .fmp_makeHandle()
    phen <- setNames(
        map(contexts, function(.) .fmp_makeSe(traits = traits)),
        contexts
    )
    QtlDataset(
        studyName = "study1",
        genotypes = gh,
        phenotypes = phen,
        genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0)
    )
}

.fmp_makeSumstatsGr <- function(
    snp_ids = sprintf("chr1:%d:A:G", 100L * (1:5))
) {
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = length(snp_ids)),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = snp_ids,
        A1 = rep("A", length(snp_ids)),
        A2 = rep("G", length(snp_ids)),
        Z = rnorm(length(snp_ids)),
        N = rep(1000L, length(snp_ids))
    )
    gr
}

.fmp_makeQtlSumStats <- function(qc = TRUE, traitPos = NULL) {
    QtlSumStats(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        entry = list(.fmp_makeSumstatsGr()),
        genome = "hg19",
        ldReference = .fmp_makeHandle(),
        qcOptions = if (qc) list(step1 = "ok") else list(),
        traitPos = traitPos
    )
}

# A trait position for the summary statistics above: a MultiStudyQtlDataset
# lays its `cisWindow` there too, and 1000 bp around it keeps every variant.
.fmp_traitPos <- function() {
    GenomicRanges::GRanges("chr1:300")
}

.fmp_makeGwasSumStats <- function(qc = TRUE, studyName = "G1") {
    GwasSumStats(
        studyName = studyName,
        entry = list(.fmp_makeSumstatsGr()),
        genome = "hg19",
        ldReference = .fmp_makeHandle(),
        qcOptions = if (qc) list(step1 = "ok") else list()
    )
}

# The fixtures most tests read, built once: every constructor validates its
# whole collection, so rebuilding them in each test was a large share of this
# file's run time. Tests only read them (a change makes a local copy).
.fmp_qdA <- .fmp_makeQtlDataset(contexts = "brain", traits = "ENSG_A")
.fmp_qdAB <- .fmp_makeQtlDataset()
.fmp_qss <- .fmp_makeQtlSumStats()
.fmp_gss <- .fmp_makeGwasSumStats()
.fmp_qdTwoCtxA <- .fmp_makeQtlDataset(
    contexts = c("brain", "liver"),
    traits = "ENSG_A"
)
.fmp_qdTwoCtxAB <- .fmp_makeQtlDataset(contexts = c("brain", "liver"))
# One individual-level study plus embedded summary statistics placed at
# .fmp_traitPos(), so a `cisWindow` can be laid on both.
.fmp_msqd <- MultiStudyQtlDataset(
    qtlDatasets = list(study1 = .fmp_qdA),
    sumStats = .fmp_makeQtlSumStats(traitPos = .fmp_traitPos())
)

# The result class each summary-statistics input comes back as.
.fmp_resultClass <- c(
    GwasSumStats = "GwasFineMappingResult",
    QtlSumStats = "QtlFineMappingResult"
)

# One susie fit of a summary-statistics input under `susieRssParam`.
.fmp_fitRss <- function(input, susieRssParam = SusieRssParam()) {
    suppressMessages(fineMappingPipeline(
        input,
        methods = "susie",
        initializeWithSusieInf = FALSE,
        susieRssParam = susieRssParam
    ))
}

# Mocks for the SuSiE fitters + post-processor. Return tiny payloads keyed
# only by the token so post-process knows what to wrap. The `userArgs`
# parameter (per-method kwargs merged in by .fmMergeUserArgs) is recorded
# when a `capture` environment is given and otherwise ignored -- the mocks
# don't simulate downstream susie behaviour.
.fmp_mockFitIndiv <- function(capture = NULL) {
    function(
        X,
        y,
        token,
        chainFromInf = NULL,
        coverage = 0.95,
        userArgs = NULL
    ) {
        if (!is.null(capture)) {
            capture$userArgs <- userArgs
        }
        list(token = token, X_cols = ncol(X))
    }
}

.fmp_mockFitRss <- function() {
    function(
        z,
        R,
        n,
        token,
        chainFromInf = NULL,
        coverage = 0.95,
        userArgs = NULL,
        rFinite = NULL,
        rMismatch = "none",
        rssControl = NULL,
        ...
    ) {
        list(token = token, n_variants = length(z))
    }
}

# RSS fitter mock that carries susieR-style finite-sample R diagnostics so the
# SER-fallback branch in .fmFitRssBlock can be exercised end-to-end. `flag`
# controls R_reliability_flag; `capture` (an environment) records the finite/EB
# args the pipeline forwarded so tests can assert on them.
.fmp_mockFitRssDiag <- function(flag = FALSE, capture = NULL) {
    function(
        z,
        R,
        n,
        token,
        chainFromInf = NULL,
        coverage = 0.95,
        userArgs = NULL,
        rFinite = NULL,
        rMismatch = "none",
        rssControl = NULL,
        ...
    ) {
        if (!is.null(capture)) {
            capture$rFinite <- rFinite
            capture$rMismatch <- rMismatch
            capture$rssControl <- rssControl
            capture$userArgs <- userArgs
        }
        fit <- list(token = token, n_variants = length(z), multiTag = "multi")
        fit$R_finite_diagnostics <- list(
            R_reliability_flag = flag,
            ser_model = list(token = token, serTag = "ser")
        )
        fit
    }
}

.fmp_mockFitSer <- function() {
    function(z, n, coverage = 0.95, userArgs = NULL) {
        list(token = "ser", n_variants = length(z))
    }
}

.fmp_mockPostprocess <- function() {
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
        # Capture the requesting method on the FineMappingRow so the test can
        # verify the right dispatch happened.
        if (is.matrix(dataX)) {
            vids <- colnames(dataX)
        } else {
            vids <- names(dataY)
            if (is.null(vids) && is.list(dataY) && !is.null(dataY$z)) {
                vids <- names(dataY$z)
            }
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

# ===========================================================================
# .fmNormalizeMethods
# ===========================================================================

test_that(".fmNormalizeMethods: rejects NULL / empty / non-character/list", {
    expect_error(
        pecotmr:::.fmNormalizeMethods(NULL, inputKind = "QtlDataset"),
        "non-empty character"
    )
    expect_error(
        pecotmr:::.fmNormalizeMethods(character(0), inputKind = "QtlDataset"),
        "non-empty character"
    )
    expect_error(
        pecotmr:::.fmNormalizeMethods(42L, inputKind = "QtlDataset"),
        "must be a character vector, a named list"
    )
})

test_that(".fmNormalizeMethods: char-vector form deduplicates + seeds susie L defaults", {
    res <- pecotmr:::.fmNormalizeMethods(
        c("susie", "susie", "susieInf"),
        inputKind = "QtlDataset"
    )
    expect_equal(res$tokens, c("susie", "susieInf"))
    expect_equal(names(res$methodArgs), c("susie", "susieInf"))
    # SuSiE-family tokens get the pipeline L default (pecotmr owns these, not the
    # CLI wrappers). L_greedy defaults to NULL (greedy off), so it is left unset.
    expect_equal(res$methodArgs$susie$L, 10L)
    expect_null(res$methodArgs$susie$L_greedy)
    expect_equal(res$methodArgs$susieInf$L, 10L)
    # Non-susie-family tokens are left untouched.
    expect_length(
        pecotmr:::.fmNormalizeMethods(
            "mvsusie",
            inputKind = "QtlDataset"
        )$methodArgs$mvsusie,
        0L
    )
})

test_that(".fmNormalizeMethods: a methods record keeps kwargs + fills susie L", {
    res <- pecotmr:::.fmNormalizeMethods(
        FineMappingMethodsParam(
            qtlDatasetMethods = list(
                susie = list(L = 1, refine = FALSE),
                susieInf = list()
            )
        ),
        inputKind = "QtlDataset"
    )
    expect_equal(res$tokens, c("susie", "susieInf"))
    expect_equal(res$methodArgs$susie$L, 1) # explicit kwarg wins
    expect_false(res$methodArgs$susie$refine)
    expect_null(res$methodArgs$susie$L_greedy) # greedy off by default (Lgreedy = NULL)
    expect_equal(res$methodArgs$susieInf$L, 10L) # both filled
})

test_that(".fmNormalizeMethods: L / Lgreedy args override the susie defaults", {
    res <- pecotmr:::.fmNormalizeMethods(
        c("susie"),
        maxNumSingleEffects = 30L,
        maxNumSingleEffectsGreedy = 7L,
        inputKind = "QtlDataset"
    )
    expect_equal(res$methodArgs$susie$L, 30L)
    expect_equal(res$methodArgs$susie$L_greedy, 7L)
})

test_that(".fmNormalizeMethods: a bare named list is routed by input class", {
    # It used to be refused, on the grounds that options must be checked
    # against the engine receiving them. They still are -- but this call
    # knows its input class, so the overrides go into that path's slot and
    # are checked there. The constructor, which knows no path, still
    # refuses a plain list under `methods`.
    ind <- pecotmr:::.fmNormalizeMethods(
        list(susie = list(L = 1)),
        inputKind = "QtlDataset"
    )
    expect_equal(ind$tokens, "susie")
    expect_equal(ind$methodArgs$susie$L, 1)
    # And an argument belonging to the other path is caught, which neither
    # the old refusal nor the union-validated constructor managed.
    expect_error(
        pecotmr:::.fmNormalizeMethods(
            list(susie = list(z = 1)),
            inputKind = "QtlDataset"
        ),
        "unknown argument"
    )
    expect_equal(
        pecotmr:::.fmNormalizeMethods(
            list(susie = list(z = 1)),
            inputKind = "QtlSumStats"
        )$methodArgs$susie$z,
        1
    )
})

test_that("FineMappingMethodsParam rejects unnamed and non-list entries", {
    expect_error(
        FineMappingMethodsParam(qtlDatasetMethods = list(list(L = 1), list())),
        "must be named for its method"
    )
    expect_error(
        FineMappingMethodsParam(qtlDatasetMethods = list(susie = 42)),
        "must be that method's Options record or a list of its arguments"
    )
})

test_that(".fmMergeUserArgs: user overrides win over the base arguments", {
    # base sets convergence_method = "pip"; user overrides to "objective"
    out <- pecotmr:::.fmMergeUserArgs(
        list(z = 1:3, convergence_method = "pip"),
        userArgs = list(convergence_method = "objective", L = 1)
    )
    expect_equal(out$convergence_method, "objective")
    expect_equal(out$L, 1)
    expect_identical(out$z, 1:3)
    expect_identical(pecotmr:::.fmMergeUserArgs(list(a = 1)), list(a = 1))
})

# ===========================================================================
# .fmCheckMethodCapabilities
# ===========================================================================

test_that(".fmCheckMethodCapabilities: unknown token errors with full menu", {
    expect_error(
        pecotmr:::.fmCheckMethodCapabilities("bogus", "QtlDataset"),
        "unknown method token"
    )
})

test_that(".fmCheckMethodCapabilities: mrmash always rejected", {
    expect_error(
        pecotmr:::.fmCheckMethodCapabilities("mrmash", "QtlDataset"),
        "TWAS-weight-oriented"
    )
})

test_that(".fmCheckMethodCapabilities: fsusie on QtlSumStats rejected (no sumstatImpl)", {
    expect_error(
        pecotmr:::.fmCheckMethodCapabilities("fsusie", "QtlSumStats"),
        "individual-only"
    )
})

test_that(".fmCheckMethodCapabilities: mvsusie on GwasSumStats rejected", {
    expect_error(
        pecotmr:::.fmCheckMethodCapabilities("mvsusie", "GwasSumStats"),
        "not supported on GwasSumStats"
    )
})

# ===========================================================================
# .fmResolveSusieChain
# ===========================================================================

test_that(".fmResolveSusieChain: chains susie from susieInf when both are requested", {
    res <- pecotmr:::.fmResolveSusieChain(
        c("susieInf", "susie"),
        initializeWithSusieInf = TRUE
    )
    expect_true(res$chainSusie)
    expect_true(res$runInf)
    expect_true(res$keepInf)
})

test_that(".fmResolveSusieChain: keeps susieInf when explicitly requested", {
    res <- pecotmr:::.fmResolveSusieChain(
        c("susieInf", "susie"),
        initializeWithSusieInf = FALSE
    )
    expect_true(res$runInf)
    expect_true(res$keepInf)
})

test_that(".fmResolveSusieChain: no chain or susieInf run for susie alone", {
    res <- pecotmr:::.fmResolveSusieChain(
        c("susie"),
        initializeWithSusieInf = FALSE
    )
    expect_false(res$chainSusie)
    expect_false(res$runInf)
})

# ===========================================================================
# .fmCacheLookup / .fmCacheLookupGwas
# ===========================================================================

test_that(".fmCacheLookup: NULL fineMappingResult returns NULL", {
    expect_null(pecotmr:::.fmCacheLookup(NULL, "s1", "c1", "t1", "susie"))
})

test_that(".fmCacheLookup: returns matching entry by 4-tuple", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(token = "susie"),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    fmr <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    hit <- pecotmr:::.fmCacheLookup(fmr, "s1", "c1", "t1", "susie")
    # The entry is a derived view now, rebuilt from the element on access.
    expect_identical(
        pecotmr:::.fmrPartsVariantIds(hit),
        pecotmr:::.fmrPartsVariantIds(e)
    )
    expect_identical(
        pecotmr:::.fmrPartsSusieFit(hit),
        pecotmr:::.fmrPartsSusieFit(e)
    )
    expect_null(pecotmr:::.fmCacheLookup(fmr, "ghost", "c1", "t1", "susie"))
})

test_that(".fmCacheLookupGwas: matches on study/method/blockId", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(token = "susie"),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    # GwasFineMappingResult names a row with no blockId by its own span; the
    # lookup keys on that span (multi-block FMRs disambiguate per-block fits
    # by it).
    fmr <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e)
    )
    key <- pecotmr:::.rtlRangeKeys(fmr)[[1L]]
    hitGwas <- pecotmr:::.fmCacheLookupGwas(fmr, "g1", "susie", key)
    expect_identical(
        pecotmr:::.fmrPartsVariantIds(hitGwas),
        pecotmr:::.fmrPartsVariantIds(e)
    )
    expect_identical(
        pecotmr:::.fmrPartsSusieFit(hitGwas),
        pecotmr:::.fmrPartsSusieFit(e)
    )
    expect_null(pecotmr:::.fmCacheLookupGwas(fmr, "ghost", "susie", key))
})

test_that(".fmCacheLookupGwas: a one-block cache answers for that block only", {
    # A cache holding one block of a study is no fit for its other blocks, so
    # the range is compared even when only one row matches (study, method).
    fmr <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(.sc_makeFineMappingRow(3))
    )
    key <- pecotmr:::.rtlRangeKeys(fmr)[[1L]]
    expect_false(is.null(
        pecotmr:::.fmCacheLookupGwas(fmr, "g1", "susie", key)
    ))
    expect_null(
        pecotmr:::.fmCacheLookupGwas(fmr, "g1", "susie", "chr1_10100_10300")
    )
})

test_that(".fmCacheLookupGwas: the range key disambiguates multi-block FMRs", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(3, offset = 10000L)
    fmr <- GwasFineMappingResult(
        studyName = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    keys <- pecotmr:::.rtlRangeKeys(fmr)
    hit <- pecotmr:::.fmCacheLookupGwas(fmr, "g1", "susie", keys[[2L]])
    expect_identical(
        pecotmr:::.fmrPartsVariantIds(hit),
        pecotmr:::.fmrPartsVariantIds(e2)
    )
    # A key naming no block is a miss, not a silent first-row match.
    expect_null(
        pecotmr:::.fmCacheLookupGwas(fmr, "g1", "susie", "chr9_1_2")
    )
})

test_that(".fmCacheLookup: non-QtlFineMappingResult input returns NULL", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(token = "susie"),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    gwasFmr <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e)
    )
    expect_null(pecotmr:::.fmCacheLookup(gwasFmr, "g1", "c1", "t1", "susie"))
})

# ===========================================================================
# .fmBuildQtlResult / .fmBuildGwasResult — empty-entries errors
# ===========================================================================

test_that(".fmBuildQtlResult: empty entries errors", {
    expect_error(
        pecotmr:::.fmBuildQtlResult(
            character(0),
            character(0),
            character(0),
            character(0),
            list()
        ),
        "no \\(study, context, trait, method\\) tuples"
    )
})

test_that(".fmBuildGwasResult: empty entries errors", {
    expect_error(
        pecotmr:::.fmBuildGwasResult(character(0), character(0), list()),
        "no \\(study, method\\) tuples"
    )
})

# ===========================================================================
# .rbindFineMappingResult — class-check branches
# ===========================================================================

test_that(".rbindFineMappingResult: rejects non-FineMappingResultBase input", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(token = "susie"),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    fmr <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    expect_error(
        pecotmr:::.rbindFineMappingResult(fmr, "not_an_fmr"),
        "Must inherit from class 'FineMappingResultBase'"
    )
    expect_error(
        pecotmr:::.rbindFineMappingResult("not_an_fmr", fmr),
        "Must inherit from class 'FineMappingResultBase'"
    )
})

test_that(".rbindFineMappingResult: rejects mixed Qtl/Gwas inputs", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(token = "susie"),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    qtlFmr <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    gwasFmr <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e)
    )
    expect_error(
        pecotmr:::.rbindFineMappingResult(qtlFmr, gwasFmr),
        "inputs must be the same concrete class"
    )
})

test_that(".rbindFineMappingResult: concatenates two GwasFineMappingResult collections", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(token = "susie"),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    a <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        entry = list(e)
    )
    b <- GwasFineMappingResult(
        studyName = "g2",
        method = "susie",
        entry = list(e)
    )
    out <- pecotmr:::.rbindFineMappingResult(a, b)
    expect_s4_class(out, "GwasFineMappingResult")
    expect_equal(nrow(out), 2L)
})

test_that("combineFineMappingResults: row-binds same-class collections; rejects mixed class/empty", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(token = "susie"),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    g1 <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        blockId = "r1",
        entry = list(e)
    )
    # A second block has to occupy a different span: rows are identified by
    # (study, method, range), so two rows over the same variants would be the
    # same block regardless of their labels.
    g2 <- GwasFineMappingResult(
        studyName = "g1",
        method = "susie",
        blockId = "r2",
        entry = list(.sc_makeFineMappingRow(3, offset = 10000L))
    )
    out <- combineFineMappingResults(g1, g2) # variadic
    expect_s4_class(out, "GwasFineMappingResult")
    expect_equal(nrow(out), 2L)
    expect_equal(nrow(combineFineMappingResults(list(g1, g2))), 2L) # list form
    q <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    expect_error(combineFineMappingResults(g1, q), "same concrete class")
    expect_error(combineFineMappingResults(), "nothing to combine")
})

# ===========================================================================
# .fmExtractZn
# ===========================================================================

test_that(".fmExtractZn: errors on missing SNP / Z / N columns", {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100))
    expect_error(pecotmr:::.fmExtractZn(gr, "x"), "no SNP mcol")
    S4Vectors::mcols(gr)$SNP <- "chr1:100:A:G"
    expect_error(pecotmr:::.fmExtractZn(gr, "x"), "no Z mcol")
    S4Vectors::mcols(gr)$Z <- 1.0
    expect_error(pecotmr:::.fmExtractZn(gr, "x"), "no N mcol")
})

# ===========================================================================
# fineMappingPipeline(QtlDataset)
# ===========================================================================

test_that("fineMappingPipeline(QtlDataset): univariate rows record traitPos", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            .fmp_qdAB,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    # 1 context x 2 traits x 1 method = 2 rows.
    expect_equal(nrow(res), 2L)
    expect_setequal(methodNames(res), "susie")
    tp <- traitPosition(res)
    expect_s4_class(tp, "GRanges")
    expect_length(tp, nrow(res))
    expect_true(all(as.character(GenomicRanges::seqnames(tp)) == "chr1"))
    # TSS (start) is the trait's own rowRanges start from the QtlDataset
    # (.fmp_makeSe: ENSG_A -> 1000, ENSG_B -> 2000, width 500), NOT the
    # cisWindow-expanded fit region. Map by the row's trait column.
    tss <- setNames(GenomicRanges::start(tp), res$trait)
    ent <- setNames(GenomicRanges::end(tp), res$trait)
    expect_equal(tss[["ENSG_A"]], 1000L)
    expect_equal(tss[["ENSG_B"]], 2000L)
    expect_equal(ent[["ENSG_A"]], 1499L)
    expect_equal(ent[["ENSG_B"]], 2499L)
    # genomicRegion() reports the REALIZED variant span now, not the nominal cis
    # window (traitPos +/- cisWindow). The window is gone because it had no
    # correct update rule under subsetting; the span is in sync with the
    # variants by construction. It still sits inside the window the fit was
    # drawn from, which is what this checks.
    reg <- genomicRegion(res)
    rs <- setNames(GenomicRanges::start(reg), res$trait)
    re <- setNames(GenomicRanges::end(reg), res$trait)
    expect_gte(rs[["ENSG_A"]], 1L) # traitPos 1000 - cisWindow, clamped
    expect_lte(re[["ENSG_A"]], 2499L) # traitPos 1499 + cisWindow
    expect_gte(rs[["ENSG_B"]], 1000L) # traitPos 2000 - cisWindow
    expect_lte(re[["ENSG_B"]], 3499L) # traitPos 2499 + cisWindow
})

test_that(".fmAfForX: returns directional effect-allele af aligned to colnames(X)", {
    qd <- .fmp_qdA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .package = "pecotmr"
    )
    region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 100000L))
    # The helper aligns af to dimnames; values come from af over the same
    # selection. Columns deliberately reordered to test name-based alignment.
    X <- matrix(
        0,
        nrow = 5L,
        ncol = 3L,
        dimnames = list(
            paste0("s", 1:5),
            c("chr1:300:A:G", "chr1:100:A:G", "chr1:200:A:G")
        )
    )
    noFilter <- pecotmr:::.qtlResolveFilter()
    af <- pecotmr:::.fmAfForX(qd, X, genotypeFilter = noFilter, region = region)
    expect_length(af, 3L)
    expect_false(anyNA(af)) # region matched -> every fitted variant has an af
    blk <- pecotmr:::.qtlExtractBlock(
        qd,
        window = region,
        samples = rownames(X),
        filter = noFilter
    )
    expect_equal(af, unname(setNames(blk$af, blk$variantIds)[colnames(X)]))
})

test_that(".fmAfForX: returns NULL for an empty block or a non-QtlDataset source", {
    qd <- .fmp_qdA
    emptyX <- matrix(numeric(0), nrow = 0L, ncol = 0L)
    expect_null(pecotmr:::.fmAfForX(qd, emptyX))
    X <- matrix(
        0,
        nrow = 2L,
        ncol = 1L,
        dimnames = list(c("s1", "s2"), "chr1:100:A:G")
    )
    expect_null(pecotmr:::.fmAfForX(list(not = "a dataset"), X))
})

test_that("fineMappingPipeline(QtlDataset): threads directional af into postprocess", {
    # Regression for af = NA in credibleSets: the individual-level univariate path must
    # forward a non-NULL, directional effect-allele frequency to the
    # post-processor (which writes it into the topLoci `af` column).
    qd <- .fmp_qdA
    captured <- new.env(parent = emptyenv())
    captured$af <- "UNSET"
    captured$cols <- NULL
    recordingPostprocess <- function(
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
        captured$af <- af
        captured$cols <- colnames(dataX)
        vids <- colnames(dataX)
        FineMappingRow(
            variantIds = vids,
            susieFit = list(method = method),
            topLoci = data.frame(
                variant_id = vids,
                pip = seq(0.9, by = -0.1, length.out = length(vids)),
                af = af,
                stringsAsFactors = FALSE
            )
        )
    }
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = recordingPostprocess,
        .package = "pecotmr"
    )
    suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        )
    )
    # af forwarded (not left at the NULL default), one value per fitted variant.
    expect_false(identical(captured$af, "UNSET"))
    expect_false(is.null(captured$af))
    expect_length(captured$af, length(captured$cols))
    expect_true(any(!is.na(captured$af)))
    expect_true(all(captured$af >= 0 & captured$af <= 1, na.rm = TRUE))
})

test_that("fineMappingPipeline(QtlDataset): seed argument is accepted and runs", {
    qd <- .fmp_qdA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE,
            seed = 42L
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(nrow(res), 1L)
})

test_that(".fmTopPcScores: samples x min(nPCs, traits) scores, else NULL", {
    set.seed(7)
    n <- 30L
    Y <- matrix(
        rnorm(n * 3L),
        nrow = n,
        ncol = 3L,
        dimnames = list(paste0("s", seq_len(n)), c("ta", "tb", "tc"))
    )
    fn <- function(...) pecotmr:::.fmTopPcScores(...)
    # (a) 3 traits, nPCs >= traits -> 3 columns named topPC1..topPC3, rows = samples.
    sc <- fn(Y, 10L)
    expect_true(is.matrix(sc))
    expect_equal(ncol(sc), 3L)
    expect_equal(colnames(sc), c("topPC1", "topPC2", "topPC3"))
    expect_equal(nrow(sc), n)
    expect_equal(rownames(sc), rownames(Y))
    # (b) nPCs caps the number of returned columns.
    sc2 <- fn(Y, 2L)
    expect_equal(ncol(sc2), 2L)
    expect_equal(colnames(sc2), c("topPC1", "topPC2"))
    # (c) single-column Y -> NULL (PCA undefined for a single trait).
    expect_null(fn(Y[, 1L, drop = FALSE], 10L))
    # (d) a zero-variance trait is dropped; k reflects only the usable traits.
    Yzv <- cbind(Y, td = rep(1, n))
    scz <- fn(Yzv, 10L)
    expect_equal(ncol(scz), 3L) # td dropped -> still 3 usable traits
    expect_equal(colnames(scz), c("topPC1", "topPC2", "topPC3"))
    # (e) rows with any NA are dropped before PCA.
    Yna <- Y
    Yna[c(1L, 2L), 1L] <- NA
    scn <- fn(Yna, 10L)
    expect_equal(nrow(scn), n - 2L)
    expect_false(any(c("s1", "s2") %in% rownames(scn)))
    # (f) fewer than two complete rows, fewer than two traits with variance,
    # or nPCs = 0 leaves nothing to score.
    expect_null(fn(matrix(c(1, NA), 2, 2), 2L))
    expect_null(fn(cbind(a = rep(1, 4), b = rnorm(4)), 2L))
    expect_null(fn(Y, 0L))
})

test_that("fineMappingPipeline(QtlDataset, usePCA): rows keyed topPC{i}", {
    qd <- .fmp_makeQtlDataset(
        contexts = "brain",
        traits = c("ENSG_A", "ENSG_B", "ENSG_C")
    )
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE,
            usePCA = TRUE,
            nPCs = 2L
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_setequal(methodNames(res), "susie")
    pcRows <- as.character(res$trait) %in% c("topPC1", "topPC2")
    # 3 per-trait univariate susie rows + 2 top-PC rows = 5.
    expect_equal(sum(pcRows), 2L)
    expect_setequal(as.character(res$trait)[pcRows], c("topPC1", "topPC2"))
    expect_setequal(as.character(res$method)[pcRows], "susie")
})

test_that(".fmLookupMrmashFit: finds the mr.mash fit by (study, trait)", {
    mkEntry <- function(fits) {
        TwasWeightsRow(
            variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
            weights = c(0.1, 0.2),
            methodFits = fits
        )
    }
    payload <- list(
        dataDrivenPriorMatrices = list(U = list(a = diag(2))),
        w0 = c(a = 1),
        V = diag(2)
    )
    # The joint fit lives on the first mrmash row of the (study, trait) group;
    # the other context row carries fits = NULL. A non-mrmash row is ignored.
    tw <- TwasWeights(
        studyName = c("S", "S", "S"),
        context = c("c1", "c2", "c1"),
        trait = c("G", "G", "G"),
        method = c("mrmash", "mrmash", "enet"),
        entry = list(mkEntry(payload), mkEntry(NULL), mkEntry(payload))
    )
    lk <- function(...) pecotmr:::.fmLookupMrmashFit(...)
    expect_identical(lk(tw, "S", "G"), payload) # first non-NULL mrmash row
    expect_null(lk(tw, "S", "OTHER")) # no such trait
    expect_null(lk(tw, "OTHER", "G")) # no such study
    expect_null(lk(NULL, "S", "G")) # no TwasWeights supplied
})

test_that(".fmLookupMrmashCv: finds the per-fold CV payload by (study, trait)", {
    mkEntry <- function(cv) {
        TwasWeightsRow(
            variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
            weights = c(0.1, 0.2),
            cvResult = cv
        )
    }
    cv <- list(
        samplePartition = data.frame(Sample = "s1", Fold = 1L),
        foldFits = list(fold_1 = list(w0 = 1))
    )
    tw <- TwasWeights(
        studyName = c("S", "S"),
        context = c("c1", "c2"),
        trait = c("G", "G"),
        method = c("mrmash", "mrmash"),
        entry = list(mkEntry(cv), mkEntry(NULL))
    )
    lk <- function(...) pecotmr:::.fmLookupMrmashCv(...)
    expect_identical(lk(tw, "S", "G"), cv)
    expect_null(lk(tw, "S", "OTHER"))
    expect_null(lk(NULL, "S", "G"))
    # A cvResult without foldFits is not a per-fold prior payload -> NULL.
    tw2 <- TwasWeights(
        studyName = "S",
        context = "c1",
        trait = "G",
        method = "mrmash",
        entry = list(mkEntry(list(predictions = 1)))
    )
    expect_null(lk(tw2, "S", "G"))
})

test_that(".fmBuildMvsusiePriorCv: mode B reweights each fold's own fit", {
    sp <- data.frame(
        Sample = paste0("s", 1:6),
        Fold = rep(1:3, each = 2),
        stringsAsFactors = FALSE
    )
    mkFold <- function(uname, v) {
        list(
            dataDrivenPriorMatrices = list(
                U = setNames(list(diag(2)), uname),
                w = setNames(0.5, uname)
            ),
            w0 = setNames(0.5, paste0(uname, "_grid1")),
            V = diag(2) * v
        )
    }
    mvCv <- list(
        samplePartition = sp,
        foldFits = list(
            fold_1 = mkFold("A", 10),
            fold_2 = mkFold("B", 20),
            fold_3 = mkFold("C", 30)
        )
    )
    local_mocked_bindings(
        rescaleCovW0 = function(w0) setNames(1, sub("_grid1$", "", names(w0))),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        create_mixture_prior = function(...) "PRIOR",
        .package = "mvsusieR"
    )
    out <- pecotmr:::.fmBuildMvsusiePriorCv(
        mvCv,
        fullFitParts = NULL,
        conditionNames = c("c1", "c2")
    )
    expect_equal(names(out), c("1", "2", "3"))
    expect_equal(out[["1"]]$residualVariance, diag(2) * 10) # fold 1's own V
    expect_equal(out[["3"]]$residualVariance, diag(2) * 30) # fold 3's own V
})

test_that(".fmBuildMvsusiePriorCv: mode C reuses full-fit w0/V with per-fold U", {
    sp <- data.frame(
        Sample = paste0("s", 1:4),
        Fold = rep(1:2, each = 2),
        stringsAsFactors = FALSE
    )
    full <- list(
        dataDrivenPriorMatrices = list(U = list(Z = diag(2)), w = c(Z = 1)),
        w0 = c(Z_grid1 = 1),
        V = diag(2) * 99
    )
    # Fold stubs carry only U (no w0) -> mode C: override the full fit's U.
    mvCv <- list(
        samplePartition = sp,
        foldFits = list(
            fold_1 = list(
                dataDrivenPriorMatrices = list(U = list(Z = diag(2) * 2))
            ),
            fold_2 = list(
                dataDrivenPriorMatrices = list(U = list(Z = diag(2) * 3))
            )
        )
    )
    captured <- list()
    local_mocked_bindings(
        rescaleCovW0 = function(w0) setNames(1, sub("_grid1$", "", names(w0))),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        create_mixture_prior = function(...) {
            captured[[length(captured) + 1L]] <<- list(...)
            "PRIOR"
        },
        .package = "mvsusieR"
    )
    out <- pecotmr:::.fmBuildMvsusiePriorCv(
        mvCv,
        fullFitParts = full,
        conditionNames = c("c1", "c2")
    )
    expect_equal(names(out), c("1", "2"))
    expect_equal(out[["1"]]$residualVariance, diag(2) * 99) # full fit's V
    expect_equal(captured[[1]]$mixture_prior$matrices$Z, diag(2) * 2) # fold 1's U
    expect_equal(captured[[2]]$mixture_prior$matrices$Z, diag(2) * 3) # fold 2's U
})

test_that("fineMappingPipeline(QtlDataset): pipCutoffToSkip skips no-signal univariate traits", {
    qd <- .fmp_qdAB
    # Stateful screen: reject the first block (ENSG_A), keep the rest (ENSG_B).
    seen <- 0L
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .fmSerScreen = function(X, y, cutoff) {
            seen <<- seen + 1L
            seen > 1L
        },
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE,
            signalScreenParam = SignalScreenParam(pip = -1)
        )
    )
    # ENSG_A screened out, ENSG_B kept -> a single row.
    expect_equal(nrow(res), 1L)
    expect_setequal(traitNames(res), "ENSG_B")
})

test_that("fineMappingPipeline(QtlDataset): CV folds attach a cvResult", {
    skip_if_not_installed("susieR")
    qd <- .fmp_qdA
    # Real susie fits drive CV; the genotype extraction and the full-data
    # post-processor are mocked (the mock SNP ids are not chr:pos:a1:a2, which the
    # real buildTopLoci requires). CV runs its own real .fmFitSusieIndiv, so the
    # cvResult is genuine.
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE,
            crossValidationParam = CrossValidationParam(folds = 3),
            verbose = 0
        )
    )
    cv <- cvResult(
        res,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "susie"
    )
    expect_false(is.null(cv))
    expect_setequal(colnames(cv$samplePartition), c("Sample", "Fold"))
    expect_setequal(sort(unique(cv$samplePartition$Fold)), 1:3)
    expect_true("susie_predicted" %in% names(cv$prediction))
    expect_true("susie_performance" %in% names(cv$performance))
    # One out-of-fold prediction per sample (no fold leaves a sample unscored).
    expect_false(anyNA(cv$prediction[["susie_predicted"]]))
})

test_that("fineMappingPipeline(QtlDataset): no CV folds leaves cvResult NULL", {
    qd <- .fmp_qdA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        )
    )
    expect_null(cvResult(
        res,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "susie"
    ))
})

# ===========================================================================
# Cross-validation internals: .cvSamplePartition / .fmWeightsCv /
# .fmSliceCv / .fmAttachCv (unit-level counterparts to the cvFolds end-to-end
# tests above).
# ===========================================================================

test_that(".fmWeightsCv returns twasWeightsCv-shaped output keyed by snake method", {
    skip_if_not_installed("susieR")
    set.seed(42)
    n <- 60L
    p <- 12L
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(
            paste0("s", seq_len(n)),
            sprintf("chr1:%d:A:G", 100L * seq_len(p))
        )
    )
    y <- X[, 2] * 1.5 + rnorm(n, sd = 0.5)
    names(y) <- rownames(X)
    cv <- pecotmr:::.fmWeightsCv(
        X,
        y,
        tokens = "susie",
        methodArgs = list(susie = list()),
        fold = 3L,
        coverage = 0.95,
        verbose = 0
    )
    expect_named(
        cv,
        c("samplePartition", "prediction", "performance", "foldFits")
    )
    expect_setequal(colnames(cv$samplePartition), c("Sample", "Fold"))
    # Keyed by the TWAS snake method name (adapter methodKey base).
    expect_true("susie_predicted" %in% names(cv$prediction))
    expect_true("susie_performance" %in% names(cv$performance))
    pred <- cv$prediction[["susie_predicted"]]
    expect_equal(dim(pred), c(n, 1L))
    # Every sample is held out exactly once => no missing out-of-fold predictions.
    expect_false(anyNA(pred))
    perf <- cv$performance[["susie_performance"]]
    expect_equal(
        colnames(perf),
        c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
    )
    # A real causal signal should yield positive out-of-fold correlation.
    expect_gt(perf[1, "corr"], 0)
})

test_that(".fmWeightsCv reuses a supplied samplePartition verbatim", {
    skip_if_not_installed("susieR")
    set.seed(7)
    n <- 40L
    p <- 8L
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(
            paste0("s", seq_len(n)),
            sprintf("chr1:%d:A:G", 100L * seq_len(p))
        )
    )
    y <- X[, 1] + rnorm(n, sd = 0.5)
    names(y) <- rownames(X)
    part <- pecotmr:::.cvSamplePartition(rownames(X), fold = 4L)
    cv <- pecotmr:::.fmWeightsCv(
        X,
        y,
        tokens = "susie",
        methodArgs = list(susie = list()),
        fold = 4L,
        samplePartition = part,
        coverage = 0.95,
        verbose = 0
    )
    expect_identical(cv$samplePartition, part)
})

test_that(".fmSliceCv / .fmAttachCv slice one method and round-trip onto an entry", {
    full <- list(
        samplePartition = data.frame(Sample = c("s1", "s2"), Fold = c(1L, 2L)),
        prediction = list(
            susie_predicted = matrix(1, 2, 1),
            susie_inf_predicted = matrix(2, 2, 1)
        ),
        performance = list(
            susie_performance = matrix(0, 1, 6),
            susie_inf_performance = matrix(0, 1, 6)
        )
    )
    sl <- pecotmr:::.fmSliceCv(full, "susieInf")
    expect_identical(names(sl$prediction), "susie_inf_predicted")
    expect_identical(names(sl$performance), "susie_inf_performance")
    expect_identical(sl$samplePartition, full$samplePartition)

    tl <- data.frame(
        variant_id = "chr1:100:A:G",
        pip = 0.5,
        stringsAsFactors = FALSE
    )
    e <- FineMappingRow("chr1:100:A:G", list(), tl)
    e2 <- pecotmr:::.fmAttachCv(e, sl)
    expect_identical(pecotmr:::.fmrPartsCvResult(e2), sl)
})

test_that("fineMappingPipeline(QtlDataset): a TWAS-only method is rejected", {
    expect_error(
        fineMappingPipeline(.fmp_qdAB, methods = "mrmash"),
        "TWAS-weight-oriented"
    )
})

test_that("fineMappingPipeline(QtlDataset): unknown context errors", {
    expect_error(
        fineMappingPipeline(.fmp_qdAB, methods = "susie", contexts = "ghost"),
        "no study has: ghost"
    )
})

test_that("fineMappingPipeline(QtlDataset): empty trait filter errors", {
    expect_error(
        fineMappingPipeline(.fmp_qdA, methods = "susie", traitId = "ENSG_Z"),
        "no traits selected"
    )
})

test_that("fineMappingPipeline(QtlDataset): mvsusie with single trait/context rejected", {
    expect_error(
        fineMappingPipeline(.fmp_qdA, methods = "mvsusie"),
        "mvsusie requires multi-trait or multi-context"
    )
})

test_that("fineMappingPipeline(QtlDataset): fsusie with single trait rejected", {
    expect_error(
        fineMappingPipeline(.fmp_qdA, methods = "fsusie"),
        "fsusie requires multi-trait"
    )
})

# ===========================================================================
# fineMappingPipeline(QtlDataset): mvsusie + fsusie dispatch (mocked)
# ===========================================================================

# Mock mvsusieR::mvsusie / create_mixture_prior so the joint-fit branches run
# without actually fitting. Returns a stub fit object tagged with the input
# shape so the test can assert it was constructed as expected.
.fmp_mockMvsusie <- function() {
    function(X, Y, prior_variance, coverage) {
        list(token = "mvsusie", n_X_cols = ncol(X), n_Y_cols = ncol(Y))
    }
}
.fmp_mockMixturePrior <- function() {
    function(R, ...) list(R = R)
}
.fmp_mockSusiF <- function() {
    # L / cov_lev / min_purity are refused by FsusieOptions as pipeline-owned,
    # so the pipeline is their only source: record them here so a test can
    # assert they actually arrive. `...` keeps the stub tolerant of further
    # susiF arguments.
    function(X, Y, pos, L = NULL, cov_lev = NULL, min_purity = NULL, ...) {
        list(
            token = "fsusie",
            n_X_cols = ncol(X),
            n_Y_cols = ncol(Y),
            pos = pos,
            L = L,
            cov_lev = cov_lev,
            min_purity = min_purity
        )
    }
}

test_that("fineMappingPipeline(QtlDataset): mvsusie multi-trait single-context dispatch", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L)
    )
    expect_s4_class(res, "QtlFineMappingResult")
    # mvsusie multi-trait fans the joint fit out across both traits.
    expect_equal(nrow(res), 2L)
    expect_setequal(traitNames(res), c("ENSG_A", "ENSG_B"))
    expect_setequal(methodNames(res), "mvsusie")
    # Engine-built (multivariate) FM rows carry BOTH provenance columns, derived in
    # the shared .runJointCell seam: traitPos = the bare trait position, region =
    # traitPos +/- cisWindow (1000). This is the path the univariate pushRow
    # threading never touches, so it is the regression guard for the engine fix.
    tp <- traitPosition(res)
    reg <- genomicRegion(res)
    expect_s4_class(tp, "GRanges")
    expect_s4_class(reg, "GRanges")
    tpS <- setNames(GenomicRanges::start(tp), res$trait)
    regS <- setNames(GenomicRanges::start(reg), res$trait)
    regE <- setNames(GenomicRanges::end(reg), res$trait)
    expect_equal(tpS[["ENSG_A"]], 1000L)
    expect_equal(tpS[["ENSG_B"]], 2000L)
    # region is the realized variant span now, not traitPos +/- cisWindow.
    # This fixture mocks the fit, so its variants are not confined to the cis
    # window and containment is not a property to assert here -- what holds is
    # that every trait gets a well-formed span of its own.
    expect_true(all(regS <= regE))
    expect_setequal(names(regS), c("ENSG_A", "ENSG_B"))
    expect_false(identical(unname(regS), c(1L, 1000L))) # not the window
})

test_that("fineMappingPipeline(QtlDataset): mvsusie multi-context single-trait dispatch", {
    qd <- .fmp_qdTwoCtxA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L)
    )
    # Multi-context fan-out: one row per context for the shared trait.
    expect_equal(nrow(res), 2L)
    expect_setequal(contexts(res), c("brain", "liver"))
    expect_setequal(traitNames(res), "ENSG_A")
})

test_that("fineMappingPipeline(QtlDataset): pipCutoffToSkip drops null contexts before joint mvsusie", {
    qd <- .fmp_makeQtlDataset(
        contexts = c("brain", "liver", "heart"),
        traits = "ENSG_A"
    )
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        # Drop the middle context (liver); keep brain + heart.
        .fmSerScreenColumns = function(X, Y, cutoff) c(TRUE, FALSE, TRUE),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "mvsusie",
            cisWindow = 1000L,
            signalScreenParam = SignalScreenParam(pip = -1)
        )
    )
    # liver screened out -> the joint fit runs on brain + heart only.
    expect_equal(nrow(res), 2L)
    expect_setequal(contexts(res), c("brain", "heart"))
})

test_that("fineMappingPipeline(QtlDataset): pipCutoffToSkip skips mvsusie when < 2 contexts survive", {
    # susie runs alongside so the result still has rows after mvsusie is skipped.
    qd <- .fmp_qdTwoCtxA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .fmSerScreen = function(X, y, cutoff) TRUE, # keep univariate susie
        .fmSerScreenColumns = function(X, Y, cutoff) c(TRUE, FALSE), # only brain
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = c("susie", "mvsusie"),
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE,
            signalScreenParam = SignalScreenParam(pip = -1)
        )
    )
    # mvsusie skipped (only 1 context survives); susie still produced per-context.
    expect_setequal(methodNames(res), "susie")
    expect_false("mvsusie" %in% methodNames(res))
})

test_that("fineMappingPipeline(QtlDataset): mvsusie both multi falls back to per-context multi-trait", {
    qd <- .fmp_qdTwoCtxAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L)
    )
    # 2 contexts * 2 traits = 4 rows (joint fit reused per context).
    expect_equal(nrow(res), 4L)
    # Auto-detection now routes through the joint engine: each per-context group
    # is cross-trait, so every row tags its co-fit trait membership.
    expect_true("jointTraits" %in% pecotmr:::.tupleColumnNames(res))
    expect_true(all(grepl(
        "ENSG_A;ENSG_B|ENSG_B;ENSG_A",
        as.character(res$jointTraits)
    )))
})

test_that("fineMappingPipeline(QtlDataset): multi-trait auto-detection USES the data-driven mr.mash prior", {
    # Regression for the original bug: the old multi-trait path hardcoded the
    # canonical create_mixture_prior(R=ncol). Routed through the engine, it now
    # looks up a prior cross-trait mr.mash fit (study, context fixed; trait
    # match-any) and reweights from it.
    qd <- .fmp_qdAB
    fitParts <- list(
        dataDrivenPriorMatrices = list(U = list(K = diag(2)), w = c(K = 1)),
        w0 = c(K_grid1 = 1),
        V = diag(2)
    )
    mkE <- function() {
        TwasWeightsRow(
            variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
            weights = c(0.1, 0.2),
            methodFits = fitParts
        )
    }
    # What twasWeightsPipeline(jointSpecification='trait') emits: per-trait rows
    # each carrying the SHARED joint fit for (study1, brain).
    tw <- TwasWeights(
        studyName = c("study1", "study1"),
        context = c("brain", "brain"),
        trait = c("ENSG_A", "ENSG_B"),
        method = c("mrmash", "mrmash"),
        entry = list(mkE(), mkE())
    )
    sawMixturePrior <- FALSE
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        rescaleCovW0 = function(w0) c(K = 1),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = function(X, Y, prior_variance, coverage, ...) {
            list(token = "mvsusie")
        },
        create_mixture_prior = function(...) {
            if (!is.null(list(...)$mixture_prior)) {
                sawMixturePrior <<- TRUE
            }
            "PRIOR"
        },
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "mvsusie",
            cisWindow = 1000L,
            mrmashPrior = tw
        )
    )
    expect_equal(nrow(res), 2L)
    expect_true(sawMixturePrior) # data-driven prior built, not canonical
})

test_that("fineMappingPipeline(QtlDataset): multi-trait without mrmashPrior keeps the canonical prior", {
    qd <- .fmp_qdAB
    sawMixturePrior <- FALSE
    sawCanonical <- FALSE
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = function(...) {
            a <- list(...)
            if (!is.null(a$mixture_prior)) {
                sawMixturePrior <<- TRUE
            }
            if (!is.null(a$R)) {
                sawCanonical <<- TRUE
            }
            "PRIOR"
        },
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(qd, methods = "mvsusie", cisWindow = 1000L)
    )
    expect_false(sawMixturePrior)
    expect_true(sawCanonical)
})

test_that("fineMappingPipeline(QtlDataset): mvsusie resume cache short-circuits the joint fitter", {
    # All conditions of the cross-trait group are already in the prior partial
    # result -> the engine reuses the cached entries and never calls mvsusie.
    qd <- .fmp_qdAB
    cachedEntry <- function() {
        FineMappingRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
            susieFit = list(token = "mvsusie_cached"),
            topLoci = data.frame(
                variant_id = sprintf("chr1:%d:A:G", 100L * (1:3)),
                pip = c(0.9, 0.5, 0.1),
                stringsAsFactors = FALSE
            )
        )
    }
    cache <- QtlFineMappingResult(
        studyName = c("study1", "study1"),
        context = c("brain", "brain"),
        trait = c("ENSG_A", "ENSG_B"),
        method = c("mvsusie", "mvsusie"),
        entry = list(cachedEntry(), cachedEntry())
    )
    mv_calls <- 0
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = function(...) {
            mv_calls <<- mv_calls + 1L
            list()
        },
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "mvsusie",
            cisWindow = 1000L,
            fineMappingResult = cache
        )
    )
    expect_equal(mv_calls, 0L) # cache hit -> fitter never called
    expect_equal(nrow(res), 2L)
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='context' produces one joint row per trait", {
    qd <- .fmp_qdTwoCtxAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "mvsusie",
            cisWindow = 1000L,
            jointSpecification = "context"
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    # Per-context rows: each trait's 2-context joint emits 2 rows (4 total),
    # sharing the joint fit; jointContexts tags each with the co-fit membership.
    expect_equal(nrow(res), 4L)
    expect_true("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
    expect_setequal(as.character(res$context), c("brain", "liver"))
    expect_setequal(traitNames(res), c("ENSG_A", "ENSG_B"))
    expect_true(all(grepl(
        "brain;liver|liver;brain",
        as.character(res$jointContexts)
    )))
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='context' + univariate compose without double-fitting mvsusie", {
    qd <- .fmp_qdTwoCtxA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = c("susie", "mvsusie"),
            cisWindow = 1000L,
            jointSpecification = "context",
            initializeWithSusieInf = FALSE
        )
    )
    # susie -> 2 univariate rows (one per context); mvsusie joint over 2 contexts
    # -> 2 per-context rows sharing the joint fit. 4 rows total.
    expect_equal(nrow(res), 4L)
    expect_equal(sum(as.character(res$method) == "mvsusie"), 2L)
    expect_equal(sum(as.character(res$method) == "susie"), 2L)
    # Univariate rows have NA in jointContexts; both mvsusie rows carry membership.
    jc <- as.character(res$jointContexts)
    expect_equal(sum(is.na(jc)), 2L)
    expect_equal(sum(!is.na(jc)), 2L)
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='context' with only one context skips with message", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    expect_error(
        suppressMessages(
            fineMappingPipeline(
                qd,
                methods = "mvsusie",
                cisWindow = 1000L,
                jointSpecification = "context"
            )
        ),
        "no joint fits produced"
    )
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='trait' produces one joint row per context", {
    qd <- .fmp_qdTwoCtxAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "mvsusie",
            cisWindow = 1000L,
            jointSpecification = "trait"
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    # Per-trait rows: each context's 2-trait joint emits 2 rows (4 total).
    expect_equal(nrow(res), 4L)
    expect_true("jointTraits" %in% pecotmr:::.tupleColumnNames(res))
    expect_setequal(as.character(res$trait), c("ENSG_A", "ENSG_B"))
    expect_setequal(as.character(res$context), c("brain", "liver"))
})

test_that("fineMappingPipeline(QtlDataset): jointSpec='trait' with fsusie wires fsusieR::susiF", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        susiF = .fmp_mockSusiF(),
        .package = "fsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "fsusie",
            cisWindow = 1000L,
            jointSpecification = "trait"
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    # fsusie fits ONE functional model across the traits and stores it once.
    # Its curve spans the trait axis -- each trait is a POSITION on a shared
    # grid -- so unlike mvsusie, whose posteriors carry a conditions axis
    # that .fmSliceCond() slices per condition, there is nothing to
    # partition: one row per trait would hold N copies of the same curve.
    # The full membership stays visible in the `jointTraits` column.
    expect_equal(nrow(res), 1L)
    expect_true(all(as.character(res$method) == "fsusie"))
    expect_equal(as.character(res$trait), "ENSG_A")
    expect_equal(as.character(res$jointTraits), "ENSG_A;ENSG_B")
    expect_true("jointTraits" %in% pecotmr:::.tupleColumnNames(res))
})

test_that("fineMappingPipeline(QtlDataset): fsusie receives the CredibleSetParam settings", {
    # Regression: .fmOwnedEngineArgs() refuses L / cov_lev / min_purity in
    # FsusieOptions on the grounds that the pipeline supplies them, but for a
    # while nothing did -- the row described the retired fsusieWrapper() while
    # the pipeline called fitFsusie(). `--L` was silently ignored and susiF
    # ran at its own default. The mocked susiF records what it was handed.
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        susiF = .fmp_mockSusiF(),
        .package = "fsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "fsusie",
            cisWindow = 1000L,
            jointSpecification = "trait",
            credibleSetParam = CredibleSetParam(
                maxNumSingleEffects = 3L,
                coverage = 0.9,
                minAbsCorr = 0.4
            )
        )
    )
    # .fmp_mockPostprocess stores the raw fit under `payload`.
    recorded <- susieFit(res[1, ])$payload
    expect_equal(recorded$L, 3L)
    expect_equal(recorded$cov_lev, 0.9)
    expect_equal(recorded$min_purity, 0.4)
})

test_that("fineMappingPipeline(QtlDataset): fsusie multi-trait per context dispatch", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        susiF = .fmp_mockSusiF(),
        .package = "fsusieR"
    )
    res <- suppressMessages(
        fineMappingPipeline(qd, methods = "fsusie", cisWindow = 1000L)
    )
    expect_s4_class(res, "QtlFineMappingResult")
    # fsusie fits ONE functional model across the traits and stores it once.
    # Its curve spans the trait axis -- each trait is a POSITION on a shared
    # grid -- so unlike mvsusie, whose posteriors carry a conditions axis
    # that .fmSliceCond() slices per condition, there is nothing to
    # partition: one row per trait would hold N copies of the same curve.
    # The full membership stays visible in the `jointTraits` column.
    expect_equal(nrow(res), 1L)
    expect_equal(traitNames(res), "ENSG_A")
    expect_equal(as.character(res$jointTraits), "ENSG_A;ENSG_B")
    expect_setequal(methodNames(res), "fsusie")
})

# ===========================================================================
# fineMappingPipeline(MultiStudyQtlDataset)
# ===========================================================================

# Two individual-level studies of one trait in one context, built once.
.fmp_twoStudies <- MultiStudyQtlDataset(
    qtlDatasets = list(
        s1 = QtlDataset(
            studyName = "s1",
            genotypes = .fmp_makeHandle(),
            phenotypes = list(brain = .fmp_makeSe(traits = "ENSG_A")),
            genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0)
        ),
        s2 = QtlDataset(
            studyName = "s2",
            genotypes = .fmp_makeHandle(),
            phenotypes = list(brain = .fmp_makeSe(traits = "ENSG_A")),
            genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0)
        )
    )
)

test_that("fineMappingPipeline(MultiStudyQtlDataset): aggregates results across constituent QtlDatasets", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            .fmp_twoStudies,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    # One row per (study, context, trait, method) tuple.
    expect_equal(nrow(res), 2L)
    expect_setequal(studyName(res), c("s1", "s2"))
    # Pure individual-level -> ldReference should be NULL.
    expect_identical(ldReference(res), list())
})

# A study measured in brain and liver, on the synthetic genotypes.
.fmp_twoContextStudy <- function(studyName) {
    QtlDataset(
        studyName = studyName,
        genotypes = .fmp_makeHandle(),
        phenotypes = list(
            brain = .fmp_makeSe(traits = c("ENSG_A", "ENSG_B")),
            liver = .fmp_makeSe(traits = c("ENSG_A", "ENSG_B"))
        ),
        genotypeCovariates = matrix(numeric(0), nrow = 0, ncol = 0)
    )
}

test_that("fineMappingPipeline(MultiStudy): list selections stay per study", {
    # Every study used to get the selection whole: list(A = "ENSG_B") also
    # restricted B to ENSG_B, and list(B = "brain") restricted A to brain.
    mt <- MultiStudyQtlDataset(
        qtlDatasets = list(
            A = .fmp_twoContextStudy("A"),
            B = .fmp_twoContextStudy("B")
        )
    )
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        mt,
        methods = "susie",
        traitId = list(A = "ENSG_B"),
        contexts = list(B = "brain"),
        cisWindow = 1000L,
        initializeWithSusieInf = FALSE
    ))
    expect_setequal(
        paste(res$study, res$context, res$trait),
        c(
            "A brain ENSG_B",
            "A liver ENSG_B",
            "B brain ENSG_A",
            "B brain ENSG_B"
        )
    )
})

test_that("fineMappingPipeline(QtlDataset): trait chosen per context", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        .fmp_twoContextStudy("A"),
        methods = "susie",
        traitId = list(A = list(brain = "ENSG_B")),
        cisWindow = 1000L,
        initializeWithSusieInf = FALSE
    ))
    # brain keeps only ENSG_B; liver, which the list does not name, keeps both.
    expect_setequal(
        paste(res$context, res$trait),
        c("brain ENSG_B", "liver ENSG_A", "liver ENSG_B")
    )
})

test_that("fineMappingPipeline(MultiStudyQtlDataset): per-study options are not dropped", {
    # fitRetention was declared on the multi-study method and passed nowhere;
    # genotypeFilter / usePCA / nPCs were reachable only through `...`.
    capX <- new.env(parent = emptyenv())
    capSs <- new.env(parent = emptyenv())
    local_mocked_bindings(
        .fmPerStudy = function(qd, cfg) {
            capX$cfg <- cfg
            NULL
        },
        .fmSumStats = function(ss, cfg) {
            capSs$cfg <- cfg
            NULL
        },
        .package = "pecotmr"
    )
    expect_error(
        suppressMessages(fineMappingPipeline(
            .fmp_msqd,
            methods = "susie",
            cisWindow = 1000L,
            credibleSetParam = CredibleSetParam(coverage = 0.9),
            fitRetention = "full",
            genotypeFilterParam = GenotypeFilterParam(mafCutoff = 0.05),
            usePCA = TRUE,
            nPCs = 3L
        )),
        "no entries produced"
    )
    expect_equal(capX$cfg$fitRetention, "full")
    expect_equal(capX$cfg$genotypeFilterParam$mafCutoff, 0.05)
    expect_true(capX$cfg$usePCA)
    expect_equal(capX$cfg$nPCs, 3L)
    # The embedded QtlSumStats used to be post-processed at the default
    # coverage while its sibling QtlDatasets used the caller's.
    expect_equal(capSs$cfg$credibleSetParam$coverage, 0.9)
    expect_equal(capSs$cfg$fitRetention, "full")
})

test_that("fineMappingPipeline(MultiStudyQtlDataset): a misspelled option errors", {
    # The method no longer takes `...`, so a typo cannot be swallowed.
    expect_error(
        fineMappingPipeline(
            .fmp_msqd,
            methods = "susie",
            credibelSet = CredibleSetParam()
        ),
        "unused argument"
    )
})

test_that("fineMappingPipeline(MSQD): L / Lgreedy reach the per-study fit", {
    # One individual-level study and one embedded QtlSumStats, so both
    # per-study workers (.fmPerStudy and .fmSumStats) are exercised.
    capX <- new.env(parent = emptyenv())
    capRss <- new.env(parent = emptyenv())
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(capture = capX),
        .fmFitSusieRss = .fmp_mockFitRssDiag(capture = capRss),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    suppressMessages(
        fineMappingPipeline(
            .fmp_msqd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE,
            credibleSetParam = CredibleSetParam(
                maxNumSingleEffects = 4L,
                maxNumSingleEffectsGreedy = 2L
            )
        )
    )
    # The multi-study method had no L / Lgreedy formals of its own before the
    # bundle: they rode `...` into the per-study recursion, so a misspelling
    # was silently discarded and .fmMsResolveTokens -- which seeds the joint
    # phase -- never saw them at all. They are checked formals instead.
    expect_equal(capX$userArgs$L, 4L)
    expect_equal(capX$userArgs$L_greedy, 2L)
    expect_equal(capRss$userArgs$L, 4L)
    expect_equal(capRss$userArgs$L_greedy, 2L)
})

test_that("fineMappingPipeline(MultiStudyQtlDataset): jointRegions=FALSE merges per region in each study", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    regions <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(50L, 350L), c(350L, 650L))
    )
    res <- suppressMessages(fineMappingPipeline(
        .fmp_twoStudies,
        methods = "susie",
        traitId = "ENSG_A",
        region = regions,
        jointRegions = FALSE,
        initializeWithSusieInf = FALSE
    ))
    expect_s4_class(res, "QtlFineMappingResult")
    # one merged row per study (region collapsed into the entry).
    expect_equal(nrow(res), 2L)
    expect_setequal(studyName(res), c("s1", "s2"))
    fit <- susieFit(
        res,
        studyName = "s1",
        context = "brain",
        trait = "ENSG_A",
        method = "susie"
    )
    expect_equal(names(fit), c("region1", "region2"))
})

test_that("fineMappingPipeline(MSQD): embedded sumstats record ldReference", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            .fmp_msqd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_true(nrow(res) >= 2L)
    # Embedded sumstats has an LD sketch -> the merged result carries it.
    expect_s4_class(ldReference(res)[[1L]], "RangedSummarizedExperiment")
})

# ===========================================================================
# fineMappingPipeline(QtlSumStats) and fineMappingPipeline(GwasSumStats)
# ===========================================================================

test_that("fineMappingPipeline(*SumStats): runs with mocked RSS fitters", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    for (input in list(.fmp_qss, .fmp_gss)) {
        kind <- class(input)[[1L]]
        res <- .fmp_fitRss(input)
        expect_s4_class(res, .fmp_resultClass[[kind]])
        expect_equal(nrow(res), 1L, info = kind)
        expect_equal(methodNames(res), "susie", info = kind)
        # The result carries the build its summary statistics name.
        expect_identical(
            unique(unname(GenomeInfoDb::genome(res))),
            "hg19",
            info = kind
        )
    }
})

test_that("fineMappingPipeline(*SumStats): un-QCd input rejected", {
    unQcd <- list(
        .fmp_makeQtlSumStats(qc = FALSE),
        .fmp_makeGwasSumStats(qc = FALSE)
    )
    for (input in unQcd) {
        expect_error(
            fineMappingPipeline(input, methods = "susie"),
            "has no QC record",
            info = class(input)[[1L]]
        )
    }
})

test_that("fineMappingPipeline(QtlSumStats): empty selection rejected", {
    expect_error(
        fineMappingPipeline(.fmp_qss, methods = "susie", traitId = "ghost"),
        "no entries matched"
    )
    # A context no study has is caught earlier, when the selection is parsed.
    expect_error(
        fineMappingPipeline(.fmp_qss, methods = "susie", contexts = "ghost"),
        "no study has: ghost"
    )
})

# ===========================================================================
# fineMappingPipeline(GwasSumStats)
# ===========================================================================

test_that("fineMappingPipeline(GwasSumStats): a blockId gives no window", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    regionOf <- function(blockId) {
        gss <- GwasSumStats(
            studyName = "G1",
            entry = list(.fmp_makeSumstatsGr()),
            genome = "hg19",
            ldReference = .fmp_makeHandle(),
            blockId = blockId,
            qcOptions = list(step1 = "ok")
        )
        res <- suppressMessages(fineMappingPipeline(
            gss,
            methods = "susie",
            initializeWithSusieInf = FALSE
        ))
        reg <- genomicRegion(res)
        c(GenomicRanges::start(reg), GenomicRanges::end(reg))
    }
    # A blockId is a name, however it is spelled, so the window is the
    # entry's span; a window the caller means goes in `region`.
    expect_equal(regionOf("chr1_50_600"), c(100L, 500L))
    expect_equal(regionOf("block1"), c(100L, 500L))
})

test_that("fineMappingPipeline(GwasSumStats): a range keeps its origin", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    fitOver <- function(...) {
        gss <- GwasSumStats(
            studyName = "G1",
            entry = list(.fmp_makeSumstatsGr()),
            genome = "hg19",
            ldReference = .fmp_makeHandle(),
            qcOptions = list(step1 = "ok"),
            ...
        )
        suppressMessages(fineMappingPipeline(
            gss,
            methods = "susie",
            initializeWithSusieInf = FALSE
        ))
    }
    # A window named in a blockId column gives the block its range and key.
    named <- fitOver(
        region = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(1L, 2000L),
            blockId = "b1"
        )
    )
    expect_identical(named$blockId, "b1")
    expect_identical(as.character(genomicRegion(named)), "chr1:1-2000")
    expect_true(named$regionGiven)
    # Without one, the row keeps its seqname as its name, and its span bounds
    # the block, computed.
    bare <- fitOver()
    expect_identical(bare$blockId, "chr1")
    expect_identical(as.character(genomicRegion(bare)), "chr1:100-500")
    expect_false(bare$regionGiven)
})

test_that("fineMappingPipeline(GwasSumStats): threads per-variant N", {
    # Regression for top_loci$N == 1. The RSS path must forward the per-variant
    # effective N (the entry's N column, aligned to the analysed variants) to
    # the post-processor -- NOT the scalar median the SuSiE-RSS fit consumes,
    # never a list-collapsed 1.
    snp_ids <- sprintf("chr1:%d:A:G", 100L * (1:5))
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = 5L),
            width = 1L
        )
    )
    # median = 1200, all distinct
    perVarN <- c(1000L, 1100L, 1200L, 1300L, 1400L)
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = snp_ids,
        A1 = rep("A", 5L),
        A2 = rep("G", 5L),
        Z = rnorm(5L),
        N = perVarN
    )
    gss <- GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .fmp_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
    captured <- new.env(parent = emptyenv())
    captured$n <- "UNSET"
    recording <- function(
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
        n = NULL,
        ...
    ) {
        captured$n <- n
        vids <- names(dataY$z)
        if (is.null(vids)) {
            vids <- snp_ids
        }
        FineMappingRow(
            variantIds = vids,
            susieFit = list(method = method),
            topLoci = data.frame(
                variant_id = vids,
                pip = 0.5,
                stringsAsFactors = FALSE
            )
        )
    }
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = recording,
        .package = "pecotmr"
    )
    suppressMessages(
        fineMappingPipeline(
            gss,
            methods = "susie",
            initializeWithSusieInf = FALSE
        )
    )
    # per-variant N forwarded (all 5 distinct), not the scalar median, not 1.
    expect_false(identical(captured$n, "UNSET"))
    expect_equal(as.integer(captured$n), perVarN)
})

# ---- PIP-screen graceful skip: a screened region -> empty result, not error ----

# A collection whose rows QC audited as `audits` (one list per row).
.fmp_withAudits <- function(x, audits) {
    x$qcAudit <- S4Vectors::SimpleList(audits)
    x
}

# A GwasSumStats whose (single) entry was emptied by summaryStatsQc's PIP screen:
# 0-variant entry + an audit with pipScreenSkipped = TRUE (+ reason).
.fmp_makeScreenedGwas <- function(
    studyName = "G1",
    reason = "no signals above PIP threshold 0.025"
) {
    gss <- GwasSumStats(
        studyName = studyName,
        entry = list(GenomicRanges::GRanges()),
        genome = "hg19",
        ldReference = .fmp_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
    .fmp_withAudits(
        gss,
        list(list(pipScreenSkipped = TRUE, pipScreenReason = reason))
    )
}

test_that("fineMappingPipeline(GwasSumStats): a PIP-screened region yields an empty result, not an error", {
    gss <- .fmp_makeScreenedGwas()
    expect_message(
        res <- fineMappingPipeline(
            gss,
            methods = "susie",
            initializeWithSusieInf = FALSE
        ),
        "region skipped"
    )
    expect_s4_class(res, "GwasFineMappingResult")
    expect_equal(nrow(res), 0L) # graceful skip, not stop("no ... tuples")
    # An audit recording the skip but no reason gets the screen's own.
    noReason <- .fmp_withAudits(
        GwasSumStats(
            studyName = "G1",
            entry = list(GenomicRanges::GRanges()),
            genome = "hg19",
            ldReference = .fmp_makeHandle(),
            qcOptions = list(step1 = "ok")
        ),
        list(list(pipScreenSkipped = TRUE))
    )
    expect_message(
        fineMappingPipeline(
            noReason,
            methods = "susie",
            initializeWithSusieInf = FALSE
        ),
        "skipped: no signal above the PIP pre-screen cutoff"
    )
})

test_that("fineMappingPipeline(GwasSumStats): mixed screened + real keeps only the real region", {
    gss <- .fmp_withAudits(
        GwasSumStats(
            studyName = c("Gskip", "Greal"),
            entry = list(GenomicRanges::GRanges(), .fmp_makeSumstatsGr()),
            genome = "hg19",
            ldReference = .fmp_makeHandle(),
            qcOptions = list(step1 = "ok")
        ),
        list(
            list(pipScreenSkipped = TRUE, pipScreenReason = "no signal"),
            list()
        )
    )
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            gss,
            methods = "susie",
            initializeWithSusieInf = FALSE
        )
    )
    expect_equal(unique(as.character(res$study)), "Greal") # screened study absent
    expect_gt(nrow(res), 0L)
})

test_that("fineMappingPipeline(GwasSumStats): a 0-variant entry is skipped even without the flag", {
    gss <- GwasSumStats(
        studyName = "G1",
        entry = list(GenomicRanges::GRanges()),
        genome = "hg19",
        ldReference = .fmp_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
    res <- suppressMessages(
        fineMappingPipeline(
            gss,
            methods = "susie",
            initializeWithSusieInf = FALSE
        )
    )
    expect_s4_class(res, "GwasFineMappingResult")
    expect_equal(nrow(res), 0L)
})

test_that("fineMappingPipeline(QtlSumStats): a screened trait is skipped gracefully", {
    qss <- .fmp_withAudits(
        QtlSumStats(
            studyName = "Q1",
            context = "c1",
            trait = "t1",
            entry = list(GenomicRanges::GRanges()),
            genome = "hg19",
            ldReference = .fmp_makeHandle(),
            qcOptions = list(step1 = "ok")
        ),
        list(list(pipScreenSkipped = TRUE, pipScreenReason = "no signal"))
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qss,
            methods = "susie",
            initializeWithSusieInf = FALSE
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(nrow(res), 0L)
})

# ---- SER fallback (shared sumstat SuSiE-RSS path: GwasSumStats + QtlSumStats) ----
#
# Both summary-statistics inputs reach the same SuSiE-RSS block, so each test
# runs it once per input kind.

test_that("fineMappingPipeline(*SumStats): serFallback keeps a reliable fit", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRssDiag(flag = FALSE),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    for (input in list(.fmp_gss, .fmp_qss)) {
        kind <- class(input)[[1L]]
        sf <- susieFit(.fmp_fitRss(
            input,
            SusieRssParam(serFallback = TRUE, rMismatch = "eb")
        ))
        # Multi-effect fit reported (payload keeps the multi tag, not the ser
        # tag).
        expect_equal(sf$payload$multiTag, "multi", info = kind)
        expect_false(sf$serFallbackUsed, info = kind)
        expect_false(sf$R_reliability_flag, info = kind)
        # keepFullFit="fallback": non-fallback region carries no retained fit.
        expect_null(sf$multiEffectFit, info = kind)
    }
})

test_that("fineMappingPipeline(*SumStats): unreliable R falls back to SER", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRssDiag(flag = TRUE),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    for (input in list(.fmp_gss, .fmp_qss)) {
        kind <- class(input)[[1L]]
        sf <- susieFit(.fmp_fitRss(
            input,
            SusieRssParam(serFallback = TRUE, rMismatch = "eb")
        ))
        # Reported result is the single-effect ser_model.
        expect_equal(sf$payload$serTag, "ser", info = kind)
        expect_true(sf$serFallbackUsed, info = kind)
        expect_true(sf$R_reliability_flag, info = kind)
        # Fallback region retains the pre-fallback multi-effect fit.
        expect_false(is.null(sf$multiEffectFit), info = kind)
        expect_equal(sf$multiEffectFit$multiTag, "multi", info = kind)
    }
})

test_that("fineMappingPipeline(*SumStats): serFallback=FALSE changes nothing", {
    cap <- new.env(parent = emptyenv())
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRssDiag(flag = TRUE, capture = cap),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    for (input in list(.fmp_gss, .fmp_qss)) {
        kind <- class(input)[[1L]]
        sf <- susieFit(.fmp_fitRss(input))
        # With serFallback off (the default), the shared block never fires the
        # fallback nor forwards a non-default rFinite/rMismatch -- even though
        # the mock flags the fit unreliable -- so the multi-effect result is
        # kept and the fit carries NO diagnostic fields.
        expect_null(cap$rFinite, info = kind)
        expect_equal(cap$rMismatch, "none", info = kind)
        expect_equal(sf$payload$multiTag, "multi", info = kind)
        expect_false("R_reliability_flag" %in% names(sf), info = kind)
        expect_false("serFallbackUsed" %in% names(sf), info = kind)
        expect_false("multiEffectFit" %in% names(sf), info = kind)
    }
})

test_that("fineMappingPipeline(*SumStats): SusieRssParam reaches the fitter", {
    # rFinite defaults to the LD-panel sample size (40 in .fmp_makeHandle) and
    # an explicit one is forwarded verbatim; rMismatch -- eb_mix too, which
    # susieR >= 0.16.6 takes -- and the susie_rss_control() settings pass
    # through unrestricted, and with no control the fitter sees NULL.
    cap <- new.env(parent = emptyenv())
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRssDiag(capture = cap),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    control <- SusieRssControlOptions(
        check_prior = TRUE,
        mismatch_estimator = "map"
    )
    for (input in list(.fmp_gss, .fmp_qss)) {
        kind <- class(input)[[1L]]
        .fmp_fitRss(input, SusieRssParam(serFallback = TRUE, rMismatch = "eb"))
        expect_equal(cap$rFinite, 40L, info = kind)
        expect_equal(cap$rMismatch, "eb", info = kind)
        expect_null(cap$rssControl, info = kind)
        .fmp_fitRss(
            input,
            SusieRssParam(
                serFallback = TRUE,
                rFinite = 12345,
                rMismatch = "eb_mix",
                control = control
            )
        )
        expect_equal(cap$rFinite, 12345, info = kind)
        expect_equal(cap$rMismatch, "eb_mix", info = kind)
        expect_equal(
            cap$rssControl,
            list(check_prior = TRUE, mismatch_estimator = "map"),
            info = kind
        )
    }
})

test_that("fineMappingPipeline(*SumStats): keepFullFit='all' retains the fit", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRssDiag(flag = FALSE),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    keepAll <- SusieRssParam(
        serFallback = TRUE,
        rMismatch = "eb",
        keepFullFit = "all"
    )
    for (input in list(.fmp_gss, .fmp_qss)) {
        kind <- class(input)[[1L]]
        # A region that did not fall back still keeps its multi-effect fit.
        sf <- susieFit(.fmp_fitRss(input, keepAll))
        expect_false(is.null(sf$multiEffectFit), info = kind)
        expect_equal(sf$multiEffectFit$multiTag, "multi", info = kind)
    }
})

test_that("fineMappingPipeline(QtlSumStats): top-level L / Lgreedy reach the susie fit", {
    cap <- new.env(parent = emptyenv())
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRssDiag(capture = cap),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    suppressMessages(
        fineMappingPipeline(
            .fmp_qss,
            methods = "susie",
            initializeWithSusieInf = FALSE,
            credibleSetParam = CredibleSetParam(
                maxNumSingleEffects = 4L,
                maxNumSingleEffectsGreedy = 2L
            )
        )
    )
    # QtlDataset and GwasSumStats have always taken these; QtlSumStats did not,
    # so its susie fits were pinned to the L = 10 default.
    expect_equal(cap$userArgs$L, 4L)
    expect_equal(cap$userArgs$L_greedy, 2L)
})

test_that("fineMappingPipeline(QtlSumStats): per-method L beats the top-level L", {
    cap <- new.env(parent = emptyenv())
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRssDiag(capture = cap),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    suppressMessages(
        fineMappingPipeline(
            .fmp_qss,
            methods = list(susie = list(L = 7L)),
            initializeWithSusieInf = FALSE,
            credibleSetParam = CredibleSetParam(maxNumSingleEffects = 4L)
        )
    )
    expect_equal(cap$userArgs$L, 7L)
})

test_that("fineMappingPipeline(GwasSumStats): non-RSS family rejected by capability check", {
    expect_error(
        fineMappingPipeline(.fmp_gss, methods = "fsusie"),
        "not supported on GwasSumStats"
    )
})

# ---- ser (LD-free single-effect regression); sumstat-only, not GWAS-only ----

test_that("fineMappingPipeline(*SumStats): ser runs via .fmFitSusieSer", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieSer = .fmp_mockFitSer(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    for (input in list(.fmp_qss, .fmp_gss)) {
        kind <- class(input)[[1L]]
        res <- suppressMessages(fineMappingPipeline(input, methods = "ser"))
        expect_s4_class(res, .fmp_resultClass[[kind]])
        expect_equal(methodNames(res), "ser", info = kind)
    }
})

test_that("ser is sumstat-only: rejected on QtlDataset, allowed on all sumstat kinds", {
    expect_error(
        pecotmr:::.fmCheckMethodCapabilities("ser", "QtlDataset"),
        "sumstat-only"
    )
    expect_silent(pecotmr:::.fmCheckMethodCapabilities("ser", "QtlSumStats"))
    expect_silent(pecotmr:::.fmCheckMethodCapabilities("ser", "GwasSumStats"))
    expect_silent(pecotmr:::.fmCheckMethodCapabilities(
        "ser",
        "MultiStudyQtlDataset"
    ))
})

test_that(".fmNormalizeMethods does not inject L / L_greedy for ser", {
    norm <- pecotmr:::.fmNormalizeMethods(
        "ser",
        maxNumSingleEffects = 20L,
        maxNumSingleEffectsGreedy = 5L,
        inputKind = "QtlDataset"
    )
    expect_null(norm$methodArgs[["ser"]][["L"]])
    expect_null(norm$methodArgs[["ser"]][["L_greedy"]])
})

# ===========================================================================
# fineMappingPipeline(ANY)
# ===========================================================================

test_that("fineMappingPipeline(ANY): unsupported input class errors", {
    expect_error(
        fineMappingPipeline(matrix(0, 5, 5), methods = "susie"),
        "does not accept inputs of class 'matrix'"
    )
})

# ===========================================================================
# QtlSumStats: mvsusie single-context rejection and auto-detection, cache hits
# ===========================================================================

.fmp_makeMultiCtxQtlSumStats <- function() {
    # Multi-row QtlSumStats: 2 contexts x 1 trait, the shape a cross-context
    # joint fit takes.
    e1 <- .fmp_makeSumstatsGr()
    e2 <- .fmp_makeSumstatsGr()
    QtlSumStats(
        studyName = c("Q1", "Q1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        entry = list(e1, e2),
        genome = "hg19",
        ldReference = .fmp_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
}

test_that("fineMappingPipeline(QtlSumStats): mvsusie needs two contexts", {
    # Build a single-context-per-(study, trait) collection. Mvsusie requires
    # at least two contexts per (study, trait) group.
    ss <- QtlSumStats(
        studyName = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t2"),
        entry = list(.fmp_makeSumstatsGr(), .fmp_makeSumstatsGr()),
        genome = "hg19",
        ldReference = .fmp_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
    expect_error(
        fineMappingPipeline(ss, methods = "mvsusie"),
        "mvsusie requires at least two"
    )
})

test_that("fineMappingPipeline(QtlSumStats): mvsusie auto-joins contexts", {
    # Multi-context (c1, c2) single-trait collection: mvsusie without an explicit
    # jointSpecification routes through the joint engine -> per-context rows that
    # share the cross-context RSS fit, each tagged with the co-fit membership.
    ss <- .fmp_makeMultiCtxQtlSumStats()
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie_rss = function(Z, R, N, prior_variance, coverage, ...) {
            list(token = "mvsusie_rss")
        },
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    res <- suppressMessages(fineMappingPipeline(ss, methods = "mvsusie"))
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(nrow(res), 2L)
    expect_setequal(as.character(res$context), c("c1", "c2"))
    expect_true("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
    expect_true(all(grepl("c1;c2|c2;c1", as.character(res$jointContexts))))
})

test_that("fineMappingPipeline(QtlSumStats): cache hit short-circuits the RSS fitter", {
    ss <- .fmp_qss
    cachedEntry <- FineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
        susieFit = list(token = "susie_cached"),
        topLoci = data.frame(
            variant_id = sprintf("chr1:%d:A:G", 100L * (1:5)),
            pip = seq(0.9, 0.1, length.out = 5),
            stringsAsFactors = FALSE
        )
    )
    cache <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(cachedEntry),
        ldReference = .fmp_makeHandle()
    )
    rss_calls <- 0
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = function(...) {
            rss_calls <<- rss_calls + 1L
            .fmp_mockFitRss()(...)
        },
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            ss,
            methods = "susie",
            initializeWithSusieInf = FALSE,
            fineMappingResult = cache
        )
    )
    expect_equal(rss_calls, 0L)
    expect_equal(nrow(res), 1L)
})

# ===========================================================================
# GwasSumStats: cache hit by (study, method)
# ===========================================================================

test_that("fineMappingPipeline(GwasSumStats): cache hit short-circuits the RSS fitter", {
    gss <- .fmp_gss
    cachedEntry <- FineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
        susieFit = list(token = "susie_cached"),
        topLoci = data.frame(
            variant_id = sprintf("chr1:%d:A:G", 100L * (1:5)),
            pip = seq(0.9, 0.1, length.out = 5),
            stringsAsFactors = FALSE
        )
    )
    # The GwasSumStats branch keys its cache by (study, method, blockId)
    # where blockId is "{seqname}_{minPos}_{maxPos}" derived from the
    # entry's GRanges. .fmp_makeSumstatsGr() yields chr1 positions 100..500,
    # so the cache row must use blockId = "chr1_100_500" to hit.
    cache <- GwasFineMappingResult(
        studyName = "G1",
        method = "susie",
        blockId = "chr1_100_500",
        entry = list(cachedEntry),
        ldReference = .fmp_makeHandle()
    )
    rss_calls <- 0
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = function(...) {
            rss_calls <<- rss_calls + 1L
            .fmp_mockFitRss()(...)
        },
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            gss,
            methods = "susie",
            initializeWithSusieInf = FALSE,
            fineMappingResult = cache
        )
    )
    expect_equal(rss_calls, 0L)
    expect_equal(nrow(res), 1L)
})

test_that("fineMappingPipeline(GwasSumStats): wrong-shape cache (QtlFineMappingResult) is ignored", {
    # When `fineMappingResult` is a QtlFineMappingResult the GwasSumStats
    # method's cache-lookup branch should treat it as a cache miss and
    # still invoke the RSS fitter.
    gss <- .fmp_gss
    cachedEntry <- FineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
        susieFit = list(token = "susie_cached"),
        topLoci = data.frame(
            variant_id = sprintf("chr1:%d:A:G", 100L * (1:5)),
            pip = rep(0.5, 5),
            stringsAsFactors = FALSE
        )
    )
    wrongCache <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(cachedEntry)
    )
    rss_calls <- 0
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = function(...) {
            rss_calls <<- rss_calls + 1L
            .fmp_mockFitRss()(...)
        },
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            gss,
            methods = "susie",
            initializeWithSusieInf = FALSE,
            fineMappingResult = wrongCache
        )
    )
    expect_equal(rss_calls, 1L)
    expect_equal(nrow(res), 1L)
})

test_that("fineMappingPipeline(QtlDataset): cache hit avoids the fitter", {
    qd <- .fmp_qdA
    # Build a cache that already has the (study1, brain, ENSG_A, susie) row.
    cachedEntry <- FineMappingRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
        susieFit = list(token = "susie_cached"),
        topLoci = data.frame(
            variant_id = sprintf("chr1:%d:A:G", 100L * (1:3)),
            pip = c(0.9, 0.5, 0.1),
            stringsAsFactors = FALSE
        )
    )
    cache <- QtlFineMappingResult(
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "susie",
        entry = list(cachedEntry)
    )
    fitter_calls <- 0
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = function(...) {
            fitter_calls <<- fitter_calls + 1L
            .fmp_mockFitIndiv()(...)
        },
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE,
            fineMappingResult = cache
        )
    )
    expect_equal(fitter_calls, 0L)
    expect_equal(nrow(res), 1L)
})

# ===========================================================================
# Residualization settings propagation
# ===========================================================================
# The accessors take the ResidualizationParam bundle itself, so there is no
# translation to test -- only that the bundle travels intact.

test_that(".fmResidGeno / .fmResidPheno forward the bundle to the accessors", {
    capturedGeno <- NULL
    capturedPheno <- NULL
    fakeGeno <- function(x, ...) {
        capturedGeno <<- list(...)
        matrix(0, 0, 0)
    }
    fakePheno <- function(x, ...) {
        capturedPheno <<- list(...)
        list(c1 = matrix(0, 0, 0))
    }
    local_mocked_bindings(
        .qtlResidualizedGenotypes = fakeGeno,
        .qtlResidualizedPhenotypes = fakePheno,
        .package = "pecotmr"
    )
    res <- ResidualizationParam(
        phenotypeCovariates = "age",
        residualizePhenotype = FALSE,
        residualizeGenotype = TRUE
    )
    # Passed, not scraped: the settings reach the accessors because the
    # caller hands them over, with no enclosing frame involved.
    pecotmr:::.fmResidGeno(
        NULL,
        contexts = "c1",
        residualizationParam = res,
        genotypeFilter = pecotmr:::.qtlResolveFilter()
    )
    pecotmr:::.fmResidPheno(NULL, contexts = "c1", residualizationParam = res)
    for (captured in list(capturedGeno, capturedPheno)) {
        got <- captured$residualizationParam
        expect_s4_class(got, "ResidualizationParam")
        expect_equal(got$phenotypeCovariates, "age")
        expect_false(got$residualizePhenotype)
        expect_true(got$residualizeGenotype)
    }
})

test_that("an absent bundle becomes the constructor's own defaults", {
    captured <- NULL
    local_mocked_bindings(
        .qtlResidualizedPhenotypes = function(x, ...) {
            captured <<- list(...)
            list(c1 = matrix(0, 0, 0))
        },
        .package = "pecotmr"
    )
    pecotmr:::.fmResidPheno(NULL, contexts = "c1", residualizationParam = NULL)
    expect_s4_class(captured$residualizationParam, "ResidualizationParam")
    # Nothing pinned: every flag is unset, which the accessors read as "on".
    expect_length(captured$residualizationParam, 0L)
})

# ===========================================================================
# Multi-region / jointRegions (P3)
# ===========================================================================

test_that(".fmCsIdx / .fmRelabelCs parse and renumber <method>_<idx> labels", {
    expect_equal(
        pecotmr:::.fmCsIdx(c("susie_0", "susie_1", "susie_2")),
        c(0L, 1L, 2L)
    )
    expect_equal(
        pecotmr:::.fmRelabelCs(c("susie_0", "susie_1", "susie_2"), 3L),
        c("susie_0", "susie_4", "susie_5")
    )
    # offset 0 is a no-op; the "_0" not-in-CS sentinel is always preserved.
    expect_equal(
        pecotmr:::.fmRelabelCs(c("susie_0", "susie_1"), 0L),
        c("susie_0", "susie_1")
    )
})

test_that(".fmMergeEntries concatenates variants, renumbers CS, lists susieFit", {
    # Ids must encode coordinates: the merge now emits a row payload, whose
    # element is built from the range and alleles the id renders.
    mk <- function(vids, cs95, fit) {
        FineMappingRow(
            variantIds = vids,
            susieFit = fit,
            topLoci = data.frame(
                variant_id = vids,
                pip = rep(0.5, length(vids)),
                cs_95 = cs95,
                stringsAsFactors = FALSE
            )
        )
    }
    v1 <- c("chr1:100:A:G", "chr1:200:A:G")
    v2 <- c("chr1:300:A:G", "chr1:400:A:G")
    e1 <- mk(v1, c("susie_1", "susie_0"), list(tag = "f1"))
    e2 <- mk(v2, c("susie_1", "susie_1"), list(tag = "f2"))
    m <- pecotmr:::.fmMergeEntries(list(e1, e2))
    tl <- pecotmr:::.fmrPartsTopLoci(m)
    expect_equal(pecotmr:::.fmrPartsVariantIds(m), c(v1, v2))
    expect_equal(as.character(tl$variant_id), c(v1, v2))
    # e1's max CS index is 1, so e2's "susie_1" is renumbered to "susie_2".
    expect_equal(tl$cs_95, c("susie_1", "susie_0", "susie_2", "susie_2"))
    fit <- pecotmr:::.fmrPartsSusieFit(m)
    expect_true(is.list(fit))
    expect_equal(names(fit), c("region1", "region2"))
})

test_that(".fmMergeEntries: no entries -> NULL; one comes back unchanged", {
    expect_null(pecotmr:::.fmMergeEntries(list(NULL, NULL)))
    # Ids must encode coordinates now: the payload builds its element from the
    # range and alleles the id renders.
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(),
        topLoci = data.frame(variant_id = "chr1:100:A:G", pip = 0.5)
    )
    expect_identical(pecotmr:::.fmMergeEntries(list(e)), e)
})

test_that(".fmMergeEntries renumbers past a region with no cs column", {
    # A region whose fit found no set reports its top loci without a cs_95
    # column; it adds no sets, so the next region's numbering continues.
    mk <- function(vids, cs95 = NULL) {
        tl <- data.frame(variant_id = vids, pip = rep(0.5, length(vids)))
        if (!is.null(cs95)) {
            tl$cs_95 <- cs95
        }
        FineMappingRow(variantIds = vids, susieFit = NULL, topLoci = tl)
    }
    m <- pecotmr:::.fmMergeEntries(list(
        mk("chr1:100:A:G", "susie_1"),
        mk("chr1:200:A:G"),
        mk("chr1:300:A:G", "susie_1")
    ))
    tl <- pecotmr:::.fmrPartsTopLoci(m)
    expect_equal(tl$cs_95, c("susie_1", NA, "susie_2"))
})

test_that("fineMappingPipeline(QtlDataset): region is a GRanges", {
    # Given in R, a window is a GRanges; a string or a table is refused rather
    # than read, so its base is never in doubt.
    expect_error(
        fineMappingPipeline(
            .fmp_qdAB,
            methods = "susie",
            region = "chr1:1-600"
        ),
        "`region` must be a GRanges (got character)",
        fixed = TRUE
    )
})

test_that("fineMappingPipeline(QtlDataset): region + cisWindow is rejected", {
    qd <- .fmp_qdAB
    expect_error(
        fineMappingPipeline(
            qd,
            methods = "susie",
            region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 600)),
            cisWindow = 1000L
        ),
        "either `region` or `cisWindow`"
    )
})

test_that("fineMappingPipeline(QtlDataset): jointRegions=FALSE merges per-region fits", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    regions <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(50L, 350L), c(350L, 650L))
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "susie",
        traitId = "ENSG_A",
        region = regions,
        jointRegions = FALSE,
        initializeWithSusieInf = FALSE
    ))
    expect_s4_class(res, "QtlFineMappingResult")
    # 1 ctx x 1 trait x 1 method -> a single merged row, not one per region.
    expect_equal(nrow(res), 1L)
    fit <- susieFit(
        res,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "susie"
    )
    expect_equal(names(fit), c("region1", "region2")) # per-region fit list
})

test_that("fineMappingPipeline(QtlDataset): jointRegions=TRUE fits one concatenated block", {
    qd <- .fmp_qdA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    regions <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(50L, 350L), c(350L, 650L))
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "susie",
        traitId = "ENSG_A",
        region = regions,
        jointRegions = TRUE,
        initializeWithSusieInf = FALSE
    ))
    expect_equal(nrow(res), 1L)
    fit <- susieFit(
        res,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "susie"
    )
    # one concatenated fit -> the mock fit object, not a per-region list.
    expect_false(identical(names(fit), c("region1", "region2")))
})

test_that("fineMappingPipeline(QtlDataset): mvsusie jointRegions=FALSE merges per-region fits", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    regions <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(50L, 350L), c(350L, 650L))
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "mvsusie",
        traitId = c("ENSG_A", "ENSG_B"),
        region = regions,
        jointRegions = FALSE
    ))
    expect_equal(nrow(res), 2L) # joint fit fanned out to both traits
    fit <- susieFit(
        res,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "mvsusie"
    )
    expect_equal(names(fit), c("region1", "region2")) # per-region merged fit list
})

test_that("fineMappingPipeline(QtlDataset): fsusie jointRegions=FALSE merges per-region fits", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(susiF = .fmp_mockSusiF(), .package = "fsusieR")
    regions <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(50L, 350L), c(350L, 650L))
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "fsusie",
        traitId = c("ENSG_A", "ENSG_B"),
        region = regions,
        jointRegions = FALSE
    ))
    # One row for the joint fsusie fit (see above); the per-REGION merge is
    # what `names(fit)` checks, and it is unaffected.
    expect_equal(nrow(res), 1L)
    fit <- susieFit(
        res,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "fsusie"
    )
    expect_equal(names(fit), c("region1", "region2"))
})

test_that("fineMappingPipeline(QtlDataset): jointSpec + jointRegions=FALSE merges per region", {
    # One trait: the cross-context joint and the per-region merge are what is
    # tested, and each further trait only repeats them.
    qd <- .fmp_qdTwoCtxA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    local_mocked_bindings(
        mvsusie = .fmp_mockMvsusie(),
        create_mixture_prior = .fmp_mockMixturePrior(),
        .package = "mvsusieR"
    )
    regions <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(50L, 350L), c(350L, 650L))
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "mvsusie",
        traitId = "ENSG_A",
        region = regions,
        jointRegions = FALSE,
        jointSpecification = "context"
    ))
    expect_s4_class(res, "QtlFineMappingResult")
    # cross-context joint -> one row per context, each merged across the 2
    # regions.
    expect_equal(nrow(res), 2L)
    expect_setequal(as.character(res$context), c("brain", "liver"))
    fit <- susieFit(
        res,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "mvsusie"
    )
    expect_equal(names(fit), c("region1", "region2")) # merged across regions
})

# ===========================================================================
# Additional coverage: tractable pure helpers. (The real-fit internals --
# .fmPostprocessOne, mvsusie/fsusie fold weights, susieInf chaining, the
# top-PC cross-trait path -- are driven by the mocked-fit pipeline tests and
# the actual fitting is left to the external solvers.)
# ===========================================================================

test_that(".fmCheckMethodCapabilities: empty token list is a no-op", {
    expect_null(pecotmr:::.fmCheckMethodCapabilities(
        character(0),
        "QtlDataset"
    ))
})

test_that(".fmCacheLookupGwas: NULL / non-GwasFineMappingResult -> NULL", {
    expect_null(pecotmr:::.fmCacheLookupGwas(NULL, "G1", "susie", "chr1:1-100"))
    fmr <- QtlFineMappingResult(
        studyName = "S",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(
            "chr1:100:A:G",
            list(),
            data.frame(variant_id = "chr1:100:A:G", pip = 0.9)
        ))
    )
    expect_null(pecotmr:::.fmCacheLookupGwas(fmr, "S", "susie", "chr1:1-100"))
})

test_that(".fmBuildMvsusiePriorCv: NULL CV -> NULL; NULL fold fits are skipped", {
    expect_null(pecotmr:::.fmBuildMvsusiePriorCv(NULL, NULL, c("c1", "c2")))
    local_mocked_bindings(
        .buildMvsusieReweightedPrior = function(...) "PRIOR",
        .package = "pecotmr"
    )
    mvCv <- list(
        samplePartition = data.frame(
            Sample = paste0("s", 1:4),
            Fold = c(1, 1, 2, 2)
        ),
        foldFits = list(fold_1 = list(w0 = 1), fold_2 = NULL)
    )
    out <- pecotmr:::.fmBuildMvsusiePriorCv(
        mvCv,
        list(w0 = 1, V = diag(2)),
        c("c1", "c2")
    )
    expect_equal(out[["1"]], "PRIOR")
    expect_null(out[["2"]]) # a NULL fold fit is skipped
})

test_that(".fmFoldFitAt finds a fold's fit by name, else by position", {
    fn <- pecotmr:::.fmFoldFitAt
    named <- list(fold_2 = "fit2", fold_1 = "fit1")
    expect_equal(fn(1L, c(1, 2), named), "fit1")
    # Fits not named fold_<id> are taken in fold order, and a fold beyond
    # them has none.
    unnamed <- list("first", "second")
    expect_equal(fn(2L, c(1, 2), unnamed), "second")
    expect_null(fn(3L, c(1, 2, 3), unnamed))
})

test_that(".fmTwasMethodKey: bare token without adapter returned unchanged", {
    expect_equal(pecotmr:::.fmTwasMethodKey("lasso"), "lasso") # no adapter
    expect_equal(pecotmr:::.fmTwasMethodKey("susie"), "susie") # adapter -> stripped
})

test_that(".fmSliceCv: NULL cv or missing predicted key -> NULL", {
    expect_null(pecotmr:::.fmSliceCv(NULL, "susie"))
    cv <- list(
        samplePartition = NULL,
        prediction = list(enet_predicted = matrix(0, 1, 1)),
        performance = list(enet_performance = matrix(0, 1, 6))
    )
    expect_null(pecotmr:::.fmSliceCv(cv, "susie")) # no susie_predicted key
    expect_true(
        "enet_predicted" %in% names(pecotmr:::.fmSliceCv(cv, "enet")$prediction)
    )
})

test_that(".fmAttachCv: NULL entry or NULL cvResult returns the entry unchanged", {
    e <- FineMappingRow(
        "chr1:100:A:G",
        list(),
        data.frame(variant_id = "chr1:100:A:G", pip = 0.9)
    )
    expect_identical(pecotmr:::.fmAttachCv(e, NULL), e)
    expect_null(pecotmr:::.fmAttachCv(NULL, list(x = 1)))
    expect_equal(
        pecotmr:::.fmrPartsCvResult(
            pecotmr:::.fmAttachCv(e, list(samplePartition = 1))
        ),
        list(samplePartition = 1)
    )
})

# ---- jointSpec dispatch branches in the QtlSumStats / MultiStudy methods -----
# (dispatchers mocked; the real joint fitting is covered in
# test_jointSpecification.R)

test_that("fineMappingPipeline(QtlSumStats): mvsusie-only jointSpec returns the joint result", {
    ss <- .fmp_qss
    jr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "mvsusie",
        entry = list(FineMappingRow(
            "chr1:100:A:G",
            list(),
            data.frame(variant_id = "chr1:100:A:G", pip = 0.9)
        ))
    )
    local_mocked_bindings(
        .fmDispatchJointSpecsQtlSumStats = function(...) jr,
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            ss,
            methods = "mvsusie",
            jointSpecification = "context"
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_setequal(as.character(res$method), "mvsusie")
})

test_that("fineMappingPipeline(QtlSumStats): joint spec plus per-entry susie", {
    # The joint spec consumes mvsusie; susie, asked for alongside it, still
    # runs per entry, and the two results are combined.
    ss <- .fmp_qss
    jr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "mvsusie",
        entry = list(FineMappingRow(
            "chr1:100:A:G",
            list(),
            data.frame(variant_id = "chr1:100:A:G", pip = 0.9)
        ))
    )
    local_mocked_bindings(
        .fmDispatchJointSpecsQtlSumStats = function(...) jr,
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        ss,
        methods = c("mvsusie", "susie"),
        jointSpecification = "context"
    ))
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(sum(res$method == "mvsusie"), 1L)
    expect_equal(sum(res$method == "susie"), nrow(ss))
})

test_that("fineMappingPipeline(QtlSumStats): no auto-joint fit is an error", {
    # Without a joint spec mvsusie is fit by the automatic cross-context
    # joint; with no other method to run, a joint that produces nothing
    # leaves no result at all.
    ss <- .fmp_makeMultiCtxQtlSumStats()
    local_mocked_bindings(
        .fmDispatchJointSpecsQtlSumStats = function(...) NULL,
        .package = "pecotmr"
    )
    expect_error(
        suppressMessages(fineMappingPipeline(ss, methods = "mvsusie")),
        "no entries produced a result"
    )
})

test_that("fineMappingPipeline(QtlSumStats): mvsusie-only jointSpec with no fits -> error", {
    local_mocked_bindings(
        .fmDispatchJointSpecsQtlSumStats = function(...) NULL,
        .package = "pecotmr"
    )
    expect_error(
        suppressMessages(fineMappingPipeline(
            .fmp_qss,
            methods = "mvsusie",
            jointSpecification = "context"
        )),
        "no joint fits produced"
    )
})

test_that("MultiStudyQtlDataset forwards the mvsusie prior settings", {
    # `mrmashPrior` and `dataDrivenPriorWeightsCutoff` configure the mvsusie
    # data-driven reweighted prior. They used to reach ONLY jointCfg (the
    # cross-study joint phase); the per-study config omitted them, so every
    # per-study recursion silently refit mr.mash and used the 1e-10 cutoff.
    mt <- .fmp_msqd
    seen <- new.env(parent = emptyenv())
    seen$perStudy <- NULL
    seen$sumStats <- NULL
    local_mocked_bindings(
        .fmPerStudy = function(qd, cfg) {
            seen$perStudy <- cfg
            NULL
        },
        .fmSumStats = function(ss, cfg) {
            seen$sumStats <- cfg
            NULL
        },
        .package = "pecotmr"
    )
    # Both workers are stubbed, so the driver rbinds nothing -- the configs
    # they were handed are what this test is about.
    try(
        suppressMessages(suppressWarnings(fineMappingPipeline(
            mt,
            methods = "susie",
            verbose = 0,
            mrmashPrior = "cache-sentinel",
            dataDrivenPriorWeightsCutoff = 0.5
        ))),
        silent = TRUE
    )
    for (cfg in list(seen$perStudy, seen$sumStats)) {
        expect_equal(cfg$mrmashPrior, "cache-sentinel")
        expect_equal(cfg$dataDrivenPriorWeightsCutoff, 0.5)
    }
})

test_that("fineMappingPipeline(MSQD) lays cisWindow on its sumstats", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        .fmp_msqd,
        methods = "susie",
        cisWindow = 100L,
        initializeWithSusieInf = FALSE
    ))
    fromSumStats <- which(res$study == "Q1")
    expect_identical(
        as.character(genomicRegion(res)[fromSumStats]),
        "chr1:200-400"
    )
    # A summary-statistics study without a trait position is an error, as it
    # is on its own.
    bare <- MultiStudyQtlDataset(
        qtlDatasets = list(study1 = .fmp_qdA),
        sumStats = .fmp_qss
    )
    expect_error(
        fineMappingPipeline(bare, methods = "susie", cisWindow = 100L),
        "fineMappingPipeline: `cisWindow` is laid around each tuple's trait",
        fixed = TRUE
    )
})

test_that("fineMappingPipeline(MultiStudyQtlDataset): region + cisWindow is rejected", {
    expect_error(
        fineMappingPipeline(
            .fmp_msqd,
            methods = "mvsusie",
            region = GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100)),
            cisWindow = 1000L
        ),
        "specify either"
    )
})

test_that("fineMappingPipeline(MultiStudyQtlDataset): mvsusie-only jointSpec returns the joint result", {
    mt <- .fmp_msqd
    jr <- QtlFineMappingResult(
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "mvsusie",
        entry = list(FineMappingRow(
            "chr1:100:A:G",
            list(),
            data.frame(variant_id = "chr1:100:A:G", pip = 0.9)
        ))
    )
    local_mocked_bindings(
        .fmDispatchJointSpecsMultiStudy = function(...) jr,
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            mt,
            methods = "mvsusie",
            jointSpecification = "context"
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_setequal(as.character(res$method), "mvsusie")
})

test_that("fineMappingPipeline(MultiStudy): joint spec, rest per component", {
    # The joint spec consumes mvsusie; susie, asked for alongside it, still
    # runs on every component -- the study and the embedded sumstats.
    mt <- .fmp_msqd
    jr <- QtlFineMappingResult(
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "mvsusie",
        entry = list(FineMappingRow(
            "chr1:100:A:G",
            list(),
            data.frame(variant_id = "chr1:100:A:G", pip = 0.9)
        ))
    )
    local_mocked_bindings(
        .fmDispatchJointSpecsMultiStudy = function(...) jr,
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        mt,
        methods = c("mvsusie", "susie"),
        jointSpecification = "context",
        cisWindow = 1000L,
        initializeWithSusieInf = FALSE
    ))
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(sum(res$method == "mvsusie"), 1L)
    expect_gt(sum(res$method == "susie"), 1L)
})

test_that("fineMappingPipeline(MultiStudy): an unselected component sits out", {
    # One context vector applies to every component; one lacking it keeps
    # none of its tuples and contributes no rows.
    mt <- .fmp_msqd
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    fit <- function(contexts) {
        suppressMessages(fineMappingPipeline(
            mt,
            methods = "susie",
            contexts = contexts,
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        ))
    }
    expect_warning(
        sumstatsOnly <- fit("c1"),
        "study 'study1' is missing requested context\\(s\\): c1"
    )
    expect_setequal(as.character(sumstatsOnly$study), "Q1")
    expect_warning(
        studyOnly <- fit("brain"),
        "study 'Q1' is missing requested context\\(s\\): brain"
    )
    expect_setequal(as.character(studyOnly$study), "study1")
})

test_that("fineMappingPipeline(MultiStudy): trait selection skips a study", {
    # A trait vector applies to every component too: one carrying none of
    # the traits sits out, and the run is refused only when none carries any.
    mt <- .fmp_msqd
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    fit <- function(traitId) {
        suppressMessages(fineMappingPipeline(
            mt,
            methods = "susie",
            traitId = traitId,
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        ))
    }
    expect_setequal(as.character(fit("t1")$study), "Q1")
    expect_setequal(as.character(fit("ENSG_A")$study), "study1")
    expect_error(
        fit("absent"),
        "fineMappingPipeline(MultiStudyQtlDataset): no traits selected.",
        fixed = TRUE
    )
})

test_that("fineMappingPipeline(MultiStudy): joint rows record scale override", {
    # The study's own slot scales its residuals; the per-call override says
    # otherwise and wins, and the joint rows record what was used.
    mt <- .fmp_msqd
    jr <- QtlFineMappingResult(
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "mvsusie",
        entry = list(FineMappingRow(
            "chr1:100:A:G",
            list(),
            data.frame(variant_id = "chr1:100:A:G", pip = 0.9)
        ))
    )
    local_mocked_bindings(
        .fmDispatchJointSpecsMultiStudy = function(...) jr,
        .package = "pecotmr"
    )
    fit <- function(residualizationParam) {
        suppressMessages(fineMappingPipeline(
            mt,
            methods = "mvsusie",
            jointSpecification = "context",
            residualizationParam = residualizationParam
        ))
    }
    expect_true(fit(ResidualizationParam())$standardized)
    expect_false(
        fit(ResidualizationParam(scaleResiduals = FALSE))$standardized
    )
})

# =============================================================================
# Coverage top-ups (mock the fitters, exercise the orchestration)
# =============================================================================

test_that(".fmRelabelCs returns non-matching labels unchanged", {
    fn <- pecotmr:::.fmRelabelCs
    # "nomatch" doesn't match ^(.*)_([0-9]+)$ (length(parts) != 3) -> returned
    # as-is; "susie_0" keeps the not-in-CS sentinel; "susie_1" shifts by offset.
    out <- fn(c("susie_1", "nomatch", "susie_0"), offset = 2L)
    expect_equal(out, c("susie_3", "nomatch", "susie_0"))
})

test_that(".fmWeightsCv + .fmFoldWeights cover the mvSuSiE CV path (mocked fitter)", {
    # Mock one level below the orchestration: fitMvsusie/mvsusieWeights return
    # canned outputs (sized to the per-fold training columns), so the real
    # .fmFoldWeights mvsusie branch + .fmWeightsCv fold loop run at ~no cost.
    local_mocked_bindings(
        fitMvsusie = function(X, Y, ...) {
            list(vn = colnames(X), R = ncol(as.matrix(Y)))
        },
        mvsusieWeights = function(mvsusieFit = NULL, ...) {
            matrix(
                0.01,
                length(mvsusieFit$vn),
                mvsusieFit$R,
                dimnames = list(mvsusieFit$vn, NULL)
            )
        },
        .package = "pecotmr"
    )
    set.seed(1)
    n <- 30L
    p <- 5L
    R <- 2L
    X <- matrix(
        rbinom(n * p, 2, 0.4),
        n,
        p,
        dimnames = list(
            paste0("s", seq_len(n)),
            sprintf("chr1:%d:A:G", 100L * seq_len(p))
        )
    )
    Y <- matrix(rnorm(n * R), n, R, dimnames = list(rownames(X), c("c1", "c2")))
    cv <- pecotmr:::.fmWeightsCv(
        X,
        Y,
        tokens = "mvsusie",
        methodArgs = list(mvsusie = list()),
        fold = 3L,
        coverage = 0.95,
        verbose = 0
    )
    expect_named(
        cv,
        c("samplePartition", "prediction", "performance", "foldFits")
    )
    expect_true("mvsusie_performance" %in% names(cv$performance))
    expect_equal(dim(cv$prediction[["mvsusie_predicted"]]), c(n, R))
})

test_that(".fmFoldWeights covers the fSuSiE branch (mocked fitter)", {
    local_mocked_bindings(
        fitFsusie = function(...) list(),
        fsusieWeights = function(fsusieFit = NULL, variantIds = NULL, ...) {
            matrix(
                0.02,
                length(variantIds),
                1L,
                dimnames = list(variantIds, NULL)
            )
        },
        .package = "pecotmr"
    )
    set.seed(2)
    n <- 30L
    p <- 5L
    X <- matrix(
        rbinom(n * p, 2, 0.4),
        n,
        p,
        dimnames = list(
            paste0("s", seq_len(n)),
            sprintf("chr1:%d:A:G", 100L * seq_len(p))
        )
    )
    Y <- matrix(rnorm(n * 4L), n, 4L, dimnames = list(rownames(X), NULL))
    W <- pecotmr:::.fmFoldWeights(
        "fsusie",
        X,
        Y,
        coverage = 0.95,
        userArgs = list(),
        pos = seq_len(p)
    )
    expect_true(is.matrix(W))
    expect_equal(rownames(W), colnames(X))
})

test_that(".fmPostprocessOne wraps a fit into a FineMappingRow", {
    local_mocked_bindings(
        postprocessFinemappingFits = function(...) list(x = 1),
        formatFinemappingOutput = function(post, primaryMethod, ...) {
            list(
                finemappingEntry = FineMappingRow(
                    "chr1:100:A:G",
                    list(),
                    data.frame(variant_id = "chr1:100:A:G", pip = 0.5)
                )
            )
        },
        .package = "pecotmr"
    )
    ent <- pecotmr:::.fmPostprocessOne(
        fit = list(),
        method = "susie",
        dataX = matrix(0, 2, 1, dimnames = list(NULL, "chr1:100:A:G")),
        dataY = c(1, 2),
        credibleSetParam = CredibleSetParam(
            coverage = 0.95,
            secondaryCoverage = 0.7,
            signalCutoff = 0.1,
            minAbsCorr = 0.5
        )
    )
    expect_s4_class(ent, "FineMappingRow")
})

test_that(".fmPostprocessOne errors when output carries no FineMappingRow", {
    local_mocked_bindings(
        postprocessFinemappingFits = function(...) list(),
        formatFinemappingOutput = function(...) list(finemappingEntry = "nope"),
        .package = "pecotmr"
    )
    expect_error(
        pecotmr:::.fmPostprocessOne(
            list(),
            "susie",
            matrix(0, 2, 1),
            c(1, 2),
            0.95,
            0.7,
            0.1,
            0.5
        ),
        "fine-mapping row"
    )
})

test_that(".fmWeightsCv: fold prior else mvPrior; weightless folds give NA", {
    # Three folds of ten. Fold 1 has its own prior in mvPriorCv and weights
    # (all 1) on the test variants; fold 2 has its own prior but its fit
    # returns no weights; fold 3 has no prior of its own, so it takes the
    # full-data mvPrior, and its weights share no variant with the test
    # genotypes. Each fold fit records the prior it was handed.
    set.seed(3)
    X <- matrix(
        rbinom(60, 2, 0.4),
        30,
        2,
        dimnames = list(str_c("s", 1:30), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    Y <- matrix(rnorm(60), 30, 2, dimnames = list(rownames(X), c("c1", "c2")))
    sp <- data.frame(Sample = rownames(X), Fold = rep(1:3, each = 10))
    fullPrior <- list(priorVariance = diag(2) * 3)
    foldPriors <- list(
        "1" = list(priorVariance = diag(2)),
        "2" = list(priorVariance = diag(2) * 2)
    )
    local_mocked_bindings(
        .fmFoldWeights = function(
            token,
            Xtr,
            Ytr,
            coverage,
            userArgs,
            pos,
            mvPrior
        ) {
            trained <- sp$Fold[is_in(sp$Sample, rownames(Xtr))]
            heldOut <- setdiff(sp$Fold, trained)
            if (heldOut == 2L) {
                return(NULL)
            }
            vars <- if (heldOut == 1L) colnames(Xtr) else "not_a_variant"
            W <- matrix(1, length(vars), ncol(Ytr), dimnames = list(vars, NULL))
            `attr<-`(W, "fit", list(prior = mvPrior))
        },
        .package = "pecotmr"
    )
    cv <- pecotmr:::.fmWeightsCv(
        X,
        Y,
        tokens = "mvsusie",
        methodArgs = list(mvsusie = list()),
        fold = 3L,
        samplePartition = sp,
        verbose = 0,
        mvPrior = fullPrior,
        mvPriorCv = foldPriors
    )
    # A fold takes its own prior when mvPriorCv has one, else mvPrior.
    expect_identical(cv$foldFits$fold_1$mvsusie$prior, foldPriors[["1"]])
    expect_identical(cv$foldFits$fold_3$mvsusie$prior, fullPrior)
    # A fold whose fit returns no weights keeps no fit either.
    expect_length(cv$foldFits$fold_2, 0L)
    # Only fold 1 predicts: X %*% W, the dosage row sums. NULL weights (fold
    # 2) and weights on no test variant (fold 3) leave their samples NA ...
    pred <- cv$prediction[["mvsusie_predicted"]]
    f1 <- sp$Sample[sp$Fold == 1L]
    expect_equal(
        pred[f1, ],
        cbind(c1 = rowSums(X[f1, ]), c2 = rowSums(X[f1, ]))
    )
    expect_true(all(is.na(pred[!is_in(rownames(pred), f1), ])))
    # ... and those samples are left out of the score, not counted.
    expect_equal(
        cv$performance[["mvsusie_performance"]][, "corr"],
        drop(cor(Y[f1, ], rowSums(X[f1, ])))
    )
})

test_that("fineMappingPipeline(*SumStats): susieInf RSS chain (mocked)", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    for (input in list(.fmp_qss, .fmp_gss)) {
        kind <- class(input)[[1L]]
        res <- suppressMessages(
            fineMappingPipeline(
                input,
                methods = "susieInf",
                initializeWithSusieInf = FALSE,
                verbose = 1
            )
        )
        expect_s4_class(res, .fmp_resultClass[[kind]])
        expect_equal(methodNames(res), "susieInf", info = kind)
    }
})

test_that(".fmFoldWeights covers mvPrior residual var, missing rownames, unknown token", {
    local_mocked_bindings(
        fitMvsusie = function(X, Y, ...) {
            list(vn = colnames(X), R = ncol(as.matrix(Y)))
        },
        # weights WITHOUT rownames -> .fmFoldWeights sets them from Xtr
        mvsusieWeights = function(mvsusieFit = NULL, ...) {
            matrix(0.01, length(mvsusieFit$vn), mvsusieFit$R)
        },
        .package = "pecotmr"
    )
    X <- matrix(
        rbinom(40, 2, 0.4),
        20,
        2,
        dimnames = list(paste0("s", 1:20), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    Y <- matrix(rnorm(40), 20, 2)
    # an mvPrior carrying a residualVariance hands it to the fit
    W <- pecotmr:::.fmFoldWeights(
        "mvsusie",
        X,
        Y,
        coverage = 0.95,
        userArgs = list(),
        pos = NULL,
        mvPrior = list(priorVariance = diag(2), residualVariance = diag(2))
    )
    expect_equal(rownames(W), colnames(X))
    # an unknown token falls through to NULL
    expect_null(pecotmr:::.fmFoldWeights("bogus", X, Y, 0.95, list(), NULL))
})

test_that(".fmWeightsCv returns NULL for empty tokens", {
    X <- matrix(
        0,
        10,
        2,
        dimnames = list(paste0("s", 1:10), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    expect_null(pecotmr:::.fmWeightsCv(
        X,
        matrix(0, 10, 1),
        tokens = character(0),
        methodArgs = list(),
        fold = 2L
    ))
})

test_that(".fmWeightsCv fills Y rownames and reports per-fold fit failures", {
    local_mocked_bindings(
        .fmFoldWeights = function(...) stop("boom"),
        .package = "pecotmr"
    )
    X <- matrix(
        rbinom(40, 2, 0.4),
        20,
        2,
        dimnames = list(paste0("s", 1:20), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    Y <- matrix(rnorm(20), 20, 1) # no rownames -> filled from X
    expect_message(
        cv <- suppressWarnings(pecotmr:::.fmWeightsCv(
            X,
            Y,
            tokens = "susie",
            methodArgs = list(susie = list()),
            fold = 2L,
            verbose = 1
        )),
        "CV fold .* failed"
    )
})

test_that("fineMappingPipeline CV retains a lean per-fold fit", {
    f <- .twcvFixture
    expect_false(is.null(f$cv$foldFits))
    expect_equal(names(f$cv$foldFits), c("fold_1", "fold_2", "fold_3"))
    # lean: only the fields the weight extractors read
    fit1 <- f$cv$foldFits[["fold_1"]][["susie"]]
    expect_true(all(c("pip", "alpha", "mu") %in% names(fit1)))
    expect_false("lbf_variable" %in% names(fit1))
    # and it slices per method onto the row payload
    expect_equal(names(f$slice$foldFits), c("fold_1", "fold_2", "fold_3"))
})

# --- method-level branches (drive the pipeline methods with mocked fitters) ---

test_that(".fmAfForX returns NULL when af yields nothing", {
    qd <- .fmp_qdA
    local_mocked_bindings(
        .qtlExtractBlock = function(...) {
            list(af = numeric(0), variantIds = character(0))
        },
        .package = "pecotmr"
    )
    X <- matrix(
        0,
        3,
        2,
        dimnames = list(paste0("s", 1:3), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    expect_null(pecotmr:::.fmAfForX(
        qd,
        X,
        genotypeFilter = pecotmr:::.qtlResolveFilter(),
        traitId = "ENSG_A"
    ))
})

test_that("fineMappingPipeline(QtlDataset): a valid context is honored", {
    qd <- .fmp_qdTwoCtxA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "susie",
            contexts = "brain",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    # The context left out is not fit.
    expect_setequal(contexts(res), "brain")
})

test_that("fineMappingPipeline(QtlDataset): region selects traits by rowRanges overlap", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    # ENSG_A spans 1000-1499 and ENSG_B 2000-2499: the region reaches only A.
    region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 1600L))
    res <- suppressMessages(
        fineMappingPipeline(
            .fmp_qdAB,
            methods = "susie",
            region = region,
            initializeWithSusieInf = FALSE
        )
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_setequal(traitNames(res), "ENSG_A")
})

test_that("fineMappingPipeline(QtlDataset): too few shared samples errors", {
    qd <- .fmp_qdA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        # residualized X carries a sample absent from Y -> < 2 shared samples
        .fmResidGeno = function(x, ...) {
            matrix(
                0,
                1L,
                2L,
                dimnames = list(
                    "ghost_sample",
                    c("chr1:100:A:G", "chr1:200:A:G")
                )
            )
        },
        .package = "pecotmr"
    )
    expect_error(
        suppressMessages(fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        )),
        "too few shared samples"
    )
})

test_that("fineMappingPipeline(QtlDataset): errors when no tuple produces a result", {
    qd <- .fmp_qdA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmSerScreen = function(...) FALSE, # screen out every block
        .package = "pecotmr"
    )
    expect_error(
        suppressMessages(fineMappingPipeline(
            qd,
            methods = "susie",
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        )),
        "no .*tuples"
    )
})

test_that("fineMappingPipeline(QtlDataset): usePCA skips single-trait contexts", {
    qd <- .fmp_qdA
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .fmSerScreen = function(...) TRUE,
        .package = "pecotmr"
    )
    res <- suppressMessages(
        fineMappingPipeline(
            qd,
            methods = "susie",
            usePCA = TRUE,
            nPCs = 1L,
            cisWindow = 1000L,
            initializeWithSusieInf = FALSE
        )
    )
    # single-trait context -> PC loop hits `length(traits) < 2L` next; only the
    # univariate susie row survives, no topPC pseudo-trait.
    expect_false(any(grepl("PC", as.character(res$trait))))
})

test_that("fineMappingPipeline(MultiStudyQtlDataset): jointSpec with no intersecting scope errors", {
    mt <- .fmp_msqd
    # The joint engine yields nothing for this scope -> the no-joint-fits stop.
    local_mocked_bindings(
        .fmDispatchJointSpecsMultiStudy = function(...) NULL,
        .package = "pecotmr"
    )
    expect_error(
        suppressMessages(fineMappingPipeline(
            mt,
            methods = "mvsusie",
            jointSpecification = "context"
        )),
        "no joint fits produced"
    )
})

# --- usePCA sub-branches. The tests using methods = "mvsusie" skip the
# univariate dispatch (and its SER screen); the PCA path still runs susie, and
# the mocked mvsusie joint dispatch keeps the result non-empty.
.fmp_jr <- function() {
    QtlFineMappingResult(
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "mvsusie",
        entry = list(FineMappingRow(
            "chr1:100:A:G",
            list(),
            data.frame(variant_id = "chr1:100:A:G", pip = 0.9)
        ))
    )
}

test_that("fineMappingPipeline(QtlDataset): usePCA skips a context whose PCA yields no scores", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .fmTopPcScores = function(...) NULL, # PCA yields no scores
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "susie",
        usePCA = TRUE,
        nPCs = 1L,
        cisWindow = 1000L,
        initializeWithSusieInf = FALSE
    ))
    expect_false(any(grepl("PC", as.character(res$trait))))
})

test_that("fineMappingPipeline(QtlDataset): usePCA reuses a cached PC entry", {
    qd <- .fmp_qdAB
    cachedFMR <- QtlFineMappingResult(
        studyName = "study1",
        context = "brain",
        trait = "topPC1",
        method = "susie",
        entry = list(FineMappingRow(
            "chr1:100:A:G",
            list(),
            data.frame(variant_id = "chr1:100:A:G", pip = 0.9)
        ))
    )
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .fmTopPcScores = function(Y, nPCs) {
            matrix(
                rnorm(nrow(Y)),
                nrow(Y),
                1L,
                dimnames = list(rownames(Y), "topPC1")
            )
        },
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "susie",
        usePCA = TRUE,
        nPCs = 1L,
        cisWindow = 1000L,
        initializeWithSusieInf = FALSE,
        fineMappingResult = cachedFMR
    ))
    expect_true(any(as.character(res$trait) == "topPC1")) # from the cache
})

test_that("fineMappingPipeline(QtlDataset): usePCA + region uses the region genotype block", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieIndiv = .fmp_mockFitIndiv(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .fmSerScreen = function(...) TRUE,
        .fmDispatchJointSpecsQtlDataset = function(...) NULL,
        .package = "pecotmr"
    )
    region <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 3000L))
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "mvsusie",
        usePCA = TRUE,
        nPCs = 1L,
        region = region
    ))
    expect_true(any(grepl("PC", as.character(res$trait))))
})

test_that("fineMappingPipeline(QtlDataset): usePCA skips a PC block with too few shared samples", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmResidGeno = function(x, ...) {
            matrix(
                0,
                1L,
                2L,
                dimnames = list(
                    "ghost_sample",
                    c("chr1:100:A:G", "chr1:200:A:G")
                )
            )
        }, # < 2 samples shared with the PCs
        .fmDispatchJointSpecsQtlDataset = function(...) .fmp_jr(),
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "mvsusie",
        usePCA = TRUE,
        nPCs = 1L,
        cisWindow = 1000L
    ))
    expect_false(any(grepl("PC", as.character(res$trait))))
})

test_that("fineMappingPipeline(QtlDataset): usePCA skips a PC block screened out by SER", {
    qd <- .fmp_qdAB
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmSerScreen = function(...) FALSE, # the PC block is screened out
        .fmDispatchJointSpecsQtlDataset = function(...) .fmp_jr(),
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "mvsusie",
        usePCA = TRUE,
        nPCs = 1L,
        cisWindow = 1000L
    ))
    expect_false(any(grepl("PC", as.character(res$trait))))
})

# ---------------------------------------------------------------------------
# RSS panel filters: mafCutoff / macCutoff / imissCutoff on QtlSumStats and
# GwasSumStats input.
#
# There is no genotype matrix on the RSS path, so the cutoffs are measured
# against the LD reference panel. Filtering happens in .fmExtractZn, before z
# and the LD matrix are built, so the two cannot fall out of alignment.
# ---------------------------------------------------------------------------

# The bundled QTL summary statistics, QC'd once: summaryStatsQc() is the slow
# step, and every test below only reads the result.
.rssf_ss <- local({
    data(qtlSumStatsExample, envir = environment())
    suppressMessages(summaryStatsQc(qtlSumStatsExample))
})

# @noRd
.rssf_n <- function(ss, ...) {
    # `...` is the panel filter: every caller passes only its fields.
    sum(lengths(suppressWarnings(suppressMessages(
        fineMappingPipeline(
            ss,
            methods = "susie",
            panelFilterParam = PanelFilterParam(...)
        )
    ))))
}

test_that("fineMappingPipeline(QtlSumStats) filters nothing by default", {
    ss <- .rssf_ss
    expect_equal(.rssf_n(ss), sum(lengths(ss)))
})

test_that("fineMappingPipeline RSS MAF cutoffs match .panelVariantFilter", {
    # The pipeline must drop exactly the variants the shared filter names --
    # no more (which would mean something else is filtering) and no fewer
    # (which would mean the cutoff silently did nothing) -- so panel-rare
    # variants go, and a tighter cutoff drops more.
    ss <- .rssf_ss
    ids <- variantIds(ss)
    kept <- sum(lengths(ss))
    for (cut in c(0.05, 0.2)) {
        expected <- length(.panelVariantFilter(
            ldReference(ss)[[1L]],
            ids,
            PanelFilterParam(mafCutoff = cut)
        ))
        n <- .rssf_n(ss, mafCutoff = cut)
        expect_equal(n, expected, label = str_c("mafCutoff ", cut))
        expect_lt(n, kept, label = str_c("mafCutoff ", cut))
        kept <- n
    }
})

test_that("fineMappingPipeline RSS treats MAC as a MAF equivalent", {
    ss <- .rssf_ss
    nSamp <- ncol(ldReference(ss)[[1L]])
    expect_equal(
        .rssf_n(ss, macCutoff = 0.1 * 2 * nSamp),
        .rssf_n(ss, mafCutoff = 0.1)
    )
})

test_that("fineMappingPipeline RSS honours a missingness cutoff", {
    ss <- .rssf_ss
    expect_lt(.rssf_n(ss, imissCutoff = 0), sum(lengths(ss)))
    expect_equal(.rssf_n(ss, imissCutoff = 1), sum(lengths(ss)))
})

test_that(".fmExtractZn keeps z, ids and n aligned after filtering", {
    ss <- .rssf_ss
    entry <- pecotmr:::.collectionEntry(ss, 1L)
    sketch <- ldReference(ss)[[1L]]
    unfiltered <- pecotmr:::.fmExtractZn(entry, "t")
    filtered <- suppressMessages(pecotmr:::.fmExtractZn(
        entry,
        "t",
        ldReference = sketch,
        cutoffs = list(mafCutoff = 0.2, macCutoff = 0, imissCutoff = 1)
    ))
    expect_lt(length(filtered$variantIds), length(unfiltered$variantIds))
    expect_equal(length(filtered$z), length(filtered$variantIds))
    # The surviving z values are the originals, not a reindexed subset.
    keep <- is_in(unfiltered$variantIds, filtered$variantIds)
    expect_equal(filtered$z, unfiltered$z[keep])
    expect_equal(filtered$variantIds, unfiltered$variantIds[keep])
})

test_that(".fmLeanFoldFit passes a NULL fit straight through", {
    expect_null(pecotmr:::.fmLeanFoldFit(NULL, "susie"))
})

test_that(".fmPerStudy / .fmSumStats answer NULL when no method fits the kind", {
    # No individual / sumstat method in the list means there is nothing for
    # that route to run, which is a NULL rather than an empty result object.
    expect_null(pecotmr:::.fmPerStudy(NULL, list(methods = list())))
    expect_null(pecotmr:::.fmSumStats(NULL, list(methods = list())))
})

test_that("a method takes its own constructor for the slot's own path", {
    skip_if_not_installed("susieR")
    skip_if_not_installed("mvsusieR")
    expect_s4_class(
        FineMappingMethodsParam(
            qtlSumStatsMethods = list(
                susie = SusieRssOptions(),
                mvsusie = MvsusieRssOptions()
            )
        ),
        "FineMappingMethodsParam"
    )
    expect_error(
        FineMappingMethodsParam(
            qtlDatasetMethods = list(susie = MvsusieOptions())
        ),
        "uses the 'susie' constructor"
    )
    # susieInf shares susie's ENGINE -- both run susieR::susie on
    # individual data -- so only the token's own constructor pair tells
    # them apart.
    expect_error(
        FineMappingMethodsParam(
            qtlDatasetMethods = list(susieInf = SusieOptions())
        ),
        "uses the 'susieInf' constructor"
    )
    # A plain list declares no engine, so it is judged by its NAMES against
    # the path the slot names -- `z` belongs to susie_rss.
    expect_s4_class(
        FineMappingMethodsParam(qtlSumStatsMethods = list(susie = list(z = 1))),
        "FineMappingMethodsParam"
    )
    expect_error(
        FineMappingMethodsParam(qtlDatasetMethods = list(susie = list(z = 1))),
        "unknown argument"
    )
})

test_that("FineMappingMethodsParam takes lists or constructors per entry", {
    skip_if_not_installed("susieR")
    a <- FineMappingMethodsParam(
        qtlDatasetMethods = list(
            susie = list(L = 5),
            susieInf = SusieInfOptions()
        )
    )
    expect_s4_class(a, "FineMappingMethodsParam")
    expect_setequal(names(a$qtlDatasetMethods), c("susie", "susieInf"))
    expect_equal(a$qtlDatasetMethods$susie$L, 5)
    expect_error(
        FineMappingMethodsParam(qtlDatasetMethods = list(susie = list(LL = 5))),
        "unknown argument"
    )
    expect_error(
        FineMappingMethodsParam(
            qtlDatasetMethods = list(susie = SusieInfOptions())
        ),
        "built with the constructor for 'susieInf'"
    )
    expect_error(
        FineMappingMethodsParam(qtlDatasetMethods = list(nope = list())),
        "not a method this pipeline runs"
    )
})

test_that("argument checking narrows to the engine the input class reaches", {
    skip_if_not_installed("susieR")
    skip_if_not_installed("mvsusieR")
    # A constructor cannot know the input class, so it checks against both
    # entry points. The pipeline knows, so it checks against the one.
    expect_equal(
        pecotmr:::.fmMethodCalleeFor("susie", "QtlDataset"),
        "susieR::susie"
    )
    expect_equal(
        pecotmr:::.fmMethodCalleeFor("susie", "QtlSumStats"),
        "susieR::susie_rss"
    )
})

test_that("an argument for the other entry point is caught per input class", {
    skip_if_not_installed("susieR")
    # `z` is a formal of susie_rss but not of susie, so it is legal on a
    # summary-statistics run and a mistake on an individual-level one. The
    # union check the constructor applies cannot tell these apart.
    norm <- pecotmr:::.fmNormalizeMethods(
        FineMappingMethodsParam(qtlSumStatsMethods = list(susie = list(z = 1))),
        inputKind = "QtlSumStats"
    )
    expect_error(
        pecotmr:::.fmCheckMethodArgsForInput(norm$methodArgs, "QtlDataset"),
        "method 'susie': unknown argument\\(s\\) z"
    )
    expect_silent(
        pecotmr:::.fmCheckMethodArgsForInput(norm$methodArgs, "QtlSumStats")
    )
})

test_that("mvsusie becomes checkable on the individual-level path", {
    skip_if_not_installed("mvsusieR")
    norm <- pecotmr:::.fmNormalizeMethods(
        FineMappingMethodsParam(
            qtlSumStatsMethods = list(mvsusie = list(nope = 1))
        ),
        inputKind = "QtlSumStats"
    )
    # mvsusieR::mvsusie_rss takes `...`, so nothing can be rejected on the
    # RSS path; mvsusieR::mvsusie enumerates its formals, so the QtlDataset
    # path rejects the same argument.
    expect_error(
        pecotmr:::.fmCheckMethodArgsForInput(norm$methodArgs, "QtlDataset"),
        "method 'mvsusie': unknown argument\\(s\\) nope"
    )
    expect_silent(
        pecotmr:::.fmCheckMethodArgsForInput(norm$methodArgs, "QtlSumStats")
    )
})

test_that("a method with no engine for the input is left to the components", {
    # A multi-study run checks against the summary-statistics engine, which
    # fsusie lacks; its arguments are checked by the per-study recursion,
    # whose engine is susiF.
    expect_silent(
        pecotmr:::.fmCheckMethodArgsForInput(
            list(fsusie = list(anyName = 1)),
            "MultiStudyQtlDataset"
        )
    )
})

test_that("pecotmr's own seeded arguments are not checked against the engine", {
    skip_if_not_installed("susieR")
    # L_greedy is a pecotmr concept seeded onto susie-family tokens; it is not
    # a susieR formal and must not be reported as a user mistake.
    norm <- pecotmr:::.fmNormalizeMethods(
        "susieInf",
        maxNumSingleEffects = 5L,
        maxNumSingleEffectsGreedy = 2L,
        inputKind = "QtlDataset"
    )
    expect_true(is_in("L_greedy", names(norm$methodArgs$susieInf)))
    expect_silent(
        pecotmr:::.fmCheckMethodArgsForInput(norm$methodArgs, "QtlDataset")
    )
})

test_that(".fmRssArgsResolved fills rFinite without leaving the class", {
    # `rFinite` defaults to NULL meaning "use the LD panel's N", which is not
    # knowable when the caller builds the record. It is resolved here, and
    # the result must still be an SusieRssParam: the previous code spliced the
    # value in with list_assign(as.list(...)), handing everything downstream
    # a bare list at exactly the point the value became concrete.
    with_mocked_bindings(
        {
            # Nothing needs an N, so the record comes back untouched.
            plain <- SusieRssParam()
            expect_identical(
                pecotmr:::.fmRssArgsResolved(plain, "sketch"),
                plain
            )
            # serFallback needs one, so it is taken from the panel.
            got <- pecotmr:::.fmRssArgsResolved(
                SusieRssParam(serFallback = TRUE),
                "sketch"
            )
            expect_s4_class(got, "SusieRssParam")
            expect_equal(got$rFinite, 4321L)
            expect_true(got$serFallback)
            expect_equal(got$keepFullFit, "fallback")
            # An rMismatch policy needs one too.
            expect_equal(
                pecotmr:::.fmRssArgsResolved(
                    SusieRssParam(rMismatch = "drop"),
                    "sketch"
                )$rFinite,
                4321L
            )
            # A value the caller gave is never overwritten.
            expect_equal(
                pecotmr:::.fmRssArgsResolved(
                    SusieRssParam(rFinite = 99, serFallback = TRUE),
                    "sketch"
                )$rFinite,
                99
            )
        },
        .ldRefNSamples = function(x) 4321L,
        .package = "pecotmr"
    )
})

test_that("each input path offers only the methods it implements", {
    expect_setequal(
        pecotmr:::.fmMethodsFor("QtlDataset"),
        c("susie", "susieInf", "susieAsh", "mvsusie", "fsusie")
    )
    expect_setequal(
        pecotmr:::.fmMethodsFor("QtlSumStats"),
        c("susie", "susieInf", "susieAsh", "ser", "mvsusie")
    )
    # fsusie has no summary-statistics implementation; ser no individual one.
    expect_null(pecotmr:::.fmMethodCalleeFor("fsusie", "QtlSumStats"))
    expect_null(pecotmr:::.fmMethodCalleeFor("ser", "QtlDataset"))
})

test_that("CredibleSetParam owns the single-effect bounds", {
    # They bound how many credible sets can exist. They reach the engine by a
    # different route from the rest of the bundle -- seeded onto each
    # SuSiE-family token as the engine's own `L` / `L_greedy` rather than
    # read by postprocessFinemappingFits() -- but that is a plumbing
    # difference, not a different setting.
    expect_equal(CredibleSetParam()$maxNumSingleEffects, 10L)
    expect_false(
        "maxNumSingleEffectsGreedy" %in% names(CredibleSetParam())
    )
    expect_equal(
        CredibleSetParam(
            maxNumSingleEffectsGreedy = 3L
        )$maxNumSingleEffectsGreedy,
        3L
    )
    # initializeWithSusieInf stands alone: it selects a chained
    # initialisation between methods rather than describing credible sets.
    expect_false(
        "initializeWithSusieInf" %in% names(formals(CredibleSetParam))
    )
    expect_false("fitStructureArgs" %in% getNamespaceExports("pecotmr"))
})

test_that("CredibleSetParam owns the topLoci per-CS columns", {
    # perCsColumns used to be two booleans on fitRetentionArgs -- `fullFit`
    # and `fullFitAlphaOnly`, the second a documented no-op unless the first.
    # They govern the topLoci TABLE, next to includeAllCs which decides the
    # column labels, not the stored fit.
    expect_equal(CredibleSetParam()$perCsColumns, "none")
    expect_equal(CredibleSetParam(perCsColumns = "full")$perCsColumns, "full")
    expect_error(CredibleSetParam(perCsColumns = "yes"), "must be one of")
    # fitRetentionArgs() is retired: one axis is not a bundle.
    expect_false("fitRetentionArgs" %in% getNamespaceExports("pecotmr"))
})

test_that("SignalScreenParam resolves each metric to a pipeline spec", {
    # pip keeps the legacy bare-numeric spelling; the other three resolve to
    # list(metric, cutoff). .asScreen() canonicalizes either.
    expect_null(pecotmr:::.screenResolve(SignalScreenParam()))
    expect_equal(pecotmr:::.screenResolve(SignalScreenParam(pip = 0.5)), 0.5)
    # A per-context pip vector passes through with its names.
    expect_equal(
        pecotmr:::.screenResolve(SignalScreenParam(pip = c(brain = 0.3))),
        c(brain = 0.3)
    )
    expect_equal(
        pecotmr:::.screenResolve(SignalScreenParam(absZ = 5)),
        list(metric = "absZ", cutoff = 5)
    )
    expect_equal(
        pecotmr:::.screenResolve(SignalScreenParam(bf = 100)),
        list(metric = "bf", cutoff = 100)
    )
    expect_equal(
        pecotmr:::.screenResolve(SignalScreenParam(logBf = 3)),
        list(metric = "logBf", cutoff = 3)
    )
})

test_that("SignalScreenParam rejects two metrics and negative absZ / bf", {
    expect_error(
        SignalScreenParam(pip = 0.5, absZ = 5),
        "only one screening metric"
    )
    # 0 is the long-standing "off" spelling, so it is not a second metric.
    expect_equal(
        pecotmr:::.screenResolve(SignalScreenParam(pip = 0, absZ = 5)),
        list(metric = "absZ", cutoff = 5)
    )
    expect_error(SignalScreenParam(absZ = -1), "must be > 0")
    expect_error(SignalScreenParam(bf = -1), "must be > 0")
    # pip < 0 is the adaptive 3 / nVariants convention, and logBf is a log.
    expect_equal(pecotmr:::.screenResolve(SignalScreenParam(pip = -1)), -1)
    expect_equal(
        pecotmr:::.screenResolve(SignalScreenParam(logBf = -2)),
        list(metric = "logBf", cutoff = -2)
    )
})

test_that(".fmSerScreen / .fmScreenActive / .fmSerScreenColumns", {
    set.seed(2)
    X <- matrix(
        rnorm(40),
        20,
        2,
        dimnames = list(paste0("s", 1:20), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    y <- rnorm(20)
    expect_true(pecotmr:::.fmSerScreen(X, y, 0)) # disabled
    expect_true(pecotmr:::.fmSerScreen(X, c(1, rep(NA, 19)), 0.5)) # < 2 obs
    expect_type(pecotmr:::.fmSerScreen(X, y, 0.5), "logical") # real susie fit
    local_mocked_bindings(
        susie = function(...) stop("boom"),
        .package = "susieR"
    )
    expect_true(pecotmr:::.fmSerScreen(X, y, 0.5)) # fit fails -> keep
    expect_false(pecotmr:::.fmScreenActive(0))
    expect_true(pecotmr:::.fmScreenActive(0.5))
    expect_length(
        pecotmr:::.fmSerScreenColumns(X, matrix(rnorm(40), 20, 2), 0),
        2L
    )
})

test_that("every recursion forwards each cfg field its callee accepts", {
    ns <- asNamespace("pecotmr")
    cases <- list(
        list(
            site = ".fmPerStudy",
            builder = ".fmPipelineMultiStudy",
            callee = ".fmPipelineQtlDataset"
        ),
        list(
            site = ".fmSumStats",
            builder = ".fmPipelineMultiStudy",
            callee = ".fmPipelineQtlSumStats"
        ),
        list(
            site = ".twasPerStudy",
            builder = ".twasMsCfg",
            callee = ".twasPipelineQtlDataset"
        ),
        list(
            site = ".twasSumStats",
            builder = ".twasMsCfg",
            callee = ".twasPipelineQtlSumStats"
        )
    )
    for (cs in cases) {
        fields <- .cfgFieldNames(cs$builder)
        accepted <- setdiff(names(formals(get(cs$callee, envir = ns))), "...")
        passed <- .argNamesPassed(cs$site)
        # A field the record carries AND the callee accepts must be passed.
        dropped <- setdiff(intersect(fields, accepted), passed)
        expect_equal(
            dropped,
            character(0),
            info = paste0(
                cs$site,
                " never forwards: ",
                paste(dropped, collapse = ", "),
                " -- the callee accepts them and the record carries them, so ",
                "they silently fall back to ",
                cs$callee,
                "'s defaults."
            )
        )
    }
})

test_that("MultiStudyQtlDataset takes every component argument", {
    # Same contract as the TWAS twin: every argument a component method takes
    # is a MultiStudy argument, and the record handed to the components
    # carries it (the test above checks record -> component).
    ms <- .methodFormalNames("fineMappingPipeline", "MultiStudyQtlDataset")
    comp <- union(
        .methodFormalNames("fineMappingPipeline", "QtlDataset"),
        .methodFormalNames("fineMappingPipeline", "QtlSumStats")
    )
    expect_equal(setdiff(comp, ms), character(0))
    expect_equal(setdiff(ms, comp), character(0))
    expect_equal(
        setdiff(
            ms,
            c(.cfgFieldNames(".fmPipelineMultiStudy"), "jointSpecification")
        ),
        character(0)
    )
})

test_that("QtlDataset: a CV-only run keeps the CV and skips the full fit", {
    skip_if_not_installed("susieR")
    qd <- .fmp_qdA
    # The full-data post-processor stops the test if it is reached: a
    # cross-validation-only run must not fit, or summarize, the full data.
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmPostprocessOne = function(...) stop("the full-data fit ran"),
        .package = "pecotmr"
    )
    res <- suppressMessages(fineMappingPipeline(
        qd,
        methods = "susie",
        cisWindow = 1000L,
        initializeWithSusieInf = FALSE,
        crossValidationParam = CrossValidationParam(
            folds = 3,
            fitFullData = FALSE
        ),
        verbose = 0
    ))
    expect_equal(nrow(res), 1L)
    expect_equal(unname(lengths(res)), 0L)
    expect_null(susieFit(
        res,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "susie"
    ))
    cv <- cvResult(
        res,
        studyName = "study1",
        context = "brain",
        trait = "ENSG_A",
        method = "susie"
    )
    expect_true("susie_predicted" %in% names(cv$prediction))
    # Fold-level results cannot be published, so the export refuses.
    expect_error(writeBed(res, tempfile(fileext = ".bed")), "fitFullData")
})

test_that(".fmMergeEntries keeps a merged cross-validation-only row fit-free", {
    # A list of per-region NULL fits would read as a fit; a merged row with no
    # full-data fit in any region keeps NULL, and each region's CV.
    cv <- list(samplePartition = data.frame(Sample = "s1", Fold = 1L))
    merged <- pecotmr:::.fmMergeEntries(list(
        pecotmr:::.fmCvOnlyRow(cv),
        pecotmr:::.fmCvOnlyRow(cv)
    ))
    expect_null(pecotmr:::.fmrPartsSusieFit(merged))
    expect_length(pecotmr:::.fmrPartsVariantIds(merged), 0L)
    expect_named(cvResult(merged), c("region1", "region2"))
})

test_that(".fmCacheForRun reuses only rows that carry what the run asks for", {
    vid <- "chr1:100:A:G"
    tl <- tibble::tibble(variant_id = vid, pip = 0.5)
    cv <- list(samplePartition = data.frame(Sample = "s1", Fold = 1L))
    cache <- QtlFineMappingResult(
        studyName = rep("s", 3L),
        context = rep("c", 3L),
        trait = c("fullCv", "full", "cvOnly"),
        method = rep("susie", 3L),
        entry = list(
            FineMappingRow(vid, list(a = 1), tl, cvResult = cv),
            FineMappingRow(vid, list(a = 1), tl),
            pecotmr:::.fmCvOnlyRow(cv)
        )
    )
    kept <- function(cvp) {
        out <- pecotmr:::.fmCacheForRun(cache, pecotmr:::.cvResolve(cvp))
        as.character(out$trait)
    }
    expect_equal(kept(CrossValidationParam()), c("fullCv", "full"))
    expect_equal(kept(CrossValidationParam(folds = 3)), "fullCv")
    expect_equal(
        kept(CrossValidationParam(folds = 3, fitFullData = FALSE)),
        c("fullCv", "cvOnly")
    )
})

test_that("CredibleSetParam fields read and replace, L is fixed", {
    cs <- CredibleSetParam(coverage = 0.95)
    expect_equal(csCoverage(cs), 0.95)
    expect_equal(maxNumSingleEffects(cs), 10L)
    expect_null(medianAbsCorr(cs))

    # Copy semantics: replacement rebinds a modified copy, so another name
    # holding the original still sees the old value.
    cs2 <- cs
    cs2$coverage <- 0.9
    expect_equal(csCoverage(cs2), 0.9)
    expect_equal(csCoverage(cs), 0.95)
    expect_s4_class(cs2, "CredibleSetParam")

    # The slot's declared type still holds through a replacement.
    expect_error(cs$coverage <- "wide", "coverage")

    # L / Lgreedy are deliberately fixed: they are seeded onto each
    # SuSiE-family token's own arguments, so changing one needs a refit, not
    # a re-summary. A replacement would hand back a record that looks
    # applied, so `$<-` refuses them.
    expect_error(cs$maxNumSingleEffects <- 5L, "needs a refit")
    expect_error(cs$maxNumSingleEffectsGreedy <- 2L, "needs a refit")
    expect_true(is.function(maxNumSingleEffectsGreedy))
})

test_that("every CredibleSetParam field has a getter", {
    # The getter is the field name itself, except `coverage`, whose bare name
    # is an IRanges generic -- so that one is named here, not derived.
    renamed <- c(coverage = "csCoverage")
    ns <- asNamespace("pecotmr")
    for (f in names(formals(CredibleSetParam))) {
        accessor <- if (is_in(f, names(renamed))) renamed[[f]] else f
        expect_true(
            exists(accessor, envir = ns, inherits = FALSE),
            info = paste("no accessor for", f)
        )
    }
})

test_that("CredibleSetParam getters read and $<- replaces every field", {
    # Dispatching each getter is what exercises the setGeneric bodies in
    # AllGenerics.R; existence alone never reaches them.
    cs <- CredibleSetParam()
    expect_equal(secondaryCoverage(cs), c(0.7, 0.5))
    expect_equal(signalCutoff(cs), 0.025)
    expect_equal(minAbsCorr(cs), 0.8)
    expect_false(includeAllCs(cs))
    expect_equal(perCsColumns(cs), "none")
    expect_null(maxNumSingleEffectsGreedy(cs))

    # Replacement on a copy; `cs` itself keeps its original values.
    edited <- cs
    edited$secondaryCoverage <- 0.6
    edited$signalCutoff <- 0.1
    edited$minAbsCorr <- 0.5
    edited$medianAbsCorr <- 0.4
    edited$includeAllCs <- TRUE
    edited$perCsColumns <- "alpha"
    expect_equal(secondaryCoverage(edited), 0.6)
    expect_equal(signalCutoff(edited), 0.1)
    expect_equal(minAbsCorr(edited), 0.5)
    expect_equal(medianAbsCorr(edited), 0.4)
    expect_true(includeAllCs(edited))
    expect_equal(perCsColumns(edited), "alpha")
    expect_equal(secondaryCoverage(cs), c(0.7, 0.5))
    expect_false(includeAllCs(cs))

    # Replacement goes through the constructor, so a wrong-typed value is
    # refused exactly as CredibleSetParam() refuses it.
    expect_error(cs$includeAllCs <- "yes", "includeAllCs")
})

# ---------------------------------------------------------------------------
# CredibleSetParam validity
# ---------------------------------------------------------------------------

test_that("CredibleSetParam validity accepts every legitimate construction", {
    expect_s4_class(CredibleSetParam(), "CredibleSetParam")
    expect_s4_class(
        CredibleSetParam(coverage = 0.9, includeAllCs = TRUE),
        "CredibleSetParam"
    )
    expect_s4_class(
        CredibleSetParam(
            maxNumSingleEffects = 4L,
            maxNumSingleEffectsGreedy = 2L
        ),
        "CredibleSetParam"
    )
    expect_s4_class(
        CredibleSetParam(perCsColumns = "alpha"),
        "CredibleSetParam"
    )
    # No secondary coverages is a legitimate request, not a missing value.
    expect_s4_class(
        CredibleSetParam(secondaryCoverage = numeric(0)),
        "CredibleSetParam"
    )
    # medianAbsCorr is optional.
    expect_null(medianAbsCorr(CredibleSetParam()))
})

test_that("CredibleSetParam validity catches range and length errors", {
    # The slot type admits any numeric, so without validity each of these
    # was accepted -- out of range, wrong length, and empty alike.
    expect_error(CredibleSetParam(coverage = 1.5), "coverage")
    expect_error(CredibleSetParam(coverage = c(0.9, 0.95)), "coverage")
    expect_error(CredibleSetParam(coverage = numeric(0)), "coverage")
    expect_error(CredibleSetParam(minAbsCorr = -1), "minAbsCorr")
    expect_error(CredibleSetParam(signalCutoff = 2), "signalCutoff")
    expect_error(
        CredibleSetParam(secondaryCoverage = c(0.7, 3)),
        "secondaryCoverage"
    )
    expect_error(
        CredibleSetParam(maxNumSingleEffects = 0),
        "maxNumSingleEffects"
    )
    expect_error(
        CredibleSetParam(maxNumSingleEffects = 2.5),
        "maxNumSingleEffects"
    )
    expect_error(
        CredibleSetParam(maxNumSingleEffectsGreedy = 0),
        "maxNumSingleEffectsGreedy"
    )
})

test_that("CredibleSetParam replacements are held to the same checks", {
    cs <- CredibleSetParam()
    # A replacement is a rebuild through CredibleSetParam(), so it is refused
    # exactly where the constructor would refuse the value.
    expect_error(cs$coverage <- 1.5, "coverage")
    expect_error(cs$coverage <- c(0.9, 0.95), "coverage")
    expect_error(cs$minAbsCorr <- -1, "minAbsCorr")
    expect_error(cs$signalCutoff <- 2, "signalCutoff")
    expect_error(cs$secondaryCoverage <- c(0.7, 3), "secondaryCoverage")
    # perCsColumns was the sharpest asymmetry of the old setters: the
    # constructor refused a bad level via arg_match() while they accepted it.
    expect_error(cs$perCsColumns <- "bogus", "perCsColumns")
    # Valid replacements still go through.
    ok <- cs
    ok$coverage <- 0.9
    ok$minAbsCorr <- 0.5
    ok$medianAbsCorr <- 0.4
    expect_equal(csCoverage(ok), 0.9)
    expect_equal(minAbsCorr(ok), 0.5)
    expect_equal(medianAbsCorr(ok), 0.4)
})

test_that("CredibleSetParam validity rejects non-scalar counts and flags", {
    # Reached by poking a VALID object: new() straight from the prototype
    # fails on `coverage` first, and the constructor guards these before
    # validity ever sees them.
    p <- CredibleSetParam()
    q <- p
    q@maxNumSingleEffects <- numeric(0)
    expect_error(validObject(q), "maxNumSingleEffects must be a single")
    q <- p
    q@maxNumSingleEffects <- NA_real_
    expect_error(validObject(q), "maxNumSingleEffects must be a single")
    q <- p
    q@includeAllCs <- NA
    expect_error(validObject(q), "includeAllCs must be a single")
    q <- p
    q@includeAllCs <- c(TRUE, FALSE)
    expect_error(validObject(q), "includeAllCs must be a single")
    # A nullable slot left NULL is fine -- nullability is the slot type's
    # job, so the value checks simply pass a NULL through.
    q <- p
    q@maxNumSingleEffectsGreedy <- NULL
    expect_true(validObject(q))
})

test_that("SerOptions and FsusieOptions refuse what their wrappers own", {
    # The susie_ser and susiF arms of .fmOwnedEngineArgs(): ser takes only
    # the coverage half of csHome, and susiF owns the three settings the
    # pipeline actually supplies to it.
    expect_error(SerOptions(z = 1), "supplied by pecotmr")
    expect_error(SerOptions(n = 100L), "supplied by pecotmr")
    expect_error(SerOptions(coverage = 0.9), "csCoverage")
    expect_error(FsusieOptions(L = 2), "maxNumSingleEffects")
    expect_error(FsusieOptions(min_purity = 0.5), "minAbsCorr")
    expect_error(FsusieOptions(cov_lev = 0.9), "csCoverage")
    # ... and ONLY those. max_scale / max_SNP_EM / prior were refused while
    # this row still described the retired fsusieWrapper(); nothing supplies
    # them on the pipeline path, so refusing them closed the only way to set
    # them at all.
    expect_s4_class(FsusieOptions(max_scale = 6), "MethodOptions")
    expect_s4_class(FsusieOptions(max_SNP_EM = 20), "MethodOptions")
    expect_s4_class(FsusieOptions(prior = "mixture_normal"), "MethodOptions")
    # ser does NOT own L -- susie_ser is single-effect by construction, so
    # the csHome entry it takes is coverage alone.
    expect_s4_class(SerOptions(prior_variance = 1), "MethodOptions")
})

test_that(".fmOwnedEngineArgs owns nothing for an unknown engine", {
    expect_equal(
        pecotmr:::.fmOwnedEngineArgs("susieR::not_a_real_engine"),
        character(0)
    )
})

test_that(".fmBlockWindow falls back to the trait's own cis window", {
    data(qtlDatasetExample)
    explicit <- GenomicRanges::GRanges("chr1:1-999")
    # An explicit block range wins untouched.
    expect_identical(
        .fmBlockWindow(explicit, "ctx", "t", list(data = qtlDatasetExample)),
        explicit
    )
    # Without one, the cis path's window is the trait's position +/- cisWindow
    # -- the same interval .fmResidGenoBlock() slices X over. Widening the
    # cis window must widen the recorded window.
    cfgNarrow <- list(data = qtlDatasetExample, cisWindow = 1000)
    cfgWide <- list(data = qtlDatasetExample, cisWindow = 500000)
    ctx <- contexts(qtlDatasetExample)[[1L]]
    tid <- rownames(molecularTraits(qtlDatasetExample, context = ctx)[[ctx]])[[
        1L
    ]]
    narrow <- .fmBlockWindow(NULL, ctx, tid, cfgNarrow)
    wide <- .fmBlockWindow(NULL, ctx, tid, cfgWide)
    expect_lt(GenomicRanges::start(wide), GenomicRanges::start(narrow))
    expect_gt(GenomicRanges::end(wide), GenomicRanges::end(narrow))
    expect_identical(
        GenomicRanges::width(wide) - GenomicRanges::width(narrow),
        2L * (500000L - 1000L)
    )
})

test_that("a GWAS row carries its block name and window into its fit", {
    data(gwasSumStatsS4Example)
    gr <- .collectionEntry(gwasSumStatsS4Example, 1L)
    derived <- .fmGwasRangeKey(gr)

    # (a) No windows: GwasSumStats() names the row by its seqname, which the
    # fit keeps, and its window is the span it was built with.
    bare <- GwasSumStats(studyName = "s", entry = list(gr), genome = "hg19")
    block <- .fmGwasBlock(bare, 1L, gr, .rowRegions(bare))
    expect_identical(block$id, "chr22")
    expect_false(block$given)
    expect_identical(as.character(block$region), as.character(range(gr)))

    # (b) A range the caller gave -- here a `region` window keyed by name --
    # IS the block boundary and keeps its key. Adjacent blocks leave gaps, so
    # it is wider than the variants that survived QC.
    named <- GwasSumStats(
        studyName = "s",
        entry = list(gr),
        genome = "hg19",
        region = GenomicRanges::GRanges(
            "chr22",
            IRanges::IRanges(1L, 50000000L),
            blockId = "b7"
        )
    )
    block <- .fmGwasBlock(named, 1L, gr, .rowRegions(named))
    expect_identical(block$id, "b7")
    expect_true(block$given)
    expect_identical(as.character(block$region), "chr22:1-50000000")

    # (c) A name is kept as it is; only a row with none, from a collection
    # saved before the column existed, is named by its span.
    expect_identical(
        .fmGwasBlockName("chr22_1_50000000", gr),
        "chr22_1_50000000"
    )
    expect_identical(.fmGwasBlockName(NA_character_, gr), derived)

    # (d) The resume cache keys on the DERIVED range either way, so supplying a
    # wider boundary cannot make the cache miss.
    expect_identical(.fmGwasRangeKey(gr), derived)
})

# ---------------------------------------------------------------------------
# Provenance columns (dataType / dataSource / standardized)
#
# These three are what writeBed() / writeGwasVcf() publish as data_type,
# data_source and weight_standardized. Only the first is an outside fact; the
# other two follow from the branch the pipeline took, so they are derived and
# the summary-statistics methods refuse to be told.
# ---------------------------------------------------------------------------

# @noRd
.fmprov_cols <- function(res) {
    as.data.frame(S4Vectors::mcols(res))[,
        intersect(
            c("dataType", "dataSource", "standardized"),
            colnames(S4Vectors::mcols(res))
        ),
        drop = FALSE
    ]
}

test_that("fineMappingPipeline(QtlDataset) records what it is given", {
    data(qtlDatasetExample, envir = environment())
    res <- suppressMessages(fineMappingPipeline(
        qtlDatasetExample,
        methods = "susie",
        cisWindow = 1e6,
        dataType = "bulk_rnaseq"
    ))
    cols <- .fmprov_cols(res)
    expect_equal(unique(cols$dataType), "bulk_rnaseq")
    # Derived, not supplied: this entry point fits individual-level genotypes.
    expect_equal(unique(cols$dataSource), "individual")
    # Also derived: the genotypes were standardized on the way in, so there
    # is nothing for a caller to declare and no way for it to be FALSE.
    expect_true(all(cols$standardized))
})

test_that("fineMappingPipeline(QtlDataset) omits dataType when not supplied", {
    data(qtlDatasetExample, envir = environment())
    res <- suppressMessages(fineMappingPipeline(
        qtlDatasetExample,
        methods = "susie",
        cisWindow = 1e6
    ))
    # An absent label leaves the column out rather than filling it with NA --
    # an all-NA column is dropped on export anyway, and its absence is what
    # the reader sees either way.
    expect_false("dataType" %in% colnames(S4Vectors::mcols(res)))
    cols <- .fmprov_cols(res)
    expect_equal(unique(cols$dataSource), "individual")
    expect_true(all(cols$standardized))
})

test_that("fineMappingPipeline(QtlSumStats) records sumstat / standardized", {
    res <- suppressWarnings(suppressMessages(fineMappingPipeline(
        .rssf_ss,
        methods = "susie",
        dataType = "bulk_rnaseq"
    )))
    cols <- .fmprov_cols(res)
    expect_equal(unique(cols$dataType), "bulk_rnaseq")
    expect_equal(unique(cols$dataSource), "sumstat")
    # The RSS likelihood is written in standardized units, so the posterior
    # effects are on that scale whatever the caller believes.
    expect_true(all(cols$standardized))
})

test_that("no pipeline takes a weight-scale declaration any more", {
    # It was an argument while the caller's genotypes went to the fitter
    # untouched. They are standardized on the way in now, so the flag is
    # derived and a declaration could only contradict it.
    expect_error(
        fineMappingPipeline(
            .rssf_ss,
            methods = "susie",
            weightStandardized = FALSE
        ),
        "unused argument"
    )
    data(qtlDatasetExample, envir = environment())
    expect_error(
        fineMappingPipeline(
            qtlDatasetExample,
            methods = "susie",
            cisWindow = 1e6,
            weightStandardized = FALSE
        ),
        "unused argument"
    )
})

test_that(".fmMsDataSource leaves a mixed collection's joint rows NA", {
    expect_true(is.na(pecotmr:::.fmMsDataSource(.fmp_msqd)))
})

test_that(".fmMsDataSource gives an individual-level collection's joint rows", {
    # A collection always holds individual-level studies (its validity
    # requires one), so without summary statistics its joint rows are those.
    mt <- MultiStudyQtlDataset(
        qtlDatasets = list(
            A = .fmp_twoContextStudy("A"),
            B = .fmp_twoContextStudy("B")
        )
    )
    expect_equal(pecotmr:::.fmMsDataSource(mt), "individual")
})

test_that("individual-level genotypes reach the fitter standardized", {
    # Through residualizedGenotypes(), which owns the scaling.
    # The single preparation point for every individual-level fit in the
    # package: fine-mapping's per-tuple and joint paths, the TWAS joint
    # engine, and the top-PC construction all come through it.
    data(qtlDatasetExample, envir = environment())
    X <- pecotmr:::.fmResidGeno(
        qtlDatasetExample,
        ResidualizationParam(),
        genotypeFilter = pecotmr:::.qtlResolveFilter(),
        contexts = "brain",
        traitId = "ENSG_example",
        cisWindow = 1e6
    )
    expect_true(max(abs(colMeans(X))) < 1e-10)
    sds <- matrixStats::colSds(X)
    expect_true(all(abs(sds[sds > 0] - 1) < 1e-10))
})

test_that("the genotype scale comes from the accessor, not the pipeline", {
    # ONE knob: ResidualizationParam(scaleResiduals=), TRUE when unset. The
    # pipelines add no mechanism of their own -- they only RECORD what the
    # residualizer resolved, so a wrong record and a wrong scale cannot come
    # apart.
    data(qtlDatasetExample, envir = environment())
    recorded <- function(..., dataset = qtlDatasetExample) {
        res <- suppressMessages(fineMappingPipeline(
            dataset,
            methods = "susie",
            cisWindow = 1e6,
            residualizationParam = ResidualizationParam(...)
        ))
        unique(as.logical(S4Vectors::mcols(res)$standardized))
    }
    expect_true(recorded())
    expect_true(recorded(scaleResiduals = TRUE))
    expect_false(recorded(scaleResiduals = FALSE))
})

test_that("the pipelines add no second standardization mechanism", {
    # The accessor owns it. A parallel step inside .fmResidGeno() would have
    # covered only its own callers -- colocboostPipeline() calls
    # residualizedGenotypes() directly and would have been left out.
    expect_false(exists(
        ".fmStandardizeGenotypes",
        envir = asNamespace("pecotmr")
    ))
    expect_false(is_in(
        "standardizeGenotypes",
        names(formals(pecotmr:::.fmPipelineQtlDataset))
    ))
})

test_that("the hatch does not reach a summary-statistics fit", {
    # There are no genotypes to scale, and the RSS likelihood is in
    # standardized units whatever anyone asks, so the flag stays TRUE.
    res <- suppressWarnings(suppressMessages(
        fineMappingPipeline(
            .rssf_ss,
            methods = "susie",
            residualizationParam = ResidualizationParam(scaleResiduals = FALSE)
        )
    ))
    expect_true(all(as.logical(S4Vectors::mcols(res)$standardized)))
})

# ---- joint fits receive the caller's settings ------------------------------

test_that("QtlDataset: joint fits get every caller setting", {
    # The joint engine reads residualization from its `args`. The dispatcher
    # filed it only in the pipeline config, which nothing reads, so joint fits
    # ran on default residualization while the result recorded the caller's
    # scaling as its provenance.
    cap <- new.env(parent = emptyenv())
    .captureJointEngine(cap)
    rp <- ResidualizationParam(scaleResiduals = FALSE)
    expect_error(
        suppressMessages(fineMappingPipeline(
            .fmp_qdTwoCtxAB,
            methods = "mvsusie",
            traitId = "ENSG_A",
            cisWindow = 1000L,
            jointSpecification = "context",
            residualizationParam = rp,
            crossValidationParam = CrossValidationParam(folds = 3),
            signalScreenParam = SignalScreenParam(pip = 0.1),
            fitRetention = "full",
            seed = 7L
        )),
        "no joint fits produced"
    )
    ind <- .capturedCall(cap, "individual")
    expect_identical(ind$args$residualizationParam, rp)
    expect_equal(ind$config$crossValidationParam$folds, 3)
    expect_equal(ind$config$fitRetention, "full")
    expect_equal(ind$config$seed, 7L)
    expect_equal(
        ind$args$pipCutoffToSkip,
        pecotmr:::.screenResolve(SignalScreenParam(pip = 0.1))
    )
})

test_that("MultiStudyQtlDataset: per-study joint fits get every setting", {
    # The per-study joint dispatch re-listed its own subset of the settings
    # and dropped cross-validation, residualization, the signal screen,
    # fitRetention, the seed and the genotype-filter overrides.
    cap <- new.env(parent = emptyenv())
    .captureJointEngine(cap)
    rp <- ResidualizationParam(scaleResiduals = FALSE)
    expect_error(
        suppressMessages(fineMappingPipeline(
            .fmp_msqd,
            methods = "mvsusie",
            cisWindow = 1000L,
            jointSpecification = "context",
            residualizationParam = rp,
            crossValidationParam = CrossValidationParam(folds = 3),
            signalScreenParam = SignalScreenParam(pip = 0.1),
            genotypeFilterParam = GenotypeFilterParam(mafCutoff = 0.05),
            fitRetention = "full",
            seed = 7L
        )),
        "no joint fits produced"
    )
    ind <- .capturedCall(cap, "individual")
    expect_identical(ind$args$residualizationParam, rp)
    expect_equal(ind$config$crossValidationParam$folds, 3)
    expect_equal(ind$config$fitRetention, "full")
    expect_equal(ind$config$seed, 7L)
    expect_equal(
        ind$args$pipCutoffToSkip,
        pecotmr:::.screenResolve(SignalScreenParam(pip = 0.1))
    )
    expect_equal(ind$args$genotypeFilter$mafCutoff, 0.05)
})

test_that("MultiStudyQtlDataset: cross-study joint fits get their settings", {
    # The cross-study dispatch omitted fitRetention, which has no default
    # further down, so any study-axis spec failed before fitting anything.
    cap <- new.env(parent = emptyenv())
    .captureJointEngine(cap)
    panel <- PanelFilterParam(mafCutoff = 0.05)
    expect_error(
        suppressMessages(fineMappingPipeline(
            .fmp_msqd,
            methods = "mvsusie",
            jointSpecification = "study",
            fitRetention = "full",
            panelFilterParam = panel
        )),
        "no joint fits produced"
    )
    ss <- .capturedCall(cap, "sumstats")
    expect_equal(ss$config$fitRetention, "full")
    expect_equal(ss$args$cutoffs, pecotmr:::.panelCutoffs(panel))
})

# =============================================================================
# The LD reference dispatch: an LdMatrix or a GenotypeMixture as ldReference
# =============================================================================

# The bundled GWAS collection, with its panel's LD copied into an LdMatrix
# and used as the reference instead.
.fmp_withLdMatrix <- function(gss, nSamples = ncol(ldReference(gss)[[1L]])) {
    panel <- ldReference(gss)[[1L]]
    ld <- LdMatrix(
        unname(computeLd(panel)),
        variants = SummarizedExperiment::rowRanges(panel),
        nSamples = nSamples
    )
    GwasSumStats(
        studyName = studyName(gss),
        entry = list(.collectionEntry(gss, 1L)),
        genome = "hg19",
        ldReference = ld,
        qcOptions = qcOptions(gss)
    )
}

# A real susie fit. L = 4 still finds the example's three credible sets, and
# fits cheaper than the default L = 10; every route compared uses it.
# @noRd
.fmp_fitSusie <- function(x) {
    suppressWarnings(suppressMessages(fineMappingPipeline(
        x,
        methods = "susie",
        initializeWithSusieInf = FALSE,
        credibleSetParam = CredibleSetParam(maxNumSingleEffects = 4L)
    )))
}

test_that("fine-mapping with an LdMatrix reference matches its panel", {
    # Every step reads LD through the reference dispatch, so a matrix holding
    # the panel's LD gives the panel's fit.
    data(gwasSumStatsS4Example)
    gss <- gwasSumStatsS4Example
    withLd <- .fmp_withLdMatrix(gss)
    viaMatrix <- .fmp_fitSusie(withLd)
    # The panel's fit reads the panel's factor; the matrix's reads the matrix.
    viaPanel <- .fmp_fitSusie(gss)
    expect_equal(pip(viaMatrix), pip(viaPanel))
    expect_equal(credibleSets(viaMatrix), credibleSets(viaPanel))
    expect_s4_class(ldReference(viaMatrix)[[1L]], "LdMatrix")
})

test_that("an RSS block reads a genotype reference's factor, not its LD", {
    data(gwasSumStatsS4Example)
    gss <- gwasSumStatsS4Example
    panel <- ldReference(gss)[[1L]]
    ids <- names(SummarizedExperiment::rowRanges(panel))[1:25]
    plain <- .fmRssLd(panel, ids, "susie", FALSE)
    expect_null(plain$R)
    expect_equal(
        crossprod(plain$X) / (nrow(plain$X) - 1),
        .ldFromRef(panel, ids, label = "t")
    )
    # susieInf (and a chain initialised from it) and susieAsh need the matrix.
    for (toRun in list(c("susie", "susieInf"), "susieAsh")) {
        both <- .fmRssLd(panel, ids, toRun, TRUE)
        expect_equal(both$R, .ldFromRef(panel, ids, label = "t"))
        expect_identical(both$X, plain$X)
    }
    ldOnly <- .fmRssLd(
        ldReference(.fmp_withLdMatrix(gss))[[1L]],
        ids,
        "susie",
        FALSE
    )
    expect_null(ldOnly$X)
    expect_equal(ldOnly$R, .ldFromRef(panel, ids, label = "t"))
})

test_that("fine-mapping runs on a GenotypeMixture reference", {
    data(gwasSumStatsS4Example)
    gss <- gwasSumStatsS4Example
    panel <- ldReference(gss)[[1L]]
    mix <- GenotypeMixture(
        list(a = panel[, 1:80], b = panel[, 81:165]),
        mixtureWeights = c(0.5, 0.5)
    )
    withMix <- GwasSumStats(
        studyName = studyName(gss),
        entry = list(.collectionEntry(gss, 1L)),
        genome = "hg19",
        ldReference = mix,
        qcOptions = qcOptions(gss)
    )
    res <- .fmp_fitSusie(withMix)
    expect_s4_class(ldReference(res)[[1L]], "GenotypeMixture")
    expect_true(all(is.finite(unlist(pip(res)))))
})

test_that("the rFinite default needs the reference's size", {
    unknown <- LdMatrix(
        diag(2),
        variants = data.frame(chrom = "1", pos = 1:2, A1 = "A", A2 = "G")
    )
    expect_error(
        .fmResolveRFinite(NULL, TRUE, "none", unknown),
        "needs the LD reference's size"
    )
    known <- LdMatrix(
        diag(2),
        variants = data.frame(chrom = "1", pos = 1:2, A1 = "A", A2 = "G"),
        nSamples = 300
    )
    expect_identical(.fmResolveRFinite(NULL, TRUE, "none", known), 300L)
    # An explicit rFinite, or no mismatch mode, never asks.
    expect_identical(.fmResolveRFinite(50, TRUE, "none", unknown), 50)
    expect_null(.fmResolveRFinite(NULL, FALSE, "none", unknown))
})

# ===========================================================================
# Per-row LD references
# ===========================================================================

# Two QtlSumStats rows (contexts c1, c2) carrying the same z over five
# variants, each on the LD reference given for it.
.fmp_twoRefSs <- function(refs) {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(100L * (1:5), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:5)),
        A1 = "A",
        A2 = "G",
        Z = c(6, 5.5, 1, 0.5, 0.2),
        N = 1000L
    )
    QtlSumStats(
        studyName = c("S", "S"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        entry = list(gr, gr),
        genome = "hg19",
        ldReference = refs,
        qcOptions = list(step1 = "ok")
    )
}

# An LdMatrix over those five variants, with nSamples recorded.
.fmp_ldm <- function(rho) {
    ids <- sprintf("chr1:%d:A:G", 100L * (1:5))
    LdMatrix(
        `dimnames<-`(rho^abs(outer(1:5, 1:5, "-")), list(ids, ids)),
        nSamples = 1000L
    )
}

test_that("each QtlSumStats row is fine-mapped against its own LD reference", {
    weak <- .fmp_ldm(0.1)
    strong <- .fmp_ldm(0.95)
    res <- suppressMessages(suppressWarnings(
        fineMappingPipeline(
            .fmp_twoRefSs(list(weak, strong)),
            methods = "susie"
        )
    ))
    c1 <- which(res$context == "c1")
    c2 <- which(res$context == "c2")
    # The result keeps each row on the reference its fit read.
    expect_identical(pecotmr:::.rowLdRef(res, c1), weak)
    expect_identical(pecotmr:::.rowLdRef(res, c2), strong)
    # Identical z under different LD gives different fits: each row was fit
    # against its own reference, not the first one's.
    pip <- function(i) {
        pecotmr:::.fmrPartsSusieFit(pecotmr:::.fmrRowParts(res, i))$pip
    }
    expect_false(isTRUE(all.equal(pip(c1), pip(c2))))
})

test_that("mvsusie on QtlSumStats refuses contexts on two LD references", {
    ss <- .fmp_twoRefSs(list(.fmp_ldm(0.1), .fmp_ldm(0.95)))
    expect_error(
        suppressMessages(fineMappingPipeline(ss, methods = "mvsusie")),
        "a joint fit uses one LD matrix"
    )
})

# ===========================================================================
# A GWAS row spanning a whole chromosome is fine-mapped, but said
# ===========================================================================

# A study's row over `starts` on chr1.
.fmp_spanningGwas <- function(starts, studies = "G1") {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(starts, width = 1L))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", starts),
        A1 = "A",
        A2 = "G",
        Z = seq_along(starts) / 10,
        N = 1000L
    )
    GwasSumStats(
        studyName = studies,
        entry = rep(list(gr), length(studies)),
        genome = "hg19",
        qcOptions = list(step1 = "ok")
    )
}

test_that(".fmWarnChromosomeRows flags rows spanning a chromosome's worth", {
    # An LD block's worth of variants passes in silence.
    expect_silent(
        pecotmr:::.fmWarnChromosomeRows(.fmp_spanningGwas(c(100L, 5e6)))
    )
    expect_warning(
        pecotmr:::.fmWarnChromosomeRows(.fmp_spanningGwas(c(100L, 3e7))),
        "20 Mb, as a whole chromosome does \\(G1 chr1_100_30000000\\)"
    )
    # Past five rows the list is cut short.
    six <- .fmp_spanningGwas(c(100L, 3e7), studies = str_c("G", 1:6))
    expect_warning(
        pecotmr:::.fmWarnChromosomeRows(six),
        "6 GWAS row\\(s\\).*G5 chr1_100_30000000, \\.\\.\\.\\)"
    )
})

test_that("fineMappingPipeline(QtlSumStats): cisWindow sets each window", {
    local_mocked_bindings(
        extractBlockGenotypes = .fmp_mockExtractor(),
        .fmFitSusieRss = .fmp_mockFitRss(),
        .fmPostprocessOne = .fmp_mockPostprocess(),
        .package = "pecotmr"
    )
    ss <- .fmp_qss
    ss$traitPos <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(300L, width = 1L)
    )
    res <- suppressMessages(fineMappingPipeline(
        ss,
        methods = "susie",
        cisWindow = 100,
        initializeWithSusieInf = FALSE
    ))
    expect_identical(as.character(genomicRegion(res)), "chr1:200-400")
    expect_equal(length(res[[1L]]), 3L)
    # The selection is parsed before the window is laid, so naming a context
    # works alongside it.
    picked <- suppressMessages(fineMappingPipeline(
        ss,
        methods = "susie",
        contexts = "c1",
        cisWindow = 100,
        initializeWithSusieInf = FALSE
    ))
    expect_identical(as.character(genomicRegion(picked)), "chr1:200-400")
    # Without one, the result keeps the tuple's own window.
    plain <- suppressMessages(fineMappingPipeline(
        .fmp_qss,
        methods = "susie",
        initializeWithSusieInf = FALSE
    ))
    expect_identical(as.character(genomicRegion(plain)), "chr1:100-500")
    expect_error(
        fineMappingPipeline(
            .fmp_qss,
            methods = "susie",
            cisWindow = 100
        ),
        "laid around each tuple's trait position"
    )
})

test_that("fineMappingPipeline refuses a phenotype QC its fits cannot take", {
    # Resolved before any study runs: a fit cannot take a missing phenotype,
    # so naAction = "keep" is refused up front, on either input form.
    keep <- PhenotypeFilterParam(naAction = "keep")
    expect_error(
        fineMappingPipeline(
            .fmp_qdA,
            methods = "susie",
            phenotypeFilterParam = keep
        ),
        "naAction = \"keep\" is not available"
    )
    expect_error(
        fineMappingPipeline(
            .fmp_msqd,
            methods = "susie",
            phenotypeFilterParam = keep
        ),
        "naAction = \"keep\" is not available"
    )
})

# ===========================================================================
# .fmAfByVar
# ===========================================================================

test_that(".fmAfByVar reads AF, and is NULL when af was never declared", {
    mk <- function(cm) {
        .dfToEntryGranges(suppressWarnings(.resolveSumstatCols(
            .afm_df("effect_allele_frequency", 0.8),
            cm,
            "lbl"
        )))
    }
    gA <- mk(.afm_cm(af = "effect_allele_frequency"))
    expect_equal(unname(.fmAfByVar(gA, "chr1:100:A:G")), 0.8)
    gM <- mk(.afm_cm(maf = "effect_allele_frequency"))
    expect_null(.fmAfByVar(gM, "chr1:100:A:G"))
})
