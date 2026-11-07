context("causalInferencePipeline")

# ===========================================================================
# Strategy: mock extractBlockGenotypes so .ldFromRef returns a real
# LD matrix on a small panel. Everything else (twasZ, MR, p-value combine)
# runs for real on the tiny fixture.
# ===========================================================================

.cip_makeHandle <- function(snp_n = 6L, n_samples = 30L, sample_prefix = "s") {
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

.cip_mockExtractor <- function(seed = 7, n_samples = 30L) {
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

.cip_makeGwasSumstats <- function(qc = TRUE) {
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
        Z = c(2.0, -1.5, 0.5, 1.2, -0.8),
        N = rep(1000L, 5),
        MAF = rep(0.3, 5)
    )
    GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .cip_makeHandle(),
        qcOptions = if (qc) list(step1 = "ok") else list()
    )
}

.cip_makeTwasWeights <- function(
    method = "susie",
    variant_ids = sprintf("chr1:%d:A:G", 100L * (1:5))
) {
    entry <- TwasWeightsRow(
        variantIds = variant_ids,
        weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
    )
    TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = method,
        entry = list(entry),
        ldReference = .cip_makeHandle()
    )
}

.cip_makeQtlFmr <- function(
    variant_ids = sprintf("chr1:%d:A:G", 100L * (1:5))
) {
    n <- length(variant_ids)
    tl <- data.frame(
        variant_id = variant_ids,
        pip = c(0.9, 0.05, 0.5, 0.8, 0.01),
        # posterior_mean / posterior_sd carry the "fine-mapped causal effect"
        # estimates that topLoci surfaces as beta / se in its projected
        # output (the column names downstream MR / TWAS code reads).
        posterior_mean = c(0.2, 0.05, -0.1, 0.3, 0.0),
        posterior_sd = rep(0.05, n),
        stringsAsFactors = FALSE
    )
    e <- FineMappingRow(
        variantIds = variant_ids,
        susieFit = list(),
        topLoci = tl
    )
    QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e),
        ldReference = .cip_makeHandle()
    )
}

# The (Q1, c1, t1) tuple of method `m`, as the pipeline's scorers take it.
.cip_tuple <- function(m) {
    list(qStudy = "Q1", qContext = "c1", qTrait = "t1", qMethod = m)
}

# ===========================================================================
# Input-type validation
# ===========================================================================

test_that("causalInferencePipeline: rejects non-GwasSumStats input", {
    expect_error(
        causalInferencePipeline(gwasSumStats = "no"),
        "Must inherit from class 'GwasSumStats'"
    )
})

test_that("causalInferencePipeline: rejects un-QCd GwasSumStats", {
    expect_error(
        causalInferencePipeline(
            gwasSumStats = .cip_makeGwasSumstats(qc = FALSE),
            twasWeights = .cip_makeTwasWeights()
        ),
        "has no QC record"
    )
})

test_that("causalInferencePipeline needs twasWeights or fineMappingResult", {
    expect_error(
        causalInferencePipeline(gwasSumStats = .cip_makeGwasSumstats()),
        "at least one of"
    )
})

test_that("causalInferencePipeline: rejects non-TwasWeights twasWeights arg", {
    expect_error(
        causalInferencePipeline(
            gwasSumStats = .cip_makeGwasSumstats(),
            twasWeights = "not a TwasWeights"
        ),
        "Must inherit from class 'TwasWeights'"
    )
})

test_that("causalInferencePipeline refuses a GWAS result in the QTL slot", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.1,
            stringsAsFactors = FALSE
        )
    )
    gfmr <- GwasFineMappingResult(
        studyName = "G1",
        method = "susie",
        entry = list(e)
    )
    expect_error(
        causalInferencePipeline(
            gwasSumStats = .cip_makeGwasSumstats(),
            fineMappingResult = gfmr
        ),
        "does not accept GWAS-side fine"
    )
})

# ===========================================================================
# .cipRequireMatchingLdRefs branches
# ===========================================================================

test_that(".cipRequireMatchingLdRefs: NULL twas-side ldReference is allowed", {
    # Build a TwasWeights without an ldReference.
    twNoLd <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(TwasWeightsRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
            weights = rep(0.1, 5)
        )),
        ldReference = NULL
    )
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    out <- causalInferencePipeline(
        gwasSumStats = .cip_makeGwasSumstats(),
        twasWeights = twNoLd
    )
    expect_s4_class(out, "GRanges")
})

test_that(".cipRequireMatchingLdRefs: sample set mismatch errors", {
    otherPanel <- .cip_makeHandle(sample_prefix = "other")
    twOther <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(TwasWeightsRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
            weights = rep(0.1, 5)
        )),
        ldReference = otherPanel
    )
    expect_error(
        causalInferencePipeline(
            gwasSumStats = .cip_makeGwasSumstats(),
            twasWeights = twOther
        ),
        "different sample sets"
    )
})

test_that(".cipRequireMatchingLdRefs: differently trimmed panels pass", {
    # Independent QC of the QTL and GWAS sides leaves one LD reference as two
    # overlapping-but-unequal panels, which the check reports rather than
    # refuses.
    tw <- .cip_makeTwasWeights()
    tw <- pecotmr:::.withRowLdRefs(tw, .cip_makeHandle(snp_n = 6L))
    gwas <- pecotmr:::.withRowLdRefs(
        .cip_makeGwasSumstats(),
        .cip_makeHandle(snp_n = 7L)
    )
    expect_warning(
        expect_null(pecotmr:::.cipRequireMatchingLdRefs(
            tw,
            gwas,
            label = "twasWeights"
        )),
        "share 6 variant"
    )
})

test_that(".cipRequireMatchingLdRefs: rows that never meet may differ", {
    # The check is per pair of rows scored together -- those whose regions
    # overlap -- so a weight row elsewhere on the genome is not held to the
    # GWAS rows' panel.
    elsewhere <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(TwasWeightsRow(
            variantIds = sprintf("chr2:%d:A:G", 100L * (1:5)),
            weights = rep(0.1, 5)
        )),
        ldReference = .cip_makeHandle(sample_prefix = "other")
    )
    expect_null(pecotmr:::.cipRequireMatchingLdRefs(
        elsewhere,
        .cip_makeGwasSumstats(),
        label = "twasWeights"
    ))
})

# ===========================================================================
# Happy path: TwasWeights only
# ===========================================================================

test_that("causalInferencePipeline: returns GRanges with TWAS Z per tuple", {
    tw <- .cip_makeTwasWeights()
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    out <- causalInferencePipeline(
        gwasSumStats = .cip_makeGwasSumstats(),
        twasWeights = tw
    )
    expect_s4_class(out, "GRanges")
    expect_equal(length(out), 1L)
    mc <- S4Vectors::mcols(out)
    expect_equal(as.character(mc$qtlStudy[[1L]]), "Q1")
    expect_equal(as.character(mc$gwasStudy[[1L]]), "G1")
    expect_true(is.finite(mc$twasZ[[1L]]))
    expect_true(is.finite(mc$twasPval[[1L]]))
    # No FMR -> MR fields stay NA.
    expect_true(is.na(mc$waldRatio[[1L]]))
    expect_true(is.na(mc$mrPval[[1L]]))
})

# ===========================================================================
# Happy path: TwasWeights + matching QtlFineMappingResult enables MR
# ===========================================================================

test_that("causalInferencePipeline fills MR fields given both inputs", {
    tw <- .cip_makeTwasWeights()
    fmr <- .cip_makeQtlFmr()
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    out <- causalInferencePipeline(
        gwasSumStats = .cip_makeGwasSumstats(),
        twasWeights = tw,
        fineMappingResult = fmr,
        mrParam = MrParam(pipCutoff = 0.5)
    )
    mc <- S4Vectors::mcols(out)
    expect_equal(as.character(mc$analysis), c("twas", "mr"))
    # MR uses PIP > 0.5 variants from the FMR (v1, v3, v4), on its own row.
    mr <- which(mc$analysis == "mr")
    expect_true(is.finite(mc$waldRatio[[mr]]))
    expect_true(is.finite(mc$mrPval[[mr]]))
    expect_gt(mc$nIV[[mr]], 0L)
    expect_true(is.na(mc$twasZ[[mr]]))
    expect_true(is.na(mc$waldRatio[[which(mc$analysis == "twas")]]))
})

# ===========================================================================
# MR comes from the fine-mapping result alone
#
# It used to attach to the TWAS row of the same method name, so it ran only
# for `susie`, vanished when weight selection kept another method, and never
# ran for a gene the TWAS weights lacked.
# ===========================================================================

# Weights for gene `trait` of methods `methods`, each with CV r-squared `rsq`.
.cip_methodsWeights <- function(methods, rsq, trait = "t1") {
    TwasWeights(
        studyName = rep("Q1", length(methods)),
        context = rep("c1", length(methods)),
        trait = rep(trait, length(methods)),
        method = methods,
        entry = map(rsq, .cip_rsqRow),
        ldReference = .cip_makeHandle()
    )
}

.cip_rsqRow <- function(rsq) {
    TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
        weights = c(0.2, 0.1, -0.1, 0.2, 0.1) * rsq,
        cvResult = list(metrics = c(rsq = rsq, pval = 0.001))
    )
}

# The rows of result `out` of analysis `a`, as a data frame.
.cip_rowsOf <- function(out, a) {
    mc <- as.data.frame(S4Vectors::mcols(out))
    mc[mc$analysis == a, ]
}

test_that("MR runs for every fine-mapping row, whatever TWAS chose", {
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    gw <- .cip_makeGwasSumstats()
    fmr <- .cip_makeQtlFmr()
    alone <- .cip_rowsOf(
        causalInferencePipeline(gw, fineMappingResult = fmr),
        "mr"
    )
    expect_equal(nrow(alone), 1L)
    # TWAS weights with no fine-mapping method, and ones whose best method is
    # not the fine-mapping one: the MR row is the fine-mapping result's.
    expectFmrMr <- function(tw) {
        out <- causalInferencePipeline(
            gw,
            twasWeights = tw,
            fineMappingResult = fmr
        )
        expect_equal(.cip_rowsOf(out, "twas")$method, "lasso")
        mr <- .cip_rowsOf(out, "mr")
        expect_equal(mr$method, "susie")
        expect_equal(mr$waldRatio, alone$waldRatio)
    }
    expectFmrMr(.cip_methodsWeights("lasso", 0.5))
    expectFmrMr(.cip_methodsWeights(c("susie", "lasso"), c(0.2, 0.5)))
    # A gene the TWAS weights do not hold still gets its MR.
    other <- causalInferencePipeline(
        gw,
        twasWeights = .cip_methodsWeights("lasso", 0.5, trait = "t9"),
        fineMappingResult = fmr
    )
    expect_equal(.cip_rowsOf(other, "twas")$trait, "t9")
    expect_equal(.cip_rowsOf(other, "mr")$trait, "t1")
})

test_that("MR never reads a fine-mapping method among the TWAS weights", {
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    gw <- .cip_makeGwasSumstats()
    fmr <- .cip_makeQtlFmr()
    # TWAS susie weights unlike the fine-mapping result's.
    tw <- .cip_methodsWeights("susie", 0.5)
    both <- causalInferencePipeline(
        gw,
        twasWeights = tw,
        fineMappingResult = fmr
    )
    alone <- causalInferencePipeline(gw, fineMappingResult = fmr)
    expect_equal(
        .cip_rowsOf(both, "mr")$waldRatio,
        .cip_rowsOf(alone, "mr")$waldRatio
    )
    # Without a fine-mapping result, those TWAS weights give no MR at all.
    expect_equal(
        nrow(.cip_rowsOf(causalInferencePipeline(gw, twasWeights = tw), "mr")),
        0L
    )
})

test_that("all, or named methods, score those methods; best one per gene", {
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    gw <- .cip_makeGwasSumstats()
    tw <- .cip_methodsWeights(c("susie", "lasso", "enet"), c(0.2, 0.5, 0.3))
    methodsFor <- function(choice) {
        .cip_rowsOf(
            causalInferencePipeline(
                gw,
                twasWeights = tw,
                weightSelectionParam = choice
            ),
            "twas"
        )$method
    }
    expect_equal(methodsFor(WeightSelectionParam()), "lasso")
    expect_equal(
        methodsFor(WeightSelectionParam(methods = "all")),
        c("susie", "lasso", "enet")
    )
    expect_equal(
        methodsFor(WeightSelectionParam(methods = c("enet", "susie"))),
        c("susie", "enet")
    )
    # The cutoffs bound every choice.
    expect_equal(
        methodsFor(WeightSelectionParam(methods = "all", cutoff = 0.25)),
        c("lasso", "enet")
    )
    expect_error(
        methodsFor(WeightSelectionParam(methods = "mrash")),
        "names mrash, which the weights do not hold"
    )
})

test_that("an MR row is reported even when no TWAS row could be scored", {
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    # Weights on one variant: too little overlap for a TWAS Z.
    tw <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "lasso",
        entry = list(TwasWeightsRow(
            variantIds = "chr1:100:A:G",
            weights = 0.3
        )),
        ldReference = .cip_makeHandle()
    )
    out <- suppressWarnings(causalInferencePipeline(
        .cip_makeGwasSumstats(),
        twasWeights = tw,
        fineMappingResult = .cip_makeQtlFmr(),
        mrParam = MrParam(pipCutoff = 0.5, pvalCutoff = 0.01)
    ))
    expect_equal(as.character(out$analysis), "mr")
    # With no TWAS p-value to read, the gate does not withhold it.
    expect_true(is.finite(out$waldRatio))
})

test_that("the MR gate reads the gene's reported TWAS p-values", {
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    gw <- .cip_makeGwasSumstats()
    fmr <- .cip_makeQtlFmr()
    tw <- .cip_methodsWeights("lasso", 0.5)
    twasP <- .cip_rowsOf(
        causalInferencePipeline(gw, twasWeights = tw),
        "twas"
    )$twasPval
    mrAt <- function(cutoff, weights = tw) {
        .cip_rowsOf(
            causalInferencePipeline(
                gw,
                twasWeights = weights,
                fineMappingResult = fmr,
                mrParam = MrParam(pvalCutoff = cutoff)
            ),
            "mr"
        )$waldRatio
    }
    # Withheld when the TWAS p is above the cutoff, the row still reported.
    expect_true(is.na(mrAt(twasP / 2)))
    expect_true(is.finite(mrAt(twasP * 2)))
    # A gene with no TWAS row is not gated.
    expect_true(is.finite(mrAt(1e-12, .cip_methodsWeights("lasso", 0.5, "t9"))))
})

# ===========================================================================
# Fine-mapping-only path: weights come from topLoci$betahat
# ===========================================================================

test_that("FMR-only causal inference takes weights from topLoci$betahat", {
    fmr <- .cip_makeQtlFmr()
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    out <- causalInferencePipeline(
        gwasSumStats = .cip_makeGwasSumstats(),
        fineMappingResult = fmr
    )
    mc <- S4Vectors::mcols(out)
    expect_true(is.finite(mc$twasZ[[1L]]))
    # The fine-mapping result gives both the TWAS and the MR row.
    expect_equal(as.character(mc$analysis), c("twas", "mr"))
})

# ===========================================================================
# combineMethods integration
# ===========================================================================

test_that("causalInferencePipeline: combineMethods appends combined rows", {
    tw1 <- TwasWeights(
        studyName = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("lasso", "enet"),
        entry = list(
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
                weights = rep(0.1, 5)
            ),
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
                weights = rep(0.05, 5)
            )
        ),
        ldReference = .cip_makeHandle()
    )
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    out <- causalInferencePipeline(
        gwasSumStats = .cip_makeGwasSumstats(),
        twasWeights = tw1,
        combineMethods = "acat"
    )
    # 2 per-tuple rows + 1 combined row = 3.
    expect_equal(length(out), 3L)
    mc <- S4Vectors::mcols(out)
    expect_true(any(grepl("^combined\\.", as.character(mc$method))))
})

# ===========================================================================
# .cipZToBeta / .cipZToSe fallbacks
# ===========================================================================

test_that(".cipZToBeta: falls back to z when maf/n are NA", {
    res <- pecotmr:::.cipZToBeta(z = c(1, 2), maf = NA, n = NA)
    expect_equal(res, c(1, 2))
})

test_that(".cipZToSe: falls back to vector of 1 when maf/n are NA", {
    res <- pecotmr:::.cipZToSe(z = c(1, 2), maf = NA, n = NA)
    expect_equal(res, c(1, 1))
})

test_that("the cutoffs bound every choice; a gene with no CV keeps all", {
    mkEntry <- function(rsq, pval = 0.01) {
        TwasWeightsRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
            weights = rep(0.1, 3),
            cvResult = list(metrics = c(corr = 0.1, rsq = rsq, pval = pval))
        )
    }
    tw <- TwasWeights(
        studyName = rep("S", 4),
        context = rep("c1", 4),
        trait = c("G", "G", "G", "G2"),
        method = c("susie", "enet", "lasso", "susie"),
        entry = list(mkEntry(0.20), mkEntry(0.05), mkEntry(0.50), mkEntry(0.01))
    )
    keep <- function(x, ...) {
        choice <- WeightSelectionParam(...)
        pecotmr:::.selectWeightCells(
            pecotmr:::.weightCells(x, choice),
            choice,
            "f()"
        )
    }
    # rsq gate: G keeps susie (.20) and lasso (.50); enet (.05) and G2 (.01)
    # are out.
    f1 <- keep(tw, methods = "all", cutoff = 0.1)
    expect_equal(f1$method, c("susie", "lasso"))
    expect_true(all(f1$trait == "G"))
    expect_equal(nrow(keep(tw, methods = "all", cutoff = 0.99)), 0L)
    # Within the gates "best" takes each gene's highest rsq.
    expect_equal(keep(tw, cutoff = 0.1)$method, "lasso")
    best <- keep(tw)
    expect_equal(paste(best$trait, best$method), c("G lasso", "G2 susie"))
    # Named methods, within the gates too.
    expect_equal(
        keep(tw, methods = c("enet", "susie"), cutoff = 0.1)$method,
        "susie"
    )
    # pval gate: a high-rsq method with a bad CV p-value is excluded.
    tw2 <- TwasWeights(
        studyName = rep("S", 2),
        context = rep("c1", 2),
        trait = rep("G", 2),
        method = c("susie", "lasso"),
        entry = list(mkEntry(0.50, pval = 0.20), mkEntry(0.20, pval = 0.01))
    )
    expect_equal(keep(tw2, pvalCutoff = 0.05)$method, "lasso")
    # No usable cvResult: nothing to gate or rank by, so every method stays.
    twss <- TwasWeights(
        studyName = rep("S", 2),
        context = rep("c1", 2),
        trait = rep("G", 2),
        method = c("susie", "lasso"),
        entry = list(
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = rep(0.1, 3)
            ),
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = rep(0.1, 3)
            )
        )
    )
    expect_equal(nrow(keep(twss, cutoff = 0.1)), 2L)
    expect_equal(nrow(keep(twss)), 2L)
})

test_that("best passes over a method whose weights are all zero", {
    cells <- tibble(
        study = "S",
        context = "c1",
        trait = "G",
        method = c("susie", "lasso"),
        rsq = c(0.5, 0.2),
        pval = 0.01,
        zero = c(TRUE, FALSE)
    )
    best <- WeightSelectionParam()
    expect_equal(
        pecotmr:::.selectWeightCells(cells, best, "f()")$method,
        "lasso"
    )
    # When every method's weights are zero, the highest rsq is kept anyway.
    allZero <- mutate(cells, zero = TRUE)
    expect_equal(
        pecotmr:::.selectWeightCells(allZero, best, "f()")$method,
        "susie"
    )
    # A tie goes to the first.
    tied <- mutate(cells, rsq = 0.3, zero = FALSE)
    expect_equal(
        pecotmr:::.selectWeightCells(tied, best, "f()")$method,
        "susie"
    )
})

test_that("causal inference picks each group's max-rsq method (rsqCutoff)", {
    tw <- TwasWeights(
        studyName = rep("Q1", 2),
        context = rep("c1", 2),
        trait = rep("t1", 2),
        method = c("susie", "lasso"),
        entry = list(
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
                weights = c(0.1, 0.05, -0.2, 0.3, 0.0),
                cvResult = list(metrics = c(rsq = 0.2, pval = 0.001))
            ),
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
                weights = c(0.2, 0.1, -0.1, 0.2, 0.1),
                cvResult = list(metrics = c(rsq = 0.5, pval = 0.001))
            )
        ),
        ldReference = .cip_makeHandle()
    )
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    out <- causalInferencePipeline(
        gwasSumStats = .cip_makeGwasSumstats(),
        twasWeights = tw,
        weightSelectionParam = WeightSelectionParam(cutoff = 0.1)
    )
    expect_equal(as.character(S4Vectors::mcols(out)$method), "lasso")
})

test_that("best passes over all-zero weights, whose TWAS Z would be NaN", {
    tw <- TwasWeights(
        studyName = rep("Q1", 2),
        context = rep("c1", 2),
        trait = rep("t1", 2),
        method = c("susie", "lasso"),
        entry = list(
            # top rsq but all-zero weights: not a method to score
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
                weights = rep(0, 5),
                cvResult = list(metrics = c(rsq = 0.9, pval = 0.001))
            ),
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
                weights = c(0.2, 0.1, -0.1, 0.2, 0.1),
                cvResult = list(metrics = c(rsq = 0.5, pval = 0.001))
            )
        ),
        ldReference = .cip_makeHandle()
    )
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    out <- causalInferencePipeline(
        gwasSumStats = .cip_makeGwasSumstats(),
        twasWeights = tw,
        weightSelectionParam = WeightSelectionParam(cutoff = 0.1)
    )
    expect_equal(as.character(S4Vectors::mcols(out)$method), "lasso")
    expect_true(is.finite(S4Vectors::mcols(out)$twasZ))
})

test_that("causalInferencePipeline: rsqPvalCutoff drops high-CV-pval methods", {
    tw <- TwasWeights(
        studyName = rep("Q1", 2),
        context = rep("c1", 2),
        trait = rep("t1", 2),
        method = c("susie", "lasso"),
        entry = list(
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
                weights = c(0.2, 0.1, -0.1, 0.2, 0.1),
                cvResult = list(metrics = c(rsq = 0.9, pval = 0.5))
            ),
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
                weights = c(0.2, 0.1, -0.1, 0.2, 0.1),
                cvResult = list(metrics = c(rsq = 0.5, pval = 0.001))
            )
        ),
        ldReference = .cip_makeHandle()
    )
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    out <- causalInferencePipeline(
        gwasSumStats = .cip_makeGwasSumstats(),
        twasWeights = tw,
        weightSelectionParam = WeightSelectionParam(
            cutoff = 0.1,
            pvalCutoff = 0.05
        )
    )
    expect_equal(as.character(S4Vectors::mcols(out)$method), "lasso")
})


context("twasZ")

# ===========================================================================
# Fixture builder: a small genotype + LD matrix.
# ===========================================================================

.tz_makeLd <- function(n = 100, p = 8, seed = 7) {
    set.seed(seed)
    X <- matrix(rbinom(n * p, 2, runif(p, 0.2, 0.8)), nrow = n, ncol = p)
    vid <- sprintf("chr1:%d:A:G", 100L * seq_len(p))
    colnames(X) <- vid
    af <- colMeans(X) / 2
    Xstd <- sweep(X, 2, 2 * af)
    Xstd <- sweep(Xstd, 2, sqrt(2 * af * (1 - af)), "/")
    R <- crossprod(Xstd) / (n - 1)
    rownames(R) <- colnames(R) <- vid
    list(X = X, Xstd = Xstd, R = R, vid = vid, n = n, p = p)
}

# ===========================================================================
# twasZ: input coercion and validity
# ===========================================================================

test_that("twasZ: numeric vector input is coerced to a single-method matrix", {
    d <- .tz_makeLd()
    w <- rnorm(d$p)
    names(w) <- d$vid
    z <- rnorm(d$p)
    res <- pecotmr:::twasZ(w, z, R = d$R)
    expect_equal(dim(res$Z), c(1L, 2L))
    expect_equal(colnames(res$Z), c("Z", "pval"))
    expect_equal(rownames(res$Z), "method1")
})

test_that("twasZ: matrix without colnames gets method1/method2 names", {
    d <- .tz_makeLd()
    W <- cbind(rnorm(d$p), rnorm(d$p))
    rownames(W) <- d$vid
    z <- rnorm(d$p)
    res <- pecotmr:::twasZ(W, z, R = d$R)
    expect_equal(rownames(res$Z), c("method1", "method2"))
})

test_that("twasZ: non-matrix non-numeric weights errors out", {
    expect_error(
        pecotmr:::twasZ(list(a = 1), z = 1, R = matrix(1)),
        "must be a numeric vector or a matrix"
    )
})

test_that("twasZ: length mismatch between weights and z errors", {
    expect_error(
        pecotmr:::twasZ(c(0.1, 0.2, 0.3), z = c(1, 2)),
        "weights.*Must have exactly 2 rows"
    )
})

test_that("twasZ: missing R/X/V triplet errors", {
    expect_error(
        pecotmr:::twasZ(c(0.1, 0.2), z = c(1, 2)),
        "provide R, X, or the .V, D, nSketch. SVD triplet"
    )
})

# ===========================================================================
# twasZ: R path
# ===========================================================================

test_that("twasZ: R path with named alignment matches the manual formula", {
    d <- .tz_makeLd()
    set.seed(1)
    z <- rnorm(d$p)
    names(z) <- d$vid
    w <- rnorm(d$p)
    names(w) <- d$vid
    res <- pecotmr:::twasZ(w, z, R = d$R)
    expected_stat <- sum(w * z)
    expected_denom <- sqrt(as.numeric(crossprod(w, d$R %*% w)))
    expect_equal(
        as.numeric(res$Z[, "Z"]),
        expected_stat / expected_denom,
        tolerance = 1e-10
    )
})

test_that("twasZ: R path realigns by rowname order", {
    d <- .tz_makeLd()
    z <- rnorm(d$p)
    w <- rnorm(d$p)
    names(w) <- d$vid
    # Permute R rows/columns; results must be identical to the un-permuted R.
    perm <- sample(seq_len(d$p))
    R_perm <- d$R[perm, perm]
    res_orig <- pecotmr:::twasZ(w, z, R = d$R)
    res_perm <- pecotmr:::twasZ(w, z, R = R_perm)
    expect_equal(res_perm$Z, res_orig$Z, tolerance = 1e-10)
})

test_that("twasZ: R path errors when R is missing rows named in weights", {
    d <- .tz_makeLd()
    w <- rnorm(d$p)
    names(w) <- d$vid
    R_short <- d$R[seq_len(d$p - 1), seq_len(d$p - 1)]
    expect_error(
        pecotmr:::twasZ(w, rnorm(d$p), R = R_short),
        "R is missing rows for"
    )
})

test_that("twasZ: R path positional alignment errors on dim mismatch", {
    d <- .tz_makeLd()
    w <- rnorm(d$p) # unnamed -> positional alignment
    R_short <- unname(d$R[seq_len(d$p - 1), seq_len(d$p - 1)])
    expect_error(
        pecotmr:::twasZ(w, rnorm(d$p), R = R_short),
        "R \\(positional alignment\\).*Must have exactly 8 rows"
    )
})

# ===========================================================================
# twasZ: X path (computes R via computeLd)
# ===========================================================================

test_that("twasZ: X path matches R from computeLd without forming it", {
    d <- .tz_makeLd()
    set.seed(2)
    w <- rnorm(d$p)
    names(w) <- d$vid
    z <- rnorm(d$p)
    # The X path computes W'RW from X's standardized dosages, R being
    # computeLd(X)'s cor(X). Compare against that matrix passed explicitly.
    R_from_X <- cor(d$X)
    rownames(R_from_X) <- colnames(R_from_X) <- d$vid
    res_X <- pecotmr:::twasZ(w, z, X = d$X)
    res_Rfx <- pecotmr:::twasZ(w, z, R = R_from_X)
    expect_equal(res_X$Z, res_Rfx$Z, tolerance = 1e-10)
})

test_that("twasZ: X path uses computeLd's LD for missing, constant columns", {
    d <- .tz_makeLd()
    X <- d$X
    X[, 3] <- 1L
    X[c(4, 9), 5] <- NA
    set.seed(4)
    w <- setNames(rnorm(d$p), d$vid)
    z <- rnorm(d$p)
    expect_equal(
        pecotmr:::twasZ(w, z, X = X)$Z,
        pecotmr:::twasZ(w, z, R = computeLd(X))$Z
    )
})

test_that("twasZ: X path aligns columns to the weights by name", {
    d <- .tz_makeLd()
    set.seed(5)
    w <- setNames(rnorm(d$p), d$vid)
    z <- rnorm(d$p)
    shuffled <- d$X[, rev(seq_len(d$p))]
    expect_equal(
        pecotmr:::twasZ(w, z, X = shuffled)$Z,
        pecotmr:::twasZ(w, z, X = d$X)$Z
    )
    expect_error(
        pecotmr:::twasZ(w, z, X = d$X[, -1]),
        "X is missing columns for 1 variant"
    )
})

# ===========================================================================
# twasZ: SVD path coverage (V-row missing error)
# ===========================================================================

test_that("twasZ: SVD path errors when V is missing rows named in weights", {
    d <- .tz_makeLd()
    s <- svd(d$Xstd)
    rownames(s$v) <- d$vid[seq_len(nrow(s$v))]
    w <- rnorm(d$p)
    names(w) <- c(d$vid[seq_len(d$p - 1)], "ghost")
    expect_error(
        pecotmr:::twasZ(w, rnorm(d$p), V = s$v, D = s$d, nSketch = d$n),
        "V is missing rows for"
    )
})

test_that("twasZ: SVD path positional alignment errors on dim mismatch", {
    d <- .tz_makeLd()
    s <- svd(d$Xstd)
    # Pass an unnamed V with the wrong nrow.
    w <- rnorm(d$p)
    V_short <- s$v[seq_len(d$p - 1), , drop = FALSE]
    expect_error(
        pecotmr:::twasZ(w, rnorm(d$p), V = V_short, D = s$d, nSketch = d$n),
        "V \\(positional alignment\\).*Must have exactly 8 rows"
    )
})

# ===========================================================================
# twasZ: combineMethods integration
# ===========================================================================

test_that("twasZ: combineMethods K=1 returns the per-tuple p-value unchanged", {
    d <- .tz_makeLd()
    w <- rnorm(d$p)
    names(w) <- d$vid
    z <- rnorm(d$p)
    res <- pecotmr:::twasZ(
        w,
        z,
        R = d$R,
        combineMethods = c("acat", "bonferroni")
    )
    expect_false(is.null(res$combined))
    expect_equal(res$combined$results$acat$pval, res$Z[1, "pval"])
    expect_equal(res$combined$results$bonferroni$pval, res$Z[1, "pval"])
    expect_equal(res$combined$input$nValid, 1L)
})

test_that("twasZ passes K>=2 methods to combinePValues with correlation", {
    d <- .tz_makeLd()
    set.seed(11)
    W <- cbind(method_a = rnorm(d$p), method_b = rnorm(d$p))
    rownames(W) <- d$vid
    z <- rnorm(d$p)
    res <- pecotmr:::twasZ(W, z, R = d$R, combineMethods = "acat")
    expect_false(is.null(res$combined))
    # ACAT combines two p-values; result must be in (0, 1).
    expect_true(
        res$combined$results$acat$pval >= 0 &&
            res$combined$results$acat$pval <= 1
    )
    expect_equal(res$combined$input$nValid, 2L)
})
# ===========================================================================
# twasZ: hand-computed values on the R path (positional alignment)
# ===========================================================================

test_that("twasZ: R path gives w'z / sqrt(w'Rw) on hand-computed cases", {
    # Each case: weights, z, R and the Z worked out by hand.
    cases <- list(
        # 0.7 * 2.5 / sqrt(0.7^2) = 2.5
        singleVariant = list(0.7, 2.5, matrix(1), 2.5),
        # 5 * 3 / sqrt(5)
        equalWeights = list(rep(1, 5), rep(3, 5), diag(5), 15 / sqrt(5)),
        # (-1 - 0.3) / sqrt(0.25 + 0.09): the sign follows the weights.
        negativeWeights = list(
            c(-0.5, -0.3),
            c(2, 1),
            diag(2),
            -1.3 / sqrt(0.34)
        ),
        # w'Rw = sum(w)^2 = 1 when every pair is perfectly correlated.
        perfectlyCorrelated = list(rep(0.25, 4), rep(2, 4), matrix(1, 4, 4), 2),
        # A single non-zero weight reads that variant's z alone.
        sparseWeights = list(c(0, 0, 1, 0, 0), 1:5, diag(5), 3),
        # 1e6 * 1e-6 / sqrt(1e12)
        largeWeightTinyZ = list(c(1e6, 0, 0), c(1e-6, 5, 5), diag(3), 1e-6),
        # Z = 0, so the p-value is 1.
        zeroZ = list(c(0.5, 0.3, 0.2), c(0, 0, 0), diag(3), 0),
        # Z = 2e6 / sqrt(2), so the p-value underflows to 0.
        largeZ = list(c(1, 1), c(1e6, 1e6), diag(2), 2e6 / sqrt(2))
    )
    for (nm in names(cases)) {
        cs <- cases[[nm]]
        res <- twasZ(cs[[1L]], cs[[2L]], R = cs[[3L]])$Z
        expect_equal(colnames(res), c("Z", "pval"), info = nm)
        expect_equal(unname(res[, "Z"]), cs[[4L]], tolerance = 1e-10, info = nm)
        expect_equal(
            unname(res[, "pval"]),
            2 * pnorm(-abs(cs[[4L]])),
            tolerance = 1e-10,
            info = nm
        )
    }
})

test_that("twasZ: all-zero weights produce NaN z-score", {
    weights <- c(0, 0, 0)
    z <- c(1.5, -0.5, 2.0)
    R <- diag(3)
    result <- twasZ(weights, z, R = R)
    # stat = 0, denom = 0, so zscore = 0/0 = NaN
    expect_true(is.nan(as.numeric(result$Z[, "Z"])))
})


test_that("twasZ: computing R from X matches providing R directly", {
    set.seed(123)
    n <- 20
    p <- 5
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- paste0("SNP", seq_len(p))
    R <- cor(X)
    weights <- rnorm(p)
    z <- rnorm(p)
    result_with_R <- twasZ(weights, z, R = R)
    result_with_X <- twasZ(weights, z, X = X)
    expect_equal(
        as.numeric(result_with_R$Z[, "Z"]),
        as.numeric(result_with_X$Z[, "Z"]),
        tolerance = 1e-6
    )
    expect_equal(
        as.numeric(result_with_R$Z[, "pval"]),
        as.numeric(result_with_X$Z[, "pval"]),
        tolerance = 1e-6
    )
})

# ===========================================================================
# twasZ: matrix weights path -- multiple methods/conditions
# ===========================================================================

test_that("twasZ: matrix weights produce one Z row per column", {
    set.seed(10)
    p <- 5
    k <- 3
    weights <- matrix(rnorm(p * k), nrow = p, ncol = k)
    rownames(weights) <- paste0("SNP", seq_len(p))
    colnames(weights) <- paste0("Cond", seq_len(k))
    z <- rnorm(p)
    R <- diag(p)
    rownames(R) <- colnames(R) <- paste0("SNP", seq_len(p))
    result <- twasZ(weights, z, R = R)
    expect_equal(nrow(result$Z), k)
    expect_equal(rownames(result$Z), paste0("Cond", seq_len(k)))
    expect_equal(colnames(result$Z), c("Z", "pval"))
    # combineMethods omitted -> combined is NULL
    expect_null(result$combined)
})


# ===========================================================================
# twasZ: the SVD path agrees with the R path
# ===========================================================================

test_that("twasZ: SVD path matches R path, full rank and n < p", {
    # Sketch sizes (samples, variants, seed): full rank, and fewer samples
    # than variants.
    shapes <- list(fullRank = c(100, 20, 42), rankDeficient = c(15, 30, 123))
    for (nm in names(shapes)) {
        n <- shapes[[nm]][[1L]]
        p <- shapes[[nm]][[2L]]
        set.seed(shapes[[nm]][[3L]])
        X <- matrix(rbinom(n * p, 2, runif(p, 0.15, 0.85)), nrow = n, ncol = p)
        af <- colMeans(X) / 2
        # A monomorphic column has no standardized dosage.
        keep <- af > 0 & af < 1
        X <- X[, keep, drop = FALSE]
        af <- af[keep]
        Xstd <- sweep(X, 2, 2 * af)
        Xstd <- sweep(Xstd, 2, sqrt(2 * af * (1 - af)), "/")
        R <- crossprod(Xstd) / (n - 1)
        s <- svd(Xstd)
        weights <- rnorm(ncol(X))
        z <- rnorm(ncol(X))
        expect_equal(
            twasZ(weights, z, V = s$v, D = s$d, nSketch = n)$Z,
            twasZ(weights, z, R = R)$Z,
            tolerance = 1e-10,
            info = nm
        )
    }
})


# ===========================================================================
# Direct unit tests for the MR / metric helpers (mock topLoci to feed
# controlled topLoci frames; gwasDf is a plain data.frame as produced by
# as.data.frame upstream).
# ===========================================================================

.cip_gwasDf <- function(
    vids = sprintf("chr1:%d:A:G", 100L * (1:5)),
    z = rep_len(c(2, -1.5, 1.8, 0.5, -2.2), length(vids))
) {
    data.frame(
        variant_id = vids,
        z = z,
        N = rep(1000, length(vids)),
        maf = rep(0.3, length(vids)),
        chrom = "chr1",
        pos = seq(100L, by = 100L, length.out = length(vids)),
        stringsAsFactors = FALSE
    )
}

test_that(".cipCalcI2: degenerate Q / single group -> 0; normal Q -> clamped", {
    expect_equal(pecotmr:::.cipCalcI2(0, 3L), 0) # Q <= 1e-3
    expect_equal(pecotmr:::.cipCalcI2(NA_real_, 3L), 0) # non-finite Q
    expect_equal(pecotmr:::.cipCalcI2(10, 1L), 0) # nGroups <= 1
    expect_equal(pecotmr:::.cipCalcI2(10, 3L), (10 - 2) / 10)
    expect_equal(pecotmr:::.cipCalcI2(1.5, 5L), 0) # negative -> clamp to 0
})

test_that(".cipComputeMr: IVW Wald-ratio over PIP-passing instruments", {
    tl <- data.frame(
        variant_id = sprintf("chr1:%d:A:G", 100L * (1:4)),
        pip = c(0.9, 0.8, 0.1, 0.7),
        beta = c(0.3, -0.2, 0.5, 0.25),
        se = rep(0.05, 4),
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(
        .cipTopLociTable = function(x) tl,
        .package = "pecotmr"
    )
    res <- pecotmr:::.cipComputeMr(
        NULL,
        .cip_gwasDf(sprintf("chr1:%d:A:G", 100L * (1:4)), c(2, -1.5, 1, 1.8)),
        pipCutoff = 0.5
    )
    expect_true(is.finite(res$waldRatio))
    expect_true(is.finite(res$mrPval))
    expect_equal(res$nIV, 3L) # v1, v2, v4 pass + overlap
})

# ---------------------------------------------------------------------------
# Exposure / outcome scale agreement. Without maf + N the outcome falls back
# to z with unit se; the exposure has to follow, or the Wald ratio is a ratio
# of incommensurable units, inflated by 1/seX (answers like -1094 +/- 1105).
# ---------------------------------------------------------------------------

test_that(".cipComputeMr: the exposure standardizes when maf/N are absent", {
    tl <- data.frame(
        variant_id = "chr1:100:A:G",
        pip = 0.9,
        beta = 0.3,
        se = 0.01,
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(
        .cipTopLociTable = function(x) tl,
        .package = "pecotmr"
    )
    g <- .cip_gwasDf("chr1:100:A:G", 2)
    noScale <- g[, setdiff(names(g), c("N", "maf"))]
    res <- pecotmr:::.cipComputeMr(NULL, noScale, pipCutoff = 0.5)
    # A single instrument, so IVW hands back its own ratio: z / (beta / se).
    # Before the fix this was z / beta -- 100x larger for this instrument.
    expect_equal(res$nIV, 1L)
    expect_equal(res$waldRatio, 2 / (0.3 / 0.01))
    # With maf and N the outcome is a real effect size, so the raw exposure
    # stands and the ratio stays interpretable.
    keep <- pecotmr:::.cipComputeMr(NULL, g, pipCutoff = 0.5)
    expect_equal(
        keep$waldRatio,
        pecotmr:::.cipZToBeta(2, 0.3, 1000) / 0.3
    )
})

test_that(".cipCsAwareAlign: exposure scale follows the GWAS, both ways", {
    raw <- list(
        cs = c(1L, 1L),
        pip = c(0.6, 0.3),
        bhatX = c(0.3, -0.2),
        sbhatX = c(0.05, 0.04),
        vids = sprintf("chr1:%d:A:G", c(100L, 200L))
    )
    g <- .cip_gwasDf(raw$vids, c(2, -1.5))
    al <- pecotmr:::.cipCsAwareAlign(raw, g, alleleFlip = TRUE)
    # maf + N present: both sides are real effect sizes, so the exposure
    # keeps its own SE rather than being forced to 1. This is the half that
    # was wrong before -- csAware standardized unconditionally.
    expect_equal(al$bhatX, raw$bhatX)
    expect_equal(al$sbhatX, raw$sbhatX)
    # Without them the outcome is z / 1, so the exposure standardizes.
    noScale <- g[, setdiff(names(g), c("N", "maf"))]
    alNo <- pecotmr:::.cipCsAwareAlign(raw, noScale, alleleFlip = TRUE)
    expect_equal(alNo$bhatX, raw$bhatX / raw$sbhatX)
    expect_equal(alNo$sbhatX, rep(1, 2))
})

# ---------------------------------------------------------------------------
# Phase 3: tuple-based (chrom/pos/allele) matching of QTL exposure to GWAS.
# Positional ids exercise matchVariants' tuple path (the v1.. fixtures above
# use the rsID string-fallback). Previously these joins used intersect()/match
# on the raw id string, so a chr-prefix or allele-swap difference silently
# dropped the variant.
# ---------------------------------------------------------------------------

test_that(".cipComputeMr matches across a chr-prefix difference", {
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:A:G"),
        pip = c(0.9, 0.9),
        beta = c(0.3, 0.2),
        se = c(0.05, 0.05),
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(
        .cipTopLociTable = function(x) tl,
        .package = "pecotmr"
    )
    # GWAS carries the same variants without the "chr" prefix.
    g <- .cip_gwasDf(c("1:100:A:G", "1:200:A:G"), z = c(2, 1.5))
    res <- pecotmr:::.cipComputeMr(NULL, g, pipCutoff = 0.5)
    expect_equal(res$nIV, 2L) # matched despite chr-prefix mismatch
    expect_true(is.finite(res$waldRatio))
})

test_that(".cipComputeMr keeps an allele-swapped GWAS variant, sign flipped", {
    tl <- data.frame(
        variant_id = "chr1:100:A:G",
        pip = 0.9,
        beta = 0.3,
        se = 0.05,
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(
        .cipTopLociTable = function(x) tl,
        .package = "pecotmr"
    )
    rSame <- pecotmr:::.cipComputeMr(
        NULL,
        .cip_gwasDf("chr1:100:A:G", z = 2.0),
        pipCutoff = 0.5
    )
    # Same locus, ref/alt swapped on the GWAS side (same |z|): previously
    # dropped (id-string mismatch); now matched with sign -1, negating the Wald
    # ratio.
    rSwap <- pecotmr:::.cipComputeMr(
        NULL,
        .cip_gwasDf("chr1:100:G:A", z = 2.0),
        pipCutoff = 0.5
    )
    expect_equal(rSwap$nIV, 1L)
    expect_equal(rSwap$waldRatio, -rSame$waldRatio)
})

test_that(".cipComputeMr(alleleFlip = FALSE) drops a swapped GWAS variant", {
    tl <- data.frame(
        variant_id = "chr1:100:A:G",
        pip = 0.9,
        beta = 0.3,
        se = 0.05,
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(
        .cipTopLociTable = function(x) tl,
        .package = "pecotmr"
    )
    # GWAS variant is the ref/alt swap; with flipping disabled it is not
    # matched.
    rSwap <- pecotmr:::.cipComputeMr(
        NULL,
        .cip_gwasDf("chr1:100:G:A", z = 2.0),
        pipCutoff = 0.5,
        alleleFlip = FALSE
    )
    expect_equal(rSwap$nIV, 0L)
    expect_true(is.na(rSwap$waldRatio))
})

test_that(".cipComputeMrCsAware matches across a chr-prefix difference", {
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:A:G"),
        cs = c(1L, 1L),
        pip = c(0.6, 0.4),
        beta = c(0.3, 0.2),
        se = c(0.05, 0.05),
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(
        .cipTopLociTable = function(x) tl,
        .package = "pecotmr"
    )
    g <- .cip_gwasDf(c("1:100:A:G", "1:200:A:G"), z = c(2, 1.5))
    res <- pecotmr:::.cipComputeMrCsAware(NULL, g, cpipCutoff = 0.5)
    expect_equal(res$nCs, 1L)
    expect_true(is.finite(res$waldRatio))
})

test_that(".cipComputeTwasZ gives a Z from a panel's factor or an LdMatrix", {
    data(gwasSumStatsS4Example)
    panel <- ldReference(gwasSumStatsS4Example)[[1L]]
    ids <- names(SummarizedExperiment::rowRanges(panel))[c(3, 8, 15, 21, 30)]
    ld <- LdMatrix(
        unname(computeLd(panel)),
        variants = SummarizedExperiment::rowRanges(panel)
    )
    g <- .cip_gwasDf(ids, z = c(2, -1.5, 0.8, 1.1, -0.3))
    w <- c(0.2, -0.1, 0.3, 0.05, -0.2)
    viaPanel <- pecotmr:::.cipComputeTwasZ(w, ids, g, panel)
    viaMatrix <- pecotmr:::.cipComputeTwasZ(w, ids, g, ld)
    expect_true(is.finite(viaPanel$Z))
    expect_equal(viaPanel, viaMatrix)
})

test_that(".cipComputeTwasZ reconciles a chr-prefix mismatch", {
    idsGwas <- c("1:100:A:G", "1:200:A:G", "1:300:A:G") # GWAS / panel coding
    idsWt <- c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G") # QTL weight ids
    h <- .cip_makeHandle(snp_n = 3L)
    h@snpInfo$SNP <- idsGwas
    h@snpInfo$A2 <- rep("A", 3)
    h@snpInfo$A1 <- rep("G", 3)
    g <- .cip_gwasDf(idsGwas, z = c(2, -1.5, 0.8))
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    res <- pecotmr:::.cipComputeTwasZ(c(0.2, -0.1, 0.3), idsWt, g, h)
    expect_false(is.null(res)) # intersect() would have returned NULL
    expect_true(is.finite(res$Z))
})

test_that(".cipComputeMr: NA on empty / no-IV / no-overlap / zero-beta", {
    g <- .cip_gwasDf(sprintf("chr1:%d:A:G", 100L * (1:4)))
    isNa <- function(r) {
        expect_true(is.na(r$waldRatio))
        expect_equal(r$nIV, 0L)
    }
    local_mocked_bindings(
        .cipTopLociTable = function(x, ...) data.frame(),
        .package = "pecotmr"
    )
    isNa(pecotmr:::.cipComputeMr(NULL, g, 0.5)) # empty
    local_mocked_bindings(
        .cipTopLociTable = function(x) {
            # none pass
            data.frame(
                variant_id = sprintf("chr1:%d:A:G", 100L * (1:4)),
                pip = rep(0.1, 4),
                beta = rep(0.2, 4),
                se = rep(0.05, 4)
            )
        },
        .package = "pecotmr"
    )
    isNa(pecotmr:::.cipComputeMr(NULL, g, 0.5))
    local_mocked_bindings(
        .cipTopLociTable = function(x) {
            # no overlap
            data.frame(
                variant_id = paste0("z", 1:4),
                pip = rep(0.9, 4),
                beta = rep(0.2, 4),
                se = rep(0.05, 4)
            )
        },
        .package = "pecotmr"
    )
    isNa(pecotmr:::.cipComputeMr(NULL, g, 0.5))
    local_mocked_bindings(
        .cipTopLociTable = function(x) {
            # zero beta
            data.frame(
                variant_id = sprintf("chr1:%d:A:G", 100L * (1:4)),
                pip = rep(0.9, 4),
                beta = rep(0, 4),
                se = rep(0.05, 4)
            )
        },
        .package = "pecotmr"
    )
    isNa(pecotmr:::.cipComputeMr(NULL, g, 0.5))
})

test_that(".cipComputeMrCsAware: CS-aware composite Wald + heterogeneity", {
    tl <- data.frame(
        variant_id = sprintf("chr1:%d:A:G", 100L * (1:4)),
        cs = c(1L, 1L, 2L, 2L),
        pip = c(0.5, 0.4, 0.6, 0.3),
        beta = c(0.3, 0.2, -0.25, -0.15),
        se = rep(0.05, 4),
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(
        .cipTopLociTable = function(x) tl,
        .package = "pecotmr"
    )
    res <- pecotmr:::.cipComputeMrCsAware(
        NULL,
        .cip_gwasDf(
            sprintf("chr1:%d:A:G", 100L * (1:4)),
            c(2, -1.5, 1.8, -1.2)
        ),
        cpipCutoff = 0.5
    )
    expect_true(is.finite(res$waldRatio))
    expect_true(is.finite(res$I2))
    expect_equal(res$nCs, 2L) # both CSs clear cpip
})

test_that(".cipComputeMrCsAware finds the CS column by the ^cs fallback", {
    # The shape a stored result has: no bare `cs`, and cs_95 holding
    # "<method>_<idx>" labels rather than integers.
    tl <- data.frame(
        variant_id = sprintf("chr1:%d:A:G", 100L * (1:2)),
        cs_95 = c("susie_1", "susie_1"),
        pip = c(0.5, 0.4),
        beta = c(0.3, 0.2),
        se = rep(0.05, 2),
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(
        .cipTopLociTable = function(x) tl,
        .package = "pecotmr"
    )
    res <- pecotmr:::.cipComputeMrCsAware(
        NULL,
        .cip_gwasDf(sprintf("chr1:%d:A:G", 100L * (1:2)), c(2, -1.5)),
        cpipCutoff = 0.5
    )
    expect_true(is.finite(res$waldRatio))
    expect_equal(res$nCs, 1L)
})

test_that(".cipComputeMrCsAware runs on a stored fine-mapping result", {
    # The labels used to be read with as.integer(), which turned every one
    # into NA: csAware returned an all-NA result on any real result, and the
    # tests above missed it because they all build the table by hand.
    data(qtlFineMappingExample)
    tl <- pecotmr:::.cipTopLociTable(qtlFineMappingExample)
    gwas <- data.frame(
        variant_id = tl$variant_id,
        z = seq(2, 4, length.out = nrow(tl)),
        N = 50000,
        maf = pmin(pmax(tl$af, 0.05), 0.5)
    )
    res <- pecotmr:::.cipComputeMrCsAware(
        qtlFineMappingExample,
        gwas,
        cpipCutoff = 0.5
    )
    # 15 variants: 14 in the one credible set, 1 labelled susie_0.
    expect_equal(res$nCs, 1L)
    expect_equal(res$nIV, 14L)
    expect_true(is.finite(res$waldRatio))
})

test_that("the MR paths abort naming the exposure columns they lack", {
    g <- .cip_gwasDf(sprintf("chr1:%d:A:G", 100L * (1:4)))
    local_mocked_bindings(
        .cipTopLociTable = function(x, ...) {
            data.frame(
                variant_id = sprintf("chr1:%d:A:G", 100L * (1:4)),
                pip = rep(0.9, 4)
            )
        },
        .package = "pecotmr"
    )
    # Absent beta/se is a property of the whole table, so it would fail
    # identically for every tuple -- it is named, not silently skipped.
    expect_error(
        pecotmr:::.cipComputeMr(NULL, g, 0.5),
        "Wald-ratio path needs pip, beta, se.*has no beta, se"
    )
    expect_error(
        pecotmr:::.cipComputeMrCsAware(NULL, g, 0.5),
        "credible-set-aware path needs.*has no beta, se, cs\\."
    )
    # A table carrying everything but `cs` is usable for the Wald ratio and
    # only the CS-aware path objects.
    local_mocked_bindings(
        .cipTopLociTable = function(x, ...) {
            data.frame(
                variant_id = sprintf("chr1:%d:A:G", 100L * (1:4)),
                pip = rep(0.9, 4),
                beta = rep(0.2, 4),
                se = rep(0.05, 4)
            )
        },
        .package = "pecotmr"
    )
    expect_true(is.finite(pecotmr:::.cipComputeMr(NULL, g, 0.5)$waldRatio))
    expect_error(
        pecotmr:::.cipComputeMrCsAware(NULL, g, 0.5),
        "this result has no cs\\."
    )
})

test_that(".cipComputeMrCsAware: NA on empty, no CS, no overlap, low cpip", {
    g <- .cip_gwasDf(sprintf("chr1:%d:A:G", 100L * (1:4)))
    isNa <- function(r) expect_true(is.na(r$waldRatio))
    local_mocked_bindings(
        .cipTopLociTable = function(x, ...) data.frame(),
        .package = "pecotmr"
    )
    isNa(pecotmr:::.cipComputeMrCsAware(NULL, g, 0.5)) # empty
    local_mocked_bindings(
        .cipTopLociTable = function(x) {
            # cs all 0
            data.frame(
                variant_id = sprintf("chr1:%d:A:G", 100L * (1:4)),
                cs = rep(0L, 4),
                pip = rep(0.5, 4),
                beta = rep(0.2, 4),
                se = rep(0.05, 4)
            )
        },
        .package = "pecotmr"
    )
    isNa(pecotmr:::.cipComputeMrCsAware(NULL, g, 0.5))
    local_mocked_bindings(
        .cipTopLociTable = function(x) {
            # no overlap
            data.frame(
                variant_id = paste0("z", 1:4),
                cs = rep(1L, 4),
                pip = rep(0.5, 4),
                beta = rep(0.2, 4),
                se = rep(0.05, 4)
            )
        },
        .package = "pecotmr"
    )
    isNa(pecotmr:::.cipComputeMrCsAware(NULL, g, 0.5))
    local_mocked_bindings(
        .cipTopLociTable = function(x) {
            # cpip < cutoff
            data.frame(
                variant_id = sprintf("chr1:%d:A:G", 100L * (1:4)),
                cs = 1:4,
                pip = rep(0.1, 4),
                beta = rep(0.2, 4),
                se = rep(0.05, 4)
            )
        },
        .package = "pecotmr"
    )
    isNa(pecotmr:::.cipComputeMrCsAware(NULL, g, cpipCutoff = 0.5))
})

test_that(".cipExtractWeights scores one column of weights, not several", {
    ids <- sprintf("chr1:%d:A:G", 100L * (1:5))
    w <- c(0.1, 0.05, -0.2, 0.3, 0.0)
    weightsOf <- function(weights) {
        TwasWeights(
            studyName = "Q1",
            context = "c1",
            trait = "t1",
            method = "susie",
            entry = list(TwasWeightsRow(variantIds = ids, weights = weights)),
            ldReference = .cip_makeHandle()
        )
    }
    one <- pecotmr:::.cipExtractWeights(
        weightsOf(matrix(w, ncol = 1L, dimnames = list(NULL, "susie"))),
        .cip_tuple("susie")
    )
    expect_equal(one$weights, w)
    expect_warning(
        several <- pecotmr:::.cipExtractWeights(
            weightsOf(cbind(w, -w)),
            .cip_tuple("susie")
        ),
        "a matrix with 2 columns"
    )
    expect_null(several)
})

test_that(".cipExtractWeights: NULL on missing tuples and bad topLoci", {
    tw <- .cip_makeTwasWeights(method = "susie")
    fmr <- .cip_makeQtlFmr()
    # TwasWeights / FMR: tuple absent.
    expect_null(pecotmr:::.cipExtractWeights(
        tw,
        .cip_tuple("absent")
    ))
    expect_null(pecotmr:::.cipExtractWeights(
        fmr,
        .cip_tuple("absent")
    ))
    # FMR path: tuple present but topLoci empty / no beta column / all-NA.
    local_mocked_bindings(
        .fmrRowTopLoci = function(parts, ...) data.frame(),
        .package = "pecotmr"
    )
    expect_warning(
        out <- pecotmr:::.cipExtractWeights(
            fmr,
            .cip_tuple("susie")
        ),
        "no weights resolved"
    )
    expect_null(out)
    local_mocked_bindings(
        .fmrRowTopLoci = function(parts, ...) {
            data.frame(variant_id = "chr1:100:A:G")
        },
        .package = "pecotmr"
    )
    expect_warning(
        out <- pecotmr:::.cipExtractWeights(
            fmr,
            .cip_tuple("susie")
        ),
        "no weights resolved"
    )
    expect_null(out)
    local_mocked_bindings(
        .fmrRowTopLoci = function(parts, ...) {
            data.frame(variant_id = NA_character_, beta = NA_real_)
        },
        .package = "pecotmr"
    )
    expect_warning(
        out <- pecotmr:::.cipExtractWeights(
            fmr,
            .cip_tuple("susie")
        ),
        "no weights resolved"
    )
    expect_null(out)
})

test_that(".cipCombineAcrossMethods: one-method groups give no combined rows", {
    gr <- pecotmr:::.cipDfToGranges(data.frame(
        qtlStudy = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        gwasStudy = "G1",
        analysis = "twas",
        twasZ = 2.0,
        twasPval = 0.05,
        waldRatio = NA_real_,
        waldRatioSe = NA_real_,
        mrPval = NA_real_,
        nIV = NA_integer_,
        Q = NA_real_,
        I2 = NA_real_,
        nCs = NA_integer_,
        chrom = "chr1",
        startPos = 100L,
        endPos = 200L,
        stringsAsFactors = FALSE
    ))
    expect_identical(pecotmr:::.cipCombineAcrossMethods(gr, c("acat")), gr)
})


# ===========================================================================
# TWAS variant-drop reporting (spec 6.8): the intersection is unbiased but
# only loses power, so the caller has to be told how much it lost.
# ===========================================================================

.cip_gwasDf <- function(ids) {
    data.frame(
        variant_id = ids,
        chrom = rep("chr1", length(ids)),
        pos = seq(100L, by = 100L, length.out = length(ids)),
        z = seq(2, by = -0.5, length.out = length(ids)),
        stringsAsFactors = FALSE
    )
}

test_that(".cipWarnVariantDrop: reports the drop count and the survivors", {
    expect_warning(
        pecotmr:::.cipWarnVariantDrop(5L, 3L, "test pair"),
        "dropped 2 of 5"
    )
    expect_warning(
        pecotmr:::.cipWarnVariantDrop(5L, 3L, "test pair"),
        "remaining 3"
    )
    expect_warning(
        pecotmr:::.cipWarnVariantDrop(5L, 3L, "QTL x GWAS pair"),
        "QTL x GWAS pair"
    )
})

test_that(".cipComputeTwasZ: below 2 matches warns instead of silent NULL", {
    ids <- paste0("chr1:", 100 * (1:5), ":A:G")
    gdf <- .cip_gwasDf(ids[1]) # a single shared variant
    expect_warning(
        out <- pecotmr:::.cipComputeTwasZ(
            weights = rep(0.1, 5),
            variantIds = ids,
            gwasDf = gdf,
            gwasLd = NULL,
            label = "test pair"
        ),
        "matched only 1 of 5"
    )
    expect_null(out)
})

test_that(".cipWarnVariantDrop: a full overlap is silent", {
    expect_silent(
        pecotmr:::.cipWarnVariantDrop(3L, 3L, "test pair")
    )
})

test_that(".cipPairLabel: names both sides of the pair", {
    tuple <- list(
        qStudy = "Q1",
        qContext = "c1",
        qTrait = "t1",
        qMethod = "susie"
    )
    label <- pecotmr:::.cipPairLabel(tuple, "G1")
    expect_true(grepl("study='Q1'", label, fixed = TRUE))
    expect_true(grepl("trait='t1'", label, fixed = TRUE))
    expect_true(grepl("study='G1'", label, fixed = TRUE))
})


# ===========================================================================
# Column indexing, CV metric lookup and V alignment
# ===========================================================================

test_that("an absent GWAS column fills NA rather than dropping rows", {
    # The result table has one row per requested index either way, so a
    # missing optional column must not shorten it.
    expect_equal(pecotmr:::.cipGwasCol(NULL, c(1L, 2L)), c(NA_real_, NA_real_))
    expect_equal(pecotmr:::.cipGwasCol(c(9, 8, 7), c(1L, 3L)), c(9, 7))
})

test_that("a tuple with no CV record scores NA", {
    # NA, not 0: no cross-validation was run, which is different from a fit
    # that cross-validated badly.
    noCv <- TwasWeightsRow(variantIds = "chr1:100:A:G", weights = 0.1)
    expect_true(is.na(pecotmr:::.weightCvValue(noCv, "rsq")))
})

test_that("V is aligned to the weight names, and refuses to guess", {
    V <- matrix(1:6, 3L, 2L, dimnames = list(c("a", "b", "c"), NULL))
    # Named alignment reorders rather than assuming position.
    aligned <- pecotmr:::.twasZAlignV(V, c("c", "a"), 2L)
    expect_equal(dim(aligned), c(2L, 2L))
    expect_equal(unname(aligned[, 1L]), c(3L, 1L))
    expect_equal(rownames(aligned), c("c", "a"))
    # A weight variant with no row in V is an error, not a silent drop.
    expect_error(
        pecotmr:::.twasZAlignV(V, "zz", 1L),
        "V is missing rows for 1 variant"
    )
})


# ===========================================================================
# Method selection and the CS-aware MR composite
# ===========================================================================

test_that("selecting the best method over an empty frame is a no-op", {
    cells0 <- tibble(
        study = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        rsq = numeric(0),
        pval = numeric(0),
        zero = logical(0)
    )
    expect_equal(
        nrow(pecotmr:::.selectWeightCells(
            cells0,
            WeightSelectionParam(),
            "f()"
        )),
        0L
    )
})

test_that("a credible set below the cumulative-PIP cutoff adds nothing", {
    # cpip is the SUM of member PIPs; a set that does not reach the cutoff is
    # dropped rather than contributing a poorly-supported Wald ratio.
    inst <- tibble(
        cs = c("L1", "L1"),
        pip = c(0.1, 0.1),
        bhatX = c(2, 2),
        sbhatX = c(0.1, 0.1),
        bhatY = c(1, 2),
        sbhatY = c(0.1, 0.1)
    )
    expect_null(pecotmr:::.cipCsComposite("L1", inst, 0.5))
    # A set id that matches nothing has cpip 0 and is likewise dropped.
    expect_null(pecotmr:::.cipCsComposite("L9", inst, 0.5))
})

test_that("a zero exposure effect yields no composite rather than infinity", {
    # by/bx is the Wald ratio; bx == 0 makes it non-finite, so the set is
    # dropped instead of propagating Inf into the IVW pooling.
    inst <- tibble(
        cs = c("L1", "L1"),
        pip = c(0.5, 0.5),
        bhatX = c(0, 0),
        sbhatX = c(0.1, 0.1),
        bhatY = c(1, 2),
        sbhatY = c(0.1, 0.1)
    )
    expect_null(pecotmr:::.cipCsComposite("L1", inst, 0.5))
})

test_that("a well-supported credible set yields a finite composite", {
    inst <- tibble(
        cs = c("L1", "L1"),
        pip = c(0.5, 0.5),
        bhatX = c(2, 2),
        sbhatX = c(0.1, 0.1),
        bhatY = c(1, 2),
        sbhatY = c(0.1, 0.1)
    )
    out <- pecotmr:::.cipCsComposite("L1", inst, 0.5)
    expect_true(is.finite(out$bhat))
    expect_gt(out$sbhat, 0)
})

# ---------------------------------------------------------------------------
# Work-list and scoring guards: the pipeline refuses to return an empty result
# silently, and a tuple that cannot be scored drops out rather than erroring.
# ---------------------------------------------------------------------------

test_that(".cipChosenCells rejects empty weight collections", {
    local_mocked_bindings(
        .weightCells = function(...) tibble(method = character(0)),
        .package = "pecotmr"
    )
    expect_error(
        pecotmr:::.cipChosenCells(NULL, WeightSelectionParam()),
        "no QTL tuples to score"
    )
})

test_that(".cipChosenCells errors when no method clears the CV cutoffs", {
    tw <- .cip_methodsWeights("susie", 0.01)
    # Naming the cutoffs in the message is what makes this actionable.
    expect_error(
        pecotmr:::.cipChosenCells(
            tw,
            WeightSelectionParam(cutoff = 0.5, pvalCutoff = 0.05)
        ),
        "cutoff = 0.5, pvalCutoff = 0.05"
    )
})

test_that(".cipRun errors when no tuple produced a result", {
    local_mocked_bindings(
        .cipValidateInputs = function(...) invisible(NULL),
        .cipCheckLdRefs = function(...) NULL,
        .cipChosenCells = function(...) {
            tibble(study = "s", context = "c", trait = "t", method = "m")
        },
        .cipTwasTuple = function(...) NULL,
        .cipCheckBlockReach = function(...) NULL,
        .rowSpans = function(...) NULL,
        .package = "pecotmr"
    )
    expect_error(
        pecotmr:::.cipRun(
            gwasSumStats = NULL,
            twasWeights = NULL,
            fineMappingResult = NULL,
            combineMethods = NULL,
            weightSelectionParam = WeightSelectionParam(),
            alleleFlip = FALSE,
            mrParam = MrParam(method = "csAware", cpipCutoff = 0)
        ),
        "no \\(qtl, gwas\\) tuples produced a result"
    )
})

test_that(".cipTwasTuple skips a tuple with no weights", {
    local_mocked_bindings(
        .cipExtractWeights = function(...) NULL,
        .package = "pecotmr"
    )
    cells <- tibble(study = "s", context = "c", trait = "t", method = "m")
    expect_null(pecotmr:::.cipTwasTuple(1L, cells, NULL))
})

test_that("causalInferencePipeline refuses a tuple past given blocks", {
    # The block window ends before chr1:500, which a weight sits on.
    gws <- suppressWarnings(GwasSumStats(
        studyName = "G1",
        entry = list(sumStats(.cip_makeGwasSumstats())),
        genome = "hg19",
        ldReference = .cip_makeHandle(),
        qcOptions = list(step1 = "ok"),
        region = GenomicRanges::GRanges("chr1:1-450")
    ))
    expect_error(
        causalInferencePipeline(
            gwasSumStats = gws,
            twasWeights = .cip_makeTwasWeights()
        ),
        "1 variant(s) past the end of block chr1:1-450 of GWAS study 'G1'",
        fixed = TRUE
    )
})

test_that(".cipCheckBlockReach warns of a tuple running past computed spans", {
    gws <- GwasSumStats(
        studyName = "G1",
        entry = list(sumStats(.cip_makeGwasSumstats())[1:3]),
        genome = "hg19"
    )
    tuple <- list(
        tuple = list(
            qStudy = "Q1",
            qContext = "c1",
            qTrait = "t1",
            qMethod = "m"
        ),
        variants = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(100L * (1:5), width = 1L)
        )
    )
    expect_warning(
        pecotmr:::.cipCheckBlockReach(list(tuple), gws),
        str_c(
            "QTL (study='Q1', context='c1', trait='t1', method='m'): ",
            "2 variant(s) past the end of block chr1 (chr1:100-300)"
        ),
        fixed = TRUE
    )
})

test_that(".cipScoreTwasPair skips a pair with no TWAS z", {
    local_mocked_bindings(
        .cipRowsLdRef = function(...) NULL,
        .ssStitchElements = function(...) GenomicRanges::GRanges(),
        .entryToSumstatDf = function(...) data.frame(SNP = "v1", Z = 1),
        .cipComputeTwasZ = function(...) NULL,
        .cipPairLabel = function(...) "lab",
        .package = "pecotmr"
    )
    expect_null(
        pecotmr:::.cipScoreTwasPair(
            1L,
            "G1",
            list(tuple = list(), weightsInfo = list()),
            gwasSumStats = NULL,
            window = NULL,
            alleleFlip = FALSE
        )
    )
})

test_that(".cipRunMr routes to the CS-aware estimator when asked", {
    local_mocked_bindings(
        .cipComputeMrCsAware = function(...) list(SENTINEL = TRUE),
        .package = "pecotmr"
    )
    out <- pecotmr:::.cipRunMr(
        fmrEntry = "notNull",
        gdf = NULL,
        genePvals = 0.001,
        alleleFlip = FALSE,
        mrParam = MrParam(method = "csAware", cpipCutoff = 0.5, pvalCutoff = 1)
    )
    expect_true(isTRUE(out$SENTINEL))
})

test_that(".cipRunMr gates on the gene's TWAS p-values, if it has any", {
    local_mocked_bindings(
        .cipComputeMr = function(...) list(SENTINEL = TRUE),
        .package = "pecotmr"
    )
    ran <- function(genePvals, mrParam = MrParam(pvalCutoff = 0.05)) {
        isTRUE(
            pecotmr:::.cipRunMr("fmr", NULL, genePvals, FALSE, mrParam)$SENTINEL
        )
    }
    expect_true(ran(c(0.2, 0.01)))
    expect_false(ran(c(0.2, NA)))
    # No TWAS row for the gene: nothing to gate on.
    expect_true(ran(NULL))
    # A cutoff of 1 gates nothing.
    expect_true(ran(0.9, MrParam()))
})

# One-row weights for (Q1, c1, t1, susie) whose cross-validation is `cv`.
.cip_cvWeights <- function(cv) {
    TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(TwasWeightsRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0),
            cvResult = cv
        ))
    )
}

# The (Q1, c1, t1, susie) cell's row of `.cip_cvWeights(cv)`.
.cip_cvPart <- function(cv) {
    pecotmr:::.twrCellParts(.cip_cvWeights(cv), 1L, "susie")
}

test_that(".weightCvValue reads metrics wrapped or bare, NA when absent", {
    # A bare metrics list (no $metrics wrapper) is read directly.
    bare <- .cip_cvPart(list(rsq = 0.42, pval = 0.01))
    expect_equal(pecotmr:::.weightCvValue(bare, "rsq"), 0.42)
    expect_true(is.na(pecotmr:::.weightCvValue(bare, "zzz")))
    wrapped <- .cip_cvPart(list(metrics = c(rsq = 0.3, pval = 0.02)))
    expect_equal(pecotmr:::.weightCvValue(wrapped, "pval"), 0.02)
    # The first candidate the record holds is the one read.
    expect_equal(
        pecotmr:::.weightCvValue(wrapped, c("adj_rsq_pval", "pval")),
        0.02
    )
    # A one-row metrics matrix, as a fine-mapping CV keeps it.
    asMatrix <- .cip_cvPart(list(
        metrics = matrix(
            c(0.3, 0.02),
            1L,
            dimnames = list(NULL, c("rsq", "pval"))
        )
    ))
    expect_equal(pecotmr:::.weightCvValue(asMatrix, "rsq"), 0.3)
})

test_that(".weightCvValue returns NA when the CV result cannot be read", {
    expect_true(is.na(pecotmr:::.weightCvValue(.cip_cvPart(NULL), "rsq")))
    local_mocked_bindings(
        cvResult = function(...) stop("unreadable"),
        .package = "pecotmr"
    )
    expect_true(is.na(pecotmr:::.weightCvValue(.cip_cvPart(NULL), "rsq")))
})

test_that("CV selection warns once about narrowed weights; scoring does not", {
    ids <- sprintf("chr1:%d:A:G", 100L * (1:5))
    tw <- TwasWeights(
        studyName = rep("Q1", 2),
        context = rep("c1", 2),
        trait = rep("t1", 2),
        method = c("susie", "lasso"),
        entry = list(
            TwasWeightsRow(
                variantIds = ids,
                weights = c(0.1, 0.05, -0.2, 0.3, 0.0),
                methodFits = list(variantIds = ids),
                cvResult = list(metrics = c(rsq = 0.2, pval = 0.001))
            ),
            TwasWeightsRow(
                variantIds = ids,
                weights = c(0.2, 0.1, -0.1, 0.2, 0.1),
                methodFits = list(variantIds = ids),
                cvResult = list(metrics = c(rsq = 0.5, pval = 0.001))
            )
        ),
        ldReference = .cip_makeHandle()
    )
    narrowed <- endoapply(tw, utils::head, 4L)
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    selected <- .fitCoverageMessages(causalInferencePipeline(
        gwasSumStats = .cip_makeGwasSumstats(),
        twasWeights = narrowed,
        weightSelectionParam = WeightSelectionParam(cutoff = 0.1)
    ))
    expect_length(selected, 1L)
    expect_match(selected, "causalInferencePipeline\\(\\): 2 row\\(s\\)")
    # A choice that reads no CV metric -- every method, no cutoff -- reads
    # only the weights, which the row holds.
    expect_length(
        .fitCoverageMessages(causalInferencePipeline(
            gwasSumStats = .cip_makeGwasSumstats(),
            twasWeights = narrowed,
            weightSelectionParam = WeightSelectionParam(methods = "all")
        )),
        0L
    )
})

test_that("twasZ: nSketch / combineMethods guards fire", {
    expect_error(
        twasZ(matrix(1, 1, 1), 1, nSketch = 0),
        "nSketch.*Must be >= 1"
    )
    expect_error(
        twasZ(matrix(1, 1, 1), 1, combineMethods = 1L),
        "combineMethods.*Must be of type 'character'"
    )
})

# ---------------------------------------------------------------------------
# Empty-input guards. These exist because list_c() errors on an empty list,
# so the guard is the contract: nothing in, empty of the right type out.
# ---------------------------------------------------------------------------

test_that("the concat helper answers an empty list for no pieces", {
    expect_identical(pecotmr:::.cipConcat(list()), list())
})

test_that("MrParam refuses the PIP cutoff its method does not read", {
    # Each method reads exactly one of the two. Setting the other used to be
    # silently ignored.
    expect_error(
        MrParam(method = "csAware", pipCutoff = 0.9),
        "not read by method"
    )
    expect_error(
        MrParam(method = "ivwPerVariant", cpipCutoff = 0.9),
        "not read by method"
    )
    # Leaving the inapplicable one at its default is not a request.
    expect_s4_class(
        MrParam(method = "csAware", cpipCutoff = 0.9),
        "MrParam"
    )
    expect_s4_class(MrParam(pipCutoff = 0.9), "MrParam")
})

test_that("WeightSelectionParam and MrParam do not overlap", {
    expect_length(
        intersect(
            names(formals(WeightSelectionParam)),
            names(formals(MrParam))
        ),
        1L # both carry a pvalCutoff, gating different things
    )
    expect_equal(WeightSelectionParam()$cutoff, 0)
    expect_equal(MrParam()$method, "ivwPerVariant")
})

test_that("WeightSelectionParam takes best, all, or method names", {
    expect_equal(WeightSelectionParam()$methods, "best")
    expect_equal(
        WeightSelectionParam(methods = c("lasso", "mrash"))$methods,
        c("lasso", "mrash")
    )
    expect_error(
        WeightSelectionParam(methods = c("best", "lasso")),
        "\"best\" and \"all\" stand alone"
    )
    expect_error(
        WeightSelectionParam(methods = character(0)),
        "Must have length >= 1"
    )
    expect_error(
        WeightSelectionParam(methods = c("lasso", "lasso")),
        "duplicated"
    )
})

test_that("combineMethods combines TWAS rows only, each where its gene is", {
    local_mocked_bindings(
        .handleDosage = .cip_mockExtractor(),
        .package = "pecotmr"
    )
    ids <- sprintf("chr1:%d:A:G", 100L * (1:5))
    tw <- TwasWeights(
        studyName = rep("Q1", 4L),
        context = rep("c1", 4L),
        trait = c("t1", "t1", "t2", "t2"),
        method = c("susie", "lasso", "susie", "lasso"),
        entry = list(
            TwasWeightsRow(ids, c(0.2, 0.1, -0.1, 0.2, 0.1)),
            TwasWeightsRow(ids, c(0.1, 0.2, 0.1, -0.2, 0.1)),
            TwasWeightsRow(ids[2:5], c(0.1, -0.1, 0.2, 0.1)),
            TwasWeightsRow(ids[2:5], c(0.2, 0.1, 0.1, -0.1))
        ),
        ldReference = .cip_makeHandle()
    )
    out <- causalInferencePipeline(
        .cip_makeGwasSumstats(),
        twasWeights = tw,
        fineMappingResult = .cip_makeQtlFmr(),
        weightSelectionParam = WeightSelectionParam(methods = "all"),
        combineMethods = "acat"
    )
    combined <- out[startsWith(as.character(out$method), "combined.")]
    expect_equal(as.character(combined$trait), c("t1", "t2"))
    expect_true(all(combined$analysis == "twas"))
    expect_equal(GenomicRanges::start(combined), c(100L, 200L))
    expect_equal(sum(out$analysis == "mr"), 1L)
})

test_that(".cipRowsLdRef refuses overlapped rows on two references", {
    # The rows a tuple overlaps are scored together against one LD matrix.
    block <- function(start) {
        gr <- GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(seq(start, by = 100L, length.out = 3L), width = 1L)
        )
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
            SNP = sprintf(
                "chr1:%d:A:G",
                seq(start, by = 100L, length.out = 3L)
            ),
            A1 = "A",
            A2 = "G",
            Z = c(1, 2, 3),
            N = 1000L
        )
        gr
    }
    two <- GwasSumStats(
        studyName = c("G1", "G1"),
        entry = list(block(100L), block(5000L)),
        genome = "hg19",
        ldReference = list(
            .cip_makeHandle(),
            .cip_makeHandle(sample_prefix = "z")
        ),
        blockId = c("b1", "b2"),
        qcOptions = list(step1 = "ok")
    )
    expect_error(
        pecotmr:::.cipRowsLdRef(two, 1:2, "G1"),
        "GWAS study 'G1' a QTL tuple overlaps use 2 different LD"
    )
    # A tuple overlapping one of them reads that row's reference alone.
    expect_identical(
        pecotmr:::.cipRowsLdRef(two, 2L, "G1"),
        ldReference(two)[[2L]]
    )
})

test_that("a gene is scored once per GWAS study, over the rows it overlaps", {
    # One study split into two LD blocks on chr1 plus a chr2 row on another
    # panel: the chr1 gene overlaps only block A, so it is scored once,
    # and the chr2 row's reference never enters into it.
    idsA <- sprintf("chr1:%d:A:G", 100L * (1:5))
    idsB <- sprintf("chr1:%d:A:G", 100L * (6:10))
    ids2 <- sprintf("chr2:%d:A:G", 100L * (1:5))
    ldm <- function(ids) {
        n <- length(ids)
        LdMatrix(`dimnames<-`(diag(n) * 0.5 + 0.5, list(ids, ids)))
    }
    entry <- function(ids, chrom, pos) {
        gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(pos, width = 1L))
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
            SNP = ids,
            A1 = "A",
            A2 = "G",
            Z = seq(4, 0.5, length.out = length(ids)),
            N = 1000L
        )
        gr
    }
    chr1 <- ldm(c(idsA, idsB))
    gws <- GwasSumStats(
        studyName = c("G1", "G1", "G1"),
        entry = list(
            entry(idsA, "chr1", 100L * (1:5)),
            entry(idsB, "chr1", 100L * (6:10)),
            entry(ids2, "chr2", 100L * (1:5))
        ),
        genome = "hg19",
        ldReference = list(chr1, chr1, ldm(ids2)),
        blockId = c("blockA", "blockB", "chr2"),
        qcOptions = list(step1 = "ok")
    )
    tw <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "gene1",
        method = "lasso",
        entry = list(TwasWeightsRow(
            variantIds = idsA,
            weights = c(0.3, 0.2, 0.1, 0.05, 0)
        )),
        ldReference = chr1
    )
    out <- suppressMessages(suppressWarnings(
        causalInferencePipeline(gwasSumStats = gws, twasWeights = tw)
    ))
    expect_equal(length(out), 1L)
    expect_equal(as.character(S4Vectors::mcols(out)$gwasStudy), "G1")
})

test_that("a TWAS z against an LdMatrix reads its packed LD", {
    data(eqtlRegionExample, envir = environment())
    X <- eqtlRegionExample$X[, 1:8]
    ldm <- LdMatrix(cor(X), nSamples = 100)
    ids <- colnames(X)[c(3L, 1L, 7L)]
    w <- c(0.2, -0.1, 0.05)
    z <- c(1.2, -0.4, 2.1)
    expect_equal(
        pecotmr:::.cipTwasZ(w, z, ldm, ids),
        twasZ(weights = w, z = z, R = pecotmr:::.ldFromRef(ldm, ids, "test"))
    )
})
