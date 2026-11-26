context("genePrioritization")

# ===========================================================================
# Fixture: genes on chr1 and chr6, five variants each, neighbours sharing a
# variant, so their TWAS z-scores correlate. On chr6 the first eight sit in
# the HLA region (TSS at 30 Mb) and the rest outside it (50 Mb). AR(1) LD
# within each chromosome. Two genes of `features` have no weights.
# ===========================================================================

.gp_nGenes <- 32L
.gp_nHla <- 8L
.gp_nVariants <- 4L * .gp_nGenes + 1L

.gp_layout <- function() {
    chr1 <- 1000L + 100L * seq_len(.gp_nVariants)
    # Genes 1-8 of chr6 start within the first 32 variants: the HLA ones.
    hla <- seq_len(.gp_nVariants) <= 4L * .gp_nHla
    chr6 <- ifelse(hla, 30e6, 50e6) + 100L * seq_len(.gp_nVariants)
    data.frame(
        chrom = rep(c("chr1", "chr6"), each = .gp_nVariants),
        pos = as.integer(c(chr1, chr6)),
        A1 = "T",
        A2 = "C"
    ) |>
        transform(SNP = paste(chrom, pos, "C", "T", sep = ":"))
}

.gp_R <- function() {
    block <- 0.5^abs(outer(
        seq_len(.gp_nVariants),
        seq_len(.gp_nVariants),
        "-"
    ))
    as.matrix(Matrix::bdiag(block, block))
}

.gp_gwas <- function(layout = .gp_layout()) {
    set.seed(21)
    gr <- GenomicRanges::GRanges(
        layout$chrom,
        IRanges::IRanges(layout$pos, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = layout$SNP,
        A1 = layout$A1,
        A2 = layout$A2,
        Z = rnorm(nrow(layout), sd = 2)
    )
    GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = LdMatrix(.gp_R(), variants = layout),
        qcOptions = list(qc = "done")
    )
}

# Gene k of a chromosome takes variants 4k + 1..5 of it.
.gp_geneVariants <- function(chromosome, k) {
    offset <- if (chromosome == "chr1") 0L else .gp_nVariants
    offset + 4L * (k - 1L) + 1:5
}

.gp_genes <- function() {
    data.frame(
        gene = c(
            sprintf("a%02d", seq_len(.gp_nGenes)),
            sprintf("b%02d", seq_len(.gp_nGenes))
        ),
        chrom = rep(c("chr1", "chr6"), each = .gp_nGenes),
        k = rep(seq_len(.gp_nGenes), 2)
    )
}

.gp_twasWeights <- function(layout = .gp_layout()) {
    genes <- .gp_genes()
    set.seed(5)
    entries <- map(seq_len(nrow(genes)), function(i) {
        idx <- .gp_geneVariants(genes$chrom[[i]], genes$k[[i]])
        TwasWeightsRow(
            layout$SNP[idx],
            c(rnorm(4), 0),
            cvResult = list(metrics = c(rsq = runif(1), pval = 0.01))
        )
    })
    n <- nrow(genes)
    TwasWeights(
        studyName = rep("Q1", n),
        context = rep("c1", n),
        trait = genes$gene,
        method = rep("ridge", n),
        entry = entries
    )
}

.gp_features <- function(layout = .gp_layout()) {
    genes <- .gp_genes()
    ids <- c(genes$gene, "c01", "c02")
    starts <- c(
        map_int(seq_len(nrow(genes)), function(i) {
            layout$pos[[.gp_geneVariants(genes$chrom[[i]], genes$k[[i]])[[1L]]]]
        }),
        2000L,
        3000L
    )
    set.seed(8)
    SummarizedExperiment::SummarizedExperiment(
        list(
            features = matrix(
                rbinom(length(ids) * 6L, 1L, 0.4),
                length(ids),
                dimnames = list(ids, sprintf("f%d", 1:6))
            )
        ),
        rowRanges = GenomicRanges::GRanges(
            c(genes$chrom, "chr1", "chr6"),
            IRanges::IRanges(starts, width = 300L),
            names = ids
        )
    )
}

.gp_buildInputs <- function() {
    gwas <- .gp_gwas()
    tw <- .gp_twasWeights()
    list(
        gwas = gwas,
        tw = tw,
        res = causalInferencePipeline(
            gwas,
            twasWeights = tw,
            weightSelectionParam = WeightSelectionParam(methods = "all")
        ),
        features = .gp_features()
    )
}

# Built once for the file: the inputs are seeded, so every test reads the
# same ones, and rebuilding them per test cost about 6 s a test.
.gp_shared <- .gp_buildInputs()

.gp_inputs <- function() {
    .gp_shared
}

# The shared inputs' TWAS studies and gene correlation (.gcCompute()),
# computed once. A test of what is fitted to them stands these in for the
# pipeline's own (.gp_useSharedCorrelation()), which cost about 3.4 s a call.
.gp_sharedStudies <- pecotmr:::.gcCompute(
    .gp_shared$res,
    .gp_shared$gwas,
    .gp_shared$tw,
    NULL,
    "twas",
    NULL,
    NULL,
    MrParam(),
    TRUE
)

# The shared inputs' variance components on squared TWAS weights, and their
# studies, likewise computed once (about 19 s together); NULL without
# CompQuadForm, whose tests skip.
.gp_sharedVc <- if (requireNamespace("CompQuadForm", quietly = TRUE)) {
    geneVarianceComponentTest(
        .gp_shared$gwas,
        twasWeights = .gp_shared$tw,
        variantWeights = "squaredWeight"
    )
}
.gp_sharedVcStudies <- if (!is.null(.gp_sharedVc)) {
    pecotmr:::.gcCompute(
        .gp_sharedVc,
        .gp_shared$gwas,
        .gp_shared$tw,
        NULL,
        "varianceComponent",
        NULL,
        NULL,
        MrParam(),
        TRUE
    )
}

# For the calling test: the shared inputs' TWAS rows, or their variance
# components, with default selection, get their studies without
# recomputing them.
.gp_useSharedCorrelation <- function(env = parent.frame()) {
    local_mocked_bindings(
        .gcCompute = function(
            x,
            gwasSumStats,
            twasWeights,
            fmr,
            analysis,
            ...
        ) {
            if (analysis == "varianceComponent") {
                .gp_sharedVcStudies
            } else {
                .gp_sharedStudies
            }
        },
        .package = "pecotmr",
        .env = env
    )
}

# A stand-in fit with a closed form, X'y on the whitened data, so the
# plumbing around it can be checked exactly.
.gp_crossFit <- function(fitMethod, whitened, options) {
    as.numeric(crossprod(whitened$X, whitened$y))
}

test_that("each chromosome is scored by the fit that held it out", {
    inp <- .gp_inputs()
    .gp_useSharedCorrelation()
    local_mocked_bindings(.gpFit = .gp_crossFit, .package = "pecotmr")
    out <- genePrioritization(
        inp$res,
        inp$features,
        inp$gwas,
        twasWeights = inp$tw
    )
    expect_named(out, "G1")
    scores <- out$G1$scores
    expect_s4_class(scores, "GRanges")
    expect_equal(names(scores), rownames(inp$features))
    # The response: the one-sided z of each gene's TWAS z.
    twas <- inp$res[inp$res$analysis == "twas"]
    y <- set_names(pecotmr:::.zOneSided(twas$twasZ), twas$trait)
    n <- 2L * .gp_nGenes
    expect_equal(scores$y[seq_len(n)], unname(y[names(scores)[seq_len(n)]]))
    expect_true(all(is.na(scores$y[n + 1:2])))
    # The expected fit for chr1: the chr6 genes outside HLA, whitened by
    # the one-sided correlation with the covariates projected out of y.
    corr <- geneCorrelation(inp$res, inp$gwas, twasWeights = inp$tw)$G1
    X <- SummarizedExperiment::assay(inp$features)
    W <- vapply(
        seq_len(n),
        function(i) {
            genes <- .gp_genes()
            idx <- .gp_geneVariants(genes$chrom[[i]], genes$k[[i]])
            w <- as.numeric(
                pecotmr:::.twrRowResolveWeights(
                    pecotmr:::.twrCellParts(inp$tw, i, "ridge")
                )$weights
            )
            replace(numeric(2L * .gp_nVariants), idx, w)
        },
        numeric(2L * .gp_nVariants)
    )
    variance <- set_names(
        diag(crossprod(W, .gp_R() %*% W)),
        .gp_genes()$gene
    )
    rsq <- map_dbl(seq_len(n), function(i) {
        pecotmr:::.cvMetricValue(
            pecotmr:::.weightCvMetrics(
                pecotmr:::.weightCellPart(
                    i,
                    pecotmr:::.weightCellIds(inp$tw),
                    inp$tw
                )
            ),
            "rsq"
        )
    })
    covariates <- cbind(
        nVariants = 4,
        logNVariants = log(4),
        weightVariance = variance,
        logWeightVariance = log(variance),
        cvRsq = rsq
    )
    train <- sprintf("b%02d", (.gp_nHla + 1L):.gp_nGenes)
    whitened <- pecotmr:::.glsWhiten(
        y[train],
        X[train, ],
        corr,
        covariates = covariates
    )
    coef <- as.numeric(crossprod(whitened$X, whitened$y))
    expect_equal(unname(out$G1$coefficients[, "chr1"]), coef)
    chr1 <- names(scores)[
        as.character(GenomicRanges::seqnames(scores)) == "chr1"
    ]
    expect_equal(
        scores[chr1]$score,
        as.numeric(X[chr1, ] %*% coef)
    )
    expect_equal(unique(scores[chr1]$chromosome), "chr1")
    expect_equal(colnames(out$G1$coefficients), c("chr1", "chr6"))
})

test_that("HLA genes are scored but fitted on only when asked", {
    inp <- .gp_inputs()
    .gp_useSharedCorrelation()
    local_mocked_bindings(.gpFit = .gp_crossFit, .package = "pecotmr")
    out <- genePrioritization(
        inp$res,
        inp$features,
        inp$gwas,
        twasWeights = inp$tw,
        covariates = FALSE
    )
    kept <- genePrioritization(
        inp$res,
        inp$features,
        inp$gwas,
        twasWeights = inp$tw,
        covariates = FALSE,
        excludeHla = FALSE
    )
    # chr6 is held out of its own fit, so only chr1's coefficients change.
    expect_equal(
        out$G1$coefficients[, "chr6"],
        kept$G1$coefficients[, "chr6"]
    )
    expect_false(isTRUE(all.equal(
        out$G1$coefficients[, "chr1"],
        kept$G1$coefficients[, "chr1"]
    )))
    hla <- sprintf("b%02d", seq_len(.gp_nHla))
    expect_true(all(is.finite(out$G1$scores[hla]$score)))
})

# genePrioritization() fitted by `fitMethod`, which needs `engine`: every
# fold's effects and every gene's score finite.
# A stand-in DPR fit: the effects X'w of its covariate w. It keeps the
# methodArgs default, which DprOptions() reads its constructor from.
.gp_dprCovariate <- function(X, y, methodArgs = DprOptions()) {
    as.numeric(crossprod(X, methodArgs$w))
}

.gp_expectFitRuns <- function(fitMethod, engine) {
    skip_if_not_installed(engine)
    inp <- .gp_inputs()
    .gp_useSharedCorrelation()
    out <- suppressWarnings(suppressMessages(capture.output(
        res <- genePrioritization(
            inp$res,
            inp$features,
            inp$gwas,
            twasWeights = inp$tw,
            fitMethod = fitMethod,
            seed = 1
        )
    )))
    expect_true(all(is.finite(res$G1$coefficients)))
    expect_true(all(is.finite(res$G1$scores$score)))
}

test_that("ridge fits without an intercept", {
    .gp_expectFitRuns("ridge", "glmnet")
})

test_that("lasso fits without an intercept", {
    .gp_expectFitRuns("lasso", "glmnet")
})

test_that("enet fits without an intercept", {
    .gp_expectFitRuns("enet", "glmnet")
})

test_that("l0learn fits without an intercept", {
    .gp_expectFitRuns("l0learn", "L0Learn")
})

test_that("abess fits without an intercept", {
    .gp_expectFitRuns("abess", "abess")
})

test_that("mrash fits without an intercept", {
    .gp_expectFitRuns("mrash", "susieR")
})

test_that("sparsevb fits without an intercept", {
    .gp_expectFitRuns("sparsevb", "sparsevb")
})

test_that("susie fits without an intercept", {
    .gp_expectFitRuns("susie", "susieR")
})

test_that("dprVb fits without an intercept", {
    # RcppDPR 0.1.10's VB fit corrupts the heap on Linux for some data, at any
    # size tried here (random 56 x 6 and 120 x 60 inputs included); glibc then
    # aborts a later test.
    skip("RcppDPR corrupts the heap on Linux (upstream bug)")
    .gp_expectFitRuns("dprVb", "RcppDPR")
})

test_that("the fit's options turn its intercept off, and only that", {
    .gp_cs <- CredibleSetParam()
    ridge <- pecotmr:::.gpOptions("ridge", GlmnetOptions(nfolds = 3), .gp_cs)
    expect_identical(S4Vectors::metadata(ridge)$constructor, "GlmnetOptions")
    expect_false(ridge$intercept)
    expect_equal(ridge$nfolds, 3)
    expect_false(pecotmr:::.gpOptions("abess", NULL, .gp_cs)$fit.intercept)
    expect_error(
        pecotmr:::.gpOptions("ridge", GlmnetOptions(intercept = TRUE), .gp_cs),
        "cannot set intercept"
    )
    expect_error(
        pecotmr:::.gpOptions("l0learn", GlmnetOptions(), .gp_cs),
        "must be built with L0learnOptions\\(\\), not GlmnetOptions\\(\\)"
    )
    expect_error(
        pecotmr:::.gpOptions("dprVb", DprOptions(w = matrix(1, 3, 1)), .gp_cs),
        "cannot set w"
    )
    expect_error(
        pecotmr:::.gpOptions("ridge", list(nfolds = 3), .gp_cs),
        "must be built with GlmnetOptions"
    )
})

test_that("credibleSetParam sets SuSiE's L, and only SuSiE's", {
    susie <- pecotmr:::.gpOptions(
        "susie",
        SusieOptions(max_iter = 50),
        CredibleSetParam(maxNumSingleEffects = 25L)
    )
    expect_equal(susie$L, 25L)
    expect_equal(susie$max_iter, 50)
    expect_false(susie$intercept)
    expect_null(susie$L_greedy)
    greedy <- pecotmr:::.gpOptions(
        "susie",
        NULL,
        CredibleSetParam(maxNumSingleEffectsGreedy = 5L)
    )
    expect_equal(greedy$L, 10L)
    expect_equal(greedy$L_greedy, 5L)
    inp <- .gp_inputs()
    .gp_useSharedCorrelation()
    local_mocked_bindings(
        .fmFitSusieIndiv = function(X, y, token, userArgs = NULL, ...) {
            abort(str_c("fitted with L = ", userArgs$L))
        },
        .package = "pecotmr"
    )
    expect_error(
        genePrioritization(
            inp$res,
            inp$features,
            inp$gwas,
            twasWeights = inp$tw,
            fitMethod = "susie",
            credibleSetParam = CredibleSetParam(maxNumSingleEffects = 25L)
        ),
        "fitted with L = 25"
    )
    expect_error(
        genePrioritization(
            inp$res,
            inp$features,
            inp$gwas,
            twasWeights = inp$tw,
            credibleSetParam = CredibleSetParam(maxNumSingleEffects = 25L)
        ),
        "sets SuSiE's L, and the ridge fit has none"
    )
})

test_that("MR covariates count the instruments; missing ones are left out", {
    study <- list(
        genes = list(
            list(gene = "a", nIV = 2),
            list(gene = "b", nIV = 5)
        ),
        variance = c(1, 1)
    )
    expect_equal(
        pecotmr:::.gpCovariates(study, "mr", NULL),
        cbind(nInstruments = c(a = 2, b = 5), logNInstruments = log(c(2, 5)))
    )
    twas <- list(
        genes = list(
            list(gene = "a", coef = c(1, 0, 2), cellKey = "ka"),
            list(gene = "b", coef = c(1, 1, 1), cellKey = "kb")
        ),
        variance = c(2, 3)
    )
    expect_message(
        table <- pecotmr:::.gpCovariates(
            twas,
            "twas",
            c(ka = 0.2, kb = NA)
        ),
        "cvRsq are missing for some genes"
    )
    expect_equal(
        colnames(table),
        c(
            "nVariants",
            "logNVariants",
            "weightVariance",
            "logWeightVariance"
        )
    )
    expect_equal(unname(table[, "nVariants"]), c(2, 3))
    # With every covariate missing for some gene, there are none.
    unknown <- list_modify(
        study,
        genes = map(study$genes, list_modify, nIV = NA_real_)
    )
    expect_message(
        none <- pecotmr:::.gpCovariates(unknown, "mr", NULL),
        "nInstruments, logNInstruments are missing"
    )
    expect_null(none)
})

test_that("a DPR fit takes the whitened intercept as its covariate", {
    options <- pecotmr:::.gpOptions("dprVb", NULL, CredibleSetParam())
    # Each stand-in returns X'w, so its effects show the w it was given.
    local_mocked_bindings(
        dprVbWeights = .gp_dprCovariate,
        dprGibbsWeights = .gp_dprCovariate,
        dprAdaptiveGibbsWeights = .gp_dprCovariate,
        .package = "pecotmr"
    )
    whitened <- list(
        X = matrix(c(1, 0, 2, 1, 1, 0), 3),
        y = c(0.5, -1, 2),
        intercept = matrix(c(0.6, 0.3, -0.2))
    )
    for (fitMethod in c("dprVb", "dprGibbs", "dprAdaptiveGibbs")) {
        expect_equal(
            pecotmr:::.gpFit(fitMethod, whitened, options),
            as.numeric(crossprod(whitened$X, whitened$intercept))
        )
    }
})

test_that("genePrioritization refuses features it cannot read", {
    inp <- .gp_inputs()
    .gp_useSharedCorrelation()
    expect_error(
        genePrioritization(
            inp$res,
            SummarizedExperiment::assay(inp$features),
            inp$gwas,
            twasWeights = inp$tw
        ),
        "RangedSummarizedExperiment"
    )
    unnamed <- inp$features
    rownames(unnamed) <- NULL
    expect_error(
        genePrioritization(inp$res, unnamed, inp$gwas, twasWeights = inp$tw),
        "row names"
    )
    bare <- SummarizedExperiment::SummarizedExperiment(
        rowRanges = SummarizedExperiment::rowRanges(inp$features)
    )
    expect_error(
        genePrioritization(inp$res, bare, inp$gwas, twasWeights = inp$tw),
        "`features` has no assay of features"
    )
    # Only chr1 genes in features: holding chr1 out leaves nothing to fit.
    chr1 <- inp$features[seq_len(.gp_nGenes), ]
    expect_error(
        suppressMessages(genePrioritization(
            inp$res,
            chr1,
            inp$gwas,
            twasWeights = inp$tw
        )),
        "holding out chr1 leaves 0 gene"
    )
})

# ===========================================================================
# geneCorrelation
# ===========================================================================

# ===========================================================================
# Fixture: 30 variants on chr1 with AR(1) LD in one LdMatrix, a GWAS split
# into two rows (variants 1-15 and 16-30), and four genes: g1 and g2 in the
# first row, g3 across both, g4 in the second. g2's weights name the variants
# with their alleles swapped, so they enter in the GWAS coding negated.
# ===========================================================================

.gc_p <- 30L
.gc_R <- 0.6^abs(outer(seq_len(.gc_p), seq_len(.gc_p), "-"))
.gc_variants <- .testLdVariants(.gc_p, "chr1", 1000L)
.gc_ids <- .gc_variants$SNP
.gc_flipped <- sprintf("chr1:%d:T:C", .gc_variants$pos)

.gc_genes <- list(
    g3 = list(idx = 12:20, w = c(0.3, -0.1, 0.2, 0.4, 0, 0.1, -0.2, 0.3, 0.1)),
    g1 = list(idx = 1:6, w = c(0.4, 0.2, 0, 0.1, 0.3, -0.2)),
    g4 = list(idx = 22:28, w = c(0.2, 0.2, 0.1, -0.3, 0.4, 0.1, 0.2)),
    g2 = list(idx = 4:10, w = c(0.1, 0.3, 0.2, -0.1, 0.2, 0.3, 0.1))
)

.gc_z <- function() {
    set.seed(11)
    rnorm(.gc_p, sd = 2)
}

.gc_gwas <- function(z = .gc_z(), study = "G1") {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(.gc_variants$pos, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = .gc_ids,
        A1 = "T",
        A2 = "C",
        Z = z,
        N = 5000L,
        MAF = .gc_variants$AF
    )
    GwasSumStats(
        studyName = study,
        entry = list(gr),
        genome = "hg19",
        ldReference = LdMatrix(.gc_R, variants = .gc_variants, nSamples = 500),
        qcOptions = list(qc = "done"),
        region = GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(c(1000L, 2500L), c(2400L, 3900L))
        )
    )
}

.gc_geneIds <- function(name) {
    gene <- .gc_genes[[name]]
    if (name == "g2") .gc_flipped[gene$idx] else .gc_ids[gene$idx]
}

.gc_twasWeights <- function(method = "lasso") {
    n <- length(.gc_genes)
    TwasWeights(
        studyName = rep("Q1", n),
        context = rep("c1", n),
        trait = names(.gc_genes),
        method = rep(method, n),
        entry = map2(
            map(names(.gc_genes), .gc_geneIds),
            map(.gc_genes, "w"),
            TwasWeightsRow
        )
    )
}

.gc_run <- function(gwas, tw) {
    causalInferencePipeline(
        gwas,
        twasWeights = tw,
        weightSelectionParam = WeightSelectionParam(methods = "all")
    )
}

# The genes' weights in the GWAS coding as columns over all variants, in
# genomic order g1..g4.
.gc_weightMatrix <- function() {
    genes <- c("g1", "g2", "g3", "g4")
    W <- vapply(
        genes,
        function(g) {
            gene <- .gc_genes[[g]]
            sign <- if (g == "g2") -1 else 1
            replace(numeric(.gc_p), gene$idx, sign * gene$w)
        },
        numeric(.gc_p)
    )
    W
}

# rho = w_g' R w_h / sqrt(...), zero for genes that share no GWAS row.
.gc_expected <- function(W) {
    C <- crossprod(W, .gc_R %*% W)
    rho <- C / sqrt(outer(diag(C), diag(C)))
    # g1, g2 in row 1; g3 in both; g4 in row 2.
    rho["g1", "g4"] <- rho["g4", "g1"] <- 0
    rho["g2", "g4"] <- rho["g4", "g2"] <- 0
    rho
}

test_that("TWAS z correlation is w'Rw normalized, within shared GWAS rows", {
    gwas <- .gc_gwas()
    tw <- .gc_twasWeights()
    res <- .gc_run(gwas, tw)
    out <- geneCorrelation(res, gwas, twasWeights = tw, scale = "signed")
    expect_named(out, "G1")
    m <- out$G1
    expect_s4_class(m, "dsCMatrix")
    # Genomic order, whatever order the weights came in.
    expect_equal(rownames(m), c("g1", "g2", "g3", "g4"))
    expect_equal(as.matrix(m), .gc_expected(.gc_weightMatrix()))
    # Genes that share no row are structural zeros, not stored values.
    expect_equal(m["g1", "g4"], 0)
    expect_false(m["g1", "g3"] == 0)
})

test_that("each scale transforms the z-scores' correlation", {
    gwas <- .gc_gwas()
    tw <- .gc_twasWeights()
    res <- .gc_run(gwas, tw)
    signed <- as.matrix(
        geneCorrelation(res, gwas, twasWeights = tw, scale = "signed")$G1
    )
    squared <- geneCorrelation(res, gwas, twasWeights = tw, scale = "squared")
    expect_equal(as.matrix(squared$G1), signed^2)
    absolute <- geneCorrelation(
        res,
        gwas,
        twasWeights = tw,
        scale = "absolute"
    )
    r <- signed["g1", "g2"]
    expect_equal(
        absolute$G1["g1", "g2"],
        (sqrt(1 - r^2) + r * asin(r) - 1) / (pi / 2 - 1)
    )
    oneSided <- geneCorrelation(res, gwas, twasWeights = tw)
    expect_equal(
        oneSided$G1["g1", "g2"],
        pecotmr:::.oneSidedZCorrelation(r)
    )
    expect_equal(unname(diag(as.matrix(oneSided$G1))), rep(1, 4))
})

test_that("the |z| and one-sided correlations match simulated z-scores", {
    set.seed(3)
    n <- 2e6
    z1 <- rnorm(n)
    rhos <- c(-0.9, 0.3, 0.7, 0.99)
    simulated <- vapply(
        rhos,
        function(r) {
            z2 <- r * z1 + sqrt(1 - r^2) * rnorm(n)
            c(
                abs = cor(abs(z1), abs(z2)),
                oneSided = cor(
                    pecotmr:::.zOneSided(z1),
                    pecotmr:::.zOneSided(z2)
                )
            )
        },
        numeric(2)
    )
    # Monte Carlo error at n = 2e6 is about 7e-4.
    expect_equal(
        pecotmr:::.gcScaleCorrelation(rhos, "absolute"),
        simulated["abs", ],
        tolerance = 3e-3
    )
    expect_equal(
        pecotmr:::.oneSidedZCorrelation(rhos),
        simulated["oneSided", ],
        tolerance = 3e-3
    )
    expect_equal(pecotmr:::.oneSidedZCorrelation(c(0, 1, -1)), c(0, 1, 1))
})

test_that("the one-sided z is Phi^-1(1 - p) of the two-sided p-value", {
    z <- c(-6, -2, -0.3, 0.01, 1, 3, 9)
    expect_equal(
        pecotmr:::.zOneSided(z),
        qnorm(2 * pnorm(-abs(z)), lower.tail = FALSE)
    )
})

test_that("one matrix per GWAS study", {
    gwas <- combineGwasSumStats(
        .gc_gwas(study = "G1"),
        .gc_gwas(study = "G2")
    )
    tw <- .gc_twasWeights()
    res <- .gc_run(gwas, tw)
    out <- geneCorrelation(res, gwas, twasWeights = tw, scale = "signed")
    expect_named(out, c("G1", "G2"))
    expect_equal(as.matrix(out$G2), .gc_expected(.gc_weightMatrix()))
})

test_that("a gene with two rows needs method to choose one", {
    gwas <- .gc_gwas()
    tw <- combineTwasWeights(
        .gc_twasWeights("lasso"),
        .gc_twasWeights("enet")
    )
    res <- .gc_run(gwas, tw)
    expect_error(
        geneCorrelation(res, gwas, twasWeights = tw),
        "more than one twas row against one GWAS study"
    )
    out <- geneCorrelation(
        res,
        gwas,
        twasWeights = tw,
        method = "enet",
        scale = "signed"
    )
    expect_equal(as.matrix(out$G1), .gc_expected(.gc_weightMatrix()))
    expect_error(
        geneCorrelation(res, gwas, twasWeights = tw, method = "susie"),
        "no twas row with a statistic for method 'susie'"
    )
})

test_that("inputs other than x's are refused", {
    gwas <- .gc_gwas()
    tw <- .gc_twasWeights()
    res <- .gc_run(gwas, tw)
    set.seed(99)
    other <- .gc_gwas(z = rnorm(.gc_p))
    expect_error(
        geneCorrelation(res, other, twasWeights = tw),
        "does not match the one these inputs give"
    )
    # Without allele flips g2 matches nothing, so it has no weights at all.
    expect_error(
        geneCorrelation(res, gwas, twasWeights = tw, alleleFlip = FALSE),
        "trait=g2"
    )
    expect_error(
        geneCorrelation(res, gwas),
        "TWAS rows need the weights"
    )
    expect_error(
        geneCorrelation(res, gwas, twasWeights = tw, analysis = "mr"),
        "MR rows need the `fineMappingResult`"
    )
    expect_error(
        geneCorrelation(GenomicRanges::GRanges(), gwas, twasWeights = tw),
        "is not a causalInferencePipeline\\(\\) result"
    )
    # Weights between the two GWAS rows meet none of them.
    gap <- TwasWeights(
        "Q1",
        "c1",
        "g1",
        "lasso",
        list(TwasWeightsRow(c("chr1:2410:C:T", "chr1:2420:C:T"), c(1, 1)))
    )
    expect_error(
        geneCorrelation(res, gwas, twasWeights = gap),
        "trait=g1, method=lasso\\) x GWAS 'G1' does not match"
    )
    expect_error(
        geneCorrelation(res[0], gwas, twasWeights = tw),
        "`x` holds no twas row with a statistic\\.$"
    )
})

# ---------------------------------------------------------------------------
# MR: a fine-mapping row per gene, its instruments the top loci with pip > 0.5
# in one credible set.
# ---------------------------------------------------------------------------

.gc_mrGenes <- list(
    g1 = list(idx = c(2L, 4L, 5L), beta = c(0.3, -0.2, 0.25)),
    g2 = list(idx = c(6L, 8L), beta = c(0.2, 0.15)),
    g3 = list(idx = c(14L, 17L), beta = c(-0.3, 0.2))
)

# `effectScale` scales the instruments' effects; `set` labels their credible
# set (NA for none).
.gc_fineMapping <- function(effectScale = 1, set = "susie_1") {
    rows <- map(.gc_mrGenes, function(gene) {
        n <- length(gene$idx)
        FineMappingRow(
            variantIds = .gc_ids[gene$idx],
            susieFit = list(),
            topLoci = data.frame(
                variant_id = .gc_ids[gene$idx],
                pip = seq(0.9, 0.6, length.out = n),
                posterior_mean = effectScale * gene$beta,
                posterior_sd = rep(0.05, n),
                cs_0.95 = rep(set, n)
            )
        )
    })
    n <- length(rows)
    QtlFineMappingResult(
        studyName = rep("Q1", n),
        context = rep("c1", n),
        trait = names(rows),
        method = rep("susie", n),
        entry = unname(rows)
    )
}

# The outcome's se and the exposure for gene `g`'s instruments.
.gc_mrInputs <- function(gene, z) {
    maf <- .gc_variants$AF[gene$idx]
    zg <- z[gene$idx]
    seY <- 1 / sqrt(2 * maf * (1 - maf) * (5000 + zg^2))
    list(
        z = zg,
        seY = seY,
        betaY = zg * seY,
        betaX = gene$beta,
        seX = rep(0.05, length(gene$idx)),
        pip = seq(0.9, 0.6, length.out = length(gene$idx))
    )
}

# Per-variant IVW: a = w seY / betaX / sum(w), w = 1 / rSe^2.
.gc_ivwColumn <- function(gene, z) {
    x <- .gc_mrInputs(gene, z)
    rSe <- sqrt(
        (x$seY / x$betaX)^2 + (x$betaY * x$seX / x$betaX^2)^2
    )
    w <- 1 / rSe^2
    replace(numeric(.gc_p), gene$idx, w * x$seY / x$betaX / sum(w))
}

# One credible set: its PIP-weighted ratio, a = pNorm seY / betaX.
.gc_csColumn <- function(gene, z) {
    x <- .gc_mrInputs(gene, z)
    replace(
        numeric(.gc_p),
        gene$idx,
        x$pip / sum(x$pip) * x$seY / x$betaX
    )
}

.gc_mrExpected <- function(column, z) {
    A <- vapply(.gc_mrGenes, column, numeric(.gc_p), z = z)
    C <- crossprod(A, .gc_R %*% A)
    C / sqrt(outer(diag(C), diag(C)))
}

test_that("MR correlation uses the IVW estimate's linear weights", {
    gwas <- .gc_gwas()
    z <- .gc_z()
    fmr <- .gc_fineMapping()
    res <- causalInferencePipeline(gwas, fineMappingResult = fmr)
    out <- geneCorrelation(
        res,
        gwas,
        fineMappingResult = fmr,
        analysis = "mr",
        scale = "signed"
    )
    expect_equal(as.matrix(out$G1), .gc_mrExpected(.gc_ivwColumn, z))
    cs <- MrParam(method = "csAware")
    resCs <- causalInferencePipeline(
        gwas,
        fineMappingResult = fmr,
        mrParam = cs
    )
    outCs <- geneCorrelation(
        resCs,
        gwas,
        fineMappingResult = fmr,
        analysis = "mr",
        scale = "signed",
        mrParam = cs
    )
    expect_equal(as.matrix(outCs$G1), .gc_mrExpected(.gc_csColumn, z))
    # The other MR method's weights give other estimates.
    expect_error(
        geneCorrelation(
            res,
            gwas,
            fineMappingResult = fmr,
            analysis = "mr",
            mrParam = cs
        ),
        "does not match"
    )
})

test_that("MR inputs that give a gene no estimate are refused", {
    gwas <- .gc_gwas()
    fmr <- .gc_fineMapping()
    res <- causalInferencePipeline(gwas, fineMappingResult = fmr)
    cs <- MrParam(method = "csAware")
    resCs <- causalInferencePipeline(
        gwas,
        fineMappingResult = fmr,
        mrParam = cs
    )
    mismatch <- "trait=g1, method=susie\\) x GWAS 'G1' does not match"
    # No instrument passes the cutoff.
    expect_error(
        geneCorrelation(
            res,
            gwas,
            fineMappingResult = fmr,
            analysis = "mr",
            mrParam = MrParam(pipCutoff = 0.95)
        ),
        mismatch
    )
    # Instruments with no effect have no finite ratio, per variant or per set.
    null <- .gc_fineMapping(effectScale = 0)
    expect_error(
        geneCorrelation(res, gwas, fineMappingResult = null, analysis = "mr"),
        mismatch
    )
    expect_error(
        geneCorrelation(
            resCs,
            gwas,
            fineMappingResult = null,
            analysis = "mr",
            mrParam = cs
        ),
        mismatch
    )
    # Instruments in no credible set.
    expect_error(
        geneCorrelation(
            resCs,
            gwas,
            fineMappingResult = .gc_fineMapping(set = NA_character_),
            analysis = "mr",
            mrParam = cs
        ),
        mismatch
    )
})

test_that("a variant in two credible sets sums its coefficients", {
    out <- pecotmr:::.gcByVariant(c("a", "b", "a"), c(1, 2, 3), c(5, 6, 5))
    expect_equal(out$variantIds, c("a", "b"))
    expect_equal(out$coef, c(4, 2))
    expect_equal(out$z, c(5, 6))
})

# ---------------------------------------------------------------------------
# Whitening
# ---------------------------------------------------------------------------

.gc_whitenFixture <- function() {
    genes <- c("a", "b", "c", "d", "e")
    sigma <- Matrix::sparseMatrix(
        i = c(1:5, 1, 1, 4),
        j = c(1:5, 2, 3, 5),
        x = c(rep(1, 5), 0.5, -0.3, 0.4),
        dims = c(5, 5),
        dimnames = list(genes, genes),
        symmetric = TRUE
    )
    set.seed(5)
    list(
        sigma = sigma,
        y = set_names(rnorm(5, sd = 2), genes),
        X = matrix(
            rbinom(15, 1, 0.5),
            5,
            dimnames = list(rev(genes), c("f1", "f2", "f3"))
        ),
        covariates = matrix(rnorm(5), 5, dimnames = list(genes, "size"))
    )
}

test_that(".glsWhiten whitens by the ridged correlation as PoPS does", {
    fx <- .gc_whitenFixture()
    out <- pecotmr:::.glsWhiten(
        fx$y,
        fx$X,
        fx$sigma,
        covariates = fx$covariates
    )
    dense <- as.matrix(fx$sigma)
    ridge <- abs(min(0, eigen(dense)$values)) +
        0.05 +
        0.9 * max(0, var(fx$y) - 1)
    Linv <- solve(t(chol(dense + diag(ridge, 5))))
    one <- Linv %*% rep(1, 5)
    yW <- Linv %*% fx$y
    XW <- Linv %*% fx$X[names(fx$y), ]
    nuisance <- cbind(one, Linv %*% fx$covariates)
    expect_equal(
        out$y,
        set_names(as.numeric(lm.fit(nuisance, yW)$residuals), names(fx$y))
    )
    expect_equal(
        out$X,
        `dimnames<-`(
            lm.fit(one, XW)$residuals,
            list(names(fx$y), colnames(fx$X))
        )
    )
})

test_that(".glsWhiten splits the correlation into its independent blocks", {
    fx <- .gc_whitenFixture()
    blocks <- pecotmr:::.glsBlocks(fx$sigma)
    expect_setequal(map(blocks, sort), list(1:3, 4:5))
    plain <- pecotmr:::.glsWhiten(fx$y, fx$X, fx$sigma, regularize = FALSE)
    Linv <- solve(t(chol(as.matrix(fx$sigma))))
    one <- Linv %*% rep(1, 5)
    expect_equal(
        unname(plain$y),
        as.numeric(lm.fit(one, Linv %*% fx$y)$residuals)
    )
})

test_that(".glsWhiten needs every gene of y in the correlation and X", {
    fx <- .gc_whitenFixture()
    y <- c(fx$y, z = 1)
    expect_error(
        pecotmr:::.glsWhiten(y, fx$X, fx$sigma),
        "1 gene\\(s\\) of `y` are not in the gene correlation: z"
    )
    expect_error(
        pecotmr:::.glsWhiten(fx$y, fx$X[1:4, ], fx$sigma),
        "not in `X`"
    )
    expect_error(
        pecotmr:::.glsWhiten(unname(fx$y), fx$X, fx$sigma),
        "must be named by gene"
    )
})

# ---------------------------------------------------------------------------
# Variance components
# ---------------------------------------------------------------------------

# rho = a_g'(R o R)a_h / sqrt(...) for weight columns `A`, zero for genes
# that share no GWAS row.
.gc_quadraticExpected <- function(A) {
    C <- crossprod(A, .gc_R^2 %*% A)
    C / sqrt(outer(diag(C), diag(C)))
}

# .gc_fineMapping() with g2's row holding no variants, so no PIPs.
.gc_noPipFineMapping <- function() {
    rows <- map(c("g1", "g3"), function(name) {
        ids <- .gc_ids[.gc_mrGenes[[name]]$idx]
        FineMappingRow(ids, list(), data.frame(variant_id = ids, pip = 0.9))
    })
    empty <- FineMappingRow(
        character(0),
        list(),
        data.frame(variant_id = character(0), pip = numeric(0))
    )
    QtlFineMappingResult(
        studyName = rep("Q1", 3),
        context = rep("c1", 3),
        trait = c("g1", "g2", "g3"),
        method = rep("susie", 3),
        entry = list(rows[[1L]], empty, rows[[2L]])
    )
}

test_that("variance-component correlation is a'(R o R)a normalized", {
    skip_if_not_installed("CompQuadForm")
    gwas <- .gc_gwas()
    tw <- .gc_twasWeights()
    vc <- geneVarianceComponentTest(
        gwas,
        twasWeights = tw,
        variantWeights = "squaredWeight"
    )
    out <- geneCorrelation(
        vc,
        gwas,
        twasWeights = tw,
        analysis = "varianceComponent"
    )
    expected <- .gc_quadraticExpected(.gc_weightMatrix()^2)
    expected["g1", "g4"] <- expected["g4", "g1"] <- 0
    expected["g2", "g4"] <- expected["g4", "g2"] <- 0
    expect_equal(as.matrix(out$G1), expected)
    # The probit's correlation is taken as the same one, as MAGMA takes it.
    probit <- geneCorrelation(
        vc,
        gwas,
        twasWeights = tw,
        analysis = "varianceComponent",
        scale = "probit"
    )
    expect_equal(as.matrix(probit$G1), expected)
})

test_that("PIP-weighted correlation reads the fine-mapping windows", {
    skip_if_not_installed("CompQuadForm")
    gwas <- .gc_gwas()
    fmr <- .gc_fineMapping()
    vc <- geneVarianceComponentTest(gwas, fineMappingResult = fmr)
    out <- geneCorrelation(
        vc,
        gwas,
        fineMappingResult = fmr,
        analysis = "varianceComponent"
    )
    A <- vapply(
        .gc_mrGenes,
        function(gene) {
            pip <- seq(0.9, 0.6, length.out = length(gene$idx))
            replace(numeric(.gc_p), gene$idx, pip)
        },
        numeric(.gc_p)
    )
    expect_equal(as.matrix(out$G1), .gc_quadraticExpected(A))
    expect_error(
        geneCorrelation(vc, gwas, analysis = "varianceComponent"),
        "\"pip\" rows need the `fineMappingResult`"
    )
})

test_that("each analysis takes its own scales and its own result", {
    skip_if_not_installed("CompQuadForm")
    gwas <- .gc_gwas()
    tw <- .gc_twasWeights()
    vc <- geneVarianceComponentTest(
        gwas,
        twasWeights = tw,
        variantWeights = "squaredWeight"
    )
    expect_error(
        geneCorrelation(
            vc,
            gwas,
            twasWeights = tw,
            analysis = "varianceComponent",
            scale = "oneSided"
        ),
        "`scale` for varianceComponent statistics must be one of"
    )
    res <- .gc_run(gwas, tw)
    expect_error(
        geneCorrelation(res, gwas, twasWeights = tw, scale = "probit"),
        "`scale` for twas statistics must be one of"
    )
    expect_error(
        geneCorrelation(
            res,
            gwas,
            twasWeights = tw,
            analysis = "varianceComponent"
        ),
        "not a geneVarianceComponentTest\\(\\) result"
    )
    set.seed(99)
    other <- .gc_gwas(z = rnorm(.gc_p))
    expect_error(
        geneCorrelation(
            vc,
            other,
            twasWeights = tw,
            analysis = "varianceComponent"
        ),
        "does not match the one these inputs give"
    )
})

test_that("variance-component inputs other than x's are refused", {
    skip_if_not_installed("CompQuadForm")
    gwas <- .gc_gwas()
    tw <- .gc_twasWeights()
    fmr <- .gc_fineMapping()
    squared <- geneVarianceComponentTest(
        gwas,
        twasWeights = tw,
        variantWeights = "squaredWeight"
    )
    pips <- geneVarianceComponentTest(gwas, fineMappingResult = fmr)
    expect_error(
        geneCorrelation(squared, gwas, analysis = "varianceComponent"),
        "\"squaredWeight\" rows need the `twasWeights` \\(or the"
    )
    expect_error(
        geneCorrelation(
            c(squared[squared$trait == "g4"], pips),
            gwas,
            twasWeights = tw,
            fineMappingResult = fmr,
            analysis = "varianceComponent"
        ),
        "more than one kind of variant weight"
    )
    # g2's weights name its variants with their alleles swapped.
    expect_error(
        geneCorrelation(
            squared,
            gwas,
            twasWeights = tw,
            analysis = "varianceComponent",
            alleleFlip = FALSE
        ),
        "trait=g2, method=lasso\\) x GWAS 'G1' does not match"
    )
    expect_warning(
        expect_error(
            geneCorrelation(
                pips,
                gwas,
                fineMappingResult = .gc_noPipFineMapping(),
                analysis = "varianceComponent"
            ),
            "trait=g2, method=susie\\) x GWAS 'G1' does not match"
        ),
        "trait='g2', method='susie'\\) holds no PIPs"
    )
})

test_that("variance components are prioritized on either response", {
    skip_if_not_installed("CompQuadForm")
    inp <- .gp_inputs()
    .gp_useSharedCorrelation()
    local_mocked_bindings(.gpFit = .gp_crossFit, .package = "pecotmr")
    vc <- .gp_sharedVc
    for (scale in c("standardized", "probit")) {
        out <- genePrioritization(
            vc,
            inp$features,
            inp$gwas,
            twasWeights = inp$tw,
            analysis = "varianceComponent",
            scale = scale
        )
        column <- if (scale == "probit") "vcZ" else "vcStandardized"
        y <- set_names(S4Vectors::mcols(vc)[[column]], vc$trait)
        scores <- out$G1$scores
        scored <- names(scores)[!is.na(scores$y)]
        expect_setequal(scored, vc$trait)
        expect_equal(scores[scored]$y, unname(y[scored]))
    }
})

test_that("variance-component covariates count variants and eigenvalues", {
    study <- list(
        genes = list(
            list(
                gene = "a",
                coef = c(1, 0.5),
                vcMean = 1.5,
                variantWeights = "squaredWeight",
                cellKey = "ka"
            ),
            list(
                gene = "b",
                coef = c(1, 1, 1),
                vcMean = 3,
                variantWeights = "squaredWeight",
                cellKey = "kb"
            )
        ),
        variance = c(1.25, 4)
    )
    expect_equal(
        pecotmr:::.gpCovariates(
            study,
            "varianceComponent",
            c(ka = 0.3, kb = 0.1)
        ),
        cbind(
            nVariants = c(a = 2, b = 3),
            logNVariants = log(c(2, 3)),
            effectiveDf = c(1.8, 2.25),
            logEffectiveDf = log(c(1.8, 2.25)),
            cvRsq = c(0.3, 0.1)
        )
    )
    pips <- list_modify(
        study,
        genes = map(study$genes, list_modify, variantWeights = "pip")
    )
    expect_equal(
        colnames(pecotmr:::.gpCovariates(pips, "varianceComponent", NULL)),
        c("nVariants", "logNVariants", "effectiveDf", "logEffectiveDf")
    )
})

# ---------------------------------------------------------------------------
# Gene-set enrichment
# ---------------------------------------------------------------------------

test_that("a set's effect is its GLS slope with covariates out of both", {
    fx <- .gc_whitenFixture()
    sets <- cbind(
        s1 = c(1, 0, 1, 1, 0),
        s2 = c(0.5, 2, -1, 0, 1),
        all = 1
    )
    rownames(sets) <- names(fx$y)
    fit <- pecotmr:::.gseNullFit(fx$y, fx$sigma, fx$covariates)
    out <- pecotmr:::.gseChunk(1L, sets, fit)
    dense <- as.matrix(fx$sigma)
    ridge <- abs(min(0, eigen(dense)$values)) +
        0.05 +
        0.9 * max(0, var(fx$y) - 1)
    Linv <- solve(t(chol(dense + diag(ridge, 5))))
    nuisance <- Linv %*% cbind(1, fx$covariates)
    for (s in c("s1", "s2")) {
        lm <- summary(stats::lm(
            Linv %*% fx$y ~ 0 + nuisance + I(Linv %*% sets[, s])
        ))$coefficients
        row <- out[out$set == s, ]
        expect_equal(row$beta, lm[3, "Estimate"])
        expect_equal(row$se, lm[3, "Std. Error"])
        expect_equal(row$t, lm[3, "t value"])
        expect_equal(row$df, 2L)
        expect_equal(row$pval, pt(lm[3, "t value"], 2, lower.tail = FALSE))
    }
    expect_equal(out$nGenes, c(3L, 4L, 5L))
    # A set of every gene is the intercept: nothing is left to estimate.
    expect_true(is.na(out$beta[out$set == "all"]))
    expect_true(is.na(out$pval[out$set == "all"]))
    # Five genes, the intercept and three covariates: no df left for a set.
    set.seed(6)
    crowded <- matrix(rnorm(15), 5, dimnames = list(names(fx$y), NULL))
    expect_error(
        pecotmr:::.gseNullFit(fx$y, fx$sigma, crowded),
        "5 gene\\(s\\) with a statistic leave no residual degrees of freedom"
    )
})

test_that("geneSetEnrichment tests each set against each study", {
    skip_if_not_installed("CompQuadForm")
    inp <- .gp_inputs()
    .gp_useSharedCorrelation()
    out <- geneSetEnrichment(
        inp$res,
        inp$features,
        inp$gwas,
        twasWeights = inp$tw
    )
    expect_named(out, "G1")
    table <- out$G1
    expect_named(
        table,
        c("set", "nGenes", "beta", "se", "t", "df", "pval")
    )
    expect_equal(table$set, colnames(SummarizedExperiment::assay(inp$features)))
    # Every gene with a statistic, less the set and the nuisance's rank: the
    # intercept and three covariates, as every gene has four nonzero weights
    # (nVariants and its log are the intercept again).
    expect_equal(unique(table$df), 2L * .gp_nGenes - 4L - 1L)
    quadratic <- geneSetEnrichment(
        .gp_sharedVc,
        inp$features,
        inp$gwas,
        twasWeights = inp$tw,
        analysis = "varianceComponent",
        scale = "probit"
    )
    expect_true(all(is.finite(quadratic$G1$pval)))
    unnamed <- inp$features
    colnames(unnamed) <- NULL
    expect_error(
        geneSetEnrichment(inp$res, unnamed, inp$gwas, twasWeights = inp$tw),
        "must name its sets"
    )
})

# ===========================================================================
# geneVarianceComponentTest
# ===========================================================================

# ===========================================================================
# Fixture: 12 variants on chr1 with AR(1) LD in one LdMatrix, one GWAS row,
# and two genes. g1's weights sit on variants 1-6; g2's on 5-12, named with
# their alleles swapped, so they match the GWAS flipped.
# ===========================================================================

.vc_p <- 12L
.vc_R <- 0.5^abs(outer(seq_len(.vc_p), seq_len(.vc_p), "-"))
.vc_variants <- .testLdVariants(.vc_p, "chr1", 1000L)
.vc_ids <- .vc_variants$SNP
.vc_flipped <- sprintf("chr1:%d:T:C", .vc_variants$pos)
.vc_z <- c(3, -2.5, 1, 0.4, -1.2, 2, 0.3, -0.8, 1.5, -0.2, 0.9, 2.2)

.vc_genes <- list(
    g1 = list(idx = 1:6, w = c(0.4, -0.2, 0, 0.3, 0.1, -0.5)),
    g2 = list(idx = 5:12, w = c(0.2, 0.1, -0.3, 0.4, 0, 0.2, 0.1, -0.1))
)

.vc_gwas <- function(z = .vc_z, qcOptions = list(qc = "done")) {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(.vc_variants$pos, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = .vc_ids,
        A1 = "T",
        A2 = "C",
        Z = z
    )
    GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = LdMatrix(.vc_R, variants = .vc_variants),
        qcOptions = qcOptions
    )
}

.vc_geneIds <- function(name) {
    idx <- .vc_genes[[name]]$idx
    if (name == "g2") .vc_flipped[idx] else .vc_ids[idx]
}

.vc_twasWeights <- function() {
    TwasWeights(
        studyName = rep("Q1", 2),
        context = rep("c1", 2),
        trait = names(.vc_genes),
        method = rep("lasso", 2),
        entry = map2(
            map(names(.vc_genes), .vc_geneIds),
            map(.vc_genes, "w"),
            TwasWeightsRow
        )
    )
}

# PIPs over each gene's whole window, the squared weights rescaled.
.vc_pips <- function(name) {
    w2 <- .vc_genes[[name]]$w^2
    w2 / sum(w2) + 0.01
}

.vc_fineMapping <- function() {
    rows <- map(names(.vc_genes), function(name) {
        FineMappingRow(
            variantIds = .vc_geneIds(name),
            susieFit = list(),
            topLoci = data.frame(
                variant_id = .vc_geneIds(name),
                pip = .vc_pips(name)
            )
        )
    })
    QtlFineMappingResult(
        studyName = rep("Q1", 2),
        context = rep("c1", 2),
        trait = names(.vc_genes),
        method = rep("susie", 2),
        entry = rows
    )
}

# Q, its null eigenvalues and p-value for weights `a` on variants `idx`.
.vc_expected <- function(a, idx) {
    keep <- a > 0
    a <- a[keep]
    idx <- idx[keep]
    root <- sqrt(a)
    lambda <- eigen(
        outer(root, root) * .vc_R[idx, idx],
        symmetric = TRUE
    )$values
    q <- sum(a * .vc_z[idx]^2)
    # Farebrother's to 1e-15: the test's Davies inversion, to an absolute
    # 1e-8, is compared with it at a matching tolerance.
    p <- CompQuadForm::farebrother(q, lambda, eps = 1e-15)$Qq
    list(
        nVariants = length(a),
        vcStatistic = q,
        vcMean = sum(lambda),
        vcSd = sqrt(2 * sum(lambda^2)),
        vcStandardized = (q - sum(lambda)) / sqrt(2 * sum(lambda^2)),
        vcPval = p,
        vcZ = qnorm(p, lower.tail = FALSE)
    )
}

.vc_columns <- function(res, gene) {
    mc <- as.data.frame(S4Vectors::mcols(res))
    as.list(mc[mc$trait == gene, pecotmr:::.VC_STATISTICS])
}

# The variance components of `tw`'s squared weights on GWAS `gr` with LD
# reference `ldReference`, and their correlation.
.vc_runOn <- function(ldReference, gr, tw) {
    gwas <- GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg38",
        ldReference = ldReference,
        qcOptions = list(qc = "done")
    )
    vc <- geneVarianceComponentTest(
        gwas,
        twasWeights = tw,
        variantWeights = "squaredWeight"
    )
    correlation <- geneCorrelation(
        vc,
        gwas,
        twasWeights = tw,
        analysis = "varianceComponent"
    )
    list(
        statistics = as.data.frame(
            S4Vectors::mcols(vc)
        )[pecotmr:::.VC_STATISTICS],
        correlation = as.matrix(correlation$G1)
    )
}

test_that("squared weights give Q = sum(w^2 z^2) and its exact null", {
    skip_if_not_installed("CompQuadForm")
    res <- geneVarianceComponentTest(
        .vc_gwas(),
        twasWeights = .vc_twasWeights(),
        variantWeights = "squaredWeight"
    )
    expect_s4_class(res, "GRanges")
    expect_equal(res$trait, c("g1", "g2"))
    expect_equal(unique(res$analysis), "varianceComponent")
    expect_equal(unique(res$variantWeights), "squaredWeight")
    expect_true(all(is.na(res$pipFloor)))
    # Zero weights do not enter; g2's allele swap leaves its squares alone.
    for (gene in names(.vc_genes)) {
        g <- .vc_genes[[gene]]
        expect_equal(
            .vc_columns(res, gene),
            .vc_expected(g$w^2, g$idx),
            tolerance = 1e-5
        )
    }
    expect_equal(GenomeInfoDb::genome(res)[["chr1"]], "hg19")
})

test_that("PIPs weight every variant of the fine-mapping window", {
    skip_if_not_installed("CompQuadForm")
    fmr <- .vc_fineMapping()
    res <- geneVarianceComponentTest(.vc_gwas(), fineMappingResult = fmr)
    expect_equal(unique(res$variantWeights), "pip")
    for (gene in names(.vc_genes)) {
        expect_equal(
            .vc_columns(res, gene),
            .vc_expected(.vc_pips(gene), .vc_genes[[gene]]$idx),
            tolerance = 1e-5
        )
    }
    floored <- geneVarianceComponentTest(
        .vc_gwas(),
        fineMappingResult = fmr,
        variantWeights = "flooredPip",
        pipFloor = 0.05
    )
    expect_equal(unique(floored$pipFloor), 0.05)
    expect_equal(
        .vc_columns(floored, "g1"),
        .vc_expected(.vc_pips("g1") + 0.05, .vc_genes$g1$idx),
        tolerance = 1e-5
    )
})

test_that("the arguments are checked against the weights they choose", {
    gwas <- .vc_gwas()
    tw <- .vc_twasWeights()
    fmr <- .vc_fineMapping()
    expect_error(
        geneVarianceComponentTest(gwas, twasWeights = tw),
        "are PIPs, read from a `fineMappingResult`"
    )
    expect_error(
        geneVarianceComponentTest(
            gwas,
            twasWeights = tw,
            fineMappingResult = fmr
        ),
        "`twasWeights` is not read"
    )
    expect_error(
        geneVarianceComponentTest(gwas, fineMappingResult = fmr, pipFloor = 1),
        "applies to `variantWeights = \"flooredPip\"` only"
    )
    expect_error(
        geneVarianceComponentTest(
            gwas,
            fineMappingResult = fmr,
            variantWeights = "flooredPip"
        ),
        "needs `pipFloor`"
    )
    expect_error(
        geneVarianceComponentTest(
            gwas,
            fineMappingResult = fmr,
            variantWeights = "flooredPip",
            pipFloor = 0
        ),
        "must be positive"
    )
    expect_error(
        geneVarianceComponentTest(gwas, variantWeights = "squaredWeight"),
        "squared weights need"
    )
    expect_error(
        geneVarianceComponentTest(
            .vc_gwas(qcOptions = list()),
            fineMappingResult = fmr
        ),
        "gwasSumStats has no QC record"
    )
})

test_that("weight sets with no PIPs, or no one weight vector, are skipped", {
    skip_if_not_installed("CompQuadForm")
    gwas <- .gc_gwas()
    expect_warning(
        pips <- geneVarianceComponentTest(
            gwas,
            fineMappingResult = .gc_noPipFineMapping()
        ),
        "trait='g2', method='susie'\\) holds no PIPs; skipping"
    )
    expect_equal(pips$trait, c("g1", "g3"))
    twoColumns <- TwasWeights(
        rep("Q1", 2),
        rep("c1", 2),
        c("g1", "g2"),
        rep("lasso", 2),
        list(
            TwasWeightsRow(.gc_ids[1:3], matrix(1:6 / 10, 3)),
            TwasWeightsRow(.gc_ids[4:6], c(0.1, 0.2, 0.3))
        )
    )
    expect_warning(
        squared <- geneVarianceComponentTest(
            gwas,
            twasWeights = twoColumns,
            variantWeights = "squaredWeight"
        ),
        "trait=g1, method=lasso\\) are a matrix with 2 columns"
    )
    expect_equal(squared$trait, "g2")
})

test_that("a genotype panel's LD gives what its LdMatrix gives", {
    skip_if_not_installed("CompQuadForm")
    skip_if_not_installed("PlinkMatrix")
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    panel <- panel[seq_len(.vc_p), ]
    ranges <- SummarizedExperiment::rowRanges(panel)
    ids <- names(ranges)
    gr <- GenomicRanges::GRanges(
        GenomicRanges::seqnames(ranges),
        IRanges::IRanges(GenomicRanges::start(ranges), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = ids,
        A1 = ranges$A1,
        A2 = ranges$A2,
        Z = .vc_z
    )
    tw <- TwasWeights(
        rep("Q1", 2),
        rep("c1", 2),
        names(.vc_genes),
        rep("lasso", 2),
        unname(map(.vc_genes, function(g) TwasWeightsRow(ids[g$idx], g$w)))
    )
    ld <- LdMatrix(unname(computeLd(panel)), variants = ranges)
    runs <- map(list(panel = panel, ld = ld), .vc_runOn, gr = gr, tw = tw)
    expect_equal(
        runs$panel$statistics,
        runs$ld$statistics,
        tolerance = 1e-6
    )
    expect_equal(runs$panel$correlation, runs$ld$correlation)
})

test_that("weights matching no GWAS variant test nothing, and say so", {
    tw <- TwasWeights(
        "Q1",
        "c1",
        "g9",
        "lasso",
        list(TwasWeightsRow(c("chr1:1000:A:C", "chr1:1100:A:C"), c(1, 1)))
    )
    expect_message(
        res <- geneVarianceComponentTest(
            .vc_gwas(),
            twasWeights = tw,
            variantWeights = "squaredWeight"
        ),
        "nothing was tested"
    )
    expect_length(res, 0L)
    expect_true(all(
        c("variantWeights", pecotmr:::.VC_STATISTICS) %in%
            names(S4Vectors::mcols(res))
    ))
})

test_that("a genotype factor's eigenvalues are the LD matrix's", {
    set.seed(4)
    X <- scale(matrix(rbinom(40 * 9, 2, 0.3), 40, 9))
    R <- crossprod(X) / 39
    a <- runif(9)
    fromMatrix <- pecotmr:::.vcEigenvalues(list(R = R), a)
    # More samples than variants, then fewer.
    expect_equal(pecotmr:::.vcEigenvalues(list(ldFactor = X), a), fromMatrix)
    few <- X[1:6, ]
    expect_equal(
        pecotmr:::.vcEigenvalues(list(ldFactor = few), a),
        pecotmr:::.vcEigenvalues(list(R = crossprod(few) / 5), a)
    )
    expect_length(pecotmr:::.vcEigenvalues(list(ldFactor = few), a), 6L)
})

test_that("the tail inverts numerically, and by saddlepoint past 1e-12", {
    skip_if_not_installed("CompQuadForm")
    skip_if_not_installed("survey")
    lambda <- rep(1, 5)
    expect_equal(pecotmr:::.quadFormUpperTail(0, lambda), 1)
    expect_equal(pecotmr:::.quadFormUpperTail(3, numeric(0)), 1)
    # Davies' inversion for p of 1e-5 or more (q = 2, 10), to an absolute
    # 1e-8; Farebrother's below it (q = 40, p = 1.5e-7).
    bulk <- c(2, 10, 40)
    expect_equal(
        map_dbl(bulk, pecotmr:::.quadFormUpperTail, lambda = lambda),
        pchisq(bulk, 5, lower.tail = FALSE),
        tolerance = 1e-6
    )
    # Past 1e-12 the saddlepoint keeps its relative error to a few percent
    # for five equal weights, where Farebrother returns 0.
    tail <- c(80, 150, 400)
    expect_equal(
        map_dbl(tail, pecotmr:::.quadFormUpperTail, lambda = lambda) /
            pchisq(tail, 5, lower.tail = FALSE),
        rep(1, 3),
        tolerance = 0.05
    )
})

test_that("squared cross forms are B'(R o R)B over several stripes", {
    set.seed(2)
    p <- 300L
    R <- 0.7^abs(outer(seq_len(p), seq_len(p), "-"))
    B <- matrix(runif(p * 3), p, 3)
    expect_equal(
        pecotmr:::.ldSquaredCrossForms(list(R = R), B),
        crossprod(B, R^2 %*% B)
    )
})
