# Tests for the GenotypeMixture class (R/GenotypeMixture.R) and the LD
# reference dispatch it plugs into.

# Two "cohorts": the two halves of one real genotype panel's samples.
.gmPanels <- function() {
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    list(a = panel[, 1:40], b = panel[, 41:100], full = panel)
}

test_that("GenotypeMixture aligns its panels to the variants they share", {
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.3, 0.7))
    expect_s4_class(mix, "GenotypeMixture")
    expect_s4_class(mix, "MultiAssayExperiment")
    expect_equal(pecotmr:::.ldRefNVariants(mix), nrow(p$full))
    expect_equal(mixtureWeights(mix), c(a = 0.3, b = 0.7))
    comps <- as.list(MultiAssayExperiment::experiments(mix))
    expect_named(comps, c("a", "b"))
    expect_equal(map_int(comps, ncol), c(a = 40L, b = 60L))
    expect_equal(map_int(comps, nrow), c(a = nrow(p$full), b = nrow(p$full)))
    # The canonical ids the shared ranges render, not the panels' labels.
    expect_identical(
        variantIds(mix),
        pecotmr:::.grVariantIds(SummarizedExperiment::rowRanges(p$full))
    )
    expect_true(all(str_detect(variantIds(mix), "^chr21:[0-9]+:[ACGT]+:")))
    # Each panel records its allele coding against the first's; these two
    # halves of one panel code every variant alike.
    expect_true(all(pecotmr:::.gmSign(mix) == 1))
    expect_equal(dim(pecotmr:::.gmSign(mix)), c(nrow(p$full), 2L))
    expect_true(validObject(mix))
    out <- capture.output(show(mix))
    expect_true(any(grepl("from 2 genotype panels", out)))
    expect_true(any(grepl("panels a, b", out)))
})

test_that("a mixture names its samples by cohort", {
    # Cohorts can reuse sample ids (every sketch does), so the primary names
    # carry the cohort; the panels keep their own column names.
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$a), mixtureWeights = c(0.5, 0.5))
    sm <- MultiAssayExperiment::sampleMap(mix)
    ids <- colnames(p$a)
    expect_equal(
        as.character(sm$primary),
        c(str_c("a_", ids), str_c("b_", ids))
    )
    expect_equal(as.character(sm$colname), c(ids, ids))
    cd <- MultiAssayExperiment::colData(mix)
    expect_equal(rownames(cd), as.character(sm$primary))
    expect_equal(cd$cohort, rep(c("a", "b"), each = length(ids)))
})

test_that("mixture weights are matched to panels by name", {
    p <- .gmPanels()
    byName <- GenotypeMixture(
        list(a = p$a, b = p$b),
        mixtureWeights = c(b = 0.7, a = 0.3)
    )
    expect_equal(mixtureWeights(byName), c(a = 0.3, b = 0.7))
    expect_error(
        GenotypeMixture(
            list(a = p$a, b = p$b),
            mixtureWeights = c(a = 0.3, c = 0.7)
        ),
        "must name each panel once: a, b"
    )
})

test_that("a mixture's LD is the weighted average of its panels' LD", {
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.3, 0.7))
    idx <- c(5L, 1L, 20L, 7L)
    ldA <- unname(computeLd(p$a, snpIdx = idx))
    ldB <- unname(computeLd(p$b, snpIdx = idx))
    want <- 0.3 * ldA + 0.7 * ldB
    expect_equal(pecotmr:::.ldRefLdAt(mix, idx), want)
    ids <- names(SummarizedExperiment::rowRanges(p$full))[idx]
    expect_equal(unname(pecotmr:::.ldFromRef(mix, ids, label = "t")), want)
})

test_that("a mixture's factor reproduces its LD exactly", {
    # Each panel standardized on its own and weighted: the stacked factor's
    # crossproduct is (N - 1) times the mixture's LD.
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.3, 0.7))
    idx <- 1:12
    f <- pecotmr:::.ldRefFactor(mix, idx)
    expect_equal(nrow(f), 100L)
    expect_equal(
        unname(crossprod(f)) / (nrow(f) - 1),
        pecotmr:::.ldRefLdAt(mix, idx)
    )
})

test_that("a panel-monomorphic variant is padded with a zero-mean pair", {
    # The pair adds the panel's weight to the variant's own variance and
    # nothing to its mean or its LD with anything else.
    block <- cbind(c(1, -1, 0.5, -0.5), 0, c(-1, 1, 0, 0))
    rows <- pecotmr:::.gmMonomorphicRows(block, 0.3)
    expect_equal(dim(rows), c(2L, 3L))
    expect_equal(colSums(rows), rep(0, 3))
    expect_equal(crossprod(rows), diag(c(0, 0.3, 0)))
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.3, 0.7))
    f <- pecotmr:::.ldRefFactor(mix, 1:12)
    expect_equal(unname(colMeans(f)), rep(0, 12))
})

test_that("a mixture answers sizes and statistics from its panels", {
    p <- .gmPanels()
    w <- c(0.3, 0.7)
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = w)
    # Noise n: one panel whose LD is as noisy; count n: the pooled samples.
    expect_equal(
        pecotmr:::.ldRefNoiseN(mix),
        as.integer(round(1 / sum(w^2 / c(40, 60))))
    )
    expect_equal(pecotmr:::.ldRefRealNSamples(mix), 100L)
    # A variant's MAF is its smallest in any panel, its missingness the
    # largest, so a cutoff drops a variant that fails anywhere.
    st <- pecotmr:::.ldRefVariantStats(mix, 1:10, TRUE, "t")
    sa <- pecotmr:::.ldRefVariantStats(p$a, 1:10, TRUE, "t")
    sb <- pecotmr:::.ldRefVariantStats(p$b, 1:10, TRUE, "t")
    expect_equal(st$maf, pmin(sa$maf, sb$maf))
    expect_equal(st$missRate, pmax(sa$missRate, sb$missRate))
    expect_equal(pecotmr:::.ldRefMacSamples(mix, 5), 40L)
    # The genotype variance is the weighted within-panel average.
    vA <- pecotmr:::.ldRefVariance(p$a, 1:5)
    vB <- pecotmr:::.ldRefVariance(p$b, 1:5)
    expect_equal(pecotmr:::.ldRefVariance(mix, 1:5), w[[1]] * vA + w[[2]] * vB)
    # Real panels record no frequencies, so neither does their mixture.
    expect_null(pecotmr:::.ldRefAf(mix))
})

test_that("subsetting a mixture's rows subsets every panel with it", {
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.5, 0.5))
    sub <- mix[c(3L, 1L), ]
    expect_s4_class(sub, "GenotypeMixture")
    expect_equal(pecotmr:::.ldRefNVariants(sub), 2L)
    comps <- as.list(MultiAssayExperiment::experiments(sub))
    expect_equal(map_int(comps, nrow), c(a = 2L, b = 2L))
    expect_equal(dim(pecotmr:::.gmSign(sub)), c(2L, 2L))
    expect_equal(
        pecotmr:::.ldRefLdAt(sub, 1:2),
        pecotmr:::.ldRefLdAt(mix, c(3L, 1L))
    )
    expect_true(validObject(sub))
    # By variant id and by range, every panel is cut to the same variants.
    ids <- names(SummarizedExperiment::rowRanges(p$full))[c(3L, 1L)]
    expect_equal(variantIds(mix[ids, ]), variantIds(sub))
    gr <- SummarizedExperiment::rowRanges(p$full)[2:4]
    expect_equal(pecotmr:::.ldRefNVariants(mix[gr, ]), 3L)
    # A range is turned into positions like any other subscript: one naming
    # no variant empties the panels as positions do, seeds and all.
    none <- mix[GenomicRanges::GRanges("chr21", IRanges::IRanges(1L, 1L)), ]
    expect_equal(pecotmr:::.ldRefNVariants(none), 0L)
    expect_true(validObject(none))
    expect_true(all(map_lgl(
        map(
            as.list(MultiAssayExperiment::experiments(none)),
            pecotmr:::.panelSeed
        ),
        is.null
    )))
    # One index per panel could cut the panels to different variants, so
    # such subscripts are refused, as is anything else that names no rows.
    refused <- "subset by position, logical mask, variant id or GRanges"
    expect_error(mix[list(a = 1:3, b = 1:2), ], refused)
    expect_error(mix[IRanges::IntegerList(a = 1:3, b = 1:2), ], refused)
    expect_error(mix[factor("x"), ], refused)
    # Emptying drops every panel's variants, and their seeds with them, but
    # keeps the samples the sampleMap names and what each panel was read
    # from.
    empty <- pecotmr:::.emptyLdRef(mix)
    expect_equal(pecotmr:::.ldRefNVariants(empty), 0L)
    expect_true(validObject(empty))
    emptied <- as.list(MultiAssayExperiment::experiments(empty))
    expect_equal(map_int(emptied, ncol), c(a = 40L, b = 60L))
    expect_true(all(map_lgl(map(emptied, pecotmr:::.panelSeed), is.null)))
    expect_equal(
        map_chr(emptied, pecotmr:::.ldRefFormat),
        c(a = "plink1", b = "plink1")
    )
})

test_that("subsetting a mixture's samples narrows its panels", {
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.5, 0.5))
    keep <- c(
        str_c("a_", colnames(p$a)[1:10]),
        str_c("b_", colnames(p$b)[1:20])
    )
    sub <- mix[, keep]
    comps <- as.list(MultiAssayExperiment::experiments(sub))
    expect_equal(map_int(comps, ncol), c(a = 10L, b = 20L))
    want <- 0.5 *
        unname(computeLd(p$a[, 1:10], snpIdx = 1:5)) +
        0.5 * unname(computeLd(p$b[, 1:20], snpIdx = 1:5))
    expect_equal(pecotmr:::.ldRefLdAt(sub, 1:5), want)
})

test_that("a mixture's panels cannot be dropped", {
    # The weights describe every panel, so dropping one would leave them
    # describing a mixture that is no longer there.
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.5, 0.5))
    expect_error(mix[,, "a"], "panels cannot be subset")
})

test_that("a mixture is checked against other references by its panels", {
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.5, 0.5))
    # Two QC'd subsets of one mixture pass, warning only that they overlap
    # in part.
    rlang::reset_warning_verbosity("pecotmrLdReferenceOverlap-pipe")
    expect_warning(
        pecotmr:::.requireMatchingLdRefs(mix[1:200], mix[100:349], "pipe"),
        "share 101 variant"
    )
    reweighted <- GenotypeMixture(
        list(a = p$a, b = p$b),
        mixtureWeights = c(0.4, 0.6)
    )
    expect_error(
        pecotmr:::.requireMatchingLdRefs(mix, reweighted, "pipe"),
        "average different panels or weight them differently"
    )
    expect_error(
        pecotmr:::.requireMatchingLdRefs(mix, p$full, "pipe"),
        "is a GenotypeMixture and the other's a genotype panel"
    )
})

test_that("raiss imputes from a mixture through the SVD route", {
    p <- .gmPanels()
    mix <- GenotypeMixture(
        list(a = p$a[1:30, ], b = p$b[1:30, ]),
        mixtureWeights = c(0.4, 0.6)
    )
    set.seed(2)
    ref <- pecotmr:::.refVariantsFromLdRef(mix)
    known <- ref[seq_len(30L) %% 3L != 0L, ]
    known$z <- rnorm(nrow(known))
    viaGenotypes <- raiss(
        knownZscores = known,
        genotypes = mix,
        svdTol = 1e-12,
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    # Both routes drop the same weak eigen-directions (`rcond`), so at the
    # defaults they are the same imputation.
    viaLd <- raiss(
        ref,
        known,
        ldMatrix = pecotmr:::.ldRefLdAt(mix, seq_len(30L)),
        r2Threshold = 0,
        minimumLd = 0,
        verbose = FALSE
    )
    expect_equal(
        viaGenotypes$resultNofilter$z,
        viaLd$resultNofilter$z,
        tolerance = 1e-8
    )
})

test_that("GenotypeMixture refuses bad inputs", {
    p <- .gmPanels()
    expect_error(
        GenotypeMixture(list(a = p$a), mixtureWeights = 1),
        "at least two genotype panels"
    )
    expect_error(
        GenotypeMixture(list(p$a, p$b), mixtureWeights = c(0.5, 0.5)),
        "must be a named list"
    )
    expect_error(
        GenotypeMixture(list(a = p$a, a = p$b), mixtureWeights = c(0.5, 0.5)),
        "must be a named list"
    )
    ld <- LdMatrix(
        diag(2),
        variants = data.frame(chrom = "1", pos = 1:2, A1 = "A", A2 = "G")
    )
    expect_error(
        GenotypeMixture(list(a = p$a, b = ld), mixtureWeights = c(0.5, 0.5)),
        "must be a genotype panel"
    )
    expect_error(
        GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.5, 0.6)),
        "summing to 1"
    )
    expect_error(
        GenotypeMixture(
            list(a = p$a[1:10, ], b = p$b[11:20, ]),
            mixtureWeights = c(0.5, 0.5)
        ),
        "share no variant to combine over"
    )
})

test_that("an experiment added with c() leaves the mixture invalid", {
    # The weights describe the panels the mixture was built over, so c()
    # cannot extend it: whatever is added, validity refuses the result.
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.5, 0.5))
    se <- SummarizedExperiment::SummarizedExperiment(
        matrix(1, 2L, 2L, dimnames = list(c("r1", "r2"), colnames(p$a)[1:2]))
    )
    expect_error(c(mix, extra = se), "at least two genotype panels")
    expect_error(
        c(mix, extra = p$b[1:5, ]),
        "every panel must hold exactly the mixture's variants"
    )
    expect_error(c(mix, extra = p$b), "must carry `mixtureSign`")
    expect_error(c(mix, extra = p$b), "one weight per panel, named after it")
    unnamed <- mix
    S4Vectors::metadata(unnamed)$mixtureWeights <- c(0.5, 0.5)
    expect_error(validObject(unnamed), "one weight per panel, named after it")
})

test_that("a mixture's allele frequency is its panels' weighted average", {
    # A stochastic sketch records its A1 frequencies; the second half gets
    # its own, so the weights matter.
    sk <- readGenotypes(test_path("test_data", "test_harmonize_regions.pgen"))
    half <- seq_len(ncol(sk) / 2)
    a <- sk[, half]
    b <- sk[, -half]
    rd <- SummarizedExperiment::rowData(b)
    rd$AF <- seq(0.1, 0.8, length.out = nrow(b))
    SummarizedExperiment::rowData(b) <- rd
    mix <- GenotypeMixture(list(a = a, b = b), mixtureWeights = c(0.25, 0.75))
    expect_equal(
        pecotmr:::.ldRefAf(mix),
        0.25 * pecotmr:::.ldRefAf(a) + 0.75 * rd$AF
    )
    # A panel coding a variant the other way round contributes the
    # frequency of the first panel's A1: the complement of its own.
    expect_equal(
        pecotmr:::.gmAlignedAf(
            2L,
            list(c(0.1, 0.2), c(0.3, 0.4)),
            cbind(c(1, 1), c(1, -1))
        ),
        c(0.3, 0.6)
    )
})

test_that("an LD statistic reads a mixture's blocks as LdMatrix ones", {
    p <- .gmPanels()
    mix <- GenotypeMixture(list(a = p$a, b = p$b), mixtureWeights = c(0.3, 0.7))
    gr <- SummarizedExperiment::rowRanges(p$full)
    pos <- GenomicRanges::start(gr)
    cuts <- list(1:100, 101:250)
    blocks <- GenomicRanges::GRanges(
        "chr21",
        IRanges::IRanges(
            c(pos[[1L]], pos[[101L]]),
            c(pos[[100L]], pos[[250L]])
        )
    )
    lazy <- buildLdScore(mix, ldBlocks = blocks, genome = "hg19")
    expect_identical(ldReference(lazy), mix)
    ld <- map(cuts, function(i) {
        LdMatrix(
            pecotmr:::.ldRefLdAt(mix, i),
            variants = gr[i],
            nSamples = ldPanelNSamples(lazy)
        )
    })
    stored <- buildLdScore(ld, genome = "hg19")
    expect_equal(
        S4Vectors::mcols(lazy)$ldScores,
        S4Vectors::mcols(stored)$ldScores,
        tolerance = 1e-10
    )
})
