# Fixtures shared by the LD-statistic tests (LdEigen, LdScore, LdStatistic),
# the heritability tests (h2EstimationWrappers, H2Estimate, h2Annotations)
# and test_genotypeHandle.R.

# =============================================================================
# Test data helpers
# =============================================================================

makeTestLdBlocks <- function() {
    blocks_gr <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(1, 5001), end = c(5000, 10000))
    )
    blocks_gr
}

makeTestSnpInfo <- function(n = 10) {
    data.frame(
        SNP = paste0("rs", seq_len(n)),
        CHR = rep("chr1", n),
        BP = seq(1000, by = 100, length.out = n),
        A1 = rep("A", n),
        A2 = rep("G", n),
        stringsAsFactors = FALSE
    )
}


# =============================================================================
# LdMatrix fixtures for buildLdEigen() / buildLdScore()
# =============================================================================

# One block's variants as LdMatrix() takes them: SNP ids, alleles and an AF
# column, at 100 bp spacing.
.testLdVariants <- function(n, chrom, startBp) {
    bp <- seq(startBp, by = 100L, length.out = n)
    data.frame(
        chrom = chrom,
        pos = as.integer(bp),
        SNP = paste0(chrom, ":", bp, ":C:T"),
        A1 = "T",
        A2 = "C",
        AF = seq(0.2, 0.8, length.out = n)
    )
}

.testArBlock <- function(n, rho) {
    rho^abs(outer(seq_len(n), seq_len(n), "-"))
}

# One LD block: an LdMatrix with AR(1) LD.
makeTestLdMatrix <- function(
    n = 6L,
    chrom = "chr1",
    startBp = 1000L,
    rho = 0.5,
    nSamples = 500L
) {
    LdMatrix(
        .testArBlock(n, rho),
        variants = .testLdVariants(n, chrom, startBp),
        nSamples = nSamples
    )
}

# Several consecutive LD blocks along one chromosome, one LdMatrix each.
makeTestLdMatrixBlocks <- function(
    sizes = c(4L, 3L),
    chrom = "chr1",
    startBp = 1000L,
    rho = 0.5,
    nSamples = 500L
) {
    starts <- startBp + 100L * c(0L, cumsum(head(sizes, -1L)))
    map2(
        sizes,
        starts,
        function(n, s) makeTestLdMatrix(n, chrom, s, rho, nSamples)
    )
}

# A genotype panel of 100 samples cut into three LD blocks of 60, 140 and 100
# variants, with the same blocks' LD as LdMatrix objects: a statistic read
# from the panel should agree with one built from the matrices. The middle
# block has more variants than the panel has samples, so its LD has 40 zero
# eigenvalues that the SVD of the panel's factor does not list.
.h2PanelFixture <- function() {
    skip_if_not_installed("PlinkMatrix")
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    gr <- SummarizedExperiment::rowRanges(panel)
    pos <- GenomicRanges::start(gr)
    cuts <- list(1:60, 61:200, 201:300)
    blocks <- GenomicRanges::GRanges(
        "chr21",
        IRanges::IRanges(
            start = map_int(cuts, function(i) pos[[min(i)]]),
            end = map_int(cuts, function(i) pos[[max(i)]])
        )
    )
    ld <- map(cuts, function(i) {
        LdMatrix(
            unname(computeLd(panel, snpIdx = i)),
            variants = gr[i],
            nSamples = ncol(panel)
        )
    })
    list(panel = panel, blocks = blocks, ld = ld)
}

# One genome-wide LDER estimate, stratified (with an enrichment table and a
# jackknife matrix) unless `with_enrichment = FALSE`.
makeTestH2Estimate <- function(with_enrichment = TRUE) {
    enrich <- if (with_enrichment) {
        data.frame(
            annotation = c("annot1", "annot2"),
            tau = c(1e-7, 2e-7),
            tauSe = c(5e-8, 6e-8),
            enrichment = c(2.0, 3.0),
            enrichmentSe = c(0.5, 0.7),
            enrichmentP = c(0.01, 0.001),
            propH2 = c(0.3, 0.5),
            propSnps = c(0.15, 0.17),
            stringsAsFactors = FALSE
        )
    } else {
        NULL
    }

    set.seed(55)
    annotationJackknifeCoefs <- if (with_enrichment) {
        matrix(
            rnorm(20),
            nrow = 10,
            ncol = 2,
            dimnames = list(NULL, c("annot1", "annot2"))
        )
    } else {
        NULL
    }

    .h2EstimateNew(
        list(GenomicRanges::GRanges()),
        list(
            study = "test",
            context = NA_character_,
            trait = NA_character_,
            h2 = 0.3,
            h2_se = 0.05,
            h2_lower = NA_real_,
            h2_upper = NA_real_,
            h2_liability = NA_real_,
            h2_liability_se = NA_real_,
            h2_liability_lower = NA_real_,
            h2_liability_upper = NA_real_,
            intercept = 1.01,
            intercept_se = 0.02,
            p_value = NA_real_,
            n_snps = 100L,
            n_components = NA_integer_,
            scale = NA_real_,
            reference_overlap = NA_real_
        ),
        list(
            enrichment = list(enrich),
            jackknife = list(annotationJackknifeCoefs),
            score_stats = list(NULL)
        ),
        GenomicRanges::GRanges(),
        .h2GenomeWideSettings("lder", with_enrichment, inSample = FALSE)
    )
}

# The published settings of a genome-wide estimator by `method`, for a
# stratified fit on `annotations` against `ref`: stratified as their design
# is, and LDER's taper from the reference's kind, as the estimators' own
# defaults take them.
.h2MethodSettings <- function(method, annotations, ref) {
    .h2GenomeWideSettings(
        method,
        .h2Stratified(annotations),
        identical(method, "lder") && inSample(ref)
    )
}

# =============================================================================
# GWAS summary statistics over a reference's variants (shared by the
# heritability and fine-mapping annotation tests)
# =============================================================================

makeTestGwasGenotypeHandle <- function() {
    new(
        "GenotypeHandle",
        path = "/tmp/test",
        format = "plink2",
        snpInfo = data.frame(),
        nSamples = 0L,
        sampleIds = character(),
        pgenPtr = NULL
    )
}

.dfToGwasGr <- function(df) {
    gr <- GenomicRanges::GRanges(
        seqnames = df$CHR,
        ranges = IRanges::IRanges(start = df$BP, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = df$SNP,
        A1 = df$A1,
        A2 = df$A2,
        Z = df$Z,
        N = df$N
    )
    gr
}

# ---- multi-study fixtures (estimateLocalH2 and estimateLocalRg) -------------

# Two studies of qtlDatasetExample's 165 donors (82 and 83), contexts a and b
# with polygenic effects shared across the studies (a with b at 0.6).
.twoStudyFixture <- function() {
    data(qtlDatasetExample, envir = environment())
    panel <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    brain <- pecotmr:::.qtlContextSe(qtlDatasetExample, "brain")
    G <- scale(t(as.matrix(SummarizedExperiment::assay(panel))))
    G[!is.finite(G)] <- 0
    set.seed(4L)
    b <- matrix(rnorm(2L * ncol(G)), ncol = 2L) %*%
        chol(matrix(c(1, 0.6, 0.6, 1), 2L))
    e <- matrix(rnorm(2L * nrow(G)), ncol = 2L)
    y <- scale(G %*% b) * sqrt(0.5) + e * sqrt(0.5)
    contexts <- stats::setNames(nm = c("a", "b"))
    full <- map(
        c(s1 = "s1", s2 = "s2"),
        .studyOf,
        panel = panel,
        phenotypes = map(contexts, .contextOf, y = y, brain = brain)
    )
    donors <- rownames(G)
    MultiStudyQtlDataset(
        qtlDatasets = list(
            s1 = full$s1[, donors[1:82], ],
            s2 = full$s2[, donors[83:165], ]
        )
    )
}


# A study of `phenotypes` on the `panel` genotypes.
.studyOf <- function(name, panel, phenotypes) {
    QtlDataset(studyName = name, genotypes = panel, phenotypes = phenotypes)
}


# One context's phenotypes: column cx of `y` on the brain context's trait.
.contextOf <- function(cx, y, brain) {
    se <- brain[, rownames(y)]
    SummarizedExperiment::assay(se) <- matrix(
        y[, match(cx, c("a", "b"))],
        1L,
        dimnames = dimnames(se)
    )
    se
}


# A summary-statistics study "s3" of context a's trait, on study s1's donors
# with their LD as its reference, beside a MultiStudyQtlDataset's studies;
# its trait anchored `shift` bases from the studies'.
.summaryStudy <- function(ms, shift = 0L) {
    qd <- as.list(ms)$s1
    panel <- pecotmr:::.qtlGenotypeSe(qd)
    rr <- SummarizedExperiment::rowRanges(panel)
    a <- pecotmr:::.qtlContextSe(qd, "a")
    y <- SummarizedExperiment::assay(a)[1L, ]
    G <- t(as.matrix(SummarizedExperiment::assay(panel)))
    donors <- intersect(names(y)[!is.na(y)], rownames(G))
    X <- scale(G[donors, ])
    X[!is.finite(X)] <- 0
    r <- suppressWarnings(as.vector(stats::cor(X, y[donors])))
    r <- replace(r, !is.finite(r), 0)
    ids <- paste(
        GenomicRanges::seqnames(rr),
        GenomicRanges::start(rr),
        rr$A2,
        rr$A1,
        sep = ":"
    )
    gr <- GenomicRanges::GRanges(
        GenomicRanges::seqnames(rr),
        IRanges::IRanges(GenomicRanges::start(rr), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = ids,
        A1 = rr$A1,
        A2 = rr$A2,
        Z = pecotmr:::.corToZ(r, length(donors)),
        N = length(donors)
    )
    R <- stats::cov2cor(
        crossprod(X) + diag(as.numeric(colSums(X^2) == 0), ncol(X))
    )
    QtlSumStats(
        studyName = "s3",
        context = "a",
        trait = rownames(a),
        entry = list(gr),
        genome = "hg19",
        ldReference = LdMatrix(
            `dimnames<-`(R, list(ids, ids)),
            nSamples = length(donors)
        ),
        traitPos = IRanges::shift(SummarizedExperiment::rowRanges(a), shift),
        qcOptions = list(qc = "ok")
    )
}
