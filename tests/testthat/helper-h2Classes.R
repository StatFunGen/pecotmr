# Fixtures shared by the LD-statistic tests (LdEigen, LdScore, LdStatistic),
# the heritability tests (h2EstimationWrappers, H2Estimate) and
# test_genotypeHandle.R.

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

    new(
        "H2Estimate",
        estimate = 0.3,
        estimateSe = 0.05,
        intercept = 1.01,
        interceptSe = 0.02,
        localBlocks = NULL,
        enrichment = enrich,
        annotationJackknifeCoefs = annotationJackknifeCoefs,
        scoreStats = NULL,
        method = "lder",
        nSnps = 100L,
        traitName = "test"
    )
}
