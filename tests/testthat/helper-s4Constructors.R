# ===========================================================================
# Shared test helpers
# ===========================================================================

.sc_makeGenotypeHandle <- function(snp_n = 5L) {
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
        nSamples = 100L,
        sampleIds = paste0("s", seq_len(100)),
        pgenPtr = NULL
    )
}

# `offset` shifts the variant positions so two fixtures can occupy DISTINCT
# genomic ranges. Row identity is now (study, method, range), so two entries
# built at the same positions are the same block by definition -- tests that
# need two blocks have to say where each one sits.
.sc_makeTopLoci <- function(n = 3, offset = 0L) {
    data.frame(
        variant_id = paste0("chr1:", offset + 100 * seq_len(n), ":A:G"),
        chrom = rep("1", n),
        pos = as.integer(offset + 100 * seq_len(n)),
        A1 = rep("G", n),
        A2 = rep("A", n),
        N = rep(1000, n),
        af = rep(0.87, n), # directional effect-allele freq (> 0.5)
        marginal_beta = rep(0.1, n),
        marginal_se = rep(0.05, n),
        marginal_z = rep(2.0, n),
        marginal_p = rep(0.05, n),
        pip = seq(0.9, by = -0.1, length.out = n),
        posterior_mean = rep(0.05, n),
        posterior_sd = rep(0.02, n),
        cs_95 = paste0("susie_", c(1L, 1L, 0L)[seq_len(n)]),
        stringsAsFactors = FALSE
    )
}

.sc_makeFineMappingRow <- function(n = 3, offset = 0L) {
    FineMappingRow(
        variantIds = paste0("chr1:", offset + 100 * seq_len(n), ":A:G"),
        susieFit = list(fake = TRUE),
        topLoci = .sc_makeTopLoci(n, offset = offset)
    )
}

# ===========================================================================
# FineMappingRow
# ===========================================================================
.sc_makeSe <- function(traits = c("ENSG1", "ENSG2"), n_samples = 10) {
    rng <- GenomicRanges::GRanges(
        seqnames = rep("chr1", length(traits)),
        ranges = IRanges::IRanges(
            start = seq(1000L, by = 1000L, length.out = length(traits)),
            width = 500L
        )
    )
    names(rng) <- traits
    expr <- matrix(
        rnorm(length(traits) * n_samples),
        nrow = length(traits),
        ncol = n_samples,
        dimnames = list(traits, paste0("s", seq_len(n_samples)))
    )
    cd <- S4Vectors::DataFrame(
        sex = rep(c("M", "F"), length.out = n_samples),
        row.names = paste0("s", seq_len(n_samples))
    )
    SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = expr),
        rowRanges = rng,
        colData = cd
    )
}

.qtlMakeGenotypeHandle <- function() {
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

.qtlMakeEntryGr <- function(
    n = 5,
    chr = "chr1",
    start_at = 100L,
    step = 100L,
    with_maf = FALSE
) {
    gr <- GenomicRanges::GRanges(
        seqnames = rep(chr, n),
        ranges = IRanges::IRanges(
            start = seq(start_at, by = step, length.out = n),
            width = 1L
        )
    )
    mcols_list <- list(
        SNP = paste0("rs", seq_len(n)),
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = seq(1.0, by = 0.5, length.out = n),
        N = rep(1000L, n)
    )
    if (with_maf) {
        mcols_list$MAF <- seq(0.05, by = 0.01, length.out = n)
    }
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mcols_list)
    gr
}

.qtlMakeOne <- function(
    studyName = "study1",
    context = "Cortex",
    trait = "ENSG001",
    n = 5,
    varY = NA_real_,
    with_maf = FALSE,
    qcOptions = list(),
    genome = "hg19"
) {
    QtlSumStats(
        studyName = studyName,
        context = context,
        trait = trait,
        entry = list(.qtlMakeEntryGr(n, with_maf = with_maf)),
        genome = genome,
        ldReference = .qtlMakeGenotypeHandle(),
        varY = varY,
        qcOptions = qcOptions
    )
}

.cbp_makeHandle <- function(snp_n = 6L, n_samples = 30L, sample_prefix = "s") {
    new(
        "GenotypeHandle",
        path = "/tmp/cb",
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

.afm_df <- function(freqcol, val) {
    d <- data.frame(
        chrom = "chr1",
        pos = 100L,
        variant_id = "chr1:100:A:G",
        A1 = "A",
        A2 = "G",
        z = 2.5,
        n_sample = 1000,
        stringsAsFactors = FALSE
    )
    d[[freqcol]] <- val
    d
}

.afm_cm <- function(...) {
    c(
        chrom = "chrom",
        pos = "pos",
        variant_id = "variant_id",
        A1 = "A1",
        A2 = "A2",
        z = "z",
        n_sample = "n_sample",
        ...
    )
}

# A deterministic per-gene QtlSumStats: one row per gene carrying the regional
# summary as row columns (p_beta / beta shapes / n_variants), each entry a
# per-variant GRanges (P / af / tss_distance / tes_distance, plus the `qvalue`
# the qvalue SNP significance query reads, with v1 tracking the gene signal).
# Deterministic so qvalue::qvalue sees a stable p-value distribution.
#   emptyIdx   gene indices whose entry is an empty GRanges (no p-values), for
#              the "skip genes with no variants" branch; the p_beta
#              distribution is unchanged, so qvalue::qvalue still estimates
#              pi0 on its default path.
.qapFixture <- function(nVar = 50L, nVarFilt = 30L, emptyIdx = integer(0)) {
    # 6 strong signal p_beta + a deterministic uniform null (ppoints) so
    # qvalue::qvalue estimates pi0 on its default path (no fallback needed).
    # 18 null genes keep its smoother well clear of failing (it does by 10).
    pBeta <- c(10^(-c(8, 7, 6, 5, 4, 3)), stats::ppoints(18))
    G <- length(pBeta)
    entries <- map(seq_len(G), function(i) {
        if (i %in% emptyIdx) {
            return(GenomicRanges::GRanges())
        } # no variants
        k <- 4L
        pv <- c(pBeta[i] / 5, 0.2, 0.5, 0.8) # min tracks the gene signal
        gr <- GenomicRanges::GRanges(
            "chr1",
            IRanges::IRanges(
                start = seq(1000L, by = 50L, length.out = k),
                width = 1L
            )
        )
        mc <- S4Vectors::DataFrame(
            SNP = paste0("g", i, "_v", seq_len(k)),
            A1 = "A",
            A2 = "G",
            P = pv,
            af = c(0.30, 0.20, 0.005, 0.40), # v3 fails MAF 0.01
            tss_distance = c(0L, 500L, 900000L, 2000000L), # v4 outside 1Mb cis
            tes_distance = c(0L, 500L, 900000L, 2000000L),
            qvalue = c(pBeta[i] / 5, 0.30, 0.60, 0.90)
        )
        S4Vectors::mcols(gr) <- mc
        gr
    })
    QtlSumStats(
        study = rep("s", G),
        context = rep("brain", G),
        trait = paste0("g", seq_len(G)),
        entry = entries,
        genome = "hg19",
        extraCols = list(
            n_variants = rep(nVar, G),
            n_variants_filtered = rep(nVarFilt, G),
            p_beta = pBeta,
            beta_shape1 = rep(1, G),
            beta_shape2 = rep(200, G)
        )
    )
}

# The default fixture is 24 genes of GRanges + mcols; building it, like every
# postprocess / significantQtls call, costs time in proportion to the gene
# count, so the tests share one instance, built on first use. S4 objects are
# copy-on-modify, so that is safe; tests needing new entries still call
# .qapFixture() with arguments.
delayedAssign(".qapDefaultFixture", .qapFixture())
