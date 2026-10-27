context("ctwasPipeline")

# ===========================================================================
# Strategy: ctwas::ctwas_sumstats does the heavy work. We mock it to a
# function that just returns its inputs back, so we can verify how the
# pipeline assembles z_snp / weights / region_info / LD loader inputs.
# ===========================================================================

# 12 = 6 variants x the 2 blocks the default GwasSumStats fixture builds. The
# LD-sketch compatibility check compares panels between the GWAS and weight
# sides, so every fixture has to draw on the same one.
.ctp_makeHandle <- function(snp_n = 12L, n_samples = 30L) {
    # Use a per-process tempfile so .ctwasLdPanelKey's file.exists check
    # succeeds against the fixture handle (real LD-sketch payloads exist
    # by construction; mock fixtures need an equivalent on-disk anchor).
    stem <- file.path(tempdir(), "ctp_sketch")
    if (!file.exists(paste0(stem, ".pgen"))) {
        file.create(paste0(stem, ".pgen"))
    }
    positions <- seq(100L, by = 100L, length.out = snp_n)
    # SNP IDs follow the canonical chr:pos:A2:A1 layout so allele
    # harmonization inside .ctwasBuildWeights / .ctwasHarmonizeWeights can
    # parse them via parseVariantId().
    snpIds <- sprintf("chr1:%d:G:A", positions)
    new(
        "GenotypeHandle",
        path = stem,
        format = "plink2",
        snpInfo = data.frame(
            SNP = snpIds,
            CHR = rep("chr1", snp_n),
            BP = positions,
            A1 = rep("A", snp_n),
            A2 = rep("G", snp_n),
            stringsAsFactors = FALSE
        ),
        nSamples = n_samples,
        sampleIds = paste0("s", seq_len(n_samples)),
        pgenPtr = NULL
    )
}

# A panel with its correlation already built, the shape callers saw before
# the LD became on-demand. Tests that compare against that matrix want it
# materialized.
.ctp_fullPanel <- function(h) {
    p <- pecotmr:::.ctwasPanelFor(h)
    c(p, list(R = pecotmr:::.ctwasPanelLd(p)))
}

# A genotype panel's variants snpInfo-style (SNP, CHR, BP, A1, A2): the
# example dataset's PLINK1 panel reads through no handle to ask instead.
.ctp_panelSnpInfo <- function(panel) {
    pecotmr:::.rangesSnpInfo(SummarizedExperiment::rowRanges(panel))
}

# Canonical SNP IDs the fixtures use. .ctp_makeHandle() emits these in
# chr:pos:A2:A1 form; tests reference them by index via .ctp_snpId(i).
.ctp_snpId <- function(i) sprintf("chr1:%d:G:A", 100L * i)

.ctp_mockExtractor <- function(seed = 5, n_samples = 30L) {
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

# Variants 1..6 for block 1, 7..12 for block 2, and so on. LD blocks partition
# the genome, and a GwasSumStats requires (study, range) to be unique, so two
# blocks cannot carry the same variants.
.ctp_blockVariants <- function(b) {
    idx <- seq_len(6L) + 6L * (b - 1L)
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = 100L * idx, width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = map_chr(idx, .ctp_snpId),
        A1 = rep("A", 6),
        A2 = rep("G", 6),
        Z = rnorm(6),
        N = rep(1000L, 6)
    )
    gr
}

# One collection, one element per LD block (§4.7), replacing the named list of
# single-block GwasSumStats keyed by region_id.
.ctp_makeGwasSumstats <- function(qc = TRUE, blockIds = c("block1", "block2")) {
    n <- length(blockIds)
    GwasSumStats(
        studyName = rep("G1", n),
        entry = map(seq_len(n), .ctp_blockVariants),
        genome = "hg19",
        ldReference = .ctp_makeHandle(snp_n = 6L * n),
        blockId = blockIds,
        qcOptions = if (qc) list(step1 = "ok") else list()
    )
}

# `variantIdx` picks which of the fixture's 12 variants carry the gene's
# weights; it has to name five of them, one per weight. The default sits
# inside block 1 -- pass 7:11 for a gene that lives in block 2.
.ctp_makeTwasWeights <- function(variantIdx = 1:5) {
    e <- TwasWeightsRow(
        variantIds = map_chr(variantIdx, .ctp_snpId),
        weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
    )
    TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e),
        ldReference = .ctp_makeHandle()
    )
}

# The same per-gene weights expressed as a QtlFineMappingResult weight source:
# the topLoci posterior_mean carries the weight vector (what resolveWeights
# reads). Mirrors .ctp_makeTwasWeights so the two sources are comparable.
.ctp_makeFmrWeightSource <- function() {
    vids <- map_chr(1:5, .ctp_snpId)
    w <- c(0.1, 0.05, -0.2, 0.3, 0.0)
    tl <- data.frame(
        variant_id = vids,
        chrom = rep("1", 5),
        pos = 100L * (1:5),
        A1 = rep("A", 5),
        A2 = rep("G", 5),
        N = rep(100, 5),
        af = rep(0.3, 5),
        pip = rep(0.5, 5),
        posterior_mean = w,
        posterior_sd = rep(0.1, 5),
        cs_95 = rep("susie_1", 5),
        stringsAsFactors = FALSE
    )
    fe <- FineMappingRow(variantIds = vids, susieFit = list(), topLoci = tl)
    QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fe),
        ldReference = .ctp_makeHandle()
    )
}

# ===========================================================================
# Input validation
# ----------------------------------------------------------------------------
# The top of ctwasPipeline requires the (non-CRAN) `ctwas` package; without
# it the function errors out before any input-validation branch fires. Skip
# entry-point tests when ctwas isn't installed, but exercise the input-
# building helpers directly (they don't gate on ctwas).
# ===========================================================================

# Helper: minimal two-block input set for the multi-block API tests.
.ctp_makeMultiBlockInputs <- function(qc = TRUE) {
    ss <- .ctp_makeGwasSumstats(qc = qc)
    list(
        gwasSumStats = ss,
        # ctwas fine-maps one region at a time, so each block's gene has to
        # carry weights on variants that block actually has LD for: block 1
        # holds variants 1..6, block 2 holds 7..12.
        twasWeights = list(
            block1 = .ctp_makeTwasWeights(),
            block2 = .ctp_makeTwasWeights(7:11)
        )
    )
}

test_that("ctwasPipeline: rejects a per-region list of GwasSumStats", {
    skip_if_not_installed("ctwas")
    # The inverse of the old contract (§4.7): one collection whose elements are
    # blocks, not a named list of single-block collections.
    ss <- .ctp_makeGwasSumstats()
    expect_error(
        ctwasPipeline(
            gwasSumStats = list(block1 = ss, block2 = ss),
            twasWeights = .ctp_makeTwasWeights(),
            method = "susie"
        ),
        "must be a GwasSumStats whose elements are LD blocks"
    )
})

test_that("ctwasPipeline: rejects a single-block named list", {
    skip_if_not_installed("ctwas")
    expect_error(
        ctwasPipeline(
            gwasSumStats = .ctp_makeGwasSumstats(blockIds = "block1"),
            twasWeights = list(block1 = .ctp_makeTwasWeights()),
            method = "susie"
        ),
        "at least two LD blocks"
    )
})

test_that("ctwasPipeline: rejects un-QCd GwasSumStats in any region", {
    skip_if_not_installed("ctwas")
    ss_qc <- .ctp_makeGwasSumstats(qc = TRUE)
    ss_noqc <- .ctp_makeGwasSumstats(qc = FALSE)
    tw <- .ctp_makeTwasWeights()
    expect_error(
        ctwasPipeline(
            gwasSumStats = ss_noqc,
            twasWeights = list(block1 = tw, block2 = tw),
            method = "susie"
        ),
        "has no QC record"
    )
})

test_that("ctwasPipeline: rejects twasWeights keys not present in gwasSumStats", {
    skip_if_not_installed("ctwas")
    ss <- .ctp_makeGwasSumstats()
    tw <- .ctp_makeTwasWeights()
    expect_error(
        ctwasPipeline(
            gwasSumStats = .ctp_makeGwasSumstats(
                blockIds = c("blockA", "blockB")
            ),
            twasWeights = list(blockA = tw, blockC = tw),
            method = "susie"
        ),
        "key.*not present in.*gwasSumStats"
    )
})

test_that("assembleCtwasInputs: allows twasWeights keys to be a subset of gwasSumStats", {
    ss <- .ctp_makeGwasSumstats()
    tw <- .ctp_makeTwasWeights()
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    # Two blocks supply zSnp; only block1 supplies TwasWeights, whose panel
    # shares block1's six variants of its twelve -- a partial overlap,
    # reported once per session.
    rlang::reset_warning_verbosity("pecotmrLdReferenceOverlap-ctwasPipeline")
    expect_warning(
        inputs <- assembleCtwasInputs(
            gwasSumStats = ss,
            twasWeights = list(block1 = tw),
            method = "susie"
        ),
        "share 6 variant"
    )
    # Both blocks appear in region_info / snp_map (SNP-only block2 contributes
    # its zSnp), but the weights list only has block1-keyed entries.
    expect_setequal(inputs$region_info$region_id, c("block1", "block2"))
    expect_setequal(names(inputs$snp_map), c("block1", "block2"))
    expect_true(all(grepl("^block1\\|", names(inputs$weights))))
})

test_that("assembleCtwasInputs: accepts a QtlFineMappingResult weight source (topLoci posterior effect)", {
    ss <- .ctp_makeGwasSumstats()
    # The FMR topLoci posterior effect is on the STANDARDIZED scale, so it matches
    # a standardized TwasWeights carrying the same effect vector (both skip cTWAS's
    # variance scaling). Compare against that, not the default (unstandardized) one.
    tw_std <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(TwasWeightsRow(
            variantIds = map_chr(1:5, .ctp_snpId),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0),
            weightStandardized = TRUE
        )),
        ldReference = .ctp_makeHandle()
    )
    fmr <- .ctp_makeFmrWeightSource()
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    itw <- assembleCtwasInputs(
        gwasSumStats = ss,
        twasWeights = list(block1 = tw_std, block2 = tw_std),
        method = "susie"
    )
    ifmr <- assembleCtwasInputs(
        gwasSumStats = ss,
        twasWeights = list(block1 = fmr, block2 = fmr),
        method = "susie"
    )
    # Same gene keys, and the resolved weights match the standardized TwasWeights.
    expect_equal(names(ifmr$weights), names(itw$weights))
    expect_true(length(ifmr$weights) > 0L)
    expect_equal(ifmr$weights[[1L]], itw$weights[[1L]])
})

test_that("assembleCtwasInputs: boundary gene fits per-region, spans all", {
    # Build two blocks with NON-OVERLAPPING GWAS variants. Block 1 covers
    # v1..v3, block 2 covers v4..v6. The gene's weight spans v2..v5 — i.e.
    # crosses the block boundary. With a per-block filter the gene would
    # lose v4..v5 (block-2 variants); with a global-union filter all four
    # weight variants survive.
    mkBlockGss <- function(studyName, snpIds, qc = TRUE) {
        gr <- GenomicRanges::GRanges(
            seqnames = "chr1",
            ranges = IRanges::IRanges(
                start = as.integer(gsub(".*:([0-9]+):.*", "\\1", snpIds)),
                width = 1L
            )
        )
        S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
            SNP = snpIds,
            A1 = rep("A", length(snpIds)),
            A2 = rep("G", length(snpIds)),
            Z = rnorm(length(snpIds)),
            N = rep(1000L, length(snpIds))
        )
        GwasSumStats(
            studyName = studyName,
            entry = list(gr),
            genome = "hg19",
            ldReference = .ctp_makeHandle(),
            qcOptions = if (qc) list(step1 = "ok") else list()
        )
    }
    ss1 <- mkBlockGss("G1", map_chr(1:3, .ctp_snpId))
    ss2 <- mkBlockGss("G2", map_chr(4:6, .ctp_snpId))
    # Cross-boundary weights: v2..v5 (4 variants spanning both blocks).
    crossEntry <- TwasWeightsRow(
        variantIds = map_chr(2:5, .ctp_snpId),
        weights = c(0.1, 0.2, 0.3, 0.4)
    )
    tw <- TwasWeights(
        studyName = "G1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(crossEntry),
        ldReference = .ctp_makeHandle()
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    both <- c(ss1, ss2)
    mcols(both)$blockId <- c("block1", "block2")
    inputs <- assembleCtwasInputs(
        gwasSumStats = both,
        twasWeights = list(block1 = tw),
        method = "susie"
    )
    # ctwas fine-maps one region at a time, so the FITTED vector holds only the
    # variants this block has LD for -- handing susie_rss the other two produces
    # a non-finite ELBO. The gene's SPAN still covers all four, which is what
    # ctwas::get_boundary_genes reads to route the gene to merge_regions; a span
    # clipped to the block would hide the boundary gene entirely.
    entry <- inputs$weights[[1L]]
    wgt <- entry$wgt
    expect_equal(nrow(wgt), 2L)
    expect_setequal(rownames(wgt), map_chr(2:3, .ctp_snpId))
    expect_length(entry$wgt_var, 1L)
    expect_gt(entry$wgt_var, 0)
    expect_equal(entry$n_wgt, 2L)
    expect_equal(entry$p0, 200L)
    expect_equal(entry$p1, 500L)
})

test_that("ctwasPipeline: a bare TwasWeights is a FLAT source, placed by region", {
    skip_if_not_installed("ctwas")
    # A bare TwasWeights is now accepted as a flat weight source and placed into
    # LD blocks by region; the fixture carries no region, so placement errors
    # (rather than the old 'must be a NAMED LIST' shape rejection).
    expect_error(
        ctwasPipeline(
            gwasSumStats = .ctp_makeGwasSumstats(),
            twasWeights = .ctp_makeTwasWeights(),
            method = "susie"
        ),
        "no `traitPos` provenance"
    )
})

test_that("ctwasPipeline: rejects non-GRanges twasZ", {
    skip_if_not_installed("ctwas")
    inp <- .ctp_makeMultiBlockInputs()
    expect_error(
        ctwasPipeline(
            gwasSumStats = inp$gwasSumStats,
            twasWeights = inp$twasWeights,
            twasZ = "not a GRanges",
            method = "susie"
        ),
        "must be a GRanges"
    )
})

test_that("CtwasPriorParam rejects an unknown varStructure", {
    # The check moved to the constructor, so it fires where it is written
    # rather than part-way into a run -- and needs no ctwas install.
    expect_error(
        CtwasPriorParam(varStructure = "bogus"),
        "`varStructure` must be one of"
    )
    expect_equal(CtwasPriorParam()$varStructure, "shared_type")
    expect_equal(
        CtwasPriorParam(varStructure = "independent")$varStructure,
        "independent"
    )
})

# ===========================================================================
# .ctwasRequireMatchingLdRefs
# ===========================================================================

test_that(".ctwasRequireMatchingLdRefs: NULL twas-side handle is allowed", {
    twNoLd <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(TwasWeightsRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
            weights = rep(0.1, 5)
        )),
        ldReference = NULL
    )
    expect_silent(pecotmr:::.ctwasRequireMatchingLdRefs(
        twNoLd,
        .ctp_makeGwasSumstats()
    ))
})

# One weight row on `handle`, for the region-level match checks.
.ctp_weightsOn <- function(handle) {
    TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(TwasWeightsRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
            weights = rep(0.1, 3)
        )),
        ldReference = handle
    )
}

test_that(".ctwasRequireMatchingLdRefs: differently trimmed panels ok", {
    # One LD reference QC'd separately on the two sides leaves each sketch
    # trimmed to its own surviving variants, so the panels overlap without
    # being identical. That is the normal case, not an error.
    #
    # The overlap note is rate-limited to once per session and the assembly
    # tests legitimately trip it first (each block's GWAS panel is narrowed to
    # its own variants while the weights carry the whole one), so clear that
    # budget here rather than depending on which files ran before this.
    rlang::reset_warning_verbosity("pecotmrLdReferenceOverlap-ctwasPipeline")
    tw <- .ctp_weightsOn(.ctp_makeHandle(snp_n = 5L))
    gwas <- pecotmr:::.withRowLdRefs(
        .ctp_makeGwasSumstats(),
        .ctp_makeHandle(snp_n = 6L)
    )
    expect_warning(
        expect_null(pecotmr:::.ctwasRequireMatchingLdRefs(tw, gwas)),
        "share 5 variant"
    )
})

test_that(".ctwasRequireMatchingLdRefs: disjoint panels error", {
    gwasLd <- .ctp_makeHandle(snp_n = 3L)
    si <- snpInfo(gwasLd)
    si$BP <- si$BP + 1e6L
    gwasLd@snpInfo <- si
    expect_error(
        pecotmr:::.ctwasRequireMatchingLdRefs(
            .ctp_weightsOn(.ctp_makeHandle(snp_n = 3L)),
            pecotmr:::.withRowLdRefs(.ctp_makeGwasSumstats(), gwasLd)
        ),
        "share no variant"
    )
})

# ===========================================================================
# Input-building helpers
# ===========================================================================

test_that(".ctwasBuildZSnp: produces a flat data.frame keyed by SNP/study", {
    ss <- .ctp_makeGwasSumstats(blockIds = "block1")
    df <- pecotmr:::.ctwasBuildZSnp(ss, map_chr(1:6, .ctp_snpId))
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 6L)
    expect_setequal(
        colnames(df),
        c("id", "chrom", "pos", "A1", "A2", "z", "study")
    )
    expect_setequal(df$id, map_chr(1:6, .ctp_snpId))
    expect_setequal(unique(df$study), "G1")
})

test_that(".ctwasBuildZSnp: negates z for a panel-flipped variant", {
    # z arrives in the GWAS frame while R/variance/weights are in the panel's.
    # An exact-string join dropped a swapped spelling silently -- present on
    # both sides, so nothing reported it missing and the gene just lost an
    # instrument.
    ss <- .ctp_makeGwasSumstats(blockIds = "block1")
    gwasIds <- map_chr(1:6, .ctp_snpId)
    panelIds <- gwasIds
    panelIds[2] <- "chr1:200:A:G" # fixture spells this one G:A
    plain <- pecotmr:::.ctwasBuildZSnp(ss, gwasIds)
    flipped <- pecotmr:::.ctwasBuildZSnp(ss, panelIds)
    i <- match("chr1:200:G:A", plain$id)
    expect_equal(flipped$id[i], "chr1:200:A:G")
    expect_equal(flipped$z[i], -plain$z[i])
    expect_equal(flipped$A1[i], plain$A2[i])
    expect_equal(flipped$A2[i], plain$A1[i])
    expect_equal(flipped$id[-i], plain$id[-i])
    expect_equal(flipped$z[-i], plain$z[-i])
})

test_that(".ctwasBuildZSnp: leaves variants absent from the panel untouched", {
    ss <- .ctp_makeGwasSumstats(blockIds = "block1")
    gwasIds <- map_chr(1:6, .ctp_snpId)
    plain <- pecotmr:::.ctwasBuildZSnp(ss, gwasIds)
    none <- pecotmr:::.ctwasBuildZSnp(ss, "chr9:999:A:T")
    expect_equal(none$id, plain$id)
    expect_equal(none$z, plain$z)
})

test_that(".ctwasFilterVariants: a CS member the fine-mapping spells flipped still survives the cap", {
    # The weights arrive in the panel frame, the fine-mapping auxiliaries in
    # the fine-mapping result's. An exact-string join lost the credible-set
    # member's must-keep protection, so the cap dropped the variant the
    # fine-mapping was most confident about and kept the largest |w| instead.
    vids <- c("chr1:100:G:A", "chr1:200:G:A", "chr1:300:G:A")
    w <- c(0.1, 0.2, 0.9)
    aux <- list(
        pip = NULL,
        csMembers = list("chr1:100:A:G"), # flipped spelling of vids[1]
        csPurity = 0.9
    )
    out <- pecotmr:::.ctwasFilterVariants(
        vids,
        w,
        finemapAux = aux,
        twasWeightCutoff = 0,
        csMinCor = 0.8,
        minPipCutoff = 0,
        maxNumVariants = 1
    )
    expect_equal(out$vids, "chr1:100:G:A")
})

test_that(".ctwasFilterVariants: PIPs spelled flipped still drive the cap", {
    vids <- c("chr1:100:G:A", "chr1:200:G:A", "chr1:300:G:A")
    w <- c(0.1, 0.2, 0.9)
    aux <- list(
        pip = c("chr1:100:A:G" = 0.95, "chr1:300:A:G" = 0.01),
        csMembers = list(),
        csPurity = numeric(0)
    )
    out <- pecotmr:::.ctwasFilterVariants(
        vids,
        w,
        finemapAux = aux,
        twasWeightCutoff = 0,
        csMinCor = 0,
        minPipCutoff = 0,
        maxNumVariants = 1
    )
    # PIP 0.95 wins over the larger |w| whose PIP is 0.01.
    expect_equal(out$vids, "chr1:100:G:A")
})

test_that(".ctwasSnpInfoForGwasBlock: keeps rows the GWAS spells flipped", {
    # Regression: `is_in()` on the id string retained none of these.
    ss <- .ctp_makeGwasSumstats(blockIds = "block1")
    flippedPanel <- data.frame(
        chrom = rep(1L, 6),
        id = sprintf("chr1:%d:A:G", 100L * (1:6)),
        pos = 100L * (1:6),
        alt = rep("A", 6),
        ref = rep("G", 6),
        stringsAsFactors = FALSE
    )
    keep <- pecotmr:::.ctwasSnpInfoForGwasBlock(ss, flippedPanel)
    expect_equal(nrow(keep), 6L)
    expect_equal(keep$id, flippedPanel$id)
})

test_that(".ctwasBuildSingleRegionInfo: a given region bounds the block", {
    # The LD block a GWAS was split by is its window, however much narrower
    # its variants are.
    gss <- GwasSumStats(
        studyName = "G1",
        entry = list(.ctp_blockVariants(1L)),
        genome = "hg19",
        region = GenomicRanges::GRanges("chr1:1-1000")
    )
    ri <- pecotmr:::.ctwasBuildSingleRegionInfo("chr1_1_1000", gss)
    expect_equal(ri$start, 1L)
    expect_equal(ri$stop, 1001L)
})

test_that(".ctwasBuildSingleRegionInfo: uses the block entry span, not the wider shared LD sketch", {
    # Regression: many LD blocks can share one whole-chromosome LD payload, so
    # the sketch span (here BP 100-1200) is NOT the block's span. The entry
    # here covers only 200-400; region bounds must follow the entry, otherwise
    # every block collapses to the whole-chromosome span and every SNP is
    # assigned to every region (inflating SNP group_size and crushing the gene
    # PIP).
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = c(200L, 300L, 400L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = c("a", "b", "c"),
        A1 = "A",
        A2 = "G",
        Z = 0,
        N = 1000L
    )
    gss <- GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ctp_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
    ri <- pecotmr:::.ctwasBuildSingleRegionInfo("blockX", gss)
    expect_equal(ri$region_id, "blockX")
    expect_equal(ri$chrom, 1L)
    # cTWAS reads [start, stop), so stop is one past the window's last base.
    expect_equal(ri$start, 200L) # entry min, not sketch min (100)
    expect_equal(ri$stop, 401L) # one past entry max, not sketch max (1200)
})

test_that(".ctwasBuildSingleRegionInfo: multi-chromosome block entry errors", {
    gr <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = c("a", "b", "c"),
        A1 = "A",
        A2 = "G",
        Z = 0,
        N = 1000L
    )
    gss <- GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ctp_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
    expect_error(
        pecotmr:::.ctwasBuildSingleRegionInfo("block1", gss),
        "spans multiple chromosomes"
    )
})

test_that(".ctwasSnpInfoForBlock: returns ctwas-required columns chrom/id/pos/alt/ref", {
    df <- pecotmr:::.ctwasSnpInfoForBlock(.ctp_makeHandle())
    # ctwas's read_snp_info_files asserts these exact column names
    expect_setequal(colnames(df), c("chrom", "id", "pos", "alt", "ref"))
})

test_that(".ctwasSnpInfoForBlock: repairs a tag-allele panel id from A1/A2", {
    # The id becomes the R dimnames and `panelSnps`, both compared against
    # harmonized (reference-allele) ids, so a tag id would drop the variant.
    h <- .ctp_makeHandle(snp_n = 3L)
    h@snpInfo$SNP[2L] <- "chr1:200:INS:A"
    df <- pecotmr:::.ctwasSnpInfoForBlock(h)
    expect_equal(df$id[2L], "chr1:200:G:A")
    expect_equal(df$ref[2L], "G")
    expect_equal(df$alt[2L], "A")
    # the LD matrix and variance vector are keyed the same way
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    panel <- .ctp_fullPanel(h)
    expect_true("chr1:200:G:A" %in% rownames(panel$R))
    expect_true("chr1:200:G:A" %in% names(panel$variance))
})

test_that(".ctwasLdPanelKey: returns the file behind an existing sketch", {
    handle <- .ctp_makeHandle()
    key <- pecotmr:::.ctwasLdPanelKey(handle)
    expect_true(file.exists(key))
    expect_equal(key, paste0(path(handle), ".pgen"))
})

test_that(".ctwasLdPanelKey: a VCF panel's key is its file", {
    skip_if_not_installed("VariantAnnotation")
    vcf <- .writeTestVcf(.vcfGtOnly(.vcfOverlapRecords()), c("s1", "s2", "s3"))
    panel <- suppressWarnings(readGenotypes(vcf))
    expect_equal(pecotmr:::.ctwasLdPanelKey(panel), normalizePath(vcf))
})

test_that(".ctwasLdPanelKey: errors when no candidate file exists", {
    ghost <- new(
        "GenotypeHandle",
        path = file.path(tempdir(), "does_not_exist_for_test.pgen_stem"),
        format = "plink2",
        snpInfo = data.frame(
            SNP = "chr1:100:A:G",
            CHR = "chr1",
            BP = 100L,
            A1 = "A",
            A2 = "G",
            stringsAsFactors = FALSE
        ),
        nSamples = 1L,
        sampleIds = "s1",
        pgenPtr = NULL
    )
    expect_error(
        pecotmr:::.ctwasLdPanelKey(ghost),
        "could not find the file behind the LD panel"
    )
})

test_that(".ctwasLdPanelKey: a mixture checks each panel; memory needs none", {
    panel <- readGenotypes(test_path("test_data", "test_variants.bed"))
    mix <- GenotypeMixture(
        list(a = panel[, 1:40], b = panel[, 41:100]),
        mixtureWeights = c(0.5, 0.5)
    )
    expect_null(pecotmr:::.ctwasLdPanelKey(mix))
    # One panel whose file has gone fails the whole mixture.
    dir <- withr::local_tempdir()
    file.copy(
        test_path(
            "test_data",
            str_c("test_variants", c(".bed", ".bim", ".fam"))
        ),
        dir
    )
    copy <- readGenotypes(file.path(dir, "test_variants.bed"))
    gone <- GenotypeMixture(
        list(a = copy[, 1:40], b = panel[, 41:100]),
        mixtureWeights = c(0.5, 0.5)
    )
    file.remove(file.path(dir, "test_variants.bed"))
    expect_error(
        pecotmr:::.ctwasLdPanelKey(gone),
        "could not find the file behind the LD panel"
    )
    # A panel held in memory reads no file, so there is none to find.
    inMemory <- .ctp_mockExtractor()(.ctp_makeHandle(snp_n = 6L), 1:6)
    expect_null(pecotmr:::.ctwasLdPanelKey(inMemory))
})

# ---------------------------------------------------------------------------
# Per-block LD scoping. Two bugs met here: combineGwasSumStats() UNIONS the
# per-piece sketches, and `[` carries that collection-level slot into every
# block, so each one referenced the whole chromosome; and the panel cache was
# keyed on the backing genotype file, which every block on a chromosome
# shares, so only block 1's panel was ever computed and the rest read it back.
# ---------------------------------------------------------------------------

test_that(".ctwasPickBlock: narrows the collection sketch to the block", {
    ss <- .ctp_makeGwasSumstats()
    # The collection's sketch spans both blocks, as combineGwasSumStats
    # leaves it.
    expect_equal(length(pecotmr:::.ldRefRanges(ldReference(ss)[[1L]])), 12L)
    byBlock <- pecotmr:::.ctwasGwasByBlock(ss)
    sizes <- map_int(byBlock, function(b) {
        length(pecotmr:::.ldRefRanges(ldReference(b)[[1L]]))
    })
    expect_equal(unname(sizes), c(6L, 6L))
    # Narrowed to DIFFERENT variants: block 1 holds 1..6, block 2 holds 7..12.
    pos <- map(byBlock, function(b) {
        GenomicRanges::start(pecotmr:::.ldRefRanges(ldReference(b)[[1L]]))
    })
    expect_equal(pos[[1L]], 100L * (1:6))
    expect_equal(pos[[2L]], 100L * (7:12))
})

test_that(".ctwasRegionLdTokens: one distinct, existing token per region", {
    ss <- .ctp_makeGwasSumstats()
    byBlock <- pecotmr:::.ctwasGwasByBlock(ss)
    sketches <- map(byBlock, function(b) ldReference(b)[[1L]])
    tokens <- pecotmr:::.ctwasRegionLdTokens(names(byBlock), sketches)
    expect_length(unique(unname(tokens)), 2L)
    expect_true(all(file.exists(tokens)))
    expect_equal(names(tokens), c("block1", "block2"))
})

test_that(".ctwasTokenSlug: sanitizes ids and keeps duplicates apart", {
    slug <- pecotmr:::.ctwasTokenSlug(
        c("chr19:1-500000", "chr19:500001-1000000", "chr19:1-500000")
    )
    expect_false(any(str_detect(slug, ":")))
    # The index prefix separates two ids that sanitize to the same stem.
    expect_length(unique(slug), 3L)
})

test_that("assembleCtwasInputs: each region gets its own block-sized panel", {
    inp <- .ctp_makeMultiBlockInputs()
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- assembleCtwasInputs(
        gwasSumStats = inp$gwasSumStats,
        twasWeights = inp$twasWeights,
        method = "susie"
    )
    expect_length(unique(out$LD_map$LD_file), 2L)
    panels <- pecotmr:::.ctwasCachedPanels(out)
    expect_length(panels, 2L)
    # Block-sized (6x6), not the 12-variant whole-panel LD.
    expect_equal(dim(pecotmr:::.ctwasPanelLd(panels[[1L]])), c(6L, 6L))
    expect_equal(dim(pecotmr:::.ctwasPanelLd(panels[[2L]])), c(6L, 6L))
    # And genuinely different LD, not block 1's twice.
    expect_false(identical(panels[[1L]]$snpInfo$id, panels[[2L]]$snpInfo$id))
})

test_that(".ctwasResolveLdPaths: re-mints sentinels lost with the session", {
    inp <- .ctp_makeMultiBlockInputs()
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- assembleCtwasInputs(
        gwasSumStats = inp$gwasSumStats,
        twasWeights = inp$twasWeights,
        method = "susie"
    )
    before <- pecotmr:::.ctwasCachedPanels(out)
    # The sentinels live in the session tempdir; a payload reloaded in a new
    # session finds them gone, and ctwas asserts file.exists() on each.
    unlink(out$LD_map$LD_file)
    restored <- pecotmr:::.ctwasResolveLdPaths(out)
    tokens <- restored$LD_map$LD_file
    expect_true(all(file.exists(tokens)))
    expect_length(unique(tokens), 2L)
    # Each region still reaches its OWN panel through the re-keyed loader.
    loader <- restored$LD_loader_fun
    expect_equal(loader(tokens[[1L]]), pecotmr:::.ctwasPanelLd(before[[1L]]))
    expect_equal(loader(tokens[[2L]]), pecotmr:::.ctwasPanelLd(before[[2L]]))
})

test_that(".ctwasResolveMethod: no-ensemble TwasWeights need `method`", {
    e <- TwasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
        weights = c(0.1, 0.2, 0.3)
    )
    tw <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "mrash",
        entry = list(e),
        ldReference = .ctp_makeHandle()
    )
    # One method present is still not a default: only "ensemble" is.
    expect_error(
        pecotmr:::.ctwasResolveMethod(list(r1 = tw)),
        str_c(
            "pass `method` to choose the weights\\. These TwasWeights carry ",
            "no 'ensemble' weights to default to \\(available: mrash\\)"
        )
    )
})

# ---------------------------------------------------------------------------
# .ctwasResolveMethod on the shapes a pipeline hands it
# ---------------------------------------------------------------------------
.ctp_multiMethodTw <- function(methods = c("mrash", "susie")) {
    TwasWeights(
        studyName = rep("Q1", length(methods)),
        context = rep("c1", length(methods)),
        trait = rep("t1", length(methods)),
        method = methods,
        entry = map(methods, function(m) {
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.1, 0.2, 0.3)
            )
        }),
        ldReference = .ctp_makeHandle()
    )
}

test_that(".ctwasResolveMethod reads one flat source or a list of them", {
    # A flat TwasWeights, and a pre-bucketed list with an empty block.
    expect_equal(
        pecotmr:::.ctwasResolveMethod(.ctp_multiMethodTw(c(
            "mrash",
            "ensemble"
        ))),
        "ensemble"
    )
    expect_equal(
        pecotmr:::.ctwasResolveMethod(
            list(r1 = .ctp_multiMethodTw(c("mrash", "susie")), r2 = NULL),
            method = "susie"
        ),
        "susie"
    )
    expect_error(
        pecotmr:::.ctwasResolveMethod(
            list(r1 = .ctp_multiMethodTw("mrash")),
            method = "lasso"
        ),
        "method 'lasso' not present in the weight sources"
    )
})

test_that(".ctwasResolveMethod: a fine-mapping weight source has no default", {
    # A fine-mapping result's weights come from one of its methods, and none
    # of them is an "ensemble" -- even one named that way would be a guess.
    src <- .ctp_makeFmrWeightSource()
    expect_error(
        pecotmr:::.ctwasResolveMethod(src),
        "one of its fine-mapping methods \\(available: susie\\)"
    )
    expect_equal(pecotmr:::.ctwasResolveMethod(src, "susie"), "susie")
})

# ---------------------------------------------------------------------------
# .ctwasGwasStudy — single-disease invariant
# ---------------------------------------------------------------------------
test_that(".ctwasGwasStudy: unique study, errors on mixed studies", {
    # It reads the collection's own `$study` column, one value per block, so a
    # data.frame stand-in with one row per block suffices (and dodges the S4
    # `$<-` coercion barrier on GwasSumStats).
    mk <- function(...) data.frame(study = c(...))
    expect_equal(pecotmr:::.ctwasGwasStudy(mk("D1", "D1")), "D1")
    expect_error(
        pecotmr:::.ctwasGwasStudy(mk("D1", "D2")),
        "multiple GWAS studies"
    )
    expect_true(is.na(pecotmr:::.ctwasGwasStudy(mk(NA_character_))))
})

# ---------------------------------------------------------------------------
# Phase 5: internal LD-block placement (by start(region) == cTWAS p0 rule)
# ---------------------------------------------------------------------------
# Three LD blocks split from one GWAS by their windows, as cTWAS takes them:
# one single-row slice per block, named by its blockId. The variants sit well
# inside each window (chr1:500, chr1:1500, chr22:800).
.ctp_windowGrid <- function() {
    gr <- GenomicRanges::GRanges(
        c("chr1", "chr1", "chr22"),
        IRanges::IRanges(start = c(500L, 1500L, 800L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = c("chr1:500:G:A", "chr1:1500:G:A", "chr22:800:G:A"),
        A1 = "A",
        A2 = "G",
        Z = 0,
        N = 1000L
    )
    blocks <- GenomicRanges::GRanges(
        c("chr1:1-1000", "chr1:1001-2000", "chr22:1-1000")
    )
    names(blocks) <- c("b1", "b2", "b3")
    pecotmr:::.ctwasGwasByBlock(GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        region = blocks
    ))
}

test_that(".ctwasPlaceByAnchor homes genes by start(region)", {
    region <- GenomicRanges::GRanges(
        c("chr1", "chr1", "chr22", "chr2"),
        IRanges::IRanges(
            start = c(500, 1500, 800, 100),
            end = c(600, 1600, 900, 200)
        )
    )
    home <- pecotmr:::.ctwasPlaceByAnchor(region, .ctp_windowGrid())
    # chr22:800 -> the chr22 block; chr2 has no block -> NA.
    expect_equal(home, c("b1", "b2", "b3", NA))
})

test_that(".ctwasPlaceByAnchor: a block holds its whole window", {
    # chr1:900 lies between block 1's variant (500) and block 2's (1500): the
    # window, not the variants, decides. The blocks' names are only names, and
    # the last base of a window (1000) is its own, the next (1001) the next
    # block's.
    region <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = c(900, 1000, 1001), width = 1)
    )
    expect_equal(
        pecotmr:::.ctwasPlaceByAnchor(region, .ctp_windowGrid()),
        c("b1", "b1", "b2")
    )
})

test_that(".ctwasIsPreBucketed / .ctwasCombineWeightSources: dispatch flat vs pre-bucketed", {
    tw <- .ctp_makeTwasWeights()
    grid <- list(block1 = 1, block2 = 2)
    expect_true(pecotmr:::.ctwasIsPreBucketed(list(block1 = tw))) # named list
    expect_false(pecotmr:::.ctwasIsPreBucketed(tw)) # flat S4
    expect_false(pecotmr:::.ctwasIsPreBucketed(list(tw))) # unnamed
    expect_s4_class(pecotmr:::.ctwasCombineWeightSources(tw), "TwasWeights")
    expect_s4_class(
        pecotmr:::.ctwasCombineWeightSources(list(tw)),
        "TwasWeights"
    )
})

test_that(".ctwasBucketWeights: places a flat 2-gene source into its home blocks", {
    mkE <- function() {
        TwasWeightsRow(
            variantIds = map_chr(1:5, .ctp_snpId),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
        )
    }
    # Placement anchors on the GENE position (traitPos), not the weight
    # variants: both genes carry the same five variants here, so only the
    # anchor can tell them apart. block1 spans 100-600, block2 spans 700-1200,
    # so gA @chr1:100 homes to block1 and gB @chr1:800 to block2.
    tw <- TwasWeights(
        studyName = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("gA", "gB"),
        method = c("susie", "susie"),
        entry = list(mkE(), mkE()),
        traitPos = GenomicRanges::GRanges(
            c("chr1", "chr1"),
            IRanges::IRanges(c(100, 800), c(150, 850))
        ),
        ldReference = .ctp_makeHandle()
    )
    grid <- pecotmr:::.ctwasGwasByBlock(.ctp_makeGwasSumstats())
    bucketed <- pecotmr:::.ctwasBucketWeights(tw, grid)
    expect_named(bucketed, c("block1", "block2"), ignore.order = TRUE)
    expect_equal(as.character(bucketed[["block1"]]$trait), "gA")
    expect_equal(as.character(bucketed[["block2"]]$trait), "gB")
    expect_length(ldReference(bucketed[["block1"]]), 1L)
})

test_that(".ctwasBucketWeights: errors when source lacks traitPos", {
    expect_error(
        pecotmr:::.ctwasBucketWeights(
            .ctp_makeTwasWeights(),
            list(
                chr1_1_350 = .ctp_makeGwasSumstats(),
                chr1_350_700 = .ctp_makeGwasSumstats()
            )
        ),
        "no `traitPos` provenance"
    )
})

# ---------------------------------------------------------------------------
# .ctwasParseGeneIds — id -> identity components
# ---------------------------------------------------------------------------
test_that(".ctwasParseGeneIds: splits region|study|context|trait|method", {
    p <- pecotmr:::.ctwasParseGeneIds(c(
        "blk|Q1|brain|gA|susie",
        "blk2|Q1|liver|gB|lasso"
    ))
    expect_equal(p$rid, c("blk", "blk2"))
    expect_equal(p$study, c("Q1", "Q1"))
    expect_equal(p$context, c("brain", "liver"))
    expect_equal(p$trait, c("gA", "gB"))
    expect_equal(p$method, c("susie", "lasso"))
    expect_error(pecotmr:::.ctwasParseGeneIds("too|few|fields"), "malformed")
})

# ---------------------------------------------------------------------------
# .ctwasRunToRows — decompose a cTWAS run into per-context CtwasResult rows
# ---------------------------------------------------------------------------
.ctp_runResult <- function(
    ids,
    pips,
    prior,
    snpIds = character(0),
    region = data.frame(region_id = "blk")
) {
    geneFm <- data.frame(
        id = ids,
        type = "gene",
        susie_pip = pips,
        stringsAsFactors = FALSE
    )
    snpFm <- if (length(snpIds)) {
        data.frame(
            id = snpIds,
            type = "SNP",
            susie_pip = seq(0.01, by = 0.01, length.out = length(snpIds)),
            stringsAsFactors = FALSE
        )
    } else {
        NULL
    }
    fm <- rbind(geneFm, snpFm)
    list(
        weights = setNames(as.list(seq_along(ids)), ids),
        finemap_res = fm,
        susie_alpha_res = cbind(fm, susie_alpha = 0.5),
        param = list(group_prior = prior),
        region_info = region
    )
}

test_that(".ctwasRunToRows: single-context run -> one row, no jointContexts", {
    run <- .ctp_runResult(
        c("blk|Q1|c1|gA|susie", "blk|Q1|c1|gB|susie"),
        c(0.9, 0.2),
        c(c1 = 0.01, SNP = 1e-4)
    )
    rows <- pecotmr:::.ctwasRunToRows(run, gwasStudy = "D1", method = "susie")
    expect_length(rows, 1L)
    expect_equal(rows[[1L]]$context, "c1")
    expect_equal(rows[[1L]]$study, "Q1")
    expect_equal(rows[[1L]]$gwasStudy, "D1")
    expect_true(is.na(rows[[1L]]$jointContexts))
    expect_equal(nrow(ctwasPosteriors(rows[[1L]]$entry)), 2L)
    expect_equal(
        ctwasGroupPriors(rows[[1L]]$entry)$group_prior[["c1"]],
        0.01
    )
})

test_that(".ctwasRunToRows: multi-context run -> per-context rows sharing jointContexts + param", {
    ids <- c(
        "blk|Q1|brain|gA|susie",
        "blk|Q1|brain|gB|susie",
        "blk|Q1|liver|gA|susie",
        "blk|Q1|liver|gB|susie"
    )
    run <- .ctp_runResult(
        ids,
        c(0.9, 0.1, 0.8, 0.2),
        c(brain = 0.01, liver = 0.02, SNP = 1e-4)
    )
    rows <- pecotmr:::.ctwasRunToRows(run, gwasStudy = "D1", method = "susie")
    expect_length(rows, 2L)
    expect_setequal(
        map_chr(rows, function(r) r$context),
        c("brain", "liver")
    )
    expect_true(all(
        map_chr(rows, function(r) r$jointContexts) == "brain,liver"
    ))
    # Each per-context row keeps only its own genes but shares the joint param.
    expect_equal(nrow(ctwasPosteriors(rows[[1L]]$entry)), 2L)
    expect_named(
        ctwasGroupPriors(rows[[1L]]$entry)$group_prior,
        c("brain", "liver", "SNP")
    )
})

test_that(".ctwasRunToRows: rejects multi-context runs with unshared gene sets", {
    run <- .ctp_runResult(
        c("blk|Q1|brain|gA|susie", "blk|Q1|liver|gC|susie"),
        c(0.9, 0.8),
        c(brain = 0.01, liver = 0.02, SNP = 1e-4)
    )
    expect_error(
        pecotmr:::.ctwasRunToRows(run, "D1", "susie"),
        "SAME gene set in every context"
    )
})

test_that(".ctwasRunToRows: empty finemap_res still yields the modeled row", {
    run <- list(
        weights = setNames(list(1), "blk|Q1|c1|gA|susie"),
        finemap_res = NULL,
        param = list(),
        region_info = NULL
    )
    rows <- pecotmr:::.ctwasRunToRows(run, "D1", "susie")
    expect_length(rows, 1L)
    expect_null(ctwasPosteriors(rows[[1L]]$entry))
    expect_null(susieAlpha(rows[[1L]]$entry))
})

test_that(".ctwasRunToRows: susieAlpha subset mirrors finemap per context", {
    run <- .ctp_runResult(
        c("blk|Q1|c1|gA|susie", "blk|Q1|c1|gB|susie"),
        c(0.9, 0.2),
        c(c1 = 0.01, SNP = 1e-4)
    )
    rows <- pecotmr:::.ctwasRunToRows(run, "D1", "susie")
    sa <- susieAlpha(rows[[1L]]$entry)
    expect_equal(nrow(sa), 2L)
    expect_true("susie_alpha" %in% names(sa))
})

test_that(".ctwasRunToRows: keepSnps adds a dedicated SNP row; default drops SNPs", {
    ids <- c("blk|Q1|c1|gA|susie", "blk|Q1|c1|gB|susie")
    run <- .ctp_runResult(
        ids,
        c(0.9, 0.2),
        c(c1 = 0.01, SNP = 1e-4),
        snpIds = c("chr1:1:G:A", "chr1:2:C:T")
    )
    # default: SNP background dropped, only the gene (context) row.
    rowsOff <- pecotmr:::.ctwasRunToRows(run, "D1", "susie")
    expect_length(rowsOff, 1L)
    expect_true(all(ctwasPosteriors(rowsOff[[1L]]$entry)$type == "gene"))
    # keepSnps: one extra study = context = "SNP" row carrying the SNP background.
    rowsOn <- pecotmr:::.ctwasRunToRows(run, "D1", "susie", keepSnps = TRUE)
    expect_length(rowsOn, 2L)
    snpRow <- rowsOn[[2L]]
    expect_equal(snpRow$study, "SNP")
    expect_equal(snpRow$context, "SNP")
    # The payloads are ranged now: SNP rows carry the coordinates their own
    # variant ids encode, so the count is length() and the columns are mcols.
    fm <- ctwasPosteriors(snpRow$entry)
    expect_s4_class(fm, "GRanges")
    expect_length(fm, 2L)
    expect_true(all(S4Vectors::mcols(fm)$type == "SNP"))
    expect_length(susieAlpha(snpRow$entry), 2L)
})

test_that(".ctwasFilterMethod: subsets rows to the requested method", {
    tw <- TwasWeights(
        studyName = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("mrash", "susie"),
        entry = list(
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.1, 0.2, 0.3)
            ),
            TwasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.4, 0.5, 0.6)
            )
        ),
        ldReference = .ctp_makeHandle()
    )
    twSub <- pecotmr:::.ctwasFilterMethod(tw, "susie")
    expect_equal(nrow(twSub), 1L)
    expect_equal(as.character(twSub$method), "susie")
})

# Build an ldPanel fixture (the materialized shape, R included)
# shape) for the 6-SNP toy panel from .ctp_makeHandle().
.ctp_makeLdPanel <- function(snp_n = 6L) {
    h <- .ctp_makeHandle(snp_n = snp_n)
    snpInfo <- pecotmr:::.ctwasSnpInfoForBlock(h)
    R <- diag(1, snp_n)
    dimnames(R) <- list(snpInfo$id, snpInfo$id)
    # Unit dosage variance — sqrt(1) = 1, so the variance-scaling step in
    # .ctwasBuildWeights is a no-op for this fixture.
    variance <- setNames(rep(1, snp_n), snpInfo$id)
    list(R = R, snpInfo = snpInfo, variance = variance)
}

test_that(".ctwasBuildWeights: only raw weights are scaled by sqrt(variance)", {
    panel <- .ctp_makeLdPanel()
    # Replace the default unit variance with non-trivial values; raw
    # weights should be multiplied by sqrt(variance) before reaching the
    # final wgt matrix, and standardized weights bypass the scaling.
    panel$variance <- setNames(c(0.5, 1, 2, 4, 8, 16), panel$snpInfo$id)
    ids5 <- map_chr(1:5, .ctp_snpId)
    rawW <- c(0.1, 0.2, 0.3, 0.4, 0.5)
    for (standardized in c(FALSE, TRUE)) {
        tw <- TwasWeights(
            studyName = "Q1",
            context = "c1",
            trait = "t1",
            method = "susie",
            entry = list(TwasWeightsRow(
                variantIds = ids5,
                weights = rawW,
                weightStandardized = standardized
            )),
            ldReference = .ctp_makeHandle()
        )
        wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
        expected <- if (standardized) {
            rawW
        } else {
            unname(rawW * sqrt(panel$variance[ids5]))
        }
        expect_equal(
            as.numeric(wl[[1L]]$wgt),
            expected,
            tolerance = 1e-12,
            info = str_c("weightStandardized = ", standardized)
        )
    }
})

# Build an LD-panel fixture with realistic chr:pos:A2:A1 variant IDs so
# `.ctwasHarmonizeWeights` (which calls parseVariantId + harmonizeAlleles)
# can do real allele matching.
.ctp_makeAllelePanel <- function() {
    ids <- c("1:100:C:T", "1:200:G:A", "1:300:A:G", "1:400:T:C")
    snpInfo <- data.frame(
        chrom = 1L,
        id = ids,
        pos = c(100L, 200L, 300L, 400L),
        alt = c("T", "A", "G", "C"), # A1 (effect)
        ref = c("C", "G", "A", "T"), # A2 (other)
        stringsAsFactors = FALSE
    )
    R <- diag(1, 4)
    dimnames(R) <- list(ids, ids)
    list(R = R, snpInfo = snpInfo, variance = setNames(rep(1, 4), ids))
}

test_that(".ctwasHarmonizeWeights: sign-flips weights for swapped A1/A2", {
    panel <- .ctp_makeAllelePanel()
    refVariants <- data.frame(
        chrom = panel$snpInfo$chrom,
        pos = panel$snpInfo$pos,
        A2 = panel$snpInfo$ref,
        A1 = panel$snpInfo$alt,
        variant_id = panel$snpInfo$id,
        stringsAsFactors = FALSE
    )
    # Variant 1: alleles match panel ("1:100:C:T" — same A2/A1 ordering)
    # Variant 2: A1/A2 swapped vs panel ("1:200:A:G" — flips relative to "1:200:G:A")
    # Output variant_id values are rebuilt via formatVariantId(), which
    # always emits a `chr` prefix.
    res <- pecotmr:::.ctwasHarmonizeWeights(
        origVids = c("1:100:C:T", "1:200:A:G"),
        origW = c(0.5, 0.3),
        refVariants = refVariants
    )
    expect_equal(nrow(res), 2L)
    # Variant 1 keeps its sign; variant 2 should be sign-flipped to -0.3.
    matches <- match(c("chr1:100:C:T", "chr1:200:G:A"), res$variant_id)
    expect_equal(res$w[matches[[1L]]], 0.5, tolerance = 1e-12)
    expect_equal(res$w[matches[[2L]]], -0.3, tolerance = 1e-12)
})

test_that(".ctwasHarmonizeWeights: drops variants not present in the panel", {
    panel <- .ctp_makeAllelePanel()
    refVariants <- data.frame(
        chrom = panel$snpInfo$chrom,
        pos = panel$snpInfo$pos,
        A2 = panel$snpInfo$ref,
        A1 = panel$snpInfo$alt,
        variant_id = panel$snpInfo$id,
        stringsAsFactors = FALSE
    )
    res <- pecotmr:::.ctwasHarmonizeWeights(
        origVids = c("1:100:C:T", "1:999:A:T"), # 1:999 not in panel
        origW = c(0.5, 0.3),
        refVariants = refVariants
    )
    expect_equal(nrow(res), 1L)
    expect_equal(res$variant_id, "chr1:100:C:T")
})

test_that(".ctwasIsSusieFit: recognizes the susie intermediate shape", {
    fits <- list(
        alpha = matrix(1 / 3, 2, 3),
        mu = matrix(0, 2, 3),
        X_column_scale_factors = rep(1, 3)
    )
    expect_true(pecotmr:::.ctwasIsSusieFit(fits))
    expect_false(pecotmr:::.ctwasIsSusieFit(NULL))
    # lbf_variable alone no longer qualifies: renormalization needs alpha.
    expect_false(pecotmr:::.ctwasIsSusieFit(list(
        lbf_variable = matrix(0, 2, 3),
        mu = matrix(0, 2, 3),
        X_column_scale_factors = rep(1, 3)
    )))
})

test_that(".ctwasRenormalizeSusieWeights: signFlip carries over to mu", {
    # All variants kept; harmonized weights have opposite sign on v2.
    origVids <- sprintf("chr1:%d:A:G", 100L * (1:3))
    origW <- c(0.1, 0.2, 0.3)
    # Strong lbf concentrated on a single effect / single variant per row,
    # so the recomputed weight directly mirrors mu (alpha ≈ identity rows).
    lbf <- rbind(c(100, -100, -100), c(-100, 100, -100))
    mu <- rbind(c(1, 2, 3), c(4, 5, 6))
    fits <- list(
        alpha = lbfToAlpha(lbf),
        mu = mu,
        X_column_scale_factors = rep(1, 3)
    )
    harmW_noflip <- c(0.1, 0.2, 0.3) # all positive
    harmW_v2flip <- c(0.1, -0.2, 0.3) # v2 flipped
    outNoFlip <- pecotmr:::.ctwasRenormalizeSusieWeights(
        fits,
        origVids,
        origW,
        keptIdx = 1:3,
        harmonizedW = harmW_noflip
    )
    outV2flip <- pecotmr:::.ctwasRenormalizeSusieWeights(
        fits,
        origVids,
        origW,
        keptIdx = 1:3,
        harmonizedW = harmW_v2flip
    )
    # v1, v3 should be unchanged between the two; v2 should flip sign.
    expect_equal(outNoFlip[[1L]], outV2flip[[1L]], tolerance = 1e-9)
    expect_equal(outNoFlip[[3L]], outV2flip[[3L]], tolerance = 1e-9)
    expect_equal(outNoFlip[[2L]], -outV2flip[[2L]], tolerance = 1e-9)
})

test_that(".ctwasRenormalizeSusieWeights: honours a non-uniform prior", {
    # The stored alpha already carries the fit's prior, so restricting and
    # renormalizing it is exact. Rebuilding alpha from lbf_variable -- what
    # this used to do -- substitutes a uniform prior and lands elsewhere.
    origVids <- sprintf("chr1:%d:A:G", 100L * (1:4))
    origW <- c(0.1, 0.2, 0.3, 0.4)
    lbf <- rbind(c(0.5, -0.2, 1.1, 0.3), c(-1.0, 0.8, 0.2, 0.4))
    prior <- c(0.6, 0.2, 0.1, 0.1)
    w <- exp(sweep(lbf, 1L, apply(lbf, 1L, max), `-`)) * rep(prior, each = 2L)
    alpha <- w / rowSums(w)
    mu <- rbind(c(1, 2, 3, 4), c(1, 2, 3, 4))
    fits <- list(
        alpha = alpha,
        lbf_variable = lbf,
        mu = mu,
        X_column_scale_factors = rep(1, 4)
    )
    keptIdx <- c(1L, 2L, 3L)
    out <- pecotmr:::.ctwasRenormalizeSusieWeights(
        fits,
        origVids = origVids,
        origW = origW,
        keptIdx = keptIdx,
        harmonizedW = origW[keptIdx]
    )
    expected <- alpha[, keptIdx, drop = FALSE]
    expected <- expected / rowSums(expected)
    expect_equal(out, as.numeric(colSums(expected * mu[, keptIdx])))

    # The old uniform-prior route disagrees on the same fit.
    uniform <- lbfToAlpha(lbf[, keptIdx, drop = FALSE])
    expect_gt(
        max(abs(as.numeric(colSums(uniform * mu[, keptIdx])) - out)),
        1e-6
    )
})

test_that(".ctwasSnpInfoForGwasBlock: restricts panel snpInfo to block GWAS variants", {
    # Block 1's GWAS rows (variants 1..6) against the whole 12-variant panel.
    ss <- .ctp_makeGwasSumstats(blockIds = "block1")
    panelInfo <- data.frame(
        chrom = 1L,
        id = map_chr(1:12, .ctp_snpId), # whole panel
        pos = seq(100L, by = 100L, length.out = 12L),
        alt = "A",
        ref = "G",
        stringsAsFactors = FALSE
    )
    blockInfo <- pecotmr:::.ctwasSnpInfoForGwasBlock(ss, panelInfo)
    # Restricted to the variant IDs on the GwasSumStats entry, in panel order.
    expect_equal(blockInfo$id, map_chr(1:6, .ctp_snpId))
})

test_that(".ctwasBuildWeights: keys per-tuple weights, adds gene metadata", {
    tw <- .ctp_makeTwasWeights()
    panel <- .ctp_makeLdPanel()
    ids5 <- map_chr(1:5, .ctp_snpId)
    wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
    expect_equal(length(wl), 1L)
    expect_equal(names(wl), "Q1|c1|t1|susie")
    expect_equal(wl[[1L]]$study, "Q1")
    expect_equal(wl[[1L]]$context, "c1")
    expect_equal(wl[[1L]]$gene_name, "t1")
    # wgt is a variants x 1 matrix with rownames = SNP IDs
    expect_true(is.matrix(wl[[1L]]$wgt))
    expect_equal(dim(wl[[1L]]$wgt), c(5L, 1L))
    expect_equal(rownames(wl[[1L]]$wgt), ids5)
    # wgt_var is the gene's w'Rw, here from the panel's own R; no weight LD
    # matrix is carried.
    w <- wl[[1L]]$wgt[, 1L]
    expect_equal(
        wl[[1L]]$wgt_var,
        as.numeric(crossprod(w, panel$R[ids5, ids5] %*% w))
    )
    expect_null(wl[[1L]]$R_wgt)
    expect_equal(wl[[1L]]$n_wgt, 5L)
})

test_that(".ctwasBuildWeights: drops variants not present in the LD panel", {
    ids3 <- map_chr(1:3, .ctp_snpId)
    missing <- c("chr1:99900:G:A", "chr1:99910:G:A") # not in panel
    tw <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(TwasWeightsRow(
            variantIds = c(ids3, missing),
            weights = c(0.1, 0.2, 0.3, 0.4, 0.5)
        )),
        ldReference = .ctp_makeHandle()
    )
    panel <- .ctp_makeLdPanel()
    wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
    expect_equal(nrow(wl[[1L]]$wgt), 3L)
    expect_equal(rownames(wl[[1L]]$wgt), ids3)
    expect_equal(wl[[1L]]$n_wgt, 3L)
})

test_that(".ctwasBuildWeights: intersects with gwasSnpIds when supplied", {
    # The LD panel covers ids 1..6; the per-block GWAS sumstats covers only
    # ids 1, 2, 4 (a subset). Weight variants that live in the panel but
    # outside the block (id 3 here) must be dropped, otherwise ctwas's
    # compute_gene_z asserts the weight variant is missing from z_snp.
    ids5 <- map_chr(1:5, .ctp_snpId)
    blockIds <- map_chr(c(1, 2, 4), .ctp_snpId)
    tw <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(TwasWeightsRow(
            variantIds = ids5,
            weights = c(0.1, 0.2, 0.3, 0.4, 0.5)
        )),
        ldReference = .ctp_makeHandle()
    )
    panel <- .ctp_makeLdPanel()
    wl <- pecotmr:::.ctwasBuildWeights(
        tw,
        panel,
        gwasSnpIds = blockIds
    )
    expect_equal(rownames(wl[[1L]]$wgt), blockIds)
    expect_equal(wl[[1L]]$n_wgt, 3L)
})

test_that(".ctwasPanelFor: cheap fields eagerly, correlation on request", {
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- pecotmr:::.ctwasPanelFor(.ctp_makeHandle(snp_n = 6L))
    ids6 <- map_chr(1:6, .ctp_snpId)
    # The n x n correlation is deliberately NOT built here: ctwas pulls a
    # region's LD only while fine-mapping, and only for the regions that
    # survive screening.
    expect_named(out, c("snpInfo", "variance", "sketch"))
    expect_null(out$R)
    expect_setequal(
        colnames(out$snpInfo),
        c("chrom", "id", "pos", "alt", "ref")
    )
    expect_named(out$variance, ids6)
    # It is built, with the right shape, when something asks.
    R <- pecotmr:::.ctwasPanelLd(out)
    expect_true(is.matrix(R))
    expect_equal(dim(R), c(6L, 6L))
    expect_equal(rownames(R), ids6)
})

test_that(".ctwasPanelWeightLd: refuses a variant absent from the panel", {
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    panel <- pecotmr:::.ctwasPanelFor(.ctp_makeHandle(snp_n = 6L))
    vids <- c(.ctp_snpId(2L), "chr22:999999999:A:G")
    # A missing id would give an NA column that silently poisons the weight
    # variance, so it is refused by name and count.
    expect_error(
        pecotmr:::.ctwasPanelWeightLd(panel, vids),
        "1 weight variant\\(s\\) absent from the block's LD panel"
    )
})

test_that(".ctwasPanelWeightLd: a panel's factor, or the R it carries", {
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    panel <- pecotmr:::.ctwasPanelFor(.ctp_makeHandle(snp_n = 6L))
    vids <- map_chr(c(4L, 1L, 3L), .ctp_snpId)
    full <- pecotmr:::.ctwasPanelLd(panel)
    # A sketch gives its factor: the LD over vids without forming it.
    viaFactor <- pecotmr:::.ctwasPanelWeightLd(panel, vids)
    expect_null(viaFactor$R)
    expect_identical(colnames(viaFactor$F), vids)
    expect_equal(
        crossprod(viaFactor$F) / (nrow(viaFactor$F) - 1),
        full[vids, vids]
    )
    # A panel that already carries R is sliced, not re-read.
    withR <- pecotmr:::.ctwasPanelWeightLd(c(panel, list(R = full)), vids)
    expect_null(withR$F)
    expect_equal(withR$R, full[vids, vids])
    # An LdMatrix has no factor to give: its matrix is read over vids.
    ldPanel <- pecotmr:::.ctwasPanelFor(LdMatrix(full, nSamples = 30))
    viaMatrix <- pecotmr:::.ctwasPanelWeightLd(ldPanel, vids)
    expect_null(viaMatrix$F)
    expect_equal(viaMatrix$R, full[vids, vids])
})

test_that(".ctwasLdMapRegionIds: falls back to positions without region_id", {
    # The granular entry points accept hand-built LD_map stubs, which need
    # not carry region_id -- and a region_id of the wrong length cannot be
    # paired with the resolved blocks either.
    expect_equal(
        pecotmr:::.ctwasLdMapRegionIds(data.frame(x = 1:3), 3L),
        c("1", "2", "3")
    )
    expect_equal(
        pecotmr:::.ctwasLdMapRegionIds(data.frame(region_id = "a"), 3L),
        c("1", "2", "3")
    )
    expect_equal(
        pecotmr:::.ctwasLdMapRegionIds(
            data.frame(region_id = c("a", "b")),
            2L
        ),
        c("a", "b")
    )
})

test_that(".ctwasGlobalPanelInfo: NULL when the collection has no sketch", {
    # The shared fixture carries a sketch, so build one without: the global
    # snpInfo table is only derivable when there is a panel behind it.
    gs <- GwasSumStats(
        studyName = "G1",
        entry = list(.ctp_blockVariants(1L)),
        genome = "hg19",
        blockId = "block1",
        qcOptions = list(step1 = "ok")
    )
    expect_identical(ldReference(gs), list())
    # .ctwasRegionGwasLd reports the missing-sketch case with its own
    # message, so this helper must stay quiet and hand back NULL.
    expect_null(pecotmr:::.ctwasGlobalPanelInfo(gs))
})

test_that(".ctwasWeightVariances: w'Rw from the factor matches the matrix", {
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    panel <- pecotmr:::.ctwasPanelFor(.ctp_makeHandle(snp_n = 6L))
    full <- pecotmr:::.ctwasPanelLd(panel)
    prepared <- list(
        list(vids = map_chr(c(1L, 2L), .ctp_snpId), w = c(0.5, -0.2)),
        list(vids = map_chr(c(2L, 5L, 6L), .ctp_snpId), w = c(0.1, 0.3, -0.4))
    )
    viaFactor <- pecotmr:::.ctwasWeightVariances(panel, prepared)
    viaMatrix <- pecotmr:::.ctwasWeightVariances(
        c(panel, list(R = full)),
        prepared
    )
    expect_equal(viaFactor, viaMatrix)
    expect_equal(
        viaMatrix[[2L]],
        as.numeric(crossprod(
            prepared[[2L]]$w,
            full[prepared[[2L]]$vids, prepared[[2L]]$vids] %*% prepared[[2L]]$w
        ))
    )
    expect_identical(pecotmr:::.ctwasWeightVariances(panel, list()), numeric(0))
})

test_that(".ctwasGeneZ gives ctwas::compute_gene_z's gene z-scores", {
    skip_if_not_installed("ctwas")
    ids <- map_chr(1:3, .ctp_snpId)
    R <- matrix(c(1, 0.4, 0.1, 0.4, 1, 0.3, 0.1, 0.3, 1), 3, 3)
    dimnames(R) <- list(ids, ids)
    wgt <- matrix(c(0.5, -0.2, 0.3), ncol = 1L, dimnames = list(ids, "wgt"))
    entry <- list(wgt = wgt, type = "c1", context = "c1")
    zSnp <- data.frame(id = ids, A1 = "A", A2 = "G", z = c(2, -1, 0.5))
    ours <- pecotmr:::.ctwasGeneZ(
        list(
            g1 = list_assign(
                entry,
                wgt_var = as.numeric(crossprod(wgt, R %*% wgt))
            )
        ),
        zSnp
    )
    theirs <- suppressMessages(ctwas::compute_gene_z(
        zSnp,
        list(g1 = list_assign(entry, R_wgt = R))
    ))
    expect_equal(ours, theirs)
    expect_error(
        pecotmr:::.ctwasGeneZ(
            list(g1 = list_assign(entry, wgt_var = 1)),
            zSnp[-1L, ]
        ),
        "1 weight variant\\(s\\) of gene 'g1' have no GWAS z-score"
    )
})

test_that(".ctwasGeneZ: none without weights; an infinite z is dropped", {
    expect_null(pecotmr:::.ctwasGeneZ(list(), data.frame()))
    # A gene whose weights carry no variance in the panel has an infinite
    # z, which ctwas's own compute_gene_z() drops too.
    ids <- map_chr(1:2, .ctp_snpId)
    wgt <- matrix(c(0.5, -0.2), ncol = 1L, dimnames = list(ids, "wgt"))
    entry <- list(wgt = wgt, type = "c1", context = "c1", wgt_var = 1)
    zSnp <- data.frame(id = ids, A1 = "A", A2 = "G", z = c(2, -1))
    expect_warning(
        out <- pecotmr:::.ctwasGeneZ(
            list(g1 = entry, g2 = list_assign(entry, wgt_var = 0)),
            zSnp
        ),
        "dropped 1 gene\\(s\\) with an infinite z-score"
    )
    expect_equal(out$id, "g1")
    expect_equal(out$z, 0.5 * 2 + 0.2)
})

test_that(".ctwasBuildZGene: builds z_gene from a TWAS-Z GRanges", {
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        qtlStudy = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        twasZ = 1.5
    )
    df <- pecotmr:::.ctwasBuildZGene(gr)
    expect_equal(nrow(df), 1L)
    expect_setequal(
        colnames(df),
        c("id", "z", "type", "context", "gene_name", "study", "method")
    )
    expect_equal(df$id, "Q1|c1|t1|susie")
})

# ===========================================================================
# LD loader / SNP-info loader closures
# ===========================================================================

test_that(".ctwasMultiBlockLdLoader: dispatches by LD_file token", {
    RA <- matrix(
        runif(4),
        2,
        2,
        dimnames = list(
            c("chr1:100:A:G", "chr1:200:A:G"),
            c("chr1:100:A:G", "chr1:200:A:G")
        )
    )
    RB <- matrix(
        runif(4),
        2,
        2,
        dimnames = list(
            c("chr1:300:A:G", "chr1:400:A:G"),
            c("chr1:300:A:G", "chr1:400:A:G")
        )
    )
    loader <- pecotmr:::.ctwasMultiBlockLdLoader(
        list(tokenA = list(R = RA), tokenB = list(R = RB))
    )
    expect_identical(loader("tokenA"), RA)
    expect_identical(loader("tokenB"), RB)
    expect_error(loader("unknown_token"), "no cached panel")
})

test_that(".ctwasMultiBlockSnpInfoLoader: dispatches by LD_file token", {
    infoA <- data.frame(
        chrom = 1L,
        id = c("chr1:100:A:G", "chr1:200:A:G"),
        pos = c(100L, 200L),
        alt = "A",
        ref = "G",
        stringsAsFactors = FALSE
    )
    infoB <- data.frame(
        chrom = 1L,
        id = c("chr1:300:A:G", "chr1:400:A:G"),
        pos = c(300L, 400L),
        alt = "C",
        ref = "T",
        stringsAsFactors = FALSE
    )
    loader <- pecotmr:::.ctwasMultiBlockSnpInfoLoader(
        list(tokenA = list(snpInfo = infoA), tokenB = list(snpInfo = infoB))
    )
    expect_identical(loader("tokenA"), infoA)
    expect_identical(loader("tokenB"), infoB)
    expect_error(loader("unknown_token"), "no cached panel")
})

# ===========================================================================
# Input-assembly shape checks via assembleCtwasInputs (no ctwas engine needed)
# ===========================================================================

test_that("assembleCtwasInputs: assembles the documented input shape for ctwas", {
    inp <- .ctp_makeMultiBlockInputs()
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    inputs <- assembleCtwasInputs(
        gwasSumStats = inp$gwasSumStats,
        twasWeights = inp$twasWeights,
        method = "susie"
    )
    # Two regions, two LD_map rows.
    expect_equal(nrow(inputs$region_info), 2L)
    expect_setequal(inputs$region_info$region_id, c("block1", "block2"))
    # zSnp is the concatenation of both blocks' Z columns.
    expect_equal(nrow(inputs$z_snp), 12L)
    # snp_map keyed by region_id.
    expect_setequal(names(inputs$snp_map), c("block1", "block2"))
    # Per-region weights keys are prefixed with the region_id.
    expect_true(all(grepl("^(block1|block2)\\|", names(inputs$weights))))
    # LD_map carries the same number of rows as regions.
    expect_equal(nrow(inputs$LD_map), 2L)
    # LD / snpInfo loader closures are present.
    expect_true(is.function(inputs$LD_loader_fun))
    expect_true(is.function(inputs$snpinfo_loader_fun))
})

test_that("assembleCtwasInputs: forwards a twasZ argument as z_gene", {
    inp <- .ctp_makeMultiBlockInputs()
    twasZ <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
    S4Vectors::mcols(twasZ) <- S4Vectors::DataFrame(
        qtlStudy = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        twasZ = 1.5
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    inputs <- assembleCtwasInputs(
        gwasSumStats = inp$gwasSumStats,
        twasWeights = inp$twasWeights,
        twasZ = twasZ,
        method = "susie"
    )
    expect_equal(inputs$z_gene$id, "Q1|c1|t1|susie")
})

# ===========================================================================
# Step-wise dispatch: estCtwasGroupPriors → screenCtwasRegions → finemapCtwasRegions
# ===========================================================================

test_that("ctwasPipeline: dispatches assemble → est → screen → finemap and accumulates state", {
    skip_if_not_installed("ctwas")
    inp <- .ctp_makeMultiBlockInputs()
    capturedAssemble <- list()
    capturedEst <- list()
    capturedScreen <- list()
    capturedFinemap <- list()
    local_mocked_bindings(
        assemble_region_data = function(...) {
            capturedAssemble <<- list(...)
            list(block1 = list(stub = TRUE))
        },
        get_boundary_genes = function(...) {
            data.frame(id = "t1", n_regions = 2L)
        },
        est_param = function(...) {
            capturedEst <<- list(...)
            list(
                group_prior = c(g = 0.1, SNP = 0.0001),
                group_prior_var = c(g = 5, SNP = 5)
            )
        },
        screen_regions = function(...) {
            capturedScreen <<- list(...)
            list(
                screened_region_data = list(block1 = list(stub = TRUE)),
                screen_res_meta = "mocked"
            )
        },
        finemap_regions = function(...) {
            capturedFinemap <<- list(...)
            # Return per-gene finemap rows keyed by our cTWAS gene ids
            # (region|study|context|trait|method) so the decomposition into
            # CtwasResult rows finds a matching finemap payload.
            list(
                finemap_res = data.frame(
                    id = c("block1|Q1|c1|t1|susie", "block2|Q1|c1|t1|susie"),
                    type = c("gene", "gene"),
                    context = c("c1", "c1"),
                    susie_pip = c(0.9, 0.8),
                    stringsAsFactors = FALSE
                ),
                susie_alpha_res = data.frame(
                    id = c("block1|Q1|c1|t1|susie", "block2|Q1|c1|t1|susie"),
                    type = c("gene", "gene"),
                    susie_alpha = c(0.9, 0.8),
                    stringsAsFactors = FALSE
                )
            )
        },
        .package = "ctwas"
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- ctwasPipeline(
        gwasSumStats = inp$gwasSumStats,
        twasWeights = inp$twasWeights,
        method = "susie"
    )
    # Each ctwas step was invoked exactly once.
    expect_true(length(capturedAssemble) > 0L)
    expect_true(length(capturedEst) > 0L)
    expect_true(length(capturedScreen) > 0L)
    expect_true(length(capturedFinemap) > 0L)
    # est sees the assembled region_data.
    expect_named(capturedEst$region_data, "block1")
    # screen receives the param estimates as group_prior / group_prior_var.
    expect_equal(unname(capturedScreen$group_prior), c(0.1, 0.0001))
    # finemap consumes screen's screened_region_data.
    expect_named(capturedFinemap$region_data, "block1")
    # Output is a CtwasResult: one row per (run, region), and this run's
    # genes, placed without a region_id, share one chromosome.
    expect_s4_class(out, "CtwasResult")
    expect_equal(nrow(out), 1L)
    expect_equal(methodNames(out), "susie")
    expect_equal(studyName(out), "Q1")
    expect_equal(contexts(out), "c1")
    expect_equal(as.character(out$gwasStudy), "G1")
    # The run's jointly-estimated param is carried on the row.
    expect_equal(unname(out$groupPriors[[1L]]$group_prior), c(0.1, 0.0001))
    # ctwasPosteriors aggregates the per-gene rows, tagged with run identity.
    fm <- ctwasPosteriors(out)
    expect_s4_class(fm, "GRanges")
    expect_equal(length(fm), 2L)
    expect_true(all(
        c("gwasStudy", "study", "context", "method", "susie_pip") %in%
            names(S4Vectors::mcols(fm))
    ))
    expect_setequal(fm$id, c("block1|Q1|c1|t1|susie", "block2|Q1|c1|t1|susie"))
    # susieAlpha aggregates the fuller per-effect table, likewise tagged.
    sa <- susieAlpha(out)
    expect_equal(length(sa), 2L)
    expect_true(all(
        c("gwasStudy", "context", "susie_alpha") %in%
            names(S4Vectors::mcols(sa))
    ))
})

test_that("ctwasPipeline: keepSnps retains the SNP background as a dedicated row", {
    skip_if_not_installed("ctwas")
    inp <- .ctp_makeMultiBlockInputs()
    local_mocked_bindings(
        assemble_region_data = function(...) list(block1 = list(stub = TRUE)),
        get_boundary_genes = function(...) {
            data.frame(id = "t1", n_regions = 2L)
        },
        est_param = function(...) {
            list(
                group_prior = c(c1 = 0.1, SNP = 1e-4),
                group_prior_var = c(c1 = 5, SNP = 5)
            )
        },
        screen_regions = function(...) {
            list(
                screened_region_data = list(block1 = list(stub = TRUE)),
                screen_res_meta = "m"
            )
        },
        finemap_regions = function(...) {
            list(
                finemap_res = data.frame(
                    id = c("block1|Q1|c1|t1|susie", "chr1:100:G:A"),
                    type = c("gene", "SNP"),
                    susie_pip = c(0.9, 0.02),
                    stringsAsFactors = FALSE
                ),
                susie_alpha_res = data.frame(
                    id = c("block1|Q1|c1|t1|susie", "chr1:100:G:A"),
                    type = c("gene", "SNP"),
                    susie_alpha = c(0.9, 0.02),
                    stringsAsFactors = FALSE
                )
            )
        },
        .package = "ctwas"
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- ctwasPipeline(
        inp$gwasSumStats,
        inp$twasWeights,
        keepSnps = TRUE,
        method = "susie"
    )
    # A gene (c1) row plus the dedicated SNP row.
    expect_equal(nrow(out), 2L)
    expect_setequal(contexts(out), c("c1", "SNP"))
    fm <- ctwasPosteriors(out)
    expect_setequal(fm$type, c("gene", "SNP"))
    expect_true("chr1:100:G:A" %in% fm$id)
})

# ===========================================================================
# mergeCtwasBoundaryRegions + ctwasPipeline(mergeBoundary = TRUE)
# ===========================================================================
.ctp_finemapResult <- function(hasLd = TRUE) {
    fm <- data.frame(
        id = "block1|Q1|c1|t1|susie",
        type = "gene",
        susie_pip = 0.90,
        stringsAsFactors = FALSE
    )
    base <- list(
        finemap_res = fm,
        susie_alpha_res = cbind(fm, susie_alpha = 0.90),
        region_data = list(block1 = list(stub = TRUE)),
        region_info = data.frame(
            region_id = "block1",
            chrom = 1,
            start = 1,
            stop = 1000
        ),
        z_snp = data.frame(id = "chr1:100:G:A", z = 1.0),
        z_gene = data.frame(id = "block1|Q1|c1|t1|susie", z = 2.0),
        # Placed by chrom / p0 / p1, as ctwas's own weights are.
        weights = setNames(
            list(list(wgt = 1, chrom = 1, p0 = 100L, p1 = 100L)),
            "block1|Q1|c1|t1|susie"
        ),
        snp_map = list(block1 = data.frame(id = "chr1:100:G:A")),
        param = list(
            group_prior = c(c1 = 0.1, SNP = 1e-4),
            group_prior_var = c(c1 = 5, SNP = 5)
        ),
        LD_map = data.frame(
            region_id = "block1",
            LD_file = "f",
            SNP_file = "f",
            stringsAsFactors = FALSE
        )
    )
    if (hasLd) {
        base$LD_loader_fun <- function(...) NULL
        base$snpinfo_loader_fun <- function(...) NULL
    }
    base
}

.ctp_mergeReturn <- function() {
    list(
        updated_finemap_res = data.frame(
            id = "block1|Q1|c1|t1|susie",
            type = "gene",
            susie_pip = 0.99,
            stringsAsFactors = FALSE
        ),
        updated_susie_alpha_res = data.frame(
            id = "block1|Q1|c1|t1|susie",
            susie_alpha = 0.99
        ),
        updated_region_data = list(merged = TRUE),
        updated_region_info = data.frame(region_id = "block1_merged"),
        updated_LD_map = data.frame(region_id = "block1_merged"),
        updated_snp_map = list(block1_merged = 1)
    )
}

test_that("mergeCtwasBoundaryRegions passes no unnamed argument to ctwas", {
    skip_if_not_installed("ctwas")
    # Regression: methodArgs was merged with c(), which appends a MethodOptions
    # (a SimpleList) as ONE opaque element instead of splicing its entries.
    # exec() then handed that object to ctwas positionally, binding it to
    # `combine_PIPs` -- so even the default CtwasOptions() broke the real call.
    captured <- NULL
    local_mocked_bindings(
        postprocess_region_merging = function(...) {
            captured <<- list(...)
            .ctp_mergeReturn()
        },
        .package = "ctwas"
    )
    mergeCtwasBoundaryRegions(.ctp_finemapResult(hasLd = TRUE))
    expect_false(any(names(captured) == ""))
    expect_false("combine_PIPs" %in% names(captured))

    # ...and a user option still arrives spliced, by its own name
    mergeCtwasBoundaryRegions(
        .ctp_finemapResult(hasLd = TRUE),
        methodArgs = CtwasOptions(
            merge = CtwasMergeOptions(min_abs_corr = 0.1)
        )
    )
    expect_equal(captured$min_abs_corr, 0.1)
    expect_false(any(names(captured) == ""))
})

test_that("mergeCtwasBoundaryRegions: no-LD path uses postprocess_region_merging_noLD", {
    skip_if_not_installed("ctwas")
    usedNoLd <- FALSE
    local_mocked_bindings(
        postprocess_region_merging = function(...) {
            .ctp_mergeReturn()
        },
        postprocess_region_merging_noLD = function(...) {
            usedNoLd <<- TRUE
            .ctp_mergeReturn()
        },
        .package = "ctwas"
    )
    out <- mergeCtwasBoundaryRegions(.ctp_finemapResult(hasLd = FALSE))
    expect_true(usedNoLd)
    expect_equal(out$finemap_res$susie_pip, 0.99)
})

test_that("ctwasPipeline: mergeBoundary = TRUE re-fine-maps and flows into CtwasResult", {
    skip_if_not_installed("ctwas")
    inp <- .ctp_makeMultiBlockInputs()
    merged <- FALSE
    local_mocked_bindings(
        assemble_region_data = function(...) list(block1 = list(stub = TRUE)),
        get_boundary_genes = function(...) {
            data.frame(id = "t1", n_regions = 2L)
        },
        est_param = function(...) {
            list(
                group_prior = c(c1 = 0.1, SNP = 1e-4),
                group_prior_var = c(c1 = 5, SNP = 5)
            )
        },
        screen_regions = function(...) {
            list(
                screened_region_data = list(block1 = list(stub = TRUE)),
                screen_res_meta = "m"
            )
        },
        finemap_regions = function(...) {
            list(
                finemap_res = data.frame(
                    id = "block1|Q1|c1|t1|susie",
                    type = "gene",
                    susie_pip = 0.90,
                    stringsAsFactors = FALSE
                ),
                susie_alpha_res = data.frame(
                    id = "block1|Q1|c1|t1|susie",
                    susie_alpha = 0.90
                )
            )
        },
        postprocess_region_merging = function(...) {
            merged <<- TRUE
            list(
                updated_finemap_res = data.frame(
                    id = "block1|Q1|c1|t1|susie",
                    type = "gene",
                    susie_pip = 0.99,
                    stringsAsFactors = FALSE
                ),
                updated_susie_alpha_res = data.frame(
                    id = "block1|Q1|c1|t1|susie",
                    susie_alpha = 0.99
                )
            )
        },
        .package = "ctwas"
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- ctwasPipeline(
        inp$gwasSumStats,
        inp$twasWeights,
        boundaryMergeParam = BoundaryMergeParam(enabled = TRUE),
        method = "susie"
    )
    expect_true(merged) # merging ran
    expect_s4_class(out, "CtwasResult")
    # The post-merge PIP (0.99) is what the decomposition carries.
    expect_equal(ctwasPosteriors(out)$susie_pip, 0.99)
})

# ===========================================================================
# asCtwasResult — structure a granular finemap result (the wrapper path)
# ===========================================================================
test_that("asCtwasResult: structures a granular finemap result into a CtwasResult", {
    fmr <- .ctp_finemapResult()
    fmr$z_snp$study <- "D1" # GWAS study read from z_snp
    cr <- asCtwasResult(fmr)
    expect_s4_class(cr, "CtwasResult")
    expect_equal(nrow(cr), 1L)
    expect_equal(as.character(cr$gwasStudy), "D1")
    expect_equal(studyName(cr), "Q1") # derived from the gene ids
    expect_equal(contexts(cr), "c1")
    expect_equal(methodNames(cr), "susie")
    expect_equal(ctwasPosteriors(cr)$susie_pip, 0.90)
    expect_false(is.null(susieAlpha(cr)))
})

test_that("asCtwasResult: keepSnps adds the dedicated SNP row", {
    fmr <- .ctp_finemapResult()
    fmr$z_snp$study <- "D1"
    fmr$finemap_res <- rbind(
        fmr$finemap_res,
        data.frame(
            id = "chr1:100:G:A",
            type = "SNP",
            susie_pip = 0.02,
            stringsAsFactors = FALSE
        )
    )
    cr <- asCtwasResult(fmr, keepSnps = TRUE)
    expect_setequal(contexts(cr), c("c1", "SNP"))
})

test_that("asCtwasResult: errors when the weights mix methods", {
    fmr <- .ctp_finemapResult()
    fmr$weights <- setNames(
        list(1, 2),
        c("block1|Q1|c1|t1|susie", "block1|Q1|c1|t2|lasso")
    )
    expect_error(asCtwasResult(fmr), "mixes weight methods")
})

test_that("estCtwasGroupPriors: fallbackToPrefit recovers from accurate-EM NaN divergence", {
    skip_if_not_installed("ctwas")
    inp <- .ctp_makeMultiBlockInputs()
    # Mock est_param to throw the documented NaN error, and fit_EM to
    # produce a stub prefit result. Verify estCtwasGroupPriors catches the
    # NaN error AND that the returned param is the prefit estimate.
    local_mocked_bindings(
        assemble_region_data = function(...) {
            list(
                block1 = list(
                    gid = "t1",
                    sid = c("s1", "s2"),
                    stub = TRUE
                )
            )
        },
        get_boundary_genes = function(...) {
            data.frame(id = "t1", n_regions = 2L)
        },
        compute_gene_z = function(...) data.frame(id = "t1", z = 1.0),
        est_param = function(...) {
            stop("Estimated group_prior_var contains NAs!")
        },
        # fit_EM is internal to ctwas — mock via the same `ctwas` namespace.
        fit_EM = function(region_data, ...) {
            list(
                group_prior = c(g = 0.05, SNP = 1e-4),
                group_prior_var = c(g = 4.0, SNP = 5.0),
                group_size = c(g = 1, SNP = 100)
            )
        },
        .package = "ctwas"
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    # Assembled once; the two estimation calls below differ only in the
    # fallback setting.
    inputs <- assembleCtwasInputs(
        inp$gwasSumStats,
        inp$twasWeights,
        method = "susie"
    )
    # Without fallback: the NaN error propagates.
    expect_error(
        estCtwasGroupPriors(
            inputs,
            ctwasPriorParam = CtwasPriorParam(fallbackToPrefit = FALSE)
        ),
        "contains NAs"
    )
    # With fallback: prefit estimates are returned as the param.
    # .ctwasFitPrefitEm thin-scales the SNP group_prior (mirroring ctwas's
    # est_param), so the mocked SNP prior 1e-4 emerges as 1e-4 * thin
    # (default thin = 0.1) → 1e-5. The group_prior_var is not thinned.
    est <- estCtwasGroupPriors(
        inputs,
        ctwasPriorParam = CtwasPriorParam(fallbackToPrefit = TRUE)
    )
    expect_equal(unname(est$param$group_prior), c(0.05, 1e-5))
    expect_equal(unname(est$param$group_prior_var), c(4.0, 5.0))
})

test_that("estCtwasGroupPriors fallback drops degenerate regions before fit_EM", {
    # Regression for the ctwas >= 0.6.0 breakage: the prefit fallback used to hand
    # ALL regions to ctwas::fit_EM, so a degenerate region (empty gid/sid, whose
    # `sid` ctwas::extract_region_data now requires) crashed with
    # "regiondata$sid ... target is NULL". The fallback must mirror est_param's
    # min_var / min_gene skip and fit only the qualifying regions.
    skip_if_not_installed("ctwas")
    inp <- .ctp_makeMultiBlockInputs()
    seen <- NULL
    local_mocked_bindings(
        assemble_region_data = function(...) {
            list(
                good = list(gid = "t1", sid = c("s1", "s2")), # 1 gene + 2 SNPs -> kept
                degenerate = list(gid = character(0), sid = NULL)
            )
        }, # no variables    -> dropped
        get_boundary_genes = function(...) {
            data.frame(id = "t1", n_regions = 2L)
        },
        compute_gene_z = function(...) data.frame(id = "t1", z = 1.0),
        est_param = function(...) stop("No regions selected!"),
        fit_EM = function(region_data, ...) {
            seen <<- names(region_data)
            list(
                group_prior = c(g = 0.05, SNP = 1e-4),
                group_prior_var = c(g = 4.0, SNP = 5.0),
                group_size = c(g = 1, SNP = 100)
            )
        },
        .package = "ctwas"
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    est <- estCtwasGroupPriors(
        assembleCtwasInputs(
            inp$gwasSumStats,
            inp$twasWeights,
            method = "susie"
        ),
        ctwasPriorParam = CtwasPriorParam(fallbackToPrefit = TRUE)
    )
    # only the qualifying region reached fit_EM; the degenerate region was filtered
    expect_equal(seen, "good")
    # ...but every region is still accounted for in the returned p_single_effect
    expect_setequal(
        est$param$p_single_effect$region_id,
        c("good", "degenerate")
    )
})

test_that("(real ctwas) prefit fallback skips a degenerate region fit_EM would reject", {
    # No-mock guard for the ctwas >= 0.6.0 contract. Runs the REAL ctwas::fit_EM
    # (via .ctwasFitPrefitEm) on a genuine assemble_region_data fixture with one
    # valid region (1 gene, 108 SNPs) and one degenerate region (0 genes/SNPs,
    # unset `sid`). Handing the degenerate region to fit_EM crashes ctwas >= 0.6.0
    # in extract_region_data ("regiondata$sid ... target is NULL"); the fallback
    # must filter it. Unlike the mocked tests above, this exercises the real engine,
    # so it would catch a FUTURE ctwas contract change (which the mocks cannot).
    skip_if_not_installed("ctwas")
    region_data <- readRDS(test_path(
        "test_data",
        "ctwas_region_data_degenerate.rds"
    ))
    expect_length(region_data, 2L)
    # Three iterations keep the test quick, short of ctwas's convergence.
    expect_warning(
        res <- .ctwasFitPrefitEm(
            region_data,
            niter = 3L,
            groupPriorVarStructure = "shared_all",
            thin = 1,
            threads = 1L
        ),
        "EM did not converge in 3 iterations"
    )
    # the prefit EM ran on the valid region only and returns finite real group priors
    expect_true("SNP" %in% names(res$group_prior))
    expect_true(all(is.finite(res$group_prior)))
    # every region (valid + degenerate) is still listed in p_single_effect
    expect_setequal(res$p_single_effect$region_id, names(region_data))
})

test_that("estCtwasGroupPriors / screenCtwasRegions / finemapCtwasRegions can be called independently", {
    skip_if_not_installed("ctwas")
    inp <- .ctp_makeMultiBlockInputs()
    local_mocked_bindings(
        assemble_region_data = function(...) list(block1 = list(stub = TRUE)),
        get_boundary_genes = function(...) {
            data.frame(id = "t1", n_regions = 2L)
        },
        est_param = function(...) {
            list(
                group_prior = c(g = 0.05, SNP = 1e-4),
                group_prior_var = c(g = 4, SNP = 5)
            )
        },
        screen_regions = function(...) {
            list(
                screened_region_data = list(block1 = list(stub = TRUE))
            )
        },
        finemap_regions = function(...) {
            list(
                finemap_res = data.frame(id = "t1"),
                susie_alpha_res = data.frame(pip = 0.5)
            )
        },
        .package = "ctwas"
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    # Step 1
    inputs <- assembleCtwasInputs(
        inp$gwasSumStats,
        inp$twasWeights,
        method = "susie"
    )
    expect_true("region_info" %in% names(inputs))
    expect_true("LD_loader_fun" %in% names(inputs))
    # Step 2
    est <- estCtwasGroupPriors(inputs)
    expect_true("region_data" %in% names(est))
    expect_true("param" %in% names(est))
    # User can OVERRIDE the estimated priors before screen/finemap — this is
    # the escape hatch for NaN-on-iter-2 EM divergence.
    est$param$group_prior <- c(g = 0.2, SNP = 1e-4)
    est$param$group_prior_var <- c(g = 4.5, SNP = 5)
    # Step 3
    screened <- screenCtwasRegions(est)
    expect_true("screened_region_data" %in% names(screened))
    # Step 4
    final <- finemapCtwasRegions(screened)
    expect_setequal(
        names(final),
        c(
            "z_gene",
            "param",
            "finemap_res",
            "susie_alpha_res",
            "region_data",
            "boundary_genes",
            "screen_res",
            "region_info",
            "z_snp",
            "weights",
            "snp_map",
            "LD_map",
            "LD_loader_fun",
            "snpinfo_loader_fun"
        )
    )
})

# ===========================================================================
# Real-engine end-to-end: drives ctwas::ctwas_sumstats with the bundled
# example PLINK panel + synthetic TwasWeights. Exercises the LD-loader
# and snp-info-loader closure bodies as actually invoked by ctwas
# (mocked tests only construct the closures, never invoke them).
# ===========================================================================

test_that("ctwasPipeline: real-engine end-to-end on the bundled example panel", {
    skip_if_not_installed("ctwas")
    data(gwasSumStatsS4Example)
    data(qtlDatasetExample)
    gss <- gwasSumStatsS4Example
    qd <- qtlDatasetExample
    gh <- pecotmr:::.qtlGenotypeSe(qd)

    # Two 5-variant synthetic genes from the bundled panel, one anchored in
    # each LD block below. cTWAS's EM runs per region and cannot fit a block
    # with no genes in it, so a real partition needs a gene on both sides --
    # the previous per-region list dodged this by putting the SAME GWAS block
    # and the SAME gene under two region ids.
    geneOne <- TwasWeightsRow(
        variantIds = .ctp_panelSnpInfo(gh)$SNP[1:5],
        weights = c(0.1, -0.2, 0.05, 0.0, 0.3)
    )
    geneTwo <- TwasWeightsRow(
        variantIds = .ctp_panelSnpInfo(gh)$SNP[9:13],
        weights = c(0.2, 0.1, -0.15, 0.05, 0.0)
    )
    # A FLAT weight source is placed into blocks by `traitPos`, so each gene
    # needs the span it sits in.
    geneSpan <- function(i) {
        bp <- .ctp_panelSnpInfo(gh)$BP[i]
        GenomicRanges::GRanges(
            str_c("chr", .ctp_panelSnpInfo(gh)$CHR[i][[1L]]),
            IRanges::IRanges(min(bp), max(bp))
        )
    }
    tw <- TwasWeights(
        studyName = rep("study1", 2L),
        context = rep("brain", 2L),
        trait = c("ENSG_example", "ENSG_example2"),
        method = rep("susie", 2L),
        entry = list(geneOne, geneTwo),
        traitPos = c(geneSpan(1:5), geneSpan(9:13)),
        ldReference = gh
    )

    # ctwasPipeline needs >= 2 blocks for the joint EM to have something to
    # estimate against. The GWAS variants are split across two blocks by
    # position; the same TWAS weights are offered to both. The bundled toy
    # panel is still too sparse for the convergence checks, so we relax the
    # gates.
    gssRanges <- range(unlist(gss))
    mid <- start(gssRanges) + floor(IRanges::width(gssRanges) / 2)
    twoBlocks <- GenomicRanges::GRanges(
        as.character(seqnames(gssRanges)),
        IRanges::IRanges(
            c(start(gssRanges), mid + 1L),
            c(mid, end(gssRanges))
        )
    )
    names(twoBlocks) <- c("blockA", "blockB")
    gssTwoBlocks <- GwasSumStats(
        studyName = studyName(gss),
        entry = list(unlist(gss)),
        genome = unique(unname(GenomeInfoDb::genome(gss))),
        ldReference = ldReference(gss),
        qcOptions = qcOptions(gss),
        region = twoBlocks
    )
    res <- suppressMessages(suppressWarnings(
        ctwasPipeline(
            gwasSumStats = gssTwoBlocks,
            twasWeights = tw,
            ctwasPriorParam = CtwasPriorParam(niter = 5L, niterPrefit = 2L),
            # Toy panel: relax the production filters that gate out tiny
            # inputs. `filter_L = FALSE` used to sit here too; no ctwas
            # function has ever had such a formal, so .ctwasInvoke dropped it
            # silently. CtwasOptions() now rejects it.
            methodArgs = CtwasOptions(
                regionData = CtwasRegionDataOptions(min_group_size = 1L),
                estParam = CtwasEstParamOptions(
                    min_group_size = 1L,
                    min_p_single_effect = 0
                )
            ),
            method = "susie"
        )
    ))

    # The one-shot pipeline returns a CtwasResult: a single (brain, study1,
    # susie) row (single-context run, so no jointContexts).
    expect_s4_class(res, "CtwasResult")
    expect_equal(nrow(res), 1L)
    expect_equal(contexts(res), "brain")
    expect_equal(studyName(res), "study1")
    expect_equal(methodNames(res), "susie")
    expect_false("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
    # The gene we passed in was fine-mapped and shows up in the aggregated
    # finemap table, tagged with the run identity.
    fm <- ctwasPosteriors(res)
    expect_gt(length(fm), 0L)
    expect_true(all(
        c("gwasStudy", "study", "context", "method", "susie_pip") %in%
            names(S4Vectors::mcols(fm))
    ))
    # Either gene may survive ctwas's region screen; what matters here is that
    # a fine-mapped gene carries the full run identity in its id.
    # ids are "<blockId>|study|context|trait|method". Either gene may survive
    # ctwas's region screen, so this pins the identity tagging rather than
    # which gene came through.
    expect_true(any(grepl(
        "\\|study1\\|brain\\|ENSG_example2?\\|susie$",
        fm$id
    )))
    expect_true(all(fm$context == "brain"))
    # The fuller per-effect table is retained and reconstructable.
    sa <- susieAlpha(res)
    expect_gt(length(sa), 0L)
    expect_true(all(
        c("susie_pip", "susie_alpha", "region_id") %in%
            names(S4Vectors::mcols(sa))
    ))
})

# ===========================================================================
# .ctwasFilterVariants
# ===========================================================================
# The filter has four knobs:
#   1. twasWeightCutoff — drop |w| < cutoff
#   2. csMinCor         — high-purity CS rescue (must-keep)
#   3. minPipCutoff     — high-PIP rescue (must-keep)
#   4. maxNumVariants   — per-gene cap, prioritized by PIP then |w|

test_that(".ctwasFilterVariants: twasWeightCutoff drops low-magnitude variants", {
    vids <- sprintf("chr1:%d:A:G", 100L * (1:6))
    w <- c(0.5, 0.001, 0.3, 0.0005, -0.4, 0)
    out <- pecotmr:::.ctwasFilterVariants(
        vids = vids,
        w = w,
        finemapAux = NULL,
        twasWeightCutoff = 0.01,
        csMinCor = 0.8,
        minPipCutoff = 0,
        maxNumVariants = Inf
    )
    # Survivors: v1 (0.5), v3 (0.3), v5 (-0.4) — three with |w| >= 0.01
    expect_setequal(out$vids, c("chr1:100:A:G", "chr1:300:A:G", "chr1:500:A:G"))
})

test_that(".ctwasFilterVariants: maxNumVariants caps by |w| when no PIP", {
    vids <- sprintf("chr1:%d:A:G", 100L * (1:5))
    w <- c(0.1, 0.5, 0.2, 0.4, 0.05)
    out <- pecotmr:::.ctwasFilterVariants(
        vids = vids,
        w = w,
        finemapAux = NULL,
        twasWeightCutoff = 0,
        csMinCor = 0.8,
        minPipCutoff = 0,
        maxNumVariants = 3
    )
    # Top 3 by |w|: v2 (0.5), v4 (0.4), v3 (0.2)
    expect_setequal(out$vids, c("chr1:200:A:G", "chr1:400:A:G", "chr1:300:A:G"))
})

test_that(".ctwasFilterVariants: minPipCutoff rescues high-PIP variants from cap", {
    vids <- sprintf("chr1:%d:A:G", 100L * (1:5))
    w <- c(0.5, 0.4, 0.3, 0.2, 0.1)
    finemapAux <- list(
        pip = setNames(c(0.01, 0.02, 0.8, 0.01, 0.95), vids),
        csMembers = list(),
        csPurity = numeric(0)
    )
    out <- pecotmr:::.ctwasFilterVariants(
        vids = vids,
        w = w,
        finemapAux = finemapAux,
        twasWeightCutoff = 0,
        csMinCor = 0.8,
        minPipCutoff = 0.5,
        maxNumVariants = 2
    )
    # Must-keep (PIP > 0.5): v3, v5. Cap is 2 → both kept.
    expect_setequal(out$vids, c("chr1:300:A:G", "chr1:500:A:G"))
})

test_that(".ctwasFilterVariants: csMinCor rescues high-purity CS members from cap", {
    vids <- sprintf("chr1:%d:A:G", 100L * (1:6))
    w <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
    finemapAux <- list(
        pip = setNames(rep(0, length(vids)), vids),
        csMembers = list(
            c("chr1:300:A:G", "chr1:600:A:G"),
            c("chr1:200:A:G", "chr1:400:A:G")
        ),
        csPurity = c(0.9, 0.5)
    ) # CS 1 (v3, v6) is high-purity
    out <- pecotmr:::.ctwasFilterVariants(
        vids = vids,
        w = w,
        finemapAux = finemapAux,
        twasWeightCutoff = 0,
        csMinCor = 0.8,
        minPipCutoff = 0,
        maxNumVariants = 3
    )
    # Must-keep from high-purity CS: v3, v6. Remaining slot filled by
    # next-highest |w| that isn't must-keep: v1 (0.5).
    expect_setequal(out$vids, c("chr1:300:A:G", "chr1:600:A:G", "chr1:100:A:G"))
})

test_that(".ctwasFilterVariants: returns NULL when no variants survive", {
    vids <- sprintf("chr1:%d:A:G", 100L * (1:3))
    w <- c(0.001, 0.0005, 0.002)
    out <- pecotmr:::.ctwasFilterVariants(
        vids = vids,
        w = w,
        finemapAux = NULL,
        twasWeightCutoff = 0.5,
        csMinCor = 0.8,
        minPipCutoff = 0,
        maxNumVariants = Inf
    )
    expect_null(out)
})

test_that(".ctwasBuildWeights: the variant cap and weight cutoff apply", {
    data(qtlDatasetExample)
    gh <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    vids <- .ctp_panelSnpInfo(gh)$SNP[1:5]
    ldPanel <- .ctp_fullPanel(gh)
    # Each case: the weights, the filter setting, and the variants kept.
    cases <- list(
        # Top 3 by |w| from c(0.1, -0.2, 0.05, 0.3, 0.15): 0.3, -0.2, 0.15
        maxNumVariants = list(
            c(0.1, -0.2, 0.05, 0.3, 0.15),
            list(maxNumVariants = 3L),
            c(4L, 2L, 5L)
        ),
        # v1 (0.005) and v3 (0.001) fall below the 0.01 cutoff
        twasWeightCutoff = list(
            c(0.005, 0.2, 0.001, 0.3, 0.1),
            list(twasWeightCutoff = 0.01),
            c(2L, 4L, 5L)
        )
    )
    for (nm in names(cases)) {
        cs <- cases[[nm]]
        tw <- TwasWeights(
            studyName = "study1",
            context = "brain",
            trait = "ENSG_example",
            method = "susie",
            entry = list(TwasWeightsRow(variantIds = vids, weights = cs[[1L]])),
            ldReference = gh
        )
        wl <- exec(pecotmr:::.ctwasBuildWeights, tw, ldPanel, !!!cs[[2L]])
        expect_equal(wl[[1L]]$n_wgt, 3L, info = nm)
        expect_equal(nrow(wl[[1L]]$wgt), 3L, info = nm)
        expect_equal(
            sort(rownames(wl[[1L]]$wgt)),
            sort(vids[cs[[3L]]]),
            info = nm
        )
    }
})

# ===========================================================================
# mergeCtwasBoundaryRegions (step 4: boundary-gene region merging)
# ===========================================================================

test_that("mergeCtwasBoundaryRegions: no first-pass finemap_res returns unchanged", {
    skip_if_not_installed("ctwas")
    fmr <- list(finemap_res = NULL, region_data = "rd")
    expect_message(
        expect_identical(mergeCtwasBoundaryRegions(fmr), fmr),
        "no first-pass finemap"
    )
    fmr0 <- list(finemap_res = data.frame()[0, ], region_data = "rd")
    expect_identical(suppressMessages(mergeCtwasBoundaryRegions(fmr0)), fmr0)
})

test_that("mergeCtwasBoundaryRegions: LD path forwards carried state + splices updated_*", {
    skip_if_not_installed("ctwas")
    captured <- NULL
    local_mocked_bindings(
        postprocess_region_merging = function(...) {
            captured <<- list(...)
            list(
                updated_finemap_res = data.frame(id = "g", susie_pip = 0.9),
                updated_susie_alpha_res = "ua_new",
                updated_region_data = "rd_new",
                updated_region_info = "ri_new",
                updated_LD_map = "ld_new",
                updated_snp_map = "sm_new",
                selected_boundary_genes = data.frame(id = "g")
            )
        },
        .package = "ctwas"
    )
    fmr <- list(
        finemap_res = data.frame(id = "g", type = "gene", susie_pip = 0.6),
        susie_alpha_res = "ua0",
        region_data = "rd",
        region_info = "ri",
        z_snp = "zs",
        z_gene = data.frame(id = "g"),
        weights = "w",
        snp_map = "sm",
        LD_map = "ld",
        LD_loader_fun = function() NULL,
        snpinfo_loader_fun = function() NULL,
        param = list(group_prior = 0.1, group_prior_var = 5)
    )
    out <- mergeCtwasBoundaryRegions(fmr, pipThresh = 0.5, maxSNP = 100)
    # dispatched to the LD path carrying the loaders + first-pass state
    expect_true(all(
        c("LD_map", "LD_loader_fun", "snpinfo_loader_fun") %in% names(captured)
    ))
    expect_true(is.function(captured$LD_loader_fun))
    expect_equal(captured$pip_thresh, 0.5)
    expect_equal(captured$maxSNP, 100)
    expect_identical(captured$region_data, "rd")
    expect_equal(captured$group_prior, 0.1)
    # updated_* spliced back into the result
    expect_equal(out$finemap_res$susie_pip, 0.9)
    expect_identical(out$region_data, "rd_new")
    expect_identical(out$region_info, "ri_new")
    expect_identical(out$LD_map, "ld_new")
    expect_identical(out$snp_map, "sm_new")
    expect_identical(out$susie_alpha_res, "ua_new")
    expect_identical(out$merge_res$selected_boundary_genes$id, "g")
})

# ===========================================================================
# assembleCtwasInputs: remaining input-validation branches
# ---------------------------------------------------------------------------
# These fire before any heavy panel work; they go through the ctwas-gated
# entry point, so skip when ctwas is absent (otherwise the requireNamespace
# guard would surface a different error).
# ===========================================================================

test_that("assembleCtwasInputs: an unnamed twasWeights list is a FLAT source (needs region)", {
    skip_if_not_installed("ctwas")
    ss <- .ctp_makeGwasSumstats()
    tw <- .ctp_makeTwasWeights() # no region provenance
    # An unnamed list is now treated as a flat weight source to place by region;
    # without region provenance that placement can't happen.
    expect_error(
        assembleCtwasInputs(
            gwasSumStats = ss,
            twasWeights = list(tw),
            method = "susie"
        ), # unnamed -> flat
        "no `traitPos` provenance"
    )
})

test_that("assembleCtwasInputs: rejects a non-GwasSumStats gwas input", {
    skip_if_not_installed("ctwas")
    tw <- .ctp_makeTwasWeights()
    expect_error(
        assembleCtwasInputs(
            gwasSumStats = "not a GwasSumStats",
            twasWeights = list(block1 = tw),
            method = "susie"
        ),
        "must be a GwasSumStats whose elements are LD blocks"
    )
})

test_that("assembleCtwasInputs: rejects a weight source that is neither TwasWeights nor QtlFineMappingResult", {
    skip_if_not_installed("ctwas")
    ss <- .ctp_makeGwasSumstats()
    expect_error(
        assembleCtwasInputs(
            gwasSumStats = ss,
            twasWeights = list(block1 = "not a weight source"),
            method = "susie"
        ),
        "must be a TwasWeights or QtlFineMappingResult"
    )
})

test_that("assembleCtwasInputs: the rescue source is a QtlFineMappingResult", {
    skip_if_not_installed("ctwas")
    ss <- .ctp_makeGwasSumstats()
    tw <- .ctp_makeTwasWeights()
    expect_error(
        assembleCtwasInputs(
            gwasSumStats = ss,
            twasWeights = list(block1 = tw),
            fineMappingResult = "not a fine-mapping result",
            method = "susie"
        ),
        "must be a QtlFineMappingResult"
    )
    # A GWAS result is keyed by a study and method no gene shares, so it can
    # only ever have matched a gene by a coincidence of names.
    tl <- data.frame(variant_id = "chr1:100:A:G", pip = 0.5)
    gwasFmr <- GwasFineMappingResult(
        studyName = "Q1",
        method = "susie",
        entry = list(FineMappingRow(tl$variant_id, NULL, tl))
    )
    expect_error(
        assembleCtwasInputs(
            gwasSumStats = ss,
            twasWeights = list(block1 = tw),
            fineMappingResult = gwasFmr,
            method = "susie"
        ),
        "the genes' own fine-mapping"
    )
})

test_that("assembleCtwasInputs: skips a block whose TwasWeights lacks the resolved method", {
    skip_if_not_installed("ctwas")
    ss <- .ctp_makeGwasSumstats()
    ids5 <- map_chr(1:5, .ctp_snpId)
    mkTw <- function(m) {
        TwasWeights(
            studyName = "Q1",
            context = "c1",
            trait = "t1",
            method = m,
            entry = list(TwasWeightsRow(
                variantIds = ids5,
                weights = c(0.1, 0.2, 0.3, 0.4, 0.5)
            )),
            ldReference = .ctp_makeHandle()
        )
    }
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    # Resolved method is "ensemble" (present in block1); block2 carries only
    # "mrash", so .ctwasFilterMethod returns NULL and the block is skipped.
    inputs <- assembleCtwasInputs(
        gwasSumStats = ss,
        twasWeights = list(block1 = mkTw("ensemble"), block2 = mkTw("mrash"))
    )
    expect_length(inputs$weights, 1L)
    expect_true(all(grepl("^block1\\|", names(inputs$weights))))
    expect_false(any(grepl("^block2\\|", names(inputs$weights))))
})

# ===========================================================================
# .ctwasFilterMethod edge branch
# ===========================================================================

test_that(".ctwasFilterMethod: returns NULL when no row matches the method", {
    tw <- .ctp_makeTwasWeights() # method == "susie"
    expect_null(pecotmr:::.ctwasFilterMethod(tw, "mrash"))
})

# ===========================================================================
# .ctwasHarmonizeWeights early NULL returns
# ===========================================================================

test_that(".ctwasHarmonizeWeights: returns NULL when there is nothing to parse", {
    panel <- .ctp_makeAllelePanel()
    refVariants <- data.frame(
        chrom = panel$snpInfo$chrom,
        pos = panel$snpInfo$pos,
        A2 = panel$snpInfo$ref,
        A1 = panel$snpInfo$alt,
        variant_id = panel$snpInfo$id,
        stringsAsFactors = FALSE
    )
    expect_null(pecotmr:::.ctwasHarmonizeWeights(
        origVids = character(0),
        origW = numeric(0),
        refVariants = refVariants
    ))
})

test_that(".ctwasHarmonizeWeights: returns NULL when harmonizeAlleles fails", {
    local_mocked_bindings(
        harmonizeAlleles = function(...) NULL,
        .package = "pecotmr"
    )
    panel <- .ctp_makeAllelePanel()
    refVariants <- data.frame(
        chrom = panel$snpInfo$chrom,
        pos = panel$snpInfo$pos,
        A2 = panel$snpInfo$ref,
        A1 = panel$snpInfo$alt,
        variant_id = panel$snpInfo$id,
        stringsAsFactors = FALSE
    )
    expect_null(pecotmr:::.ctwasHarmonizeWeights(
        origVids = "1:100:C:T",
        origW = 0.5,
        refVariants = refVariants
    ))
})

# ===========================================================================
# .ctwasBuildWeights: per-gene skip branches + variance / renorm paths
# ===========================================================================

test_that(".ctwasBuildWeights: skips a gene when harmonization yields nothing", {
    local_mocked_bindings(
        .ctwasHarmonizeWeights = function(...) NULL,
        .package = "pecotmr"
    )
    expect_length(
        pecotmr:::.ctwasBuildWeights(
            .ctp_makeTwasWeights(),
            .ctp_makeLdPanel()
        ),
        0L
    )
})

test_that(".ctwasBuildWeights: skips a gene when no variant survives gwasSnpIds intersect", {
    expect_length(
        pecotmr:::.ctwasBuildWeights(
            .ctp_makeTwasWeights(),
            .ctp_makeLdPanel(),
            gwasSnpIds = .ctp_snpId(6)
        ), # gene covers ids 1..5 only
        0L
    )
})

test_that(".ctwasBuildWeights: SuSiE renormalization fires when variants are dropped", {
    panel <- .ctp_makeLdPanel()
    ids4 <- map_chr(1:4, .ctp_snpId)
    bogus <- "chr1:99900:G:A" # absent from the 6-SNP panel
    # Fit dims line up with the 5 original variants. Two single effects,
    # each concentrated on one of the first two variants; lbfToAlpha
    # softmaxes to ~identity rows, so the renormalized weight over the 4
    # kept columns mirrors mu: w[v1]=mu[1,1]=1, w[v2]=mu[2,2]=7, rest ~0.
    fits <- list(
        alpha = lbfToAlpha(rbind(
            c(100, -100, -100, -100, -100),
            c(-100, 100, -100, -100, -100)
        )),
        mu = rbind(c(1, 2, 3, 4, 5), c(6, 7, 8, 9, 10)),
        X_column_scale_factors = rep(1, 5L)
    )
    ent <- TwasWeightsRow(
        variantIds = c(ids4, bogus),
        weights = c(0.1, 0.2, 0.3, 0.4, 0.5),
        methodFits = fits
    )
    tw <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(ent),
        ldReference = .ctp_makeHandle()
    )
    wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
    expect_equal(wl[[1L]]$n_wgt, 4L)
    expect_equal(
        unname(wl[[1L]]$wgt[ids4, 1L]),
        c(1, 7, 0, 0),
        tolerance = 1e-6
    )
})

test_that(".ctwasBuildWeights: errors when the LD panel lacks a kept variant's variance", {
    panel <- .ctp_makeLdPanel()
    panel$variance <- panel$variance[1:4] # drop variance for ids 5, 6
    expect_error(
        pecotmr:::.ctwasBuildWeights(.ctp_makeTwasWeights(), panel),
        "missing genotype variance"
    )
})

test_that(".ctwasBuildWeights: skips a gene when the filter removes every variant", {
    expect_length(
        pecotmr:::.ctwasBuildWeights(
            .ctp_makeTwasWeights(),
            .ctp_makeLdPanel(),
            twasWeightCutoff = 1.0
        ), # |w| max 0.3 -> all dropped
        0L
    )
})

# ===========================================================================
# .ctwasGetFinemapAux — PIP + credible-set membership / purity extraction
# ===========================================================================

test_that(".ctwasGetFinemapAux: parses pip + cs_95 membership + purity", {
    tl <- data.frame(
        variant_id = c(
            "chr1:100:G:A",
            "chr1:200:G:A",
            "chr1:300:G:A",
            "chr1:400:G:A"
        ),
        pip = c(0.3, 0.6, 0.8, 0.02),
        cs_95 = c("susie_1", "susie_1", "susie_2", "susie_0"),
        cs_95_purity = c(0.95, 0.95, 0.7, NA),
        stringsAsFactors = FALSE
    )
    fmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(tl$variant_id, NULL, tl))
    )
    aux <- pecotmr:::.ctwasGetFinemapAux(fmr, "Q1", "c1", "t1")
    expect_equal(aux$pip[["chr1:200:G:A"]], 0.6)
    expect_length(aux$csMembers, 2L)
    expect_setequal(aux$csMembers[[1L]], c("chr1:100:G:A", "chr1:200:G:A"))
    expect_setequal(aux$csMembers[[2L]], "chr1:300:G:A")
    expect_equal(aux$csPurity, c(0.95, 0.7))
})

test_that(".ctwasGetFinemapAux: cs_95 without a purity column yields NA purity", {
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        pip = c(0.1, 0.5, 0.9),
        cs_95 = c("susie_1", "susie_1", "susie_0"),
        stringsAsFactors = FALSE
    )
    fmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(tl$variant_id, NULL, tl))
    )
    aux <- pecotmr:::.ctwasGetFinemapAux(fmr, "Q1", "c1", "t1")
    expect_length(aux$csMembers, 1L)
    expect_setequal(aux$csMembers[[1L]], c("chr1:100:A:G", "chr1:200:A:G"))
    expect_true(all(is.na(aux$csPurity)))
})

test_that(".ctwasGetFinemapAux: no cs_95 column yields empty CS membership", {
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:A:G"),
        pip = c(0.2, 0.8),
        stringsAsFactors = FALSE
    )
    fmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(tl$variant_id, NULL, tl))
    )
    aux <- pecotmr:::.ctwasGetFinemapAux(fmr, "Q1", "c1", "t1")
    expect_length(aux$csMembers, 0L)
    expect_length(aux$csPurity, 0L)
    expect_equal(aux$pip[["chr1:200:A:G"]], 0.8)
})

test_that(".ctwasGetFinemapAux: a cs_95 column naming no set yields none", {
    # A fit that found no credible set still reports its top loci, each one
    # "in no set" (index 0).
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:A:G"),
        pip = c(0.2, 0.3),
        cs_95 = c("susie_0", "susie_0"),
        stringsAsFactors = FALSE
    )
    fmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(tl$variant_id, NULL, tl))
    )
    aux <- pecotmr:::.ctwasGetFinemapAux(fmr, "Q1", "c1", "t1")
    expect_length(aux$csMembers, 0L)
    expect_identical(aux$csPurity, numeric(0))
    expect_equal(aux$pip[["chr1:200:A:G"]], 0.3)
})

test_that(".ctwasRelabelFinemapAux: ids absent from the weights are kept", {
    # Fine-mapping can cover variants the gene's weights do not; those keep
    # their own labels.
    aux <- list(
        pip = c("chr1:900:G:A" = 0.4),
        csMembers = list("chr1:900:G:A"),
        csPurity = 0.9
    )
    expect_identical(
        pecotmr:::.ctwasRelabelFinemapAux(aux, map_chr(1:3, .ctp_snpId)),
        aux
    )
})

test_that(".ctwasFixMolecularId: rows without a gene id are left alone", {
    # SNP-background rows carry variant ids, which have no trait to restore.
    snps <- data.frame(
        id = map_chr(1:2, .ctp_snpId),
        molecular_id = c("x", "y")
    )
    expect_identical(pecotmr:::.ctwasFixMolecularId(snps), snps)
})

test_that(".ctwasGetFinemapAux: NULL input, no-match tuple, and empty topLoci all return NULL", {
    expect_null(pecotmr:::.ctwasGetFinemapAux(NULL, "Q1", "c1", "t1"))
    tl <- data.frame(
        variant_id = "chr1:100:A:G",
        pip = 0.5,
        stringsAsFactors = FALSE
    )
    fmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(tl$variant_id, NULL, tl))
    )
    # No matching (study, context, trait, method) tuple -> NULL.
    expect_null(pecotmr:::.ctwasGetFinemapAux(fmr, "NOPE", "c1", "t1"))
    # Matching tuple but an empty topLoci -> NULL.
    fmrEmpty <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(character(0), NULL, data.frame()))
    )
    expect_null(pecotmr:::.ctwasGetFinemapAux(fmrEmpty, "Q1", "c1", "t1"))
})

# ===========================================================================
# .ctwasFilterVariants / .ctwasSnpInfoForGwasBlock early returns
# ===========================================================================

test_that(".ctwasFilterVariants: returns NULL for an empty variant set", {
    expect_null(pecotmr:::.ctwasFilterVariants(
        vids = character(0),
        w = numeric(0),
        finemapAux = NULL,
        twasWeightCutoff = 0,
        csMinCor = 0.8,
        minPipCutoff = 0,
        maxNumVariants = Inf
    ))
})

test_that(".ctwasSnpInfoForGwasBlock: returns an empty frame when the block has no SNP ids", {
    # GRanges entry with no mcols -> no SNP column -> blockIds is empty.
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100L, width = 1L))
    gss <- GwasSumStats(
        studyName = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .ctp_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
    panelInfo <- data.frame(
        chrom = 1L,
        id = "chr1:100:G:A",
        pos = 100L,
        alt = "A",
        ref = "G",
        stringsAsFactors = FALSE
    )
    out <- pecotmr:::.ctwasSnpInfoForGwasBlock(gss, panelInfo)
    expect_equal(nrow(out), 0L)
    expect_setequal(colnames(out), colnames(panelInfo))
})

# ===========================================================================
# Coverage: internal error / edge branches
# ===========================================================================
test_that(".ctwasCombineWeightSources: type + emptiness guards; FMR list", {
    expect_error(
        pecotmr:::.ctwasCombineWeightSources(42),
        "must be a TwasWeights"
    )
    expect_error(
        pecotmr:::.ctwasCombineWeightSources(list(NULL)),
        "empty weight source"
    )
    fmr <- .ctp_makeFmrWeightSource() # a QtlFineMappingResult
    expect_s4_class(
        pecotmr:::.ctwasCombineWeightSources(list(fmr)),
        "QtlFineMappingResult"
    )
})

test_that(".ctwasGwasStudy / .ctwasGwasStudyFromZSnp: NA when no study present", {
    expect_true(is.na(pecotmr:::.ctwasGwasStudy(list(
        a = data.frame(study = NA_character_)
    ))))
    expect_true(is.na(pecotmr:::.ctwasGwasStudyFromZSnp(NULL)))
    expect_true(is.na(pecotmr:::.ctwasGwasStudyFromZSnp(data.frame(
        study = character(0)
    ))))
})

test_that(".ctwasSubsetById / .ctwasSubsetSnp: NULL on empty / no type column", {
    df <- data.frame(
        id = c("a", "b"),
        type = c("gene", "SNP"),
        stringsAsFactors = FALSE
    )
    expect_null(pecotmr:::.ctwasSubsetById(df, "zzz"))
    expect_null(pecotmr:::.ctwasSubsetById(NULL, "a"))
    expect_null(pecotmr:::.ctwasSubsetSnp(data.frame(id = "a"))) # no type column
    expect_equal(nrow(pecotmr:::.ctwasSubsetSnp(df)), 1L)
})

test_that(".ctwasRowsToResult / .ctwasMethodFromWeights: empty-input guards", {
    expect_error(pecotmr:::.ctwasRowsToResult(list()), "no genes were modeled")
    expect_error(pecotmr:::.ctwasMethodFromWeights(NULL), "carries no weights")
    expect_error(
        pecotmr:::.ctwasMethodFromWeights(setNames(
            list(1, 2),
            c("blk|Q1|c1|gA|susie", "blk|Q1|c1|gB|lasso")
        )),
        "mixes weight methods"
    )
})

test_that(".ctwasRunToRows: empty weights -> no rows; a context mixing studies errors", {
    expect_equal(
        pecotmr:::.ctwasRunToRows(
            list(weights = setNames(list(), character(0))),
            "D1",
            "susie"
        ),
        list()
    )
    run <- list(
        weights = setNames(
            list(1, 2),
            c("blk|Q1|c1|gA|susie", "blk|Q2|c1|gB|susie")
        ),
        finemap_res = NULL,
        susie_alpha_res = NULL,
        param = list(),
        region_info = NULL
    )
    expect_error(
        pecotmr:::.ctwasRunToRows(run, "D1", "susie"),
        "mixes multiple QTL studies"
    )
})

test_that(".ctwasBucketWeights: unplaced genes warn + drop; empty blocks skipped", {
    mkE <- function() {
        TwasWeightsRow(
            variantIds = map_chr(1:5, .ctp_snpId),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
        )
    }
    tw <- TwasWeights(
        studyName = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("gA", "gX"),
        method = c("susie", "susie"),
        entry = list(mkE(), mkE()),
        traitPos = GenomicRanges::GRanges(
            c("chr1", "chr9"),
            IRanges::IRanges(c(100, 100), c(150, 150))
        ),
        ldReference = .ctp_makeHandle()
    )
    grid <- list(
        chr1_1_350 = .ctp_makeGwasSumstats(),
        chr1_350_700 = .ctp_makeGwasSumstats()
    )
    b <- expect_warning(
        pecotmr:::.ctwasBucketWeights(tw, grid),
        "fell in no LD block"
    )
    expect_named(b, "chr1_1_350") # 2nd block empty -> skipped
    expect_equal(as.character(b[["chr1_1_350"]]$trait), "gA")
})

test_that(".ctwasBucketWeights: errors when no gene lands in any block", {
    mkE <- function() {
        TwasWeightsRow(
            variantIds = map_chr(1:5, .ctp_snpId),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
        )
    }
    twOff <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "gX",
        method = "susie",
        entry = list(mkE()),
        traitPos = GenomicRanges::GRanges("chr9", IRanges::IRanges(100, 150)),
        ldReference = .ctp_makeHandle()
    )
    grid <- list(
        chr1_1_350 = .ctp_makeGwasSumstats(),
        chr1_350_700 = .ctp_makeGwasSumstats()
    )
    expect_error(
        suppressWarnings(pecotmr:::.ctwasBucketWeights(twOff, grid)),
        "no gene placed"
    )
})

test_that("assembleCtwasInputs: single-block + NULL twasWeights guards", {
    skip_if_not_installed("ctwas")
    expect_error(
        assembleCtwasInputs(
            .ctp_makeGwasSumstats(blockIds = "block1"),
            .ctp_makeTwasWeights(),
            method = "susie"
        ),
        "at least two LD blocks"
    )
    expect_error(
        assembleCtwasInputs(.ctp_makeGwasSumstats(), NULL),
        "twasWeights` is required"
    )
})

test_that("weight-source / method / gwas-study guards reject degenerate input", {
    expect_error(
        pecotmr:::.ctwasRequireNamedLists(.ctp_makeGwasSumstats(), 42),
        "must be a TwasWeights"
    )
    expect_error(pecotmr:::.ctwasResolveMethod(list()), "no method entries")
    expect_error(
        pecotmr:::.ctwasGwasStudyFromZSnp(data.frame(study = c("D1", "D2"))),
        "multiple GWAS studies"
    )
})

test_that("finemapCtwasRegions: an empty screened-region set returns a NULL finemap", {
    skip_if_not_installed("ctwas")
    screenStub <- list(
        screened_region_data = list(),
        z_gene = NULL,
        region_data = list(),
        boundary_genes = NULL,
        screen_res = NULL,
        param = list(group_prior = c(g = 0.1), group_prior_var = c(g = 5)),
        region_info = data.frame(),
        z_snp = data.frame(),
        weights = list(),
        snp_map = list(),
        LD_map = data.frame(),
        LD_loader_fun = NULL,
        snpinfo_loader_fun = NULL
    )
    out <- finemapCtwasRegions(screenStub)
    expect_null(out$finemap_res)
    expect_null(out$susie_alpha_res)
})

# =============================================================================
# Bundled example payloads: portable LD tokens
#
# `LD_map$LD_file` is both something ctwas asserts exists on disk and the key
# pecotmr dispatches on into the cached LD panels. Those two duties are why
# the token cannot be the backing genotype file: every block on a chromosome
# shares it, so each block was handed the FIRST block's panel. Each region
# now carries its own sentinel, and a payload read back in a later session
# finds those gone -- .ctwasResolveLdPaths() re-mints one per region and
# re-keys the loader cache onto them.
# =============================================================================

test_that("bundled cTWAS payloads carry one resolvable token per region", {
    data(ctwasInputsExample)
    data(ctwasEstExample)
    data(ctwasFinemapExample)

    for (payload in list(
        ctwasInputsExample,
        ctwasEstExample,
        ctwasFinemapExample
    )) {
        tokens <- as.character(payload$LD_map$LD_file)
        # One per region, never shared: the bug was every region reading
        # block 1's panel back out of a single key.
        expect_equal(length(unique(tokens)), nrow(payload$LD_map))
        expect_equal(as.character(payload$LD_map$SNP_file), tokens)
        # Whatever state the sentinels are in on this machine, the resolver
        # leaves every region with a token that exists and its own panel.
        resolved <- pecotmr:::.ctwasResolveLdPaths(payload)
        live <- as.character(resolved$LD_map$LD_file)
        expect_true(all(file.exists(live)))
        expect_equal(length(unique(live)), nrow(payload$LD_map))
        loader <- resolved$LD_loader_fun
        panels <- pecotmr:::.ctwasCachedPanels(resolved)
        # Both loader closures are re-keyed onto the new tokens; rewriting
        # LD_map alone would leave them unable to find their panels.
        for (i in seq_along(live)) {
            expect_equal(
                loader(live[[i]]),
                pecotmr:::.ctwasPanelLd(panels[[i]])
            )
            expect_false(is.null(resolved$snpinfo_loader_fun(live[[i]])))
        }
    }
})

test_that("an ordinary LD path is left alone", {
    # A token still on disk -- the caller's own file, or this session's
    # sentinel -- needs nothing re-minted.
    path <- withr::local_tempfile(fileext = ".ld")
    file.create(path)
    payload <- list(
        LD_map = data.frame(region_id = "r1", LD_file = path, SNP_file = path)
    )
    expect_identical(pecotmr:::.ctwasResolveLdPaths(payload), payload)
    noLdMap <- list(z_snp = NULL)
    expect_equal(pecotmr:::.ctwasResolveLdPaths(noLdMap), noLdMap)
    # Hand-built payloads reach the granular steps with an LD_map stub that is
    # not a table at all; it has to survive untouched rather than error.
    stub <- list(LD_map = "ld", LD_loader_fun = function() NULL)
    expect_equal(pecotmr:::.ctwasResolveLdPaths(stub), stub)
    expect_equal(
        pecotmr:::.ctwasResolveLdPaths(list(LD_map = list(other = 1))),
        list(LD_map = list(other = 1))
    )
})

test_that("finemapCtwasRegions runs from the bundled est payload", {
    skip_if_not_installed("ctwas")
    data(ctwasEstExample)
    # The step that the stale absolute paths broke.
    screened <- suppressMessages(
        screenCtwasRegions(
            ctwasEstExample,
            methodArgs = CtwasOptions(
                screen = CtwasScreenOptions(min_nonSNP_PIP = 0)
            )
        )
    )
    out <- suppressMessages(finemapCtwasRegions(screened))
    expect_gt(nrow(out$finemap_res), 0L)
})


# =============================================================================
# ctwasPipeline shares the LD-sketch validator
#
# ctwasPipeline cannot use `.ldFromRef()` itself: it needs the whole panel
# rather than a matched subset, and returns per-variant variance alongside R.
# It must still fail the same way when the LD reference is missing, rather
# than surfacing an S4 dispatch error from deep inside.
# =============================================================================

test_that("ctwas LD assembly reports a missing LD sketch like everything else", {
    expect_error(
        pecotmr:::.ctwasPanelFor(NULL),
        "carries no ldReference"
    )
    expect_error(
        pecotmr:::.ctwasPanelFor("not a panel"),
        "must be a genotype panel"
    )
})

test_that("the shared validator gives ctwas and .ldFromRef one message", {
    shared <- tryCatch(
        pecotmr:::.ldFromRef(NULL, "chr1:1:A:G", label = "ctwasPipeline"),
        error = conditionMessage
    )
    own <- tryCatch(
        pecotmr:::.ctwasPanelFor(NULL),
        error = conditionMessage
    )
    expect_equal(own, shared)
})


# ===========================================================================
# Degenerate inputs to the ctwas assembly helpers
# ===========================================================================

test_that("gene coordinates are NULL when there are no weights", {
    expect_null(pecotmr:::.ctwasGeneCoords(NULL))
    expect_null(pecotmr:::.ctwasGeneCoords(list()))
})

test_that("z harmonization is a no-op when either side is empty", {
    # Both exits return the input untouched: an empty GWAS has nothing to
    # relabel, and an empty panel has nothing to relabel it to.
    empty <- tibble(
        id = character(0),
        chrom = integer(0),
        pos = integer(0),
        A1 = character(0),
        A2 = character(0),
        z = numeric(0),
        study = character(0)
    )
    expect_equal(
        nrow(pecotmr:::.ctwasHarmonizeZToPanel(empty, "chr1:1:A:G")),
        0L
    )
    one <- tibble(
        id = "chr1:1:A:G",
        chrom = 1L,
        pos = 1L,
        A1 = "G",
        A2 = "A",
        z = 1.5,
        study = "S"
    )
    expect_identical(
        pecotmr:::.ctwasHarmonizeZToPanel(one, character(0)),
        one
    )
    # A panel that shares no variant leaves the ids alone.
    expect_equal(
        pecotmr:::.ctwasHarmonizeZToPanel(one, "chr9:9:C:T")$id,
        "chr1:1:A:G"
    )
})

test_that("alpha renormalization refuses a row with no finite mass", {
    # log(0) across a whole row leaves no maximum to subtract, so the
    # renormalized alpha would be NaN rather than a distribution.
    alpha <- matrix(c(0, 0, 0.5, 0.5), nrow = 2L, byrow = TRUE)
    expect_null(pecotmr:::.ctwasRenormAlpha(alpha, 1:2))
    ok <- matrix(c(0.25, 0.75, 0.5, 0.5), nrow = 2L, byrow = TRUE)
    out <- pecotmr:::.ctwasRenormAlpha(ok, 1:2)
    expect_equal(unname(rowSums(out)), c(1, 1))
})


test_that("the SNP background row is omitted when no SNP rows were fitted", {
    # ctwas tags its non-gene background with type == "SNP"; a run with none
    # contributes no row rather than an empty one.
    expect_null(
        pecotmr:::.ctwasSnpRow("G", "susie", "", NULL, NULL, NULL)
    )
})

test_that("subsetting to SNP rows answers NULL when there are none", {
    expect_null(pecotmr:::.ctwasSubsetSnp(NULL))
    geneOnly <- data.frame(
        id = "a",
        type = "gene",
        stringsAsFactors = FALSE
    )
    expect_null(pecotmr:::.ctwasSubsetSnp(geneOnly))
    mixed <- data.frame(
        id = c("a", "b"),
        type = c("SNP", "gene"),
        stringsAsFactors = FALSE
    )
    expect_equal(nrow(pecotmr:::.ctwasSubsetSnp(mixed)), 1L)
})


# ===========================================================================
# Assembly shortcuts and block placement
# ===========================================================================

test_that("a precomputed z_gene is used rather than recomputed", {
    # compute_gene_z is the expensive step; supplying it must short-circuit.
    expect_equal(
        pecotmr:::.ctwasEnsureZGene(list(z_gene = "PRECOMPUTED"), 1L),
        "PRECOMPUTED"
    )
})

test_that("boundary genes are NULL for a single-region run", {
    # Boundary adjustment only means something across two or more regions.
    expect_null(pecotmr:::.ctwasBoundaryGenes(
        list(region_info = data.frame(a = 1)),
        1L,
        list()
    ))
})

test_that("relabelling ids is a no-op on an empty set and flips when matched", {
    expect_equal(
        pecotmr:::.ctwasRelabelIds(character(0), "chr1:1:A:G"),
        character(0)
    )
    # A swapped spelling is the same variant, relabelled to the panel frame.
    expect_equal(
        pecotmr:::.ctwasRelabelIds("chr1:1:G:A", "chr1:1:A:G"),
        "chr1:1:A:G"
    )
})

test_that("re-keying LD loaders is a no-op when no panels are cached", {
    payload <- list(a = 1)
    expect_identical(
        pecotmr:::.ctwasRekeyLdLoaders(payload, list(x = "y")),
        payload
    )
})

test_that("a variant is placed in the block whose half-open window holds it", {
    # [start, end): an unplaceable or out-of-range anchor is NA rather than
    # being snapped to the nearest block.
    f <- pecotmr:::.ctwasBlockIdForVariant
    expect_equal(f(1L, 50L, "1", "1", 1L, 100L, "B1"), "B1")
    expect_true(is.na(f(1L, 500L, "1", "1", 1L, 100L, "B1")))
    expect_true(is.na(f(1L, NA_integer_, "1", "1", 1L, 100L, "B1")))
    # The window excludes its end.
    expect_true(is.na(f(1L, 100L, "1", "1", 1L, 100L, "B1")))
})


# ===========================================================================
# Row coordinates, trait positions and the loader panel cache
# ===========================================================================

test_that("row coordinates prefer the gene cache, then the variant id", {
    # Genes are keyed by name and have no parseable coordinates; variants
    # carry theirs in the id. An id that is neither is NULL rather than a
    # fabricated position.
    f <- pecotmr:::.ctwasRowCoord
    cached <- f(1L, "G1", list(G1 = list(chrom = "chr9", start = 5L, end = 9L)))
    expect_equal(cached$chrom, "chr9")
    parsed <- f(1L, "chr1:100:A:G", NULL)
    expect_equal(parsed$chrom, "chr1")
    expect_equal(parsed$start, 100L)
    # A single-base variant spans one position.
    expect_equal(parsed$end, 100L)
    expect_null(f(1L, "rsNOTPARSEABLE", NULL))
})

test_that("ranged results drop only the rows they cannot place, saying so", {
    df <- data.frame(
        id = c("chr1:100:A:G", "geneX", "geneY", "rsNOTPARSEABLE"),
        susie_pip = c(0.1, 0.2, 0.3, 0.4)
    )
    gc <- list(geneX = list(chrom = "chr1", start = 500L, end = 900L))
    expect_warning(
        out <- pecotmr:::.ctwasRangePayload(df, gc),
        "2 result row\\(s\\) could not be placed on the genome"
    )
    expect_s4_class(out, "GRanges")
    expect_equal(out$id, c("chr1:100:A:G", "geneX"))
    expect_equal(start(out), c(100L, 500L))
    expect_equal(end(out), c(100L, 900L))
})

test_that("a table with no placeable row stays whole, without a warning", {
    df <- data.frame(id = c("geneY", "rsNOTPARSEABLE"), susie_pip = c(0.3, 0.4))
    expect_no_warning(out <- pecotmr:::.ctwasRangePayload(df, NULL))
    expect_identical(out, df)
    expect_null(pecotmr:::.ctwasRangePayload(NULL, NULL))
})

test_that("the unplaced-row warning names at most five ids", {
    cnd <- tryCatch(
        pecotmr:::.ctwasWarnUnplaced(sprintf("gene%d", 1:7)),
        warning = identity
    )
    expect_match(conditionMessage(cnd), "7 result row\\(s\\)")
    expect_equal(stringr::str_count(conditionMessage(cnd), "gene[0-9]"), 5L)
    expect_null(pecotmr:::.ctwasWarnUnplaced(character(0)))
})

test_that("trait positions are NULL when the weights do not carry them", {
    data(twasWeightsExample)
    expect_false(
        is_in("traitPos", pecotmr:::.tupleColumnNames(twasWeightsExample))
    )
    expect_null(pecotmr:::.ctwasTraitPosAt(twasWeightsExample, 1L))
})

test_that("the panel cache is NULL for loaders this package did not build", {
    # Re-keying only makes sense over our own closure; a foreign or
    # non-function loader is left alone rather than reached into.
    expect_null(pecotmr:::.ctwasCachedPanels(list(LD_loader_fun = "nope")))
    expect_null(
        pecotmr:::.ctwasCachedPanels(list(LD_loader_fun = function(x) x))
    )
})


# ===========================================================================
# Input validation
# ===========================================================================

test_that("block ids must be present, non-empty and unique", {
    # blockId becomes the region id, so a missing or empty one leaves an
    # element unkeyable and a repeated one would have two elements silently
    # overwrite each other.
    gs <- .ctp_makeGwasSumstats(blockIds = c("block1", "block2"))
    expect_silent(pecotmr:::.ctwasValidateGwasEntries(gs))

    missing <- gs
    S4Vectors::mcols(missing)$blockId <- c("b", NA_character_)
    expect_error(
        pecotmr:::.ctwasValidateGwasEntries(missing),
        "empty or missing `blockId`"
    )

    blank <- gs
    S4Vectors::mcols(blank)$blockId <- c("", "b2")
    expect_error(
        pecotmr:::.ctwasValidateGwasEntries(blank),
        "empty or missing `blockId`"
    )

    dup <- gs
    S4Vectors::mcols(dup)$blockId <- c("dup", "dup")
    expect_error(
        pecotmr:::.ctwasValidateGwasEntries(dup),
        "block ids must be unique; repeated: dup"
    )
})

test_that("each weight entry must be a TwasWeights or QtlFineMappingResult", {
    # The per-gene weight source is one of those two; anything else would
    # fail later inside the assembly with a far less specific message.
    expect_error(
        pecotmr:::.ctwasValidateWeightEntries(list(r1 = "not a TwasWeights")),
        "twasWeights\\[\\['r1'\\]\\] must be a TwasWeights or QtlFineMappingResult"
    )
})

# ---------------------------------------------------------------------------
# Input guards and coordinate fallbacks: each rejects a shape the pipeline
# cannot proceed on, or substitutes a placeholder for an unplaced feature.
# ---------------------------------------------------------------------------

test_that(".ctwasValidateGwasList requires a blockId column", {
    skip_if_not_installed("ctwas")
    gss <- .ctp_makeGwasSumstats()
    # A GwasSumStats built by a current constructor always carries blockId;
    # dropping it stands in for an object deserialized from an older one.
    S4Vectors::mcols(gss)$blockId <- NULL
    expect_false(is_in("blockId", colnames(gss)))
    expect_error(
        pecotmr:::.ctwasValidateGwasList(gss),
        "has no `blockId` column"
    )
})

test_that(".ctwasResolveAndValidateWeights requires a weight source", {
    expect_error(
        pecotmr:::.ctwasResolveAndValidateWeights(NULL, NULL),
        "`twasWeights` is required"
    )
    expect_error(
        pecotmr:::.ctwasResolveAndValidateWeights(),
        "`twasWeights` is required"
    )
})

test_that(".ctwasPrefitRegionFilter errors when nothing survives", {
    regionData <- list(
        r1 = list(gid = character(0), sid = "s1"),
        r2 = list(gid = character(0), sid = character(0))
    )
    # Default minGene = 1 drops both regions: neither carries a gene.
    expect_error(
        pecotmr:::.ctwasPrefitRegionFilter(regionData, list()),
        "No regions selected!"
    )
    kept <- pecotmr:::.ctwasPrefitRegionFilter(
        list(r1 = list(gid = "g1", sid = c("s1", "s2"))),
        list()
    )
    expect_equal(names(kept), "r1")
})

test_that(".ctwasPrefitRegionFilter honours min_var = 0 and min_gene = 0", {
    # The est_param options reach the prefit skip too: min_var = 0 keeps a
    # region of any size, and min_gene = 0 one without a gene.
    regionData <- list(
        r1 = list(gid = "g1", sid = character(0)),
        r2 = list(gid = character(0), sid = "s1")
    )
    expect_named(
        pecotmr:::.ctwasPrefitRegionFilter(
            regionData,
            list(min_var = 0, min_gene = 0)
        ),
        c("r1", "r2")
    )
    expect_named(
        pecotmr:::.ctwasPrefitRegionFilter(regionData, list(min_var = 0)),
        "r1"
    )
})

test_that(".ctwasCoordField falls back for an unplaced feature", {
    expect_equal(pecotmr:::.ctwasCoordField(NULL, "chrom"), "chrUnplaced")
    expect_equal(pecotmr:::.ctwasCoordField2(NULL, "start"), 1L)
    co <- list(chrom = "chr7", start = 42L)
    expect_equal(pecotmr:::.ctwasCoordField(co, "chrom"), "chr7")
    expect_equal(pecotmr:::.ctwasCoordField2(co, "start"), 42L)
})

test_that(".ctwasBuildSingleRegionInfo reports an empty block accurately", {
    emptyGr <- GenomicRanges::GRanges()
    S4Vectors::mcols(emptyGr)$variant_id <- character(0)
    S4Vectors::mcols(emptyGr)$Z <- numeric(0)
    S4Vectors::mcols(emptyGr)$N <- integer(0)
    gss <- GwasSumStats(
        studyName = "G1",
        entry = list(emptyGr),
        genome = "hg19",
        blockId = "b1"
    )
    # Emptiness is checked before the window -- a block empty from the start
    # has only a placeholder for one.
    expect_error(
        pecotmr:::.ctwasBuildSingleRegionInfo("b1", gss),
        "has no variants for cTWAS to fit"
    )
})

# ---------------------------------------------------------------------------
# Per-gene weight assembly: each skip below returns NULL so the gene drops out
# of the run, rather than contributing a mis-sliced weight vector.
# ---------------------------------------------------------------------------

.ctp_noLdGwas <- function() {
    mkGr <- function(chr, pos) {
        g <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos, width = 1))
        S4Vectors::mcols(g)$variant_id <- str_c(chr, ":", pos, ":A:G")
        S4Vectors::mcols(g)$Z <- rep(1, length(pos))
        S4Vectors::mcols(g)$N <- rep(100L, length(pos))
        g
    }
    GwasSumStats(
        studyName = c("G1", "G1"),
        entry = list(mkGr("chr1", c(100L, 200L)), mkGr("chr1", c(300L, 400L))),
        genome = "hg19",
        blockId = c("b1", "b2"),
        qcOptions = list(step1 = "ok")
    )
}

test_that(".ctwasFirstPass requires an LD reference on every region", {
    byBlock <- pecotmr:::.ctwasGwasByBlock(.ctp_noLdGwas())
    expect_identical(ldReference(byBlock[[1L]]), list())
    expect_error(
        pecotmr:::.ctwasFirstPass("b1", byBlock, list(b1 = NULL)),
        "region 'b1' carries no ldReference"
    )
})

test_that("ctwasPipeline aborts when no genes could be modeled", {
    local_mocked_bindings(
        .ctwasRequireNamedLists = function(...) invisible(NULL),
        .ctwasResolveMethod = function(...) "susie",
        .ctwasGwasStudy = function(...) "G1",
        .ctwasRunMethod = function(...) list(),
        .package = "pecotmr"
    )
    expect_error(
        ctwasPipeline(gwasSumStats = .ctp_noLdGwas(), twasWeights = list()),
        "no genes were modeled"
    )
})

test_that(".ctwasAlignGeneWeights returns NULL when the row has no variants", {
    row <- TwasWeightsRow(variantIds = character(0), weights = numeric(0))
    expect_null(pecotmr:::.ctwasAlignGeneWeights(row, NULL, NULL))
})

test_that(".ctwasTraitPosAt is bounded by the traitPos column length", {
    row <- FineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:C:T"),
        susieFit = list(pip = c(0.8, 0.2)),
        topLoci = data.frame(
            variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
            pip = c(0.8, 0.2),
            stringsAsFactors = FALSE
        )
    )
    res <- GwasFineMappingResult(
        studyName = "G1",
        method = "susie",
        entry = list(row),
        traitPos = GenomicRanges::GRanges("chr1", IRanges::IRanges(500, 900))
    )
    expect_s4_class(pecotmr:::.ctwasTraitPosAt(res, 1L), "GRanges")
    # One row, so row 2 has no anchor to report.
    expect_null(pecotmr:::.ctwasTraitPosAt(res, 2L))
})

test_that(".ctwasRenormalizeSusieWeights skips fits it cannot slice", {
    alpha <- matrix(c(0.5, 0.5, 0.3, 0.7), nrow = 2L, byrow = TRUE)
    mu <- matrix(1, nrow = 2L, ncol = 2L)
    good <- list(alpha = alpha, mu = mu, X_column_scale_factors = c(1, 1))
    origVids <- c("v1", "v2")
    origW <- c(0.5, 0.5)
    keptIdx <- c(1L, 2L)
    harmonizedW <- c(0.5, 0.5)
    run <- function(fits) {
        pecotmr:::.ctwasRenormalizeSusieWeights(
            fits,
            origVids,
            origW,
            keptIdx,
            harmonizedW
        )
    }
    # Any missing susie field: nothing to renormalize from.
    expect_null(run(list(mu = mu, X_column_scale_factors = c(1, 1))))
    expect_null(run(list(alpha = alpha, X_column_scale_factors = c(1, 1))))
    expect_null(run(list(alpha = alpha, mu = mu)))
    # susieInf / susieAsh carry an infinitesimal term that alpha and mu alone
    # cannot reproduce, so those fits take the plain subset path instead.
    expect_null(run(c(good, list(theta = 1))))
    expect_null(run(c(good, list(omega_weights = 1))))
    # Fit-vs-entry dimension mismatch (e.g. a null_weight fit's extra column).
    expect_null(run(list(
        alpha = matrix(0.25, nrow = 2L, ncol = 3L),
        mu = matrix(1, nrow = 2L, ncol = 3L),
        X_column_scale_factors = c(1, 1, 1)
    )))
    # An all-zero alpha has no finite row maximum in log space.
    expect_null(run(list(
        alpha = matrix(0, nrow = 2L, ncol = 2L),
        mu = mu,
        X_column_scale_factors = c(1, 1)
    )))
    expect_equal(run(good), c(0.8, 1.2))
})

test_that(".ctwasResolveAndValidateWeights requires region_id names", {
    # A partially-named list is as unusable as an unnamed one: the blank key
    # cannot address a region.
    buckets <- list(
        unnamed = list("a", "b"),
        blank = stats::setNames(list("a", "b"), c("r1", ""))
    )
    current <- NULL
    local_mocked_bindings(
        .ctwasResolveWeightBuckets = function(twasWeights, gwasSumStats) {
            current
        },
        .package = "pecotmr"
    )
    for (nm in names(buckets)) {
        current <- buckets[[nm]]
        expect_error(
            pecotmr:::.ctwasResolveAndValidateWeights(list(x = 1), NULL),
            "must resolve to a named list keyed by region_id",
            info = nm
        )
    }
})

test_that("asCtwasResult: keepSnps must be a flag", {
    expect_error(asCtwasResult(NULL, keepSnps = NA), "keepSnps.*May not be NA")
    expect_error(
        asCtwasResult(NULL, keepSnps = 1L),
        "keepSnps.*logical flag"
    )
})

test_that("the ctwas concat helpers answer an empty vector of their own type", {
    expect_identical(pecotmr:::.ctwasConcat(list()), list())
    expect_identical(pecotmr:::.ctwasConcatChr(list()), character(0))
})

test_that("ctwas entry points refuse a bare list for methodArgs", {
    expect_error(
        screenCtwasRegions(list(), methodArgs = list(min_nonSNP_PIP = 0)),
        "must be built with CtwasOptions\\(\\)"
    )
})


test_that("ctwas reads variance and LD from an LdMatrix block", {
    ids <- c("chr22:100:A:G", "chr22:200:C:T", "chr22:300:G:A")
    R <- matrix(c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1), 3L, 3L)
    vars <- data.frame(
        chrom = "chr22",
        pos = c(100L, 200L, 300L),
        SNP = ids,
        A1 = c("G", "T", "A"),
        A2 = c("A", "C", "G")
    )
    bare <- LdMatrix(R, variants = vars, nSamples = 100)
    # No AF: no genotype variance, which only non-standardized weights need.
    panel <- .ctwasPanelFor(bare)
    expect_true(all(is.na(panel$variance)))
    expect_equal(unname(.ctwasPanelLd(panel)), R)
    withAf <- LdMatrix(
        R,
        variants = transform(vars, AF = c(0.1, 0.3, 0.5)),
        nSamples = 100
    )
    expect_equal(
        unname(.ctwasPanelFor(withAf)$variance),
        2 * c(0.1, 0.3, 0.5) * (1 - c(0.1, 0.3, 0.5))
    )
    # An LdMatrix has no file payload to check.
    expect_null(.ctwasLdPanelKey(bare))
})

test_that(".ctwasRegionGwasLd refuses a region whose rows use two references", {
    # A region's LD is one matrix, so its GWAS rows must share a reference.
    gwas <- .ctp_makeGwasSumstats(blockIds = c("blockA", "blockA"))
    gwas <- pecotmr:::.withRowLdRefs(
        gwas,
        list(.ctp_makeHandle(snp_n = 12L), .ctp_makeHandle(snp_n = 6L))
    )
    expect_error(
        pecotmr:::.ctwasRegionGwasLd(1L, list(gwas), list()),
        "2 different LD references; a region's LD is one matrix"
    )
})

# ===========================================================================
# The rescues read the fine-mapping method they are told to
# ===========================================================================

# A gene's fine-mapping under two methods that disagree about variant 3: susie
# puts its PIP at 0.9, susieInf at 0.01.
.ctp_twoMethodRescue <- function(vids) {
    fit <- function(pip3, rest) {
        pip <- rep(rest, 5L)
        pip[[3L]] <- pip3
        FineMappingRow(
            variantIds = vids,
            susieFit = list(),
            topLoci = data.frame(variant_id = vids, pip = pip)
        )
    }
    QtlFineMappingResult(
        studyName = rep("study1", 2L),
        context = rep("brain", 2L),
        trait = rep("ENSG_example", 2L),
        method = c("susie", "susieInf"),
        entry = list(fit(0.9, 0.05), fit(0.01, 0.2))
    )
}

test_that(".ctwasRescueSource narrows the result to the named method", {
    data(qtlDatasetExample)
    vids <- .ctp_panelSnpInfo(pecotmr:::.qtlGenotypeSe(qtlDatasetExample))$SNP
    fmr <- .ctp_twoMethodRescue(vids[1:5])
    expect_null(pecotmr:::.ctwasRescueSource(NULL, NULL))
    narrowed <- pecotmr:::.ctwasRescueSource(fmr, "susieInf")
    expect_identical(as.character(narrowed$method), "susieInf")
    expect_error(
        pecotmr:::.ctwasRescueSource(fmr, NULL),
        "must name the method of it .*available: susie, susieInf"
    )
    expect_error(
        pecotmr:::.ctwasRescueSource(fmr, "ensemble"),
        "'ensemble' is not a method of `fineMappingResult`"
    )
    expect_error(
        pecotmr:::.ctwasRescueSource(NULL, "susie"),
        "names a method of `fineMappingResult`, which was not given"
    )
})

test_that("an ensemble gene is rescued by the fine-mapping method named", {
    data(qtlDatasetExample)
    gh <- pecotmr:::.qtlGenotypeSe(qtlDatasetExample)
    vids <- .ctp_panelSnpInfo(gh)$SNP[1:5]
    # Variant 3 carries the smallest weight, so the cap drops it unless a
    # rescue keeps it.
    tw <- TwasWeights(
        studyName = "study1",
        context = "brain",
        trait = "ENSG_example",
        method = "ensemble",
        entry = list(TwasWeightsRow(
            variantIds = vids,
            weights = c(0.1, -0.2, 0.05, 0.3, 0.15)
        )),
        ldReference = gh
    )
    fmr <- .ctp_twoMethodRescue(vids)
    kept <- function(fineMappingMethod) {
        wl <- pecotmr:::.ctwasBuildWeights(
            tw,
            .ctp_fullPanel(gh),
            fineMappingResult = pecotmr:::.ctwasRescueSource(
                fmr,
                fineMappingMethod
            ),
            minPipCutoff = 0.5,
            maxNumVariants = 3L
        )
        rownames(wl[[1L]]$wgt)
    }
    expect_true(vids[[3L]] %in% kept("susie"))
    expect_false(vids[[3L]] %in% kept("susieInf"))
})

# ===========================================================================
# Narrowed rows: the SuSiE fits that renormalize weights
# ===========================================================================

# One gene's SuSiE weights on variants `variantIdx`, fitted on all five and
# narrowed to the first four; `fits` is the stored fit.
.ctp_narrowedSusieWeights <- function(
    variantIdx = 1:5,
    method = "susie",
    fits = list(
        alpha = matrix(0.2, 1L, 5L),
        mu = matrix(0.1, 1L, 5L),
        X_column_scale_factors = rep(1, 5L)
    )
) {
    tw <- TwasWeights(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = method,
        entry = list(TwasWeightsRow(
            variantIds = map_chr(variantIdx, .ctp_snpId),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0),
            methodFits = fits
        )),
        ldReference = .ctp_makeHandle()
    )
    endoapply(tw, utils::head, 4L)
}

test_that("the narrowed rows checked are SuSiE-fitted rows of the method", {
    narrowed <- .ctp_narrowedSusieWeights()
    cnd <- tryCatch(
        pecotmr:::.ctwasCheckFitCoverage(
            list(block1 = narrowed, block2 = .ctp_narrowedSusieWeights(7:11)),
            "susie",
            "f()"
        ),
        warning = identity
    )
    expect_s3_class(cnd, "pecotmrFitCoverage")
    expect_match(conditionMessage(cnd), "f\\(\\): 2 row\\(s\\)")
    # Another method's rows are not read.
    expect_no_warning(
        pecotmr:::.ctwasCheckFitCoverage(narrowed, "lasso", "f()")
    )
    # Nor is a fit that is not SuSiE's, which renormalizes nothing.
    lasso <- .ctp_narrowedSusieWeights(
        method = "lasso",
        fits = list(beta = rep(0.1, 5L))
    )
    expect_null(pecotmr:::.ctwasSusieFitRows(lasso, "lasso"))
    expect_no_warning(pecotmr:::.ctwasCheckFitCoverage(lasso, "lasso", "f()"))
    # A fine-mapping weight source gives its topLoci effects, not its fit.
    vids <- map_chr(1:5, .ctp_snpId)
    tl <- data.frame(variant_id = vids, pip = rep(0.2, 5L))
    fmr <- QtlFineMappingResult(
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(vids, list(pip = tl$pip), tl))
    )
    expect_no_warning(pecotmr:::.ctwasCheckFitCoverage(
        list(block1 = endoapply(fmr, utils::head, 4L)),
        "susie",
        "f()"
    ))
})

test_that("assembleCtwasInputs warns about narrowed SuSiE weights", {
    skip_if_not_installed("ctwas")
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    cnd <- .fitCoverageCondition(assembleCtwasInputs(
        gwasSumStats = .ctp_makeGwasSumstats(),
        twasWeights = list(
            block1 = .ctp_narrowedSusieWeights(),
            block2 = .ctp_narrowedSusieWeights(7:11)
        ),
        method = "susie"
    ))
    expect_s3_class(cnd, "pecotmrFitCoverage")
    expect_match(conditionMessage(cnd), "assembleCtwasInputs\\(\\): 2 row")
})

# A stand-in for assembleCtwasInputs() that warns as it does about a narrowed
# row, then stops the run.
.ctp_warnThenStop <- function(...) {
    warn("assembleCtwasInputs(): narrowed", class = "pecotmrFitCoverage")
    abort("assembled")
}

test_that("ctwasPipeline warns once, under its own name", {
    skip_if_not_installed("ctwas")
    local_mocked_bindings(
        assembleCtwasInputs = .ctp_warnThenStop,
        .package = "pecotmr"
    )
    msgs <- character(0)
    expect_error(
        withCallingHandlers(
            ctwasPipeline(
                gwasSumStats = .ctp_makeGwasSumstats(),
                twasWeights = .ctp_narrowedSusieWeights(),
                method = "susie"
            ),
            pecotmrFitCoverage = function(w) {
                msgs <<- c(msgs, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        "assembled"
    )
    expect_length(msgs, 1L)
    expect_match(msgs, "ctwasPipeline\\(\\): 1 row")
})
