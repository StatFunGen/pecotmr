context("ctwasPipeline")

# ===========================================================================
# Strategy: ctwas::ctwas_sumstats does the heavy work. We mock it to a
# function that just returns its inputs back, so we can verify how the
# pipeline assembles z_snp / weights / region_info / LD loader inputs.
# ===========================================================================

# 12 = 6 variants x the 2 blocks the default GwasSumStats fixture builds. The
# LD-sketch identity check compares panels between the GWAS and weight sides,
# so every fixture has to draw on the same one.
.ctp_makeHandle <- function(snp_n = 12L, n_samples = 30L) {
    # Use a per-process tempfile so .ctwasLdPanelKey's file.exists check
    # succeeds against the fixture handle (real LD-sketch payloads exist
    # by construction; mock fixtures need an equivalent on-disk anchor).
    gdsPath <- file.path(tempdir(), "ctp_sketch.gds")
    if (!file.exists(gdsPath)) {
        file.create(gdsPath)
    }
    positions <- seq(100L, by = 100L, length.out = snp_n)
    # SNP IDs follow the canonical chr:pos:A2:A1 layout so allele
    # harmonization inside .ctwasBuildWeights / .ctwasHarmonizeWeights can
    # parse them via parseVariantId().
    snpIds <- sprintf("chr1:%d:G:A", positions)
    new(
        "GenotypeHandle",
        path = gdsPath,
        format = "gds",
        snpInfo = data.frame(
            SNP = snpIds,
            CHR = rep("1", snp_n),
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

# Canonical SNP IDs the fixtures use. .ctp_makeHandle() emits these in
# chr:pos:A2:A1 form; tests reference them by index via .ctp_snpId(i).
.ctp_snpId <- function(i) sprintf("chr1:%d:G:A", 100L * i)

.ctp_mockExtractor <- function(seed = 5, n_samples = 30L) {
    function(handle, snpIdx, meanImpute = TRUE) {
        set.seed(seed)
        panel <- matrix(
            rbinom(n_samples * nrow(handle@snpInfo), 2, 0.3),
            nrow = n_samples,
            ncol = nrow(handle@snpInfo),
            dimnames = list(handle@sampleIds, handle@snpInfo$SNP)
        )
        sub <- panel[, snpIdx, drop = FALSE]
        rr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", handle@snpInfo$CHR[snpIdx]),
            ranges = IRanges::IRanges(
                start = handle@snpInfo$BP[snpIdx],
                width = 1L
            )
        )
        S4Vectors::mcols(rr) <- S4Vectors::DataFrame(
            SNP = handle@snpInfo$SNP[snpIdx],
            A1 = handle@snpInfo$A1[snpIdx],
            A2 = handle@snpInfo$A2[snpIdx]
        )
        cd <- S4Vectors::DataFrame(
            sampleId = handle@sampleIds,
            row.names = handle@sampleIds
        )
        dosage <- t(sub)
        rownames(dosage) <- handle@snpInfo$SNP[snpIdx]
        colnames(dosage) <- handle@sampleIds
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
        SNP = vapply(idx, .ctp_snpId, character(1)),
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
        study = rep("G1", n),
        entry = map(seq_len(n), .ctp_blockVariants),
        genome = "hg19",
        ldSketch = .ctp_makeHandle(snp_n = 6L * n),
        blockId = blockIds,
        qcInfo = if (qc) list(step1 = "ok") else list()
    )
}

# `variantIdx` picks which of the fixture's 12 variants carry the gene's
# weights; it has to name five of them, one per weight. The default sits
# inside block 1 -- pass 7:11 for a gene that lives in block 2.
.ctp_makeTwasWeights <- function(variantIdx = 1:5) {
    e <- twasWeightsRow(
        variantIds = vapply(variantIdx, .ctp_snpId, character(1)),
        weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
    )
    TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e),
        ldSketch = .ctp_makeHandle()
    )
}

# The same per-gene weights expressed as a QtlFineMappingResult weight source:
# the topLoci posterior_mean carries the weight vector (what resolveWeights
# reads). Mirrors .ctp_makeTwasWeights so the two sources are comparable.
.ctp_makeFmrWeightSource <- function() {
    vids <- vapply(1:5, .ctp_snpId, character(1))
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
    fe <- fineMappingRow(variantIds = vids, susieFit = list(), topLoci = tl)
    QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fe),
        ldSketch = .ctp_makeHandle()
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
            twasWeights = .ctp_makeTwasWeights()
        ),
        "must be a GwasSumStats whose elements are LD blocks"
    )
})

test_that("ctwasPipeline: rejects a single-block named list", {
    skip_if_not_installed("ctwas")
    expect_error(
        ctwasPipeline(
            gwasSumStats = .ctp_makeGwasSumstats(blockIds = "block1"),
            twasWeights = list(block1 = .ctp_makeTwasWeights())
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
            twasWeights = list(block1 = tw, block2 = tw)
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
            twasWeights = list(blockA = tw, blockC = tw)
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
    # Two blocks supply zSnp; only block1 supplies TwasWeights.
    inputs <- assembleCtwasInputs(
        gwasSumStats = ss,
        twasWeights = list(block1 = tw)
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
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(twasWeightsRow(
            variantIds = vapply(1:5, .ctp_snpId, character(1)),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0),
            standardized = TRUE
        )),
        ldSketch = .ctp_makeHandle()
    )
    fmr <- .ctp_makeFmrWeightSource()
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    itw <- assembleCtwasInputs(
        gwasSumStats = ss,
        twasWeights = list(block1 = tw_std, block2 = tw_std)
    )
    ifmr <- assembleCtwasInputs(
        gwasSumStats = ss,
        twasWeights = list(block1 = fmr, block2 = fmr)
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
    mkBlockGss <- function(study, snpIds, qc = TRUE) {
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
            study = study,
            entry = list(gr),
            genome = "hg19",
            ldSketch = .ctp_makeHandle(),
            qcInfo = if (qc) list(step1 = "ok") else list()
        )
    }
    ss1 <- mkBlockGss("G1", vapply(1:3, .ctp_snpId, character(1)))
    ss2 <- mkBlockGss("G2", vapply(4:6, .ctp_snpId, character(1)))
    # Cross-boundary weights: v2..v5 (4 variants spanning both blocks).
    crossEntry <- twasWeightsRow(
        variantIds = vapply(2:5, .ctp_snpId, character(1)),
        weights = c(0.1, 0.2, 0.3, 0.4)
    )
    tw <- TwasWeights(
        study = "G1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(crossEntry),
        ldSketch = .ctp_makeHandle()
    )
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    both <- c(ss1, ss2)
    mcols(both)$blockId <- c("block1", "block2")
    inputs <- assembleCtwasInputs(
        gwasSumStats = both,
        twasWeights = list(block1 = tw)
    )
    # ctwas fine-maps one region at a time, so the FITTED vector holds only the
    # variants this block has LD for -- handing susie_rss the other two produces
    # a non-finite ELBO. The gene's SPAN still covers all four, which is what
    # ctwas::get_boundary_genes reads to route the gene to merge_regions; a span
    # clipped to the block would hide the boundary gene entirely.
    entry <- inputs$weights[[1L]]
    wgt <- entry$wgt
    expect_equal(nrow(wgt), 2L)
    expect_setequal(rownames(wgt), vapply(2:3, .ctp_snpId, character(1)))
    expect_equal(nrow(entry$R_wgt), 2L)
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
            twasWeights = .ctp_makeTwasWeights()
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
            twasZ = "not a GRanges"
        ),
        "must be a GRanges"
    )
})

test_that("ctwasPipeline: rejects unknown groupPriorVarStructure value", {
    skip_if_not_installed("ctwas")
    inp <- .ctp_makeMultiBlockInputs()
    expect_error(
        ctwasPipeline(
            gwasSumStats = inp$gwasSumStats,
            twasWeights = inp$twasWeights,
            groupPriorVarStructure = "bogus"
        ),
        "groupPriorVarStructure"
    )
})

# ===========================================================================
# .ctwasRequireMatchingLdSketches
# ===========================================================================

test_that(".ctwasRequireMatchingLdSketches: NULL twas-side handle is allowed", {
    twNoLd <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(twasWeightsRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:5)),
            weights = rep(0.1, 5)
        )),
        ldSketch = NULL
    )
    expect_silent(pecotmr:::.ctwasRequireMatchingLdSketches(
        twLd = NULL,
        gwasLd = .ctp_makeHandle()
    ))
})

test_that(".ctwasRequireMatchingLdSketches: panel-size mismatch errors", {
    twLd <- .ctp_makeHandle(snp_n = 5L)
    gwasLd <- .ctp_makeHandle(snp_n = 6L)
    expect_error(
        pecotmr:::.ctwasRequireMatchingLdSketches(twLd, gwasLd),
        "ldSketch panels differ in size"
    )
})

# ===========================================================================
# Input-building helpers
# ===========================================================================

test_that(".ctwasBuildZSnp: produces a flat data.frame keyed by SNP/study", {
    ss <- .ctp_makeGwasSumstats(blockIds = "block1")
    df <- pecotmr:::.ctwasBuildZSnp(ss)
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 6L)
    expect_setequal(
        colnames(df),
        c("id", "chrom", "pos", "A1", "A2", "z", "study")
    )
    expect_setequal(df$id, vapply(1:6, .ctp_snpId, character(1)))
    expect_setequal(unique(df$study), "G1")
})

test_that(".ctwasBuildSingleRegionInfo: pulls chrom + bp span from the GWAS block entry", {
    # Bounds come from the block's GWAS variants (the GwasSumStats entry), NOT
    # the LD sketch — many blocks can share one whole-chromosome LD payload.
    ri <- pecotmr:::.ctwasBuildSingleRegionInfo(
        "block1",
        .ctp_makeGwasSumstats(blockIds = "block1")
    )
    expect_equal(ri$region_id, "block1")
    expect_equal(ri$chrom, 1L)
    expect_equal(ri$start, 100L)
    expect_equal(ri$stop, 600L)
})

test_that(".ctwasBuildSingleRegionInfo: uses the block entry span, not the wider shared LD sketch", {
    # Regression: many LD blocks can share one whole-chromosome LD payload, so
    # the sketch span (here BP 100-600) is NOT the block's span. The entry here
    # covers only 200-400; region bounds must follow the entry, otherwise every
    # block collapses to the whole-chromosome span and every SNP is assigned to
    # every region (inflating SNP group_size and crushing the gene PIP).
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
        study = "G1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ctp_makeHandle(),
        qcInfo = list(step1 = "ok")
    )
    ri <- pecotmr:::.ctwasBuildSingleRegionInfo("blockX", gss)
    expect_equal(ri$start, 200L) # entry min, not sketch min (100)
    expect_equal(ri$stop, 400L) # entry max, not sketch max (600)
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
        study = "G1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ctp_makeHandle(),
        qcInfo = list(step1 = "ok")
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

test_that(".ctwasLdPanelKey: returns the on-disk path for an existing GDS sketch", {
    handle <- .ctp_makeHandle()
    key <- pecotmr:::.ctwasLdPanelKey(handle)
    expect_true(file.exists(key))
    expect_equal(key, getPath(handle))
})

test_that(".ctwasLdPanelKey: errors when no candidate file exists", {
    ghost <- new(
        "GenotypeHandle",
        path = file.path(tempdir(), "does_not_exist_for_test.pgen_stem"),
        format = "plink2",
        snpInfo = data.frame(
            SNP = "chr1:100:A:G",
            CHR = "1",
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
        "could not derive an existing LD-file token"
    )
})

test_that(".ctwasResolveMethod: caller-supplied method wins when present", {
    e <- twasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
        weights = c(0.1, 0.2, 0.3)
    )
    tw <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "mrash",
        entry = list(e),
        ldSketch = .ctp_makeHandle()
    )
    expect_equal(pecotmr:::.ctwasResolveMethod(list(r1 = tw), "mrash"), "mrash")
})

test_that(".ctwasResolveMethod: caller-supplied unknown method errors", {
    e <- twasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
        weights = c(0.1, 0.2, 0.3)
    )
    tw <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "mrash",
        entry = list(e),
        ldSketch = .ctp_makeHandle()
    )
    expect_error(
        pecotmr:::.ctwasResolveMethod(list(r1 = tw), "bogus"),
        "not present in TwasWeights"
    )
})

test_that(".ctwasResolveMethod: defaults to ensemble when present among multiple", {
    mkTw <- function(m) {
        e <- twasWeightsRow(
            variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
            weights = c(0.1, 0.2, 0.3)
        )
        TwasWeights(
            study = "Q1",
            context = "c1",
            trait = "t1",
            method = m,
            entry = list(e),
            ldSketch = .ctp_makeHandle()
        )
    }
    # Build a multi-method TwasWeights by stitching two methods together.
    tw <- TwasWeights(
        study = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("mrash", "ensemble"),
        entry = list(
            twasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.1, 0.2, 0.3)
            ),
            twasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.4, 0.5, 0.6)
            )
        ),
        ldSketch = .ctp_makeHandle()
    )
    expect_equal(pecotmr:::.ctwasResolveMethod(list(r1 = tw)), "ensemble")
})

test_that(".ctwasResolveMethod: single method auto-picked when only one available", {
    e <- twasWeightsRow(
        variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
        weights = c(0.1, 0.2, 0.3)
    )
    tw <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "mrash",
        entry = list(e),
        ldSketch = .ctp_makeHandle()
    )
    expect_equal(pecotmr:::.ctwasResolveMethod(list(r1 = tw)), "mrash")
})

test_that(".ctwasResolveMethod: multi-method + no ensemble + no caller method errors", {
    tw <- TwasWeights(
        study = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("mrash", "susie"),
        entry = list(
            twasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.1, 0.2, 0.3)
            ),
            twasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.4, 0.5, 0.6)
            )
        ),
        ldSketch = .ctp_makeHandle()
    )
    expect_error(
        pecotmr:::.ctwasResolveMethod(list(r1 = tw)),
        "Supply a `method` argument"
    )
})

# ---------------------------------------------------------------------------
# .ctwasResolveMethods (plural) — the pipeline's fan-out resolver
# ---------------------------------------------------------------------------
.ctp_multiMethodTw <- function(methods = c("mrash", "susie")) {
    TwasWeights(
        study = rep("Q1", length(methods)),
        context = rep("c1", length(methods)),
        trait = rep("t1", length(methods)),
        method = methods,
        entry = lapply(methods, function(m) {
            twasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.1, 0.2, 0.3)
            )
        }),
        ldSketch = .ctp_makeHandle()
    )
}

test_that(".ctwasResolveMethods: unspecified + multiple + no ensemble -> ALL methods", {
    tw <- .ctp_multiMethodTw(c("mrash", "susie"))
    expect_setequal(
        pecotmr:::.ctwasResolveMethods(list(r1 = tw)),
        c("mrash", "susie")
    )
})

test_that(".ctwasResolveMethods: ensemble wins, single auto-picks, explicit restricts", {
    expect_equal(
        pecotmr:::.ctwasResolveMethods(
            list(r1 = .ctp_multiMethodTw(c("mrash", "ensemble")))
        ),
        "ensemble"
    )
    expect_equal(
        pecotmr:::.ctwasResolveMethods(
            list(r1 = .ctp_multiMethodTw("mrash"))
        ),
        "mrash"
    )
    expect_equal(
        pecotmr:::.ctwasResolveMethods(
            list(r1 = .ctp_multiMethodTw(c("mrash", "susie")), r2 = NULL),
            method = "susie"
        ),
        "susie"
    )
    expect_error(
        pecotmr:::.ctwasResolveMethods(
            list(r1 = .ctp_multiMethodTw("mrash")),
            method = "lasso"
        ),
        "not present"
    )
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
test_that(".ctwasPlaceByAnchor: homes each gene by start(region), half-open, chr-agnostic", {
    region <- GenomicRanges::GRanges(
        c("chr1", "chr1", "22", "chr2"),
        IRanges::IRanges(
            start = c(500, 1500, 800, 100),
            end = c(600, 1600, 900, 200)
        )
    )
    grid <- setNames(
        vector("list", 3),
        c("chr1_1_1000", "chr1_1000_2000", "chr22_1_1000")
    )
    home <- pecotmr:::.ctwasPlaceByAnchor(region, grid)
    # chr1:500 -> [1,1000); chr1:1500 -> [1000,2000); "22":800 -> chr22 block
    # (chr prefix normalized); chr2 -> no block -> NA.
    expect_equal(home, c("chr1_1_1000", "chr1_1000_2000", "chr22_1_1000", NA))
})

test_that(".ctwasPlaceByAnchor: block interval is half-open [start, stop)", {
    region <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = c(1000, 999), end = c(1000, 999))
    )
    grid <- setNames(vector("list", 2), c("chr1_1_1000", "chr1_1000_2000"))
    # start==1000 falls in the SECOND block (>= 1000), 999 in the first.
    expect_equal(
        pecotmr:::.ctwasPlaceByAnchor(region, grid),
        c("chr1_1000_2000", "chr1_1_1000")
    )
})

test_that(".ctwasIsPreBucketed / .ctwasCombineWeightSources: dispatch flat vs pre-bucketed", {
    tw <- .ctp_makeTwasWeights()
    grid <- list(block1 = 1, block2 = 2)
    expect_true(pecotmr:::.ctwasIsPreBucketed(list(block1 = tw), grid)) # named list
    expect_false(pecotmr:::.ctwasIsPreBucketed(tw, grid)) # flat S4
    expect_false(pecotmr:::.ctwasIsPreBucketed(list(tw), grid)) # unnamed
    expect_s4_class(pecotmr:::.ctwasCombineWeightSources(tw), "TwasWeights")
    expect_s4_class(
        pecotmr:::.ctwasCombineWeightSources(list(tw)),
        "TwasWeights"
    )
})

test_that(".ctwasBucketWeights: places a flat 2-gene source into its home blocks", {
    mkE <- function() {
        twasWeightsRow(
            variantIds = vapply(1:5, .ctp_snpId, character(1)),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
        )
    }
    # Placement anchors on the GENE position (traitPos), not the weight
    # variants: both genes carry the same five variants here, so only the
    # anchor can tell them apart. block1 spans 100-600, block2 spans 700-1200,
    # so gA @chr1:100 homes to block1 and gB @chr1:800 to block2.
    tw <- TwasWeights(
        study = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("gA", "gB"),
        method = c("susie", "susie"),
        entry = list(mkE(), mkE()),
        traitPos = GenomicRanges::GRanges(
            c("chr1", "chr1"),
            IRanges::IRanges(c(100, 800), c(150, 850))
        ),
        ldSketch = .ctp_makeHandle()
    )
    grid <- pecotmr:::.ctwasGwasByBlock(.ctp_makeGwasSumstats())
    bucketed <- pecotmr:::.ctwasBucketWeights(tw, grid)
    expect_named(bucketed, c("block1", "block2"), ignore.order = TRUE)
    expect_equal(as.character(bucketed[["block1"]]$trait), "gA")
    expect_equal(as.character(bucketed[["block2"]]$trait), "gB")
    expect_false(is.null(getLdSketch(bucketed[["block1"]])))
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
    expect_equal(nrow(getFinemap(rows[[1L]]$entry)), 2L)
    expect_equal(getCtwasParam(rows[[1L]]$entry)$group_prior[["c1"]], 0.01)
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
        vapply(rows, function(r) r$context, ""),
        c("brain", "liver")
    )
    expect_true(all(
        vapply(rows, function(r) r$jointContexts, "") == "brain,liver"
    ))
    # Each per-context row keeps only its own genes but shares the joint param.
    expect_equal(nrow(getFinemap(rows[[1L]]$entry)), 2L)
    expect_named(
        getCtwasParam(rows[[1L]]$entry)$group_prior,
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
    expect_null(getFinemap(rows[[1L]]$entry))
    expect_null(getSusieAlpha(rows[[1L]]$entry))
})

test_that(".ctwasRunToRows: susieAlpha subset mirrors finemap per context", {
    run <- .ctp_runResult(
        c("blk|Q1|c1|gA|susie", "blk|Q1|c1|gB|susie"),
        c(0.9, 0.2),
        c(c1 = 0.01, SNP = 1e-4)
    )
    rows <- pecotmr:::.ctwasRunToRows(run, "D1", "susie")
    sa <- getSusieAlpha(rows[[1L]]$entry)
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
    expect_true(all(getFinemap(rowsOff[[1L]]$entry)$type == "gene"))
    # keepSnps: one extra study = context = "SNP" row carrying the SNP background.
    rowsOn <- pecotmr:::.ctwasRunToRows(run, "D1", "susie", keepSnps = TRUE)
    expect_length(rowsOn, 2L)
    snpRow <- rowsOn[[2L]]
    expect_equal(snpRow$study, "SNP")
    expect_equal(snpRow$context, "SNP")
    # The payloads are ranged now: SNP rows carry the coordinates their own
    # variant ids encode, so the count is length() and the columns are mcols.
    fm <- getFinemap(snpRow$entry)
    expect_s4_class(fm, "GRanges")
    expect_length(fm, 2L)
    expect_true(all(S4Vectors::mcols(fm)$type == "SNP"))
    expect_length(getSusieAlpha(snpRow$entry), 2L)
})

test_that(".ctwasFilterMethod: subsets rows to the requested method", {
    tw <- TwasWeights(
        study = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("mrash", "susie"),
        entry = list(
            twasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.1, 0.2, 0.3)
            ),
            twasWeightsRow(
                variantIds = sprintf("chr1:%d:A:G", 100L * (1:3)),
                weights = c(0.4, 0.5, 0.6)
            )
        ),
        ldSketch = .ctp_makeHandle()
    )
    twSub <- pecotmr:::.ctwasFilterMethod(tw, "susie")
    expect_equal(nrow(twSub), 1L)
    expect_equal(as.character(twSub$method), "susie")
})

# Build an ldPanel fixture (matches .ctwasComputeFullPanelLd's return
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

test_that(".ctwasBuildWeights: scales non-standardized weights by sqrt(variance)", {
    panel <- .ctp_makeLdPanel()
    # Replace the default unit variance with non-trivial values; raw
    # weights should be multiplied by sqrt(variance) before reaching the
    # final wgt matrix.
    panel$variance <- setNames(c(0.5, 1, 2, 4, 8, 16), panel$snpInfo$id)
    ids5 <- vapply(1:5, .ctp_snpId, character(1))
    rawW <- c(0.1, 0.2, 0.3, 0.4, 0.5)
    tw <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(twasWeightsRow(variantIds = ids5, weights = rawW)),
        ldSketch = .ctp_makeHandle()
    )
    wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
    expected <- unname(rawW * sqrt(panel$variance[ids5]))
    expect_equal(as.numeric(wl[[1L]]$wgt), expected, tolerance = 1e-12)
})

test_that(".ctwasBuildWeights: standardized weights bypass variance scaling", {
    panel <- .ctp_makeLdPanel()
    panel$variance <- setNames(c(0.5, 1, 2, 4, 8, 16), panel$snpInfo$id)
    ids5 <- vapply(1:5, .ctp_snpId, character(1))
    rawW <- c(0.1, 0.2, 0.3, 0.4, 0.5)
    tw <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(twasWeightsRow(
            variantIds = ids5,
            weights = rawW,
            standardized = TRUE
        )),
        ldSketch = .ctp_makeHandle()
    )
    wl <- pecotmr:::.ctwasBuildWeights(tw, panel)
    expect_equal(as.numeric(wl[[1L]]$wgt), rawW, tolerance = 1e-12)
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

test_that(".ctwasRenormalizeSusieWeights: alpha renorm + colSums", {
    # Original fit covers 4 variants; drop variant 4, renormalize over {1,2,3}.
    origVids <- sprintf("chr1:%d:A:G", 100L * (1:4))
    origW <- c(0.1, 0.2, 0.3, 0.4)
    # Toy lbf_variable: 2 effects, 4 variants. Constructed so the kept
    # subset yields easily-predictable softmax weights.
    lbf <- rbind(
        c(0, 0, 0, 100), # effect 1: only v4 has signal
        c(10, 0, -10, 0)
    ) # effect 2: v1 dominates
    mu <- matrix(c(1, 2, 3, 4, 1, 2, 3, 4), nrow = 2, byrow = TRUE)
    xCol <- rep(1, 4)
    fits <- list(
        alpha = lbfToAlpha(lbf),
        mu = mu,
        X_column_scale_factors = xCol
    )
    out <- pecotmr:::.ctwasRenormalizeSusieWeights(
        fits,
        origVids = origVids,
        origW = origW,
        keptIdx = c(1L, 2L, 3L),
        harmonizedW = origW[1:3]
    )
    # Effect 1 lbf over {v1,v2,v3} = c(0,0,0) -> uniform alpha = 1/3
    # Effect 2 lbf over {v1,v2,v3} = c(10,0,-10) -> v1 ≈ 1
    # weight[v1] = (1/3)*1 + ~1*1 = ~1.33
    expect_equal(length(out), 3L)
    expect_true(out[[1L]] > out[[2L]] && out[[1L]] > out[[3L]])
})

test_that(".ctwasRenormalizeSusieWeights: returns NULL on fit/entry dimension mismatch", {
    fits <- list(
        alpha = matrix(1 / 5, 2, 5),
        mu = matrix(0, 2, 5),
        X_column_scale_factors = rep(1, 5)
    )
    out <- pecotmr:::.ctwasRenormalizeSusieWeights(
        fits,
        # entry says 3, fit covers 5 -> mismatch
        origVids = sprintf("chr1:%d:A:G", 100L * (1:3)),
        origW = rep(0.1, 3),
        keptIdx = 1:3,
        harmonizedW = rep(0.1, 3)
    )
    expect_null(out)
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

test_that(".ctwasRenormalizeSusieWeights: skips Omega-weighted susieInf fits", {
    # coef.susie for susieInf / susieAsh is colSums(alpha * mu)/scale +
    # theta/scale. Recomputing from alpha and mu alone would drop theta, so
    # these fits fall back to the plain subset-w-and-R path instead.
    origVids <- sprintf("chr1:%d:A:G", 100L * (1:4))
    origW <- c(0.1, 0.2, 0.3, 0.4)
    base <- list(
        alpha = lbfToAlpha(rbind(c(2, 0, -2, 0), c(0, 2, 0, -2))),
        mu = rbind(c(1, 2, 3, 4), c(1, 2, 3, 4)),
        X_column_scale_factors = rep(1, 4)
    )
    args <- list(
        origVids = origVids,
        origW = origW,
        keptIdx = c(1L, 2L, 3L),
        harmonizedW = origW[1:3]
    )
    expect_type(
        do.call(
            pecotmr:::.ctwasRenormalizeSusieWeights,
            c(list(base), args)
        ),
        "double"
    )
    withTheta <- c(base, list(theta = rep(0.1, 4)))
    expect_null(do.call(
        pecotmr:::.ctwasRenormalizeSusieWeights,
        c(list(withTheta), args)
    ))
    withOmega <- c(base, list(omega_weights = rep(1, 4)))
    expect_null(do.call(
        pecotmr:::.ctwasRenormalizeSusieWeights,
        c(list(withOmega), args)
    ))
})

test_that(".ctwasSnpInfoForGwasBlock: restricts panel snpInfo to block GWAS variants", {
    ss <- .ctp_makeGwasSumstats()
    panelInfo <- data.frame(
        chrom = 1L,
        id = vapply(1:6, .ctp_snpId, character(1)), # whole panel
        pos = seq(100L, by = 100L, length.out = 6L),
        alt = "A",
        ref = "G",
        stringsAsFactors = FALSE
    )
    blockInfo <- pecotmr:::.ctwasSnpInfoForGwasBlock(ss, panelInfo)
    # Restricted to the variant IDs on the GwasSumStats entry.
    expect_true(all(
        blockInfo$id %in%
            as.character(S4Vectors::mcols(ss[[1L]])$SNP)
    ))
    expect_true(nrow(blockInfo) <= nrow(panelInfo))
})

test_that(".ctwasBuildWeights: keys per-tuple weights, adds gene metadata", {
    tw <- .ctp_makeTwasWeights()
    panel <- .ctp_makeLdPanel()
    ids5 <- vapply(1:5, .ctp_snpId, character(1))
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
    # R_wgt is a 5x5 slice of the cached panel R
    expect_true(is.matrix(wl[[1L]]$R_wgt))
    expect_equal(dim(wl[[1L]]$R_wgt), c(5L, 5L))
    expect_equal(rownames(wl[[1L]]$R_wgt), ids5)
    expect_equal(wl[[1L]]$n_wgt, 5L)
    # And it is literally a slice of the panel R (no recompute path).
    expect_equal(wl[[1L]]$R_wgt, panel$R[ids5, ids5])
})

test_that(".ctwasBuildWeights: drops variants not present in the LD panel", {
    ids3 <- vapply(1:3, .ctp_snpId, character(1))
    missing <- c("chr1:99900:G:A", "chr1:99910:G:A") # not in panel
    tw <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(twasWeightsRow(
            variantIds = c(ids3, missing),
            weights = c(0.1, 0.2, 0.3, 0.4, 0.5)
        )),
        ldSketch = .ctp_makeHandle()
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
    ids5 <- vapply(1:5, .ctp_snpId, character(1))
    blockIds <- vapply(c(1, 2, 4), .ctp_snpId, character(1))
    tw <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(twasWeightsRow(
            variantIds = ids5,
            weights = c(0.1, 0.2, 0.3, 0.4, 0.5)
        )),
        ldSketch = .ctp_makeHandle()
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

test_that(".ctwasComputeFullPanelLd: extracts once + returns cached R + snpInfo + variance", {
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    out <- pecotmr:::.ctwasComputeFullPanelLd(.ctp_makeHandle(snp_n = 6L))
    ids6 <- vapply(1:6, .ctp_snpId, character(1))
    expect_named(out, c("R", "snpInfo", "variance"))
    expect_true(is.matrix(out$R))
    expect_equal(dim(out$R), c(6L, 6L))
    expect_equal(rownames(out$R), ids6)
    expect_setequal(
        colnames(out$snpInfo),
        c("chrom", "id", "pos", "alt", "ref")
    )
    expect_named(out$variance, ids6)
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
        twasWeights = inp$twasWeights
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
        twasZ = twasZ
    )
    expect_equal(inputs$z_gene$id, "Q1|c1|t1|susie")
})

# ===========================================================================
# Step-wise dispatch: estCtwasParam → screenCtwasRegions → finemapCtwasRegions
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
        twasWeights = inp$twasWeights
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
    # Output is a CtwasResult: one (gwasStudy, study, context, method) row.
    expect_s4_class(out, "CtwasResult")
    expect_equal(nrow(out), 1L)
    expect_equal(getMethodNames(out), "susie")
    expect_equal(getStudy(out), "Q1")
    expect_equal(getContexts(out), "c1")
    expect_equal(as.character(out$gwasStudy), "G1")
    # The run's jointly-estimated param is carried on the row.
    expect_equal(
        unname(getCtwasParam(pecotmr:::.collectionEntry(out, 1L))$group_prior),
        c(0.1, 0.0001)
    )
    # getFinemap aggregates the per-gene rows, tagged with run identity.
    fm <- getFinemap(out)
    expect_equal(nrow(fm), 2L)
    expect_true(all(
        c("gwasStudy", "study", "context", "method", "susie_pip") %in% names(fm)
    ))
    expect_setequal(fm$id, c("block1|Q1|c1|t1|susie", "block2|Q1|c1|t1|susie"))
    # getSusieAlpha aggregates the fuller per-effect table, likewise tagged.
    sa <- getSusieAlpha(out)
    expect_equal(nrow(sa), 2L)
    expect_true(all(c("gwasStudy", "context", "susie_alpha") %in% names(sa)))
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
    out <- ctwasPipeline(inp$gwasSumStats, inp$twasWeights, keepSnps = TRUE)
    # A gene (c1) row plus the dedicated SNP row.
    expect_equal(nrow(out), 2L)
    expect_setequal(getContexts(out), c("c1", "SNP"))
    fm <- getFinemap(out)
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
        weights = setNames(list(list(wgt = 1)), "block1|Q1|c1|t1|susie"),
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

test_that("mergeCtwasBoundaryRegions: LD path splices updated_* back + merge_res", {
    skip_if_not_installed("ctwas")
    captured <- NULL
    local_mocked_bindings(
        postprocess_region_merging = function(...) {
            captured <<- list(...)
            .ctp_mergeReturn()
        },
        .package = "ctwas"
    )
    out <- mergeCtwasBoundaryRegions(
        .ctp_finemapResult(hasLd = TRUE),
        pipThresh = 0.5
    )
    # LD path passed the loader closures + pip_thresh through.
    expect_equal(captured$pip_thresh, 0.5)
    expect_true(is.function(captured$LD_loader_fun))
    # updated_* spliced onto the result; merge_res carries the full postprocess out.
    expect_equal(out$finemap_res$susie_pip, 0.99)
    expect_equal(out$susie_alpha_res$susie_alpha, 0.99)
    expect_equal(out$region_info$region_id, "block1_merged")
    expect_true(!is.null(out$merge_res))
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

test_that("mergeCtwasBoundaryRegions: empty first-pass finemap returns unchanged", {
    skip_if_not_installed("ctwas")
    fmr <- .ctp_finemapResult()
    fmr$finemap_res <- fmr$finemap_res[0, ]
    expect_message(
        out <- mergeCtwasBoundaryRegions(fmr),
        "no first-pass finemap"
    )
    expect_null(out$merge_res)
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
        mergeBoundary = TRUE
    )
    expect_true(merged) # merging ran
    expect_s4_class(out, "CtwasResult")
    # The post-merge PIP (0.99) is what the decomposition carries.
    expect_equal(getFinemap(out)$susie_pip, 0.99)
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
    expect_equal(getStudy(cr), "Q1") # derived from the gene ids
    expect_equal(getContexts(cr), "c1")
    expect_equal(getMethodNames(cr), "susie")
    expect_equal(getFinemap(cr)$susie_pip, 0.90)
    expect_false(is.null(getSusieAlpha(cr)))
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
    expect_setequal(getContexts(cr), c("c1", "SNP"))
})

test_that("asCtwasResult: errors when the weights mix methods", {
    fmr <- .ctp_finemapResult()
    fmr$weights <- setNames(
        list(1, 2),
        c("block1|Q1|c1|t1|susie", "block1|Q1|c1|t2|lasso")
    )
    expect_error(asCtwasResult(fmr), "mixes weight methods")
})

test_that("estCtwasParam: fallbackToPrefit recovers from accurate-EM NaN divergence", {
    skip_if_not_installed("ctwas")
    inp <- .ctp_makeMultiBlockInputs()
    # Mock est_param to throw the documented NaN error, and fit_EM to
    # produce a stub prefit result. Verify estCtwasParam catches the
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
    # Without fallback: the NaN error propagates.
    expect_error(
        estCtwasParam(
            assembleCtwasInputs(inp$gwasSumStats, inp$twasWeights),
            fallbackToPrefit = FALSE
        ),
        "contains NAs"
    )
    # With fallback: prefit estimates are returned as the param.
    # .ctwasFitPrefitEm thin-scales the SNP group_prior (mirroring ctwas's
    # est_param), so the mocked SNP prior 1e-4 emerges as 1e-4 * thin
    # (default thin = 0.1) → 1e-5. The group_prior_var is not thinned.
    est <- estCtwasParam(
        assembleCtwasInputs(inp$gwasSumStats, inp$twasWeights),
        fallbackToPrefit = TRUE
    )
    expect_equal(unname(est$param$group_prior), c(0.05, 1e-5))
    expect_equal(unname(est$param$group_prior_var), c(4.0, 5.0))
})

test_that("estCtwasParam fallback drops degenerate regions before fit_EM", {
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
    est <- estCtwasParam(
        assembleCtwasInputs(inp$gwasSumStats, inp$twasWeights),
        fallbackToPrefit = TRUE
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
    res <- .ctwasFitPrefitEm(
        region_data,
        niter = 3L,
        groupPriorVarStructure = "shared_all",
        thin = 1,
        ncore = 1L
    )
    # the prefit EM ran on the valid region only and returns finite real group priors
    expect_true("SNP" %in% names(res$group_prior))
    expect_true(all(is.finite(res$group_prior)))
    # every region (valid + degenerate) is still listed in p_single_effect
    expect_setequal(res$p_single_effect$region_id, names(region_data))
})

test_that("estCtwasParam / screenCtwasRegions / finemapCtwasRegions can be called independently", {
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
    inputs <- assembleCtwasInputs(inp$gwasSumStats, inp$twasWeights)
    expect_true("region_info" %in% names(inputs))
    expect_true("LD_loader_fun" %in% names(inputs))
    # Step 2
    est <- estCtwasParam(inputs)
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
    gh <- qd@genotypes

    # Two 5-variant synthetic genes from the bundled panel, one anchored in
    # each LD block below. cTWAS's EM runs per region and cannot fit a block
    # with no genes in it, so a real partition needs a gene on both sides --
    # the previous per-region list dodged this by putting the SAME GWAS block
    # and the SAME gene under two region ids.
    geneOne <- twasWeightsRow(
        variantIds = gh@snpInfo$SNP[1:5],
        weights = c(0.1, -0.2, 0.05, 0.0, 0.3)
    )
    geneTwo <- twasWeightsRow(
        variantIds = gh@snpInfo$SNP[9:13],
        weights = c(0.2, 0.1, -0.15, 0.05, 0.0)
    )
    # A FLAT weight source is placed into blocks by `traitPos`, so each gene
    # needs the span it sits in.
    geneSpan <- function(i) {
        bp <- gh@snpInfo$BP[i]
        GenomicRanges::GRanges(
            str_c("chr", gh@snpInfo$CHR[i][[1L]]),
            IRanges::IRanges(min(bp), max(bp))
        )
    }
    tw <- TwasWeights(
        study = rep("study1", 2L),
        context = rep("brain", 2L),
        trait = c("ENSG_example", "ENSG_example2"),
        method = rep("susie", 2L),
        entry = list(geneOne, geneTwo),
        traitPos = c(geneSpan(1:5), geneSpan(9:13)),
        ldSketch = gh
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
        study = getStudy(gss),
        entry = list(unlist(gss)),
        genome = getGenome(gss),
        ldSketch = getLdSketch(gss),
        qcInfo = getQcInfo(gss),
        ldBlocks = twoBlocks
    )
    res <- suppressMessages(suppressWarnings(
        ctwasPipeline(
            gwasSumStats = gssTwoBlocks,
            twasWeights = tw,
            niter = 5L,
            niterPrefit = 2L,
            # Toy panel: relax the production filters that gate out tiny inputs.
            min_group_size = 1L,
            min_p_single_effect = 0,
            filter_L = FALSE
        )
    ))

    # The one-shot pipeline returns a CtwasResult: a single (brain, study1,
    # susie) row (single-context run, so no jointContexts).
    expect_s4_class(res, "CtwasResult")
    expect_equal(nrow(res), 1L)
    expect_equal(getContexts(res), "brain")
    expect_equal(getStudy(res), "study1")
    expect_equal(getMethodNames(res), "susie")
    expect_false("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
    # The gene we passed in was fine-mapped and shows up in the aggregated
    # finemap table, tagged with the run identity.
    fm <- getFinemap(res)
    expect_true(!is.null(fm) && nrow(fm) > 0L)
    expect_true(all(
        c("gwasStudy", "study", "context", "method", "susie_pip") %in% names(fm)
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
    sa <- getSusieAlpha(res)
    expect_true(!is.null(sa) && nrow(sa) > 0L)
    expect_true(all(c("susie_pip", "susie_alpha", "region_id") %in% names(sa)))
})

# ===========================================================================
# .ctwasFilterVariants — ported from R/ctwasWrapper.R::trimCtwasVariants
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

test_that(".ctwasBuildWeights: maxNumVariants caps the per-gene weight matrix", {
    data(qtlDatasetExample)
    qd <- qtlDatasetExample
    gh <- qd@genotypes
    vids <- gh@snpInfo$SNP[1:5]
    ent <- twasWeightsRow(
        variantIds = vids,
        weights = c(0.1, -0.2, 0.05, 0.3, 0.15)
    )
    tw <- TwasWeights(
        study = "study1",
        context = "brain",
        trait = "ENSG_example",
        method = "susie",
        entry = list(ent),
        ldSketch = gh
    )
    ldPanel <- pecotmr:::.ctwasComputeFullPanelLd(gh)
    wl <- pecotmr:::.ctwasBuildWeights(tw, ldPanel, maxNumVariants = 3L)
    expect_equal(wl[[1L]]$n_wgt, 3L)
    expect_equal(nrow(wl[[1L]]$wgt), 3L)
    # Top 3 by |w| from c(0.1, -0.2, 0.05, 0.3, 0.15): 0.3, -0.2, 0.15
    expect_setequal(rownames(wl[[1L]]$wgt), vids[c(4L, 2L, 5L)])
})

test_that(".ctwasBuildWeights: twasWeightCutoff drops low-magnitude variants", {
    data(qtlDatasetExample)
    qd <- qtlDatasetExample
    gh <- qd@genotypes
    vids <- gh@snpInfo$SNP[1:5]
    ent <- twasWeightsRow(
        variantIds = vids,
        # v1 (0.005) and v3 (0.001) will be dropped at cutoff 0.01
        weights = c(0.005, 0.2, 0.001, 0.3, 0.1)
    )
    tw <- TwasWeights(
        study = "study1",
        context = "brain",
        trait = "ENSG_example",
        method = "susie",
        entry = list(ent),
        ldSketch = gh
    )
    ldPanel <- pecotmr:::.ctwasComputeFullPanelLd(gh)
    wl <- pecotmr:::.ctwasBuildWeights(tw, ldPanel, twasWeightCutoff = 0.01)
    expect_equal(wl[[1L]]$n_wgt, 3L)
    expect_setequal(rownames(wl[[1L]]$wgt), vids[c(2L, 4L, 5L)])
})

# ===========================================================================
# mergeCtwasBoundaryRegions (step 4: boundary-gene region merging)
# ===========================================================================

test_that("mergeCtwasBoundaryRegions: no first-pass finemap_res returns unchanged", {
    fmr <- list(finemap_res = NULL, region_data = "rd")
    expect_identical(mergeCtwasBoundaryRegions(fmr), fmr)
    fmr0 <- list(finemap_res = data.frame()[0, ], region_data = "rd")
    expect_identical(mergeCtwasBoundaryRegions(fmr0), fmr0)
})

test_that("mergeCtwasBoundaryRegions: LD path forwards carried state + splices updated_*", {
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

test_that("mergeCtwasBoundaryRegions: no-LD path used when LD loaders are absent", {
    called <- NULL
    local_mocked_bindings(
        postprocess_region_merging_noLD = function(...) {
            called <<- "noLD"
            list(
                updated_finemap_res = data.frame(id = "g"),
                updated_susie_alpha_res = NULL
            )
        },
        .package = "ctwas"
    )
    fmr <- list(
        finemap_res = data.frame(id = "g", type = "gene", susie_pip = 0.6),
        susie_alpha_res = NULL,
        region_data = "rd",
        region_info = "ri",
        z_snp = "zs",
        z_gene = data.frame(id = "g"),
        weights = "w",
        snp_map = "sm",
        LD_map = NULL,
        LD_loader_fun = NULL,
        snpinfo_loader_fun = NULL,
        param = list(group_prior = 0.1, group_prior_var = 5)
    )
    out <- mergeCtwasBoundaryRegions(fmr)
    expect_equal(called, "noLD")
})

# ===========================================================================
# assembleCtwasInputs: remaining input-validation branches
# ---------------------------------------------------------------------------
# These fire before any heavy panel work; they go through the ctwas-gated
# entry point, so skip when ctwas is absent (otherwise the requireNamespace
# guard would surface a different error).
# ===========================================================================

test_that("assembleCtwasInputs: rejects an unnamed gwasSumStats list", {
    skip_if_not_installed("ctwas")
    ss <- .ctp_makeGwasSumstats()
    tw <- .ctp_makeTwasWeights()
    expect_error(
        assembleCtwasInputs(
            gwasSumStats = list(ss, ss), # unnamed
            twasWeights = list(block1 = tw)
        ),
        "must be a GwasSumStats whose elements are LD blocks"
    )
})

test_that("assembleCtwasInputs: an unnamed twasWeights list is a FLAT source (needs region)", {
    skip_if_not_installed("ctwas")
    ss <- .ctp_makeGwasSumstats()
    tw <- .ctp_makeTwasWeights() # no region provenance
    # An unnamed list is now treated as a flat weight source to place by region;
    # without region provenance that placement can't happen.
    expect_error(
        assembleCtwasInputs(
            gwasSumStats = ss,
            twasWeights = list(tw)
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
            twasWeights = list(block1 = tw)
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
            twasWeights = list(block1 = "not a weight source")
        ),
        "must be a TwasWeights or QtlFineMappingResult"
    )
})

test_that("assembleCtwasInputs: rejects a non-FineMappingResultBase fineMappingResult", {
    skip_if_not_installed("ctwas")
    ss <- .ctp_makeGwasSumstats()
    tw <- .ctp_makeTwasWeights()
    expect_error(
        assembleCtwasInputs(
            gwasSumStats = ss,
            twasWeights = list(block1 = tw),
            fineMappingResult = "not a FineMappingResult"
        ),
        "must be a FineMappingResultBase"
    )
})

test_that("assembleCtwasInputs: skips a block whose TwasWeights lacks the resolved method", {
    skip_if_not_installed("ctwas")
    ss <- .ctp_makeGwasSumstats()
    ids5 <- vapply(1:5, .ctp_snpId, character(1))
    mkTw <- function(m) {
        TwasWeights(
            study = "Q1",
            context = "c1",
            trait = "t1",
            method = m,
            entry = list(twasWeightsRow(
                variantIds = ids5,
                weights = c(0.1, 0.2, 0.3, 0.4, 0.5)
            )),
            ldSketch = .ctp_makeHandle()
        )
    }
    local_mocked_bindings(
        extractBlockGenotypes = .ctp_mockExtractor(),
        .package = "pecotmr"
    )
    # Resolved method is "ensemble" (present in block1); block2 carries only
    # "mrash", so .ctwasFilterMethod returns NULL and the block is skipped
    # in the second pass (the `if (is.null(twMethod)) next` branch).
    inputs <- assembleCtwasInputs(
        gwasSumStats = ss,
        twasWeights = list(block1 = mkTw("ensemble"), block2 = mkTw("mrash"))
    )
    expect_length(inputs$weights, 1L)
    expect_true(all(grepl("^block1\\|", names(inputs$weights))))
    expect_false(any(grepl("^block2\\|", names(inputs$weights))))
})

# ===========================================================================
# .ctwasResolveMethod / .ctwasFilterMethod edge branches
# ===========================================================================

test_that(".ctwasResolveMethod: errors when no method entries exist", {
    expect_error(pecotmr:::.ctwasResolveMethod(list()), "no method entries")
})

test_that(".ctwasFilterMethod: returns NULL when no row matches the method", {
    tw <- .ctp_makeTwasWeights() # method == "susie"
    expect_null(pecotmr:::.ctwasFilterMethod(tw, "mrash"))
})

# ===========================================================================
# .ctwasHarmonizeWeights / .ctwasRenormalizeSusieWeights early NULL returns
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

test_that(".ctwasRenormalizeSusieWeights: returns NULL when fit components are NULL", {
    fits <- list(
        lbf_variable = NULL,
        mu = matrix(0, 2, 3),
        X_column_scale_factors = rep(1, 3)
    )
    expect_null(pecotmr:::.ctwasRenormalizeSusieWeights(
        fits,
        origVids = sprintf("chr1:%d:A:G", 100L * (1:3)),
        origW = rep(0.1, 3),
        keptIdx = 1:3,
        harmonizedW = rep(0.1, 3)
    ))
})

# ===========================================================================
# .ctwasBuildWeights: per-gene skip branches + variance / renorm paths
# ===========================================================================

test_that("TwasWeightsRow: rejects mismatched variantIds/weights lengths", {
    # The weights are per-variant mcols on the entry's element now, so a
    # length mismatch cannot be represented. It is refused at construction
    # rather than producing a gene that gets silently dropped downstream.
    expect_error(
        twasWeightsRow(
            variantIds = vapply(1:5, .ctp_snpId, character(1)),
            weights = c(0.1, 0.2, 0.3)
        ),
        "length\\(weights\\) is 3 but 5 variants were supplied"
    )
})

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
    ids4 <- vapply(1:4, .ctp_snpId, character(1))
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
    ent <- twasWeightsRow(
        variantIds = c(ids4, bogus),
        weights = c(0.1, 0.2, 0.3, 0.4, 0.5),
        fits = fits
    )
    tw <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(ent),
        ldSketch = .ctp_makeHandle()
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
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fineMappingRow(tl$variant_id, NULL, tl))
    )
    aux <- pecotmr:::.ctwasGetFinemapAux(fmr, "Q1", "c1", "t1", "susie")
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
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fineMappingRow(tl$variant_id, NULL, tl))
    )
    aux <- pecotmr:::.ctwasGetFinemapAux(fmr, "Q1", "c1", "t1", "susie")
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
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fineMappingRow(tl$variant_id, NULL, tl))
    )
    aux <- pecotmr:::.ctwasGetFinemapAux(fmr, "Q1", "c1", "t1", "susie")
    expect_length(aux$csMembers, 0L)
    expect_length(aux$csPurity, 0L)
    expect_equal(aux$pip[["chr1:200:A:G"]], 0.8)
})

test_that(".ctwasGetFinemapAux: NULL input, no-match tuple, and empty topLoci all return NULL", {
    expect_null(pecotmr:::.ctwasGetFinemapAux(NULL, "Q1", "c1", "t1", "susie"))
    tl <- data.frame(
        variant_id = "chr1:100:A:G",
        pip = 0.5,
        stringsAsFactors = FALSE
    )
    fmr <- QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fineMappingRow(tl$variant_id, NULL, tl))
    )
    # No matching (study, context, trait, method) tuple -> NULL.
    expect_null(pecotmr:::.ctwasGetFinemapAux(fmr, "NOPE", "c1", "t1", "susie"))
    # Matching tuple but an empty topLoci -> NULL.
    fmrEmpty <- QtlFineMappingResult(
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(fineMappingRow(character(0), NULL, data.frame()))
    )
    expect_null(pecotmr:::.ctwasGetFinemapAux(
        fmrEmpty,
        "Q1",
        "c1",
        "t1",
        "susie"
    ))
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
        study = "G1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .ctp_makeHandle(),
        qcInfo = list(step1 = "ok")
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
        twasWeightsRow(
            variantIds = vapply(1:5, .ctp_snpId, character(1)),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
        )
    }
    tw <- TwasWeights(
        study = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("gA", "gX"),
        method = c("susie", "susie"),
        entry = list(mkE(), mkE()),
        traitPos = GenomicRanges::GRanges(
            c("chr1", "chr9"),
            IRanges::IRanges(c(100, 100), c(150, 150))
        ),
        ldSketch = .ctp_makeHandle()
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
        twasWeightsRow(
            variantIds = vapply(1:5, .ctp_snpId, character(1)),
            weights = c(0.1, 0.05, -0.2, 0.3, 0.0)
        )
    }
    twOff <- TwasWeights(
        study = "Q1",
        context = "c1",
        trait = "gX",
        method = "susie",
        entry = list(mkE()),
        traitPos = GenomicRanges::GRanges("chr9", IRanges::IRanges(100, 150)),
        ldSketch = .ctp_makeHandle()
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
            .ctp_makeTwasWeights()
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
    expect_error(pecotmr:::.ctwasResolveMethods(list()), "no method entries")
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
