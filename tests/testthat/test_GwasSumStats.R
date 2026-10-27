# Tests for R/gwasSumStats.R (the GwasSumStats class). makeGwasSumStatsFromDf(),
# below, is a test helper over GwasSumStats().

makeTestSumstatsDf <- function(n = 50) {
    set.seed(42)
    data.frame(
        SNP = paste0("rs", seq_len(n)),
        CHR = rep("1", n),
        BP = seq(1000, by = 100, length.out = n),
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = rnorm(n),
        N = rep(10000, n),
        stringsAsFactors = FALSE
    )
}

# Bridge helper: turn the legacy per-study data.frame shape into a
# single-row GwasSumStats collection using the new API. Keeps the bulk of
# the per-accessor tests below readable.
.testGenotypeHandle <- function() {
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

.dfToSumstatsGr <- function(df) {
    chrs <- as.character(df$CHR)
    if (!all(grepl("^chr", chrs))) {
        chrs <- paste0("chr", sub("^chr", "", chrs))
    }
    gr <- GenomicRanges::GRanges(
        seqnames = chrs,
        ranges = IRanges::IRanges(start = as.integer(df$BP), width = 1L)
    )
    mc <- df[, setdiff(colnames(df), c("CHR", "BP")), drop = FALSE]
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(mc)
    gr
}

makeGwasSumStatsFromDf <- function(
    df,
    traitName = "test",
    genome = "hg19",
    varY = NA_real_
) {
    if (is.null(varY)) {
        varY <- NA_real_
    }
    GwasSumStats(
        studyName = traitName,
        entry = list(.dfToSumstatsGr(df)),
        genome = genome,
        ldReference = .testGenotypeHandle(),
        varY = varY
    )
}

# A handle over a real (fake-path) variant set, so tests can build the
# per-block LD panels a block-parallel pipeline produces: same genotype source,
# each narrowed to its own variants.
.blockGenotypeHandle <- function(bp, path = "/tmp/test.pgen") {
    new(
        "GenotypeHandle",
        path = path,
        format = "plink2",
        snpInfo = data.frame(
            SNP = paste0("chr1:", bp, ":A:G"),
            CHR = rep("chr1", length(bp)),
            BP = as.integer(bp),
            A1 = rep("G", length(bp)),
            A2 = rep("A", length(bp)),
            fileIdx = seq_along(bp),
            stringsAsFactors = FALSE
        ),
        nSamples = 10L,
        sampleIds = paste0("s", seq_len(10)),
        pgenPtr = NULL
    )
}

# One LD block's worth of GWAS summary statistics: `n` variants starting at
# `start`, carrying a QC record and an LD panel narrowed to just those
# variants (what a per-block pipeline step writes out).
makeGwasBlock <- function(
    blockId,
    start,
    n = 5L,
    genome = "hg19",
    qcOptions = list(mafCutoff = 0.01),
    ldReference = .blockGenotypeHandle(bp)
) {
    bp <- seq(start, by = 100L, length.out = n)
    df <- data.frame(
        SNP = paste0("chr1:", bp, ":A:G"),
        CHR = rep("1", n),
        BP = bp,
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = seq_len(n) / 10,
        N = rep(1000L, n),
        stringsAsFactors = FALSE
    )
    obj <- GwasSumStats(
        studyName = "trait1",
        entry = list(.dfToSumstatsGr(df)),
        genome = genome,
        ldReference = ldReference,
        blockId = blockId,
        qcOptions = list(timestamp = NA_character_, options = qcOptions)
    )
    obj$qcAudit <- S4Vectors::SimpleList(list(list(
        variantsIn = n,
        block = blockId
    )))
    obj
}

test_that("GwasSumStats() records the study, genome build and variants", {
    df <- makeTestSumstatsDf(20)
    obj <- makeGwasSumStatsFromDf(df, traitName = "height", genome = "hg38")

    expect_s4_class(obj, "GwasSumStats")
    expect_true(methods::validObject(obj))
    # One chromosome, so one row, keyed on the study it was given.
    expect_identical(obj$study, "height")
    expect_equal(unname(GenomeInfoDb::genome(obj)), "hg38")
    expect_identical(sumStats(obj)$SNP, df$SNP)
})


test_that("GwasSumStats() reads an entry's chromosomes in UCSC style", {
    # Built here rather than through makeGwasSumStatsFromDf(), which adds
    # the "chr" prefix itself before the constructor sees the names.
    gr <- GenomicRanges::GRanges(
        c("1", "1", "CHR2"),
        IRanges::IRanges(c(100L, 200L, 300L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = c("rs1", "rs2", "rs3"),
        A1 = "A",
        A2 = "G",
        Z = c(1, 2, 3),
        N = 100L
    )
    obj <- GwasSumStats(studyName = "t1", entry = list(gr), genome = "hg19")
    expect_identical(
        as.character(GenomicRanges::seqnames(sumStats(obj))),
        c("chr1", "chr1", "chr2")
    )
    # The per-chromosome rows are keyed on the renamed chromosomes.
    expect_identical(obj$blockId, c("chr1", "chr2"))
    # Names already in UCSC style are kept, not prefixed a second time.
    again <- GwasSumStats(
        studyName = "t1",
        entry = list(sumStats(obj)),
        genome = "hg19"
    )
    expect_identical(GenomeInfoDb::seqlevels(again), c("chr1", "chr2"))
})


test_that("z() returns correct Z vector", {
    set.seed(99)
    df <- makeTestSumstatsDf(5)
    obj <- makeGwasSumStatsFromDf(df)
    z <- z(obj)
    expect_type(z, "double")
    expect_equal(length(z), 5)
})


test_that("sampleSize() returns correct N vector", {
    df <- makeTestSumstatsDf(5)
    obj <- makeGwasSumStatsFromDf(df)
    n <- sampleSize(obj)
    expect_equal(length(n), 5)
    expect_true(all(n == 10000))
})


test_that("maf() returns MAF when present, NULL when absent", {
    df <- makeTestSumstatsDf(5)
    obj_no_maf <- makeGwasSumStatsFromDf(df)
    expect_null(maf(obj_no_maf))

    df$MAF <- runif(5, 0.01, 0.5)
    obj_with_maf <- makeGwasSumStatsFromDf(df)
    maf <- maf(obj_with_maf)
    expect_type(maf, "double")
    expect_equal(length(maf), 5)
})


test_that("nSnps() returns correct count", {
    df <- makeTestSumstatsDf(30)
    obj <- makeGwasSumStatsFromDf(df)
    expect_equal(nSnps(obj), 30)
})


test_that("keepSeqlevels() narrows a GwasSumStats to a chromosome", {
    df <- makeTestSumstatsDf(10)
    df$CHR <- c(rep("1", 6), rep("2", 4))
    obj <- makeGwasSumStatsFromDf(df)

    chr1 <- GenomeInfoDb::keepSeqlevels(obj, "chr1", pruning.mode = "coarse")
    expect_s4_class(chr1, "GwasSumStats")
    expect_equal(nSnps(chr1), 6)
    chr2 <- GenomeInfoDb::keepSeqlevels(obj, "chr2", pruning.mode = "coarse")
    expect_equal(nSnps(chr2), 4)
})


test_that("keepSeqlevels() preserves study-level nCase/nControl/nSample", {
    # The per-study case/control counts and total N are row columns, so a
    # chromosome subset carries them with the rows it keeps.
    gr <- GenomicRanges::GRanges(
        c("chr1", "chr1", "chr2"),
        IRanges::IRanges(start = c(100L, 200L, 300L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = paste0("rs", 1:3),
        A1 = rep("A", 3),
        A2 = rep("G", 3),
        Z = c(1.0, -0.5, 2.0),
        N = rep(100L, 3)
    )
    obj <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .sh_makeGenotypeHandle(),
        nCase = 5000,
        nControl = 15000,
        nSample = 20000
    )

    sub <- GenomeInfoDb::keepSeqlevels(obj, "chr1", pruning.mode = "coarse")
    expect_equal(nSnps(sub), 2)
    expect_equal(as.numeric(sub$nCase), 5000)
    expect_equal(as.numeric(sub$nControl), 15000)
    expect_equal(as.numeric(sub$nSample), 20000)
})


test_that("varY() returns var_y and NULL cases", {
    df <- makeTestSumstatsDf(5)

    obj_null <- makeGwasSumStatsFromDf(df, varY = NULL)
    expect_null(varY(obj_null))

    obj_vy <- makeGwasSumStatsFromDf(df, varY = 4.5)
    expect_equal(varY(obj_vy), 4.5)
})


test_that("as.data.frame() is the standard coercion, one row per variant", {
    df_in <- makeTestSumstatsDf(15)
    obj <- makeGwasSumStatsFromDf(df_in)
    # The inherited GRangesList coercion: the element ranges plus their
    # stored mcols, under the stored names.
    std <- as.data.frame(obj)
    expect_true(all(c("group", "seqnames", "start", "Z") %in% names(std)))
    expect_equal(nrow(std), 15)
    # The standardized layout the pipelines use is built internally.
    df_out <- pecotmr:::.entryToSumstatDf(sumStats(obj))
    expect_true(is.data.frame(df_out))
    expect_true(all(
        c("variant_id", "chrom", "pos", "A1", "A2", "z", "N") %in%
            names(df_out)
    ))
    expect_equal(nrow(df_out), 15)
    expect_equal(df_out$variant_id, df_in$SNP)
    # position should round-trip
    expect_equal(df_out$pos, as.integer(df_in$BP))
})


# =============================================================================
# show()
# =============================================================================

test_that("show.GwasSumStats prints nrow and genome build", {
    ss <- GwasSumStats(
        studyName = c("g1", "g2"),
        entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
        genome = "hg19",
        ldReference = .sh_makeGenotypeHandle()
    )
    out <- capture.output(show(ss))
    expect_true(any(grepl(
        "GwasSumStats: 2 studies in 2 rows, genome build hg19",
        out
    )))
    expect_true(any(grepl("LD reference: plink2 @ /tmp/test", out)))
})

test_that("show.GwasSumStats counts a study split across chromosomes once", {
    gr <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(c(100L, 100L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = c("chr1:100:A:G", "chr2:100:A:G"),
        A1 = "A",
        A2 = "G",
        Z = c(1, 2),
        N = 100L
    )
    ss <- GwasSumStats(studyName = "g1", entry = list(gr), genome = "hg19")
    expect_equal(nrow(ss), 2L)
    out <- capture.output(show(ss))
    # One header line, not one per chromosome the study spans.
    expect_identical(
        grep("^GwasSumStats:", out, value = TRUE),
        "GwasSumStats: 1 studies in 2 rows, genome build hg19"
    )
})

# =============================================================================
# GwasSumStats() constructor validation + accessor error branches
# (.sh_makeQtlSumstatsGr / .sh_makeGenotypeHandle come from
# helper-showMethods.R)
# =============================================================================

test_that("GwasSumStats() errors when required args are missing", {
    expect_error(GwasSumStats(studyName = "g1"), "are all required")
})

test_that("GwasSumStats() errors when genome is not a single string", {
    expect_error(
        GwasSumStats(
            studyName = "g1",
            entry = list(.sh_makeQtlSumstatsGr()),
            genome = c("hg19", "hg38"),
            ldReference = .sh_makeGenotypeHandle()
        ),
        "single character string"
    )
})

test_that("GwasSumStats() errors when entry is not a list", {
    expect_error(
        GwasSumStats(
            studyName = "g1",
            entry = "not_a_list",
            genome = "hg19",
            ldReference = .sh_makeGenotypeHandle()
        ),
        "must be a list"
    )
})

test_that("GwasSumStats() errors when length(entry) != length(study)", {
    expect_error(
        GwasSumStats(
            studyName = c("g1", "g2"),
            entry = list(.sh_makeQtlSumstatsGr()),
            genome = "hg19",
            ldReference = .sh_makeGenotypeHandle()
        ),
        "must equal length"
    )
})

test_that("GwasSumStats() errors when a per-study column has a bad length", {
    expect_error(
        GwasSumStats(
            studyName = c("g1", "g2"),
            entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
            genome = "hg19",
            ldReference = .sh_makeGenotypeHandle(),
            nCase = c(1, 2, 3)
        ),
        "must have length 1 or length"
    )
})

test_that("GwasSumStats() attaches extra per-study columns via ...", {
    obj <- GwasSumStats(
        studyName = c("g1", "g2"),
        entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
        genome = "hg19",
        ldReference = .sh_makeGenotypeHandle(),
        extraCols = list(cohort = c("UKB", "FinnGen"))
    )
    expect_equal(as.character(obj$cohort), c("UKB", "FinnGen"))
    # Not variadic: the base class's GRangesList(...) means ELEMENTS,
    # so a bare `cohort =` here would read as adding an entry.
    expect_error(
        GwasSumStats(
            studyName = "g1",
            entry = list(),
            genome = "hg19",
            cohort = "UKB"
        ),
        "unused argument"
    )
})

test_that("sumStats() errors on an empty GwasSumStats", {
    empty <- GwasSumStats(
        studyName = character(0),
        entry = list(),
        genome = "hg19",
        ldReference = .sh_makeGenotypeHandle(),
        varY = numeric(0)
    )
    expect_equal(nrow(empty), 0L)
    expect_error(sumStats(empty), "has no rows")
})

test_that("sumStats() on a multi-study GwasSumStats needs a study selector", {
    two <- GwasSumStats(
        studyName = c("g1", "g2"),
        entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
        genome = "hg19",
        ldReference = .sh_makeGenotypeHandle()
    )
    expect_error(sumStats(two), "studies. Pass")
    expect_error(sumStats(two, studyName = "ghost"), "Unknown study")
})

test_that("GwasSumStats: ldReference is optional (NULL when LD-free)", {
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(.sh_makeQtlSumstatsGr()),
        genome = "hg19"
    ) # ldReference omitted -> LD-free
    expect_identical(ldReference(ss), list())
    expect_true(is.na(ss$ldReference))
    expect_output(show(ss), "none \\(LD-free\\)")
})

test_that("GwasSumStats: a non-GenotypeHandle ldReference is rejected", {
    expect_error(
        GwasSumStats(
            studyName = "g1",
            entry = list(.sh_makeQtlSumstatsGr()),
            genome = "hg19",
            ldReference = "not_a_handle"
        ),
        "must be a genotype panel"
    )
})


# blockId + region (§4.2). A genome-wide GWAS split by seqname gives one
# element per chromosome, which is too coarse for cTWAS: its EM estimates
# parameters across LD blocks. `blockId` is always present so downstream code
# can key regions without asking how the collection was built.

# @noRd
.gss_variants <- function(chrom, pos) {
    g <- GenomicRanges::GRanges(chrom, IRanges::IRanges(pos, width = 1L))
    mcols(g) <- S4Vectors::DataFrame(
        SNP = str_c("rs", seq_along(pos)),
        A1 = "A",
        A2 = "G",
        Z = seq_along(pos) / 10,
        N = 100L
    )
    g
}

# @noRd
.gss_blocks <- function() {
    b <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 800L), c(500L, 1200L))
    )
    names(b) <- c("b1", "b2")
    b
}

test_that("GwasSumStats records blockId as the seqname without a manifest", {
    g <- .gss_variants(c("chr1", "chr1", "chr2"), c(100L, 900L, 300L))
    x <- GwasSumStats(studyName = "t1", entry = list(g), genome = "hg38")
    expect_equal(length(x), 2L)
    expect_equal(x$blockId, c("chr1", "chr2"))
    expect_true(is_in("blockId", colnames(x)))
})

test_that("GwasSumStats splits by LD block when a manifest is supplied", {
    g <- .gss_variants("chr1", c(100L, 250L, 900L))
    x <- GwasSumStats(
        studyName = "t1",
        entry = list(g),
        genome = "hg38",
        region = .gss_blocks()
    )
    # Same chromosome, so the seqname split would have left ONE element; the
    # manifest is what produces per-block granularity.
    expect_equal(length(x), 2L)
    expect_equal(x$blockId, c("b1", "b2"))
    expect_equal(unname(lengths(x)), c(2L, 1L))
    expect_true(validObject(x, test = TRUE) == TRUE)
})

test_that("GwasSumStats replicates the study row against each block", {
    g <- .gss_variants("chr1", c(100L, 900L))
    x <- GwasSumStats(
        studyName = "t1",
        entry = list(g),
        genome = "hg38",
        varY = 2,
        region = .gss_blocks()
    )
    expect_equal(x$study, c("t1", "t1"))
    expect_equal(x$varY, c(2, 2))
})

test_that("GwasSumStats warns about variants outside every LD block", {
    g <- .gss_variants("chr1", c(100L, 5000L))
    expect_warning(
        x <- GwasSumStats(
            studyName = "t1",
            entry = list(g),
            genome = "hg38",
            region = .gss_blocks()
        ),
        "1 variant"
    )
    expect_equal(sum(lengths(x)), 1L)
})

test_that("GwasSumStats keeps blockId aligned across two studies", {
    x <- GwasSumStats(
        studyName = c("t1", "t2"),
        entry = list(
            .gss_variants("chr1", c(100L, 900L)),
            .gss_variants("chr1", 900L)
        ),
        genome = "hg38",
        region = .gss_blocks()
    )
    expect_equal(x$study, c("t1", "t1", "t2"))
    expect_equal(x$blockId, c("b1", "b2", "b2"))
})

test_that("GwasSumStats blockId survives subsetting", {
    x <- GwasSumStats(
        studyName = "t1",
        entry = list(.gss_variants("chr1", c(100L, 900L))),
        genome = "hg38",
        region = .gss_blocks()
    )
    expect_equal(x[2L]$blockId, "b2")
})

# ===========================================================================
# Block-key resolution: region vs blockId
#
# Both arguments answer "how is this entry split into elements", so supplying
# both is ambiguous rather than redundant. Neither error branch had a test.
# ===========================================================================

.gss_entry <- function(n = 3L) {
    gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", n),
        ranges = IRanges::IRanges(
            seq(100L, by = 100L, length.out = n),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        variant_id = str_c(
            "chr1:",
            seq(100L, by = 100L, length.out = n),
            ":A:G"
        ),
        SNP = str_c("rs", seq_len(n)),
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = seq_len(n) / 2,
        N = rep(100L, n)
    )
    gr
}

test_that("supplying both region and blockId is refused", {
    expect_error(
        GwasSumStats(
            studyName = "g1",
            entry = list(.gss_entry()),
            genome = "hg19",
            region = GenomicRanges::GRanges(
                "chr1",
                IRanges::IRanges(1L, 10000L)
            ),
            blockId = "b1"
        ),
        "not both"
    )
})

test_that("blockId must carry one value per entry", {
    expect_error(
        GwasSumStats(
            studyName = c("g1", "g2"),
            entry = list(.gss_entry(), .gss_entry()),
            genome = "hg19",
            blockId = "only-one"
        ),
        "one value per `entry`"
    )
})

test_that("a supplied blockId survives onto the elements", {
    obj <- GwasSumStats(
        studyName = "g1",
        entry = list(.gss_entry()),
        genome = "hg19",
        blockId = "chr1:1-10000"
    )
    expect_equal(unique(as.character(mcols(obj)$blockId)), "chr1:1-10000")
})


# === combineGwasSumStats() ===

test_that("combineGwasSumStats() row-binds per-block pieces in order", {
    a <- makeGwasBlock("chr1_1_1000", 100L, n = 5L)
    b <- makeGwasBlock("chr1_1001_2000", 1100L, n = 3L)
    out <- combineGwasSumStats(a, b)

    expect_s4_class(out, "GwasSumStats")
    expect_equal(nrow(out), nrow(a) + nrow(b))
    expect_equal(as.character(out$blockId), c("chr1_1_1000", "chr1_1001_2000"))
    expect_identical(as.list(out), c(as.list(a), as.list(b)))
    expect_equal(unname(GenomeInfoDb::genome(out)), "hg19")
})


test_that("combineGwasSumStats() accepts a list and passes a lone input through", {
    a <- makeGwasBlock("chr1_1_1000", 100L)
    b <- makeGwasBlock("chr1_1001_2000", 1100L)
    expect_identical(
        combineGwasSumStats(list(a, b)),
        combineGwasSumStats(a, b)
    )
    expect_identical(combineGwasSumStats(a), a)
})


test_that("combineGwasSumStats() concatenates the per-element QC audit", {
    a <- makeGwasBlock("chr1_1_1000", 100L, n = 5L)
    b <- makeGwasBlock("chr1_1001_2000", 1100L, n = 3L)
    out <- combineGwasSumStats(a, b)

    audit <- out$qcAudit
    expect_length(audit, nrow(out))
    # element i's audit still describes element i
    expect_equal(audit[[1L]]$block, "chr1_1_1000")
    expect_equal(audit[[2L]]$block, "chr1_1001_2000")
    expect_equal(qcOptions(out)$options, qcOptions(a)$options)
})


test_that("combineGwasSumStats() unions the per-block LD panels", {
    a <- makeGwasBlock("chr1_1_1000", 100L, n = 5L)
    b <- makeGwasBlock("chr1_1001_2000", 1100L, n = 3L)
    out <- combineGwasSumStats(a, b)

    # Each block's panel was narrowed to its own variants, but one panel's
    # files back both, so the bound collection holds them once, unioned.
    refs <- ldReference(out)
    expect_length(refs, 1L)
    expect_equal(nrow(refs[[1L]]), 8L)
    expect_equal(
        rownames(refs[[1L]]),
        c(rownames(ldReference(a)[[1L]]), rownames(ldReference(b)[[1L]]))
    )
    expect_equal(unique(out$ldReference), names(refs))
})


test_that("combineGwasSumStats() honours an explicit ldReference", {
    a <- makeGwasBlock("chr1_1_1000", 100L)
    b <- makeGwasBlock("chr1_1001_2000", 1100L)
    out <- combineGwasSumStats(a, b, ldReference = ldReference(a))
    expect_identical(ldReference(out), ldReference(a))
})


test_that("combineGwasSumStats() rejects collection-level disagreement", {
    a <- makeGwasBlock("chr1_1_1000", 100L)

    expect_error(
        combineGwasSumStats(
            a,
            makeGwasBlock("chr1_1001_2000", 1100L, genome = "hg38")
        ),
        "share one genome build"
    )
    expect_error(
        combineGwasSumStats(
            a,
            makeGwasBlock(
                "chr1_1001_2000",
                1100L,
                qcOptions = list(mafCutoff = 0.05)
            )
        ),
        "different summaryStatsQc\\(\\) options"
    )
    noQc <- makeGwasBlock("chr1_1001_2000", 1100L)
    noQc@qcOptions <- list()
    expect_error(combineGwasSumStats(a, noQc), "carry no QC record")
})


test_that("combineGwasSumStats() keeps each row on its own LD reference", {
    a <- makeGwasBlock("chr1_1_1000", 100L)
    other <- makeGwasBlock(
        "chr1_1001_2000",
        1100L,
        ldReference = .blockGenotypeHandle(
            seq(1100L, by = 100L, length.out = 5L),
            path = "/tmp/other.pgen"
        )
    )
    out <- combineGwasSumStats(a, other)
    # Panels read from different files cannot be unioned, so both are kept.
    expect_length(ldReference(out), 2L)
    expect_identical(
        pecotmr:::.rowLdRefs(out),
        list(ldReference(a)[[1L]], ldReference(other)[[1L]])
    )
    # An LD-free part's rows stay LD-free beside them.
    noLd <- makeGwasBlock("chr1_2001_3000", 2100L, ldReference = NULL)
    mixed <- combineGwasSumStats(a, noLd)
    expect_equal(is.na(mixed$ldReference), c(FALSE, TRUE))
    expect_length(ldReference(mixed), 1L)
})


test_that("combineGwasSumStats() validates its inputs", {
    expect_error(combineGwasSumStats(list()), "nothing to combine")
    expect_error(
        combineGwasSumStats(makeGwasBlock("chr1_1_1000", 100L), 1L),
        "every input must be a GwasSumStats"
    )
})

test_that("sumStats(GwasSumStats) refuses the QtlSumStats selectors", {
    # The shared SumStatsBase accessors pass the union of both classes'
    # selectors, so a GWAS sees context / trait, and sumStats() also takes
    # annotateSignificance. It names them in order to reject them, rather
    # than absorbing them in `...`.
    gss <- GwasSumStats(
        studyName = "g1",
        entry = list(.gss_entry()),
        genome = "hg19"
    )
    expect_error(sumStats(gss, context = "brain"), "QtlSumStats selector")
    expect_error(sumStats(gss, trait = "G1"), "QtlSumStats selector")
    expect_error(
        sumStats(gss, annotateSignificance = "qvalue"),
        "QtlSumStats selector"
    )
    # and the same through an accessor that forwards them
    expect_error(z(gss, context = "brain"), "QtlSumStats selector")
    # the supported selectors still work
    expect_s4_class(sumStats(gss), "GRanges")
    expect_type(z(gss), "double")
})

test_that("varY reads one value for a study split across rows", {
    gr <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(c(100L, 100L), width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = c("chr1:100:A:G", "chr2:100:A:G"),
        A1 = "A",
        A2 = "G",
        Z = c(1, 2),
        N = 100L
    )
    ss <- GwasSumStats(
        studyName = "g1",
        entry = list(gr),
        genome = "hg19",
        varY = 4
    )
    expect_equal(nrow(ss), 2L)
    expect_equal(varY(ss), 4)
    # Rows of one study that disagree have no single answer.
    split <- GwasSumStats(
        studyName = c("g1", "g1"),
        entry = list(gr[1L], gr[2L]),
        genome = "hg19",
        varY = c(1, 2),
        blockId = c("a", "b")
    )
    expect_error(varY(split), "the 2 rows of study 'g1' disagree on varY")
})

test_that("GwasSumStats records each row's window and its origin", {
    data(gwasSumStatsS4Example)
    gr <- sumStats(gwasSumStatsS4Example)
    cut <- GenomicRanges::start(gr)[[100L]]
    # A `region` gives each window its interval, however it is keyed.
    blocks <- GenomicRanges::GRanges(
        "chr22",
        IRanges::IRanges(c(1L, cut + 1L), c(cut, 60000000L))
    )
    names(blocks) <- c("b1", "b2")
    x <- GwasSumStats(
        studyName = "s",
        entry = list(gr),
        genome = "hg19",
        region = blocks
    )
    expect_identical(x$blockId, c("b1", "b2"))
    expect_identical(x$regionGiven, c(TRUE, TRUE))
    expect_identical(
        as.character(x$region),
        c(str_c("chr22:1-", cut), str_c("chr22:", cut + 1L, "-60000000"))
    )
    # Read from a BED file, a window's start is 0-based and the fourth column
    # names it. The variants past its one window are dropped, with a warning.
    bedPath <- withr::local_tempfile(fileext = ".bed")
    writeLines(str_c("chr22\t0\t", cut, "\tb1"), bedPath)
    expect_warning(
        bed <- GwasSumStats(
            studyName = "s",
            entry = list(gr),
            genome = "hg19",
            region = bedPath
        ),
        "100 variant(s) lie outside every `region` window",
        fixed = TRUE
    )
    expect_identical(bed$blockId, "b1")
    expect_identical(as.character(bed$region), str_c("chr22:1-", cut))
    # Given in R, a window is a GRanges: a table is refused, not read.
    expect_error(
        GwasSumStats(
            studyName = "s",
            entry = list(gr),
            genome = "hg19",
            region = data.frame(chrom = "chr22", start = 1L, end = cut)
        ),
        "must be a GRanges or the path of a BED file"
    )
    # A blockId is a name and gives no window, however it is spelled: each
    # row's window is its variants' span.
    y <- GwasSumStats(
        studyName = c("s", "s"),
        entry = list(gr[1:100], gr[101:200]),
        genome = "hg19",
        blockId = c("chr22_1_15700000", "b2")
    )
    expect_identical(y$blockId, c("chr22_1_15700000", "b2"))
    expect_identical(y$regionGiven, c(FALSE, FALSE))
    expect_identical(
        as.character(y$region),
        c(as.character(range(gr[1:100])), as.character(range(gr[101:200])))
    )
    # One set of windows per study, or none for a study split by chromosome.
    per <- suppressWarnings(GwasSumStats(
        studyName = c("s", "t"),
        entry = list(gr, gr),
        genome = "hg19",
        region = list(blocks[1L], NULL)
    ))
    expect_identical(per$study, c("s", "t"))
    expect_identical(per$blockId, c("b1", "chr22"))
    expect_identical(per$regionGiven, c(TRUE, FALSE))
    expect_error(
        GwasSumStats(
            studyName = "s",
            entry = list(gr),
            genome = "hg19",
            region = list(blocks, blocks)
        ),
        "one element per `entry`"
    )
    # Split by chromosome, a row's window is its variants' span.
    z <- GwasSumStats(studyName = "s", entry = list(gr), genome = "hg19")
    expect_false(z$regionGiven)
    expect_identical(as.character(z$region), "chr22:14560203-15826991")
})

test_that("combineGwasSumStats() carries each row's block range through", {
    data(gwasSumStatsS4Example)
    gr <- sumStats(gwasSumStatsS4Example)
    given <- GwasSumStats(
        studyName = "s",
        entry = list(gr),
        genome = "hg19",
        region = GenomicRanges::GRanges("chr22:1-20000000")
    )
    computed <- GwasSumStats(studyName = "t", entry = list(gr), genome = "hg19")
    both <- combineGwasSumStats(given, computed)
    expect_identical(both$regionGiven, c(TRUE, FALSE))
    expect_identical(
        as.character(both$region),
        c("chr22:1-20000000", "chr22:14560203-15826991")
    )
})

# === metadata, qcOptions and a tabular extraCols ===

test_that("GwasSumStats() takes extra columns as a DataFrame or data.frame", {
    build <- function(extra) {
        GwasSumStats(
            studyName = c("g1", "g2"),
            entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
            genome = "hg19",
            ldReference = .sh_makeGenotypeHandle(),
            extraCols = extra
        )
    }
    cohorts <- c("UKB", "FinnGen")
    fromDf <- build(S4Vectors::DataFrame(cohort = cohorts, wave = 1:2))
    expect_equal(fromDf$cohort, cohorts)
    expect_equal(fromDf$wave, 1:2)
    expect_equal(build(data.frame(cohort = cohorts))$cohort, cohorts)
})

test_that("GwasSumStats() keeps the caller's metadata and QC options", {
    obj <- GwasSumStats(
        studyName = "g1",
        entry = list(.sh_makeQtlSumstatsGr()),
        genome = "hg19",
        ldReference = .sh_makeGenotypeHandle(),
        qcOptions = list(step1 = "ok"),
        metadata = list(source = "unit test")
    )
    expect_equal(S4Vectors::metadata(obj), list(source = "unit test"))
    expect_equal(qcOptions(obj), list(step1 = "ok"))
})

# ===========================================================================
# nCase and nControl
# ===========================================================================

test_that("GwasSumStats: nCase/nControl are optional columns (absent by default)", {
    gr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = 5L),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = sprintf("chr1:%d:A:G", 100L * (1:5)),
        A1 = "A",
        A2 = "G",
        Z = rnorm(5),
        N = rep(1000L, 5)
    )
    base <- list(
        study = "G1",
        entry = list(gr),
        genome = "hg19",
        ldReference = .cbp_makeHandle(),
        qcOptions = list(ok = 1)
    )
    g0 <- exec(GwasSumStats, !!!base)
    expect_false(any(
        c("nCase", "nControl") %in% colnames(S4Vectors::mcols(g0))
    ))
    g1 <- exec(GwasSumStats, !!!c(base, list(nCase = 500, nControl = 1500)))
    expect_true(all(c("nCase", "nControl") %in% colnames(S4Vectors::mcols(g1))))
    expect_equal(g1$nCase, 500)
    expect_equal(g1$nControl, 1500)
})
