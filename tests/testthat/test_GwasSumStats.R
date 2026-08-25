# Tests migrated from test_h2ClassesSumstats.R

# === Tests migrated from test_h2ClassesSumstats.R (GwasSumStats) ===

test_that("GwasSumStats(df) errors when required mcols are missing", {
    df <- data.frame(
        SNP = "rs1",
        CHR = "1",
        BP = 100,
        A1 = "A",
        stringsAsFactors = FALSE
    )
    expect_error(
        makeGwasSumStatsFromDf(df),
        "Missing required columns"
    )
})


test_that("GwasSumStats valid object passes with all required mcols", {
    set.seed(1)
    df <- data.frame(
        SNP = paste0("rs", 1:5),
        CHR = rep("1", 5),
        BP = 1:5,
        A1 = rep("A", 5),
        A2 = rep("G", 5),
        Z = rnorm(5),
        N = rep(1000, 5),
        stringsAsFactors = FALSE
    )
    obj <- makeGwasSumStatsFromDf(df)
    expect_true(methods::validObject(obj))
})


test_that("makeGwasSumStatsFromDf() constructor creates object from data.frame", {
    df <- make_test_sumstats_df(20)
    obj <- makeGwasSumStatsFromDf(df, traitName = "height", genome = "hg38")

    expect_s4_class(obj, "GwasSumStats")
    expect_equal(as.character(obj$study)[[1L]], "height")
    expect_equal(getGenome(obj), "hg38")
    expect_equal(length(getSumStats(obj)), 20)
})


test_that("makeGwasSumStatsFromDf() normalizes chr prefix", {
    df <- make_test_sumstats_df(5)
    # Input has CHR = "1" (no prefix)
    obj <- makeGwasSumStatsFromDf(df)
    chrs <- as.character(GenomicRanges::seqnames(getSumStats(obj)))
    expect_true(all(startsWith(chrs, "chr")))

    # Input already has "chr" prefix
    df2 <- df
    df2$CHR <- "chr1"
    obj2 <- makeGwasSumStatsFromDf(df2)
    chrs2 <- as.character(GenomicRanges::seqnames(getSumStats(obj2)))
    # Should not double-prefix
    expect_true(all(chrs2 == "chr1"))
    expect_false(any(grepl("^chrchr", chrs2)))
})


test_that("makeGwasSumStatsFromDf() errors on missing columns", {
    df <- data.frame(SNP = "rs1", CHR = "1", BP = 100)
    expect_error(makeGwasSumStatsFromDf(df), "Missing required columns")
})


test_that("makeGwasSumStatsFromDf() removes rows with NA in required columns", {
    df <- make_test_sumstats_df(10)
    df$Z[1] <- NA
    df$N[3] <- NA
    expect_message(
        obj <- makeGwasSumStatsFromDf(df),
        "Removed.*SNPs with missing"
    )
    expect_equal(length(getSumStats(obj)), 8)
})


test_that("getz() returns correct Z vector", {
    set.seed(99)
    df <- make_test_sumstats_df(5)
    obj <- makeGwasSumStatsFromDf(df)
    z <- getZ(obj)
    expect_type(z, "double")
    expect_equal(length(z), 5)
})


test_that("getn() returns correct N vector", {
    df <- make_test_sumstats_df(5)
    obj <- makeGwasSumStatsFromDf(df)
    n <- getN(obj)
    expect_equal(length(n), 5)
    expect_true(all(n == 10000))
})


test_that("getmaf() returns MAF when present, NULL when absent", {
    df <- make_test_sumstats_df(5)
    obj_no_maf <- makeGwasSumStatsFromDf(df)
    expect_null(getMaf(obj_no_maf))

    df$MAF <- runif(5, 0.01, 0.5)
    obj_with_maf <- makeGwasSumStatsFromDf(df)
    maf <- getMaf(obj_with_maf)
    expect_type(maf, "double")
    expect_equal(length(maf), 5)
})


test_that("nSnps() returns correct count", {
    df <- make_test_sumstats_df(30)
    obj <- makeGwasSumStatsFromDf(df)
    expect_equal(nSnps(obj), 30)
})


test_that("subsetchr() filters correctly", {
    df <- make_test_sumstats_df(10)
    df$CHR <- c(rep("1", 6), rep("2", 4))
    obj <- makeGwasSumStatsFromDf(df)

    chr1 <- subsetChr(obj, "1")
    expect_equal(nSnps(chr1), 6)

    # Also works with "chr" prefix
    chr2 <- subsetChr(obj, "chr2")
    expect_equal(nSnps(chr2), 4)
})


test_that("subsetChr() preserves study-level nCase/nControl/nSample scalars", {
    # Regression: the chromosome subset rebuilds the GwasSumStats and previously
    # dropped the optional per-study case/control counts + total N.
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
        study = "g1",
        entry = list(gr),
        genome = "hg19",
        ldSketch = .sh_makeGenotypeHandle(),
        nCase = 5000,
        nControl = 15000,
        nSample = 20000
    )

    sub <- subsetChr(obj, "chr1")
    expect_equal(nSnps(sub), 2)
    expect_equal(as.numeric(sub$nCase), 5000)
    expect_equal(as.numeric(sub$nControl), 15000)
    expect_equal(as.numeric(sub$nSample), 20000)
})


test_that("getvary() returns var_y and NULL cases", {
    df <- make_test_sumstats_df(5)

    obj_null <- makeGwasSumStatsFromDf(df, varY = NULL)
    expect_null(getVarY(obj_null))

    obj_vy <- makeGwasSumStatsFromDf(df, varY = 4.5)
    expect_equal(getVarY(obj_vy), 4.5)
})


test_that("as.data.frame.makeGwasSumStatsFromDf() round-trips", {
    df_in <- make_test_sumstats_df(15)
    obj <- makeGwasSumStatsFromDf(df_in)
    df_out <- as.data.frame(obj)

    expect_true(is.data.frame(df_out))
    expect_true(all(
        c("SNP", "CHR", "BP", "A1", "A2", "Z", "N") %in%
            names(df_out)
    ))
    expect_equal(nrow(df_out), 15)
    expect_equal(df_out$SNP, df_in$SNP)
    # BP should round-trip
    expect_equal(df_out$BP, as.integer(df_in$BP))
})


# =============================================================================
# AnnotationMatrix (h2Annotations.R)
# =============================================================================

# === Tests migrated from test_showMethods.R (GwasSumStats) ===

test_that("show.GwasSumStats prints nrow and genome build", {
    ss <- GwasSumStats(
        study = c("g1", "g2"),
        entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
        genome = "hg19",
        ldSketch = .sh_makeGenotypeHandle()
    )
    out <- capture.output(show(ss))
    expect_true(any(grepl("GwasSumStats: 2 studies, genome build hg19", out)))
    expect_true(any(grepl("LD sketch: gds @ /tmp/test.gds", out)))
})

# === Tests migrated from test_h2ClassesSumstats.R (showMethods) ===

# =============================================================================
# GwasSumStats() constructor validation + accessor error branches
# (.sh_makeQtlSumstatsGr / .sh_makeGenotypeHandle come from
# helper-showMethods.R)
# =============================================================================

test_that("GwasSumStats() errors when required args are missing", {
    expect_error(GwasSumStats(study = "g1"), "are all required")
})

test_that("GwasSumStats() errors when genome is not a single string", {
    expect_error(
        GwasSumStats(
            study = "g1",
            entry = list(.sh_makeQtlSumstatsGr()),
            genome = c("hg19", "hg38"),
            ldSketch = .sh_makeGenotypeHandle()
        ),
        "single character string"
    )
})

test_that("GwasSumStats() errors when entry is not a list", {
    expect_error(
        GwasSumStats(
            study = "g1",
            entry = "not_a_list",
            genome = "hg19",
            ldSketch = .sh_makeGenotypeHandle()
        ),
        "must be a list"
    )
})

test_that("GwasSumStats() errors when length(entry) != length(study)", {
    expect_error(
        GwasSumStats(
            study = c("g1", "g2"),
            entry = list(.sh_makeQtlSumstatsGr()),
            genome = "hg19",
            ldSketch = .sh_makeGenotypeHandle()
        ),
        "must equal length"
    )
})

test_that("GwasSumStats() errors when a per-study column has a bad length", {
    expect_error(
        GwasSumStats(
            study = c("g1", "g2"),
            entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
            genome = "hg19",
            ldSketch = .sh_makeGenotypeHandle(),
            nCase = c(1, 2, 3)
        ),
        "must have length 1 or length"
    )
})

test_that("GwasSumStats() attaches extra per-study columns via ...", {
    obj <- GwasSumStats(
        study = c("g1", "g2"),
        entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
        genome = "hg19",
        ldSketch = .sh_makeGenotypeHandle(),
        cohort = c("UKB", "FinnGen")
    )
    expect_equal(as.character(obj$cohort), c("UKB", "FinnGen"))
})

test_that("getSumStats() errors on an empty GwasSumStats", {
    empty <- GwasSumStats(
        study = character(0),
        entry = list(),
        genome = "hg19",
        ldSketch = .sh_makeGenotypeHandle(),
        varY = numeric(0)
    )
    expect_equal(nrow(empty), 0L)
    expect_error(getSumStats(empty), "has no rows")
})

test_that("getSumStats() on a multi-study GwasSumStats needs a study selector", {
    two <- GwasSumStats(
        study = c("g1", "g2"),
        entry = list(.sh_makeQtlSumstatsGr(), .sh_makeQtlSumstatsGr()),
        genome = "hg19",
        ldSketch = .sh_makeGenotypeHandle()
    )
    expect_error(getSumStats(two), "studies. Pass")
    expect_error(getSumStats(two, study = "ghost"), "Unknown study")
})

test_that("GwasSumStats: ldSketch is optional (NULL for LD-free workflows)", {
    ss <- GwasSumStats(
        study = "g1",
        entry = list(.sh_makeQtlSumstatsGr()),
        genome = "hg19"
    ) # ldSketch omitted -> NULL
    expect_null(getLdSketch(ss))
    expect_output(show(ss), "none \\(LD-free\\)")
})

test_that("GwasSumStats: a non-GenotypeHandle ldSketch is rejected", {
    expect_error(
        GwasSumStats(
            study = "g1",
            entry = list(.sh_makeQtlSumstatsGr()),
            genome = "hg19",
            ldSketch = "not_a_handle"
        ),
        "must be a genotype panel"
    )
})

# show() smoke test, moved here from test_showMethods.R so the test
# tree mirrors R/.
test_that("show(GwasSumStats) does not error", {
    ss <- makeGwasSumStatsFromDf(make_test_sumstats_df(10))
    expect_output(show(ss), "GwasSumStats")
})


# blockId + ldBlocks (§4.2). A genome-wide GWAS split by seqname gives one
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
.gss_ldBlocks <- function() {
    b <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 800L), c(500L, 1200L))
    )
    names(b) <- c("b1", "b2")
    b
}

test_that("GwasSumStats records blockId as the seqname without a manifest", {
    g <- .gss_variants(c("chr1", "chr1", "chr2"), c(100L, 900L, 300L))
    x <- GwasSumStats(study = "t1", entry = list(g), genome = "hg38")
    expect_equal(length(x), 2L)
    expect_equal(x$blockId, c("chr1", "chr2"))
    expect_true(is_in("blockId", colnames(x)))
})

test_that("GwasSumStats splits by LD block when a manifest is supplied", {
    g <- .gss_variants("chr1", c(100L, 250L, 900L))
    x <- GwasSumStats(
        study = "t1",
        entry = list(g),
        genome = "hg38",
        ldBlocks = .gss_ldBlocks()
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
        study = "t1",
        entry = list(g),
        genome = "hg38",
        varY = 2,
        ldBlocks = .gss_ldBlocks()
    )
    expect_equal(x$study, c("t1", "t1"))
    expect_equal(x$varY, c(2, 2))
})

test_that("GwasSumStats warns about variants outside every LD block", {
    g <- .gss_variants("chr1", c(100L, 5000L))
    expect_warning(
        x <- GwasSumStats(
            study = "t1",
            entry = list(g),
            genome = "hg38",
            ldBlocks = .gss_ldBlocks()
        ),
        "1 variant"
    )
    expect_equal(sum(lengths(x)), 1L)
})

test_that("GwasSumStats keeps blockId aligned across two studies", {
    x <- GwasSumStats(
        study = c("t1", "t2"),
        entry = list(
            .gss_variants("chr1", c(100L, 900L)),
            .gss_variants("chr1", 900L)
        ),
        genome = "hg38",
        ldBlocks = .gss_ldBlocks()
    )
    expect_equal(x$study, c("t1", "t1", "t2"))
    expect_equal(x$blockId, c("b1", "b2", "b2"))
})

test_that("GwasSumStats blockId survives subsetting", {
    x <- GwasSumStats(
        study = "t1",
        entry = list(.gss_variants("chr1", c(100L, 900L))),
        genome = "hg38",
        ldBlocks = .gss_ldBlocks()
    )
    expect_equal(x[2L]$blockId, "b2")
})

# ===========================================================================
# Block-key resolution: ldBlocks vs blockId
#
# Both arguments answer "how is this entry split into elements", so supplying
# both is ambiguous rather than redundant. Neither error branch had a test.
# ===========================================================================

.gss_entry <- function(n = 3L) {
    gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", n),
        ranges = IRanges::IRanges(seq(100L, by = 100L, length.out = n),
                                  width = 1L)
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        variant_id = str_c("chr1:", seq(100L, by = 100L, length.out = n),
                           ":A:G"),
        SNP = str_c("rs", seq_len(n)),
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = seq_len(n) / 2,
        N = rep(100L, n)
    )
    gr
}

test_that("supplying both ldBlocks and blockId is refused", {
    expect_error(
        GwasSumStats(
            study = "g1",
            entry = list(.gss_entry()),
            genome = "hg19",
            ldBlocks = GenomicRanges::GRanges(
                "chr1", IRanges::IRanges(1L, 10000L)
            ),
            blockId = "b1"
        ),
        "not both"
    )
})

test_that("blockId must carry one value per entry", {
    expect_error(
        GwasSumStats(
            study = c("g1", "g2"),
            entry = list(.gss_entry(), .gss_entry()),
            genome = "hg19",
            blockId = "only-one"
        ),
        "one value per `entry`"
    )
})

test_that("a supplied blockId survives onto the elements", {
    obj <- GwasSumStats(
        study = "g1",
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
    expect_equal(getGenome(out), "hg19")
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

    audit <- getQcInfo(out)$entryAudit
    expect_length(audit, nrow(out))
    # element i's audit still describes element i
    expect_equal(audit[[1L]]$block, "chr1_1_1000")
    expect_equal(audit[[2L]]$block, "chr1_1001_2000")
    expect_equal(getQcInfo(out)$options, getQcInfo(a)$options)
})


test_that("combineGwasSumStats() unions the per-block LD panels", {
    a <- makeGwasBlock("chr1_1_1000", 100L, n = 5L)
    b <- makeGwasBlock("chr1_1001_2000", 1100L, n = 3L)
    out <- combineGwasSumStats(a, b)

    # Keeping only the first block's panel would leave a two-block collection
    # whose LD reference covers one block.
    expect_equal(nrow(getLdSketch(out)), 8L)
    expect_equal(
        rownames(getLdSketch(out)),
        c(rownames(getLdSketch(a)), rownames(getLdSketch(b)))
    )
})


test_that("combineGwasSumStats() honours an explicit ldSketch", {
    a <- makeGwasBlock("chr1_1_1000", 100L)
    b <- makeGwasBlock("chr1_1001_2000", 1100L)
    out <- combineGwasSumStats(a, b, ldSketch = getLdSketch(a))
    expect_identical(getLdSketch(out), getLdSketch(a))
})


test_that("combineGwasSumStats() rejects collection-level disagreement", {
    a <- makeGwasBlock("chr1_1_1000", 100L)

    expect_error(
        combineGwasSumStats(a, makeGwasBlock("chr1_1001_2000", 1100L,
            genome = "hg38")),
        "share one genome build"
    )
    expect_error(
        combineGwasSumStats(a, makeGwasBlock("chr1_1001_2000", 1100L,
            qcOptions = list(mafCutoff = 0.05))),
        "different summaryStatsQc\\(\\) options"
    )
    noQc <- makeGwasBlock("chr1_1001_2000", 1100L)
    noQc@qcInfo <- list()
    expect_error(combineGwasSumStats(a, noQc), "carry no QC record")

    noLd <- makeGwasBlock("chr1_1001_2000", 1100L)
    noLd@ldSketch <- NULL
    expect_error(combineGwasSumStats(a, noLd), "carry no ldSketch")

    other <- makeGwasBlock("chr1_1001_2000", 1100L)
    other@ldSketch <- pecotmr:::.asLdSketch(
        .blockGenotypeHandle(seq(1100L, by = 100L, length.out = 5L),
            path = "/tmp/other.pgen")
    )
    expect_error(
        combineGwasSumStats(a, other),
        "different genotype panels"
    )
})


test_that("combineGwasSumStats() validates its inputs", {
    expect_error(combineGwasSumStats(list()), "nothing to combine")
    expect_error(
        combineGwasSumStats(makeGwasBlock("chr1_1_1000", 100L), 1L),
        "every input must be a GwasSumStats"
    )
})
