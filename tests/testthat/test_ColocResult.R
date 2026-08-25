# Tests for the ColocResult class and its four views.
#
# The per-variant layer is the point: colocPipeline used to discard
# coloc.bf_bf's $results, so every assertion below that touches SNP.PP.H4 is
# guarding the gap this class closes.

.cr_pairs <- function(
    qtlCs = 1L,
    gwasCs = 1L,
    blockId = "chr1_1_1000",
    pp4 = 0.6,
    trait = "g1"
) {
    n <- max(
        length(qtlCs),
        length(gwasCs),
        length(blockId),
        length(pp4),
        length(trait)
    )
    data.frame(
        study = "s1",
        context = "c1",
        trait = trait,
        method = "susie",
        gwasStudy = "G1",
        gwasMethod = "susie",
        blockId = blockId,
        qtlCs = as.integer(qtlCs),
        gwasCs = as.integer(gwasCs),
        nSnps = 3L,
        PP.H0.abf = (1 - pp4) / 4,
        PP.H1.abf = (1 - pp4) / 4,
        PP.H2.abf = (1 - pp4) / 4,
        PP.H3.abf = (1 - pp4) / 4,
        PP.H4.abf = pp4,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )[seq_len(n), , drop = FALSE]
}

.cr_variants <- function(n = 1L, pp = c(0.6, 0.3, 0.1)) {
    one <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A"),
        SNP.PP.H4 = pp
    )
    rep(list(one), n)
}

test_that("ColocResult: stores one element per pair with its variants", {
    x <- ColocResult(.cr_pairs(), .cr_variants())
    expect_s4_class(x, "ColocResult")
    expect_equal(nrow(x), 1L)
    expect_equal(lengths(x), 3L)
    expect_equal(x$PP.H4.abf, 0.6)
})

test_that("ColocResult: rejects an element with no per-variant layer", {
    # The whole point of the class is that SNP.PP.H4 survives; an element
    # without it must not construct.
    bad <- list(data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
        somethingElse = c(1, 2)
    ))
    expect_error(ColocResult(.cr_pairs(), bad), "SNP.PP.H4")
})

test_that("ColocResult: rejects a variants list that is not parallel", {
    expect_error(
        ColocResult(.cr_pairs(), .cr_variants(2L)),
        "parallel"
    )
})

test_that("ColocResult: an empty result still carries the column schema", {
    x <- ColocResult(.cr_pairs()[0, , drop = FALSE], list())
    expect_equal(nrow(x), 0L)
    expect_true(is_in("PP.H4.abf", colnames(x)))
    expect_equal(nrow(getColocPairs(x)), 0L)
})

test_that("ColocResult: variants keep their alleles, so ids round-trip", {
    # mcols must be APPENDED to the parsed A1/A2, not replace them: without
    # the alleles an element cannot render its own variant ids.
    x <- ColocResult(.cr_pairs(), .cr_variants())
    expect_equal(
        getColocVariants(x)$variant_id,
        c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    )
})

test_that("ColocResult: is reachable by range, not just by identity", {
    x <- ColocResult(.cr_pairs(), .cr_variants())
    win <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 250))
    expect_equal(sum(lengths(subsetRegion(x, win))), 1L)
    expect_equal(length(IRanges::subsetByOverlaps(x, win)), 1L)
})

test_that("getColocVariants: colocPp is the pair posterior times the variant", {
    x <- ColocResult(.cr_pairs(pp4 = 0.5), .cr_variants())
    v <- getColocVariants(x)
    expect_equal(nrow(v), 3L)
    expect_equal(v$colocPp, 0.5 * c(0.6, 0.3, 0.1))
})

test_that("getColocGenes: sums within a QTL CS, noisy-ORs across them", {
    # Two QTL credible sets, each tested against two blocks. Within a CS the
    # per-block results are mutually exclusive (sum); across CSs they are
    # independent signals (noisy-OR).
    pairs <- .cr_pairs(
        qtlCs = c(1L, 1L, 2L, 2L),
        blockId = c("b1", "b2", "b1", "b2"),
        pp4 = c(0.1, 0.2, 0.3, 0.4)
    )
    x <- ColocResult(pairs, .cr_variants(4L))
    genes <- getColocGenes(x)
    expect_equal(nrow(genes), 1L)
    expect_equal(genes$PP.H4, 1 - (1 - 0.3) * (1 - 0.7))
    expect_equal(genes$nQtlCs, 2L)
    expect_equal(genes$nPairs, 4L)
})

test_that("getColocGenes: clips a sum above 1 and says so", {
    pairs <- .cr_pairs(
        qtlCs = c(1L, 1L),
        blockId = c("b1", "b2"),
        pp4 = c(0.7, 0.6)
    )
    x <- ColocResult(pairs, .cr_variants(2L))
    expect_warning(genes <- getColocGenes(x), "competing")
    expect_equal(genes$PP.H4, 1)
})

test_that("getColocGenes: keeps distinct genes apart", {
    pairs <- .cr_pairs(trait = c("g1", "g2"), pp4 = c(0.5, 0.25))
    x <- ColocResult(pairs, .cr_variants(2L))
    genes <- getColocGenes(x)
    expect_equal(nrow(genes), 2L)
    expect_setequal(genes$trait, c("g1", "g2"))
})

test_that("getColocVariants(pooled): pools by the same rule as genes", {
    pairs <- .cr_pairs(
        qtlCs = c(1L, 2L),
        blockId = c("b1", "b1"),
        pp4 = c(0.5, 0.25)
    )
    x <- ColocResult(pairs, .cr_variants(2L))
    pooled <- getColocVariants(x, pooled = TRUE)
    expect_equal(nrow(pooled), 3L)
    lead <- pooled[pooled$variant_id == "chr1:100:A:G", ]
    expect_equal(lead$colocPp, 1 - (1 - 0.5 * 0.6) * (1 - 0.25 * 0.6))
})

test_that("getColocCredibleSets: takes the smallest set reaching coverage", {
    x <- ColocResult(.cr_pairs(), .cr_variants())
    expect_equal(getColocCredibleSets(x, coverage = 0.5)$csSize, 1L)
    expect_equal(getColocCredibleSets(x, coverage = 0.9)$csSize, 2L)
    expect_equal(getColocCredibleSets(x, coverage = 0.95)$csSize, 3L)
})

test_that("getColocCredibleSets: reports the lead variant and coverage", {
    x <- ColocResult(.cr_pairs(), .cr_variants())
    cs <- getColocCredibleSets(x, coverage = 0.9)
    expect_equal(cs$leadVariant, "chr1:100:A:G")
    expect_equal(cs$leadPp, 0.6)
    expect_equal(cs$csCoverage, 0.9)
})

test_that("getColocCredibleSets: minPp4 drops pairs before any LD work", {
    pairs <- .cr_pairs(qtlCs = c(1L, 2L), pp4 = c(0.8, 0.1))
    x <- ColocResult(pairs, .cr_variants(2L))
    expect_equal(nrow(getColocCredibleSets(x)), 2L)
    expect_equal(nrow(getColocCredibleSets(x, minPp4 = 0.5)), 1L)
    expect_equal(nrow(getColocCredibleSets(x, minPp4 = 0.99)), 0L)
})

test_that("getColocCredibleSets: requireMaxH4 keeps only H4-dominant pairs", {
    pairs <- .cr_pairs(qtlCs = c(1L, 2L), pp4 = c(0.8, 0.1))
    x <- ColocResult(pairs, .cr_variants(2L))
    cs <- getColocCredibleSets(x, requireMaxH4 = TRUE)
    expect_equal(nrow(cs), 1L)
    expect_equal(cs$qtlCs, 1L)
})

test_that("getColocCredibleSets: purity is NA, not passing, without LD", {
    # No LD reference means no evidence either way; reporting a passing value
    # would let a filtered view silently depend on whether a sketch happened
    # to be attached.
    x <- ColocResult(.cr_pairs(), .cr_variants())
    cs <- getColocCredibleSets(x)
    expect_true(is.na(cs$purity))
    expect_equal(nrow(getColocCredibleSets(x, minAbsCorr = 0.99)), 1L)
})

test_that("as.data.frame returns the pair view", {
    x <- ColocResult(.cr_pairs(), .cr_variants())
    df <- as.data.frame(x)
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 1L)
    expect_equal(df$PP.H4.abf, 0.6)
})

test_that(".crPivotColocResults: maps results column K to summary row K", {
    # Verified against coloc.bf_bf: SNP.PP.H4.rowK is positionally the K-th
    # summary row. Getting this backwards would silently attach the wrong
    # variant posteriors to a pair.
    results <- data.frame(
        snp = c("chr1:100:A:G", "chr1:200:C:T"),
        SNP.PP.H4.row1 = c(0.9, 0.1),
        SNP.PP.H4.row2 = c(0.2, 0.8)
    )
    out <- .crPivotColocResults(results, 2L)
    expect_equal(out[[1L]]$SNP.PP.H4, c(0.9, 0.1))
    expect_equal(out[[2L]]$SNP.PP.H4, c(0.2, 0.8))
})

test_that(".crPivotColocResults: empty results give empty variant tables", {
    out <- .crPivotColocResults(NULL, 2L)
    expect_length(out, 2L)
    expect_equal(nrow(out[[1L]]), 0L)
})

test_that(".crNoisyOr: matches the closed form", {
    expect_equal(.crNoisyOr(c(0.5, 0.5)), 0.75)
    expect_equal(.crNoisyOr(0.3), 0.3)
    expect_equal(.crNoisyOr(c(1, 0.5)), 1)
})

test_that(".crFirstAtLeast: absorbs floating-point error, not real shortfall", {
    expect_equal(.crFirstAtLeast(cumsum(c(0.6, 0.3, 0.1)), 0.9), 2L)
    expect_equal(.crFirstAtLeast(cumsum(c(0.6, 0.3, 0.1)), 0.95), 3L)
    # Never reaching the target returns everything, not an error.
    expect_equal(.crFirstAtLeast(cumsum(c(0.3, 0.2)), 0.9), 2L)
})

# ===========================================================================
# show
#
# Every other collection class has an expect_output(show(x), ...) test; these
# two did not, which left both show methods entirely unexecuted.
# ===========================================================================

test_that("show summarizes the pairs, studies and best PP.H4", {
    x <- ColocResult(.cr_pairs(), .cr_variants())
    expect_output(show(x), "ColocResult with 1 colocalized pair\\(s\\)")
    expect_output(show(x), "QTL studies")
    expect_output(show(x), "GWAS studies")
    expect_output(show(x), "variants")
    expect_output(show(x), "max PP.H4")
})

test_that("show on an empty ColocResult stops after the header", {
    x <- ColocResult(.cr_pairs()[0, , drop = FALSE], list())
    expect_output(show(x), "with 0 colocalized pair\\(s\\)")
    expect_false(any(grepl(
        "max PP.H4", capture.output(show(x)), fixed = TRUE
    )))
})

# ===========================================================================
# Credible-set purity
#
# Purity is the minimum absolute correlation among a set's members, so it
# needs an LD reference. Every test above builds a ColocResult without one,
# which left the whole extraction path dark.
# ===========================================================================

.cr_panel <- function() {
    readGenotypes(
        test_path("test_data", "test_variants"),
        format = "plink1"
    )
}

.cr_panelIds <- function(n = 3L) {
    panel <- .cr_panel()
    normalizeVariantId(
        .refVariantsFromSketch(.ldSketchHandle(panel))$variant_id
    )[seq_len(n)]
}

test_that("a singleton credible set is pure by definition", {
    # No pair to correlate; susie treats it as pure rather than unmeasurable.
    expect_equal(.crPurityOne(.cr_panelIds(1L), .cr_panel()), 1)
})

test_that("purity is the minimum absolute correlation among members", {
    ids <- .cr_panelIds(3L)
    pur <- .crPurityOne(ids, .cr_panel())
    R <- computeLd(.cr_panel(), snpIdx = 1:3)
    expect_equal(pur, min(abs(R[upper.tri(R)])), tolerance = 1e-8)
})

test_that("purity is NA when the panel cannot supply the members", {
    # Reported honestly rather than defaulted to a passing value: an
    # unmeasurable purity must not look like a pure set.
    expect_true(is.na(
        .crPurityOne(c("chr9:1:A:G", "chr9:2:C:T"), .cr_panel())
    ))
})

test_that("purity is NA for every set when there is no LD reference", {
    out <- tibble(csVariants = list(c("a", "b"), c("c", "d")))
    expect_equal(.crPuritiesFor(out, NULL), c(NA_real_, NA_real_))
})

# ===========================================================================
# Validity
# ===========================================================================

test_that("validity names missing pair columns", {
    bad <- ColocResult(.cr_pairs(), .cr_variants())
    mcols(bad)$gwasStudy <- NULL
    expect_error(methods::validObject(bad), "missing columns: gwasStudy")
})

test_that("validity names missing posterior columns", {
    bad <- ColocResult(.cr_pairs(), .cr_variants())
    mcols(bad)$PP.H4.abf <- NULL
    expect_error(
        methods::validObject(bad),
        "missing posterior columns: PP.H4.abf"
    )
})

test_that("validity names elements whose variants lack SNP.PP.H4", {
    # The expectation wraps the ASSIGNMENT: `[[<-` on the GRangesList
    # re-validates, so the object never survives long enough for a separate
    # validObject() call to see it.
    bad <- ColocResult(.cr_pairs(), .cr_variants())
    expect_error(
        mcols(bad[[1L]])$SNP.PP.H4 <- NULL,
        "have no SNP.PP.H4 column in their variant metadata"
    )
})
