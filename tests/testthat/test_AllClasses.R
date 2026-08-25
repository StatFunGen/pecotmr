context("AllClasses (virtual base classes)")

# Most slots / accessors on the concrete subclasses are exercised in their
# own test files; these tests target the *base-class* behaviors that the
# concrete subclasses inherit without overriding (getStudy on SumStatsBase,
# getQcDiagnostics body branches, and the zero-row adjustPips short-circuit
# on FineMappingResultBase).

# ===========================================================================
# Helpers
# ===========================================================================

.alc_makeHandle <- function(snp_n = 3L) {
    new(
        "GenotypeHandle",
        path = "/tmp/test.gds",
        format = "gds",
        snpInfo = data.frame(
            SNP = paste0("rs", seq_len(snp_n)),
            CHR = rep("1", snp_n),
            BP = seq(100L, by = 100L, length.out = snp_n),
            A1 = rep("A", snp_n),
            A2 = rep("G", snp_n),
            stringsAsFactors = FALSE
        ),
        nSamples = 10L,
        sampleIds = paste0("s", seq_len(10L)),
        pgenPtr = NULL
    )
}

.alc_makeGr <- function(n = 3) {
    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(
            start = seq(100L, by = 100L, length.out = n),
            width = 1L
        )
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        SNP = paste0("rs", seq_len(n)),
        A1 = rep("A", n),
        A2 = rep("G", n),
        Z = rnorm(n),
        N = rep(1000L, n)
    )
    gr
}

.alc_makeGwasSumStats <- function(qcInfo = list()) {
    GwasSumStats(
        study = "g1",
        entry = list(.alc_makeGr()),
        genome = "hg19",
        ldSketch = .alc_makeHandle(),
        qcInfo = qcInfo
    )
}

.alc_makeFmEntry <- function(n = 3) {
    tl <- data.frame(
        variant_id = paste0("chr1:", 100 * seq_len(n), ":A:G"),
        chrom = rep("1", n),
        pos = as.integer(100 * seq_len(n)),
        A1 = rep("G", n),
        A2 = rep("A", n),
        N = rep(1000, n),
        MAF = rep(0.1, n),
        marginal_beta = rep(0.1, n),
        marginal_se = rep(0.05, n),
        marginal_z = rep(2.0, n),
        marginal_p = rep(0.05, n),
        pip = seq(0.9, by = -0.1, length.out = n),
        posterior_mean = rep(0.05, n),
        posterior_sd = rep(0.02, n),
        stringsAsFactors = FALSE
    )
    fineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl
    )
}

# ===========================================================================
# SumStatsBase: getStudy (inherited by QtlSumStats / GwasSumStats)
# ===========================================================================

test_that("SumStatsBase: getStudy on a GwasSumStats returns unique study names", {
    ss <- .alc_makeGwasSumStats()
    expect_equal(getStudy(ss), "g1")
})

# ===========================================================================
# SumStatsBase: getQcDiagnostics — every branch
# ===========================================================================

test_that("SumStatsBase: getQcDiagnostics returns NULL on empty qcInfo", {
    ss <- .alc_makeGwasSumStats() # qcInfo = list() by default
    expect_null(getQcDiagnostics(ss))
})

test_that("SumStatsBase: getQcDiagnostics returns NULL when entryAudit slot is absent", {
    # qcInfo has steps but no entryAudit -> nothing to return.
    ss <- .alc_makeGwasSumStats(qcInfo = list(step1 = "ok"))
    expect_null(getQcDiagnostics(ss))
})

test_that("SumStatsBase: getQcDiagnostics returns the per-entry diagnostics by index", {
    diag1 <- data.frame(SNP = "rs1", outlier = FALSE, stringsAsFactors = FALSE)
    diag2 <- data.frame(SNP = "rs2", outlier = TRUE, stringsAsFactors = FALSE)
    qc <- list(
        entryAudit = list(
            list(ldMismatchDiagnostics = diag1),
            list(ldMismatchDiagnostics = diag2)
        )
    )
    ss <- .alc_makeGwasSumStats(qcInfo = qc)
    expect_identical(getQcDiagnostics(ss, entry = 1L), diag1)
    expect_identical(getQcDiagnostics(ss, entry = 2L), diag2)
})

test_that("SumStatsBase: getQcDiagnostics(entry = NULL) returns the populated entries only", {
    diag1 <- data.frame(SNP = "rs1", outlier = FALSE, stringsAsFactors = FALSE)
    # Entry 2's audit has no ldMismatchDiagnostics field; should be filtered.
    qc <- list(
        entryAudit = list(
            list(ldMismatchDiagnostics = diag1),
            list(other = "no diagnostics here")
        )
    )
    ss <- .alc_makeGwasSumStats(qcInfo = qc)
    out <- getQcDiagnostics(ss, entry = NULL)
    expect_type(out, "list")
    expect_equal(length(out), 1L)
    expect_named(out, "1")
    expect_identical(out[["1"]], diag1)
})

test_that("SumStatsBase: getQcDiagnostics(entry = NULL) returns NULL when no entry has diagnostics", {
    qc <- list(entryAudit = list(list(other = 1), list(other = 2)))
    ss <- .alc_makeGwasSumStats(qcInfo = qc)
    expect_null(getQcDiagnostics(ss, entry = NULL))
})

test_that("SumStatsBase: getQcDiagnostics errors on out-of-range entry", {
    qc <- list(
        entryAudit = list(list(ldMismatchDiagnostics = data.frame(z = 1)))
    )
    ss <- .alc_makeGwasSumStats(qcInfo = qc)
    expect_error(getQcDiagnostics(ss, entry = 0L), "must be a single integer")
    expect_error(getQcDiagnostics(ss, entry = 99L), "must be a single integer")
    expect_error(
        getQcDiagnostics(ss, entry = c(1L, 2L)),
        "must be a single integer"
    )
    expect_error(
        getQcDiagnostics(ss, entry = "first"),
        "must be a single integer"
    )
})

# ===========================================================================
# FineMappingResultBase: adjustPips zero-row short-circuit
# ===========================================================================

test_that("FineMappingResultBase: adjustPips on a zero-row collection returns the input unchanged", {
    e <- .alc_makeFmEntry(3)
    res <- GwasFineMappingResult(
        study = "g1",
        method = "susie",
        entry = list(e)
    )
    empty <- res[integer(0), ]
    expect_s4_class(empty, "GwasFineMappingResult")
    expect_equal(nrow(empty), 0L)
    # Should hit the `if (nrow(x) == 0L) return(x)` early-return.
    out <- adjustPips(empty, character(0))
    expect_identical(out, empty)
})

# ===========================================================================
# getTopLoci(type = "GRanges") + aggregate-view branches on FineMappingResultBase
# ===========================================================================
.alc_makeFmr2 <- function() {
    QtlFineMappingResult(
        study = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t2"),
        method = c("susie", "susie"),
        entry = list(.alc_makeFmEntry(), .alc_makeFmEntry()),
        ldSketch = .alc_makeHandle()
    )
}

test_that("getTopLoci(type='GRanges'): a pinned single entry returns a GRanges", {
    fmr <- .alc_makeFmr2()
    gr <- getTopLoci(
        fmr,
        type = "GRanges",
        study = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie"
    )
    expect_s4_class(gr, "GRanges")
})

test_that("getTopLoci(type='GRanges'): aggregating >1 entry is rejected", {
    expect_error(
        getTopLoci(.alc_makeFmr2(), type = "GRanges"),
        "requires type = 'data.frame'"
    )
})

test_that("FineMappingResultBase aggregate view: empty per-entry views -> 0-row frame", {
    # .alc_makeFmEntry has no credible sets, so getCs aggregates to nothing.
    expect_equal(nrow(getCs(.alc_makeFmr2())), 0L)
})

test_that("FineMappingResultBase aggregate view: a no-match selector re-raises the selection error", {
    expect_error(getCs(.alc_makeFmr2(), study = "nope"), "entries")
})

test_that("FineMappingResultBase aggregate view: an empty collection yields a 0-row frame", {
    empty <- QtlFineMappingResult(
        study = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        entry = list(),
        ldSketch = .alc_makeHandle()
    )
    expect_equal(nrow(getCs(empty)), 0L)
})


# ===========================================================================
# Variant reconciliation: intersectVariants + getRetainedMass
# ===========================================================================

.rc_makeFmr <- function(vids, study = "s1", L = 3L, seed = 3L) {
    set.seed(seed)
    p <- length(vids)
    lbf <- matrix(rnorm(L * p, sd = 3), L, p, dimnames = list(NULL, vids))
    alpha <- lbfToAlpha(lbf)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    e <- fineMappingRow(
        vids,
        list(pip = pip, alpha = alpha, lbf_variable = lbf, V = rep(1, L)),
        data.frame(variant_id = vids, pip = pip, stringsAsFactors = FALSE)
    )
    QtlFineMappingResult(
        study = study,
        context = "c1",
        trait = "g1",
        method = "susie",
        entry = list(e)
    )
}

.rc_vids <- function(idx) sprintf("chr1:%d:A:G", 100L * idx)

test_that("intersectVariants restricts both sides to the shared variants", {
    x <- .rc_makeFmr(.rc_vids(1:30))
    y <- .rc_makeFmr(.rc_vids(20:50), study = "s2")
    out <- intersectVariants(x, y)

    expect_setequal(names(out), c("x", "y"))
    shared <- .rc_vids(20:30)
    expect_setequal(
        getVariantIds(pecotmr:::.collectionEntry(out$x, 1L)),
        shared
    )
    expect_setequal(
        getVariantIds(pecotmr:::.collectionEntry(out$y, 1L)),
        shared
    )
    # Both sides end up scored on the SAME variant set, which is what coloc
    # needs -- a pair scored on two different sets is not comparable.
    expect_equal(
        getVariantIds(pecotmr:::.collectionEntry(out$x, 1L)),
        getVariantIds(pecotmr:::.collectionEntry(out$y, 1L))
    )
})

test_that("intersectVariants oneSided adjusts only x", {
    x <- .rc_makeFmr(.rc_vids(1:30))
    y <- .rc_makeFmr(.rc_vids(20:50), study = "s2")
    out <- intersectVariants(x, y, oneSided = TRUE)

    expect_s4_class(out, "QtlFineMappingResult")
    expect_setequal(
        getVariantIds(pecotmr:::.collectionEntry(out, 1L)),
        .rc_vids(20:30)
    )
    # y is untouched: TWAS / MR / cTWAS reconcile the QTL side to the GWAS
    # variant set and leave the GWAS side alone.
    expect_length(unlist(y, use.names = FALSE), 31L)
})

test_that("intersectVariants errors when the two share no variants", {
    # Returning two empty collections would be indistinguishable from
    # "reconciled fine, nothing colocalizes" -- a different conclusion.
    x <- .rc_makeFmr(.rc_vids(1:10))
    y <- .rc_makeFmr(sprintf("chr9:%d:A:G", 100L * 1:10), study = "s2")
    expect_error(intersectVariants(x, y), "share no variants")
})

test_that("intersectVariants matches variants allele-aware", {
    # A chr-prefix difference must not read as no-overlap.
    x <- .rc_makeFmr(.rc_vids(1:10))
    y <- .rc_makeFmr(sprintf("%d:%d:A:G", 1L, 100L * 5:14), study = "s2")
    out <- intersectVariants(x, y)
    expect_length(unlist(out$x, use.names = FALSE), 6L)
})

test_that("getRetainedMass reports per-effect surviving mass", {
    x <- .rc_makeFmr(.rc_vids(1:40))
    adj <- adjustPips(x, .rc_vids(1:20))
    mass <- getRetainedMass(adj)

    expect_equal(nrow(mass), 3L) # one row per effect
    expect_true(all(mass$retainedMass > 0 & mass$retainedMass <= 1))
    expect_true(any(mass$retainedMass < 0.999)) # something was actually lost
    expect_equal(unique(mass$nVariants), 20L)
    # Identity columns come along so the diagnostic is attributable.
    expect_true(all(c("study", "context", "trait", "method") %in% names(mass)))
})

test_that("getRetainedMass is the PRE-renormalization share", {
    # After renormalization every effect row sums to 1, so the surviving share
    # is only knowable at adjustment time. A value of 1 across the board would
    # mean the diagnostic was reading the post-normalized alpha.
    x <- .rc_makeFmr(.rc_vids(1:40))
    adj <- adjustPips(x, .rc_vids(1:20))
    mass <- getRetainedMass(adj)
    expect_false(all(abs(mass$retainedMass - 1) < 1e-9))

    fit <- getSusieFit(pecotmr:::.collectionEntry(adj, 1L))
    expect_equal(rowSums(fit$alpha), rep(1, 3L), tolerance = 1e-10)
})

test_that("getRetainedMass reports nothing for a fit never reconciled", {
    # Not a vector of 1s: nothing was dropped, so there is no surviving-share
    # to report, and implying a subset happened would be misleading.
    x <- .rc_makeFmr(.rc_vids(1:10))
    expect_equal(nrow(getRetainedMass(x)), 0L)
    expect_equal(nrow(getRetainedMass(x[0L])), 0L)
})

test_that("getRetainedMass: empty result has the populated result's schema", {
    # A caller that selects `study` must not break only when there is nothing
    # to report, so the zero-row shape has to match the non-zero-row shape.
    data(qtlFineMappingExample, envir = environment())
    data(gwasFineMappingExample, envir = environment())
    empty <- getRetainedMass(qtlFineMappingExample)
    expect_equal(nrow(empty), 0L)
    both <- intersectVariants(
        qtlFineMappingExample,
        gwasFineMappingExample
    )
    full <- getRetainedMass(both$x)
    expect_gt(nrow(full), 0L)
    expect_equal(colnames(empty), colnames(full))
    expect_equal(nrow(bind_rows(empty, full)), nrow(full))
})

test_that("getRetainedMass: identity columns follow the concrete class", {
    # A GWAS collection has no context / trait, so neither should its table.
    data(gwasFineMappingExample, envir = environment())
    cols <- colnames(getRetainedMass(gwasFineMappingExample))
    expect_true(all(is_in(c("study", "method"), cols)))
    expect_false(any(is_in(c("context", "trait"), cols)))
})


# ---------------------------------------------------------------------------
# getVariantIds on SumStatsBase.
#
# One method serves GwasSumStats and QtlSumStats: each class's own
# getSumStats() knows its selectors and owns the ambiguity error, so this
# delegates rather than re-deciding when a collection is addressable.
# ---------------------------------------------------------------------------

test_that("getVariantIds returns a single-row collection's ids", {
    data(qtlSumStatsExample, envir = environment())
    data(gwasSumStatsS4Example, envir = environment())
    qtl <- getVariantIds(qtlSumStatsExample)
    gwas <- getVariantIds(gwasSumStatsS4Example)
    expect_type(qtl, "character")
    expect_equal(length(qtl), sum(lengths(qtlSumStatsExample)))
    expect_equal(length(gwas), sum(lengths(gwasSumStatsS4Example)))
})

test_that("getVariantIds renders ids the same way the row classes do", {
    # .grVariantIds is the one renderer, so an id is the same string whichever
    # object produced it -- which is what lets these be intersected at all.
    data(qtlSumStatsExample, envir = environment())
    expect_equal(
        getVariantIds(qtlSumStatsExample),
        getSumstatDf(qtlSumStatsExample, require = "Z")$variant_id
    )
    expect_match(getVariantIds(qtlSumStatsExample)[[1L]], "^chr[^:]+:\\d+:")
})

test_that("getVariantIds selects one row of a multi-row collection", {
    data(qtlSumStatsMulticontextExample, envir = environment())
    mc <- qtlSumStatsMulticontextExample
    expect_gt(nrow(mc), 1L)
    ids <- getVariantIds(
        mc,
        study = mc$study[[1L]],
        context = "blood",
        trait = mc$trait[[1L]]
    )
    expect_equal(length(ids), lengths(mc)[[which(mc$context == "blood")]])
})

test_that("getVariantIds defers the ambiguity error to the class", {
    # Not re-implemented here: the message names the selectors that class
    # actually takes, and duplicating it would let the two drift.
    data(qtlSumStatsMulticontextExample, envir = environment())
    expect_error(
        getVariantIds(qtlSumStatsMulticontextExample),
        "Pass `study`, `context`, and `trait`"
    )
})

test_that("getVariantIds on an empty collection defers too", {
    data(qtlSumStatsExample, envir = environment())
    expect_error(getVariantIds(qtlSumStatsExample[0]), "no rows")
})

# ===========================================================================
# adjustPips: entries that share no variant with `keepVariants`
#
# The zero-row short-circuit is covered above and the all-overlap path by
# test_ColocResult; the two partial/disjoint branches were not.
# ===========================================================================

# Two entries built by .rc_makeFmr (which carries a real alpha matrix, as
# renormalization requires): one on chr1, one on chr9 so it can be made to
# miss `keepVariants` entirely.
.ac_twoEntry <- function() {
    a <- .rc_makeFmr(.rc_vids(1:10), study = "s1")
    b <- .rc_makeFmr(sprintf("chr9:%d:A:G", 100L * 1:10), study = "s2")
    QtlFineMappingResult(
        study = c("s1", "s2"),
        context = c("c1", "c1"),
        trait = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(.collectionEntries(a)[[1L]], .collectionEntries(b)[[1L]])
    )
}

test_that("adjustPips refuses a keepVariants set disjoint from every entry", {
    # Nothing to renormalize anywhere: silently returning the input would
    # hand back unadjusted PIPs that look adjusted.
    expect_error(
        adjustPips(.ac_twoEntry(), "chr22:999999:A:G"),
        "the two variant sets are disjoint"
    )
})

test_that("adjustPips drops non-overlapping entries and says so", {
    # Matches subsetRegion's rule: an element that trims to zero variants
    # goes away rather than lingering unadjusted.
    expect_message(
        out <- adjustPips(.ac_twoEntry(), .rc_vids(1:5)),
        "dropping 1 of 2 entries"
    )
    expect_equal(nrow(out), 1L)
    expect_equal(as.character(mcols(out)$study), "s1")
})
