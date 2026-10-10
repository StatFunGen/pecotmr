context("AllClasses (virtual base classes)")

# Most slots / accessors on the concrete subclasses are exercised in their
# own test files; these tests target the *base-class* behaviors that the
# concrete subclasses inherit without overriding (studyName on SumStatsBase,
# qcDiagnostics body branches, and the zero-row adjustPips short-circuit
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
        studyName = "g1",
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
# SumStatsBase: studyName (inherited by QtlSumStats / GwasSumStats)
# ===========================================================================

test_that("SumStatsBase: studyName on a GwasSumStats returns unique study names", {
    ss <- .alc_makeGwasSumStats()
    expect_equal(studyName(ss), "g1")
})

# ===========================================================================
# SumStatsBase: qcDiagnostics — every branch
# ===========================================================================

test_that("SumStatsBase: qcDiagnostics returns NULL on empty qcInfo", {
    ss <- .alc_makeGwasSumStats() # qcInfo = list() by default
    expect_null(qcDiagnostics(ss))
})

test_that("SumStatsBase: qcDiagnostics returns NULL when entryAudit slot is absent", {
    # qcInfo has steps but no entryAudit -> nothing to return.
    ss <- .alc_makeGwasSumStats(qcInfo = list(step1 = "ok"))
    expect_null(qcDiagnostics(ss))
})

test_that("SumStatsBase: qcDiagnostics returns the per-entry diagnostics by index", {
    diag1 <- data.frame(SNP = "rs1", outlier = FALSE, stringsAsFactors = FALSE)
    diag2 <- data.frame(SNP = "rs2", outlier = TRUE, stringsAsFactors = FALSE)
    qc <- list(
        entryAudit = list(
            list(ldMismatchDiagnostics = diag1),
            list(ldMismatchDiagnostics = diag2)
        )
    )
    ss <- .alc_makeGwasSumStats(qcInfo = qc)
    expect_identical(qcDiagnostics(ss, entry = 1L), diag1)
    expect_identical(qcDiagnostics(ss, entry = 2L), diag2)
})

test_that("SumStatsBase: qcDiagnostics(entry = NULL) returns the populated entries only", {
    diag1 <- data.frame(SNP = "rs1", outlier = FALSE, stringsAsFactors = FALSE)
    # Entry 2's audit has no ldMismatchDiagnostics field; should be filtered.
    qc <- list(
        entryAudit = list(
            list(ldMismatchDiagnostics = diag1),
            list(other = "no diagnostics here")
        )
    )
    ss <- .alc_makeGwasSumStats(qcInfo = qc)
    out <- qcDiagnostics(ss, entry = NULL)
    expect_type(out, "list")
    expect_equal(length(out), 1L)
    expect_named(out, "1")
    expect_identical(out[["1"]], diag1)
})

test_that("SumStatsBase: qcDiagnostics(entry = NULL) returns NULL when no entry has diagnostics", {
    qc <- list(entryAudit = list(list(other = 1), list(other = 2)))
    ss <- .alc_makeGwasSumStats(qcInfo = qc)
    expect_null(qcDiagnostics(ss, entry = NULL))
})

test_that("SumStatsBase: qcDiagnostics errors on out-of-range entry", {
    qc <- list(
        entryAudit = list(list(ldMismatchDiagnostics = data.frame(z = 1)))
    )
    ss <- .alc_makeGwasSumStats(qcInfo = qc)
    expect_error(qcDiagnostics(ss, entry = 0L), "must be a single integer")
    expect_error(qcDiagnostics(ss, entry = 99L), "must be a single integer")
    expect_error(
        qcDiagnostics(ss, entry = c(1L, 2L)),
        "must be a single integer"
    )
    expect_error(
        qcDiagnostics(ss, entry = "first"),
        "must be a single integer"
    )
})

# ===========================================================================
# FineMappingResultBase: adjustPips zero-row short-circuit
# ===========================================================================

test_that("FineMappingResultBase: adjustPips on a zero-row collection returns the input unchanged", {
    e <- .alc_makeFmEntry(3)
    res <- GwasFineMappingResult(
        studyName = "g1",
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
# topLoci() ranges + aggregate-view branches on FineMappingResultBase
# ===========================================================================
.alc_makeFmr2 <- function() {
    QtlFineMappingResult(
        studyName = c("Q1", "Q1"),
        context = c("c1", "c1"),
        trait = c("t1", "t2"),
        method = c("susie", "susie"),
        entry = list(.alc_makeFmEntry(), .alc_makeFmEntry()),
        ldSketch = .alc_makeHandle()
    )
}

test_that("topLoci(): a pinned single entry returns a GRanges", {
    fmr <- .alc_makeFmr2()
    gr <- topLoci(
        fmr,
        studyName = "Q1",
        context = "c1",
        trait = "t1",
        method = "susie"
    )
    expect_s4_class(gr, "GRanges")
})

test_that("topLoci(): aggregating >1 entry also returns a GRanges", {
    # Converting to ranges after stacking is what makes this work: the
    # row-identity columns become mcols. The old type = "GRanges" branch had
    # to refuse anything but a single pinned entry.
    gr <- topLoci(.alc_makeFmr2())
    expect_s4_class(gr, "GRanges")
    expect_true(all(
        c("study", "context", "trait", "method") %in%
            names(S4Vectors::mcols(gr))
    ))
})

test_that("topLoci(raw = TRUE) stays a table", {
    # The stored canonical table has no posterior projection to put on a
    # range, so raw is the one form that is not a GRanges.
    raw <- topLoci(.alc_makeFmr2(), raw = TRUE)
    expect_true(is.data.frame(raw))
})

test_that("FineMappingResultBase aggregate view: empty per-entry views -> 0-row frame", {
    # .alc_makeFmEntry has no credible sets, so credibleSets aggregates to nothing.
    expect_equal(nrow(credibleSets(.alc_makeFmr2())), 0L)
})

test_that("FineMappingResultBase aggregate view: a no-match selector re-raises the selection error", {
    expect_error(credibleSets(.alc_makeFmr2(), studyName = "nope"), "entries")
})

test_that("FineMappingResultBase aggregate view: an empty collection yields a 0-row frame", {
    empty <- QtlFineMappingResult(
        studyName = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        entry = list(),
        ldSketch = .alc_makeHandle()
    )
    expect_equal(nrow(credibleSets(empty)), 0L)
})


# ===========================================================================
# Variant reconciliation: intersectVariants + retainedMass
# ===========================================================================

.rc_makeFmr <- function(vids, studyName = "s1", maxNumSingleEffects = 3L, seed = 3L) {
    set.seed(seed)
    p <- length(vids)
    lbf <- matrix(rnorm(maxNumSingleEffects * p, sd = 3), maxNumSingleEffects, p, dimnames = list(NULL, vids))
    alpha <- lbfToAlpha(lbf)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    e <- fineMappingRow(
        vids,
        list(pip = pip, alpha = alpha, lbf_variable = lbf, V = rep(1, maxNumSingleEffects)),
        data.frame(variant_id = vids, pip = pip, stringsAsFactors = FALSE)
    )
    QtlFineMappingResult(
        studyName = studyName,
        context = "c1",
        trait = "g1",
        method = "susie",
        entry = list(e)
    )
}

.rc_vids <- function(idx) sprintf("chr1:%d:A:G", 100L * idx)

test_that("intersectVariants restricts both sides to the shared variants", {
    x <- .rc_makeFmr(.rc_vids(1:30))
    y <- .rc_makeFmr(.rc_vids(20:50), studyName = "s2")
    out <- intersectVariants(x, y)

    expect_setequal(names(out), c("x", "y"))
    shared <- .rc_vids(20:30)
    expect_setequal(
        variantIds(pecotmr:::.collectionEntry(out$x, 1L)),
        shared
    )
    expect_setequal(
        variantIds(pecotmr:::.collectionEntry(out$y, 1L)),
        shared
    )
    # Both sides end up scored on the SAME variant set, which is what coloc
    # needs -- a pair scored on two different sets is not comparable.
    expect_equal(
        variantIds(pecotmr:::.collectionEntry(out$x, 1L)),
        variantIds(pecotmr:::.collectionEntry(out$y, 1L))
    )
})

test_that("intersectVariants oneSided adjusts only x", {
    x <- .rc_makeFmr(.rc_vids(1:30))
    y <- .rc_makeFmr(.rc_vids(20:50), studyName = "s2")
    out <- intersectVariants(x, y, oneSided = TRUE)

    expect_s4_class(out, "QtlFineMappingResult")
    expect_setequal(
        variantIds(pecotmr:::.collectionEntry(out, 1L)),
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
    y <- .rc_makeFmr(sprintf("chr9:%d:A:G", 100L * 1:10), studyName = "s2")
    expect_error(intersectVariants(x, y), "share no variants")
})

test_that("intersectVariants matches variants allele-aware", {
    # A chr-prefix difference must not read as no-overlap.
    x <- .rc_makeFmr(.rc_vids(1:10))
    y <- .rc_makeFmr(sprintf("%d:%d:A:G", 1L, 100L * 5:14), studyName = "s2")
    out <- intersectVariants(x, y)
    expect_length(unlist(out$x, use.names = FALSE), 6L)
})

test_that("retainedMass reports per-effect surviving mass", {
    x <- .rc_makeFmr(.rc_vids(1:40))
    adj <- adjustPips(x, .rc_vids(1:20))
    mass <- retainedMass(adj)

    expect_equal(nrow(mass), 3L) # one row per effect
    expect_true(all(mass$retainedMass > 0 & mass$retainedMass <= 1))
    expect_true(any(mass$retainedMass < 0.999)) # something was actually lost
    expect_equal(unique(mass$nVariants), 20L)
    # Identity columns come along so the diagnostic is attributable.
    expect_true(all(c("study", "context", "trait", "method") %in% names(mass)))
})

test_that("retainedMass is the PRE-renormalization share", {
    # After renormalization every effect row sums to 1, so the surviving share
    # is only knowable at adjustment time. A value of 1 across the board would
    # mean the diagnostic was reading the post-normalized alpha.
    x <- .rc_makeFmr(.rc_vids(1:40))
    adj <- adjustPips(x, .rc_vids(1:20))
    mass <- retainedMass(adj)
    expect_false(all(abs(mass$retainedMass - 1) < 1e-9))

    fit <- susieFit(pecotmr:::.collectionEntry(adj, 1L))
    expect_equal(rowSums(fit$alpha), rep(1, 3L), tolerance = 1e-10)
})

test_that("retainedMass reports nothing for a fit never reconciled", {
    # Not a vector of 1s: nothing was dropped, so there is no surviving-share
    # to report, and implying a subset happened would be misleading.
    x <- .rc_makeFmr(.rc_vids(1:10))
    expect_equal(nrow(retainedMass(x)), 0L)
    expect_equal(nrow(retainedMass(x[0L])), 0L)
})

test_that("retainedMass: empty result has the populated result's schema", {
    # A caller that selects `study` must not break only when there is nothing
    # to report, so the zero-row shape has to match the non-zero-row shape.
    data(qtlFineMappingExample, envir = environment())
    data(gwasFineMappingExample, envir = environment())
    empty <- retainedMass(qtlFineMappingExample)
    expect_equal(nrow(empty), 0L)
    both <- intersectVariants(
        qtlFineMappingExample,
        gwasFineMappingExample
    )
    full <- retainedMass(both$x)
    expect_gt(nrow(full), 0L)
    expect_equal(colnames(empty), colnames(full))
    expect_equal(nrow(bind_rows(empty, full)), nrow(full))
})

test_that("retainedMass: identity columns follow the concrete class", {
    # A GWAS collection has no context / trait, so neither should its table.
    data(gwasFineMappingExample, envir = environment())
    cols <- colnames(retainedMass(gwasFineMappingExample))
    expect_true(all(is_in(c("study", "method"), cols)))
    expect_false(any(is_in(c("context", "trait"), cols)))
})


# ---------------------------------------------------------------------------
# variantIds on SumStatsBase.
#
# One method serves GwasSumStats and QtlSumStats: each class's own
# sumStats() knows its selectors and owns the ambiguity error, so this
# delegates rather than re-deciding when a collection is addressable.
# ---------------------------------------------------------------------------

test_that("variantIds returns a single-row collection's ids", {
    data(qtlSumStatsExample, envir = environment())
    data(gwasSumStatsS4Example, envir = environment())
    qtl <- variantIds(qtlSumStatsExample)
    gwas <- variantIds(gwasSumStatsS4Example)
    expect_type(qtl, "character")
    expect_equal(length(qtl), sum(lengths(qtlSumStatsExample)))
    expect_equal(length(gwas), sum(lengths(gwasSumStatsS4Example)))
})

test_that("variantIds renders ids the same way the row classes do", {
    # .grVariantIds is the one renderer, so an id is the same string whichever
    # object produced it -- which is what lets these be intersected at all.
    data(qtlSumStatsExample, envir = environment())
    expect_equal(
        variantIds(qtlSumStatsExample),
        as.data.frame(qtlSumStatsExample, require = "Z")$variant_id
    )
    expect_match(variantIds(qtlSumStatsExample)[[1L]], "^chr[^:]+:\\d+:")
})

test_that("variantIds selects one row of a multi-row collection", {
    data(qtlSumStatsMulticontextExample, envir = environment())
    mc <- qtlSumStatsMulticontextExample
    expect_gt(nrow(mc), 1L)
    ids <- variantIds(
        mc,
        studyName = mc$study[[1L]],
        context = "blood",
        trait = mc$trait[[1L]]
    )
    expect_equal(length(ids), lengths(mc)[[which(mc$context == "blood")]])
})

test_that("variantIds defers the ambiguity error to the class", {
    # Not re-implemented here: the message names the selectors that class
    # actually takes, and duplicating it would let the two drift.
    data(qtlSumStatsMulticontextExample, envir = environment())
    expect_error(
        variantIds(qtlSumStatsMulticontextExample),
        "Pass `study`, `context`, and `trait`"
    )
})

test_that("variantIds on an empty collection defers too", {
    data(qtlSumStatsExample, envir = environment())
    expect_error(variantIds(qtlSumStatsExample[0]), "no rows")
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
    a <- .rc_makeFmr(.rc_vids(1:10), studyName = "s1")
    b <- .rc_makeFmr(sprintf("chr9:%d:A:G", 100L * 1:10), studyName = "s2")
    QtlFineMappingResult(
        studyName = c("s1", "s2"),
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


# ===========================================================================
# .ssStitchElements()'s `ranges` filter
#
# sumStats(x, ranges = ) is a documented argument that nothing exercised:
# every call left it NULL and returned through the early exit, so the
# chromosome test and the overlap intersection were both cold.
# ===========================================================================

test_that("sumStats(ranges=) narrows the variants to the window", {
    data(gwasSumStatsS4Example)
    full <- sumStats(gwasSumStatsS4Example)
    span <- range(GenomicRanges::start(full))
    win <- GenomicRanges::GRanges(
        as.character(GenomicRanges::seqnames(full))[[1L]],
        IRanges::IRanges(
            span[[1L]],
            span[[1L]] + (span[[2L]] - span[[1L]]) %/% 4L
        )
    )
    sub <- sumStats(gwasSumStatsS4Example, ranges = win)
    expect_s4_class(sub, "GRanges")
    expect_lt(length(sub), length(full))
    expect_gt(length(sub), 0L)
    expect_true(all(IRanges::overlapsAny(sub, win)))
})

test_that("sumStats(ranges=) on another chromosome returns nothing", {
    # The chromosome test short-circuits before overlapsAny(), which would
    # otherwise warn about disjoint seqlevels.
    data(gwasSumStatsS4Example)
    off <- GenomicRanges::GRanges("chrZZ", IRanges::IRanges(1L, 100L))
    sub <- suppressWarnings(sumStats(gwasSumStatsS4Example, ranges = off))
    expect_s4_class(sub, "GRanges")
    expect_equal(length(sub), 0L)
})


test_that("the genome check names a missing build and a mixed one", {
    f <- pecotmr:::.sumStatsCheckGenome
    gr <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(c(100L, 200L), width = 1L)
    )
    # No build recorded anywhere.
    expect_match(f(gr), "no genome build in seqinfo")
    # Two different builds across seqlevels.
    GenomeInfoDb::genome(gr) <- c("hg19", "hg38")
    expect_match(f(gr), "names more than one genome build")
    # A single build is accepted.
    GenomeInfoDb::genome(gr) <- "hg38"
    expect_null(f(gr))
})


test_that("an empty collection has no variants", {
    setClass(
        "RcEmptyKid",
        contains = "RangedTupleList",
        representation(genome = "character", ldSketch = "ANY", qcInfo = "list")
    )
    empty <- new(
        "RcEmptyKid",
        GenomicRanges::GRangesList(),
        genome = "hg38",
        ldSketch = NULL,
        qcInfo = list()
    )
    expect_equal(nrow(empty), 0L)
    expect_equal(pecotmr:::.rcAllVariants(empty), character(0))
})

test_that("retained mass is NULL for a fit that never went through reconciliation", {
    # NULL rather than a vector of 1s: nothing was dropped, so there is no
    # mass to report.
    expect_null(pecotmr:::.rcFitRetainedMass(NULL))
    expect_null(pecotmr:::.rcFitRetainedMass("not a fit"))
    expect_equal(pecotmr:::.rcFitRetainedMass(list(retained_mass = 0.8)), 0.8)
})

# ===========================================================================
# Variant-level replacement: af<- / maf<- / nSamples<- dispatch on `value`
# ===========================================================================

# A replacement GRanges over the same variants as `x`, with the alleles in
# the given orientation and one payload column.
.alc_replacementGr <- function(x, swap = FALSE, k = NULL, payload = 0.3) {
    ids <- variantIds(x)
    p <- parseVariantId(ids)
    idx <- if (is.null(k)) seq_along(ids) else seq_len(k)
    gr <- GenomicRanges::GRanges(
        p$chrom[idx],
        IRanges::IRanges(p$pos[idx], width = 1L)
    )
    a1 <- if (swap) p$A2[idx] else p$A1[idx]
    a2 <- if (swap) p$A1[idx] else p$A2[idx]
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        A1 = a1,
        A2 = a2,
        v = rep(payload, length(idx))
    )
    gr
}

test_that("af<- numeric is positional and checks length hard", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    n <- length(variantIds(x))
    # The motivating case: the source carried no effect-allele frequency.
    expect_null(af(x))
    af(x) <- seq(0.1, 0.4, length.out = n)
    expect_equal(length(af(x)), n)
    expect_equal(af(x)[[1L]], 0.1)
    # A recycled or truncated frequency vector is silently wrong, so length
    # is refused rather than recycled.
    expect_error(af(x) <- c(0.1, 0.2), "one value per variant")
    expect_error(af(x) <- 0.3, "one value per variant")
})

test_that("af<- GRanges matches on (chrom, pos, A1, A2) regardless of order", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    gr <- .alc_replacementGr(x, payload = 0.25)
    af(x) <- gr[rev(seq_along(gr))]
    expect_equal(unique(af(x)), 0.25)
})

test_that("a flip complements af but leaves maf and nSamples alone", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    swapped <- .alc_replacementGr(x, swap = TRUE, payload = 0.3)
    # An allele swap sends af -> 1 - af. maf is min(af, 1 - af) and nSamples
    # is a count, so both are swap-invariant -- complementing them would
    # corrupt them.
    y <- x
    af(y) <- swapped
    expect_equal(unique(af(y)), 0.7)
    z <- x
    maf(z) <- swapped
    expect_equal(unique(maf(z)), 0.3)
    w <- x
    nSamples(w) <- swapped
    expect_equal(unique(nSamples(w)), 0.3)
})

test_that("af<- GRanges requires A1/A2 and exactly one payload column", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    n <- length(variantIds(x))
    p <- parseVariantId(variantIds(x))
    bare <- GenomicRanges::GRanges(
        p$chrom,
        IRanges::IRanges(p$pos, width = 1L)
    )
    S4Vectors::mcols(bare) <- S4Vectors::DataFrame(v = rep(0.3, n))
    # Without A1/A2 a flip cannot be detected, so the value cannot be
    # oriented and the write is refused rather than guessed at.
    expect_error(af(x) <- bare, "must carry A1 and A2")
    two <- .alc_replacementGr(x)
    S4Vectors::mcols(two)$extra <- 1
    expect_error(af(x) <- two, "exactly one value column")
})

test_that("af<- sets NA for entry variants the replacement omits", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    n <- length(variantIds(x))
    partial <- .alc_replacementGr(x, k = 5L, payload = 0.42)
    expect_warning(af(x) <- partial, "were set NA")
    a <- af(x)
    # Unknown, not unchanged.
    expect_equal(sum(!is.na(a)), 5L)
    expect_equal(sum(is.na(a)), n - 5L)
    expect_equal(unique(a[!is.na(a)]), 0.42)
})

test_that("af<- refuses a replacement that matches nothing", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    far <- GenomicRanges::GRanges(
        "chr9",
        IRanges::IRanges(1:3, width = 1L)
    )
    S4Vectors::mcols(far) <- S4Vectors::DataFrame(
        A1 = rep("A", 3),
        A2 = rep("G", 3),
        v = c(0.1, 0.2, 0.3)
    )
    expect_error(af(x) <- far, "no variant in the replacement matches")
})
