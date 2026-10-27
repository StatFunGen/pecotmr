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
        nSamples = 10L,
        sampleIds = paste0("s", seq_len(10L)),
        pgenPtr = NULL
    )
}

.alc_makeGr <- function(n = 3, chrom = "chr1") {
    gr <- GenomicRanges::GRanges(
        chrom,
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

.alc_makeGwasSumStats <- function(qcOptions = list()) {
    GwasSumStats(
        studyName = "g1",
        entry = list(.alc_makeGr()),
        genome = "hg19",
        ldReference = .alc_makeHandle(),
        qcOptions = qcOptions
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
    FineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl
    )
}

# ===========================================================================
# SumStatsBase: studyName (inherited by QtlSumStats / GwasSumStats)
# ===========================================================================

test_that("SumStatsBase: studyName on a GwasSumStats names each study once", {
    # g1's variants span two chromosomes, so it is split into two rows and
    # the study column repeats it; studyName() names each study once.
    ss <- GwasSumStats(
        studyName = c("g1", "g2"),
        entry = list(
            .alc_makeGr(chrom = c("chr1", "chr1", "chr2")),
            .alc_makeGr()
        ),
        genome = "hg19",
        ldReference = .alc_makeHandle()
    )
    expect_identical(ss$study, c("g1", "g1", "g2"))
    expect_identical(studyName(ss), c("g1", "g2"))
})

# ===========================================================================
# SumStatsBase: qcDiagnostics — every branch
# ===========================================================================

test_that("SumStatsBase: qcDiagnostics returns NULL without per-row audits", {
    # QC options recorded but no qcAudit column -> nothing to return.
    ss <- .alc_makeGwasSumStats(qcOptions = list(step1 = "ok"))
    expect_null(qcDiagnostics(ss))
})

# A two-study GwasSumStats carrying one QC audit per row.
.alc_twoStudyQc <- function(audits) {
    x <- GwasSumStats(
        studyName = c("g1", "g2"),
        entry = list(.alc_makeGr(), .alc_makeGr()),
        genome = "hg19",
        ldReference = .alc_makeHandle(),
        qcOptions = list(step1 = "ok")
    )
    x$qcAudit <- S4Vectors::SimpleList(audits)
    x
}

.alc_diagQc <- function(...) {
    list(...)
}

test_that("SumStatsBase: qcDiagnostics returns the selected study's table", {
    diag1 <- data.frame(SNP = "rs1", outlier = FALSE)
    diag2 <- data.frame(SNP = "rs2", outlier = TRUE)
    ss <- .alc_twoStudyQc(.alc_diagQc(
        list(ldMismatchDiagnostics = diag1),
        list(ldMismatchDiagnostics = diag2)
    ))
    out <- qcDiagnostics(ss, studyName = "g2")
    expect_s4_class(out, "DataFrame")
    expect_equal(as.data.frame(out), diag2)
    expect_equal(as.data.frame(qcDiagnostics(ss, studyName = "g1")), diag1)
    # Two studies and no selector: which one is meant cannot be known.
    expect_error(qcDiagnostics(ss), "2 studies")
})

test_that("SumStatsBase: qcDiagnostics is NULL for a selection without one", {
    ss <- .alc_twoStudyQc(.alc_diagQc(
        list(ldMismatchDiagnostics = data.frame(SNP = "rs1")),
        list(other = "no diagnostics here")
    ))
    expect_null(qcDiagnostics(ss, studyName = "g2"))
})

test_that("SumStatsBase: qcDiagnostics stacks a split study's rows", {
    data(gwasSumStatsS4Example)
    g <- gwasSumStatsS4Example
    split <- GwasSumStats(
        studyName = "trait1",
        entry = list(sumStats(g)),
        genome = "hg19",
        ldReference = ldReference(g)[[1L]],
        qcOptions = list(step1 = "ok"),
        region = GenomicRanges::GRanges(
            "chr22",
            IRanges::IRanges(c(1L, 15680321L), c(15680320L, 51304566L))
        )
    )
    split$qcAudit <- S4Vectors::SimpleList(.alc_diagQc(
        list(ldMismatchDiagnostics = data.frame(SNP = "rs1")),
        list(ldMismatchDiagnostics = data.frame(SNP = c("rs2", "rs3")))
    ))
    expect_equal(nrow(split), 2L)
    expect_equal(qcDiagnostics(split)$SNP, c("rs1", "rs2", "rs3"))
})

test_that("subsetting keeps each element's QC audit with it", {
    # Each row's audit is a column of mcols(), so it subsets with the row: a
    # subset never reports another element's diagnostics.
    diag2 <- data.frame(SNP = "rs2", outlier = TRUE)
    ss <- .alc_twoStudyQc(.alc_diagQc(
        list(ldMismatchDiagnostics = data.frame(SNP = "rs1")),
        list(ldMismatchDiagnostics = diag2)
    ))
    expect_equal(as.data.frame(qcDiagnostics(ss[2])), diag2)
    expect_length(ss[2]$qcAudit, 1L)
    expect_equal(as.data.frame(qcDiagnostics(ss[c(FALSE, TRUE)])), diag2)
    expect_equal(as.data.frame(qcDiagnostics(ss[-1])), diag2)
    # A collection without an audit subsets as before.
    expect_null(.alc_makeGwasSumStats()[1]$qcAudit)
})

# ===========================================================================
# SumStatsBase: qcInfo / qcOptions, and the saved-object migration
# ===========================================================================

test_that("qcInfo() is one row per tuple: its identity, then each count", {
    ss <- .alc_twoStudyQc(.alc_diagQc(
        list(
            variantsIn = 10L,
            contentFilters = list(nDropped = 2L),
            ldMismatchDiagnostics = data.frame(SNP = "rs1")
        ),
        list(variantsIn = 8L)
    ))
    info <- qcInfo(ss)
    expect_s4_class(info, "DataFrame")
    expect_equal(info$study, c("g1", "g2"))
    expect_equal(info$variantsIn, c(10L, 8L))
    # A nested count is named by its path; a row the step skipped is NA.
    expect_equal(info$contentFilters.nDropped, c(2L, NA))
    # A per-variant table is qcDiagnostics()'s, never a column here.
    expect_false(is_in("ldMismatchDiagnostics", colnames(info)))
})

test_that("qcInfo() before QC names the rows and records nothing", {
    info <- qcInfo(.alc_makeGwasSumStats())
    expect_equal(nrow(info), 1L)
    expect_equal(colnames(info), c("study", "blockId"))
})

test_that("an audit that is not a list contributes no counts", {
    expect_equal(pecotmr:::.qcAuditScalars(NULL), list())
    expect_equal(pecotmr:::.qcAuditScalars(list()), list())
    expect_equal(pecotmr:::.qcAuditScalars(data.frame(a = 1)), list())
    # Only scalars are counts: a vector is left out.
    expect_equal(
        pecotmr:::.qcAuditScalars(list(a = 1L, b = 1:3)),
        list(a = 1L)
    )
})

test_that("qcOptions() returns the settings, separate from the audits", {
    ss <- .alc_makeGwasSumStats(qcOptions = list(step1 = "ok"))
    expect_equal(qcOptions(ss), list(step1 = "ok"))
    expect_equal(qcOptions(.alc_makeGwasSumStats()), list())
})

test_that(".rowQcAudit reads one row's audit, NULL where there is none", {
    ss <- .alc_twoStudyQc(.alc_diagQc(list(a = 1), list(a = 2)))
    expect_equal(pecotmr:::.rowQcAudit(ss, 2L), list(a = 2))
    expect_null(pecotmr:::.rowQcAudit(ss, 3L))
    expect_null(pecotmr:::.rowQcAudit(.alc_makeGwasSumStats(), 1L))
})

# `x` as it was saved before the per-row audits moved into mcols(): the
# whole QC record in a `qcInfo` slot, and no `qcOptions` slot.
.alc_savedWithQcInfo <- function(x, qcInfo) {
    attr(x, "qcOptions") <- NULL
    attr(x, "qcInfo") <- qcInfo
    x
}

test_that("updateObject() splits a saved QC record into options and audits", {
    ss <- .alc_twoStudyQc(.alc_diagQc(list(), list()))
    ss$qcAudit <- NULL
    saved <- .alc_savedWithQcInfo(
        ss,
        list(
            step1 = "ok",
            entryAudit = list(list(variantsIn = 3L), list(variantsIn = 4L))
        )
    )
    up <- updateObject(saved)
    expect_true(validObject(up))
    expect_null(attr(up, "qcInfo", exact = TRUE))
    expect_equal(qcOptions(up), list(step1 = "ok"))
    expect_equal(qcInfo(up)$variantsIn, c(3L, 4L))
})

test_that("updateObject() keeps the options of a record with no audits", {
    ss <- .alc_makeGwasSumStats()
    up <- updateObject(.alc_savedWithQcInfo(ss, list(step1 = "ok")))
    expect_equal(qcOptions(up), list(step1 = "ok"))
    expect_null(up$qcAudit)
})

test_that("updateObject() drops audits it cannot match to rows, saying so", {
    ss <- .alc_twoStudyQc(.alc_diagQc(list(), list()))
    ss$qcAudit <- NULL
    saved <- .alc_savedWithQcInfo(
        ss,
        list(
            step1 = "ok",
            entryAudit = list(list(variantsIn = 3L))
        )
    )
    expect_warning(
        up <- updateObject(saved),
        "1 row audit\\(s\\) for 2 row\\(s\\)"
    )
    expect_equal(qcOptions(up), list(step1 = "ok"))
    expect_null(up$qcAudit)
})

test_that("updateObject() leaves a current summary-statistics object as is", {
    ss <- .alc_makeGwasSumStats(qcOptions = list(step1 = "ok"))
    up <- updateObject(ss, verbose = FALSE)
    expect_equal(qcOptions(up), list(step1 = "ok"))
    expect_equal(as.list(up), as.list(ss))
    expect_message(updateObject(ss, verbose = TRUE), "updateObject")
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
        ldReference = .alc_makeHandle()
    )
}

test_that("topLoci() has no type switch: the stored table stays internal", {
    # `raw = TRUE` once returned the stored table, a data.frame, from the same
    # accessor that otherwise returns ranges.
    expect_error(topLoci(.alc_makeFmr2(), raw = TRUE), "unused argument")
})

test_that("QtlFineMappingResult: topLoci aggregates entries when selectors do not pin one", {
    # Two entries, no selectors: the old behaviour errored ("2 entries; pass
    # study/context/trait/method"). Now it stacks the per-variant tables,
    # prefixed with the row identity so variants stay attributable.
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    agg <- topLoci(res, signalCutoff = 0)
    expect_s4_class(agg, "GRanges")
    expect_equal(length(agg), 5L)
    expect_true(all(
        c("study", "context", "trait", "blockId", "method") %in%
            names(S4Vectors::mcols(agg))
    ))
    expect_equal(agg$context, c("c1", "c1", "c1", "c2", "c2"))
    expect_true("variant_id" %in% names(S4Vectors::mcols(agg)))
    # QTL results key on context/trait, so blockId is NA-filled.
    expect_true(all(is.na(agg$blockId)))
    # the added per-variant `method` must not duplicate the identity column
    expect_equal(sum(names(S4Vectors::mcols(agg)) == "method"), 1L)
})

test_that("QtlFineMappingResult: topLoci with a full tuple keeps the bare (id-free) table", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    tl <- topLoci(
        res,
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        signalCutoff = 0
    )
    expect_s4_class(tl, "GRanges")
    expect_equal(length(tl), 3L)
    # single-entry fast path, unchanged
    expect_false(is_in("context", names(S4Vectors::mcols(tl))))
})

test_that("QtlFineMappingResult: topLoci aggregates only the matching subset", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    sub <- topLoci(res, context = "c2", signalCutoff = 0)
    expect_equal(length(sub), 2L)
    expect_equal(unique(sub$context), "c2")
})

test_that("QtlFineMappingResult: topLoci on a single-row collection still carries identity columns", {
    # A no-selector call takes the aggregate path even for one row, so the PIP
    # table has a uniform shape (identity columns present) whatever the row count
    # -- callers don't special-case single-entry results.
    row <- .sc_makeFineMappingRow(4)
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(row)
    )
    tl <- topLoci(res, signalCutoff = 0)
    expect_equal(length(tl), 4L)
    # The row's own (projected) table, variant for variant.
    expect_equal(
        tl$variant_id,
        .fmrRowTopLoci(row, signalCutoff = 0)$variant_id
    )
    expect_true(all(
        c("study", "context", "trait", "method", "variant_id") %in%
            names(S4Vectors::mcols(tl))
    ))
    expect_equal(unique(tl$context), "c1")
})

test_that("aggregate view: empty per-entry views give empty ranges", {
    # .alc_makeFmEntry has no credible sets, so credibleSets aggregates to nothing.
    expect_length(credibleSets(.alc_makeFmr2()), 0L)
})

test_that("FineMappingResultBase aggregate view: a no-match selector re-raises the selection error", {
    expect_error(credibleSets(.alc_makeFmr2(), studyName = "nope"), "entries")
})

test_that("aggregate view: an empty collection gives empty ranges", {
    empty <- QtlFineMappingResult(
        studyName = character(0),
        context = character(0),
        trait = character(0),
        method = character(0),
        entry = list(),
        ldReference = .alc_makeHandle()
    )
    expect_length(credibleSets(empty), 0L)
})

# ===========================================================================
# Variant reconciliation: intersectVariants + retainedMass
# ===========================================================================

.rc_makeFmr <- function(
    vids,
    studyName = "s1",
    maxNumSingleEffects = 3L,
    seed = 3L
) {
    set.seed(seed)
    p <- length(vids)
    lbf <- matrix(
        rnorm(maxNumSingleEffects * p, sd = 3),
        maxNumSingleEffects,
        p,
        dimnames = list(NULL, vids)
    )
    alpha <- lbfToAlpha(lbf)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    e <- FineMappingRow(
        vids,
        list(
            pip = pip,
            alpha = alpha,
            lbf_variable = lbf,
            V = rep(1, maxNumSingleEffects)
        ),
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

test_that("retainedMass is each effect's pre-renormalization surviving mass", {
    x <- .rc_makeFmr(.rc_vids(1:40))
    adj <- adjustPips(x, .rc_vids(1:20))
    mass <- retainedMass(adj)

    expect_equal(nrow(mass), 3L) # one row per effect
    expect_true(all(mass$retainedMass > 0 & mass$retainedMass <= 1))
    expect_true(any(mass$retainedMass < 0.999)) # something was actually lost
    expect_equal(unique(mass$nVariants), 20L)
    # Identity columns come along so the diagnostic is attributable.
    expect_true(all(c("study", "context", "trait", "method") %in% names(mass)))

    # After renormalization every effect row sums to 1, so the surviving share
    # is only knowable at adjustment time. A value of 1 across the board would
    # mean the diagnostic was reading the post-normalized alpha.
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
    expect_s4_class(full, "DataFrame")
    expect_equal(nrow(rbind(empty, full)), nrow(full))
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
        pecotmr:::.entryToSumstatDf(
            sumStats(qtlSumStatsExample),
            require = "Z"
        )$variant_id
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
        "Pass `studyName`, `context`, and `trait`"
    )
})

test_that("variantIds on an empty collection defers too", {
    data(qtlSumStatsExample, envir = environment())
    expect_error(variantIds(qtlSumStatsExample[0]), "no rows")
})

# ===========================================================================
# Fit coverage: a row narrowed after its fit
#
# Narrowing a row's variants leaves its fit as it was made. Every reader of
# the fit says so, and nothing is rescaled.
# ===========================================================================

# `x` with every row narrowed to its first `k` variants, its fits untouched.
.alc_narrowed <- function(x, k) {
    endoapply(x, utils::head, k)
}

test_that("a fine-mapping row keeps its fit's variants when narrowed", {
    x <- .rc_makeFmr(.rc_vids(1:10))
    expect_equal(as.character(x$fitVariants[[1L]]), .rc_vids(1:10))
    narrowed <- .alc_narrowed(x, 4L)
    expect_length(narrowed[[1L]], 4L)
    expect_equal(as.character(narrowed$fitVariants[[1L]]), .rc_vids(1:10))
})

test_that("reading a narrowed row's fit warns and returns it unchanged", {
    x <- .rc_makeFmr(.rc_vids(1:10))
    narrowed <- .alc_narrowed(x, 4L)
    expect_warning(
        fit <- susieFit(narrowed),
        "made with 10 variants, 6 of which this row no longer holds"
    )
    expect_identical(fit, susieFit(x))
    expect_warning(
        credibleSetSummary(narrowed),
        "credibleSetSummary\\(\\): 1 row"
    )
    # The rows that hold every variant of their fit read silently.
    expect_no_warning(susieFit(x))
    expect_no_warning(lbf(x))
    expect_no_warning(credibleSetSummary(x))
})

test_that("lbf() of a narrowed row keeps the held variants' values as fitted", {
    x <- .rc_makeFmr(.rc_vids(1:10))
    expect_warning(out <- lbf(.alc_narrowed(x, 4L)), "lbf\\(\\): 1 row")
    full <- lbf(x)
    expect_equal(out$variant_id, .rc_vids(1:4))
    expect_equal(
        as.data.frame(mcols(out)),
        as.data.frame(mcols(full))[1:4, , drop = FALSE]
    )
})

test_that("the fSuSiE readers warn on a narrowed row too", {
    data(fsusieFineMappingExample, envir = environment())
    narrowed <- .alc_narrowed(fsusieFineMappingExample, 3L)
    expect_warning(fsusieCredibleBand(narrowed), "fsusieCredibleBand")
    expect_warning(fsusieTraitEffects(narrowed), "fsusieTraitEffects")
    expect_warning(fsusieAffectedRegions(narrowed), "fsusieAffectedRegions")
})

test_that("intersectVariants restricts a narrowed row's fit by id", {
    # The route the coverage warning points to. The kept variants are found
    # in the fit by id, so a row narrowed before reconciling restricts the
    # fit's own columns rather than failing on the width mismatch.
    x <- .rc_makeFmr(.rc_vids(1:30))
    y <- .rc_makeFmr(.rc_vids(5:50), studyName = "s2")
    narrowed <- endoapply(x, utils::tail, 20L)
    out <- intersectVariants(narrowed, y, oneSided = TRUE)
    expect_no_warning(fit <- susieFit(out))
    expect_equal(colnames(fit$alpha), .rc_vids(11:30))
    kept <- susieFit(x)$alpha[, .rc_vids(11:30)]
    expect_equal(unname(fit$alpha), unname(kept / rowSums(kept)))
    expect_equal(as.character(out$fitVariants[[1L]]), .rc_vids(11:30))
})

test_that("updateObject() backfills the variants each fit was made with", {
    x <- .rc_makeFmr(.rc_vids(1:10))
    saved <- x
    saved$fitVariants <- NULL
    expect_equal(
        as.character(updateObject(saved)$fitVariants[[1L]]),
        .rc_vids(1:10)
    )
    # A row narrowed before it was saved: its fit still names all ten.
    savedNarrowed <- .alc_narrowed(x, 4L)
    savedNarrowed$fitVariants <- NULL
    up <- updateObject(savedNarrowed)
    expect_equal(as.character(up$fitVariants[[1L]]), .rc_vids(1:10))
    expect_warning(susieFit(up), "no longer holds")
    # A current collection is left as it is.
    expect_identical(updateObject(x)$fitVariants, x$fitVariants)
    expect_message(updateObject(x, verbose = TRUE), "updateObject")
})

test_that("updateObject(): an unnamed fit records its row, no fit none", {
    vids <- .rc_vids(1:3)
    tl <- data.frame(variant_id = vids, pip = c(0.5, 0.3, 0.2))
    saved <- QtlFineMappingResult(
        studyName = c("s1", "s2"),
        context = c("c1", "c1"),
        trait = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(
            FineMappingRow(vids, list(pip = c(0.5, 0.3, 0.2)), tl),
            FineMappingRow(vids, NULL, tl)
        )
    )
    saved$fitVariants <- NULL
    up <- updateObject(saved)
    expect_equal(as.character(up$fitVariants[[1L]]), vids)
    expect_equal(as.character(up$fitVariants[[2L]]), character(0))
})

test_that(".fitNamedVariants reads a fit's names, never susieR's null column", {
    f <- pecotmr:::.fitNamedVariants
    expect_equal(f(NULL), character(0))
    expect_equal(f(list()), character(0))
    # A null_weight fit: one more effect-matrix column, named "null".
    alpha <- matrix(0.5, 1L, 3L, dimnames = list(NULL, c("v1", "v2", "null")))
    expect_equal(f(list(alpha = alpha)), c("v1", "v2"))
    expect_equal(
        f(list(pip = c(v1 = 0.1, v2 = 0.2), alpha = alpha)),
        c("v1", "v2")
    )
    lbfMat <- matrix(0, 1L, 2L, dimnames = list(NULL, c("a", "b")))
    expect_equal(f(list(lbf_variable = lbfMat)), c("a", "b"))
    # fSuSiE names its variables on csd_X.
    expect_equal(f(list(csd_X = c(x1 = 1, x2 = 2))), c("x1", "x2"))
})

# ===========================================================================
# adjustPips: entries that share no variant with `keepVariants`
#
# The zero-row short-circuit is covered above and the all-overlap path by
# test_fineMappingRow.R; these are the partial and disjoint branches.
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
    # As with range subsetting: an element that trims to zero variants
    # goes away rather than lingering unadjusted.
    expect_message(
        out <- adjustPips(.ac_twoEntry(), .rc_vids(1:5)),
        "dropping 1 of 2 entries"
    )
    expect_equal(nrow(out), 1L)
    expect_equal(as.character(mcols(out)$study), "s1")
})

test_that("adjustPips on a FineMappingResultBase collection renormalizes each entry", {
    vidsA <- paste0("chr1:", 1:6, ":A:G")
    vidsB <- paste0("chr1:", 3:8, ":A:G")
    entryA <- .makeAdjustEntry(vidsA)
    entryB <- .makeAdjustEntry(vidsB)
    fmr <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(entryA, entryB)
    )
    # Keep only variants shared by both entries' raw sets.
    keep <- intersect(vidsA, vidsB)
    adj <- adjustPips(fmr, keep)
    expect_s4_class(adj, "QtlFineMappingResult")
    expect_equal(nrow(adj), 2L)
    expect_equal(variantIds(pecotmr:::.collectionEntry(adj, 1L)), keep)
    expect_equal(variantIds(pecotmr:::.collectionEntry(adj, 2L)), keep)
})

# ===========================================================================
# .ssStitchElements()'s `window` filter
#
# The pipelines narrow a tuple to the span a pair's weights reach; the
# public sumStats() returns the whole tuple.
# ===========================================================================

test_that(".ssStitchElements() narrows the variants to a window", {
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
    sub <- pecotmr:::.ssStitchElements(gwasSumStatsS4Example, 1L, win)
    expect_s4_class(sub, "GRanges")
    expect_lt(length(sub), length(full))
    expect_gt(length(sub), 0L)
    expect_true(all(IRanges::overlapsAny(sub, win)))
})

test_that(".ssStitchElements() on another chromosome returns nothing", {
    # The chromosome test short-circuits before overlapsAny(), which would
    # otherwise warn about disjoint seqlevels.
    data(gwasSumStatsS4Example)
    off <- GenomicRanges::GRanges("chrZZ", IRanges::IRanges(1L, 100L))
    sub <- pecotmr:::.ssStitchElements(gwasSumStatsS4Example, 1L, off)
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
        representation(
            genome = "character",
            ldReference = "ANY",
            qcOptions = "list"
        )
    )
    empty <- new(
        "RcEmptyKid",
        GenomicRanges::GRangesList(),
        genome = "hg38",
        ldReference = NULL,
        qcOptions = list()
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
# Variant-level replacement: af<- / maf<- / sampleSize<- dispatch on `value`
# ===========================================================================

# A fillFromReference() reference over the same variants as `x`, with the
# alleles in the given orientation and one value column.
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

test_that("maf<- and sampleSize<- numeric are positional too", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    n <- length(variantIds(x))
    # The numeric (positional) half of the dispatch, for the two columns that
    # are swap-invariant and so have no allele-matching branch to exercise.
    maf(x) <- seq(0.05, 0.2, length.out = n)
    expect_equal(length(maf(x)), n)
    expect_equal(maf(x)[[1L]], 0.05)
    sampleSize(x) <- rep(1234, n)
    expect_equal(unique(sampleSize(x)), 1234)
    expect_error(maf(x) <- c(0.1, 0.2), "one value per variant")
    expect_error(sampleSize(x) <- 1L, "one value per variant")
})

test_that("maf returns NULL when MAF column is absent, vector when present", {
    obj_no_maf <- .qtlMakeOne(with_maf = FALSE)
    expect_null(maf(obj_no_maf))

    obj_maf <- .qtlMakeOne(n = 4, with_maf = TRUE)
    expect_equal(length(maf(obj_maf)), 4L)
})

test_that("fillFromReference warns about variants it could not place", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    gr <- .alc_replacementGr(x, payload = 0.4)
    # Two variants shifted far enough to match nothing. Silently dropping
    # them would hide a mis-joined reference table, so it warns and says
    # how many were ignored.
    extra <- GenomicRanges::shift(gr[1:2], 10000000L)
    expect_warning(x <- fillFromReference(x, c(gr, extra)), "were ignored")
    expect_equal(unique(af(x)), 0.4)
})

test_that("the variant-level setters work on a GwasSumStats as well", {
    data(gwasSumStatsS4Example)
    x <- gwasSumStatsS4Example
    n <- length(variantIds(x))
    # A GWAS collection is keyed by study alone, so the row lookup takes the
    # non-QtlSumStats branch and the QTL-only selectors are refused.
    af(x) <- rep(0.3, n)
    expect_equal(unique(af(x)), 0.3)
    expect_error(
        af(x, context = "brain") <- rep(0.3, n),
        "context"
    )
})

test_that("fillFromReference matches on (chrom, pos, A1, A2) in any order", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    gr <- .alc_replacementGr(x, payload = 0.25)
    x <- fillFromReference(x, gr[rev(seq_along(gr))])
    expect_equal(unique(af(x)), 0.25)
})

test_that("a flip complements AF but leaves MAF and N alone", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    swapped <- .alc_replacementGr(x, swap = TRUE, payload = 0.3)
    # An allele swap sends af -> 1 - af. maf is min(af, 1 - af) and the
    # sample size is a count, so both are swap-invariant -- complementing
    # them would corrupt them.
    expect_equal(unique(af(fillFromReference(x, swapped, "AF"))), 0.7)
    expect_equal(unique(maf(fillFromReference(x, swapped, "MAF"))), 0.3)
    expect_equal(unique(sampleSize(fillFromReference(x, swapped, "N"))), 0.3)
})

test_that("fillFromReference requires A1/A2 and says which column it read", {
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
    expect_error(fillFromReference(x, bare), "must carry A1 and A2")
    # Two value columns, neither named like the target: which one is meant
    # cannot be known, so it is refused rather than guessed.
    two <- .alc_replacementGr(x)
    S4Vectors::mcols(two)$extra <- 1
    expect_error(
        fillFromReference(x, two),
        "has no `AF` column.*the one column besides A1/A2 \\(got 2: v, extra"
    )
})

test_that("fillFromReference reads the column named like the target", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    # A panel carries several per-variant columns at once; the one named
    # like the target is read, whatever else rides along.
    panel <- .alc_replacementGr(x, payload = 0.1)
    S4Vectors::mcols(panel)$AF <- 0.6
    S4Vectors::mcols(panel)$MAF <- 0.4
    expect_equal(unique(af(fillFromReference(x, panel, "AF"))), 0.6)
    expect_equal(unique(maf(fillFromReference(x, panel, "MAF"))), 0.4)
})

test_that("fillFromReference sets NA for variants the reference omits", {
    data(qtlSumStatsExample)
    x <- qtlSumStatsExample
    n <- length(variantIds(x))
    partial <- .alc_replacementGr(x, k = 5L, payload = 0.42)
    expect_warning(x <- fillFromReference(x, partial), "were set NA")
    a <- af(x)
    # Unknown, not unchanged.
    expect_equal(sum(!is.na(a)), 5L)
    expect_equal(sum(is.na(a)), n - 5L)
    expect_equal(unique(a[!is.na(a)]), 0.42)
})

test_that("fillFromReference refuses a reference that matches nothing", {
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
    expect_error(fillFromReference(x, far), "no variant in `reference` matches")
})

test_that("fillFromReference checks what it is handed", {
    data(qtlSumStatsExample)
    data(gwasSumStatsS4Example)
    gr <- .alc_replacementGr(qtlSumStatsExample)
    expect_error(fillFromReference(sumStats(qtlSumStatsExample), gr), "x")
    expect_error(
        fillFromReference(qtlSumStatsExample, as.data.frame(gr)),
        "reference"
    )
    expect_error(
        fillFromReference(qtlSumStatsExample, gr, column = "BETA"),
        "must be one of"
    )
    # A GWAS collection is keyed by study alone, so a QTL selector is refused
    # exactly as the positional setter refuses it.
    expect_error(
        fillFromReference(gwasSumStatsS4Example, gr, context = "brain"),
        "context"
    )
})

# ===========================================================================
# A GWAS study split across rows: one study, several elements
# ===========================================================================

# gwasSumStatsS4Example's one study split into two LD blocks of 100 variants.
.alc_blockSplit <- function() {
    data(gwasSumStatsS4Example)
    g <- gwasSumStatsS4Example
    GwasSumStats(
        studyName = "trait1",
        entry = list(sumStats(g)),
        genome = "hg19",
        ldReference = ldReference(g)[[1L]],
        qcOptions = qcOptions(g),
        region = GenomicRanges::GRanges(
            "chr22",
            IRanges::IRanges(c(1L, 15680321L), c(15680320L, 51304566L))
        )
    )
}

test_that("a positional setter writes across every row of a split study", {
    x <- .alc_blockSplit()
    expect_equal(nrow(x), 2L)
    v <- seq(0.01, 0.5, length.out = 200L)
    af(x, studyName = "trait1") <- v
    expect_equal(nrow(x), 2L)
    # sumStats() order is the rows' order, so each row takes its own piece.
    expect_equal(af(x), v)
    expect_equal(S4Vectors::mcols(x[[2L]])$AF, v[101:200])
    expect_error(af(x) <- v[-1L], "one value per variant \\(200")
})

test_that("fillFromReference places values across a split study's rows", {
    x <- .alc_blockSplit()
    gr <- sumStats(x)
    rep <- gr
    S4Vectors::mcols(rep) <- S4Vectors::DataFrame(
        A1 = S4Vectors::mcols(gr)$A1,
        A2 = S4Vectors::mcols(gr)$A2,
        v = seq_len(200L) / 1000
    )
    # Every reference variant is in one row or the other, so neither the
    # unused nor the missed warning fires.
    expect_silent(x <- fillFromReference(x, rep[rev(seq_along(rep))], "MAF"))
    expect_equal(maf(x), seq_len(200L) / 1000)
    expect_silent(x <- fillFromReference(x, rep, "N"))
    expect_equal(sampleSize(x)[c(1L, 200L)], c(0.001, 0.2))
})

test_that("emptied rows of one study do not collide", {
    x <- .alc_blockSplit()
    # A PIP screen or filter() leaves each block an empty row; each keeps
    # its own block, so two of them are not duplicates.
    emptied <- combineGwasSumStats(
        GwasSumStats(
            studyName = "trait1",
            entry = list(sumStats(x)[0L]),
            genome = "hg19",
            blockId = "b1"
        ),
        GwasSumStats(
            studyName = "trait1",
            entry = list(sumStats(x)[0L]),
            genome = "hg19",
            blockId = "b2"
        )
    )
    expect_equal(unname(lengths(emptied)), c(0L, 0L))
    # The same block twice is still a duplicate.
    expect_error(
        combineGwasSumStats(emptied[1L], emptied[1L]),
        "\\(study, range\\) must be unique"
    )
    # Without a blockId, an empty row is keyed by its position; a non-empty
    # one keeps its span key.
    bare <- GenomicRanges::GRangesList(
        GenomicRanges::GRanges(),
        GenomicRanges::GRanges(),
        GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 2L))
    )
    S4Vectors::mcols(bare) <- S4Vectors::DataFrame(study = rep("a", 3L))
    keys <- pecotmr:::.keyEmptyElements(bare, c("", "", "chr1:1-2"))
    expect_false(anyDuplicated(keys) > 0L)
    expect_identical(keys[[3L]], "chr1:1-2")
})

test_that("a show() line names one genome build however many chromosomes", {
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
    x <- GwasSumStats(studyName = "g1", entry = list(gr), genome = "hg38")
    expect_identical(pecotmr:::.showGenomeBuild(x), "hg38")
    expect_identical(
        pecotmr:::.showGenomeBuild(GenomicRanges::GRangesList()),
        "NA"
    )
})

test_that(".inputsBuild carries the one build the inputs name", {
    hg19 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, width = 1L))
    GenomeInfoDb::genome(hg19) <- "hg19"
    unset <- GenomicRanges::GRanges("chr2", IRanges::IRanges(1L, width = 1L))
    f <- pecotmr:::.inputsBuild
    expect_identical(f(list(hg19, unset, NULL), "f"), "hg19")
    expect_identical(f(list(unset, NULL), "f"), NA_character_)
    # A plain list is read element by element.
    expect_identical(f(list(list(unset, hg19)), "f"), "hg19")
    hg38 <- GenomeInfoDb::`genome<-`(unset, value = "hg38")
    expect_error(
        f(list(hg19, hg38), "colocPipeline"),
        "colocPipeline: the inputs name different genome builds (hg19, hg38)",
        fixed = TRUE
    )
})

test_that(".namedBuilds reads a QtlDataset's contexts, a MultiStudy's parts", {
    tagged <- .sc_makeSe()
    GenomeInfoDb::genome(tagged) <- "hg38"
    qd <- QtlDataset(
        studyName = "s1",
        genotypes = .sc_makeGenotypeHandle(),
        phenotypes = list(brain = tagged)
    )
    expect_identical(pecotmr:::.namedBuilds(qd), "hg38")
    ms <- MultiStudyQtlDataset(
        qtlDatasets = list(s1 = qd),
        sumStats = qtlSumStatsExample
    )
    expect_setequal(pecotmr:::.namedBuilds(ms), c("hg38", "hg19"))
    # An unranged experiment places nothing, so it names no build.
    bare <- SummarizedExperiment::SummarizedExperiment(list(a = matrix(1)))
    expect_identical(pecotmr:::.namedBuilds(bare), character(0))
})

test_that(".oneBuild and .withGenome leave an unknown build unset", {
    mixed <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(1L, width = 1L)
    )
    GenomeInfoDb::genome(mixed) <- c("hg19", "hg38")
    expect_identical(pecotmr:::.oneBuild(mixed), NA_character_)
    expect_identical(pecotmr:::.withGenome(mixed, NA_character_), mixed)
    expect_null(pecotmr:::.withGenome(NULL, "hg19"))
    tagged <- pecotmr:::.withGenome(mixed, "hg19")
    expect_identical(unname(GenomeInfoDb::genome(tagged)), c("hg19", "hg19"))
})

test_that("a result's ranged views carry its genome build", {
    data(qtlFineMappingExample)
    data(qtlFineMappingLbfExample)
    x <- GenomeInfoDb::`genome<-`(qtlFineMappingExample, value = "hg19")
    build <- function(gr) unique(unname(GenomeInfoDb::genome(gr)))
    expect_identical(build(topLoci(x)), "hg19")
    expect_identical(build(credibleSets(x)), "hg19")
    expect_identical(build(marginalEffects(x)), "hg19")
    expect_identical(build(lbf(qtlFineMappingLbfExample)), NA_character_)
    expect_identical(
        build(lbf(GenomeInfoDb::`genome<-`(
            qtlFineMappingLbfExample,
            value = "hg19"
        ))),
        "hg19"
    )
})

# ===========================================================================
# Collection aggregates: lbf, credibleSetSummary and the fSuSiE accessors
# ===========================================================================

test_that("lbf aggregates over a collection, carrying identity + NA-filling ragged effects", {
    vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    mkEntry <- function(maxNumSingleEffects) {
        lbf <- matrix(
            seq_len(maxNumSingleEffects * length(vids)),
            nrow = maxNumSingleEffects,
            byrow = TRUE
        ) # L x p
        fit <- structure(
            list(lbf_variable = lbf, pip = rep(0.3, length(vids))),
            class = "susie"
        )
        FineMappingRow(
            variantIds = vids,
            susieFit = fit,
            topLoci = data.frame(
                variant_id = vids,
                pip = rep(0.3, length(vids)),
                stringsAsFactors = FALSE
            )
        )
    }
    # Two entries with different effect counts (L = 2 vs 3) exercise the NA-fill of
    # the ragged lbf_L<k> columns described in the collection method's contract.
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("brain", "blood"),
        trait = c("g", "g"),
        method = c("susie", "susie"),
        entry = list(mkEntry(2L), mkEntry(3L))
    )
    w <- lbf(res)
    expect_s4_class(w, "GRanges")
    expect_true(all(
        c("study", "context", "trait", "method", "variant_id") %in%
            names(S4Vectors::mcols(w))
    ))
    expect_equal(length(w), 6L) # 2 entries x 3 variants
    expect_setequal(unique(w$context), c("brain", "blood"))
    expect_true("lbf_L3" %in% names(S4Vectors::mcols(w)))
    expect_true(all(is.na(w$lbf_L3[w$context == "brain"]))) # L = 2 entry: no effect 3
    expect_true(all(!is.na(w$lbf_L3[w$context == "blood"]))) # L = 3 entry: present
})

test_that("credibleSetSummary aggregates across a collection with entry identity", {
    mkEntry <- function() {
        vn <- c("chr1:100:A:G", "chr1:200:C:T")
        tl <- data.frame(
            variant_id = vn,
            pip = c(0.7, 0.2),
            logBF = c(2.5, 0.3),
            cs_95 = c("susie_1", "susie_1"),
            cs_95_purity = c(0.9, 0.9),
            stringsAsFactors = FALSE
        )
        FineMappingRow(variantIds = vn, susieFit = list(x = 1), topLoci = tl)
    }
    res <- QtlFineMappingResult(
        studyName = c("s", "s"),
        context = c("brain", "blood"),
        trait = c("g", "g"),
        method = c("susie", "susie"),
        entry = list(mkEntry(), mkEntry())
    )
    s <- credibleSetSummary(res)
    expect_true(all(
        c("study", "context", "trait", "method", "cs") %in% names(s)
    ))
    expect_equal(nrow(s), 2L) # one CS per entry
    expect_setequal(unique(s$context), c("brain", "blood"))
})

test_that("fsusieCredibleBand on a result with no fSuSiE fit is empty ranges", {
    data(qtlFineMappingExample, envir = environment())
    cb <- fsusieCredibleBand(qtlFineMappingExample)
    expect_s4_class(cb, "GRanges")
    expect_length(cb, 0L)
})

test_that("fsusieCredibleBand + fsusieAffectedRegions aggregate across a collection", {
    skip_if_not_installed("fsusieR")
    skip_if_not_installed("wavethresh")
    fit <- .fsa_fits$none # one fit, read by both rows
    res <- QtlFineMappingResult(
        studyName = c("s", "s"),
        context = c("brain", "blood"),
        trait = c("g", "g"),
        method = c("fsusie", "fsusie"),
        entry = list(.fsa_entry(fit), .fsa_entry(fit))
    )

    cb <- fsusieCredibleBand(res)
    expect_true(all(
        c(
            "study",
            "context",
            "trait",
            "method",
            "cs",
            "effect",
            "lower",
            "upper"
        ) %in%
            names(S4Vectors::mcols(cb))
    ))
    expect_s4_class(cb, "GRanges")
    expect_gt(length(cb), 0)
    expect_setequal(unique(cb$context), c("brain", "blood"))

    gr <- fsusieAffectedRegions(res)
    expect_s4_class(gr, "GRanges")
    expect_gt(length(gr), 0)
    expect_true(all(
        c("study", "context", "trait", "cs", "purity", "direction") %in%
            names(S4Vectors::mcols(gr))
    ))
    expect_setequal(unique(S4Vectors::mcols(gr)$context), c("brain", "blood"))
})

test_that("fsusieAffectedRegions on a collection of non-fSuSiE entries is an empty GRanges", {
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(pip = 0.5),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    res <- QtlFineMappingResult(
        studyName = "s",
        context = "c",
        trait = "t",
        method = "susie",
        entry = list(e)
    )
    expect_equal(length(fsusieAffectedRegions(res)), 0L)
})
