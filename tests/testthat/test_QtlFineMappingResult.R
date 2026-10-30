# Tests for R/QtlFineMappingResult.R

# === Tests migrated from test_s4Constructors.R (QtlFineMappingResult) ===

test_that("QtlFineMappingResult: builds a collection keyed by 4-tuple", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(3)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_equal(nrow(res), 2L)
    expect_identical(ldReference(res), list())
})


test_that("QtlFineMappingResult: validity does not recurse on key subset (#546)", {
    # The validity method builds a key-column data.frame to check tuple
    # uniqueness. Doing so via `object[, keyCols]` preserves the
    # QtlFineMappingResult class while dropping the required `entry` column;
    # older S4Vectors revalidates that intermediate and fails with
    # "missing columns: entry". Guard the reporter's exact scenario.
    e <- .sc_makeFineMappingRow(3)
    res <- QtlFineMappingResult(
        studyName = "s",
        context = "c",
        trait = "t",
        method = "qsusie",
        entry = list(e),
        ldReference = NULL
    )
    expect_s4_class(res, "QtlFineMappingResult")
    expect_true(validObject(res))

    # A bare key-column subset is itself an invalid QtlFineMappingResult;
    # validity must never re-run over it.
    sub <- S4Vectors::mcols(res)[,
        c("study", "context", "trait", "method"),
        drop = FALSE
    ]
    expect_false("entry" %in% names(S4Vectors::mcols(sub)))
})


test_that("QtlFineMappingResult: stores an LD sketch when supplied", {
    e <- .sc_makeFineMappingRow(3)
    gh <- .sc_makeGenotypeHandle()
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e),
        ldReference = gh
    )
    # The reference is a panel; the handle it reads through is its seed.
    refs <- ldReference(res)
    expect_named(refs, res$ldReference)
    expect_s4_class(refs[[1L]], "RangedSummarizedExperiment")
    expect_identical(pecotmr:::.ldRefHandle(refs[[1L]]), gh)
})


test_that("QtlFineMappingResult: errors on length mismatch", {
    e <- .sc_makeFineMappingRow(3)
    expect_error(
        QtlFineMappingResult(
            studyName = c("s1", "s2"),
            context = c("c1"),
            trait = c("t1"),
            method = c("susie"),
            entry = list(e)
        ),
        "same length"
    )
})


test_that("QtlFineMappingResult: validity rejects duplicate 4-tuples", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(3)
    expect_error(
        QtlFineMappingResult(
            studyName = c("s1", "s1"),
            context = c("c1", "c1"),
            trait = c("t1", "t1"),
            method = c("susie", "susie"),
            entry = list(e1, e2)
        ),
        "uniqueness violated"
    )
})


test_that("QtlFineMappingResult: fineMappingResult returns selected entry", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(3)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    picked <- fineMappingResult(
        res,
        studyName = "s1",
        context = "c2",
        trait = "t1",
        method = "susie"
    )
    # fineMappingResult() returns the single-row COLLECTION for that tuple
    # -- there is no detached entry object -- so compare what the row carries.
    expect_s4_class(picked, "QtlFineMappingResult")
    expect_equal(nrow(picked), 1L)
    expect_identical(variantIds(picked), variantIds(e2))
    expect_identical(susieFit(picked), susieFit(e2))
    expect_equal(pip(picked), .fmrRowPip(e2))
})


test_that("QtlFineMappingResult: fineMappingResult errors on missing tuple", {
    e <- .sc_makeFineMappingRow(3)
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    expect_error(
        fineMappingResult(
            res,
            studyName = "ghost",
            context = "c1",
            trait = "t1",
            method = "susie"
        ),
        "No entry for"
    )
})


test_that("QtlFineMappingResult: single-row collection allows omitting selectors", {
    e <- .sc_makeFineMappingRow(3)
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    picked <- fineMappingResult(res)
    expect_s4_class(picked, "QtlFineMappingResult")
    expect_equal(nrow(picked), 1L)
    expect_identical(variantIds(picked), .fmrPartsVariantIds(e))
    expect_identical(susieFit(picked), .fmrPartsSusieFit(e))
})


test_that("QtlFineMappingResult: joint columns absent by default", {
    e <- .sc_makeFineMappingRow(3)
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    expect_false("jointStudies" %in% pecotmr:::.tupleColumnNames(res))
    expect_false("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
    expect_false("jointTraits" %in% pecotmr:::.tupleColumnNames(res))
})


test_that("QtlFineMappingResult: accepts jointContexts column", {
    e <- .sc_makeFineMappingRow(3)
    # Univariate susie at c1 + the c1 slice of an mvsusie joint over (c1, c2):
    # both real context c1, distinguished by method and the jointContexts tag.
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("susie", "mvsusie"),
        entry = list(e, e),
        jointContexts = c(NA_character_, "c1;c2")
    )
    expect_true("jointContexts" %in% pecotmr:::.tupleColumnNames(res))
    expect_identical(res$jointContexts, c(NA_character_, "c1;c2"))
})


test_that("QtlFineMappingResult: jointStudies + jointTraits combine cleanly", {
    e <- .sc_makeFineMappingRow(3)
    # Three per-context rows sharing the real (s1, c1, t1) tuple, each a slice of
    # a different joint fit: univariate; a cross-study+trait mvsusieRss; a
    # cross-context mvsusie. The joint* tags carry each fit's co-fit membership.
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1", "s1"),
        context = c("c1", "c1", "c1"),
        trait = c("t1", "t1", "t1"),
        method = c("susie", "mvsusieRss", "mvsusie"),
        entry = list(e, e, e),
        jointStudies = c(NA_character_, "s1;s2", NA_character_),
        jointContexts = c(NA_character_, NA_character_, "c1;c2"),
        jointTraits = c(NA_character_, "t1;t2", NA_character_)
    )
    expect_equal(nrow(res), 3L)
    expect_identical(res$jointStudies, c(NA_character_, "s1;s2", NA_character_))
    expect_identical(res$jointTraits, c(NA_character_, "t1;t2", NA_character_))
})


test_that("QtlFineMappingResult: uniqueness distinguishes joint members", {
    e <- .sc_makeFineMappingRow(3)
    # Real scenario: context c1 participates in two different mvsusie joint fits
    # -- one over (c1, c2), one over (c1, c3) -- producing two c1 rows with the
    # same 4-tuple, kept distinct only by their jointContexts membership.
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c1"),
        trait = c("t1", "t1"),
        method = c("mvsusie", "mvsusie"),
        entry = list(e, e),
        jointContexts = c("c1;c2", "c1;c3")
    )
    expect_equal(nrow(res), 2L)
    # same 4-tuple AND same jointContexts -> duplicate
    expect_error(
        QtlFineMappingResult(
            studyName = c("s1", "s1"),
            context = c("c1", "c1"),
            trait = c("t1", "t1"),
            method = c("mvsusie", "mvsusie"),
            entry = list(e, e),
            jointContexts = c("c1;c2", "c1;c2")
        ),
        "uniqueness violated"
    )
})


test_that("QtlFineMappingResult: length-mismatched joint vector errors", {
    e <- .sc_makeFineMappingRow(3)
    expect_error(
        QtlFineMappingResult(
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "susie",
            entry = list(e),
            jointContexts = c("c1;c2", "c1;c3")
        ),
        "same length"
    )
})

# === Tests migrated from test_showMethods.R (QtlFineMappingResult) ===

test_that("show.QtlFineMappingResult prints entry/study/context/trait/method counts", {
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susieRss"),
        entry = list(.sh_makeFmEntry(), .sh_makeFmEntry())
    )
    out <- capture.output(show(res))
    expect_true(any(grepl("QtlFineMappingResult: 2 entries", out)))
    expect_true(any(grepl("1 studies.*2 contexts.*1 traits.*2 methods", out)))
    expect_true(any(grepl("LD reference: NULL", out)))
})


test_that("show.QtlFineMappingResult reports the ldReference source", {
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(.sh_makeFmEntry()),
        ldReference = .sh_makeGenotypeHandle()
    )
    out <- capture.output(show(res))
    expect_true(any(grepl("LD reference: plink2 @ /tmp/test", out)))
})


# === Tests migrated from test_collectionAccessors.R (QtlFineMappingResult) ===

test_that("QtlFineMappingResult: pip returns named pip vector for selected tuple", {
    e1 <- .ca_makeFmEntry(3)
    e2 <- .ca_makeFmEntry(4)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    pip <- pip(
        res,
        studyName = "s1",
        context = "c2",
        trait = "t1",
        method = "susie"
    )
    expect_equal(length(pip), 4L)
    expect_equal(names(pip)[1L], "chr1:100:A:G")
})


test_that("QtlFineMappingResult: pip() is one named numeric vector", {
    e <- .ca_makeFmEntry(3)
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    out <- pip(res)
    expect_type(out, "double")
    expect_length(out, 3L)
    expect_error(pip(res, returnList = TRUE), "unused")
})


test_that("QtlFineMappingResult: credibleSets aggregates every entry's credible sets with identity columns", {
    # e1 (n=3) has cs_95 = susie_1/susie_1/susie_0 -> 2 CS members; e2 (n=2) -> 2.
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    cs <- credibleSets(res)
    expect_equal(length(cs), 4L)
    expect_true(all(
        c("study", "context", "trait", "blockId", "method", "variant_id") %in%
            names(S4Vectors::mcols(cs))
    ))
    expect_equal(cs$context, c("c1", "c1", "c2", "c2"))
    # full tuple still returns the bare per-entry credible-set table
    bare <- credibleSets(
        res,
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie"
    )
    expect_false("context" %in% names(S4Vectors::mcols(bare)))
    expect_equal(length(bare), 2L)
})

test_that("QtlFineMappingResult: marginalEffects aggregates every entry with identity columns", {
    e1 <- .sc_makeFineMappingRow(3)
    e2 <- .sc_makeFineMappingRow(2)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    me <- marginalEffects(res)
    expect_equal(length(me), 5L) # all topLoci rows (no CS filter)
    expect_true(all(
        c("study", "context", "trait", "method") %in%
            names(S4Vectors::mcols(me))
    ))
    expect_equal(unique(me$context), c("c1", "c2"))
})


test_that("QtlFineMappingResult: susieFit reads the entry's trimmedFit", {
    e <- .ca_makeFmEntry(3)
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    expect_equal(susieFit(res), list(payload = "fit_n=3"))
})


test_that("QtlFineMappingResult: variantIds reads the entry's variantIds", {
    e <- .ca_makeFmEntry(3)
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(e)
    )
    expect_equal(length(variantIds(res)), 3L)
})


test_that("QtlFineMappingResult: studyName/contexts/traitNames/methodNames are unique", {
    e <- .ca_makeFmEntry(3)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1", "s2"),
        context = c("c1", "c2", "c1"),
        trait = c("t1", "t1", "t1"),
        method = c("susie", "susieRss", "susie"),
        entry = list(e, e, e)
    )
    expect_setequal(studyName(res), c("s1", "s2"))
    expect_setequal(contexts(res), c("c1", "c2"))
    expect_equal(traitNames(res), "t1")
    expect_setequal(methodNames(res), c("susie", "susieRss"))
})


test_that("cvResult works at the QtlFineMappingResult collection level", {
    tl <- data.frame(
        variant_id = "chr1:100:A:G",
        pip = 0.5,
        stringsAsFactors = FALSE
    )
    cv <- list(
        samplePartition = data.frame(Sample = "s1", Fold = 1L),
        prediction = list(susie_predicted = matrix(0, 1, 1)),
        performance = list(susie_performance = matrix(0, 1, 6))
    )
    e <- FineMappingRow("chr1:100:A:G", list(), tl, cvResult = cv)
    fmr <- QtlFineMappingResult(
        studyName = "S",
        context = "C",
        trait = "T",
        method = "susie",
        entry = list(e)
    )
    expect_identical(
        cvResult(
            fmr,
            studyName = "S",
            context = "C",
            trait = "T",
            method = "susie"
        ),
        cv
    )
})


test_that("QtlFineMappingResult: marginalEffects with tuple selectors", {
    e1 <- .ca_makeFmEntry(3)
    e2 <- .ca_makeFmEntry(4)
    res <- QtlFineMappingResult(
        studyName = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("t1", "t1"),
        method = c("susie", "susie"),
        entry = list(e1, e2)
    )
    # Collection-level selection picks the (s1, c2, t1, susie) entry, then
    # delegates to the entry-level getMarginalEffects.
    me <- marginalEffects(
        res,
        studyName = "s1",
        context = "c2",
        trait = "t1",
        method = "susie"
    )
    expect_s4_class(me, "GRanges")
    expect_equal(length(me), 4L)
    expect_true(all(
        c("variant_id", "beta", "se", "z", "p") %in% names(S4Vectors::mcols(me))
    ))
})

# ===========================================================================
# Validity: the branches that name what is wrong
#
# Every test above builds a valid collection, so validity's early returns
# were never executed. They are reachable by editing mcols on a built object,
# which is what a careless downstream edit would do.
# ===========================================================================

.qfmr_res <- function() {
    QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(.sc_makeFineMappingRow(3))
    )
}

test_that("validity names a missing identity column", {
    bad <- .qfmr_res()
    mcols(bad)$trait <- NULL
    expect_error(
        methods::validObject(bad),
        "missing elements \\{'trait'\\}"
    )
})

test_that("validity names a missing entry payload column", {
    bad <- .qfmr_res()
    mcols(bad)$susieFit <- NULL
    expect_error(
        methods::validObject(bad),
        "missing entry payload columns: .*missing elements \\{'susieFit'\\}"
    )
})

test_that("a joint column must be character, not a factor or number", {
    # The joint columns hold semicolon-joined member identities; a factor
    # would compare and paste differently and silently corrupt the tuple key.
    bad <- .qfmr_res()
    mcols(bad)$jointStudies <- 1L
    expect_error(
        methods::validObject(bad),
        "'jointStudies' column must be character \\(got integer\\)"
    )
})

test_that("a QtlFineMappingResult refuses a blockId", {
    data(qtlFineMappingExample)
    x <- qtlFineMappingExample
    # Block selection belongs to GwasFineMappingResult. The shared
    # FineMappingResultBase accessors pass `blockId` to every selector, and
    # this one refuses it rather than return the study/context/trait/method
    # match as if it had been honoured.
    expect_error(susieFit(x, blockId = "b1"), "has no LD blocks")
    expect_error(variantIds(x, blockId = "b1"), "has no LD blocks")
    # The supported selectors still work.
    expect_type(susieFit(x, method = "susie"), "list")
})

test_that("fineMappingResult(QtlFineMappingResult) refuses a blockId", {
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(.sc_makeFineMappingRow(3))
    )
    # A QTL fit is one row; `blockId` names a GWAS block, and NA is unset.
    expect_equal(nrow(fineMappingResult(res, blockId = NA)), 1L)
    expect_error(
        fineMappingResult(res, blockId = "chr1_1_500"),
        "`blockId` names an LD block of a GWAS fit"
    )
})

test_that("a QTL fit across chromosomes is refused by name", {
    ids <- c("chr1:100:A:G", "chr2:100:A:G")
    trans <- FineMappingRow(
        variantIds = ids,
        susieFit = list(),
        topLoci = data.frame(variant_id = ids, pip = c(0.6, 0.4))
    )
    expect_error(
        QtlFineMappingResult(
            studyName = c("s1", "s1"),
            context = c("c1", "c1"),
            trait = c("t1", "t2"),
            method = c("susie", "susie"),
            entry = list(.sc_makeFineMappingRow(3), trans)
        ),
        "`entry` element\\(s\\) 2 span more than one chromosome"
    )
})

test_that("QtlFineMappingResult() keeps the caller's metadata", {
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(.sc_makeFineMappingRow(3)),
        metadata = list(run = "r1")
    )
    expect_equal(S4Vectors::metadata(res), list(run = "r1"))
    expect_error(
        QtlFineMappingResult(
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "susie",
            entry = list(.sc_makeFineMappingRow(3)),
            metadata = "r1"
        ),
        "metadata"
    )
})

test_that("rows record no fit variants: a row is its fit's variables", {
    row <- .sc_makeFineMappingRow(3)
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(row)
    )
    expect_false("fitVariants" %in% colnames(S4Vectors::mcols(res)))
    expect_identical(susieFit(res), susieFit(row))
})

test_that("cvResult() of a narrowed row warns and returns it unchanged", {
    # A row is narrowed when it holds fewer variants than its fit has, so
    # it is the fit that tells.
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A"),
        pip = c(0.5, 0.3, 0.2)
    )
    cv <- list(prediction = list(susie = 1:4))
    res <- QtlFineMappingResult(
        studyName = "s1",
        context = "c1",
        trait = "t1",
        method = "susie",
        entry = list(FineMappingRow(
            tl$variant_id,
            list(pip = tl$pip),
            tl,
            cvResult = cv
        ))
    )
    narrowed <- endoapply(res, utils::head, 2L)
    expect_warning(
        out <- cvResult(narrowed),
        "cvResult\\(\\): the stored fit was made with 3 variants"
    )
    expect_identical(out, cv)
    expect_no_warning(cvResult(res))
})

# ===========================================================================
# Constructor guards
# ===========================================================================

test_that("QtlFineMappingResult: rejects rows that are not row payloads", {
    expect_error(
        QtlFineMappingResult(
            studyName = "s1",
            context = "c1",
            trait = "t1",
            method = "susie",
            entry = list("not_an_entry")
        ),
        "must be a fine-mapping row"
    )
})
