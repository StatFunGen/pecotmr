# Tests for the ColocBoostResult class and its views.
#
# The class exists because ColocBoost's output cannot be expressed in
# ColocResult's schema: its unit is a SET of colocalized outcomes rather than a
# pair, and it computes no PP.H0-PP.H4 decomposition at all. Several tests
# below guard exactly that boundary.

test_that("ColocBoostResult: one element per confidence set", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_s4_class(x, "ColocBoostResult")
    expect_equal(nrow(x), 1L)
    expect_equal(unname(lengths(x)), 1L)
    expect_true(x$isColocalized)
})

test_that("ColocBoostResult: a set holds an arbitrary number of outcomes", {
    # This is the reason the class is not ColocResult: coloc is pairwise, a
    # CoS is not. Three outcomes must live in ONE row.
    x <- ColocBoostResult(
        list(.cbr_fake(outcomes = list(c("t1", "t2", "t3")))),
        "xqtl_coloc",
        outcomeInfo = .cbr_info(c("t1", "t2", "t3"))
    )
    expect_equal(nrow(x), 1L)
    expect_equal(lengths(x$outcomes), 3L)
    expect_equal(colocPairs(x)$nOutcomes, 3L)
})

test_that("ColocBoostResult: shares a base with ColocResult", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_true(is(x, "ColocResultBase"))
    expect_true(is(x, "RangedTupleList"))
    # ... but is NOT a ColocResult: the schemas are not interchangeable.
    expect_false(is(x, "ColocResult"))
})

test_that("ColocBoostResult: rejects an element with no vcp layer", {
    # The per-variant layer is the point of the class; an element without it
    # must not construct. Built directly so the validity method is what
    # rejects it, rather than an earlier shape error.
    gr <- .variantIdsToGRanges("chr1:200:C:T", what = "variant_id")
    grl <- GenomicRanges::GRangesList(list(gr))
    md <- S4Vectors::DataFrame(
        cosId = "c1",
        analysis = "xqtl_coloc",
        gwasStudy = NA_character_,
        purity = 1,
        cosNpc = 1,
        minNpcOutcome = 1,
        nVariables = 1L,
        topVariable = "chr1:200:C:T",
        topVariableVcp = 1,
        focalOutcome = NA_character_,
        isColocalized = TRUE
    )
    md$outcomes <- IRanges::CharacterList(list("t1"))
    mcols(grl) <- md
    expect_error(
        new("ColocBoostResult", grl, outcomeInfo = .cbr_info()),
        "vcp"
    )
})

test_that("ColocBoostResult: rejects outcomes missing from outcomeInfo", {
    # A name that does not resolve would silently drop rows from the identity
    # join, which reads as "fewer traits colocalized" rather than as a bug.
    expect_error(
        ColocBoostResult(
            list(.cbr_fake(outcomes = list(c("t1", "ghost")))),
            "xqtl_coloc",
            outcomeInfo = .cbr_info()
        ),
        "outcomeInfo"
    )
})

test_that("ColocBoostResult: an empty result keeps the column schema", {
    x <- ColocBoostResult(list(), character(0), outcomeInfo = .cbr_info())
    expect_equal(nrow(x), 0L)
    expect_true(all(is_in(
        c("cosId", "analysis", "cosNpc", "isColocalized"),
        colnames(x)
    )))
    expect_equal(nrow(colocPairs(x)), 0L)
})

test_that("ColocBoostResult: NULL runs are skipped, not errors", {
    # .cbRun returns NULL for a failed analysis.
    x <- ColocBoostResult(
        list(NULL, .cbr_fake()),
        c("xqtl_coloc", "joint_gwas"),
        outcomeInfo = .cbr_info()
    )
    expect_equal(nrow(x), 1L)
    expect_equal(x$analysis, "joint_gwas")
})

test_that("ColocBoostResult: purity comes off the matrix diagonal", {
    # cos_purity is a square matrix over SETS; indexing it by id (which is what
    # it looks like it wants) silently yields NA.
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_equal(x$purity, 1)
})

test_that("ColocBoostResult: focalOutcome is character, never logical", {
    # colocboost returns `character` when a focal outcome is set and
    # `logical FALSE` when it is not; unnormalized, binding two analyses
    # coerces the column.
    noFocal <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    withFocal <- ColocBoostResult(
        list(.cbr_fake(focal = "t1")),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_type(noFocal$focalOutcome, "character")
    expect_true(is.na(noFocal$focalOutcome))
    expect_equal(withFocal$focalOutcome, "t1")
})

test_that("ColocBoostResult: variants keep alleles so ids round-trip", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_equal(colocVariants(x)$variant_id, "chr1:200:C:T")
})

test_that("ColocBoostResult: reachable by range, not just identity", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    win <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 250))
    expect_equal(length(IRanges::subsetByOverlaps(x, win)), 1L)
    expect_equal(
        sum(lengths(endoapply(x, IRanges::subsetByOverlaps, ranges = win))),
        1L
    )
})

test_that("colocBoostOutcomes: joins each outcome to its identity", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    out <- colocBoostOutcomes(x)
    expect_equal(nrow(out), 2L)
    expect_setequal(out$outcome, c("t1", "t2"))
    expect_setequal(out$context, c("ctx1", "ctx2"))
    expect_true(all(out$trait == "GENE1"))
})

test_that("regionVcp spans the region, not just set members", {
    x <- ColocBoostResult(
        list(.cbr_fake(nRegion = 6L)),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_equal(length(regionVcp(x)), 6L)
    expect_equal(sum(lengths(x)), 1L)
})

test_that("colocPairs is the set-level view; as.data.frame is standard", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    sets <- colocPairs(x)
    expect_s4_class(sets, "DataFrame")
    expect_equal(sets$cosId, "cos1:y1_y2")
    expect_equal(sets$outcomes, "t1; t2")
    # as.data.frame() is the inherited GRangesList coercion: one row per
    # (set, variant), placed by its range.
    df <- as.data.frame(x)
    expect_true(all(c("group", "seqnames", "start") %in% names(df)))
    expect_equal(nrow(df), sum(lengths(x)))
})

test_that(".cbrFocalOutcome maps the no-focal sentinel to NA", {
    expect_true(is.na(.cbrFocalOutcome(FALSE)))
    expect_true(is.na(.cbrFocalOutcome(NULL)))
    expect_equal(.cbrFocalOutcome("t3"), "t3")
})

test_that(".cbrPurity returns NA rather than guessing", {
    expect_true(is.na(.cbrPurity(list(), "x")))
    m <- matrix(0.7, 1, 1, dimnames = list("s1", "s1"))
    expect_equal(.cbrPurity(list(min_abs_cor = m), "s1"), 0.7)
    expect_true(is.na(.cbrPurity(list(min_abs_cor = m), "absent")))
})

# ===========================================================================
# show
# ===========================================================================

test_that("show summarizes sets, analyses and the best cos_npc", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_output(show(x), "ColocBoostResult with 1 confidence set\\(s\\)")
    expect_output(show(x), "colocalized :")
    expect_output(show(x), "outcome-only:")
    expect_output(show(x), "analyses    :")
    expect_output(show(x), "max cos_npc :")
})

test_that("show on an empty ColocBoostResult stops after the header", {
    x <- ColocBoostResult(list(), character(0), outcomeInfo = .cbr_info())
    expect_output(show(x), "with 0 confidence set\\(s\\)")
    expect_false(any(grepl(
        "max cos_npc",
        capture.output(show(x)),
        fixed = TRUE
    )))
})

# ===========================================================================
# Outcome-specific (uncolocalized) sets
#
# colocboost can report a credible set that belongs to ONE outcome -- "there
# is signal here, but it does not colocalize". The class keeps those as
# elements so that stays distinguishable from "no signal at all", but no
# fixture built one, leaving the whole ucos_details path unexecuted.
# ===========================================================================

.cbrWithUcos <- function(ucosIds = "ucos1:y1", outcome = "t1") {
    res <- .cbr_fake()
    vars <- list("chr1:300:C:T")
    purity <- matrix(
        1,
        nrow = length(ucosIds),
        ncol = length(ucosIds),
        dimnames = list(ucosIds, ucosIds)
    )
    res$ucos_details <- list(
        ucos = list(
            ucos_index = set_names(list(3L), ucosIds),
            ucos_variables = set_names(vars, ucosIds)
        ),
        ucos_outcomes = list(
            outcome_name = set_names(list(outcome), ucosIds)
        ),
        ucos_purity = list(min_abs_cor = purity),
        ucos_top_variables = data.frame(
            top_index = 3L,
            top_variables = unlist(vars),
            row.names = ucosIds,
            stringsAsFactors = FALSE
        )
    )
    res
}

test_that("an outcome-specific set becomes its own uncolocalized element", {
    x <- ColocBoostResult(
        list(.cbrWithUcos()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_equal(nrow(x), 2L)
    expect_equal(sum(x$isColocalized), 1L)
    expect_equal(sum(!x$isColocalized), 1L)
})

test_that("an uncolocalized set carries no colocalization statistics", {
    # cosNpc / vcp describe agreement BETWEEN outcomes, so a single-outcome
    # set has none to report -- NA rather than a misleading zero.
    x <- ColocBoostResult(
        list(.cbrWithUcos()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    md <- mcols(x, use.names = FALSE)
    u <- which(!md$isColocalized)
    expect_true(is.na(md$cosNpc[u]))
    expect_true(is.na(md$minNpcOutcome[u]))
    expect_true(all(is.na(mcols(x[[u]])$vcp)))
    expect_equal(md$nVariables[u], 1L)
})

test_that("show counts colocalized and outcome-only sets separately", {
    x <- ColocBoostResult(
        list(.cbrWithUcos()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_output(show(x), "colocalized : 1")
    expect_output(show(x), "outcome-only: 1")
})

# ===========================================================================
# Pooled variants, and the fallbacks for a result missing its summary tables
# ===========================================================================

test_that("colocPooledVariants collapses to one row per variant", {
    # Noisy-OR only: ColocBoost's sets are not per-block slices of one signal,
    # so there is no mutually-exclusive sum stage the way coloc has.
    x <- ColocBoostResult(
        list(.cbr_fake(), .cbr_fake()),
        c("xqtl_coloc", "gwas_only"),
        outcomeInfo = .cbr_info()
    )
    long <- colocVariants(x)
    pooled <- colocPooledVariants(x)
    expect_lte(length(pooled), length(long))
    pooledCols <- as.data.frame(S4Vectors::mcols(pooled))
    expect_equal(
        anyDuplicated(pooledCols[, c("analysis", "variant_id")]),
        0L
    )
    expect_true(all(pooled$vcp >= 0 & pooled$vcp <= 1))
})

test_that("a set with no cos_summary row falls back to NA descriptors", {
    # colocboost need not report a summary row for every set; the class must
    # still place the set rather than drop it.
    res <- .cbr_fake()
    res$cos_summary <- res$cos_summary[0, , drop = FALSE]
    x <- ColocBoostResult(list(res), "xqtl_coloc", outcomeInfo = .cbr_info())
    md <- mcols(x, use.names = FALSE)
    expect_equal(nrow(x), 1L)
    # topVariable comes from cos_top_variables, which is still present; the
    # summary is what supplies the vcp and the focal-outcome flag.
    expect_true(is.na(md$topVariableVcp))
    expect_true(is.na(md$focalOutcome))
})

test_that(".cbrTopVariable is NA when the id is absent from the table", {
    tv <- data.frame(
        top_index = 1L,
        top_variables = "chr1:100:C:T",
        row.names = "cos1",
        stringsAsFactors = FALSE
    )
    expect_equal(.cbrTopVariable(tv, "cos1"), "chr1:100:C:T")
    expect_true(is.na(.cbrTopVariable(tv, "nope")))
    expect_true(is.na(.cbrTopVariable(NULL, "cos1")))
})

test_that("validity names a missing identity column", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    bad <- x
    mcols(bad)$analysis <- NULL
    expect_error(
        methods::validObject(bad),
        "missing elements \\{'analysis'\\}"
    )
})


# ===========================================================================
# Empty-input paths
#
# The views and row builders were only ever driven by populated results, so
# their empty answers were never produced.
# ===========================================================================

test_that("uncolocalized rows are empty when there is nothing to report", {
    # Two separate exits: no ucos_details at all, and details carrying an
    # empty index.
    expect_equal(.cbrUncolocalizedRows(list(), "xqtl", NA_character_), list())
    expect_equal(
        .cbrUncolocalizedRows(
            list(ucos_details = list(ucos = list(ucos_index = list()))),
            "xqtl",
            NA_character_
        ),
        list()
    )
})

test_that("a summary lookup that matches no set yields an NA row", {
    # NA rather than zero rows: the set exists, its summary just is not there.
    out <- .cbrSummaryFor(
        list(cos_summary = data.frame(cos_id = "other")),
        "missing_id"
    )
    expect_equal(nrow(out), 1L)
    expect_true(is.na(out$top_variable))
    expect_true(is.na(out$top_variable_vcp))
})

test_that("the views answer an empty collection with empty tables", {
    data(colocboostResultExample)
    empty <- colocboostResultExample[0]
    expect_equal(length(empty), 0L)
    expect_equal(nrow(.cbrLongVariants(empty)), 0L)
    expect_equal(nrow(colocBoostOutcomes(empty)), 0L)
    expect_length(colocPooledVariants(empty), 0L)
    expect_s4_class(colocBoostOutcomes(empty), "DataFrame")
    expect_s4_class(colocPairs(colocboostResultExample), "DataFrame")
})

test_that("outcome info defaults to the empty frame and is not joined", {
    # With no outcomeInfo the constructor substitutes the empty frame, and the
    # outcomes view returns before attempting a join against it.
    cb <- ColocBoostResult(list(.cbr_fake()), analysis = "xqtl")
    expect_equal(nrow(cb@outcomeInfo), 0L)
    out <- colocBoostOutcomes(cb)
    expect_gt(nrow(out), 0L)
    expect_false(is_in("outcomeTrait", colnames(out)))
})


test_that("validity names the outcomeInfo columns that are missing", {
    cb <- ColocBoostResult(list(.cbr_fake()), analysis = "xqtl")
    cb@outcomeInfo <- data.frame(name = "t1")
    expect_match(
        .cbrCheckOutcomeInfo(cb),
        paste0(
            "outcomeInfo is missing columns: .*",
            "missing elements \\{'study','context','trait','dataForm'\\}"
        )
    )
})

# A small LdMatrix to stand in for one side's reference.
.cbr_ldm <- function(ids = c("chr1:100:C:T", "chr1:200:A:G")) {
    LdMatrix(matrix(c(1, 0.2, 0.2, 1), 2, 2, dimnames = list(ids, ids)))
}

test_that("ColocBoostResult keeps a named LD reference per side", {
    refs <- list(
        qtl = .cbr_ldm(),
        gwas = .cbr_ldm(c("chr1:300:G:A", "chr1:400:T:C"))
    )
    info <- mutate(
        .cbr_info(),
        dataForm = "sumstats",
        ldReference = c("qtl", "gwas")
    )
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "joint_gwas",
        outcomeInfo = info,
        ldReference = refs
    )
    expect_identical(ldReference(x), refs)
    out <- colocBoostOutcomes(x)
    expect_equal(
        out$ldReference[match(c("t1", "t2"), out$outcome)],
        c("qtl", "gwas")
    )
})

test_that("ColocBoostResult with no LD reference holds an empty list", {
    # An outcomeInfo without the column gets it, NA throughout.
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_identical(ldReference(x), list())
    expect_true(all(is.na(colocBoostOutcomes(x)$ldReference)))
})

test_that("ColocBoostResult refuses a lone LD reference", {
    # Its name is what outcomeInfo keys on, so it is not made up.
    expect_error(
        ColocBoostResult(
            list(.cbr_fake()),
            "xqtl_coloc",
            outcomeInfo = .cbr_info(),
            ldReference = .cbr_ldm()
        ),
        "named list of LD references"
    )
})

test_that("ColocBoostResult refuses an outcome keyed to a missing reference", {
    info <- mutate(.cbr_info(), ldReference = c("qtl", "gwas"))
    expect_error(
        ColocBoostResult(
            list(.cbr_fake()),
            "joint_gwas",
            outcomeInfo = info,
            ldReference = list(qtl = .cbr_ldm())
        ),
        "not in ldReference: gwas"
    )
})

test_that("validity wants uniquely named LD references", {
    # The shared check every class with per-row references validates with.
    problems <- function(refs) .ldRefDictProblems(refs, character(0))
    expect_match(problems(list(.cbr_ldm())), "unique, non-empty names")
    expect_match(
        problems(list(a = .cbr_ldm(), a = .cbr_ldm())),
        "unique, non-empty names"
    )
    expect_match(
        problems(list(a = .cbr_ldm(), b = "panel")),
        "not LD references: b"
    )
    expect_null(problems(list(a = .cbr_ldm())))
})

test_that("the region vcp skips runs that cannot supply one", {
    # Two ways a run is unusable: no vcp at all, and a vcp whose names are
    # missing so the variants cannot be placed.
    unnamed <- c(0.5, 0.5)
    expect_null(.cbrRegionVcp(list(list(vcp = unnamed))))
    expect_null(.cbrRegionVcp(list(list(vcp = NULL))))
    # ...and a later run that does carry one is still found.
    gr <- .cbrRegionVcp(list(
        list(vcp = NULL),
        list(vcp = set_names(0.7, "chr1:100:C:T"))
    ))
    expect_s4_class(gr, "GRanges")
    expect_equal(length(gr), 1L)
})

test_that("updateObject keys an old result's sumstat outcomes to its ref", {
    ref <- .cbr_ldm()
    info <- mutate(.cbr_info(), dataForm = c("sumstats", "individual"))
    x <- ColocBoostResult(list(.cbr_fake()), "joint_gwas", outcomeInfo = info)
    # As an earlier version saved it: one reference, no outcome keys.
    old <- `slot<-`(x, "ldReference", check = FALSE, value = ref)
    old <- `slot<-`(
        old,
        "outcomeInfo",
        value = select(info, -any_of("ldReference"))
    )
    new <- updateObject(old)
    expect_true(methods::validObject(new))
    expect_identical(ldReference(new), list(ld1 = ref))
    expect_identical(new@outcomeInfo$ldReference, c("ld1", NA))
    # A current result is left as it is.
    expect_identical(updateObject(x), x)
})

# ===========================================================================
# computingTime(): one row per analysis that ran
# ===========================================================================

test_that("computingTime() tabulates each analysis: its runs and seconds", {
    timing <- list(
        Analysis = list(
            xqtl_coloc = as.difftime(90, units = "secs"),
            joint_gwas = NULL,
            separate_gwas = list(
                total = as.difftime(2, units = "mins"),
                n_studies = 3L
            )
        )
    )
    x <- ColocBoostResult(
        list(),
        character(0),
        outcomeInfo = .cbr_info(),
        computingTime = timing
    )
    tab <- computingTime(x)
    expect_s4_class(tab, "DataFrame")
    # An analysis that did not run is no row.
    expect_equal(tab$analysis, c("xqtl_coloc", "separate_gwas"))
    expect_equal(tab$nRuns, c(1L, 3L))
    expect_equal(tab$seconds, c(90, 120))
})

test_that("computingTime() has no rows when nothing was timed", {
    tab <- computingTime(
        ColocBoostResult(list(), character(0), outcomeInfo = .cbr_info())
    )
    expect_equal(nrow(tab), 0L)
    expect_equal(colnames(tab), c("analysis", "nRuns", "seconds"))
})

test_that("a timing as a bare number, or without a total, still reads", {
    expect_equal(pecotmr:::.cbrTimingSeconds(12.5), 12.5)
    expect_true(is.na(pecotmr:::.cbrTimingSeconds(list(n_studies = 2L))))
    expect_true(is.na(pecotmr:::.cbrTimingRuns(list(total = 3))))
})

test_that("ColocBoostResult() keeps the caller's metadata", {
    x <- ColocBoostResult(
        list(),
        character(0),
        outcomeInfo = .cbr_info(),
        metadata = list(run = "r1")
    )
    expect_equal(S4Vectors::metadata(x), list(run = "r1"))
})
