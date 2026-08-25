# Tests for the ColocBoostResult class and its views.
#
# The class exists because ColocBoost's output cannot be expressed in
# ColocResult's schema: its unit is a SET of colocalized outcomes rather than a
# pair, and it computes no PP.H0-PP.H4 decomposition at all. Several tests
# below guard exactly that boundary.

# A minimal stand-in for a colocboost object, shaped like the real one:
# cos_details keyed by cos_id, purity as a SQUARE matrix over sets, and
# top variables as a data frame with set ids for rownames.
.cbr_fake <- function(
    ids = "cos1:y1_y2",
    outcomes = list(c("t1", "t2")),
    members = list(2L),
    variants = list("chr1:200:C:T"),
    npc = 0.9,
    nRegion = 4L,
    focal = FALSE
) {
    regionIds <- str_c("chr1:", seq_len(nRegion) * 100L, ":C:T")
    vcp <- set_names(seq_len(nRegion) / (nRegion * 2), regionIds)
    purity <- matrix(
        1,
        nrow = length(ids),
        ncol = length(ids),
        dimnames = list(ids, ids)
    )
    structure(
        list(
            cos_summary = tibble(
                cos_id = ids,
                focal_outcome = focal,
                top_variable = map_chr(variants, 1L),
                top_variable_vcp = rep(0.8, length(ids))
            ),
            vcp = vcp,
            cos_details = list(
                cos = list(
                    cos_index = set_names(members, ids),
                    cos_variables = set_names(variants, ids)
                ),
                cos_outcomes = list(
                    outcome_name = set_names(outcomes, ids)
                ),
                cos_vcp = set_names(
                    rep(list(as.numeric(vcp)), length(ids)),
                    ids
                ),
                cos_npc = set_names(rep(npc, length(ids)), ids),
                cos_min_npc_outcome = set_names(rep(npc, length(ids)), ids),
                cos_purity = list(min_abs_cor = purity),
                cos_top_variables = data.frame(
                    top_index = unlist(members),
                    top_variables = unlist(variants),
                    row.names = ids,
                    stringsAsFactors = FALSE
                )
            )
        ),
        class = "colocboost"
    )
}

.cbr_info <- function(names = c("t1", "t2")) {
    data.frame(
        name = names,
        context = str_c("ctx", seq_along(names)),
        trait = "GENE1",
        study = "study1",
        dataForm = "individual",
        stringsAsFactors = FALSE
    )
}

test_that("ColocBoostResult: one element per confidence set", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_s4_class(x, "ColocBoostResult")
    expect_equal(nrow(x), 1L)
    expect_equal(lengths(x), 1L)
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
    expect_equal(getColocPairs(x)$nOutcomes, 3L)
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
    expect_equal(nrow(getColocPairs(x)), 0L)
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
    expect_equal(getColocVariants(x)$variant_id, "chr1:200:C:T")
})

test_that("ColocBoostResult: reachable by range, not just identity", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    win <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 250))
    expect_equal(length(IRanges::subsetByOverlaps(x, win)), 1L)
    expect_equal(sum(lengths(subsetRegion(x, win))), 1L)
})

test_that("getColocBoostOutcomes: joins each outcome to its identity", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    out <- getColocBoostOutcomes(x)
    expect_equal(nrow(out), 2L)
    expect_setequal(out$outcome, c("t1", "t2"))
    expect_setequal(out$context, c("ctx1", "ctx2"))
    expect_true(all(out$trait == "GENE1"))
})

test_that("getRegionVcp spans the region, not just set members", {
    x <- ColocBoostResult(
        list(.cbr_fake(nRegion = 6L)),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_equal(length(getRegionVcp(x)), 6L)
    expect_equal(sum(lengths(x)), 1L)
})

test_that("as.data.frame returns the set-level view", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    df <- as.data.frame(x)
    expect_s3_class(df, "data.frame")
    expect_equal(df$cosId, "cos1:y1_y2")
    expect_equal(df$outcomes, "t1; t2")
})

test_that("separate_gwas runs are keyed by gwasStudy", {
    # The three analyses are not parallel: separate_gwas is a named list one
    # level deeper, so without gwasStudy its sets would be indistinguishable.
    raw <- list(
        xqtl_coloc = NULL,
        joint_gwas = NULL,
        separate_gwas = list(G1 = .cbr_fake(), G2 = .cbr_fake()),
        computing_time = list()
    )
    x <- .cbToResultObject(raw, .cbr_info())
    expect_equal(nrow(x), 2L)
    expect_equal(x$analysis, rep("separate_gwas", 2L))
    expect_setequal(x$gwasStudy, c("G1", "G2"))
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
        "max cos_npc", capture.output(show(x)), fixed = TRUE
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

.cbr_with_ucos <- function(ucosIds = "ucos1:y1", outcome = "t1") {
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
        list(.cbr_with_ucos()),
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
        list(.cbr_with_ucos()),
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
        list(.cbr_with_ucos()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    expect_output(show(x), "colocalized : 1")
    expect_output(show(x), "outcome-only: 1")
})

# ===========================================================================
# Pooled variants, and the fallbacks for a result missing its summary tables
# ===========================================================================

test_that("getColocVariants(pooled = TRUE) collapses to one row per variant", {
    # Noisy-OR only: ColocBoost's sets are not per-block slices of one signal,
    # so there is no mutually-exclusive sum stage the way coloc has.
    x <- ColocBoostResult(
        list(.cbr_fake(), .cbr_fake()),
        c("xqtl_coloc", "gwas_only"),
        outcomeInfo = .cbr_info()
    )
    long <- getColocVariants(x)
    pooled <- getColocVariants(x, pooled = TRUE)
    expect_lte(nrow(pooled), nrow(long))
    expect_equal(
        anyDuplicated(pooled[, c("analysis", "variant_id")]),
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

test_that("validity names missing identity and outcomeInfo columns", {
    x <- ColocBoostResult(
        list(.cbr_fake()),
        "xqtl_coloc",
        outcomeInfo = .cbr_info()
    )
    bad <- x
    mcols(bad)$analysis <- NULL
    expect_error(methods::validObject(bad), "missing columns: analysis")
})
