context("CtwasResult")

# CtwasResultEntry (the per-row payload class) is tested in
# test_CtwasResultEntry.R; this file covers the CtwasResult collection.

.cr_entry <- function(ids = c("g1", "g2"), pip = c(0.9, 0.1), prior = NULL) {
    # Ranged, as the pipeline stores them: each id placed at a position.
    at <- IRanges::IRanges(100L * seq_along(ids), width = 1L)
    CtwasResultEntry(
        posteriors = GenomicRanges::GRanges(
            "chr1",
            at,
            id = ids,
            susie_pip = pip
        ),
        susieAlpha = GenomicRanges::GRanges(
            "chr1",
            at,
            id = ids,
            susie_alpha = pip
        ),
        groupPriors = prior,
        regionInfo = data.frame(region_id = "r1", stringsAsFactors = FALSE)
    )
}

test_that("CtwasResult rows are keyed by (gwasStudy, study, context, method)", {
    cr <- CtwasResult(
        gwasStudy = c("D1", "D1"),
        studyName = c("Q1", "Q1"),
        context = c("brain", "liver"),
        method = c("susie", "susie"),
        entry = list(.cr_entry(), .cr_entry(ids = "g3", pip = 0.5))
    )
    expect_s4_class(cr, "CtwasResult")
    expect_equal(nrow(cr), 2L)
    expect_equal(methodNames(cr), "susie")
    expect_setequal(contexts(cr), c("brain", "liver"))
    expect_equal(studyName(cr), "Q1")
})

test_that("CtwasResult: c() appends runs, their tables kept", {
    run <- function(context, entry) {
        CtwasResult(
            gwasStudy = "D1",
            studyName = "Q1",
            context = context,
            method = "susie",
            entry = list(entry)
        )
    }
    brain <- run("brain", .cr_entry())
    liver <- run("liver", .cr_entry(ids = "g3", pip = 0.5))
    both <- c(brain, liver)
    expect_s4_class(both, "CtwasResult")
    expect_identical(as.character(both$context), c("brain", "liver"))
    expect_identical(both$regionInfo, c(brain$regionInfo, liver$regionInfo))
    expect_identical(length(ctwasPosteriors(both)), 3L)
})

test_that("ctwasPosteriors aggregates rows tagged with run identity", {
    cr <- CtwasResult(
        gwasStudy = c("D1", "D1"),
        studyName = c("Q1", "Q1"),
        context = c("brain", "liver"),
        method = c("susie", "susie"),
        entry = list(
            .cr_entry(ids = c("g1", "g2")),
            .cr_entry(ids = "g3", pip = 0.5)
        )
    )
    fm <- ctwasPosteriors(cr)
    expect_s4_class(fm, "GRanges")
    expect_true(all(
        c("gwasStudy", "study", "context", "method", "id", "susie_pip") %in%
            names(S4Vectors::mcols(fm))
    ))
    expect_equal(length(fm), 3L) # 2 + 1 gene rows
    expect_equal(fm$context, c("brain", "brain", "liver"))
})

test_that("susieAlpha aggregates the per-effect table with identity", {
    cr <- CtwasResult(
        gwasStudy = c("D1", "D1"),
        studyName = c("Q1", "Q1"),
        context = c("brain", "liver"),
        method = c("susie", "susie"),
        entry = list(
            .cr_entry(ids = c("g1", "g2")),
            .cr_entry(ids = "g3", pip = 0.5)
        )
    )
    sa <- susieAlpha(cr)
    expect_equal(length(sa), 3L)
    expect_true(all(
        c("gwasStudy", "study", "context", "method", "susie_alpha") %in%
            names(S4Vectors::mcols(sa))
    ))
    expect_equal(sa$context, c("brain", "brain", "liver"))
})

test_that("CtwasResult: uniqueness includes joint columns", {
    # same (gwasStudy, study, context, method) collides without joint columns
    expect_error(
        CtwasResult(
            gwasStudy = c("D1", "D1"),
            studyName = c("Q1", "Q1"),
            context = c("brain", "brain"),
            method = c("susie", "susie"),
            entry = list(.cr_entry(), .cr_entry())
        ),
        "uniqueness violated"
    )
    # a single-context row and a multi-context (jointContexts) row over the same
    # context are distinct once jointContexts joins the key
    cr <- CtwasResult(
        gwasStudy = c("D1", "D1"),
        studyName = c("Q1", "Q1"),
        context = c("brain", "brain"),
        method = c("susie", "susie"),
        entry = list(.cr_entry(), .cr_entry()),
        jointContexts = c(NA, "brain,liver")
    )
    expect_equal(nrow(cr), 2L)
})

test_that("CtwasResult: rejects a non-CtwasResultEntry payload", {
    expect_error(
        CtwasResult(
            gwasStudy = "D1",
            studyName = "Q1",
            context = "brain",
            method = "susie",
            entry = list("not an entry")
        ),
        "must be a CtwasResultEntry"
    )
})

test_that("CtwasResult: constructor rejects mismatched core-vector lengths", {
    expect_error(
        CtwasResult(
            gwasStudy = c("D1", "D1"),
            studyName = "Q1",
            context = "brain",
            method = "susie",
            entry = list(.cr_entry())
        ),
        "same length"
    )
})

test_that("CtwasResult: a joint column of the wrong length is rejected", {
    expect_error(
        CtwasResult(
            gwasStudy = "D1",
            studyName = "Q1",
            context = "brain",
            method = "susie",
            entry = list(.cr_entry()),
            jointContexts = c("a", "b")
        ),
        "same length as"
    )
})

test_that("CtwasResult: the views are empty ranges when no row has a payload", {
    empty <- CtwasResult(
        gwasStudy = character(0),
        studyName = character(0),
        context = character(0),
        method = character(0),
        entry = list()
    )
    expect_length(ctwasPosteriors(empty), 0L)
    expect_s4_class(susieAlpha(empty), "GRanges")
    nullE <- CtwasResult(
        gwasStudy = "D1",
        studyName = "Q1",
        context = "brain",
        method = "susie",
        entry = list(CtwasResultEntry(posteriors = NULL, susieAlpha = NULL))
    )
    expect_length(ctwasPosteriors(nullE), 0L)
    expect_length(susieAlpha(nullE), 0L)
})

test_that("CtwasResult: show() prints a one-line-per-run summary", {
    cr <- CtwasResult(
        gwasStudy = "D1",
        studyName = "Q1",
        context = "brain",
        method = "susie",
        entry = list(.cr_entry())
    )
    expect_output(show(cr), "CtwasResult: 1 run")
    expect_output(show(cr), "susie")
})

# ===========================================================================
# Validity: the branches that name what is wrong
#
# A CtwasResult's identity and payload columns are mcols, so a downstream edit
# can drop or retype one. Validity guards that, but every test above built a
# valid object.
# ===========================================================================

.ctr_res <- function() {
    CtwasResult(
        gwasStudy = "D1",
        studyName = "Q1",
        context = "brain",
        method = "susie",
        entry = list(.cr_entry())
    )
}

test_that("validity names a missing required column", {
    bad <- .ctr_res()
    bad$context <- NULL
    expect_error(
        methods::validObject(bad),
        "missing elements \\{'context'\\}"
    )
})

test_that("a joint column must be character", {
    bad <- .ctr_res()
    S4Vectors::mcols(bad)$jointStudies <- 1L
    expect_error(
        methods::validObject(bad),
        "'jointStudies' column must be character \\(got integer\\)"
    )
})


test_that("a payload with no coordinates is refused at construction", {
    # An entry built by hand from a plain table has nowhere to place its rows
    # in a ranged collection, so the constructor says so instead of dropping
    # them.
    tableOnly <- CtwasResultEntry(
        posteriors = data.frame(id = "g1", susie_pip = 0.5)
    )
    expect_error(
        CtwasResult(
            gwasStudy = "D1",
            studyName = "Q1",
            context = "brain",
            method = "susie",
            entry = list(tableOnly)
        ),
        "run 1's posteriors carry no genomic coordinates"
    )
    # The entry itself still reads back as given.
    expect_s3_class(ctwasPosteriors(tableOnly), "data.frame")
})

# ===========================================================================
# One element per cTWAS region
# ===========================================================================

# A run over two cTWAS regions on two chromosomes, with its region metadata.
.cr_twoRegionEntry <- function() {
    regions <- c("chr1_1_1000", "chr2_1_1000")
    pos <- GenomicRanges::GRanges(
        c("chr1", "chr1", "chr2"),
        IRanges::IRanges(c(100L, 200L, 300L), width = 1L),
        id = c("g1", "g2", "g3"),
        region_id = regions[c(1L, 1L, 2L)],
        susie_pip = c(0.9, 0.1, 0.5)
    )
    sa <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(c(100L, 300L), width = 1L),
        id = c("g1", "g3"),
        region_id = regions,
        susie_alpha = c(0.8, 0.4)
    )
    ri <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(1L, 999L),
        region_id = regions
    )
    CtwasResultEntry(
        posteriors = pos,
        susieAlpha = sa,
        groupPriors = list(group_prior = c(gene = 0.01)),
        regionInfo = ri
    )
}

test_that("a run splits into one element per cTWAS region", {
    cr <- CtwasResult(
        gwasStudy = "D1",
        studyName = "Q1",
        context = "brain",
        method = "susie",
        entry = list(.cr_twoRegionEntry())
    )
    expect_s4_class(cr, "RangedTupleList")
    expect_equal(length(cr), 2L)
    expect_equal(cr$blockId, c("chr1_1_1000", "chr2_1_1000"))
    expect_equal(unname(lengths(cr)), c(2L, 1L))
    # Each element sits on one chromosome, so the range machinery applies.
    expect_equal(unname(lengths(range(cr))), c(1L, 1L))
    # The region's own payloads ride with it; the run's priors repeat.
    expect_equal(cr$susieAlpha[[2L]]$id, "g3")
    expect_equal(cr$regionInfo[[1L]]$region_id, "chr1_1_1000")
    expect_identical(cr$groupPriors[[1L]], cr$groupPriors[[2L]])
    # Named by the run identity plus the region.
    expect_equal(names(cr)[[1L]], "Q1|brain|susie|D1|chr1_1_1000")
    # The stacked views put every region back together.
    expect_equal(ctwasPosteriors(cr)$id, c("g1", "g2", "g3"))
    expect_equal(susieAlpha(cr)$id, c("g1", "g3"))
    # And a region narrows like any ranged collection.
    one <- IRanges::subsetByOverlaps(
        cr,
        GenomicRanges::GRanges("chr2", IRanges::IRanges(1L, 1000L))
    )
    expect_equal(one$blockId, "chr2_1_1000")
    expect_s4_class(one, "CtwasResult")
})

# ===========================================================================
# One set of group priors per run
# ===========================================================================

test_that("a run's rows must carry one set of group priors", {
    run <- function(chrom, prior, method = "susie") {
        gr <- GenomicRanges::GRanges(
            chrom,
            IRanges::IRanges(100L, width = 1L),
            id = "g1",
            susie_pip = 0.9
        )
        CtwasResult(
            gwasStudy = "D1",
            studyName = "Q1",
            context = "brain",
            method = method,
            entry = list(CtwasResultEntry(
                posteriors = gr,
                groupPriors = list(group_prior = c(gene = prior))
            ))
        )
    }
    # Fitted chromosome by chromosome: one identity, two fits.
    expect_error(
        c(run("chr1", 0.01), run("chr2", 0.05)),
        "\\(D1, Q1, brain, susie\\) carry different group priors"
    )
    # Under distinct identities they are two runs.
    both <- c(run("chr1", 0.01), run("chr2", 0.05, method = "lasso"))
    expect_equal(methodNames(both), c("susie", "lasso"))
    # One fit's regions taken apart and put back together are still one run.
    whole <- CtwasResult(
        gwasStudy = "D1",
        studyName = "Q1",
        context = "brain",
        method = "susie",
        entry = list(.cr_twoRegionEntry())
    )
    expect_equal(length(c(whole[1L], whole[2L])), 2L)
})

test_that("ctwasGroupPriors() tabulates each run's priors once, by group", {
    prior <- list(
        group_prior = c(SNP = 1e-4, brain = 0.01),
        group_prior_var = c(SNP = 10, brain = 20),
        group_size = c(SNP = 5000, brain = 30),
        converged = TRUE
    )
    cr <- CtwasResult(
        gwasStudy = c("D1", "D1"),
        studyName = c("Q1", "Q1"),
        context = c("brain", "liver"),
        method = c("susie", "susie"),
        entry = list(.cr_entry(prior = prior), .cr_entry(ids = "g3", pip = 0.5))
    )
    tab <- ctwasGroupPriors(cr)
    expect_s3_class(tab, "tbl_df")
    expect_named(
        tab,
        c(
            "gwasStudy",
            "study",
            "context",
            "method",
            "group",
            "group_study",
            "group_context",
            "group_type",
            "group_size",
            "group_prior",
            "group_prior_var",
            "log_enrichment",
            "log_enrichment_se",
            "enrichment_p",
            "group_pve",
            "prop_heritability"
        )
    )
    # The liver run recorded no priors, so it adds no rows.
    expect_equal(tab$context, c("brain", "brain"))
    expect_equal(tab$group, c("SNP", "brain"))
    expect_equal(tab$group_size, c(5000, 30))
    expect_equal(tab$group_prior, c(1e-4, 0.01))
    expect_equal(tab$group_prior_var, c(10, 20))
    # A hand-built run recorded neither the groups' identities nor a summary.
    expect_true(all(is.na(tab$group_study)))
    expect_true(all(is.na(tab$group_pve)))
    # No run with priors: no rows, the same columns.
    expect_equal(ctwasGroupPriors(cr[2L]), tab[0L, ])
    # A run over two regions is listed once; what it did not record is NA.
    two <- CtwasResult(
        gwasStudy = "D1",
        studyName = "Q1",
        context = "brain",
        method = "susie",
        entry = list(.cr_twoRegionEntry())
    )
    expect_equal(length(two), 2L)
    expect_equal(
        ctwasGroupPriors(two)[, c("group", "group_prior", "group_prior_var")],
        tibble::tibble(
            group = "gene",
            group_prior = 0.01,
            group_prior_var = NA_real_
        )
    )
})

test_that("ctwasGroupPriors() lists a joint run's shared priors per context", {
    prior <- list(group_prior = c(SNP = 1e-4, brain = 0.01, liver = 0.02))
    cr <- CtwasResult(
        gwasStudy = c("D1", "D1"),
        studyName = c("Q1", "Q1"),
        context = c("brain", "liver"),
        method = c("susie", "susie"),
        entry = list(
            .cr_entry(prior = prior),
            .cr_entry(ids = "g3", pip = 0.5, prior = prior)
        ),
        jointContexts = c("brain,liver", "brain,liver")
    )
    tab <- ctwasGroupPriors(cr)
    expect_equal(tab$context, rep(c("brain", "liver"), each = 3L))
    expect_equal(unique(tab$jointContexts), "brain,liver")
    expect_equal(tab$group_prior, rep(c(1e-4, 0.01, 0.02), 2L))
})

test_that("updateObject rebuilds a CtwasResult saved as a table of runs", {
    # The earlier shape: a DFrame of whole runs, one entry per row.
    old <- S4Vectors::DataFrame(
        gwasStudy = "D1",
        study = "Q1",
        context = "brain",
        method = "susie",
        entry = S4Vectors::SimpleList(list(.cr_twoRegionEntry()))
    )
    attr(old, "class") <- methods::getClass("CtwasResult")@className
    up <- updateObject(old)
    expect_s4_class(up, "CtwasResult")
    expect_true(methods::validObject(up))
    expect_equal(length(up), 2L)
    # A current object needs nothing.
    expect_identical(updateObject(up), up)
    expect_message(
        updateObject(up, verbose = TRUE),
        "updateObject\\(object = 'CtwasResult'\\)"
    )
})

test_that("CtwasResult() keeps the caller's metadata", {
    cr <- CtwasResult(
        gwasStudy = "D1",
        studyName = "Q1",
        context = "brain",
        method = "susie",
        entry = list(.cr_entry()),
        metadata = list(run = "r1")
    )
    expect_equal(S4Vectors::metadata(cr), list(run = "r1"))
})
