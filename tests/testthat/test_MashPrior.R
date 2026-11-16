# Tests for R/MashPrior.R — the data-driven (mash) prior bundle: full-data prior
# (fullFit) + per-fold cross-validated priors (cvFits) consumed by
# twasWeightsPipeline (mr.mash) and, transitively, fineMappingPipeline
# (mvSuSiE).

test_that("MashPrior: construct + accessors (full + cv)", {
    priors <- list(context = mashPriorPayload())
    sp <- data.frame(
        Sample = paste0("s", 1:6),
        Fold = rep(1:3, each = 2),
        stringsAsFactors = FALSE
    )
    perFold <- list(priors, priors, priors)
    mp <- MashPrior(
        fullFit = priors,
        cvFits = list(samplePartition = sp, perFoldFits = perFold)
    )
    expect_s4_class(mp, "MashPrior")
    expect_identical(fullFit(mp), priors)
    expect_length(cvFits(mp)$perFoldFits, 3L)
    expect_identical(cvFits(mp)$samplePartition, sp)
    expect_output(show(mp), "fullFit: context")
})

test_that("MashPrior: full-only and cv-only bundles", {
    priors <- list(context = mashPriorPayload())
    mpFull <- MashPrior(fullFit = priors)
    expect_null(cvFits(mpFull))
    mpCv <- MashPrior(cvFits = list(perFoldFits = list(priors, priors)))
    expect_null(fullFit(mpCv))
    expect_length(cvFits(mpCv)$perFoldFits, 2L)
})

test_that("MashPrior: each prior must say what its rows are", {
    prior <- mashPriorPayload()
    # The bare list(U, w) of a single prior is not a set of named priors.
    expect_error(
        MashPrior(fullFit = list(prior)),
        "named list with one prior per joint spec"
    )
    expect_error(
        MashPrior(fullFit = list(context = prior[c("U", "w")])),
        "must carry U, alpha, inputScale, axes, conditions"
    )
    expect_error(
        MashPrior(
            fullFit = list(context = list_assign(prior, conditions = "x"))
        ),
        "columns study, context, trait, n"
    )
    expect_error(
        MashPrior(
            fullFit = list(
                context = list_assign(prior, U = list(a = diag(3)))
            )
        ),
        "list of 2 x 2 matrices"
    )
    # A residual correlation is optional (earlier mashPipeline() output has
    # none), but when present it is one row and column per condition.
    expect_s4_class(
        MashPrior(
            fullFit = list(
                context = list_assign(prior, residualCorrelation = diag(2))
            )
        ),
        "MashPrior"
    )
    expect_error(
        MashPrior(
            fullFit = list(
                context = list_assign(prior, residualCorrelation = diag(3))
            )
        ),
        "residualCorrelation must be a 2 x 2 matrix"
    )
    expect_error(
        MashPrior(
            fullFit = list(context = list_assign(prior, axes = "tissue"))
        ),
        "axes must name study, context or trait"
    )
    expect_error(
        MashPrior(
            fullFit = list(context = list_assign(prior, inputScale = "auto"))
        ),
        "inputScale"
    )
    # The same rules hold for each fold's priors, which name the fold.
    expect_error(
        MashPrior(cvFits = list(perFoldFits = list(list(prior)))),
        "perFoldFits\\[\\[1\\]\\]"
    )
})

test_that("MashPrior: validity rejects malformed bundles", {
    sp <- data.frame(
        Sample = paste0("s", 1:6),
        Fold = rep(1:3, each = 2),
        stringsAsFactors = FALSE
    )
    # Empty bundle.
    expect_error(MashPrior(), "at least one")
    # perFoldFits count must match the partition's fold count.
    priors <- list(context = mashPriorPayload())
    expect_error(
        MashPrior(
            cvFits = list(samplePartition = sp, perFoldFits = list(priors))
        ),
        "fold"
    )
    # perFoldFits must be a non-empty list.
    expect_error(MashPrior(cvFits = list(perFoldFits = "nope")))
    # samplePartition must carry Sample + Fold columns.
    expect_error(
        MashPrior(
            cvFits = list(
                samplePartition = data.frame(x = 1),
                perFoldFits = list(priors)
            )
        ),
        "Sample"
    )
})

test_that("MashPrior: show reports 'cvFits: none' for a full-only bundle", {
    mp <- MashPrior(
        fullFit = list(
            context = mashPriorPayload(),
            trait = mashPriorPayload()
        )
    )
    out <- capture.output(show(mp))
    expect_true(any(grepl("cvFits: none", out)))
    expect_true(any(grepl("fullFit: context, trait", out)))
})


test_that("the cvFits check rejects a list with no perFoldFits", {
    f <- pecotmr:::.mashPriorCheckCvFits
    msg <- "`cvFits` must be a list with a `perFoldFits` element"
    expect_equal(f(list(somethingElse = 1)), msg)
    expect_equal(f("not a list"), msg)
    # NULL is the "no CV was run" case and is accepted.
    expect_null(f(NULL))
})


test_that("the partition check stops when perFoldFits is not a list", {
    # The fold-count comparison below it would be meaningless, and the
    # perFoldFits check already reports the type problem.
    expect_null(pecotmr:::.mashPriorCheckPartition(list(
        perFoldFits = "not a list",
        samplePartition = data.frame(Sample = "a", Fold = 1L)
    )))
})

test_that("MashPrior: cvFits must be a list or NULL", {
    priors <- list(context = mashPriorPayload())
    expect_error(
        MashPrior(fullFit = priors, cvFits = "nope"),
        "cvFits.*Must be of type 'list'"
    )
    expect_silent(MashPrior(fullFit = priors))
})
