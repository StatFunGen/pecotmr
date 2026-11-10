# Tests for R/JointGroup.R — the uniform contract for the joint-analysis engine:
# the JointGroup hierarchy (conditions-table identity model), the
# JointDispatchCell wiring row, and the pipeline marker classes. Construction is
# validated, so a malformed group / mistyped cell fails loudly at the source.

.jg_cond <- function(studyName = "S", context = c("c1", "c2"), trait = "G") {
    data.frame(
        study = studyName,
        context = context,
        trait = trait,
        stringsAsFactors = FALSE
    )
}

test_that("JointGroup subclasses construct from a conditions table", {
    X <- matrix(
        0,
        10,
        2,
        dimnames = list(paste0("s", 1:10), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    Y <- matrix(0, 10, 2, dimnames = list(paste0("s", 1:10), c("c1", "c2")))
    g <- new("IndividualJointGroup", conditions = .jg_cond(), X = X, Y = Y)
    expect_s4_class(g, "JointGroup")
    expect_s4_class(g, "IndividualJointGroup")
    expect_equal(nrow(.jgConditions(g)), 2L)

    Z <- matrix(0, 3, 2, dimnames = list(paste0("v", 1:3), NULL))
    sg <- new(
        "SumStatsJointGroup",
        conditions = .jg_cond(),
        Z = Z,
        ldReference = NULL,
        N = c(100, 120)
    )
    expect_s4_class(sg, "JointGroup")
    expect_s4_class(sg, "SumStatsJointGroup")
})

test_that(".jgKeepConditions subsets each condition's column with it", {
    X <- matrix(0, 4, 1, dimnames = list(paste0("s", 1:4), "chr1:100:A:G"))
    Y <- matrix(
        as.numeric(1:12),
        4,
        3,
        dimnames = list(paste0("s", 1:4), c("c1", "c2", "c3"))
    )
    g <- new(
        "IndividualJointGroup",
        conditions = .jg_cond(context = c("c1", "c2", "c3")),
        X = X,
        Y = Y,
        traitPos = c(10, 20, 30)
    )
    kept <- .jgKeepConditions(g, c(TRUE, FALSE, TRUE))
    expect_equal(.jgConditions(kept)$context, c("c1", "c3"))
    expect_equal(.jgY(kept), Y[, c(1L, 3L)])
    expect_equal(.jgTraitPos(kept), c(10, 30))
    expect_identical(.jgKeepConditions(g, rep(TRUE, 3L)), g)
    expect_null(.jgKeepConditions(g, rep(FALSE, 3L)))
    sg <- new(
        "SumStatsJointGroup",
        conditions = .jg_cond(),
        Z = matrix(
            as.numeric(1:6),
            3,
            2,
            dimnames = list(paste0("v", 1:3), NULL)
        ),
        ldReference = NULL,
        N = c(100, 120)
    )
    keptSs <- .jgKeepConditions(sg, c(FALSE, TRUE))
    expect_equal(.jgZ(keptSs), .jgZ(sg)[, 2L, drop = FALSE])
    expect_equal(.jgN(keptSs), 120)
})

test_that("JointGroup validity rejects malformed groups", {
    X <- matrix(
        0,
        10,
        2,
        dimnames = list(paste0("s", 1:10), c("chr1:100:A:G", "chr1:200:A:G"))
    )
    Y <- matrix(0, 10, 2, dimnames = list(paste0("s", 1:10), c("c1", "c2")))
    # 1 condition is valid (the univariate cell); 0 conditions is not.
    expect_s4_class(
        new(
            "IndividualJointGroup",
            conditions = .jg_cond(context = "c1"),
            X = X[, 1, drop = FALSE],
            Y = Y[, 1, drop = FALSE]
        ),
        "IndividualJointGroup"
    )
    expect_error(
        new(
            "IndividualJointGroup",
            conditions = data.frame(
                study = character(0),
                context = character(0),
                trait = character(0)
            ),
            X = X,
            Y = Y[, 0, drop = FALSE]
        ),
        ">= 1 condition"
    )
    # Missing identity column.
    expect_error(
        new(
            "IndividualJointGroup",
            conditions = data.frame(study = "S", context = c("c1", "c2")),
            X = X,
            Y = Y
        ),
        "must have columns"
    )
    # X/Y row mismatch.
    expect_error(
        new(
            "IndividualJointGroup",
            conditions = .jg_cond(),
            X = X,
            Y = Y[1:5, , drop = FALSE]
        ),
        "dimension"
    )
    # ncol(Y) must equal nrow(conditions).
    expect_error(
        new(
            "IndividualJointGroup",
            conditions = .jg_cond(context = c("c1", "c2", "c3")),
            X = X,
            Y = Y
        ),
        "ncol\\(Y\\)"
    )
    # Z must name its variants: they are what the LD matrix is derived over,
    # so an unnamed Z would leave the group unable to say what it covers.
    expect_error(
        new(
            "SumStatsJointGroup",
            conditions = .jg_cond(),
            Z = matrix(0, 3, 2),
            ldReference = NULL,
            N = 1
        ),
        "variant ids as rownames"
    )
})

test_that("JointDispatchCell + pipeline markers validate at construction", {
    cell <- new(
        "JointDispatchCell",
        pattern = "context",
        dataForm = "individual",
        enumerate = function(data, scope, args) list(),
        minGroup = 2L
    )
    expect_s4_class(cell, "JointDispatchCell")
    expect_error(
        new(
            "JointDispatchCell",
            pattern = "context",
            dataForm = "bogus",
            enumerate = function() NULL,
            minGroup = 2L
        ),
        "dataForm"
    )
    expect_error(
        new(
            "JointDispatchCell",
            pattern = "context",
            dataForm = "individual",
            enumerate = function() NULL,
            minGroup = 0L
        ),
        "minGroup"
    )

    fm <- new("FmJointPipeline", config = list(coverage = 0.95))
    tw <- new("TwasJointPipeline", config = list(fitRetention = "slim"))
    expect_s4_class(fm, "JointPipeline")
    expect_s4_class(tw, "JointPipeline")
    expect_equal(fm@config$coverage, 0.95)
})

test_that("a summary-statistics group may hold one LD reference per study", {
    ids <- sprintf("chr1:%d:A:G", c(100L, 200L))
    ldm <- function(rho) {
        LdMatrix(`dimnames<-`(matrix(c(1, rho, rho, 1), 2L), list(ids, ids)))
    }
    cond <- .jg_cond(context = c("c1", "c1"))
    cond$study <- c("A", "B")
    refs <- list(ldm(0.2), ldm(0.8))
    Z <- matrix(c(1, 2, NA, 3), 2L, 2L, dimnames = list(ids, c("A", "B")))
    sg <- new(
        "SumStatsJointGroup",
        conditions = cond,
        Z = Z,
        ldReference = refs,
        N = c(100, 120)
    )
    expect_identical(pecotmr:::.jgLdRefAt(sg, 2L), refs[[2L]])
    # A fit taking one LD matrix has none to take.
    expect_error(
        pecotmr:::.jgOneLdRef(sg, "mrmash"),
        "mrmash: a joint fit uses one LD matrix"
    )
    kept <- .jgKeepConditions(sg, c(FALSE, TRUE))
    expect_identical(pecotmr:::.jgLdRef(kept), refs[2L])
    expect_error(
        new(
            "SumStatsJointGroup",
            conditions = cond,
            Z = Z,
            ldReference = refs[1L],
            N = c(100, 120)
        ),
        "one LD reference per condition"
    )
    one <- methods::initialize(sg, ldReference = refs[[1L]])
    expect_identical(pecotmr:::.jgOneLdRef(one, "mrmash"), refs[[1L]])
    expect_identical(pecotmr:::.jgLdRefAt(one, 2L), refs[[1L]])
})
