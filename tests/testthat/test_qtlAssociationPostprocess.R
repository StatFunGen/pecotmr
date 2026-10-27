context("qtlAssociationPostprocess")

# The shared fixture without the named row columns, for the required-column /
# missing-input error paths: dropping a column is a cheap mcols edit, where a
# rebuild costs as much as the test itself.
.qapWithout <- function(cols) {
    x <- .qapDefaultFixture
    kept <- setdiff(colnames(S4Vectors::mcols(x)), cols)
    S4Vectors::mcols(x) <- S4Vectors::mcols(x)[, kept, drop = FALSE]
    x
}

test_that("qtlAssociationPostprocess enriches with package-computed columns", {
    x <- .qapDefaultFixture
    r <- qtlAssociationPostprocess(
        x,
        fdrThreshold = 0.1,
        mafCutoff = 0.01,
        cisWindow = 1e6
    )

    expect_s4_class(r, "QtlSumStats")
    expect_false(is.null(qcOptions(r)$associationPostprocess)) # recipe stashed

    # Bonferroni original == min over variants of p.adjust(P, "bonferroni", n).
    expP <- map_dbl(
        seq_len(nrow(x)),
        function(i) {
            min(stats::p.adjust(
                S4Vectors::mcols(x[[i]])$P,
                "bonferroni",
                n = 50
            ))
        }
    )
    expect_equal(as.numeric(r$p_bonferroni_min_original), expP)
    expect_equal(
        as.numeric(r$fdr_bonferroni_min_original),
        stats::p.adjust(expP, "fdr")
    )
    expect_true(all(
        r$q_bonferroni_min_original >= 0 &
            r$q_bonferroni_min_original <= 1
    ))

    # Filtered flavour present and uses the smaller test count (n_variants_filtered).
    expect_true(all(
        c(
            "p_bonferroni_min_filtered",
            "fdr_bonferroni_min_filtered",
            "q_bonferroni_min_filtered"
        ) %in%
            colnames(S4Vectors::mcols(r))
    ))
})

test_that("q-values come from qvalue::qvalue, not a hand-rolled fallback", {
    skip_if_not_installed("qvalue")
    x <- .qapDefaultFixture
    r <- qtlAssociationPostprocess(x, methods = "permutation")
    expect_equal(
        as.numeric(r$q_beta),
        qvalue::qvalue(as.numeric(x$p_beta))$qvalues
    )
    expect_equal(
        as.numeric(r$fdr_beta),
        stats::p.adjust(as.numeric(x$p_beta), "fdr")
    )
})

test_that("permutation nominal threshold == stats::qbeta of the empirical cutoff", {
    x <- .qapDefaultFixture
    r <- qtlAssociationPostprocess(
        x,
        fdrThreshold = 0.1,
        methods = "permutation"
    )
    qb <- as.numeric(r$q_beta)
    pb <- as.numeric(x$p_beta)
    lb <- sort(pb[qb <= 0.1])
    ub <- sort(pb[qb > 0.1])
    expect_gt(length(lb), 0)
    thr <- if (length(ub) > 0) (lb[length(lb)] + ub[1]) / 2 else lb[length(lb)]
    expect_equal(
        as.numeric(r$p_nominal_threshold),
        stats::qbeta(thr, rep(1, nrow(x)), rep(200, nrow(x)))
    )
})

test_that("significantQtls (bonferroni) matches the derived threshold rule", {
    x <- .qapDefaultFixture
    r <- qtlAssociationPostprocess(x, mafCutoff = 0.01, cisWindow = 1e6)
    sig <- significantQtls(r, "bonferroni_original", threshold = 0.5)
    expect_s4_class(sig, "GRanges")
    expect_true("trait" %in% names(S4Vectors::mcols(sig)))

    # Reproduce the rule: variants with P*n_variants <= max(p_bonf_min over FDR-sig genes).
    fdr <- as.numeric(r$fdr_bonferroni_min_original)
    sigGenes <- which(fdr < 0.5)
    expect_gt(length(sigGenes), 0)
    varThr <- max(as.numeric(r$p_bonferroni_min_original)[sigGenes])
    expN <- sum(map_int(
        seq_len(nrow(r)),
        function(i) {
            sum(pmin(1, S4Vectors::mcols(r[[i]])$P * 50) <= varThr)
        }
    ))
    expect_equal(length(sig), expN)
})

test_that("significance accessors require a postprocessed object", {
    x <- .qapDefaultFixture # not postprocessed
    expect_error(
        significantQtls(x, "bonferroni_original"),
        "qtlAssociationPostprocess"
    )
})

# ---- qvalue-native edge-case retries (the ONLY substitute is qvalue itself) --
# The fixture carries p_beta but no q_beta, so permutation postprocessing routes
# q_beta through .qapSafeQvalue -> qvalue::qvalue. Mock qvalue to raise each
# qvalue-native error on the primary (no-lambda) call and succeed on the retry
# (lambda supplied), then assert the retry value flows through unchanged.

test_that(".qapSafeQvalue retries both qvalue-native degenerate cases", {
    skip_if_not_installed("qvalue")
    # 'missing or infinite' retries with lambda = 0, 'pi0 <= 0' with a
    # bootstrap pi0.
    for (err in c(
        "missing or infinite values in the smoother",
        "ERROR: The estimated pi0 <= 0."
    )) {
        testthat::local_mocked_bindings(
            qvalue = function(p, lambda, ...) {
                if (missing(lambda)) {
                    stop(err)
                }
                list(qvalues = rep(0.111, length(p)))
            },
            .package = "qvalue"
        )
        r <- qtlAssociationPostprocess(
            .qapDefaultFixture,
            methods = "permutation"
        )
        expect_true(all(as.numeric(r$q_beta) == 0.111), info = err)
    }
})

test_that(".qapSafeQvalue re-raises a non-native qvalue error (no hand-rolled q)", {
    skip_if_not_installed("qvalue")
    testthat::local_mocked_bindings(
        qvalue = function(p, ...) stop("some unrelated numerical failure"),
        .package = "qvalue"
    )
    expect_error(
        qtlAssociationPostprocess(.qapDefaultFixture, methods = "permutation"),
        "qvalue::qvalue failed"
    )
})

test_that("permutation nominal threshold is NA for every gene when none pass FDR", {
    # An impossibly small FDR threshold leaves no q_beta-significant gene, so the
    # empirical p_beta bracketing is empty and every nominal threshold is NA.
    r <- qtlAssociationPostprocess(
        .qapDefaultFixture,
        fdrThreshold = 1e-300,
        methods = "permutation"
    )
    expect_true(all(is.na(as.numeric(r$p_nominal_threshold))))
})

# ---- significantQtls: one branch of the significance mask per method -------

test_that("significantQtls (permutation) uses each gene's nominal threshold", {
    r <- qtlAssociationPostprocess(
        .qapDefaultFixture,
        fdrThreshold = 0.1,
        methods = "permutation"
    )
    sig <- significantQtls(r, "permutation") # threshold defaults to recipe
    expect_s4_class(sig, "GRanges")
    expect_true("trait" %in% names(S4Vectors::mcols(sig)))
    # Reproduce: per gene, variants with P < p_nominal_threshold[gene].
    thr <- as.numeric(r$p_nominal_threshold)
    expN <- sum(map_int(
        seq_len(nrow(r)),
        function(i) {
            if (is.na(thr[i])) {
                return(0L)
            }
            sum(S4Vectors::mcols(r[[i]])$P < thr[i])
        }
    ))
    expect_gt(expN, 0)
    expect_equal(length(sig), expN)
})

test_that("significantQtls (permutation) skips genes with an NA threshold", {
    # Gene 1's NA beta_shape1 makes its qbeta nominal threshold NA, so the mask
    # loop must `next` past it while still processing the non-NA genes.
    x <- .qapDefaultFixture
    S4Vectors::mcols(x)$beta_shape1[1L] <- NA_real_
    r <- qtlAssociationPostprocess(
        x,
        fdrThreshold = 0.1,
        methods = "permutation"
    )
    thr <- as.numeric(r$p_nominal_threshold)
    expect_true(is.na(thr[1])) # NA-threshold gene present ...
    expect_true(any(!is.na(thr))) # ... alongside non-NA genes
    sig <- significantQtls(r, "permutation")
    expect_s4_class(sig, "GRanges")
    # Gene 1 contributes nothing; total is the sum over the non-NA-threshold genes.
    expN <- sum(map_int(
        seq_len(nrow(r)),
        function(i) {
            if (is.na(thr[i])) {
                return(0L)
            }
            sum(S4Vectors::mcols(r[[i]])$P < thr[i])
        }
    ))
    expect_equal(length(sig), expN)
})

test_that("significantQtls (permutation) errors without a nominal threshold", {
    # No beta shapes -> postprocess produces no p_nominal_threshold column.
    r <- qtlAssociationPostprocess(
        .qapWithout(c("beta_shape1", "beta_shape2")),
        methods = "permutation"
    )
    expect_null(r$p_nominal_threshold)
    expect_error(significantQtls(r, "permutation"), "p_nominal_threshold")
})

test_that("significantQtls (bonferroni) errors when its columns are absent", {
    r <- qtlAssociationPostprocess(.qapDefaultFixture, methods = "permutation") # no bonferroni
    expect_error(significantQtls(r, "bonferroni_original"), "columns absent")
})

test_that("significantQtls (qvalue) selects variants by their qvalue mcol", {
    skip_if_not_installed("qvalue")
    r <- qtlAssociationPostprocess(
        .qapDefaultFixture,
        fdrThreshold = 0.1,
        methods = "permutation"
    )
    sig <- significantQtls(r, "qvalue", threshold = 0.1)
    expect_s4_class(sig, "GRanges")
    # Reproduce: within genes whose event q_beta < 0.1, variants with qvalue < 0.1.
    qb <- as.numeric(r$q_beta)
    sigGenes <- which(qb < 0.1)
    expect_gt(length(sigGenes), 0)
    expN <- sum(map_int(
        sigGenes,
        function(i) {
            sum(S4Vectors::mcols(r[[i]])$qvalue < 0.1)
        }
    ))
    expect_gt(expN, 0)
    expect_equal(length(sig), expN)
})

test_that("significantQtls (qvalue) errors when no event q column exists", {
    # Drop p_beta so permutation adds no q_beta, and skip bonferroni so there is
    # no q_bonferroni_min_original either -- yet the object is still postprocessed.
    r <- qtlAssociationPostprocess(
        .qapWithout("p_beta"),
        methods = "permutation"
    )
    expect_null(r$q_beta)
    expect_null(r$q_bonferroni_min_original)
    expect_error(significantQtls(r, "qvalue"), "event q-value column")
})

test_that("significantQtls returns an empty GRanges when nothing is significant", {
    r <- qtlAssociationPostprocess(
        .qapDefaultFixture,
        mafCutoff = 0.01,
        cisWindow = 1e6
    )
    sig <- significantQtls(r, "bonferroni_original", threshold = 1e-300)
    expect_s4_class(sig, "GRanges")
    expect_length(sig, 0)
})

# ---- qtlAssociationPostprocess: required-column + empty-entry guards ----------

test_that("qtlAssociationPostprocess requires n_variants for Bonferroni", {
    expect_error(
        qtlAssociationPostprocess(
            .qapWithout("n_variants"),
            methods = "bonferroni"
        ),
        "required for the Bonferroni"
    )
})

test_that("qtlAssociationPostprocess requires n_variants_filtered when filtering", {
    expect_error(
        qtlAssociationPostprocess(
            .qapWithout("n_variants_filtered"),
            methods = "bonferroni",
            mafCutoff = 0.01
        ),
        "n_variants_filtered"
    )
})

test_that("qtlAssociationPostprocess skips genes whose entry has no variants", {
    r <- qtlAssociationPostprocess(
        .qapFixture(emptyIdx = 1L),
        methods = "bonferroni"
    )
    pmin1 <- as.numeric(r$p_bonferroni_min_original)
    expect_true(is.na(pmin1[1])) # empty gene -> NA Bonferroni min
    expect_true(all(!is.na(pmin1[-1]))) # every other gene is finite
})

test_that("FILTERED Bonferroni drops MAF/cis-failing variants + uses n_variants_filtered", {
    # Fixture entries: af = (.30, .20, .005, .40), tss/tes_distance = (0, 500, 9e5, 2e6).
    # With mafCutoff 0.01 + cisWindow 1e6 the filtered set keeps only v1, v2
    # (v3 fails MAF, v4 is outside the cis window); the filtered count is 30.
    x <- .qapDefaultFixture
    r <- qtlAssociationPostprocess(
        x,
        mafCutoff = 0.01,
        cisWindow = 1e6,
        methods = "bonferroni"
    )
    P <- map(seq_len(nrow(x)), function(i) S4Vectors::mcols(x[[i]])$P)
    expFilt <- map_dbl(P, function(p) min(pmin(1, p[1:2] * 30))) # v1,v2 @ n=30
    expOrig <- map_dbl(P, function(p) min(pmin(1, p * 50))) # all @ n=50
    expect_equal(as.numeric(r$p_bonferroni_min_filtered), expFilt)
    expect_equal(as.numeric(r$p_bonferroni_min_original), expOrig)
    # A smaller test count makes the filtered flavour no less significant.
    expect_true(all(
        r$p_bonferroni_min_filtered <= r$p_bonferroni_min_original + 1e-12
    ))
})

test_that("significantQtls(bonferroni_filtered) applies the derived rule on the filtered set", {
    x <- .qapDefaultFixture
    r <- qtlAssociationPostprocess(
        x,
        mafCutoff = 0.01,
        cisWindow = 1e6,
        methods = "bonferroni"
    )
    sig <- significantQtls(r, "bonferroni_filtered", threshold = 0.5)
    expect_s4_class(sig, "GRanges")
    fdr <- as.numeric(r$fdr_bonferroni_min_filtered)
    sigGenes <- which(!is.na(fdr) & fdr < 0.5)
    expect_gt(length(sigGenes), 0) # signal genes pass
    varThr <- max(as.numeric(r$p_bonferroni_min_filtered)[sigGenes])
    # only MAF/cis-passing variants (v1,v2) with P*n_filtered <= threshold qualify
    expN <- sum(map_int(
        seq_len(nrow(r)),
        function(i) {
            p <- S4Vectors::mcols(r[[i]])$P[1:2]
            sum(pmin(1, p * 30) <= varThr)
        }
    ))
    expect_gt(expN, 0)
    expect_equal(length(sig), expN)
    # Every returned variant satisfies the filtered keep rule (MAF and the
    # cis window).
    af <- S4Vectors::mcols(sig)$af
    expect_true(all(pmin(af, 1 - af) > 0.01))
    expect_true(all(
        S4Vectors::mcols(sig)$tss_distance >= -1e6 &
            S4Vectors::mcols(sig)$tes_distance <= 1e6
    ))
})

test_that(".qapSignificanceMask returns empty masks for an unknown method", {
    # Both public callers wrap this in arg_match(), so the fall-through is
    # only reachable by calling the helper directly. It is tested here so the
    # contract is explicit: an unrecognized method yields the untouched empty
    # masks rather than erroring or silently returning the wrong ones.
    g <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(100, 200), width = 1)
    )
    S4Vectors::mcols(g)$variant_id <- c("chr1:100:A:G", "chr1:200:C:T")
    S4Vectors::mcols(g)$Z <- c(1, 2)
    S4Vectors::mcols(g)$N <- c(10L, 10L)
    qss <- QtlSumStats(
        studyName = "s1",
        context = "c1",
        trait = "g1",
        entry = list(g),
        genome = "hg19"
    )
    local_mocked_bindings(
        qcOptions = function(x) {
            list(
                associationPostprocess = list(
                    fdrThreshold = 0.05,
                    pvalueCol = "pval_nominal"
                )
            )
        },
        .qapEmptyMask = function(i, x, pcol) str_c("EMPTY", i),
        .package = "pecotmr"
    )
    expect_equal(
        pecotmr:::.qapSignificanceMask(qss, "not_a_method", 0.05),
        list("EMPTY1")
    )
})

test_that(".qapPermutationCols reuses a supplied q_beta", {
    x <- list(
        p_beta = c(0.01, 0.2),
        q_beta = c(0.02, 0.3),
        beta_shape1 = c(1, 1),
        beta_shape2 = c(1, 1),
        pval_nominal = c(0.01, 0.2)
    )
    cols <- pecotmr:::.qapPermutationCols(x, 0.05)
    # q_beta is already present, so it is not recomputed or re-added.
    expect_setequal(names(cols), c("fdr_beta", "p_nominal_threshold"))
})

test_that("QvalueOptions refuses the p-value vector", {
    expect_error(QvalueOptions(p = 0.1), "the p-value vector")
    # qvalue::qvalue() takes dots, so nothing else can be rejected.
    expect_s4_class(QvalueOptions(nosuch = 1), "MethodOptions")
    expect_equal(
        QvalueOptions(pi0.method = "bootstrap")$pi0.method,
        "bootstrap"
    )
})

test_that(".qapSafeQvalue forwards qvalueArgs to qvalue::qvalue", {
    skip_if_not_installed("qvalue")
    seen <- NULL
    real <- qvalue::qvalue
    set.seed(4)
    p <- c(runif(200, 0, 0.05), runif(800))
    with_mocked_bindings(
        pecotmr:::.qapSafeQvalue(p, QvalueOptions(pfdr = TRUE)),
        qvalue = function(...) {
            seen <<- list(...)
            real(...)
        },
        .package = "qvalue"
    )
    expect_true(seen$pfdr)
})

test_that("the qvalue degenerate-case retries keep the caller's options", {
    skip_if_not_installed("qvalue")
    seen <- list()
    real <- qvalue::qvalue
    set.seed(5)
    p <- runif(200)
    res <- with_mocked_bindings(
        pecotmr:::.qapSafeQvalue(p, QvalueOptions(pfdr = TRUE)),
        qvalue = function(...) {
            args <- list(...)
            seen[[length(seen) + 1L]] <<- args
            # Force the first attempt into qvalue's own "pi0 <= 0" case.
            if (is.null(args$pi0.method)) {
                stop("ERROR: The estimated pi0 <= 0.")
            }
            real(...)
        },
        .package = "qvalue"
    )
    expect_length(seen, 2L)
    # The retry sets lambda/pi0.method but must not drop what the caller set.
    expect_equal(seen[[2]]$pi0.method, "bootstrap")
    expect_true(seen[[2]]$pfdr)
    expect_length(res, length(p))
})
