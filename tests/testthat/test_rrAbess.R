context("regularized_regression - abess")

# The least-squares fit of each condition's phenotype on `support`, on the
# samples observing it, as abess refits its chosen subset.
.abessRefit <- function(X, Y, support) {
    sapply(seq_len(ncol(Y)), function(k) {
        rows <- !is.na(Y[, k])
        coef(lm(Y[rows, k] ~ X[rows, support, drop = FALSE]))[-1L]
    })
}

# ---- l0learnWeights on abess ----
test_that("l0learnWeights on abess fits a best subset of one phenotype", {
    skip_if_not_installed("abess")
    d <- .rrwXy(n = 100, p = 10)
    set.seed(1)
    w <- l0learnWeights(d$X, d$y, backend = "abess", fitRetention = "slim")
    expect_equal(dim(w), c(10L, 1L))
    support <- which(w[, 1L] != 0)
    expect_equal(support, 1:2)
    expect_equal(
        w[support, 1L],
        .abessRefit(d$X, cbind(d$y), support)[, 1L],
        tolerance = 1e-6,
        check.attributes = FALSE
    )
    expect_named(
        attr(w, "fit"),
        c("coef", "intercept", "supportSize", "cvMean", "variantIds")
    )
    expect_equal(attr(w, "fit")$supportSize, 2L)
})

test_that("abess fits conditions on the same samples with mgaussian", {
    skip_if_not_installed("abess")
    d <- .rrwMulti(n = 80, p = 8, K = 3)
    set.seed(1)
    # abess's mgaussian family stops on an unnamed phenotype matrix.
    w <- l0learnWeights(d$X, unname(d$Y), backend = "abess")
    expect_equal(dim(w), c(8L, 3L))
    support <- which(rowSums(w != 0) > 0)
    expect_true(all(rowSums(w != 0) %in% c(0, 3)))
    expect_equal(
        unname(w[support, , drop = FALSE]),
        unname(.abessRefit(d$X, d$Y, support)),
        tolerance = 1e-6
    )
})

test_that("abess fits conditions on samples of their own by groups", {
    skip_if_not_installed("abess")
    d <- .rrwMulti(n = 120, p = 8, K = 2)
    Y <- d$Y
    Y[1:60, 1L] <- NA
    Y[61:120, 2L] <- NA
    set.seed(1)
    w <- l0learnWeights(d$X, Y, backend = "abess", fitRetention = "slim")
    support <- which(rowSums(w != 0) > 0)
    expect_true(length(support) > 0L)
    expect_true(all(rowSums(w != 0) %in% c(0, 2)))
    expect_equal(
        unname(w[support, , drop = FALSE]),
        unname(.abessRefit(d$X, Y, support)),
        tolerance = 1e-6
    )
    expect_equal(attr(w, "fit")$supportSize, length(support))
})

test_that("abess leaves R's random numbers as its folds left them", {
    skip_if_not_installed("abess")
    d <- .rrwXy(n = 60, p = 6)
    set.seed(9)
    invisible(sample(rep_len(seq_len(5L), 60L)))
    expected <- runif(1L)
    set.seed(9)
    invisible(abessWeights(d$X, d$y))
    expect_equal(runif(1L), expected)
})

test_that("abess keeps its whole fit without the data its call holds", {
    skip_if_not_installed("abess")
    d <- .rrwXy(n = 60, p = 6)
    set.seed(1)
    full <- attr(
        l0learnWeights(d$X, d$y, backend = "abess", fitRetention = "full"),
        "fit"
    )
    expect_s3_class(full, "abess")
    expect_identical(full$call$x, quote(x))
    expect_identical(full$call$y, quote(y))
})

test_that("abess fits the L0 penalty, and L0Learn one phenotype", {
    skip_if_not_installed("abess")
    skip_if_not_installed("L0Learn")
    d <- .rrwMulti(n = 40, p = 5, K = 2)
    expect_error(
        l0learnWeights(d$X, d$Y[, 1L], penalty = "L0L2", backend = "abess"),
        "L0L1 and L0L2 take backend = \"l0learn\""
    )
    expect_error(l0learnWeights(d$X, d$Y), "backend = \"abess\"")
})

# ---- AbessOptions ----
test_that("AbessOptions refuses what pecotmr sets and what abess mishandles", {
    skip_if_not_installed("abess")
    expect_error(AbessOptions(num.threads = 2), "one thread")
    expect_error(AbessOptions(weight = rep(1, 3)), "square roots")
    expect_error(AbessOptions(init.active.set = 1), "an index of its own")
    expect_error(AbessOptions(foldid = 1:3), "deals each sample")
    expect_error(AbessOptions(tune.type = "gic"), "cross-validation")
    expect_error(AbessOptions(supportsize = 2), "unknown argument")
    expect_equal(AbessOptions(support.size = 0:3)$support.size, 0:3)
    expect_equal(pecotmr:::.optionsInputType(AbessOptions()), "QtlDataset")
})
