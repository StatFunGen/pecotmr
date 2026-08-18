context("credibleSetSummary")

test_that("getCredibleSetSummary returns one row per CS with size/purity/V/logBF/lead", {
    vn <- paste0("chr1:", (1:5) * 100, ":A:G")
    L <- 2L
    P <- 5L
    alpha <- matrix(0.1, L, P)
    alpha[1, 1:2] <- c(0.5, 0.4)
    alpha[2, 3:4] <- c(0.5, 0.4)
    lbf <- matrix(
        c(3.0, 2.5, 0.1, 0.1, 0.0, 0.1, 0.1, 2.0, 1.8, 0.0),
        L,
        P,
        byrow = TRUE
    )
    fit <- list(
        alpha = alpha,
        mu = matrix(0.3, L, P),
        mu2 = matrix(1.2, L, P),
        pip = c(0.6, 0.5, 0.4, 0.4, 0.02),
        V = c(0.8, 0.5),
        lbf_variable = lbf,
        sets = list(
            purity = data.frame(
                min.abs.corr = c(0.9, 0.3),
                mean.abs.corr = c(0.95, 0.4),
                row.names = c("L1", "L2")
            )
        )
    )
    class(fit) <- "susie"
    cst <- list(
        list(
            sets = list(
                cs = list(L1 = c(1L, 2L), L2 = c(3L, 4L)),
                purity = data.frame(min.abs.corr = c(0.9, 0.3))
            )
        ),
        list(sets = list(cs = list())),
        list(sets = list(cs = list()))
    )
    attr(cst, "coverage") <- c(0.95, 0.70, 0.50)
    tl <- buildTopLoci(fit, cst, variantNames = vn, method = "susie")
    e <- FineMappingEntry(variantIds = vn, susieFit = fit, topLoci = tl)

    s <- getCredibleSetSummary(e)
    s <- s[order(s$cs), , drop = FALSE]
    expect_equal(nrow(s), 2L)
    expect_equal(s$cs, c("susie_1", "susie_2"))
    expect_equal(s$effect_id, c("L1", "L2"))
    expect_equal(s$n_variants, c(2L, 2L))
    expect_equal(s$purity_min, c(0.9, 0.3))
    expect_equal(s$purity_mean, c(0.95, 0.4))
    expect_equal(s$V, c(0.8, 0.5))
    # cs_log10bf = max member logBF (per-variant max single-effect lbf)
    expect_equal(s$cs_log10bf, c(3.0, 2.0))
    # cs_log_bf = true per-effect single-effect logBF = logSumExp(lbf[L, ]) - log(P)
    expect_equal(s$cs_log_bf, c(1.9597, 1.2031), tolerance = 1e-3)
    # cs_pip = summed member PIP (inclusion mass captured)
    expect_equal(s$cs_pip, c(1.1, 0.8), tolerance = 1e-6)
    # cs_mean_effect = NA for univariate susie (no conditional_effect column)
    expect_true(all(is.na(s$cs_mean_effect)))
    # lead = highest-PIP member
    expect_equal(s$lead_variant, c("chr1:100:A:G", "chr1:300:A:G"))
})

test_that("getCredibleSetSummary is empty when there are no credible sets", {
    vn <- c("chr1:100:A:G", "chr1:200:C:T")
    fit <- list(
        alpha = matrix(c(0.2, 0.1), 1, 2),
        mu = matrix(0.1, 1, 2),
        mu2 = matrix(1, 1, 2),
        pip = c(0.05, 0.04)
    )
    class(fit) <- "susie"
    cst <- list(
        list(sets = list(cs = list())),
        list(sets = list(cs = list())),
        list(sets = list(cs = list()))
    )
    attr(cst, "coverage") <- c(0.95, 0.70, 0.50)
    tl <- buildTopLoci(fit, cst, variantNames = vn, method = "susie")
    e <- FineMappingEntry(variantIds = vn, susieFit = fit, topLoci = tl)
    expect_equal(nrow(getCredibleSetSummary(e)), 0L)
})

test_that("getCredibleSetSummary is empty when the requested coverage column is absent", {
    # topLoci carries only cs_95 / cs_70 / cs_50; asking for coverage 0.99 finds no
    # cs_99 column -> the .csSummaryFit guard returns the empty summary.
    vn <- c("chr1:100:A:G", "chr1:200:C:T")
    tl <- data.frame(
        variant_id = vn,
        pip = c(0.7, 0.2),
        logBF = c(2.5, 0.3),
        cs_95 = c("susie_1", "susie_1"),
        cs_95_purity = c(0.9, 0.9),
        stringsAsFactors = FALSE
    )
    e <- FineMappingEntry(variantIds = vn, susieFit = list(x = 1), topLoci = tl)
    expect_equal(nrow(getCredibleSetSummary(e, coverage = 0.99)), 0L)
})

test_that("getCredibleSetSummary reports NA cs_log_bf when the fit carries no lbf matrix", {
    # A trimmed fit drops lbf_variable: the per-effect single-effect logBF degrades
    # to NA while the max-member cs_log10bf still comes from the topLoci.
    vn <- c("chr1:100:A:G", "chr1:200:C:T")
    tl <- data.frame(
        variant_id = vn,
        pip = c(0.7, 0.2),
        logBF = c(2.5, 0.3),
        cs_95 = c("susie_1", "susie_1"),
        cs_95_purity = c(0.88, 0.88),
        stringsAsFactors = FALSE
    )
    e <- FineMappingEntry(
        variantIds = vn,
        susieFit = list(trimmed = TRUE),
        topLoci = tl
    )
    s <- getCredibleSetSummary(e)
    expect_equal(nrow(s), 1L)
    expect_true(is.na(s$cs_log_bf))
    expect_equal(s$cs_log10bf, 2.5)
})

test_that("getCredibleSetSummary reports NA cs_log_bf when the effect's lbf row is all non-finite", {
    # effect 1's lbf row is entirely -Inf, so the finite subset is empty and the
    # logSumExp path short-circuits to NA.
    vn <- c("chr1:100:A:G", "chr1:200:C:T")
    lbf <- matrix(c(-Inf, -Inf, 2.0, 0.1), nrow = 2, byrow = TRUE)
    tl <- data.frame(
        variant_id = vn,
        pip = c(0.7, 0.2),
        logBF = c(1.0, 0.3),
        cs_95 = c("susie_1", "susie_1"),
        cs_95_purity = c(0.9, 0.9),
        stringsAsFactors = FALSE
    )
    fit <- structure(list(lbf_variable = lbf), class = "susie")
    e <- FineMappingEntry(variantIds = vn, susieFit = fit, topLoci = tl)
    expect_true(is.na(getCredibleSetSummary(e)$cs_log_bf))
})

test_that("getCredibleSetSummary aggregates across a collection with entry identity", {
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
        FineMappingEntry(variantIds = vn, susieFit = list(x = 1), topLoci = tl)
    }
    res <- QtlFineMappingResult(
        study = c("s", "s"),
        context = c("brain", "blood"),
        trait = c("g", "g"),
        method = c("susie", "susie"),
        entry = list(mkEntry(), mkEntry())
    )
    s <- getCredibleSetSummary(res)
    expect_true(all(
        c("study", "context", "trait", "method", "cs") %in% names(s)
    ))
    expect_equal(nrow(s), 2L) # one CS per entry
    expect_setequal(unique(s$context), c("brain", "blood"))
})
