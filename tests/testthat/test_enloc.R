test_that("the hypothesis probabilities are a partition", {
    set.seed(1)
    p <- 40L
    r <- enlocClusterProbabilities(
        gwasLog10Bf = rnorm(p, 1, 1),
        qtlPip = {
            v <- runif(p)
            v / sum(v) * 0.6
        },
        pi1e = 0.02,
        pi1ne = 0.001
    )
    expect_equal(r$ppH0 + r$ppH1 + r$ppH2 + r$ppH3 + r$ppH4, 1)
    expect_true(all(
        c(r$ppH0, r$ppH1, r$ppH2, r$ppH3, r$ppH4) >= 0
    ))
})

test_that("RCP and LCP are the coloc hypotheses under other names", {
    set.seed(2)
    p <- 25L
    r <- enlocClusterProbabilities(
        gwasLog10Bf = rnorm(p),
        qtlPip = {
            v <- runif(p)
            v / sum(v) * 0.4
        },
        pi1e = 0.05,
        pi1ne = 0.002
    )
    # fastenloc sums the grid's diagonal for RCP and the whole grid for LCP;
    # those are exactly "one shared causal variant" and "both traits causal".
    expect_equal(r$rcp, r$ppH4)
    expect_equal(r$lcp, r$ppH3 + r$ppH4)
    expect_equal(sum(r$scp), r$rcp)
})

test_that("enrichment raises the shared-variant probability", {
    set.seed(3)
    p <- 30L
    bf <- rnorm(p, 1, 1)
    pip <- runif(p)
    pip <- pip / sum(pip) * 0.5
    none <- enlocClusterProbabilities(bf, pip, pi1e = 0.001, pi1ne = 0.001)
    some <- enlocClusterProbabilities(bf, pip, pi1e = 0.05, pi1ne = 0.001)
    # a1 = 0 means pi1e == pi1ne: being an eQTL says nothing about GWAS
    # causality, which is the no-enrichment case.
    expect_gt(some$rcp, none$rcp)
})

test_that("the priors match fastenloc's set_enrich_params", {
    a0 <- -6
    a1 <- 3
    pQtl <- 0.01
    pr <- pecotmr:::.enlocPriorsFromEnrichment(a0, a1, pQtl)
    expect_equal(pr$pi1e, stats::plogis(a0 + a1))
    expect_equal(pr$pi1ne, stats::plogis(a0))
    # p1 = (1 - P_eqtl) * exp(a0)/(1 + exp(a0))
    expect_equal(pr$p1, (1 - pQtl) * exp(a0) / (1 + exp(a0)))
    # p2 = P_eqtl / (1 + exp(a0 + a1))
    expect_equal(pr$p2, pQtl / (1 + exp(a0 + a1)))
    # p12 = P_eqtl * exp(a0 + a1)/(1 + exp(a0 + a1))
    expect_equal(pr$p12, pQtl * exp(a0 + a1) / (1 + exp(a0 + a1)))
})

test_that("the enrichment parameters round-trip through the priors", {
    a0 <- -6
    a1 <- 3
    pr <- pecotmr:::.enlocPriorsFromEnrichment(a0, a1, pQtl = 0.01)
    back <- pecotmr:::.enlocEnrichmentFromPriors(pr$p1, pr$p2, pr$p12)
    expect_equal(back$a0, a0)
    expect_equal(back$a1, a1)
})

test_that("a1 is capped before the priors are formed", {
    uncapped <- pecotmr:::.enlocPriorsFromEnrichment(-6, 5, 0.01)
    capped <- pecotmr:::.enlocPriorsFromEnrichment(-6, 5, 0.01, capA1 = 2)
    expect_equal(capped$pi1e, stats::plogis(-6 + 2))
    expect_lt(capped$p12, uncapped$p12)
})

test_that("the SNP-level floor is applied only when asked", {
    bf <- c(0, 0, 0)
    pip <- c(1e-12, 1e-12, 1e-12)
    without <- enlocClusterProbabilities(
        bf,
        pip,
        pi1e = 0.02,
        pi1ne = 0.001,
        p12 = 1e-5
    )
    with <- enlocClusterProbabilities(
        bf,
        pip,
        pi1e = 0.02,
        pi1ne = 0.001,
        p12 = 1e-5,
        applyFloor = TRUE
    )
    # fastenloc floors each SNP-level term at p12 when the full variant set
    # is present, so a vanishing term cannot read as certainty against
    # colocalisation.
    expect_gt(with$rcp, without$rcp)
})

test_that("recovered GWAS Bayes factors follow the DAP-1 contrast", {
    pip <- c(0.6, 0.2, 0.0)
    bf <- pecotmr:::.enlocGwasLog10Bf(pip, pGwas = 1e-4, clusterPip = 0.8)
    # A zero PIP is floored rather than giving -Inf.
    expect_true(all(is.finite(bf)))
    # Monotone in the PIP, since the other terms are shared.
    expect_true(bf[[1]] > bf[[2]])
    expect_true(bf[[2]] > bf[[3]])
})

test_that("mismatched input lengths are rejected", {
    expect_error(
        enlocClusterProbabilities(c(1, 2), c(0.1), pi1e = 0.02, pi1ne = 0.001),
        "same variants"
    )
})
