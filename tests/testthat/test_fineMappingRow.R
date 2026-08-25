# Tests for FineMappingRow (S4 class)

# Helper for adjustPips tests, migrated from test_dataStructures.R

.makeAdjustEntry <- function(vids, L = 2L) {
    p <- length(vids)
    set.seed(11L)
    lbf <- matrix(rnorm(L * p), nrow = L, ncol = p)
    colnames(lbf) <- vids
    alpha <- lbfToAlpha(lbf)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    fineMappingRow(
        variantIds = vids,
        susieFit = list(
            pip = pip,
            alpha = alpha,
            lbf_variable = lbf,
            mu = matrix(0, L, p),
            X_column_scale_factors = rep(1, p)
        ),
        topLoci = data.frame(
            variant_id = vids,
            pip = pip,
            betahat = rep(0, p),
            sebetahat = rep(1, p),
            stringsAsFactors = FALSE
        )
    )
}


# ===========================================================================
# Tests migrated from test_dataStructures.R (getTopLoci, adjustPips)
# ===========================================================================

test_that(".fmrRowTopLoci(type='GRanges') converts topLoci data.frame to GRanges", {
    tl <- data.frame(
        variant_id = c("1:100:A:G", "1:200:C:T"),
        pip = c(0.9, 0.1),
        betahat = c(0.5, -0.2),
        sebetahat = c(0.1, 0.2),
        cs = c(1L, 0L),
        method = "susie",
        stringsAsFactors = FALSE
    )
    ent <- fineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl
    )
    gr <- .fmrRowTopLoci(ent, type = "GRanges")
    expect_s4_class(gr, "GRanges")
    expect_equal(length(gr), 2)
    expect_equal(S4Vectors::mcols(gr)$pip, c(0.9, 0.1))
})


test_that(".fmrRowTopLoci(type='GRanges') handles empty input", {
    ent <- fineMappingRow(
        variantIds = character(0),
        susieFit = list(),
        topLoci = data.frame()
    )
    gr <- .fmrRowTopLoci(ent, type = "GRanges")
    expect_s4_class(gr, "GRanges")
    expect_equal(length(gr), 0)
})


test_that("getTopLoci defaults to data.frame", {
    tl <- data.frame(
        variant_id = "1:100:A:G",
        pip = 0.9,
        betahat = 0.5,
        sebetahat = 0.1,
        cs = 1L,
        stringsAsFactors = FALSE
    )
    ent <- fineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl
    )
    expect_s3_class(.fmrRowTopLoci(ent), "data.frame")
})

# =============================================================================
# extractBlockGenotypes returns RSE
# =============================================================================

test_that("adjustPips renormalizes PIPs on a kept FineMappingRow subset", {
    vids <- paste0("chr1:", 1:6, ":A:G")
    entry <- .makeAdjustEntry(vids)
    keep <- vids[2:5]
    adj <- adjustPips(entry, keep)
    expect_s4_class(adj, "FineMappingRow")
    expect_equal(getVariantIds(adj), keep)
    expect_equal(ncol(getSusieFit(adj)$lbf_variable), 4)
    # Renormalized: each effect's alpha row sums to 1 (when row has any signal)
    expect_true(all(abs(rowSums(getSusieFit(adj)$alpha) - 1) < 1e-10))
    # PIPs match topLoci
    expect_equal(.fmrPartsTopLoci(adj)$pip, getSusieFit(adj)$pip)
    # PIPs change under renormalization
    origPips <- .fmrRowPip(entry)
    expect_false(identical(unname(origPips[keep]), getSusieFit(adj)$pip))
})


test_that("adjustPips errors when the intersection is empty", {
    vids <- paste0("chr1:", 1:4, ":A:G")
    entry <- .makeAdjustEntry(vids)
    expect_error(
        adjustPips(entry, paste0("chr2:", 1:4, ":A:G")),
        "intersection.*empty"
    )
})


test_that("adjustPips tolerates a chr-prefix difference between entry and keepVariants", {
    vids <- paste0("chr1:", 1:6, ":A:G")
    entry <- .makeAdjustEntry(vids)
    keep <- paste0("1:", 2:5, ":A:G") # same variants, no "chr" prefix
    adj <- adjustPips(entry, keep)
    expect_s4_class(adj, "FineMappingRow")
    expect_equal(getVariantIds(adj), vids[2:5]) # entry keeps its own labels
    expect_equal(ncol(getSusieFit(adj)$lbf_variable), 4)
})


test_that("adjustPips on a FineMappingResultBase collection renormalizes each entry", {
    vidsA <- paste0("chr1:", 1:6, ":A:G")
    vidsB <- paste0("chr1:", 3:8, ":A:G")
    entryA <- .makeAdjustEntry(vidsA)
    entryB <- .makeAdjustEntry(vidsB)
    fmr <- QtlFineMappingResult(
        study = c("s1", "s1"),
        context = c("c1", "c2"),
        trait = c("g1", "g1"),
        method = c("susie", "susie"),
        entry = list(entryA, entryB)
    )
    # Keep only variants shared by both entries' raw sets.
    keep <- intersect(vidsA, vidsB)
    adj <- adjustPips(fmr, keep)
    expect_s4_class(adj, "QtlFineMappingResult")
    expect_equal(nrow(adj), 2L)
    expect_equal(getVariantIds(pecotmr:::.collectionEntry(adj, 1L)), keep)
    expect_equal(getVariantIds(pecotmr:::.collectionEntry(adj, 2L)), keep)
})


# === Tests migrated from test_s4Constructors.R (FineMappingRow) ===

test_that("FineMappingRow: constructor stores slots and accessors return them", {
    tl <- .sc_makeTopLoci(3)
    tl$variant_id <- sprintf("chr9:%d:A:G", 100L * (1:3))
    entry <- fineMappingRow(
        variantIds = sprintf("chr9:%d:A:G", 100L * (1:3)),
        susieFit = list(payload = 1L),
        topLoci = tl
    )
    expect_s4_class(entry, "FineMappingRow")
    expect_equal(
        .fmrPartsVariantIds(entry),
        sprintf("chr9:%d:A:G", 100L * (1:3))
    )
    expect_equal(.fmrPartsSusieFit(entry), list(payload = 1L))
    # getTopLoci returns the projected posterior view, not the raw slot
    out <- .fmrRowTopLoci(entry, signalCutoff = 0)
    expect_equal(out$variant_id, sprintf("chr9:%d:A:G", 100L * (1:3)))
})


test_that("FineMappingRow: getPip returns named pip vector keyed by variant_id", {
    entry <- .sc_makeFineMappingRow(3)
    pip <- .fmrRowPip(entry)
    expect_equal(length(pip), 3L)
    expect_equal(names(pip), paste0("chr1:", 100 * 1:3, ":A:G"))
})

test_that("FineMappingRow: resolveWeights returns topLoci posterior effect aligned to variant_id", {
    entry <- .sc_makeFineMappingRow(3) # topLoci posterior_mean = 0.05 for all
    wr <- .fmrRowResolveWeights(entry)
    expect_equal(wr$variantIds, paste0("chr1:", 100 * 1:3, ":A:G"))
    expect_equal(wr$weights, rep(0.05, 3))
    expect_equal(length(wr$variantIds), length(wr$weights))
    # empty topLoci -> empty pair
    empty <- fineMappingRow(
        variantIds = character(0),
        susieFit = list(),
        topLoci = data.frame(
            variant_id = character(0),
            pip = numeric(0),
            stringsAsFactors = FALSE
        )
    )
    expect_length(.fmrRowResolveWeights(empty)$variantIds, 0L)
})


test_that("FineMappingRow: getPip returns numeric(0) when topLoci is empty", {
    entry <- fineMappingRow(
        variantIds = character(0),
        susieFit = list(),
        topLoci = data.frame(
            variant_id = character(0),
            pip = numeric(0),
            stringsAsFactors = FALSE
        )
    )
    expect_equal(.fmrRowPip(entry), numeric(0))
})


test_that("FineMappingRow: getCs filters to rows in any credible set", {
    entry <- .sc_makeFineMappingRow(3) # last row has cs_95 = "susie_0"
    res <- .fmrRowCs(entry)
    expect_equal(nrow(res), 2L)
})


test_that("FineMappingRow: getCs/getTopLoci surface the directional af column", {
    # Regression: the posterior view must carry the topLoci `af` (effect-allele
    # frequency) through to .fmrRowCs() / .fmrRowTopLoci(), not drop it to NA. The
    # value is directional (0.87 > 0.5), so a folded MAF would be a bug.
    entry <- .sc_makeFineMappingRow(3)
    cs <- .fmrRowCs(entry)
    expect_true("af" %in% names(cs))
    expect_equal(unname(cs$af), rep(0.87, nrow(cs)))
    expect_false(anyNA(cs$af))

    tl <- .fmrRowTopLoci(entry, signalCutoff = 0)
    expect_true("af" %in% names(tl))
    expect_equal(unname(tl$af), rep(0.87, nrow(tl)))
})


test_that("FineMappingRow: validity errors when topLoci is missing required cols", {
    expect_error(
        fineMappingRow(
            variantIds = "chr1:100:A:G",
            susieFit = list(),
            topLoci = data.frame(other = 1, stringsAsFactors = FALSE)
        ),
        "topLoci missing required columns"
    )
})


# === cvResult slot (cross-validation payload) ===

test_that("FineMappingRow stores and returns a cvResult payload", {
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:A:G"),
        pip = c(0.8, 0.2),
        stringsAsFactors = FALSE
    )
    cv <- list(
        samplePartition = data.frame(Sample = c("s1", "s2"), Fold = c(1L, 2L)),
        prediction = list(susie_predicted = matrix(0, 2, 1)),
        performance = list(susie_performance = matrix(0, 1, 6))
    )
    e <- fineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl,
        cvResult = cv
    )
    expect_identical(.fmrPartsCvResult(e), cv)
})

test_that("FineMappingRow cvResult defaults to NULL and rejects non-list", {
    tl <- data.frame(
        variant_id = "chr1:100:A:G",
        pip = 0.5,
        stringsAsFactors = FALSE
    )
    e <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(),
        topLoci = tl
    )
    expect_null(.fmrPartsCvResult(e))
    expect_error(
        fineMappingRow(
            variantIds = "chr1:100:A:G",
            susieFit = list(),
            topLoci = tl,
            cvResult = 1:3
        ),
        "cvResult must be NULL or a list"
    )
})

# ===========================================================================
# TwasWeightsRow
# ===========================================================================

test_that("QtlFineMappingResult: rejects rows that are not row payloads", {
    expect_error(
        QtlFineMappingResult(
            study = "s1",
            context = "c1",
            trait = "t1",
            method = "susie",
            entry = list("not_an_entry")
        ),
        "must be a fine-mapping row"
    )
})


# === Tests migrated from test_showMethods.R (FineMappingRow) ===

test_that("show.FineMappingRow reports variant count and CS count", {
    e_with_cs <- .sh_makeFmEntry(n = 3, with_cs = TRUE) # 2 distinct cs > 0
    out <- capture.output(show(e_with_cs))
    expect_true(any(grepl(
        "FineMappingRow: 3 variants.*1 credible sets",
        out
    )))

    # No cs column -> 0 credible sets reported.
    tl <- data.frame(
        variant_id = c("chr9:100:A:G", "chr9:200:A:G"),
        pip = c(0.1, 0.2),
        stringsAsFactors = FALSE
    )
    e_no_cs <- fineMappingRow(
        variantIds = c("chr9:100:A:G", "chr9:200:A:G"),
        susieFit = list(),
        topLoci = tl
    )
    out_no <- capture.output(show(e_no_cs))
    expect_true(any(grepl("0 credible sets", out_no)))
})


# === getMarginalEffects maxPval filter ===

test_that("FineMappingRow: getMarginalEffects applies the maxPval filter", {
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:A:G", "chr1:300:A:G"),
        pip = c(0.9, 0.5, 0.1),
        marginal_p = c(0.001, 0.5, NA_real_),
        stringsAsFactors = FALSE
    )
    entry <- fineMappingRow(
        variantIds = tl$variant_id,
        susieFit = list(),
        topLoci = tl
    )
    out <- .fmrRowMarginalEffects(entry, maxPval = 0.01)
    # Drops the p = 0.5 row and the NA-p row; keeps only v1.
    expect_equal(nrow(out), 1L)
    expect_equal(out$variant_id, "chr1:100:A:G")
})


# === getCs empty / cs-less topLoci projections ===

test_that("FineMappingRow: getCs returns empty posterior view when topLoci is empty", {
    entry <- fineMappingRow(
        variantIds = character(0),
        susieFit = list(),
        topLoci = data.frame(
            variant_id = character(0),
            pip = numeric(0),
            stringsAsFactors = FALSE
        )
    )
    res <- .fmrRowCs(entry)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 0L)
    expect_true(all(c("variant_id", "pip") %in% names(res)))
})


test_that("FineMappingRow: getCs returns empty posterior view when the cs column is absent", {
    tl <- data.frame(
        variant_id = c("chr9:100:A:G", "chr9:200:A:G"),
        pip = c(0.1, 0.2),
        stringsAsFactors = FALSE
    )
    entry <- fineMappingRow(
        variantIds = c("chr9:100:A:G", "chr9:200:A:G"),
        susieFit = list(),
        topLoci = tl
    )
    res <- .fmrRowCs(entry) # no cs_95 column -> empty view
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 0L)
})


# === adjustPips: missing alpha + mu2-driven posterior recompute ===

test_that("FineMappingRow: adjustPips errors when susieFit lacks alpha", {
    tl <- data.frame(
        variant_id = c("chr1:100:A:G", "chr1:200:A:G"),
        pip = c(0.6, 0.4),
        stringsAsFactors = FALSE
    )
    entry <- fineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
        susieFit = list(),
        topLoci = tl
    )
    expect_error(
        adjustPips(entry, c("chr1:100:A:G", "chr1:200:A:G")),
        "no `alpha` matrix"
    )
})


test_that("FineMappingRow: adjustPips subsets mu2 and recomputes posterior_sd", {
    vids <- paste0("chr1:", 1:5, ":A:G")
    p <- length(vids)
    L <- 2L
    set.seed(101L)
    lbf <- matrix(rnorm(L * p), nrow = L, ncol = p)
    colnames(lbf) <- vids
    alpha <- lbfToAlpha(lbf)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    mu <- matrix(rnorm(L * p), L, p)
    mu2 <- mu^2 + 1 # plausible second moment (>= mean^2)
    entry <- fineMappingRow(
        variantIds = vids,
        susieFit = list(
            pip = pip,
            alpha = alpha,
            lbf_variable = lbf,
            mu = mu,
            mu2 = mu2,
            X_column_scale_factors = rep(1, p)
        ),
        topLoci = data.frame(
            variant_id = vids,
            pip = pip,
            stringsAsFactors = FALSE
        )
    )
    keep <- vids[2:4]
    adj <- adjustPips(entry, keep)
    expect_s4_class(adj, "FineMappingRow")
    # mu2 carried through the variant subsetting alongside lbf/mu.
    expect_equal(ncol(getSusieFit(adj)$mu2), 3L)
    # posterior_mean / posterior_sd recomputed from the subset alpha/mu/mu2.
    expect_equal(nrow(.fmrPartsTopLoci(adj)), 3L)
    expect_true("posterior_sd" %in% names(.fmrPartsTopLoci(adj)))
    expect_equal(length(.fmrPartsTopLoci(adj)$posterior_sd), 3L)
    expect_true(all(.fmrPartsTopLoci(adj)$posterior_sd >= 0))
})

test_that("getCs / getTopLoci minPurity filters CS by purity, independent of coverage + pip", {
    # Two 0.95 credible sets: CS1 (v1,v2) pure (0.9), CS2 (v3,v4) impure (0.3); v5 non-CS.
    vn <- paste0("chr1:", (1:5) * 100, ":A:G")
    L <- 2L
    P <- 5L
    alpha <- matrix(0.1, L, P)
    alpha[1, 1:2] <- c(0.5, 0.4)
    alpha[2, 3:4] <- c(0.5, 0.4)
    fit <- list(
        alpha = alpha,
        mu = matrix(0.3, L, P),
        mu2 = matrix(1.2, L, P),
        pip = c(0.6, 0.5, 0.4, 0.4, 0.02)
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
    e <- fineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)

    # getCs: minPurity is orthogonal to coverage -> keeps only the pure CS members
    expect_equal(nrow(.fmrRowCs(e, coverage = 0.95)), 4L)
    expect_equal(
        .fmrRowCs(e, coverage = 0.95, minPurity = 0.8)$variant_id,
        vn[1:2]
    )
    # getTopLoci: minPurity is orthogonal to the pip cutoff -> drops impure-CS
    # variants (v3,v4), keeps pure-CS (v1,v2) and non-CS (v5)
    expect_setequal(
        .fmrRowTopLoci(e, signalCutoff = 0, minPurity = 0.8)$variant_id,
        vn[c(1, 2, 5)]
    )
    expect_equal(nrow(.fmrRowTopLoci(e, signalCutoff = 0)), 5L) # default (NULL) unchanged
})


# ===========================================================================
# adjustPips: renormalize the stored alpha (prior-correct), reject
# Omega-weighted fits, handle null_weight, and recompute credible sets
# ===========================================================================

# A fit whose alpha was built from a deliberately NON-uniform prior, so the
# stored-alpha route and the rebuild-from-lbf route disagree.
.makeNonUniformPriorEntry <- function(vids, prior, L = 2L, seed = 404L) {
    p <- length(vids)
    set.seed(seed)
    lbf <- matrix(rnorm(L * p, sd = 2), L, p, dimnames = list(NULL, vids))
    w <- exp(sweep(lbf, 1L, apply(lbf, 1L, max), `-`)) * rep(prior, each = L)
    alpha <- w / rowSums(w)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    list(
        lbf = lbf,
        alpha = alpha,
        entry = fineMappingRow(
            variantIds = vids,
            susieFit = list(
                pip = pip,
                alpha = alpha,
                lbf_variable = lbf,
                pi = prior
            ),
            topLoci = data.frame(
                variant_id = vids,
                pip = pip,
                stringsAsFactors = FALSE
            )
        )
    )
}

test_that("adjustPips renormalizes stored alpha, honouring a non-uniform prior", {
    vids <- paste0("chr1:", 1:6, ":A:G")
    prior <- c(0.5, 0.2, 0.1, 0.1, 0.05, 0.05)
    built <- .makeNonUniformPriorEntry(vids, prior)
    keepIdx <- 2:5
    adj <- adjustPips(built$entry, vids[keepIdx])

    expected <- built$alpha[, keepIdx, drop = FALSE]
    expected <- expected / rowSums(expected)
    expect_equal(unname(getSusieFit(adj)$alpha), unname(expected))
    expect_equal(
        getSusieFit(adj)$pip,
        as.numeric(1 - apply(1 - expected, 2, prod))
    )
    # The old route -- rebuilding alpha from lbf_variable with a uniform prior
    # -- gets a materially different answer on the same fit.
    uniform <- lbfToAlpha(built$lbf[, keepIdx, drop = FALSE])
    expect_gt(max(abs(unname(uniform) - unname(expected))), 1e-6)
    # `pi` is carried along, subset to the retained variants.
    expect_equal(getSusieFit(adj)[["pi"]], prior[keepIdx])
})

test_that("adjustPips does not mistake `pip` for the prior `pi`", {
    # `$` on a list prefix-matches when the exact name is absent, so a `fit$pi`
    # read would silently return `fit$pip` and mis-size the subset.
    vids <- paste0("chr1:", 1:5, ":A:G")
    entry <- .makeAdjustEntry(vids)
    expect_null(.fmrPartsSusieFit(entry)[["pi"]])
    adj <- adjustPips(entry, vids[2:4])
    expect_null(getSusieFit(adj)[["pi"]])
    expect_length(getSusieFit(adj)[["pip"]], 3L)
})

test_that("adjustPips refuses an Omega-weighted (susieInf / susieAsh) fit", {
    vids <- paste0("chr1:", 1:4, ":A:G")
    entry <- .makeAdjustEntry(vids)
    fit <- .fmrPartsSusieFit(entry)
    fit$theta <- rep(0.1, length(vids))
    weighted <- fineMappingRow(
        variantIds = vids,
        susieFit = fit,
        topLoci = .fmrRowTopLoci(entry, raw = TRUE)
    )
    expect_error(adjustPips(weighted, vids[1:3]), "Omega-weighted")
})

test_that("adjustPips keeps the null_weight column in the renormalization", {
    vids <- paste0("chr1:", 1:5, ":A:G")
    p <- length(vids)
    L <- 2L
    set.seed(77L)
    # susieR appends the null as alpha column p + 1 and leaves pip length p.
    raw <- matrix(runif(L * (p + 1)), L, p + 1)
    alpha <- raw / rowSums(raw)
    pip <- as.numeric(1 - apply(1 - alpha[, seq_len(p), drop = FALSE], 2, prod))
    entry <- fineMappingRow(
        variantIds = vids,
        susieFit = list(pip = pip, alpha = alpha, null_index = p + 1L),
        topLoci = data.frame(
            variant_id = vids,
            pip = pip,
            stringsAsFactors = FALSE
        )
    )
    keepIdx <- c(1L, 3L, 5L)
    adj <- adjustPips(entry, vids[keepIdx])
    adjAlpha <- getSusieFit(adj)$alpha

    expect_equal(ncol(adjAlpha), length(keepIdx) + 1L)
    expect_equal(getSusieFit(adj)[["null_index"]], length(keepIdx) + 1L)
    expect_equal(rowSums(adjAlpha), rep(1, L))
    expect_length(getSusieFit(adj)[["pip"]], length(keepIdx))
    # Null mass survives: the retained variants alone do not sum to 1.
    expect_lt(max(rowSums(adjAlpha[, seq_along(keepIdx), drop = FALSE])), 1)
    expected <- alpha[, c(keepIdx, p + 1L), drop = FALSE]
    expect_equal(unname(adjAlpha), unname(expected / rowSums(expected)))
})

test_that("adjustPips errors when alpha's width matches neither variant count", {
    vids <- paste0("chr1:", 1:4, ":A:G")
    entry <- .makeAdjustEntry(vids)
    fit <- .fmrPartsSusieFit(entry)
    fit$alpha <- cbind(fit$alpha, fit$alpha)
    broken <- fineMappingRow(
        variantIds = vids,
        susieFit = fit,
        topLoci = .fmrRowTopLoci(entry, raw = TRUE)
    )
    expect_error(adjustPips(broken, vids[1:3]), "8 columns but the entry has 4")
})

test_that("adjustPips recomputes credible sets instead of remapping indices", {
    vids <- paste0("chr1:", 1:6, ":A:G")
    p <- length(vids)
    L <- 2L
    set.seed(303L)
    lbf <- matrix(rnorm(L * p, sd = 3), L, p, dimnames = list(NULL, vids))
    alpha <- lbfToAlpha(lbf)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    entry <- fineMappingRow(
        variantIds = vids,
        susieFit = list(
            pip = pip,
            alpha = alpha,
            lbf_variable = lbf,
            V = rep(1, L),
            # Deliberately stale: indices into the pre-subset variant set.
            sets = list(
                cs = list(L1 = c(5L, 6L), L2 = 1L),
                cs_index = 1:2,
                requested_coverage = 0.95,
                purity = data.frame(min.abs.corr = c(0.9, 0.8))
            )
        ),
        topLoci = data.frame(
            variant_id = vids,
            pip = pip,
            method = "susie",
            cs_95 = c("susie_1", rep("susie_0", p - 1L)),
            cs_95_purity = rep(0.9, p),
            stringsAsFactors = FALSE
        )
    )
    keepIdx <- 1:3
    adj <- adjustPips(entry, vids[keepIdx])
    sets <- getSusieFit(adj)$sets

    # Every recomputed index addresses a retained variant.
    expect_true(all(unlist(sets$cs) %in% seq_along(keepIdx)))
    expect_equal(sets$requested_coverage, 0.95)
    # Purity is LD-dependent: dropped rather than carried forward stale.
    expect_null(sets$purity)
    expect_true(all(is.na(.fmrRowTopLoci(adj, raw = TRUE)$cs_95_purity)))
    # Labels are rebuilt from the recomputed sets, still "<method>_<k>".
    labels <- .fmrRowTopLoci(adj, raw = TRUE)$cs_95
    expect_true(all(str_detect(labels, "^susie_[0-9]+$")))
    expect_equal(length(labels), length(keepIdx))
})

test_that("adjustPips subsets mvsusie-shaped coef rows and 3-D clfsr", {
    vids <- paste0("chr1:", 1:5, ":A:G")
    p <- length(vids)
    L <- 2L
    R <- 3L
    entry <- .makeAdjustEntry(vids)
    fit <- .fmrPartsSusieFit(entry)
    fit$coef <- matrix(seq_len(p * R), nrow = p, ncol = R)
    fit$clfsr <- array(seq_len(L * p * R), dim = c(L, p, R))
    mv <- fineMappingRow(
        variantIds = vids,
        susieFit = fit,
        topLoci = .fmrRowTopLoci(entry, raw = TRUE)
    )
    keepIdx <- c(2L, 4L)
    adj <- adjustPips(mv, vids[keepIdx])

    expect_equal(getSusieFit(adj)$coef, fit$coef[keepIdx, , drop = FALSE])
    expect_equal(dim(getSusieFit(adj)$clfsr), c(L, length(keepIdx), R))
    expect_equal(
        getSusieFit(adj)$clfsr,
        fit$clfsr[, keepIdx, , drop = FALSE]
    )
})

test_that(".fmrRowTopLoci(raw = TRUE) returns the stored table verbatim", {
    vids <- paste0("chr1:", 1:4, ":A:G")
    entry <- .makeAdjustEntry(vids)
    raw <- .fmrRowTopLoci(entry, raw = TRUE)
    expect_equal(nrow(raw), length(vids))
    expect_true(all(c("betahat", "sebetahat") %in% names(raw)))
    # The default view filters by pip and projects; raw does neither.
    expect_false(identical(names(raw), names(.fmrRowTopLoci(entry))))
})

# ===========================================================================
# getLbf (merged from test_getLbf.R when R/getLbf.R was folded
# into R/FineMappingRow.R + R/AllClasses.R)
# ===========================================================================

test_that("getLbf returns the wide variant x effect lbf matrix", {
    vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    lbf <- matrix(c(3.0, 2.5, 0.1, 0.1, 0.1, 2.0), nrow = 2, byrow = TRUE) # L=2 x p=3
    fit <- list(
        alpha = matrix(0.2, 2, 3),
        mu = matrix(0.1, 2, 3),
        mu2 = matrix(1, 2, 3),
        pip = c(0.5, 0.4, 0.3),
        lbf_variable = lbf
    )
    class(fit) <- "susie"
    e <- fineMappingRow(
        variantIds = vids,
        susieFit = fit,
        topLoci = data.frame(
            variant_id = vids,
            pip = fit$pip,
            stringsAsFactors = FALSE
        )
    )
    w <- .fmrRowLbf(e)
    expect_equal(w$variant_id, vids)
    expect_true(all(c("lbf_L1", "lbf_L2") %in% names(w)))
    expect_equal(w$lbf_L1, lbf[1, ]) # effect 1's lbf across variants
    expect_equal(w$lbf_L2, lbf[2, ]) # effect 2's lbf across variants
})

test_that("getLbf is empty when the fit carries no lbf matrix", {
    e <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(pip = 0.5),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    expect_equal(nrow(.fmrRowLbf(e)), 0L)
})

test_that("getLbf aggregates over a collection, carrying identity + NA-filling ragged effects", {
    vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    mkEntry <- function(L) {
        lbf <- matrix(seq_len(L * length(vids)), nrow = L, byrow = TRUE) # L x p
        fit <- structure(
            list(lbf_variable = lbf, pip = rep(0.3, length(vids))),
            class = "susie"
        )
        fineMappingRow(
            variantIds = vids,
            susieFit = fit,
            topLoci = data.frame(
                variant_id = vids,
                pip = rep(0.3, length(vids)),
                stringsAsFactors = FALSE
            )
        )
    }
    # Two entries with different effect counts (L = 2 vs 3) exercise the NA-fill of
    # the ragged lbf_L<k> columns described in the collection method's contract.
    res <- QtlFineMappingResult(
        study = c("s1", "s1"),
        context = c("brain", "blood"),
        trait = c("g", "g"),
        method = c("susie", "susie"),
        entry = list(mkEntry(2L), mkEntry(3L))
    )
    w <- getLbf(res)
    expect_true(all(
        c("study", "context", "trait", "method", "variant_id") %in% names(w)
    ))
    expect_equal(nrow(w), 6L) # 2 entries x 3 variants
    expect_setequal(unique(w$context), c("brain", "blood"))
    expect_true("lbf_L3" %in% names(w))
    expect_true(all(is.na(w$lbf_L3[w$context == "brain"]))) # L = 2 entry: no effect 3
    expect_true(all(!is.na(w$lbf_L3[w$context == "blood"]))) # L = 3 entry: present
})

# ===========================================================================
# getCredibleSetSummary (merged from test_credibleSetSummary.R
# when R/credibleSetSummary.R was folded in)
# ===========================================================================

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
    e <- fineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)

    s <- .fmrRowCredibleSetSummary(e)
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
    e <- fineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)
    expect_equal(nrow(.fmrRowCredibleSetSummary(e)), 0L)
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
    e <- fineMappingRow(variantIds = vn, susieFit = list(x = 1), topLoci = tl)
    expect_equal(nrow(.fmrRowCredibleSetSummary(e, coverage = 0.99)), 0L)
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
    e <- fineMappingRow(
        variantIds = vn,
        susieFit = list(trimmed = TRUE),
        topLoci = tl
    )
    s <- .fmrRowCredibleSetSummary(e)
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
    e <- fineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)
    expect_true(is.na(.fmrRowCredibleSetSummary(e)$cs_log_bf))
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
        fineMappingRow(variantIds = vn, susieFit = list(x = 1), topLoci = tl)
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

# ===========================================================================
# fsusieCredibleBand / fsusieAffectedRegions (merged from
# test_fsusieAccessors.R when R/fsusieAccessors.R was folded in)
# ===========================================================================

.fsa_makeFit <- function() {
    set.seed(1)
    n <- 150L
    p <- 24L
    J <- 16L
    X <- matrix(
        rnorm(n * p),
        n,
        p,
        dimnames = list(
            paste0("s", seq_len(n)),
            paste0("chr1:", (seq_len(p)) * 100, ":A:G")
        )
    )
    b1 <- sin(seq(0, 2 * pi, length.out = J))
    b2 <- cos(seq(0, pi, length.out = J))
    Y <- X[, 3] %o% b1 + X[, 10] %o% b2 + matrix(rnorm(n * J, sd = 0.3), n, J)
    colnames(Y) <- paste0("f", seq_len(J))
    suppressWarnings(fsusieR::susiF(
        X = X,
        Y = Y,
        pos = seq_len(J),
        L = 5,
        post_processing = "none",
        verbose = FALSE
    ))
}
.fsa_entry <- function(fit, purityVal = 0.85) {
    vids <- names(fit$csd_X)
    pipv <- if (!is.null(fit$pip)) {
        as.numeric(fit$pip)
    } else {
        rep(0.1, length(vids))
    }
    # Canonical CS labels + purity from the fit's CS membership, mirroring what
    # fsusieWrapper writes into a real FMR's topLoci (cs_95 / cs_95_purity).
    cs95 <- rep("fsusie_0", length(vids))
    for (l in seq_along(fit$cs)) {
        idx <- fit$cs[[l]]
        if (is.numeric(idx)) {
            idx <- as.integer(idx)
        }
        cs95[idx] <- paste0("fsusie_", l)
    }
    tl <- data.frame(
        variant_id = vids,
        pip = pipv,
        cs_95 = cs95,
        cs_95_purity = ifelse(cs95 != "fsusie_0", purityVal, 0),
        stringsAsFactors = FALSE
    )
    fineMappingRow(variantIds = vids, susieFit = fit, topLoci = tl)
}

test_that("fsusieCredibleBand returns a long effect + band table (lower <= effect <= upper)", {
    skip_if_not_installed("fsusieR")
    skip_if_not_installed("wavethresh")
    cb <- .fmrRowFsusieCredibleBand(.fsa_entry(.fsa_makeFit()))
    expect_named(cb, c("cs", "chrom", "pos", "effect", "lower", "upper"))
    expect_gt(nrow(cb), 0)
    expect_true(all(
        cb$lower <= cb$effect + 1e-9 & cb$effect <= cb$upper + 1e-9
    ))
    expect_true(all(grepl("^fsusie_", cb$cs)))
})

test_that("fsusieAffectedRegions returns a GRanges with cs / purity / direction", {
    skip_if_not_installed("fsusieR")
    skip_if_not_installed("wavethresh")
    gr <- .fmrRowFsusieAffectedRegions(.fsa_entry(
        .fsa_makeFit(),
        purityVal = 0.85
    ))
    expect_s4_class(gr, "GRanges")
    expect_gt(length(gr), 0)
    expect_true(all(
        c("cs", "purity", "direction") %in% names(S4Vectors::mcols(gr))
    ))
    expect_true(all(S4Vectors::mcols(gr)$direction %in% c("pos", "neg", NA)))
    # purity sourced from the entry's cs_95_purity (matched by CS variant membership),
    # NOT the (nonexistent) fit$purity slot
    expect_true(all(S4Vectors::mcols(gr)$purity == 0.85))
    expect_true(all(grepl("^fsusie_", S4Vectors::mcols(gr)$cs)))
})

test_that("fsusie accessors degrade to empty for a non-fSuSiE / trimmed fit (no wavelet slots)", {
    e <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(pip = 0.5),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    expect_equal(nrow(.fmrRowFsusieCredibleBand(e)), 0L)
    expect_equal(length(.fmrRowFsusieAffectedRegions(e)), 0L)
})

test_that("fsusieCredibleBand + fsusieAffectedRegions aggregate across a collection", {
    skip_if_not_installed("fsusieR")
    skip_if_not_installed("wavethresh")
    fit <- .fsa_makeFit() # deterministic; reuse for both rows
    res <- QtlFineMappingResult(
        study = c("s", "s"),
        context = c("brain", "blood"),
        trait = c("g", "g"),
        method = c("fsusie", "fsusie"),
        entry = list(.fsa_entry(fit), .fsa_entry(fit))
    )

    cb <- fsusieCredibleBand(res)
    expect_true(all(
        c(
            "study",
            "context",
            "trait",
            "method",
            "cs",
            "effect",
            "lower",
            "upper"
        ) %in%
            names(cb)
    ))
    expect_gt(nrow(cb), 0)
    expect_setequal(unique(cb$context), c("brain", "blood"))

    gr <- fsusieAffectedRegions(res)
    expect_s4_class(gr, "GRanges")
    expect_gt(length(gr), 0)
    expect_true(all(
        c("study", "context", "trait", "cs", "purity", "direction") %in%
            names(S4Vectors::mcols(gr))
    ))
    expect_setequal(unique(S4Vectors::mcols(gr)$context), c("brain", "blood"))
})

test_that("fsusieAffectedRegions on a collection of non-fSuSiE entries is an empty GRanges", {
    e <- fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(pip = 0.5),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    res <- QtlFineMappingResult(
        study = "s",
        context = "c",
        trait = "t",
        method = "susie",
        entry = list(e)
    )
    expect_equal(length(fsusieAffectedRegions(res)), 0L)
})

test_that(".fsusieChrom returns NA when the fit has no named variants", {
    expect_true(is.na(pecotmr:::.fsusieChrom(list(csd_X = c(1, 2, 3))))) # unnamed
    expect_true(is.na(pecotmr:::.fsusieChrom(list()))) # no csd_X
})

test_that(".fsusieCsMapFromTopLoci falls back to index labels when topLoci lacks CS columns", {
    fit <- list(
        cs = list(1L, 2L),
        csd_X = stats::setNames(1:3, paste0("chr1:", 1:3 * 100, ":A:G"))
    )
    m <- pecotmr:::.fsusieCsMapFromTopLoci(fit, topLoci = NULL)
    expect_equal(unname(m$label), c("fsusie_1", "fsusie_2"))
    expect_true(all(is.na(m$purity)))
    # a topLoci missing the required cs_95 / cs_95_purity columns takes the same path
    m2 <- pecotmr:::.fsusieCsMapFromTopLoci(
        fit,
        topLoci = data.frame(variant_id = "x")
    )
    expect_equal(unname(m2$label), c("fsusie_1", "fsusie_2"))
})

# A minimal object that passes .isFsusieFit (the wavelet-slot gate) so the band /
# affected-region degradation paths can be exercised with the (untrimmed) band
# computation replaced by a mock returning a controlled, degenerate fit.
.fsa_fakeFit <- function() {
    list(fitted_wc2 = 1, fitted_func = list(1), outing_grid = 1:2, alpha = 1)
}
.fsa_bareEntry <- function() {
    fineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = .fsa_fakeFit(),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
}

test_that("fsusieCredibleBand skips NULL effect bands and degrades to empty", {
    skip_if_not_installed("fsusieR")
    testthat::local_mocked_bindings(
        .fsusiePopulateCredibleBand = function(fit) {
            list(
                outing_grid = c(100, 200),
                csd_X = stats::setNames(1:2, c("chr1:100:A:G", "chr1:200:A:G")),
                cred_band = list(NULL, NULL),
                fitted_func = list(NULL, NULL)
            )
        },
        .package = "pecotmr"
    )
    expect_equal(nrow(.fmrRowFsusieCredibleBand(.fsa_bareEntry())), 0L)
})

test_that("fsusieAffectedRegions degrades to empty when affected_reg finds nothing", {
    skip_if_not_installed("fsusieR")
    testthat::local_mocked_bindings(
        .fsusiePopulateCredibleBand = function(fit) {
            list(
                outing_grid = c(100, 200),
                csd_X = stats::setNames(1:2, c("chr1:100:A:G", "chr1:200:A:G")),
                fitted_func = list(c(1, 1))
            )
        },
        .package = "pecotmr"
    )
    testthat::local_mocked_bindings(
        affected_reg = function(...) NULL,
        .package = "fsusieR"
    )
    expect_equal(length(.fmrRowFsusieAffectedRegions(.fsa_bareEntry())), 0L)
})

test_that("fsusieAffectedRegions yields NA direction for a region outside the grid", {
    skip_if_not_installed("fsusieR")
    testthat::local_mocked_bindings(
        .fsusiePopulateCredibleBand = function(fit) {
            list(
                outing_grid = c(100, 200),
                csd_X = stats::setNames(1:2, c("chr1:100:A:G", "chr1:200:A:G")),
                fitted_func = list(c(1, 1))
            )
        },
        .package = "pecotmr"
    )
    # Region [500, 600] lies outside the grid [100, 200]: no grid points fall in it,
    # so the per-region effect-direction summary is NA.
    testthat::local_mocked_bindings(
        affected_reg = function(...) {
            data.frame(CS = 1L, Start = 500, End = 600)
        },
        .package = "fsusieR"
    )
    gr <- .fmrRowFsusieAffectedRegions(.fsa_bareEntry())
    expect_s4_class(gr, "GRanges")
    expect_equal(length(gr), 1L)
    expect_true(is.na(S4Vectors::mcols(gr)$direction))
})

test_that("getCs warns when the coverage has no purity column to filter on", {
    # minPurity is orthogonal to coverage, so asking for it at a coverage that
    # produced no credible sets leaves nothing to filter. Skipping silently
    # would look like "the filter ran and kept everything".
    vn <- str_c("chr1:", (1:5) * 100, ":A:G")
    alpha <- matrix(0.1, 2L, 5L)
    alpha[1, 1:2] <- c(0.5, 0.4)
    alpha[2, 3:4] <- c(0.5, 0.4)
    fit <- list(
        alpha = alpha,
        mu = matrix(0.3, 2L, 5L),
        mu2 = matrix(1.2, 2L, 5L),
        pip = c(0.6, 0.5, 0.4, 0.4, 0.02)
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
    # Dropped by hand: buildTopLoci ALWAYS emits cs_<cov>_purity, even when
    # the fitter reported no purity table, so this branch is unreachable for
    # anything the package builds itself. It guards a topLoci supplied from
    # outside -- and a coverage with no sets returns earlier, so that route
    # cannot reach it either.
    tl$cs_95_purity <- NULL
    e <- fineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)

    expect_warning(
        .fmrRowCs(e, coverage = 0.95, minPurity = 0.8),
        "no purity column"
    )
})

# ===========================================================================
# Construction guards: topLoci must line up with the variants
#
# topLoci's columns are taken POSITIONALLY, so a table of the right length in
# the wrong order silently mis-assigns every column. These guards catch that;
# none had a test.
# ===========================================================================

.fmr_vids <- function(n = 3L) str_c("chr1:", 100L * seq_len(n), ":A:G")

test_that("topLoci must have one row per supplied variant", {
    vids <- .fmr_vids(3L)
    expect_error(
        fineMappingRow(
            variantIds = vids,
            susieFit = list(fake = TRUE),
            topLoci = data.frame(
                variant_id = vids[1:2],
                pip = c(0.4, 0.3),
                stringsAsFactors = FALSE
            )
        ),
        "aligned row-for-row"
    )
})

test_that("topLoci must list the variants in the same order", {
    # Same set, different order: length checks pass, so only an explicit
    # order check catches it.
    vids <- .fmr_vids(3L)
    expect_error(
        fineMappingRow(
            variantIds = vids,
            susieFit = list(fake = TRUE),
            topLoci = data.frame(
                variant_id = rev(vids),
                pip = c(0.4, 0.3, 0.2),
                stringsAsFactors = FALSE
            )
        ),
        "first differ at position 1"
    )
})

# ===========================================================================
# adjustPips: fits this code cannot slice safely
# ===========================================================================

.fmr_susieFit <- function(p = 5L, L = 2L, seed = 7L) {
    set.seed(seed)
    lbf <- matrix(rnorm(L * p, sd = 3), L, p)
    alpha <- lbfToAlpha(lbf)
    list(
        alpha = alpha,
        mu = matrix(0.3, L, p),
        mu2 = matrix(1.2, L, p),
        lbf_variable = lbf,
        V = rep(1, L),
        pip = as.numeric(1 - apply(1 - alpha, 2, prod))
    )
}

test_that("adjustPips rejects a null_index that is not the last column", {
    # The null column carries the "no signal here" mass; slicing assumes it
    # sits last, so any other position is a layout this code cannot honour.
    fit <- .fmr_susieFit()
    expect_equal(.adjustPipsNullIdx(fit, 5L), 5L)
    fit$null_index <- 2L
    expect_error(
        .adjustPipsNullIdx(fit, 5L),
        "null column must be the last of 5 alpha columns"
    )
})

test_that("adjustPips refuses effects with no mass on the retained variants", {
    # Renormalizing zero mass is not "a small number" -- it is undefined, and
    # returning it would mix adjusted and meaningless posteriors.
    alpha <- matrix(c(1, 0, 0, 0, 0, 1), nrow = 2L, byrow = TRUE)
    expect_error(
        .adjustPipsRenormAlpha(alpha, cols = 2L),
        "carry no .*posterior mass on the retained variants"
    )
})

# ===========================================================================
# Credible-set summary completeness (ported from PR #577)
#
# The summary is built from the fit's own sets$cs when it carries them, so
# every set the fit found gets a row -- including one whose members were all
# relabelled into a smaller overlapping set in the per-variant column.
# ===========================================================================

test_that("every fit credible set gets a row, with its true membership", {
    vn <- str_c("chr1:", (1:5) * 100, ":A:G")
    alpha <- matrix(0, 2, 5, dimnames = list(c("L1", "L2"), vn))
    alpha[1, 1:3] <- 0.3
    alpha[2, 2] <- 0.9
    fit <- list(
        pip = c(0.3, 0.9, 0.3, 0.01, 0.01),
        alpha = alpha,
        V = c(0.1, 0.2),
        lbf_variable = matrix(
            c(1, 1, 1, 0, 0, 0, 2, 0, 0, 0), 2, 5,
            byrow = TRUE, dimnames = list(c("L1", "L2"), vn)
        ),
        sets = list(
            # L2 = {2} sits inside L1 = {1,2,3}: v2 tags to the smaller L2 in
            # the per-variant column, so L1 loses it there but not in the fit.
            cs = list(L1 = c(1L, 2L, 3L), L2 = c(2L)),
            purity = data.frame(
                min.abs.corr = c(0.8, 1.0),
                mean.abs.corr = c(0.85, 1.0),
                row.names = c("L1", "L2")
            ),
            requested_coverage = 0.95
        )
    )
    class(fit) <- "susie"
    tl <- data.frame(
        variant_id = vn,
        pip = fit$pip,
        logBF = c(1, 2, 1, 0, 0),
        cs_95 = c("susie_1", "susie_2", "susie_1", "susie_0", "susie_0"),
        cs_95_purity = c(0.8, 1.0, 0.8, 0, 0),
        stringsAsFactors = FALSE
    )
    s <- .csSummaryFit(tl, fit, 0.95)
    s <- s[order(s$effect_id), , drop = FALSE]
    expect_equal(nrow(s), 2L)
    expect_equal(s$cs, c("susie_1", "susie_2"))
    expect_equal(s$effect_id, c("L1", "L2"))
    # L1 keeps its full size 3 even though only v1/v3 carry its column tag.
    expect_equal(s$n_variants, c(3L, 1L))
    # cs_pip is summed over the TRUE membership, so v2 counts toward L1.
    expect_equal(s$cs_pip, c(1.5, 0.9), tolerance = 1e-9)
})

test_that("the summary carries the fit's true (gapped) effect index", {
    # sets$cs is ordered L2, L1 (e.g. by purity): the size-3 set must stay L2,
    # not be renumbered to position 1.
    vn <- str_c("chr1:", (1:4) * 100, ":A:G")
    fit <- list(
        pip = c(0.3, 0.3, 0.3, 0.6),
        alpha = matrix(0, 2, 4, dimnames = list(c("L2", "L1"), vn)),
        V = c(0.1, 0.2),
        sets = list(
            cs = list(L2 = c(1L, 2L, 3L), L1 = c(4L)),
            purity = data.frame(
                min.abs.corr = c(0.9, 0.5),
                row.names = c("L2", "L1")
            ),
            requested_coverage = 0.95
        )
    )
    class(fit) <- "susie"
    tl <- data.frame(
        variant_id = vn, pip = fit$pip, logBF = c(1, 1, 1, 2),
        cs_95 = c("susie_2", "susie_2", "susie_2", "susie_1"),
        stringsAsFactors = FALSE
    )
    s <- .csSummaryFit(tl, fit, 0.95)
    big <- s[s$n_variants == 3L, , drop = FALSE]
    expect_equal(big$effect_id, "L2")
    expect_equal(big$cs, "susie_2")
})

test_that("a fit with no stored sets falls back to the per-variant column", {
    # Minimal / hand-built fits carry no sets$cs; the column is then the only
    # source of membership there is.
    vn <- str_c("chr1:", (1:3) * 100, ":A:G")
    fit <- list(pip = c(0.5, 0.4, 0.1))
    class(fit) <- "susie"
    tl <- data.frame(
        variant_id = vn, pip = fit$pip, logBF = c(1, 1, 0),
        cs_95 = c("susie_1", "susie_1", "susie_0"),
        stringsAsFactors = FALSE
    )
    s <- .csSummaryFit(tl, fit, 0.95)
    expect_equal(nrow(s), 1L)
    expect_equal(s$n_variants, 2L)
    expect_equal(s$cs, "susie_1")
})
