# Tests for FineMappingRow (S4 class)

# ===========================================================================
# A row as GRanges, and adjustPips on a kept subset
# ===========================================================================

test_that("as(row, 'GRanges') is the element a collection stores", {
    ids <- c("chr1:100:A:G", "chr1:200:C:T")
    row <- FineMappingRow(
        variantIds = ids,
        susieFit = list(),
        topLoci = data.frame(variant_id = ids, pip = c(0.7, 0.2))
    )
    gr <- as(row, "GRanges")
    expect_s4_class(gr, "GRanges")
    expect_identical(gr$pip, c(0.7, 0.2))
    x <- QtlFineMappingResult(
        studyName = "s",
        context = "c",
        trait = "t",
        method = "susie",
        entry = list(row)
    )
    expect_identical(unname(x[[1L]]), unname(gr))
})

# =============================================================================
# adjustPips: renormalize PIPs on a kept subset of variants
# =============================================================================

test_that("fold fits become fold weights whatever each fold kept", {
    noFolds <- list(metrics = c(rsq = 0.2))
    expect_identical(pecotmr:::.cvPayloadFoldWeights(noFolds, "v1"), noFolds)
    expect_null(pecotmr:::.foldFitWeights(NULL))
    coef <- matrix(c(0.5, -0.2), 2L, 1L, dimnames = list(c("v1", "v2"), "y"))
    expect_identical(
        pecotmr:::.foldFitWeights(list(coef = coef)),
        c(v1 = 0.5, v2 = -0.2)
    )
    expect_identical(
        pecotmr:::.alignFoldWeights(NULL, c("v1", "v2")),
        c(NA_real_, NA_real_)
    )
})

test_that("adjustPips renormalizes PIPs on a kept FineMappingRow subset", {
    vids <- paste0("chr1:", 1:6, ":A:G")
    entry <- .makeAdjustEntry(vids)
    keep <- vids[2:5]
    adj <- adjustPips(entry, keep)
    expect_s4_class(adj, "FineMappingRow")
    expect_equal(variantIds(adj), keep)
    expect_equal(ncol(susieFit(adj)$lbf_variable), 4)
    # Renormalized: each effect's alpha row sums to 1 (when row has any signal)
    expect_true(all(abs(rowSums(susieFit(adj)$alpha) - 1) < 1e-10))
    # PIPs match topLoci
    expect_equal(.fmrPartsTopLoci(adj)$pip, susieFit(adj)$pip)
    # PIPs change under renormalization
    origPips <- .fmrRowPip(entry)
    expect_false(identical(unname(origPips[keep]), susieFit(adj)$pip))
})

test_that("adjustPips errors when the intersection is empty", {
    vids <- paste0("chr1:", 1:4, ":A:G")
    entry <- .makeAdjustEntry(vids)
    expect_error(
        adjustPips(entry, paste0("chr2:", 1:4, ":A:G")),
        "intersection.*empty"
    )
})

test_that("adjustPips tolerates a chr-prefix difference with keepVariants", {
    vids <- paste0("chr1:", 1:6, ":A:G")
    entry <- .makeAdjustEntry(vids)
    keep <- paste0("1:", 2:5, ":A:G") # same variants, no "chr" prefix
    adj <- adjustPips(entry, keep)
    expect_s4_class(adj, "FineMappingRow")
    expect_equal(variantIds(adj), vids[2:5]) # entry keeps its own labels
    expect_equal(ncol(susieFit(adj)$lbf_variable), 4)
})

# === Tests migrated from test_s4Constructors.R (FineMappingRow) ===

test_that("FineMappingRow stores its slots and the accessors return them", {
    tl <- .sc_makeTopLoci(3)
    tl$variant_id <- sprintf("chr9:%d:A:G", 100L * (1:3))
    entry <- FineMappingRow(
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
    # topLoci returns the projected posterior view, not the raw slot
    out <- .fmrRowTopLoci(entry, signalCutoff = 0)
    expect_equal(out$variant_id, sprintf("chr9:%d:A:G", 100L * (1:3)))
})

test_that("FineMappingRow: pip returns named pip vector keyed by variant_id", {
    entry <- .sc_makeFineMappingRow(3)
    pip <- .fmrRowPip(entry)
    expect_equal(length(pip), 3L)
    expect_equal(names(pip), paste0("chr1:", 100 * 1:3, ":A:G"))
})

test_that("resolveWeights gives the row's posterior effect by variant_id", {
    entry <- .sc_makeFineMappingRow(3) # topLoci posterior_mean = 0.05 for all
    wr <- .fmrRowResolveWeights(entry)
    expect_equal(wr$variantIds, paste0("chr1:", 100 * 1:3, ":A:G"))
    expect_equal(wr$weights, rep(0.05, 3))
    expect_equal(length(wr$variantIds), length(wr$weights))
    # empty topLoci -> empty pair
    empty <- FineMappingRow(
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

test_that("FineMappingRow: pip returns numeric(0) when topLoci is empty", {
    entry <- FineMappingRow(
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

test_that("FineMappingRow: credibleSets filters to rows in any credible set", {
    entry <- .sc_makeFineMappingRow(3) # last row has cs_95 = "susie_0"
    res <- .fmrRowCs(entry)
    expect_equal(nrow(res), 2L)
})

test_that("a row's credibleSets and topLoci carry the directional af", {
    # Regression: the posterior view must carry the topLoci `af` (effect-allele
    # frequency) through to .fmrRowCs() / .fmrRowTopLoci(), not drop it to NA.
    # The value is directional (0.87 > 0.5), so a folded MAF would be a bug.
    entry <- .sc_makeFineMappingRow(3)
    cs <- .fmrRowCs(entry)
    expect_true("af" %in% names(cs))
    expect_equal(unname(cs$af), rep(0.87, nrow(cs)))
    expect_false(anyNA(cs$af))

    tl <- .fmrRowTopLoci(entry, signalCutoff = 0)
    expect_true("af" %in% names(tl))
    expect_equal(unname(tl$af), rep(0.87, nrow(tl)))
})

test_that("FineMappingRow validity errors on topLoci missing a column", {
    expect_error(
        FineMappingRow(
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
    e <- FineMappingRow(
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
    e <- FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = list(),
        topLoci = tl
    )
    expect_null(.fmrPartsCvResult(e))
    expect_error(
        FineMappingRow(
            variantIds = "chr1:100:A:G",
            susieFit = list(),
            topLoci = tl,
            cvResult = 1:3
        ),
        "cvResult.*Must be of type 'list' \\(or 'NULL'\\)"
    )
})

# ===========================================================================
# show(), credibleSets and adjustPips on one row
# ===========================================================================

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
    e_no_cs <- FineMappingRow(
        variantIds = c("chr9:100:A:G", "chr9:200:A:G"),
        susieFit = list(),
        topLoci = tl
    )
    out_no <- capture.output(show(e_no_cs))
    expect_true(any(grepl("0 credible sets", out_no)))

    # An empty row has no variants and so no sets.
    empty <- FineMappingRow(
        variantIds = character(0),
        susieFit = list(pip = numeric(0)),
        topLoci = data.frame(
            variant_id = character(0),
            pip = numeric(0),
            stringsAsFactors = FALSE
        )
    )
    expect_output(show(empty), "0 variants, 0 credible sets")
})

# === credibleSets empty / cs-less topLoci projections ===

test_that("credibleSets of a row with empty topLoci is an empty view", {
    entry <- FineMappingRow(
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

test_that("credibleSets of a row without its cs column is an empty view", {
    tl <- data.frame(
        variant_id = c("chr9:100:A:G", "chr9:200:A:G"),
        pip = c(0.1, 0.2),
        stringsAsFactors = FALSE
    )
    entry <- FineMappingRow(
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
    entry <- FineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:A:G"),
        susieFit = list(),
        topLoci = tl
    )
    expect_error(
        adjustPips(entry, c("chr1:100:A:G", "chr1:200:A:G")),
        "no `alpha` matrix"
    )
})

test_that("adjustPips on a row subsets mu2 and recomputes posterior_sd", {
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
    entry <- FineMappingRow(
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
    expect_equal(ncol(susieFit(adj)$mu2), 3L)
    # posterior_mean / posterior_sd recomputed from the subset alpha/mu/mu2,
    # the single effects' variances adding.
    tl <- .fmrPartsTopLoci(adj)
    sub <- susieFit(adj)
    expect_equal(nrow(tl), 3L)
    expect_equal(tl$posterior_mean, as.numeric(colSums(sub$alpha * sub$mu)))
    expect_equal(
        tl$posterior_sd,
        as.numeric(sqrt(colSums(
            sub$alpha * sub$mu2 - (sub$alpha * sub$mu)^2
        )))
    )
})

test_that("adjustPips reports the posterior on topLoci's scale", {
    # mu and mu2 are stored for the standardized X, and buildTopLoci() divides
    # by X_column_scale_factors. The recomputation did not, so adjusted rows
    # switched to standardized units.
    vids <- paste0("chr1:", 1:4, ":A:G")
    alpha <- matrix(c(0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.7, 0.1), 2L, byrow = TRUE)
    mu <- matrix(c(1, 2, 3, 4, -1, -2, -3, -4), 2L, byrow = TRUE)
    mu2 <- mu^2 + 0.5
    sf <- c(2, 4, 5, 10)
    entry <- FineMappingRow(
        variantIds = vids,
        susieFit = list(
            pip = as.numeric(1 - apply(1 - alpha, 2, prod)),
            alpha = alpha,
            mu = mu,
            mu2 = mu2,
            X_column_scale_factors = sf
        ),
        topLoci = data.frame(variant_id = vids, pip = 0.5)
    )
    tl <- .fmrPartsTopLoci(adjustPips(entry, vids[1:3]))
    sub <- susieFit(adjustPips(entry, vids[1:3]))
    expect_equal(
        tl$posterior_mean,
        as.numeric(colSums(sub$alpha * sub$mu) / sf[1:3])
    )
    expect_equal(
        tl$posterior_sd,
        as.numeric(
            sqrt(colSums(sub$alpha * sub$mu2 - (sub$alpha * sub$mu)^2)) /
                sf[1:3]
        )
    )
})

test_that("minPurity filters CSs by purity alone, not coverage or pip", {
    # Two 0.95 credible sets: CS1 (v1,v2) pure (0.9), CS2 (v3,v4) impure (0.3);
    # v5 non-CS.
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
    e <- FineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)

    # credibleSets: minPurity is orthogonal to coverage -> keeps only the pure
    # CS members
    expect_equal(nrow(.fmrRowCs(e, coverage = 0.95)), 4L)
    expect_equal(
        .fmrRowCs(e, coverage = 0.95, minPurity = 0.8)$variant_id,
        vn[1:2]
    )
    # topLoci: minPurity is orthogonal to the pip cutoff -> drops impure-CS
    # variants (v3,v4), keeps pure-CS (v1,v2) and non-CS (v5)
    expect_setequal(
        .fmrRowTopLoci(e, signalCutoff = 0, minPurity = 0.8)$variant_id,
        vn[c(1, 2, 5)]
    )
    # default (NULL) unchanged
    expect_equal(nrow(.fmrRowTopLoci(e, signalCutoff = 0)), 5L)
})

# ===========================================================================
# adjustPips: renormalize the stored alpha (prior-correct), reject
# Omega-weighted fits, handle null_weight, and recompute credible sets
# ===========================================================================

# A fit whose alpha was built from a deliberately NON-uniform prior, so the
# stored-alpha route and the rebuild-from-lbf route disagree.
.makeNonUniformPriorEntry <- function(
    vids,
    prior,
    maxNumSingleEffects = 2L,
    seed = 404L
) {
    p <- length(vids)
    set.seed(seed)
    lbf <- matrix(
        rnorm(maxNumSingleEffects * p, sd = 2),
        maxNumSingleEffects,
        p,
        dimnames = list(NULL, vids)
    )
    w <- exp(sweep(lbf, 1L, apply(lbf, 1L, max), `-`)) *
        rep(prior, each = maxNumSingleEffects)
    alpha <- w / rowSums(w)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    list(
        lbf = lbf,
        alpha = alpha,
        entry = FineMappingRow(
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

test_that("adjustPips renormalizes alpha, honouring a non-uniform prior", {
    vids <- paste0("chr1:", 1:6, ":A:G")
    prior <- c(0.5, 0.2, 0.1, 0.1, 0.05, 0.05)
    built <- .makeNonUniformPriorEntry(vids, prior)
    keepIdx <- 2:5
    adj <- adjustPips(built$entry, vids[keepIdx])

    expected <- built$alpha[, keepIdx, drop = FALSE]
    expected <- expected / rowSums(expected)
    expect_equal(unname(susieFit(adj)$alpha), unname(expected))
    expect_equal(
        susieFit(adj)$pip,
        as.numeric(1 - apply(1 - expected, 2, prod))
    )
    # The old route -- rebuilding alpha from lbf_variable with a uniform prior
    # -- gets a materially different answer on the same fit.
    uniform <- lbfToAlpha(built$lbf[, keepIdx, drop = FALSE])
    expect_gt(max(abs(unname(uniform) - unname(expected))), 1e-6)
    # `pi` is carried along, subset to the retained variants.
    expect_equal(susieFit(adj)[["pi"]], prior[keepIdx])
})

test_that("adjustPips does not mistake `pip` for the prior `pi`", {
    # `$` on a list prefix-matches when the exact name is absent, so a `fit$pi`
    # read would silently return `fit$pip` and mis-size the subset.
    vids <- paste0("chr1:", 1:5, ":A:G")
    entry <- .makeAdjustEntry(vids)
    expect_null(.fmrPartsSusieFit(entry)[["pi"]])
    adj <- adjustPips(entry, vids[2:4])
    expect_null(susieFit(adj)[["pi"]])
    expect_length(susieFit(adj)[["pip"]], 3L)
})

test_that("adjustPips refuses an Omega-weighted (susieInf / susieAsh) fit", {
    vids <- paste0("chr1:", 1:4, ":A:G")
    entry <- .makeAdjustEntry(vids)
    fit <- .fmrPartsSusieFit(entry)
    fit$theta <- rep(0.1, length(vids))
    weighted <- FineMappingRow(
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
    entry <- FineMappingRow(
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
    adjAlpha <- susieFit(adj)$alpha

    expect_equal(ncol(adjAlpha), length(keepIdx) + 1L)
    expect_equal(susieFit(adj)[["null_index"]], length(keepIdx) + 1L)
    expect_equal(rowSums(adjAlpha), rep(1, L))
    expect_length(susieFit(adj)[["pip"]], length(keepIdx))
    # Null mass survives: the retained variants alone do not sum to 1.
    expect_lt(max(rowSums(adjAlpha[, seq_along(keepIdx), drop = FALSE])), 1)
    expected <- alpha[, c(keepIdx, p + 1L), drop = FALSE]
    expect_equal(unname(adjAlpha), unname(expected / rowSums(expected)))
})

test_that("adjustPips errors when alpha's width matches no variant count", {
    vids <- paste0("chr1:", 1:4, ":A:G")
    entry <- .makeAdjustEntry(vids)
    fit <- .fmrPartsSusieFit(entry)
    fit$alpha <- cbind(fit$alpha, fit$alpha)
    broken <- FineMappingRow(
        variantIds = vids,
        susieFit = fit,
        topLoci = .fmrRowTopLoci(entry, raw = TRUE)
    )
    expect_error(adjustPips(broken, vids[1:3]), "8 columns but the entry has 4")
})

test_that("adjustPips restricts a narrowed row's fit by the fit's own ids", {
    # The row holds variants 3:6 of the six its fit was made with. The kept
    # ones are found among the fit's columns by the names the fit gives them,
    # not by their position in the row, which would pick the fit's columns
    # 1:2.
    vids <- paste0("chr1:", 1:6, ":A:G")
    entry <- .makeAdjustEntry(vids)
    narrowed <- methods::new(
        "FineMappingRow",
        variants = variants(entry)[3:6],
        susieFit = susieFit(entry)
    )
    adj <- adjustPips(narrowed, vids[3:4])
    alpha <- susieFit(entry)$alpha[, 3:4]
    expect_equal(variantIds(adj), vids[3:4])
    expect_equal(unname(susieFit(adj)$alpha), unname(alpha / rowSums(alpha)))
    expect_length(susieFit(adj)$pip, 2L)
})

test_that("adjustPips refuses a row holding a variant its fit never had", {
    fitted <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    expect_equal(
        pecotmr:::.adjustPipsFitIdx(fitted[c(3L, 1L)], fitted),
        c(3L, 1L)
    )
    expect_error(
        pecotmr:::.adjustPipsFitIdx(c(fitted[[1L]], "chr1:900:A:G"), fitted),
        "1 of the row's variants are not among the 3 its fit was made with"
    )
})

test_that("adjustPips recomputes credible sets instead of remapping indices", {
    vids <- paste0("chr1:", 1:6, ":A:G")
    p <- length(vids)
    L <- 2L
    set.seed(303L)
    lbf <- matrix(rnorm(L * p, sd = 3), L, p, dimnames = list(NULL, vids))
    alpha <- lbfToAlpha(lbf)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    entry <- FineMappingRow(
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
    sets <- susieFit(adj)$sets

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
    mv <- FineMappingRow(
        variantIds = vids,
        susieFit = fit,
        topLoci = .fmrRowTopLoci(entry, raw = TRUE)
    )
    keepIdx <- c(2L, 4L)
    adj <- adjustPips(mv, vids[keepIdx])

    expect_equal(susieFit(adj)$coef, fit$coef[keepIdx, , drop = FALSE])
    expect_equal(dim(susieFit(adj)$clfsr), c(L, length(keepIdx), R))
    expect_equal(
        susieFit(adj)$clfsr,
        fit$clfsr[, keepIdx, , drop = FALSE]
    )
})

# ===========================================================================
# lbf (merged from test_getLbf.R when R/getLbf.R was folded
# into R/FineMappingRow.R + R/AllClasses.R)
# ===========================================================================

test_that("lbf returns the wide variant x effect lbf matrix", {
    vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    # L=2 x p=3
    lbf <- matrix(c(3.0, 2.5, 0.1, 0.1, 0.1, 2.0), nrow = 2, byrow = TRUE)
    fit <- list(
        alpha = matrix(0.2, 2, 3),
        mu = matrix(0.1, 2, 3),
        mu2 = matrix(1, 2, 3),
        pip = c(0.5, 0.4, 0.3),
        lbf_variable = lbf
    )
    class(fit) <- "susie"
    e <- FineMappingRow(
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

test_that("lbf is empty when the fit carries no lbf matrix", {
    e <- FineMappingRow(
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

# ===========================================================================
# credibleSetSummary (merged from test_credibleSetSummary.R
# when R/credibleSetSummary.R was folded in)
# ===========================================================================

test_that("credibleSetSummary gives each CS's size, purity, V, logBF, lead", {
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
    e <- FineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)

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
    # cs_log_bf = true per-effect single-effect logBF = logSumExp(lbf[L, ]) -
    # log(P)
    expect_equal(s$cs_log_bf, c(1.9597, 1.2031), tolerance = 1e-3)
    # cs_pip = summed member PIP (inclusion mass captured)
    expect_equal(s$cs_pip, c(1.1, 0.8), tolerance = 1e-6)
    # cs_mean_effect = mean(posterior_mean / pip) over the members. It was NA
    # for every univariate susie fit while it read a `conditional_effect`
    # column that only multi-condition fits ever got -- an artefact of where
    # the column came from, not of the quantity.
    expect_true(all(is.finite(s$cs_mean_effect)))
    # lead = highest-PIP member
    expect_equal(s$lead_variant, c("chr1:100:A:G", "chr1:300:A:G"))
})

test_that("credibleSetSummary is empty when there are no credible sets", {
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
    e <- FineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)
    expect_equal(nrow(.fmrRowCredibleSetSummary(e)), 0L)
})

test_that("credibleSetSummary is empty without the requested coverage column", {
    # topLoci carries only cs_95 / cs_70 / cs_50; asking for coverage 0.99 finds
    # no cs_99 column -> the .csSummaryFit guard returns the empty summary.
    vn <- c("chr1:100:A:G", "chr1:200:C:T")
    tl <- data.frame(
        variant_id = vn,
        pip = c(0.7, 0.2),
        logBF = c(2.5, 0.3),
        cs_95 = c("susie_1", "susie_1"),
        cs_95_purity = c(0.9, 0.9),
        stringsAsFactors = FALSE
    )
    e <- FineMappingRow(variantIds = vn, susieFit = list(x = 1), topLoci = tl)
    expect_equal(nrow(.fmrRowCredibleSetSummary(e, coverage = 0.99)), 0L)
})

test_that("credibleSetSummary: NA cs_log_bf when the fit has no lbf matrix", {
    # A trimmed fit drops lbf_variable: the per-effect single-effect logBF
    # degrades to NA while the max-member cs_log10bf still comes from the
    # topLoci.
    vn <- c("chr1:100:A:G", "chr1:200:C:T")
    tl <- data.frame(
        variant_id = vn,
        pip = c(0.7, 0.2),
        logBF = c(2.5, 0.3),
        cs_95 = c("susie_1", "susie_1"),
        cs_95_purity = c(0.88, 0.88),
        stringsAsFactors = FALSE
    )
    e <- FineMappingRow(
        variantIds = vn,
        susieFit = list(trimmed = TRUE),
        topLoci = tl
    )
    s <- .fmrRowCredibleSetSummary(e)
    expect_equal(nrow(s), 1L)
    expect_true(is.na(s$cs_log_bf))
    expect_equal(s$cs_log10bf, 2.5)
})

test_that("credibleSetSummary: NA cs_log_bf for an all-non-finite lbf row", {
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
    e <- FineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)
    expect_true(is.na(.fmrRowCredibleSetSummary(e)$cs_log_bf))
})

# ===========================================================================
# fsusieCredibleBand / fsusieAffectedRegions (merged from
# test_fsusieAccessors.R when R/fsusieAccessors.R was folded in)
# ===========================================================================

test_that("fsusieCredibleBand gives a long table, lower <= effect <= upper", {
    skip_if_not_installed("fsusieR")
    skip_if_not_installed("wavethresh")
    cb <- .fmrRowFsusieCredibleBand(.fsa_entry(.fsa_fits$none))
    expect_named(
        cb,
        c("cs", "chrom", "pos", "effect", "lower", "upper", "halfwidth")
    )
    expect_gt(nrow(cb), 0)
    expect_true(all(
        cb$lower <= cb$effect + 1e-9 & cb$effect <= cb$upper + 1e-9
    ))
    expect_true(all(grepl("^fsusie_", cb$cs)))
})

test_that("trimmed fSuSiE keeps curves, bands and coloc inputs", {
    skip_if_not_installed("fsusieR")
    skip_if_not_installed("wavethresh")
    for (mode in c("none", "TI")) {
        fit <- .fsa_fits[[mode]]
        expectedFit <- pecotmr:::.fsusiePopulateCredibleBand(fit)
        tables <- list(list(sets = list(cs = list())))
        trimmed <- pecotmr:::trimFinemappingFit(
            fit,
            seq_along(fit$alpha),
            "fsusie",
            tables
        )
        entry <- unserialize(serialize(.fsa_entry(trimmed), NULL))
        saved <- entry@susieFit
        expect_null(saved$fitted_wc)
        expect_null(saved$fitted_wc2)
        expect_identical(saved$fitted_func, fit$fitted_func)
        expect_identical(saved$cred_band, expectedFit$cred_band)
        expect_equal(
            .fmrRowFsusieCredibleBand(entry),
            .fmrRowFsusieCredibleBand(.fsa_entry(fit))
        )
        expect_equal(
            .fmrRowFsusieAffectedRegions(entry),
            .fmrRowFsusieAffectedRegions(.fsa_entry(fit))
        )
        # The five positional filter arguments were folded into
        # ColocLbfFilterParam(); these are its defaults, pinned explicitly so
        # a later default change cannot quietly alter what this asserts.
        lbf <- pecotmr:::.colocExtractLbfFromEntry(
            entry,
            ColocLbfFilterParam(
                filterLbfCs = FALSE,
                secondary = NULL,
                concentration = 0.5,
                priorTol = 1e-9
            )
        )
        expected <- pecotmr:::.asLbfMatrix(fit)
        colnames(expected) <- variantIds(entry)
        expect_equal(lbf$lbf, expected)
        colocResult <- coloc::coloc.bf_bf(lbf$lbf, lbf$lbf)
        expect_true("PP.H4.abf" %in% names(colocResult$summary))
    }
})

test_that("fsusieAffectedRegions gives a GRanges of cs, purity, direction", {
    skip_if_not_installed("fsusieR")
    skip_if_not_installed("wavethresh")
    gr <- .fmrRowFsusieAffectedRegions(.fsa_entry(
        .fsa_fits$none,
        purityVal = 0.85
    ))
    expect_s4_class(gr, "GRanges")
    expect_gt(length(gr), 0)
    expect_true(all(
        c("cs", "purity", "direction") %in% names(S4Vectors::mcols(gr))
    ))
    expect_true(all(S4Vectors::mcols(gr)$direction %in% c("pos", "neg", NA)))
    # purity sourced from the entry's cs_95_purity (matched by CS variant
    # membership), NOT the (nonexistent) fit$purity slot
    expect_true(all(S4Vectors::mcols(gr)$purity == 0.85))
    expect_true(all(grepl("^fsusie_", S4Vectors::mcols(gr)$cs)))
})

test_that("fsusie accessors are empty for a fit without wavelet slots", {
    e <- FineMappingRow(
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

# ---- fsusieTraitEffects ---------------------------------------------

# The pipeline records trait_positions on the stored fit (jointEngine.R) and
# keeps the matching names in the row's jointTraits column; rebuild that pairing
# here, since .fsa_makeFit() returns the raw susiF fit.
.fsa_jointFit <- function(fit, pos) {
    fit$trait_positions <- pos
    fit
}

test_that("fsusieTraitEffects evaluates the curve at each joint trait", {
    skip_if_not_installed("fsusieR")
    skip_if_not_installed("wavethresh")
    fit <- .fsa_fits$none
    nTrait <- length(fit$outing_grid)
    # Positions BETWEEN grid points, so the interpolation is doing real work.
    pos <- as.numeric(seq_len(nTrait)) + 0.5
    names <- paste0("gene", seq_len(nTrait))
    res <- QtlFineMappingResult(
        studyName = "s",
        context = "brain",
        trait = "g",
        method = "fsusie",
        entry = list(.fsa_entry(.fsa_jointFit(fit, pos))),
        jointTraits = paste(names, collapse = ";")
    )

    te <- fsusieTraitEffects(res)
    expect_true(all(
        c(
            "study",
            "context",
            "trait",
            "method",
            "cs",
            "chrom",
            "jointTrait",
            "traitPos",
            "effect"
        ) %in%
            names(te)
    ))
    nCs <- length(fit$fitted_func)
    expect_equal(nrow(te), nCs * nTrait)
    expect_equal(unique(te$jointTrait), names)
    expect_equal(sort(unique(te$traitPos)), sort(pos))

    # The effect is the fitted curve interpolated onto the trait positions,
    # with constant extrapolation past the grid ends.
    first <- te[te$cs == "fsusie_1", ]
    expect_equal(
        first$effect,
        stats::approx(
            as.numeric(fit$outing_grid),
            as.numeric(fit$fitted_func[[1L]]),
            xout = pos,
            rule = 2,
            ties = "ordered"
        )$y
    )
})

test_that("fsusieTraitEffects refuses a fit whose positions, names disagree", {
    skip_if_not_installed("fsusieR")
    skip_if_not_installed("wavethresh")
    fit <- .fsa_fits$none
    res <- QtlFineMappingResult(
        studyName = "s",
        context = "brain",
        trait = "g",
        method = "fsusie",
        entry = list(.fsa_entry(.fsa_jointFit(
            fit,
            seq_along(fit$outing_grid)
        ))),
        jointTraits = "only;two"
    )
    expect_error(fsusieTraitEffects(res), "trait position")
})

test_that("fsusieTraitEffects returns zero rows for a non-fsusie fit", {
    res <- QtlFineMappingResult(
        studyName = "s",
        context = "brain",
        trait = "g",
        method = "susie",
        entry = list(FineMappingRow(
            variantIds = "chr1:100:A:G",
            susieFit = list(pip = 0.5),
            topLoci = data.frame(variant_id = "chr1:100:A:G", pip = 0.5)
        ))
    )
    te <- fsusieTraitEffects(res)
    expect_equal(nrow(te), 0L)
    expect_true(all(c("cs", "jointTrait", "traitPos", "effect") %in% names(te)))
})

test_that(".fsusieChrom returns NA when the fit has no named variants", {
    # unnamed
    expect_true(is.na(pecotmr:::.fsusieChrom(list(csd_X = c(1, 2, 3)))))
    expect_true(is.na(pecotmr:::.fsusieChrom(list()))) # no csd_X
})

test_that(".fsusieCsMapFromTopLoci uses index labels without CS columns", {
    fit <- list(
        cs = list(1L, 2L),
        csd_X = stats::setNames(1:3, paste0("chr1:", 1:3 * 100, ":A:G"))
    )
    m <- pecotmr:::.fsusieCsMapFromTopLoci(fit, topLoci = NULL)
    expect_equal(unname(m$label), c("fsusie_1", "fsusie_2"))
    expect_true(all(is.na(m$purity)))
    # a topLoci missing the required cs_95 / cs_95_purity columns takes the same
    # path
    m2 <- pecotmr:::.fsusieCsMapFromTopLoci(
        fit,
        topLoci = data.frame(variant_id = "x")
    )
    expect_equal(unname(m2$label), c("fsusie_1", "fsusie_2"))
})

# A minimal object that passes .isFsusieFit (the wavelet-slot gate) so the band
# / affected-region degradation paths can be exercised with the (untrimmed) band
# computation replaced by a mock returning a controlled, degenerate fit.
.fsa_fakeFit <- function() {
    list(fitted_wc2 = 1, fitted_func = list(1), outing_grid = 1:2, alpha = 1)
}
.fsa_bareEntry <- function() {
    FineMappingRow(
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

test_that("fsusieAffectedRegions is empty when affected_reg finds nothing", {
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

test_that("fsusieAffectedRegions: NA direction for a region off the grid", {
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
    # Region [500, 600] lies outside the grid [100, 200]: no grid points fall in
    # it, so the per-region effect-direction summary is NA.
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

test_that("credibleSets warns when the coverage has no purity column", {
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
    e <- FineMappingRow(variantIds = vn, susieFit = fit, topLoci = tl)

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
        FineMappingRow(
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
        FineMappingRow(
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
# Fold weights: the row holds their variant ids once
# ===========================================================================

.fmr_cv <- function(ids, folds = 2L) {
    fw <- matrix(
        seq_len(length(ids) * folds) / 10,
        length(ids),
        folds,
        dimnames = list(ids, str_c("fold_", seq_len(folds)))
    )
    list(
        samplePartition = data.frame(Sample = "s1", Fold = 1L),
        foldWeights = fw
    )
}

.fmr_tl <- function(ids) {
    data.frame(variant_id = ids, pip = rep(0.1, length(ids)))
}

test_that("a row stores its own fold weights unnamed; cvResult() names them", {
    vids <- .fmr_vids(3L)
    cv <- .fmr_cv(vids)
    row <- FineMappingRow(vids, NULL, .fmr_tl(vids), cvResult = cv)
    expect_null(rownames(pecotmr:::.fmRowCvStored(row)$foldWeights))
    expect_identical(cvResult(row), cv)
    # Weights on other variants keep their ids: a cross-validation-only row
    # has no variants to name them by, and a row in another order cannot
    # name them positionally.
    cvOnly <- pecotmr:::.fmCvOnlyRow(cv)
    expect_identical(pecotmr:::.fmRowCvStored(cvOnly), cv)
    expect_identical(cvResult(cvOnly), cv)
    swapped <- FineMappingRow(rev(vids), NULL, .fmr_tl(rev(vids)), cv)
    expect_identical(cvResult(swapped), cv)
    # A payload of any other shape passes through.
    other <- FineMappingRow(vids, NULL, .fmr_tl(vids), list(rsq = 0.5))
    expect_identical(cvResult(other), list(rsq = 0.5))
})

test_that("a merged row names each region's fold weights by its variants", {
    a <- .fmr_vids(3L)[1:2]
    b <- .fmr_vids(3L)[3L]
    rows <- list(
        FineMappingRow(a, NULL, .fmr_tl(a), cvResult = .fmr_cv(a)),
        FineMappingRow(b, NULL, .fmr_tl(b), cvResult = .fmr_cv(b))
    )
    merged <- pecotmr:::.fmMergeEntries(rows)
    stored <- pecotmr:::.fmRowCvStored(merged)
    expect_named(stored, c("region1", "region2"))
    expect_null(rownames(stored$region1$foldWeights))
    expect_null(rownames(stored$region2$foldWeights))
    expect_identical(
        cvResult(merged),
        list(region1 = .fmr_cv(a), region2 = .fmr_cv(b))
    )
    # A region without CV leaves the others unable to say whose variants
    # they cover, so they keep their ids.
    rows[[2L]] <- FineMappingRow(b, NULL, .fmr_tl(b))
    partial <- pecotmr:::.fmMergeEntries(rows)
    expect_identical(
        pecotmr:::.fmRowCvStored(partial),
        list(region1 = .fmr_cv(a), region2 = NULL)
    )
    expect_identical(cvResult(partial), pecotmr:::.fmRowCvStored(partial))
})

# ===========================================================================
# topLoci in a collection: stored compact, handed out plain
# ===========================================================================

# A topLoci table in the order the fine-mapping post-processor writes it.
.fmr_fullTl <- function(ids, p = NULL) {
    z <- seq(-2, 2, length.out = length(ids))
    tibble::tibble(
        variant_id = ids,
        N = rep(100, length(ids)),
        marginal_z = z,
        marginal_z_supplied = rep(FALSE, length(ids)),
        marginal_p = p %||% pecotmr:::.fmDerivedP(z),
        pip = seq(0.1, 0.5, length.out = length(ids)),
        cs_95 = c("susie_1", rep("susie_0", length(ids) - 1L)),
        method = rep("susie", length(ids))
    )
}

.fmr_tlRes <- function(tl) {
    ids <- tl$variant_id
    QtlFineMappingResult(
        studyName = "S",
        context = "C",
        trait = "T",
        method = "susie",
        entry = list(FineMappingRow(ids, list(pip = tl$pip), tl))
    )
}

test_that("a collection stores topLoci compactly and hands it back plain", {
    tl <- .fmr_fullTl(.fmr_vids(6L))
    res <- .fmr_tlRes(tl)
    stored <- mcols(unlist(res, use.names = FALSE))
    # marginal_p is the p marginal_z gives, and the pip the fit's, so both
    # are left out ...
    expect_false(any(is_in(c("marginal_p", "pip"), colnames(stored))))
    # ... and the repeating columns are run-length encoded, the rest plain.
    expect_s4_class(stored$N, "Rle")
    expect_s4_class(stored$cs_95, "Rle")
    expect_s4_class(stored$method, "Rle")
    expect_true(is.numeric(stored$marginal_z))
    # The row handed out holds them plain, marginal_p back in its place.
    back <- pecotmr:::.fmrPartsTopLoci(pecotmr:::.collectionEntry(res, 1L))
    expect_identical(back[names(tl)], tl)
    expect_identical(
        names(back)[1:8],
        c(
            "variant_id",
            "chrom",
            "pos",
            "A1",
            "A2",
            "N",
            "marginal_z",
            "marginal_z_supplied"
        )
    )
    expect_equal(match("marginal_p", names(back)), 9L)
})

test_that("a supplied marginal_p that differs from the derived one is kept", {
    ids <- .fmr_vids(4L)
    p <- c(0.01, 0.2, 0.5, 0.9)
    res <- .fmr_tlRes(.fmr_fullTl(ids, p = p))
    stored <- mcols(unlist(res, use.names = FALSE))
    expect_equal(as.vector(stored$marginal_p), p)
    back <- pecotmr:::.fmrPartsTopLoci(pecotmr:::.collectionEntry(res, 1L))
    expect_equal(back$marginal_p, p)
})

test_that("marginal_p out of its place, or without marginal_z, is kept", {
    md <- S4Vectors::DataFrame(
        marginal_p = 0.5,
        marginal_z = 0,
        marginal_z_supplied = FALSE
    )
    expect_identical(pecotmr:::.fmDropDerivedP(md), md)
    noZ <- S4Vectors::DataFrame(marginal_p = 0.5)
    expect_identical(pecotmr:::.fmDropDerivedP(noZ), noZ)
    expect_true(is.na(pecotmr:::.fmDerivedPPosition("pip")))
    # Without marginal_z_supplied, marginal_p sits right after marginal_z.
    expect_equal(
        pecotmr:::.fmDerivedPPosition(c("N", "marginal_z", "pip")),
        3L
    )
})

test_that("updateObject() stores a saved collection's topLoci compactly", {
    res <- .fmr_tlRes(.fmr_fullTl(.fmr_vids(6L)))
    # As saved before collections stored topLoci compactly.
    flat <- unlist(res, use.names = FALSE)
    plain <- pecotmr:::.fmPlainElement(flat)
    saved <- res
    attr(saved, "unlistData") <- plain
    expect_false(methods::is(mcols(unlist(saved))$N, "Rle"))
    up <- updateObject(saved)
    expect_identical(up, res)
    expect_identical(updateObject(up), up)
})

# ===========================================================================
# adjustPips: fits this code cannot slice safely
# ===========================================================================

.fmr_susieFit <- function(p = 5L, maxNumSingleEffects = 2L, seed = 7L) {
    set.seed(seed)
    lbf <- matrix(
        rnorm(maxNumSingleEffects * p, sd = 3),
        maxNumSingleEffects,
        p
    )
    alpha <- lbfToAlpha(lbf)
    list(
        alpha = alpha,
        mu = matrix(0.3, maxNumSingleEffects, p),
        mu2 = matrix(1.2, maxNumSingleEffects, p),
        lbf_variable = lbf,
        V = rep(1, maxNumSingleEffects),
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
            c(1, 1, 1, 0, 0, 0, 2, 0, 0, 0),
            2,
            5,
            byrow = TRUE,
            dimnames = list(c("L1", "L2"), vn)
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
        variant_id = vn,
        pip = fit$pip,
        logBF = c(1, 1, 1, 2),
        cs_95 = c("susie_2", "susie_2", "susie_2", "susie_1"),
        stringsAsFactors = FALSE
    )
    s <- .csSummaryFit(tl, fit, 0.95)
    big <- s[s$n_variants == 3L, , drop = FALSE]
    expect_equal(big$effect_id, "L2")
    expect_equal(big$cs, "susie_2")
})

test_that("a trimmed fit's summary reads each set's own effect", {
    # Trimming keeps the non-null effects, so effects 1 and 3 of three sit on
    # rows 1 and 2. Looking V and the log BF up by effect index read past the
    # end for L3, and would have read L3's row for an L2.
    vn <- str_c("chr1:", (1:4) * 100, ":A:G")
    fit <- list(
        pip = c(0.6, 0.3, 0.05, 0.05),
        alpha = matrix(
            c(0.7, 0.1, 0.1, 0.1, 0.1, 0.7, 0.1, 0.1),
            2L,
            byrow = TRUE,
            dimnames = list(NULL, vn)
        ),
        lbf_variable = matrix(c(3, 0, 0, 0, 0, 2, 0, 0), 2L, byrow = TRUE),
        V = c(0.4, 0.2),
        effect_index = c(1L, 3L),
        sets = list(
            cs = list(L1 = 1L, L3 = 2L),
            cs_index = c(1L, 3L),
            requested_coverage = 0.95
        )
    )
    class(fit) <- "susie"
    tl <- data.frame(
        variant_id = vn,
        pip = fit$pip,
        logBF = c(3, 2, 0, 0),
        cs_95 = c("susie_1", "susie_3", "susie_0", "susie_0")
    )
    s <- .csSummaryFit(tl, fit, 0.95)
    s <- s[order(s$effect_id), , drop = FALSE]
    expect_equal(s$effect_id, c("L1", "L3"))
    expect_equal(s$V, c(0.4, 0.2))
    expect_equal(
        s$cs_log_bf,
        c(log(mean(exp(c(3, 0, 0, 0)))), log(mean(exp(c(0, 2, 0, 0)))))
    )
})

test_that("adjustPips names a trimmed fit's rebuilt sets by effect", {
    # The rebuilt sets took their L<k> names and cs_index from varAlpha's
    # rows, so on a trimmed fit holding effects 1 and 3 an L3 came back as
    # L2 -- relabelling topLoci's susie_3 as susie_2.
    varAlpha <- matrix(c(0.9, 0.1, 0.1, 0.9), 2L, byrow = TRUE)
    sets <- .adjustPipsCsAtCoverage(varAlpha, c(0.4, 0.2), 0.8, c(1L, 3L))
    expect_equal(names(sets$cs), c("L1", "L3"))
    expect_equal(sets$cs_index, c(1L, 3L))
    fit <- list(V = c(0.4, 0.2), effect_index = c(1L, 3L))
    expect_equal(
        .adjustPipsCsLabels(varAlpha, fit$V, "cs_80", "susie", 2L, c(1L, 3L)),
        c("susie_1", "susie_3")
    )
})

test_that("a fit with no stored sets falls back to the per-variant column", {
    # Minimal / hand-built fits carry no sets$cs; the column is then the only
    # source of membership there is.
    vn <- str_c("chr1:", (1:3) * 100, ":A:G")
    fit <- list(pip = c(0.5, 0.4, 0.1))
    class(fit) <- "susie"
    tl <- data.frame(
        variant_id = vn,
        pip = fit$pip,
        logBF = c(1, 1, 0),
        cs_95 = c("susie_1", "susie_1", "susie_0"),
        stringsAsFactors = FALSE
    )
    s <- .csSummaryFit(tl, fit, 0.95)
    expect_equal(nrow(s), 1L)
    expect_equal(s$n_variants, 2L)
    expect_equal(s$cs, "susie_1")
})

# ===========================================================================
# Purity filtering, coverage defaults and column projection
# ===========================================================================

test_that("purity filtering keeps everything when there is nothing to judge", {
    # No purity columns, or purity columns with no matching cs column: the
    # entry cannot be assessed, so nothing is dropped rather than everything.
    tl <- tibble(variant_id = c("a", "b"), pip = c(0.1, 0.9))
    expect_equal(pecotmr:::.fmePurityKeep(tl, 0.8), c(TRUE, TRUE))
    tl2 <- tibble(
        variant_id = c("a", "b"),
        cs_0.95_purity = c(0.9, 0.2)
    )
    expect_equal(pecotmr:::.fmePurityKeep(tl2, 0.8), c(TRUE, TRUE))
})

test_that("purity filtering drops only in-CS variants below the threshold", {
    # `susie_0` means "in no credible set", so that row survives regardless
    # of the purity recorded beside it.
    tl <- tibble(
        variant_id = c("a", "b"),
        cs_0.95 = c("susie_1", "susie_0"),
        cs_0.95_purity = c(0.2, 0.9)
    )
    expect_equal(pecotmr:::.fmePurityKeep(tl, 0.8), c(FALSE, TRUE))
})

test_that("an unusable coverage falls back to 0.95", {
    expect_equal(pecotmr:::.adjustPipsCoverage(NULL), 0.95)
})

test_that("projecting an absent column yields typed NAs of the right length", {
    # One NA per row, of the requested type, so the column binds against the
    # rest of the table instead of collapsing it.
    tl <- tibble(variant_id = c("a", "b"), pip = c(0.1, 0.9))
    expect_equal(
        pecotmr:::.tlCol(tl, "nope", "character"),
        rep(NA_character_, 2L)
    )
    expect_equal(pecotmr:::.tlCol(tl, "nope", "numeric"), rep(NA_real_, 2L))
    expect_equal(pecotmr:::.tlCol(tl, "nope", "integer"), rep(NA_integer_, 2L))
    expect_equal(pecotmr:::.tlCol(tl, "pip", "numeric"), c(0.1, 0.9))
})

# ===========================================================================
# Credible-set summary assembly and PIP adjustment
# ===========================================================================

test_that("an absent or empty topLoci summarises to the empty CS table", {
    expect_equal(nrow(pecotmr:::.csSummaryFit(NULL, list(), 0.95)), 0L)
    expect_equal(
        nrow(pecotmr:::.csSummaryFit(tibble(a = integer(0)), list(), 0.95)),
        0L
    )
})

test_that("credible sets are resolved from the secondary store by coverage", {
    # A trimmed fit keeps non-primary coverages in sets_secondary, keyed
    # cs_<pct>_; asking for a coverage it does not hold is NULL, not the
    # primary sets, which would silently answer the wrong question.
    sec <- list(cs_70_x = list(sets = list(cs = list(L1 = 1:2))))
    fit <- list(sets = list(), sets_secondary = sec)
    expect_equal(
        pecotmr:::.csSetsForCoverage(fit, 0.70)$cs$L1,
        1:2
    )
    expect_null(pecotmr:::.csSetsForCoverage(fit, 0.95))
})

test_that("an untrimmed fit with no recorded coverage uses its primary sets", {
    # No requested_coverage and no secondary store: the primary sets are the
    # only ones there are, so they answer for any coverage.
    fit <- list(sets = list(cs = list(L1 = 1:2)))
    expect_equal(pecotmr:::.csSetsForCoverage(fit, 0.95)$cs$L1, 1:2)
})

test_that("mean conditional effect is the posterior mean over the PIP", {
    # NA when either column the quotient needs is absent.
    expect_true(is.na(pecotmr:::.csMeanEffect(tibble(a = 1))))
    expect_true(is.na(pecotmr:::.csMeanEffect(tibble(pip = c(0.5, 0.5)))))
    expect_equal(
        pecotmr:::.csMeanEffect(
            tibble(posterior_mean = c(0.5, 1.5), pip = c(0.5, 0.5))
        ),
        2
    )
    # A zero PIP has no conditional effect to speak of, so it drops out
    # rather than contributing an infinity.
    expect_equal(
        pecotmr:::.csMeanEffect(
            tibble(posterior_mean = c(0.5, 1), pip = c(0.5, 0))
        ),
        1
    )
})

test_that("rebuilding an empty topLoci returns it unchanged", {
    tl <- tibble(variant_id = character(0))
    expect_equal(
        nrow(pecotmr:::.adjustPipsRebuildTopLoci(tl, list(), character(0))),
        0L
    )
})

test_that("posterior columns are left in place when the fit has no alpha", {
    # Without alpha there is nothing to recompute from, so the existing
    # values stay rather than being overwritten with NA.
    tl <- tibble(variant_id = "a")
    expect_equal(pecotmr:::.adjustPipsPosterior(tl, list())$variant_id, "a")
})

# ---------------------------------------------------------------------------
# adjustPips index resolution and the CS-label fallbacks: each decides how a
# fit slot lines up with the entry's variants, or what to call a set when the
# table does not say.
# ---------------------------------------------------------------------------

test_that(".csMethodTag falls back to the fit class when no label survives", {
    fit <- structure(list(), class = "susieAsh")
    expect_equal(pecotmr:::.csMethodTag(NULL, fit, "cs_95"), "susie_ash")
    # "_0" means "not in a credible set", so it is filtered out and cannot
    # supply the tag either.
    onlyNull <- data.frame(
        variant_id = "v1",
        cs_95 = "L1_0",
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.csMethodTag(onlyNull, fit, "cs_95"), "susie_ash")
    labelled <- data.frame(
        variant_id = "v1",
        cs_95 = "susie_1",
        stringsAsFactors = FALSE
    )
    expect_equal(pecotmr:::.csMethodTag(labelled, fit, "cs_95"), "susie")
})

test_that(".csSummaryRow returns NULL when the set has no rows in the table", {
    ctx <- list(
        tl = tibble::tibble(variant_id = c("v1", "v2"), pip = c(0.5, 0.5))
    )
    expect_null(pecotmr:::.csSummaryRow(list(ids = "zzz", eff = 1L), ctx))
})

test_that(".adjustPipsIdxFor accepts the null-weight extra column", {
    p <- list(nVariants = 3L, keepIdx = c(1L, 2L, 3L), cols = c(1L, 2L, 3L, 4L))
    expect_equal(pecotmr:::.adjustPipsIdxFor(3L, p, "alpha"), c(1L, 2L, 3L))
    # A null_weight fit carries one column MORE than the entry's variants.
    expect_equal(
        pecotmr:::.adjustPipsIdxFor(4L, p, "alpha"),
        c(1L, 2L, 3L, 4L)
    )
    expect_error(
        pecotmr:::.adjustPipsIdxFor(9L, p, "alpha"),
        "fit slot `alpha` spans 9 variants but the entry has 3"
    )
})

test_that(".adjustPipsCols treats a dimensionless slot as a vector", {
    p <- list(nVariants = 3L, keepIdx = c(1L, 3L), cols = c(1L, 2L, 3L))
    expect_equal(pecotmr:::.adjustPipsCols(c(10, 20, 30), p, "pip"), c(10, 30))
    expect_null(pecotmr:::.adjustPipsCols(NULL, p, "pip"))
})

test_that(".adjustPipsRebuildSecondary passes through an entry with no sets", {
    entry <- list(a = 1)
    expect_identical(
        pecotmr:::.adjustPipsRebuildSecondary(entry, NULL, NULL),
        entry
    )
})

test_that(".adjustPipsMethodLabel defaults to 'cs'", {
    expect_equal(pecotmr:::.adjustPipsMethodLabel(data.frame(x = 1)), "cs")
    # A method column of all NA names nothing either.
    expect_equal(
        pecotmr:::.adjustPipsMethodLabel(
            data.frame(method = NA_character_, stringsAsFactors = FALSE)
        ),
        "cs"
    )
    expect_equal(
        pecotmr:::.adjustPipsMethodLabel(
            data.frame(method = "susie", stringsAsFactors = FALSE)
        ),
        "susie"
    )
})

test_that(".csLog10Bf and .csMeanEffect return NA, never -Inf or NaN", {
    # Both previously disagreed with their own missing-column branch: an empty
    # or all-NA column yielded -Inf / NaN from max() / mean() while a missing
    # column yielded NA_real_. These values reach output tables.
    expect_identical(pecotmr:::.csLog10Bf(list()), NA_real_)
    expect_identical(
        pecotmr:::.csLog10Bf(list(logBF = numeric(0))),
        NA_real_
    )
    expect_identical(
        pecotmr:::.csLog10Bf(list(logBF = c(NA_real_, NA_real_))),
        NA_real_
    )
    expect_identical(pecotmr:::.csLog10Bf(list(logBF = c(1, 3, 2))), 3)

    expect_identical(pecotmr:::.csMeanEffect(list()), NA_real_)
    expect_identical(
        pecotmr:::.csMeanEffect(
            list(posterior_mean = numeric(0), pip = numeric(0))
        ),
        NA_real_
    )
    expect_identical(
        pecotmr:::.csMeanEffect(
            list(posterior_mean = c(NA, NA), pip = c(0.5, 0.5))
        ),
        NA_real_
    )
    expect_identical(
        pecotmr:::.csMeanEffect(
            list(posterior_mean = c(0.5, 1.5), pip = c(0.5, 0.5))
        ),
        2
    )
})

test_that("fsusieCredibleBand() carries the band half-width", {
    data(fsusieFineMappingExample)
    band <- as.data.frame(fsusieCredibleBand(fsusieFineMappingExample))
    expect_true(is_in("halfwidth", names(band)))
    expect_gt(nrow(band), 0L)
    expect_equal(band$halfwidth, (band$upper - band$lower) / 2)
    # A band half-width is a magnitude, never negative.
    expect_true(all(band$halfwidth >= 0))
})

test_that("an empty credible band still declares halfwidth", {
    expect_true(is_in("halfwidth", names(.emptyCredibleBand())))
    expect_identical(nrow(.emptyCredibleBand()), 0L)
})

# ---------------------------------------------------------------------------
# fsusie trait-effect edge paths, and the two projections either side of them
# ---------------------------------------------------------------------------

# The minimum a fit needs to pass `.isFsusieFit()` and carry trait positions:
# one curve per effect over `outing_grid`, plus the two fields that identify a
# susiF fit. A `NULL` element of `fitted_func` is what an effect that produced
# no curve looks like.
# @noRd
.fte_stubFit <- function(curves, pos) {
    list(
        fitted_func = curves,
        outing_grid = c(1, 2, 3),
        fitted_wc2 = 1,
        alpha = 1,
        trait_positions = pos
    )
}

# @noRd
.fte_row <- function(fit) {
    FineMappingRow(
        variantIds = "chr1:100:A:G",
        susieFit = fit,
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
}

test_that("an effect with no fitted curve contributes no trait effects", {
    res <- QtlFineMappingResult(
        studyName = "s",
        context = "brain",
        trait = "g",
        method = "fsusie",
        entry = list(.fte_row(.fte_stubFit(list(NULL, NULL), 1.5))),
        jointTraits = "g1"
    )
    te <- fsusieTraitEffects(res)
    # Every effect skipped, so there is nothing to report -- but the empty
    # table keeps the schema, so a caller can bind or filter it exactly as it
    # would a populated one.
    expect_equal(nrow(te), 0L)
    expect_true(all(
        c("cs", "chrom", "jointTrait", "traitPos", "effect") %in% names(te)
    ))
})

test_that("a fit with no jointTraits names its curves after the row", {
    res <- QtlFineMappingResult(
        studyName = "s",
        context = "brain",
        trait = "g",
        method = "fsusie",
        entry = list(.fte_row(.fte_stubFit(list(c(1, 2, 3)), 1.5)))
    )
    te <- fsusieTraitEffects(res)
    # The joint member names ride on `jointTraits`; a fit over a single trait
    # has none, so the row's own trait is the name of its curve.
    expect_equal(nrow(te), 1L)
    expect_equal(te$jointTrait, "g")
    # The curve is linear on the grid, so the midpoint interpolates exactly.
    expect_equal(te$effect, 1.5)
})

test_that("a GWAS fsusie fit has no trait axis to name its curves", {
    res <- GwasFineMappingResult(
        studyName = "s",
        method = "fsusie",
        entry = list(.fte_row(.fte_stubFit(list(c(1, 2, 3)), 1.5)))
    )
    # A GWAS result is not keyed by trait at all, so there is no name to fall
    # back to and the fit's one trait position has nothing to pair with.
    expect_error(
        fsusieTraitEffects(res),
        "1 trait position\\(s\\) but 0 trait name\\(s\\)"
    )
})

test_that("a credible-set label matching no fit effect claims none", {
    tl <- data.frame(
        variant_id = "chr1:100:A:G",
        cs_95 = "fsusie_1",
        stringsAsFactors = FALSE
    )
    # The label's variants are not any effect's membership, so the label has
    # no effect to attach to -- NA rather than silently taking the first.
    expect_true(is.na(.fsusieFirstMatchingEffect(
        "fsusie_1",
        fitSets = list("chr1:200:A:G", "chr1:300:A:G"),
        topLoci = tl
    )))
    expect_equal(
        .fsusieFirstMatchingEffect(
            "fsusie_1",
            fitSets = list("chr1:300:A:G", "chr1:100:A:G"),
            topLoci = tl
        ),
        2L
    )
})

test_that("a NULL posterior table projects to an empty GRanges", {
    gr <- .variantTableGRanges(NULL)
    expect_s4_class(gr, "GRanges")
    expect_length(gr, 0L)
})

test_that("a fit whose prior variances all vanished labels nothing", {
    vids <- str_c("chr1:", (1:4) * 100, ":A:G")
    p <- length(vids)
    L <- 2L
    set.seed(3L)
    lbf <- matrix(stats::rnorm(L * p), nrow = L, ncol = p)
    colnames(lbf) <- vids
    alpha <- lbfToAlpha(lbf)
    pip <- as.numeric(1 - apply(1 - alpha, 2, prod))
    entry <- FineMappingRow(
        variantIds = vids,
        susieFit = list(
            pip = pip,
            alpha = alpha,
            lbf_variable = lbf,
            # Every prior variance driven to zero: no effect clears the
            # tolerance, so there is no credible set left to label with.
            V = rep(0, L)
        ),
        topLoci = data.frame(
            variant_id = vids,
            pip = pip,
            method = "susie",
            cs_95 = c("susie_1", rep("susie_0", p - 1L)),
            stringsAsFactors = FALSE
        )
    )
    adj <- adjustPips(entry, vids[1:3])
    expect_equal(.fmrRowTopLoci(adj, raw = TRUE)$cs_95, rep("susie_0", 3L))
})

# ===========================================================================
# A row stores no per-variant pip column when its fit reproduces it: every
# view restores it from the fit. These tests rebuild each bundled example --
# saved with the column -- through FineMappingRow() and check that nothing a
# caller can see changes.
# ===========================================================================

.fpd_hasStoredPip <- function(x) {
    "pip" %in% .fmStoredColumns(x)
}

.fpd_csColumns <- c(
    "cs_95",
    "cs_70",
    "cs_50",
    "cs_95_purity",
    "cs_70_purity",
    "cs_50_purity"
)

test_that("a rebuilt row stores no pip its fit reproduces", {
    x <- .fmLoadExample("qtlFineMappingExample")
    expect_true(.fpd_hasStoredPip(x))
    expect_false(.fpd_hasStoredPip(.fmRebuildExample(x)))
})

test_that("a rebuilt row stores no credible-set column its fit reproduces", {
    x <- .fmLoadExample("qtlFineMappingExample")
    expect_true(all(.fpd_csColumns %in% .fmStoredColumns(x)))
    new <- .fmRebuildExample(x)
    expect_false(any(.fpd_csColumns %in% .fmStoredColumns(new)))
    # The row as GRanges derives them back, in buildTopLoci()'s order: the
    # order the example's own row hands them out in.
    restored <- as(.fmExampleRow(x, 1L), "GRanges")
    held <- as(pecotmr:::.fmrRowParts(x, 1L), "GRanges")
    expect_identical(
        colnames(S4Vectors::mcols(restored)),
        colnames(S4Vectors::mcols(held))
    )
})

test_that("credible-set columns the fit cannot give back are stored", {
    inputs <- .fmExampleRowInputs(.fmLoadExample("qtlFineMappingExample"), 1L)
    stored <- function(inputs) {
        row <- do.call(FineMappingRow, inputs)
        q <- QtlFineMappingResult("s", "c", "t", "susie", list(row))
        intersect(.fpd_csColumns, .fmStoredColumns(q))
    }
    # Without the secondary credible sets the 0.7 and 0.5 columns have no
    # source, and the block is kept whole.
    noSecondary <- inputs
    noSecondary$susieFit$sets_secondary <- NULL
    expect_identical(stored(noSecondary), .fpd_csColumns)
    # Nor is a relabelled column re-derived over.
    relabelled <- inputs
    relabelled$topLoci$cs_70 <- sub("_[0-9]+$", "_9", relabelled$topLoci$cs_70)
    expect_identical(stored(relabelled), .fpd_csColumns)
    row <- do.call(FineMappingRow, relabelled)
    expect_identical(
        S4Vectors::mcols(as(row, "GRanges"))$cs_70,
        relabelled$topLoci$cs_70
    )
    # A table with only some of the block is not laid out as the fit writes
    # it.
    partial <- inputs
    partial$topLoci$cs_50_purity <- NULL
    expect_identical(stored(partial), setdiff(.fpd_csColumns, "cs_50_purity"))
})

test_that("every view of a rebuilt example is unchanged", {
    for (name in .fmExampleNames()) {
        old <- .fmLoadExample(name)
        new <- .fmRebuildExample(old)
        same <- function(f) {
            expect_equal(f(new), f(old), info = name)
        }
        same(function(x) topLoci(x))
        same(function(x) topLoci(x, signalCutoff = 0))
        same(function(x) topLoci(x, signalCutoff = 0.1, minPurity = 0.5))
        for (coverage in c(0.95, 0.7, 0.5)) {
            same(function(x) credibleSets(x, coverage = coverage))
        }
        same(function(x) credibleSetSummary(x))
        same(function(x) marginalEffects(x))
        same(flattenTupleRanges)
        same(function(x) as(x, "GRanges"))
    }
})

test_that("pip() of a rebuilt row comes from its fit", {
    old <- .fmLoadExample("gwasFineMappingExample")
    new <- .fmRebuildExample(old)
    expect_equal(pip(new), pip(old))
    expect_false(anyNA(pip(new)))
})

test_that("the BED writer writes the same table for a rebuilt example", {
    for (name in c("qtlFineMappingExample", "gwasFineMappingExample")) {
        old <- .fmLoadExample(name)
        paths <- c(withr::local_tempfile(), withr::local_tempfile())
        writeBed(old, paths[[1L]], signalCutoff = 0)
        writeBed(.fmRebuildExample(old), paths[[2L]], signalCutoff = 0)
        expect_identical(
            readLines(paths[[2L]]),
            readLines(paths[[1L]]),
            info = name
        )
    }
})

test_that("a flatten / nest round trip keeps the pip column dropped", {
    new <- .fmRebuildExample(.fmLoadExample("qtlFineMappingExample"))
    flat <- flattenTupleRanges(new)
    expect_true("pip" %in% colnames(S4Vectors::mcols(flat)))
    expect_true(all(.fpd_csColumns %in% colnames(S4Vectors::mcols(flat))))
    back <- nestTupleRanges(flat, new)
    expect_false(.fpd_hasStoredPip(back))
    expect_false(any(.fpd_csColumns %in% .fmStoredColumns(back)))
    expect_equal(
        topLoci(back, signalCutoff = 0),
        topLoci(new, signalCutoff = 0)
    )
    # A verb that changes the PIPs keeps its own values.
    halved <- flat
    S4Vectors::mcols(halved)$pip <- S4Vectors::mcols(halved)$pip / 2
    edited <- nestTupleRanges(halved, new)
    expect_true(.fpd_hasStoredPip(edited))
    expect_equal(
        S4Vectors::mcols(topLoci(edited, signalCutoff = 0))$pip,
        S4Vectors::mcols(topLoci(new, signalCutoff = 0))$pip / 2
    )
})

test_that("a pip its fit does not give is stored", {
    inputs <- .fmExampleRowInputs(.fmLoadExample("qtlFineMappingExample"), 1L)
    inputs$topLoci$pip <- inputs$topLoci$pip / 2
    edited <- do.call(FineMappingRow, inputs)
    expect_equal(
        S4Vectors::mcols(as(edited, "GRanges"))$pip,
        inputs$topLoci$pip
    )
    # No fit at all: the table's pip is the only one there is.
    bare <- FineMappingRow(
        variantIds = c("chr1:100:A:G", "chr1:200:C:T"),
        susieFit = list(),
        topLoci = data.frame(
            variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
            pip = c(0.9, 0.1)
        )
    )
    expect_equal(S4Vectors::mcols(as(bare, "GRanges"))$pip, c(0.9, 0.1))
})

test_that("a collection joining saved and rebuilt rows reads every pip", {
    old <- .fmLoadExample("qtlFineMappingLbfExample")
    new <- .fmRebuildExample(old)
    mixed <- c(old[1L], new[2L])
    expect_equal(
        topLoci(mixed, signalCutoff = 0),
        topLoci(old, signalCutoff = 0)
    )
})

# ===========================================================================
# Restoring the columns a row leaves to its fit: where nothing is restored
# ===========================================================================

.fmr_ids <- c("chr1:100:A:G", "chr1:200:C:T")

# A row's variants as stored, from a fit that gives none of its columns
# back, so the table's pip and credible-set columns all stay.
.fmr_storedVariants <- function() {
    variants(FineMappingRow(
        .fmr_ids,
        list(),
        data.frame(
            variant_id = .fmr_ids,
            pip = c(0.4, 0.6),
            cs_95 = c("susie_1", "susie_0"),
            cs_95_purity = c(0.9, 0)
        )
    ))
}

test_that("the pip and credible-set slots follow their anchor columns", {
    expect_identical(pecotmr:::.fmRowPipSlot(c("variant_id", "beta")), 0L)
    expect_identical(pecotmr:::.fmRowPipSlot(c("A2", "marginal_p", "x")), 2L)
    expect_identical(pecotmr:::.fmRowCsSlot(c("A2", "logBF", "x")), 2L)
    expect_identical(pecotmr:::.fmRowCsSlot(c("A2", "pip", "x")), 2L)
    expect_identical(pecotmr:::.fmRowCsSlot(c("A2", "x")), 1L)
})

test_that("credible-set labels give a method tag only when they share one", {
    tag <- pecotmr:::.fmRowCsLabelTag
    expect_identical(
        tag(S4Vectors::DataFrame(cs_95 = c("susie_1", "susie_0"))),
        "susie"
    )
    expect_null(tag(S4Vectors::DataFrame(cs_95 = c("susie_1", "fsusie_1"))))
    expect_null(tag(S4Vectors::DataFrame(cs_95 = c("_1", "_2"))))
})

test_that("a fit gives no credible-set columns it cannot place or label", {
    fit <- list(pip = c(0.5, 0.5))
    cols <- pecotmr:::.fitCsColumns
    expect_null(cols(fit, "cs_95", NULL, 1:2, 2L))
    expect_null(cols(fit, "cs_95", "susie", NULL, 2L))
    expect_null(cols(fit, "cs_95", "susie", c(1L, NA), 2L))
    # Nor without the credible sets of the coverage asked for.
    expect_null(cols(fit, "cs_95", "susie", 1:2, 2L))
})

test_that("restoring leaves what a row stores, or its fit cannot give", {
    gr <- .fmr_storedVariants()
    row <- methods::new(
        "FineMappingRow",
        variants = gr,
        susieFit = list(pip = c(0.4, 0.6)),
        derivedColumns = c("cs_95", "method=susie")
    )
    # Every credible-set column it records is stored.
    expect_identical(pecotmr:::.fmRowRestoreFitCs(gr, row), gr)
    # Without them, a fit with no credible sets gives nothing back.
    bare <- gr
    S4Vectors::mcols(bare) <- S4Vectors::mcols(gr)[, c("A1", "A2", "pip")]
    expect_identical(pecotmr:::.fmRowRestoreFitCs(bare, row), bare)
    # A fit with no PIPs gives no pip back.
    noPip <- bare
    S4Vectors::mcols(noPip) <- S4Vectors::mcols(gr)[, c("A1", "A2")]
    fitless <- methods::new(
        "FineMappingRow",
        variants = noPip,
        susieFit = list(alpha = matrix(0.5, 1L, 2L))
    )
    expect_identical(pecotmr:::.fmRowRestoreFitPip(noPip, fitless), noPip)
})

test_that(".fmRowFillColumn fills only a column left all NA", {
    md <- S4Vectors::DataFrame(a = c(1, 2), b = c(NA_real_, NA_real_))
    derived <- list(a = c(9, 9), b = c(3, 4))
    expect_identical(pecotmr:::.fmRowFillColumn(md, "a", derived), md)
    filled <- pecotmr:::.fmRowFillColumn(md, "b", derived)
    expect_identical(filled$b, c(3, 4))
    expect_identical(filled$a, md$a)
})

test_that("a row saved before rows recorded derived columns derives none", {
    row <- FineMappingRow(
        .fmr_ids,
        list(pip = c(0.4, 0.6)),
        data.frame(variant_id = .fmr_ids, pip = c(0.4, 0.6))
    )
    old <- row
    attr(old, "derivedColumns") <- NULL
    expect_false(methods::.hasSlot(old, "derivedColumns"))
    expect_identical(pecotmr:::.fmRowDerivedColumns(old), character(0))
})
