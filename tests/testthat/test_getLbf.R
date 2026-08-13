context("getLbf")

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
    e <- FineMappingEntry(
        variantIds = vids,
        susieFit = fit,
        topLoci = data.frame(
            variant_id = vids,
            pip = fit$pip,
            stringsAsFactors = FALSE
        )
    )
    w <- getLbf(e)
    expect_equal(w$variant_id, vids)
    expect_true(all(c("lbf_L1", "lbf_L2") %in% names(w)))
    expect_equal(w$lbf_L1, lbf[1, ]) # effect 1's lbf across variants
    expect_equal(w$lbf_L2, lbf[2, ]) # effect 2's lbf across variants
})

test_that("getLbf is empty when the fit carries no lbf matrix", {
    e <- FineMappingEntry(
        variantIds = "chr1:100:A:G",
        susieFit = list(pip = 0.5),
        topLoci = data.frame(
            variant_id = "chr1:100:A:G",
            pip = 0.5,
            stringsAsFactors = FALSE
        )
    )
    expect_equal(nrow(getLbf(e)), 0L)
})

test_that("getLbf aggregates over a collection, carrying identity + NA-filling ragged effects", {
    vids <- c("chr1:100:A:G", "chr1:200:C:T", "chr1:300:G:A")
    mkEntry <- function(L) {
        lbf <- matrix(seq_len(L * length(vids)), nrow = L, byrow = TRUE) # L x p
        fit <- structure(
            list(lbf_variable = lbf, pip = rep(0.3, length(vids))),
            class = "susie"
        )
        FineMappingEntry(
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
