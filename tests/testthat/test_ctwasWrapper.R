context("ctwasWrapper")

# =============================================================================
# Tests for the ctwas engine interface (R/ctwasWrapper.R): CtwasOptions and
# the screen step. The pipeline orchestration tests, and the ctwas:: helpers
# that are exercised only through a full run, stay in test_ctwasPipeline.R.
# =============================================================================

test_that("CtwasOptions now reaches compute_gene_z too", {
    skip_if_not_installed("ctwas")
    # compute_gene_z was the one ctwas step methodArgs could not configure:
    # it was called directly rather than through .ctwasInvoke, so its
    # `logfile` was unreachable and the name was not in the accepted union.
    expect_true("ctwas::compute_gene_z" %in% pecotmr:::.ctwasCallees())
    expect_s4_class(CtwasOptions(logfile = "gene-z.log"), "MethodOptions")
    seen <- NULL
    local_mocked_bindings(
        # The mock must carry the REAL formals: .ctwasInvoke filters the
        # bundle to the callee's explicit formals, so a `function(...)` mock
        # would drop `logfile` and hide the thing being tested.
        compute_gene_z = function(z_snp, weights, ncore = 1L, logfile = NULL) {
            seen <<- list(ncore = ncore, logfile = logfile)
            data.frame(id = "t1", z = 1.0)
        },
        .package = "ctwas"
    )
    inputs <- list(
        z_snp = data.frame(id = "s1", z = 1),
        weights = list(t1 = list(wgt = 1)),
        z_gene = NULL
    )
    pecotmr:::.ctwasEnsureZGene(
        inputs,
        numThreads = 1L,
        extra = CtwasOptions(logfile = "gene-z.log")
    )
    expect_equal(seen$logfile, "gene-z.log")
    expect_equal(seen$ncore, 1L)
    # A supplied z_gene still short-circuits without calling the engine.
    seen <- NULL
    inputs$z_gene <- data.frame(id = "t1", z = 2.0)
    pecotmr:::.ctwasEnsureZGene(inputs, numThreads = 1L)
    expect_null(seen)
})

test_that("screenCtwasRegions no longer advertises an L it cannot use", {
    skip_if_not_installed("ctwas")
    # ctwas::screen_regions has no `L` formal -- screening is always SER --
    # so the parameter could never be honoured. It was documented as
    # "Unused. Retained for call-site compatibility"; a knob that cannot
    # work is better removed than explained.
    expect_false("L" %in% names(formals(screenCtwasRegions)))
    expect_false("L" %in% names(formals(ctwas::screen_regions)))
    # finemapCtwasRegions is where L genuinely applies, and still takes it.
    expect_true("L" %in% names(formals(finemapCtwasRegions)))
})

test_that("CtwasOptions validates against what the ctwas steps accept together", {
    a <- CtwasOptions(min_group_size = 2L, min_gene = 1L)
    expect_s4_class(a, "MethodOptions")
    expect_equal(a$min_group_size, 2L)
    # A name valid for only ONE of the steps is still accepted here; the
    # per-step filtering in .ctwasInvoke routes it to the step that takes it.
    expect_true(is_in(
        "min_nonSNP_PIP",
        names(CtwasOptions(min_nonSNP_PIP = 0))
    ))
    expect_true(is_in("min_gene", names(CtwasOptions(min_gene = 1L))))
})

test_that("CtwasOptions rejects a name no ctwas step accepts", {
    # `filter_L` is the real case that prompted this: it is not a formal of
    # any ctwas function, so it was silently dropped for years.
    expect_error(
        CtwasOptions(filter_L = FALSE),
        "unknown argument\\(s\\) filter_L"
    )
    expect_error(CtwasOptions(min_genee = 1L), "unknown argument")
})

test_that("CtwasOptions does not accept pecotmr's own pipeline parameters", {
    # These ARE real ctwas formals, so the union check would wave them
    # through -- but .ctwasInvoke drops any name the pipeline already
    # supplies, so a value set here was silently discarded. The error names
    # the setting that actually takes effect.
    expect_error(CtwasOptions(thin = 0.1), "CtwasPriorParam\\(thin =\\)")
    expect_error(CtwasOptions(niter = 10L), "CtwasPriorParam\\(niter =\\)")
    expect_error(
        CtwasOptions(group_prior_var_structure = "shared_all"),
        "CtwasPriorParam\\(varStructure =\\)"
    )
    expect_error(CtwasOptions(L = 5L), "the pipeline's own")
    expect_error(CtwasOptions(ncore = 4L), "the pipeline's own")
    expect_error(CtwasOptions(maxSNP = 10L), "BoundaryMergeParam\\(maxSnp =\\)")
})
