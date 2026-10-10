context("ctwasWrapper")

# =============================================================================
# Tests for the ctwas engine interface (R/ctwasWrapper.R): CtwasOptions and
# the screen step. The pipeline orchestration tests, and the ctwas:: helpers
# that are exercised only through a full run, stay in test_ctwasPipeline.R.
# =============================================================================

test_that("the geneZ bundle reaches compute_gene_z", {
    skip_if_not_installed("ctwas")
    # compute_gene_z was the one ctwas step methodArgs could not configure:
    # it was called directly rather than through .ctwasInvoke, so its
    # `logfile` was unreachable. It now has its own bundle.
    expect_s4_class(CtwasGeneZOptions(logfile = "gene-z.log"), "MethodOptions")
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
        extra = CtwasGeneZOptions(logfile = "gene-z.log")
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
    expect_true(
        "maxNumSingleEffects" %in% names(formals(finemapCtwasRegions))
    )
})

test_that("CtwasOptions nests one bundle per ctwas step", {
    skip_if_not_installed("ctwas")
    a <- CtwasOptions(
        estParam = CtwasEstParamOptions(min_group_size = 2L),
        screen = CtwasScreenOptions(min_gene = 1L)
    )
    expect_s4_class(a, "MethodOptions")
    expect_equal(
        names(a),
        c(
            "geneZ", "regionData", "estParam", "screen", "expand",
            "finemap", "merge", "boundaryGenes"
        )
    )
    expect_equal(a$estParam$min_group_size, 2L)
    expect_equal(a$screen$min_gene, 1L)
    # A bare list is refused the same way every other bundle refuses one.
    expect_error(CtwasOptions(screen = list(min_gene = 1L)), "CtwasScreen")
})

test_that("each step's bundle accepts only that step's names", {
    skip_if_not_installed("ctwas")
    # This is what the flat bundle could not do: it was checked against the
    # UNION of every step's formals, so of 55 names 29 were accepted by
    # exactly one step -- and a name aimed at the wrong step was taken and
    # then silently dropped. Now it is a construction-time error.
    expect_s4_class(CtwasScreenOptions(min_nonSNP_PIP = 0), "MethodOptions")
    expect_error(
        CtwasFinemapOptions(min_nonSNP_PIP = 0),
        "unknown argument"
    )
    expect_s4_class(CtwasFinemapOptions(min_abs_corr = 0.1), "MethodOptions")
    expect_error(CtwasScreenOptions(min_abs_corr = 0.1), "unknown argument")
})

test_that("the region-assembly step is no longer undeclared", {
    skip_if_not_installed("ctwas")
    # assemble_region_data was missing from the flat bundle's callee set even
    # though the pipeline calls it, so four of its real options were
    # rejected as unknown.
    for (nm in c("trim_by", "thin_by", "adjust_boundary_genes", "seed")) {
        expect_s4_class(
            do.call(CtwasRegionDataOptions, setNames(list(1L), nm)),
            "MethodOptions"
        )
    }
})

test_that("the estParam bundle covers the prefit-EM fallback too", {
    skip_if_not_installed("ctwas")
    # The pipeline falls back from est_param() to ctwas's internal fit_EM(),
    # whose formals are NOT a subset -- so the bundle is checked against
    # both and .ctwasInvoke drops whichever the running step lacks.
    expect_s4_class(
        CtwasEstParamOptions(warn_converge_fail = FALSE), "MethodOptions"
    )
    expect_s4_class(CtwasEstParamOptions(EM_tol = 1e-4), "MethodOptions")
})

test_that("the merge bundle also accepts its fine-mapping rerun's names", {
    skip_if_not_installed("ctwas")
    # Both postprocess_region_merging*() forward `...` into a fine-mapping
    # rerun, so finemap_regions' formals are legitimate here -- twelve names
    # reach the rerun and nothing else. That is why the merge call must not
    # be filtered down to the merge function's own formals.
    expect_s4_class(CtwasMergeOptions(combine_PIPs = TRUE), "MethodOptions")
    expect_s4_class(CtwasMergeOptions(min_abs_corr = 0.1), "MethodOptions")
    expect_s4_class(CtwasMergeOptions(coverage = 0.9), "MethodOptions")
})

test_that("a step bundle rejects a name no ctwas step accepts", {
    skip_if_not_installed("ctwas")
    # `filter_L` is the real case that prompted this: it is not a formal of
    # any ctwas function, so it was silently dropped for years.
    expect_error(
        CtwasEstParamOptions(filter_L = FALSE),
        "unknown argument\\(s\\) filter_L"
    )
    expect_error(CtwasScreenOptions(min_genee = 1L), "unknown argument")
})

test_that("step bundles refuse pecotmr's own pipeline parameters", {
    skip_if_not_installed("ctwas")
    # These ARE real ctwas formals, so a formals check alone would wave them
    # through -- but the pipeline supplies them, so a value set here was
    # silently discarded. The error names the setting that takes effect, and
    # it is now raised by the step that would have dropped it.
    expect_error(
        CtwasRegionDataOptions(thin = 0.1), "CtwasPriorParam\\(thin =\\)"
    )
    expect_error(
        CtwasEstParamOptions(niter = 10L), "CtwasPriorParam\\(niter =\\)"
    )
    expect_error(
        CtwasEstParamOptions(group_prior_var_structure = "shared_all"),
        "CtwasPriorParam\\(varStructure =\\)"
    )
    expect_error(CtwasFinemapOptions(L = 5L), "the pipeline's own")
    expect_error(CtwasGeneZOptions(ncore = 4L), "the pipeline's own")
    expect_error(
        CtwasMergeOptions(maxSNP = 10L), "BoundaryMergeParam\\(maxSnp =\\)"
    )
})
