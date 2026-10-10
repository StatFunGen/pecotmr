# =============================================================================
# cTWAS engine interface
# -----------------------------------------------------------------------------
# Every call into the `ctwas` package lives here; ctwasPipeline.R orchestrates
# the three steps and holds the Param bundles. The split is mechanical: a
# function belongs here if it calls `ctwas::`, or is the Options bag that
# forwards arguments to one of those calls.
# =============================================================================

# z_gene for assemble_region_data (which requires it non-NULL): use the caller's
# inputs$z_gene, else compute via ctwas::compute_gene_z (mirrors
# ctwas_sumstats). Routed through .ctwasInvoke so `methodArgs` reaches it --
# this was the one ctwas step CtwasOptions() could not configure.
# @noRd
.ctwasEnsureZGene <- function(inputs, numThreads, extra = list()) {
    if (!is.null(inputs$z_gene)) {
        return(inputs$z_gene)
    }
    .ctwasInvoke(
        ctwas::compute_gene_z,
        list(
            z_snp = inputs$z_snp,
            weights = inputs$weights,
            ncore = numThreads
        ),
        extra = extra
    )
}

# assemble_region_data -> the per-region list (keyed by region_id). It does NOT
# echo z_gene or return boundary genes; both are recovered separately.
# @noRd
.ctwasAssembleRegionData <- function(inputs, zGene, thin, numThreads, extra) {
    .ctwasInvoke(
        ctwas::assemble_region_data,
        list(
            region_info = inputs$region_info,
            z_snp = inputs$z_snp,
            z_gene = zGene,
            weights = inputs$weights,
            snp_map = inputs$snp_map,
            thin = thin,
            ncore = numThreads
        ),
        extra = extra
    )
}

# Boundary genes (computed internally by assemble_region_data for adjustment but
# never returned) recovered via ctwas::get_boundary_genes; NULL for one region.
# @noRd
.ctwasBoundaryGenes <- function(inputs, numThreads, extra) {
    if (nrow(inputs$region_info) <= 1L) {
        return(NULL)
    }
    .ctwasInvoke(
        ctwas::get_boundary_genes,
        list(
            region_info = inputs$region_info,
            weights = inputs$weights,
            ncore = numThreads
        ),
        extra = extra
    )
}

# The accurate EM (ctwas::est_param) invocation, split out so
# .ctwasEstParamOrFallback can run it directly (no fallback) or inside a
# tryCatch
# (with fallback) without duplicating the argument assembly.
# @noRd
.ctwasEstParamAccurate <- function(
    regionData,
    niterPrefit,
    niter,
    groupPriorVarStructure,
    numThreads,
    extra
) {
    .ctwasInvoke(
        ctwas::est_param,
        list(
            region_data = regionData,
            niter_prefit = as.integer(niterPrefit),
            niter = as.integer(niter),
            group_prior_var_structure = groupPriorVarStructure,
            ncore = numThreads
        ),
        extra = extra
    )
}

#' Screen cTWAS regions
#'
#' @description Step 3 of the three-step \code{\link{ctwasPipeline}}: runs
#'   \code{ctwas::screen_regions} on the \code{\link{estCtwasGroupPriors}}
#'   result and returns the screened-region set. Use this entry point to
#'   substitute hand-tuned priors for the ones estimated in step 2 (e.g. when
#'   the accurate EM diverges to NaN and you want to recover the prefit
#'   values).
#'
#' @param estResult A list returned by \code{\link{estCtwasGroupPriors}}.
#' @param numThreads Number of cores.
#' @param methodArgs Additional arguments forwarded to ctwas, built
#'   with \code{\link{CtwasOptions}}. Names are checked against what
#'   the ctwas steps accept between them.
#' @return The \code{estResult} list augmented with \code{screen_res} (the full
#'   ctwas output) and \code{screened_region_data}.
#' @importFrom purrr map compact
#' @examples
#' data(ctwasEstExample)
#' screenCtwasRegions(ctwasEstExample)
#' @export
screenCtwasRegions <- function(
    estResult,
    numThreads = 1L,
    methodArgs = CtwasOptions()
) {
    .assertMethodOptions(methodArgs, "CtwasOptions", "methodArgs")
    if (!requireNamespace("ctwas", quietly = TRUE)) {
        abort("Package 'ctwas' is required for screenCtwasRegions.")
    }
    estResult <- .ctwasResolveLdPaths(estResult)
    # ctwas::screen_regions requires thin = 1 region_data; expand the
    # thinned set first when assemble_region_data was called with thin < 1
    # (matches ctwas_sumstats's own expand-before-screen step).
    thinVals <- compact(map(estResult$region_data, "thin"))
    needsExpand <- length(thinVals) > 0L && min(list_c(thinVals)) < 1
    regionDataForScreen <- if (needsExpand) {
        .ctwasInvoke(
            ctwas::expand_region_data,
            list(
                region_data = estResult$region_data,
                snp_map = estResult$snp_map,
                z_snp = estResult$z_snp,
                ncore = as.integer(numThreads)
            ),
            extra = methodArgs$expand
        )
    } else {
        estResult$region_data
    }
    screenRes <- .ctwasInvoke(
        ctwas::screen_regions,
        list(
            region_data = regionDataForScreen,
            group_prior = estResult$param$group_prior,
            group_prior_var = estResult$param$group_prior_var,
            ncore = as.integer(numThreads)
        ),
        extra = methodArgs$screen
    )
    c(
        estResult,
        list(
            screen_res = screenRes,
            screened_region_data = screenRes$screened_region_data
        )
    )
}

# Fine-map the screened regions, or hand back the empty shape when the
# screen kept none -- ctwas::finemap_regions() has no meaningful answer for
# an empty region set, and the caller's result still needs both slots.
# @noRd
.ctwasFinemapOrEmpty <- function(
    screenResult,
    maxNumSingleEffects,
    numThreads,
    methodArgs
) {
    rd <- screenResult$screened_region_data
    if (length(rd) == 0L) {
        return(list(finemap_res = NULL, susie_alpha_res = NULL))
    }
    .ctwasInvoke(
        ctwas::finemap_regions,
        list(
            region_data = rd,
            LD_map = screenResult$LD_map,
            weights = screenResult$weights,
            group_prior = screenResult$param$group_prior,
            group_prior_var = screenResult$param$group_prior_var,
            L = as.integer(maxNumSingleEffects),
            LD_format = "custom",
            LD_loader_fun = screenResult$LD_loader_fun,
            snpinfo_loader_fun = screenResult$snpinfo_loader_fun,
            ncore = as.integer(numThreads)
        ),
        extra = methodArgs$finemap
    )
}

# Pick the LD vs no-LD region-merging fn + args. ctwas's postprocess_*()
# forward `...` into finemap_regions, so the LD loader closures must ride in the
# explicit arg list (not filtered through .ctwasInvoke).
# @noRd
.ctwasMergeDispatch <- function(finemapResult, common, maxNumSingleEffects) {
    if (is.null(finemapResult$LD_loader_fun)) {
        return(list(
            fn = ctwas::postprocess_region_merging_noLD,
            args = common
        ))
    }
    args <- c(
        common,
        list(
            LD_map = finemapResult$LD_map,
            L = as.integer(maxNumSingleEffects),
            LD_format = "custom",
            LD_loader_fun = finemapResult$LD_loader_fun,
            snpinfo_loader_fun = finemapResult$snpinfo_loader_fun
        )
    )
    list(fn = ctwas::postprocess_region_merging, args = args)
}

# Every ctwas step bundle is built the same way: refuse what the pipeline
# supplies, then check the rest against that step's live formals.
#
# `filtered = TRUE` because several ctwas steps take `...`, which would
# otherwise make the record accept any name -- exactly what the per-step
# split exists to prevent. For the steps .ctwasInvoke() drives it is literally
# true (it intersects the bundle with the running callee's formals); for the
# merge step nothing filters, but every name the bundle accepts is a formal
# of either the merge call or the fine-mapping rerun it forwards to, so the
# union is the right set either way.
# @noRd
.ctwasStepOptions <- function(extra, owned, label, engine, callees) {
    .configRefuseOwned(extra, owned, label)
    .newMethodOptions(
        callees,
        defaults = list(),
        extra = extra,
        label = label,
        engine = engine,
        filtered = TRUE
    )
}

# What the pipeline itself supplies at each step's call site -- taken from the
# base `list()` handed to .ctwasInvoke(), which is exactly the set that gets
# deduped away. `ncore` is the pipeline's `numThreads` everywhere.
# @noRd
.ctwasOwnedGeneZ <- function() {
    c(
        z_snp = "supplied from the data by the pipeline",
        weights = "supplied from the data by the pipeline",
        ncore = "the pipeline's own `numThreads`"
    )
}

#' @title Arguments For ctwas's Gene Z-Score Step
#' @description Arguments for \code{ctwas::compute_gene_z()}, which derives
#'   gene-level z-scores from the SNP z-scores and the TWAS weights.
#' @param ... Any \code{ctwas::compute_gene_z()} argument, under ctwas's own
#'   names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasGeneZOptions(logfile = "gene-z.log")
#' @export
CtwasGeneZOptions <- function(...) {
    .ctwasStepOptions(
        list(...), .ctwasOwnedGeneZ(), "CtwasGeneZOptions", "ctwasGeneZ",
        "ctwas::compute_gene_z"
    )
}

#' @title Arguments For ctwas's Region-Assembly Step
#' @description Arguments for \code{ctwas::assemble_region_data()}, which
#'   builds the per-region list the EM and fine-mapping steps consume.
#'
#'   This step was absent from the flat bundle's callee set, so four of its
#'   real options --- \code{trim_by}, \code{thin_by},
#'   \code{adjust_boundary_genes}, \code{seed} --- were rejected as unknown
#'   even though the pipeline does call it.
#' @param ... Any \code{ctwas::assemble_region_data()} argument, under
#'   ctwas's own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasRegionDataOptions(trim_by = "z", seed = 1L)
#' @export
CtwasRegionDataOptions <- function(...) {
    .ctwasStepOptions(
        list(...),
        c(
            region_info = "supplied from the data by the pipeline",
            z_snp = "supplied from the data by the pipeline",
            z_gene = "computed by the pipeline from the TWAS weights",
            weights = "supplied from the data by the pipeline",
            snp_map = "supplied from the data by the pipeline",
            thin = "CtwasPriorParam(thin =)",
            ncore = "the pipeline's own `numThreads`"
        ),
        "CtwasRegionDataOptions", "ctwasRegionData",
        "ctwas::assemble_region_data"
    )
}

#' @title Arguments For ctwas's Group-Prior EM
#' @description Arguments for the EM that estimates ctwas's group priors.
#'   Two callees, because the pipeline falls back from
#'   \code{ctwas::est_param()} to ctwas's internal \code{fit_EM()} when the
#'   accurate EM diverges; \code{fit_EM()} is not a subset of
#'   \code{est_param()} (\code{groups}, \code{types}, \code{contexts},
#'   \code{warn_converge_fail} are its own), so the bundle is checked against
#'   both and \code{.ctwasInvoke()} drops whichever the running step does not
#'   take.
#' @param ... Any argument of either EM entry point, under ctwas's own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasEstParamOptions(min_group_size = 1, EM_tol = 1e-4)
#' @export
CtwasEstParamOptions <- function(...) {
    .ctwasStepOptions(
        list(...),
        c(
            region_data = "supplied from the data by the pipeline",
            niter_prefit = "CtwasPriorParam(niterPrefit =)",
            niter = "CtwasPriorParam(niter =)",
            group_prior_var_structure = "CtwasPriorParam(varStructure =)",
            ncore = "the pipeline's own `numThreads`"
        ),
        "CtwasEstParamOptions", "ctwasEstParam",
        c("ctwas::est_param", "ctwas::fit_EM")
    )
}

#' @title Arguments For ctwas's Region-Screening Step
#' @description Arguments for \code{ctwas::screen_regions()}.
#' @param ... Any \code{ctwas::screen_regions()} argument, under ctwas's own
#'   names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasScreenOptions(min_nonSNP_PIP = 0, min_gene = 1)
#' @export
CtwasScreenOptions <- function(...) {
    .ctwasStepOptions(
        list(...),
        c(
            region_data = "supplied from the data by the pipeline",
            group_prior = "estimated by the EM step",
            group_prior_var = "estimated by the EM step",
            ncore = "the pipeline's own `numThreads`"
        ),
        "CtwasScreenOptions", "ctwasScreen", "ctwas::screen_regions"
    )
}

#' @title Arguments For ctwas's Region-Expansion Step
#' @description Arguments for \code{ctwas::expand_region_data()}, run before
#'   screening when the region data was assembled with \code{thin < 1}.
#' @param ... Any \code{ctwas::expand_region_data()} argument, under ctwas's
#'   own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasExpandOptions(maxSNP = 1000L)
#' @export
CtwasExpandOptions <- function(...) {
    .ctwasStepOptions(
        list(...),
        c(
            region_data = "supplied from the data by the pipeline",
            snp_map = "supplied from the data by the pipeline",
            z_snp = "supplied from the data by the pipeline",
            ncore = "the pipeline's own `numThreads`"
        ),
        "CtwasExpandOptions", "ctwasExpand", "ctwas::expand_region_data"
    )
}

#' @title Arguments For ctwas's Fine-Mapping Step
#' @description Arguments for \code{ctwas::finemap_regions()}.
#' @param ... Any \code{ctwas::finemap_regions()} argument, under ctwas's own
#'   names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasFinemapOptions(min_abs_corr = 0.1, include_cs = TRUE)
#' @export
CtwasFinemapOptions <- function(...) {
    .ctwasStepOptions(
        list(...),
        c(
            region_data = "supplied from the data by the pipeline",
            LD_map = "supplied from the data by the pipeline",
            weights = "supplied from the data by the pipeline",
            group_prior = "estimated by the EM step",
            group_prior_var = "estimated by the EM step",
            L = "the pipeline's own `maxNumSingleEffects`",
            LD_format = "fixed by the pipeline's LD loaders",
            LD_loader_fun = "fixed by the pipeline's LD loaders",
            snpinfo_loader_fun = "fixed by the pipeline's LD loaders",
            ncore = "the pipeline's own `numThreads`"
        ),
        "CtwasFinemapOptions", "ctwasFinemap", "ctwas::finemap_regions"
    )
}

#' @title Arguments For ctwas's Region-Merging Postprocess
#' @description Arguments for \code{ctwas::postprocess_region_merging()} (or
#'   its \code{_noLD} twin, whose formals are a strict subset).
#'
#'   The accepted set is deliberately WIDER than those functions' own
#'   formals: both of them forward \code{...} into a fine-mapping rerun for
#'   the merged regions, so \code{ctwas::finemap_regions()}'s formals are
#'   legitimate here too --- twelve names (\code{min_abs_corr},
#'   \code{include_cs}, \code{coverage}, ...) reach the rerun and nothing
#'   else. That is why this step's call must NOT be filtered down to the
#'   merge function's explicit formals.
#' @param ... Any argument of the merging postprocess, or of the fine-mapping
#'   rerun it forwards to, under ctwas's own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasMergeOptions(combine_PIPs = TRUE, min_abs_corr = 0.1)
#' @export
CtwasMergeOptions <- function(...) {
    .ctwasStepOptions(
        list(...),
        c(
            region_info = "supplied from the data by the pipeline",
            region_data = "supplied from the data by the pipeline",
            z_snp = "supplied from the data by the pipeline",
            z_gene = "computed by the pipeline from the TWAS weights",
            weights = "supplied from the data by the pipeline",
            snp_map = "supplied from the data by the pipeline",
            finemap_res = "the first-pass fine-mapping result",
            susie_alpha_res = "the first-pass fine-mapping result",
            group_prior = "estimated by the EM step",
            group_prior_var = "estimated by the EM step",
            pip_thresh = "the wrapper's `pipThresh`",
            filter_cs = "the wrapper's `filterCs`",
            maxSNP = "BoundaryMergeParam(maxSnp =)",
            L = "the pipeline's own `maxNumSingleEffects`",
            LD_map = "supplied from the data by the pipeline",
            LD_format = "fixed by the pipeline's LD loaders",
            LD_loader_fun = "fixed by the pipeline's LD loaders",
            snpinfo_loader_fun = "fixed by the pipeline's LD loaders",
            ncore = "the pipeline's own `numThreads`"
        ),
        "CtwasMergeOptions", "ctwasMerge",
        c(
            "ctwas::postprocess_region_merging",
            "ctwas::postprocess_region_merging_noLD",
            "ctwas::finemap_regions"
        )
    )
}

#' @title Arguments For ctwas's Boundary-Gene Step
#' @description Arguments for \code{ctwas::get_boundary_genes()}, which
#'   recovers the boundary genes \code{assemble_region_data()} computes but
#'   does not return.
#' @param ... Any \code{ctwas::get_boundary_genes()} argument, under ctwas's
#'   own names.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasBoundaryGenesOptions(show_mapping = TRUE)
#' @export
CtwasBoundaryGenesOptions <- function(...) {
    .ctwasStepOptions(
        list(...),
        c(
            region_info = "supplied from the data by the pipeline",
            weights = "supplied from the data by the pipeline",
            ncore = "the pipeline's own `numThreads`"
        ),
        "CtwasBoundaryGenesOptions", "ctwasBoundaryGenes",
        "ctwas::get_boundary_genes"
    )
}

#' @title Arguments For The ctwas Fitting Steps
#' @description ctwas is not one function but nine, run in sequence, so each
#'   step's arguments travel in their own bundle and are checked against that
#'   step's own live formals --- the same shape as
#'   \code{\link{CovUdrOptions}}.
#'
#'   This replaces a single flat bundle that was checked against the
#'   \emph{union} of every step's formals. That mirrored ctwas's own
#'   \code{ctwas_sumstats}, which takes forty arguments in one signature, but
#'   it could not tell which step a name was meant for: of the 55 names in
#'   the union, 29 are accepted by exactly one step, so a name aimed at the
#'   wrong step was accepted and then silently dropped. Per-step bundles make
#'   that a construction-time error instead.
#' @param geneZ Arguments for \code{ctwas::compute_gene_z()}, built with
#'   \code{\link{CtwasGeneZOptions}}.
#' @param regionData Arguments for \code{ctwas::assemble_region_data()},
#'   built with \code{\link{CtwasRegionDataOptions}}.
#' @param estParam Arguments for the group-prior EM, built with
#'   \code{\link{CtwasEstParamOptions}}.
#' @param screen Arguments for \code{ctwas::screen_regions()}, built with
#'   \code{\link{CtwasScreenOptions}}.
#' @param expand Arguments for \code{ctwas::expand_region_data()}, built with
#'   \code{\link{CtwasExpandOptions}}.
#' @param finemap Arguments for \code{ctwas::finemap_regions()}, built with
#'   \code{\link{CtwasFinemapOptions}}.
#' @param merge Arguments for the region-merging postprocess, built with
#'   \code{\link{CtwasMergeOptions}}.
#' @param boundaryGenes Arguments for \code{ctwas::get_boundary_genes()},
#'   built with \code{\link{CtwasBoundaryGenesOptions}}.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasOptions(
#'     estParam = CtwasEstParamOptions(min_group_size = 1),
#'     screen = CtwasScreenOptions(min_gene = 1)
#' )
#' @export
CtwasOptions <- function(
    geneZ = CtwasGeneZOptions(),
    regionData = CtwasRegionDataOptions(),
    estParam = CtwasEstParamOptions(),
    screen = CtwasScreenOptions(),
    expand = CtwasExpandOptions(),
    finemap = CtwasFinemapOptions(),
    merge = CtwasMergeOptions(),
    boundaryGenes = CtwasBoundaryGenesOptions()
) {
    .assertMethodOptions(geneZ, "CtwasGeneZOptions", "geneZ")
    .assertMethodOptions(regionData, "CtwasRegionDataOptions", "regionData")
    .assertMethodOptions(estParam, "CtwasEstParamOptions", "estParam")
    .assertMethodOptions(screen, "CtwasScreenOptions", "screen")
    .assertMethodOptions(expand, "CtwasExpandOptions", "expand")
    .assertMethodOptions(finemap, "CtwasFinemapOptions", "finemap")
    .assertMethodOptions(merge, "CtwasMergeOptions", "merge")
    .assertMethodOptions(
        boundaryGenes, "CtwasBoundaryGenesOptions", "boundaryGenes"
    )
    .newMethodOptions(
        NULL,
        defaults = list(
            geneZ = geneZ,
            regionData = regionData,
            estParam = estParam,
            screen = screen,
            expand = expand,
            finemap = finemap,
            merge = merge,
            boundaryGenes = boundaryGenes
        ),
        extra = list(),
        label = "CtwasOptions",
        engine = "ctwas"
    )
}
