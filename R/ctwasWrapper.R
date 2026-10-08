# =============================================================================
# cTWAS engine interface
# -----------------------------------------------------------------------------
# Every call into the `ctwas` package lives here; ctwasPipeline.R orchestrates
# the three steps and holds the Param bundles. The split is mechanical: a
# function belongs here if it calls `ctwas::`, or is the Options bag that
# forwards arguments to one of those calls.
# =============================================================================

# Every ctwas function the pipeline fans `CtwasOptions` out to. The bundle is
# validated against the union of their explicit formals, which is sound
# because .ctwasInvoke filters it per callee (see there).
# @noRd
.ctwasCallees <- function() {
    c(
        "ctwas::compute_gene_z",
        "ctwas::est_param",
        "ctwas::screen_regions",
        "ctwas::finemap_regions",
        "ctwas::expand_region_data",
        "ctwas::postprocess_region_merging",
        "ctwas::postprocess_region_merging_noLD",
        "ctwas::get_boundary_genes"
    )
}

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
            extra = methodArgs
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
        extra = methodArgs
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
.ctwasFinemapOrEmpty <- function(screenResult, L, numThreads, methodArgs) {
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
            L = as.integer(L),
            LD_format = "custom",
            LD_loader_fun = screenResult$LD_loader_fun,
            snpinfo_loader_fun = screenResult$snpinfo_loader_fun,
            ncore = as.integer(numThreads)
        ),
        extra = methodArgs
    )
}

# Pick the LD vs no-LD region-merging fn + args. ctwas's postprocess_*()
# forward `...` into finemap_regions, so the LD loader closures must ride in the
# explicit arg list (not filtered through .ctwasInvoke).
# @noRd
.ctwasMergeDispatch <- function(finemapResult, common, L) {
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
            L = as.integer(L),
            LD_format = "custom",
            LD_loader_fun = finemapResult$LD_loader_fun,
            snpinfo_loader_fun = finemapResult$snpinfo_loader_fun
        )
    )
    list(fn = ctwas::postprocess_region_merging, args = args)
}

#' @title Arguments For The ctwas Fitting Steps
#' @description Options forwarded to ctwas. \code{ctwasPipeline} runs seven
#'   ctwas functions in sequence and hands this one bundle to each, forwarding
#'   only the arguments that step actually accepts -- so a setting meant for
#'   the screening step does not disturb the fitting step. Names are therefore
#'   checked against what those seven functions accept \emph{between them}.
#'
#'   The interface is deliberately flat, matching ctwas's own
#'   \code{ctwas_sumstats}, which takes forty arguments in one signature and
#'   distributes them to its internal steps itself.
#' @param ... Any argument accepted by one of the ctwas steps, under ctwas's
#'   own names (\code{niter}, \code{min_gene}, \code{min_nonSNP_PIP},
#'   \code{numThreads}, ...). Note that settings pecotmr exposes as its own
#'   parameters --- \code{thin}, \code{L}, \code{numThreads} on the pipeline ---
#'   are passed there, not here.
#' @return A \code{\link{MethodOptions}} object.
#' @examples
#' CtwasOptions(min_group_size = 1, min_gene = 1)
#' @export
CtwasOptions <- function(...) {
    extra <- list(...)
    .ctwasRefusePipelineOwned(extra)
    .newMethodOptions(
        .ctwasCallees(),
        defaults = list(),
        extra = extra,
        label = "CtwasOptions",
        engine = "ctwas",
        # .ctwasInvoke intersects the bundle with each step's explicit
        # formals, so the union check is sound even though three of the seven
        # take `...`.
        filtered = TRUE
    )
}
