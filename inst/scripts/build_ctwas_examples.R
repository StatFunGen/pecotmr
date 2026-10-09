# =============================================================================
# Rebuild the bundled cTWAS example objects.
# -----------------------------------------------------------------------------
# ctwasInputsExample / ctwasEstExample / ctwasFinemapExample are the payloads
# of the granular chain
#
#   assembleCtwasInputs -> estCtwasGroupPriors -> screenCtwasRegions
#     -> finemapCtwasRegions
#
# captured after each step, so a user can enter the chain anywhere without
# re-running what came before.
#
# The payloads carry one `LD_map$LD_file` token PER REGION, which ctwas
# asserts exists on disk and which pecotmr also uses as the dispatch key into
# the cached per-region LD panels. Those two duties are why the token cannot
# simply be the backing genotype file: every block on a chromosome shares that
# file, so each block was handed the FIRST block's LD panel. The tokens are
# per-region sentinels under the session tempdir instead, and are therefore
# gone when a saved object is read back -- `.ctwasResolveLdPaths()` re-mints
# one per region on the way in and re-keys the loader cache onto them, so no
# path from the building machine is ever baked into the .rda.
#
# Run from the package root:
#   Rscript inst/scripts/build_ctwas_examples.R
# =============================================================================

devtools::load_all(".", quiet = TRUE)

# A bundled stem as the portable reference the serialized objects carry, so
# no .rda holds a path from the machine that built it.
asBundledResource <- function(stem) {
    rel <- sub(
        paste0("^.*", .Platform$file.sep, "extdata", .Platform$file.sep),
        "",
        stem
    )
    paste0("pecotmr://extdata/", rel)
}

# A GenotypeHandle for `stem`, carrying the portable reference rather than
# the concrete path it had to be read from.
.portableHandle <- function(stem) {
    handle <- pecotmr:::.ldSketchHandle(readGenotypes(plink2Prefix = stem))
    handle@path <- asBundledResource(stem)
    handle
}

# -----------------------------------------------------------------------------
# 1. Inputs: the bundled chr22 LD panel, GWAS sumstats, and TWAS weights.
# -----------------------------------------------------------------------------
ldStem <- file.path(
    system.file("extdata", "ld_reference", "chr22", package = "pecotmr"),
    "protocol_example.LD.chr22"
)
gwasTsv <- system.file(
    "extdata",
    "manifests",
    "protocol_example.twas.gwas_sumstats.chr22.tsv.gz",
    package = "pecotmr"
)

# cTWAS needs at least two LD blocks: its EM cannot converge on one region.
blocks <- GenomicRanges::GRanges(
    "chr22",
    IRanges::IRanges(c(10000000, 15000001), c(15000000, 19000000)),
    blockId = c("chr22_1", "chr22_2")
)

gwasSumStats <- loadGwasSumStatsFromManifest(
    manifest = data.frame(study = "gwas1", sumStatsPath = gwasTsv),
    genome = "hg38",
    ldSketch = ldStem,
    region = "chr22:10000000-19000000",
    ldBlocks = blocks
)
gwasByRegion <- summaryStatsQc(
    gwasSumStats,
    panelFilterArgs = PanelFilterParam(mafCutoff = 0.0025)
)

data(ctwasWeightsExample)

# -----------------------------------------------------------------------------
# 1b. Point the weight example's LD sketch at the bundled panel.
# -----------------------------------------------------------------------------
# ctwasWeightsExample was originally built against a local xqtl-protocol
# checkout, so it carried that machine's absolute genotype path and a
# 1975-variant sketch read from the untrimmed source. All 168 of its weight
# variants lie inside the bundled 175-variant chr22 panel, and the two span
# exactly the same range, so the bundled copy is the right panel rather than
# a substitute for it.
#
# Construction needs a real path (it reads .pvar to populate snpInfo), so
# the handle is built from the concrete stem and its path rewritten to the
# portable "pecotmr://extdata/<stem>" form afterwards -- the same two-step
# every other bundled example uses. .resolveGenotypeResourcePath() turns it
# back into a real path at read time, on whatever machine has the package.
ctwasWeightsExample <- methods::initialize(
    ctwasWeightsExample,
    ldSketch = pecotmr:::.asLdSketch(.portableHandle(ldStem))
)

stopifnot(
    getPath(pecotmr:::.ldSketchHandle(
        getLdSketch(ctwasWeightsExample)
    )) == asBundledResource(ldStem)
)

# -----------------------------------------------------------------------------
# 2. Run the chain, capturing each payload.
# -----------------------------------------------------------------------------
ctwasInputsExample <- assembleCtwasInputs(
    gwasSumStats = gwasByRegion,
    twasWeights = list(ctwasWeightsExample)
)

ctwasEstExample <- estCtwasGroupPriors(
    ctwasInputsExample,
    ctwasPriorArgs = CtwasPriorParam(
        thin = 1,
        niterPrefit = 3,
        niter = 10,
        fallbackToPrefit = TRUE
    ),
    methodArgs = CtwasOptions(min_group_size = 1, min_p_single_effect = 0)
)

# The toy GWAS carries no genome-wide-significant signal, so the default
# screen (min_nonSNP_PIP = 0.5) selects nothing and the example would be an
# empty result. Keep every region so the finemap payload is populated.
screened <- screenCtwasRegions(
    ctwasEstExample,
    methodArgs = CtwasOptions(min_nonSNP_PIP = 0)
)
ctwasFinemapExample <- finemapCtwasRegions(screened)

# -----------------------------------------------------------------------------
# 3. Materialize each region's LD and let its sketch go.
# -----------------------------------------------------------------------------
# The pipeline keeps a sketch per panel so a region's n x n correlation is
# built only if ctwas actually fine-maps that region. A bundled example has
# the opposite requirement: it has to work in an installation with no access
# to the machine that built it, and a serialized sketch carries that
# machine's genotype path. These panels are small, so materialize the
# correlation and drop the sketch -- .ctwasPanelLd() uses a panel's own R
# when it has one, which is the path every reader of these objects takes.
materializePanel <- function(p) {
    list(
        R = pecotmr:::.ctwasPanelLd(p),
        snpInfo = p$snpInfo,
        variance = p$variance
    )
}

materializeLd <- function(payload) {
    panels <- pecotmr:::.ctwasCachedPanels(payload)
    if (is.null(panels)) {
        return(payload)
    }
    built <- lapply(panels, materializePanel)
    names(built) <- names(panels)
    # Both loaders dispatch on the LD_file token, so they are rebuilt over
    # the materialized cache exactly as .ctwasAssembleResult() builds them.
    payload$LD_loader_fun <- pecotmr:::.ctwasMultiBlockLdLoader(built)
    payload$snpinfo_loader_fun <-
        pecotmr:::.ctwasMultiBlockSnpInfoLoader(built)
    payload
}

ctwasInputsExample <- materializeLd(ctwasInputsExample)
ctwasEstExample <- materializeLd(ctwasEstExample)
ctwasFinemapExample <- materializeLd(ctwasFinemapExample)

stopifnot(
    all(!vapply(
        pecotmr:::.ctwasCachedPanels(ctwasInputsExample),
        function(p) is.null(p$R) || !is.null(p$sketch),
        logical(1)
    ))
)

# -----------------------------------------------------------------------------
# 4. Check the token invariant.
# -----------------------------------------------------------------------------
# Every region now carries its OWN LD_file token -- a sentinel under the
# session tempdir -- because ctwas asserts file.exists() on each token AND
# dispatches LD_loader_fun on that same string, so pointing them all at the
# one shared genotype file handed every block the FIRST block's LD. The
# sentinels are gone by the time a saved object is read back, which is
# precisely the case .ctwasResolveLdPaths() handles: it re-mints one per
# region and re-keys the loader cache onto the new tokens. So nothing is
# rewritten here. The old "pecotmr://extdata/<path>" step only worked while
# the token was a real file under extdata and would mangle a sentinel path.
#
# What is worth asserting is the invariant the bug violated: one distinct
# token per region.
stopifnot(
    anyDuplicated(ctwasInputsExample$LD_map$LD_file) == 0L,
    anyDuplicated(ctwasEstExample$LD_map$LD_file) == 0L,
    anyDuplicated(ctwasFinemapExample$LD_map$LD_file) == 0L
)

# -----------------------------------------------------------------------------
# 5. Verify each payload still drives the step that consumes it, from the
#    portable form -- the check the previous objects would have failed.
# -----------------------------------------------------------------------------
# Drop the live sentinels first. Without this the verification below runs
# against tokens that still exist in this session and so would not exercise
# the re-mint path every reader of the saved object actually takes.
unlink(unique(c(
    ctwasInputsExample$LD_map$LD_file,
    ctwasEstExample$LD_map$LD_file,
    ctwasFinemapExample$LD_map$LD_file
)))

invisible(estCtwasGroupPriors(
    ctwasInputsExample,
    ctwasPriorArgs = CtwasPriorParam(
        thin = 1,
        niterPrefit = 3,
        niter = 10,
        fallbackToPrefit = TRUE
    ),
    methodArgs = CtwasOptions(min_group_size = 1, min_p_single_effect = 0)
))
invisible(finemapCtwasRegions(
    screenCtwasRegions(
        ctwasEstExample,
        methodArgs = CtwasOptions(min_nonSNP_PIP = 0)
    )
))
invisible(asCtwasResult(ctwasFinemapExample))
invisible(mergeCtwasBoundaryRegions(ctwasFinemapExample))

# -----------------------------------------------------------------------------
# 6. Save.
# -----------------------------------------------------------------------------
usethis::use_data(ctwasWeightsExample, overwrite = TRUE, compress = "xz")
usethis::use_data(ctwasInputsExample, overwrite = TRUE, compress = "xz")
usethis::use_data(ctwasEstExample, overwrite = TRUE, compress = "xz")
usethis::use_data(ctwasFinemapExample, overwrite = TRUE, compress = "xz")
