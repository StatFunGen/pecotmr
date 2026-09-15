# =============================================================================
# Rebuild the bundled cTWAS example objects.
# -----------------------------------------------------------------------------
# ctwasInputsExample / ctwasEstExample / ctwasFinemapExample are the payloads
# of the granular chain
#
#   assembleCtwasInputs -> estCtwasParam -> screenCtwasRegions
#     -> finemapCtwasRegions
#
# captured after each step, so a user can enter the chain anywhere without
# re-running what came before.
#
# The payloads carry an `LD_map$LD_file` token that ctwas asserts exists on
# disk and that pecotmr uses as the dispatch key into the cached per-region LD
# panels. Baking an absolute path into the .rda makes that token point at the
# machine that built it -- which is what left the previous objects unusable
# past `screenCtwasRegions` (finemapCtwasRegions failed on
# `all(file.exists(LD_matrix_files))`). Step 3 rewrites the token to the
# portable "pecotmr://extdata/<path>" form the genotype handles already use;
# `.ctwasResolveLdPaths()` resolves it at extraction time on whatever machine
# the package is installed.
#
# Run from the package root:
#   Rscript inst/scripts/build_ctwas_examples.R
# =============================================================================

devtools::load_all(".", quiet = TRUE)

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
gwasByRegion <- summaryStatsQc(gwasSumStats, mafCutoff = 0.0025)

data(ctwasWeightsExample)

# -----------------------------------------------------------------------------
# 2. Run the chain, capturing each payload.
# -----------------------------------------------------------------------------
ctwasInputsExample <- assembleCtwasInputs(
    gwasSumStats = gwasByRegion,
    twasWeights = list(ctwasWeightsExample)
)

ctwasEstExample <- estCtwasParam(
    ctwasInputsExample,
    thin = 1,
    niterPrefit = 3,
    niter = 10,
    min_group_size = 1,
    min_p_single_effect = 0,
    fallbackToPrefit = TRUE
)

# The toy GWAS carries no genome-wide-significant signal, so the default
# screen (min_nonSNP_PIP = 0.5) selects nothing and the example would be an
# empty result. Keep every region so the finemap payload is populated.
screened <- screenCtwasRegions(ctwasEstExample, min_nonSNP_PIP = 0)
ctwasFinemapExample <- finemapCtwasRegions(screened)

# -----------------------------------------------------------------------------
# 3. Make the LD token portable.
# -----------------------------------------------------------------------------
asResource <- function(p) {
    rel <- sub(
        paste0("^.*", .Platform$file.sep, "extdata", .Platform$file.sep),
        "",
        p
    )
    paste0("pecotmr://extdata/", rel)
}

repointLd <- function(payload) {
    if (is.null(payload$LD_map)) {
        return(payload)
    }
    stored <- as.character(payload$LD_map$LD_file)
    portable <- vapply(stored, asResource, character(1), USE.NAMES = FALSE)
    keyMap <- setNames(portable, stored)
    payload$LD_map$LD_file <- portable
    payload$LD_map$SNP_file <- portable
    # The loader closures dispatch on the same token, so their cache has to be
    # re-keyed in step with it or the resolved token finds no panel.
    pecotmr:::.ctwasRekeyLdLoaders(payload, keyMap)
}

ctwasInputsExample <- repointLd(ctwasInputsExample)
ctwasEstExample <- repointLd(ctwasEstExample)
ctwasFinemapExample <- repointLd(ctwasFinemapExample)

stopifnot(
    all(grepl("^pecotmr://", ctwasInputsExample$LD_map$LD_file)),
    all(grepl("^pecotmr://", ctwasEstExample$LD_map$LD_file)),
    all(grepl("^pecotmr://", ctwasFinemapExample$LD_map$LD_file))
)

# -----------------------------------------------------------------------------
# 4. Verify each payload still drives the step that consumes it, from the
#    portable form -- the check the previous objects would have failed.
# -----------------------------------------------------------------------------
invisible(estCtwasParam(
    ctwasInputsExample,
    thin = 1,
    niterPrefit = 3,
    niter = 10,
    min_group_size = 1,
    min_p_single_effect = 0,
    fallbackToPrefit = TRUE
))
invisible(finemapCtwasRegions(
    screenCtwasRegions(ctwasEstExample, min_nonSNP_PIP = 0)
))
invisible(asCtwasResult(ctwasFinemapExample))
invisible(mergeCtwasBoundaryRegions(ctwasFinemapExample))

# -----------------------------------------------------------------------------
# 5. Save.
# -----------------------------------------------------------------------------
usethis::use_data(ctwasInputsExample, overwrite = TRUE, compress = "xz")
usethis::use_data(ctwasEstExample, overwrite = TRUE, compress = "xz")
usethis::use_data(ctwasFinemapExample, overwrite = TRUE, compress = "xz")
