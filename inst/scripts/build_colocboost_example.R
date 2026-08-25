#!/usr/bin/env Rscript
#
# build_colocboost_example.R
#
# Builds data/colocboostResultExample.rda.
#
# ColocBoostResult() takes `colocboost` objects, which are far too large and
# too intricate to hand-write in an @examples block -- the minimal stand-in
# used by the unit tests runs to ~50 lines, and hand-rolling one would teach a
# shape nobody should construct themselves. So the fixture comes from a real
# ColocBoost run.
#
# Source: the xqtl-protocol MWE, which is NOT part of this repository:
#   <MWE>/output/colocboost_sim/colocboost/
#     colocboost_sim.chr22_GENESIM1.cb_xqtl.rds
#
# That run is the only one in the MWE or the protocol test fixtures with a
# non-empty cos_summary -- the others found no colocalization, which would
# make for a fixture that cannot demonstrate anything. It carries one
# confidence set colocalizing two outcomes, which is exactly the property
# that distinguishes ColocBoostResult from the pairwise ColocResult.
#
# Usage
#   GITHUB_WORKSPACE=$PWD pixi run --environment r45 Rscript \
#     inst/scripts/build_colocboost_example.R [<path to MWE root>]

suppressMessages(devtools::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
mweRoot <- if (length(args) > 0L) {
    args[[1L]]
} else {
    "~/Downloads/fungen_xqtl/xqtl-protocol"
}

src <- path.expand(file.path(
    mweRoot,
    "output/colocboost_sim/colocboost",
    "colocboost_sim.chr22_GENESIM1.cb_xqtl.rds"
))

if (!file.exists(src)) {
    stop(
        "source not found: ",
        src,
        "\nPass the xqtl-protocol MWE root as the first argument."
    )
}

raw <- readRDS(src)
cb <- if (methods::is(raw, "colocboost")) raw else raw[[1L]]
if (NROW(cb$cos_summary) == 0L) {
    stop("source run found no colocalization; fixture would be empty")
}

# Map the two outcome names onto the identity tuple, as the pipeline does.
outcomeInfo <- data.frame(
    name = c("eQTL_GENESIM1", "psiQTL_GENESIM1"),
    study = "protocol_example",
    context = c("eQTL", "psiQTL"),
    trait = "GENESIM1",
    dataForm = "individual",
    stringsAsFactors = FALSE
)

colocboostResultExample <- ColocBoostResult(
    list(cb),
    "xqtl_coloc",
    outcomeInfo = outcomeInfo
)
validObject(colocboostResultExample)

stopifnot(
    nrow(colocboostResultExample) == 1L,
    length(colocboostResultExample$outcomes[[1L]]) == 2L
)

save(
    colocboostResultExample,
    file = file.path("data", "colocboostResultExample.rda"),
    compress = "xz"
)

cat(sprintf(
    "  colocboostResultExample  %d row, %d outcomes, %d KB\n",
    nrow(colocboostResultExample),
    length(colocboostResultExample$outcomes[[1L]]),
    round(file.size("data/colocboostResultExample.rda") / 1024)
))
